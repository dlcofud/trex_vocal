"""
음원-필터 모델로 티라노사우루스 발성 합성
보고서 표 1의 범위에서 '중앙값(대표값)'을 뽑아 네 가지 발성 조합을 소리로 만든다.

- 음원: 막을 진동하는 줄로 근사  F0 = 1/(2L(1+eps)) * sqrt(sigma/rho),  sigma = a*exp(b*eps)
- 열린 발성 필터: 한쪽이 막힌 관, 홀수 배 공명  Fn = (2n-1) c / 4(L + 0.61 r)
- 닫힌 발성 필터: 헬름홀츠 공명  f = c/(2pi) sqrt(A / (V * L_eff)), 몸 벽을 통한 저역 통과
- 출력: 물리 스케일 원본 + 사람이 들을 수 있게 모든 주파수를 4배(2옥타브) 올린 청취용
"""
import json
import numpy as np
from scipy.io import wavfile
from scipy.signal import lfilter, butter, sosfilt

SR = 48000
RNG = np.random.default_rng(66)

# ---------- 표 1 대표값 ----------
M = 7750.0                     # kg, 6000-9500 균등의 중앙
rho_tissue = 1040.0            # kg/m^3
a = 2.4e3 / np.sqrt(1 + (2.9 / 2.4) ** 2)   # 로그정규 중앙값 (평균 2.4, SD 2.9 kPa)
b = 16.2
eps_min, eps_max = 0.03, 0.25  # 최대 늘임 15-35%의 중앙
L0 = 3.28 * M ** 0.4 / 1000    # m, 보정 배수 1.0
c = 353.5                      # m/s, 351-356

def f0(eps):
    sigma = a * np.exp(b * eps)          # trex_vocal_sim.ipynb 의 f0()와 같은 식
    return 1 / (2 * L0 * (1 + eps)) * np.sqrt(sigma / rho_tissue)

F0_MIN, F0_MAX = f0(eps_min), f0(eps_max)

# 통로
L_larynx = 1.35; r_tube = 0.055
L_syrinx = L_larynx + 1.2 + 0.3          # + 목 기관 + 가슴 기관
def tube_f1(L): return c / (4 * (L + 0.61 * r_tube))
F1_LARYNX, F1_SYRINX = tube_f1(L_larynx), tube_f1(L_syrinx)

# 헬름홀츠 주머니 (몸 부피의 0.1-1%, 기하 중앙 ≈ 0.316%)
V_body = M / 1000.0
V = V_body * 10 ** -2.5
r_neck, l_neck = 0.04, 0.125
A = np.pi * r_neck ** 2
L_eff = l_neck + 1.7 * r_neck
F_HELM = c / (2 * np.pi) * np.sqrt(A / (V * L_eff))

# ---------- 음원 ----------
def f0_contour(dur, kind):
    t = np.linspace(0, 1, int(dur * SR))
    lo, hi = F0_MIN, F0_MAX
    if kind == "open":   # 올라갔다 내려오는 울음: 하한 근처 → 범위 중간 → 하한
        peak = np.sqrt(lo * hi) * 1.6
        shape = np.sin(np.pi * t ** 0.7) ** 1.5
        f = lo * 1.3 + (peak - lo * 1.3) * shape
    else:                # 닫힌 발성: 낮게 깔리며 천천히 처지는 붐
        f = (lo * 2.2) * (1 - 0.25 * t)
    # 지터: 느린 무작위 흔들림 (~2%)
    wob = lfilter([1], [1, -0.9995], RNG.normal(0, 1, len(t)))
    wob = wob / (np.abs(wob).max() + 1e-9) * 0.02
    return f * (1 + wob)

def glottal_source(f, oq=0.55):
    """Rosenberg 펄스열. 순시 주파수 f[n]을 따라 위상 누적."""
    phase = np.cumsum(f) / SR
    ph = phase % 1.0
    tp = oq * 0.65; tn = oq - tp
    g = np.where(ph < tp, 0.5 * (1 - np.cos(np.pi * ph / tp)),
         np.where(ph < oq, np.cos(np.pi / 2 * (ph - tp) / tn), 0.0))
    # 쉬머: 주기마다 진폭 ±8%
    cyc = np.floor(phase).astype(int)
    amp = 1 + 0.08 * RNG.normal(0, 1, cyc.max() + 1)
    g = g * amp[cyc]
    dg = np.diff(g, prepend=0)                 # 흐름 미분 = 음원 파형
    # 열린 구간에 섞이는 기류 잡음
    noise = RNG.normal(0, 1, len(g)) * g * 0.004
    return dg / np.abs(dg).max() + noise

# ---------- 필터 ----------
def resonator(x, fc, bw):
    r = np.exp(-np.pi * bw / SR)
    th = 2 * np.pi * fc / SR
    a1, a2 = -2 * r * np.cos(th), r * r
    gain = 1 + a1 + a2
    return lfilter([gain], [1, a1, a2], x)

def bp_resonator(x, fc, bw):
    """피크 이득 1의 대역통과 공명기 (병렬 합성용)"""
    r = np.exp(-np.pi * bw / SR)
    th = 2 * np.pi * fc / SR
    return lfilter([(1 - r * r) / 2, 0, -(1 - r * r) / 2], [1, -2 * r * np.cos(th), r * r], x)

def open_tract(x, F1, n_res=12):
    """한쪽 막힌 관: F1, 3F1, 5F1 ... 병렬 공명 + 입 방사(+6dB/oct)"""
    y = np.zeros_like(x)
    for n in range(1, n_res + 1):
        fc = (2 * n - 1) * F1
        if fc > SR / 2.5: break
        bw = 0.12 * fc + 8          # 벽 손실로 대역폭 증가
        y += bp_resonator(x, fc, bw) / n ** 0.6
    return np.diff(y, prepend=0) * SR / (2 * np.pi * F1)  # 열린 입의 방사

def closed_tract(x, Fh):
    """헬름홀츠 단일 공명 + 몸통·피부를 통한 방사(저역 통과)"""
    y = resonator(x, Fh, 0.25 * Fh)
    sos = butter(4, min(Fh * 4, SR / 2.5), btype="low", fs=SR, output="sos")
    return sosfilt(sos, y)

# ---------- 합성 ----------
def envelope(n, att=0.35, rel=0.9):
    t = np.arange(n) / SR; T = n / SR
    e = np.minimum(1, t / att) * np.minimum(1, (T - t) / rel)
    return np.clip(e, 0, 1) ** 1.5

def synth(source_kind, mouth, k=1.0, dur=4.5):
    """k: 주파수 배율 (1=물리값, 4=청취용 2옥타브 위)"""
    f = f0_contour(dur, mouth) * k
    x = glottal_source(f)
    if mouth == "open":
        F1 = (F1_LARYNX if source_kind == "larynx" else F1_SYRINX) * k
        y = open_tract(x, F1)
    else:
        y = closed_tract(x, F_HELM * k)
    sos = butter(2, 2.0, btype="high", fs=SR, output="sos")  # DC 제거
    y = sosfilt(sos, y) * envelope(len(y))
    pad = np.zeros(int(0.3 * SR))
    y = np.concatenate([pad, y, pad])
    return y / np.abs(y).max() * 0.89          # -1 dBFS

if __name__ == "__main__":
    combos = [("larynx", "open", "후두+열린"), ("larynx", "closed", "후두+닫힌"),
              ("syrinx", "open", "울대+열린"), ("syrinx", "closed", "울대+닫힌")]
    # 닫힌 발성은 이 모델에서 후두형과 울대형이 같은 소리이므로 한 파일로 저장
    combos = [c for c in combos if not (c[0] == "syrinx" and c[1] == "closed")]
    for src, mouth, label in combos:
        for k, tag in [(1, "original"), (4, "listen_x4")]:
            y = synth(src, mouth, k)
            fn = f"trex_closed_{tag}.wav" if mouth == "closed" else f"trex_{src}_{mouth}_{tag}.wav"
            wavfile.write(fn, SR, (y * 32767).astype(np.int16))
    params = dict(F0_min=F0_MIN, F0_max=F0_MAX, F1_larynx_open=F1_LARYNX,
                  F1_syrinx_open=F1_SYRINX, F_helmholtz_closed=F_HELM,
                  a_median_Pa=a, L0_mm=L0 * 1000, V_sac_L=V * 1000)
    print(json.dumps({k: round(v, 2) for k, v in params.items()}, indent=1))
