# ch07 §데이터에서는 이렇게 보여요 mRNA가 잠깐 4배가 될 때 반감기 4 h·46 h 단백질의 시간 지연  (redraw of figures/ch07_lag.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

# ── chapter's R model, exactly ───────────────────────────────────────────────
t = np.arange(0, 121)                                  # 0..120 h
m = np.where((t >= 12) & (t < 36), 4.0, 1.0)           # mRNA: 12-36 h만 4배 (가상)

def protein(half_life):
    kd = np.log(2) / half_life
    p = np.ones(len(t))
    for i in range(1, len(t)):
        p[i] = p[i - 1] + kd * (m[i - 1] - p[i - 1])
    return p

p4, p46 = protein(4), protein(46)

assert round(p46[20], 2) == 1.34, p46[20]
assert round(p46[50], 2) == 1.74, p46[50]
assert round(p46[120], 2) == 1.26, p46[120]
assert round(p4[20], 2) == 3.35, p4[20]
assert round(p4[50], 2) == 1.21, p4[50]

# ── colours: three parallel series ───────────────────────────────────────────
C_M, C_4, C_46 = CAT[3], CAT[1], CAT[0]

fig, ax = plt.subplots(figsize=(7.6, 3.9))
ax.set_xlim(0, 120)
ax.set_ylim(0, 4.6)
for y in range(1, 5):
    ax.axhline(y, color=LINE, lw=0.8, zorder=0)

ax.step(t, m, where="post", color=C_M, lw=2.4, zorder=3)
ax.plot(t, p4, color=C_4, lw=2.4, zorder=4)
ax.plot(t, p46, color=C_46, lw=2.6, zorder=5)

# time points
for x in (20, 50):
    ax.axvline(x, color=TEXT2, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.text(x, 4.5, f"{x} h", ha="center", va="center", fontsize=12,
            fontweight="bold", color=TEXT)

def dot(x, y, c):
    ax.plot([x], [y], "o", ms=7, mfc=c, mec="white", mew=1.4, zorder=7)

# 20 h values
dot(20, m[20], C_M);   ax.text(21.5, 3.84, f"{m[20]:.0f}", ha="left", va="center", fontsize=11, color=TEXT)
dot(20, p4[20], C_4);  ax.text(21.5, 3.08, f"{p4[20]:.2f}", ha="left", va="center", fontsize=11, color=TEXT)
dot(20, p46[20], C_46); ax.text(21.5, 1.13, f"{p46[20]:.2f}", ha="left", va="center", fontsize=11, color=TEXT)

# 50 h values
dot(50, m[50], C_M);   ax.text(51.2, 0.8, f"{m[50]:.0f}", ha="left", va="center", fontsize=11, color=TEXT)
dot(50, p4[50], C_4);  ax.text(51.2, 1.38, f"{p4[50]:.2f}", ha="left", va="center", fontsize=11, color=TEXT)
dot(50, p46[50], C_46)
ax.annotate(f"{p46[50]:.2f}배 (mRNA는 1배)", xy=(50, p46[50]), xytext=(53, 2.25),
            fontsize=12, fontweight="bold", color=ON_ACCENT, ha="left", va="center",
            bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.3", fc=ACCENT, ec="none"),
            arrowprops=dict(arrowstyle="-", color=ACCENT_INK, lw=1.0,
                            shrinkA=0, shrinkB=4), zorder=8)

# direct series labels (kept clear of the dashed time lines)
ax.text(29, 4.22, "mRNA", ha="center", va="center", fontsize=12, fontweight="bold", color=C_M)
ax.text(40, 3.05, "단백질 (반감기 4 h)", ha="left", va="center", fontsize=12, bbox=dict(fc="white", ec="none", pad=1.5), zorder=6,
        fontweight="bold", color=C_4)
ax.text(119, 1.62, "단백질 (반감기 46 h)", ha="right", va="center", fontsize=12,
        fontweight="bold", color=C_46)

ax.set_xticks(range(0, 121, 24))
ax.set_yticks(range(0, 5))
ax.tick_params(labelsize=11)
ax.set_xlabel("시간 (h)", fontsize=12)
ax.set_ylabel("평소 대비 배수", fontsize=12)

save(fig, os.path.join(HERE, "../../figures/ch07_lag.png"))
