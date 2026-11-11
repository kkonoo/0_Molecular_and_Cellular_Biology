# ch04 §4 재조합으로 haplotype이 잘게 나뉘어도 가까운 변이는 함께 남고, GWAS에서 LD로 묶인 이웃도 함께 유의하게 나오는 모습  (new figure)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.7, 6.5
fig, ax = new(W, H)

# ───────────────────────── ① 재조합 ─────────────────────────
label(ax, 0.1, 6.2, "① 재조합", fs=12, w="bold", ha="left")

BX0, BX1, BH = 1.55, 6.0, 0.28
def X1(f):
    return BX0 + f * (BX1 - BX0)

rows = [
    (5.55, "1세대", [(0, 0.62, 0), (0.62, 1, 1)]),
    (4.95, "몇 세대 후", [(0, 0.22, 0), (0.22, 0.36, 2), (0.36, 0.62, 0),
                       (0.62, 0.80, 1), (0.80, 1, 3)]),
    (4.35, "여러 세대 후", [(0, 0.08, 2), (0.08, 0.17, 0), (0.17, 0.26, 3), (0.26, 0.34, 1),
                         (0.34, 0.37, 2), (0.37, 0.51, 0), (0.51, 0.60, 1), (0.60, 0.66, 3),
                         (0.66, 0.78, 2), (0.78, 0.86, 0), (0.86, 0.93, 1), (0.93, 1, 3)]),
]

for y, name, segs in rows:
    label(ax, 0.1, y, name, fs=11, c=TEXT2, ha="left")
    for a, b, ci in segs:
        ax.add_patch(mp.Rectangle((X1(a), y - BH / 2), X1(b) - X1(a), BH,
                                  fc=CAT[ci], ec="white", lw=1.5, zorder=2))

label(ax, X1(0.2), 5.55, "조상 1", fs=11, w="bold", c=on(CAT[0]))
label(ax, X1(0.81), 5.55, "조상 2", fs=11, w="bold", c=on(CAT[1]))

arrow(ax, 1.35, 5.62, 1.35, 4.25, color="#9AA1AB", lw=1.4, ms=10)

vf = [0.41, 0.46]
for f in vf:
    ax.plot([X1(f), X1(f)], [5.76, 4.12], ls=(0, (2, 2)), color=TEXT2, lw=0.9, zorder=1)
    for y, _, _ in rows:
        ax.add_patch(mp.Circle((X1(f), y), 0.065, fc="white", ec=TEXT, lw=1.2, zorder=4))
label(ax, (X1(vf[0]) + X1(vf[1])) / 2, 5.92, "가까운 두 변이", fs=11, w="bold")

label(ax, 6.2, 5.55, "감수분열마다\n조각 교환", fs=11, c=TEXT2, ha="left", linespacing=1.3)
label(ax, 6.2, 4.35, "조각은 짧아져도\n두 변이는 함께", fs=11, c=TEXT, ha="left", linespacing=1.3)

# ───────────────────────── ② GWAS에서의 뜻 ─────────────────────────
label(ax, 0.1, 3.55, "② GWAS에서의 뜻", fs=12, w="bold", ha="left")

PX0, PX1, PY0, PH, UMAX = 1.15, 6.25, 0.55, 2.55, 12.0
def X2(f):
    return PX0 + f * (PX1 - PX0)
def Y2(u):
    return PY0 + u / UMAX * PH

ax.plot([PX0 - 0.1, PX0 - 0.1], [PY0 - 0.05, PY0 + PH + 0.1], color="#9AA1AB", lw=1.0, zorder=1)
ax.plot([PX0 - 0.1, PX1 + 0.1], [PY0 - 0.05, PY0 - 0.05], color="#9AA1AB", lw=1.0, zorder=1)
label(ax, PX0 - 0.35, PY0 + PH / 2, r"$-\log_{10}\,P$", fs=12, c=TEXT, rotation=90)
label(ax, (PX0 + PX1) / 2, PY0 - 0.3, "염색체 위치", fs=12, c=TEXT)

THR = 5.0
ax.plot([PX0 - 0.1, PX1 + 0.1], [Y2(THR), Y2(THR)], ls=(0, (4, 3)), color=TEXT2, lw=1.0, zorder=1)
label(ax, PX1 + 0.1, Y2(THR) + 0.13, "유의 기준", fs=11, c=TEXT2, ha="right")

# illustrative variants
rng = np.random.default_rng(7)
CF, LF = 0.47, 0.53                      # causal, lead positions
B0, B1 = 0.39, 0.59                      # LD block around the causal variant

# linked neighbours: evenly spread positions, height falls off with distance from the causal one
xb = np.linspace(B0, B1, 24) + rng.uniform(-0.003, 0.003, 24)
xb = xb[(np.abs(xb - CF) > 0.015) & (np.abs(xb - LF) > 0.015)]
d = np.abs(xb - CF)
ub = 5.8 + 3.2 * (1 - d / 0.12) + rng.normal(0, 0.35, xb.size)
ub = np.clip(ub, 5.6, 9.0)

# unlinked variants everywhere along the region
xo = np.sort(rng.uniform(0.02, 0.98, 110))
uo = np.clip(rng.exponential(0.6, xo.size), 0.05, 3.2)

grey = mix(TEXT2, 0.45)
ax.scatter(X2(xo), Y2(uo), s=22, c=grey, lw=0, zorder=3)
ax.scatter(X2(xb), Y2(ub), s=26, c=PRIMARY, lw=0, zorder=3)

cx, cy = X2(CF), Y2(9.4)
lx, ly = X2(LF), Y2(11.0)
ax.scatter([lx], [ly], s=110, c=PRIMARY, edgecolors=TEXT, linewidths=1.6, zorder=6)
ax.scatter([cx], [cy], s=140, c=ACCENT, edgecolors=ACCENT_INK, linewidths=1.3, zorder=6)

ax.plot([cx - 0.08, cx - 0.32], [cy + 0.06, cy + 0.36], color=ACCENT_INK, lw=1.0, zorder=5)
label(ax, cx - 0.36, cy + 0.42, "원인 변이", fs=12, w="bold", ha="right")
ax.plot([lx + 0.1, lx + 0.3], [ly + 0.03, ly + 0.1], color=TEXT2, lw=1.0, zorder=5)
label(ax, lx + 0.34, ly + 0.12, "lead SNP (P 최소)", fs=12, w="bold", ha="left")

# neighbours label → leftmost point of the tower
i0 = np.argmin(xb)
tx, ty = X2(xb[i0]), Y2(ub[i0])
lab_y = Y2(8.0)
label(ax, 1.25, lab_y, "LD로 묶인 이웃 변이", fs=11, c=PRIMARY, w="bold", ha="left")
ax.plot([2.78, tx - 0.07], [lab_y - 0.04, ty + 0.02], color=PRIMARY, lw=0.9, zorder=5)

label(ax, X2(0.67), Y2(7.2), "후보 좁히기\n= fine-mapping", fs=11, c=TEXT2, ha="left",
      linespacing=1.3)

label(ax, X2(0.12), Y2(3.8), "묶이지 않은 변이", fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch04_ld.png"))
