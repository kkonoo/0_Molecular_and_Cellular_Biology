# ch05 §3 enhancer가 DNA 고리를 통해 promoter 곁으로 오는 모습  (redraw of figures/ch05_regulation.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W = 8.0
fig, ax = new(W, 5.65)

DNA_C = "#8C939D"
GENE_N = "#B8BEC6"            # 이웃 유전자 (회색)
ENH_C = CAT[1]                # enhancer (teal)
TF2_C = CAT[2]                # 두 번째 TF (green)
PROM_C = PRIMARY              # promoter · 표적 유전자 · Pol II (blue)
GTF_C = STEPS[1]
RNA_C = shade(CAT[3], 0.15)

def dna(x1, x2, y):
    ax.plot([x1, x2], [y, y], color=DNA_C, lw=3, solid_capstyle="round", zorder=2)

def gene_block(x1, x2, y, fc, h=0.24):
    box(ax, x1, y - h / 2, x2 - x1, h, fc=fc, ec="none", r=0.03, z=3)

# ── ① 한 줄로 펼친 DNA ─────────────────────────────────────────────
box(ax, 0, 4.0, W, 1.62, fc=FAINT, ec="none", r=0.1, z=0)
label(ax, 0.18, 5.4, "① 한 줄로 펼친 DNA", fs=12, w="bold", ha="left")

y1 = 4.75
dna(0.3, 7.8, y1)
box(ax, 0.6, y1 - 0.12, 1.0, 0.24, fc=mix(ENH_C, 0.18), ec=ENH_C, lw=1.4, r=0.04, z=3)
label(ax, 1.1, 5.05, "enhancer", fs=12, w="bold", c=ENH_C)
for x in (2.3, 3.0, 3.7):
    gene_block(x, x + 0.3, y1, GENE_N)
label(ax, 3.15, 5.05, "이웃 유전자", fs=11, c=TEXT2)
box(ax, 4.6, y1 - 0.12, 1.0, 0.24, fc=mix(PROM_C, 0.15), ec=PROM_C, lw=1.4, r=0.04, z=3)
label(ax, 5.1, 5.05, "promoter", fs=12, w="bold", c=PROM_C)
ax.plot([5.6, 5.6], [y1 + 0.12, 5.1], color=TEXT, lw=1.6, zorder=4)
arrow(ax, 5.6, 5.1, 5.92, 5.1, color=TEXT, lw=1.6, ms=11, z=4)
label(ax, 5.97, 5.1, "TSS", fs=12, w="bold", ha="left")
for a, b in ((6.25, 6.65), (6.85, 7.05), (7.25, 7.65)):
    gene_block(a, b, y1, PROM_C)
label(ax, 6.95, 5.05, "표적 유전자", fs=12, w="bold", c=PROM_C)
yb = 4.47
ax.plot([1.6, 1.6, 4.6, 4.6], [yb + 0.09, yb, yb, yb + 0.09], color=TEXT2, lw=1.1, zorder=2)
label(ax, 3.1, 4.24, "멀리 떨어져 있을 수 있음", fs=11, c=TEXT2)

# ── ② DNA 고리 ─────────────────────────────────────────────────────
box(ax, 0, 0.6, W, 3.25, fc=FAINT, ec="none", r=0.1, z=0)
label(ax, 0.18, 3.62, "② DNA가 고리를 만들면", fs=12, w="bold", ha="left")

y2 = 2.4
dna(0.3, 2.9, y2)
dna(3.2, 7.8, y2)

def bez(p0, p1, p2, p3, t):
    t = np.asarray(t)[:, None]
    p0, p1, p2, p3 = map(np.array, (p0, p1, p2, p3))
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3

yl = 0.95                                   # 고리 바닥
L = ((2.9, y2), (2.9, 1.65), (2.0, yl), (3.05, yl))
R = ((3.05, yl), (4.1, yl), (3.2, 1.65), (3.2, y2))
tt = np.linspace(0, 1, 120)
for seg in (L, R):
    p = bez(*seg, tt)
    ax.plot(p[:, 0], p[:, 1], color=DNA_C, lw=3, solid_capstyle="round", zorder=2)
for seg, (a, b) in ((L, (0.62, 0.74)), (L, (0.86, 0.98)), (R, (0.04, 0.16)), (R, (0.28, 0.40))):
    p = bez(*seg, np.linspace(a, b, 20))
    ax.plot(p[:, 0], p[:, 1], color=GENE_N, lw=9, solid_capstyle="butt", zorder=3)
label(ax, 2.2, 1.12, "이웃 유전자", fs=11, c=TEXT2, ha="right")
label(ax, 3.95, 1.5, "DNA 고리 (loop)", fs=11, c=TEXT2, ha="left")

# enhancer + TF
box(ax, 1.9, y2 - 0.12, 1.0, 0.24, fc=mix(ENH_C, 0.18), ec=ENH_C, lw=1.4, r=0.04, z=3)
label(ax, 2.3, 2.03, "enhancer", fs=12, w="bold", c=ENH_C)
ax.add_patch(mp.Circle((2.2, 2.74), 0.22, fc=ENH_C, ec="white", lw=1.2, zorder=5))
ax.add_patch(mp.Circle((2.62, 2.74), 0.22, fc=TF2_C, ec="white", lw=1.2, zorder=5))
ax.plot([1.72, 2.0], [3.12, 2.86], color=TEXT2, lw=1, zorder=4)
label(ax, 1.68, 3.15, "특이 전사인자 (TF)", fs=12, w="bold", ha="right")

# promoter + 일반 전사인자 + Pol II
box(ax, 3.2, y2 - 0.12, 1.0, 0.24, fc=mix(PROM_C, 0.15), ec=PROM_C, lw=1.4, r=0.04, z=3)
label(ax, 3.7, 2.03, "promoter", fs=12, w="bold", c=PROM_C)
for x in (3.36, 3.56):
    ax.add_patch(mp.Circle((x, 2.63), 0.1, fc=GTF_C, ec="white", lw=1, zorder=5))
ax.add_patch(mp.Circle((3.95, 2.82), 0.32, fc=PROM_C, ec="white", lw=1.2, zorder=6))
label(ax, 3.95, 2.82, "Pol II", fs=12, w="bold", c="white", z=7)
ax.plot([4.28, 4.64], [2.92, 3.05], color=TEXT2, lw=1, zorder=4)
label(ax, 4.68, 3.1, "RNA 중합효소 II (Pol II)", fs=12, w="bold", c=PROM_C, ha="left")
label(ax, 4.68, 2.84, "+ 일반 전사인자", fs=11, c=TEXT2, ha="left")
ys = np.linspace(3.14, 3.58, 60)
xs = 3.97 + 0.06 * np.sin((ys - 3.14) * 28)
ax.plot(xs, ys, color=RNA_C, lw=2, zorder=6)
label(ax, 4.12, 3.52, "새로 만들어지는 RNA", fs=11, c=RNA_C, ha="left")
ax.plot([4.2, 4.2], [y2 - 0.12, 1.98], color=TEXT, lw=1.6, zorder=4)
arrow(ax, 4.2, 1.98, 4.5, 1.98, color=TEXT, lw=1.6, ms=11, z=4)
label(ax, 4.55, 1.98, "TSS", fs=12, w="bold", ha="left")
for a, b in ((5.2, 5.7), (6.0, 6.25), (6.55, 7.0), (7.25, 7.65)):
    gene_block(a, b, y2, PROM_C)
label(ax, 6.4, 2.03, "표적 유전자", fs=12, w="bold", c=PROM_C)

# 한 가지 핵심: 서로 가까이 모임 (DNA·상자 뒤에 깔리도록 낮은 zorder)
ax.add_patch(mp.Ellipse((3.1, 2.72), 2.9, 1.15, fc="none", ec=ACCENT_INK, lw=1.6,
                        ls=(0, (4, 3)), zorder=1))
label(ax, 2.6, 3.45, "서로 가까이 모임", fs=12, w="bold", c=ACCENT_INK)

save(fig, os.path.join(HERE, "../../figures/ch05_regulation.png"))
