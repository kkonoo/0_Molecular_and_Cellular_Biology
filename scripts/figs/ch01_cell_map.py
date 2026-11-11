# ch01 §2 진핵세포 모식도: 소기관, 구획별 RNA, 세포막 확대  (redraw of figures/ch01_cell_map.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.0, 4.8)

PRE = PRIMARY          # pre-mRNA exon (핵)
MAT = CAT[1]           # 성숙 mRNA (세포질)
MTR = CAT[2]           # 미토콘드리아 / mtDNA 전사체
DNA_C = TEXT2          # DNA (핵 DNA, mtDNA)
INTRON = "#B0B6BF"
RIBO = "#4B5563"
ER_C = "#C5CBD3"
GOLGI_C = "#AEB5BF"

# ── cell ────────────────────────────────────────────────────────────────────
CX, CY, CW, CH = 2.55, 2.4, 5.0, 4.5
ax.add_patch(mp.Ellipse((CX, CY), CW, CH, fc=FAINT, ec=TEXT2, lw=2, zorder=1))

# ── nucleus with double envelope and pores ─────────────────────────────────
NX, NY, NR = 2.0, 2.45, 0.95
ax.add_patch(mp.Circle((NX, NY), NR, fc=STEPS[3], ec="none", zorder=2))
pores = [10, 55, 110, 165, 215, 285, 330]
gap = 5
for i, a in enumerate(pores):
    a2 = pores[(i + 1) % len(pores)] + (360 if i == len(pores) - 1 else 0)
    for r in (NR, NR - 0.06):
        ax.add_patch(mp.Arc((NX, NY), 2 * r, 2 * r, theta1=a + gap, theta2=a2 - gap,
                            color=STEPS[0], lw=1.8, zorder=3))

# DNA squiggle
xs = np.linspace(1.25, 2.75, 200)
ys = 2.66 + 0.1 * np.sin(xs * 9) + 0.05 * np.sin(xs * 23)
ax.plot(xs, ys, color=DNA_C, lw=1.4, zorder=4, solid_capstyle="round")
label(ax, 2.6, 2.98, "DNA", fs=11, c=TEXT2, ha="right")
label(ax, NX - 0.25, 3.12, "핵", fs=13, w="bold")

def pre_mrna(x0, x1, y, n=4, ew=0.16):
    ax.plot([x0, x1], [y, y], color=INTRON, lw=1.5, zorder=4)
    for cx in np.linspace(x0 + ew / 2, x1 - ew / 2, n):
        box(ax, cx - ew / 2, y - 0.045, ew, 0.09, fc=PRE, ec="none", r=0.03, z=5)

pre_mrna(1.35, 2.65, 2.3)
pre_mrna(1.4, 2.6, 2.02)
label(ax, NX, 1.77, "pre-mRNA", fs=11, w="bold", c=PRE)

# 핵공 label (left pore at 165°)
px, py = NX + NR * np.cos(np.radians(165)), NY + NR * np.sin(np.radians(165))
ax.plot([0.8, px - 0.02], [py, py], color=TEXT2, lw=1, zorder=4)
label(ax, 0.76, py, "핵공", fs=12, w="bold", ha="right")

# ── mature mRNA ─────────────────────────────────────────────────────────────
def mrna(x, y, w=0.42):
    box(ax, x - w / 2, y - 0.04, w, 0.08, fc=MAT, ec="none", r=0.04, z=4)
    for dx in (-0.09, 0.09):
        ax.add_patch(mp.Circle((x + dx, y + 0.085), 0.045, fc=RIBO, ec="none", zorder=5))

mrna(3.25, 3.9)
mrna(0.95, 1.25)
mrna(2.05, 0.95)
label(ax, 1.55, 0.68, "세포질", fs=12, w="bold", c=TEXT2)

# export arrow through pore at 55°
ex, ey = NX + NR * np.cos(np.radians(55)), NY + NR * np.sin(np.radians(55))
arrow(ax, ex + 0.03, ey + 0.04, 2.98, 3.8, color=TEXT2, lw=1.4, ms=11)
label(ax, 2.98, 3.52, "핵공으로 나감", fs=11, c=TEXT2, ha="left")

# ── rough ER ────────────────────────────────────────────────────────────────
for x0 in (3.3, 3.57, 3.84):
    yy = np.linspace(1.45, 3.25, 120)
    xx = x0 + 0.05 * np.sin(yy * 7 + x0 * 3)
    ax.plot(xx, yy, color=ER_C, lw=5, zorder=3, solid_capstyle="round")
    for k, yr in enumerate(np.arange(1.55, 3.2, 0.22)):
        xr = x0 + 0.05 * np.sin(yr * 7 + x0 * 3)
        ax.add_patch(mp.Circle((xr + (0.075 if k % 2 == 0 else -0.075), yr), 0.025,
                               fc=RIBO, ec="none", zorder=4))
label(ax, 4.0, 2.75, "거친 ER", fs=12, w="bold", ha="left")
label(ax, 4.0, 2.5, "리보솜이 붙음", fs=11, c=TEXT2, ha="left")

# ── Golgi ───────────────────────────────────────────────────────────────────
for i in range(3):
    r = 0.45 + 0.13 * i
    a = mp.Arc((3.4, 0.35), 2 * r, 2 * r, theta1=55, theta2=125, color=GOLGI_C, lw=4.5, zorder=3)
    a.set_capstyle("round")
    ax.add_patch(a)
label(ax, 3.4, 0.55, "골지체", fs=12, w="bold")

# ── mitochondria ────────────────────────────────────────────────────────────
def mito(cx, cy, ang, L, W, cristae=True):
    ax.add_patch(mp.Ellipse((cx, cy), L, W, angle=ang, fc=mix(MTR, 0.15), ec=MTR, lw=1.6, zorder=3))
    c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
    P = lambda u, v: (cx + u * c - v * s, cy + u * s + v * c)
    if cristae:
        us = np.linspace(-0.36 * L, 0.36 * L, 9)
        pts = [P(u, (0.28 if k % 2 == 0 else -0.28) * W) for k, u in enumerate(us)]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=mix(MTR, 0.45), lw=1, zorder=4)
    rx, ry = P(-0.2 * L, 0)
    ax.add_patch(mp.Circle((rx, ry), 0.05, fc="white", ec=DNA_C, lw=1.3, zorder=5))
    a1, a2 = P(0.0, 0), P(0.24 * L, 0)
    ax.plot([a1[0], a2[0]], [a1[1], a2[1]], color=MTR, lw=4, zorder=5, solid_capstyle="round")

mito(1.05, 3.8, 20, 0.85, 0.4)
mito(4.42, 1.7, 65, 0.75, 0.36)
label(ax, 0.86, 3.3, "미토콘드리아", fs=12, w="bold")

# ── lysosome ────────────────────────────────────────────────────────────────
ax.add_patch(mp.Circle((1.95, 4.1), 0.16, fc=mix(TEXT2, 0.12), ec=TEXT2, lw=1.4, zorder=3))
for dx, dy in [(-0.05, 0.04), (0.05, 0.02), (0.0, -0.06)]:
    ax.add_patch(mp.Circle((1.95 + dx, 4.1 + dy), 0.022, fc=TEXT2, ec="none", zorder=4))
label(ax, 2.17, 4.1, "리소좀", fs=12, w="bold", ha="left")

# 세포막 label
ax.plot([4.6, 4.37], [0.5, 0.81], color=TEXT2, lw=1, zorder=4)
label(ax, 4.62, 0.42, "세포막", fs=12, w="bold", ha="left")

# ── ① 구획별 RNA ────────────────────────────────────────────────────────────
PX, PW = 5.35, 2.6
box(ax, PX, 2.15, PW, 2.45, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, PX + 0.15, 4.35, "① 구획별 RNA", fs=12, w="bold", ha="left")

ix0, ix1, tx = 5.5, 6.02, 6.2
# 핵
y = 3.85
ax.plot([ix0, ix1], [y, y], color=INTRON, lw=1.5, zorder=4)
for cx in (5.56, 5.76, 5.96):
    box(ax, cx - 0.06, y - 0.045, 0.12, 0.09, fc=PRE, ec="none", r=0.03, z=5)
label(ax, tx, y + 0.12, "핵", fs=12, w="bold", ha="left")
label(ax, tx, y - 0.13, "pre-mRNA (intron 포함)", fs=11, c=TEXT2, ha="left")
# 세포질
y = 3.2
box(ax, ix0, y - 0.04, ix1 - ix0, 0.08, fc=MAT, ec="none", r=0.04, z=4)
for dx in (-0.09, 0.09):
    ax.add_patch(mp.Circle(((ix0 + ix1) / 2 + dx, y + 0.085), 0.045, fc=RIBO, ec="none", zorder=5))
label(ax, tx, y + 0.12, "세포질", fs=12, w="bold", ha="left")
label(ax, tx, y - 0.13, "성숙 mRNA (intron 제거)", fs=11, c=TEXT2, ha="left")
# 미토콘드리아
y = 2.55
mito((ix0 + ix1) / 2, y, 0, 0.54, 0.28, cristae=False)
label(ax, tx, y + 0.12, "미토콘드리아", fs=12, w="bold", ha="left")
label(ax, tx, y - 0.13, "mtDNA 전사체 (MT-)", fs=11, c=TEXT2, ha="left")

# ── ② 세포막 확대 ───────────────────────────────────────────────────────────
box(ax, PX, 0.15, PW, 1.75, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, PX + 0.15, 1.65, "② 세포막 확대", fs=12, w="bold", ha="left")
HEAD, TAIL = "#9AA1AB", "#B0B6BF"
ytop, ybot = 1.32, 0.78
heads_x = list(np.arange(5.55, 6.31, 0.15)) + list(np.arange(7.02, 7.8, 0.15))
for hx in heads_x:
    for yh, sgn in ((ytop, -1), (ybot, 1)):
        for dx in (-0.02, 0.02):
            ax.plot([hx + dx, hx + dx], [yh, yh + sgn * 0.2], color=TAIL, lw=1.1, zorder=3)
        ax.add_patch(mp.Circle((hx, yh), 0.055, fc=HEAD, ec="none", zorder=4))
box(ax, 6.36, 0.62, 0.6, 0.86, fc="white", ec=TEXT2, lw=1.4, r=0.12, z=5)
label(ax, 6.66, 1.05, "막\n단백질", fs=11, w="bold", linespacing=1.15, z=6)
label(ax, PX + PW / 2, 0.4, "인지질 이중층", fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch01_cell_map.png"))
