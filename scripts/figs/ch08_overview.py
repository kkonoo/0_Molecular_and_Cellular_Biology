# ch08 §1 신호전달의 일반 구조: 리간드 → 수용체 → kinase cascade → 전사인자 → 표적 유전자, 음성 되먹임  (redraw of figures/ch08_overview.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *
import matplotlib.patches as mp
import numpy as np

W, H = 8.4, 5.25
fig, ax = new(W, H)

ARW = "#6B7280"
MEMB = "#9AA1AB"
FB = CAT[2]                 # 음성 되먹임 / 억제 단백질
GENE = CAT[1]               # 표적 유전자 / mRNA

# ── 세포 (막 = 위쪽 테두리) ─────────────────────────────────────────
box(ax, 0.2, 0.2, 7.95, 4.1, fc="white", ec=MEMB, lw=2.2, r=0.3, z=1)
label(ax, 8.05, 4.95, "세포 밖", c=TEXT2, ha="right")
label(ax, 8.05, 4.45, "세포막", c=TEXT2, ha="right")
label(ax, 7.95, 4.1, "세포 안", c=TEXT2, ha="right")

# ── ① 신호 (리간드) ─────────────────────────────────────────────────
for (x, y) in [(1.7, 4.84), (2.35, 5.0), (2.75, 4.78)]:
    ax.add_patch(mp.Circle((x, y), 0.13, fc=PRIMARY, ec="none", zorder=4))
num(ax, 3.2, 4.9, 1)
label(ax, 3.42, 4.9, "신호 (리간드)", fs=12, w="bold", ha="left")

# ── ② 수용체 ────────────────────────────────────────────────────────
for x0 in (1.52, 1.72):
    box(ax, x0, 3.95, 0.16, 0.75, fc=STEPS[2], ec=PRIMARY, lw=1.2, r=0.07, z=3)
num(ax, 0.55, 4.0, 2)
label(ax, 0.77, 4.0, "수용체", fs=12, w="bold", ha="left")

# 2차 전달자
arrow(ax, 1.95, 4.05, 2.3, 4.05, color=ARW, lw=1.3, ms=10)
for (x, y) in [(2.45, 4.1), (2.62, 3.98), (2.78, 4.1)]:
    ax.add_patch(mp.Circle((x, y), 0.045, fc=CAT[3], ec="none", zorder=4))
label(ax, 2.95, 4.04, "2차 전달자 (cAMP, Ca²⁺ 등)", c=TEXT2, ha="left")

# ── ③ 전달: kinase cascade ──────────────────────────────────────────
arrow(ax, 1.7, 3.92, 1.7, 3.17, color=ARW)
ky, kh, kw = 2.75, 0.4, 1.15
kx = [1.7, 3.4, 5.1]
for i, cx in enumerate(kx):
    box(ax, cx - kw / 2, ky, kw, kh, fc=STEPS[3], ec=PRIMARY, lw=1.2, z=3)
    label(ax, cx, ky + kh / 2, f"kinase {i + 1}", fs=12, w="bold")
for a, b in [(0, 1), (1, 2)]:
    x1, x2 = kx[a] + kw / 2 + 0.03, kx[b] - kw / 2 - 0.03
    arrow(ax, x1, ky + kh / 2, x2, ky + kh / 2, color=ARW)
    label(ax, (x1 + x2) / 2, ky + kh / 2 + 0.18, "+p", c=PRIMARY, w="bold")
num(ax, 0.55, 2.4, 3)
label(ax, 0.77, 2.4, "전달: kinase cascade", fs=12, w="bold", ha="left")

# kinase 3 → 전사인자
tfx, tfy, tfh = 5.1, 1.95, 0.4
arrow(ax, 5.1, ky - 0.03, 5.1, tfy + tfh + 0.03, color=ARW)
label(ax, 5.22, (ky + tfy + tfh) / 2, "+p", c=PRIMARY, w="bold", ha="left")
box(ax, tfx - kw / 2, tfy, kw, tfh, fc=STEPS[3], ec=PRIMARY, lw=1.2, z=3)
label(ax, tfx, tfy + tfh / 2, "전사인자", fs=12, w="bold")

# ── 핵 ──────────────────────────────────────────────────────────────
ax.add_patch(mp.Ellipse((5.6, 1.0), 4.6, 1.4, fc=mix(GENE, 0.07), ec=GENE, lw=1.5, zorder=2))
label(ax, 3.75, 1.0, "핵", fs=12, c=TEXT2)
for y in (0.97, 1.05):
    ax.plot([4.1, 7.3], [y, y], color=MEMB, lw=1.4, zorder=3, solid_capstyle="round")
# 전사인자가 DNA에 붙음
arrow(ax, 4.9, tfy - 0.03, 4.78, 1.26, color=ARW)
ax.add_patch(mp.Ellipse((4.75, 1.12), 0.42, 0.22, fc=PRIMARY, ec="none", zorder=4))
# 표적 유전자
box(ax, 5.2, 0.83, 1.4, 0.36, fc=GENE, ec="none", r=0.05, z=4)
label(ax, 5.9, 1.01, "표적 유전자", fs=12, w="bold", c="white")
num(ax, 4.45, 0.6, 4)
label(ax, 4.67, 0.6, "반응: 표적 유전자 발현 변화", fs=12, w="bold", ha="left")

# ── mRNA → 번역 → 억제 단백질 ──────────────────────────────────────
t = np.linspace(0, 1, 120)
mx = 6.3 + 0.07 * np.sin(t * 4 * np.pi)
my = 1.19 + t * 1.11
ax.plot(mx, my, color=GENE, lw=2.2, zorder=5, solid_capstyle="round")
arrow(ax, 6.3, 2.28, 6.3, 2.46, color=GENE, lw=2.2, ms=12, z=5)
label(ax, 6.5, 2.25, "mRNA", fs=12, w="bold", c=GENE, ha="left")
pill(ax, 6.5, 1.92, "RNA-seq이 세는 것", fc=ACCENT, tc=ON_ACCENT, fs=11, ha="left")

arrow(ax, 6.3, 2.52, 6.3, 3.37, color=ARW, lw=1.3, ls="--", ms=10)
label(ax, 6.42, 2.95, "번역", c=TEXT2, ha="left")

box(ax, 5.4, 3.4, 2.5, 0.4, fc=mix(FB, 0.12), ec=FB, lw=1.3, z=3)
label(ax, 6.65, 3.6, "억제 단백질 (SOCS, DUSP 등)", fs=12, w="bold")

# 음성 되먹임: 억제 단백질 ⊣ kinase 2
ax.plot([5.4, 3.4, 3.4], [3.6, 3.6, 3.21], color=FB, lw=2.0, zorder=3,
        solid_joinstyle="miter")
ax.plot([3.25, 3.55], [3.21, 3.21], color=FB, lw=3.0, zorder=3, solid_capstyle="butt")
label(ax, 4.4, 3.77, "음성 되먹임", fs=12, w="bold", c=shade(FB, 0.25))

save(fig, os.path.join(HERE, "../../figures/ch08_overview.png"))
