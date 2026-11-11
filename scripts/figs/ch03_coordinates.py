# ch03 §coordinates 같은 염기 10개 위의 1-based(글자 번호)와 0-based(경계 번호) 좌표계  (redraw of figures/ch03_coordinates.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(7.2, 3.4)

seq = "ACGTTGCAAG"
x0, cw = 1.85, 0.5          # left edge of first cell, cell width
y0, ch = 1.6, 0.52          # cell row bottom, height
C1 = CAT[0]                 # 1-based (blue)
C0 = CAT[1]                 # 0-based (teal)
EDGE = "#C5CAD2"

# ── base cells ──────────────────────────────────────────────
for i, b in enumerate(seq):
    hit = 2 <= i <= 5       # GTTG = letters 3..6
    fc = ACCENT if hit else "white"
    ax.add_patch(mp.Rectangle((x0 + i * cw, y0), cw, ch, fc=fc, ec=EDGE, lw=1.0, zorder=2))
    label(ax, x0 + i * cw + cw / 2, y0 + ch / 2, b, fs=15, w="bold",
          c=ON_ACCENT if hit else TEXT)

# ── 1-based: numbers on letters (top) ───────────────────────
ytop = y0 + ch + 0.2
for i in range(10):
    label(ax, x0 + i * cw + cw / 2, ytop, str(i + 1), fs=12, w="bold", c=C1)

# ── 0-based: numbers on boundaries (bottom) ─────────────────
for i in range(11):
    x = x0 + i * cw
    ax.plot([x, x], [y0 - 0.17, y0], color=C0, lw=1.6, zorder=3, solid_capstyle="butt")
    label(ax, x, y0 - 0.36, str(i), fs=12, w="bold", c=C0)

# ── brackets ────────────────────────────────────────────────
def bracket(xa, xb, y, tip, color):
    ax.plot([xa, xa, xb, xb], [y + tip, y, y, y + tip], color=color, lw=1.5, zorder=3)

xa, xb = x0 + 2 * cw, x0 + 6 * cw
xm = (xa + xb) / 2

yb = ytop + 0.25
bracket(xa + 0.04, xb - 0.04, yb, -0.08, C1)
label(ax, xm, yb + 0.24, "start 3, end 6  →  [3, 6]", fs=12, w="bold", c=C1)
label(ax, xm, yb + 0.52, "길이 6 − 3 + 1 = 4", fs=11, c=TEXT2)

yb2 = y0 - 0.62
bracket(xa, xb, yb2, 0.08, C0)
label(ax, xm, yb2 - 0.24, "start 2, end 6  →  [2, 6)", fs=12, w="bold", c=C0)
label(ax, xm, yb2 - 0.52, "길이 6 − 2 = 4", fs=11, c=TEXT2)

# ── row labels (left) ───────────────────────────────────────
xl = x0 - 0.25
label(ax, xl, ytop + 0.1, "1-based", fs=12, w="bold", c=C1, ha="right")
label(ax, xl, ytop - 0.15, "GTF·VCF·SAM", fs=11, c=TEXT2, ha="right")
label(ax, xl, y0 - 0.26, "0-based", fs=12, w="bold", c=C0, ha="right")
label(ax, xl, y0 - 0.51, "BED", fs=11, c=TEXT2, ha="right")

save(fig, os.path.join(HERE, "../../figures/ch03_coordinates.png"))
