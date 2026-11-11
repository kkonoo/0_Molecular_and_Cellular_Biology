# ch12 §3-3 paired-end: one fragment read from both ends; bulk RNA-seq vs 10x 3′  (new figure)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(6.9, 3.9)

ADAPTER = "#D5D9DF"
CDNA = mix(CAT[0], 0.28)
BC = mix(CAT[1], 0.32)
READ = "#4B5563"

X_A1, X_I0, X_I1, X_A2 = 1.5, 2.25, 5.85, 6.6
FH = 0.36


def seg(x0, x1, y, fc, text):
    ax.add_patch(mp.Rectangle((x0, y), x1 - x0, FH, fc=fc, ec="white", lw=1.5, zorder=2))
    label(ax, (x0 + x1) / 2, y + FH / 2, text, fs=11, c=TEXT)


def reads(y, r1_end, r2_end):
    ya = y + FH + 0.24
    arrow(ax, X_I0, ya, r1_end, ya, color=READ, lw=2.2, ms=14)
    arrow(ax, X_I1, ya, r2_end, ya, color=READ, lw=2.2, ms=14)
    label(ax, X_I0, ya + 0.24, "Read 1", fs=12, w="bold", ha="left")
    label(ax, X_I1, ya + 0.24, "Read 2", fs=12, w="bold", ha="right")


# ── row 1: bulk RNA-seq ─────────────────────────────────────────────
y1 = 2.75
label(ax, 0.1, y1 + FH / 2, "bulk RNA-seq", fs=12, w="bold", ha="left")
seg(X_A1, X_I0, y1, ADAPTER, "어댑터")
seg(X_I0, X_I1, y1, CDNA, "cDNA")
seg(X_I1, X_A2, y1, ADAPTER, "어댑터")
reads(y1, 3.85, 4.25)
# bracket over the whole insert: both ends known
yb = y1 - 0.14
ax.plot([X_I0, X_I0, X_I1, X_I1], [yb + 0.08, yb, yb, yb + 0.08], color=TEXT2, lw=1.3, zorder=3)
label(ax, (X_I0 + X_I1) / 2, yb - 0.22, "조각의 양 끝 위치를 앎", fs=12, c=TEXT2)

# divider
ax.plot([0.1, X_A2], [2.02, 2.02], color=LINE, lw=1.0, zorder=1)

# ── row 2: 10x 3′ ───────────────────────────────────────────────────
y2 = 0.95
X_BC = 3.6
# "10x" drawn in two pieces so the font does not turn "10x" into "10×"
t10 = label(ax, 0.1, y2 + FH / 2, "10", fs=12, w="bold", ha="left")
fig.canvas.draw()
w10 = t10.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
label(ax, 0.1 + w10 + 0.005, y2 + FH / 2, "x 3′", fs=12, w="bold", ha="left")
seg(X_A1, X_I0, y2, ADAPTER, "어댑터")
seg(X_I0, X_BC, y2, BC, "바코드 + UMI")
seg(X_BC, X_I1, y2, CDNA, "cDNA")
seg(X_I1, X_A2, y2, ADAPTER, "어댑터")
reads(y2, X_BC, 4.25)
label(ax, (X_I0 + X_I1) / 2, y2 - 0.3, "두 read의 역할이 다름", fs=12, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch12_paired_end.png"))
