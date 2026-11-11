# ch06 §4 isoform A(exon 1-2-3-4)와 B(exon 1-2-4), read가 어느 isoform에 맞는지  (redraw of figures/ch06_isoforms.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.7, 5.4
fig, ax = new(W, H)

COL_A = CAT[0]          # isoform A, A에서만 가능한 read
COL_B = CAT[1]          # isoform B, B에서만 가능한 read
COL_AMB = "#A3A9B2"     # 둘 다 가능한 read
DNA_FILL = "#C9CDD3"

EX = {1: (1.75, 2.55), 2: (3.20, 3.90), 3: (4.60, 5.30), 4: (6.00, 7.35)}
LX = 0.05               # left label column

# ── 유전자 (DNA) ──────────────────────────────────────────
y_dna, bh = 4.65, 0.24
ax.plot([1.55, 7.55], [y_dna, y_dna], color="#9AA1AB", lw=1.6, zorder=1)
for k, (a, b) in EX.items():
    box(ax, a, y_dna - bh / 2, b - a, bh, fc=DNA_FILL, ec="none", r=0.03, z=2)
    label(ax, (a + b) / 2, y_dna + 0.32, f"exon {k}", fs=11, c=TEXT2, w="bold")
label(ax, LX, y_dna, "유전자 (DNA)", fs=12, w="bold", ha="left")

# ── isoforms ──────────────────────────────────────────────
def isoform(y, exons, col):
    for k in exons:
        a, b = EX[k]
        box(ax, a, y - bh / 2, b - a, bh, fc=col, ec="none", r=0.03, z=2)
    for k1, k2 in zip(exons[:-1], exons[1:]):
        x1, x2 = EX[k1][1], EX[k2][0]
        yt = y + bh / 2
        ax.plot([x1, (x1 + x2) / 2, x2], [yt, yt + 0.16, yt], color=col, lw=1.2, zorder=1)

y_a, y_b = 3.85, 3.15
isoform(y_a, [1, 2, 3, 4], COL_A)
isoform(y_b, [1, 2, 4], COL_B)
label(ax, LX, y_a + 0.1, "isoform A", fs=12, w="bold", c=COL_A, ha="left")
label(ax, LX, y_a - 0.17, "exon 3 포함", fs=11, c=TEXT2, ha="left")
label(ax, LX, y_b + 0.1, "isoform B", fs=12, w="bold", c=COL_B, ha="left")
label(ax, LX, y_b - 0.17, "exon 3 건너뜀", fs=11, c=TEXT2, ha="left")

ax.plot([LX, 7.55], [2.72, 2.72], color=LINE, lw=1.0, zorder=1)

# ── reads ─────────────────────────────────────────────────
label(ax, LX, 2.47, "read", fs=12, w="bold", ha="left")
label(ax, LX + 0.47, 2.465, "(정렬 위치)", fs=11, c=TEXT2, ha="left")

y_top, y_bot = 2.32, 0.42
for k, (a, b) in EX.items():          # faint exon columns for alignment
    box(ax, a, y_bot, b - a, y_top - y_bot, fc=FAINT, ec="none", r=0.02, z=0)

rh = 0.09
def read(y, segs, col):
    for a, b in segs:
        box(ax, a, y - rh / 2, b - a, rh, fc=col, ec="none", r=0.02, z=3)
    for (a1, b1), (a2, b2) in zip(segs[:-1], segs[1:]):
        ax.plot([b1, a2], [y, y], color=col, lw=1.0, zorder=2)

# 둘 다 가능 (exon 1, 2, 4 안)
for y, segs_list in [(2.12, [[(1.85, 2.30)], [(3.28, 3.72)], [(6.10, 6.55)], [(6.85, 7.28)]]),
                     (1.92, [[(2.05, 2.48)], [(3.45, 3.86)], [(6.50, 6.95)]])]:
    for s in segs_list:
        read(y, s, COL_AMB)
# A에서만 (exon 3, exon 2-3 · 3-4 경계)
read(1.55, [(4.75, 5.20)], COL_A)
read(1.35, [(3.55, 3.90), (4.60, 4.85)], COL_A)
read(1.15, [(5.00, 5.30), (6.00, 6.28)], COL_A)
# B에서만 (exon 2-4 경계)
read(0.80, [(3.55, 3.90), (6.00, 6.18)], COL_B)
read(0.60, [(3.75, 3.90), (6.00, 6.35)], COL_B)

label(ax, LX, 2.02 + 0.11, "둘 다 가능", fs=12, w="bold", c=TEXT2, ha="left")
label(ax, LX, 2.02 - 0.16, "→ 모호", fs=11, c=TEXT2, ha="left")
label(ax, LX, 1.35, "A에서만", fs=12, w="bold", c=COL_A, ha="left")
label(ax, LX, 0.70, "B에서만", fs=12, w="bold", c=COL_B, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch06_isoforms.png"))
