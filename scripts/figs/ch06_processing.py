# ch06 §1 pre-mRNA 가공: DNA → pre-mRNA → cap·splicing·3' 끝/poly(A) → 성숙 mRNA  (redraw of figures/ch06_processing.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *
import numpy as np

UTR, CDS = STEPS[1], PRIMARY          # 연파랑 = UTR, 파랑 = CDS
CAP = CAT[1]                          # 5' cap (청록)
PA = CAT[2]                           # poly(A) (초록)
PA_FILL, PA_INK = mix(PA, 0.22), shade(PA, 0.3)
BAR = 0.22

fig, ax = new(8.4, 5.7)

def seg(x0, x1, y, fc, z=3):
    ax.add_patch(mp.Rectangle((x0, y - BAR / 2), x1 - x0, BAR, fc=fc, ec="none", zorder=z))

# exon layout on the gene: (start, end, [(a, b, colour), ...])
EXONS = [
    (2.0, 3.2, [(2.0, 2.6, UTR), (2.6, 3.2, CDS)]),
    (4.2, 5.0, [(4.2, 5.0, CDS)]),
    (6.0, 7.2, [(6.0, 6.5, CDS), (6.5, 7.2, UTR)]),
]

def exons(y):
    for _, _, parts in EXONS:
        for a, b, c in parts:
            seg(a, b, y, c)

def cap(x, y):
    ax.add_patch(mp.Circle((x, y), 0.1, fc=CAP, ec="none", zorder=4))

# ── row labels ─────────────────────────────────────────────
def rowlab(y, main, sub=None):
    label(ax, 0.1, y + (0.1 if sub else 0), main, fs=12, w="bold", ha="left")
    if sub:
        label(ax, 0.1, y - 0.17, sub, fs=11, c=TEXT2, ha="left")

Y_DNA, Y_PRE, Y_PROC, Y_MAT = 5.0, 3.6, 1.9, 0.75
rowlab(Y_DNA, "유전자 (DNA)")
rowlab(Y_PRE, "pre-mRNA", "막 전사된 RNA")
rowlab(Y_PROC, "가공 단계")
rowlab(Y_MAT, "성숙 mRNA", "세포질로 나감")

# ── ① gene (DNA) ───────────────────────────────────────────
ax.plot([1.75, 7.6], [Y_DNA, Y_DNA], color="#9AA1AB", lw=3, zorder=1, solid_capstyle="butt")
exons(Y_DNA)
ax.plot([2.0, 2.0, 2.22], [Y_DNA + BAR / 2, Y_DNA + 0.42, Y_DNA + 0.42], color=TEXT, lw=1.5, zorder=4)
ax.annotate("", xy=(2.3, Y_DNA + 0.42), xytext=(2.15, Y_DNA + 0.42),
            arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.5, mutation_scale=12))
label(ax, 2.36, Y_DNA + 0.42, "TSS", fs=12, w="bold", ha="left")
for x, t in [(3.7, "intron 1"), (5.5, "intron 2")]:
    label(ax, x, Y_DNA + 0.24, t, fs=11, c=TEXT2)
for (a, b, _), t in zip(EXONS, ["exon 1", "exon 2", "exon 3"]):
    label(ax, (a + b) / 2, Y_DNA - 0.3, t, fs=12, w="bold")

# arrow: 전사
arrow(ax, 4.6, 4.5, 4.6, 4.0)
label(ax, 4.72, 4.25, "전사 (Pol II)", fs=11, c=TEXT2, ha="left")

# ── pre-mRNA ───────────────────────────────────────────────
ax.plot([2.0, 7.5], [Y_PRE, Y_PRE], color="#9AA1AB", lw=1.6, zorder=1)
exons(Y_PRE)
label(ax, 1.85, Y_PRE, "5'", fs=12, w="bold", ha="right")
label(ax, 7.6, Y_PRE, "3'", fs=12, w="bold", ha="left")

# arrow: 가공
arrow(ax, 4.6, 3.35, 4.6, 2.98)
label(ax, 4.72, 3.17, "가공 (핵 안)", fs=11, c=TEXT2, ha="left")

# ── processing ─────────────────────────────────────────────
cap(1.88, Y_PROC)
exons(Y_PROC)
# introns loop out (dashed arcs)
for a, b in [(3.2, 4.2), (5.0, 6.0)]:
    t = np.linspace(0, np.pi, 60)
    xs = (a + b) / 2 - (b - a) / 2 * np.cos(t)
    ys = Y_PROC + BAR / 2 + 0.42 * np.sin(t)
    ax.plot(xs, ys, color=TEXT2, lw=1.5, ls=(0, (3, 2.5)), zorder=2)
    label(ax, a + 0.12, Y_PROC - 0.3, "GU", fs=11, w="bold")
    label(ax, b - 0.12, Y_PROC - 0.3, "AG", fs=11, w="bold")
# 3' end cut + poly(A)
ax.plot([7.2, 7.42], [Y_PROC, Y_PROC], color="#9AA1AB", lw=1.6, zorder=1)
for s in (1, -1):
    ax.plot([7.27, 7.37], [Y_PROC - 0.13 * s, Y_PROC + 0.13 * s], color=TEXT, lw=1.6, zorder=4)
label(ax, 7.32, Y_PROC - 0.3, "3' 끝 자르기", fs=11, c=TEXT2)
box(ax, 7.52, Y_PROC - BAR / 2, 0.62, BAR, fc=PA_FILL, ec=PA, lw=1.0, r=BAR / 2, z=3)
label(ax, 7.83, Y_PROC, "AAAA", fs=11, w="bold", c=PA_INK)

# step markers
Y_STEP = 2.72
for x, n, t in [(1.82, 1, "5' cap"), (4.12, 2, "splicing"), (7.0, 3, "poly(A) 꼬리")]:
    num(ax, x, Y_STEP, n, fc=TEXT, r=0.12, fs=11)
    label(ax, x + 0.18, Y_STEP, t, fs=12, w="bold", ha="left")

# arrow: to mature
arrow(ax, 4.6, 1.42, 4.6, 1.05)

# ── mature mRNA ────────────────────────────────────────────
cap(1.88, Y_MAT)
for a, b, c in [(2.0, 2.6, UTR), (2.6, 3.2, CDS), (3.2, 4.0, CDS), (4.0, 4.5, CDS), (4.5, 5.2, UTR)]:
    seg(a, b, Y_MAT, c)
for x in (3.2, 4.0):     # exon–exon junctions
    ax.plot([x, x], [Y_MAT - BAR / 2, Y_MAT + BAR / 2], color="white", lw=2, zorder=4)
box(ax, 5.25, Y_MAT - BAR / 2, 1.45, BAR, fc=PA_FILL, ec=PA, lw=1.0, r=BAR / 2, z=3)
label(ax, 5.975, Y_MAT, "AAAAAAA … A", fs=11, w="bold", c=PA_INK)

# brackets
Yb = Y_MAT - 0.28
for a, b, t, c, w in [(2.0, 2.58, "5' UTR", TEXT, "bold"),
                      (2.62, 4.48, "단백질 암호 구간 (CDS)", PRIMARY, "bold"),
                      (4.52, 5.2, "3' UTR", TEXT, "bold"),
                      (5.25, 6.7, "poly(A) 꼬리", PA_INK, "bold")]:
    ax.plot([a, a, b, b], [Yb + 0.06, Yb, Yb, Yb + 0.06], color=TEXT2, lw=1.0, zorder=2)
    label(ax, (a + b) / 2, Yb - 0.2, t, fs=11, w=w, c=c)

save(fig, os.path.join(HERE, "../../figures/ch06_processing.png"))
