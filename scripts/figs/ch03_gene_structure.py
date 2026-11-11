# ch03 §2 유전자 하나의 구조: 유전체 위 부품 → pre-mRNA → 성숙한 mRNA  (redraw of figures/ch03_gene_structure.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.0, 5.6)

CDS_C = PRIMARY                 # 번역되는 부분
UTR_C = STEPS[2]                # 번역 안 되는 부분 (옅은 파랑)
REG_C = CAT[1]                  # 조절 서열 (청록)
INTRON_C = mix(TEXT2, 0.35)     # RNA 속 intron
BH = 0.34                       # box height


def seg(x0, x1, y, kind):
    if kind == "cds":
        box(ax, x0, y - BH / 2, x1 - x0, BH, fc=CDS_C, ec=CDS_C, lw=1.0, r=0.02, z=3)
    else:
        box(ax, x0, y - BH / 2, x1 - x0, BH, fc=UTR_C, ec=CDS_C, lw=1.0, r=0.02, z=3)


# ── 위: 유전체 위의 유전자 ─────────────────────────────────────────
yg = 4.55
label(ax, 0.15, yg, "5'", fs=11, c=TEXT2)
label(ax, 7.8, yg, "3'", fs=11, c=TEXT2)
ax.plot([0.3, 1.5], [yg, yg], color=TEXT2, lw=1.4, zorder=1)
ax.plot([1.72, 7.65], [yg, yg], color=TEXT2, lw=1.4, zorder=1)
for xb in (1.54, 1.66):                      # 멀리 떨어짐 표시 //
    ax.plot([xb - 0.05, xb + 0.05], [yg - 0.13, yg + 0.13], color=TEXT2, lw=1.4, zorder=1)

# 조절 서열
box(ax, 0.35, yg - BH / 2, 1.0, BH, fc=REG_C, ec=REG_C, r=0.06, z=3)
label(ax, 0.85, yg, "enhancer", fs=11, w="bold", c="white")
label(ax, 0.85, yg - 0.42, "멀리 있을 수 있음", fs=11, c=TEXT2)
box(ax, 1.85, yg - BH / 2, 1.0, BH, fc=REG_C, ec=REG_C, r=0.06, z=3)
label(ax, 2.35, yg, "promoter", fs=11, w="bold", c="white")

# 전사되는 부분
x_tss, x_start, x_e1, x_e2a, x_e2b, x_e3a, x_stop, x_end = 2.85, 3.25, 4.0, 4.8, 5.55, 6.25, 6.85, 7.45
seg(x_tss, x_start, yg, "utr")
seg(x_start, x_e1, yg, "cds")
seg(x_e2a, x_e2b, yg, "cds")
seg(x_e3a, x_stop, yg, "cds")
seg(x_stop, x_end, yg, "utr")

# TSS
ax.plot([x_tss, x_tss], [yg + BH / 2, yg + 0.62], color=TEXT, lw=1.6, zorder=4)
arrow(ax, x_tss - 0.008, yg + 0.62, 3.35, yg + 0.62, color=TEXT, lw=1.6, ms=12, z=4)
label(ax, 3.43, yg + 0.62, "TSS", fs=12, w="bold", ha="left")

# polyA 신호
ax.add_patch(mp.Polygon([(7.16, yg + 0.33), (7.34, yg + 0.33), (7.25, yg + BH / 2 + 0.01)],
                        closed=True, fc=TEXT, ec="none", zorder=5))
label(ax, 7.25, yg + 0.5, "polyA 신호", fs=11, w="bold")

# start / stop codon
for xc, t in ((x_start, "start codon"), (x_stop, "stop codon")):
    ax.plot([xc, xc], [yg - BH / 2, yg - 0.3], color=TEXT, lw=1.0, zorder=4)
    label(ax, xc, yg - 0.42, t, fs=11, c=TEXT)

# exon / intron 이름 (exon 1·3은 UTR까지 포함한 전체 구간 가운데)
ye = yg - 0.75
label(ax, (x_tss + x_e1) / 2, ye, "exon 1", fs=12, w="bold", c=PRIMARY)
label(ax, (x_e2a + x_e2b) / 2, ye, "exon 2", fs=12, w="bold", c=PRIMARY)
label(ax, (x_e3a + x_end) / 2, ye, "exon 3", fs=12, w="bold", c=PRIMARY)
label(ax, (x_e1 + x_e2a) / 2, ye, "intron 1", fs=12, c=TEXT2)
label(ax, (x_e2b + x_e3a) / 2, ye, "intron 2", fs=12, c=TEXT2)

# ── 전사 ───────────────────────────────────────────────────────────
xa = 4.4
arrow(ax, xa, ye - 0.2, xa, 3.22, color="#9AA1AB", lw=1.6)
label(ax, xa + 0.13, (ye - 0.2 + 3.22) / 2, "전사 (5장)", fs=11, ha="left", c=TEXT)

# ── pre-mRNA ───────────────────────────────────────────────────────
yp = 3.0
label(ax, 2.65, yp + 0.08, "pre-mRNA", fs=12, w="bold", ha="right")
label(ax, 2.65, yp - 0.2, "intron 포함", fs=11, c=TEXT2, ha="right")
for x0, x1 in ((x_e1, x_e2a), (x_e2b, x_e3a)):
    box(ax, x0, yp - 0.06, x1 - x0, 0.12, fc=INTRON_C, ec="none", r=0.0, z=2)
seg(x_tss, x_start, yp, "utr")
seg(x_start, x_e1, yp, "cds")
seg(x_e2a, x_e2b, yp, "cds")
seg(x_e3a, x_stop, yp, "cds")
seg(x_stop, x_end, yp, "utr")

# ── RNA 가공 ──────────────────────────────────────────────────────
arrow(ax, xa, yp - 0.25, xa, 2.0, color="#9AA1AB", lw=1.6)
label(ax, xa + 0.13, (yp - 0.25 + 2.0) / 2, "splicing · poly(A) 꼬리 붙이기 (6장)", fs=11, ha="left", c=TEXT)

# ── 성숙한 mRNA ───────────────────────────────────────────────────
ym = 1.75
label(ax, 2.65, ym, "성숙한 mRNA", fs=12, w="bold", ha="right")
m = [2.85, 3.25, 4.0, 4.75, 5.35, 5.95]
seg(m[0], m[1], ym, "utr")
seg(m[1], m[4], ym, "cds")
seg(m[4], m[5], ym, "utr")
for xj in (m[2], m[3]):                      # exon 이음매
    ax.plot([xj, xj], [ym - BH / 2 + 0.02, ym + BH / 2 - 0.02], color="white", lw=2.0, zorder=4)
label(ax, (m[0] + m[1]) / 2, ym + 0.36, "5'UTR", fs=11, c=TEXT2)
label(ax, (m[4] + m[5]) / 2, ym + 0.36, "3'UTR", fs=11, c=TEXT2)
label(ax, m[5] + 0.08, ym, "AAAA…", fs=12, w="bold", c=PRIMARY, ha="left")
label(ax, m[5] + 0.85, ym, "poly(A) 꼬리", fs=11, c=TEXT2, ha="left")

# CDS 묶음
yb = ym - BH / 2 - 0.13
ax.plot([m[1], m[1], m[4], m[4]], [yb + 0.07, yb, yb, yb + 0.07], color=PRIMARY, lw=1.3, zorder=3)
label(ax, (m[1] + m[4]) / 2, yb - 0.22, "CDS: 번역되는 부분 → 단백질 (7장)", fs=12, w="bold", c=PRIMARY)

save(fig, os.path.join(HERE, "../../figures/ch03_gene_structure.png"))
