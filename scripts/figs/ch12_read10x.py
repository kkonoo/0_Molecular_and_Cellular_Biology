# ch12 §read10x 10x 3' library 조각 하나의 구조와 R1·R2·I1·I2 FASTQ 파일의 관계  (redraw of figures/ch12_read10x.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.7, 4.35
fig, ax = new(W, H)

C_BC, C_UMI, C_CDNA, C_IDX = CAT[0], CAT[3], CAT[1], CAT[2]
MONO = "DejaVu Sans Mono"

# ── library 조각 막대 ──────────────────────────────────────────────
segs = [  # (name, sub, width, fill, edge)
    ("어댑터", None, 0.80, LINE, "none"),
    ("세포 바코드", "16 bp", 1.30, C_BC, "none"),
    ("UMI", "12 bp", 0.95, C_UMI, "none"),
    ("poly(T)", None, 0.90, FAINT, "#C9CED6"),
    ("cDNA", None, 2.35, C_CDNA, "none"),
    ("어댑터", None, 0.80, LINE, "none"),
]
gap = 0.04
y0, bh = 3.0, 0.62
x = 0.3
edges = []
for name, sub, w, fc, ec in segs:
    box(ax, x, y0, w, bh, fc=fc, ec=ec, lw=1.0, r=0.04)
    tc = on(fc)
    if sub:
        label(ax, x + w / 2, y0 + bh / 2 + 0.11, name, fs=12, w="bold", c=tc)
        label(ax, x + w / 2, y0 + bh / 2 - 0.13, sub, fs=11, c=tc)
    else:
        label(ax, x + w / 2, y0 + bh / 2, name, fs=12, w="bold", c=tc)
    edges.append((x, x + w))
    x += w + gap

adL, bc, umi, pt, cd, adR = edges

# ── Read 1 (위, 오른쪽으로) ─────────────────────────────────────────
yr1 = 3.85
arrow(ax, bc[0], yr1, umi[1], yr1, color=C_BC, lw=1.8, ms=13)
ax.plot([umi[1] + 0.03, pt[1]], [yr1, yr1], color=C_BC, lw=1.6, ls=(0, (3, 2.5)), zorder=3)
label(ax, bc[0], yr1 + 0.27, "Read 1 · 16 + 12 = 28 사이클 (v3·v4)", fs=12, w="bold",
      c=C_BC, ha="left")

# ── Read 2 (아래, 왼쪽으로) ─────────────────────────────────────────
yr2 = 2.6
arrow(ax, cd[1] - 0.05, yr2, cd[0] + 0.35, yr2, color=C_CDNA, lw=1.8, ms=13)
label(ax, (cd[0] + cd[1]) / 2 + 0.15, yr2 - 0.3, "Read 2", fs=12, w="bold", c=C_CDNA)

# ── index read ───────────────────────────────────────────────────────
for (a, b), nm in [(adL, "I2"), (adR, "I1")]:
    cx = (a + b) / 2
    pill(ax, cx, yr2, nm, fc=C_IDX, fs=12, h=0.3, padx=0.18)
    label(ax, cx, yr2 - 0.3, "시료 index", fs=11, c=TEXT2)

# ── FASTQ 파일 ───────────────────────────────────────────────────────
box(ax, 0.3, 0.15, adR[1] - 0.3, 1.6, fc=FAINT, ec="none", r=0.1, z=1)
rows = [
    ("sample_S1_L001_R1_001.fastq.gz", "세포 바코드 + UMI", C_BC),
    ("sample_S1_L001_R2_001.fastq.gz", "cDNA → 정렬해서 유전자를 정함", C_CDNA),
    ("sample_S1_L001_I1_001.fastq.gz", "시료 index (I2 파일도 같은 꼴)", C_IDX),
]
for i, (fn, what, c) in enumerate(rows):
    yy = 1.42 - i * 0.45
    ink = shade(c, 0.25) if c == C_UMI else c
    ax.text(0.55, yy, fn, fontsize=11, fontfamily=MONO, fontweight="bold", color=ink,
            ha="left", va="center", zorder=5)
    label(ax, 4.45, yy, what, fs=11.5, c=TEXT, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch12_read10x.png"))
