# ch00 §3 count matrix 일부: 행 = 유전자, 열 = 세포 바코드, 칸 = UMI 수, 마지막 줄 = 열 합  (redraw of figures/ch00_matrix.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

genes = ["CD3E", "MS4A1", "LYZ", "NKG7", "FOXP3"]
barcodes = ["AAACCTGA…", "AAACCTGC…", "AAACGGGA…",
            "AAAGATGT…", "AAAGCAAG…", "AACACGTC…"]
M = [[3, 5, 0, 0, 0, 2],
     [0, 0, 4, 0, 0, 0],
     [0, 1, 0, 18, 0, 0],
     [0, 0, 0, 0, 7, 1],
     [0, 0, 0, 0, 0, 1]]
sums = ["3,105", "4,422", "2,890", "6,317", "3,560", "2,974"]
HL = (0, 1)   # CD3E × 둘째 세포

W, H = 8.0, 4.3
fig, ax = new(W, H)

x0, cw = 1.15, 1.1          # matrix left edge, column width
top, rh = 3.45, 0.42        # matrix top, row height
nR, nC = len(genes), len(barcodes)
xr = x0 + nC * cw

# axis labels
label(ax, (x0 + xr) / 2, 4.10, "열 = 세포 (세포 바코드)", fs=12, w="bold", c=TEXT2)
label(ax, x0 - 0.12, 3.70, "행 = 유전자", fs=12, w="bold", c=TEXT2, ha="right")

# barcodes (shortened so neighbouring columns keep a clear gap)
for j, b in enumerate(barcodes):
    label(ax, x0 + (j + 0.5) * cw, 3.70, b, fs=11, c=TEXT2)

# cells
for i, g in enumerate(genes):
    y = top - (i + 1) * rh
    label(ax, x0 - 0.12, y + rh / 2, g, fs=12, ha="right")
    for j in range(nC):
        v = M[i][j]
        x = x0 + j * cw
        if (i, j) == HL:
            fc, tc, wt = ACCENT, ON_ACCENT, "bold"
        elif v > 0:
            fc, tc, wt = STEPS[3], TEXT, "bold"
        else:
            fc, tc, wt = "white", "#B4BAC3", "regular"
        ax.add_patch(mp.Rectangle((x, y), cw, rh, fc=fc, ec=LINE, lw=1.0, zorder=2))
        label(ax, x + cw / 2, y + rh / 2, str(v), fs=12, w=wt, c=tc)

# highlighted cell outline on top of grid
hx, hy = x0 + HL[1] * cw, top - (HL[0] + 1) * rh
ax.add_patch(mp.Rectangle((hx, hy), cw, rh, fc="none", ec=ACCENT_INK, lw=1.6, zorder=4))

# "more rows" dots
ybot = top - nR * rh
for j in range(nC):
    cx = x0 + (j + 0.5) * cw
    ax.scatter([cx] * 3, [ybot - 0.10, ybot - 0.18, ybot - 0.26], s=4, c=TEXT2, zorder=5)
ax.scatter([x0 - 0.45] * 3, [ybot - 0.10, ybot - 0.18, ybot - 0.26], s=4, c=TEXT2, zorder=5)

# column sums
yl = ybot - 0.40
ax.plot([x0, xr], [yl, yl], color=TEXT2, lw=1.2, zorder=3)
ys = yl - 0.28
label(ax, x0 - 0.12, ys, "열 합 (총 UMI)", fs=12, w="bold", c=PRIMARY, ha="right")
for j, s in enumerate(sums):
    label(ax, x0 + (j + 0.5) * cw, ys, s, fs=12, w="bold", c=PRIMARY)

# note for the highlighted cell
yn = ys - 0.48
box(ax, x0, yn - 0.11, 0.30, 0.22, fc=ACCENT, ec=ACCENT_INK, lw=1.2, r=0.03)
label(ax, x0 + 0.42, yn, "둘째 세포에서 잡힌 CD3E mRNA 5분자 = UMI 5개", fs=12, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch00_matrix.png"))
