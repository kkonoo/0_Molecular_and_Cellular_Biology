# ch07 §1-2) 리보솜이 mRNA 코돈을 읽고 tRNA가 아미노산을 가져오는 번역 모식도  (redraw of figures/ch07_translation.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *
import matplotlib.patches as mp

W, H = 7.6, 5.0
fig, ax = new(W, H)

C_MRNA = PRIMARY
C_TRNA = CAT[1]
C_AA = CAT[2]
C_AA_TXT = shade(CAT[2], 0.2)
RIBO_FC, RIBO_EC = STEPS[3], STEPS[2]

# ── mRNA line ───────────────────────────────────────────────────────────────
y_line = 1.9
x0, x1 = 0.4, 7.3
ax.plot([x0, x1], [y_line, y_line], color=C_MRNA, lw=3, solid_capstyle="round", zorder=4)
label(ax, x0 - 0.18, y_line, "5'", fs=13, w="bold", c=C_MRNA)
label(ax, x1 + 0.18, y_line, "3'", fs=13, w="bold", c=C_MRNA)
label(ax, 0.45, 2.45, "mRNA", fs=13, w="bold", c=C_MRNA, ha="left")

cw = 0.68
cds0 = 1.6
codons = ["AUG", "GUG", "AUU", "GGC", "UGG", "AAA", "UAA"]
aas = ["Met", "Val", "Ile", "Gly", "Trp", "Lys", "종결"]
cx = [cds0 + cw / 2 + cw * i for i in range(7)]
cds1 = cds0 + cw * 7

for k in range(8):
    xb = cds0 + cw * k
    ax.plot([xb, xb], [y_line - 0.07, y_line + 0.07], color=C_MRNA, lw=1.4, zorder=4)

y_cod = 2.08
for c, cod in zip(cx, codons):
    for j, ch in enumerate(cod):
        label(ax, c + (j - 1) * 0.19, y_cod, ch, fs=13, w="bold", c=C_MRNA)

label(ax, (x0 + cds0) / 2, y_cod, "5' UTR", fs=11, c=TEXT2)
label(ax, (cds1 + x1) / 2, y_cod, "3' UTR", fs=11, c=TEXT2)
label(ax, cx[0], 2.42, "시작 코돈", fs=12, w="bold", c=TEXT2)
label(ax, cx[6], 2.42, "종결 코돈", fs=12, w="bold", c=TEXT2)

# ── ribosome ────────────────────────────────────────────────────────────────
rc = (cx[3] + cx[4]) / 2
ax.add_patch(mp.Ellipse((rc, 3.05), 2.6, 1.6, fc=RIBO_FC, ec=RIBO_EC, lw=1.2, zorder=1))
ax.add_patch(mp.Ellipse((rc, 1.55), 2.6, 0.65, fc=RIBO_FC, ec=RIBO_EC, lw=1.2, zorder=1))

# leader labels (right)
ax.plot([rc + 1.3 * np.cos(np.radians(30)), 5.72], [3.05 + 0.8 * 0.5, 3.6], color="#9AA1AB", lw=0.9, zorder=2)
label(ax, 5.76, 3.6, "리보솜 큰 소단위", fs=12, w="bold", ha="left")
ax.plot([rc + 1.3 * np.cos(np.radians(-40)), 5.72], [1.55 + 0.325 * np.sin(np.radians(-40)), 1.25],
        color="#9AA1AB", lw=0.9, zorder=2)
label(ax, 5.76, 1.25, "리보솜 작은 소단위", fs=12, w="bold", ha="left")

# ── tRNAs ───────────────────────────────────────────────────────────────────
anti = {3: "CCG", 4: "ACC"}
for i, ac in anti.items():
    c = cx[i]
    box(ax, c - 0.29, 2.3, 0.58, 0.85, fc=mix(C_TRNA, 0.15), ec=C_TRNA, lw=1.2, r=0.08, z=3)
    label(ax, c, 2.9, "tRNA", fs=11, c=shade(C_TRNA, 0.1))
    for j, ch in enumerate(ac):
        label(ax, c + (j - 1) * 0.17, 2.47, ch, fs=12, w="bold", c=shade(C_TRNA, 0.1))
    ax.plot([c, c], [3.15, 3.25], color=C_TRNA, lw=1.6, zorder=3)

# anticodon label (left)
ax.plot([2.95, cx[3] - 0.3], [2.62, 2.5], color="#9AA1AB", lw=0.9, zorder=2)
label(ax, 2.9, 2.62, "안티코돈", fs=12, w="bold", c=shade(C_TRNA, 0.1), ha="right")

# ── amino-acid chain ────────────────────────────────────────────────────────
r = 0.2
chain = [("Met", 2.84, 4.55), ("Val", 3.22, 4.2), ("Ile", 3.6, 3.85), ("Gly", cx[3], 3.45)]
for (a, xa, ya), (b, xb, yb) in zip(chain[:-1], chain[1:]):
    ax.plot([xa, xb], [ya, yb], color=TEXT2, lw=2, zorder=3.5)
for name, xa, ya in chain + [("Trp", cx[4], 3.45)]:
    ax.add_patch(mp.Circle((xa, ya), r, fc=C_AA, ec="white", lw=1.2, zorder=4))
    label(ax, xa, ya, name, fs=11, w="bold", c="white", z=6)

ax.plot([cx[3] + r, cx[4] - r], [3.45, 3.45], color=TEXT2, lw=1.4, ls=(0, (2, 2)), zorder=3.5)
label(ax, rc, 3.75, "다음 결합", fs=11, c=TEXT2)

label(ax, 2.56, 4.64, "아미노산 사슬", fs=12, w="bold", ha="right")
label(ax, 2.56, 4.4, "N-말단", fs=11, c=TEXT2, ha="right")

# ── translation read-out ────────────────────────────────────────────────────
y_aa = 0.8
label(ax, (x0 + cds0) / 2, y_aa, "번역", fs=11, c=TEXT2)
for c, a in zip(cx, aas):
    col = TEXT if a == "종결" else C_AA_TXT
    label(ax, c, y_aa, a, fs=12, w="bold", c=col)

# movement direction
arrow(ax, 2.6, 0.35, 4.55, 0.35, color=TEXT2, lw=1.6)
label(ax, 4.7, 0.35, "리보솜 이동: 한 코돈씩 5'→3'", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch07_translation.png"))
