# ch02 §3 뉴클레오타이드 하나(5탄당 1′–5′, 인산, 염기)와 역평행 이중 가닥 사다리  (redraw of figures/ch02_nucleotide.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.3, 4.5
fig, ax = new(W, H)

GRAY = "#9AA1AB"
P_FC = TEXT2                      # phosphate (same in both panels)
SUGAR_FC, SUGAR_EC = FAINT, GRAY  # sugar (same in both panels)
BASE_C = {"A": CAT[0], "G": CAT[1], "C": CAT[2], "T": CAT[3]}

# ── panel markers ──────────────────────────────────────────────────────────
label(ax, 0.3, 4.3, "① 뉴클레오타이드 하나", fs=12, w="bold", ha="left")
label(ax, 4.0, 4.3, "② 역평행 이중 가닥", fs=12, w="bold", ha="left")
ax.plot([3.8, 3.8], [0.1, 4.3], color=LINE, lw=1.0, zorder=1)

# ── left: one nucleotide ───────────────────────────────────────────────────
cx, cy, r = 1.75, 2.55, 0.62
ang = np.deg2rad([90, 18, -54, -126, 162])      # O, 1′, 2′, 3′, 4′
vx = cx + r * np.cos(ang); vy = cy + r * np.sin(ang)
ax.add_patch(mp.Polygon(list(zip(vx, vy)), closed=True, fc=SUGAR_FC, ec=SUGAR_EC,
                        lw=1.6, zorder=2))
label(ax, cx, cy - 0.05, "5탄당", fs=12, w="bold")
label(ax, vx[0], vy[0] + 0.14, "O", fs=11, c=TEXT2)

c1 = (vx[1], vy[1]); c2 = (vx[2], vy[2]); c3 = (vx[3], vy[3]); c4 = (vx[4], vy[4])
for (x, y) in [c1, c2, c3, c4]:
    ax.add_patch(mp.Circle((x, y), 0.04, fc=TEXT2, ec="none", zorder=4))

# 4′ → 5′ carbon → phosphate
c5 = (c4[0] - 0.27, c4[1] + 0.52)
P = (c5[0] - 0.3, c5[1] + 0.5)
ax.plot([c4[0], c5[0], P[0]], [c4[1], c5[1], P[1]], color=TEXT2, lw=1.4, zorder=1)
ax.add_patch(mp.Circle(c5, 0.04, fc=TEXT2, ec="none", zorder=4))
ax.add_patch(mp.Circle(P, 0.22, fc=P_FC, ec="none", zorder=4))
label(ax, P[0], P[1], "P", fs=12, w="bold", c="white")
label(ax, P[0] + 0.3, P[1] + 0.03, "인산", fs=12, w="bold", ha="left")

label(ax, c5[0] + 0.12, c5[1] + 0.02, "5′", fs=12, w="bold", ha="left")
label(ax, c4[0] - 0.1, c4[1] - 0.06, "4′", fs=12, w="bold", ha="right")
label(ax, c1[0] + 0.02, c1[1] - 0.2, "1′", fs=12, w="bold", ha="left")
label(ax, c2[0] + 0.12, c2[1] + 0.06, "2′", fs=12, w="bold", ha="left")
label(ax, c3[0] - 0.12, c3[1] + 0.06, "3′", fs=12, w="bold", ha="right")

# 3′-OH and 2′ group
ax.plot([c3[0], c3[0] - 0.12], [c3[1], c3[1] - 0.36], color=TEXT2, lw=1.4, zorder=1)
label(ax, c3[0] - 0.14, c3[1] - 0.5, "OH", fs=12, w="bold")
ax.plot([c2[0], c2[0] + 0.12], [c2[1], c2[1] - 0.36], color=TEXT2, lw=1.4, zorder=1)
label(ax, c2[0] + 0.16, c2[1] - 0.5, "H 또는 OH", fs=12, w="bold")
label(ax, c2[0] + 0.16, c2[1] - 0.78, "DNA: H · RNA: OH", fs=11, c=TEXT2)

# 1′ → base
bx, bw, bh = 2.75, 0.85, 0.44
ax.plot([c1[0], bx], [c1[1], c1[1]], color=TEXT2, lw=1.4, zorder=1)
box(ax, bx, c1[1] - bh / 2, bw, bh, fc=FAINT, ec=GRAY, lw=1.2, r=0.08)
label(ax, bx + bw / 2, c1[1], "염기", fs=12, w="bold")
letters = [("A", BASE_C["A"]), ("G", BASE_C["G"]), ("C", BASE_C["C"]), ("T/U", BASE_C["T"])]
xs = [bx + 0.0, bx + 0.22, bx + 0.44, bx + 0.72]
for (t, c), x in zip(letters, xs):
    label(ax, x, c1[1] - 0.42, t, fs=12, w="bold", c=shade(c, 0.15) if t == "T/U" else c)

# ── right: antiparallel ladder ─────────────────────────────────────────────
xl, xr = 4.45, 7.75
ys = [3.55, 2.95, 2.35, 1.75, 1.15]
left_seq = "AGCTG"
pair = {"A": "T", "T": "A", "G": "C", "C": "G"}
nh = {"A": 2, "T": 2, "G": 3, "C": 3}
sq, bh2 = 0.22, 0.34
bl0, bl1 = xl + 0.18, 5.65
br0, br1 = 6.55, xr - 0.18

for x in (xl, xr):
    ax.plot([x, x], [ys[0] + 0.2, ys[-1] - 0.2], color=GRAY, lw=1.6, zorder=1)
    for i in range(len(ys) - 1):
        ax.add_patch(mp.Circle((x, (ys[i] + ys[i + 1]) / 2), 0.06, fc=P_FC, ec="none", zorder=4))

for y, b in zip(ys, left_seq):
    p = pair[b]
    for x in (xl, xr):
        box(ax, x - sq / 2, y - sq / 2, sq, sq, fc=SUGAR_FC, ec=SUGAR_EC, lw=1.4, r=0.04, z=4)
    box(ax, bl0, y - bh2 / 2, bl1 - bl0, bh2, fc=BASE_C[b], ec="none", r=0.05, z=3)
    label(ax, (bl0 + bl1) / 2, y, b, fs=12, w="bold", c=on(BASE_C[b]))
    box(ax, br0, y - bh2 / 2, br1 - br0, bh2, fc=BASE_C[p], ec="none", r=0.05, z=3)
    label(ax, (br0 + br1) / 2, y, p, fs=12, w="bold", c=on(BASE_C[p]))
    n = nh[b]
    offs = [-0.065, 0.065] if n == 2 else [-0.1, 0.0, 0.1]
    for o in offs:
        ax.plot([bl1 + 0.05, br0 - 0.05], [y + o, y + o], color=TEXT2, lw=1.3,
                ls=(0, (3, 2)), zorder=2)

# strand ends and directions
top, bot = ys[0] + 0.42, ys[-1] - 0.42
label(ax, xl, top, "5′", fs=13, w="bold")
label(ax, xl, bot, "3′", fs=13, w="bold")
label(ax, xr, top, "3′", fs=13, w="bold")
label(ax, xr, bot, "5′", fs=13, w="bold")
arrow(ax, 4.12, ys[0] + 0.15, 4.12, ys[-1] - 0.15, color=TEXT2, lw=1.6)
arrow(ax, 8.08, ys[-1] - 0.15, 8.08, ys[0] + 0.15, color=TEXT2, lw=1.6)

# direct key under the ladder
ky1, ky2 = 0.38, 0.05
kx = 4.3
box(ax, kx - 0.08, ky1 - 0.08, 0.16, 0.16, fc=SUGAR_FC, ec=SUGAR_EC, lw=1.4, r=0.03, z=4)
ax.add_patch(mp.Circle((kx + 0.2, ky1), 0.06, fc=P_FC, ec="none", zorder=4))
label(ax, kx + 0.4, ky1, "당–인산 뼈대: 공유결합", fs=11, ha="left")
ax.plot([kx - 0.08, kx + 0.27], [ky2, ky2], color=TEXT2, lw=1.3, ls=(0, (3, 2)))
label(ax, kx + 0.4, ky2, "수소결합: A–T 2개, G–C 3개", fs=11, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch02_nucleotide.png"))
