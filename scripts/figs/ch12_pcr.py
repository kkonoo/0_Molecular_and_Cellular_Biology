# ch12 §2 PCR: one cycle = three steps, molecules double per cycle; efficiency gap grows over 20 cycles  (redraw of figures/ch12_pcr.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.3, 4.1
fig, ax = new(W, H)

# ── left: three steps of one cycle ──────────────────────────────────────────
steps = [("① 변성", "두 가닥 분리"), ("② 결합", "primer가 붙음"), ("③ 신장", "새 가닥 합성")]
cw, ch, cx0, cy0 = 1.33, 0.40, 0.05, 3.58
for i, (name, sub) in enumerate(steps):
    x = cx0 + i * cw
    chevron(ax, x, cy0, cw - 0.04, ch, fc=STEPS[2], first=(i == 0))
    xc = x + (cw - 0.04) / 2 + (0.0 if i == 0 else 0.06)
    label(ax, xc, cy0 + ch / 2, name, fs=12, w="bold")
    label(ax, xc, cy0 - 0.2, sub, fs=11, c=TEXT2)

# ── left: molecules per cycle ───────────────────────────────────────────────
ORIG, NEW = PRIMARY, CAT[3]          # original strand / newly made strand
O, N = ORIG, NEW
rows = [
    (0, [(O, O)]),
    (1, [(O, N), (N, O)]),
    (2, [(O, N), (N, N), (N, N), (N, O)]),
    (3, [(O, N)] + [(N, N)] * 6 + [(N, O)]),
]
ys = [2.85, 2.25, 1.65, 1.05]
mw, step, mx0 = 0.26, 0.32, 0.95
bh = 0.065

def strand_pair(x, yc, top, bot):
    box(ax, x, yc + 0.025, mw, bh, fc=top, ec="none", r=0.02)
    box(ax, x, yc - 0.025 - bh, mw, bh, fc=bot, ec="none", r=0.02)

for (n, mols), y in zip(rows, ys):
    label(ax, 0.05, y, f"{n} 사이클", fs=12, w="bold", ha="left")
    for j, (t, b) in enumerate(mols):
        strand_pair(mx0 + j * step, y, t, b)
    label(ax, 4.0, y, f"{2 ** n}개", fs=12, w="bold", ha="right")

# key for strand colours + rule
ky = 0.45
box(ax, 0.05, ky - bh / 2, 0.26, bh, fc=ORIG, ec="none", r=0.02)
label(ax, 0.39, ky, "처음 가닥", fs=11, ha="left")
box(ax, 1.35, ky - bh / 2, 0.26, bh, fc=NEW, ec="none", r=0.02)
label(ax, 1.69, ky, "새 가닥", fs=11, ha="left")
label(ax, 4.0, ky, "n 사이클 → 2ⁿ개", fs=11, w="bold", ha="right")

# ── right: molecules vs cycle, log scale ───────────────────────────────────
px0, py0, pw, ph = 4.75, 0.62, 2.6, 3.3
pax = fig.add_axes([px0 / W, py0 / H, pw / W, ph / H])
pax.set_facecolor("none")
YL = (0.7, 2.5e6)
n = np.arange(0, 21)
pax.plot(n, 2.0 ** n, color=TEXT2, lw=1.8, ls=(0, (4, 3)), zorder=3)
pax.plot(n, 1.95 ** n, color=CAT[1], lw=2.2, zorder=4)
pax.plot(n, 1.80 ** n, color=CAT[2], lw=2.2, zorder=4)
pax.set_yscale("log")
pax.set_xlim(0, 20)
pax.set_ylim(*YL)
pax.set_yticks([10 ** k for k in range(7)])
pax.set_yticklabels(["1", "10"] + [f"$10^{k}$" for k in range(2, 7)], fontsize=11)
pax.minorticks_off()
pax.set_xticks([0, 5, 10, 15, 20])
pax.tick_params(labelsize=11)
pax.grid(axis="y", color=LINE, lw=0.8, zorder=0)
pax.set_xlabel("사이클 수", fontsize=11)
pax.set_ylabel("분자 수 (로그 눈금)", fontsize=11)

pax.text(0.04, 0.95, "효율 100%: 2²⁰ ≈ 105만", transform=pax.transAxes, fontsize=11,
         fontweight="bold", color=TEXT2, va="top")
pax.text(0.04, 0.86, "95%: 약 63만", transform=pax.transAxes, fontsize=11,
         fontweight="bold", color=CAT[1], va="top")
pax.text(0.04, 0.77, "80%: 약 13만", transform=pax.transAxes, fontsize=11,
         fontweight="bold", color=CAT[2], va="top")

# bracket + the one key point: ~5× gap after 20 cycles (drawn in inch coords on the canvas)
def yin(v):
    lo, hi = np.log10(YL[0]), np.log10(YL[1])
    return py0 + ph * (np.log10(v) - lo) / (hi - lo)

ya, yb = yin(1.8 ** 20), yin(1.95 ** 20)
bx = px0 + pw + 0.12
ax.plot([bx - 0.06, bx, bx, bx - 0.06], [ya, ya, yb, yb], color=TEXT, lw=1.3, zorder=5)
pill(ax, bx + 0.1, (ya + yb) / 2, "약 5배", fc=ACCENT, tc=ON_ACCENT, fs=12, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch12_pcr.png"))
