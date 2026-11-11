# ch00 §4.1 세포 속 mRNA가 포획·PCR·시퀀싱·UMI 합치기를 거쳐 count가 되는 과정  (redraw of figures/ch00_counting.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.3, 4.45
fig, ax = new(W, H)

GA, GB, GC = CAT[0], CAT[1], CAT[2]
EDGE = "#9AA1AB"

# ── panel geometry ──────────────────────────────────────────────────────────
pw, ph, gap = 1.38, 2.3, 0.3
x0, yb = 0.03, 1.65
cy = yb + ph / 2
xs = [x0 + i * (pw + gap) for i in range(5)]
heads = ["① 세포 속 mRNA", "② 포획 + 이름표", "③ PCR로 불리기", "④ 시퀀싱 → read", "⑤ UMI로 합치기"]

for i, (x, t) in enumerate(zip(xs, heads)):
    box(ax, x, yb, pw, ph, fc=FAINT, ec=LINE, lw=1.0, r=0.1, z=1)
    label(ax, x + pw / 2, yb + ph + 0.2, t, fs=12, w="bold")
    if i < 4:
        arrow(ax, x + pw + 0.06, cy, x + pw + gap - 0.06, cy, ms=11)


def mol(xc, y, c, w=0.36):
    box(ax, xc - w / 2, y - 0.045, w, 0.09, fc=c, ec="none", r=0.045, z=4)


def tagged(xc, y, u, c):
    """UMI tag + molecule, centred at xc."""
    xl = xc - 0.42
    box(ax, xl, y - 0.095, 0.3, 0.19, fc="white", ec=EDGE, lw=0.8, r=0.04, z=4)
    label(ax, xl + 0.15, y - 0.005, u, fs=11, c=TEXT2)
    box(ax, xl + 0.36, y - 0.045, 0.48, 0.09, fc=c, ec="none", r=0.045, z=4)


# ① cell with 10 mRNA molecules: A 6, B 3, C 1
c1 = xs[0] + pw / 2
ax.add_patch(mp.Ellipse((c1, cy), 1.22, 1.95, fc="white", ec=EDGE, ls="--", lw=1.0, zorder=2))
order = [GA, GA, GB, GA, GC, GA, GB, GA, GB, GA]
for k, c in enumerate(order):
    dy = 0.72 - k * 0.16
    dx = -0.18 if k % 2 == 0 else 0.18
    mol(c1 + dx, cy + dy, c)

# ② captured: u1, u2 from A, u3 from B; C not captured
c2 = xs[1] + pw / 2
for y, u, c in [(cy + 0.6, "u1", GA), (cy + 0.15, "u2", GA), (cy - 0.3, "u3", GB)]:
    tagged(c2, y, u, c)
label(ax, c2, yb + 0.25, "C는 안 잡힘", fs=11, c=TEXT2)

# ③ PCR copies: u1 ×3, u2 ×2, u3 ×4
c3 = xs[2] + pw / 2
rows = [("u1", GA)] * 3 + [("u2", GA)] * 2 + [("u3", GB)] * 4
for k, (u, c) in enumerate(rows):
    tagged(c3, cy + 0.86 - k * 0.215, u, c)

# ④ reads
c4 = xs[3] + pw / 2
for y, t, c in [(cy + 0.55, "A: read 5개", GA), (cy + 0.1, "B: read 4개", GB),
                (cy - 0.35, "C: read 0개", GC)]:
    label(ax, c4, y, t, fs=12, w="bold", c=c)
label(ax, c4, yb + 0.25, "UMI도 함께 읽힘", fs=11, c=TEXT2)

# ⑤ UMI counts
c5 = xs[4] + pw / 2
for y, t, c in [(cy + 0.55, "A: UMI 2개", GA), (cy + 0.1, "B: UMI 1개", GB),
                (cy - 0.35, "C: UMI 0개", GC)]:
    label(ax, c5, y, t, fs=12, w="bold", c=c)
label(ax, c5, yb + 0.25, "→ count matrix", fs=11, c=TEXT2)

# ── table ───────────────────────────────────────────────────────────────────
lx = 3.3                      # right edge of row-label column
vx = [3.95, 5.0, 6.05]        # value column centres
tl, tr = 1.5, 6.6
yh = 1.27
for x, g, c in zip(vx, "ABC", [GA, GB, GC]):
    label(ax, x, yh, f"유전자 {g}", fs=12, w="bold", c=c)
ax.plot([tl, tr], [yh - 0.18, yh - 0.18], color=EDGE, lw=1.0, zorder=3)

table = [("세포 속 실제 분자 수", [6, 3, 1]),
         ("read 수", [5, 4, 0]),
         ("UMI 수 (count)", [2, 1, 0])]
for r, (name, vals) in enumerate(table):
    y = 0.9 - r * 0.32
    last = r == len(table) - 1
    if last:   # the one take-away: count matrix holds UMI counts
        box(ax, tl, y - 0.15, tr - tl, 0.3, fc=ACCENT, ec="none", r=0.06, z=1)
    tc = ON_ACCENT if last else TEXT
    wt = "bold" if last else "regular"
    label(ax, lx, y, name, fs=12, w=wt, c=tc, ha="right")
    for x, v in zip(vx, vals):
        label(ax, x, y, str(v), fs=12, w=wt, c=tc)

save(fig, os.path.join(HERE, "../../figures/ch00_counting.png"))
