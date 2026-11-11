# ch09 §1 세포주기 네 단계(G1·S·G2·M)의 시간 비율, G0 이탈, 검문점 세 곳, 단계별 DNA 양  (redraw of figures/ch09_cell_cycle.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(7.6, 6.1)

cx, cy = 3.4, 3.0
R_OUT, RW = 1.9, 0.42
R_IN = R_OUT - RW
R_MID = R_OUT - RW / 2


def pt(theta, r):
    """Compass angle (degrees, clockwise from top) → (x, y)."""
    t = np.radians(theta)
    return cx + r * np.sin(t), cy + r * np.cos(t)


# phases: name, hours, colour  (24 h example: 11 / 8 / 4 / 1)
phases = [("G1", 11, CAT[0]), ("S", 8, CAT[1]), ("G2", 4, CAT[2]), ("M", 1, CAT[3])]
ink = {"G1": CAT[0], "S": CAT[1], "G2": CAT[2], "M": shade(CAT[3], 0.35)}
letter_at = {"M": 350.6}   # keep the M letter clear of the in-M checkpoint tick

bounds = {}
th = 0.0
for name, hrs, col in phases:
    span = hrs / 24 * 360
    a, b = th, th + span
    bounds[name] = (a, b)
    ax.add_patch(mp.Wedge((cx, cy), R_OUT, 90 - b, 90 - a, width=RW,
                          fc=col, ec="white", lw=2.2, zorder=2))
    xm, ym = pt(letter_at.get(name, (a + b) / 2), R_MID)
    label(ax, xm, ym, name, fs=13, w="bold", c=on(col), z=4)
    th = b

# DNA amount inside the ring (phase colour)
for name, txt, r in [("G1", "2C", 1.12), ("S", "2C → 4C", 1.05), ("G2", "4C", 1.05)]:
    a, b = bounds[name]
    x, y = pt((a + b) / 2, r)
    label(ax, x, y, txt, fs=12, w="bold", c=ink[name])

# centre
label(ax, cx, cy + 0.13, "세포주기", fs=13, w="bold")
label(ax, cx, cy - 0.16, "24시간 주기 예시", fs=11, c=TEXT2)

# phase descriptions (direct labels)
def desc(x, y, title, sub, c, ha):
    label(ax, x, y, title, fs=12, w="bold", c=c, ha=ha)
    label(ax, x, y - 0.25, sub, fs=11, c=TEXT2, ha=ha)

desc(5.5, 3.42, "G1 · 자라며 분열 결정", "약 11시간", ink["G1"], "left")
desc(1.62, 1.42, "S · DNA 복제", "약 8시간", ink["S"], "right")
desc(1.62, 4.62, "G2 · 분열 준비", "약 4시간", ink["G2"], "right")
desc(3.1, 5.82, "M · 딸세포 2개로", "약 1시간", ink["M"], "center")

# checkpoints: dark ticks across the ring + numbered circles
CK_R = 2.27
checkpoints = [(1, 165.0, 165.0), (2, 345.0, 332.0), (3, 357.2, 8.0)]  # n, tick angle, circle angle
for n, t_tick, t_circ in checkpoints:
    x1, y1 = pt(t_tick, R_IN - 0.06)
    x2, y2 = pt(t_tick, R_OUT + 0.07)
    ax.plot([x1, x2], [y1, y2], color=TEXT, lw=2.4, solid_capstyle="butt", zorder=5)
    xc, yc = pt(t_circ, CK_R)
    ax.plot([x2, xc], [y2, yc], color=TEXT, lw=1.0, zorder=4)
    num(ax, xc, yc, n, fc=TEXT, r=0.14, fs=11, z=6)

# checkpoint key (numbers → names)
for i, (n, name) in enumerate([(1, "G1/S 검문점"), (2, "G2/M 검문점"), (3, "M기 안 검문점")]):
    y = 1.65 - i * 0.4
    num(ax, 5.62, y, n, fc=TEXT, r=0.14, fs=11)
    label(ax, 5.85, y, name, fs=11, ha="left")

# G0: leave the cycle from G1 and come back
g0x, g0y = 5.5, 5.12
ax.add_patch(mp.Ellipse((g0x, g0y), 1.15, 0.64, fc=FAINT, ec=TEXT2, lw=1.1, ls="--", zorder=3))
label(ax, g0x, g0y + 0.11, "G0", fs=12, w="bold")
label(ax, g0x, g0y - 0.14, "분열을 멈춤", fs=11, c=TEXT2)
xo, yo = pt(26, R_OUT + 0.06)
arrow(ax, xo, yo, 4.93, 5.08, color=TEXT2, lw=1.4, rad=-0.15)
xb, yb = pt(42, R_OUT + 0.06)
arrow(ax, 5.25, 4.82, xb, yb, color=TEXT2, lw=1.2, ls="--", rad=-0.15)
label(ax, 6.15, 5.24, "신호가 오면", fs=11, c=TEXT2, ha="left")
label(ax, 6.15, 4.99, "돌아오기도 함", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch09_cell_cycle.png"))
