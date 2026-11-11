# ch09 §3 apoptosis의 외재·내재 경로가 실행 caspase로 모이는 흐름도  (redraw of figures/ch09_apoptosis.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

EXT = CAT[1]   # 외재 경로 (teal)
INT = CAT[0]   # 내재 경로 (blue)
EXE = CAT[2]   # 공통 실행 단계 (green)
ARW = "#9AA1AB"

fig, ax = new(8.3, 5.25)


def node(x, y, w, h, title, sub, color, solid=False, fc=None, ec=None):
    if solid:
        f, e = color, "none"
        tc = sc = on(color)
    else:
        f = fc or mix(color, 0.12)
        e = ec or color
        tc, sc = TEXT, TEXT2
    box(ax, x, y, w, h, fc=f, ec=e, lw=1.3, r=0.08, z=2)
    cx, cy = x + w / 2, y + h / 2
    label(ax, cx, cy + 0.12, title, fs=12.5, w="bold", c=tc)
    label(ax, cx, cy - 0.15, sub, fs=11, c=sc)


H = 0.62
ROW1, ROW2, ROW3 = 4.0, 2.95, 1.9

# ── panels ──
box(ax, 0.0, 1.7, 2.4, 3.45, fc=FAINT, ec="none", r=0.12, z=0)
box(ax, 2.55, 1.7, 5.7, 3.45, fc=FAINT, ec="none", r=0.12, z=0)
label(ax, 0.15, 4.88, "외재 경로", fs=13, w="bold", c=EXT, ha="left")
label(ax, 1.05, 4.87, "사멸 수용체 경로", fs=11, c=TEXT2, ha="left")
label(ax, 2.72, 4.88, "내재 경로", fs=13, w="bold", c=INT, ha="left")
label(ax, 3.62, 4.87, "미토콘드리아 경로", fs=11, c=TEXT2, ha="left")

# ── extrinsic column ──
LX, LW = 0.15, 2.1
LC = LX + LW / 2
node(LX, ROW1, LW, H, "사멸 리간드", "FAS 리간드 · TNF", EXT)
node(LX, ROW2, LW, H, "사멸 수용체", "FAS · TNFR1 (+ FADD)", EXT)
node(LX, ROW3, LW, H, "caspase-8", "개시 caspase", EXT, solid=True)
arrow(ax, LC, ROW1, LC, ROW2 + H, color=ARW)
arrow(ax, LC, ROW2, LC, ROW3 + H, color=ARW)

# ── intrinsic ──
SX, SW = 2.72, 1.45
BX, BW = 4.37, 1.9
RX, RW = 6.47, 1.63
RC = RX + RW / 2
node(SX, ROW1, SW, H, "스트레스", "DNA 손상 → p53", INT)
node(BX, ROW1, BW, H, "BH3-only", "PUMA·NOXA·BIM·BID", INT)
node(RX, ROW1, RW, H, "BAX · BAK", "외막에 구멍", INT)
node(RX, ROW2, RW, H, "cytochrome c", "세포질로 방출", INT)
node(RX, ROW3, RW, H, "caspase-9", "APAF1과 함께 · 개시", INT, solid=True)
node(BX, ROW2, BW, H, "BCL-2·BCL-XL·MCL-1", "BAX·BAK을 붙잡음", INT,
     fc="white", ec=STEPS[1])

yc = ROW1 + H / 2
arrow(ax, SX + SW, yc, BX, yc, color=ARW)
arrow(ax, BX + BW, yc, RX, yc, color=ARW)
arrow(ax, RC, ROW1, RC, ROW2 + H, color=ARW)
arrow(ax, RC, ROW2, RC, ROW3 + H, color=ARW)

# inhibition: BCL-2 family ⊣ BAX·BAK  (bar stops just short of the BAX·BAK box)
p0 = np.array([6.1, ROW2 + H])
p1 = np.array([6.56, ROW1 - 0.1])
ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=TEXT2, lw=1.8, zorder=3, solid_capstyle="butt")
d = (p1 - p0) / np.linalg.norm(p1 - p0)
perp = np.array([-d[1], d[0]]) * 0.1
ax.plot([p1[0] - perp[0], p1[0] + perp[0]], [p1[1] - perp[1], p1[1] + perp[1]],
        color=TEXT2, lw=2.2, zorder=3, solid_capstyle="butt")
label(ax, 6.12, 3.79, "억제", fs=11, c=TEXT2, ha="right")

# caspase-8 → BID (amplification into intrinsic path)
arrow(ax, LX + LW, ROW3 + 0.4, BX + 0.23, ROW1, color=EXT, lw=1.6, ls=(0, (4, 3)))
label(ax, 2.78, 2.5, "BID를 잘라 증폭", fs=11, c=EXT, ha="left", va="center")

# ── shared execution ──
EX, EW = 2.72, 2.1
EC = EX + EW / 2
DX, DW = 5.2, 2.9
BY, BH = 0.3, 0.65
node(EX, BY, EW, BH, "실행 caspase", "caspase-3 · caspase-7", EXE, solid=True)
node(DX, BY, DW, BH, "단백질 절단 → 세포 해체", "식세포가 먹음 · 염증 거의 없음", EXE)
ymid = BY + BH / 2
arrow(ax, EX + EW, ymid, DX, ymid, color=ARW)

# caspase-8 → 실행 caspase
ax.plot([LC, LC], [ROW3, ymid], color=ARW, lw=1.6, zorder=3, solid_capstyle="butt")
arrow(ax, LC, ymid, EX, ymid, color=ARW)
# caspase-9 → 실행 caspase
yj = 1.25
ax.plot([RC, RC, EC], [ROW3, yj, yj], color=ARW, lw=1.6, zorder=3,
        solid_joinstyle="miter", solid_capstyle="butt")
arrow(ax, EC, yj, EC, BY + BH, color=ARW)

save(fig, os.path.join(HERE, "../../figures/ch09_apoptosis.png"))
