# ch10 §2-2 분화 능력의 다섯 단계와 iPSC 재프로그래밍  (redraw of figures/ch10_potency.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(7.4, 4.4)

rows = [
    ("전능성", "totipotent", "수정란", "몸 전체 + 태반 쪽 조직"),
    ("만능성", "pluripotent", "내세포괴, ESC, iPSC", "세 배엽의 모든 세포\n(태반은 못 만듦)"),
    ("다능성", "multipotent", "조혈줄기세포 등 성체줄기세포", "한 계통 안의 여러 세포\n(예: 모든 혈액세포)"),
    ("단능성", "unipotent", "", "한 종류만"),
    ("분화된 세포", "differentiated", "B세포, T세포 등", "자기 일을 함\n대개 다른 종류로 바뀌지 않음"),
]
bar_w = [2.2, 1.8, 1.35, 0.85, 0.3]
cols = ramp(5)

y_top = 3.95
dy = 0.8
bh = 0.7
x0, x1 = 0.6, 6.6
bar_x = 2.05
ys = []

for i, (ko, en, ex, rng) in enumerate(rows):
    y = y_top - i * dy
    ys.append(y)
    box(ax, x0, y - bh / 2, x1 - x0, bh, fc=FAINT, ec="none", r=0.08, z=1)
    # name
    label(ax, 0.78, y + 0.12, ko, fs=13, w="bold", ha="left")
    label(ax, 0.78, y - 0.15, en, fs=11, c=TEXT2, ha="left")
    # conceptual bar
    by = y + 0.1 if ex else y
    box(ax, bar_x, by - 0.11, bar_w[i], 0.22, fc=cols[i], ec="none", r=0.05, z=2)
    if ex:
        label(ax, bar_x, y - 0.19, ex, fs=11, c=TEXT2, ha="left")
    # what it can make
    label(ax, 4.5, y, rng, fs=11, ha="left", linespacing=1.3)

top = ys[0] + bh / 2
bot = ys[-1] - bh / 2

# left: differentiation direction
arrow(ax, 0.4, top, 0.4, bot, color=TEXT2, lw=1.6, ms=13)
label(ax, 0.2, (top + bot) / 2, "분화 → 선택지가 줄어듦", fs=11, c=TEXT2, rotation=90)

# right: reprogramming (iPSC) — the one reverse path
xr = 6.8
ax.plot([x1 + 0.03, xr, xr], [ys[4], ys[4], ys[1]], color=ACCENT_INK, lw=1.6,
        ls=(0, (4, 3)), zorder=3, solid_capstyle="butt")
arrow(ax, xr, ys[1], x1 + 0.04, ys[1], color=ACCENT_INK, lw=1.6, ms=13)
label(ax, xr + 0.27, (ys[1] + ys[4]) / 2, "재프로그래밍: iPSC\n(전사인자 4개)", fs=11,
      w="bold", c=ACCENT_INK, rotation=90, linespacing=1.3)

save(fig, os.path.join(HERE, "../../figures/ch10_potency.png"))
