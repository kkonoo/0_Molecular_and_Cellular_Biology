# ch02 §reverse complement 좌표 101–115의 + 가닥·− 가닥과 complement / reverse complement  (redraw of figures/ch02_strands.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

plus = "ATGGCGTACGTTAGC"
pair = {"A": "T", "T": "A", "G": "C", "C": "G"}
minus = "".join(pair[b] for b in plus)          # complement만: 3'→5' 순서
rc = minus[::-1]                                # reverse complement: 5'→3'

# 네 염기 = 나란한 항목 → CAT 네 색, 그림 전체에서 같은 글자 = 같은 색
BASE = {"A": CAT[0], "C": CAT[1], "G": CAT[2], "T": CAT[3]}

W, H = 7.8, 4.45
fig, ax = new(W, H)

X0 = 2.2          # 첫 칸 왼쪽
P = 0.36          # 칸 간격
S = 0.30          # 칸 크기
XE = X0 + 14 * P + S

def cx(i):
    return X0 + i * P + S / 2

def row(y, seq, solid, end_l, end_r):
    for i, b in enumerate(seq):
        c = BASE[b]
        fc = c if solid else mix(c, 0.25)
        box(ax, cx(i) - S / 2, y - S / 2, S, S, fc=fc, ec="none", r=0.05)
        label(ax, cx(i), y, b, fs=12, w="bold", c=on(fc) if solid else TEXT)
    label(ax, X0 - 0.08, y, end_l, fs=12, w="bold", ha="right")
    label(ax, XE + 0.08, y, end_r, fs=12, w="bold", ha="left")

def rowlabel(y, main, sub=None):
    label(ax, 0.1, y + 0.1, main, fs=12, w="bold", ha="left")
    if sub:
        label(ax, 0.1, y - 0.15, sub, fs=11, c=TEXT2, ha="left")

# 레퍼런스 좌표
yc = 4.2
label(ax, 0.1, yc, "레퍼런스 좌표", fs=11, c=TEXT2, ha="left")
for i in range(0, 15, 2):
    label(ax, cx(i), yc, str(101 + i), fs=11, c=TEXT2)

# + 가닥 / − 가닥
yp, ym = 3.7, 3.1
for i in range(15):
    ax.plot([cx(i), cx(i)], [yp - S / 2 - 0.02, ym + S / 2 + 0.02], color="#B8BEC7",
            lw=0.9, ls=(0, (2, 2)), zorder=1)
row(yp, plus, True, "5′", "3′")
row(ym, minus, False, "3′", "5′")
rowlabel(yp, "+ 가닥", "FASTA에 저장되는 쪽")
rowlabel(ym, "− 가닥", "+ 가닥에서 계산")

# − 가닥 유전자의 방향 (오른쪽 → 왼쪽)
ya = 2.62
arrow(ax, XE, ya, X0, ya, color=TEXT2, lw=1.8, ms=14)
label(ax, (X0 + XE) / 2, ya - 0.25, "strand = − 인 유전자의 5′→3′ 방향", fs=11, c=TEXT2)

# 구분선
ax.plot([0.1, XE + 0.3], [2.05, 2.05], color=LINE, lw=1.2, zorder=1)

# complement만 / reverse complement
y1, y2 = 1.55, 0.85
row(y1, minus, False, "3′", "5′")
row(y2, rc, True, "5′", "3′")
rowlabel(y1, "complement만", "− 가닥을 3′→5′로")
label(ax, 0.1, y2 + 0.1, "reverse complement", fs=12, w="bold", ha="left")
pill(ax, 0.1, y2 - 0.2, "− 가닥을 5′→3′로", fc=ACCENT, tc=ON_ACCENT, fs=11, h=0.28, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch02_strands.png"))
