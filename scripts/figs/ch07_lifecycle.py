# ch07 §2 단백질의 일생 다섯 단계(합성→접힘→수식→위치→분해)와 각 단계가 데이터에서 보이는지  (redraw of figures/ch07_lifecycle.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.15, 4.75
fig, ax = new(W, H)

n = 5
gap = 0.22
cw = (W - (n - 1) * gap) / n
cols = ramp(n)
xs = [i * (cw + gap) for i in range(n)]
cxs = [x + cw / 2 for x in xs]

card_y0, card_y1 = 2.15, 4.62
titles = ["① 합성", "② 접힘", "③ 수식 (PTM)", "④ 위치", "⑤ 분해"]
descs = ["리보솜이 mRNA를\n읽어 사슬을 만듦",
         "3차원 모양으로\n접힘",
         "화학 그룹을\n붙이거나 뗌",
         "서열 속 신호를\n따라 제자리로",
         "유비퀴틴 표지 →\n프로테아좀이 분해"]

for i in range(n):
    box(ax, xs[i], card_y0, cw, card_y1 - card_y0, fc="white", ec=LINE, lw=1.2, r=0.1)
    pill(ax, cxs[i], 4.36, titles[i], fc=cols[i], fs=12, h=0.32)
    label(ax, cxs[i], 2.5, descs[i], fs=11, c=TEXT, linespacing=1.35)
    if i < n - 1:
        xr = xs[i] + cw
        arrow(ax, xr + 0.03, 3.45, xr + gap - 0.03, 3.45, color="#9AA1AB", lw=1.6, ms=11)

# ① 합성: ribosome on mRNA, chain coming out
c, cx = cols[0], cxs[0]
rx = cx + 0.12
ax.add_patch(mp.Ellipse((rx, 3.36), 0.74, 0.42, fc=mix(c, 0.25), ec=c, lw=1.2, zorder=3))
ax.add_patch(mp.Ellipse((rx, 2.98), 0.52, 0.18, fc=mix(c, 0.25), ec=c, lw=1.2, zorder=3))
ax.plot([cx - 0.66, cx + 0.66], [3.08, 3.08], color=TEXT2, lw=2.2, zorder=4, solid_capstyle="round")
label(ax, cx - 0.66, 2.9, "mRNA", fs=11, c=TEXT2, ha="left")
bx = np.linspace(rx + 0.03, rx - 0.47, 4)
by = np.linspace(3.5, 3.98, 4)
ax.plot(bx, by, color=c, lw=1.4, zorder=5)
for x, y in zip(bx, by):
    ax.add_patch(mp.Circle((x, y), 0.058, fc=c, ec="none", zorder=6))

# ② 접힘: chain → folded blob, chaperone underneath
c, cx = cols[1], cxs[1]
wx = np.linspace(cx - 0.66, cx - 0.2, 80)
ax.plot(wx, 3.55 + 0.07 * np.sin((wx - wx[0]) * 26), color=c, lw=2.2, zorder=4)
arrow(ax, cx - 0.16, 3.55, cx + 0.0, 3.55, color=TEXT2, lw=1.3, ms=9)
th = np.linspace(0, 2 * np.pi, 200)
rr = 0.23 + 0.04 * np.sin(5 * th)
ax.add_patch(mp.Polygon(np.c_[cx + 0.36 + rr * np.cos(th), 3.55 + rr * np.sin(th)],
                        closed=True, fc=mix(c, 0.3), ec=c, lw=1.6, zorder=4))
ax.add_patch(mp.Arc((cx + 0.2, 3.55), 1.0, 0.8, theta1=205, theta2=335,
                    color=mix(c, 0.6), lw=6, zorder=3))
label(ax, cx + 0.2, 2.98, "샤페론", fs=11, c=shade(c, 0.1), w="bold")

# ③ 수식: protein with tags
c, cx = cols[2], cxs[2]
ax.add_patch(mp.Ellipse((cx, 3.45), 0.66, 0.42, fc=mix(c, 0.3), ec=c, lw=1.6, zorder=4))
for (dx, dy, t) in [(-0.46, 0.4, "P"), (0.46, 0.4, "당"), (-0.46, -0.42, "Ub"), (0.46, -0.42, "Ac")]:
    ax.plot([cx + dx * 0.45, cx + dx], [3.45 + dy * 0.4, 3.45 + dy], color=c, lw=1.4, zorder=3)
    ax.add_patch(mp.Circle((cx + dx, 3.45 + dy), 0.16, fc="white", ec=c, lw=1.6, zorder=5))
    label(ax, cx + dx, 3.45 + dy, t, fs=11, w="bold", c=c, z=6)

# ④ 위치: cell with nucleus, membrane protein, secreted dots
c, cx = cols[3], cxs[3]
ax.add_patch(mp.Circle((cx - 0.08, 3.45), 0.5, fc="white", ec=c, lw=2, zorder=3))
ax.add_patch(mp.Circle((cx - 0.2, 3.36), 0.18, fc=mix(c, 0.25), ec=c, lw=1.4, zorder=4))
label(ax, cx - 0.2, 3.36, "핵", fs=11, w="bold", c=shade(c, 0.2))
box(ax, cx - 0.15, 3.86, 0.14, 0.16, fc=c, ec="none", r=0.03, z=5)
box(ax, cx + 0.33, 3.6, 0.16, 0.13, fc=c, ec="none", r=0.03, z=5)
for (x, y) in [(cx + 0.15, 3.3), (cx + 0.05, 3.1)]:
    ax.add_patch(mp.Circle((x, y), 0.045, fc=c, ec="none", zorder=5))
for (x, y) in [(cx + 0.5, 3.12), (cx + 0.6, 2.98), (cx + 0.46, 2.92)]:
    ax.add_patch(mp.Circle((x, y), 0.045, fc=c, ec="none", zorder=5))

# ⑤ 분해: ubiquitin-tagged chain → proteasome
c, cx = cols[4], cxs[4]
px = np.linspace(cx - 0.63, cx - 0.3, 4)
ax.plot(px, [3.4] * 4, color=c, lw=1.4, zorder=4)
for x in px:
    ax.add_patch(mp.Circle((x, 3.4), 0.055, fc=c, ec="none", zorder=5))
for (x, y) in [(cx - 0.46, 3.58), (cx - 0.46, 3.74), (cx - 0.46, 3.9)]:
    ax.add_patch(mp.Circle((x, y), 0.06, fc="white", ec=c, lw=1.4, zorder=5))
ax.plot([cx - 0.46, cx - 0.46], [3.4, 3.84], color=c, lw=1.2, zorder=4)
arrow(ax, cx - 0.2, 3.45, cx + 0.02, 3.45, color=TEXT2, lw=1.3, ms=9)
for k in range(4):
    box(ax, cx + 0.08, 3.08 + k * 0.19, 0.5, 0.15, fc=mix(c, 0.3), ec=c, lw=1.2, r=0.05, z=4)

# data row
label(ax, 0.02, 1.95, "데이터에서", fs=12, w="bold", c=TEXT2, ha="left")
ry0, rh = 1.0, 0.74
box(ax, xs[0], ry0, cw, rh, fc=FAINT, ec=LINE, r=0.1)
label(ax, cxs[0], ry0 + rh / 2, "RNA-seq이 세는 것\n= mRNA 수", fs=11, linespacing=1.35)
box(ax, xs[1], ry0, cw, rh, fc=FAINT, ec=LINE, r=0.1)
label(ax, cxs[1], ry0 + rh / 2, "HSP 유전자 발현\n(HSPA1A 등)", fs=11, linespacing=1.35)
x3 = xs[2]
box(ax, x3, ry0, xs[4] + cw - x3, rh, fc=ACCENT, ec="none", r=0.1)
label(ax, (x3 + xs[4] + cw) / 2, ry0 + rh / 2, "PTM·위치·분해는 mRNA에 흔적 없음",
      fs=13, w="bold", c=ON_ACCENT)

save(fig, os.path.join(HERE, "../../figures/ch07_lifecycle.png"))
