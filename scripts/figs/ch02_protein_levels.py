# ch02 §4-2 단백질 구조의 네 단계(1차~4차)와 도메인  (redraw of figures/ch02_protein_levels.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.6, 4.45
fig, ax = new(W, H)

GREY = "#9AA1AB"
HELIX, STRAND = CAT[0], CAT[1]          # α나선 = 파랑, β병풍 = 청록 (그림 전체에서 같은 뜻)
DOM1, DOM2 = CAT[2], CAT[3]             # 도메인 1·2 (나란한 항목)


def smooth(pts, n=24):
    """Catmull-Rom spline through pts."""
    P = np.array(pts, float)
    P = np.vstack([P[0], P, P[-1]])
    out = []
    for i in range(1, len(P) - 2):
        p0, p1, p2, p3 = P[i - 1], P[i], P[i + 1], P[i + 2]
        for t in np.linspace(0, 1, n, endpoint=False):
            out.append(0.5 * (2 * p1 + (-p0 + p2) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t * t
                              + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(P[-2])
    return np.array(out)


def strand(x1, x2, y, w=0.09, hw=0.2, hl=0.13, z=4):
    ax.add_patch(mp.FancyArrow(x1, y, x2 - x1, 0, width=w, head_width=hw, head_length=hl,
                               length_includes_head=True, fc=STRAND, ec="none", zorder=z))


# ── 위: 네 단계 패널 ─────────────────────────────────────────────────────────
pw, gap, py0, py1 = 1.66, 0.32, 2.2, 4.4
heads = ["1차 · 서열", "2차 · 국소 구조", "3차 · 접힘", "4차 · 조립"]
xs = [i * (pw + gap) for i in range(4)]
for i, x in enumerate(xs):
    box(ax, x, py0, pw, py1 - py0, fc=FAINT, ec="none", r=0.1, z=1)
    label(ax, x + pw / 2, 4.12, heads[i], fs=12.5, w="bold")
    if i < 3:
        arrow(ax, x + pw + 0.06, 3.1, x + pw + gap - 0.06, 3.1, color=GREY, lw=1.5, ms=11)

# 1차: 아미노산 구슬 사슬 (N → C)
cx = xs[0] + pw / 2 + 0.09
sp, rb = 0.34, 0.13
cols = [cx + k * sp for k in (-1.5, -0.5, 0.5, 1.5)]
rows = [3.55, 3.1, 2.65]
beads = [(cols[0], rows[0], "M"), (cols[1], rows[0], "K"), (cols[2], rows[0], "T"), (cols[3], rows[0], "A"),
         (cols[3], rows[1], "Y"), (cols[2], rows[1], "I"), (cols[1], rows[1], "A"), (cols[0], rows[1], "K"),
         (cols[0], rows[2], "Q"), (cols[1], rows[2], "R")]
bx = [b[0] for b in beads]; by = [b[1] for b in beads]
ax.plot(bx, by, color=GREY, lw=1.6, zorder=2, solid_capstyle="round")
for x, y, ch in beads:
    ax.add_patch(mp.Circle((x, y), rb, fc="white", ec=TEXT2, lw=1.2, zorder=3))
    label(ax, x, y, ch, fs=11, w="bold")
label(ax, cols[0] - 0.25, rows[0], "N", fs=12, w="bold")
label(ax, cols[1] + 0.26, rows[2], "C", fs=12, w="bold")

# 2차: α나선 + β병풍
x0 = xs[1] + 0.14
t = np.linspace(0, 10 * np.pi, 500)
hx = x0 + 0.06 + 0.027 * t + 0.06 * np.sin(t)
hy = 3.6 + 0.16 * np.cos(t)
ax.plot(hx, hy, color=HELIX, lw=2.2, zorder=3, solid_capstyle="round")
label(ax, x0 + 0.98, 3.6, "α나선", fs=11.5, w="bold", c=HELIX, ha="left")
ys = [2.95, 2.68, 2.41]
for k, y in enumerate(ys):
    if k % 2 == 0:
        strand(x0, x0 + 0.88, y)
    else:
        strand(x0 + 0.88, x0, y)
for y_a, y_b in zip(ys[:-1], ys[1:]):
    for xx in (x0 + 0.2, x0 + 0.42, x0 + 0.64):
        ax.plot([xx, xx], [y_a - 0.055, y_b + 0.055], color=GREY, lw=1.1, ls=(0, (1, 1.5)), zorder=2)
label(ax, x0 + 0.98, 2.68, "β병풍", fs=11.5, w="bold", c=STRAND, ha="left")

# 3차: 사슬 하나가 접힌 모양
c3x, c3y = xs[2] + pw / 2 + 0.01, 3.1
th = np.linspace(0, 2 * np.pi, 200)
rr = 1 + 0.05 * np.sin(3 * th + 0.6)
ax.add_patch(mp.Polygon(np.c_[c3x + 0.75 * rr * np.cos(th), c3y + 0.72 * rr * np.sin(th)],
                        closed=True, fc="white", ec=GREY, lw=1.2, ls=(0, (3, 2)), zorder=1.5))
for seg in ([(4.18, 3.25), (4.2, 3.45), (4.35, 3.55)],
            [(5.05, 3.55), (5.3, 3.52), (5.42, 3.25), (5.4, 2.85), (5.2, 2.7)],
            [(4.4, 2.7), (4.22, 2.75), (4.2, 2.92), (4.42, 3.0)]):
    s = smooth(seg)
    ax.plot(s[:, 0], s[:, 1], color=GREY, lw=1.8, zorder=2, solid_capstyle="round")
box(ax, 4.33, 3.45, 0.74, 0.2, fc=HELIX, ec="none", r=0.1, z=4)
strand(5.2, 4.4, 2.7)
strand(4.42, 5.12, 3.0)

# 4차: 접힌 사슬 4개가 모임
c4x, c4y = xs[3] + pw / 2, 3.1
for dx in (-0.31, 0.31):
    for dy in (0.31, -0.31):
        x, y = c4x + dx, c4y + dy
        ax.add_patch(mp.Circle((x, y), 0.3, fc="white", ec=GREY, lw=1.2, ls=(0, (3, 2)), zorder=2))
        box(ax, x - 0.14, y + 0.04, 0.28, 0.09, fc=HELIX, ec="none", r=0.045, z=4)
        ax.add_patch(mp.FancyArrow(x - 0.14, y - 0.1, 0.28, 0, width=0.05, head_width=0.12,
                                   head_length=0.08, length_includes_head=True,
                                   fc=STRAND, ec="none", zorder=4))

# ── 아래: 도메인 ─────────────────────────────────────────────────────────────
pb0, pb1 = 0.5, 1.88
box(ax, 0, pb0, W, pb1 - pb0, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, 0.22, 1.62, "도메인", fs=12.5, w="bold", ha="left")

by0 = 1.02
ax.plot([0.5, 4.55], [by0, by0], color=GREY, lw=6, solid_capstyle="round", zorder=2)
label(ax, 0.28, by0, "N", fs=12, w="bold")
label(ax, 4.78, by0, "C", fs=12, w="bold")
box(ax, 1.0, by0 - 0.19, 1.35, 0.38, fc=DOM1, ec="none", r=0.08, z=3)
label(ax, 1.675, by0, "도메인 1", fs=11.5, w="bold", c=on(DOM1))
box(ax, 2.95, by0 - 0.19, 1.2, 0.38, fc=DOM2, ec="none", r=0.08, z=3)
label(ax, 3.55, by0, "도메인 2", fs=11.5, w="bold", c=on(DOM2))

arrow(ax, 5.0, by0, 5.45, by0, color=GREY, lw=1.5, ms=11)
label(ax, 5.22, by0 + 0.22, "접히면", fs=11, c=TEXT2)

ax.add_patch(mp.Ellipse((6.02, by0), 0.9, 0.78, fc=DOM1, ec="none", zorder=3))
label(ax, 6.02, by0, "도메인 1", fs=11.5, w="bold", c=on(DOM1))
lx = np.linspace(6.45, 6.76, 60)
ly = by0 - 0.02 + 0.06 * np.sin((lx - 6.45) / 0.31 * 3 * np.pi)
ax.plot(lx, ly, color=GREY, lw=1.8, zorder=2)
ax.add_patch(mp.Ellipse((7.13, by0 - 0.04), 0.82, 0.7, fc=DOM2, ec="none", zorder=3))
label(ax, 7.13, by0 - 0.04, "도메인 2", fs=11.5, w="bold", c=on(DOM2))

save(fig, os.path.join(HERE, "../../figures/ch02_protein_levels.png"))
