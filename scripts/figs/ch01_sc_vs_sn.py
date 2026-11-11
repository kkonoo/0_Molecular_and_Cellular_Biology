# ch01 §4 scRNA-seq(세포 통째) vs snRNA-seq(핵만): 잡는 범위와 read 구성  (redraw of figures/ch01_sc_vs_sn.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

EXON, INTRON, MITO = CAT[0], CAT[1], CAT[2]
CAPT = TEXT2

fig, ax = new(7.7, 5.1)


def nucleus(cx, cy, r=0.62):
    ax.add_patch(mp.Circle((cx, cy), r, fc=STEPS[3], ec=STEPS[0], lw=1.6, zorder=3))
    # pre-mRNA: exon blocks joined by intron segments
    for yy, xs in [(cy + 0.15, [-0.38, -0.12, 0.14, 0.38]), (cy - 0.18, [-0.3, -0.02, 0.28])]:
        ax.plot([cx + xs[0], cx + xs[-1]], [yy, yy], color=INTRON, lw=2.2, zorder=4,
                solid_capstyle="butt")
        for dx in xs:
            box(ax, cx + dx - 0.055, yy - 0.055, 0.11, 0.11, fc=EXON, ec="none", r=0.02, z=5)


def mrna(x1, x2, y):
    ax.plot([x1, x2], [y, y], color=EXON, lw=5, solid_capstyle="round", zorder=3)


def mito(cx, cy, ang):
    ax.add_patch(mp.Ellipse((cx, cy), 0.52, 0.27, angle=ang, fc=mix(MITO, 0.22),
                            ec=MITO, lw=1.4, zorder=3))
    t = np.deg2rad(ang)
    dx, dy = 0.1 * np.cos(t), 0.1 * np.sin(t)
    ax.plot([cx - dx, cx + dx], [cy - dy, cy + dy], color=MITO, lw=3.2,
            solid_capstyle="round", zorder=4)


def read_bar(x0, y0, parts, w=3.2, h=0.32):
    label(ax, x0, y0 + h + 0.16, "read가 떨어지는 곳", fs=11, c=TEXT2, ha="left")
    x = x0
    for name, frac, col in parts:
        seg = w * frac
        box(ax, x, y0, seg - 0.03, h, fc=col, ec="none", r=0.03, z=3)
        label(ax, x + (seg - 0.03) / 2, y0 + h / 2, name, fs=11, w="bold", c=on(col))
        x += seg


# ── ① scRNA-seq ───────────────────────────────────────────────
box(ax, 0.1, 2.35, 3.6, 2.65, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, 1.9, 4.74, "scRNA-seq · 세포 통째", fs=13, w="bold")
ax.add_patch(mp.Ellipse((1.9, 3.42), 3.0, 1.9, fc="white", ec=STEPS[0], lw=1.6, zorder=2))
nucleus(2.05, 3.4)
for x1, x2, y in [(0.72, 1.02, 3.5), (0.8, 1.1, 3.0), (1.15, 1.45, 2.72),
                  (2.75, 3.05, 3.85), (2.8, 3.1, 3.12), (2.05, 2.35, 2.6)]:
    mrna(x1, x2, y)
mito(1.15, 4.0, 20)
mito(2.8, 2.82, -25)
box(ax, 0.28, 2.42, 3.24, 2.04, fc="none", ec=CAPT, lw=1.4, r=0.1, z=6, ls=(0, (4, 3)))
label(ax, 3.44, 2.56, "잡는 범위", fs=11, w="bold", c=CAPT, ha="right")

read_bar(0.3, 1.62, [("exon", 0.66, EXON), ("intron", 0.18, INTRON), ("MT-", 0.16, MITO)])
label(ax, 0.3, 1.33, "exon read 위주", fs=11, ha="left")
label(ax, 0.3, 1.03, "MT- read 있음 · 높으면 손상 의심", fs=11, ha="left")

# ── ② snRNA-seq ───────────────────────────────────────────────
box(ax, 3.9, 2.35, 3.7, 2.65, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, 5.75, 4.74, "snRNA-seq · 핵만", fs=13, w="bold")
nucleus(5.2, 3.45)
box(ax, 4.3, 2.48, 1.8, 1.95, fc="none", ec=CAPT, lw=1.4, r=0.1, z=6, ls=(0, (4, 3)))
label(ax, 6.02, 2.62, "잡는 범위", fs=11, w="bold", c=CAPT, ha="right")
label(ax, 6.25, 3.45, "세포질과\n미토콘드리아는\n벗겨 냄", fs=11, c=TEXT2, ha="left",
      linespacing=1.4)

read_bar(4.1, 1.62, [("exon", 0.52, EXON), ("intron", 0.48, INTRON)])
label(ax, 4.1, 1.33, "intron read 많음 (pre-mRNA)", fs=11, ha="left")
label(ax, 4.1, 1.03, "MT- read 거의 없어야 정상", fs=11, ha="left")
pill(ax, 4.1, 0.68, "많으면 세포질이 덜 벗겨진 신호", fc=ACCENT, tc=ON_ACCENT, fs=11,
     ha="left", h=0.32)

save(fig, os.path.join(HERE, "../../figures/ch01_sc_vs_sn.png"))
