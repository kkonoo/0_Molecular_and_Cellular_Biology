# ch03 §3 DNA 압축 다섯 단계와 닫힌·열린 염색질  (redraw of figures/ch03_chromatin.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

DNA = CAT[1]          # DNA strand: teal everywhere
HIS = PRIMARY         # nucleosome / chromatin: blue everywhere
TF = CAT[2]           # transcription factor: green

W, H = 8.4, 4.95
fig, ax = new(W, H)


def nucleosome(x, y, r=0.16, wrap=True):
    ax.add_patch(mp.Circle((x, y), r, fc=HIS, ec="white", lw=0.8, zorder=4))
    if wrap:
        ax.add_patch(mp.Ellipse((x, y), 2.5 * r, 0.9 * r, angle=22, fc="none",
                                ec=DNA, lw=1.6, zorder=5))


# ── top row: five compaction stages ──────────────────────────────────────────
pw, gap, x0 = 1.48, 0.22, 0.06
yb, yt = 2.8, 4.85
yc = 3.78
names = ["DNA 이중나선", "뉴클레오솜 사슬", "접힌 염색질", "고리 · TAD", "분열기 염색체"]
sizes = ["지름 약 2 nm", "약 11 nm", "5–24 nm 굵기", "Mb 규모", "가장 압축"]

for i in range(5):
    px = x0 + i * (pw + gap)
    box(ax, px, yb, pw, yt - yb, fc=FAINT, ec="none", r=0.1, z=1)
    label(ax, px + pw / 2, yt - 0.25, names[i], fs=12, w="bold")
    label(ax, px + pw / 2, yb + 0.24, sizes[i], fs=11, c=TEXT2)
    if i < 4:
        arrow(ax, px + pw + 0.04, yc, px + pw + gap - 0.04, yc, ms=10, lw=1.4)

# 1 double helix
px = x0
xs = np.linspace(px + 0.14, px + pw - 0.14, 300)
ph = 2 * np.pi * (xs - xs[0]) / 0.6
y1 = yc + 0.22 * np.sin(ph)
y2 = yc - 0.22 * np.sin(ph)
for xr in np.arange(xs[0] + 0.03, xs[-1], 0.06):
    a = 0.22 * np.sin(2 * np.pi * (xr - xs[0]) / 0.6)
    ax.plot([xr, xr], [yc - a, yc + a], color=mix(DNA, 0.35), lw=1, zorder=2)
ax.plot(xs, y1, color=DNA, lw=2, zorder=3)
ax.plot(xs, y2, color=DNA, lw=2, zorder=3)

# 2 nucleosome chain (beads on a string)
px = x0 + (pw + gap)
xs = np.linspace(px + 0.1, px + pw - 0.1, 200)
ys = yc + 0.1 * np.sin(2 * np.pi * (xs - xs[0]) / 0.9)
ax.plot(xs, ys, color=DNA, lw=1.8, zorder=3)
for xn in [px + 0.32, px + 0.74, px + 1.16]:
    nucleosome(xn, np.interp(xn, xs, ys))

# 3 folded chromatin: irregular clumped chain, centred in the panel
px = x0 + 2 * (pw + gap)
cx = px + pw / 2
rng = np.random.default_rng(4)
pts = [np.array([0.0, 0.0])]
for _ in range(27):
    for _try in range(40):
        ang = rng.uniform(0, 2 * np.pi)
        p = pts[-1] + 0.16 * np.array([np.cos(ang), np.sin(ang)])
        if np.hypot(p[0], p[1]) < 0.5:
            break
    pts.append(p)
pts = np.array(pts)
pts -= (pts.max(0) + pts.min(0)) / 2
pts += np.array([cx, yc])
ax.plot(pts[:, 0], pts[:, 1], color=DNA, lw=1.5, zorder=3)
for p in pts:
    ax.add_patch(mp.Circle(p, 0.085, fc=HIS, ec="white", lw=0.7, zorder=4))

# 4 loops / TAD
px = x0 + 3 * (pw + gap)
ybase = yc - 0.4
ax.plot([px + 0.1, px + pw - 0.1], [ybase, ybase], color=DNA, lw=1.8, zorder=3)
s = np.linspace(0, 2 * np.pi, 200)
for xm, a, b in [(px + 0.45, 0.26, 0.4), (px + 1.05, 0.23, 0.33)]:
    ax.plot(xm - a * np.sin(s), ybase + b * (1 - np.cos(s)), color=DNA, lw=1.8, zorder=3)

# 5 mitotic chromosome
px = x0 + 4 * (pw + gap)
cx = px + pw / 2
for sgn in (1, -1):
    ax.plot([cx - 0.24 * sgn, cx + 0.24 * sgn], [yc - 0.46, yc + 0.46], color=HIS, lw=15,
            solid_capstyle="round", zorder=3)

# ── bottom row: closed vs open chromatin ────────────────────────────────────
bw = (5 * pw + 4 * gap - gap) / 2
by0, by1 = 0.45, 2.55
ynote = 0.76

# closed (heterochromatin)
px = x0
box(ax, px, by0, bw, by1 - by0, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, px + 0.2, by1 - 0.28, "닫힌 염색질 (heterochromatin)", fs=12.5, w="bold", ha="left")
sp, r = 0.32, 0.15
rows = [(1.80, 10, 0.0), (1.52, 9, sp / 2), (1.24, 10, 0.0)]
path = []
for k, (yy, n, off) in enumerate(rows):
    xl = px + bw / 2 - (10 - 1) * sp / 2 + off
    xr = [xl + j * sp for j in range(n)]
    if k % 2:
        xr = xr[::-1]
    for j, xx in enumerate(xr):
        path.append((xx, yy + (0.07 if j % 2 else -0.07)))
        ax.add_patch(mp.Circle((xx, yy), r, fc=HIS, ec="white", lw=0.8, zorder=4))
path = np.array(path)
ax.plot(path[:, 0], path[:, 1], color=DNA, lw=1.5, zorder=5)
label(ax, px + 0.2, ynote, "빽빽함 · 발현 억제", fs=11, c=TEXT2, ha="left")

# open chromatin
px = x0 + bw + gap
box(ax, px, by0, bw, by1 - by0, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, px + 0.2, by1 - 0.28, "열린 염색질 (open chromatin)", fs=12.5, w="bold", ha="left")
yl = 1.45
ax.plot([px + 0.2, px + bw - 0.2], [yl, yl], color=DNA, lw=1.8, zorder=3)
xg = px + 2.1
box(ax, xg - 0.48, yl - 0.11, 0.96, 0.22, fc=ACCENT, ec="none", r=0.06, z=2)
for xn in [px + 0.45, px + 0.9, px + 1.35, px + 2.85, px + 3.3]:
    nucleosome(xn, yl)
box(ax, xg - 0.21, yl + 0.01, 0.42, 0.3, fc=TF, ec="none", r=0.13, z=6)
label(ax, xg, yl + 0.52, "전사인자", fs=12, w="bold")
label(ax, xg, yl - 0.36, "뉴클레오솜이 비킨 짧은 구간", fs=11)
label(ax, px + 0.2, ynote, "ATAC-seq이 이 구간을 읽어 peak로 보여 줌", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch03_chromatin.png"))
