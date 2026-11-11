# ch12 §3-2 Illumina 합성하며 읽기: cluster, 사이클별 사진, read, 한 사이클의 네 단계  (redraw of figures/ch12_sbs.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.35, 4.0)

BASE = {"A": CAT[0], "C": CAT[1], "G": CAT[2], "T": CAT[3]}
reads = ["ACGTG", "GGTAC", "TACCA"]          # cluster 1, 2, 3 (설명용)
NEUTRAL = TEXT2
PHOTO = "#2A2F38"

head_y = 3.45
row_y = [2.97, 2.575, 2.18]

# ── 왼쪽: flow cell 위 cluster ───────────────────────────────────────────
box(ax, 0.1, 1.62, 1.95, 2.05, fc=FAINT, ec=LINE, r=0.1, z=1)
label(ax, 1.075, head_y, "cluster", fs=12, w="bold")
ax.add_patch(mp.Rectangle((0.25, 2.05), 1.65, 0.07, fc="#B8BDC5", ec="none", zorder=2))
label(ax, 1.075, 1.85, "flow cell 표면", fs=11, c=TEXT2)
rng = np.random.default_rng(3)
for k, cx in enumerate([0.55, 1.075, 1.6]):
    for dx in np.linspace(-0.17, 0.17, 6):
        hgt = rng.uniform(0.5, 0.78)
        ax.plot([cx + dx, cx + dx], [2.12, 2.12 + hgt], color=NEUTRAL, lw=1.6,
                solid_capstyle="round", zorder=3)
    num(ax, cx, 3.13, k + 1, fc=NEUTRAL, r=0.13)

# ── 가운데: 사이클별 사진 ────────────────────────────────────────────────
px0, pw, pg = 2.6, 0.72, 0.12
for i, y in enumerate(row_y):
    num(ax, 2.36, y, i + 1, fc=NEUTRAL, r=0.13)
for c in range(5):
    x = px0 + c * (pw + pg)
    box(ax, x, 1.92, pw, 1.3, fc=PHOTO, ec="none", r=0.06, z=2)
    label(ax, x + pw / 2, head_y, f"사이클 {c + 1}", fs=11, w="bold")
    for i, y in enumerate(row_y):
        b = reads[i][c]
        ax.add_patch(mp.Circle((x + pw / 2, y), 0.13, fc=BASE[b], ec="none", zorder=4))
        label(ax, x + pw / 2, y, b, fs=11, w="bold", c=on(BASE[b]), z=6)
px1 = px0 + 5 * pw + 4 * pg
label(ax, (px0 + px1) / 2, 1.68, "read 길이 = 사이클 수", fs=11, c=TEXT2)

# ── 오른쪽: cluster별 read ───────────────────────────────────────────────
arrow(ax, px1 + 0.1, row_y[1], px1 + 0.4, row_y[1], color="#9AA1AB", lw=1.6)
rx0, tw, ts = 7.15, 0.2, 0.225
label(ax, rx0 + (4 * ts + tw) / 2, head_y, "read", fs=12, w="bold")
for i, y in enumerate(row_y):
    for j, b in enumerate(reads[i]):
        x = rx0 + j * ts
        box(ax, x, y - 0.14, tw, 0.28, fc=BASE[b], ec="none", r=0.04, z=3)
        label(ax, x + tw / 2, y, b, fs=11, w="bold", c=on(BASE[b]))
pill(ax, 7.25, 1.68, "cluster 하나 = read 하나", fc=ACCENT, tc=ON_ACCENT, fs=11)

# ── 아래: 한 사이클의 네 단계 ────────────────────────────────────────────
label(ax, 0.1, 1.22, "한 사이클", fs=12, w="bold", ha="left")
steps = ["형광 염기 1개 붙음", "사진 찍기", "형광·막힌 끝 떼기", "다음 사이클"]
sx0, sx1, gap, tip = 0.1, 8.23, 0.06, 0.18
sw = (sx1 - sx0 - 3 * gap) / 4
for k, s in enumerate(steps):
    x = sx0 + k * (sw + gap)
    chevron(ax, x, 0.5, sw, 0.5, fc=STEPS[3], tip=tip, first=(k == 0))
    nx = x + (0.24 if k == 0 else tip + 0.2)
    num(ax, nx, 0.75, k + 1, fc=NEUTRAL, r=0.13)
    label(ax, nx + 0.22, 0.75, s, fs=11, w="bold", c=TEXT, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch12_sbs.png"))
