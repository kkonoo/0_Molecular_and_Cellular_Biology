# ch03 §3 DNA 압축 다섯 단계와 닫힌/열린 염색질  (redraw of figures/ch03_chromatin.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

DNA = CAT[1]          # DNA strand: teal everywhere
NUC = PRIMARY         # nucleosome / histone: blue everywhere
TF = CAT[2]           # transcription factor: green
DNA2 = mix(DNA, 0.5)  # second strand of the helix
RUNG = "#C4C9D1"

W, H = 8.2, 4.35
fig, ax = new(W, H)

# ── top row: five compaction stages ─────────────────────────────────────────
gap = 0.26
pw = (W - 4 * gap) / 5
ty0, th = 2.30, 1.95
cy = 3.24
heads = ["DNA 이중나선", "뉴클레오솜 사슬", "접힌 염색질", "고리 · TAD", "분열기 염색체"]
notes = ["지름 약 2 nm", "약 11 nm", "불규칙하게 뭉침", "Mb 규모", "가장 압축"]

for i in range(5):
    x0 = i * (pw + gap)
    box(ax, x0, ty0, pw, th, fc=FAINT, ec="none", r=0.1, z=1)
    label(ax, x0 + pw / 2, ty0 + th - 0.27, heads[i], fs=12, w="bold")
    label(ax, x0 + pw / 2, ty0 + 0.24, notes[i], fs=11, c=TEXT2)
    if i < 4:
        arrow(ax, x0 + pw + 0.05, cy, x0 + pw + gap - 0.05, cy, ms=11, lw=1.5)

# 1) double helix
cx = 0 * (pw + gap) + pw / 2
xs = np.linspace(cx - 0.58, cx + 0.58, 300)
L = xs[-1] - xs[0]
ph = 2 * np.pi * 2 * (xs - xs[0]) / L
y1 = cy + 0.2 * np.sin(ph)
y2 = cy + 0.2 * np.sin(ph + np.pi)
for xr in np.linspace(cx - 0.54, cx + 0.54, 15):
    p = 2 * np.pi * 2 * (xr - xs[0]) / L
    ax.plot([xr, xr], [cy + 0.2 * np.sin(p), cy + 0.2 * np.sin(p + np.pi)],
            color=RUNG, lw=1.0, zorder=2)
ax.plot(xs, y2, color=DNA2, lw=2.4, zorder=3, solid_capstyle="round")
ax.plot(xs, y1, color=DNA, lw=2.4, zorder=4, solid_capstyle="round")

# 2) nucleosome chain (beads on a string)
cx = 1 * (pw + gap) + pw / 2
xs = np.linspace(cx - 0.64, cx + 0.64, 200)
ax.plot(xs, cy + 0.07 * np.sin(2 * np.pi * (xs - xs[0]) / 0.42), color=DNA, lw=2.0, zorder=2)
for dx, dy in [(-0.42, 0.03), (0.0, -0.03), (0.42, 0.03)]:
    bx, by = cx + dx, cy + dy
    ax.add_patch(mp.Circle((bx, by), 0.16, fc=NUC, ec="none", zorder=3))
    t = np.linspace(-1, 1, 30)
    ax.plot(bx + 0.17 * t, by - 0.07 + 0.14 * (t + 1) / 2 + 0.05 * (1 - t ** 2),
            color=DNA, lw=2.0, zorder=4, solid_capstyle="round")

# 3) folded chromatin: irregular clump of a chain
cx = 2 * (pw + gap) + pw / 2
rng = np.random.default_rng(7)
pts = [np.array([cx - 0.12, cy - 0.08])]
while len(pts) < 24:
    a = rng.uniform(0, 2 * np.pi)
    q = pts[-1] + 0.13 * np.array([np.cos(a), np.sin(a)])
    if np.hypot(*(q - [cx, cy])) < 0.38 and min(np.hypot(*(q - p)) for p in pts) > 0.09:
        pts.append(q)
pts = np.array(pts)
ax.plot(pts[:, 0], pts[:, 1], color=DNA, lw=1.6, zorder=2)
for p in pts:
    ax.add_patch(mp.Circle(p, 0.075, fc=NUC, ec="white", lw=0.6, zorder=3))

# 4) loops / TAD
cx = 3 * (pw + gap) + pw / 2
base = cy - 0.34
ax.plot([cx - 0.64, cx + 0.64], [base, base], color=DNA, lw=2.0, zorder=2)
t = np.linspace(0, 1, 200)
for xb, a, hgt in [(cx - 0.3, 0.24, 0.70), (cx + 0.32, 0.22, 0.62)]:
    ax.plot(xb + a * np.sin(2 * np.pi * t), base + hgt * (1 - np.cos(2 * np.pi * t)) / 2,
            color=DNA, lw=2.0, zorder=2)
    ax.add_patch(mp.Circle((xb, base), 0.06, fc="white", ec=TEXT2, lw=1.4, zorder=3))

# 5) mitotic chromosome (smaller, clear of the note)
cx = 4 * (pw + gap) + pw / 2
for sgn in (1, -1):
    ax.plot([cx - 0.22 * sgn, cx + 0.22 * sgn], [cy - 0.34, cy + 0.38], color=NUC, lw=15,
            solid_capstyle="round", zorder=3)
ax.add_patch(mp.Circle((cx, cy + 0.02), 0.055, fc=STEPS[1], ec="none", zorder=4))

# ── bottom row: closed vs open chromatin ────────────────────────────────────
by0, bh = 0.05, 2.0
bw = (W - 0.2) / 2
# closed
box(ax, 0, by0, bw, bh, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, 0.2, by0 + bh - 0.24, "닫힌 염색질 (heterochromatin)", fs=12, w="bold", ha="left")
sp, r = 0.36, 0.16
rows = [1.38, 1.07, 0.76]
centres = []
for k, yy in enumerate(rows):
    n = 10 if k % 2 == 0 else 9
    x_first = bw / 2 - sp * 4.5 + (sp / 2 if k % 2 else 0)
    row = [(x_first + j * sp, yy) for j in range(n)]
    centres.append(row if k % 2 == 0 else row[::-1])
path = np.array([p for row in centres for p in row])
ax.plot(path[:, 0], path[:, 1] + 0.05, color=DNA, lw=1.6, zorder=2)
for row in centres:
    for p in row:
        ax.add_patch(mp.Circle(p, r, fc=NUC, ec="white", lw=0.8, zorder=3))
label(ax, bw / 2, 0.32, "빽빽함 · 발현 억제", fs=11, c=TEXT2)

# open
ox = W - bw
box(ax, ox, by0, bw, bh, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, ox + 0.2, by0 + bh - 0.24, "열린 염색질 (open chromatin)", fs=12, w="bold", ha="left")
dy = 1.0
ax.plot([ox + 0.18, W - 0.18], [dy, dy], color=DNA, lw=2.2, zorder=3)
for bx in [ox + 0.42, ox + 0.84, ox + 1.26, W - 0.84, W - 0.42]:
    ax.add_patch(mp.Circle((bx, dy), 0.17, fc=NUC, ec="none", zorder=4))
g0, g1 = ox + 1.55, W - 1.13
box(ax, g0, dy - 0.15, g1 - g0, 0.3, fc=ACCENT, ec="none", r=0.06, z=2)
tfx = (g0 + g1) / 2
ax.add_patch(mp.Ellipse((tfx, dy + 0.12), 0.34, 0.26, fc=TF, ec="none", zorder=4))
label(ax, tfx, dy + 0.42, "전사인자", fs=12, w="bold")
label(ax, tfx, dy - 0.36, "뉴클레오솜이 비킨 구간", fs=11)
label(ax, tfx, 0.32, "→ ATAC-seq peak", fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch03_chromatin.png"))
