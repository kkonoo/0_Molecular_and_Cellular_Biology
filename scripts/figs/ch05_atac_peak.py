# ch05 §data-box ATAC-seq: 열린 구간 → fragment 쌓임 → coverage peak(chromStart, chromEnd, summit)  (redraw of figures/ch05_atac_peak.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.1, 5.2
fig, ax = new(W, H)

X0, X1 = 1.7, 7.0           # track span
OS, OE = 3.75, 5.05         # open region
NUC_C, NUC_E = STEPS[2], STEPS[0]
TF_C = CAT[1]
TN5_C = CAT[2]
READ_C = PRIMARY

# faint band marking the open region through all rows
ax.add_patch(mp.Rectangle((OS, 0.85), OE - OS, 3.9, fc=FAINT, ec="none", zorder=0))

# ── row 1: chromatin ─────────────────────────────────────────────────────────
yl = 4.3
ax.plot([X0, X1], [yl, yl], color="#9AA1AB", lw=1.6, zorder=1)
left_c = [1.9 + 0.33 * i for i in range(6)]
right_c = [5.2 + 0.33 * i for i in range(6)]
for cx in left_c + right_c:
    ax.add_patch(mp.Circle((cx, yl), 0.13, fc=NUC_C, ec=NUC_E, lw=1.0, zorder=2))
ax.add_patch(mp.Ellipse((4.4, yl), 0.3, 0.2, fc=TF_C, ec="none", zorder=3))
label(ax, 4.4, 4.03, "TF", fs=11, w="bold", c=TF_C)
label(ax, (left_c[0] + left_c[-1]) / 2, 4.03, "뉴클레오솜", fs=11, c=TEXT2)

tn5 = list(np.linspace(3.84, 4.24, 4)) + list(np.linspace(4.56, 4.98, 5))
tn5 += [left_c[1] + 0.165, left_c[4] + 0.165, right_c[1] + 0.165, right_c[4] + 0.165]
ax.plot(tn5, [yl + 0.24] * len(tn5), ls="none", marker="v", ms=7, color=TN5_C, zorder=4)

# region brackets
def bracket(x1, x2, text, bold=False):
    yb = 4.75
    ax.plot([x1, x1, x2, x2], [yb - 0.07, yb, yb, yb - 0.07], color=TEXT2, lw=1.0, zorder=2)
    label(ax, (x1 + x2) / 2, yb + 0.17, text, fs=12, w="bold" if bold else "regular",
          c=TEXT if bold else TEXT2)
bracket(1.78, 3.67, "닫힌 구간")
bracket(OS + 0.04, OE - 0.04, "열린 구간", bold=True)
bracket(5.13, 6.95, "닫힌 구간")

# ── row 2: fragments ─────────────────────────────────────────────────────────
rng = np.random.default_rng(3)
frags = []
for _ in range(20):
    c = np.clip(rng.normal(4.4, 0.27), 3.8, 5.0)
    L = rng.uniform(0.12, 0.32)
    frags.append((c - L / 2, c + L / 2))
frags += [(1.95, 2.35), (2.95, 3.2), (5.35, 5.75), (6.3, 6.55), (6.6, 6.9), (2.45, 2.7)]
frags.sort()
rows_end = []
placed = []
for s, e in frags:
    for i, re_ in enumerate(rows_end):
        if s > re_ + 0.05:
            rows_end[i] = e
            placed.append((s, e, i))
            break
    else:
        rows_end.append(e)
        placed.append((s, e, len(rows_end) - 1))
nrows = len(rows_end)
ytop = 3.6
dy = min(0.1, 0.9 / max(nrows - 1, 1))      # stack stays between 3.6 and 2.7
for s, e, i in placed:
    ax.plot([s, e], [ytop - i * dy] * 2, color=READ_C, lw=3, solid_capstyle="butt", zorder=3)

# ── row 3: coverage computed from the fragments above ───────────────────────
yb, hmax = 0.85, 1.3
xg = np.arange(X0, X1 + 1e-9, 0.005)
cov = np.zeros_like(xg)
for s, e, _ in placed:
    cov += (xg >= s) & (xg <= e)
k = np.exp(-0.5 * (np.arange(-24, 25) / 8.0) ** 2)
cov = np.convolve(cov, k / k.sum(), mode="same")
y = yb + cov / cov.max() * hmax
im = int(np.argmax(cov))
thr = 0.3 * cov.max()
a = im
while a > 0 and cov[a - 1] >= thr:
    a -= 1
b = im
while b < len(cov) - 1 and cov[b + 1] >= thr:
    b += 1
xs, xe, xm = xg[a], xg[b], xg[im]

ax.fill_between(xg, yb, y, color=mix(PRIMARY, 0.18), lw=0, zorder=2)
sel = (xg >= xs) & (xg <= xe)
ax.fill_between(xg[sel], yb, y[sel], color=mix(PRIMARY, 0.5), lw=0, zorder=2)
ax.plot(xg, y, color=PRIMARY, lw=1.4, zorder=3)
ax.plot([X0, X1], [yb, yb], color=TEXT2, lw=1.2, zorder=3)

# summit
ax.plot([xm, xm], [yb, y[im] + 0.05], color=TEXT, lw=1.2, ls=(0, (3, 2)), zorder=4)
label(ax, xm, y[im] + 0.2, "summit", fs=12, w="bold")
# chromStart / chromEnd
for xx in (xs, xe):
    ax.plot([xx, xx], [yb - 0.18, yb + 0.02], color=TEXT, lw=1.2, zorder=4)
label(ax, xs + 0.03, yb - 0.3, "chromStart", fs=11, ha="right")
label(ax, xe - 0.03, yb - 0.3, "chromEnd", fs=11, ha="left")
label(ax, xe + 0.12, yb + 0.6, "peak", fs=12, w="bold", c=PRIMARY, ha="left")
label(ax, X1, 0.2, "유전체 좌표 →", fs=11, c=TEXT2, ha="right")

# ── left labels ─────────────────────────────────────────────────────────────
label(ax, 0.05, 4.38, "① 염색질", fs=13, w="bold", ha="left")
ax.plot([0.14], [4.08], ls="none", marker="v", ms=7, color=TN5_C)
label(ax, 0.28, 4.08, "Tn5 자리", fs=11, c=TN5_C, ha="left")
label(ax, 0.05, 3.3, "② 읽힌 fragment", fs=13, w="bold", ha="left")
label(ax, 0.05, 3.02, "막대 = fragment", fs=11, c=TEXT2, ha="left")
label(ax, 0.05, 1.55, "③ 쌓인 수", fs=13, w="bold", ha="left")
label(ax, 0.05, 1.27, "(coverage)", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch05_atac_peak.png"))
