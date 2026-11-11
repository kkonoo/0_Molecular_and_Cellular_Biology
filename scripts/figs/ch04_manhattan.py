# ch04 §4 Manhattan plot: towers on chr 2, 6, 11 cross 5×10⁻⁸, chr 15 stays below  (redraw of figures/ch04_manhattan.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

rng = np.random.default_rng(4)

# illustrative chromosome lengths (Mb), only used for relative widths
L = np.array([248, 242, 198, 190, 181, 171, 159, 145, 138, 134, 135, 133,
              114, 107, 102, 90, 83, 80, 59, 64, 47, 51], dtype=float)
GAP = 14.0
starts = np.concatenate([[0], np.cumsum(L + GAP)[:-1]])
centers = starts + L / 2

# towers: chromosome -> (relative position, peak height, number of LD neighbours)
TOWERS = {2: (0.35, 11.0, 45), 6: (0.18, 9.0, 35), 11: (0.62, 8.0, 35), 15: (0.55, 6.2, 35)}

N_BG = 2860
n_per = np.round(N_BG * L / L.sum()).astype(int)

xs, ys, cs = [], [], []
col = [PRIMARY, CAT[3]]
lead = None
for i in range(22):
    chrom = i + 1
    n = n_per[i]
    pos = rng.uniform(0, L[i], n)
    y = -np.log10(rng.uniform(0, 1, n))
    big = y > 4.2
    y[big] = rng.uniform(0, 2, big.sum())          # keep background below ~4
    xs.append(starts[i] + pos); ys.append(y); cs += [col[i % 2]] * n
    if chrom in TOWERS:
        rel, peak, k = TOWERS[chrom]
        c = rel * L[i]
        tp = c + rng.normal(0, 1.2, k)
        ty = peak * rng.uniform(0, 1, k) ** 1.4
        ty = np.minimum(ty, peak * 0.93)
        tp = np.append(tp, c); ty = np.append(ty, peak)   # lead SNP at the peak
        xs.append(starts[i] + tp); ys.append(ty); cs += [col[i % 2]] * (k + 1)
        if chrom == 2:
            lead = (starts[i] + c, peak)

x = np.concatenate(xs); y = np.concatenate(ys)

fig, ax = plt.subplots(figsize=(7.2, 3.7))
ax.scatter(x, y, s=7, c=cs, linewidths=0, zorder=3)

xmax = starts[-1] + L[-1] + 6
ax.set_xlim(-6, xmax)
ax.set_ylim(0, 12.6)
ax.set_yticks(range(0, 13, 2))
ax.yaxis.grid(True, color=LINE, lw=0.8, zorder=0)
ax.set_axisbelow(True)

# thresholds
gw = -np.log10(5e-8)
ax.axhline(gw, color=ACCENT_INK, lw=1.8, ls=(0, (5, 3)), zorder=4)
ax.text(xmax - 8, gw + 0.25, "5×10⁻⁸ (−log₁₀ = 7.3)", ha="right", va="bottom",
        fontsize=12, fontweight="bold", color=ACCENT_INK, zorder=5)
ax.axhline(5, color=TEXT2, lw=1.0, ls=(0, (2, 2.5)), zorder=4)
ax.text(xmax - 8, 5.2, "1×10⁻⁵ (참고선)", ha="right", va="bottom",
        fontsize=11, color=TEXT2, zorder=5)

# lead SNP of the chromosome-2 tower
lx, ly = lead
ax.scatter([lx], [ly], s=70, facecolors="none", edgecolors=TEXT, linewidths=1.6, zorder=6)
ax.annotate("lead SNP", xy=(lx + 9, ly), xytext=(lx + 150, ly - 0.3),
            fontsize=12, fontweight="bold", color=TEXT, va="center",
            arrowprops=dict(arrowstyle="-", color=TEXT2, lw=1.0), zorder=6)

# chromosome ticks (skip the narrowest ones to avoid overlap)
show = [c for c in range(1, 23) if c not in (17, 19, 21)]
ax.set_xticks([centers[c - 1] for c in show])
ax.set_xticklabels([str(c) for c in show], fontsize=11)
ax.tick_params(axis="x", length=0, pad=5)
ax.tick_params(axis="y", labelsize=11)
ax.set_xlabel("염색체 (각 염색체 안에서는 위치 순서)", fontsize=12, labelpad=8)
ax.set_ylabel("−log₁₀(P)", fontsize=12)

fig.tight_layout()
save(fig, os.path.join(HERE, "../../figures/ch04_manhattan.png"))
