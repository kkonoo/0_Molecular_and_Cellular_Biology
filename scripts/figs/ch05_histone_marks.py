# ch05 §4.3 한 유전체 구간에서 히스톤 표지 6가지가 주로 나타나는 자리 (genome-browser 형식)  (new figure)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.4, 4.95
fig, ax = new(W, H)

X0, X1 = 1.45, 8.3                 # track area
LAB_X = 1.3                        # right edge of track labels
TRACK_H, STEP = 0.42, 0.58
GENE_Y = 3.85
REG_Y = 4.58

# ── positions of the regions (inches) ────────────────────────────────────────
ENH_A, ENH_S, TSS = 1.9, 2.98, 4.05
GENE_END = 5.85
K27_A, K27_B = 6.1, 7.2
K9_A, K9_B = 7.38, 8.1

# faint vertical guides behind the narrow regions
for a, b in [(ENH_A - 0.27, ENH_A + 0.27), (ENH_S - 0.27, ENH_S + 0.27), (TSS - 0.22, TSS + 0.27)]:
    ax.add_patch(mp.Rectangle((a, 0.1), b - a, 4.0, fc=FAINT, ec="none", zorder=0))

# ── region labels on top ─────────────────────────────────────────────────────
label(ax, ENH_A, REG_Y, "활성\nenhancer", fs=12, w="bold", linespacing=1.15)
pill(ax, ENH_S, REG_Y, "대기\nenhancer", fc=ACCENT, tc=ON_ACCENT, fs=12, h=0.56, padx=0.14)
label(ax, TSS, REG_Y, "promoter\n(TSS)", fs=12, w="bold", linespacing=1.15)
label(ax, 5.12, REG_Y, "발현 유전자", fs=12, w="bold")
label(ax, (K27_A + K27_B) / 2, REG_Y, "억제 구간\n(Polycomb)", fs=12, w="bold", linespacing=1.15)
label(ax, (K9_A + K9_B) / 2, REG_Y, "이질염색질\n반복 서열", fs=12, w="bold", linespacing=1.15)

# ── gene track ───────────────────────────────────────────────────────────────
label(ax, LAB_X, GENE_Y, "유전자", fs=12, w="bold", ha="right")
ax.plot([X0, X1], [GENE_Y, GENE_Y], color=TEXT2, lw=1.0, zorder=2)
for a, b in [(TSS, 4.25), (4.55, 4.68), (4.95, 5.08), (5.35, 5.5), (5.7, GENE_END)]:
    ax.add_patch(mp.Rectangle((a, GENE_Y - 0.1), b - a, 0.2, fc=TEXT, ec="none", zorder=3))
# TSS bent arrow
ax.plot([TSS, TSS], [GENE_Y + 0.1, GENE_Y + 0.27], color=TEXT, lw=1.4, zorder=3)
arrow(ax, TSS, GENE_Y + 0.27, TSS + 0.26, GENE_Y + 0.27, color=TEXT, lw=1.4, ms=10)
# enhancer elements (open boxes)
for c in (ENH_A, ENH_S):
    ax.add_patch(mp.Rectangle((c - 0.11, GENE_Y - 0.07), 0.22, 0.14, fc="white", ec=TEXT2,
                              lw=1.2, zorder=3))
# repeats
for xr in np.arange(K9_A + 0.05, K9_B - 0.05, 0.12):
    ax.add_patch(mp.Rectangle((xr, GENE_Y - 0.07), 0.07, 0.14, fc="#9AA1AB", ec="none", zorder=3))

# ── signal shapes (illustrative) ─────────────────────────────────────────────
x = np.linspace(X0, X1, 1600)
def g(mu, sd):
    return np.exp(-0.5 * ((x - mu) / sd) ** 2)
def block(a, b, w=0.03):
    return 1 / (1 + np.exp(-(x - a) / w)) - 1 / (1 + np.exp(-(x - b) / w))
rng = np.random.default_rng(3)
def noise(amp=0.035):
    n = np.convolve(rng.random(x.size), np.ones(25) / 25, mode="same")
    return amp * n

tracks = [
    ("H3K4me3",  CAT[0],    0.9 * (g(TSS - 0.07, 0.07) + g(TSS + 0.13, 0.08))),
    ("H3K4me1",  CAT[3],    0.75 * (g(ENH_A - 0.1, 0.08) + g(ENH_A + 0.1, 0.08)
                                    + g(ENH_S - 0.1, 0.08) + g(ENH_S + 0.1, 0.08))),
    ("H3K27ac",  CAT[1],    0.95 * g(ENH_A, 0.09) + 0.9 * g(TSS + 0.04, 0.09)),
    ("H3K36me3", CAT[2],    0.55 * block(4.3, GENE_END, 0.06) * (1 + 0.08 * np.sin(x * 13))),
    ("H3K27me3", "#6B7280", 0.65 * block(K27_A, K27_B, 0.04) * (1 + 0.1 * np.sin(x * 9))),
    ("H3K9me3",  "#4B5563", 0.65 * block(K9_A, K9_B, 0.04) * (1 + 0.1 * np.sin(x * 11 + 1))),
]

for i, (name, col, s) in enumerate(tracks):
    yb = 0.15 + (5 - i) * STEP
    s = np.clip(s + noise(), 0, 1)
    ax.plot([X0, X1], [yb, yb], color=LINE, lw=1.0, zorder=1)
    ax.fill_between(x, yb, yb + s * TRACK_H, color=col, lw=0, zorder=2)
    label(ax, LAB_X, yb + 0.16, name, fs=12, w="bold", ha="right")
    if name == "H3K27ac":
        label(ax, ENH_S, yb + 0.14, "없음", fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch05_histone_marks.png"))
