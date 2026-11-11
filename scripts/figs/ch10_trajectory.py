# ch10 §4 연속 분화의 궤적: pseudotime으로 색칠한 세포와 pseudotime에 따른 marker 발현  (redraw of figures/ch10_trajectory.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

rng = np.random.default_rng(10)

# ── 가상의 세포: 줄기(trunk) → 갈림점 → 계통 A / 계통 B ───────────────────────
T_BP = 0.45                                   # 갈림점의 pseudotime (설명용)
root = np.array([0.0, 0.0])
bp = np.array([2.0, 0.0])
endA = np.array([4.0, 1.35])
endB = np.array([4.0, -1.2])

def along(p0, p1, s):
    return p0[None, :] + s[:, None] * (p1 - p0)[None, :]

n_tr, n_br = 70, 55
s_tr = np.sort(rng.uniform(0, 1, n_tr))
s_a = np.sort(rng.uniform(0, 1, n_br))
s_b = np.sort(rng.uniform(0, 1, n_br))
xy = np.vstack([along(root, bp, s_tr), along(bp, endA, s_a), along(bp, endB, s_b)])
xy += rng.normal(0, 0.11, xy.shape)
pt = np.concatenate([s_tr * T_BP, T_BP + s_a * (1 - T_BP), T_BP + s_b * (1 - T_BP)])

cmap = LinearSegmentedColormap.from_list("pt", [STEPS[1], PRIMARY, shade(PRIMARY, 0.3)])

COL_PRO = CAT[3]          # 전구 marker
COL_A = CAT[1]            # 계통 A (왼쪽 라벨과 오른쪽 곡선에 같은 색)
COL_B = CAT[2]            # 계통 B

fig = plt.figure(figsize=(8.2, 3.6))

# ── ① 왼쪽: UMAP 모식도 ──────────────────────────────────────────────────────
ax1 = fig.add_axes([0.0, 0.06, 0.47, 0.82])
ax1.set_xlim(-0.6, 4.6); ax1.set_ylim(-1.95, 2.0); ax1.axis("off")

ax1.scatter([root[0]], [root[1]], s=820, c=ACCENT, ec="none", zorder=1)   # 시작점 (단 하나의 강조)
ax1.scatter(xy[:, 0], xy[:, 1], c=pt, cmap=cmap, vmin=0, vmax=1, s=30,
            ec="white", lw=0.5, zorder=3)
ax1.scatter([bp[0]], [bp[1]], s=320, facecolors="none", edgecolors=TEXT2,
            lw=1.2, linestyles="--", zorder=4)

ax1.text(0.0, -0.58, "시작점 (root)", ha="center", va="center", fontsize=12,
         fontweight="bold", color=TEXT, zorder=5)
ax1.text(0.0, -0.9, "사전 지식으로 지정", ha="center", va="center", fontsize=11,
         color=TEXT2, zorder=5)
ax1.text(1.5, 0.62, "갈림점", ha="center", va="center", fontsize=12,
         fontweight="bold", color=TEXT, zorder=5)
ax1.plot([1.72, 1.9], [0.48, 0.2], color=TEXT2, lw=1.0, zorder=4)
ax1.text(endA[0] - 0.05, endA[1] + 0.45, "계통 A", ha="center", va="center",
         fontsize=12, fontweight="bold", color=COL_A, zorder=5)
ax1.text(endB[0] - 0.05, endB[1] - 0.48, "계통 B", ha="center", va="center",
         fontsize=12, fontweight="bold", color=COL_B, zorder=5)

# pseudotime 척도: 단색 칸으로 나눈 막대
nb, x0, x1, yb, hb = 8, -0.45, 1.45, -1.72, 0.16
wseg = (x1 - x0) / nb
for i in range(nb):
    ax1.add_patch(mp.Rectangle((x0 + i * wseg, yb), wseg, hb,
                               fc=to_hex(cmap((i + 0.5) / nb)), ec="white", lw=0.6, zorder=3))
ax1.text(x0, yb + hb + 0.17, "pseudotime", ha="left", va="center", fontsize=11, color=TEXT)
ax1.text(x0, yb - 0.17, "0", ha="center", va="center", fontsize=11, color=TEXT2)
ax1.text(x1, yb - 0.17, "1", ha="center", va="center", fontsize=11, color=TEXT2)

fig.text(0.01, 0.95, "① 세포를 길 위에 늘어세움", fontsize=12, fontweight="bold", color=TEXT)

# ── ② 오른쪽: pseudotime에 따른 발현 ─────────────────────────────────────────
def sig(z):
    return 1 / (1 + np.exp(-z))

t = np.linspace(0, 1, 300)
y_pro = 0.04 + 0.94 * (1 - sig((t - T_BP) / 0.06))
y_a = 0.06 + 0.90 * sig((t - 0.68) / 0.06)
y_b = 0.03 + 0.72 * sig((t - 0.75) / 0.06)

ax2 = fig.add_axes([0.565, 0.2, 0.30, 0.66])
ax2.set_xlim(0, 1); ax2.set_ylim(-0.04, 1.16)
ax2.axvline(T_BP, color="#9AA1AB", lw=1.0, ls="--", zorder=1)
ax2.text(T_BP + 0.025, 1.1, "갈림점", ha="left", va="center", fontsize=11, color=TEXT2)
ax2.plot(t, y_pro, color=COL_PRO, lw=2.6, zorder=3)
ax2.plot(t, y_a, color=COL_A, lw=2.6, zorder=3)
ax2.plot(t, y_b, color=COL_B, lw=2.6, ls=(0, (4, 2.2)), zorder=3)

ax2.text(1.03, y_a[-1], "계통 A marker", ha="left", va="center", fontsize=11,
         fontweight="bold", color=COL_A, clip_on=False)
ax2.text(1.03, y_b[-1], "계통 B marker", ha="left", va="center", fontsize=11,
         fontweight="bold", color=COL_B, clip_on=False)
ax2.text(1.03, y_pro[-1] + 0.02, "전구 marker", ha="left", va="center", fontsize=11,
         fontweight="bold", color=shade(COL_PRO, 0.35), clip_on=False)

ax2.set_xticks([0, 1]); ax2.set_xticklabels(["0", "1"], fontsize=11)
ax2.set_yticks([])
ax2.set_xlabel("pseudotime (실제 시간 아님)", fontsize=11, labelpad=4)
ax2.set_ylabel("발현", fontsize=11, rotation=0, labelpad=16, va="center")

fig.text(0.53, 0.95, "② pseudotime에 따라 바뀌는 유전자", fontsize=12, fontweight="bold", color=TEXT)

save(fig, os.path.join(HERE, "../../figures/ch10_trajectory.png"))
