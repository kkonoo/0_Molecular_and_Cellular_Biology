# ch09 §2 S.Score–G2M.Score 평면에서 Phase가 나뉘는 영역과 CC.Difference·분열 방향  (redraw of figures/ch09_score_plane.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

# Phase colours: parallel items -> CAT
C_G1, C_S, C_G2M = CAT[0], CAT[1], CAT[2]
LO, HI = -0.4, 0.8

fig, ax = plt.subplots(figsize=(5.6, 5.4))
ax.set_xlim(LO, HI); ax.set_ylim(LO, HI)
ax.set_aspect("equal")

# ── regions (Seurat rule: both < 0 → G1; S > G2M → S; G2M > S → G2M) ──
ax.fill([LO, 0, 0, LO], [LO, LO, 0, 0], color=mix(C_G1, 0.13), lw=0, zorder=0)
ax.fill([0, HI, HI, 0], [LO, LO, HI, 0], color=mix(C_S, 0.13), lw=0, zorder=0)
ax.fill([LO, 0, HI, LO], [0, 0, HI, HI], color=mix(C_G2M, 0.13), lw=0, zorder=0)

# zero lines and the S = G2M diagonal
ax.axhline(0, color="#B8BEC6", lw=0.9, ls=":", zorder=1)
ax.axvline(0, color="#B8BEC6", lw=0.9, ls=":", zorder=1)
ax.plot([0, HI], [0, HI], color=TEXT2, lw=1.2, ls=(0, (5, 3)), zorder=1)
label(ax, 0.665, 0.735, "S = G2M", fs=11, c=TEXT2, rotation=45)

# ── illustrative cells (fixed seed), coloured by the rule itself ──
rng = np.random.default_rng(9)
g1 = np.column_stack([rng.normal(-0.11, 0.055, 70), rng.normal(-0.11, 0.055, 70)])
s = np.column_stack([rng.normal(0.36, 0.13, 45), rng.normal(0.04, 0.05, 45)])
g2m = np.column_stack([rng.normal(0.06, 0.08, 45), rng.normal(0.45, 0.10, 45)])
pts = np.clip(np.vstack([g1, s, g2m]), LO + 0.03, HI - 0.03)
x, y = pts[:, 0], pts[:, 1]
phase = np.where((x < 0) & (y < 0), "G1", np.where(x > y, "S", "G2M"))
for ph, c in [("G1", C_G1), ("S", C_S), ("G2M", C_G2M)]:
    m = phase == ph
    ax.scatter(x[m], y[m], s=26, color=c, edgecolor="white", lw=0.6, zorder=3)

# ── region names + rule as direct labels ──
label(ax, -0.37, 0.77, "G2M", fs=14, w="bold", c=C_G2M, ha="left", va="top")
label(ax, -0.37, 0.685, "G2M.Score가 더 큼", fs=11, c=TEXT, ha="left", va="top")
label(ax, 0.77, -0.24, "S", fs=14, w="bold", c=C_S, ha="right", va="center")
label(ax, 0.77, -0.33, "S.Score가 더 큼", fs=11, c=TEXT, ha="right", va="center")
label(ax, -0.37, -0.255, "G1", fs=14, w="bold", c=C_G1, ha="left", va="center")
label(ax, -0.37, -0.345, "두 점수 모두 < 0", fs=11, c=TEXT, ha="left", va="center")

# ── CC.Difference direction: perpendicular to the diagonal ──
ax.annotate("", xy=(0.66, 0.36), xytext=(0.36, 0.66),
            arrowprops=dict(arrowstyle="<|-|>", color=TEXT, lw=1.6, mutation_scale=14,
                            shrinkA=0, shrinkB=0), zorder=5)
label(ax, 0.66, 0.29, "CC.Difference\n= S − G2M", fs=11, c=TEXT, ha="center", va="top")

# ── dividing direction: along the diagonal (dashed) ──
ax.plot([-0.22, 0.14], [-0.22, 0.14], color=TEXT, lw=1.6, ls=(0, (4, 3)), zorder=5)
ax.annotate("", xy=(0.185, 0.185), xytext=(0.13, 0.13),
            arrowprops=dict(arrowstyle="-|>", color=TEXT, lw=1.6, mutation_scale=14,
                            shrinkA=0, shrinkB=0), zorder=5)
label(ax, -0.04, 0.11, "둘 다 높아짐\n= 분열 중인가", fs=11, c=TEXT, ha="right", va="center")

# ── axes ──
ticks = [-0.4, 0, 0.4, 0.8]
ax.set_xticks(ticks); ax.set_yticks(ticks)
ax.set_xticklabels(["-0.4", "0", "0.4", "0.8"], fontsize=11)
ax.set_yticklabels(["-0.4", "0", "0.4", "0.8"], fontsize=11)
ax.set_xlabel("S.Score", fontsize=12, fontweight="bold", color=C_S)
ax.set_ylabel("G2M.Score", fontsize=12, fontweight="bold", color=C_G2M)

fig.tight_layout()
save(fig, os.path.join(HERE, "../../figures/ch09_score_plane.png"))
