# ch10 §4 같은 가상의 세포를 해상도 낮음·중간·높음으로 클러스터링한 결과 비교  (redraw of figures/ch10_resolution.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

rng = np.random.default_rng(7)

# ── 가상의 UMAP 좌표 (패널 안 0..1 단위), 세 패널에 같은 점을 씀 ─────────────────
def blob(cx, cy, sx, sy, n):
    return np.column_stack([rng.normal(cx, sx, n), rng.normal(cy, sy, n)])

def streak(cx, cy, dx, dy, n, noise):
    t = rng.uniform(-1, 1, n)
    return np.column_stack([cx + t * dx + rng.normal(0, noise, n),
                            cy + t * dy + rng.normal(0, noise, n)])

groups = {
    "cd4n":   streak(0.28, 0.68, 0.075, 0.10, 70, 0.025),
    "cd4m":   blob(0.43, 0.80, 0.040, 0.035, 40),
    "prolif": blob(0.52, 0.94, 0.018, 0.015, 12),
    "cd8":    blob(0.13, 0.45, 0.035, 0.035, 35),
    "nk":     blob(0.14, 0.27, 0.035, 0.030, 25),
    "b":      blob(0.82, 0.82, 0.045, 0.040, 40),
    "cd14":   blob(0.56, 0.26, 0.060, 0.040, 60),
    "fcg":    blob(0.82, 0.16, 0.035, 0.030, 25),
    "dc":     blob(0.91, 0.48, 0.014, 0.014, 10),
    "dbl":    blob(0.46, 0.60, 0.014, 0.014, 8),
}
for k in groups:
    groups[k] = np.clip(groups[k], 0.02, 0.98)

# ── 클러스터 색: 같은 집단은 세 패널에서 같은 색 계열 ──────────────────────────
GREY_B = "#7D8590"
COL = {
    "cd4n": PRIMARY, "cd4m": mix(PRIMARY, 0.45), "prolif": shade(PRIMARY, 0.5),
    "cd8": CAT[1], "nk": mix(CAT[1], 0.45),
    "b": GREY_B,
    "cd14": CAT[2], "fcg": mix(CAT[2], 0.45), "dc": shade(CAT[2], 0.45),
    "dbl": "#BFC5CD",
}
low = {k: PRIMARY for k in ["cd4n", "cd4m", "prolif", "cd8", "nk"]}
low.update({"b": GREY_B})
low.update({k: CAT[2] for k in ["cd14", "fcg", "dc", "dbl"]})

mid = {k: PRIMARY for k in ["cd4n", "cd4m", "prolif"]}
mid.update({"cd8": CAT[1], "nk": CAT[1], "b": GREY_B,
            "cd14": CAT[2], "dbl": CAT[2], "fcg": COL["fcg"], "dc": COL["dc"]})

high = dict(COL)

# 라벨: (텍스트, u, v, ha)
lab_low = [("T·NK", 0.33, 0.50, "center"), ("B", 0.69, 0.82, "right"),
           ("골수계", 0.52, 0.06, "center")]
lab_mid = [("CD4 T", 0.33, 0.50, "center"), ("CD8 T·NK", 0.17, 0.08, "center"),
           ("B", 0.69, 0.82, "right"), ("CD14 단핵구", 0.56, 0.06, "center"),
           ("FCGR3A", 0.82, 0.33, "center"), ("DC", 0.85, 0.48, "right")]
lab_high = [("증식 상태", 0.44, 0.95, "right"), ("doublet", 0.51, 0.60, "left"),
            ("DC (희귀)", 0.85, 0.46, "right")]

panels = [("해상도 낮음 · 3개", low, lab_low),
          ("중간 · 6개", mid, lab_mid),
          ("높음 · 10개", high, lab_high)]

PW, PH, GAP, M = 2.65, 2.45, 0.15, 0.10
W = 3 * PW + 2 * GAP
fig, ax = new(W, PH + 0.45)

for i, (head, cmap, labs) in enumerate(panels):
    x0, y0 = i * (PW + GAP), 0.0
    box(ax, x0, y0, PW, PH, fc=FAINT, ec=LINE, lw=1.0, r=0.10, z=1)
    label(ax, x0 + PW / 2, PH + 0.22, head, fs=12, w="bold")

    def to_in(u, v):
        return x0 + M + u * (PW - 2 * M), y0 + M + v * (PH - 2 * M)

    for k, pts in groups.items():
        xs, ys = to_in(pts[:, 0], pts[:, 1])
        ax.scatter(xs, ys, s=11, c=cmap[k], edgecolors="white", linewidths=0.3, zorder=3)

    for text, u, v, ha in labs:
        x, y = to_in(u, v)
        label(ax, x, y, text, fs=11, w="bold", ha=ha)

save(fig, os.path.join(HERE, "../../figures/ch10_resolution.png"))
