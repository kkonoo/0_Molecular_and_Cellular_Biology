# ch05 §TF 모티프 PFM 개수 표와 같은 표로 그린 sequence logo  (redraw of figures/ch05_motif_logo.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Affine2D

# ── data: chapter R code (cnt), rows = A,C,G,T, columns = positions 1..6 ──────────
BASES = ["A", "C", "G", "T"]
cols = [[10, 3, 4, 3], [1, 0, 18, 1], [19, 0, 0, 1],
        [2, 1, 0, 17], [8, 2, 8, 2], [6, 2, 6, 6]]
cnt = np.array(cols, dtype=float).T            # 4 × 6
p = cnt / cnt.sum(axis=0)                      # column frequencies (sum 20)
with np.errstate(divide="ignore", invalid="ignore"):
    plogp = np.where(p > 0, p * np.log2(p), 0.0)
ic = 2 - (-plogp.sum(axis=0))                  # bits, max log2(4) = 2
hts = p * ic                                   # letter height = freq × IC

BC = {b: CAT[i] for i, b in enumerate(BASES)}  # four bases = parallel items
HILITE = mix(TEXT2, 0.16)                      # neutral tint: "most common base" cells

fig, ax = new(8.3, 3.6)

# ── ① count table ────────────────────────────────────────────────────────────
label(ax, 0.15, 3.35, "① 위치별 개수 (PFM)", fs=13, w="bold", ha="left")
cw, ch = 0.5, 0.45
tx0 = 0.65
row_y = [2.45, 2.0, 1.55, 1.1]
label(ax, 0.40, 2.9, "위치", fs=11, c=TEXT2)
for j in range(6):
    label(ax, tx0 + cw * (j + 0.5), 2.9, str(j + 1), fs=12, w="bold")
for i, b in enumerate(BASES):
    y = row_y[i]
    label(ax, 0.40, y, b, fs=13, w="bold", c=BC[b])
    for j in range(6):
        v = int(cnt[i, j])
        top = v == cnt[:, j].max()
        fc = HILITE if top else "white"
        ax.add_patch(mp.Rectangle((tx0 + cw * j, y - ch / 2), cw, ch,
                                  fc=fc, ec=LINE, lw=1.0, zorder=2))
        label(ax, tx0 + cw * (j + 0.5), y, str(v), fs=12,
              w="bold" if top else "regular", c=TEXT if top else TEXT2)
label(ax, 0.15, 0.58, "열마다 합 = 20 (서열 수)", fs=11, c=TEXT2, ha="left")
label(ax, 0.15, 0.30, "색칠한 칸 = 그 위치에서 가장 흔한 염기", fs=11, c=TEXT2, ha="left")

# ── ② sequence logo ──────────────────────────────────────────────────────────
label(ax, 4.30, 3.35, "② sequence logo", fs=13, w="bold", ha="left")
base_y = 1.1 - ch / 2          # align logo baseline with table bottom
unit = 1.0                     # inches per bit
axx = 4.78
lx0, lw, pitch = 4.92, 0.48, 0.56

ax.plot([axx, axx], [base_y, base_y + 2 * unit], color="#9AA1AB", lw=1.0, zorder=2)
ax.plot([axx, lx0 + 5 * pitch + lw + 0.04], [base_y, base_y], color="#9AA1AB", lw=1.0, zorder=2)
for t in (0, 1, 2):
    yy = base_y + t * unit
    ax.plot([axx - 0.06, axx], [yy, yy], color="#9AA1AB", lw=1.0, zorder=2)
    label(ax, axx - 0.11, yy, str(t), fs=11, c=TEXT2, ha="right")
    if t > 0:
        ax.plot([axx, lx0 + 5 * pitch + lw + 0.04], [yy, yy], color=LINE, lw=0.8, zorder=1)
label(ax, 4.40, base_y + unit, "정보량 (bits)", fs=11, c=TEXT2, rotation=90)

fp = FontProperties(family=["Pretendard", "DejaVu Sans"], weight="bold")
glyphs = {}
for b in BASES:
    tp = TextPath((0, 0), b, size=1, prop=fp)
    glyphs[b] = (tp, tp.get_extents())

def letter(b, x, y, w, h):
    tp, ext = glyphs[b]
    tr = (Affine2D().translate(-ext.x0, -ext.y0)
          .scale(w / ext.width, h / ext.height).translate(x, y))
    ax.add_patch(mp.PathPatch(tr.transform_path(tp), fc=BC[b], ec="none", zorder=3))

for j in range(6):
    x = lx0 + pitch * j
    y = base_y
    for i in np.argsort(hts[:, j]):            # smallest at bottom, largest on top
        h = hts[i, j] * unit
        if h < 0.004:
            continue
        letter(BASES[i], x, y, lw, h)
        y += h
    label(ax, x + lw / 2, base_y - 0.24, str(j + 1), fs=12, w="bold")
    label(ax, x + lw / 2, base_y - 0.52, f"{ic[j]:.2f}", fs=11, c=TEXT2)
label(ax, axx - 0.11, base_y - 0.24, "위치", fs=11, c=TEXT2, ha="right")
label(ax, axx - 0.11, base_y - 0.52, "bits", fs=11, c=TEXT2, ha="right")

save(fig, os.path.join(HERE, "../../figures/ch05_motif_logo.png"))
