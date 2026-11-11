# ch08 §2 KEGG pathway 기호(노드 모양·선 종류)와 MAPK ERK 줄기 예  (redraw of figures/ch08_kegg_legend.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.2, 5.8
fig, ax = new(W, H)

LC = "#4B5563"                       # relation lines
GFILL, GEDGE = STEPS[2], PRIMARY     # gene product (rectangle) — same everywhere
CFILL, CEDGE = mix(CAT[1], 0.15), CAT[1]   # compound
MFILL, MEDGE = mix(CAT[2], 0.15), CAT[2]   # map link


def gene_box(x, yc, w, h, text, fs=12):
    box(ax, x, yc - h / 2, w, h, fc=GFILL, ec=GEDGE, lw=1.3, r=0.02, z=3)
    label(ax, x + w / 2, yc, text, fs=fs, w="bold")


def edge(x1, x2, y, kind, mod=None):
    if kind == "act":
        arrow(ax, x1, y, x2, y, color=LC, lw=1.4, ms=11)
    elif kind == "inh":
        ax.plot([x1, x2 - 0.03], [y, y], color=LC, lw=1.4, zorder=3)
        ax.plot([x2 - 0.03] * 2, [y - 0.1, y + 0.1], color=LC, lw=2.6,
                solid_capstyle="butt", zorder=3)
    elif kind == "ind":
        ax.plot([x1, x2 - 0.1], [y, y], color=LC, lw=1.4, ls=(0, (3, 2)), zorder=3)
        arrow(ax, x2 - 0.12, y, x2, y, color=LC, lw=1.4, ms=11)
    elif kind == "bind":
        ax.plot([x1, x2], [y, y], color=LC, lw=1.4, zorder=3)
    if mod:
        label(ax, (x1 + x2) / 2, y + 0.15, mod, fs=11, w="bold")


# ── ① 노드 모양 ───────────────────────────────────────────────────────────
ytop, ybot = 5.4, 2.35
label(ax, 0.1, 5.62, "① 노드 모양", fs=12, w="bold", ha="left")
box(ax, 0.1, ybot, 3.45, ytop - ybot, fc=FAINT, ec="none", r=0.08, z=1)

sx, sw, sh = 0.25, 0.9, 0.36
tx = sx + sw + 0.15
rows = [4.95, 4.23, 3.51, 2.79]

# gene product
yc = rows[0]
gene_box(sx, yc, sw, sh, "EGFR")
label(ax, tx, yc + 0.11, "유전자 산물", fs=12, ha="left")
label(ax, tx, yc - 0.14, "주로 단백질", fs=11, c=TEXT2, ha="left")

# one box = several genes (the one must-see point)
yc = rows[1]
gene_box(sx, yc, sw, sh, "RAS")
pill(ax, tx, yc + 0.12, "상자 하나 = 유전자 여러 개", fc=ACCENT, tc=ON_ACCENT,
     fs=11, padx=0.1, h=0.28, ha="left")
label(ax, tx, yc - 0.17, "HRAS, KRAS, NRAS 등", fs=11, c=TEXT2, ha="left")

# compound
yc = rows[2]
ax.add_patch(mp.Circle((sx + sw / 2, yc), 0.12, fc=CFILL, ec=CEDGE, lw=1.3, zorder=3))
label(ax, tx, yc + 0.11, "화합물", fs=12, ha="left")
label(ax, tx, yc - 0.14, "2차 전달자, 대사물질 등", fs=11, c=TEXT2, ha="left")

# map link
yc = rows[3]
box(ax, sx, yc - sh / 2, sw, sh, fc=MFILL, ec=MEDGE, lw=1.3, r=sh / 2, z=3)
label(ax, sx + sw / 2, yc, "다른 지도", fs=11, w="bold")
label(ax, tx, yc, "다른 경로 지도로 연결", fs=12, ha="left")

# ── ② 선 (관계) ───────────────────────────────────────────────────────────
px = 3.75
label(ax, px, 5.62, "② 선 (관계)", fs=12, w="bold", ha="left")
box(ax, px, ybot, 8.1 - px, ytop - ybot, fc=FAINT, ec="none", r=0.08, z=1)

ax_, bw, bh, gap = 3.88, 0.36, 0.3, 0.75
bx = ax_ + bw + gap
kx, ex = 5.45, 6.38
edges = [
    ("act", None, "활성화", "activation"),
    ("inh", None, "억제", "inhibition"),
    ("ind", None, "간접 효과", "indirect effect"),
    ("bind", None, "결합", "binding/association"),
    ("act", "+p", "인산화", "phosphorylation"),
    ("inh", "−p", "탈인산화", "dephosphorylation"),
    ("act", "+u", "유비퀴틴화", "ubiquitination"),
]
for i, (kind, mod, ko, en) in enumerate(edges):
    y = ytop - 0.27 - i * 0.42
    gene_box(ax_, y, bw, bh, "A", fs=11)
    gene_box(bx, y, bw, bh, "B", fs=11)
    edge(ax_ + bw, bx, y, kind, mod)
    label(ax, kx, y, ko, fs=12, ha="left")
    label(ax, ex, y, en, fs=11, c=TEXT2, ha="left")

# ── ③ 예: MAPK 지도의 ERK 줄기 ────────────────────────────────────────────
label(ax, 0.1, 2.02, "③ 예: MAPK 지도(hsa04010)의 ERK 줄기", fs=12, w="bold", ha="left")
names = ["EGFR", "GRB2", "SOS", "RAS", "RAF", "MEK", "ERK", "ELK1"]
links = [("bind", None), ("act", None), ("act", None), ("act", None),
         ("act", "+p"), ("act", "+p"), ("act", "+p")]
nw, step, x0, yr = 0.72, 1.04, 0.1, 1.45
for i, n in enumerate(names):
    gene_box(x0 + i * step, yr, nw, sh, n)
for i, (kind, mod) in enumerate(links):
    edge(x0 + i * step + nw, x0 + (i + 1) * step, yr, kind, mod)

# DUSP: negative feedback by dephosphorylation of ERK
ex_c = x0 + 6 * step + nw / 2
yd = 0.5
gene_box(ex_c - nw / 2, yd, nw, sh, "DUSP")
ytopd, yerk = yd + sh / 2, yr - sh / 2
ax.plot([ex_c, ex_c], [ytopd, yerk - 0.03], color=LC, lw=1.4, zorder=3)
ax.plot([ex_c - 0.1, ex_c + 0.1], [yerk - 0.03] * 2, color=LC, lw=2.6,
        solid_capstyle="butt", zorder=3)
label(ax, ex_c + 0.1, (ytopd + yerk) / 2, "−p", fs=11, w="bold", ha="left")
label(ax, ex_c - nw / 2 - 0.15, yd, "음성 되먹임", fs=12, c=TEXT2, ha="right")

save(fig, os.path.join(HERE, "../../figures/ch08_kegg_legend.png"))
