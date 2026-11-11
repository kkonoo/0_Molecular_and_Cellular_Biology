# A2 §2 Ensembl ID를 부분으로 나눠 읽기 + 유전자·전사체·단백질 층마다 다른 ID  (redraw of figures/A2_gene_ids.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

MONO = "DejaVu Sans Mono"
W, H = 7.6, 5.85
fig, ax = new(W, H)

C_GENE, C_TX, C_PROT = CAT[0], CAT[1], CAT[2]   # 유전자 · 전사체 · 단백질 (두 패널 공통)
GREY = "#9AA1AB"

# ── panel backgrounds ───────────────────────────────────────────────────────
box(ax, 0.1, 3.08, W - 0.2, 2.67, fc=FAINT, ec="none", r=0.12, z=0)
box(ax, 0.1, 0.1, W - 0.2, 2.85, fc=FAINT, ec="none", r=0.12, z=0)

# ── ① Ensembl ID 읽기 ──────────────────────────────────────────────────────
label(ax, 0.32, 5.5, "① Ensembl ID 읽기", fs=12.5, w="bold", ha="left")

segs = [  # text, above, below, fill, edge, text colour, above colour
    ("ENS", "Ensembl", "", "white", GREY, TEXT, TEXT),
    ("MUS", "종 코드", "마우스", "white", GREY, TEXT, TEXT),
    ("G", "객체 종류", "유전자", mix(C_GENE, 0.15), C_GENE, shade(C_GENE, 0.25), C_GENE),
    ("00000017167", "번호", "", "white", GREY, TEXT, TEXT),
    (".6", "버전", "6번째", ACCENT, ACCENT, ON_ACCENT, ACCENT_INK),
]
cw, pad, gap = 0.167, 0.26, 0.1
widths = [len(s[0]) * cw + pad for s in segs]
x = (W - (sum(widths) + gap * (len(segs) - 1))) / 2
x_start = x
yb, hb = 4.25, 0.56
for (txt, above, below, fc, ec, tc, ac), w in zip(segs, widths):
    box(ax, x, yb, w, hb, fc=fc, ec=ec, lw=1.6, r=0.07, z=2)
    ax.text(x + w / 2, yb + hb / 2, txt, fontsize=20, fontweight="bold", family=MONO,
            color=tc, ha="center", va="center", zorder=4)
    label(ax, x + w / 2, yb + hb + 0.2, above, fs=12, w="bold", c=ac)
    if below:
        ax.plot([x + w / 2, x + w / 2], [yb - 0.03, yb - 0.15], color=GREY, lw=1, zorder=2)
        label(ax, x + w / 2, yb - 0.3, below, fs=11, c=TEXT2)
    x += w + gap

# 객체 종류 범례 (색 = 아래 패널의 층 색)
yl = 3.62
label(ax, x_start, yl, "객체 종류", fs=11, c=TEXT2, ha="left")
xx = x_start + 0.82
for t, c in [("G 유전자", C_GENE), ("T 전사체", C_TX), ("P 단백질", C_PROT), ("E 엑손", None)]:
    fc = mix(c, 0.15) if c else "white"
    tc = shade(c, 0.25) if c else TEXT2
    w = pill(ax, xx, yl, t, fc=fc, tc=tc, fs=11, h=0.28, ha="left")
    xx += w + 0.1
label(ax, x_start, 3.3, "사람은 종 코드 없음 → ENSG00000141510 (TP53)", fs=11, c=TEXT2, ha="left")

# ── ② 층마다 다른 ID ──────────────────────────────────────────────────────
label(ax, 0.32, 2.7, "② 층마다 다른 ID", fs=12.5, w="bold", ha="left")
label(ax, W - 0.3, 2.7, "예: 사람 TP53", fs=11, c=TEXT2, ha="right")

cols = [(1.7, 1.9, "유전자", C_GENE), (3.95, 1.45, "전사체", C_TX), (5.9, 1.45, "단백질", C_PROT)]
for cx, cwid, name, c in cols:
    label(ax, cx + cwid / 2, 2.38, name, fs=12, w="bold", c=c)
label(ax, cols[1][0] + cols[1][1] / 2, 2.16, "여럿일 수 있음", fs=11, c=TEXT2)

rows = [  # row label, y, gene cell, tx cell, protein cell (text, mono?)
    ("Ensembl", 1.78, ("ENSG00000141510", True), ("ENST…", True), ("ENSP…", True)),
    ("NCBI", 1.14, ("Gene ID 7157", False), ("NM_… · NR_…", True), ("NP_…", True)),
    ("HGNC · UniProt", 0.5, ("HGNC:11998", False), None, ("UniProt accession", False)),
]
ch = 0.42

def cell(cx, cwid, yc, c, content, stack=False):
    if stack:
        for d in (2, 1):
            box(ax, cx + 0.06 * d, yc - ch / 2 + 0.045 * d, cwid, ch,
                fc="white", ec=mix(c, 0.6), lw=1.0, r=0.06, z=2)
    box(ax, cx, yc - ch / 2, cwid, ch, fc=mix(c, 0.12), ec=c, lw=1.4, r=0.06, z=3)
    txt, mono = content
    ax.text(cx + cwid / 2, yc, txt, fontsize=12 if mono else 11.5,
            family=MONO if mono else None, color=shade(c, 0.3),
            ha="center", va="center", zorder=4)

for rl, yc, g, t, p in rows:
    label(ax, 0.35, yc, rl, fs=12, w="bold", ha="left")
    (gx, gw, _, gc), (tx, tw, _, tc), (px, pw, _, pc) = cols
    cell(gx, gw, yc, gc, g)
    cell(px, pw, yc, pc, p)
    if t:
        cell(tx, tw, yc, tc, t, stack=True)
        arrow(ax, gx + gw + 0.04, yc, tx - 0.04, yc, color=GREY, lw=1.4, ms=11)
        arrow(ax, tx + tw + 0.16, yc, px - 0.04, yc, color=GREY, lw=1.4, ms=11)
    else:
        arrow(ax, gx + gw + 0.04, yc, px - 0.04, yc, color=GREY, lw=1.3, ms=11, ls="--")

save(fig, os.path.join(HERE, "../../figures/A2_gene_ids.png"))
