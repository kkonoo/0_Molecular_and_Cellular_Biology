# ch04 §2-1 레퍼런스와 비교한 SNV·삽입·결실(VCF 표기)과 구간 중복·결실(CNV)  (redraw of figures/ch04_variant_types.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.3, 5.75
fig, ax = new(W, H)

X0 = 1.35          # first sequence cell left edge
CW = 0.30          # cell width
BW, BH = 0.26, 0.32
XL = 0.1           # row-label left edge
XPOS, XREF, XALT = 6.15, 6.85, 7.65

ANCHOR = CAT[1]                 # teal: VCF가 함께 적는 앞 글자
CHG_FC, CHG_EC = STEPS[2], PRIMARY   # 바뀐/더해진/빠진 글자

def cx(i):
    return X0 + CW * i + CW / 2

def seq_row(y, cells):
    """cells: list of (letter, kind) with kind in {'', 'chg', 'del', 'anc'}"""
    for i, (ch, kind) in enumerate(cells):
        x = cx(i)
        c = TEXT
        if kind == "chg":
            box(ax, x - BW / 2, y - BH / 2, BW, BH, fc=CHG_FC, ec=CHG_EC, lw=1.3, r=0.04)
            c = PRIMARY
        elif kind == "del":
            box(ax, x - BW / 2, y - BH / 2, BW, BH, fc="white", ec=CHG_EC, lw=1.2, r=0.04, ls=(0, (2, 1.5)))
            c = PRIMARY
        elif kind == "anc":
            box(ax, x - BW / 2, y - BH / 2, BW, BH, fc=mix(ANCHOR, 0.12), ec=ANCHOR, lw=1.6, r=0.04)
            c = ANCHOR
        label(ax, x, y, ch, fs=13, w="bold", c=c)

def row_label(y, main, sub=None):
    if sub:
        label(ax, XL, y + 0.09, main, fs=12, w="bold", ha="left")
        label(ax, XL, y - 0.15, sub, fs=11, c=TEXT2, ha="left")
    else:
        label(ax, XL, y, main, fs=12, w="bold", ha="left")

def vcf(y, pos, ref, alt):
    for x, t in [(XPOS, pos), (XREF, ref), (XALT, alt)]:
        label(ax, x, y, t, fs=13, w="bold")

ref = "CTGACCTAGGCA"

# ── ① 짧은 변이 ─────────────────────────────────────────────
label(ax, XL, 5.55, "① 짧은 변이 (50 bp 미만)", fs=12, w="bold", ha="left")
label(ax, XPOS - 0.25, 5.55, "VCF", fs=12, w="bold", ha="left")

yh = 5.15
label(ax, X0 - 0.1, yh, "위치", fs=11, c=TEXT2, ha="right")
for i in range(12):
    label(ax, cx(i), yh, str(i + 1), fs=11, c=TEXT2)
for x, t in [(XPOS, "POS"), (XREF, "REF"), (XALT, "ALT")]:
    label(ax, x, yh, t, fs=11, c=TEXT2)

y = 4.72
row_label(y, "레퍼런스")
seq_row(y, [(c, "") for c in ref])

ax.plot([XL, W - 0.1], [4.38, 4.38], color=LINE, lw=1.0)

# SNV: 위치 5 C→T
y = 4.0
row_label(y, "SNV", "한 글자 치환")
cells = [(c, "") for c in ref]
cells[4] = ("T", "chg")
seq_row(y, cells)
vcf(y, "5", "C", "T")

# 삽입: 위치 6 C 뒤에 GA
y = 3.38
row_label(y, "삽입", "insertion")
cells = [(c, "") for c in ref[:5]] + [("C", "anc"), ("G", "chg"), ("A", "chg")] + [(c, "") for c in ref[6:]]
seq_row(y, cells)
vcf(y, "6", "C", "CGA")

# 결실: 위치 7–8 TA 빠짐
y = 2.76
row_label(y, "결실", "deletion")
cells = [(c, "") for c in ref[:5]] + [("C", "anc"), ("–", "del"), ("–", "del")] + [(c, "") for c in ref[8:]]
seq_row(y, cells)
vcf(y, "6", "CTA", "C")

# 앞 글자 표시 설명
yl = 2.3
box(ax, X0 + 0.02, yl - 0.11, 0.22, 0.22, fc=mix(ANCHOR, 0.12), ec=ANCHOR, lw=1.6, r=0.04)
label(ax, X0 + 0.36, yl, "VCF가 함께 적는 앞 글자 (POS = 이 글자의 위치)", fs=11, c=ANCHOR, ha="left")

ax.plot([XL, W - 0.1], [1.98, 1.98], color=LINE, lw=1.0)

# ── ② 구조변이 ──────────────────────────────────────────────
label(ax, XL, 1.7, "② 구조변이 (SV, 50 bp 이상)", fs=12, w="bold", ha="left")

SW, SH, SG = 0.6, 0.32, 0.66
def seg_row(y, names):
    for i, n in enumerate(names):
        x = X0 + SG * i
        fc = PRIMARY if n == "B" else STEPS[2]
        chevron(ax, x, y - SH / 2, SW, SH, fc, tip=0.12, first=True)
        label(ax, x + (SW - 0.12) / 2 + 0.02, y, n, fs=12, w="bold", c=on(fc))

XA = X0 + SG * 5 + 0.35
y = 1.25
row_label(y, "레퍼런스")
seg_row(y, ["A", "B", "C", "D"])
label(ax, XA, y, "중복·결실 = 복사 수 변이 (CNV)", fs=12, w="bold", ha="left")

y = 0.75
row_label(y, "중복")
seg_row(y, ["A", "B", "B", "C", "D"])
label(ax, XA, y, "B 2개 → 복사 수 증가", fs=11, c=TEXT2, ha="left")

y = 0.25
row_label(y, "결실")
seg_row(y, ["A", "C", "D"])
label(ax, XA, y, "B 없음 → 복사 수 감소", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch04_variant_types.png"))
