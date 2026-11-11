# ch04 §3-2 Punnett 표: 상염색체 열성(Aa×Aa)과 X연관 열성(보인자 어머니×발현 없는 아버지)  (redraw of figures/ch04_inheritance.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.4, 3.6
fig, ax = new(W, H)

CW, CH = 1.05, 0.78          # cell size
GAP = 0.04
FILL = {"none": "white", "carrier": STEPS[2], "affected": PRIMARY}


def X(s):
    """Superscript genotype in the regular text font."""
    return r"$\mathregular{" + s + "}$"


def cell(x, y, geno, note, kind):
    fc = FILL[kind]
    ec = LINE if kind == "none" else "none"
    box(ax, x + GAP, y + GAP, CW - 2 * GAP, CH - 2 * GAP, fc=fc, ec=ec, lw=1.0, r=0.06)
    tc = on(fc)
    if note:
        label(ax, x + CW / 2, y + CH * 0.62, geno, fs=15, w="bold", c=tc)
        label(ax, x + CW / 2, y + CH * 0.27, note, fs=11, c=tc if kind == "affected" else TEXT2)
    else:
        label(ax, x + CW / 2, y + CH / 2, geno, fs=15, w="bold", c=tc)


def panel(x0, marker, cross, top_who, left_who, cols, rows, cells, bottom):
    gx = x0 + 1.05                 # grid left
    gtop = 2.2                     # grid top
    label(ax, x0, 3.42, marker, fs=13, w="bold", ha="left")
    label(ax, x0, 3.08, cross, fs=11, c=TEXT2, ha="left")
    # column header
    label(ax, gx + CW, 2.7, top_who, fs=11, c=TEXT2)
    for j, a in enumerate(cols):
        label(ax, gx + CW * (j + 0.5), 2.42, a, fs=14, w="bold")
    # row header (horizontal)
    label(ax, x0, gtop - CH, left_who, fs=11, c=TEXT2, ha="left")
    for i, a in enumerate(rows):
        label(ax, gx - 0.25, gtop - CH * (i + 0.5), a, fs=14, w="bold")
    # cells
    for i in range(2):
        for j in range(2):
            g, n, k = cells[i][j]
            cell(gx + CW * j, gtop - CH * (i + 1), g, n, k)
    label(ax, x0, 0.36, bottom, fs=12, ha="left")


# ① 상염색체 열성
panel(0.2, "① 상염색체 열성", "보인자 부모  Aa × Aa", "부모 2", "부모 1",
      ["A", "a"], ["A", "a"],
      [[("AA", "", "none"), ("Aa", "보인자", "carrier")],
       [("Aa", "보인자", "carrier"), ("aa", "발현", "affected")]],
      "aa만 발현 → 자녀의 1/4")

# ② X연관 열성
panel(3.95, "② X연관 열성", "보인자 어머니 " + X("X^AX^a") + " × 아버지 " + X("X^AY"),
      "아버지", "어머니",
      [X("X^A"), "Y"], [X("X^A"), X("X^a")],
      [[(X("X^AX^A"), "딸", "none"), (X("X^AY"), "아들", "none")],
       [(X("X^AX^a"), "딸 · 보인자", "carrier"), (X("X^aY"), "아들 · 발현", "affected")]],
      X("X^a") + "를 받은 아들만 발현")

save(fig, os.path.join(HERE, "../../figures/ch04_inheritance.png"))
