# ch06 §4 full-length vs 3' 방식: 같은 유전자에서 read coverage가 쌓이는 모양  (redraw of figures/ch06_coverage.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.7, 4.85
fig, ax = new(W, H)

# ── geometry (inches) ───────────────────────────────────────────────────────
X0, X1 = 1.9, 7.30                       # gene line
EXONS = [(1.98, 2.61), (3.30, 3.83), (4.54, 5.07), (5.76, 7.22)]
GENE_Y = 3.70                            # exon box centre
EH = 0.24
FL_BASE, SP_BASE = 2.40, 0.95            # track baselines
C_FL, C_3P = CAT[1], CAT[2]              # parallel methods: teal / green

# ── exon 2 band (the one must-see point) ────────────────────────────────────
bx0, bx1 = 3.15, 3.99
box(ax, bx0, SP_BASE - 0.20, bx1 - bx0, 4.30 - (SP_BASE - 0.20),
    fc=mix(ACCENT, 0.30), ec="none", r=0.06, z=0)
pill(ax, (bx0 + bx1) / 2, 4.55, "isoform마다 다른 exon", fc=ACCENT, tc=ON_ACCENT, fs=11)

# ── gene track ──────────────────────────────────────────────────────────────
ax.plot([X0, X1], [GENE_Y, GENE_Y], color="#9AA1AB", lw=2.0, zorder=1, solid_capstyle="butt")
for i, (a, b) in enumerate(EXONS):
    box(ax, a, GENE_Y - EH / 2, b - a, EH, fc=PRIMARY, ec="none", r=0.03, z=3)
    label(ax, (a + b) / 2, GENE_Y + 0.38, f"exon {i + 1}", fs=11, w="bold", c=TEXT)
# poly(A) tail
box(ax, 7.24, GENE_Y - 0.08, 0.24, 0.16, fc=STEPS[1], ec=PRIMARY, lw=0.8, r=0.07, z=3)
label(ax, 7.36, GENE_Y - 0.33, "poly(A)", fs=11, c=TEXT2)

# ── coverage tracks ─────────────────────────────────────────────────────────
x = np.linspace(X0, X1, 1600)
inexon = np.zeros_like(x, dtype=bool)
for a, b in EXONS:
    inexon |= (x >= a) & (x <= b)

rng = np.random.default_rng(6)
wobble = np.zeros_like(x)
for k in range(4):
    f = rng.uniform(3, 9)
    wobble += rng.uniform(0.008, 0.02) * np.sin(2 * np.pi * f * x + rng.uniform(0, 6.28))
y_fl = np.where(inexon, 0.52 + wobble, 0.0)

y_3p = 0.90 * np.exp(-((x - 6.90) / 0.21) ** 2)
y_3p = np.where((x >= EXONS[3][0]) & (x <= EXONS[3][1]), y_3p, 0.0)
y_3p[y_3p < 0.004] = 0.0

for base, y, c in [(FL_BASE, y_fl, C_FL), (SP_BASE, y_3p, C_3P)]:
    ax.fill_between(x, base, base + y, color=mix(c, 0.22), lw=0, zorder=2)
    ax.plot(x, base + y, color=c, lw=1.4, zorder=3)
    ax.plot([X0, X1], [base, base], color=c, lw=1.2, zorder=3)

# 3' 방식: exon 2 아래 read 거의 없음
label(ax, (bx0 + bx1) / 2, SP_BASE + 0.32, "read\n거의 없음", fs=11, w="bold", c=ON_ACCENT,
      linespacing=1.2)

# ── left labels ─────────────────────────────────────────────────────────────
LX = 0.10
label(ax, LX, GENE_Y + 0.12, "유전자", fs=13, w="bold", c=TEXT, ha="left")
label(ax, LX, GENE_Y - 0.17, "5' → 3'", fs=11, c=TEXT2, ha="left")
label(ax, LX, FL_BASE + 0.40, "full-length", fs=13, w="bold", c=C_FL, ha="left")
label(ax, LX, FL_BASE + 0.11, "예: Smart-seq2", fs=11, c=TEXT2, ha="left")
label(ax, LX, SP_BASE + 0.40, "3' 방식", fs=13, w="bold", c=C_3P, ha="left")
label(ax, LX, SP_BASE + 0.11, "예: 10x Genomics 3'", fs=11, c=TEXT2, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch06_coverage.png"))
