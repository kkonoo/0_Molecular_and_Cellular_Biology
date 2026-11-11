# ch12 §5 CRISPR–Cas9: target + PAM + guide RNA + cut site, NHEJ indel → frameshift, CRISPRi/a  (redraw of figures/ch12_crispr.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.4, 5.9)

C_TGT = CAT[0]     # target / guide RNA (same sequence → same colour)
C_PAM = CAT[1]     # PAM

# ── ① target, PAM, guide, cut ─────────────────────────────────────────────
label(ax, 0.1, 5.74, "① 표적 찾기와 절단", fs=12, w="bold", ha="left")

top = "CTAGCTGAAGCCAACCTGGAGAACGGTTC"
comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
bot = "".join(comp[b] for b in top)
x0, s = 0.75, 0.25
yt, yb = 5.0, 4.55
X = lambda i: x0 + i * s

# Cas9 (flat faint shape)
box(ax, X(0.5), 3.6, X(27.5) - X(0.5), 1.9, fc=FAINT, ec=LINE, lw=1.0, r=0.3, z=1)
label(ax, X(0.5) + 0.12, 5.32, "Cas9", fs=12, w="bold", c=TEXT2, ha="left")

# target 20 bp (i = 3..22) and PAM (i = 23..25)
box(ax, X(2.5), yt - 0.17, 20 * s, 0.34, fc=mix(C_TGT, 0.18), ec="none", r=0.04, z=2)
label(ax, X(12.5), 5.32, "표적 20 bp", fs=11, w="bold", c=C_TGT)
box(ax, X(22.5), yt - 0.17, 3 * s, 0.34, fc=mix(C_PAM, 0.25), ec="none", r=0.04, z=2)
label(ax, X(24), 5.32, "PAM (NGG)", fs=11, w="bold", c=C_PAM)

for i, (a, b) in enumerate(zip(top, bot)):
    label(ax, X(i), yt, a, fs=12, w="bold", c=TEXT)
    label(ax, X(i), yb, b, fs=12, c=TEXT2)
label(ax, X(0) - 0.33, yt, "5′", fs=11, c=TEXT2)
label(ax, X(28) + 0.33, yt, "3′", fs=11, c=TEXT2)
label(ax, X(0) - 0.33, yb, "3′", fs=11, c=TEXT2)
label(ax, X(28) + 0.33, yb, "5′", fs=11, c=TEXT2)

# guide RNA pairing with the bottom strand under the target
yg = 4.1
for i in range(3, 23):
    ax.plot([X(i), X(i)], [yg + 0.03, yb - 0.13], color=C_TGT, lw=1.0, zorder=3)
ax.plot([X(3), X(25)], [yg, yg], color=C_TGT, lw=3.0, solid_capstyle="round", zorder=4)
r = 0.18
th = np.linspace(np.pi / 2, -np.pi * 1.15, 60)
ax.plot(X(25) + r * np.cos(th), yg - r + r * np.sin(th), color=C_TGT, lw=3.0,
        solid_capstyle="round", zorder=4)
label(ax, X(3) - 0.2, yg, "5′", fs=11, c=C_TGT)
label(ax, X(3), 3.8, "guide RNA: 20염기가 표적과 짝", fs=11, w="bold", c=C_TGT, ha="left")

# cut site: between i = 19 and 20 (3 bp before PAM), both strands
xc = X(19.5)
ax.plot([xc, xc], [5.22, 3.48], color=TEXT, lw=1.6, ls=(0, (3, 2)), zorder=6)
label(ax, xc, 3.33, "절단: PAM 3–4 bp 앞, 두 가닥 모두", fs=11, w="bold", c=TEXT)

# ── ② NHEJ → indel → frameshift ───────────────────────────────────────────
label(ax, 0.1, 2.92, "② 세포가 이어 붙임 (NHEJ) → indel → frameshift", fs=12, w="bold", ha="left")
rows = [("원래", ["GCT", "GAA", "GCC", "AAC", "CTG", "GAG", "AAC"]),
        ("1 bp 결실", ["GCT", "GAA", "GCC", "AAC", "CTG", "GAA", "ACG"]),
        ("1 bp 삽입", ["GCT", "GAA", "GCC", "AAC", "CTG", "GAT", "GAA"])]
cw, cg, ch = 0.58, 0.08, 0.38
cx0 = 1.1
for k, (name, cods) in enumerate(rows):
    y = 2.38 - k * 0.55
    label(ax, 0.1, y, name, fs=12, w="bold", ha="left")
    for j, c in enumerate(cods):
        shifted = k > 0 and j >= 5
        fc = ACCENT if shifted else "white"
        ec = "none" if shifted else "#C9CED6"
        box(ax, cx0 + j * (cw + cg), y - ch / 2, cw, ch, fc=fc, ec=ec, lw=1.0, r=0.06, z=2)
        label(ax, cx0 + j * (cw + cg) + cw / 2, y, c, fs=12, w="bold",
              c=ON_ACCENT if shifted else TEXT)
label(ax, 0.1, 0.6, "절단 뒤 codon 틀이 밀림 → 기능을 잃은 단백질 → knockout",
      fs=11, c=TEXT2, ha="left")

# ── CRISPRi / CRISPRa (no cutting) ─────────────────────────────────────────
px, pw = 5.95, 2.4
label(ax, px, 2.92, "자르지 않는 CRISPR", fs=12, w="bold", ha="left")
box(ax, px, 0.36, pw, 2.31, fc=FAINT, ec="none", r=0.1, z=1)
label(ax, px + 0.14, 2.42, "dCas9", fs=12, w="bold", ha="left")
label(ax, px + 0.14, 2.15, "자르는 기능을 없앤 Cas9", fs=11, c=TEXT2, ha="left")

wi = pill(ax, px + 0.14, 1.78, "CRISPRi", fc=CAT[2], ha="left")
wa = pill(ax, px + 0.14, 1.2, "CRISPRa", fc=CAT[3], ha="left")
tx = px + 0.14 + max(wi, wa) + 0.12
label(ax, tx, 1.78, "전사 낮춤", fs=11, w="bold", ha="left")
label(ax, px + 0.14, 1.51, "+ 억제 단백질", fs=11, c=TEXT2, ha="left")
label(ax, tx, 1.2, "전사 높임", fs=11, w="bold", ha="left")
label(ax, px + 0.14, 0.93, "+ 활성 단백질", fs=11, c=TEXT2, ha="left")

ax.plot([px + 0.14, px + pw - 0.14], [0.78, 0.78], color="#D5D9DF", lw=1.0, zorder=2)
label(ax, px + 0.14, 0.6, "DNA 서열은 그대로", fs=11, w="bold", ha="left")

save(fig, os.path.join(HERE, "../../figures/ch12_crispr.png"))
