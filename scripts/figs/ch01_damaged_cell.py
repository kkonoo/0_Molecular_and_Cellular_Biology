# ch01 §4 intact vs membrane-damaged cell: cytoplasmic mRNA leaks out, mitochondrial RNA stays, MT- % rises  (redraw of figures/ch01_damaged_cell.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 8.1, 5.3
fig, ax = new(W, H)

C_MRNA = CAT[0]          # cytoplasmic mRNA / other genes (blue)
C_MT = CAT[1]            # mitochondria / MT- transcripts (teal)
C_MEM = TEXT2            # membranes
A, B = 1.4, 1.1          # cell ellipse semi-axes
CY = 3.6

def mrna(x, y):
    # a patch (not a line) so the cropping in save() sees it
    box(ax, x - 0.15, y - 0.045, 0.30, 0.09, fc=C_MRNA, ec="none", r=0.045, z=4)

def mito(x, y, ang):
    ax.add_patch(mp.Ellipse((x, y), 0.62, 0.32, angle=ang, fc=mix(C_MT, 0.12),
                            ec=C_MT, lw=1.5, zorder=3))
    u = np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
    for s in (-0.09, 0.09):
        c = np.array([x, y]) + s * u
        p0, p1 = c - 0.04 * u, c + 0.04 * u
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=C_MT, lw=3.2,
                solid_capstyle="round", zorder=4)

def membrane(cx, cy, t0, t1, lw=2.0):
    t = np.radians(np.linspace(t0, t1, 120))
    ax.plot(cx + A * np.cos(t), cy + B * np.sin(t), color=C_MEM, lw=lw,
            solid_capstyle="round", zorder=2)

def nucleus(cx, cy):
    ax.add_patch(mp.Circle((cx, cy), 0.38, fc=FAINT, ec=C_MEM, lw=1.4, zorder=3))
    label(ax, cx, cy, "핵", fs=11, c=TEXT2)

MITOS = [(0.75, 0.55, -20), (-0.8, 0.2, 25), (0.7, -0.5, -30)]
MRNAS = [(-0.5, 0.85), (0.1, 0.9), (-0.95, -0.35), (-0.3, -0.75),
         (0.25, -0.85), (1.05, 0.05), (-0.35, 0.45), (0.55, 0.05)]

# ── intact cell ───────────────────────────────────────────────────
X1 = 2.5
ax.add_patch(mp.Ellipse((X1, CY), 2 * A, 2 * B, fc="white", ec="none", zorder=1))
membrane(X1, CY, 0, 360)
nucleus(X1, CY)
for dx, dy, ang in MITOS:
    mito(X1 + dx, CY + dy, ang)
for dx, dy in MRNAS:
    mrna(X1 + dx, CY + dy)

# direct labels (left of intact cell)
label(ax, 1.05, 4.2, "미토콘드리아", fs=11, w="bold", c=C_MT, ha="right")
ax.plot([1.08, 1.48], [4.17, 3.86], color=C_MT, lw=1.0, zorder=3)
label(ax, 1.05, 2.75, "mRNA", fs=11, w="bold", c=C_MRNA, ha="right")
ax.plot([1.08, 1.38], [2.79, 3.17], color=C_MRNA, lw=1.0, zorder=3)

label(ax, X1, 1.95, "온전한 세포", fs=12.5, w="bold")
label(ax, X1, 1.66, "막이 mRNA를 가둠", fs=11, c=TEXT2)

# ── damaged cell ──────────────────────────────────────────────────
X2 = 6.0
holes = [25, 95, 155, 210, 275, 335]
gap = 9
for i, h0 in enumerate(holes):
    h1 = holes[(i + 1) % len(holes)]
    if h1 < h0:
        h1 += 360
    membrane(X2, CY, h0 + gap, h1 - gap)
nucleus(X2, CY)
for dx, dy, ang in MITOS:
    mito(X2 + dx, CY + dy, ang)
for dx, dy in [(0.55, 0.05), (-0.3, -0.75)]:
    mrna(X2 + dx, CY + dy)
for h in holes:
    c, s = np.cos(np.radians(h)), np.sin(np.radians(h))
    arrow(ax, X2 + 1.05 * A * c, CY + 1.05 * B * s, X2 + 1.17 * A * c, CY + 1.17 * B * s,
          color="#9AA1AB", lw=1.2, ms=9)
    mrna(X2 + 1.28 * A * c, CY + 1.28 * B * s)

label(ax, X2, 1.95, "막이 손상된 세포", fs=12.5, w="bold")
label(ax, X2, 1.66, "mRNA는 빠져나가고 미토콘드리아는 남음", fs=11, c=TEXT2)

# ── UMI bars (illustrative numbers from the chapter's R example) ─
X0, SCALE, BH = 1.7, 3.6 / 5000, 0.3
rows = [("온전한 세포", 5000, 200, 1.12), ("손상된 세포", 800, 200, 0.55)]
for name, total, mt, y in rows:
    label(ax, 0.3, y, name, fs=12, w="bold", ha="left")
    other_w = (total - mt) * SCALE
    mt_w = mt * SCALE
    box(ax, X0, y - BH / 2, other_w, BH, fc=C_MRNA, ec="none", r=0.02)
    box(ax, X0 + other_w, y - BH / 2, mt_w, BH, fc=C_MT, ec="none", r=0.02)
    end = X0 + other_w + mt_w
    t = label(ax, end + 0.12, y, f"{total:,} UMI 중 MT- {mt}", fs=11, ha="left")
    fig.canvas.draw()
    x1 = t.get_window_extent().x1 / fig.dpi
    if total == 5000:
        label(ax, x1 + 0.12, y, "4%", fs=12, w="bold", ha="left")
        label(ax, X0 + 0.12, y, "그 밖의 유전자", fs=11, c="white", ha="left")
        label(ax, end - mt_w / 2, y + BH / 2 + 0.13, "MT-", fs=11, w="bold", c=C_MT)
    else:
        pill(ax, x1 + 0.12, y, "25%", fc=ACCENT, tc=ON_ACCENT, fs=12, ha="left", h=0.3)

save(fig, os.path.join(HERE, "../../figures/ch01_damaged_cell.png"))
