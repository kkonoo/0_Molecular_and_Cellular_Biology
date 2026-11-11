# ch11 §3 같은 가상 조직 위에 놓은 Visium spot · Visium HD bin · 이미징 기반 분자 위치  (redraw of figures/ch11_spot.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *
import numpy as np
import matplotlib.patches as mp

# ── geometry ────────────────────────────────────────────────────────────────
FIELD = 150.0                 # µm shown in each panel
P, GAP = 2.55, 0.25           # panel size / gap (in)
S = P / FIELD                 # in per µm
X0 = [0.1 + i * (P + GAP) for i in range(3)]
PY = 1.25                     # panel bottom
W, H = 8.35, PY + P + 0.62
fig, ax = new(W, H)

COL = {"E": CAT[0], "F": CAT[1], "I": CAT[2]}     # 상피 / 섬유모 / 면역
RAD = {"E": 8.5, "F": 7.5, "I": 5.0}
GRID_C = "#7C838D"

# ── one cell layout (fixed seed), shared by all panels ──────────────────────
C = np.array([75.0, 72.0]); RS = 27.5             # spot centre, radius (55 µm diameter)
inner = [("E", -17, 12), ("E", 1, 18), ("E", -19, -6),
         ("F", 18, 10), ("F", -4, 0), ("F", 8, -14), ("F", -10, -20),
         ("I", 20, -5)]                           # 8 cells centred inside the spot: 3 / 4 / 1
cells = [(t, C[0] + dx, C[1] + dy, RAD[t]) for t, dx, dy in inner]
rng = np.random.default_rng(11)

def ctype(x, y):
    t = (x + (FIELD - y)) / (2 * FIELD)
    return "E" if t < 0.45 else ("F" if t < 0.72 else "I")

for _ in range(8000):
    x, y = rng.uniform(-2, FIELD + 2, 2)
    t = ctype(x, y)
    r = RAD[t] * rng.uniform(0.93, 1.07)
    if np.hypot(x - C[0], y - C[1]) < RS + 2.0:
        continue
    gap = 7.0 if t == "I" else 2.0
    if any(np.hypot(x - a, y - b) < r + rr + gap for _, a, b, rr in cells):
        continue
    cells.append((t, x, y, r))

# ── helpers: circles clipped to a rectangle (so nothing spills past the panel) ─
def clip(pts, lo, hi):
    def inter(p, q, a, v):
        return p + (q - p) * ((v - p[a]) / (q[a] - p[a]))
    for a in (0, 1):
        for v, side in ((lo[a], 1), (hi[a], -1)):
            if not pts:
                return pts
            out = []
            for i in range(len(pts)):
                cur, prev = pts[i], pts[i - 1]
                cin = (cur[a] - v) * side >= 0
                pin = (prev[a] - v) * side >= 0
                if cin:
                    if not pin:
                        out.append(inter(prev, cur, a, v))
                    out.append(cur)
                elif pin:
                    out.append(inter(prev, cur, a, v))
            pts = out
    return pts

def fill_circle(m, lo, hi, cx, cy, r, fc, z):
    th = np.linspace(0, 2 * np.pi, 90, endpoint=False)
    pts = [np.array([cx + r * np.cos(t), cy + r * np.sin(t)]) for t in th]
    poly = clip(pts, lo, hi)
    if len(poly) >= 3:
        ax.add_patch(mp.Polygon([m(p) for p in poly], closed=True, fc=fc, ec="none", zorder=z))

def ring(m, lo, hi, cx, cy, r, color, lw, z, ls="-"):
    th = np.linspace(0, 2 * np.pi, 240)
    x = cx + r * np.cos(th); y = cy + r * np.sin(th)
    out = (x < lo[0]) | (x > hi[0]) | (y < lo[1]) | (y > hi[1])
    x[out] = np.nan; y[out] = np.nan
    pts = [m((a, b)) for a, b in zip(x, y)]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=lw, ls=ls,
            zorder=z, solid_capstyle="butt")

def pmap(k):
    return lambda p: (X0[k] + p[0] * S, PY + p[1] * S)

LO, HI = (0, 0), (FIELD, FIELD)

def solid_cells(m, lo, hi, z=2):
    for t, x, y, r in cells:
        c = COL[t]
        fill_circle(m, lo, hi, x, y, r, mix(c, 0.22), z)
        fill_circle(m, lo, hi, x, y, r * (0.45 if t == "I" else 0.38), mix(c, 0.65), z + 0.1)
        ring(m, lo, hi, x, y, r, c, 1.0, z + 0.2)

# ── panel frames and headers ───────────────────────────────────────────────
heads = [("(가) Visium", "55 µm spot"),
         ("(나) Visium HD", "2 µm 칸 → 8 µm bin"),
         ("(다) 이미징 기반", "분자 하나의 위치")]
for k in range(3):
    box(ax, X0[k], PY, P, P, fc="white", ec="none", r=0.04, z=0.5)
    box(ax, X0[k], PY, P, P, fc="none", ec="#B8BEC7", lw=1.0, r=0.04, z=20)
    label(ax, X0[k], PY + P + 0.43, heads[k][0], fs=12, w="bold", ha="left")
    label(ax, X0[k], PY + P + 0.17, heads[k][1], fs=11, c=TEXT2, ha="left")

# ── (가) Visium spot ───────────────────────────────────────────────────────
m0 = pmap(0)
solid_cells(m0, LO, HI)
for a in range(6):                                # neighbouring spots, 100 µm apart
    nx = C[0] + 100 * np.cos(np.pi / 3 * a); ny = C[1] + 100 * np.sin(np.pi / 3 * a)
    ring(m0, LO, HI, nx, ny, RS, TEXT2, 1.3, 6, ls=(0, (4, 3)))
ring(m0, LO, HI, C[0], C[1], RS, TEXT, 2.4, 7)
# scale bar
box(ax, X0[0] + 0.07, PY + 0.06, 1.02, 0.4, fc="white", ec="none", r=0.05, z=8)
ax.plot([X0[0] + 0.17, X0[0] + 0.17 + 50 * S], [PY + 0.15] * 2, color=TEXT, lw=2.2,
        solid_capstyle="butt", zorder=9)
label(ax, X0[0] + 0.17 + 25 * S, PY + 0.32, "50 µm", fs=11, z=9)

# notes under (가)
label(ax, X0[0], PY - 0.25, "지름 55 µm · 중심 간격 100 µm", fs=11, ha="left")
pill(ax, X0[0], PY - 0.57, "spot count = 세포 8개의 합", fc=ACCENT, tc=ON_ACCENT, ha="left")
segs = [("상피 3", 3, CAT[0]), ("섬유모 4", 4, CAT[1]), ("면역 1", 1, CAT[2])]
xb = X0[0]
for i, (txt, n, c) in enumerate(segs):
    w = P * n / 8
    ax.add_patch(mp.Rectangle((xb, PY - 0.95), w, 0.17, fc=c, ec="white", lw=1.2, zorder=3))
    if i < 2:
        label(ax, xb + w / 2, PY - 1.12, txt, fs=11, w="bold", c=c)
    else:
        label(ax, xb + w, PY - 1.12, txt, fs=11, w="bold", c=c, ha="right")
    xb += w

# ── (나) Visium HD grid ────────────────────────────────────────────────────
m1 = pmap(1)
solid_cells(m1, LO, HI)
for v in np.arange(8, FIELD, 8):
    ax.plot([X0[1] + v * S] * 2, [PY, PY + P], color=GRID_C, lw=0.45, zorder=6)
    ax.plot([X0[1], X0[1] + P], [PY + v * S] * 2, color=GRID_C, lw=0.45, zorder=6)
B0 = (64.0, 72.0)                                 # highlighted 8 µm bin (spans two cells)
bx, by = m1(B0)
ax.add_patch(mp.Rectangle((bx, by), 8 * S, 8 * S, fill=False, ec=TEXT, lw=1.8, zorder=9))
# magnified bin
ibx, iby, ibw, ibh = X0[1] + P - 1.08, PY + P - 1.12, 1.0, 1.04
box(ax, ibx, iby, ibw, ibh, fc="white", ec="#9AA1AB", lw=1.0, r=0.05, z=10)
IG = 0.72; IS = IG / 8
ix0, iy0 = ibx + (ibw - IG) / 2, iby + 0.27
mi = lambda p: (ix0 + (p[0] - B0[0]) * IS, iy0 + (p[1] - B0[1]) * IS)
blo, bhi = B0, (B0[0] + 8, B0[1] + 8)
for t, x, y, r in cells:
    if np.hypot(max(abs(x - (B0[0] + 4)) - 4, 0), max(abs(y - (B0[1] + 4)) - 4, 0)) < r:
        c = COL[t]
        fill_circle(mi, blo, bhi, x, y, r, mix(c, 0.22), 11)
        fill_circle(mi, blo, bhi, x, y, r * 0.38, mix(c, 0.65), 11.1)
        ring(mi, blo, bhi, x, y, r, c, 1.4, 11.2)
for v in (2, 4, 6):
    ax.plot([ix0 + v * IS] * 2, [iy0, iy0 + IG], color=GRID_C, lw=0.9, zorder=12)
    ax.plot([ix0, ix0 + IG], [iy0 + v * IS] * 2, color=GRID_C, lw=0.9, zorder=12)
ax.add_patch(mp.Rectangle((ix0, iy0), IG, IG, fill=False, ec=TEXT, lw=1.8, zorder=12))
label(ax, ibx + ibw / 2, iby + 0.13, "2 µm × 16칸", fs=11, z=12)
ax.plot([bx + 8 * S, ix0], [by + 8 * S, iy0], color=TEXT2, lw=0.8, zorder=9)
ax.plot([bx, ix0], [by + 8 * S, iy0 + IG], color=TEXT2, lw=0.8, zorder=9)

label(ax, X0[1], PY - 0.25, "8 µm bin = 2 µm 칸 4 × 4", fs=11, ha="left")
label(ax, X0[1], PY - 0.53, "bin 경계 ≠ 세포 경계", fs=11, ha="left")

# ── (다) imaging-based: molecules + segmentation ───────────────────────────
m2 = pmap(2)
NMOL = {"E": 11, "F": 9, "I": 6}
for t, x, y, r in cells:
    c = COL[t]
    fill_circle(m2, LO, HI, x, y, r, mix(c, 0.10), 2)
    fill_circle(m2, LO, HI, x, y, r * (0.45 if t == "I" else 0.38), mix(c, 0.22), 2.1)
    ring(m2, LO, HI, x, y, r * 1.08, c, 1.0, 3, ls=(0, (2.5, 1.8)))
    a = rng.uniform(0, 2 * np.pi, NMOL[t]); d = r * 0.88 * np.sqrt(rng.uniform(0, 1, NMOL[t]))
    px, py = x + d * np.cos(a), y + d * np.sin(a)
    ok = (px > 1) & (px < FIELD - 1) & (py > 1) & (py < FIELD - 1)
    pts = [m2((u, v)) for u, v in zip(px[ok], py[ok])]
    if pts:
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=7, color=c, lw=0, zorder=4)
grey = []
while len(grey) < 40:
    u, v = rng.uniform(2, FIELD - 2, 2)
    if all(np.hypot(u - a, v - b) > rr * 1.15 for _, a, b, rr in cells):
        grey.append(m2((u, v)))
ax.scatter([p[0] for p in grey], [p[1] for p in grey], s=7, color="#9AA1AB", lw=0, zorder=4)

gx = X0[2]
for i, c in enumerate([CAT[0], CAT[1], CAT[2]]):
    ax.scatter([gx + 0.06 + 0.09 * i], [PY - 0.25], s=14, color=c, lw=0, zorder=4)
label(ax, gx + 0.36, PY - 0.25, "점 하나 = 분자 하나", fs=11, ha="left")
ax.add_patch(mp.Circle((gx + 0.15, PY - 0.53), 0.09, fill=False, ec=TEXT2, lw=1.1,
                       ls=(0, (2.5, 1.8)), zorder=4))
label(ax, gx + 0.36, PY - 0.53, "추정 세포 경계", fs=11, ha="left")
ax.scatter([gx + 0.15], [PY - 0.81], s=14, color="#9AA1AB", lw=0, zorder=4)
label(ax, gx + 0.36, PY - 0.81, "소속이 애매한 분자", fs=11, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch11_spot.png"))
