# ch10 §4 조혈 계통: HSC → CMP/CLP → MEP/GMP → 혈액 세포 → 조직에서 마지막 분화, PBMC 범위  (redraw of figures/ch10_hematopoiesis.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.2, 6.6)
OX, OY = 0.15, 0.25
def X(x): return x + OX
def Y(y): return y + OY

M = CAT[0]   # 골수계 (myeloid)
L = CAT[2]   # 림프계 (lymphoid)
P = CAT[1]   # PBMC 테두리
BH = 0.38    # 세포 상자 높이
NH = 0.62    # 두 줄 상자 높이
LW = 1.4

def cell(xl, yc, w, name, ec, h=BH, sub=None, fc="white", tc=TEXT, sc=TEXT2):
    box(ax, X(xl), Y(yc - h / 2), w, h, fc=fc, ec=ec, lw=LW, r=0.07, z=2)
    if sub:
        label(ax, X(xl + w / 2), Y(yc + 0.11), name, fs=12, w="bold", c=tc)
        label(ax, X(xl + w / 2), Y(yc - 0.14), sub, fs=11, c=sc)
    else:
        label(ax, X(xl + w / 2), Y(yc), name, fs=12, w="bold", c=tc)

def path(pts, color, ls="-", head=True):
    xs = [X(p[0]) for p in pts]
    ys = [Y(p[1]) for p in pts]
    ax.plot(xs, ys, color=color, lw=LW, ls=ls, zorder=1.5,
            solid_capstyle="butt", solid_joinstyle="miter", dash_capstyle="butt")
    if head:
        arrow(ax, xs[-1] - 0.1, ys[-1], xs[-1], ys[-1], color=color, lw=LW, ms=11, z=3)

# ── 배경 구역 ────────────────────────────────────────────────
box(ax, X(0.0), Y(-0.1), 4.62, 6.3, fc=FAINT, ec="none", r=0.12, z=0)
box(ax, X(4.70), Y(-0.1), 1.75, 6.3, fc=FAINT, ec="none", r=0.12, z=0)
box(ax, X(6.55), Y(-0.1), 1.40, 6.3, fc=FAINT, ec="none", r=0.12, z=0)
label(ax, X(2.31), Y(5.93), "골수", fs=12, w="bold", c=TEXT2)
label(ax, X(5.575), Y(5.93), "혈액", fs=12, w="bold", c=TEXT2)
label(ax, X(7.25), Y(5.93), "조직", fs=12, w="bold", c=TEXT2)

# ── PBMC 범위 ───────────────────────────────────────────────
box(ax, X(4.95), Y(0.32), 1.45, 2.73, fc=mix(P, 0.08), ec=P, lw=1.6, r=0.1, z=1)
label(ax, X(5.675), Y(0.12), "PBMC", fs=12, w="bold", c=P)

# ── 줄기세포·전구세포 ─────────────────────────────────────────
cell(0.10, 2.70, 1.08, "HSC", TEXT, h=NH, sub="조혈줄기세포", fc=TEXT, tc="white", sc="white")
cell(1.45, 4.25, 1.10, "CMP", M, h=NH, sub="골수계 전구")
cell(1.45, 1.15, 1.10, "CLP", L, h=NH, sub="림프계 전구")
cell(2.90, 5.20, 1.65, "MEP", M, h=NH, sub="거핵세포·적혈구 전구")
cell(2.90, 3.30, 1.65, "GMP", M, h=NH, sub="과립구·대식세포 전구")

# HSC → CMP / CLP
ax.plot([X(1.18), X(1.30)], [Y(2.70), Y(2.70)], color="#9AA1AB", lw=LW, zorder=1.5,
        solid_capstyle="butt")
path([(1.30, 2.70), (1.30, 4.25), (1.45, 4.25)], M)
path([(1.30, 2.70), (1.30, 1.15), (1.45, 1.15)], L)

# CMP → MEP / GMP
path([(2.55, 4.25), (2.72, 4.25), (2.72, 5.20), (2.90, 5.20)], M)
path([(2.72, 4.25), (2.72, 3.30), (2.90, 3.30)], M)

# ── 혈액 속 세포 ─────────────────────────────────────────────
BX, BW = 5.05, 1.25
mep_cells = [("적혈구", 5.45), ("혈소판", 4.95)]
gmp_cells = [("호중구", 4.35), ("호산구", 3.85), ("호염기구", 3.35),
             ("단핵구", 2.75), ("수지상세포", 2.25)]
clp_cells = [("NK세포", 1.65), ("B세포", 1.15), ("T세포", 0.65)]

for name, y in mep_cells:
    cell(BX, y, BW, name, M)
    path([(4.55, 5.20), (4.85, 5.20), (4.85, y), (BX, y)], M)
for name, y in gmp_cells:
    cell(BX, y, BW, name, M)
    path([(4.55, 3.30), (4.85, 3.30), (4.85, y), (BX, y)], M)
for name, y in clp_cells:
    cell(BX, y, BW, name, L)

# CLP → NK / B / T(흉선 경유)
path([(2.55, 1.15), (BX, 1.15)], L)
path([(4.85, 1.15), (4.85, 1.65), (BX, 1.65)], L)
path([(3.15, 1.15), (3.15, 0.65), (BX, 0.65)], L)
pill(ax, X(3.95), Y(0.65), "흉선", fc=mix(L, 0.18), tc=shade(L, 0.35), fs=11, h=0.3)

# ── 조직에서 마지막 분화 ──────────────────────────────────────
TX, TW = 6.70, 1.10
cell(TX, 4.65, TW, "비만세포", M)
cell(TX, 2.75, TW, "대식세포", M)
cell(TX, 1.15, TW, "형질세포", L)

path([(4.85, 4.35), (4.85, 4.65), (TX, 4.65)], M, ls=(0, (4, 3)))
path([(BX + BW, 2.75), (TX, 2.75)], M)
path([(BX + BW, 1.15), (TX, 1.15)], L)

save(fig, os.path.join(HERE, "../../figures/ch10_hematopoiesis.png"))
