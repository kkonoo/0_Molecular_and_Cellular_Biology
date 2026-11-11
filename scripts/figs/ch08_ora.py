# ch08 §data-box ORA: universe 안 DE 유전자와 경로 유전자의 겹침(벤 그림)과 같은 숫자의 2×2 표  (redraw of figures/ch08_ora.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.2, 3.2)

C_DE = CAT[0]     # DE 유전자
C_PW = CAT[1]     # 경로 유전자

# ── 왼쪽: universe 안의 벤 그림 ─────────────────────────────────────────────
box(ax, 0.1, 0.3, 3.7, 2.75, fc=FAINT, ec=LINE, r=0.1, z=1)
label(ax, 0.3, 2.82, "universe  N = 12,000", fs=12, w="bold", ha="left")
label(ax, 0.3, 2.55, "검정한(발현된) 유전자 전체", fs=11, c=TEXT2, ha="left")

cy, r = 1.35, 0.9
xl, xr = 1.45, 2.45
cl = mp.Circle((xl, cy), r, fc=mix(C_DE, 0.14), ec="none", zorder=2)
cr = mp.Circle((xr, cy), r, fc=mix(C_PW, 0.14), ec="none", zorder=2)
ax.add_patch(cl); ax.add_patch(cr)
# 겹침 영역 = 오른쪽 원을 왼쪽 원으로 잘라 accent로 칠함
ov = mp.Circle((xr, cy), r, fc=ACCENT, ec="none", zorder=3)
ax.add_patch(ov)
ov.set_clip_path(cl)
# 테두리
ax.add_patch(mp.Circle((xl, cy), r, fc="none", ec=C_DE, lw=1.8, zorder=4))
ax.add_patch(mp.Circle((xr, cy), r, fc="none", ec=C_PW, lw=1.8, zorder=4))

label(ax, 1.05, 1.5, "DE 유전자", fs=12, w="bold", c=C_DE)
label(ax, 1.05, 1.2, "n = 300", fs=12, c=C_DE)
label(ax, 2.85, 1.5, "경로 유전자", fs=12, w="bold", c=C_PW)
label(ax, 2.85, 1.2, "M = 200", fs=12, c=C_PW)
label(ax, 1.95, 1.5, "겹침", fs=11, c=ON_ACCENT)
label(ax, 1.95, 1.2, "k = 12", fs=12, w="bold", c=ON_ACCENT)

# ── 오른쪽: 2×2 표 ──────────────────────────────────────────────────────────
x0 = 4.15
cw0, cw = 0.8, 1.05
xs = [x0 + cw0 + cw * (i + 0.5) for i in range(3)]       # 열 중심
xlab = x0 + cw0 / 2
x_end = x0 + cw0 + 3 * cw

y_head, y_r1, y_r2, y_tot = 2.8, 2.2, 1.6, 0.9
lines_y = [2.5, 1.9]
for yy in lines_y:
    ax.plot([x0, x_end], [yy, yy], color=LINE, lw=1.2, zorder=2)
ax.plot([x0, x_end], [1.3, 1.3], color="#9AA1AB", lw=1.4, zorder=2)
ax.plot([x0, x_end], [3.08, 3.08], color="#9AA1AB", lw=1.4, zorder=2)
ax.plot([x0, x_end], [0.45, 0.45], color="#9AA1AB", lw=1.4, zorder=2)

# 열 머리
label(ax, xs[0], y_head, "DE 유전자", fs=12, w="bold", c=C_DE)
label(ax, xs[1], y_head, "DE 아님", fs=12, w="bold", c=TEXT2)
label(ax, xs[2], y_head, "합", fs=12, w="bold")
# 행 머리
label(ax, xlab, y_r1, "경로 안", fs=12, w="bold", c=C_PW)
label(ax, xlab, y_r2, "경로 밖", fs=12, w="bold", c=TEXT2)
label(ax, xlab, y_tot, "합", fs=12, w="bold")

dy = 0.13
# 1행: 경로 안
box(ax, xs[0] - cw / 2 + 0.08, 1.98, cw - 0.16, 0.44, fc=ACCENT, ec="none", r=0.06, z=3)
label(ax, xs[0], y_r1 + dy, "12", fs=12, w="bold", c=ON_ACCENT)
label(ax, xs[0], y_r1 - dy, "k", fs=11, c=ON_ACCENT)
label(ax, xs[1], y_r1, "188", fs=12)
label(ax, xs[2], y_r1 + dy, "200", fs=12, w="bold")
label(ax, xs[2], y_r1 - dy, "M", fs=11, c=TEXT2)
# 2행: 경로 밖
label(ax, xs[0], y_r2, "288", fs=12)
label(ax, xs[1], y_r2, "11,512", fs=12)
label(ax, xs[2], y_r2, "11,800", fs=12, w="bold")
# 합 행
label(ax, xs[0], y_tot + dy, "300", fs=12, w="bold")
label(ax, xs[0], y_tot - dy - 0.02, "n", fs=11, c=TEXT2)
label(ax, xs[1], y_tot + dy, "11,700", fs=12, w="bold")
label(ax, xs[2], y_tot + dy, "12,000", fs=12, w="bold")
label(ax, xs[2], y_tot - dy - 0.02, "N", fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch08_ora.png"))
