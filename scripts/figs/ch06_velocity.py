# ch06 §RNA velocity: u/s 생성·splicing·분해 단계와 (s, u) 평면의 정상상태 선 위/아래 세포  (redraw of figures/ch06_velocity.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.3, 4.45)

C_ON, C_OFF = CAT[2], CAT[1]          # 켜지는 중 / 꺼지는 중
EXON = PRIMARY

# ── panel markers & divider ────────────────────────────────────────────────
label(ax, 0.15, 4.28, "① 단계", fs=12, w="bold", ha="left")
label(ax, 4.15, 4.28, "② 세포마다 (s, u)", fs=12, w="bold", ha="left")
ax.plot([3.95, 3.95], [0.45, 4.2], color=LINE, lw=1.0, zorder=1)

# ── ① steps ────────────────────────────────────────────────────────────────
box(ax, 0.15, 1.75, 2.3, 2.2, fc=FAINT, ec="#C9CED6", lw=1.0, r=0.12, z=1, ls="--")
label(ax, 0.28, 3.8, "핵", fs=11, c=TEXT2, ha="left")

def exon(x, y, w=0.4, h=0.16):
    box(ax, x, y - h / 2, w, h, fc=EXON, ec="none", r=0.03, z=3)

# DNA
label(ax, 1.3, 3.66, "유전자 (DNA)", fs=12, w="bold")
ax.plot([0.5, 2.1], [3.38, 3.38], color=TEXT2, lw=1.8, zorder=2)
exon(0.75, 3.38); exon(1.45, 3.38)

# α 전사
arrow(ax, 1.3, 3.22, 1.3, 2.76, color=TEXT2)
label(ax, 1.42, 2.99, "α 전사", fs=12, w="bold", ha="left")

# unspliced
ax.plot([1.15, 1.45], [2.55, 2.55], color=TEXT2, lw=2.6, zorder=2)
exon(0.75, 2.55); exon(1.45, 2.55)
label(ax, 1.3, 2.25, "unspliced (u)", fs=12, w="bold")
label(ax, 1.3, 2.02, "intron이 남은 RNA", fs=11, c=TEXT2)

# β splicing
arrow(ax, 1.97, 2.45, 2.75, 1.93, color=TEXT2)
label(ax, 2.55, 2.5, "β splicing", fs=12, w="bold", ha="left")

# spliced
exon(2.6, 1.75, w=0.35); exon(2.97, 1.75, w=0.35)
box(ax, 3.32, 1.75 - 0.06, 0.24, 0.12, fc="white", ec="#9AA1AB", lw=1.0, r=0.06, z=3)
label(ax, 3.05, 1.45, "spliced (s)", fs=12, w="bold")
label(ax, 3.05, 1.22, "성숙 mRNA", fs=11, c=TEXT2)

# γ 분해
arrow(ax, 3.05, 1.05, 3.05, 0.62, color=TEXT2)
label(ax, 3.17, 0.84, "γ 분해", fs=12, w="bold", ha="left")

# equations
label(ax, 0.2, 1.32, "du/dt = α − β·u", fs=12, ha="left")
label(ax, 0.2, 0.97, "ds/dt = β·u − γ·s", fs=12, ha="left")

# ── ② phase plane (model: du/dt = α − βu, ds/dt = βu − γs) ─────────────────
alpha, beta, gamma = 1.0, 1.0, 0.4
s_ss, u_ss = alpha / gamma, alpha / beta           # 정상상태 (2.5, 1)
x0, y0 = 4.55, 0.8
smax, umax = 2.8, 1.15
sx, sy = 3.3 / smax, 2.8 / umax
X = lambda s: x0 + np.asarray(s) * sx
Y = lambda u: y0 + np.asarray(u) * sy

# axes
arrow(ax, x0, y0, x0 + 3.45, y0, color="#9AA1AB", lw=1.2, ms=10)
arrow(ax, x0, y0, x0, y0 + 3.0, color="#9AA1AB", lw=1.2, ms=10)
label(ax, x0 + 3.45, y0 - 0.25, "spliced (s)", fs=12, w="bold", ha="right")
label(ax, x0 - 0.08, y0 + 3.15, "unspliced (u)", fs=12, w="bold", ha="left")

# steady-state line u = (γ/β) s
sl = np.array([0, smax])
ax.plot(X(sl), Y(gamma / beta * sl), color=TEXT2, lw=1.3, ls=(0, (5, 4)), zorder=2)
label(ax, X(smax), Y(gamma / beta * smax) + 0.17, "정상상태 선", fs=11, c=TEXT2, ha="right")

# analytic solutions
def induction(t):      # from (0, 0), α on
    u = u_ss * (1 - np.exp(-beta * t))
    B = alpha / (beta - gamma)
    s = s_ss + (-s_ss - B) * np.exp(-gamma * t) + B * np.exp(-beta * t)
    return s, u

def repression(t):     # from steady state, α = 0
    u = u_ss * np.exp(-beta * t)
    B = beta * u_ss / (gamma - beta)
    s = (s_ss - B) * np.exp(-gamma * t) + B * np.exp(-beta * t)
    return s, u

tt = np.linspace(0, 12, 400)
rng = np.random.default_rng(6)
# (function, colour, sample times, arrow times, side of the curve for the arrow: +1 left of travel)
for f, col, tp, ta, side in [
        (induction, C_ON, [0.25, 0.5, 0.8, 1.2, 1.6, 2.1, 2.7, 3.4, 4.3, 5.5, 7.5], (0.95, 1.45), +1),
        (repression, C_OFF, [0.3, 0.6, 0.95, 1.35, 1.8, 2.4, 3.1, 4.0, 5.2, 7.0, 10.0], (1.5, 2.3), -1)]:
    s, u = f(tt)
    ax.plot(X(s), Y(u), color=mix(col, 0.55), lw=1.2, ls=(0, (2, 2.5)), zorder=2)
    sp, up = f(np.array(tp))
    sp = sp + rng.normal(0, 0.03, len(tp))
    up = up + rng.normal(0, 0.015, len(tp))
    ax.scatter(X(sp), Y(np.clip(up, 0.005, None)), s=34, color=col, ec="white", lw=0.6, zorder=4)
    # direction arrow drawn just beside the curve, so it does not sit on the dots or the axis
    (sa, ua), (sb, ub) = f(np.array(ta[0])), f(np.array(ta[1]))
    xa, ya, xb, yb = float(X(sa)), float(Y(ua)), float(X(sb)), float(Y(ub))
    dx, dy = xb - xa, yb - ya
    n = np.hypot(dx, dy)
    ox, oy = -dy / n * 0.16 * side, dx / n * 0.16 * side
    arrow(ax, xa + ox, ya + oy, xb + ox, yb + oy, color=col, lw=1.6, ms=12, z=5)

# steady-state point
ax.scatter([X(s_ss)], [Y(u_ss)], s=34, color=C_ON, ec="white", lw=0.6, zorder=4)

# direct labels (between curve and line)
label(ax, float(X(1.25)), float(Y(0.76)), "켜지는 중", fs=12, w="bold", c=C_ON)
label(ax, float(X(1.75)), float(Y(0.42)), "꺼지는 중", fs=12, w="bold", c=C_OFF)

save(fig, os.path.join(HERE, "../../figures/ch06_velocity.png"))
