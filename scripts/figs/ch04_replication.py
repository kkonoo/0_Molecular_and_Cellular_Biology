# ch04 §1 반보존적 복제: 부모 DNA → 복제 중 → 1회 → 2회 복제 후 가닥 구성  (redraw of figures/ch04_replication.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

ORIG = PRIMARY          # 원래 가닥: 파랑
NEW = CAT[2]            # 새로 만든 가닥: 초록
POL = TEXT2             # DNA 중합효소: 회색
RUNG = "#B8BEC7"
LW = 4.2
GAP = 0.16              # 두 가닥 사이 간격 (in)

fig, ax = new(8.3, 4.3)


def strand(x0, y0, x1, y1, c):
    ax.plot([x0, x1], [y0, y1], color=c, lw=LW, solid_capstyle="round", zorder=3)


def rungs(p, q, n_off, step=0.085):
    """p→q 선분을 따라, 법선 n_off 방향으로 GAP만큼 떨어진 가닥과 잇는 가로대."""
    p, q = np.array(p), np.array(q)
    L = np.linalg.norm(q - p)
    d = (q - p) / L
    for t in np.arange(0.06, L - 0.04, step):
        a = p + d * t
        b = a + n_off * GAP
        ax.plot([a[0], b[0]], [a[1], b[1]], color=RUNG, lw=0.9, zorder=2)


def duplex(x0, x1, yc, top, bot):
    rungs((x0, yc + GAP / 2), (x1, yc + GAP / 2), np.array([0, -1.0]))
    strand(x0, yc + GAP / 2, x1, yc + GAP / 2, top)
    strand(x0, yc - GAP / 2, x1, yc - GAP / 2, bot)


YC = 2.2
HEAD_Y = 4.05

# ── 단계 1: 부모 DNA ──────────────────────────────────────────
S1 = (0.25, 1.6)
label(ax, sum(S1) / 2, HEAD_Y, "부모 DNA", fs=12.5, w="bold")
duplex(S1[0], S1[1], YC, ORIG, ORIG)
label(ax, sum(S1) / 2, YC - 0.4, "원래 가닥 2개", fs=11, c=TEXT2)

arrow(ax, 1.8, YC, 2.2, YC)

# ── 단계 2: 복제 중 (복제 분기점) ──────────────────────────────
X0, XF = 2.35, 3.0                # 아직 안 풀린 부분
label(ax, 3.2, HEAD_Y, "복제 중", fs=12.5, w="bold")
duplex(X0, XF, YC, ORIG, ORIG)

for s in (+1, -1):                # +1 = 위쪽 팔, -1 = 아래쪽 팔
    p0 = np.array([XF, YC + s * GAP / 2])
    pb = np.array([3.3, YC + s * 0.42])
    p1 = np.array([4.05, YC + s * 0.66])
    strand(*p0, *pb, ORIG)
    strand(*pb, *p1, ORIG)
    d = (p1 - pb) / np.linalg.norm(p1 - pb)
    n_in = np.array([d[1], -d[0]]) * s          # 두 팔 사이(안쪽)를 향하는 법선
    a = pb + d * 0.12 + n_in * GAP
    b = p1 + n_in * GAP
    rungs(pb + d * 0.12, p1, n_in)
    strand(*a, *b, NEW)
    ax.add_patch(mp.Circle(tuple(a), 0.1, fc=POL, ec="white", lw=1.2, zorder=5))
    if s == -1:
        pol_xy = a

label(ax, 3.62, YC + 1.05, "원래 가닥", fs=11.5, w="bold", c=ORIG)
label(ax, 3.82, YC, "새 가닥", fs=11.5, w="bold", c=NEW)
label(ax, 3.05, YC - 1.05, "DNA 중합효소", fs=11, c=POL)
ax.plot([pol_xy[0], 3.15], [pol_xy[1] - 0.1, YC - 0.92], color=POL, lw=0.9, zorder=4)

arrow(ax, 4.3, YC, 4.7, YC)

# ── 단계 3: 1회 복제 후 ──────────────────────────────────────
S3 = (4.9, 6.2)
label(ax, sum(S3) / 2, HEAD_Y, "1회 복제 후", fs=12.5, w="bold")
duplex(S3[0], S3[1], 2.75, ORIG, NEW)
duplex(S3[0], S3[1], 1.65, NEW, ORIG)
label(ax, sum(S3) / 2, 1.27, "원래 1 + 새 1", fs=11, c=TEXT2)

arrow(ax, 6.4, YC, 6.8, YC)

# ── 단계 4: 2회 복제 후 ──────────────────────────────────────
S4 = (6.95, 8.15)
label(ax, sum(S4) / 2, HEAD_Y, "2회 복제 후", fs=12.5, w="bold")
duplex(S4[0], S4[1], 3.45, ORIG, NEW)
duplex(S4[0], S4[1], 2.9, NEW, NEW)
duplex(S4[0], S4[1], 1.6, NEW, NEW)
duplex(S4[0], S4[1], 1.05, NEW, ORIG)
label(ax, sum(S4) / 2, 0.72, "원래 가닥 있는 것 2개\n새 가닥만 2개", fs=11, c=TEXT2,
      va="top", linespacing=1.35)

save(fig, os.path.join(HERE, "../../figures/ch04_replication.png"))
