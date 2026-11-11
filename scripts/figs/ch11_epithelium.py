# ch11 §2 상피 조직 단면: 네 가지 접합, 꼭대기면·바닥옆면, 기저막과 결합조직 ECM  (redraw of figures/ch11_epithelium.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

fig, ax = new(8.2, 5.35)

EDGE = "#9AA1AB"
CELL = STEPS[3]
NUC = STEPS[2]
C_TJ, C_AJ, C_DS, C_GJ = CAT[0], CAT[1], CAT[2], CAT[3]

# ── geometry ───────────────────────────────────────────────────────────────
CW, GAP = 1.44, 0.12
X0 = 0.3
CY0, CY1 = 1.42, 4.5
cells = [X0 + i * (CW + GAP) for i in range(3)]
borders = [cells[0] + CW + GAP / 2, cells[1] + CW + GAP / 2]
Y_TJ, Y_AJ, Y_DS, Y_GJ = 4.22, 3.70, 3.18, 2.66
Y_DS2 = 1.75
BM_Y = 1.32                      # centre of the basement-membrane band

# ── cells, microvilli, nuclei ──────────────────────────────────────────────
for i, x0 in enumerate(cells):
    xs = np.arange(x0 + 0.12, x0 + CW - 0.08, 0.15)
    for x in xs:
        box(ax, x - 0.04, CY1 - 0.04, 0.08, 0.18, fc=CELL, ec=EDGE, lw=1.0, r=0.04, z=1)
    box(ax, x0, CY0, CW, CY1 - CY0, fc=CELL, ec=EDGE, lw=1.2, r=0.1, z=2)
    cx = x0 + CW / 2
    ax.add_patch(mp.Ellipse((cx, 2.2), 0.7, 0.85, fc=NUC, ec=EDGE, lw=1.0, zorder=3))
    if i == 0:
        label(ax, cx, 2.2, "핵", fs=11, c=TEXT2)

# ── junction glyphs ────────────────────────────────────────────────────────
def wavy(x1, x2, y, color, amp=0.025):
    xx = np.linspace(x1, x2, 60)
    ax.plot(xx, y + amp * np.sin((xx - x1) * 28), color=color, lw=0.9, zorder=3,
            solid_capstyle="round")

def desmosome(c, y, reach=0.55):
    for s in (-1, 1):
        ax.add_patch(mp.Rectangle((c + s * 0.07 - (0.06 if s < 0 else 0), y - 0.1), 0.06, 0.2,
                                  fc=C_DS, ec="none", zorder=5))
        for dy in (-0.05, 0.05):
            wavy(c + s * 0.13, c + s * reach, y + dy, mix(C_DS, 0.75))
    ax.add_patch(mp.Rectangle((c - 0.07, y - 0.06), 0.14, 0.12, fc=mix(C_DS, 0.35),
                              ec="none", zorder=4))

for c in borders:
    for dy in (-0.06, 0, 0.06):
        ax.plot([c - 0.14, c + 0.14], [Y_TJ + dy] * 2, color=C_TJ, lw=2.2, zorder=5,
                solid_capstyle="butt")
    ax.add_patch(mp.Rectangle((c - 0.11, Y_AJ - 0.08), 0.22, 0.16, fc=C_AJ, ec="none", zorder=5))
    for s in (-1, 1):
        ax.plot([c + s * 0.11, c + s * 0.42], [Y_AJ] * 2, color=C_AJ, lw=1.1,
                ls=(0, (3, 2)), zorder=3)
    desmosome(c, Y_DS)
    desmosome(c, Y_DS2, reach=0.42)
    for dy in (-0.07, 0, 0.07):
        ax.add_patch(mp.Rectangle((c - 0.12, Y_GJ + dy - 0.018), 0.24, 0.036, fc=C_GJ,
                                  ec="none", zorder=5))

for y in (Y_TJ, Y_AJ, Y_DS, Y_GJ):
    ax.plot([borders[1] + 0.62, 5.0], [y, y], color="#B8BDC5", lw=0.9, ls=":", zorder=4)

# ── basement membrane and connective-tissue ECM ───────────────────────────
XR = cells[2] + CW + 0.1
ax.add_patch(mp.Rectangle((0.2, BM_Y - 0.06), XR - 0.2, 0.12, fc=TEXT2, ec="none", zorder=2))
box(ax, 0.2, 0.1, XR - 0.2, 1.05, fc=FAINT, ec="none", r=0.08, z=1)
for i, y in enumerate([0.45, 0.6, 0.75, 0.9]):
    xx = np.linspace(0.4, 3.0, 80)
    ax.plot(xx, y + 0.035 * np.sin(xx * 2.4 + i * 1.3), color=EDGE, lw=1.6, zorder=2,
            solid_capstyle="round")
label(ax, 0.4, 0.25, "콜라겐 섬유", fs=11, c=TEXT2, ha="left")

fx, fy, fl = 3.85, 0.62, 0.6
t = np.linspace(-1, 1, 60)
hh = 0.16 * (1 - t ** 2) ** 0.8
pts = list(zip(fx + fl * t, fy + hh)) + list(zip(fx + fl * t[::-1], fy - hh[::-1]))
ax.add_patch(mp.Polygon(pts, closed=True, fc=CELL, ec=EDGE, lw=1.1, zorder=3))
ax.add_patch(mp.Ellipse((fx, fy), 0.36, 0.14, fc=NUC, ec=EDGE, lw=0.9, zorder=4))
label(ax, fx, 0.96, "섬유모세포", fs=11, c=TEXT2)

# ── polarity labels ───────────────────────────────────────────────────────
label(ax, 0.3, 5.12, "바깥 · 관 안쪽 (lumen)", fs=11, c=TEXT2, ha="left")
label(ax, 0.3, 4.86, "꼭대기면 (apical)", fs=12, w="bold", ha="left")

# ── right label column ─────────────────────────────────────────────────────
TX = 5.32

def entry(y, name, eng, note, sw=None, swh=0.16):
    """y = centre of the first (name) line."""
    if sw is not None:
        box(ax, 5.06, y - swh / 2, 0.16, swh, fc=sw, ec="none", r=0.03, z=4)
    t1 = label(ax, TX, y, name, fs=12, w="bold", ha="left")
    fig.canvas.draw()
    wd = t1.get_window_extent().width / fig.dpi
    if eng:
        label(ax, TX + wd + 0.08, y, eng, fs=11, c=TEXT2, ha="left")
    label(ax, TX, y - 0.24, note, fs=11, c=TEXT2, ha="left")

entry(Y_TJ, "밀착연접", "tight junction", "틈을 막음 · claudin, occludin", C_TJ)
entry(Y_AJ, "부착연접", "adherens junction", "E-cadherin · actin과 연결", C_AJ)
entry(Y_DS, "데스모솜", "desmosome", "중간섬유와 연결 · 힘을 견딤", C_DS)
entry(Y_GJ, "간극연접", "gap junction", "connexin 통로 · 작은 분자 이동", C_GJ)

arrow(ax, 5.24, 2.12, 4.96, 2.12, color=TEXT2, lw=1.2, ms=10)
entry(2.12, "바닥옆면", "basolateral", "이웃 세포 · 기저막 쪽")

entry(BM_Y, "기저막", "basement membrane", "라미닌 · IV형 콜라겐 · nidogen 등", TEXT2, swh=0.12)
entry(0.72, "결합조직의 ECM", "", "콜라겐 · 피브로넥틴 · 프로테오글리칸")

save(fig, os.path.join(HERE, "../../figures/ch11_epithelium.png"))
