# ch11 §3 분비 신호·접촉 신호(Notch)와 scRNA-seq에서 보이는 리간드·수용체 mRNA  (redraw of figures/ch11_signaling.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

LIG = PRIMARY          # 리간드 (모든 패널 공통)
REC = CAT[1]           # 수용체 (모든 패널 공통)
MEM = "#9AA1AB"        # 세포막 윤곽

fig, ax = new(8.2, 6.85)

# ── (가) 분비 신호 ────────────────────────────────────────────────
label(ax, 0.0, 6.62, "(가) 분비 신호", fs=12, w="bold", ha="left")
box(ax, 0, 3.85, 3.95, 2.55, fc=FAINT, ec="none", z=1)

cy, a, b = 5.1, 0.65, 1.0
RX = 3.15
for cx in (0.8, RX):
    ax.add_patch(mp.Ellipse((cx, cy), 2 * a, 2 * b, fc="white", ec=MEM, lw=1.5, zorder=2))
label(ax, 0.8, 5.27, "송신 세포", fs=12, w="bold")
label(ax, 0.8, 4.95, "리간드 분비", fs=11, c=TEXT2)
label(ax, RX + 0.05, 5.27, "수신 세포", fs=12, w="bold")
label(ax, RX + 0.05, 4.95, "수용체", fs=11, c=REC, w="bold")

# 수용체 (수신 세포 왼쪽 막 위)
for ry in (4.72, 5.1, 5.48):
    xm = RX - a * np.sqrt(1 - ((ry - cy) / b) ** 2)
    ax.plot([xm, xm - 0.12], [ry, ry], color=REC, lw=2.6, zorder=4, solid_capstyle="round")
    ax.plot([xm - 0.12, xm - 0.23], [ry, ry + 0.08], color=REC, lw=2.6, zorder=4, solid_capstyle="round")
    ax.plot([xm - 0.12, xm - 0.23], [ry, ry - 0.08], color=REC, lw=2.6, zorder=4, solid_capstyle="round")
    if ry == 5.1:   # 결합한 리간드 (수용체 오목한 곳에)
        ax.add_patch(mp.Circle((xm - 0.27, ry), 0.07, fc=LIG, ec="none", zorder=5))

# 떠다니는 리간드
for (x, y) in [(1.68, 5.62), (1.93, 5.28), (1.72, 4.9), (2.13, 5.78), (2.0, 4.66), (1.58, 5.2)]:
    ax.add_patch(mp.Circle((x, y), 0.07, fc=LIG, ec="none", zorder=4))
arrow(ax, 1.52, 4.38, 2.55, 4.38, color=LIG, lw=1.6)
label(ax, 2.03, 4.13, "리간드", fs=11, c=LIG, w="bold")

# ── (나) 접촉 신호 (Notch) ────────────────────────────────────────
label(ax, 4.15, 6.62, "(나) 접촉 신호 (Notch)", fs=12, w="bold", ha="left")
box(ax, 4.15, 3.85, 3.95, 2.55, fc=FAINT, ec="none", z=1)

box(ax, 4.3, 4.0, 1.65, 2.25, fc="white", ec=MEM, lw=1.5, r=0.15, z=2)
box(ax, 6.05, 4.0, 1.9, 2.25, fc="white", ec=MEM, lw=1.5, r=0.15, z=2)

label(ax, 5.05, 6.0, "송신 세포", fs=12, w="bold")
label(ax, 5.05, 5.7, "막에 박힌 리간드", fs=11, c=LIG, w="bold")
label(ax, 5.05, 5.45, "(Delta, Jagged)", fs=11, c=TEXT2)
label(ax, 7.0, 6.0, "수신 세포", fs=12, w="bold")
label(ax, 7.0, 5.7, "Notch 수용체", fs=11, c=REC, w="bold")

# 위: 리간드–Notch 결합
ax.add_patch(mp.Rectangle((5.77, 5.09), 0.24, 0.12, fc=LIG, ec="none", zorder=4))
ax.add_patch(mp.Rectangle((6.01, 5.07), 0.26, 0.16, fc=REC, ec="none", zorder=4))
# 아래: 잘린 뒤 NICD가 떨어져 핵으로
ax.add_patch(mp.Rectangle((5.77, 4.49), 0.24, 0.12, fc=LIG, ec="none", zorder=4))
ax.add_patch(mp.Rectangle((6.01, 4.47), 0.14, 0.16, fc=REC, ec="none", zorder=4))
box(ax, 6.38, 4.46, 0.4, 0.18, fc=REC, ec="none", r=0.05, z=4)
label(ax, 6.58, 4.22, "NICD", fs=11, c=REC, w="bold")
ax.add_patch(mp.Ellipse((7.5, 4.55), 0.7, 0.55, fc=STEPS[2], ec=STEPS[0], lw=1.2, zorder=3))
label(ax, 7.5, 4.55, "핵", fs=11, w="bold")
arrow(ax, 6.84, 4.55, 7.13, 4.55, color=REC, lw=1.6)

# ── (다) cell type별 mRNA ─────────────────────────────────────────
label(ax, 0.0, 3.5, "(다) cell type별 mRNA", fs=12, w="bold", ha="left")
box(ax, 0, 0.2, 8.1, 3.1, fc=FAINT, ec="none", z=1)

base = 0.95
label(ax, 0.35, 2.98, "평균 발현", fs=11, c=TEXT2, ha="left")
ax.plot([0.35, 3.95], [base, base], color=MEM, lw=1.2, zorder=3)
groups = [("type A", "송신 후보", 1.75, 0.12),
          ("type B", "수신 후보", 0.15, 1.4),
          ("type C", "", 0.4, 0.46)]          # 가상의 값
gcs = [0.95, 2.15, 3.35]
for (name, sub, hl, hr), gc in zip(groups, gcs):
    ax.add_patch(mp.Rectangle((gc - 0.34, base), 0.32, hl, fc=LIG, ec="none", zorder=4))
    ax.add_patch(mp.Rectangle((gc + 0.02, base), 0.32, hr, fc=REC, ec="none", zorder=4))
    label(ax, gc, 0.72, name, fs=12, w="bold")
    if sub:
        label(ax, gc, 0.46, sub, fs=11, c=TEXT2)
label(ax, 1.0, base + 1.75 - 0.08, "리간드 L mRNA", fs=11, c=LIG, w="bold", ha="left")
label(ax, 2.58, base + 1.4 - 0.1, "수용체 R mRNA", fs=11, c=REC, w="bold", ha="left")

# 가운데: 후보 고르기
MC = 5.07
label(ax, MC, 2.2, "A의 L  +  B의 R", fs=12, w="bold")
arrow(ax, MC, 1.98, MC, 1.6, color=MEM, lw=1.6)
pill(ax, MC, 1.36, "A → B 상호작용 후보", fc=ACCENT, tc=ON_ACCENT, fs=12, h=0.38)

# 오른쪽: 데이터에 없는 것
box(ax, 6.2, 0.95, 1.75, 2.0, fc="white", ec=LINE, lw=1.0, r=0.08, z=2)
label(ax, 6.35, 2.68, "데이터에 없는 것", fs=12, w="bold", ha="left")
for i, t in enumerate(["· 단백질 양", "· 실제로 분비됐나", "· 두 세포가 가까운가", "· 정말 결합했나"]):
    label(ax, 6.35, 2.3 - 0.34 * i, t, fs=11, ha="left")

save(fig, os.path.join(HERE, "../../figures/ch11_signaling.png"))
