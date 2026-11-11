# A3 §3 GEO와 SRA 번호의 대응: 연구·시료·run 단위로 이어지는 모습  (redraw of figures/A3_geo_sra.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.6, 4.95
fig, ax = new(W, H)

GEO_C, SRA_C = CAT[0], CAT[1]
GREY = "#9AA1AB"

# ── panels ─────────────────────────────────────────────
box(ax, 0.10, 0.10, 3.20, 4.75, fc=FAINT, ec="none", r=0.12, z=0)
box(ax, 4.30, 0.10, 3.20, 4.75, fc=FAINT, ec="none", r=0.12, z=0)

label(ax, 0.35, 4.56, "GEO", fs=13, w="bold", c=GEO_C, ha="left")
label(ax, 0.92, 4.56, "시료 정보·처리된 데이터", fs=11, c=TEXT2, ha="left")
label(ax, 4.55, 4.56, "SRA", fs=13, w="bold", c=SRA_C, ha="left")
label(ax, 5.10, 4.56, "시퀀서가 낸 raw read", fs=11, c=TEXT2, ha="left")


def card(x, y, w, h, code, desc, c, stack=False, ls="-"):
    if stack:
        for dx, dy in [(0.12, 0.10), (0.06, 0.05)]:
            box(ax, x + dx, y + dy, w, h, fc="white", ec=mix(c, 0.35), lw=1.0, z=1)
    box(ax, x, y, w, h, fc="white", ec=c, lw=1.5, z=2, ls=ls)
    label(ax, x + 0.16, y + h - 0.22, code, fs=13, w="bold", c=c, ha="left")
    if desc:
        label(ax, x + 0.16, y + 0.20, desc, fs=11, c=TEXT2, ha="left")


BW, BH = 2.40, 0.70
LX, RX = 0.35, 4.55
Y_STUDY, Y_SAMPLE, Y_THIRD, Y_FILE = 3.60, 2.45, 1.30, 0.30

# ── GEO column ─────────────────────────────────────────
card(LX, Y_STUDY, BW, BH, "GSE…", "시리즈 = 연구 하나", GEO_C)
card(LX, Y_SAMPLE, BW, BH, "GSM…", "시료 하나의 측정", GEO_C, stack=True)
card(LX, Y_THIRD, BW, BH, "GPL…", "플랫폼 (기기·칩)", GEO_C)
card(LX, Y_FILE, BW, BH, "보충 파일", "저자가 올린 행렬·추가 결과", GEO_C, ls="--")

cx = LX + BW / 2
arrow(ax, cx, Y_STUDY - 0.04, cx, Y_SAMPLE + BH + 0.14, color=GREY)
ax.plot([cx, cx], [Y_SAMPLE - 0.04, Y_THIRD + BH + 0.04], color=GREY, lw=1.4,
        ls=(0, (2, 2.5)), zorder=1)
label(ax, cx + 0.14, (Y_SAMPLE + Y_THIRD + BH) / 2, "어떤 기기로 쟀나", fs=11, c=TEXT2, ha="left")

# ── SRA column ─────────────────────────────────────────
card(RX, Y_STUDY, BW, BH, "SRP…", "study = 연구 하나", SRA_C)
card(RX, Y_SAMPLE, BW, BH, "SRX…", "experiment", SRA_C, stack=True)
GAP = 0.20
SW = (BW - GAP) / 2
srr_x = [RX, RX + SW + GAP]
for x in srr_x:
    card(x, Y_THIRD, SW, BH, "SRR…", "run", SRA_C)
card(RX, Y_FILE, BW, BH, "FASTQ", "raw read 파일", SRA_C, ls="--")

rcx = RX + BW / 2
arrow(ax, rcx, Y_STUDY - 0.04, rcx, Y_SAMPLE + BH + 0.14, color=GREY)
for x in srr_x:
    arrow(ax, rcx, Y_SAMPLE - 0.04, x + SW / 2, Y_THIRD + BH + 0.04, color=GREY)
    arrow(ax, x + SW / 2, Y_THIRD - 0.04, x + SW / 2, Y_FILE + BH + 0.04, color=GREY)
label(ax, srr_x[1] + SW / 2 + 0.18, Y_THIRD + BH + 0.22, "1개 이상", fs=11, c=TEXT2, ha="left")

# ── GEO ↔ SRA links ───────────────────────────────────
for y, txt, x0 in [(Y_STUDY + BH / 2, "같은 연구", LX + BW + 0.08),
                   (Y_SAMPLE + BH / 2, "같은 시료", LX + BW + 0.22)]:
    ax.add_patch(mp.FancyArrowPatch((x0, y), (RX - 0.08, y), arrowstyle="<|-|>",
                                    mutation_scale=11, color=TEXT2, lw=1.3,
                                    ls=(0, (4, 3)), shrinkA=0, shrinkB=0, zorder=3))
    label(ax, 3.80, y + 0.17, txt, fs=11, c=TEXT)

save(fig, os.path.join(HERE, "../../figures/A3_geo_sra.png"))
