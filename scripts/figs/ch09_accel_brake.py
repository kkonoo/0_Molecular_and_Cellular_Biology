# ch09 §4 가속기(종양유전자)와 브레이크(종양억제유전자): 두 벌 중 몇 벌이 고장 나야 효과가 나는지  (redraw of figures/ch09_accel_brake.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.55, 4.25
fig, ax = new(W, H)

MUT = "#4B5563"          # 변이 allele — 두 줄 모두 같은 색
BAR = "#E6E9EE"          # allele 막대 (유전자 바깥 부분)
COLS = [2.88, 4.63, 6.38]
CW = 1.66                # 칸 너비
HEADS = ["두 벌 모두 정상", "한 벌에 변이", "두 벌 모두 변이"]

ROWS = [
    dict(y0=2.2, c=CAT[0], name="가속기", kind="종양유전자", en="(oncogene)",
         genes="KRAS · MYC",
         cells=[("신호 있을 때만 켜짐", None, False),
                ("계속 켜짐", "한 벌로 충분", True),
                ("계속 켜짐", None, True)]),
    dict(y0=0.35, c=CAT[1], name="브레이크", kind="종양억제유전자", en="(tumor suppressor gene)",
         genes="TP53 · RB1",
         cells=[("브레이크 작동", None, False),
                ("아직 작동", "첫 번째 hit", False),
                ("브레이크 없음", "두 번째 hit", True)]),
]
RH = 1.65

# column headers
for cx, h in zip(COLS, HEADS):
    label(ax, cx, ROWS[0]["y0"] + RH + 0.22, h, fs=12, w="bold")

def allele(cx, y, fc, mutant):
    box(ax, cx - 0.68, y - 0.07, 1.36, 0.14, fc=BAR, ec="none", r=0.07, z=2)
    box(ax, cx + 0.02, y - 0.13, 0.46, 0.26, fc=MUT if mutant else fc, ec="none", r=0.04, z=3)
    if mutant:
        label(ax, cx + 0.25, y, "변이", fs=11, w="bold", c="white", z=4)

for row in ROWS:
    y0, c = row["y0"], row["c"]
    box(ax, 0.05, y0, W - 0.15, RH, fc=FAINT, ec="none", r=0.1, z=0)

    # left labels
    label(ax, 0.25, y0 + 1.3, row["name"], fs=13, w="bold", c=c, ha="left")
    label(ax, 0.25, y0 + 0.97, row["kind"], fs=11, ha="left")
    label(ax, 0.25, y0 + 0.72, row["en"], fs=11, c=TEXT2, ha="left")
    label(ax, 0.25, y0 + 0.36, row["genes"], fs=11, w="bold", ha="left")

    for i, (cx, (main, sub, effect)) in enumerate(zip(COLS, row["cells"])):
        if effect:   # 효과가 난 칸: 줄 색의 옅은 면
            box(ax, cx - CW / 2, y0 + 0.1, CW, RH - 0.2, fc=mix(c, 0.13), ec="none", r=0.08, z=1)
        allele(cx, y0 + 1.28, c, i >= 1)
        allele(cx, y0 + 0.93, c, i >= 2)
        label(ax, cx, y0 + 0.54, main, fs=12, w="bold" if effect else "regular",
              c=TEXT if effect else TEXT2)
        if sub:
            label(ax, cx, y0 + 0.27, sub, fs=11, c=TEXT2)

save(fig, os.path.join(HERE, "../../figures/ch09_accel_brake.png"))
