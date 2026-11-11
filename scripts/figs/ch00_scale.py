# ch00 §2 크기의 사다리: 1 nm–100 µm 로그 눈금 위의 DNA·단백질·리보솜·미토콘드리아·핵·림프구·Visium spot  (redraw of figures/ch00_scale.svg)
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from medfig import *

W, H = 7.4, 3.0
fig, ax = new(W, H)

X0, DEC = 0.8, 1.2          # x of 1 nm, inches per decade
Y = 1.3                     # axis height

def xp(nm):
    return X0 + DEC * np.log10(nm)

AXC = "#9AA1AB"
BIO = PRIMARY               # 세포 안의 구조·세포
SPOT = CAT[1]               # 실험 장치(Visium spot)

# ── axis ────────────────────────────────────────────────────────────────
ax.plot([xp(1), xp(1e5)], [Y, Y], color=AXC, lw=1.4, zorder=2, solid_capstyle="butt")
for d in range(6):
    x = xp(10 ** d)
    ax.plot([x, x], [Y - 0.08, Y + 0.08], color=AXC, lw=1.4, zorder=2)
    if d < 5:
        for m in range(2, 10):
            xm = xp(m * 10 ** d)
            ax.plot([xm, xm], [Y - 0.04, Y + 0.04], color=AXC, lw=0.9, zorder=2)
ticks = ["1 nm", "10 nm", "100 nm", "1 µm", "10 µm", "100 µm"]
for d, t in enumerate(ticks):
    label(ax, xp(10 ** d), Y - 0.14, t, fs=11, c=TEXT2, va="top")

# ── objects ─────────────────────────────────────────────────────────────
LEAD = dict(color=AXC, lw=1.0, ls=(0, (2, 2)), zorder=1)

def point(nm, c):
    ax.plot([xp(nm)], [Y], marker="o", ms=9, mfc=c, mec="white", mew=1.2, zorder=4)
    return xp(nm)

def span(a, b, c):
    ax.plot([xp(a), xp(b)], [Y, Y], color=c, lw=7, solid_capstyle="round", zorder=4)
    return (xp(a) + xp(b)) / 2

def above(x, L, name, size, c=TEXT):
    ax.plot([x, x], [Y + 0.1, L - 0.13], **LEAD)
    label(ax, x, L, size, fs=11, c=TEXT2)
    label(ax, x, L + 0.26, name, fs=12, w="bold", c=c)

def below(x, L, name, size, c=TEXT):
    ax.plot([x, x], [Y - 0.1, L + 0.14], **LEAD)
    label(ax, x, L, name, fs=12, w="bold", c=c)
    label(ax, x, L - 0.26, size, fs=11, c=TEXT2)

# DNA 이중나선: 지름 약 2 nm
x = point(2, BIO);            above(x, 2.30, "DNA 이중나선", "지름 약 2 nm")
# 단백질 하나: 약 4 nm
x = point(4, BIO);            below(x, 0.62, "단백질 하나", "약 4 nm")
# 리보솜: 20–30 nm
x = span(20, 30, BIO);        above(x, 1.78, "리보솜", "20–30 nm")
# 미토콘드리아: 길이 1–2 µm (leader from the left part of the bar, clear of the "1 µm" tick label)
span(1000, 2000, BIO);        below(xp(1180), 0.62, "미토콘드리아", "길이 1–2 µm")
# 핵: 약 5 µm
x = point(5000, BIO);         below(x, 0.62, "핵", "약 5 µm")
# 림프구: 약 7 µm
x = point(7000, BIO);         above(x, 1.78, "림프구", "약 7 µm")
# Visium spot: 지름 55 µm
x = point(55000, SPOT);       above(x, 1.78, "Visium spot", "지름 55 µm", c=SPOT)

save(fig, os.path.join(HERE, "../../figures/ch00_scale.png"))
