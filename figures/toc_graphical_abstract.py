"""Graphical abstract (ACS/JPR Table-of-Contents graphic).

Schematic only — no data, no specific results (ACS TOC rules). Sized to the ACS
maximum of 3.25 x 1.75 in; sans-serif type, nothing below 6.5 pt; colours are the
paper's own LOD/LOQ/ULOQ code from `src.style`, darkened to clear WCAG contrast.

Writes output/TOC_graphical_abstract.{svg,png} (PNG at 300 dpi).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.fom_io import load_config, repo_path

# --- palette -------------------------------------------------------------------
# hues from src/style.py (PRIMARY navy, talusbio blue, FOM_COLORS); the FOM hues
# are darkened here so text/markers clear WCAG 4.5:1 (text) / 3:1 (non-text).
NAVY = "#0C015B"
BLUE = "#0086BB"
GREY = "#5A5A66"
LOD, LOQ, ULOQ = "#C0491F", "#1F7A5E", "#A83E68"     # 4.98 / 5.25 / 5.91 : 1
S_A, S_A_TXT = "#64C0CA", "#1B6F7B"                  # method A (style.LIGHT) + label
S_B, S_B_TXT = "#0086BB", "#00658D"                  # method B (style.DARK)  + label
BAND = {"noise": "#EFEFF4", "linear": "#E4F2F8", "sat": "#EFEFF4"}

W, H = 3.25, 1.75                                    # ACS TOC maximum, inches
FS_HEAD, FS_SUB, FS_TINY = 9.0, 7.0, 6.5

matplotlib.rcParams.update({
    "svg.fonttype": "none",                          # keep text as text in the SVG
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.6,
})

fx = lambda v: v / W                                 # inches -> figure fraction
fy = lambda v: v / H

fig = plt.figure(figsize=(W, H))
fig.patch.set_facecolor("white")

AX_B, AX_T = 0.24, 1.30                              # panel bottom / top, inches
axL = fig.add_axes([fx(0.24), fy(AX_B), fx(1.17), fy(AX_T - AX_B)])
axR = fig.add_axes([fx(1.88), fy(AX_B), fx(1.20), fy(AX_T - AX_B)])

for ax in (axL, axR):
    ax.set_xlim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GREY)

# --- left panel: one peptide, trilinear fit ------------------------------------
XB1, XB2, NOISE, TOP = 0.22, 0.62, 0.09, 0.88
SLOPE = (TOP - NOISE) / (XB2 - XB1)
tri = lambda x: np.clip(NOISE + SLOPE * (np.asarray(x, float) - XB1), NOISE, TOP)

axL.set_ylim(0, 1.28)
axL.axvspan(0, XB1, color=BAND["noise"], lw=0)
axL.axvspan(XB1, XB2, color=BAND["linear"], lw=0)
axL.axvspan(XB2, 1, color=BAND["sat"], lw=0)

rng = np.random.default_rng(7)
levels = np.array([0.00, 0.05, 0.10, 0.16, 0.24, 0.33, 0.43, 0.54, 0.66, 0.80, 0.95])
px = np.repeat(levels, 3)
py = tri(px) + rng.normal(0, 1, px.size) * (0.022 + 0.055 * tri(px))
axL.scatter(px, np.clip(py, 0.01, 1.02), s=3.2, color=BLUE, alpha=0.75, lw=0, zorder=3)

xs = np.linspace(0, 1, 200)
axL.plot(xs, tri(xs), color=NAVY, lw=1.1, zorder=4, solid_capstyle="round")

for x, c in ((0.28, LOD), (0.40, LOQ), (XB2, ULOQ)):
    axL.vlines(x, 0, tri(x), color=c, lw=0.7, ls=(0, (2.2, 1.6)), zorder=5)
    axL.plot([x], [tri(x)], "o", ms=3.0, mfc=c, mec="white", mew=0.6, zorder=6)

for i, (lab, c) in enumerate((("LOD", LOD), ("LOQ", LOQ), ("ULOQ", ULOQ))):
    y = 1.15 - 0.145 * i
    axL.plot([0.035, 0.105], [y, y], color=c, lw=1.4, solid_capstyle="butt", zorder=6)
    axL.text(0.135, y, lab, color=c, fontsize=FS_TINY, va="center", ha="left", zorder=6)

axL.set_ylabel("signal", fontsize=FS_SUB, color=GREY, labelpad=2)
axL.set_xlabel("concentration", fontsize=FS_SUB, color=GREY, labelpad=2)

# --- right panel: LOQ distributions compare two methods ------------------------
gauss = lambda x, mu, sd: np.exp(-0.5 * ((x - mu) / sd) ** 2)
gx = np.linspace(0, 1, 400)
gB = gauss(gx, 0.34, 0.115) * 0.80          # lower LOQ  -> better method
gA = gauss(gx, 0.63, 0.135) * 0.68          # higher LOQ

axR.set_ylim(0, 1.28)
for g, face, edge in ((gA, S_A, S_A), (gB, S_B, S_B)):
    axR.fill_between(gx, 0, g, color=face, alpha=0.30, lw=0)
    axR.plot(gx, g, color=edge, lw=1.1, zorder=3)

axR.annotate("", xy=(0.34, 1.10), xytext=(0.63, 1.10),
             arrowprops=dict(arrowstyle="-|>,head_width=0.10,head_length=0.22",
                             color=NAVY, lw=0.7, shrinkA=0, shrinkB=0), zorder=5)
axR.text(0.485, 1.15, "lower LOQ", ha="center", va="baseline",
         fontsize=FS_TINY, color=NAVY)

axR.text(0.34, 0.88, "method B", color=S_B_TXT, fontsize=FS_TINY, ha="center",
         va="baseline", zorder=6)
axR.text(0.63, 0.76, "method A", color=S_A_TXT, fontsize=FS_TINY, ha="center",
         va="baseline", zorder=6)

axR.set_ylabel("peptides", fontsize=FS_SUB, color=GREY, labelpad=2)
axR.set_xlabel("LOQ", fontsize=FS_SUB, color=GREY, labelpad=2)

# --- annotation ----------------------------------------------------------------
fig.text(0.5, fy(1.625), "Bounding the quantifiable range",
         ha="center", va="baseline", fontsize=FS_HEAD, color=NAVY, weight="bold")

fig.text(fx(0.24 + 1.17 / 2), fy(1.485), "one peptide", ha="center", va="baseline",
         fontsize=FS_SUB, color=GREY, style="italic")
fig.text(fx(1.88 + 1.20 / 2), fy(1.485), "every peptide", ha="center", va="baseline",
         fontsize=FS_SUB, color=GREY, style="italic")  # population, two methods

for x, lab in ((0.11, "noise"), (0.42, "linear"), (0.81, "saturation")):
    fig.text(fx(0.24 + x * 1.17), fy(1.365), lab, ha="center", va="baseline",
             fontsize=FS_TINY, color=GREY)
fig.text(fx(1.88 + 1.20 / 2), fy(1.365), "gradient · instrument · software",
         ha="center", va="baseline", fontsize=FS_TINY, color=GREY)

fig.add_artist(FancyArrow(fx(1.46), fy(0.77), fx(0.22), 0, transform=fig.transFigure,
                          width=fy(0.012), head_width=fy(0.055), head_length=fx(0.08),
                          length_includes_head=True, color=NAVY, lw=0))

cfg = load_config()
for ext, kw in (("svg", {}), ("png", {"dpi": 300})):
    out = repo_path(cfg["output"], f"TOC_graphical_abstract.{ext}")
    fig.savefig(out, facecolor="white", **kw)
    print("wrote", out)
