"""Example ULOQ calibration curves from each dataset (Exploris DIA, IL15 PRM,
Bruker Ultra II DIA) on LINEAR axes, with the raw measured points and the trilinear
auto fit + LOD/ULOQ markers. No log scaling. NEEDS RAW DATA (config.yaml -> raw:) and
the Ultra II FOM (data/figuresofmerit/main/bruker_ultraII.csv) for peptide selection.
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL
from src.fom_io import load_config, load_tool, load_fom, raw, repo_path
from src.prep_curves_pep import reshape

PRIMARY, ACCENT = set_style()
cl = load_tool()
cfg = load_config()


def fit_pep(sub, mnp):
    sub = sub.sort_values("curvepoint")
    x = sub.curvepoint.to_numpy(float); y = sub.area.to_numpy(float)
    res, _ = cl.fit_by_lmfit_yang(x, y, "auto")
    a = res.params["a"].value; b = res.params["b"].value
    c = res.params["c"].value; c_high = res.params["c_high"].value
    LOD, sd = cl.calculate_lod(np.array([0.0, c, a, b]), sub, 2, mnp, 1, x, "auto")
    ULOQ = cl.calculate_uloq(a, b, c_high, sub, 2, 2)
    return x, y, a, b, c, c_high, LOD, ULOQ


def panel(ax, x, y, a, b, c, c_high, LOD, ULOQ, title, xlabel):
    xs = np.linspace(0, x.max(), 400)
    ys = np.minimum(np.maximum(a * xs + b, c), c_high)
    ax.axvspan(LOD, ULOQ, color=ACCENT, alpha=0.12, lw=0)
    ax.plot(xs, ys, "-", color=PAL[0], lw=2.2, zorder=3)
    if np.isfinite(c_high):
        ax.axhline(c_high, color="0.6", ls=":", lw=1.3)
    ax.scatter(x, y, s=55, color=PRIMARY, zorder=5, edgecolor="white", linewidth=0.7)
    ax.axvline(LOD, color=PAL[1], lw=1.8)
    ax.axvline(ULOQ, color=PAL[3], lw=1.8)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(xlabel); ax.set_ylabel("peak area")
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.set_xlim(-x.max() * 0.02, x.max() * 1.02)
    ax.text(0.97, 0.05, "LOD %.3g\nULOQ %.3g" % (LOD, ULOQ), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.85))


exp = cl.read_input(raw("exploris_elib", cfg), repo_path("data/maps/exploris_filenamemap.csv"))
il15 = cl.read_input(raw("il15_elib", cfg), repo_path("data/maps/il15_filenamemap.csv"))
brk = reshape(raw("ultraII_curves_pep", cfg))          # melted (peptide, curvepoint, area), human-only
brk_fom = load_fom("main", "bruker_ultraII", cfg)

# pick two clean, abundant Bruker peptides with a finite LOD and ULOQ, excluding the
# sequences already shown in the Exploris row
used = {"VFLENVIR", "YPIEHGIITNWDDMEK", "AGLQFPVGR"}
brk_u = brk_fom[np.isfinite(pd.to_numeric(brk_fom.LOD, errors="coerce"))
                & np.isfinite(pd.to_numeric(brk_fom.ULOQ, errors="coerce"))].copy()
brk_u = brk_u[~brk_u.peptide.str.replace(r"_[\d.]+$", "", regex=True).isin(used)]
brk_u["ymax"] = brk_u.peptide.map(brk.groupby("peptide")["area"].max())
brk_peps = list(brk_u.sort_values("ymax", ascending=False).peptide.head(2))

ROWS = [
    ("Exploris 480 DIA", exp, 2, "% analyte", ["VFLENVIR", "YPIEHGIITNWDDMEK"]),
    ("IL15/IL2 PRM", il15, 0, "% IL", ["RVAEDDEDDDVDTK", "VAEDDEDDDVDTKK"]),
    ("Bruker Ultra II DIA", brk, 2, "fraction", brk_peps),
]

fig, axes = plt.subplots(3, 2, figsize=(11, 12))
for r, (label, df, mnp, xlabel, peps) in enumerate(ROWS):
    for cc, pep in enumerate(peps):
        x, y, a, b, c, c_high, LOD, ULOQ = fit_pep(df[df.peptide == pep], mnp)
        disp = re.sub(r"_(\d+)\.0$", r" +\1", pep)  # "SEQ_2.0" -> "SEQ +2"
        short = disp if len(disp) <= 22 else disp[:21] + "…"
        panel(axes[r, cc], x, y, a, b, c, c_high, LOD, ULOQ, f"{short}  ({label})", xlabel)

plt.tight_layout()
out = repo_path(cfg["output"], "ULOQ_examples_linear_raw.png")
plt.savefig(out, dpi=200, bbox_inches="tight")
print("Bruker peptides used:", brk_peps)
print("wrote", out)
