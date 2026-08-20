"""Diagnostic - all 26 SIS/AQUA peptide-response curves on the Ultra, log-log, one
column per spike-level group (A..E). The question this figure answers: does the top of
any curve flatten into a saturation plateau (which would give a ULOQ), or is it still
climbing? Each panel carries a dashed slope-1 reference (ideal linear response) anchored
at the second curve point, and the measured top step ratio -- nominal is 3.0 for the
3-fold dilution, ~1.0 would mean a genuine plateau.

NEEDS RAW DATA (config.yaml -> raw: sis_ultra_export / sis_ultra_sky).
Optional: --fom <figuresofmerit.csv> overlays LOD/LOQ/ULOQ markers.
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL, FOM_COLORS
from src.fom_io import load_config, raw, repo_path
from src.prep_skyline_sis import reshape

set_style()
cfg = load_config()

ap = argparse.ArgumentParser()
ap.add_argument("--fom", default=None, help="optional figuresofmerit.csv for LOD/LOQ/ULOQ markers")
ap.add_argument("--out", default=None)
args = ap.parse_args()

tidy = reshape(raw("sis_ultra_export", cfg), raw("sis_ultra_sky", cfg))
tidy = tidy[tidy["area"] > 0]

fom = None
if args.fom and os.path.exists(args.fom):
    fom = pd.read_csv(args.fom).set_index("peptide")

GROUPS = ["group_A", "group_B", "group_C", "group_D", "group_E"]
GCOLOR = dict(zip(GROUPS, [PAL[4], PAL[3], PAL[2], PAL[1], PAL[0]]))
nrow = max(tidy.groupby("group")["peptide"].nunique())

fig, axes = plt.subplots(nrow, len(GROUPS), figsize=(4.0 * len(GROUPS), 2.7 * nrow),
                         sharex=False, sharey=False)

for j, grp in enumerate(GROUPS):
    sub = tidy[tidy["group"] == grp]
    mult = sub["multiplier"].iloc[0]
    peps = sorted(sub["peptide"].unique())
    for i in range(nrow):
        ax = axes[i, j]
        if i >= len(peps):
            ax.axis("off")
            continue
        pep = peps[i]
        p = sub[sub["peptide"] == pep]
        med = p.groupby("eff_conc")["area"].median().sort_index()

        ax.scatter(p["eff_conc"], p["area"], s=26, color=GCOLOR[grp],
                   alpha=0.85, zorder=3, edgecolor="none")
        ax.plot(med.index, med.values, color=GCOLOR[grp], lw=1.4, alpha=0.75, zorder=2)

        # slope-1 (ideal linear) reference anchored at the 2nd point
        if len(med) >= 2:
            x0, y0 = med.index[1], med.values[1]
            xr = np.array([med.index[0], med.index[-1]], dtype=float)
            ax.plot(xr, y0 * (xr / x0), ls="--", lw=1.1, color="0.45", zorder=1)

        # top step ratio: the saturation tell
        ratio = med.values[-1] / med.values[-2] if len(med) >= 2 and med.values[-2] > 0 else np.nan
        ax.text(0.04, 0.93, f"top step {ratio:.2f}x", transform=ax.transAxes,
                fontsize=9, va="top", ha="left",
                bbox=dict(fc="white", ec="0.8", lw=0.6, alpha=0.9, pad=1.8))

        if fom is not None and pep in fom.index:
            for key in ("LOD", "LOQ", "ULOQ"):
                v = pd.to_numeric(pd.Series([fom.loc[pep, key]]), errors="coerce").iloc[0]
                if np.isfinite(v):
                    ax.axvline(v, color=FOM_COLORS[key], lw=1.3, ls=":", zorder=1)

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(pep, fontsize=9, pad=3)
        ax.tick_params(labelsize=8)
        ax.grid(True, which="major", alpha=0.25)
        if i == 0:
            ax.annotate(f"{grp.replace('group_','group ')}  (x{mult:,.0f})",
                        xy=(0.5, 1.30), xycoords="axes fraction", ha="center",
                        fontsize=13, fontweight="bold", color=GCOLOR[grp])

handles = [Line2D([], [], ls="--", color="0.45", label="ideal linear (slope 1)")]
if fom is not None:
    handles += [Line2D([], [], ls=":", color=FOM_COLORS[k], label=k) for k in ("LOD", "LOQ", "ULOQ")]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=11)

fig.suptitle("SIS peptide-response curves, timsTOF Ultra — no top-end plateau "
             "(top step stays well above 1, so no ULOQ)", fontsize=15, y=0.995)
fig.supxlabel("effective concentration (analyte conc × per-peptide multiplier)", fontsize=12, y=0.035)
fig.supylabel("heavy Total Area Fragment", fontsize=12)
fig.tight_layout(rect=(0.012, 0.055, 1, 0.975))

out = args.out or repo_path(cfg["output"], "diag_sis_curves.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=180)
print("wrote", out)
