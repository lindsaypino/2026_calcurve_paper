"""Supplement - SIS Ultra (heavy-only) good/bad/ugly example calibration curves.

A three-panel triptych illustrating the range of SIS peptide-response curve quality
on the timsTOF Ultra: a clean quantifiable curve, a detectable-but-not-quantifiable
curve, and a pure-background peptide with no usable response. Each panel is annotated
with the fit R^2 (log and linear space), the median within-level replicate CV, and the
LOD/LOQ/ULOQ from the committed figures of merit. No SIS peptide reaches saturation, so
ULOQ is inf throughout (see docs/diagnostics.md).

Runs from committed artifacts only (output/sis_run/sis_melted.csv + fom_mnp0_msp2.csv);
the curve shape is refit with the pinned tool's fit_by_lmfit_yang. No raw MS data needed.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL, FOM_COLORS
from src.fom_io import load_config, load_tool, repo_path

PRIMARY, ACCENT = set_style()
cfg = load_config()
cl = load_tool()

SIS_DIR = repo_path("output", "sis_run")
melted = pd.read_csv(os.path.join(SIS_DIR, "sis_melted.csv"))
fom = pd.read_csv(os.path.join(SIS_DIR, "fom_mnp0_msp2.csv")).set_index("peptide")

C_FIT = PAL[0]                       # blue trilinear fit
C_LOD, C_LOQ = FOM_COLORS["LOD"], FOM_COLORS["LOQ"]   # orange / teal

# good / bad / ugly, chosen from the per-peptide fit-quality audit (R^2, CV)
PICKS = [
    ("LDFAVEDGENAQK", "Good", "clean trilinear, low scatter"),
    ("SGGLLQLWK", "Bad", "detectable, not quantifiable"),
    ("FVINYDYPNSSEDYIHR", "Ugly", "pure background, no response"),
]


def num(v):
    return pd.to_numeric(v, errors="coerce")


def r2(y, yhat):
    y = np.asarray(y, float); yhat = np.asarray(yhat, float)
    ss_res = np.sum((y - yhat) ** 2); ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan


fig, axes = plt.subplots(1, 3, figsize=(14.5, 5))
for ax, (pep, tier, note) in zip(axes, PICKS):
    sub = melted[melted.peptide == pep].sort_values("curvepoint")
    x = sub.curvepoint.to_numpy(float); y = sub.area.to_numpy(float)
    row = fom.loc[pep]
    LOD, LOQ, ULOQ = num(row.LOD), num(row.LOQ), num(row.ULOQ)

    cv = sub.groupby("curvepoint")["area"].agg(
        lambda s: s.std(ddof=1) / s.mean() if s.mean() else np.nan)
    med_cv = np.nanmedian(cv.values)

    res, _ = cl.fit_by_lmfit_yang(x, y, "auto")
    a = res.params["a"].value; b = res.params["b"].value
    c = res.params["c"].value; c_high = res.params["c_high"].value
    yhat = np.clip(a * x + b, c, c_high if np.isfinite(c_high) else np.inf)
    r2_lin = r2(y, yhat)
    pos = (y > 0) & (yhat > 0)
    r2_log = r2(np.log10(y[pos]), np.log10(yhat[pos])) if pos.sum() > 2 else np.nan

    xs = np.logspace(np.log10(x.min()), np.log10(x.max()), 300)
    ys = np.maximum(a * xs + b, c)
    if np.isfinite(c_high):
        ys = np.minimum(ys, c_high)
    ys = np.clip(ys, max(y.min() * 0.3, 1), y.max() * 3)
    ax.plot(xs, ys, "-", color=C_FIT, lw=2.4, zorder=3)
    ax.scatter(x, y, s=48, color=PRIMARY, zorder=5, edgecolor="white", linewidth=0.7)
    if np.isfinite(LOD):
        ax.axvline(LOD, color=C_LOD, lw=1.8)
    if np.isfinite(LOQ):
        ax.axvline(LOQ, color=C_LOQ, lw=1.8, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(f"{tier} — {pep}\n{note}", fontsize=11)
    ax.set_xlabel("concentration (per-peptide units)")
    ax.set_ylabel("heavy peak area")
    ax.grid(True, alpha=0.3)
    stats = (u"R²(log) %.2f\nR²(lin) %.2f\nmed CV %.0f%%\n"
             "LOD %s\nLOQ %s\nULOQ %s") % (
        r2_log, r2_lin, 100 * med_cv,
        ("%.3g" % LOD) if np.isfinite(LOD) else "inf",
        ("%.3g" % LOQ) if np.isfinite(LOQ) else "inf",
        ("%.3g" % ULOQ) if np.isfinite(ULOQ) else "inf")
    ax.text(0.03, 0.97, stats, transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.8", alpha=0.9))

handles = [Line2D([], [], color=PRIMARY, marker="o", ls="", label="replicates (n=3/level)"),
           Line2D([], [], color=C_FIT, lw=2.4, label="trilinear auto fit"),
           Line2D([], [], color=C_LOD, lw=1.8, label="LOD"),
           Line2D([], [], color=C_LOQ, lw=1.8, ls="--", label="LOQ")]
fig.legend(handles=handles, frameon=False, ncol=4, loc="lower center",
           bbox_to_anchor=(0.5, -0.02))
fig.suptitle("SIS Ultra (heavy-only) — good / bad / ugly example calibration curves")
fig.tight_layout(rect=(0, 0.05, 1, 0.96))

out = repo_path(cfg["output"], "SUPP_sis_curve_examples.png")
fig.savefig(out, dpi=300)
print("wrote", out)
