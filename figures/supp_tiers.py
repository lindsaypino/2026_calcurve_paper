"""Supplement - under min_noise_points=0, split LOD/LOQ into two tiers per DIA-NN
dataset: resolved noise floor (stndev_noise finite) vs detection-floored (stndev_noise
NaN; LOD pinned to the peptide's lowest detected dilution -> the 'lumps'). Runs from
committed FOM (supp_mnp0)."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, LIGHT, DARK
from src.fom_io import load_config, load_fom, repo_path

set_style()
cfg = load_config()
DS = [("ultra", "Ultra"), ("60spd", "60 SPD"), ("100spd", "100 SPD")]


def tiers(d, fom):
    v = pd.to_numeric(d[fom], errors="coerce")
    sd = pd.to_numeric(d["stndev_noise"], errors="coerce")
    fin = np.isfinite(v) & (v > 0) & (v <= 1.0)
    return np.log10(v[fin & np.isfinite(sd)]), np.log10(v[fin & sd.isna()])


fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for j, (name, nice) in enumerate(DS):
    d = load_fom("mnp0", f"bruker_{name}", cfg)
    for i, fom in enumerate(["LOD", "LOQ"]):
        ax = axes[i, j]
        res, det = tiers(d, fom)
        if det.size > 2:
            sns.kdeplot(det, ax=ax, color=DARK, fill=True, alpha=0.3, linewidth=2,
                        label=f"detection-floored (n={det.size:,})")
        if res.size > 2:
            sns.kdeplot(res, ax=ax, color=LIGHT, fill=True, alpha=0.3, linewidth=2,
                        label=f"resolved noise floor (n={res.size:,})")
        ax.set_xlabel(f"{fom} (log)"); ax.set_xlim(-3, 0)
        ax.set_xticks(range(-3, 1)); ax.set_xticklabels([f"$10^{{{t}}}$" for t in range(-3, 1)])
        ax.grid(True, alpha=0.3); ax.set_title(nice)

# one shared legend above the panels (generic labels — per-panel n varies), so it
# never overlaps a panel's data or ticks
from matplotlib.patches import Patch
legend_handles = [
    Patch(facecolor=DARK, alpha=0.3, edgecolor=DARK, linewidth=2, label="detection-floored"),
    Patch(facecolor=LIGHT, alpha=0.3, edgecolor=LIGHT, linewidth=2, label="resolved noise floor"),
]
fig.legend(handles=legend_handles, loc="upper center", ncol=2, frameon=False,
           bbox_to_anchor=(0.5, 1.02))
plt.tight_layout(rect=(0, 0, 1, 0.95))
out = repo_path(cfg["output"], "SUPP_lod_loq_tiers.png")
plt.savefig(out, dpi=1000, bbox_inches="tight")
print("wrote", out)
