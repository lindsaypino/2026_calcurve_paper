"""Supplement - old (default min_noise_points=2) vs new (min_noise_points=0) LOD/LOQ
for the DIA-NN datasets. Shows the settings correction. Runs from committed FOM
(legacy_mnp2 = old code, supp_mnp0 = new code)."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, LIGHT, DARK
from src.fom_io import load_config, load_fom, finite, repo_path

set_style()
cfg = load_config()
DS = [("ultra", "Ultra"), ("60spd", "60 SPD"), ("100spd", "100 SPD")]


def flog(v):
    v = v[(v > 0) & (v <= 1.0)]
    return np.log10(v)


fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for j, (name, nice) in enumerate(DS):
    old = load_fom("legacy", f"bruker_{name}", cfg)
    new = load_fom("mnp0", f"bruker_{name}", cfg)
    for i, fom in enumerate(["LOD", "LOQ"]):
        ax = axes[i, j]
        for d, c, lab in [(old, LIGHT, "old (mnp=2)"), (new, DARK, "new (mnp=0)")]:
            lg = flog(finite(d[fom]))
            n = len(finite(d[fom]))
            if lg.size > 2:
                sns.kdeplot(lg, ax=ax, color=c, fill=True, alpha=0.3, linewidth=2, label=f"{lab}: {n:,}")
        ax.set_xlabel(f"{fom} (log)"); ax.set_xlim(-3, 0)
        ax.set_xticks(range(-3, 1)); ax.set_xticklabels([f"$10^{{{t}}}$" for t in range(-3, 1)])
        ax.grid(True, alpha=0.3); ax.set_title(nice)
        if i == 0 and j == 2:
            ax.legend()
plt.tight_layout()
out = repo_path(cfg["output"], "SUPP_rerun_before_after.png")
plt.savefig(out, dpi=1000, bbox_inches="tight")
print("wrote", out)
