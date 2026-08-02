"""Fig 1B - LOD / LOQ / ULOQ density distributions for the Exploris DIA dataset
(the 'auto model works across instrumentation' panel). Runs from committed FOM."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL
from src.fom_io import load_config, load_fom, finite, repo_path

set_style()
cfg = load_config()
COLORS = {"LOD": PAL[0], "LOQ": PAL[1], "ULOQ": PAL[2]}
MAXC = 100.0  # Exploris curve tops out at 100%

d = load_fom("main", "exploris_dia", cfg)
fig, ax = plt.subplots(figsize=(6, 4))
for fom in ["LOD", "LOQ", "ULOQ"]:
    v = finite(d[fom]); v = v[(v > 0) & (v <= MAXC)]
    if v.size > 2:
        sns.kdeplot(np.log10(v), ax=ax, color=COLORS[fom], fill=True, alpha=0.3,
                    linewidth=2, label=f"{fom} (n={v.size:,})")
lo, hi = -3, np.log10(MAXC)
ax.set_xlim(lo, hi)
ticks = range(int(lo), int(np.ceil(hi)) + 1)
ax.set_xticks(list(ticks)); ax.set_xticklabels([f"$10^{{{t}}}$" for t in ticks])
ax.set_xlabel("concentration (log)"); ax.grid(True, alpha=0.3); ax.legend()
plt.tight_layout()
out = repo_path(cfg["output"], "FIG1B_exploris_fom_distributions.png")
plt.savefig(out, dpi=1000, bbox_inches="tight")
print("wrote", out)
