"""Fig 2 - LOD / LOQ / ULOQ triptychs for the three MMCC comparisons. Runs entirely
from committed figures-of-merit CSVs (no raw data needed).

  2A gradient : 60SPD vs 100SPD          (asms diann_report)
  2B hardware : Ultra vs Ultra II        (250320 separate-search CURVES_pep)
  2C software : 60SPD report vs 60SPD_pr (asms report vs pr_matrix)
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, kde_pair
from src.fom_io import load_config, load_fom, finite, repo_path

set_style()
cfg = load_config()
FOMS = ["LOD", "LOQ", "ULOQ"]

PANELS = {
    "FIG2A_gradient_triptych": (("60SPD", "bruker_60spd"), ("100SPD", "bruker_100spd")),
    "FIG2B_hardware_triptych": (("Ultra", "bruker_ultra"), ("Ultra II", "bruker_ultraII")),
    "FIG2C_software_triptych": (("report", "bruker_60spd"), ("pr_matrix", "bruker_60spd_pr")),
}

for stem, ((la, na), (lb, nb)) in PANELS.items():
    da, db = load_fom("main", na, cfg), load_fom("main", nb, cfg)
    fig, axes = plt.subplots(1, 3, figsize=(9, 3), sharey=True)
    for ax, fom in zip(axes, FOMS):
        kde_pair(ax, finite(da[fom]), finite(db[fom]), la, lb)
        ax.set_xlabel(fom)
        ax.set_xlim(-0.1, 1.1)
    axes[-1].legend(fontsize=8)
    plt.tight_layout()
    out = repo_path(cfg["output"], f"{stem}.png")
    plt.savefig(out, dpi=1000, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
