"""Fig 2 - LOD / LOQ / ULOQ triptychs for the three MMCC comparisons. Runs entirely
from committed figures-of-merit CSVs (no raw data needed).

  2A gradient : 60SPD vs 100SPD          (asms diann_report)
  2B hardware : Ultra vs Ultra II        (250320 separate-search CURVES_pep)
  2C software : 60SPD report vs 60SPD_pr (asms report vs pr_matrix)

Each row is four axes: the three FOM densities, then a fourth panel carrying the
same three FOMs as box-and-whisker pairs. Whiskers are the usual 1.5 x IQR;
fliers are hidden because n is in the thousands and would swamp the boxes.
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, kde_pair, LIGHT, DARK, PRIMARY
from src.fom_io import load_config, load_fom, finite, repo_path

set_style()
cfg = load_config()
FOMS = ["LOD", "LOQ", "ULOQ"]
BOX_W, BOX_OFF = 0.30, 0.18          # box width, and the +/- offset of the two series

PANELS = {
    "FIG2A_gradient_triptych": (("60SPD", "bruker_60spd"), ("100SPD", "bruker_100spd")),
    "FIG2B_hardware_triptych": (("Ultra", "bruker_ultra"), ("Ultra II", "bruker_ultraII")),
    "FIG2C_software_triptych": (("report", "bruker_60spd"), ("pr_matrix", "bruker_60spd_pr")),
}

for stem, ((la, na), (lb, nb)) in PANELS.items():
    da, db = load_fom("main", na, cfg), load_fom("main", nb, cfg)
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3),
                             gridspec_kw={"width_ratios": [1, 1, 1, 1.3]})

    # --- the three density panels (unchanged; y shared among themselves only) ---
    for ax in axes[1:3]:
        ax.sharey(axes[0])
        ax.tick_params(labelleft=False)
    for ax, fom in zip(axes[:3], FOMS):
        kde_pair(ax, finite(da[fom]), finite(db[fom]), la, lb)
        ax.set_xlabel(fom)
        ax.set_xlim(-0.1, 1.1)

    # --- fourth panel: the same three FOMs as box-and-whisker pairs -------------
    axb = axes[3]
    for i, fom in enumerate(FOMS):
        for off, d, c in ((-BOX_OFF, da, LIGHT), (BOX_OFF, db, DARK)):
            v = finite(d[fom])
            if v.size == 0:
                continue
            bp = axb.boxplot(v, positions=[i + 1 + off], widths=BOX_W,
                             patch_artist=True, showfliers=False, whis=1.5)
            bp["boxes"][0].set(facecolor=c, alpha=0.3, edgecolor=c, linewidth=1.5)
            for part in ("whiskers", "caps"):
                for art in bp[part]:
                    art.set(color=c, linewidth=1.5)
            bp["medians"][0].set(color=PRIMARY, linewidth=1.5)
    axb.set_xticks(range(1, len(FOMS) + 1))
    axb.set_xticklabels(FOMS)
    axb.set_xlim(0.5, len(FOMS) + 0.5)
    axb.set_ylim(-0.1, 1.1)
    axb.set_ylabel("concentration")
    axb.grid(True, axis="y", alpha=0.3)
    axb.legend(handles=[Patch(facecolor=c, alpha=0.3, edgecolor=c, label=lab)
                        for c, lab in ((LIGHT, la), (DARK, lb))], fontsize=8)

    plt.tight_layout()
    out = repo_path(cfg["output"], f"{stem}.png")
    plt.savefig(out, dpi=1000, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)
