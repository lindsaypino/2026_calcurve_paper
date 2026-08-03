"""Fig 1B - example MMCC auto-model curve fits across instrumentation / curve
constructions. Top row = Exhausted CD8 dimethyl-background DIA (Thermo Exploris 480);
bottom row = IL15/IL2-stimulated CD8 PRM. Refits raw curve data with the pinned tool
so panels match the production figures-of-merit. NEEDS RAW DATA (config.yaml -> raw:).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL
from src.fom_io import load_config, load_tool, raw, repo_path

PRIMARY, ACCENT = set_style()
cl = load_tool()
cfg = load_config()

C_DATA, C_FIT, C_NOISE = PRIMARY, PAL[0], "#9aa7b0"
C_LOD, C_LOQ, C_ULOQ, C_BAND = PAL[1], PAL[2], PAL[3], ACCENT
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300})

DATASETS = {
    "Exploris 480 DIA": dict(
        elib=raw("exploris_elib", cfg), cmap=repo_path("data/maps/exploris_filenamemap.csv"),
        min_noise_points=2, xlabel="% analyte",
        peptides=["VFLENVIR", "YPIEHGIITNWDDMEK", "AGLQFPVGR"]),
    "IL15/IL2 CD8 PRM": dict(
        elib=raw("il15_elib", cfg), cmap=repo_path("data/maps/il15_filenamemap.csv"),
        min_noise_points=0, xlabel="% IL",
        peptides=["RVAEDDEDDDVDTK", "TEIVPLFTNLASDEQDSVR", "VAEDDEDDDVDTKK"]),
}


def process(sub, min_noise_points):
    sub = sub.sort_values("curvepoint")
    x = sub.curvepoint.to_numpy(float); y = sub.area.to_numpy(float)
    res, _ = cl.fit_by_lmfit_yang(x, y, "auto")
    a = res.params["a"].value; b = res.params["b"].value; c = res.params["c"].value
    c_high = res.params["c_high"].value if "c_high" in res.params else np.inf
    mp = np.array([0.0, c, a, b])
    LOD, std_noise = cl.calculate_lod(mp, sub, 2, min_noise_points, 1, x, "auto")
    ULOQ = cl.calculate_uloq(a, b, c_high, sub, 2, 2)
    if np.isfinite(ULOQ) and np.isfinite(LOD) and ULOQ <= LOD:
        ULOQ = np.inf
    mp = np.append(mp, [LOD, std_noise])
    LOQ = np.inf
    if np.isfinite(LOD):
        upper = min(ULOQ, x.max()) if np.isfinite(ULOQ) else x.max()
        if upper <= LOD:
            upper = x.max()
        xi = np.linspace(LOD, upper, 100)
        boot_model = "trilinear" if np.isfinite(c_high) else "bilinear"
        bdf = cl.bootstrap_many(sub, new_x=xi, num_bootreps=200, model=boot_model)
        LOQ = cl.calculate_loq(np.append(mp, np.inf), bdf, 0.2)
    return dict(x=x, y=y, a=a, b=b, c=c, c_high=c_high, LOD=LOD, LOQ=LOQ, ULOQ=ULOQ)


def plot_panel(ax, r, title, xlabel):
    x, y = r["x"], r["y"]
    xmin_nz = np.nanmin(np.where(x > 0, x, np.nan))
    xplot = np.where(x > 0, x, xmin_nz / 3.0)
    xs = np.logspace(np.log10(xmin_nz / 3.0), np.log10(np.nanmax(x)), 400)
    ys = np.minimum(np.maximum(r["a"] * xs + r["b"], r["c"]), r["c_high"])
    lo = r["LOQ"] if np.isfinite(r["LOQ"]) else r["LOD"]
    hi = r["ULOQ"] if np.isfinite(r["ULOQ"]) else np.nanmax(x)
    if np.isfinite(lo) and hi > lo:
        ax.axvspan(lo, hi, color=C_BAND, alpha=0.07, lw=0)
    ax.plot(xs, ys, "-", color=C_FIT, lw=2.2, zorder=3)
    if np.isfinite(r["c_high"]):
        ax.axhline(r["c_high"], color=C_NOISE, ls=":", lw=1.2, zorder=1)
    ax.scatter(xplot, y, s=48, color=C_DATA, zorder=5, edgecolor="white", linewidth=0.6)
    for val, col in [(r["LOD"], C_LOD), (r["LOQ"], C_LOQ), (r["ULOQ"], C_ULOQ)]:
        if np.isfinite(val) and val > 0:
            ax.axvline(val, color=col, lw=1.6, zorder=4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel(xlabel); ax.set_ylabel("peak area")
    ax.grid(True, which="major", alpha=0.15)

    def f(v):
        return "n/a" if not np.isfinite(v) else ("%.3g" % v)
    ax.text(0.03, 0.97, f"LOD {f(r['LOD'])}\nLOQ {f(r['LOQ'])}\nULOQ {f(r['ULOQ'])}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cccccc", alpha=0.85))


fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for row, (dsname, ds) in enumerate(DATASETS.items()):
    df = cl.read_input(ds["elib"], ds["cmap"])
    for col, pep in enumerate(ds["peptides"]):
        r = process(df[df.peptide == pep], ds["min_noise_points"])
        short = pep if len(pep) <= 20 else pep[:19] + "…"
        plot_panel(axes[row, col], r, f"{short}\n{dsname}", ds["xlabel"])

handles = [
    Line2D([0], [0], color=C_FIT, lw=2.2, label="auto fit (noise–linear–saturation)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=C_DATA, markersize=8, label="curve points"),
    Line2D([0], [0], color=C_LOD, lw=1.6, label="LOD"),
    Line2D([0], [0], color=C_LOQ, lw=1.6, label="LOQ"),
    Line2D([0], [0], color=C_ULOQ, lw=1.6, label="ULOQ"),
]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=[0, 0.03, 1, 1])
out = repo_path(cfg["output"], "FIG1B_auto_model_examples.png")
fig.savefig(out, bbox_inches="tight")
print("wrote", out)
