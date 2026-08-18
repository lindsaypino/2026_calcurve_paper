"""Supplement - where the dilution points sit determines what LOD/LOQ you can measure.

Two halves:

  simulated (A, B)  the same measurement budget (42 injections) is spent on four
                    different dilution designs. Ground truth is known, so each
                    design has its own achievable LOQ and its own estimation
                    spread. Log spacing with replicates wins; a linear series
                    spends most of its runs where nothing is happening.

  real data (C, D)  the tool sample dataset, thinned to coarser designs, and read
                    on different CV grids. Dropping the dense low end moves LOD and
                    LOQ for real peptides; and reading the CV on a uniform grid --
                    what the tool does -- skips over the region where the crossing
                    actually is, so peptides get reported as having no crossing.

Runs from the pinned submodule and its committed sample data; no raw MS data
needed. Caches to output/ (gitignored); pass --force to recompute.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL, PRIMARY
from src.fom_io import load_config, repo_path, load_tool

set_style()
cfg = load_config()
calc = load_tool(cfg)

TOOL_DIR = os.path.dirname(os.path.dirname(repo_path(cfg["tool"]["calculate_loq"])))
CURVE = os.path.join(TOOL_DIR, "data", "one_protein.csv")
CMAP = os.path.join(TOOL_DIR, "data", "filename2samplegroup_map.csv")
CACHE = repo_path(cfg["output"], "_cache_curve_spacing.npz")

CV_THRESH, EPS = 0.2, np.finfo(float).eps
A, C0, ADD, PROP = 1.0e6, 3.0e3, 8.0e3, 0.25
BUDGET = 42
M_TRUTH, N_EXP, B = 800, 60, 60
GRID = np.geomspace(0.001, 1.0, 40)

LOG14 = np.array([0, .001, .003, .005, .007, .01, .03, .05, .07, .1, .3, .5, .7, 1.])
DESIGNS = {
    "log, 14 x 3\n(the real design)": (LOG14, 3),
    "linear, 14 x 3": (np.linspace(0, 1, 14), 3),
    "log, 7 x 6\n(fewer levels)": (LOG14[::2], 6),
    "top half only, 7 x 6": (np.array([0, .1, .3, .5, .7, .9, 1.]), 6),
}
DCOLORS = [PAL[0], PAL[1], PAL[2], PAL[3]]


def wts(x):
    return np.minimum(1.0 / (np.sqrt(x) + EPS), 1000)


def fit(x, y, w=None):
    res, _ = calc._fit_one_model(x, y, wts(x) if w is None else w, "bilinear")
    p = res.params
    return p["a"].value, p["b"].value, p["c"].value, p["c_high"].value


def pred(pr, grid):
    a, b, c, ch = pr
    return np.minimum(np.maximum(grid * a + b, c), ch)


def cv_of(ps, grid):
    P = np.vstack([pred(p, grid) for p in ps])
    return P.std(axis=0, ddof=1) / P.mean(axis=0)


def crossing(cv, grid, lo=None):
    """Interpolated last crossing of the CV threshold; explicit outcomes."""
    ok = np.isfinite(cv) & (grid > (lo if lo is not None else -np.inf))
    g, c = grid[ok], cv[ok]
    if len(g) < 2:
        return np.nan, "above_everywhere"
    above = np.where(c >= CV_THRESH)[0]
    if above.size == 0:
        return g[0], "no_crossing"
    i = above[-1]
    if i == len(g) - 1:
        return np.nan, "above_everywhere"
    x1, c1, x2, c2 = g[i], c[i], g[i + 1], c[i + 1]
    return (x2 if c2 == c1 else x1 + (CV_THRESH - c1) * (x2 - x1) / (c2 - c1)), "resolved"


def draw_design(levels, reps, seed):
    """One simulated experiment for a design: returns (x, y)."""
    x = np.repeat(levels, reps)
    mu = np.maximum(C0, A * x)
    sd = np.sqrt(ADD ** 2 + (PROP * mu) ** 2)
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    return x, np.maximum(mu + rng.normal(0, sd), 0.0)


# ------------------------------------------------------------------ simulated half
def simulate_designs():
    out = {}
    for k, (levels, reps) in DESIGNS.items():
        truth = []
        for i in range(M_TRUTH):
            x, y = draw_design(levels, reps, 700_000 + i)
            try:
                truth.append(fit(x, y))
            except Exception:
                pass
        true_loq, _ = crossing(cv_of(truth, GRID), GRID)
        est = []
        for e in range(N_EXP):
            x, y = draw_design(levels, reps, e)
            ps = []
            for i in range(B):
                rng = np.random.default_rng(np.random.SeedSequence([e, i, 7]))
                idx = rng.integers(0, len(x), len(x))
                if len(np.unique(y[idx])) < 2:
                    continue
                try:
                    ps.append(fit(x[idx], y[idx]))
                except Exception:
                    pass
            if len(ps) < 15:
                continue
            q, o = crossing(cv_of(ps, GRID), GRID)
            est.append(q if o == "resolved" else np.nan)
        out[k] = (true_loq, np.array(est, float))
        print(f"  design {k.splitlines()[0]:<24} true LOQ {true_loq:.5g}", flush=True)
    return out


# ------------------------------------------------------------------ real-data half
def real_data():
    df = calc.read_input(CURVE, CMAP)
    levels = np.sort(df["curvepoint"].unique())
    thin = np.array([0, .1, .3, .5, .7, 1.])          # drop the dense low end
    rows, grid_rows = [], []
    for pep in sorted(df["peptide"].dropna().unique()):
        sub_full = df[df["peptide"] == pep].sort_values(["curvepoint", "area"],
                                                        kind="mergesort")
        for label, keep in (("full", levels), ("thinned", thin)):
            s = sub_full[sub_full["curvepoint"].isin(keep)]
            x = np.asarray(s["curvepoint"], float)
            y = np.asarray(s["area"], float)
            if len(np.unique(x)) < 4:
                continue
            try:
                pr = fit(x, y)
            except Exception:
                continue
            mp = np.asarray([0.0, pr[2], pr[0], pr[1]])
            lod, _ = calc.calculate_lod(mp, s, 2.0, 2, 1, x, "auto")
            if not np.isfinite(lod) or lod <= 0:
                rows.append((pep, label, np.nan, np.nan))
                continue
            ps = []
            for i in range(60):
                rng = np.random.default_rng(np.random.SeedSequence([i, 3]))
                idx = rng.integers(0, len(x), len(x))
                if len(np.unique(y[idx])) < 2:
                    continue
                try:
                    ps.append(fit(x[idx], y[idx]))
                except Exception:
                    pass
            if len(ps) < 15:
                rows.append((pep, label, lod, np.nan))
                continue
            g = np.geomspace(lod, x.max(), 60)
            q, o = crossing(cv_of(ps, g), g, lo=lod)
            rows.append((pep, label, lod, q if o == "resolved" else np.nan))
            if label == "full":
                for gname, gg in (("uniform", np.linspace(lod, x.max(), 100)),
                                  ("log", np.geomspace(lod, x.max(), 100)),
                                  ("measured", np.unique(x[x > lod]))):
                    if len(gg) < 2:
                        continue
                    _, oo = crossing(cv_of(ps, gg), gg, lo=lod)
                    grid_rows.append((pep, gname, oo))
        print(f"  {pep}", flush=True)
    return pd.DataFrame(rows, columns=["peptide", "design", "LOD", "LOQ"]), \
        pd.DataFrame(grid_rows, columns=["peptide", "grid", "outcome"])


if os.path.exists(CACHE) and "--force" not in sys.argv:
    print("using cached results:", CACHE)
    z = np.load(CACHE, allow_pickle=True)
    sim = z["sim"].item()
    real = pd.DataFrame(z["real"], columns=["peptide", "design", "LOD", "LOQ"])
    real[["LOD", "LOQ"]] = real[["LOD", "LOQ"]].astype(float)
    grids = pd.DataFrame(z["grids"], columns=["peptide", "grid", "outcome"])
else:
    print("simulating designs...")
    sim = simulate_designs()
    print("re-fitting the sample dataset under thinned designs...")
    real, grids = real_data()
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, sim=np.array(sim, dtype=object),
                        real=real.to_numpy(dtype=object),
                        grids=grids.to_numpy(dtype=object))

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(17.5, 14.5))
gs = fig.add_gridspec(3, 4, hspace=0.52, wspace=0.30, top=0.945)

BLANK_X = 4e-4      # the 0-concentration point, drawn at the left edge of a log axis
Y_FLOOR = 1.0e3     # responses at or below zero are drawn here so they stay visible
PANEL = "ABCD"

for j, (k, (levels, reps)) in enumerate(DESIGNS.items()):
    ax = fig.add_subplot(gs[0, j])
    xd, yd = draw_design(levels, reps, 0)
    pr = fit(xd, yd)
    xs = np.geomspace(BLANK_X, 1.35, 400)
    ax.plot(xs, np.maximum(C0, A * xs), "--", color=PRIMARY, lw=2, label="true curve")
    ax.plot(xs, pred(pr, xs), "-", color=DCOLORS[j], lw=2.5, label="fitted curve")
    ax.plot(np.where(xd > 0, xd, BLANK_X), np.maximum(yd, Y_FLOOR), "o", ms=6.5,
            color=DCOLORS[j], alpha=0.75, markeredgecolor="white",
            markeredgewidth=0.6, label="simulated runs")
    t = sim[k][0]
    if np.isfinite(t):
        ax.axvline(t, color="0.35", lw=1.6, ls=":")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(BLANK_X * 0.6, 2.4)
    ax.set_ylim(Y_FLOOR * 0.6, 5e6)
    ax.set_xlabel("quantity")
    ax.set_title(PANEL[j] + "  " + " ".join(k.splitlines()), fontsize=14)
    ax.grid(True, alpha=0.3)
    if j == 0:
        ax.set_ylabel("response")
        ax.legend(loc="upper left", fontsize=10.5, frameon=False)
        ax.annotate("blank drawn at the left edge;", xy=(0.40, 0.12),
                    xycoords="axes fraction", fontsize=9.5, color="0.45")
        ax.annotate("zero responses at the axis floor", xy=(0.40, 0.04),
                    xycoords="axes fraction", fontsize=9.5, color="0.45")
    if j == 1:
        ax.annotate("dotted line = achievable LOQ,", xy=(0.34, 0.12),
                    xycoords="axes fraction", fontsize=10.5, color="0.3")
        ax.annotate("with no runs near it", xy=(0.34, 0.04),
                    xycoords="axes fraction", fontsize=10.5, color="0.3")

ax = fig.add_subplot(gs[1, :2])
for i, (k, (levels, reps)) in enumerate(DESIGNS.items()):
    xs = np.where(levels > 0, levels, BLANK_X)
    ax.plot(xs, [i] * len(levels), "o", ms=11, color=DCOLORS[i], alpha=0.9)
    ax.text(1.45, i, f"{reps} reps", va="center", fontsize=12, color=DCOLORS[i])
ax.set_xscale("log")
ax.set_yticks(range(len(DESIGNS)))
ax.set_yticklabels([" ".join(k.splitlines()) for k in DESIGNS], fontsize=12)
ax.set_xlabel("quantity (leftmost point is the blank)")
ax.set_xlim(2.5e-4, 3.2)
ax.set_ylim(-0.6, len(DESIGNS) - 0.4)
ax.invert_yaxis()
ax.set_title("E  Four ways to spend the same 42 injections", loc="left", fontsize=15)
ax.grid(True, alpha=0.3, axis="x")

ax = fig.add_subplot(gs[1, 2:])
keys = list(DESIGNS)
vals = [sim[k][1][np.isfinite(sim[k][1])] for k in keys]
bp = ax.boxplot(vals, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="white", lw=2))
for patch, c in zip(bp["boxes"], DCOLORS):
    patch.set_facecolor(c)
    patch.set_alpha(0.85)
    patch.set_edgecolor("none")
for i, k in enumerate(keys, start=1):
    t = sim[k][0]
    if np.isfinite(t):
        ax.plot([i - 0.36, i + 0.36], [t, t], color=PRIMARY, lw=2.5, ls="--", zorder=5)
ax.plot([], [], color=PRIMARY, lw=2.5, ls="--", label="achievable LOQ for that design")
ax.set_yscale("log")
ax.set_xticks(range(1, len(keys) + 1))
ax.set_xticklabels([" ".join(k.splitlines()).replace(" (the real design)", "")
                    .replace(" (fewer levels)", "") for k in keys],
                   fontsize=11, rotation=12, ha="right")
ax.set_ylabel("estimated LOQ (log)")
ax.set_title("F  The design sets both the achievable LOQ and the spread", loc="left",
             fontsize=15)
ax.legend(loc="upper left", fontsize=11.5, frameon=False)
ax.grid(True, alpha=0.3, axis="y")

ax = fig.add_subplot(gs[2, :2])
piv = real.pivot_table(index="peptide", columns="design", values="LOQ", aggfunc="first")
piv = piv.dropna()
for _, r in piv.iterrows():
    ax.plot([0, 1], [r["full"], r["thinned"]], "-o", ms=5, lw=1.2,
            color=PAL[0] if r["thinned"] > r["full"] else PAL[1], alpha=0.8)
worse = int((piv["thinned"] > piv["full"]).sum())
n_full = int(np.isfinite(real.loc[real["design"] == "full", "LOQ"]).sum())
n_thin = int(np.isfinite(real.loc[real["design"] == "thinned", "LOQ"]).sum())
lod_full = int(np.isfinite(real.loc[real["design"] == "full", "LOD"]).sum())
lod_thin = int(np.isfinite(real.loc[real["design"] == "thinned", "LOD"]).sum())
ax.set_xlim(-0.3, 1.3)
ax.set_xticks([0, 1])
ax.set_xticklabels(["full 14-level curve", "top half only"], fontsize=12)
ax.set_yscale("log")
ax.set_ylabel("LOQ (log)")
ax.set_title(f"G  Real peptides: thinning costs {n_full - n_thin} of {n_full} "
             f"their LOQ outright", loc="left", fontsize=15)
ax.grid(True, alpha=0.3, axis="y")
med = np.median(piv["thinned"] / piv["full"])
ax.annotate(f"peptides with a finite LOQ: {n_full} -> {n_thin}"
            f"   (finite LOD: {lod_full} -> {lod_thin})",
            xy=(0.5, 0.11), xycoords="axes fraction", ha="center", fontsize=12.5,
            color="0.3")
ax.annotate(f"of the {len(piv)} that survive, median {med:.1f}x higher",
            xy=(0.5, 0.04), xycoords="axes fraction", ha="center", fontsize=12.5,
            color="0.3")

ax = fig.add_subplot(gs[2, 2:])
order = ["uniform", "log", "measured"]
labels = ["uniform (the tool)", "log-spaced", "measured levels"]
res = [int(((grids["grid"] == g) & (grids["outcome"] == "resolved")).sum()) for g in order]
noc = [int(((grids["grid"] == g) & (grids["outcome"] == "no_crossing")).sum()) for g in order]
xs = np.arange(len(order))
ax.bar(xs, res, width=0.55, color=PAL[2], label="LOQ resolved")
ax.bar(xs, noc, width=0.55, bottom=res, color=PAL[3], label="reported as no crossing")
for i in range(len(order)):
    ax.text(i, res[i] / 2, str(res[i]), ha="center", va="center", color="white",
            fontsize=14)
    if noc[i]:
        ax.text(i, res[i] + noc[i] / 2, str(noc[i]), ha="center", va="center",
                color="#5A2A44", fontsize=14)
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=12)
ax.set_ylabel("peptides")
ax.set_ylim(0, max(np.array(res) + np.array(noc)) * 1.32)
ax.set_title("H  ...and where the CV is read decides how many resolve at all",
             loc="left", fontsize=15)
ax.legend(loc="upper right", fontsize=11.5, frameon=False)
ax.grid(True, alpha=0.3, axis="y")

fig.suptitle("Supplementary Figure - dilution spacing, not the bootstrap, sets what "
             "the curve can measure", fontsize=18, y=0.995)
out = repo_path(cfg["output"], "SUPP_curve_spacing")
fig.savefig(out + ".png", dpi=300, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png / .pdf")
for k in keys:
    v = sim[k][1][np.isfinite(sim[k][1])]
    print(f"  {k.splitlines()[0]:<26} true {sim[k][0]:.5g}  est median "
          f"{np.median(v):.5g}  n={len(v)}")
print(f"  real: finite LOQ {n_full} -> {n_thin} when thinned; of the {len(piv)} "
      f"surviving, {worse} worsen, median {med:.2f}x higher")
print(f"  grid outcomes resolved: {dict(zip(order, res))}")
