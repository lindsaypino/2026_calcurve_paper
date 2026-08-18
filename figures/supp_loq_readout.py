"""Supplement - how the LOQ is read off the CV curve: grid spacing and the crossing rule.

Calibrated the same way as the resampling supplement: ground truth is known, so the
true CV curve and the true LOQ are known, and every readout can be scored against
them. Each simulated experiment is bootstrapped once (case resampling, settled in the
companion supplement) and the SAME replicates are read back several ways.

  grids        uniform linspace (what the tool does), log geomspace, measured levels
               only; 100 and 400 points where meaningful
  readers      current = lowest GRID POINT above the LOD whose CV is under threshold
               interp  = interpolated last crossing, with an explicit
                         "no crossing in range" outcome

Two truth scenarios, because a readout can fail in opposite directions:

  A  a crossing exists  -> can it be found, and where?   (bias, false negatives)
  B  no crossing exists -> is one invented anyway?       (fabrication)

Scenario B is the decisive one for the current rule: it always emits a number, so it
cannot express "the CV never got that bad".

Self-contained simulation; no raw MS data. Caches to output/ (gitignored); pass
--force to recompute.
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

LEVELS = np.array([0, .001, .003, .005, .007, .01, .03, .05, .07, .1, .3, .5, .7, 1.])
REPS, EPS, CV_THRESH = 3, np.finfo(float).eps, 0.2
A, C0 = 1.0e6, 3.0e3
SCENARIOS = {"A_crossing": (1.5e3, 0.50), "B_no_crossing": (1.2e3, 0.04)}
M_TRUTH, N_EXP, B = 1200, 120, 100

X = np.repeat(LEVELS, REPS)
N = len(X)
DENSE = np.geomspace(1e-4, 1.0, 4000)
CACHE = repo_path(cfg["output"], "_cache_loq_readout.npz")

GRIDS = ["uniform_100", "uniform_400", "log_100", "log_400", "measured"]
NICE_G = {"uniform_100": "uniform\n100 pts", "uniform_400": "uniform\n400 pts",
          "log_100": "log\n100 pts", "log_400": "log\n400 pts",
          "measured": "measured\nlevels"}
C_CUR, C_INT = PAL[1], PAL[0]


def wts(x):
    return np.minimum(1.0 / (np.sqrt(x) + EPS), 1000)


W = wts(X)


def fit(y, x=X, w=W):
    res, _ = calc._fit_one_model(x, y, w, "bilinear")
    p = res.params
    return p["a"].value, p["b"].value, p["c"].value, p["c_high"].value


def pred(pr, grid):
    a, b, c, ch = pr
    return np.minimum(np.maximum(grid * a + b, c), ch)


def cv_of(ps, grid):
    P = np.vstack([pred(p, grid) for p in ps])
    return P.std(axis=0, ddof=1) / P.mean(axis=0)


def read_current(cv, grid, lod):
    ok = (grid > lod) & np.isfinite(cv)
    g, c = grid[ok], cv[ok]
    if len(g) == 0:
        return np.nan, "above_everywhere"
    good = c < CV_THRESH
    if not good.any():
        return np.nan, "above_everywhere"
    q = g[good].min()
    return q, ("at_grid_floor" if np.isclose(q, g[0], rtol=1e-12) else "resolved")


def read_interp(cv, grid, lod):
    ok = (grid > lod) & np.isfinite(cv)
    g, c = grid[ok], cv[ok]
    if len(g) < 2:
        return np.nan, "above_everywhere"
    above = np.where(c >= CV_THRESH)[0]
    if above.size == 0:
        return np.nan, "no_crossing"
    i = above[-1]
    if i == len(g) - 1:
        return np.nan, "above_everywhere"
    x1, c1, x2, c2 = g[i], c[i], g[i + 1], c[i + 1]
    q = x2 if c2 == c1 else x1 + (CV_THRESH - c1) * (x2 - x1) / (c2 - c1)
    return q, "resolved"


READERS = {"current": read_current, "interp": read_interp}


def grids_for(lod):
    return {"uniform_100": np.linspace(lod, 1.0, 100),
            "uniform_400": np.linspace(lod, 1.0, 400),
            "log_100": np.geomspace(lod, 1.0, 100),
            "log_400": np.geomspace(lod, 1.0, 400),
            "measured": np.unique(X[X > lod])}


def simulate():
    rows, extra = [], {}
    for scen, (add, prop) in SCENARIOS.items():
        mu = np.maximum(C0, A * X)
        sd = np.sqrt(add ** 2 + (prop * mu) ** 2)
        truth, lods = [], []
        for i in range(M_TRUTH):
            rng = np.random.default_rng(np.random.SeedSequence(500_000 + i))
            y = np.maximum(mu + rng.normal(0, sd), 0.0)
            try:
                pr = fit(y)
            except Exception:
                continue
            truth.append(pr)
            mp = np.asarray([0.0, pr[2], pr[0], pr[1]])
            l, _ = calc.calculate_lod(mp, pd.DataFrame({"curvepoint": X, "area": y}),
                                      2.0, 2, 1, X, "auto")
            if np.isfinite(l):
                lods.append(l)
        true_cv = cv_of(truth, DENSE)
        med_lod = float(np.median(lods)) if lods else 0.0
        tq, _ = read_interp(true_cv, DENSE, med_lod)
        if scen == "A_crossing":
            extra["true_cv"] = true_cv
            extra["true_loq"] = np.array([tq])
            extra["med_lod"] = np.array([med_lod])
        print(f"  {scen}: true LOQ {tq:.5g}, median LOD {med_lod:.5g}", flush=True)

        for e in range(N_EXP):
            rng0 = np.random.default_rng(np.random.SeedSequence(e))
            y = np.maximum(mu + rng0.normal(0, sd), 0.0)
            try:
                base = fit(y)
            except Exception:
                continue
            mp = np.asarray([0.0, base[2], base[0], base[1]])
            lod, _ = calc.calculate_lod(mp, pd.DataFrame({"curvepoint": X, "area": y}),
                                        2.0, 2, 1, X, "auto")
            if not np.isfinite(lod) or lod <= 0:
                continue
            ps = []
            for i in range(B):
                rng = np.random.default_rng(np.random.SeedSequence([e, i, 11]))
                idx = rng.integers(0, N, N)
                if len(np.unique(y[idx])) < 2:
                    continue
                try:
                    ps.append(fit(y[idx], X[idx], wts(X[idx])))
                except Exception:
                    pass
            if len(ps) < 20:
                continue
            gg = grids_for(lod)
            if scen == "A_crossing" and "ex_lod" not in extra:
                extra["ex_lod"] = np.array([lod])
                extra["ex_gu"] = gg["uniform_100"]
                extra["ex_cu"] = cv_of(ps, gg["uniform_100"])
                extra["ex_gl"] = gg["log_100"]
                extra["ex_cl"] = cv_of(ps, gg["log_100"])
            for gname, grid in gg.items():
                if len(grid) < 2:
                    continue
                cv = cv_of(ps, grid)
                for rname, reader in READERS.items():
                    q, o = reader(cv, grid, lod)
                    rows.append((scen, e, gname, rname, q, o, tq, lod))
        print(f"  {scen} done", flush=True)
    return pd.DataFrame(rows, columns=["scenario", "exp", "grid", "reader", "loq",
                                       "outcome", "true_loq", "lod"]), extra


if os.path.exists(CACHE) and "--force" not in sys.argv:
    print("using cached simulation:", CACHE)
    z = np.load(CACHE, allow_pickle=True)
    r = pd.DataFrame(z["rows"], columns=["scenario", "exp", "grid", "reader", "loq",
                                         "outcome", "true_loq", "lod"])
    for col in ("loq", "true_loq", "lod"):
        r[col] = r[col].astype(float)
    ex = {k: z[k] for k in z.files if k != "rows"}
else:
    print("simulating two truth scenarios...")
    r, ex = simulate()
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, rows=r.to_numpy(dtype=object), **ex)

TRUE_LOQ = float(ex["true_loq"][0])
MED_LOD = float(ex["med_lod"][0])
a = r[r.scenario == "A_crossing"]
b = r[r.scenario == "B_no_crossing"]

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(17.5, 9.6))
gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.28, top=0.90)

# --- A: where each grid actually samples the CV curve
ax = fig.add_subplot(gs[0, :2])
ax.plot(DENSE, ex["true_cv"], color=PRIMARY, lw=3, ls="--", label="true CV")
lod = float(ex["ex_lod"][0])
gu, cu, gl, cl = ex["ex_gu"], ex["ex_cu"], ex["ex_gl"], ex["ex_cl"]
ax.plot(gl, cl, "o", ms=4.5, color=C_INT, alpha=0.7, label="read on a log grid")
ax.plot(gu, cu, "s", ms=5.5, color=C_CUR, alpha=0.95, markeredgecolor="white",
        markeredgewidth=0.5, label="read on a uniform grid")
ax.axhline(CV_THRESH, color="0.45", lw=1.5, ls=":")
ax.axvline(TRUE_LOQ, color="0.25", lw=2)
ax.annotate("true LOQ", xy=(TRUE_LOQ * 1.1, 0.62), fontsize=12, color="0.25")
ax.axvspan(lod, gu[1], color=C_CUR, alpha=0.13)
ax.annotate("the uniform grid has no point", xy=(lod * 1.03, 0.135),
            fontsize=11, color=C_CUR)
ax.annotate("anywhere in this shaded band", xy=(lod * 1.03, 0.075),
            fontsize=11, color=C_CUR)
ax.set_xscale("log")
ax.set_xlim(lod * 0.75, 1.05)
ax.set_ylim(0, 0.75)
ax.set_xlabel("quantity")
ax.set_ylabel("bootstrap CV")
ax.set_title("A  One experiment, one set of replicates, read two ways", loc="left",
             fontsize=15)
ax.legend(loc="upper right", fontsize=11, frameon=False)
ax.grid(True, alpha=0.3)

# --- B: dependence on the arbitrary grid density
ax = fig.add_subplot(gs[0, 2])
pos, labels, colors, data = [], [], [], []
for k, (rd, fam) in enumerate([("current", "uniform"), ("current", "log"),
                               ("interp", "uniform"), ("interp", "log")]):
    p1 = a[(a.grid == fam + "_100") & (a.reader == rd)].set_index("exp")["loq"]
    p4 = a[(a.grid == fam + "_400") & (a.reader == rd)].set_index("exp")["loq"]
    j = pd.concat([p1, p4], axis=1, keys=["g1", "g4"]).dropna()
    data.append((np.abs(j.g4 - j.g1) / j.g1 * 100).to_numpy())
    pos.append(k + 1)
    labels.append(f"{rd}\n{fam}")
    colors.append(C_CUR if rd == "current" else C_INT)
bp = ax.boxplot(data, positions=pos, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="white", lw=2))
for patch, c in zip(bp["boxes"], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.85)
    patch.set_edgecolor("none")
ax.set_yscale("log")
ax.set_xticks(pos)
ax.set_xticklabels(labels, fontsize=11)
ax.set_ylabel("LOQ shift, 100 vs 400 points (%)")
ax.set_title("B  Dependence on grid density", loc="left", fontsize=15)
ax.grid(True, alpha=0.3, axis="y")

# --- C: accuracy of genuine crossings
ax = fig.add_subplot(gs[1, 0])
meds, floors = [], []
for g in GRIDS:
    s = a[(a.grid == g) & (a.reader == "current")]
    gen = s.loc[s.outcome == "resolved", "loq"].to_numpy(float)
    meds.append(np.median(gen) / TRUE_LOQ - 1 if len(gen) else np.nan)
    floors.append(int((s.outcome == "at_grid_floor").sum()))
xs = np.arange(len(GRIDS))
bars = ax.bar(xs, [m * 100 for m in meds], width=0.6,
              color=[C_INT if "log" in g else (PAL[2] if g == "measured" else C_CUR)
                     for g in GRIDS])
for i, (m, f) in enumerate(zip(meds, floors)):
    ax.text(i, m * 100 + 4, f"{f} floor", ha="center", fontsize=10, color="0.35")
ax.axhline(0, color=PRIMARY, lw=2, ls="--")
ax.set_xticks(xs)
ax.set_xticklabels([NICE_G[g] for g in GRIDS], fontsize=10)
ax.set_ylabel("bias of genuine crossings (%)")
ax.set_ylim(0, max(meds) * 100 * 1.30)
ax.set_title("C  Even read well, the LOQ is biased high", loc="left", fontsize=15)
ax.grid(True, alpha=0.3, axis="y")

# --- D: fabrication when no crossing exists
ax = fig.add_subplot(gs[1, 1])
w = 0.36
xs = np.arange(len(GRIDS))
cur_rate = [float((np.isfinite(b[(b.grid == g) & (b.reader == "current")].loq)).mean())
            for g in GRIDS]
int_rate = [float((np.isfinite(b[(b.grid == g) & (b.reader == "interp")].loq)).mean())
            for g in GRIDS]
ax.bar(xs - w / 2, [v * 100 for v in cur_rate], width=w, color=C_CUR,
       label="current rule")
ax.bar(xs + w / 2, [v * 100 for v in int_rate], width=w, color=C_INT,
       label="interpolated + explicit")
ax.set_xticks(xs)
ax.set_xticklabels([NICE_G[g] for g in GRIDS], fontsize=10)
ax.set_ylabel("experiments given an LOQ (%)")
ax.set_ylim(0, 118)
ax.set_title("D  With no true LOQ, one rule still invents one", loc="left",
             fontsize=15)
ax.legend(loc="center right", fontsize=11, frameon=False)
ax.grid(True, alpha=0.3, axis="y")

# --- E: the cost of declining
ax = fig.add_subplot(gs[1, 2])
dec = [float((a[(a.grid == g) & (a.reader == "interp")].outcome == "no_crossing").mean())
       for g in GRIDS]
ax.bar(np.arange(len(GRIDS)), [v * 100 for v in dec], width=0.6,
       color=[C_INT if "log" in g else (PAL[2] if g == "measured" else C_CUR)
              for g in GRIDS])
ax.set_xticks(np.arange(len(GRIDS)))
ax.set_xticklabels([NICE_G[g] for g in GRIDS], fontsize=10)
ax.set_ylabel("declined despite a crossing (%)")
ax.set_ylim(0, max(dec) * 100 * 1.45)
ax.set_title("E  ...at the price of silence", loc="left", fontsize=15)
ax.grid(True, alpha=0.3, axis="y")
ax.annotate("log spacing lowers this too", xy=(0.5, 0.90), xycoords="axes fraction",
            ha="center", fontsize=11, color="0.35")

fig.suptitle("Supplementary Figure - the LOQ readout: grid spacing and the crossing rule",
             fontsize=18, y=0.975)
out = repo_path(cfg["output"], "SUPP_loq_readout")
fig.savefig(out + ".png", dpi=300, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png / .pdf")
print(f"  true LOQ {TRUE_LOQ:.5g}, median LOD {MED_LOD:.5g}")
for g, m, f, cr, ir, d in zip(GRIDS, meds, floors, cur_rate, int_rate, dec):
    print(f"  {g:<12} bias {m:+.0%}  floors {f:>3}  fabricate cur {cr:.0%} "
          f"int {ir:.0%}  declined {d:.0%}")
