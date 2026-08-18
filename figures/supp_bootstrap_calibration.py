"""Supplement - why changing the bootstrap resampling scheme does not fix the LOQ.

Ground truth is known here: a segmented curve sampled on the real log-spaced design
(14 levels x 3 replicates) with combined additive + proportional error. The true
sampling variability of the fitted response comes from simulating many INDEPENDENT
experiments and taking the spread of their fits -- exactly the quantity a bootstrap
tries to estimate from a single experiment.

Four schemes are compared against it:
  case         rows resampled across the whole curve (what the tool does)
  stratified   rows resampled within each curve point
  wild         x fixed; each residual keeps its magnitude, sign randomised
  bayesian     x fixed; every row kept with a random Dirichlet weight

A scheme is calibrated when bootstrap CV / true CV is about 1. The current scheme
turns out to be the best calibrated of the four, the cleaner-looking alternatives
understate uncertainty, and the LOQ spread across experiments dwarfs the difference
between schemes -- the limit is the design, not the resampler.

Self-contained simulation; needs no raw MS data. Caches to output/ (gitignored);
pass --force to recompute.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.style import set_style, PAL, PRIMARY
from src.fom_io import load_config, repo_path, load_tool

set_style()
cfg = load_config()
calc = load_tool(cfg)

# ---------------------------------------------------------------- simulation setup
LEVELS = np.array([0, .001, .003, .005, .007, .01, .03, .05, .07, .1, .3, .5, .7, 1.])
REPS = 3
A, C0, ADD, PROP = 1.0e6, 3.0e3, 8.0e3, 0.25      # truth: slope, noise floor, error
CV_THRESH, EPS = 0.2, np.finfo(float).eps
M_TRUTH, N_EXP, B = 1500, 100, 100
SCHEMES = ("case", "stratified", "wild", "bayesian")
NICE = {"case": "case\n(current)", "stratified": "stratified",
        "wild": "wild", "bayesian": "Bayesian"}
COLOR = dict(zip(SCHEMES, [PAL[1], PAL[0], PAL[2], PAL[3]]))

X = np.repeat(LEVELS, REPS)
N = len(X)
BLOCKS = [np.arange(j * REPS, (j + 1) * REPS) for j in range(len(LEVELS))]
MU = np.maximum(C0, A * X)
SD = np.sqrt(ADD ** 2 + (PROP * MU) ** 2)
GRID = np.geomspace(0.001, 1.0, 40)
BAND = (GRID >= 0.003) & (GRID <= 0.3)
CACHE = repo_path(cfg["output"], "_cache_bootstrap_calibration.npz")


def wts(x):
    return np.minimum(1.0 / (np.sqrt(x) + EPS), 1000)


W = wts(X)


def fit(y, x=X, w=W):
    res, _ = calc._fit_one_model(x, y, w, "bilinear")
    p = res.params
    return p["a"].value, p["b"].value, p["c"].value, p["c_high"].value


def pred(pr, grid=GRID):
    a, b, c, ch = pr
    return np.minimum(np.maximum(grid * a + b, c), ch)


def cv_of(ps):
    P = np.vstack([pred(p) for p in ps])
    return P.std(axis=0, ddof=1) / P.mean(axis=0)


def loq_of(cv, grid=GRID):
    ok = np.isfinite(cv)
    g, c = grid[ok], cv[ok]
    above = np.where(c >= CV_THRESH)[0]
    if above.size == 0:
        return g[0], "no_crossing"
    i = above[-1]
    if i == len(g) - 1:
        return np.inf, "above_everywhere"
    x1, c1, x2, c2 = g[i], c[i], g[i + 1], c[i + 1]
    return (x2 if c2 == c1 else x1 + (CV_THRESH - c1) * (x2 - x1) / (c2 - c1)), "resolved"


def replicate(scheme, y, fitted, resid, rng):
    if scheme == "case":
        while True:
            idx = rng.integers(0, N, N)
            if len(np.unique(y[idx])) > 1:
                break
        return fit(y[idx], X[idx], wts(X[idx]))
    if scheme == "stratified":
        idx = np.concatenate([rng.choice(b, size=len(b), replace=True) for b in BLOCKS])
        return fit(y[idx], X[idx], wts(X[idx]))
    if scheme == "wild":
        return fit(fitted + rng.choice([-1.0, 1.0], size=N) * resid)
    w = rng.exponential(size=N)
    return fit(y, X, W * np.sqrt(w / w.mean()))


def simulate():
    truth = []
    for i in range(M_TRUTH):
        rng = np.random.default_rng(np.random.SeedSequence(900_000 + i))
        try:
            truth.append(fit(np.maximum(MU + rng.normal(0, SD), 0.0)))
        except Exception:
            pass
    true_cv = cv_of(truth)
    curves = {s: [] for s in SCHEMES}
    loqs = {s: [] for s in SCHEMES}
    for e in range(N_EXP):
        rng0 = np.random.default_rng(np.random.SeedSequence(e))
        y = np.maximum(MU + rng0.normal(0, SD), 0.0)
        try:
            base = fit(y)
        except Exception:
            continue
        fitted = pred(base, X)
        resid = y - fitted
        for s in SCHEMES:
            ps = []
            for i in range(B):
                rng = np.random.default_rng(
                    np.random.SeedSequence([e, i, abs(hash(s)) % 9973]))
                try:
                    ps.append(replicate(s, y, fitted, resid, rng))
                except Exception:
                    pass
            if len(ps) < 20:
                continue
            cv = cv_of(ps)
            curves[s].append(cv)
            q, o = loq_of(cv)
            loqs[s].append(q if o == "resolved" else np.nan)
        if (e + 1) % 25 == 0:
            print(f"  {e + 1}/{N_EXP} experiments", flush=True)
    out = {"true_cv": true_cv, "grid": GRID}
    for s in SCHEMES:
        out["cv_" + s] = np.vstack(curves[s])
        out["loq_" + s] = np.array(loqs[s], float)
    return out


if os.path.exists(CACHE) and "--force" not in sys.argv:
    print("using cached simulation:", CACHE)
    d = dict(np.load(CACHE))
else:
    print(f"simulating {N_EXP} experiments x {len(SCHEMES)} schemes x {B} replicates...")
    d = simulate()
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    np.savez_compressed(CACHE, **d)

TRUE_CV = d["true_cv"]
GRID = d["grid"]
TRUE_LOQ, _ = loq_of(TRUE_CV, GRID)
ratios = {s: np.ravel(d["cv_" + s][:, BAND] / TRUE_CV[BAND]) for s in SCHEMES}

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(17, 9.5))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.95], hspace=0.48, wspace=0.30)

for j, s in enumerate(SCHEMES):
    ax = fig.add_subplot(gs[0, j])
    cvs = d["cv_" + s]
    for row in cvs[:30]:
        ax.plot(GRID, row, color=COLOR[s], lw=0.9, alpha=0.35)
    ax.plot(GRID, np.median(cvs, axis=0), color=COLOR[s], lw=3, label="bootstrap")
    ax.plot(GRID, TRUE_CV, color=PRIMARY, lw=3, ls="--", label="true")
    ax.axhline(CV_THRESH, color="0.45", lw=1.5, ls=":")
    ax.set_xscale("log")
    ax.set_ylim(0, 0.85)
    ax.set_xlabel("quantity")
    if j == 0:
        ax.set_ylabel("CV of fitted response")
        ax.legend(loc="upper right", fontsize=12, frameon=False)
    ax.set_title(f"{'ABCD'[j]}  {NICE[s]}\nbootstrap / true = {np.median(ratios[s]):.2f}",
                 fontsize=14)
    ax.grid(True, alpha=0.3)

ax = fig.add_subplot(gs[1, :2])
bp = ax.boxplot([ratios[s] for s in SCHEMES], patch_artist=True, showfliers=False,
                widths=0.6, medianprops=dict(color="white", lw=2))
for patch, s in zip(bp["boxes"], SCHEMES):
    patch.set_facecolor(COLOR[s])
    patch.set_alpha(0.85)
    patch.set_edgecolor("none")
ax.axhline(1.0, color=PRIMARY, lw=2, ls="--")
ax.text(4.55, 1.04, "calibrated", va="bottom", color=PRIMARY, fontsize=12)
ax.set_xticks(range(1, 5))
ax.set_xticklabels([NICE[s].replace("\n", " ") for s in SCHEMES], fontsize=13)
ax.set_ylabel("bootstrap CV / true CV")
ax.set_ylim(0.35, 1.8)
ax.set_xlim(0.4, 5.4)
ax.set_title("E  Only the current scheme estimates the right uncertainty", loc="left",
             fontsize=15)
ax.grid(True, alpha=0.3, axis="y")
for i, s in enumerate(SCHEMES, start=1):
    ax.annotate(f"{np.median(ratios[s]):.2f}", xy=(i, 1.66), ha="center", fontsize=13,
                color=COLOR[s])

ax = fig.add_subplot(gs[1, 2:])
vals = [d["loq_" + s][np.isfinite(d["loq_" + s])] for s in SCHEMES]
bp = ax.boxplot(vals, patch_artist=True, showfliers=False, widths=0.6,
                medianprops=dict(color="white", lw=2))
for patch, s in zip(bp["boxes"], SCHEMES):
    patch.set_facecolor(COLOR[s])
    patch.set_alpha(0.85)
    patch.set_edgecolor("none")
ax.axhline(TRUE_LOQ, color=PRIMARY, lw=2, ls="--")
ax.text(4.55, TRUE_LOQ * 1.06, "true LOQ", va="bottom", color=PRIMARY, fontsize=12)
ax.set_xticks(range(1, 5))
ax.set_xticklabels([NICE[s].replace("\n", " ") for s in SCHEMES], fontsize=13)
ax.set_ylabel("estimated LOQ")
ax.set_xlim(0.4, 5.4)
ax.set_ylim(0, 0.056)
ax.set_title("F  ...but each scheme's spread dwarfs the gap between them", loc="left",
             fontsize=15)
ax.grid(True, alpha=0.3, axis="y")
iqrs = [np.percentile(v, 75) - np.percentile(v, 25) for v in vals]
ax.annotate(f"across replicate experiments the IQR is\n"
            f"{min(iqrs) / TRUE_LOQ:.0%}-{max(iqrs) / TRUE_LOQ:.0%} of the true LOQ itself",
            xy=(0.03, 0.90), xycoords="axes fraction", fontsize=12.5, color="0.3")

fig.suptitle("Supplementary Figure - the bootstrap scheme is not what limits the LOQ",
             fontsize=18, y=0.99)
out = repo_path(cfg["output"], "SUPP_bootstrap_calibration")
fig.savefig(out + ".png", dpi=300, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png / .pdf")
print(f"  true LOQ {TRUE_LOQ:.5g}")
for s in SCHEMES:
    v = d["loq_" + s][np.isfinite(d["loq_" + s])]
    print(f"  {s:<11} ratio {np.median(ratios[s]):.2f}  "
          f"LOQ median {np.median(v):.5g} ({np.median(v) / TRUE_LOQ - 1:+.0%})  "
          f"IQR {np.percentile(v, 75) - np.percentile(v, 25):.5g}")
