# Note: how precisely is the LOQ determined?

Working note for the methods/limitations discussion. Recorded 2026-08-11, revised
the same day — see [Correction](#correction-to-the-first-version-of-this-note) at
the end. Companion to [`fit_weighting_note.md`](fit_weighting_note.md).

> **Superseded in part (2026-08-12).** This note names the resampling scheme as the
> dominant cause of LOQ instability. Later work found that the *readout grid* matters
> more: the CV is read on a uniform `linspace` while the design is log-spaced, so the
> crossing often falls inside the interval the grid steps over. Correcting the spacing
> makes case and Bayesian resampling produce identical outcomes, which means the
> scheme differences measured here were substantially grid spacing in disguise. Read
> [`loq_grid_and_resampling_note.md`](loq_grid_and_resampling_note.md) alongside this
> one; the measurements below stand, the attribution does not.
>
> **Superseded again (2026-08-12).** This note treats stratified resampling as
> the fix. A calibration study against a known ground truth has since shown that
> stratified resampling *understates* uncertainty by about 20%, and that the
> scheme already in the tool is the best calibrated of four tested. Do not read
> this note as an argument for changing the resampler.

**Not measured on the paper's datasets.** Every number below comes from
`data/one_protein.csv`, the tool's 27-peptide sample dataset (yeast, EncyclopeDIA,
14 curve points × 3 replicates). It characterizes the *estimator*, not our results.
Re-run the recipes on the real datasets before quoting any figure in the manuscript.

## The finding

The reported LOQ is a single draw from a distribution that, for some peptides, is
wide. The tool is deterministic — `bootstrap_many` seeds replicate *i* with
`SeedSequence(i)`, so the same input always gives the same answer — but determinism
pins one draw without making it precise.

Re-running the same peptides against 8 independent families of bootstrap seeds, at
the default `--bootreps 100`:

| peptide | LOD | LOQ min | LOQ max | spread |
|---|---|---|---|---|
| GEGFMVVTATGDNTFVGR | 0.0195 | 0.0294 | 0.247 | **310%** |
| LMNGKPMK | 0.0700 | 0.126 | 0.230 | 59% |
| KVTAVVESPEGER | 0.00893 | 0.0189 | 0.0290 | 36% |
| VLEFHPFDPVSK | 0.00910 | 0.0191 | 0.0191 | 0% |
| KADTGIAVEGATDAAR | 0.00838 | 0.0284 | 0.0284 | 0% |
| TLANTAVVIR | 0.00373 | 0.0138 | 0.0138 | 0% |

The LOD and the fit coefficients are stable across the same perturbations (varying
by 1e-5 relative or less). It is specifically the bootstrap-derived LOQ that moves,
not the segmented fit.

## The dominant cause: the bootstrap discards whole curve levels

`_bootstrap_once` resamples **rows**, across the peptide's whole curve:

```python
resampled_df = df.sample(n=len(df), replace=True, random_state=rng)
```

This is case resampling, and nothing in it ties the draw to the curve's design. Each
row has a ~37% chance of going undrawn, and a concentration level survives only if at
least one of its three replicates is drawn. Measured over 5000 draws on the real 14×3
curve:

| | |
|---|---|
| replicates missing **at least one level entirely** | **50.2%** |
| mean levels missing per replicate | 0.63 |
| missing exactly two levels | 10.4% |
| the 0-concentration blank absent entirely | 5.1% |

So half of all bootstrap replicates fit a calibration curve with at least one
concentration simply not present. There is a second-order effect too:
`_initialize_params_auto` defines the noise region as the two lowest levels *present
in the resample*, so in 8.9% of replicates it treats genuine signal points as noise
and drags that replicate's noise plateau upward.

**Stratified resampling — drawing within each concentration level, preserving the
design — removes most of this.** Five seed families each, `--bootreps 100`:

| peptide | case resampling (current) | stratified |
|---|---|---|
| GEGFMVVTATGDNTFVGR | 250.8% | **0.0%** |
| LMNGKPMK | 60.9% | 42.5% |
| KVTAVVESPEGER | 37.1% | 37.1% |
| TLANTAVVIR | 0.0% | 0.0% |

The argument for stratified resampling is not just that it is more stable. The
concentration levels are a *designed* factor — chosen and pipetted, not sampled.
Case resampling treats x as random and lets the design dissolve; stratified
resampling respects it and estimates variation in the response at each fixed level,
which is what a calibration curve's CV is meant to describe.

## The residual: Monte Carlo noise

`KVTAVVESPEGER` is unchanged by stratification (37.1% either way) but collapses to
0% by 400 replicates. Nothing is wrong with the estimator there; 100 replicates is
simply too few to pin the CV curve, and a higher `--bootreps` fixes it outright.

`LMNGKPMK` is only partly helped (60.9% → 42.5%), so stratification is not a
complete answer for every peptide.

## Prevalence

Probing the same 27 peptides by permuting input row order instead of reseeding
(equivalent perturbation, done before the row order was made canonical): the LOQ
moved for 23 of 27 peptides. Binned by relative spread — 14 stable (<0.1%), 5 at
5–25%, 5 above 25%, 3 non-finite. So roughly 10 of the 24 peptides with a finite
LOQ carried more than 5% of arbitrary variation.

The two probes are not interchangeable, and disagree on individual peptides —
`VLEFHPFDPVSK` moved 46% under row permutation but 0% under reseeding. Treat the
prevalence figure as an order-of-magnitude indication, not a per-peptide result.

## Relationship to the code changes

This is *not* the row-order bug (matrix-matched_calcurves issue #16, fixed in
`ffa5118`). That bug made the LOQ depend on the order rows arrived in, which is now
canonical. Fixing it removed the visible symptom — the tool no longer gives
different answers for the same data — while leaving the underlying estimator
variance untouched. That is why this note exists: the instability is now invisible
in normal use.

**Stratified resampling is not implemented.** As of this revision the tool still
uses case resampling; the comparison above was run with a patched copy. Adopting it
would change LOQ values, so it is a decision to make *before* regenerating the
paper's figures of merit, not after.

**Submodule pin.** `tools/matrix-matched_calcurves` is currently at `ac0b951`, which
predates today's tool commits (`3d18164`, issue #15 DIA-NN densification; `ffa5118`,
issue #16 input normalization; `ffb1087`, piecewise + multiplier fixes). The
measurements above were taken against `ffa5118`. Bump the pin before reproducing.

We decided not to change the tool's *output* — reporting an interval or a
"poorly determined" flag would need a defensible definition and would widen the
output schema mid-project. Recording it here instead.

## Correction to the first version of this note

The first version attributed `GEGFMVVTATGDNTFVGR`'s spread to **non-identifiability**
— a CV curve approaching the 0.2 threshold too shallowly for the crossing to be
determined by the data, and therefore not fixable by resampling. That was wrong. The
evidence for it was that raising `--bootreps` to 1600 left the spread at ~40%, which
is true but does not imply what I took it to imply: more replicates of a *biased*
resampling scheme do not remove the bias. Stratified resampling takes the same
peptide to 0.0% spread.

The corrected reading is that the dominant cause is a flaw in the resampling scheme,
not a property of that peptide's curve. Do not describe any peptide here as
intrinsically non-identifiable on the strength of this note.

## Reproducing

Against a checkout of the tool at `ffa5118` or later, from its repo root. This
compares the two resampling schemes; drop the `stratified` branch to reproduce the
seed-family table alone.

```python
import importlib, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "bin")
calc = importlib.import_module("calculate-loq")

df = calc.read_input("data/one_protein.csv", "data/filename2samplegroup_map.csv")
sub = (df[df["peptide"] == "GEGFMVVTATGDNTFVGR"]
       .sort_values(calc.SORT_KEYS, kind="mergesort"))
x, y = np.asarray(sub["curvepoint"], float), np.asarray(sub["area"], float)

res, _ = calc.fit_by_lmfit_yang(x, y, "auto")
a, b = res.params["a"].value, res.params["b"].value
mp = np.asarray([0.0, res.params["c"].value, a, b])
lod, _ = calc.calculate_lod(mp, sub, 2.0, 2, 1, x, "auto")
grid = np.linspace(lod, max(x), 100)
model = "trilinear" if np.isfinite(res.params["c_high"].value) else "bilinear"
groups = [g for _, g in sub.groupby("curvepoint", sort=True)]

def cv_curve(family, stratified, reps=100):
    rows = []
    for i in range(reps):
        rng = np.random.default_rng(np.random.SeedSequence(family * 10_000 + i))
        while True:
            r = (pd.concat([g.sample(n=len(g), replace=True, random_state=rng)
                            for g in groups]) if stratified
                 else sub.sample(n=len(sub), replace=True, random_state=rng))
            if r["area"].nunique() > 1:
                break
        rf, _ = calc.fit_by_lmfit_yang(np.asarray(r["curvepoint"], float),
                                       np.asarray(r["area"], float), model)
        p = rf.params
        rows.append(np.minimum(np.maximum(grid * p["a"].value + p["b"].value,
                                          p["c"].value), p["c_high"].value))
    mat = np.vstack(rows)
    return mat.std(axis=0, ddof=1) / mat.mean(axis=0)

for stratified in (False, True):
    loqs = []
    for family in range(5):
        cv = cv_curve(family, stratified)
        good = (grid > lod) & (cv < 0.2)
        loqs.append(grid[good].min() if good.any() else np.inf)
    print("stratified" if stratified else "case", loqs)
```

To count how often a replicate loses a level, resample and compare
`set(r["curvepoint"].unique())` against the full set of levels.
