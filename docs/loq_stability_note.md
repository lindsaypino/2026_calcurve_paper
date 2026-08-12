# Note: how precisely is the LOQ determined?

Working note for the methods/limitations discussion. Recorded 2026-08-11.

**Not measured on the paper's datasets.** Every number below comes from
`data/one_protein.csv`, the tool's 27-peptide sample dataset (yeast, EncyclopeDIA,
14 curve points × 3 replicates). It characterizes the *estimator*, not our results.
Re-run the recipe on the real datasets before quoting any figure in the manuscript.

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

## Two distinct causes, which matter differently

**Monte Carlo noise — curable.** `KVTAVVESPEGER` spreads 36% at 100 replicates and
collapses to 0% by 400. Nothing is wrong with the estimator here; 100 replicates is
simply too few to pin the CV curve. A higher `--bootreps` fixes it outright.

**Non-identifiability — not curable by resampling.** `GEGFMVVTATGDNTFVGR` spreads
85% at 100 replicates and is still at 40% with 1600. Its bootstrap CV curve
approaches the 0.2 threshold so shallowly that the crossing point genuinely is not
determined by the data. The LOQ for this peptide ranges over nearly an order of
magnitude (0.029 to 0.247) depending only on which resamples were drawn. More
replicates buy nothing, because the limit being estimated is not sharp.

The distinction is worth making explicitly: the first is a tuning parameter, the
second is a property of the peptide's curve. Both currently report a bare number
with no indication of which situation you are in.

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

**Submodule pin.** `tools/matrix-matched_calcurves` is currently at `ac0b951`, which
predates both of today's tool commits (`3d18164`, issue #15 DIA-NN densification;
`ffa5118`, issue #16 input normalization). The measurements above were taken against
`ffa5118`. Bump the pin before reproducing.

We decided not to change the tool's output — reporting an interval or a
"poorly determined" flag would need a defensible definition and would widen the
output schema mid-project. Recording it here instead.

## Reproducing

Against a checkout of the tool at `ffa5118` or later, from its repo root:

```python
import importlib, sys
import numpy as np
sys.path.insert(0, "bin")
calc = importlib.import_module("calculate-loq")

df = calc.read_input("data/one_protein.csv", "data/filename2samplegroup_map.csv")
sub = df[df["peptide"] == "GEGFMVVTATGDNTFVGR"].sort_values(calc.SORT_KEYS, kind="mergesort")
x = np.asarray(sub["curvepoint"], float)
y = np.asarray(sub["area"], float)

res, _ = calc.fit_by_lmfit_yang(x, y, "auto")
a, b = res.params["a"].value, res.params["b"].value
mp = np.asarray([0.0, res.params["c"].value, a, b])
lod, _ = calc.calculate_lod(mp, sub, 2.0, 2, 1, x, "auto")
grid = np.linspace(lod, max(x), 100)
model = "trilinear" if np.isfinite(res.params["c_high"].value) else "bilinear"

for family in range(8):                     # shift the replicate seeds
    mat = np.vstack([calc._bootstrap_once(sub, grid, family * 10_000 + i, model)
                     for i in range(100)])
    cv = mat.std(axis=0, ddof=1) / mat.mean(axis=0)
    good = (grid > lod) & (cv < 0.2)
    print(grid[good].min() if good.any() else np.inf)
```

Vary the `100` to reproduce the convergence check.
