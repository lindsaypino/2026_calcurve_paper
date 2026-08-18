# Note: where the CV is read matters more than how it is resampled

Working note for the methods/limitations discussion. Recorded 2026-08-12. Supersedes
the "dominant cause" framing in [`loq_stability_note.md`](loq_stability_note.md);
companion to [`fit_weighting_note.md`](fit_weighting_note.md).

**Not measured on the paper's datasets.** Everything below comes from
`data/one_protein.csv`, the tool's 27-peptide sample dataset (yeast, EncyclopeDIA,
14 log-spaced levels × 3 replicates), against tool commit `54c3f36`. 26 peptides
have a finite LOD and are analysed; the other one has no LOD, so no CV grid.

**No code was changed.** These are evaluations run from patched copies. The tool
still uses case resampling on a uniform grid.

## The finding

`calculate_loq` reads the bootstrap CV on `np.linspace(LOD, max_x, 100)` — a
**uniform** grid — while the curve's concentrations are **log-spaced** (1, 0.7, 0.5,
0.3, 0.1 … 0.005, 0.003, 0.001, 0). The two disagree at exactly the end where the
LOQ lives, so the tool cannot see the bottom of its own curve.

How bad the mismatch is, per peptide:

| | |
|---|---|
| peptides whose **first** grid point above the LOD already skips ≥1 measured level | **16 of 26** |
| median measured levels skipped by that single first step | 1 (max 4) |
| fraction of the usable design skipped by one grid step | median 11%, max 33% |

`VVEILQNR` is the worst: its first eligible grid point is 0.0129, stepping over 4 of
the 12 levels that sit above its LOD. Any LOQ it reports is coarser than the
experiment that produced it.

## Consequence: "no crossing" is mostly a grid artifact

Reading the *same bootstrap replicates* at different x, with the crossing
interpolated and eligibility `x > LOD` (matching `calculate_loq`):

| where the CV is read | case: resolved / no crossing / never | bayesian |
|---|---|---|
| uniform, `linspace` (current) | 15 / 9 / 2 | 13 / 11 / 2 |
| log-spaced, `geomspace` | **19 / 5 / 2** | **19 / 5 / 2** |
| measured concentrations only | 18 / 6 / 2 | 17 / 7 / 2 |

Log spacing rescues 4 peptides under case resampling and 6 under Bayesian, and their
LOQs land roughly half the value the uniform grid was assigning them:

| peptide | LOD | uniform grid floor | log-grid LOQ |
|---|---|---|---|
| VVEILQNR | 0.0028 | 0.0129 | 0.0058 |
| TLANTAVVIR | 0.0037 | 0.0138 | 0.0084 |
| ADTGIAVEGATDAAR | 0.0063 | 0.0163 | 0.0102 |
| TVEEDHPIPEDVHENYENK | 0.0079 | 0.0179 | 0.0102 |

The crossing was inside the interval the uniform grid steps over. So most of what
earlier notes called "floor pinning" is not the data failing to determine an LOQ —
it is the readout never sampling the region where the answer is.

**And the resampling scheme largely stops mattering.** Under the uniform grid, case
and Bayesian disagree (15 vs 13 resolved); under a log grid they are *identical*
(19/5/2). The scheme differences measured earlier were substantially grid spacing in
disguise.

## Interpolating the crossing

Separately from spacing, the LOQ is currently the lowest **grid point** whose CV is
under threshold, rather than the interpolated crossing. Comparing a 100-point against
a 400-point grid on identical replicates, restricted to genuinely resolved crossings:

| definition | median shift | 90th pct | max |
|---|---|---|---|
| current (grid point) | 5.23% | 17.96% | 28.00% |
| interpolated | 0.91% | 2.87% | 5.11% |

Consistent across all four schemes tested (current 4.6–6.0%, interpolated
0.78–0.96%). The reported LOQ moves by a median 5%, up to 28%, purely from the choice
of `num=100` — an implementation detail with no scientific content.

Two things interpolation does **not** do:

- It does not recover any floor-pinned peptide. Outcomes relabel exactly 1:1 from
  "grid floor" to "no crossing in range". The gain is that the result stops
  masquerading as a measured LOQ.
- It does not improve seed-to-seed stability. Restricted to resolved crossings, the
  current definition shows a median spread of 0.0% (case) purely because grid
  snapping hides variation; interpolated, all four schemes converge on ~8% at 100
  replicates. That ~8% is the estimator's honest precision.

## Replicates buy Monte Carlo precision, not information

Worth stating explicitly in Methods, because the two are easy to conflate:

- **Monte Carlo error** — from using B finite replicates. Vanishes as B grows. This
  is the ~8% spread, and the seed-driven category flips below.
- **Statistical error** — from having 42 real measurements. Completely unaffected by
  B. Pushing B → ∞ converges to the exact bootstrap distribution, which is itself a
  plug-in estimate built from those 42 points.

The resample space is also small and discrete: stratified resampling of 3 replicates
admits only 10 distinct multisets per level. Only the Bayesian bootstrap's Dirichlet
weights vary continuously, which is a practical argument for it at this n. Note the
tool's CV is a CV of *model predictions*, so the fit pools all 42 points rather than
relying on 3 per level — but no resampling scheme creates information the design did
not collect.

## The outcome category is itself seed-dependent

Over 3 seed families, whether a peptide gets an LOQ at all can flip:

| scheme | resolved in all 3 | flips with the seed | never resolved |
|---|---|---|---|
| case | 13 | 3 | 10 |
| stratified | 11 | 2 | 13 |
| wild | 9 | 4 | 13 |
| bayesian | 12 | 3 | 11 |

Under case resampling the flippers are `KVTAVVESPEGER`, `TVEEDHPIPEDVHENYENK` and
`YGLNQMADEK`. **The single-family counts quoted above therefore carry ±2–4 peptides
of seed noise and should not be quoted as exact.** The ordering across conditions
survives; the individual integers do not. If these numbers head for the manuscript,
re-run over ~10 families and report the modal category per peptide with its
frequency.

## Which bootstrap is correct? Settled by simulation

Earlier versions of these notes assumed case resampling inflated the CV spuriously,
because it lets whole levels drop out of a replicate. A simulation with known ground
truth says the opposite. Simulating many INDEPENDENT experiments from a known curve
gives the real sampling variability of the fitted response; a bootstrap is calibrated
when its CV divided by that true CV is about 1.

| scheme | config A (prop 0.25) | config B (prop 0.15) | verdict |
|---|---|---|---|
| **case (what the tool does)** | **0.96** | **0.97** | well calibrated in both |
| stratified | 0.81 | 0.80 | understates by about 20% |
| wild | 0.86 | 0.90 | understates by 10-15% |
| Bayesian | 0.85 | 0.88 | understates by 12-15% |

The current scheme is the best calibrated of the four. Every "more principled"
alternative is anti-conservative: it would report the assay as more precise than it
is. The likely mechanism is small-sample bootstrap bias - with 3 replicates per
level, within-stratum resampling cannot express the true variance, and the wild
bootstrap inherits shrunken fitted residuals.

On the LOQ specifically no scheme wins consistently: case is +30% in config A and
+2% in B, Bayesian +7% and -4%, stratified -12% and -6%, wild -20% and -10%. The
bias is swamped by spread - the interquartile range across replicate experiments is
91-149% of the true LOQ itself. The schemes with the tightest IQR are exactly the
ones that understate uncertainty, so that tightness is anti-conservatism showing up
rather than precision.

**Conclusion: do not change the resampler.** The evaluation that led to a
`--bootstrap stratified` flag was reverted, and this is why.

Figure: `output/SUPP_bootstrap_calibration.png`, from
`figures/supp_bootstrap_calibration.py`.

## What the design does instead

Spending the same 42 injections on different dilution designs changes the
*achievable* LOQ, not merely its estimate:

| design | achievable LOQ |
|---|---|
| log-spaced, 14 levels x 3 | 0.0100 |
| log-spaced, 7 levels x 6 | 0.0112 |
| linear, 14 levels x 3 | 0.0536 |
| top half only, 7 levels x 6 | 0.0637 |

A linear series is roughly 5x worse than log spacing for the same instrument time,
because it spends most of its runs where the response is already strong and none
where the limit actually is. Trading levels for replicates (14x3 to 7x6) costs
very little.

On the sample dataset, thinning the design to drop the dense low end takes finite
LODs from 26 of 27 peptides down to 9, and finite LOQs from 21 down to 7.

Figure: `output/SUPP_curve_spacing.png`, from `figures/supp_curve_spacing.py`.

## The readout, calibrated against a known truth

The sections above compare readouts to each other. This one scores them against a
known answer, the same way the resampler was settled. Ground truth is simulated, so
the true CV curve and therefore the true LOQ are known exactly; each experiment is
bootstrapped once (case resampling) and the same replicates are read back several
ways. Two scenarios, because a readout can fail in opposite directions.

**Scenario A - a true crossing exists** (true LOQ 0.0237, median LOD 0.0070).
Counting only genuine crossings, and separating out the grid-floor values the
current rule emits:

| grid | floor values | genuine | bias of genuine |
|---|---|---|---|
| uniform, 100 pts | 41 of 116 | 75 | +93% |
| uniform, 400 pts | 37 | 79 | +69% |
| log, 100 pts | 29 | 87 | **+56%** |
| log, 400 pts | 28 | 88 | +54% |
| measured levels only | 33 | 82 | **+111%** |

Dependence on the arbitrary `num=` parameter, same replicates read at 100 vs 400
points:

| reader | uniform | log |
|---|---|---|
| current | median 15.9%, max 78.6% | median 2.6%, max 4.9% |
| interpolated | median 1.1%, max 29.0% | median 0.07%, max 0.58% |

**Scenario B - the truth has no crossing at all**, so any number reported is
fabricated:

| grid | current rule | interpolated + explicit |
|---|---|---|
| uniform, 100 pts | **100%** (median 0.0154) | 0% |
| uniform, 400 pts | **100%** (median 0.0078) | 0% |
| log, 100 pts | **100%** (median 0.0056) | 1% |
| measured levels | **100%** (median 0.0070) | 0% |

The current rule invents an LOQ in every experiment where none exists, and the
invented value is itself an artifact of grid density - it halves from 0.0154 to
0.0078 when the grid goes from 100 to 400 points. Interpolation with an explicit
outcome essentially never fabricates.

On matched pairs, where both readers report for the same experiment, interpolation
is also the less biased of the two: +67% against +93% on a uniform grid, +50%
against +53% on a log grid. The current rule only *appears* less biased in the raw
tables because its fabricated floor values sit low and drag the median down.

The cost of the explicit outcome is silence: in scenario A, where a crossing really
does exist, interpolation declines to report for 22% of experiments on a log grid
and 34% on a uniform one. Log spacing reduces the false negatives as well as the
bias, which is one more reason the two changes belong together.

**Best combination tested: log spacing plus interpolation plus an explicit
no-crossing outcome** - bias +50%, grid dependence 0.07%, fabrication 1%,
declines 22%.

Even so, the honest headline matches the resampling supplement: the best readout is
still biased about +50% high with an interquartile range near 100% of the true
value. Fixing the readout removes artifacts, it does not make the LOQ accurate. A
14 x 3 design does not pin the LOQ.

Figure: `output/SUPP_loq_readout.png`, from `figures/supp_loq_readout.py`.

### Correction: reading only at measured levels is not the best option

An earlier version of this note called reading the CV only at measured
concentrations the strictest and most principled choice, on the grounds that it
never claims quantitation at a concentration that was not run. That argument is
still true but the option performs worst: **+111% bias, the largest of the five
tested**, because with only 8-11 points above the LOD it is too coarse to locate
the crossing. Appealing in principle, the loser in practice. A log grid is the
better recommendation.

Caveats: one truth pair, `bilinear` model only, B = 100, 120 experiments per
scenario. Scenario A had to be rebuilt once - the first attempt put the true
crossing below the LOD, so there was nothing to recover.

## Ranking, if anything is ever changed

1. **Log grid spacing.** Largest effect: recovers 4-6 peptides, roughly halves their
   LOQ, and against a known truth cuts the bias of genuine crossings from +93% to
   +56% and the dependence on `num=` from a median 15.9% to 2.6%. Reading only at
   measured levels was the earlier recommendation and is worse than either
   (+111% bias) - see the correction above. A log grid suits a log-spaced design;
   a linear design would want spacing matched to it instead.
2. **Interpolated crossing**, with an explicit "no crossing in range" outcome. Removes
   the grid-resolution dependence, and stops the tool inventing an LOQ in 100% of
   the cases where none exists. Costs 22% false negatives on a log grid.
3. **Resampling scheme.** Leave it alone. Simulation against a known truth shows
   the current scheme is the best calibrated of the four tested; see above.
4. **Replicate count.** Only addresses the ~8% Monte Carlo spread, and bottoms out
   against the design.

## Reproducing

Scripts live in the session scratchpad, not in a repo: `crossing_experiment.py`
(4 schemes × 3 LOQ definitions × 2 grid sizes) and `grid_spacing.py` (uniform vs log
vs measured-level readout). Both import the tool from `bin/`, store each replicate's
fit parameters, and evaluate CV curves afterwards, so one set of fits answers several
questions. Re-create them from this note if needed; the key detail is that
eligibility must be `x > LOD`, matching `calculate_loq` — including `x == LOD` gives
every grid the same high-CV anchor point and hides the spacing effect entirely.
