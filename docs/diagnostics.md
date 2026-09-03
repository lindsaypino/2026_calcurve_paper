# Diagnostics and open questions

Internal notes from investigating the figures-of-merit machinery — findings that
aren't a paper panel but that we may want to revisit, either across the other
figures or in the pinned tool source. Line numbers refer to the pinned submodule
commit `ffb1087`.

## Settled: the 240610 SIS dataset has no ULOQ in it

The candidate "Bruker dataset with a ULOQ"
(`240610_SIS_peptide_response_Ultra_BC_refined`, a Skyline AQUA/SIS peptide-response
curve on the timsTOF Ultra) **does not contain saturation.** Ran 2026-08-03.

Design: 5 curve points, 3-fold serial dilution (0.00024692 → 0.02), triplicate = 15
runs; 26 SIS peptides in groups A–E spiked at per-peptide `concentration_multiplier`
1/10/100/1000/10000; plus 11 iRT peptides. Analyte is the **heavy** channel (light is
100% `#N/A`). All `sample_type="standard"`, no blank point. Reshaped by
[`src/prep_skyline_sis.py`](../src/prep_skyline_sis.py); see its docstring for the four
format quirks (positional peptide identity, exact column-name sniffing, the multiplier
file, concentrations living in the `.sky`).

Results, `--model auto --min_noise_points 0 --min_saturation_points 2`:

| | finite | of |
|---|---|---|
| LOD | 19 | 26 |
| LOQ | 18 | 26 |
| ULOQ | **0** | 26 |

`--min_noise_points 0` is required: the high-multiplier groups sit so far above the
noise floor that the fit resolves no noise plateau, and at the mnp=2 default only 6/26
get a finite LOD. Group A (×1) is pure noise and correctly returns `inf` — the tool
independently reproduced the analyst's own `Peptide Note` QC calls, including that
`EASGLSADSLAR` is the only usable group B peptide.

**Why no ULOQ:** the clause `n_saturated >= min_saturation_points` fails 26/26; `n_sat`
is 0 or 1, never 2. With 3-fold spacing two points share a plateau only if the top step
ratio ≈ 1.0; measured top steps are 2.2–2.8 against a nominal 3.0. The curves are
compressing slightly but still climbing at the highest level spiked.
`LGQHLATEPLGTNSWER` (1.36) is the only one that visibly bends.

Relaxing to `--min_saturation_points 1` produces 7 ULOQs, but they are artifacts, and
the tool's own plots show why (see `output/diag_sis_uloq_toolfits.png`): the fitted
ceiling lands **inside the top point's replicate scatter** in all 7, and the ULOQ falls
below the highest measured concentration in 7/7 and below the *second*-highest in 2/7.
By the CV metric the tool itself thresholds (CV of bootstrap-resampled means) the top
point passes in 7/7. Conclusion: the `min_saturation_points=2` guard is doing its job;
a real ULOQ needs new acquisition above the current top point, not a settings change.

## Reusable: the ULOQ clause audit

`calculate-loq.py` reports no reason when a peptide gets no ULOQ — acceptance is a
seven-clause conjunction (`fit_by_lmfit_yang`, the `supported = (...)` block) and any
one clause silently kills it. [`src/uloq_audit.py`](../src/uloq_audit.py) re-fits with
the tool's own `_fit_one_model` and reports all seven, plus which clause is the *sole*
blocker:

```bash
# melted.csv needs columns (peptide, curvepoint, area) — what the tool sees
# AFTER any --multiplier_file is applied
python -m src.uloq_audit <melted.csv> [min_saturation_points]
```

Validated against the SIS data: at msp=2 it reports `enough_sat_points` as the sole
blocker for exactly 7 peptides, and at msp=1 it predicts the identical 7-peptide
accepted set the tool actually produced.

## Open — worth a second look across the other figures

- [ ] **Audit `n_sat` behind the paper's existing ULOQ calls.** Fig 2's triptychs and
      Fig 1B's distributions plot ULOQ from `main/*.csv`, and
      `fig_uloq_examples_linear.py` picks example ULOQ peptides. None of those have been
      checked for how many distinct plateau points support each ULOQ. The risk is much
      lower than for the SIS curve — the Bruker A–N curves have 14 points, Exploris 15,
      IL15 9, so `n_sat` has real headroom (2–7) rather than the single viable value a
      5-point curve allows — but it is *untested*, and an example peptide resting on a
      2-point plateau would be a weak choice to feature.
      **Blocked on** populating `data/figuresofmerit/` (gitignored), or on raw data.
- [ ] **The manifest's claim that "ULOQ counts are `min_noise_points`-independent" is
      right but incomplete** — they are emphatically *not*
      `min_saturation_points`-independent, and that sensitivity has never been swept for
      the paper's datasets. Consider a supplementary sweep, or at least a Methods
      sentence.
- [ ] Decide whether the SIS work earns a supplement. The defensible framings are "the
      Ultra's linear range extends past the top of a 5-point SIS curve" and/or a clean
      demonstration of the mnp=0 rescue (6/26 → 19/26 finite LODs).
- [ ] **Track down `250528_Ultra2023_HumanYeast_100ng_Calcurves`.** Named 2026-08-20 as a
      possible replacement ULOQ dataset, but not found on this machine: not under
      `2023_bruker_ultra/data/` (29 folders, no `250528_*`), not anywhere in
      `G:/My Drive/Lindsay Pino/` to depth 6 for `*250528*` (only two unrelated
      `code-benchmark_*_20250528.ipynb` notebooks), not in `G:/My Drive/` to depth 4 for
      `*Ultra2023*`, not in `2025_calcurves_asms/data/`, not in `Downloads/`, and not
      mentioned in this repo. Nearest Drive names are `250515_CalCurve_HumanYeast_Timbaux`
      and `250526_CalCurve_HumanYeast_Desnaux_{60,100}SPD`. Two things to settle: where it
      lives (instrument PC, Bruker share, S3, PRIDE, or a Drive folder outside
      `Lindsay Pino/`), and whether it is meant to replace the SIS curve. The name reads as
      a matrix-matched human/yeast curve at 100 ng — the same design as the Timbaux/Desnaux
      curves in Fig 2, not an SIS peptide-response curve. That would be the better ULOQ
      candidate, because a matrix-matched curve can carry its top points into saturation,
      which is what the 240610 SIS curve lacked.

## Open — worth a second look in the tool source (submodule)

- [ ] **`std_sat` pushes the ULOQ the wrong way.** `calculate_uloq` (line 560) computes
      `ULOQ = (c_high - std_mult*std_sat - intercept)/slope` at line 577. When the
      plateau rests on one concentration, `std_sat` is just that concentration's
      replicate spread, so a *noisier* top point yields a *lower* ULOQ — i.e. a
      narrower claimed quantifiable range. That is backwards: added noise should widen
      uncertainty, not tighten the range. On the SIS data it drags ULOQ well below the
      onset (45 vs onset ~90 for `LGQHLATEPLGTNSWER`). Check whether it materially
      moves ULOQ on the real datasets.
- [ ] **Short curves have exactly one viable `n_sat`.** Line 414 refuses a trilinear fit
      below 5 distinct curve points, and line 438 caps the plateau via
      `2*n_sat <= n_distinct`. At `n_distinct == 5` that leaves `n_sat == 2` as the only
      value satisfying both it and the msp=2 minimum. Worth deciding whether the
      minority rule should scale differently for short curves, and worth stating the
      ≥5-point gate in Methods — it silently makes a ULOQ impossible for ≤4-point
      curves, and isn't in the README settings table.
- [ ] **Document which CV the LOQ thresholds.** Lines 723–730 plot and threshold the CV
      of bootstrap-resampled *means*, not the raw replicate CV — with n=3 that is a ~1.7×
      difference and it directly sets the LOQ. Legitimate choice, but Methods should say
      which. Related asymmetry worth a sentence: LOQ is CV-gated, ULOQ is purely
      geometric with no CV check at all.
- [ ] **Cosmetic: ULOQ is drawn but never labeled.** `build_plots` draws the ULOQ
      vertical in the top panel (line 700) but builds the legend from the *bottom*
      subplot (line 759), which re-draws only LOD and LOQ (lines 734–740). So generated
      PNGs show an unlabeled orange line and no ULOQ legend entry — worth fixing before
      any of these plots go into a supplement.

## Artifacts

Diagnostic figures land in `output/` (gitignored):

| file | script |
|---|---|
| `diag_sis_curves.png` | [`figures/diag_sis_curves.py`](../figures/diag_sis_curves.py) — all 26 curves log-log, one column per spike group, slope-1 reference + top-step ratio |
| `diag_sis_uloq_toolfits.png` | the tool's own `--plot y` output for the 7 msp=1 ULOQ peptides, tiled |

The run's inputs and figures of merit are committed under
[`output/sis_run/`](../output/sis_run/): the three tool-ready inputs, the melted frame,
and the three FOM runs (`fom_mnp0_msp2.csv` is the defensible one). The 26 per-peptide
`--plot y` PNGs are not — they are regenerable, and gitignored.

To regenerate the whole SIS chain from the raw export:

```bash
python -m src.prep_skyline_sis "<export.csv>" "<document.sky>" <out_dir>
python tools/matrix-matched_calcurves/bin/calculate-loq.py \
    <out_dir>/sis_ultra_skyline.csv <out_dir>/sis_ultra_map.csv \
    --multiplier_file <out_dir>/sis_ultra_multipliers.csv \
    --model auto --plot n --min_noise_points 0 --min_saturation_points 2 \
    --output_path <out_dir>/fom
python figures/diag_sis_curves.py --fom <out_dir>/fom/figuresofmerit.csv
```

## Reusable: never subset a report to a handful of peptides

Recorded 2026-09-03. Measured on `260521_RH4_Calcurves` — RH4, human+yeast database,
36 runs, 12 curve points x 3 replicates, DIA-NN `diann_report.tsv`. Line numbers in
this section and the next refer to submodule commit `61414c6`.

To plot six PAX3/FOXO1 precursors, the report was filtered down to just those
precursors' rows. That produces **wrong figures of merit**, silently — no note, no
warning, `n_curvepoints` still reads 12.

The six are identified in only 17 of the 36 runs, and in **none of the three blanks**.
`_complete_grid` ([`bin/calculate-loq.py:77`](../tools/matrix-matched_calcurves/bin/calculate-loq.py#L77))
builds its dense grid from the runs present in the input and, by design, refuses to
invent a column for a run that is absent — its docstring is explicit that doing so
"would fabricate data rather than restore it". The behaviour is documented; the failure
mode is not. Filtering by peptide removes runs as a side effect, and with the blanks
gone the noise plateau the LOD is fit from goes with them.

Against the full 55,717-precursor run:

| precursor | LOD, subset | LOD, full | |
|---|---|---|---|
| TTFTAEQLEELER2 | 0.0158 | 0.0496 | -68% |
| SGFPLEVSTPLGQGR2 | 0.0434 | 0.0662 | -34% |
| LPSDLDGMFIER2 | 0.0545 | 0.0743 | -27% |
| SSWWMLNPEGGK2 | `inf` | 0.0353 | flipped |
| YQETGSIRPGAIGGSKPK3 | `inf` | 0.0770 | flipped |
| ENPGMFSWEIR2 | 0.3000 | 0.3000 | unchanged |

`slope_linear` was low by 3-13% as well, so this is not confined to the noise segment.

**How to subset safely.** Add one carrier precursor identified in all runs, so the run
set is restored, then confirm the subset reproduces the full run before trusting it.
Doing that here made all six agree with the full run exactly, on every column. A
carrier is easy to find:

```bash
# precursors seen in all 36 runs
awk -F'\t' 'NR>1{k=$3; if(!(k SUBSEP $1 in s)){s[k SUBSEP $1]=1; c[k]++}}
            END{for(k in c) if(c[k]==36) print k}' diann_report.tsv | head
```

## Settled: `min_noise_points` reports the *most* sensitive peptides as `inf`

Recorded 2026-09-03. Second sighting of the effect already noted for the SIS dataset
above, on unrelated data and with the code path pinned down.

[`bin/calculate-loq.py:520-522`](../tools/matrix-matched_calcurves/bin/calculate-loq.py#L520)
forces **both** `LOD` and `stndev_noise` to `inf` when fewer than `min_noise_points`
distinct curve points fall *below* the computed LOD:

```python
elif df["curvepoint"][~mask].nunique() < min_noise_points:
    # if there's not enough below the LOD
    lod_results = [np.inf, np.inf]
```

At the mnp=2 default, a peptide sensitive enough that its LOD lands below the
second-lowest dilution is rejected for lack of noise points. In `figuresofmerit.csv`
it is indistinguishable from a peptide that never rose out of the noise — same `inf`,
empty `notes`.

Two of ten SIX1/SIX2 precursors on `260521_RH4_Calcurves` hit this. Re-running with
`--min_noise_points 1`:

| precursor | LOD default | LOD mnp=1 | LOQ default | LOQ mnp=1 |
|---|---|---|---|---|
| FLWSLPAC(UniMod4)DHLHK3 | `inf` | 1.3e-11 | `inf` | **0.0202** |
| LQQLWLK2 | `inf` | 1.1e-11 | `inf` | **0.0505** |
| AVVAFHR2 | 0.0184 | 0.0184 | `inf` | `inf` |
| TIWDGEETSYC(UniMod4)FK2 | 0.1495 | 0.1495 | `inf` | `inf` |
| ILESHQFSPHNHPK3 | 0.0118 | 0.0118 | 0.0517 | 0.0517 |

The other eight are bit-identical, which confirms the gate is what fired rather than a
different fit. The two suppressed precursors have the **best** LOQs of the ten; their
LODs sit below the lowest point tested (0.005), so the curve bottoms out above their
detection limit.

Consequence for the paper: an `inf` LOD count is not a count of undetectable peptides.
Anywhere we report finite-LOD fractions, the denominator mixes two populations —
too noisy to detect, and too sensitive for the curve's dynamic range. Separating them
needs either a more dilute bottom point or an mnp sweep.

## Settled: the search database sets the blank, and the blank sets the LOD

Recorded 2026-09-03. Same 36 RH4 acquisitions searched two ways.

Residual human signal at the 0% point, as a fraction of the 100% point:

| search database | human at blank | ratio |
|---|---|---|
| human only | 4.5e9 / 1.09e11 | 4.1% |
| human + yeast | 6.15e8 / 6.15e10 | **1.0%** |

Adding the second proteome removes about three quarters of the apparent human signal
in a sample containing no human analyte — that signal was matrix peptides with nowhere
correct to go. Human identifications at the blank are 2,395, against 157,526 at full
strength.

**The two rows are not a matched measurement.** The human-only figures are protein-group
totals carried over from a separate analysis of `260526_calcurve-rodeo_RH4_EXP61` and
were not re-measured here; the human+yeast figures are sums of `Precursor.Quantity` over
identified precursors, measured directly. The direction and rough magnitude are solid,
the exact percentages are not comparable.

Why it matters here rather than as a general search-quality point: LOD and LOQ are fit
from the blank. A database that leaves matrix peptides unassigned inflates the noise
plateau, and every LOD downstream of it is correspondingly conservative. Any
figures-of-merit comparison across datasets has to hold the database constant, or say
that it does not.
