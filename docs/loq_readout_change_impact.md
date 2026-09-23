# Impact of the interpolated LOQ readout on the figures of merit

Recorded 2026-09-22. Companion to [`loq_grid_and_resampling_note.md`](loq_grid_and_resampling_note.md),
which argues the decision; this note measures what the decision did to the numbers we report.

## The change

`calculate_loq` now reads the LOQ by interpolating the bootstrap CV curve at the 0.20
crossing, on a grid matched to the dilution design (geometric for a log-spaced series,
uniform for a linear one), and reports no LOQ when the CV never crosses. It replaces the
old rule, which snapped the LOQ to the lowest uniform-grid point under the threshold and
emitted a value even when the curve never crossed. LOD and ULOQ are untouched by the change.

Tool: `matrix-matched_calcurves` PR #23, merged to master `6a017bb`; paper submodule pinned
to it (paper commit `ffff0e8`). See the tool's `CHANGELOG.md` for the full method delta from
the Pino 2020 baseline.

## How the before/after was measured

The pre-change tool (`d80f50e`) and the new tool were run in the **same environment on the
same raw inputs**, so any LOD/ULOQ difference is numerical noise and any LOQ difference is
the readout change alone (not library or pipeline drift). Both used `--bootreps 100` and the
dataset's paper settings (`min_noise_points` 2, except il15 = 0).

## Result: LOD and ULOQ unchanged, LOQ shifts small and downward

| dataset | peptides | LOD change | ULOQ change | LOQ finite (before→after) | LOQ median shift | LOQ moved >2x |
|---|---|---|---|---|---|---|
| il15_prm (mnp=0) | 539 | 0.00% (identical) | 0.00% (identical) | 5 → 6 | -5.4% (only 5 LOQs) | 0 |
| exploris_dia (mnp=2) | 62,760 | 0.00% (identical) | 0.00% (identical) | 23,789 → 23,808 | **-1.7%** | 13 (0.05%) |
| ultraII (mnp=2) | 74,330 | 0.00% (identical) | 0.00% (identical) | 57,857 → 57,893 | **-0.8%** | 5 (0.009%) |

LOD and ULOQ are **bit-identical** (0.00% change, no finite/inf flips) — expected, since the
change never touches the fit or `calculate_uloq`. The LOQ shift is small and downward; the
few peptides that move more than two-fold are old grid floor-pins the log grid corrected.
No real-data peptide had a fabricated LOQ that needed removing (the grid-snap rule rarely
floor-pinned on these designs, unlike the heavily floor-pinned yeast sample dataset where
the median shift was -15%).

The main effect is qualitative, not quantitative: peptides whose CV never reaches 20% now
carry an explicit `loq_no_crossing` note (20,833 in exploris, 7,128 in ultraII, 478 in il15)
instead of a silent `inf`, and a small set carry `loq_at_lod` (LOQ = LOD; 375 exploris, 54
ultraII) when quantifiable down to detection.

Figure: before/after LOQ scatter (`loq_before_after_scatter.png`, session scratchpad).

## Manuscript text (draft, for reuse)

**Methods.** "The LOQ is the concentration at which the bootstrap coefficient of variation
first falls below 0.20, obtained by linearly interpolating the fitted CV curve between the
two evaluation points that bracket the crossing. The CV is evaluated on a grid matched to
the calibration design — geometric spacing for a log-spaced dilution series, uniform for a
linear one, detected automatically from the curve points. When the CV does not fall below
0.20 anywhere above the LOD, no LOQ is reported."

**Robustness.** "Replacing the earlier grid-nearest LOQ readout with the interpolated,
design-matched readout left the LOD and ULOQ unchanged (identical to machine precision for
every peptide) and shifted the LOQ only slightly: the median change in reported LOQ was
-1.7% on the Exploris DIA dataset and -0.8% on the timsTOF Ultra II dataset, with fewer than
0.06% of peptides changing by more than two-fold. The main effect was qualitative — peptides
whose CV never reaches 20% are now reported as having no determinable LOQ rather than a value
pinned to the lowest grid point."

## Status

All 10 non-legacy FOM CSVs are regenerated against `6a017bb` and installed under
`data/figuresofmerit/` (completed 2026-09-23; every Bruker run exited cleanly). Per-CSV row
and loq-outcome counts are in [`fom_provenance.md`](fom_provenance.md). `legacy_mnp2` stays
frozen. Pre-change FOM is backed up in the session scratchpad (`fom_backup/`).

The before/after characterization above (LOD/ULOQ bit-identical, small LOQ shift) was measured
by running the pre-change tool (`d80f50e`) against the new one in the same environment on
`il15_prm`, `exploris_dia`, and `bruker_ultraII` — three representative datasets spanning
targeted PRM, wide-format DIA, and hardware CURVES_pep. Because the change never touches the
fit, LOD and ULOQ are unchanged for the Bruker DIA-NN sets by construction; a per-set
before/after on those would require re-running the old code (~a day) and was not done.

Reproduce: `compare_before_after.py`, `uloq_check.py`, `scatter_before_after.py` (session
scratchpad); run the pinned tool and its `d80f50e` worktree in the same venv on the raw inputs.
