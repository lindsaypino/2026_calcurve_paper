# Where things stand — 2026-08-03

Session handoff. Delete this file once the queue below is drained.

## Read this first: nothing is committed

Working tree has 7 uncommitted changes on `master` (in sync with `origin/master`):

```
 M config.yaml              # + SIS raw paths (currently a local Downloads path)
 M docs/TODO.md             # ULOQ item closed; follow-ups + stale items folded in
 M docs/figure_manifest.md  # SIS row, diag script, min_saturation_points caveat
?? docs/diagnostics.md      # NEW - the diagnostics note + open questions
?? figures/diag_sis_curves.py
?? src/prep_skyline_sis.py
?? src/uloq_audit.py
```

Review and commit when you're ready. Worth branching rather than committing straight to
`master`. Nothing here is half-finished — every script runs — but none of it is reviewed.

## What happened this session

1. **Submodule initialized.** `tools/matrix-matched_calcurves` is checked out at the
   pinned `ac0b951`. Also installed `lmfit` + `pyteomics`, which were missing and were
   breaking `load_tool()`. They went into **system Python 3.13.3** — there is still no
   `.venv` despite the README describing one.
2. **The SIS ULOQ dataset is settled: it has no ULOQ in it.** Wrote
   `src/prep_skyline_sis.py` to make the Skyline export readable by the tool, got it
   running end to end, and confirmed the saturation clause fails 26/26 because the curves
   are still climbing at the top spike level. Not a settings problem. Full write-up in
   [`diagnostics.md`](diagnostics.md).
3. **Wrote a reusable clause audit**, `src/uloq_audit.py` — answers "why did this dataset
   get no ULOQ" for any melted `(peptide, curvepoint, area)` frame. Validated: predicts
   the tool's accepted set exactly at both msp=2 and msp=1.
4. **Filed** https://github.com/lindsaypino/matrix-matched_calcurves/issues/14 (the ULOQ
   legend bug — cosmetic, verified with an isolated repro).

## Queue, roughly in priority order

1. **Populate `data/figuresofmerit/`.** This is the top blocker and it gates two separate
   things: all four FOM-only figures currently fail with `FileNotFoundError`, *and* the
   `n_sat` audit below can't run. Drop in the supplemental tables or regenerate per the
   README's per-dataset table.
2. **Audit `n_sat` behind the existing ULOQ panels** (Fig 1B distributions, Fig 2
   triptychs, and especially the peptides `fig_uloq_examples_linear.py` chooses to
   feature). Use `python -m src.uloq_audit <melted.csv>`. Risk is lower than for the SIS
   curve — the real curves have 9–15 points vs 5 — but it is untested, and a featured
   example resting on a thin plateau is a soft target in review.
3. **Send the collaborator email.** Draft is at `output/sis_run/sis_collaborator_update.md`
   and needs only the `[Name]` filled in. It asks whether spike levels above the current
   top point are feasible — which is the only route to a real ULOQ from this experiment.
   Deliberately does *not* mention the `--min_saturation_points 1` workaround; keep it
   that way.
4. **Decide the SIS dataset's fate.** No ULOQ panel. Two defensible supplements: "the
   Ultra's linear range extends past the top of a 5-point SIS curve", and/or the mnp=0
   rescue demonstration (6/26 → 19/26 finite LODs). Both figures already exist in
   `output/`.
5. **Consider a second issue on the `std_sat` inversion** — `calculate_uloq` makes a
   noisier top point yield a *narrower* quantifiable range, which is backwards. More
   substantive than the legend bug, but wants evidence from a real 14-point curve first,
   so it's gated on item 1.
6. **Methods text** now has three extra things to state (the ≥5-distinct-point gate for
   trilinear fits; that LOQ thresholds the CV of resampled *means*, ~1.7x off raw
   replicate CV at n=3; and that LOQ is CV-gated while ULOQ is purely geometric). Already
   noted in `TODO.md`.
7. **Housekeeping carried over:** no `LICENSE` (the tool is Apache-2.0 if you want to
   match); `figure_manifest.md` lists a `fig03_loq_model.py` that doesn't exist and the
   README references a `build_fom.py` that doesn't either; `config.yaml`'s SIS paths point
   at `Downloads` and should move to the Drive root.

## Artifacts

Preserved under `output/sis_run/` (gitignored, but on disk — the session temp dir they
were generated in will be cleaned up):

| file | what |
|---|---|
| `sis_ultra_skyline.csv`, `sis_ultra_map.csv`, `sis_ultra_multipliers.csv` | tool-ready inputs from `prep_skyline_sis` |
| `sis_melted.csv` | melted frame for `src.uloq_audit` |
| `fom_mnp0_msp2.csv` | **the defensible run** — 19/26 LOD, 18/26 LOQ, 0 ULOQ |
| `fom_mnp2_msp2.csv` | paper default, for the mnp comparison (6/26 LOD) |
| `fom_mnp0_msp1.csv` | the relaxed run that yields 7 artifact ULOQs |
| `tool_plots/` | all 26 built-in `--plot y` PNGs |
| `sis_collaborator_update.md` | the email draft |
| `issue_body_clean.md` | body as filed on issue #14 |

Plus `output/diag_sis_curves.png` (26 curves, log-log, per-group) and
`output/diag_sis_uloq_toolfits.png` (the 7 ULOQ peptides' tool fits, tiled).

To regenerate the whole SIS chain from scratch:

```bash
python -m src.prep_skyline_sis "<export.csv>" "<document.sky>" <out_dir>
python tools/matrix-matched_calcurves/bin/calculate-loq.py <out_dir>/sis_ultra_skyline.csv <out_dir>/sis_ultra_map.csv --multiplier_file <out_dir>/sis_ultra_multipliers.csv --model auto --plot n --min_noise_points 0 --min_saturation_points 2 --output_path <out_dir>/fom
python figures/diag_sis_curves.py --fom <out_dir>/fom/figuresofmerit.csv
```
