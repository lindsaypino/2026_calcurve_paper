# Figure manifest

Maps each paper panel to its generator, data source, and processing settings.
All LOD/LOQ/ULOQ come from `calculate-loq.py` (`tools/matrix-matched_calcurves`,
`--model auto`). Committed FOM CSVs live in `data/figuresofmerit/`.

## Processing settings by dataset

| dataset | instrument / method | search input | min_noise_points | committed FOM |
|---|---|---|---|---|
| Exploris DIA | Thermo Exploris 480, 15-pt | EncyclopeDIA `.elib.peptides.txt` | 2 (default) | `main/exploris_dia.csv` |
| IL15/IL2 PRM | targeted PRM, 9-pt, single inj. | EncyclopeDIA `.elib.peptides.txt` | **0** (single injection) | `main/il15_prm.csv` |
| Bruker 60SPD | Ultra II, 60 SPD | DIA-NN `diann_report` (quantms) | 2 | `main/bruker_60spd.csv` |
| Bruker 100SPD | Ultra II, 100 SPD | DIA-NN `diann_report` (quantms) | 2 | `main/bruker_100spd.csv` |
| Bruker 60SPD_pr | Ultra II, 60 SPD | DIA-NN `pr_matrix` (normalized) | 2 | `main/bruker_60spd_pr.csv` |
| Bruker Ultra | timsTOF Ultra (Timbaux), 60 SPD | separate-search `CURVES_pep` | 2 | `main/bruker_ultra.csv` |
| Bruker Ultra II | timsTOF Ultra II (Desnaux), 60 SPD | separate-search `CURVES_pep` | 2 | `main/bruker_ultraII.csv` |

`supp_mnp0/` = same DIA-NN datasets re-run at `min_noise_points=0`; `legacy_mnp2/`
= the original pre-fix outputs (default mnp=2) for the before/after supplement.

## Panels

| panel | script | inputs | needs raw? |
|---|---|---|---|
| Fig 1B examples | `figures/fig01b_examples.py` | Exploris + IL15 `.elib`; refits peptides | yes |
| Fig 1B distributions | `figures/fig01b_distributions.py` | `main/exploris_dia.csv` | no |
| Fig 2A gradient | `figures/fig02_triptychs.py` | `main/bruker_{60,100}spd.csv` | no |
| Fig 2B hardware | `figures/fig02_triptychs.py` | `main/bruker_ultra{,II}.csv` | no |
| Fig 2C software | `figures/fig02_triptychs.py` | `main/bruker_60spd{,_pr}.csv` | no |
| Fig 3 RT/mz/abund | `figures/fig03_retentiontime.py` | `legacy_mnp2/bruker_ultra.csv` + ultra `diann_report` | yes |
| Fig 3 LOQ model | `figures/fig03_loq_model.py` | ultra report-derived features | yes |
| Supp before/after | `figures/supp_before_after.py` | `legacy_mnp2/` + `supp_mnp0/` | no |
| Supp tiers | `figures/supp_tiers.py` | `supp_mnp0/` | no |
| Supp bootstrap calibration | `figures/supp_bootstrap_calibration.py` | simulation only (ground truth known) | no |
| Supp curve spacing | `figures/supp_curve_spacing.py` | simulation + submodule `data/one_protein.csv` | no |
| Supp LOQ readout | `figures/supp_loq_readout.py` | simulation only (ground truth known) | no |
| ULOQ examples (linear) | `figures/fig_uloq_examples_linear.py` | Exploris/IL15 `.elib` + Ultra II `CURVES_pep` | yes |

## Provenance notes
- **Ultra vs Ultra II (2B)** are method-matched (60 SPD, 100 ng, 75 ms, 250310),
  distinguished only by instrument nickname: **Timbaux = Ultra**, **Desnaux = Ultra II**.
  Both filenames say "Ultra" (product family).
- **2A/2C** use the asms quantms pipeline; **2B** uses the 250320 separate-search
  pipeline — different experiments, so absolute counts are NOT comparable across
  panels; within-panel comparisons are.
- ULOQ counts are `min_noise_points`-independent (ULOQ derives from the trilinear
  fit, not the noise-point logic).

## Methods supplements (simulation)

Two supplements answer methods questions rather than showing results, and
need no raw MS data. Both cache their simulations to `output/` (gitignored);
pass `--force` to recompute.

- **Bootstrap calibration** - simulates experiments from a known truth and asks
  which resampling scheme recovers the real sampling variability. The scheme the
  tool already uses (case) is the best calibrated; stratified, wild and Bayesian
  resampling all understate uncertainty by 15-20%. The LOQ spread across
  replicate experiments is larger than any difference between schemes.
- **Curve spacing** - the same 42 injections spent on four dilution designs, plus
  the sample dataset thinned to a coarser design. Log spacing beats linear by
  ~5x on achievable LOQ; thinning the low end costs most peptides their LOQ.
  Also shows the readout grid effect: reading the CV on a uniform grid, as the
  tool does, reports fewer peptides as resolved than log or measured-level
  spacing.

- **LOQ readout** - scores grid spacing and the crossing rule against a known
  truth. A log grid roughly halves the bias of genuine crossings versus the
  uniform grid the tool uses (+56% vs +93%) and cuts the dependence on the
  arbitrary point count from 15.9% to 2.6%. In a scenario with no true LOQ at
  all, the current rule invents one in 100% of experiments while an
  interpolated crossing with an explicit no-crossing outcome essentially never
  does.

See `docs/loq_grid_and_resampling_note.md` for the numbers and caveats.
