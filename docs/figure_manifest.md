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
| SIS peptide response | timsTOF Ultra, 17 min, 5-pt 3-fold | Skyline `Peptide Total Area Fragment` | **0** (never reaches noise) | *diagnostic only — no ULOQ; see [`diagnostics.md`](diagnostics.md)* |

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
| ULOQ examples (linear) | `figures/fig_uloq_examples_linear.py` | Exploris/IL15 `.elib` + Ultra II `CURVES_pep` | yes |
| *diag:* SIS curves | `figures/diag_sis_curves.py` | SIS Skyline export + `.sky` | yes |

## Provenance notes
- **Ultra vs Ultra II (2B)** are method-matched (60 SPD, 100 ng, 75 ms, 250310),
  distinguished only by instrument nickname: **Timbaux = Ultra**, **Desnaux = Ultra II**.
  Both filenames say "Ultra" (product family).
- **2A/2C** use the asms quantms pipeline; **2B** uses the 250320 separate-search
  pipeline — different experiments, so absolute counts are NOT comparable across
  panels; within-panel comparisons are.
- ULOQ counts are `min_noise_points`-independent (ULOQ derives from the trilinear
  fit, not the noise-point logic). They are **not** `min_saturation_points`-independent,
  and that sensitivity has not been swept for these datasets — see
  [`diagnostics.md`](diagnostics.md) for the audit tool and the open question.
