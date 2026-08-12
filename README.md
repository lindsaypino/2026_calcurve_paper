# Matrix-matched calibration curves — paper figures

Figure-generation code for the JASMS matrix-matched calibration curves paper. The
LOD/LOQ/**ULOQ** algorithm itself lives in a separate repository, pinned here as a
submodule ([`matrix-matched_calcurves`](https://github.com/lindsaypino/matrix-matched_calcurves)).

**Tool version.** The submodule is pinned at **`ffb1087`** ("Restore piecewise's
legacy LOD edge cases; keep peptides with no multiplier", 2026-08-11). Not every
figure of merit in `data/figuresofmerit/` has been regenerated against that pin yet —
[`docs/fom_provenance.md`](docs/fom_provenance.md) records, per CSV, the exact input,
concentration map, `min_noise_points`, and producing commit, and marks the ones that
are still stale.

> **What's left to do:** see [`docs/TODO.md`](docs/TODO.md). All figure scripts run
> (`python run_all.py --with-raw`); remaining work is mostly manuscript text + the
> outstanding Bruker ULOQ dataset.

## Layout

```
config.yaml              central paths (edit `raw:` for your machine)
requirements.txt         figure-making dependencies
run_all.py               regenerate every figure
src/                     shared code
  style.py               vizta 'talk' theme + talusbio palette (the ONE style source)
  fom_io.py              load figures-of-merit CSVs + config
  prep_curves_pep.py     reshape the Bruker separate-search CURVES_pep files
figures/                 one script per panel (import src)
data/
  figuresofmerit/        LOD/LOQ/ULOQ outputs — GITIGNORED (manuscript supplemental
                         tables); reproduce locally with `build_fom.py`
    main/                mnp=2 (Exploris mnp=2, IL15 mnp=0) — the paper's main panels
    supp_mnp0/           DIA-NN datasets re-run at min_noise_points=0 (supplement)
    legacy_mnp2/         original pre-fix outputs (before/after; NOT reproducible here)
  maps/                  filename -> concentration maps (tracked; recipe inputs)
  raw/                   (gitignored) large MS data — see config.yaml
tools/matrix-matched_calcurves/   the LOD/LOQ tool (submodule, pinned)
output/                  generated PNGs (gitignored)
docs/figure_manifest.md  panel -> script -> data -> settings
```

## Setup

```bash
git clone --recurse-submodules <this repo>
# or, if already cloned:  git submodule update --init --recursive
python -m venv .venv && . .venv/Scripts/activate   # Windows; use bin/activate on *nix
pip install -r requirements.txt
```

## Regenerate figures

Figures read `data/figuresofmerit/` (gitignored). First populate it — drop in the
manuscript supplemental FOM tables, or regenerate from raw data (see below):

```bash
python run_all.py             # panels that read only FOM CSVs
python run_all.py --with-raw  # also the example-curve / Fig 3 figures (need raw MS data)
```

Once the FOM CSVs are present, most panels regenerate in seconds — no multi-GB
DIA-NN reports or multi-hour LOD/LOQ runs. Only the example-curve figures (Fig 1B
examples, Fig 3, linear ULOQ examples) refit raw data and need `config.yaml → raw:`
to point at the MS data on your machine.

## Regenerating the figures-of-merit

The FOM CSVs (`data/figuresofmerit/`, gitignored) are distributed with the manuscript
as supplemental tables, alongside the raw `diann_report` outputs. To regenerate one
from its raw input, run the pinned tool:

```bash
python tools/matrix-matched_calcurves/bin/calculate-loq.py \
    <input> <map.csv> --model auto --plot n \
    --min_noise_points <N> --min_saturation_points 2 --output_path <dir>
# then rename/move the resulting figuresofmerit.csv to data/figuresofmerit/<...>.csv
```

Per dataset (paths relative to the Google Drive roots in `config.yaml`; `A` =
`2025_calcurves_asms/data`, `B` = `2023_bruker_ultra/data`):

| output CSV | dataset directory | input file | `min_noise_points` |
|---|---|---|---|
| `main/exploris_dia.csv`   | `A/Exhausted_CD8_DimethylBackground_DIA_Exploris480` | `calcurve2_KBB.elib.peptides.txt` | 2 |
| `main/il15_prm.csv`       | `A/IL15_IL2_stimulated_CD8_Curves` | `50perCycle_Quant.elib.peptides.txt` | **0** (single injection) |
| `main/bruker_60spd.csv`   | `B/asms/60SPD`    | `diann_report_noYEAST.tsv` | 2 |
| `main/bruker_100spd.csv`  | `B/asms/100SPD`   | `diann_report_noYEAST.tsv` | 2 |
| `main/bruker_60spd_pr.csv`| `B/asms/60SPD_pr` | `diann_report.pr_matrix_noYEAST.tsv` | 2 |
| `main/bruker_ultra.csv`   | `B/250320_Timbaux_Ultra_HumanYeast_22min`   | `Timbaux_60SPD_DDM_CURVES_pep.tsv` † | 2 |
| `main/bruker_ultraII.csv` | `B/250320_Desnaux_UltraII_HumanYeast_22min` | `UltraII_60spd_100ng_CURVES_pep.tsv` † | 2 |
| `supp_mnp0/bruker_{60,100}spd.csv` | `B/asms/{60,100}SPD` | `diann_report_noYEAST.tsv` | **0** |

Maps are in `data/maps/`. **†** The two Bruker hardware datasets are hand-merged
per-point `CURVES_pep` files; reshape into a DIA-NN report + A–N map first:
`python -m src.prep_curves_pep <CURVES_pep.tsv> report.tsv map.csv`, then run the tool
on `report.tsv map.csv`.

- **Ultra vs Ultra II** are method-matched (60 SPD, 100 ng, 250310); the instrument
  identity is only the nickname — **Timbaux = Ultra**, **Desnaux = Ultra II**.
- `legacy_mnp2/` (the before/after "before") is from the ORIGINAL pre-fix tool;
  reproduce it only by checking out an earlier tool commit.
- Write tool output to **local disk**, not Google Drive — Drive rejects the tool's
  per-row `fsync` under concurrent writers (`OSError: Errno 22`).
