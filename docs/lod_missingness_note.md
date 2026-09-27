# Note: how search-engine missingness sets the LOD

Working note for the methods/limitations discussion. Recorded 2026-09-27. Companion to
[`fit_weighting_note.md`](fit_weighting_note.md) and
[`loq_stability_note.md`](loq_stability_note.md).

## The finding

The LOD depends on how the search engine reports low-concentration runs, not only on the
data. EncyclopeDIA integrates the expected retention-time / *m/z* window at every run, so it
returns a real quantity even where the peptide does not clear a detection score. DIA-NN
writes a row only where the precursor was identified, so a peptide that drops out of the
low-concentration runs has no row there at all.

The tool fills every missing peptide x run cell with zero (`_complete_grid` and
`_normalize_input`, `bin/calculate-loq.py:96-142`; the code carries its own
`# TODO: is this appropriate?` at the fill step). For EncyclopeDIA this changes little,
because the matrix is already dense. For DIA-NN it rebuilds the low end out of fabricated
zeros. That is fine when a peptide truly drops out at the bottom of an otherwise-sampled
curve, but for a peptide seen at only one or two levels it manufactures the whole noise
plateau and hands back a confident LOD.

## Evidence

Share of `area == 0` cells by dilution level, after the same read + densify step:

| level | Exploris (EncyclopeDIA) | Bruker Ultra (DIA-NN) |
|---|---|---|
| bottom (0.007-0.01) | 22-26% | 97-98% |
| 0.03-0.1 | 22-33% | 80-93% |
| 0.3 | 22% | 53% |
| 0.7 | 22% | 22% |
| top (1.0) | 22% | 6% |
| **overall** | **16.5%** | **65.9%** |

EncyclopeDIA keeps a real (nonzero) noise plateau at the low end, so the noise/linear
intersection is continuous. DIA-NN's bottom is almost all fabricated zeros.

**Why the zeros pin the LOD to a dilution level.** The LOD is the concentration where the
linear segment meets the noise level (`bin/calculate-loq.py:523`,
`LOD = (b_noise + k*std_noise - b_linear) / m_linear`). For a sparse peptide the noise
level is zero (all filled zeros) and the "linear" segment is a line anchored on those zeros
up to the one or two detected points. A line drawn through zero-valued points crosses the
(also-zero) noise level at the last all-zero dilution step. So the LOD lands on the grid,
reporting "the highest concentration at which the peptide was still undetected" rather than
a fitted detection limit. Two examples from `main/bruker_ultra` (LOD with vs without the
zero-fill):

| peptide | real detections | LOD with zero-fill | LOD without (real points only) |
|---|---|---|---|
| STAGDTHLGGEDFDNR | 1 / 10 | 0.70 | cannot fit (too few points) |
| SLGVSNFNR | 2 / 10 | 0.30 | inf (not determinable) |
| AC(UniMod4)ANPAAGSVILLENLR | 9 / 10 | 0.029 | 0.018 (no pinning) |

Strip the zeros and the sparse peptides give no LOD, which is the honest answer. The
well-sampled peptide barely moves.

## Scope

Affects the DIA-NN inputs (`bruker_ultra`, `bruker_ultraII`, `bruker_60spd`,
`bruker_100spd`): ~99% of Ultra / Ultra II peptides have a zero-variance noise block, and
~20% of their finite LODs sit on the single value 0.70 (the second-highest dilution step).
The EncyclopeDIA inputs (`exploris_dia`, `il15_prm`) do not show it. `bruker_60spd_pr`
(DIA-NN `pr_matrix`, normalized) piles a different way (~55% near 0.68) because the
normalization, not the fit, sets the low end.

The interpolated LOQ readout (matrix-matched_calcurves#21) does not touch this. LOQ is
smoother because it interpolates a CV crossing; the LOD grid-pinning comes from the input
missingness plus the zero-fill, not from a readout grid.

## Consequence for the manuscript

The detection gain reported for the Bruker DIA-NN sets is inflated by these zero-filled
peptides. On the same `bruker_ultra` input, finite LODs rise from 6,614 (no densification)
to 48,754 (current) — most of the added "detections" are peptides seen at only one or two
real levels. Any detection-count comparison across the DIA-NN sets should either state this
caveat or apply a gate.

Two fixes, both input-aware, neither yet applied:
1. **Detection gate** — require at least N real (reported, nonzero) levels before returning
   a finite LOD/LOQ. Removes the one- and two-point peptides and deflates the detection
   count to what was actually measured.
2. **Do not treat unreported DIA-NN cells as measured zeros** for the noise-plateau
   estimate (the `bin/calculate-loq.py:141` TODO). Keeps the densification where it helps
   (genuine low-end dropout) without fabricating a plateau from nothing.

## Reproducing

Against the committed FOM and the raw inputs in `config.yaml`. The `area == 0` share per
level comes from reading each input through the tool (`load_tool().read_input(path, map)`)
and grouping by `curvepoint`. The with/without-zero-fill LOD comes from refitting a peptide
on its densified grid versus its real detected points only (`fit_by_lmfit_yang` +
`calculate_lod`). EncyclopeDIA vs DIA-NN inputs are named in
[`figure_manifest.md`](figure_manifest.md).
