# TODO — what's left

Status of the paper's figures + analysis. All figure scripts run; the four
FOM-only panels reproduce pixel-identically, the three raw-data figures reproduce
the same content (see `run_all.py --with-raw`).

## Data / analysis
- [ ] **Regenerate the remaining FOMs against the current tool.** The submodule is now
      bumped to `ffb1087` (was `ac0b951`, three commits behind). The two EncyclopeDIA
      datasets are done and came back **numerically identical** — `3d18164` (DIA-NN
      densification) and `ffa5118` (canonical row order) have nothing to bite on in an
      already-dense wide-format matrix, confirmed by a control run at the old pin. The
      six Bruker DIA-NN FOMs are where the movement lives and are still stale; see
      [`fom_provenance.md`](fom_provenance.md) for the per-CSV recipe and status.
      Budget ~0.8 peptides/s per worker: the two 1.2 GB reports are ~4 h each at 6
      threads and must run one at a time on a 16 GB box.
- [ ] **Decide stratified vs case resampling before regenerating the Bruker FOMs.**
      See [`loq_stability_note.md`](loq_stability_note.md) — `ffb1087` still uses
      unstratified case resampling, and switching would invalidate everything
      regenerated under it. Under discussion as of 2026-08-12.
- [ ] **Separate Bruker dataset with a ULOQ** — locate it and run it through the
      tool; decide whether it earns a panel or a supplement.
- [ ] *(optional)* Regenerate `asms/ultra` with the **new code at mnp=2** so Fig 3's
      LOQ comes from a current-tool FOM instead of `legacy_mnp2/bruker_ultra.csv`
      (the original pre-fix output). Cosmetic — values are close.

## Manuscript text
- [ ] **Methods:** describe the algorithm changes — trilinear (noise+linear+
      saturation) fit + `--min_saturation_points` + the ULOQ figure of merit; the
      `min_noise_points=0` rescue fix (compare vs lowest *nonzero* curve point);
      never-drop-a-peptide noted rows.
- [ ] **Results / Fig 2B:** the manuscript's "≈24% increase, Ultra n=5756 / Ultra II
      n=7137" is from the ORIGINAL (pre-improvement) code. Update to the new numbers:
      detection **+36%**, quantifiable **+7–10%** (finite LOD 6,614→7,084; finite LOQ
      5,885→6,449). Clarify which metric the "increase" refers to.
- [ ] Confirm the instrument nickname mapping for Fig 2B: **Timbaux = Ultra**,
      **Desnaux = Ultra II** (both filenames say "Ultra" = product family).

## Figures
- [ ] Coauthor review of example-peptide selection and panel labels.
- [ ] Decide final home/order of the two supplements (before/after, detection-floored
      tiers) and the ULOQ linear-examples / distribution panels.

## Repo
- [ ] Push the figure-port commit.
- [ ] Add a `LICENSE`.
- [ ] Point `config.yaml → raw:` at wherever the MS data lives for anyone reproducing
      the raw-data figures (`run_all.py --with-raw`).

## Done
- [x] Fig 1B (distributions + examples), Fig 2A/2B/2C LOD/LOQ/ULOQ triptychs, Fig 3,
      before/after + tier supplements, ULOQ linear examples.
- [x] Tool changes (trilinear/ULOQ, mnp=0 rescue, noted rows) committed + pinned.
- [x] Repo scaffolded, smoke-tested, initial commit pushed.
