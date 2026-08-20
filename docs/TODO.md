# TODO — what's left

Status of the paper's figures + analysis. All figure scripts run; the four
FOM-only panels reproduce pixel-identically, the three raw-data figures reproduce
the same content (see `run_all.py --with-raw`).

## Data / analysis
- [x] **Separate Bruker dataset with a ULOQ** — located (`240610_SIS_peptide_response_
      Ultra_BC_refined`) and run through the tool. **It has no ULOQ in it**: the
      saturation clause fails 26/26 because the curves are still climbing at the top
      spike level (top step 2.2–2.8x vs nominal 3.0). Not a settings problem — needs
      acquisition above the current top point. Full write-up in
      [`docs/diagnostics.md`](diagnostics.md).
- [ ] **Decide the SIS dataset's fate** — no ULOQ panel, but two defensible supplements:
      "the Ultra's linear range extends past the top of a 5-point SIS curve", and/or the
      mnp=0 rescue demonstration (6/26 → 19/26 finite LODs). Collaborator has been asked
      whether higher spike levels are feasible.
- [ ] **Audit `n_sat` behind the existing ULOQ panels** (Fig 1B distributions, Fig 2
      triptychs, the ULOQ example peptides). Untested, and an example peptide resting on
      a thin plateau would be a weak feature. Blocked on populating
      `data/figuresofmerit/`. See [`docs/diagnostics.md`](diagnostics.md).
- [ ] *(optional)* Regenerate `asms/ultra` with the **new code at mnp=2** so Fig 3's
      LOQ comes from a current-tool FOM instead of `legacy_mnp2/bruker_ultra.csv`
      (the original pre-fix output). Cosmetic — values are close.

## Manuscript text
- [ ] **Methods:** describe the algorithm changes — trilinear (noise+linear+
      saturation) fit + `--min_saturation_points` + the ULOQ figure of merit; the
      `min_noise_points=0` rescue fix (compare vs lowest *nonzero* curve point);
      never-drop-a-peptide noted rows. Three things surfaced by the SIS work that
      Methods should also state (details in [`docs/diagnostics.md`](diagnostics.md)):
      the trilinear fit requires **≥5 distinct curve points**; the LOQ thresholds the CV
      of bootstrap-resampled *means* rather than raw replicate CV (~1.7x apart at n=3);
      and LOQ is CV-gated while ULOQ is purely geometric.
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
- [x] Push the figure-port commit.
- [ ] Add a `LICENSE`. (The pinned tool is Apache-2.0, if we want to match.)
- [ ] Point `config.yaml → raw:` at wherever the MS data lives for anyone reproducing
      the raw-data figures (`run_all.py --with-raw`). The SIS entries currently point at
      a local `Downloads` path — move to the Drive root once the data lands there.
- [ ] `docs/figure_manifest.md` lists a `figures/fig03_loq_model.py` that doesn't exist;
      the README and `.gitignore` reference a `build_fom.py` that doesn't either. Either
      port them or drop the references (`scikit-learn` in `requirements.txt` is only
      there for the missing model script).

## Done
- [x] Fig 1B (distributions + examples), Fig 2A/2B/2C LOD/LOQ/ULOQ triptychs, Fig 3,
      before/after + tier supplements, ULOQ linear examples.
- [x] Tool changes (trilinear/ULOQ, mnp=0 rescue, noted rows) committed + pinned.
- [x] Repo scaffolded, smoke-tested, initial commit pushed.
