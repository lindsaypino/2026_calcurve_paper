# Note: the fit weights, and the arbitrary constant in them

Working note for the methods discussion. Recorded 2026-08-11. Companion to
[`loq_stability_note.md`](loq_stability_note.md).

**Not measured on the paper's datasets.** Numbers below come from
`data/one_protein.csv`, the tool's 27-peptide sample dataset. They characterize the
tool's sensitivity, not our results.

## What the tool does

The segmented fit is weighted, in `fit_by_lmfit_yang`:

```python
weights = np.minimum(1.0 / (np.sqrt(x) + np.finfo(float).eps), 1000)
```

Weighting by `1/sqrt(x)` is an ordinary proportional-error model — it assumes
variance grows with the quantity injected, which is the right general shape for MS
response. That part is defensible and worth stating plainly in Methods.

The `1000` is not derived from anything. At a 0-concentration point `1/sqrt(x)`
diverges, so the cap is what keeps the blank finite, and its value sets how much the
blank replicates dominate the fit:

| curve point | weight | squared leverage vs. x=1 |
|---|---|---|
| 0 | 1000 | 1,000,000× |
| 0.001 | 31.6 | 1,000× |
| 0.005 | 14.1 | 200× |
| 0.05 | 4.5 | 20× |
| 0.5 | 1.4 | 2× |
| 1.0 | 1.0 | 1× |

Least squares minimizes weighted *squared* residuals, so a blank replicate carries
roughly a million times the leverage of a top-of-curve replicate. In effect the noise
plateau is pinned almost exactly to the mean of the blanks, and the linear segment is
fit around that. For estimating an LOD that is arguably the behavior you want — but
it is a consequence of an unexamined constant, not a stated modeling choice.

## How much does the constant actually matter?

Re-fitting all 27 peptides with the cap changed:

| cap | LODs that move | median change | max change |
|---|---|---|---|
| 1000 → 100 | 26 of 26 finite | 0.49% | 2.65% |
| 1000 → 100,000 | 24 of 26 finite | 0.01% | 16.4% |

The count of peptides with a finite LOD is unchanged in both directions. So the
results are not balanced on this constant: a 10× or 100× change moves the typical
LOD by well under 1%. One peptide moves 16% at cap=100,000, so it is not entirely
inert either.

**Suggested framing:** state the `1/sqrt(x)` weighting as the modeling choice it is,
note that zero-concentration points are capped to keep the weight finite, and say the
cap's value does not materially affect the reported limits. If a reviewer asks, the
sensitivity numbers above are the answer.

## Interaction with the densification change

Tool commit `3d18164` (issue #15) fills in zero areas for runs where a peptide was
not identified, which for DIA-NN inputs adds many rows at *low* concentration — that
is, at high weight. This is the intended fix (those runs really did measure nothing,
and the noise plateau depends on them), but it does mean DIA-NN curves now lean on
the high-weight end of this scheme harder than they did when those rows were simply
absent. The sensitivity above suggests that is not a problem; worth a sentence if the
before/after comparison is discussed.

## Reproducing

Patch the cap and re-fit, from the tool's repo root:

```python
import io, importlib.util, sys
import numpy as np

src = io.open("bin/calculate-loq.py", encoding="utf-8").read()
anchor = "weights = np.minimum(1.0 / (np.sqrt(x) + np.finfo(float).eps), 1000)"

def variant(cap, name):
    path = f"{name}.py"
    io.open(path, "w", encoding="utf-8").write(
        src.replace(anchor, anchor.replace("), 1000)", f"), {cap})")))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod

base, low = variant(1000, "cap1000"), variant(100, "cap100")
df = base.read_input("data/one_protein.csv", "data/filename2samplegroup_map.csv")

for mod in (base, low):
    sub = df[df["peptide"] == "TLANTAVVIR"].sort_values(mod.SORT_KEYS, kind="mergesort")
    x = np.asarray(sub["curvepoint"], float)
    res, _ = mod.fit_by_lmfit_yang(x, np.asarray(sub["area"], float), "auto")
    mp = np.asarray([0.0, res.params["c"].value, res.params["a"].value,
                     res.params["b"].value])
    print(mod.__name__, mod.calculate_lod(mp, sub, 2.0, 2, 1, x, "auto")[0])
```
