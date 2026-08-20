### What happens

In the per-peptide PNGs from `--plot y`, a trilinear fit with a finite ULOQ draws an
orange vertical line in the top (signal) panel, but the legend only ever lists `LOD` and
`LOQ`. The ULOQ line ends up unlabeled, so there's no way to tell from the figure alone
what the orange line is or what value it sits at.

### Why

`build_plots` draws the ULOQ marker in the **top** subplot:

```python
# bin/calculate-loq.py, ~line 645
if np.isfinite(ULOQ):
    plt.axvline(x=ULOQ, color='orange', label=('ULOQ = %.3e' % ULOQ))
```

but the legend is built *after* switching to the **bottom** (CV) subplot at
`plt.subplot(2, 1, 2)`:

```python
# ~line 705
legend = plt.legend(loc=8, bbox_to_anchor=(0, -.75, 1., .102), ncol=2)
```

`plt.legend()` collects handles from the current axes only. The bottom subplot re-draws
`LOD` and `LOQ` (~lines 680–686), which is why those two *do* show up — but ULOQ is not
re-drawn there, so it can never be picked up. The comment above the call still reads
"add legend with LOD and LOQ values", which suggests it just wasn't revisited when ULOQ
was added.

### Minimal reproduction

Independent of this codebase — shows the legend scope behaviour:

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(2, 1, 1)
plt.axvline(x=1, color="m", label="LOD = 1")
plt.axvline(x=2, color="c", label="LOQ = 2")
plt.axvline(x=3, color="orange", label="ULOQ = 3")   # only drawn here

plt.subplot(2, 1, 2)
plt.axvline(x=1, color="m", label="LOD = 1")          # re-drawn
plt.axvline(x=2, color="c", label="LOQ = 2")          # re-drawn
# ULOQ deliberately not re-drawn, matching build_plots

leg = plt.legend()
print([t.get_text() for t in leg.get_texts()])
# -> ['LOD = 1', 'LOQ = 2']
```

Observed on `ac0b951` against a 5-point Skyline calibration curve run with
`--model auto --plot y --min_noise_points 0 --min_saturation_points 1`; every peptide
with a finite ULOQ shows the unlabeled orange line.

### Suggested fix

Mirror the existing LOD/LOQ handling and re-draw ULOQ in the bottom subplot alongside
them, next to lines 680–686:

```python
if np.isfinite(ULOQ):
    plt.axvline(x=ULOQ, color='orange', label=('ULOQ = %.3e' % ULOQ))
```

That gets it into the legend and, as a side benefit, marks the ULOQ on the CV panel too —
useful for seeing whether CV actually degrades above the claimed upper limit. Bumping
`ncol=2` to `3` would keep the legend on one row.

Alternative, if you'd rather not duplicate the line: gather handles from both axes and
build a figure-level legend via `fig.legend(...)`.

### Impact

Cosmetic — no effect on any computed figure of merit; `figuresofmerit.csv` is correct.
Only matters where the generated plots are shown to someone, e.g. a supplemental figure.
