"""Audit *why* peptides do or don't get a ULOQ, clause by clause.

`calculate-loq.py` accepts a trilinear (saturating) fit only when every clause of a
seven-part test passes (see `fit_by_lmfit_yang`, the `supported = (...)` block). When a
dataset yields no ULOQ the tool says nothing about which clause blocked it — this module
re-fits each peptide with the tool's own `_fit_one_model` and reports all seven
individually, so "no ULOQ" can be attributed to data shape vs. a threshold.

Dataset-agnostic: takes a melted frame with columns (peptide, curvepoint, area), i.e.
whatever the tool itself would see *after* any `--multiplier_file` is applied.

    from src.uloq_audit import audit, summarize
    a = audit(melted, min_saturation_points=2)
    print(summarize(a))

CLI:  python -m src.uloq_audit <melted.csv> [min_saturation_points]

Caveat: this reproduces the acceptance test, not the bootstrap, so it predicts which
peptides *can* get a ULOQ. Confirm against a real `--plot n` run before quoting numbers.
"""
import sys
import numpy as np
import pandas as pd

from src.fom_io import load_tool

# clause names mirror the order of the `supported = (...)` conjunction in the tool
CLAUSES = ["c_high_finite", "below_ymax", "at_least_half_ymax", "onset_positive",
           "enough_sat_points", "plateau_is_minority", "aic_prefers_trilinear"]


def audit(melted, min_saturation_points=2, tool=None):
    """Per-peptide clause-by-clause ULOQ acceptance. Returns a DataFrame."""
    cl = tool or load_tool()
    rows = []
    for pep, sub in melted.groupby("peptide"):
        sub = sub.sort_values("curvepoint")
        x = sub["curvepoint"].to_numpy(float)
        y = np.nan_to_num(sub["area"].to_numpy(float))  # the tool fills missing areas with 0
        w = np.minimum(1.0 / (np.sqrt(x) + np.finfo(float).eps), 1000)
        n_distinct = len(np.unique(x))
        rec = {"peptide": pep, "n_distinct": n_distinct}

        # the tool refuses to even attempt a third segment below 5 distinct points
        if n_distinct < 5:
            rows.append({**rec, "attempted": False, **{c: False for c in CLAUSES}})
            continue
        try:
            res_bi, _ = cl._fit_one_model(x, y, w, "bilinear")
            res_tri, _ = cl._fit_one_model(x, y, w, "trilinear")
        except Exception as e:
            rows.append({**rec, "attempted": False, "error": type(e).__name__,
                         **{c: False for c in CLAUSES}})
            continue

        a = res_tri.params["a"].value
        b = res_tri.params["b"].value
        c_high = res_tri.params["c_high"].value
        onset = (c_high - b) / a if a > 0 else np.inf
        n_sat = len(np.unique(x[x >= onset]))
        y_max = float(np.max(y)) if len(y) else np.nan
        aic_tri = getattr(res_tri, "aic", np.inf)
        aic_bi = getattr(res_bi, "aic", np.inf)

        rows.append({
            **rec, "attempted": True, "n_sat": n_sat, "y_max": y_max, "c_high": c_high,
            "ceiling_frac_of_ymax": c_high / y_max if y_max else np.nan, "onset": onset,
            "aic_tri": aic_tri, "aic_bi": aic_bi,
            "c_high_finite": bool(np.isfinite(c_high)),
            "below_ymax": bool(c_high < y_max),
            "at_least_half_ymax": bool(c_high >= 0.5 * y_max),
            "onset_positive": bool(onset > 0),
            "enough_sat_points": bool(n_sat >= min_saturation_points),
            "plateau_is_minority": bool(2 * n_sat <= n_distinct),
            "aic_prefers_trilinear": bool(aic_tri < aic_bi),
        })
    out = pd.DataFrame(rows)
    for c in CLAUSES:
        out[c] = out[c].fillna(False).astype(bool)
    out["accepted"] = out["attempted"].fillna(False).astype(bool) & out[CLAUSES].all(axis=1)
    return out


def summarize(a, min_saturation_points=2):
    """Human-readable rollup: how often each clause is the blocker."""
    n = len(a)
    lines = [f"{n} peptides; ULOQ accepted for {int(a['accepted'].sum())}",
             f"trilinear attempted for {int(a['attempted'].sum())} "
             f"(needs >=5 distinct curve points)", "",
             "clause failures (a peptide can fail several):"]
    for c in CLAUSES:
        nf = int((~a[c]).sum())
        lines.append(f"  {c:24s} {nf:4d}/{n}  ({100 * nf / n:5.1f}%)")
    # the clause that is the *sole* blocker is the actionable one
    sole = []
    for c in CLAUSES:
        others = [o for o in CLAUSES if o != c]
        m = (~a[c]) & a[others].all(axis=1) & a["attempted"]
        if int(m.sum()):
            sole.append((c, int(m.sum())))
    lines += ["", "sole blocker (everything else passed):"]
    lines += [f"  {c:24s} {k:4d}" for c, k in sorted(sole, key=lambda t: -t[1])] or ["  (none)"]
    if "n_sat" in a:
        lines += ["", f"distinct plateau points (n_sat): "
                      f"{dict(a['n_sat'].value_counts().sort_index())}",
                  f"  min_saturation_points is {min_saturation_points}; note the tool also "
                  f"caps n_sat via 2*n_sat <= n_distinct"]
    return "\n".join(lines)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    msp = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    melted = pd.read_csv(sys.argv[1])
    a = audit(melted, min_saturation_points=msp)
    print(summarize(a, msp))
    print("\nper-peptide:")
    cols = ["peptide", "n_distinct", "n_sat", "ceiling_frac_of_ymax"] + CLAUSES + ["accepted"]
    print(a[[c for c in cols if c in a]].to_string(index=False,
                                                   float_format=lambda v: f"{v:.3f}"))
