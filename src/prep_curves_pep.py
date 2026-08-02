"""Reshape a hand-merged *_CURVES_pep.tsv (14 side-by-side per-point DIA-NN blocks,
A..N, each with its own metadata cols + replicate .d columns) into either a tidy
DIA-NN-report the LOD/LOQ tool can read, or a melted (peptide, curvepoint, area)
DataFrame for plotting. Used for the Bruker Ultra / Ultra II 60SPD separate-search data.

CLI:  python -m src.prep_curves_pep <curves_pep.tsv> <out_report.tsv> <out_map.csv>
API:  reshape(path) -> melted df ; reshape(path, human_only=False) keeps yeast.
"""
import re
import sys
import numpy as np
import pandas as pd

# standard A..N dilution scheme (fraction of undiluted), shared with one_protein/ultra
AN_CONC = {"A": 1, "B": 0.7, "C": 0.5, "D": 0.3, "E": 0.1, "F": 0.07, "G": 0.05,
           "H": 0.03, "I": 0.01, "J": 0.007, "K": 0.005, "L": 0.003, "M": 0.001, "N": 0}


def _blocks(df):
    """Walk columns, associating each sample (.d) column with its block's metadata columns."""
    cur_mod = cur_chg = cur_prot = None
    recs = []
    for c in df.columns:
        base = re.sub(r"_\d+$", "", c)
        if base == "Protein.Names":
            cur_prot = c
        elif base == "Modified.Sequence":
            cur_mod = c
        elif base == "Precursor.Charge":
            cur_chg = c
        elif c.strip().endswith(".d") or "_DIA" in c or re.search(r"_[A-N]_R\d", c):
            recs.append((c, cur_mod, cur_chg, cur_prot))
    return recs


def reshape(path, human_only=True):
    """Return a melted (peptide, curvepoint, area) DataFrame from a CURVES_pep.tsv."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    out = []
    for scol, mcol, ccol, pcol in _blocks(df):
        m = re.search(r"_([A-N])_R\d", scol) or re.search(r"_([A-N])[_.]", scol)
        if not m:
            continue
        conc = AN_CONC[m.group(1)]
        prot = df[pcol].astype(str) if pcol else pd.Series("", index=df.index)
        pep = (df[mcol].astype(str) + "_" + df[ccol].astype(str)).str.replace(":", "", regex=False)
        s = pd.DataFrame({"peptide": pep, "curvepoint": conc,
                          "area": pd.to_numeric(df[scol], errors="coerce"), "prot": prot})
        s = s[s.peptide.notna() & (s.peptide != "nan_nan")].dropna(subset=["area"])
        if human_only:
            s = s[~s.prot.str.contains("YEAST", na=False)]
        out.append(s.drop(columns=["prot"]))
    return pd.concat(out, ignore_index=True)


def to_report(path, out_report, out_map, human_only=True):
    """Write a DIA-NN-report TSV + filename->concentration map for calculate-loq.py."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    longs, conc_rows = [], []
    for scol, mcol, ccol, pcol in _blocks(df):
        m = re.search(r"_([A-N])_R\d", scol) or re.search(r"_([A-N])[_.]", scol)
        if not m:
            continue
        conc = AN_CONC[m.group(1)]
        conc_rows.append((scol, conc))
        prot = df[pcol].astype(str) if pcol else pd.Series("", index=df.index)
        sub = pd.DataFrame({
            "Precursor.Id": df[mcol].astype(str) + "_" + df[ccol].astype(str),
            "Stripped.Sequence": df[mcol].astype(str),
            "File.Name": scol,
            "Precursor.Quantity": pd.to_numeric(df[scol], errors="coerce"),
            "Protein.Names": prot,
        })
        sub = sub[sub["Precursor.Id"].notna() & (sub["Precursor.Id"] != "nan_nan")].dropna(subset=["Precursor.Quantity"])
        if human_only:
            sub = sub[~sub["Protein.Names"].str.contains("YEAST", na=False)]
        longs.append(sub.drop(columns=["Protein.Names"]))
    pd.concat(longs, ignore_index=True).to_csv(out_report, sep="\t", index=False)
    pd.DataFrame(conc_rows, columns=["filename", "concentration"]).to_csv(out_map, index=False)


if __name__ == "__main__":
    to_report(sys.argv[1], sys.argv[2], sys.argv[3])
    print("wrote", sys.argv[2], "and", sys.argv[3])
