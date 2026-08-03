"""Fig 3 - general trends: LOQ vs RT, m/z, and undiluted-A abundance for the Ultra
DIA-NN calcurve. Plain seaborn whitegrid / black-scatter style with a Spearman rho
(deliberately NOT the vizta theme). LOQ from the asms/ultra default-mnp2 FOM (stored
as legacy_mnp2/bruker_ultra.csv). NEEDS RAW DATA: the ultra diann_report (config raw).
"""
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pyteomics import mass

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.fom_io import load_config, load_fom, raw, repo_path

cfg = load_config()
fom_df = load_fom("legacy", "bruker_ultra", cfg)

# read only the columns needed (keeps the ~560MB report manageable)
df = pd.read_csv(raw("ultra_diann_report", cfg), sep="\t",
                 usecols=["Precursor.Id", "Modified.Sequence", "Precursor.Charge",
                          "RT", "File.Name", "Precursor.Quantity"])
df["peptide"] = df["Precursor.Id"]

# --- per-peptide RT ---
rt_df = df.groupby("peptide", as_index=False)["RT"].median()

# --- per-peptide m/z (pyteomics) ---
unimod_mass = {"4": 57.021464, "35": 15.994915}
PROTON = mass.nist_mass["H+"][0][0]


def mod_mass(seq):
    ids = re.findall(r"\(UniMod:(\d+)\)", seq)
    if not ids:
        return 0.0
    return np.nan if any(mid not in unimod_mass for mid in ids) else sum(unimod_mass[m] for m in ids)


def calc_mz(seq, z):
    if pd.isna(seq) or pd.isna(z):
        return np.nan
    delta = mod_mass(seq)
    if np.isnan(delta):
        return np.nan
    return (mass.calculate_mass(sequence=re.sub(r"\(UniMod:\d+\)", "", seq)) + delta + int(z) * PROTON) / int(z)


uniq = df.drop_duplicates("peptide")[["peptide", "Modified.Sequence", "Precursor.Charge"]].copy()
uniq["mz"] = uniq.apply(lambda r: calc_mz(r["Modified.Sequence"], r["Precursor.Charge"]), axis=1)
mz_df = uniq[["peptide", "mz"]]

# --- per-peptide undiluted (curve A) abundance ---
undiluted = df[df["File.Name"].str.contains("_A_", na=False)]
avg_signal = (undiluted.groupby("peptide", as_index=False)["Precursor.Quantity"].mean()
              .rename(columns={"Precursor.Quantity": "avg_precursor_qty_A"}))

# --- assemble panels in the original plain style ---
sns.set_style("whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(16, 4))


def panel(ax, merged, xcol, xlabel, logx=False):
    merged = merged[np.isfinite(merged["LOQ"])]
    ax.scatter(merged[xcol], merged["LOQ"], s=12, alpha=0.4, color="black")
    ax.set_xlabel(xlabel); ax.set_ylabel("LOQ")
    if logx:
        ax.set_xscale("log")
    rho = merged[xcol].corr(merged["LOQ"], method="spearman")
    ax.text(0.02, 0.98, f"Spearman rho = {rho:.2f}", transform=ax.transAxes, ha="left", va="top")


panel(axes[0], pd.merge(fom_df, rt_df, on="peptide", how="inner"), "RT", "RT (min)")
panel(axes[1], pd.merge(fom_df, mz_df, on="peptide", how="inner"), "mz", "m/z (pyteomics)")
panel(axes[2], pd.merge(fom_df, avg_signal, on="peptide", how="inner"),
      "avg_precursor_qty_A", "avg abundance, undiluted curve A", logx=True)

plt.tight_layout()
out = repo_path(cfg["output"], "FIG3_trends.png")
plt.savefig(out, dpi=300, bbox_inches="tight")
print("wrote", out)
