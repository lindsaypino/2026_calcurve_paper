"""Reshape a Skyline "Peptide Total Area Fragment" export from the SIS (AQUA)
peptide-response curve into the three inputs calculate-loq.py needs:

  1. a Skyline-shaped CSV with the exact columns the tool sniffs for
     (`File Name`, `Peptide Sequence`, `Total Area Fragment`),
  2. a filename -> concentration map, and
  3. a per-peptide `--multiplier_file`.

Why the extra prep. This export differs from every other dataset in the paper:

  * **No sequence column.** Skyline exported `Protein Name, File Name, Peptide Note`
    and one area column per label channel. Peptide identity is positional — the
    export is peptide-major, one contiguous block of `n_runs` rows per peptide, in
    the same order the peptides appear in the .sky document. We recover sequences
    from the .sky and *verify* the alignment against `Peptide Note` before trusting it.
  * **Two label channels.** `light` is the endogenous analyte and is entirely #N/A
    for the SIS peptides; the curve lives in `heavy` (the spiked standard). The tool
    wants a single `Total Area Fragment`, so a channel must be chosen (`--channel`).
  * **Per-peptide spike level.** The peptides are split into groups spiked at
    different absolute levels, recorded in the .sky as `concentration_multiplier`
    (1 / 10 / 100 / 1000 / 10000). A filename->concentration map cannot express a
    per-peptide x-axis, so the multipliers go in a separate `--multiplier_file`;
    the tool multiplies each peptide's curvepoints by its own multiplier.
  * **Concentrations live in the .sky**, not the filenames (which only carry S1..S5),
    as per-replicate `analyte_concentration` attributes.

iRT peptides are dropped by default: they are light-channel only, carry no
multiplier, and are retention-time standards rather than curve analytes.

CLI:  python -m src.prep_skyline_sis <export.csv> <document.sky> <out_dir>
API:  reshape(csv, sky) -> tidy (peptide, curvepoint, area, ...) DataFrame
"""
import os
import sys
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

# Skyline column names the tool's format sniffer requires, verbatim.
TOOL_COLS = ["File Name", "Peptide Sequence", "Total Area Fragment"]
IRT_PROTEIN = "irt"


def _norm_note(s):
    """Collapse Skyline note whitespace so CSV and XML copies compare equal."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    return " ".join(str(s).split()).replace("nan", "").strip()


def read_document(sky_path):
    """Peptides in .sky document order -> DataFrame(protein, peptide, multiplier, note).

    Order matters: it is the key that recovers peptide identity for the export.
    """
    root = ET.parse(sky_path).getroot()
    rows = []
    for plist in root.findall("peptide_list"):
        for pep in plist.findall("peptide"):
            mult = pep.get("concentration_multiplier")
            rows.append({
                "protein": plist.get("label_name"),
                "peptide": pep.get("modified_sequence"),
                "multiplier": float(mult) if mult else np.nan,
                "note": _norm_note(pep.findtext("note")),
            })
    return pd.DataFrame(rows)


def read_concentrations(sky_path):
    """Replicate name -> analyte_concentration, for replicates that declare one."""
    root = ET.parse(sky_path).getroot()
    return {r.get("name"): float(r.get("analyte_concentration"))
            for r in root.iter("replicate") if r.get("analyte_concentration")}


def _map_filenames(filenames, concs):
    """Match each export File Name to its .sky replicate by name suffix.

    Skyline stores the replicate as the run name without the acquisition prefix
    (e.g. `1_20ng_..._662` for `230707_AQUA_peptide_HEK_S1_20ng_..._662.d`), so the
    replicate name is a suffix of the stem. Suffix matching avoids re-deriving the
    dilution from the filename (which is also where the `ulltra` typo lives).
    """
    out = {}
    for fn in filenames:
        stem = fn[:-2] if fn.endswith(".d") else fn
        hits = [r for r in concs if stem.endswith(r)]
        if len(hits) != 1:
            raise ValueError(
                f"{'no' if not hits else 'ambiguous'} .sky replicate for run {fn!r}"
                + (f" (matched {hits})" if hits else ""))
        out[fn] = concs[hits[0]]
    return out


def reshape(csv_path, sky_path, channel="heavy", drop_irt=True):
    """Tidy the Skyline export into one row per (peptide, run) with concentrations.

    Peptide identity comes from block position; the mapping is validated against
    `Peptide Note` and raises if the export and document disagree.
    """
    d = pd.read_csv(csv_path)
    doc = read_document(sky_path)

    n_runs = d["File Name"].nunique()
    if len(d) != n_runs * len(doc):
        raise ValueError(
            f"export has {len(d)} rows; expected {len(doc)} peptides x {n_runs} runs "
            f"= {len(doc) * n_runs}. The positional identity assumption does not hold.")
    d["_block"] = d.index // n_runs

    # every block must cover the full run set exactly once, or the blocks are not runs-within-peptide
    per_block = d.groupby("_block")["File Name"].agg(["size", "nunique"])
    if not ((per_block["size"] == n_runs) & (per_block["nunique"] == n_runs)).all():
        raise ValueError("export blocks do not each contain the full run set once; "
                         "re-export sorted by peptide then replicate.")

    # validate positional identity against the notes Skyline did export
    csv_note = d.groupby("_block")["Peptide Note"].first().map(_norm_note).to_numpy()
    csv_prot = d.groupby("_block")["Protein Name"].first().to_numpy()
    bad = [(i, a, b) for i, (a, b) in enumerate(zip(csv_note, doc["note"].to_numpy())) if a != b]
    if bad or not (csv_prot == doc["protein"].to_numpy()).all():
        raise ValueError(
            "export block order does not match .sky document order — cannot recover "
            f"peptide identity positionally. First mismatches: {bad[:3]}")

    area_col = f"{channel} Total Area Fragment"
    if area_col not in d.columns:
        raise ValueError(f"no {area_col!r} column; available: {list(d.columns)}")

    conc_by_file = _map_filenames(d["File Name"].unique(), read_concentrations(sky_path))

    out = pd.DataFrame({
        "peptide": doc["peptide"].to_numpy()[d["_block"].to_numpy()],
        "protein": doc["protein"].to_numpy()[d["_block"].to_numpy()],
        "multiplier": doc["multiplier"].to_numpy()[d["_block"].to_numpy()],
        "note": doc["note"].to_numpy()[d["_block"].to_numpy()],
        "filename": d["File Name"].to_numpy(),
        "area": pd.to_numeric(d[area_col], errors="coerce").to_numpy(),
    })
    out["curvepoint"] = out["filename"].map(conc_by_file)
    out["eff_conc"] = out["curvepoint"] * out["multiplier"]
    out["group"] = out["note"].str.extract(r"(group_[A-E])")[0]

    if drop_irt:
        out = out[out["protein"] != IRT_PROTEIN].reset_index(drop=True)
    return out


def to_skyline(csv_path, sky_path, out_dir, channel="heavy", drop_irt=True, prefix="sis_ultra"):
    """Write the tool-ready export CSV, concentration map, and multiplier file.

    Returns (export_path, map_path, multiplier_path).
    """
    tidy = reshape(csv_path, sky_path, channel=channel, drop_irt=drop_irt)
    os.makedirs(out_dir, exist_ok=True)

    export = os.path.join(out_dir, f"{prefix}_skyline.csv")
    cmap = os.path.join(out_dir, f"{prefix}_map.csv")
    cmult = os.path.join(out_dir, f"{prefix}_multipliers.csv")

    # 1. the export, renamed to exactly what the tool's sniffer looks for. Areas are
    #    left as-is: the tool fills missing with 0 itself.
    tidy.rename(columns={"filename": "File Name", "peptide": "Peptide Sequence",
                         "area": "Total Area Fragment"})[TOOL_COLS].to_csv(export, index=False)

    # 2. filename -> concentration (one row per run)
    (tidy[["filename", "curvepoint"]].drop_duplicates()
        .rename(columns={"curvepoint": "concentration"})
        .sort_values(["concentration", "filename"])
        .to_csv(cmap, index=False))

    # 3. peptide -> multiplier. Peptides absent here are dropped by the tool's inner
    #    join, so a missing multiplier silently loses a peptide — fail loudly instead.
    mult = tidy[["peptide", "multiplier"]].drop_duplicates()
    if mult["multiplier"].isna().any():
        missing = mult.loc[mult["multiplier"].isna(), "peptide"].tolist()
        raise ValueError(f"peptides with no concentration_multiplier in the .sky: {missing}")
    mult.to_csv(cmult, index=False)

    return export, cmap, cmult


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    paths = to_skyline(sys.argv[1], sys.argv[2], sys.argv[3])
    for p in paths:
        print("wrote", p)
