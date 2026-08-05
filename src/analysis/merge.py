"""
merge.py — manual merge module for building feature matrices from arbitrary sources.

Supports merging per-sample files into a single features × samples matrix TSV.
Each modality defines its own groups, sample paths, and file type.

Supported types:
  bedgraph      : CpG methylation bedGraph (MethylDackel or user-provided)
                  → merged via bedtools unionbedg
  occupancy_tsv : bigWigAverageOverBed output (col0=region, col5=mean)
  wps_tsv       : wps.py output (chr/start/end/WPS/mean_WPS columns)
  delfi_tsv     : finaletoolkit delfi output (ratio_corrected column)

Output location per modality:
  cpg       → 1_process/5_merged_matrix/cpg_matrix.tsv
  occupancy → 4_fragmentomics/occupancy/occupancy_matrix.tsv
  wps       → 4_fragmentomics/wps/wps_matrix.tsv
  delfi     → 4_fragmentomics/delfi/delfi_matrix.tsv
  (any other modality name → {output_dir}/{modality}_matrix.tsv)
"""

import os
import shutil
import sys

import pandas as pd

from util import disp, recorded_run


# ── Type readers ──────────────────────────────────────────────────────────────

def _read_bedgraph(entries, out_path):
    """
    Merge CpG bedGraph files via bedtools unionbedg.
    Input files must be pre-sorted (user responsibility).
    entries: list of {"name": col_name, "path": filepath}
    """
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        sys.exit("[merge] ERROR: bedtools not found in PATH.")

    col_names = [e["name"] for e in entries]
    names_str = " ".join(col_names)
    stripped  = " ".join(
        f"<(grep -v '^track' {e['path']} | cut -f1-4)"
        for e in entries
    )
    tmp = out_path + ".tmp"
    cmd = (
        f"bash -c '"
        f"{bedtools} unionbedg -header -names {names_str} -filler NA "
        f"-i {stripped}"
        f" | cut -f1,3-"
        f" | sed s/end/pos/"
        f" > {tmp}'"
    )
    disp(f"[merge] bedtools unionbedg → {out_path}")
    ret = recorded_run(cmd, shell=True, label="merge bedGraph matrix")
    if ret.returncode != 0:
        sys.exit("[merge] ERROR: bedtools unionbedg failed.")

    matrix = pd.read_csv(tmp, sep="\t", index_col=0)
    os.rename(tmp, out_path)
    return matrix


def _read_occupancy_tsv(entries, out_path):
    """
    Merge bigWigAverageOverBed output files.
    No header; col0=region_name, col5=mean (covered bases only).
    """
    frames = {}
    for e in entries:
        df = pd.read_csv(
            e["path"], sep="\t", header=None,
            usecols=[0, 5], names=["region", e["name"]]
        )
        frames[e["name"]] = df.set_index("region")[e["name"]]
    matrix = pd.DataFrame(frames)
    matrix.to_csv(out_path, sep="\t")
    return matrix


def _read_wps_tsv(entries, out_path):
    """
    Merge wps.py output files.
    Columns: chr, start, end, WPS, mean_WPS.
    Index: chr:start-end.
    """
    frames = {}
    for e in entries:
        df = pd.read_csv(e["path"], sep="\t")
        df.index = (df["chr"] + ":" + df["start"].astype(str)
                    + "-" + df["end"].astype(str))
        frames[e["name"]] = df["mean_WPS"]
    matrix = pd.DataFrame(frames)
    matrix.to_csv(out_path, sep="\t")
    return matrix


def _read_delfi_tsv(entries, out_path):
    """
    Merge finaletoolkit delfi output files.
    Uses ratio_corrected column (falls back to ratio).
    Index: chr:start-end.
    """
    frames = {}
    for e in entries:
        with open(e["path"]) as f:
            first = f.readline()
        df = pd.read_csv(e["path"], sep="\t")
        if first.startswith("#"):
            df.columns = df.columns.str.replace("^#", "", regex=True).str.strip()
        ratio_col = "ratio_corrected" if "ratio_corrected" in df.columns else "ratio"
        df = df.dropna(subset=[ratio_col])
        df["contig"] = df["contig"].astype(str)
        df.index = (df["contig"] + ":" + df["start"].astype(str)
                    + "-" + df["end"].astype(str))
        frames[e["name"]] = df[ratio_col]
    matrix = pd.DataFrame(frames)
    matrix.to_csv(out_path, sep="\t")
    return matrix


# ── Type dispatch ─────────────────────────────────────────────────────────────

_TYPE_READERS = {
    "bedgraph":     _read_bedgraph,
    "occupancy_tsv": _read_occupancy_tsv,
    "wps_tsv":      _read_wps_tsv,
    "delfi_tsv":    _read_delfi_tsv,
}

SUPPORTED_TYPES = sorted(_TYPE_READERS.keys())


# ── Output path resolver ──────────────────────────────────────────────────────

def _out_path(modality, paths):
    """Canonical output path for each modality matrix."""
    mapping = {
        "cpg":       os.path.join(paths["cpg_matrix"], "cpg_matrix.tsv"),
        "occupancy": os.path.join(paths["occ_out"],    "occupancy_matrix.tsv"),
        "wps":       os.path.join(paths["wps_out"],    "wps_matrix.tsv"),
        "delfi":     os.path.join(paths["delfi_out"],  "delfi_matrix.tsv"),
    }
    if modality in mapping:
        return mapping[modality]
    # unknown modality → place under fragmentomics
    return os.path.join(paths["fragmentomics"], f"{modality}_matrix.tsv")


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_merge_cfg(modality, mod_cfg, cfg):
    """Validate a single modality block in the merge config."""
    errors = []
    valid_groups = list(cfg.get("samples", {}).keys())
    ga, gb = cfg["comparison"].split("_vs_", 1)

    groups = mod_cfg.get("groups", {})
    if set(groups.keys()) != {ga, gb}:
        errors.append(
            f"merge.{modality}.groups must contain exactly the two comparison "
            f"groups: {ga}, {gb}. Got: {sorted(groups.keys())}"
        )
        return errors  # no point checking further

    for grp, entries in groups.items():
        if not isinstance(entries, list) or len(entries) == 0:
            errors.append(f"merge.{modality}.groups.{grp} must be a non-empty list.")
            continue
        for e in entries:
            if not e.get("name"):
                errors.append(f"merge.{modality}.groups.{grp}: entry missing 'name'.")
            if e.get("type") not in _TYPE_READERS:
                errors.append(
                    f"merge.{modality}.groups.{grp}.{e.get('name','?')}: "
                    f"unsupported type '{e.get('type')}'. "
                    f"Supported: {SUPPORTED_TYPES}"
                )
            if not e.get("path"):
                errors.append(
                    f"merge.{modality}.groups.{grp}.{e.get('name','?')}: "
                    "missing 'path'."
                )
    return errors


# ── Main merge function ───────────────────────────────────────────────────────

def run_merge(modality, cfg, paths):
    """
    Merge all samples for one modality into a feature matrix TSV.
    Reads configuration from cfg["merge"][modality].
    Column names = sample["name"] (user-defined, must start with group name
    for diff/mesa group matching).
    Overwrites existing matrix if present.
    """
    mod_cfg = cfg.get("merge", {}).get(modality)
    if mod_cfg is None:
        sys.exit(
            f"[merge] ERROR: 'merge.{modality}' not found in cftk_init.json.\n"
            f"Add a 'merge' block after 'process' in your config."
        )

    errors = _validate_merge_cfg(modality, mod_cfg, cfg)
    if errors:
        disp(f"[merge] ERROR: invalid merge.{modality} config:")
        for e in errors:
            disp(f"  x {e}")
        sys.exit(1)

    ga, gb = cfg["comparison"].split("_vs_", 1)
    groups = mod_cfg["groups"]

    # flatten entries in group order (group_a first, group_b second)
    all_entries = []
    for grp in [ga, gb]:
        for e in groups[grp]:
            all_entries.append({
                "name":  e["name"],
                "type":  e["type"],
                "path":  e["path"],
                "group": grp,
            })

    # validate files exist
    missing = [e["path"] for e in all_entries if not os.path.exists(e["path"])]
    if missing:
        disp("[merge] ERROR: the following input files were not found:")
        for p in missing:
            disp(f"  {p}")
        sys.exit(1)

    # check all types are the same (mixing types in one modality not supported)
    types = {e["type"] for e in all_entries}
    if len(types) > 1:
        sys.exit(
            f"[merge] ERROR: all entries in merge.{modality} must have the same type. "
            f"Found: {types}"
        )
    file_type = types.pop()
    reader    = _TYPE_READERS[file_type]

    # resolve output path and create dir
    out_path = _out_path(modality, paths)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    disp(f"[merge] modality={modality}  type={file_type}  "
         f"samples={len(all_entries)}  → {out_path}")

    # run merge
    matrix = reader(all_entries, out_path)

    disp(f"[merge] {matrix.shape[0]} features × {matrix.shape[1]} samples "
         f"→ {out_path}")
    return out_path
