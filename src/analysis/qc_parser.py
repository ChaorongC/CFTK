"""
qc_parser.py — Parse existing process output files into per-sample QC metrics.

Key changes vs previous version:
  N7  parse_mbias_txt()      — Rewritten for MethylDackel --txt TSV format:
                                Strand/Read/Position/nMethylated/nUnmethylated
                                Saves structured TSV to mbias_data/ for report plots.
  N7b parse_chh_bedgraph()   — Legacy helper for historical CHH bedGraphs.

Run-order requirements:
  - markdup_dup_pct:          process step 3 with P1 patch (sambamba stderr capture)
  - median_frag_len:           cftk qc -s 2 (bamPEFragmentSize)
  All missing files → NaN (never an error).
"""

from __future__ import annotations
import os, re, zipfile
import numpy as np
import pandas as pd
from typing import Any

NaN = float("nan")


def _safe_float(s: Any) -> float:
    try:
        return float(str(s).replace(",", "").replace("%", "").strip())
    except (ValueError, TypeError):
        return NaN


def _read_lines(path: str) -> list[str]:
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return f.readlines()
    except OSError:
        return []


# ── N1: Trim Galore trimming report ──────────────────────────────────────────

def parse_trimming_report(path: str) -> dict:
    out = {"total_reads": NaN, "adapter_reads_pct": NaN,
           "quality_trimmed_bp_pct": NaN, "bases_written_pct": NaN}
    for line in _read_lines(path):
        line = line.strip()
        m = re.match(r"Total reads processed:\s+([\d,]+)", line)
        if m:
            out["total_reads"] = _safe_float(m.group(1))
        m = re.match(r"Reads with adapters:\s+[\d,]+\s+\(([\d.]+)%\)", line)
        if m:
            out["adapter_reads_pct"] = _safe_float(m.group(1))
        m = re.match(r"Quality-trimmed:\s+[\d,]+\s+bp\s+\(([\d.]+)%\)", line)
        if m:
            out["quality_trimmed_bp_pct"] = _safe_float(m.group(1))
        m = re.match(r"Total written \(filtered\):\s+[\d,]+\s+bp\s+\(([\d.]+)%\)", line)
        if m:
            out["bases_written_pct"] = _safe_float(m.group(1))
    return out


# ── N2: FastQC zip ────────────────────────────────────────────────────────────

_FASTQC_MODULES = [
    ("Per base sequence quality",    "fqc_per_base_quality"),
    ("Per sequence quality scores",  "fqc_per_seq_quality"),
    ("Per base sequence content",    "fqc_per_base_content"),
    ("Per sequence GC content",      "fqc_per_seq_gc"),
    ("Per base N content",           "fqc_per_base_n"),
    ("Sequence Length Distribution", "fqc_seq_length_dist"),
    ("Sequence Duplication Levels",  "fqc_seq_duplication"),
    ("Overrepresented sequences",    "fqc_overrepresented"),
    ("Adapter Content",              "fqc_adapter_content"),
    ("Per tile sequence quality",    "fqc_per_tile_quality"),
]
_STATUS_MAP = {"pass": 0, "warn": 1, "fail": 2}


def parse_fastqc_zip(path: str) -> dict:
    out = {"fqc_total_seqs": NaN, "fqc_pct_gc": NaN, "fqc_seq_length": NaN,
           "fqc_pct_duplicates": NaN, "fqc_pct_modules_failed": NaN}
    for _, key in _FASTQC_MODULES:
        out[key] = NaN
    if not path or not os.path.exists(path):
        return out
    try:
        with zipfile.ZipFile(path, "r") as zf:
            data_file = next((n for n in zf.namelist() if n.endswith("fastqc_data.txt")), None)
            if not data_file:
                return out
            content = zf.read(data_file).decode("utf-8", errors="replace")
    except Exception:
        return out

    statuses: list[int] = []
    for line in content.splitlines():
        if line.startswith(">>") and not line.startswith(">>END"):
            for full_name, key in _FASTQC_MODULES:
                if line.startswith(f">>{full_name}"):
                    parts = line.split("\t")
                    s = _STATUS_MAP.get(parts[-1].strip().lower() if len(parts) > 1 else "fail", 2)
                    out[key] = s
                    statuses.append(s)
                    break
        if "\t" in line and not line.startswith(">>") and not line.startswith("#"):
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) >= 2:
                k, v = parts[0], parts[1]
                if k == "Total Sequences":
                    out["fqc_total_seqs"] = _safe_float(v)
                elif k == "%GC":
                    out["fqc_pct_gc"] = _safe_float(v)
                elif k == "Sequence length":
                    m = re.search(r"(\d+)$", v)
                    out["fqc_seq_length"] = _safe_float(m.group(1)) if m else NaN
                elif k == "Total Deduplicated Percentage":
                    d = _safe_float(v)
                    out["fqc_pct_duplicates"] = round(100 - d, 2) if not np.isnan(d) else NaN
    if statuses:
        out["fqc_pct_modules_failed"] = round(100 * statuses.count(2) / len(statuses), 1)
    return out


# ── N3: samtools flagstat ─────────────────────────────────────────────────────

def parse_flagstat(path: str) -> dict:
    out = {"flagstat_total_reads": NaN, "flagstat_mapped_pct": NaN,
           "flagstat_properly_paired_pct": NaN, "flagstat_mapped_reads": NaN,
           "flagstat_properly_paired_reads": NaN}
    for line in _read_lines(path):
        m = re.match(r"(\d+) \+ \d+ in total", line)
        if m:
            out["flagstat_total_reads"] = _safe_float(m.group(1))
        # "N + 0 mapped (XX.XX%)" — capture both absolute count and pct
        m = re.search(r"(\d+) \+ \d+ mapped \((\d+\.?\d*)%", line)
        if m:
            out["flagstat_mapped_reads"] = _safe_float(m.group(1))
            out["flagstat_mapped_pct"]   = _safe_float(m.group(2))
        m = re.search(r"(\d+) \+ \d+ properly paired \((\d+\.?\d*)%", line)
        if m:
            out["flagstat_properly_paired_reads"] = _safe_float(m.group(1))
            out["flagstat_properly_paired_pct"]   = _safe_float(m.group(2))
    return out


# ── N4: samtools stats ────────────────────────────────────────────────────────

def parse_samtools_stats(path: str) -> dict:
    out = {"stats_error_rate": NaN, "stats_insert_size_mean": NaN,
           "stats_insert_size_std": NaN, "stats_raw_total_seqs": NaN}
    for line in _read_lines(path):
        if not line.startswith("SN\t"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        key = parts[1].rstrip(":").strip()
        val = _safe_float(parts[2])
        if key == "error rate":
            out["stats_error_rate"] = round(val * 100, 4) if not np.isnan(val) else NaN
        elif key == "insert size average":
            out["stats_insert_size_mean"] = val
        elif key == "insert size standard deviation":
            out["stats_insert_size_std"] = val
        elif key == "raw total sequences":
            out["stats_raw_total_seqs"] = val
    return out


# ── N5: sambamba markdup ──────────────────────────────────────────────────────

def parse_markdup_metrics(path: str) -> dict:
    """Parse sambamba markdup metrics file.

    Sambamba format (no inline percentage):
      sorted M end pairs
      found N duplicates
    → dup_pct = N / (M * 2) * 100   (M end pairs = M*2 total reads)

    Also handles formats with inline percentage for other tools.
    """
    out = {"markdup_dup_pct": NaN, "markdup_dup_reads": NaN}
    lines = _read_lines(path)
    # Collect all numbers first for sambamba format
    dup_reads = NaN
    total_pairs = NaN
    for line in lines:
        # Format 1: inline percentage (non-sambamba tools)
        m = re.search(
            r"found\s+([\d,]+)\s+duplicates\s+in\s+([\d,]+)\s+reads\s+\(([\d.]+)%\)", line)
        if m:
            out["markdup_dup_reads"] = _safe_float(m.group(1))
            out["markdup_dup_pct"]   = _safe_float(m.group(3))
            return out
        m = re.search(r"([\d,]+)\s*/\s*([\d,]+)\s*=\s*([\d.]+)%\s+duplicates", line)
        if m:
            out["markdup_dup_reads"] = _safe_float(m.group(1))
            out["markdup_dup_pct"]   = _safe_float(m.group(3))
            return out
        # Format 2: sambamba — "found N duplicates"
        m = re.search(r"found\s+([\d,]+)\s+duplicates", line)
        if m:
            dup_reads = _safe_float(m.group(1))
        # Format 2: sambamba — "sorted M end pairs"
        m = re.search(r"sorted\s+([\d,]+)\s+end\s+pairs", line)
        if m:
            total_pairs = _safe_float(m.group(1))
    # Compute from sambamba numbers if found
    if not np.isnan(dup_reads) and not np.isnan(total_pairs) and total_pairs > 0:
        total_reads = total_pairs * 2
        out["markdup_dup_reads"] = dup_reads
        out["markdup_dup_pct"]   = round(100.0 * dup_reads / total_reads, 3)
    return out


# ── N6: CpG bedGraph ──────────────────────────────────────────────────────────

def parse_cpg_bedgraph(path: str, min_depth: int = 1) -> dict:
    out = {"cpg_global_meth_pct": NaN, "cpg_covered_sites": NaN, "cpg_mean_depth": NaN}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(
            path, sep="\t", comment="#", header=None,
            names=["chr", "start", "end", "meth_pct", "count_m", "count_u"],
            dtype={"meth_pct": float, "count_m": float, "count_u": float},
            on_bad_lines="skip",
        )
        df = df[~df["chr"].astype(str).str.startswith("track")]
        df["depth"] = df["count_m"] + df["count_u"]
        df = df[df["depth"] >= min_depth]
        if df.empty:
            return out
        total_m = df["count_m"].sum()
        total_u = df["count_u"].sum()
        if total_m + total_u > 0:
            out["cpg_global_meth_pct"] = round(100.0 * total_m / (total_m + total_u), 2)
        out["cpg_covered_sites"] = int(len(df))
        out["cpg_mean_depth"]    = round(float(df["depth"].mean()), 2)
    except Exception:
        pass
    return out


# ── N7: MethylDackel mbias --txt TSV output ───────────────────────────────────
#
# Process P3a changes the mbias command to:
#   MethylDackel mbias --txt ... > {prefix}_mbias.txt 2> {prefix}_mbias_OT_OB.temp
#
# _mbias.txt now contains pure stdout from --txt flag:
#   Strand  Read  Position  nMethylated  nUnmethylated
#   OT      1     1         182437       180917
#   OT      2     1         163408       171646
#   OB      1     1         ...
#   OB      2     1         ...
#   (CpG context only)
#
# %Methylation = nMethylated / (nMethylated + nUnmethylated) * 100
#
# The file is also saved as structured TSV in mbias_data/ for the interactive
# M-bias report plot.  MultiQC can also read _mbias.txt directly.
#
# NOTE: --txt outputs CpG context only → bisulfite_conversion_rate cannot be
# computed from this file. Use parse_chh_bedgraph() (N7b) instead.

def parse_mbias_txt(path: str, save_tsv_dir: str | None = None) -> dict:
    """
    Parse MethylDackel *_mbias.txt produced with --txt flag.
    Returns CpG methylation stats and saves per-position TSV for report plots.
    bisulfite_conversion_rate is NOT computed here (use parse_chh_bedgraph).
    """
    out = {"mbias_cpg_meth_pct": NaN}

    lines = _read_lines(path)
    if not lines:
        return out

    # Detect format: --txt output starts with header "Strand\tRead\tPosition\t..."
    # Legacy format (pre-P3a) starts with "Suggested inclusion options:" or context headers
    has_txt_header = any(
        line.strip().startswith("Strand\t") or line.strip().startswith("Strand ")
        for line in lines[:5]
    )

    if not has_txt_header:
        # Legacy format or empty — cannot extract per-position data
        # Try to extract CpG mean from "Avg. CpG methylation:" if present
        for line in lines:
            m = re.search(r"Avg\.?\s*CpG\s+methylation[:\s]+([\d.]+)%?", line, re.IGNORECASE)
            if m:
                out["mbias_cpg_meth_pct"] = _safe_float(m.group(1))
        return out

    # ── Parse --txt TSV format ────────────────────────────────────────────────
    records = []
    in_data  = False
    for line in lines:
        line = line.rstrip()
        # Header line
        if re.match(r"Strand\s+Read\s+Position", line, re.IGNORECASE):
            in_data = True
            continue
        if not in_data:
            continue
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            strand       = parts[0]   # OT / OB / CTOT / CTOB
            read         = int(parts[1])
            pos          = int(parts[2])
            n_meth       = int(parts[3])
            n_unmeth     = int(parts[4])
            total        = n_meth + n_unmeth
            pct          = round(100.0 * n_meth / total, 3) if total > 0 else NaN
            records.append({
                "strand":   strand,
                "read":     read,
                "pos":      pos,
                "n_meth":   n_meth,
                "n_unmeth": n_unmeth,
                "pct_meth": pct,
            })
        except (ValueError, IndexError):
            continue

    if not records:
        return out

    df = pd.DataFrame(records)

    # Compute overall CpG methylation% (weighted across all positions/strands/reads)
    total_m = df["n_meth"].sum()
    total_u = df["n_unmeth"].sum()
    if total_m + total_u > 0:
        out["mbias_cpg_meth_pct"] = round(100.0 * total_m / (total_m + total_u), 3)

    # Save per-position TSV for interactive M-bias report plot
    if save_tsv_dir:
        os.makedirs(save_tsv_dir, exist_ok=True)
        sample  = os.path.basename(path).replace("_mbias.txt", "")
        out_tsv = os.path.join(save_tsv_dir, f"{sample}_mbias.tsv")
        df["sample"] = sample
        df[["sample", "strand", "read", "pos", "pct_meth"]].to_csv(
            out_tsv, sep="\t", index=False
        )

    return out


# ── N7b: CHH bedGraph → bisulfite_conversion_rate ────────────────────────────
#
# Process P3b adds:
#   MethylDackel extract --CHH --noCpG --minDepth 1 -o {prefix}_chh {ref} {bam}
# Output: {prefix}_chh_CHH.bedGraph
# Format: chr  start  end  pct_meth  nMethylated  nUnmethylated
#
# bisulfite_conversion_rate = 100 - weighted_mean(CHH pct_meth)
# Expected: CHH meth% ~0.3-0.8% for good conversion (>99% conversion rate)

def parse_chh_bedgraph(path: str) -> dict:
    """
    Parse {prefix}_chh_CHH.bedGraph for bisulfite conversion rate.
    Returns NaN if file doesn't exist (P3b patch not yet applied or step 4 not re-run).
    """
    out = {"bisulfite_conversion_rate": NaN, "chh_meth_pct": NaN}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(
            path, sep="\t", comment="#", header=None,
            names=["chr", "start", "end", "meth_pct", "count_m", "count_u"],
            dtype={"meth_pct": float, "count_m": float, "count_u": float},
            on_bad_lines="skip",
        )
        df = df[~df["chr"].astype(str).str.startswith("track")]
        df["depth"] = df["count_m"] + df["count_u"]
        df = df[df["depth"] > 0]
        if df.empty:
            return out
        # Weighted mean CHH methylation%
        total_m = df["count_m"].sum()
        total_u = df["count_u"].sum()
        if total_m + total_u > 0:
            chh_pct = round(100.0 * total_m / (total_m + total_u), 3)
            out["chh_meth_pct"]              = chh_pct
            out["bisulfite_conversion_rate"] = round(100.0 - chh_pct, 3)
    except Exception:
        pass
    return out


# ── N8: fragment length CSV ───────────────────────────────────────────────────

def parse_fragment_csv(path: str) -> dict:
    out = {"median_frag_len": NaN, "mononuc_peak_bp": NaN,
           "short_frag_ratio": NaN, "mono_nuc_ratio": NaN, "di_nuc_ratio": NaN}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", skiprows=1, header=None)
        df.columns = ["size", "count"] + (["scaled"] if df.shape[1] > 2 else [])
        df = df[df["size"] < 1000].copy()
        df["size"]  = df["size"].astype(int)
        df["count"] = df["count"].astype(float)
        total = df["count"].sum()
        if total == 0:
            return out
        df["cumsum"] = df["count"].cumsum()
        out["median_frag_len"] = int(df[df["cumsum"] >= total / 2.0].iloc[0]["size"])
        df_filt = df[(df["size"] >= 100) & (df["size"] <= 400)]
        if not df_filt.empty:
            out["mononuc_peak_bp"] = int(df_filt.loc[df_filt["count"].idxmax(), "size"])
        out["short_frag_ratio"] = round(df[df["size"] < 120]["count"].sum() / total * 100, 2)
        out["mono_nuc_ratio"]   = round(
            df[(df["size"] >= 140) & (df["size"] <= 180)]["count"].sum() / total * 100, 2)
        out["di_nuc_ratio"]     = round(
            df[(df["size"] >= 300) & (df["size"] <= 360)]["count"].sum() / total * 100, 2)
    except Exception:
        pass
    return out


# ── Top-level collector ───────────────────────────────────────────────────────


def parse_beta_density_qc(matrix_path: str) -> pd.DataFrame:
    """
    Compute beta-density QC metrics from cpg_matrix.tsv.

    For each sample column, computes:
      left_peak  = max histogram density in beta <= 0.15
      right_peak = max histogram density in beta >= 0.85
      mid_density = median histogram density in 0.35 <= beta <= 0.65

      M_score            = min(left_peak, right_peak) / mid_density
      peak_balance       = min(left_peak, right_peak) / max(left_peak, right_peak)
      boundary_mid_ratio = frac(beta<=0.15 or beta>=0.85) / frac(0.35<=beta<=0.65)

    PASS: M_score >= 2, peak_balance >= 0.30, boundary_mid_ratio >= 2
    FAIL: otherwise (no WARN level)

    Returns a DataFrame indexed by sample name with three columns:
      beta_M_score, beta_peak_balance, beta_boundary_mid_ratio
    """
    if not matrix_path or not os.path.exists(matrix_path):
        print(f"[qc_parser] beta_density_qc: matrix not found at {matrix_path}")
        return pd.DataFrame()

    print(f"[qc_parser] beta_density_qc: reading {matrix_path} ...")
    try:
        df = pd.read_csv(matrix_path, sep="\t", index_col=0, low_memory=False)
        df = df.apply(pd.to_numeric, errors="coerce")
    except Exception as e:
        print(f"[qc_parser] beta_density_qc: failed to read matrix — {e}")
        return pd.DataFrame()

    BINS  = 200
    RANGE = (0.0, 1.0)

    records = {}
    for col in df.columns:
        vals = df[col].dropna().values
        if len(vals) < 100:
            records[col] = {
                "beta_M_score": NaN,
                "beta_peak_balance": NaN,
                "beta_boundary_mid_ratio": NaN,
            }
            continue

        # Auto-detect scale: 0-1 or 0-100 (percentage), normalise to 0-1
        vmax_val = float(vals.max())
        vals_norm = vals / 100.0 if vmax_val > 1.5 else vals.copy()
        vals_norm = np.clip(vals_norm, 0.0, 1.0)

        # Histogram density on normalised 0-1 values
        hist, edges = np.histogram(vals_norm, bins=BINS, range=RANGE, density=True)
        centers = (edges[:-1] + edges[1:]) / 2

        left_peak   = hist[centers <= 0.15].max()  if (centers <= 0.15).any() else NaN
        right_peak  = hist[centers >= 0.85].max()  if (centers >= 0.85).any() else NaN
        mid_mask    = (centers >= 0.35) & (centers <= 0.65)
        mid_density = float(np.median(hist[mid_mask])) if mid_mask.any() else NaN

        # M-score — if mid_density is 0 (perfect bimodal, no mid CpGs),
        # treat as very large M-score (PASS). Use epsilon floor to avoid /0.
        if not np.isnan(left_peak) and not np.isnan(right_peak) and not np.isnan(mid_density):
            denom = mid_density if mid_density > 1e-10 else 1e-10
            M_score = float(min(left_peak, right_peak)) / denom
        else:
            M_score = NaN

        # Peak balance
        if not np.isnan(left_peak) and not np.isnan(right_peak) and \
                max(left_peak, right_peak) > 0:
            peak_balance = float(min(left_peak, right_peak)) / float(max(left_peak, right_peak))
        else:
            peak_balance = NaN

        # Boundary/mid ratio (fraction-based, on normalised values)
        n = len(vals_norm)
        boundary_frac = float(((vals_norm <= 0.15) | (vals_norm >= 0.85)).sum()) / n
        mid_frac      = float(((vals_norm >= 0.35) & (vals_norm <= 0.65)).sum()) / n
        if mid_frac > 0:
            boundary_mid_ratio = boundary_frac / mid_frac
        else:
            boundary_mid_ratio = NaN

        records[col] = {
            "beta_M_score":            round(M_score, 4)            if not np.isnan(M_score) else NaN,
            "beta_peak_balance":       round(peak_balance, 4)       if not np.isnan(peak_balance) else NaN,
            "beta_boundary_mid_ratio": round(boundary_mid_ratio, 4) if not np.isnan(boundary_mid_ratio) else NaN,
        }

    result = pd.DataFrame.from_dict(records, orient="index")
    result.index.name = "sample"
    print(f"[qc_parser] beta_density_qc: computed for {len(result)} samples")
    return result


def collect_qc_metrics(
    samples: list[dict],
    paths: dict,
    out_tsv: str,
    cores: int = 1,
) -> pd.DataFrame:
    """
    Collect all QC metrics for all samples → qc_summary.tsv.

    Missing data handling (all return NaN, never error):
      - markdup_dup_pct:          P1 patch + step 3 re-run required
      - mbias_cpg_meth_pct:       P3a patch + step 4 re-run required
      - median_frag_len:          cftk qc -s 2 must run before cftk qc -s 0
    """
    from joblib import Parallel, delayed

    # mbias_data lives alongside the methylation outputs, not the QC output dir
    mbias_tsv_dir = os.path.join(paths.get("methylation", ""), "mbias_data")
    if not paths.get("methylation"):
        mbias_tsv_dir = os.path.join(os.path.dirname(os.path.abspath(out_tsv)), "mbias_data")

    def _one(sample: dict) -> dict:
        name  = sample["name"]
        group = sample.get("group", "")
        row: dict = {"sample": name, "group": group}

        # N1: trimming report
        row.update(parse_trimming_report(
            os.path.join(paths["trimming"], f"{name}_R1_trimming_report.txt")))

        # N2: FastQC trimmed R1
        row.update(parse_fastqc_zip(
            os.path.join(paths["trimming"], f"{name}_R1_fastqc.zip")))

        # N3: flagstat
        row.update(parse_flagstat(
            os.path.join(paths["alignment"], f"{name}.bam.flagstat")))

        # N4: samtools stats
        row.update(parse_samtools_stats(
            os.path.join(paths["alignment"], f"{name}.bam.stats")))

        # Recalculate TRUE mapping rate.
        # bwameth BAM contains only mapped reads, so both flagstat "mapped %" AND
        # samtools-stats "raw total sequences" come from the mapped-only BAM
        # → their ratio is always ~100%. The correct denominator is the number of
        # reads that ENTERED alignment, i.e. trimmed input reads.
        #
        # Trim Galore reports "Total reads processed" per FASTQ (the R1 file).
        # For paired-end data, total input reads = total_reads (R1) * 2.
        # mapped_pct = mapped_reads / (R1_reads * 2) * 100.
        mapped_rd = row.get("flagstat_mapped_reads", NaN)
        r1_reads  = row.get("total_reads", NaN)   # R1 reads (from trimming report)
        input_reads = r1_reads * 2 if not np.isnan(r1_reads) else NaN

        if not np.isnan(input_reads) and input_reads > 0 and not np.isnan(mapped_rd):
            pct = 100.0 * mapped_rd / input_reads
            row["flagstat_mapped_pct"] = round(min(pct, 100.0), 3)

        # N5: markdup (requires P1 patch; NaN if missing)
        row.update(parse_markdup_metrics(
            os.path.join(paths["markdup"], f"{name}.markdup_metrics.txt")))

        # N6: CpG bedGraph
        row.update(parse_cpg_bedgraph(
            os.path.join(paths["methylation"], f"{name}_CpG.bedGraph")))

        # N7: mbias --txt TSV (requires P3a patch; saves per-position TSV for plots)
        row.update(parse_mbias_txt(
            os.path.join(paths["methylation"], f"{name}_mbias.txt"),
            save_tsv_dir=mbias_tsv_dir,
        ))

        # N8: fragment CSV (NaN if cftk qc -s 2 not yet run)
        frag_prefix = os.path.join(paths["qc"], "2_fragment_length", "fragment_length")
        frag_csv    = next(
            (f"{frag_prefix}.{stem}.raw.csv"
             for stem in [name, f"{name}.markdup"]
             if os.path.exists(f"{frag_prefix}.{stem}.raw.csv")),
            "",
        )
        row.update(parse_fragment_csv(frag_csv))

        return row

    rows = Parallel(n_jobs=max(1, cores), backend="threading")(
        delayed(_one)(s) for s in samples
    )
    df   = pd.DataFrame(rows)
    cols = ["sample", "group"] + [c for c in df.columns if c not in ("sample", "group")]
    df   = df[cols]

    # Beta-density QC from cpg_matrix.tsv
    meth_dir    = paths.get("methylation", "")
    work_dir    = os.path.dirname(meth_dir) if meth_dir else ""
    matrix_path = os.path.join(work_dir, "5_merged_matrix", "cpg_matrix.tsv")
    beta_df = parse_beta_density_qc(matrix_path)
    if not beta_df.empty:
        df = df.merge(
            beta_df.reset_index(),  # index = sample name
            on="sample", how="left"
        )

    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[qc_parser] qc_summary.tsv → {out_tsv}  ({len(df)} samples)")
    return df
