"""
dmr.py — DMR analysis pipeline.
Step 1: prepare metilene input (bedtools unionbedg, integrated from prepare_metilene_input.py)
Step 2: run metilene to call DMRs
Step 3: annotate DMRs via dmr_annotation.r (annotatr hg38)
"""

import os
import shutil
import sys

from util import disp, recorded_run


# ── Internal helpers ──────────────────────────────────────────────────────────

def _run(cmd, label=""):
    disp(f"[dmr] {label}")
    ret = recorded_run(cmd, shell=True, label=f"dmr {label}")
    if ret.returncode != 0:
        sys.exit(f"[dmr] ERROR: failed — {label}")


def _resolve_bedtools():
    """Find bedtools in PATH."""
    path = shutil.which("bedtools")
    if path is None:
        sys.exit("[dmr] ERROR: bedtools not found in PATH.")
    return path


def _validate_bedgraph(fp):
    if not os.path.exists(fp):
        sys.exit(f"[dmr] ERROR: bedGraph not found: {fp}")
    if not os.access(fp, os.R_OK):
        sys.exit(f"[dmr] ERROR: bedGraph not readable: {fp}")
    return os.path.abspath(fp)


# ── Step 1: prepare metilene input ───────────────────────────────────────────

def prepare_metilene_input(
    files_a: list,
    files_b: list,
    group_a: str,
    group_b: str,
    out_path: str,
):
    """
    Merge two groups of bedGraph files into a metilene-compatible input file
    using bedtools unionbedg.

    files_a / files_b : list of absolute paths to CpG.bedGraph files
    group_a / group_b : group name labels (used as column headers)
    out_path          : output file path
    """
    disp(f"[dmr] Preparing metilene input")
    disp(f"  group_a ({group_a}): {len(files_a)} sample(s)")
    disp(f"  group_b ({group_b}): {len(files_b)} sample(s)")

    # validate all input files
    files_a = [_validate_bedgraph(f) for f in files_a]
    files_b = [_validate_bedgraph(f) for f in files_b]

    bedtools = _resolve_bedtools()

    # column header names: group_a repeated, then group_b repeated
    names = " ".join([group_a] * len(files_a) + [group_b] * len(files_b))

    # strip track headers, keep 4 cols; unionbedg merges; cut chr+pos+values; rename end→pos
    stripped = " ".join(
        f"<(grep -v '^track' {f} | cut -f1-4)" for f in files_a + files_b
    )
    cmd = (
        f"bash -c '"
        f"{bedtools} unionbedg -header -names {names} -filler NA -i {stripped}"
        f" | cut -f1,3-"
        f" | sed s/end/pos/"
        f" > {out_path}'"
    )

    disp(f"[dmr] Writing metilene input → {out_path}")
    ret = recorded_run(cmd, shell=True, label="dmr annotation")
    if ret.returncode != 0:
        sys.exit("[dmr] ERROR: bedtools unionbedg pipeline failed.")

    disp(f"[dmr] Metilene input ready: {out_path}")
    return out_path


# ── Step 2 + 3: metilene + annotation ────────────────────────────────────────

def run_dmr(args):
    """
    Full DMR pipeline:
      1. Prepare metilene input from per-sample bedGraph files.
      2. Run metilene to call DMRs.
      3. Annotate DMRs via dmr_annotation.r.

    args fields (set by _cmd_dmr in cftk.py):
      output_dir       : output directory
      group_a          : group A name (label=0, from comparison)
      group_b          : group B name (label=1, from comparison)
      bedgraph_a       : list of bedGraph paths for group A
      bedgraph_b       : list of bedGraph paths for group B
      metilene_tool    : tool name (default: "metilene")
      threads          : metilene threads
      dmr_extra_args   : extra args passed to metilene
    """
    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)

    metilene_input = os.path.join(out_dir, "metilene_input.bedGraph")
    raw_bed        = os.path.join(out_dir, "dmr_raw.bed")
    ann_bed        = os.path.join(out_dir, "dmr_annotated.bed")

    # ── Step 1: prepare metilene input ───────────────────────────────────────
    prepare_metilene_input(
        files_a  = args.bedgraph_a,
        files_b  = args.bedgraph_b,
        group_a  = args.group_a,
        group_b  = args.group_b,
        out_path = metilene_input,
    )

    # ── Step 2: metilene ─────────────────────────────────────────────────────
    tool    = getattr(args, "metilene_tool", "metilene")
    threads = getattr(args, "threads", 20)
    extra   = getattr(args, "dmr_extra_args", "")

    if tool == "metilene":
        cmd = (
            f"metilene -a {args.group_a} -b {args.group_b} "
            f"-t {threads} {extra} {metilene_input} > {raw_bed}"
        )
    else:
        sys.exit(f"[dmr] unsupported tool: {tool}")

    _run(cmd, "metilene")

    # ── Step 3: R annotation ─────────────────────────────────────────────────
    r_script = os.path.join(os.path.dirname(__file__), "dmr_annotation.r")
    _run(f"Rscript {r_script} {raw_bed} {ann_bed}", "dmr_annotation.r")

    disp(f"[dmr] raw   → {raw_bed}")
    disp(f"[dmr] annot → {ann_bed}")
    return raw_bed, ann_bed
