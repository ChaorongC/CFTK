"""
process.py — Part 1 raw data processing.
Steps: 1=trimming  2=alignment  3=markdup  4=methylation

Fixes vs previous version:
  - _step1_trim(): after Trim Galore runs, ALL output files are renamed to
    a consistent R1/R2 convention:
      {name}_val_1.fq.gz        → {name}_R1.fq.gz
      {name}_val_2.fq.gz        → {name}_R2.fq.gz
      {name}_val_1_fastqc.zip   → {name}_R1_fastqc.zip
      {name}_val_2_fastqc.zip   → {name}_R2_fastqc.zip
      (trimming reports already renamed correctly)
    The "Input filename:" line inside each trimming report is also rewritten
    to the sample name so MultiQC can merge all metrics into one row per sample.

  - _run_step1(): checkpoint now recognises both old val_1/val_2 and new R1/R2
    naming so already-trimmed samples are not re-processed.

  - _run_step2(): reads trimmed FASTQs using R1/R2 names first,
    falls back to val_1/val_2 for backward compatibility.

  - _run_multiqc(): step 1 now auto-generates a replace_names.tsv that maps
    original FASTQ stems → sample R1/R2 names, so MultiQC displays unified
    sample names across trimming reports and FastQC modules.

  - P1: sambamba markdup stderr captured to {name}.markdup_metrics.txt.
"""

import hashlib
import os
import re
import shlex
import sys
from pathlib import Path
from joblib import Parallel, delayed
from util import configure_command_log, disp, recorded_run, run_command
from resource_planning import (
    detect_scheduler_allocation,
    ensure_scheduler_capacity,
    plan_parallelism,
)
from init import (
    get_all_samples,
    get_bam,
    get_sequence_dictionary_path,
    get_work_paths,
    load_config,
)


DEFAULT_TARGET_BED = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed"
)

STEPS = {
    1: "Adapter trimming",
    2: "Bisulfite alignment",
    3: "Mark duplicates",
    4: "CpG methylation calling",
}

SUPPORTED_TOOLS = {
    1: {"trim_galore", "fastp"},
    2: {"bwameth", "bismark"},
    3: {"sambamba", "picard", "samblaster"},
    4: {"methyldackel", "bismark_extractor"},
}
DEFAULT_PICARD_JAVA_MEMORY = "8g"
_JAVA_MEMORY_RE = re.compile(r"^[1-9][0-9]*[kKmMgGtT]$")


def _picard_command(program, java_memory=DEFAULT_PICARD_JAVA_MEMORY):
    """Return a Picard command prefix with a validated maximum Java heap."""
    if (
        not isinstance(java_memory, str)
        or not _JAVA_MEMORY_RE.fullmatch(java_memory)
    ):
        sys.exit(
            "[process] ERROR: picard_java_memory must use a positive JVM size "
            "such as '8g' or '4096m'."
        )
    return f"picard -Xmx{java_memory.lower()} {program}"


# ── Step implementations (single sample) ─────────────────────────────────────

def _step1_trim(sample, step_cfg, ref_data, paths, per_cores):
    if sample["input_type"] != "fastq":
        disp(f"  [trim] skip BAM input: {sample['name']}")
        return None, None

    tool  = step_cfg["tool"]
    extra = step_cfg["params"].get("extra_args", "")
    name  = sample["name"]
    out   = paths["trimming"]
    os.makedirs(out, exist_ok=True)

    # Determine extension from input
    ext = "fq.gz" if sample["r1"].endswith(".gz") else "fq"

    if tool == "trim_galore":
        cmd = (
            f"trim_galore --paired --2colour 20 --cores {per_cores} "
            f"--fastqc "
            f"-o {out} --basename {name} {extra} "
            f"{sample['r1']} {sample['r2']} || exit 1"
        )
        run_command(cmd, label=f"trim [{name}]")

        # ── Unified rename: all outputs → R1/R2 convention ────────────────
        # 1. Extract original FASTQ base names (for trimming report matching)
        r1_base = os.path.splitext(os.path.splitext(
            os.path.basename(sample["r1"]))[0])[0]
        r2_base = os.path.splitext(os.path.splitext(
            os.path.basename(sample["r2"]))[0])[0]

        rename_map = [
            # Trimming reports — Trim Galore names after the INPUT file
            # Try multiple possible patterns (fastq.gz vs fq.gz)
            (f"{r1_base}.fastq.gz_trimming_report.txt", f"{name}_R1_trimming_report.txt"),
            (f"{r2_base}.fastq.gz_trimming_report.txt", f"{name}_R2_trimming_report.txt"),
            (f"{r1_base}.fq.gz_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
            (f"{r2_base}.fq.gz_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
            (f"{r1_base}_trimming_report.txt",           f"{name}_R1_trimming_report.txt"),
            (f"{r2_base}_trimming_report.txt",           f"{name}_R2_trimming_report.txt"),
            # With --basename, Trim Galore 0.6.11 uses the output mate stem.
            (f"{name}_R1.fastq.gz_trimming_report.txt", f"{name}_R1_trimming_report.txt"),
            (f"{name}_R2.fastq.gz_trimming_report.txt", f"{name}_R2_trimming_report.txt"),
            (f"{name}_R1.fq.gz_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
            (f"{name}_R2.fq.gz_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
            (f"{name}_R1.fastq_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
            (f"{name}_R2.fastq_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
            (f"{name}_R1.fq_trimming_report.txt",       f"{name}_R1_trimming_report.txt"),
            (f"{name}_R2.fq_trimming_report.txt",       f"{name}_R2_trimming_report.txt"),
            # FastQC files — Trim Galore names after the --basename output
            (f"{name}_val_1_fastqc.zip",  f"{name}_R1_fastqc.zip"),
            (f"{name}_val_1_fastqc.html", f"{name}_R1_fastqc.html"),
            (f"{name}_val_2_fastqc.zip",  f"{name}_R2_fastqc.zip"),
            (f"{name}_val_2_fastqc.html", f"{name}_R2_fastqc.html"),
            # Trimmed FASTQ files — rename val_1/val_2 → R1/R2
            (f"{name}_val_1.{ext}", f"{name}_R1.{ext}"),
            (f"{name}_val_2.{ext}", f"{name}_R2.{ext}"),
        ]

        for old_name, new_name in rename_map:
            old_path = os.path.join(out, old_name)
            new_path = os.path.join(out, new_name)
            if old_path == new_path:
                continue
            if os.path.exists(old_path):
                if os.path.exists(new_path):
                    # new already exists (e.g. from a previous run): remove old
                    os.remove(old_path)
                else:
                    os.rename(old_path, new_path)
                    disp(f"  [trim] renamed: {old_name} → {new_name}")

        # ── Rewrite "Input filename:" in trimming reports ─────────────────
        # MultiQC reads this field (not the file name) to determine sample name.
        # Replace the original FASTQ path with the cftk sample name so MultiQC
        # merges trimming-report metrics with FastQC metrics into one row.
        for report_name, new_stem in [
            (f"{name}_R1_trimming_report.txt", f"{name}_R1"),
            (f"{name}_R2_trimming_report.txt", f"{name}_R2"),
        ]:
            report_path = os.path.join(out, report_name)
            if os.path.exists(report_path):
                try:
                    with open(report_path, "r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                    # Replace "Input filename: /any/path/to/original.fq.gz"
                    new_content = re.sub(
                        r"(?m)^(Input filename:).*$",
                        rf"\1 {new_stem}",
                        content,
                    )
                    if new_content != content:
                        with open(report_path, "w", encoding="utf-8") as f:
                            f.write(new_content)
                        disp(f"  [trim] updated Input filename → {new_stem} in {report_name}")
                except Exception as e:
                    disp(f"  [trim] WARNING: could not patch {report_name}: {e}")

    elif tool == "fastp":
        r1_out = os.path.join(out, f"{name}_R1.{ext}")
        r2_out = os.path.join(out, f"{name}_R2.{ext}")
        cmd = (
            f"fastp -i {sample['r1']} -I {sample['r2']} "
            f"-o {r1_out} -O {r2_out} -w {per_cores} {extra} || exit 1"
        )
        run_command(cmd, label=f"trim [{name}]")
    else:
        sys.exit(f"[step1] unsupported tool: {tool}. Supported: {SUPPORTED_TOOLS[1]}")

    r1_out = os.path.join(out, f"{name}_R1.{ext}")
    r2_out = os.path.join(out, f"{name}_R2.{ext}")
    return r1_out, r2_out


def _step2_align(sample, r1, r2, step_cfg, ref_data, paths, per_cores):
    tool  = step_cfg["tool"]
    extra = step_cfg["params"].get("extra_args", "")
    ref   = ref_data["genome_fa"]
    name  = sample["name"]
    out   = paths["alignment"]
    os.makedirs(out, exist_ok=True)
    bam   = os.path.join(out, f"{name}.bam")

    if tool == "bwameth":
        read_group = (
            f"@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:ILLUMINA"
        )
        cmd = (
            f"bwameth.py --reference {ref} -t {per_cores} "
            f"--read-group {shlex.quote(read_group)} {extra} {r1} {r2} | "
            f"sambamba view -t {per_cores} "
            f"-F 'not secondary_alignment and not failed_quality_control "
            f"and not supplementary and proper_pair and mapping_quality > 0' "
            f"-f bam -S -l 0 /dev/stdin | "
            f"sambamba sort -t {per_cores} -o {bam} /dev/stdin || exit 1; "
            f"samtools index -@ {per_cores} {bam} || exit 1"
        )
    elif tool == "bismark":
        cmd = (
            f"bismark --genome {os.path.dirname(ref)} -1 {r1} -2 {r2} "
            f"-p {per_cores} -o {out} {extra} || exit 1"
        )
    else:
        sys.exit(f"[step2] unsupported tool: {tool}. Supported: {SUPPORTED_TOOLS[2]}")

    run_command(cmd, label=f"align [{name}]")
    run_command(
        f"samtools flagstat -@ {per_cores} {bam} > {bam}.flagstat",
        label=f"flagstat [{name}]"
    )
    run_command(
        f"samtools stats   -@ {per_cores} {bam} > {bam}.stats",
        label=f"samtools_stats [{name}]"
    )
    return bam


def _step3_markdup(sample, bam_in, step_cfg, ref_data, paths, per_cores):
    tool    = step_cfg["tool"]
    extra   = step_cfg["params"].get("extra_args", "")
    picard_java_memory = step_cfg["params"].get(
        "picard_java_memory", DEFAULT_PICARD_JAVA_MEMORY
    )
    ref     = ref_data["genome_fa"]
    name    = sample["name"]
    out     = paths["markdup"]
    os.makedirs(out, exist_ok=True)
    bam_out = os.path.join(out, f"{name}.markdup.bam")

    if tool == "sambamba":
        # P1: capture sambamba markdup stderr → {name}.markdup_metrics.txt
        metrics_txt = os.path.join(out, f"{name}.markdup_metrics.txt")
        cmd = (
            f"sambamba markdup -t {per_cores} {extra} {bam_in} {bam_out} "
            f"2>{metrics_txt} || exit 1; "
            f"samtools index -@ {per_cores} {bam_out} || exit 1"
        )
    elif tool == "picard":
        metrics = bam_out.replace(".bam", "_metrics.txt")
        cmd = (
            f"{_picard_command('MarkDuplicates', picard_java_memory)} "
            f"I={bam_in} O={bam_out} R={ref} "
            f"M={metrics} SORTING_COLLECTION_SIZE_RATIO=0.15 "
            f"ASSUME_SORT_ORDER=coordinate OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 "
            f"MAX_RECORDS_IN_RAM=1000 {extra} || exit 1; "
            f"samtools index -@ {per_cores} {bam_out} || exit 1"
        )
    elif tool == "samblaster":
        cmd = (
            f"samblaster {extra} --addMateTags "
            f"< {bam_in} > {bam_out} || exit 1; "
            f"samtools index -@ {per_cores} {bam_out} || exit 1"
        )
    else:
        sys.exit(f"[step3] unsupported tool: {tool}. Supported: {SUPPORTED_TOOLS[3]}")

    run_command(cmd, label=f"markdup [{name}]")
    return bam_out


def _resolve_target_bed(target_bed=None):
    """Resolve an explicit target BED or the repository Twist panel default."""
    candidate = Path(target_bed).expanduser() if target_bed else DEFAULT_TARGET_BED
    if not candidate.is_file():
        source = f"requested target BED {candidate}" if target_bed else (
            f"default Twist target BED {candidate}"
        )
        sys.exit(
            f"[step3] ERROR: {source} was not found. Provide --target-bed PATH "
            "or use --skip-picard-metrics for a workflow that does not require "
            "Twist panel metrics."
        )
    return str(candidate.resolve())


def _prepare_picard_interval_list(
    target_bed,
    reference_fa,
    paths,
    picard_java_memory=DEFAULT_PICARD_JAVA_MEMORY,
):
    """Create the target interval list once for all post-markdup metrics."""
    sequence_dict = get_sequence_dictionary_path(reference_fa)
    if not os.path.isfile(sequence_dict):
        sys.exit(
            f"[step3] ERROR: Picard sequence dictionary not found: {sequence_dict}. "
            "Run 'cftk init' to prepare the reference."
        )

    stem = os.path.basename(target_bed)
    if stem.endswith(".bed"):
        stem = stem[:-4]
    digest = hashlib.sha256()
    for input_path in (target_bed, sequence_dict):
        with open(input_path, "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    profile_key = f"{stem}.{digest.hexdigest()[:12]}"
    metrics_dir = os.path.join(paths["markdup"], "picard_metrics", profile_key)
    os.makedirs(metrics_dir, exist_ok=True)
    interval_list = os.path.join(metrics_dir, f"{stem}.interval_list")
    if os.path.exists(interval_list):
        disp(f"[picard] interval list already present: {interval_list}")
        return interval_list

    cmd = (
        f"{_picard_command('BedToIntervalList', picard_java_memory)} "
        f"I={shlex.quote(target_bed)} "
        f"O={shlex.quote(interval_list)} SD={shlex.quote(sequence_dict)} || exit 1"
    )
    run_command(cmd, label="picard BedToIntervalList")
    return interval_list


def _run_picard_metrics(
    sample,
    bam_in,
    reference_fa,
    interval_list,
    paths,
    picard_java_memory=DEFAULT_PICARD_JAVA_MEMORY,
):
    """Collect Twist target and general Picard performance metrics for a BAM."""
    name = sample["name"]
    metrics_dir = os.path.dirname(interval_list)
    os.makedirs(metrics_dir, exist_ok=True)
    hs_metrics = os.path.join(metrics_dir, f"{name}.hs_metrics.txt")
    per_target = os.path.join(metrics_dir, f"{name}.per_target_coverage.txt")
    multiple_prefix = os.path.join(metrics_dir, f"{name}.multiple_metrics")
    multiple_done = f"{multiple_prefix}.done"

    if os.path.exists(hs_metrics) and os.path.exists(per_target):
        disp(f"  [picard] {name} — HsMetrics already present, skipping")
    else:
        cmd = (
            f"{_picard_command('CollectHsMetrics', picard_java_memory)} "
            f"I={shlex.quote(bam_in)} "
            f"O={shlex.quote(hs_metrics)} R={shlex.quote(reference_fa)} "
            f"BAIT_INTERVALS={shlex.quote(interval_list)} "
            f"TARGET_INTERVALS={shlex.quote(interval_list)} "
            f"PER_TARGET_COVERAGE={shlex.quote(per_target)} "
            f"MINIMUM_MAPPING_QUALITY=20 COVERAGE_CAP=1000 "
            f"NEAR_DISTANCE=500 || exit 1"
        )
        run_command(cmd, label=f"picard CollectHsMetrics [{name}]")

    if os.path.exists(multiple_done):
        disp(f"  [picard] {name} — multiple metrics already present, skipping")
    else:
        cmd = (
            f"{_picard_command('CollectMultipleMetrics', picard_java_memory)} "
            f"I={shlex.quote(bam_in)} "
            f"O={shlex.quote(multiple_prefix)} R={shlex.quote(reference_fa)} "
            f"PROGRAM=CollectGcBiasMetrics "
            f"PROGRAM=CollectInsertSizeMetrics "
            f"PROGRAM=CollectAlignmentSummaryMetrics "
            f"|| exit 1; touch {shlex.quote(multiple_done)}"
        )
        run_command(cmd, label=f"picard CollectMultipleMetrics [{name}]")


def _step4_methylation(sample, bam_in, step_cfg, ref_data, paths, per_cores):
    tool  = step_cfg["tool"]
    depth = step_cfg["params"].get("min_depth", 10)
    extra = step_cfg["params"].get("extra_args", "")
    ref   = ref_data["genome_fa"]
    name  = sample["name"]
    out   = paths["methylation"]
    os.makedirs(out, exist_ok=True)
    prefix = os.path.join(out, name)

    if tool == "methyldackel":
        mbias_txt  = f"{prefix}_mbias.txt"
        mbias_temp = f"{prefix}_mbias_OT_OB.temp"
        # P3a: mbias --txt outputs per-position TSV to stdout; OT/OB coords to stderr.
        #      Capture stdout/stderr instead of using shell redirection so the
        #      M-bias table and OT/OB bounds remain independently available.
        # run_command() is used for the extract step (no stdout capture needed).
        import re as _re

        # Step A: mbias --txt → capture stdout (TSV) and stderr (OT/OB coords)
        mbias_args = [
            "MethylDackel", "mbias", "--txt",
            "-@", str(per_cores),
            ref, bam_in, prefix,
        ]
        disp(f"  [mbias] {name} — running mbias --txt")
        mbias_proc = recorded_run(
            mbias_args, capture_output=True, text=True, label=f"mbias [{name}]"
        )

        # Write stdout (TSV) to _mbias.txt
        with open(mbias_txt, "w") as _f:
            _f.write(mbias_proc.stdout)

        # Write stderr (OT/OB coords) to _mbias_OT_OB.temp
        with open(mbias_temp, "w") as _f:
            _f.write(mbias_proc.stderr)

        if mbias_proc.returncode != 0:
            disp(f"  [mbias] ERROR (rc={mbias_proc.returncode}): {mbias_proc.stderr[:300]}")
            import sys as _sys; _sys.exit(f"[step4] mbias failed for {name}")

        # Validate: stdout must start with "Strand" header
        first_line = mbias_proc.stdout.splitlines()[0] if mbias_proc.stdout.strip() else ""
        if not first_line.startswith("Strand"):
            disp(f"  [mbias] WARNING: stdout does not look like --txt TSV "
                 f"(first line: {first_line[:80]!r})")
            disp(f"  [mbias] stderr: {mbias_proc.stderr[:200]}")

        # Extract OT/OB from stderr for the extract step
        ot_ob_m = _re.search(r"--OT\s+(\S+)\s+--OB\s+(\S+)", mbias_proc.stderr)
        if not ot_ob_m:
            sys.exit(
                f"[step4] could not parse OT/OB inclusion bounds for {name}; "
                f"inspect {mbias_temp}. Methylation extraction was not run."
            )
        ot_flag = f"--OT {ot_ob_m.group(1)}"
        ob_flag = f"--OB {ot_ob_m.group(2)}"

        # Step B: CpG extract
        cmd = (
            f"MethylDackel extract --mergeContext --minDepth {depth} "
            f"--maxVariantFrac 0.25 "
            f"-@ {per_cores} "
            f"{ot_flag} {ob_flag} "
            f"-o {prefix} {extra} {ref} {bam_in} || exit 1"
        )
    elif tool == "bismark_extractor":
        cmd = (
            f"bismark_methylation_extractor --paired-end --gzip "
            f"--CpG_only --cytosine_report "
            f"--genome_folder {os.path.dirname(ref)} "
            f"-o {out} {extra} {bam_in} || exit 1"
        )
    else:
        sys.exit(f"[step4] unsupported tool: {tool}. Supported: {SUPPORTED_TOOLS[4]}")

    run_command(cmd, label=f"methylation [{name}]")
    return f"{prefix}_CpG.bedGraph"


# ── CPG merge (auto after step4) ──────────────────────────────────────────────

def _merge_cpg(bedgraph_files, samples, paths):
    import shutil
    import pandas as pd

    out_dir = paths["cpg_matrix"]
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "cpg_matrix.tsv")

    name_map = {s["name"]: s["name"] for s in samples}
    col_names = []
    for fp in bedgraph_files:
        sample_name = os.path.basename(fp).replace("_CpG.bedGraph", "")
        col_names.append(name_map.get(sample_name, sample_name))
    tmp_path = None

    if len(bedgraph_files) == 1:
        rows = []
        with open(bedgraph_files[0], "r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped_line = line.strip()
                if not stripped_line or stripped_line.startswith(
                    ("track", "browser", "#")
                ):
                    continue
                fields = stripped_line.split("\t")
                if len(fields) < 4:
                    sys.exit(
                        f"[merge] ERROR: invalid bedGraph row at "
                        f"{bedgraph_files[0]}:{line_number}; expected at "
                        "least four tab-separated columns."
                    )
                rows.append(fields[:4])
        if not rows:
            sys.exit(
                f"[merge] ERROR: no CpG rows found in {bedgraph_files[0]}."
            )
        matrix = pd.DataFrame(
            rows, columns=["chrom", "start", "end", col_names[0]]
        )
        disp(f"[merge] single bedGraph: 1 file → {out_path}")
    else:
        bedtools = shutil.which("bedtools")
        if bedtools is None:
            sys.exit("[merge] ERROR: bedtools not found in PATH.")

        col_names_str = " ".join(col_names)
        stripped = " ".join(
            f"<(grep -v '^track' {fp} | cut -f1-4)" for fp in bedgraph_files
        )
        tmp_path = out_path + ".tmp"
        cmd = (
            f"bash -c '"
            f"{bedtools} unionbedg -header -names {col_names_str} -filler NA "
            f"-i {stripped}"
            f" > {tmp_path}'"
        )

        disp(
            f"[merge] bedtools unionbedg: {len(bedgraph_files)} files → "
            f"{out_path}"
        )
        ret = recorded_run(cmd, shell=True, label="merge CpG matrix")
        if ret.returncode != 0:
            sys.exit("[merge] ERROR: bedtools unionbedg failed.")

        matrix = pd.read_csv(tmp_path, sep="\t")
        if matrix.shape[1] < 4:
            sys.exit("[merge] ERROR: bedtools unionbedg returned an invalid table.")
    chrom_col, start_col, end_col = matrix.columns[:3]
    # MethylDackel merged CpGs span C and G. Keep the established 1-based
    # cytosine coordinate rather than using BED end, which would identify G.
    matrix.insert(1, "pos", pd.to_numeric(matrix[start_col]) + 1)
    matrix = matrix.drop(columns=[start_col, end_col])
    pos_col = "pos"
    matrix.index = matrix[chrom_col].astype(str) + "_" + matrix[pos_col].astype(str)
    matrix.index.name = "cpg_id"
    matrix = matrix.drop(columns=[chrom_col, pos_col])
    n_before = len(matrix)
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    n_dropped = n_before - len(matrix)
    if n_dropped > 0:
        disp(f"[merge] dropped {n_dropped} duplicate CpG positions")
    matrix.to_csv(out_path, sep="\t")
    if tmp_path is not None:
        os.remove(tmp_path)
    disp(f"[merge] {matrix.shape[0]} CpGs × {matrix.shape[1]} samples → {out_path}")
    return out_path


# ── Per-step single-sample runners ────────────────────────────────────────────

def _run_step1(sample, proc, ref, paths, per_cores):
    """
    Trimming checkpoint: recognises both new R1/R2 naming and legacy val_1/val_2,
    so samples trimmed before the rename fix are not re-processed.
    """
    name = sample["name"]
    ext  = "fq.gz" if sample.get("r1", "").endswith(".gz") else "fq"

    # New naming (post-fix)
    r1_new = os.path.join(paths["trimming"], f"{name}_R1.{ext}")
    r2_new = os.path.join(paths["trimming"], f"{name}_R2.{ext}")
    # Legacy naming (pre-fix)
    r1_old = os.path.join(paths["trimming"], f"{name}_val_1.{ext}")
    r2_old = os.path.join(paths["trimming"], f"{name}_val_2.{ext}")

    if os.path.exists(r1_new) and os.path.exists(r2_new):
        disp(f"  [step1] {name} — already done, skipping")
        return r1_new, r2_new

    if os.path.exists(r1_old) and os.path.exists(r2_old):
        # Legacy files exist: apply the rename fix now without re-trimming
        disp(f"  [step1] {name} — legacy val_1/val_2 found, applying rename fix")
        _apply_trim_rename_fix(sample, paths["trimming"], ext)
        return (
            r1_new if os.path.exists(r1_new) else r1_old,
            r2_new if os.path.exists(r2_new) else r2_old,
        )

    disp(f"  [step1] {name}")
    return _step1_trim(sample, proc["step1_trimming"], ref, paths, per_cores)


def _apply_trim_rename_fix(sample, out_dir: str, ext: str):
    """
    Apply the rename + trimming-report content fix to an already-trimmed sample
    without re-running Trim Galore.  Used by _run_step1 when legacy files exist.
    """
    name    = sample["name"]
    r1_base = os.path.splitext(os.path.splitext(
        os.path.basename(sample["r1"]))[0])[0]
    r2_base = os.path.splitext(os.path.splitext(
        os.path.basename(sample["r2"]))[0])[0]

    rename_map = [
        (f"{r1_base}.fastq.gz_trimming_report.txt", f"{name}_R1_trimming_report.txt"),
        (f"{r2_base}.fastq.gz_trimming_report.txt", f"{name}_R2_trimming_report.txt"),
        (f"{r1_base}.fq.gz_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
        (f"{r2_base}.fq.gz_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
        (f"{r1_base}_trimming_report.txt",           f"{name}_R1_trimming_report.txt"),
        (f"{r2_base}_trimming_report.txt",           f"{name}_R2_trimming_report.txt"),
        (f"{name}_R1.fastq.gz_trimming_report.txt", f"{name}_R1_trimming_report.txt"),
        (f"{name}_R2.fastq.gz_trimming_report.txt", f"{name}_R2_trimming_report.txt"),
        (f"{name}_R1.fq.gz_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
        (f"{name}_R2.fq.gz_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
        (f"{name}_R1.fastq_trimming_report.txt",    f"{name}_R1_trimming_report.txt"),
        (f"{name}_R2.fastq_trimming_report.txt",    f"{name}_R2_trimming_report.txt"),
        (f"{name}_R1.fq_trimming_report.txt",       f"{name}_R1_trimming_report.txt"),
        (f"{name}_R2.fq_trimming_report.txt",       f"{name}_R2_trimming_report.txt"),
        (f"{name}_val_1_fastqc.zip",  f"{name}_R1_fastqc.zip"),
        (f"{name}_val_1_fastqc.html", f"{name}_R1_fastqc.html"),
        (f"{name}_val_2_fastqc.zip",  f"{name}_R2_fastqc.zip"),
        (f"{name}_val_2_fastqc.html", f"{name}_R2_fastqc.html"),
        (f"{name}_val_1.{ext}",       f"{name}_R1.{ext}"),
        (f"{name}_val_2.{ext}",       f"{name}_R2.{ext}"),
    ]
    for old_name, new_name in rename_map:
        old_path = os.path.join(out_dir, old_name)
        new_path = os.path.join(out_dir, new_name)
        if old_path == new_path:
            continue
        if os.path.exists(old_path):
            if os.path.exists(new_path):
                os.remove(old_path)
            else:
                os.rename(old_path, new_path)
                disp(f"  [trim-fix] renamed: {old_name} → {new_name}")

    # Patch trimming report content
    for report_name, new_stem in [
        (f"{name}_R1_trimming_report.txt", f"{name}_R1"),
        (f"{name}_R2_trimming_report.txt", f"{name}_R2"),
    ]:
        report_path = os.path.join(out_dir, report_name)
        if os.path.exists(report_path):
            try:
                with open(report_path, "r", encoding="utf-8", errors="replace") as f:
                    content = f.read()
                new_content = re.sub(
                    r"(?m)^(Input filename:).*$",
                    rf"\1 {new_stem}",
                    content,
                )
                if new_content != content:
                    with open(report_path, "w", encoding="utf-8") as f:
                        f.write(new_content)
            except Exception as e:
                disp(f"  [trim-fix] WARNING: could not patch {report_name}: {e}")


def _run_step2(sample, proc, ref, paths, per_cores):
    """
    Find trimmed FASTQs using R1/R2 naming first (new), then val_1/val_2 (legacy).
    """
    name  = sample["name"]
    itype = sample["input_type"]
    disp(f"  [step2] {name}")

    if itype == "bam":
        disp(f"  [step2] skip — BAM input: {name}")
        return sample.get("bam")

    ext = "fq.gz" if sample["r1"].endswith(".gz") else "fq"

    # Resolve R1 path — prefer new R1/R2 naming, fall back to val_1/val_2
    r1 = _find_trimmed(paths["trimming"], name, "R1", "val_1", ext)
    r2 = _find_trimmed(paths["trimming"], name, "R2", "val_2", ext)

    missing = [p for p in [r1, r2] if not p or not os.path.exists(p)]
    if missing:
        sys.exit(
            f"[step2] ERROR: trimmed reads not found for '{name}'. "
            f"Expected {name}_R1.{ext} or {name}_val_1.{ext} in "
            f"{paths['trimming']}. Run step 1 first."
        )

    bam_done = os.path.join(paths["alignment"], f"{name}.bam")
    if os.path.exists(bam_done):
        disp(f"  [step2] {name} — already done, skipping")
        return bam_done
    return _step2_align(sample, r1, r2, proc["step2_alignment"],
                        ref, paths, per_cores)


def _find_trimmed(trimming_dir: str, name: str,
                  new_suffix: str, old_suffix: str, ext: str) -> str:
    """Return path of trimmed FASTQ, preferring new R1/R2 naming over val_1/val_2."""
    new_path = os.path.join(trimming_dir, f"{name}_{new_suffix}.{ext}")
    old_path = os.path.join(trimming_dir, f"{name}_{old_suffix}.{ext}")
    if os.path.exists(new_path):
        return new_path
    if os.path.exists(old_path):
        return old_path
    return new_path   # return new_path so error message is informative


def _run_step3(sample, proc, ref, paths, per_cores):
    name     = sample["name"]
    bam_done = os.path.join(paths["markdup"], f"{name}.markdup.bam")
    if os.path.exists(bam_done):
        disp(f"  [step3] {name} — already done, skipping")
        return bam_done
    disp(f"  [step3] {name}")
    bam_in = get_bam(sample, paths)
    if not (bam_in and os.path.exists(bam_in)):
        disp(f"  [step3] ERROR: BAM not found for {sample['name']}: {bam_in}")
        return None
    return _step3_markdup(sample, bam_in,
                          proc["step3_markdup"], ref, paths, per_cores)


def _run_step3_with_metrics(
    sample, proc, ref, paths, per_cores, interval_list
):
    """Run or reuse markdup, then fill any missing Picard metric outputs."""
    bam_out = _run_step3(sample, proc, ref, paths, per_cores)
    if bam_out and os.path.exists(bam_out):
        picard_java_memory = proc["step3_markdup"]["params"].get(
            "picard_java_memory", DEFAULT_PICARD_JAVA_MEMORY
        )
        _run_picard_metrics(
            sample,
            bam_out,
            ref["genome_fa"],
            interval_list,
            paths,
            picard_java_memory,
        )
    return bam_out


def _run_step4(sample, proc, ref, paths, per_cores):
    name      = sample["name"]
    prefix    = os.path.join(paths["methylation"], name)
    bg_done   = f"{prefix}_CpG.bedGraph"
    mbias_txt  = f"{prefix}_mbias.txt"
    mbias_temp = f"{prefix}_mbias_OT_OB.temp"

    # Full checkpoint: all three outputs must exist.
    # _CpG.bedGraph alone is NOT sufficient — if mbias.txt / mbias_OT_OB.temp
    # are missing or contain only legacy stderr content, we must re-run.
    def _mbias_is_legacy(path):
        """Return True if file exists but has no --txt TSV data (only stderr)."""
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return True
        try:
            with open(path) as f:
                first = f.readline()
            # --txt TSV starts with "Strand	"; legacy stderr starts with "[" or "Suggested"
            return not first.strip().startswith("Strand")
        except OSError:
            return True

    if (os.path.exists(bg_done)
            and os.path.exists(mbias_temp)
            and not _mbias_is_legacy(mbias_txt)):
        disp(f"  [step4] {name} — already done, skipping")
        return bg_done

    if os.path.exists(bg_done):
        disp(f"  [step4] {name} — CpG bedGraph exists but mbias --txt output missing, re-running")

    disp(f"  [step4] {name}")
    bam_in = get_bam(sample, paths)
    if not (bam_in and os.path.exists(bam_in)):
        disp(f"  [step4] ERROR: BAM not found for {sample['name']}: {bam_in}")
        return None
    return _step4_methylation(sample, bam_in,
                              proc["step4_methylation"], ref, paths, per_cores)


# ── Step dispatcher ───────────────────────────────────────────────────────────

_STEP_FN = {
    1: _run_step1,
    2: _run_step2,
    3: _run_step3,
    4: _run_step4,
}


# ── MultiQC ──────────────────────────────────────────────────────────────────

_MULTIQC_STEP = {
    1: "trimming",
    2: "alignment",
    # Step 4 (methylation) intentionally excluded: MethylDackel output is
    # parsed directly by qc_parser.py; MultiQC adds no value here.
}

_MULTIQC_TITLE = {
    1: "Trim Galore QC",
    2: "bwameth Alignment",
}


def _build_replace_names_tsv(samples: list, out_path: str):
    """
    Write a MultiQC replace-names TSV that maps original FASTQ stems
    to cftk sample names (e.g. S10000Nr4_R1 → Control_2_R1).

    MultiQC uses the "Input filename:" field inside trimming reports to name
    samples.  Even after we patch that field, this TSV acts as a safety net
    so MultiQC can merge all metrics (FastQC + trimming) into one row.
    """
    rows = ["Search\tReplace"]
    for s in samples:
        if s.get("input_type") != "fastq":
            continue
        name    = s["name"]
        r1_base = os.path.splitext(os.path.splitext(
            os.path.basename(s["r1"]))[0])[0]
        r2_base = os.path.splitext(os.path.splitext(
            os.path.basename(s["r2"]))[0])[0]
        rows.append(f"{r1_base}\t{name}_R1")
        rows.append(f"{r2_base}\t{name}_R2")
        # Also map the legacy val_ names in case any remain
        rows.append(f"{name}_val_1\t{name}_R1")
        rows.append(f"{name}_val_2\t{name}_R2")
    with open(out_path, "w") as f:
        f.write("\n".join(rows) + "\n")
    disp(f"[multiqc] replace_names.tsv → {out_path}")


def _run_multiqc(step, paths, samples=None):
    import shutil
    multiqc = shutil.which("multiqc")
    if not multiqc:
        disp(f"[multiqc] step {step}: multiqc not found, skipping.")
        return

    if step not in _MULTIQC_STEP:
        return
    step_dir = paths.get(_MULTIQC_STEP[step], "")
    if not step_dir or not os.path.exists(step_dir):
        disp(f"[multiqc] step {step}: directory not found ({step_dir}), skipping.")
        return

    out_dir = os.path.join(step_dir, "multiqc")
    os.makedirs(out_dir, exist_ok=True)

    title        = _MULTIQC_TITLE.get(step, f"Step {step}")
    # Step 4 (methylation): ignore large genome/BAM files; also ignore temp files
    if step == 4:
        ignore_flags = (
            "--ignore '*.fq.gz' --ignore '*.fastq.gz' "
            "--ignore '*.bam' --ignore '*.bai' "
            "--ignore '*.bedGraph' --ignore '*.temp'"
        )
    else:
        ignore_flags = (
            "--ignore '*.fq.gz' --ignore '*.fastq.gz' "
            "--ignore '*.bam' --ignore '*.bai'"
        )

    # Step 1: generate replace_names.tsv so MultiQC unifies sample names
    replace_flag = ""
    if step == 1 and samples:
        replace_tsv = os.path.join(out_dir, "replace_names.tsv")
        _build_replace_names_tsv(samples, replace_tsv)
        replace_flag = f"--replace-names {replace_tsv}"

    cmd = (
        f"{multiqc} {step_dir} "
        f"--outdir {out_dir} "
        f"--filename multiqc_report.html "
        f"--title \"{title}\" "
        f"--export "
        f"{ignore_flags} "
        f"{replace_flag} "
        f"--force "
        f"--quiet"
    )
    disp(f"[multiqc] step {step}: running MultiQC on {step_dir}")
    run_command(cmd, label=f"multiqc [step {step}]")
    html_out = os.path.join(out_dir, "multiqc_report.html")
    if os.path.exists(html_out):
        disp(f"[multiqc] step {step}: report → {html_out}")
    else:
        disp(f"[multiqc] step {step}: WARNING — multiqc_report.html not found")


# ── Main entry ────────────────────────────────────────────────────────────────

def process(args, config_path="./cftk_init.json"):
    cfg     = load_config(config_path)
    paths   = get_work_paths(cfg)
    provenance = paths.get("provenance")
    if provenance:
        configure_command_log(os.path.join(provenance, "commands.jsonl"))
    samples = get_all_samples(cfg)
    steps   = sorted(set(args.step))
    ref     = cfg["reference_data"]
    proc    = cfg["process"]

    requested_parallel = getattr(args, "parallel", None) \
                         or proc.get("parallel_samples", 1)
    step_keys = {
        1: "step1_trimming",
        2: "step2_alignment",
        3: "step3_markdup",
        4: "step4_methylation",
    }

    invalid = [s for s in steps if s not in STEPS]
    if invalid:
        sys.exit(f"[process] Invalid steps: {invalid}. Valid: {list(STEPS.keys())}")

    step_resources = {}
    scheduler = detect_scheduler_allocation()
    try:
        for step in steps:
            total = proc[step_keys[step]]["params"].get("cores", 20)
            ensure_scheduler_capacity(total, scheduler)
            step_resources[step] = plan_parallelism(
                total, requested_parallel, len(samples)
            )
    except ValueError as exc:
        sys.exit(f"[process] ERROR: invalid CPU resource plan: {exc}")

    picard_interval_list = None
    if 3 in steps and not getattr(args, "skip_picard_metrics", False):
        picard_java_memory = proc["step3_markdup"]["params"].get(
            "picard_java_memory", DEFAULT_PICARD_JAVA_MEMORY
        )
        target_bed = _resolve_target_bed(
            getattr(args, "target_bed", None) or ref.get("target_bed")
        )
        picard_interval_list = _prepare_picard_interval_list(
            target_bed, ref["genome_fa"], paths, picard_java_memory
        )

    disp(f"Samples          : {len(samples)}")
    disp(f"Steps            : {[STEPS[s] for s in steps]}")
    disp(f"Parallel samples : {requested_parallel}")
    for s in steps:
        resources = step_resources[s]
        disp(f"  Step {s} cores  : {resources['threads_per_sample']} per sample "
             f"({resources['concurrent_samples']} concurrent; total budget "
             f"{resources['total_core_budget']})")

    bedgraph_results = [None] * len(samples)

    for step in steps:
        resources = step_resources[step]
        parallel = resources["concurrent_samples"]
        sep = "=" * 50
        disp(sep)
        disp(f"[step {step}] {STEPS[step]} — {len(samples)} samples "
             f"(parallel={parallel})")
        disp(sep)

        per_cores = resources["threads_per_sample"]
        if step == 3 and picard_interval_list:
            results = Parallel(n_jobs=parallel, backend="multiprocessing")(
                delayed(_run_step3_with_metrics)(
                    s, proc, ref, paths, per_cores, picard_interval_list
                )
                for s in samples
            )
        else:
            fn = _STEP_FN[step]
            results = Parallel(n_jobs=parallel, backend="multiprocessing")(
                delayed(fn)(s, proc, ref, paths, per_cores)
                for s in samples
            )

        if step == 4:
            bedgraph_results = results

        disp(f"[step {step}] all samples complete.")
        # Pass samples to multiqc so replace_names.tsv can be generated for step 1
        _run_multiqc(step, paths, samples=samples)

    if 4 in steps:
        successful = [bg for bg in bedgraph_results
                      if bg and os.path.exists(bg)]
        if successful:
            _merge_cpg(successful, samples, paths)
        else:
            disp("[merge] WARNING: no successful bedGraph files found.")

    disp("\nAll samples processed.")
