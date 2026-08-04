"""
QC analysis — three independent steps driven by cftk_init.json.
Step 0: parse all process outputs → qc_summary.tsv + qc_scores.tsv  (NEW M1)
Step 1: methylation beta-value distribution  (requires cpg_matrix.tsv)
Step 2: fragment length distribution          (requires BAM files) — NOW PARALLEL + CHECKPOINT
Step 3: dinucleotide frequency               (requires BAM files + reference) — NOW PARALLEL + CHECKPOINT
"""

import os
import glob
import subprocess
import sys
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from util import disp


def run_qc(args):
    """
    Dispatcher: route to the correct QC computation based on args.step.
    Step 0: collect QC metrics from all process output files and score.
    Steps 1-3: unchanged purpose, 2 and 3 now parallel with checkpoints.
    """
    step = args.step
    os.makedirs(args.output_dir, exist_ok=True)

    if step == 0:
        return _run_qc_summary(args)
    elif step == 1:
        return _run_meth_distribution(args)
    elif step == 2:
        return _run_fragment_length(args)
    elif step == 3:
        return _run_dinucleotide(args)
    else:
        sys.exit(f"[qc] Invalid step: {step}. Valid: 0, 1, 2, 3.")


# ── Step 0: QC summary (NEW — M1) ────────────────────────────────────────────

def _run_qc_summary(args):
    """
    Parse all process output files and compute per-sample QC scores.

    Run-order guarantee: must run AFTER step 2 (fragment length) so that
    fragment_length.*.raw.csv files exist for median_frag_len.
    In run-all this is enforced by pipeline ordering in cftk.py
    (step 2 executes before step 0).
    When run manually (`cftk qc -s 0`), fragment data is silently omitted
    if step 2 hasn't run yet — no user prompt, no error.

    Checkpoint: skips if both output files exist AND no fragment CSV is
    newer than the existing summary (i.e., step 2 hasn't updated data).
    Pass --force to bypass.
    """
    from analysis.qc_parser import collect_qc_metrics
    from analysis.qc_scorer import score_samples

    out_dir     = args.output_dir
    summary_tsv = os.path.join(out_dir, "qc_summary.tsv")
    scores_tsv  = os.path.join(out_dir, "qc_scores.tsv")

    # Checkpoint: skip only if outputs exist AND fragment data hasn't changed
    if (
        os.path.exists(summary_tsv)
        and os.path.exists(scores_tsv)
        and not getattr(args, "force", False)
    ):
        frag_dir  = os.path.join(out_dir, "2_fragment_length")
        frag_csvs = glob.glob(os.path.join(frag_dir, "*.raw.csv"))
        sum_mtime = os.path.getmtime(summary_tsv)
        # Re-run if any fragment CSV is newer than the existing summary
        # (means step 2 completed after the last step 0 run)
        newer = [f for f in frag_csvs if os.path.getmtime(f) > sum_mtime]
        if not newer:
            disp(f"[qc step0] already done — {summary_tsv}")
            return {"summary": summary_tsv, "scores": scores_tsv}
        disp(f"[qc step0] {len(newer)} fragment CSV(s) updated since last run — re-collecting")

    samples = getattr(args, "all_samples", [])
    if not samples:
        disp("[qc step0] WARNING: no samples found in args.all_samples")
        return {}

    paths   = args.paths
    cores   = getattr(args, "cores", 1)

    # Collect raw metrics
    summary_df = collect_qc_metrics(samples, paths, summary_tsv, cores=cores)

    # Score and write qc_scores.tsv
    score_samples(summary_df, out_tsv=scores_tsv)

    return {"summary": summary_tsv, "scores": scores_tsv}


# ── Step 1: Methylation distribution ─────────────────────────────────────────

def _run_meth_distribution(args):
    """
    Load cpg_matrix.tsv and return path for visualization.
    No additional computation needed — plot_qc reads the matrix directly.
    """
    matrix = os.path.join(args.matrices_dir, "cpg_matrix.tsv")
    if not os.path.exists(matrix):
        sys.exit(
            f"[qc step1] cpg_matrix.tsv not found: {matrix}\n"
            "Run 'cftk merge --modality cpg' first."
        )
    disp(f"[qc step1] cpg matrix: {matrix}")
    return {"matrix": matrix}


# ── Step 2: Fragment length distribution — PARALLEL + CHECKPOINT (M1b) ────────

def _run_fragment_length(args):
    """
    Run bamPEFragmentSize for each BAM in parallel.
    Per-sample checkpoint: skip if raw CSV already exists.
    """
    bams    = args.infile
    out_dir = args.output_dir
    cores   = getattr(args, "cores", 1)
    workers = getattr(args, "parallel", 1) or 1
    sub_dir = os.path.join(out_dir, "2_fragment_length")
    os.makedirs(sub_dir, exist_ok=True)
    prefix  = os.path.join(sub_dir, "fragment_length")

    if not bams:
        sys.exit("[qc step2] No BAM files provided.")

    def _process_one(bam):
        if not os.path.exists(bam):
            disp(f"[qc step2] WARNING: BAM not found: {bam}")
            return None
        stem    = os.path.splitext(os.path.basename(bam))[0]
        raw_csv = f"{prefix}.{stem}.raw.csv"
        # Checkpoint: skip if already done
        if os.path.exists(raw_csv) and not getattr(args, "force", False):
            disp(f"[qc step2] {stem} — already done, skipping")
            return raw_csv
        hist_png = f"{prefix}.{stem}.hist.png"
        cmd = (
            f"bamPEFragmentSize "
            f"--outRawFragmentLengths {raw_csv} "
            f"-hist {hist_png} "
            f"-p {cores} -b {bam}"
        )
        disp(f"[qc step2] {stem}")
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode != 0:
            disp(f"[qc step2] WARNING: bamPEFragmentSize failed for {bam}")
            return None
        return raw_csv

    raw_files = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_one, bam): bam for bam in bams}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                raw_files.append(res)

    raw_files = sorted(raw_files)
    disp(f"[qc step2] {len(raw_files)} raw fragment length files ready.")
    return {"prefix": prefix, "raw_files": raw_files}


# ── Step 3: Dinucleotide frequency — PARALLEL bamtobed + CHECKPOINT (M1c) ─────

def _run_dinucleotide(args):
    """
    Extract dinucleotide frequencies around fragment centres using bedtools nuc.
    bamtobed now runs in parallel across BAMs.
    Per-sample checkpoint on fragment BED; dinucleotide files also checkpointed.
    """
    bams      = args.infile
    ref_fa    = args.ref_fa
    frag_len  = getattr(args, "fragment",  167)
    cores     = getattr(args, "cores",     1)
    workers   = getattr(args, "parallel",  1) or 1
    out_dir   = args.output_dir
    sub_dir   = os.path.join(out_dir, "3_dinucleotide_freq")
    os.makedirs(sub_dir, exist_ok=True)
    prefix    = os.path.join(sub_dir, "dinucleotide")

    if not ref_fa or not os.path.exists(ref_fa):
        sys.exit(
            "[qc step3] reference_data.genome_fa is required for dinucleotide analysis."
        )

    all_frag_file = f"{prefix}.all_fragment"
    window_file   = f"{prefix}.all_fragment.window2bp"

    # Per-BAM fragment extraction — parallel (M1c)
    # Each BAM writes to its own temp file to avoid parallel >> race condition
    def _extract_one_bam(bam):
        if not os.path.exists(bam):
            disp(f"[qc step3] WARNING: BAM not found: {bam}")
            return None
        sample    = os.path.splitext(os.path.basename(bam))[0]
        frag_done = f"{prefix}.{sample}.frags_done"
        # Checkpoint per sample
        if os.path.exists(frag_done) and not getattr(args, "force", False):
            disp(f"[qc step3] {sample} — fragments already extracted, skipping")
            return frag_done
        # Write to per-sample temp file — avoids parallel >> line-interleaving
        sample_frag_tmp = f"{prefix}.{sample}.frags_tmp"
        cmd = (
            f"bedtools bamtobed -bedpe -mate1 -i {bam} 2>/dev/null | "
            f"awk -v OFS='\\t' -v sample={sample} "
            f"'{{"
            f"start = ($2 < $5) ? $2 : $5;"
            f"end   = ($3 > $6) ? $3 : $6;"
            f"print $1, start, end, sample;"
            f"}}' | "
            f"awk -v OFS='\\t' '$3-$2=={frag_len} {{print}}' "
            f"> {sample_frag_tmp} || exit 1"
        )
        disp(f"[qc step3] bamtobed: {sample}")
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode != 0:
            disp(f"[qc step3] WARNING: bamtobed failed for {bam}")
            return None
        open(frag_done, "w").close()
        return frag_done

    # Run bamtobed in parallel (each writes its own tmp file)
    if not (os.path.exists(all_frag_file) and not getattr(args, "force", False)):
        if os.path.exists(all_frag_file):
            os.remove(all_frag_file)
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futures = {ex.submit(_extract_one_bam, bam): bam for bam in bams}
            for fut in as_completed(futures):
                fut.result()  # propagate exceptions

        # Serial concatenation after all parallel jobs finish — no race condition
        disp(f"[qc step3] concatenating per-sample fragment files → {all_frag_file}")
        with open(all_frag_file, "wb") as out_fh:
            for bam in bams:
                sample = os.path.splitext(os.path.basename(bam))[0]
                tmp_f  = f"{prefix}.{sample}.frags_tmp"
                if os.path.exists(tmp_f):
                    with open(tmp_f, "rb") as in_fh:
                        out_fh.write(in_fh.read())
                    os.remove(tmp_f)
    else:
        disp(f"[qc step3] all_fragment file exists — skipping bamtobed")

    if not os.path.exists(all_frag_file) or os.path.getsize(all_frag_file) == 0:
        sys.exit(f"[qc step3] No fragments extracted to {all_frag_file}.")

    # Build 2bp sliding window BED (checkpoint on window_file)
    if not (os.path.exists(window_file) and not getattr(args, "force", False)):
        if os.path.exists(window_file):
            os.remove(window_file)
        frags   = pd.read_table(all_frag_file, header=None,
                                names=["chrom", "start", "end", "sample"],
                                usecols=[0, 1, 2, 3],
                                on_bad_lines="skip",
                                dtype={"chrom": str, "start": int,
                                       "end": int, "sample": str})
        half_w  = (250 - frag_len) // 2
        with open(window_file, "ab") as fh:
            for _, row in frags.iterrows():
                start_w = int(row["start"]) - half_w
                pos     = np.arange(250)
                arr     = np.column_stack([
                    [str(row["chrom"])] * 250,
                    (start_w + pos).astype(str),
                    (start_w + pos + 2).astype(str),
                    pos.astype(str),
                    [str(row["sample"])] * 250,
                ])
                np.savetxt(fh, arr, delimiter="\t", fmt="%s")
        disp(f"[qc step3] window2bp file written → {window_file}")
    else:
        disp(f"[qc step3] window2bp file exists — skipping")

    # Dinucleotide counts — parallel across patterns (unchanged logic)
    dinu_list = ["AA", "AT", "TA", "TT", "GG", "GC", "CG", "CC"]

    def _nuc(pattern):
        out_file = f"{prefix}.all_fragment_{pattern}.txt"
        # Checkpoint: skip only if file exists AND is non-empty
        if (os.path.exists(out_file)
                and os.path.getsize(out_file) > 0
                and not getattr(args, "force", False)):
            disp(f"[qc step3] {pattern} — already done, skipping")
            return
        # Remove stale empty file before writing
        if os.path.exists(out_file):
            os.remove(out_file)
        # Capture stderr so we can report the actual bedtools error
        cmd = (
            f"bedtools nuc -pattern {pattern} -C -fi {ref_fa} "
            f"-bed {window_file}"
        )
        with open(out_file, "w") as out_fh:
            ret = subprocess.run(cmd, shell=True, stdout=out_fh,
                                 stderr=subprocess.PIPE, text=True)
        if ret.returncode != 0 or not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
            stderr_msg = ret.stderr.strip() if ret.stderr else "(no stderr)"
            # Clean up empty output file
            if os.path.exists(out_file) and os.path.getsize(out_file) == 0:
                os.remove(out_file)
            raise RuntimeError(
                f"[qc step3] bedtools nuc failed for pattern {pattern}.\n"
                f"  ref_fa:      {ref_fa}\n"
                f"  window_file: {window_file} "
                f"(exists={os.path.exists(window_file)}, "
                f"size={os.path.getsize(window_file) if os.path.exists(window_file) else 0})\n"
                f"  stderr: {stderr_msg}"
            )

    from joblib import Parallel, delayed
    Parallel(n_jobs=cores)(delayed(_nuc)(p) for p in dinu_list)

    disp(f"[qc step3] dinucleotide files written to {sub_dir}")
    return {"prefix": prefix, "dinu_list": dinu_list}
