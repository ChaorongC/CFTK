"""End-motif analysis — wraps finaletoolkit end-motifs CLI."""

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from util import recorded_run


def _run(cmd, label="end_motif"):
    print(f"[{label}] {cmd[:120]}", file=sys.stderr)
    ret = recorded_run(cmd, shell=True, label=label)
    if ret.returncode != 0:
        print(f"[{label}] ERROR: command failed.", file=sys.stderr)


def run_end_motif(args):
    """Run finaletoolkit end-motifs for each input BAM, in parallel."""
    out_dir = args.end_motif_out
    os.makedirs(out_dir, exist_ok=True)

    kmer    = getattr(args, "kmer",      4)
    mapq    = getattr(args, "mapq",     30)
    minf    = getattr(args, "min_frag", 100)
    maxf    = getattr(args, "max_frag", 220)
    extra   = getattr(args, "em_extra", "")
    workers = getattr(args, "parallel",  1) or 1

    def _process_one(bam):
        sample  = os.path.splitext(os.path.basename(bam))[0].replace(".markdup", "")
        out_tsv = os.path.join(out_dir, f"{sample}_{kmer}mer.tsv")
        if os.path.exists(out_tsv):
            print(f"[end_motif] {sample} — already done, skipping", file=sys.stderr)
            return out_tsv
        cmd = (
            f"finaletoolkit end-motifs "
            f"-k {kmer} -min {minf} -max {maxf} "
            f"-o {out_tsv} -q {mapq} -v {extra} "
            f"{bam} {args.genome2bit}"
        )
        _run(cmd, "end_motif")
        print(f"[end_motif] saved → {out_tsv}")
        return out_tsv

    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_one, bam): bam for bam in args.infile}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                results.append(res)
    return results
