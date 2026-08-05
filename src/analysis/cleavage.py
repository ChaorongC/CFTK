"""Cleavage profile analysis — wraps finaletoolkit cleavage-profile CLI."""

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from util import recorded_run


def _run(cmd, label="cleavage"):
    print(f"[{label}] {cmd[:120]}", file=sys.stderr)
    ret = recorded_run(cmd, shell=True, label=label)
    if ret.returncode != 0:
        print(f"[{label}] ERROR: command failed.", file=sys.stderr)


def run_cleavage(args):
    """Run finaletoolkit cleavage-profile for each input BAM, in parallel."""
    out_dir = args.cleavage_out
    os.makedirs(out_dir, exist_ok=True)

    if not args.bed:
        sys.exit("[cleavage] ERROR: reference_data.ctcf_bed is required.")

    minf    = getattr(args, "min_frag", 100)
    maxf    = getattr(args, "max_frag", 220)
    mapq    = getattr(args, "cl_mapq",  getattr(args, "mapq", 30))
    window  = getattr(args, "window",   20)
    extra   = getattr(args, "cl_extra", "")
    workers = getattr(args, "parallel", 1) or 1

    def _process_one(bam):
        sample = os.path.splitext(os.path.basename(bam))[0].replace(".markdup", "")
        out_bw = os.path.join(out_dir, f"{sample}_cleavage.bw")
        if os.path.exists(out_bw):
            print(f"[cleavage] {sample} — already done, skipping", file=sys.stderr)
            return out_bw
        cmd = (
            f"finaletoolkit cleavage-profile "
            f"-o {out_bw} -min {minf} -max {maxf} "
            f"-q {mapq} -w {window} -v {extra} "
            f"{bam} {args.bed} {args.chrom_sizes}"
        )
        _run(cmd, "cleavage")
        print(f"[cleavage] saved → {out_bw}")
        return out_bw

    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_one, bam): bam for bam in args.infile}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                results.append(res)
    return results
