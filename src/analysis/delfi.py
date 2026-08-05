"""DELFI analysis — wraps finaletoolkit delfi CLI."""

import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from util import recorded_run


def _run(cmd, label="delfi"):
    print(f"[{label}] {cmd[:120]}", file=sys.stderr)
    ret = recorded_run(cmd, shell=True, label=label)
    if ret.returncode != 0:
        print(f"[{label}] ERROR: command failed.", file=sys.stderr)


def run_delfi(args):
    """Run finaletoolkit delfi for each input BAM, in parallel."""
    out_dir = args.delfi_out
    os.makedirs(out_dir, exist_ok=True)

    mapq    = getattr(args, "delfi_mapq",   getattr(args, "mapq",    30))
    window  = getattr(args, "delfi_window", getattr(args, "window",  20))
    extra   = getattr(args, "delfi_extra",  "")
    workers = getattr(args, "parallel",     1) or 1

    def _process_one(bam):
        sample  = os.path.splitext(os.path.basename(bam))[0].replace(".markdup", "")
        out_tsv = os.path.join(out_dir, f"{sample}_delfi.tsv")
        if os.path.exists(out_tsv):
            print(f"[delfi] {sample} — already done, skipping", file=sys.stderr)
            return out_tsv
        cmd = (
            f"finaletoolkit delfi "
            f"-b {args.blacklist} -g {args.gap} "
            f"-o {out_tsv} -R -q {mapq} -w {window} -M -v {extra} "
            f"{bam} {args.chrom_sizes} {args.genome2bit} {args.bins}"
        )
        _run(cmd, "delfi")
        print(f"[delfi] saved → {out_tsv}")
        return out_tsv

    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_one, bam): bam for bam in args.infile}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                results.append(res)
    return results
