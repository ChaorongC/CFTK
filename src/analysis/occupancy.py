"""Nucleosome occupancy via DANPOS3. Parallel per-sample, auto-merges matrix."""

import os
import sys
import glob
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from util import disp, run_command


def run_occupancy(args):
    """
    DANPOS3 dpos → wigToBigWig → bigWigAverageOverBed per sample, in parallel.
    Output dir  : 4_fragmentomics/occupancy/
    Per-sample  : {out_dir}/{name}.occupancy.tsv
    Matrix      : {out_dir}/occupancy_matrix.tsv  (auto when >1 sample)
    """
    out_dir    = args.occ_out
    chrom_sz   = args.chrom_sizes
    region_bed = args.region
    tool       = getattr(args, "danpos",       "danpos.py")
    tool_cmd   = f"python {tool}" if tool.endswith(".py") else tool
    extra      = getattr(args, "danpos_extra", "--paired 1 -u 0 -c 1000000")
    workers    = getattr(args, "parallel",     1) or 1
    os.makedirs(out_dir, exist_ok=True)

    def _process_one(bam):
        if not os.path.exists(bam):
            disp(f"[occupancy] WARNING: BAM not found: {bam}")
            return None
        stem       = os.path.splitext(os.path.basename(bam))[0]
        name       = stem.replace(".markdup", "")
        occ        = os.path.join(out_dir, f"{name}.occupancy.tsv")
        if os.path.exists(occ):
            disp(f"[occupancy] {name} — already done, skipping")
            return (occ, name)
        danpos_tmp = os.path.join(out_dir, f"danpos_tmp_{name}")
        os.makedirs(danpos_tmp, exist_ok=True)
        wig = os.path.join(danpos_tmp, "pooled", f"{stem}.Fnor.smooth.wig")
        bw  = os.path.join(out_dir, f"{name}.bw")
        cmd = (
            f"{tool_cmd} dpos {bam} {extra} -o {danpos_tmp} && "
            f"wigToBigWig -clip {wig} {chrom_sz} {bw} && "
            f"bigWigAverageOverBed {bw} {region_bed} {occ} || exit 1"
        )
        run_command(cmd, label=f"occupancy [{name}]")
        disp(f"[occupancy] saved → {occ}")
        return (occ, name)

    raw = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_process_one, bam): bam for bam in args.infile}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                raw.append(res)

    results = [r for r in raw if r is not None]
    if len(results) > 1:
        _merge_occupancy(results, out_dir)
    elif len(results) == 1:
        disp("[occupancy] Single sample — skipping matrix merge.")
    return [r[0] for r in results]


def _merge_occupancy(occ_col_pairs, out_dir):
    out_path = os.path.join(out_dir, "occupancy_matrix.tsv")
    disp(f"[occupancy] merging {len(occ_col_pairs)} files → {out_path}")
    frames = {}
    for fp, col_name in occ_col_pairs:
        df = pd.read_csv(fp, sep="\t", header=None,
                         usecols=[0, 5], names=["region", col_name])
        frames[col_name] = df.set_index("region")[col_name]
    matrix = pd.DataFrame(frames)
    matrix.to_csv(out_path, sep="\t")
    disp(f"[occupancy] {matrix.shape[0]} regions × {matrix.shape[1]} samples → {out_path}")
    return out_path
