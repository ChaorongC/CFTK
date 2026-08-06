"""
visualization.py — unified visualization dispatcher.
All plot functions save both PNG (300 dpi) and PDF.

M3d: step 0 plot removed from plot_qc() — QC summary is now an interactive
     Plotly table rendered directly in report_generator._qc_table().
"""

import os


def _out(base_dir, *parts, stem):
    d = os.path.join(base_dir, *parts)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{stem}.png"), os.path.join(d, f"{stem}.pdf")


# ── QC ────────────────────────────────────────────────────────────────────────

def plot_qc(args):
    """
    Visualization dispatcher for QC steps.

    Step 0: No plot — QC summary is rendered as an interactive Plotly table
            in report_generator._qc_table() (M3d).
    Step 1: methylation β-value density
    Step 2: fragment length distribution
    Step 3: dinucleotide frequency
    """
    from visualization.plot_qc import (
        plot_methylation_distribution,
        plot_fragment_length,
        plot_dinucleotide_freq,
    )
    out_dir = getattr(args, "output_dir", "results/2_qc")
    step    = args.step

    # M3d: step 0 has no plot output — skip silently
    if step == 0:
        return

    if step == 1:
        sub_dir = os.path.join(out_dir, "1_methylation_distribution")
        os.makedirs(sub_dir, exist_ok=True)
        matrix = os.path.join(getattr(args, "matrices_dir", out_dir), "cpg_matrix.tsv")
        png, pdf = _out(sub_dir, stem="methylation_distribution")
        plot_methylation_distribution(matrix, png, pdf, args)

    elif step == 2:
        sub_dir = os.path.join(out_dir, "2_fragment_length")
        os.makedirs(sub_dir, exist_ok=True)
        prefix = os.path.join(sub_dir, "fragment_length")
        png, pdf = _out(sub_dir, stem="fragment_length")
        plot_fragment_length(prefix, png, pdf, args)

    elif step == 3:
        sub_dir = os.path.join(out_dir, "3_dinucleotide_freq")
        os.makedirs(sub_dir, exist_ok=True)
        prefix = os.path.join(sub_dir, "dinucleotide")
        png, pdf = _out(sub_dir, stem="dinucleotide_freq")
        plot_dinucleotide_freq(prefix, png, pdf, args)


# ── Differential ──────────────────────────────────────────────────────────────

def plot_differential(args):
    from visualization.plot_differential import plot_pca, plot_violin, plot_heatmap

    modality = args.modality
    out_dir  = os.path.join(args.output_dir, modality)
    coord_f  = os.path.join(out_dir, "pca_coordinates.txt")
    var_f    = os.path.join(out_dir, "pca_variance.txt")
    matrix_f = args.infile

    png, pdf = _out(out_dir, stem="pca")
    plot_pca(coord_f, var_f, png, pdf,
             feature_name=getattr(args, "feature_name", modality),
             colors=getattr(args, "colors", None))

    png, pdf = _out(out_dir, stem="violin")
    plot_violin(matrix_f, args.group_labels, png, pdf,
                feature_name=getattr(args, "feature_name", modality),
                colors=getattr(args, "colors", None))

    png, pdf = _out(out_dir, stem="heatmap")
    plot_heatmap(matrix_f, args.group_labels, png, pdf,
                 feature_name=getattr(args, "feature_name", modality),
                 colors=getattr(args, "colors", None),
                 top_n=getattr(args, "top_n", 500))


# ── DMR ───────────────────────────────────────────────────────────────────────

def plot_dmr(args):
    from visualization.plot_dmr import plot_dmr_volcano

    out_dir = getattr(args, "output_dir", "results/3_differential/dmr")
    ann_bed = os.path.join(out_dir, "dmr_annotated.bed")

    import pandas as pd
    df = pd.read_csv(ann_bed, sep="\t")

    png, pdf = _out(out_dir, stem="dmr_volcano")
    plot_dmr_volcano(
        df=df,
        group_a=args.group_a,
        group_b=args.group_b,
        png_path=png,
        pdf_path=pdf,
        top_n=getattr(args, "top_n", 20),
        q_thr=getattr(args, "q_thr", 0.05),
    )


# ── Fragmentomics ─────────────────────────────────────────────────────────────

def plot_fragmentomics(args, mode):
    from visualization.plot_fragmentomics import (
        plot_occupancy, plot_delfi, plot_end_motif, plot_cleavage, plot_wps
    )
    scope = getattr(args, "scope_info", {}) or {}
    scope_mode = scope.get("mode")

    if mode == "occupancy":
        import glob
        out_dir   = args.occ_out
        tsv_files = sorted(glob.glob(os.path.join(out_dir, "*.occupancy.tsv")))
        for tsv in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv))[0].replace(".occupancy", "")
            png, pdf = _out(out_dir, stem=f"{stem}_occupancy")
            plot_occupancy(tsv, png, pdf, scope_mode=scope_mode)

    elif mode == "delfi":
        import glob
        out_dir      = args.delfi_out
        group_labels = getattr(args, "group_labels", {})
        tsv_files    = sorted(glob.glob(os.path.join(out_dir, "*_delfi.tsv")))

        for tsv in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv))[0]
            png, pdf = _out(out_dir, stem=f"{stem}_genome")
            plot_delfi(tsv, png, pdf, scope_mode=scope_mode)

        if group_labels and len(tsv_files) > 1:
            from visualization.plot_fragmentomics import plot_delfi_group
            name_to_tsv = {
                os.path.splitext(os.path.basename(t))[0]
                  .replace("_delfi", "").replace(".markdup", ""): t
                for t in tsv_files
            }
            for grp, col_names in group_labels.items():
                grp_tsvs = [name_to_tsv[n] for n in col_names if n in name_to_tsv]
                if grp_tsvs:
                    png, pdf = _out(out_dir, stem=f"delfi_{grp}")
                    plot_delfi_group(grp_tsvs, png, pdf, label=grp)
            grp_means = {}
            for grp, col_names in group_labels.items():
                grp_tsvs = [name_to_tsv[n] for n in col_names if n in name_to_tsv]
                if grp_tsvs:
                    grp_means[grp] = grp_tsvs
            if len(grp_means) >= 2:
                from visualization.plot_fragmentomics import plot_delfi_comparison
                png, pdf = _out(out_dir, stem="delfi_comparison")
                plot_delfi_comparison(grp_means, png, pdf)

    elif mode == "end_motif":
        import glob
        out_dir      = args.end_motif_out
        group_labels = getattr(args, "group_labels", {})
        tsv_files    = sorted(glob.glob(os.path.join(out_dir, "*mer.tsv")))

        for tsv in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv))[0]
            png, pdf = _out(out_dir, stem=f"{stem}_top20")
            plot_end_motif(tsv, png, pdf, n=20)

        if group_labels and len(tsv_files) > 1:
            from visualization.plot_fragmentomics import plot_end_motif_group
            import re as _re
            name_to_tsv = {
                _re.sub(r"\.markdup|_\d+mer$",
                        "", os.path.splitext(os.path.basename(t))[0]): t
                for t in tsv_files
            }
            for grp, col_names in group_labels.items():
                grp_tsvs = [name_to_tsv[n] for n in col_names if n in name_to_tsv]
                if grp_tsvs:
                    png, pdf = _out(out_dir, stem=f"end_motif_{grp}")
                    plot_end_motif_group(grp_tsvs, png, pdf, label=grp, n=20)

    elif mode == "cleavage":
        import glob
        out_dir      = args.cleavage_out
        group_labels = getattr(args, "group_labels", {})
        bw_files     = sorted(glob.glob(os.path.join(out_dir, "*.bw")))
        bed          = getattr(args, "bed", None)
        up           = getattr(args, "upstream", 1500)
        dn           = getattr(args, "downstream", 1500)

        if bw_files and bed:
            name_to_bw = {
                os.path.splitext(os.path.basename(b))[0]
                  .replace("_cleavage", "").replace(".markdup", ""): b
                for b in bw_files
            }
            for grp, col_names in group_labels.items():
                grp_bws    = [name_to_bw[n] for n in col_names if n in name_to_bw]
                grp_labels = [n for n in col_names if n in name_to_bw]
                if grp_bws:
                    png, pdf = _out(out_dir, stem=f"cleavage_{grp}_samples")
                    plot_cleavage(grp_bws, bed, png, pdf,
                                  upstream=up, downstream=dn,
                                  labels=grp_labels)
            if len(group_labels) >= 2:
                from visualization.plot_fragmentomics import plot_cleavage_comparison
                grp_bw_dict = {}
                for grp, col_names in group_labels.items():
                    grp_bws = [name_to_bw[n] for n in col_names if n in name_to_bw]
                    if grp_bws:
                        grp_bw_dict[grp] = grp_bws
                if len(grp_bw_dict) >= 2:
                    png, pdf = _out(out_dir, stem="cleavage_comparison")
                    plot_cleavage_comparison(grp_bw_dict, bed, png, pdf,
                                             upstream=up, downstream=dn)

    elif mode == "wps":
        import glob
        out_dir   = args.wps_out
        tsv_files = sorted(glob.glob(os.path.join(out_dir, "*.wps.tsv")))
        for tsv in tsv_files:
            stem = os.path.splitext(os.path.basename(tsv))[0]
            png, pdf = _out(out_dir, stem=f"{stem}_profile")
            plot_wps(tsv, png, pdf, scope_mode=scope_mode)


# ── MESA ──────────────────────────────────────────────────────────────────────

def plot_mesa(args):
    from visualization.plot_mesa import plot_roc, plot_prob_heatmap, plot_spearman

    out_dir  = getattr(args, "output_dir", "results/5_mesa")
    pred_tsv = os.path.join(out_dir, "loocv_predictions.tsv")

    import pandas as pd
    pred_df = pd.read_csv(pred_tsv, sep="\t", index_col=0)

    png, pdf = _out(out_dir, stem="mesa_roc")
    plot_roc(pred_df, png, pdf)

    png, pdf = _out(out_dir, stem="mesa_heatmap")
    plot_prob_heatmap(pred_df, png, pdf)

    png, pdf = _out(out_dir, stem="mesa_spearman")
    plot_spearman(pred_df, png, pdf)


# ── Power ─────────────────────────────────────────────────────────────────────

def plot_power(args):
    from visualization.plot_qc import plot_power_curves

    out_dir = getattr(args, "output_dir", "results/0_power")
    tsv     = os.path.join(out_dir, "power_cumulative.tsv")

    import pandas as pd
    if not os.path.exists(tsv):
        return
    data = pd.read_csv(tsv, sep="\t", index_col=0)

    png, pdf = _out(out_dir, stem="power_cumulative")
    plot_power_curves(data, png, pdf,
                      threshold=getattr(args, "plot_threshold", 0.8))
