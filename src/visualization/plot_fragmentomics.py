"""Fragmentomics plots: occupancy, DELFI, end-motif, cleavage profile, WPS."""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
})


def _save(fig, png_path, pdf_path):
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_frag] saved → {png_path}")


def _clean_header(df):
    df.columns = df.columns.str.replace("^#", "", regex=True).str.strip()
    return df


def _scope_suffix(scope_mode):
    return f" [{scope_mode} scope]" if scope_mode else ""


# ── Nucleosome occupancy ──────────────────────────────────────────────────────

def plot_occupancy(tsv_path, png_path, pdf_path, figsize=(14, 4), scope_mode=None):
    """Plot occupancy scores from ``bigWigAverageOverBed`` output.

    The UCSC utility emits six columns (region name, size, covered, sum,
    mean0, mean).  Older CFTK outputs used a nine-column chromosome table;
    both formats remain supported.  Target-panel region names are plotted in
    input order because they do not encode chromosome coordinates.
    """
    raw = pd.read_csv(tsv_path, sep="\t", header=None)
    if raw.empty:
        print(f"[plot_occ] empty occupancy table {tsv_path}, skipping.")
        return

    if raw.shape[1] == 6:
        df = raw.copy()
        df.columns = ["region", "size", "covered", "sum", "mean0", "mean"]
        df["mean"] = pd.to_numeric(df["mean"], errors="coerce")
        df = df.dropna(subset=["mean"]).reset_index(drop=True)
        if df.empty:
            print(f"[plot_occ] no numeric occupancy values in {tsv_path}, skipping.")
        else:
            x = np.arange(len(df))
            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(x, df["mean"], color="#e67e22", linewidth=0.6, alpha=0.8)
            ax.set_xlabel("Target interval (input order)")
            ax.set_ylabel("Occupancy Score")
            ax.set_title(f"{Path(tsv_path).stem} — Nucleosome Occupancy{_scope_suffix(scope_mode)}")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
            fig.tight_layout()
            _save(fig, png_path, pdf_path)
        return

    if raw.shape[1] != 9:
        print(f"[plot_occ] unexpected {raw.shape[1]} columns in {tsv_path}, skipping.")
        return

    df = raw.copy()
    df.columns = ["chrom", "start", "end", "name", "size",
                  "mean0", "mean", "max", "summit"]
    df["mean"] = pd.to_numeric(df["mean"], errors="coerce")
    df = df.dropna(subset=["mean"]).reset_index(drop=True)
    if df.empty:
        print(f"[plot_occ] no numeric occupancy values in {tsv_path}, skipping.")
        return

    chrom_order = [f"chr{i}" for i in range(1, 23)]
    recognized = set(df["chrom"].astype(str)) & set(chrom_order)
    if not recognized:
        x = np.arange(len(df))
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(x, df["mean"], color="#e67e22", linewidth=0.6, alpha=0.8)
        ax.set_xlabel("Target interval (input order)")
        ax.set_ylabel("Occupancy Score")
        ax.set_title(f"{Path(tsv_path).stem} — Nucleosome Occupancy{_scope_suffix(scope_mode)}")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
        fig.tight_layout()
        _save(fig, png_path, pdf_path)
        return

    df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
    df = df.sort_values(["chrom", "start"]).reset_index(drop=True)

    current_pos = 0
    chrom_pos, chrom_ctr = {}, {}
    gpos_list = []
    for chrom in chrom_order:
        sub = df[df["chrom"] == chrom]
        if sub.empty:
            continue
        chrom_pos[chrom] = current_pos
        gpos_list.append(pd.Series(current_pos + np.arange(len(sub)), index=sub.index))
        chrom_ctr[chrom] = current_pos + len(sub) / 2
        current_pos += len(sub)
    df["_gpos"] = pd.concat(gpos_list).reindex(df.index)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(df["_gpos"], df["mean"], color="#e67e22", linewidth=0.6, alpha=0.8)

    for i, chrom in enumerate(chrom_order):
        if chrom not in chrom_pos:
            continue
        pos = chrom_pos[chrom]
        n = len(df[df["chrom"] == chrom])
        if i % 2 == 0:
            ax.axvspan(pos, pos + n, alpha=0.05, color="gray", zorder=0)
        if i > 0:
            ax.axvline(pos, color="lightgray", linewidth=0.5, linestyle="--", alpha=0.5)

    ax.set_xticks([chrom_ctr[c] for c in chrom_order if c in chrom_ctr])
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order if c in chrom_ctr], fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Occupancy Score")
    ax.set_title(f"{Path(tsv_path).stem} — Nucleosome Occupancy{_scope_suffix(scope_mode)}")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── DELFI ─────────────────────────────────────────────────────────────────────

def plot_delfi(tsv_path, png_path, pdf_path,
               use_corrected=True, figsize=(20, 4), ylim=None, scope_mode=None):
    """DELFI score line plot for the supplied bins."""
    with open(tsv_path) as f:
        first = f.readline()
    df = pd.read_csv(tsv_path, sep="\t")
    if first.startswith("#"):
        df.columns = df.columns.str.replace("^#", "", regex=True).str.strip()

    ratio_col = "ratio_corrected" if (use_corrected and "ratio_corrected" in df.columns) \
                else "ratio"
    df = df.dropna(subset=[ratio_col])

    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df["contig"] = pd.Categorical(df["contig"], categories=chrom_order, ordered=True)
    df = df.sort_values(["contig", "start"])

    current_pos = 0
    chrom_positions, chrom_centers = {}, {}
    gpos_list = []
    for chrom in chrom_order:
        sub = df[df["contig"] == chrom]
        if sub.empty:
            continue
        chrom_positions[chrom] = current_pos
        gpos_list.append(
            pd.Series(current_pos + np.arange(len(sub)), index=sub.index)
        )
        chrom_centers[chrom]   = current_pos + len(sub) / 2
        current_pos           += len(sub)
    df["_gpos"] = pd.concat(gpos_list).reindex(df.index)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(df["_gpos"], df[ratio_col], color="#1f77b4", linewidth=0.8, alpha=0.8)

    for i, chrom in enumerate(chrom_order):
        if chrom not in chrom_positions:
            continue
        pos = chrom_positions[chrom]
        n   = len(df[df["contig"] == chrom])
        if i % 2 == 0:
            ax.axvspan(pos, pos + n, alpha=0.05, color="gray", zorder=0)
        if i > 0:
            ax.axvline(pos, color="lightgray", linewidth=0.5, linestyle="--", alpha=0.5)

    ax.set_xticks([chrom_centers[c] for c in chrom_order if c in chrom_centers])
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order if c in chrom_centers],
                       fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("DELFI Score")
    ax.set_title(f"{Path(tsv_path).stem} — DELFI Score{_scope_suffix(scope_mode)}")
    if ylim:
        ax.set_ylim(*ylim)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── End-motif ─────────────────────────────────────────────────────────────────

def plot_end_motif(tsv_path, png_path, pdf_path, n=20, show_values=True):
    """Horizontal bar chart of top N most frequent k-mers."""
    df   = pd.read_csv(tsv_path, sep="\t", header=None, names=["kmer", "frequency"])
    data = df.nlargest(n, "frequency")

    height = max(6, n * 0.35)
    fig, ax = plt.subplots(figsize=(10, height))
    ax.barh(range(len(data)), data["frequency"],
            color="steelblue", alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data["kmer"], fontfamily="monospace", fontsize=9)
    ax.set_xlabel("Frequency")
    ax.set_ylabel("K-mer")
    ax.set_title(f"{Path(tsv_path).stem} — Top {n} End Motifs")
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)

    if show_values:
        for i, val in enumerate(data["frequency"]):
            ax.text(val, i, f" {val:.4f}", va="center", fontsize=7)

    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── Cleavage profile ──────────────────────────────────────────────────────────

def plot_cleavage(bw_paths, bed_path, png_path, pdf_path,
                  upstream=1500, downstream=1500, labels=None,
                  smooth_window=10, colors=None):
    """Aggregate cleavage profile over BED regions from bigWig files."""
    try:
        import pyBigWig
    except ImportError:
        print("[plot_cleavage] pyBigWig not installed, skipping.")
        return

    default_colors = ["darkred", "darkblue", "darkgreen", "darkorange", "purple"]
    colors = colors or default_colors
    labels = labels or [Path(b).stem for b in bw_paths]
    window = upstream + downstream

    def _extract(bw_file):
        bw      = pyBigWig.open(bw_file)
        regions = pd.read_csv(bed_path, sep="\t", header=None,
                              names=["chrom", "start", "end", "name", "score", "strand"])
        mat = []
        for _, r in regions.iterrows():
            center = (r["start"] + r["end"]) // 2
            ws, we = max(0, center - upstream), center + downstream
            try:
                vals = bw.values(str(r["chrom"]), ws, we)
                if vals and len(vals) == window:
                    mat.append(np.nan_to_num(np.array(vals)))
            except Exception:
                pass
        bw.close()
        if not mat:
            return None
        sig = np.mean(mat, axis=0)
        if smooth_window > 1:
            sig = np.convolve(sig, np.ones(smooth_window) / smooth_window, mode="same")
        return sig

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(-upstream, downstream)

    for bw_path, label, color in zip(bw_paths, labels, colors):
        sig = _extract(bw_path)
        if sig is None:
            continue
        ax.plot(x, sig * 100, linewidth=2, label=label, color=color, alpha=0.8)

    ax.axvline(0, color="black", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.set_xlabel("Distance to CTCF motif (bp)")
    ax.set_ylabel("Cleavage Proportion (%)")
    ax.set_title("Cleavage Profile")
    ax.set_xticks([-1500, -750, 0, 750, 1500])
    ax.set_xticklabels(["-1.5kb", "-750", "0", "750", "1.5kb"])
    ax.legend(frameon=True, fontsize=9)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── WPS ───────────────────────────────────────────────────────────────────────

def plot_wps(tsv_path, png_path, pdf_path, figsize=(14, 4), scope_mode=None):
    """Mean WPS per supplied region."""
    df = pd.read_csv(tsv_path, sep="\t")
    if "mean_WPS" not in df.columns:
        print(f"[plot_wps] mean_WPS column not found in {tsv_path}, skipping.")
        return

    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df["chr"]   = pd.Categorical(df["chr"], categories=chrom_order, ordered=True)
    df          = df.sort_values(["chr", "start"]).reset_index(drop=True)

    current_pos  = 0
    chrom_pos, chrom_ctr = {}, {}
    gpos_list = []
    for chrom in chrom_order:
        sub = df[df["chr"] == chrom]
        if sub.empty:
            continue
        chrom_pos[chrom] = current_pos
        gpos_list.append(
            pd.Series(current_pos + np.arange(len(sub)), index=sub.index)
        )
        chrom_ctr[chrom]  = current_pos + len(sub) / 2
        current_pos      += len(sub)
    df["_gpos"] = pd.concat(gpos_list).reindex(df.index)

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(df["_gpos"], df["mean_WPS"], color="#2c7bb6", linewidth=0.6, alpha=0.8)

    for i, chrom in enumerate(chrom_order):
        if chrom not in chrom_pos:
            continue
        pos = chrom_pos[chrom]
        n   = len(df[df["chr"] == chrom])
        if i % 2 == 0:
            ax.axvspan(pos, pos + n, alpha=0.05, color="gray", zorder=0)
        if i > 0:
            ax.axvline(pos, color="lightgray", linewidth=0.5, linestyle="--", alpha=0.5)

    ax.set_xticks([chrom_ctr[c] for c in chrom_order if c in chrom_ctr])
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order if c in chrom_ctr],
                       fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Mean WPS")
    ax.set_title(f"{Path(tsv_path).stem} — WPS{_scope_suffix(scope_mode)}")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── DELFI group / comparison plots ───────────────────────────────────────────

def _load_delfi_ratio(tsv_path, use_corrected=True):
    """Load DELFI ratio series sorted by genome position."""
    with open(tsv_path) as f:
        first = f.readline()
    df = pd.read_csv(tsv_path, sep="\t")
    if first.startswith("#"):
        df.columns = df.columns.str.replace("^#", "", regex=True).str.strip()
    ratio_col = "ratio_corrected" if (use_corrected and "ratio_corrected" in df.columns)                 else "ratio"
    df = df.dropna(subset=[ratio_col])
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    df["contig"] = pd.Categorical(df["contig"], categories=chrom_order, ordered=True)
    df = df.sort_values(["contig", "start"]).reset_index(drop=True)
    # assign genome position
    current_pos = 0
    chrom_positions, chrom_centers, gpos_list = {}, {}, []
    for chrom in chrom_order:
        sub = df[df["contig"] == chrom]
        if sub.empty:
            continue
        chrom_positions[chrom] = current_pos
        gpos_list.append(pd.Series(current_pos + np.arange(len(sub)), index=sub.index))
        chrom_centers[chrom] = current_pos + len(sub) / 2
        current_pos += len(sub)
    if gpos_list:
        df["_gpos"] = pd.concat(gpos_list).reindex(df.index)
    return df, ratio_col, chrom_order, chrom_positions, chrom_centers


def _delfi_ax_style(ax, df, chrom_order, chrom_positions):
    """Apply chromosome background, grid and labels to DELFI axes."""
    chrom_centers = {}
    for i, chrom in enumerate(chrom_order):
        if chrom not in chrom_positions:
            continue
        sub = df[df["contig"] == chrom]
        n   = len(sub)
        pos = chrom_positions[chrom]
        chrom_centers[chrom] = pos + n / 2
        if i % 2 == 0:
            ax.axvspan(pos, pos + n, alpha=0.05, color="gray", zorder=0)
        if i > 0:
            ax.axvline(pos, color="lightgray", linewidth=0.5, linestyle="--", alpha=0.5)
    ax.set_xticks([chrom_centers[c] for c in chrom_order if c in chrom_centers])
    ax.set_xticklabels([c.replace("chr", "") for c in chrom_order if c in chrom_centers],
                       fontsize=7)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("DELFI Score")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.5)


def _align_delfi_ratios(tsv_paths, use_corrected=True):
    """Inner-join sample ratios by genomic coordinates for cohort plotting."""

    coordinate_columns = ["contig", "start", "stop"]
    aligned = None
    value_columns = []
    for index, tsv in enumerate(tsv_paths):
        df, ratio_col, _chrom_order, _chrom_pos, _chrom_centers = _load_delfi_ratio(
            tsv, use_corrected
        )
        missing = [column for column in coordinate_columns if column not in df.columns]
        if missing:
            raise ValueError(
                f"DELFI table {tsv} is missing coordinate column(s): {', '.join(missing)}"
            )
        if df.duplicated(coordinate_columns).any():
            raise ValueError(f"DELFI table {tsv} contains duplicate genomic bins")
        value_column = f"sample_{index}"
        frame = df[coordinate_columns + [ratio_col]].copy()
        frame["contig"] = frame["contig"].astype("object")
        frame = frame.rename(columns={ratio_col: value_column})
        aligned = frame if aligned is None else aligned.merge(
            frame,
            on=coordinate_columns,
            how="inner",
            sort=False,
            validate="one_to_one",
        )
        value_columns.append(value_column)

    if aligned is None or aligned.empty:
        raise ValueError("DELFI samples have no shared genomic bins with valid ratios")

    chrom_order = [f"chr{i}" for i in range(1, 23)]
    aligned["contig"] = pd.Categorical(
        aligned["contig"], categories=chrom_order, ordered=True
    )
    aligned = aligned.sort_values(["contig", "start", "stop"]).reset_index(drop=True)
    current_pos = 0
    chrom_positions = {}
    gpos_list = []
    for chrom in chrom_order:
        sub = aligned[aligned["contig"] == chrom]
        if sub.empty:
            continue
        chrom_positions[chrom] = current_pos
        gpos_list.append(pd.Series(current_pos + np.arange(len(sub)), index=sub.index))
        current_pos += len(sub)
    aligned["_gpos"] = pd.concat(gpos_list).reindex(aligned.index)
    return aligned, value_columns, chrom_order, chrom_positions


def plot_delfi_group(tsv_paths, png_path, pdf_path, label="Group",
                     use_corrected=True, figsize=(20, 4), color="#1f77b4"):
    """Mean DELFI score across samples in one group."""
    ref_df, value_columns, chrom_order, chrom_positions = _align_delfi_ratios(
        tsv_paths, use_corrected
    )
    mean_ratio = ref_df[value_columns].mean(axis=1).to_numpy()

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(ref_df["_gpos"], mean_ratio, color=color, linewidth=0.8, alpha=0.8,
            label=label)
    _delfi_ax_style(ax, ref_df, chrom_order, chrom_positions)
    ax.set_title(f"DELFI Score — {label} (mean)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


def plot_delfi_comparison(grp_tsv_dict, png_path, pdf_path,
                          use_corrected=True, figsize=(20, 4)):
    """Two-group comparison: one mean line per group."""
    colors = ["#FF8C00", "#1f77b4", "#2ca02c", "#d62728"]
    fig, ax = plt.subplots(figsize=figsize)
    paths = [path for tsv_paths in grp_tsv_dict.values() for path in tsv_paths]
    ref_df, value_columns, chrom_order, chrom_positions = _align_delfi_ratios(
        paths, use_corrected
    )
    offset = 0
    for i, (grp, tsv_paths) in enumerate(grp_tsv_dict.items()):
        group_columns = value_columns[offset:offset + len(tsv_paths)]
        offset += len(tsv_paths)
        mean_ratio = ref_df[group_columns].mean(axis=1).to_numpy()
        ax.plot(ref_df["_gpos"], mean_ratio, color=colors[i % len(colors)],
                linewidth=0.8, alpha=0.8, label=grp)

    _delfi_ax_style(ax, ref_df, chrom_order, chrom_positions)
    ax.set_title("DELFI Score — Group Comparison")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── End-motif group plot ──────────────────────────────────────────────────────

def plot_end_motif_group(tsv_paths, png_path, pdf_path, label="Group",
                         n=20, color="steelblue"):
    """Mean end-motif frequency barplot for one group."""
    all_dfs = []
    for tsv in tsv_paths:
        df = pd.read_csv(tsv, sep="\t", header=None, names=["kmer", "frequency"])
        all_dfs.append(df.set_index("kmer")["frequency"])
    mean_freq = pd.concat(all_dfs, axis=1).mean(axis=1)
    data = mean_freq.nlargest(n).reset_index()
    data.columns = ["kmer", "frequency"]

    height = max(6, n * 0.35)
    fig, ax = plt.subplots(figsize=(10, height))
    ax.barh(range(len(data)), data["frequency"],
            color=color, alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(data["kmer"], fontfamily="monospace", fontsize=9)
    ax.set_xlabel("Mean Frequency")
    ax.set_ylabel("K-mer")
    ax.set_title(f"End Motif — {label} (mean top {n})")
    ax.grid(axis="x", alpha=0.3, linestyle="--")
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    for i, val in enumerate(data["frequency"]):
        ax.text(val, i, f" {val:.4f}", va="center", fontsize=7)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── Cleavage comparison plot ──────────────────────────────────────────────────

def plot_cleavage_comparison(grp_bw_dict, bed_path, png_path, pdf_path,
                              upstream=1500, downstream=1500,
                              smooth_window=10):
    """Two-group comparison: one mean cleavage line per group."""
    try:
        import pyBigWig
    except ImportError:
        print("[plot_cleavage] pyBigWig not installed, skipping.")
        return

    colors = ["#2c7bb6", "#d7191c", "#1a9641", "#fdae61"]
    window = upstream + downstream

    def _mean_signal(bw_paths):
        regions = pd.read_csv(bed_path, sep="\t", header=None,
                              names=["chrom", "start", "end", "name", "score", "strand"])
        all_sigs = []
        for bw_file in bw_paths:
            bw  = pyBigWig.open(bw_file)
            mat = []
            for _, r in regions.iterrows():
                center = (r["start"] + r["end"]) // 2
                ws, we = max(0, center - upstream), center + downstream
                try:
                    vals = bw.values(str(r["chrom"]), ws, we)
                    if vals and len(vals) == window:
                        mat.append(np.nan_to_num(np.array(vals)))
                except Exception:
                    pass
            bw.close()
            if mat:
                all_sigs.append(np.mean(mat, axis=0))
        if not all_sigs:
            return None
        sig = np.mean(all_sigs, axis=0)
        if smooth_window > 1:
            sig = np.convolve(sig, np.ones(smooth_window) / smooth_window, mode="same")
        return sig

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(-upstream, downstream)

    for i, (grp, bw_paths) in enumerate(grp_bw_dict.items()):
        sig = _mean_signal(bw_paths)
        if sig is None:
            continue
        ax.plot(x, sig * 100, linewidth=2, label=grp,
                color=colors[i % len(colors)], alpha=0.8)

    ax.axvline(0, color="black", linestyle="--", linewidth=1.5, alpha=0.7)
    ax.set_xlabel("Distance to CTCF motif (bp)")
    ax.set_ylabel("Cleavage Proportion (%)")
    ax.set_title("Cleavage Profile — Group Comparison")
    xtick_vals = [-upstream, -upstream // 2, 0, upstream // 2, upstream]
    ax.set_xticks(xtick_vals)
    ax.set_xticklabels([f"{v:+d}" if v != 0 else "0" for v in xtick_vals])
    ax.legend(frameon=True, fontsize=9)
    ax.grid(axis="y", alpha=0.3, linestyle="--")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _save(fig, png_path, pdf_path)
