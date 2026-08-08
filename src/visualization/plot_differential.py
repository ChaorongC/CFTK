"""Differential analysis plots: PCA, violin, heatmap."""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy.stats import mannwhitneyu, gaussian_kde
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "figure.dpi": 150,
})

DEFAULT_COLORS = ["#2c7bb6", "#d7191c", "gray", "#1a9641", "#fdae61"]


def _save(fig, png_path, pdf_path):
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] saved → {png_path}")


def _palette(group_names, colors):
    colors = colors or []
    return {n: (colors[i] if i < len(colors) else DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
            for i, n in enumerate(group_names)}


def _resolve_groups(df_cols, group_labels_raw):
    """
    Resolve group → list of exact column names.
    Accepts:
      dict {"GroupA": ["col1","col2"]}  ← exact lists (primary form)
      dict {"GroupA": "prefix_"}        ← startswith prefix (legacy)
      list ["GroupA","prefA",...]        ← legacy pairs
    """
    if isinstance(group_labels_raw, dict):
        groups = group_labels_raw
    else:
        it     = iter(group_labels_raw)
        groups = {name: pref for name, pref in zip(it, it)}

    resolved = {}
    for name, spec in groups.items():
        if isinstance(spec, (list, tuple)):
            # exact column list — validate existence
            missing = [s for s in spec if s not in df_cols]
            if missing:
                import sys
                sys.exit(
                    f"[plot_diff] ERROR: columns not found for group '{name}': "
                    f"{missing}. Available: {list(df_cols)[:10]}"
                )
            resolved[name] = list(spec)
        elif isinstance(spec, str):
            # prefix matching (legacy)
            matched = [c for c in df_cols if c.startswith(spec)]
            if not matched:
                import sys
                sys.exit(
                    f"[plot_diff] ERROR: no columns start with '{spec}' "
                    f"for group '{name}'. Available: {list(df_cols)[:10]}"
                )
            resolved[name] = matched
    return resolved


def _load_matrix(filepath):
    sep = "," if filepath.endswith(".csv") else "\t"
    return pd.read_csv(filepath, sep=sep, header=0, index_col=0)


# ── PCA ───────────────────────────────────────────────────────────────────────

def plot_pca(coord_txt, var_txt, png_path, pdf_path,
             feature_name="", colors=None):
    """Plot PCA scatter from saved coordinate file."""
    coord = pd.read_csv(coord_txt, sep="\t", index_col=0)
    var   = pd.read_csv(var_txt, sep="\t")

    group_names = coord["group"].unique().tolist()
    palette     = _palette(group_names, colors)
    markers     = ["s", "o", "^", "D", "v"]

    sns.set_theme("paper", style="ticks")
    mpl.rcParams.update({"axes.spines.top": True, "axes.spines.right": True})
    fig, ax = plt.subplots(figsize=(3.5, 3.5))

    for i, name in enumerate(group_names):
        sub = coord[coord["group"] == name]
        ax.scatter(sub["PC1"], sub["PC2"],
                   c=palette[name], marker=markers[i % len(markers)],
                   s=14, alpha=0.75, label=name, linewidths=0)

    pc1_var = var.loc[var["PC"] == "PC1", "variance_explained_pct"].values[0]
    pc2_var = var.loc[var["PC"] == "PC2", "variance_explained_pct"].values[0]
    ax.set_xlabel(f"PC1 ({pc1_var:.2f}%)")
    ax.set_ylabel(f"PC2 ({pc2_var:.2f}%)")

    title = f"{feature_name} — PCA" if feature_name else "PCA"
    ax.set_title(title, pad=32)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01),
              ncol=len(group_names), frameon=False, fontsize=8)

    mpl.rcParams.update({"axes.spines.top": False, "axes.spines.right": False})
    fig.tight_layout()
    _save(fig, png_path, pdf_path)


# ── Violin ────────────────────────────────────────────────────────────────────

def _kde_violin(ax, vals, pos, color, width=0.6):
    vals = vals[~np.isnan(vals)]
    if len(vals) < 2:
        return
    kde    = gaussian_kde(vals)
    bw     = kde.factor * vals.std(ddof=1)
    y_pts  = np.linspace(vals.min() - 2.5 * bw, vals.max() + 2.5 * bw, 400)
    dens   = kde(y_pts)
    half_w = (dens / dens.max()) * (width / 2)
    face   = mpl.colors.to_rgba(color, alpha=0.55)
    edge   = mpl.colors.to_rgba(color, alpha=0.90)
    ax.fill_betweenx(y_pts, pos - half_w, pos + half_w,
                     facecolor=face, edgecolor=edge, linewidth=0.8)


def plot_violin(matrix_path, group_labels_raw, png_path, pdf_path,
                feature_name="", colors=None):
    df      = _load_matrix(matrix_path)
    groups  = _resolve_groups(df.columns, group_labels_raw)
    gnames  = list(groups.keys())
    palette = _palette(gnames, colors)

    # sample-level mean across all features
    sample_means = {}
    for name, cols in groups.items():
        sub = df[cols].copy()
        sub = sub.loc[~sub.isna().all(axis=1)]
        sample_means[name] = np.nanmean(sub.values, axis=0)

    sns.set_theme("paper", style="ticks")
    mpl.rcParams.update({"axes.spines.top": True, "axes.spines.right": True})
    fig, ax = plt.subplots(figsize=(2.8, 3.8))
    fig.subplots_adjust(bottom=0.22)

    positions = list(range(1, len(gnames) + 1))
    medians   = {}
    all_vals  = np.concatenate(list(sample_means.values()))
    y_min, y_max = np.nanmin(all_vals), np.nanmax(all_vals)

    for pos, name in zip(positions, gnames):
        vals = sample_means[name]
        _kde_violin(ax, vals, pos, palette[name])
        med = np.nanmedian(vals)
        medians[name] = med
        ax.hlines(med, pos - 0.13, pos + 0.13, colors="black", linewidths=1.4)

    if len(gnames) == 2:
        a, b = sample_means[gnames[0]], sample_means[gnames[1]]
        _, pv = mannwhitneyu(a, b, alternative="two-sided")
        ax.text((positions[0] + positions[-1]) / 2,
                y_min + (y_max - y_min) * 0.97,
                f"p = {pv:.2g}", ha="center", va="top", fontsize=8)

    ax.set_xticks(positions)
    ax.set_xticklabels(gnames)
    ax.set_ylabel(f"Average {feature_name} value" if feature_name
                  else "Average value of all features")
    ax.set_title(f"Metrics {feature_name} differential analysis" if feature_name
                 else "Metrics differential analysis", pad=8)

    # median labels below x-axis
    fig.canvas.draw()
    trans    = ax.transData + fig.transFigure.inverted()
    y_fig    = ax.get_position().y0 - 0.055
    for pos, name in zip(positions, gnames):
        xf, _ = trans.transform((pos, y_min))
        fig.text(xf, y_fig, f"{medians[name]:.4f}",
                 ha="center", va="top", fontsize=7,
                 transform=fig.transFigure)

    mpl.rcParams.update({"axes.spines.top": False, "axes.spines.right": False})
    _save(fig, png_path, pdf_path)


# ── Heatmap ───────────────────────────────────────────────────────────────────

def plot_heatmap(matrix_path, group_labels_raw, png_path, pdf_path,
                 feature_name="Feature", colors=None, top_n=500,
                 vmin=-5, vmax=5):
    import warnings

    df      = _load_matrix(matrix_path)
    groups  = _resolve_groups(df.columns, group_labels_raw)
    gnames  = list(groups.keys())
    palette = _palette(gnames, colors)

    ordered_cols = [c for name in gnames for c in groups[name]]
    labels       = [name for name in gnames for _ in groups[name]]

    sub = df[ordered_cols].copy()

    # ── NaN handling: same as original visualization_v1_2.py _prepare_data ──
    # Step 1: drop all-NaN features
    sub = sub.loc[~sub.isna().all(axis=1)]
    # Step 2: drop features with NaN > 50%
    sub = sub.loc[sub.isna().mean(axis=1) <= 0.5]
    # Step 3: drop zero-variance features
    variances = sub.var(axis=1, skipna=True)
    sub = sub.loc[~(variances.isna() | (variances == 0))]

    nan_pct = sub.isna().mean().mean() * 100
    if nan_pct > 0:
        print(f"[heatmap] Remaining NaN: {nan_pct:.1f}% → imputed with column mean")

    if sub.shape[0] < 2:
        print(f"[heatmap] WARNING: only {sub.shape[0]} features remain after "
              "NaN filtering — skipping heatmap.")
        return

    # Step 4: top-N by variance
    if top_n and top_n < sub.shape[0]:
        sub = sub.loc[sub.var(axis=1, skipna=True).nlargest(top_n).index]
        print(f"[heatmap] Top {top_n} features selected")

    # Step 5: impute remaining NaN with column mean, then z-score
    # (matches original: SimpleImputer → StandardScaler pipeline on samples axis)
    pipe = make_pipeline(SimpleImputer(strategy="mean"), StandardScaler())
    # WPS matrices may legitimately contain repeated interval labels when
    # overlapping windows are retained.  Passing ``sub.T`` as a DataFrame
    # makes recent scikit-learn versions validate those labels as feature
    # names and fail on duplicates.  The labels are only metadata here; use
    # the numeric array for preprocessing so every retained row is preserved.
    X    = pipe.fit_transform(sub.T.to_numpy(dtype=float)).T
    hm   = pd.DataFrame(X, index=sub.index, columns=sub.columns)

    col_colors = pd.Series(labels, index=ordered_cols, name="Group").map(palette)
    sns.set_theme("paper", style="ticks")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        try:
            from scipy.stats import SmallSampleWarning
            warnings.filterwarnings("ignore", category=SmallSampleWarning)
        except ImportError:
            pass  # SmallSampleWarning not available in scipy < 1.11

        cg = sns.clustermap(
            hm.clip(vmin, vmax),
            row_cluster=False, col_cluster=True,
            method="centroid", col_colors=col_colors,
            cmap="vlag", yticklabels=False, xticklabels=False,
            vmin=vmin, vmax=vmax,
            cbar_pos=(0.02, 0.60, 0.015, 0.16),
            figsize=(6, 4),
        )

    cg.ax_heatmap.set_xlabel("Sample")
    cg.ax_heatmap.set_ylabel(feature_name)

    cax = cg.ax_cbar
    cax.set_ylabel("Z score", rotation=90, va="center", fontsize=8)
    cax.yaxis.set_ticks_position("left")
    cax.tick_params(axis="y", pad=2, labelsize=7)

    handles = [Patch(facecolor=palette[n], edgecolor="none", label=n) for n in gnames]
    cg.fig.legend(handles=handles, title="Group",
                  loc="upper left", bbox_to_anchor=(0.01, 0.50),
                  bbox_transform=cg.fig.transFigure,
                  frameon=False, fontsize=8, title_fontsize=8)

    cg.fig.savefig(png_path, dpi=300, bbox_inches="tight")
    cg.fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(cg.fig)
    print(f"[plot] saved → {png_path}")
