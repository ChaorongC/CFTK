"""
data_builder.py — compute interactive Plotly figure data for the HTML report.

Called from generate_report() before HTML assembly. Each function:
  1. Reads source data files from results/
  2. Builds a Plotly figure dict
  3. Saves it as JSON to data/{key}.txt
  4. Returns the figure dict (embedded as JS variable in report.html)

Data embedding strategy:
  All figure dicts are injected as window.__CFTK_DATA__ = {...} in a <script>
  block so they work with file:// protocol (no fetch() needed).

Sections handled here (interactive charts):
  2.2  Fragment Length   — line chart, all samples + group means
  2.3  Dinucleotide      — AT vs GC line chart
  3.1  PCA               — scatter per modality (reads pca_coordinates.txt)
  3.4  DMR Volcano       — scatter with gene labels
  4.2  End Motif         — per-group bar chart (keys: end_motif_{group})
  5.2  MESA ROC          — line chart per modality
  5.3  MESA Heatmap      — prediction probability heatmap
  5.4  MESA Spearman     — correlation heatmap across modalities

NOTE: the Sample Correlation heatmap (build_correlation) is retained for
      reference but is no longer built by build_all() — section 2.5 was
      removed from the report.
"""

from __future__ import annotations

import glob
import json
import os

import numpy as np
import pandas as pd

# ── Plotly colour palette (Tableau-10) ────────────────────────────────────────
_COLORS = [
    "#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac",
]

_GROUP_COLORS = ["#1e5fa0", "#c0392b", "#1a9641", "#9b59b6", "#e67e22"]


# ── helpers ───────────────────────────────────────────────────────────────────

def _source_roots(source_results: dict[str, str] | None) -> list[str]:
    """Return unique source result roots in configured sample order."""
    roots = []
    for root in (source_results or {}).values():
        root = str(root)
        if root and root not in roots:
            roots.append(root)
    return roots


def _fragment_raw_files(
    rd: str,
    source_results: dict[str, str] | None,
    sample_names: list[str],
) -> list[str]:
    """Find one fragment-length CSV per configured sample, preferring ``rd``."""
    selected: dict[str, str] = {}
    wanted = set(sample_names)
    for root in [rd] + _source_roots(source_results):
        frag_dir = os.path.join(root, "2_qc", "2_fragment_length")
        for path in sorted(glob.glob(os.path.join(frag_dir, "fragment_length.*.raw.csv"))):
            name = (os.path.basename(path).replace("fragment_length.", "")
                    .replace(".raw.csv", "").replace(".markdup", ""))
            if (not wanted or name in wanted) and name not in selected:
                selected[name] = path
    return [selected[name] for name in sample_names if name in selected] + [
        path for name, path in selected.items() if name not in sample_names
    ]

def _layout(title="", xaxis_title="", yaxis_title="", height=380, **kw) -> dict:
    base = {
        "title":  {"text": title, "font": {"size": 13, "color": "#000000"}},
        "height": height,
        "margin": {"l": 60, "r": 20, "t": 44, "b": 80},
        "legend": {"orientation": "h", "y": -0.28, "x": 0,
                   "font": {"color": "#000000"}},
        "xaxis":  {"title": {"text": xaxis_title, "font": {"color": "#000000"}},
                   "automargin": True, "color": "#000000",
                   "tickangle": -40, "tickfont": {"size": 10, "color": "#000000"}},
        "yaxis":  {"title": {"text": yaxis_title, "font": {"color": "#000000"}},
                   "automargin": True, "color": "#000000",
                   "tickfont": {"color": "#000000"}},
        "plot_bgcolor":  "#ffffff",
        "paper_bgcolor": "#ffffff",
        "font": {"family": "Times New Roman, Times, serif", "size": 11,
                 "color": "#000000"},
        "hoverlabel": {"bgcolor": "white", "font": {"size": 11}},
    }
    base.update(kw)
    return base


def _save(data_dir: str, key: str, fig: dict) -> dict:
    """Save figure dict to data/{key}.txt and return it."""
    os.makedirs(data_dir, exist_ok=True)
    path = os.path.join(data_dir, f"{key}.txt")
    with open(path, "w") as f:
        json.dump(fig, f, separators=(",", ":"))
    print(f"[data_builder] saved → {path}")
    return fig


def _missing_fig(msg: str, key: str) -> dict:
    return {
        "data": [],
        "layout": _layout(title=f"⚠ {msg}", height=120),
        "_missing": True,
        "_key": key,
    }


# ── 2.2 Fragment Length ────────────────────────────────────────────────────────

def build_fragment_length(
    rd: str,
    data_dir: str,
    group_labels: dict,
    source_results: dict[str, str] | None = None,
) -> dict:
    """
    Read all fragment_length.*.raw.csv files from 2_qc/2_fragment_length/.
    Produce a line chart with:
      - one trace per sample (initially hidden, toggled by dropdown)
      - group-mean traces (visible by default)
      - overall mean trace
    """
    sample_names = [sample for members in group_labels.values() for sample in members]
    raw_files = _fragment_raw_files(rd, source_results, sample_names)
    if not raw_files:
        return _save(data_dir, "fragment_length",
                     _missing_fig("Fragment length data not found — run cftk qc -s 2", "fragment_length"))

    # Load all samples
    sample_data: dict[str, pd.Series] = {}
    size_arr = np.arange(500)

    for fp in raw_files:
        stem = os.path.basename(fp)
        # strip prefix "fragment_length." and suffix ".raw.csv"
        name = stem.replace("fragment_length.", "").replace(".raw.csv", "")
        # also strip ".markdup" suffix common in bam names
        name = name.replace(".markdup", "")
        try:
            t = pd.read_table(fp, skiprows=1).iloc[:, :2]
            t.columns = ["Size", "Occurrences"]
            t = t[t["Size"] < 500]
            base = pd.DataFrame({"Size": np.arange(500)})
            merged = pd.merge(base, t, on="Size", how="left").sort_values("Size")["Occurrences"].fillna(0)
            vals = merged.values.astype(float)
            total = vals.sum()
            ratio = pd.Series(100 * vals / total if total > 0 else vals, index=size_arr)
            sample_data[name] = ratio
        except Exception:
            pass

    if not sample_data:
        return _save(data_dir, "fragment_length",
                     _missing_fig("Could not parse fragment length CSV files", "fragment_length"))

    x_vals = list(map(int, size_arr))
    traces = []

    # Per-sample traces (hidden by default, shown only when user toggles)
    for i, (sname, ratio) in enumerate(sample_data.items()):
        s_peak = int(size_arr[ratio.values.argmax()])
        traces.append({
            "type": "scatter", "mode": "lines",
            "name": f"{sname} (peak={s_peak}bp)",
            "x": x_vals, "y": [round(v, 4) for v in ratio.values],
            "line": {"color": _COLORS[i % len(_COLORS)], "width": 1.2},
            "opacity": 0.5,
            "visible": "legendonly",
            "legendgroup": "samples",
            "hovertemplate": f"{sname}<br>Length: %{{x}} bp<br>%{{y:.3f}}%<extra></extra>",
        })

    # Group mean traces
    all_groups = {s: g for g, members in group_labels.items() for s in members}
    for gi, (grp, members) in enumerate(group_labels.items()):
        grp_ratios = [v for k, v in sample_data.items() if k in members]
        if not grp_ratios:
            continue
        mean_ratio = pd.concat(grp_ratios, axis=1).mean(axis=1)
        peak = int(size_arr[mean_ratio.values.argmax()])
        traces.append({
            "type": "scatter", "mode": "lines",
            "name": f"{grp} (mean, peak={peak}bp)",
            "x": x_vals, "y": [round(v, 4) for v in mean_ratio.values],
            "line": {"color": _GROUP_COLORS[gi % len(_GROUP_COLORS)], "width": 2.5},
            "hovertemplate": f"{grp} mean<br>%{{x}} bp<br>%{{y:.3f}}%<extra></extra>",
        })

    # Overall mean
    all_ratio = pd.concat(list(sample_data.values()), axis=1).mean(axis=1)
    overall_peak = int(size_arr[all_ratio.values.argmax()])
    traces.append({
        "type": "scatter", "mode": "lines",
        "name": f"All samples (mean, peak={overall_peak}bp)",
        "x": x_vals, "y": [round(v, 4) for v in all_ratio.values],
        "line": {"color": "#555", "width": 2, "dash": "dot"},
        "hovertemplate": f"Overall mean<br>%{{x}} bp<br>%{{y:.3f}}%<extra></extra>",
    })

    # Build per-trace peak metadata for JS-driven dynamic peak line
    # Each trace carries customdata[0] = peak position
    for t in traces:
        name = t["name"]
        # Find peak from name string e.g. "Control (mean, peak=161bp)"
        import re as _re
        m = _re.search(r"peak=(\d+)bp", name)
        peak_val = int(m.group(1)) if m else overall_peak
        t["customdata"] = [[peak_val]] * len(t["x"])

    fig = {
        "data": traces,
        "layout": {
            **_layout(
                title="Fragment Length Distribution",
                xaxis_title="Fragment length (bp)",
                yaxis_title="% fragments",
                height=420,
            ),
            "xaxis": {"title": "Fragment length (bp)", "range": [50, 250],
                      "automargin": True, "tickfont": {"size": 10}},
            "yaxis": {"title": "% fragments", "automargin": True},
            "margin": {"l": 60, "r": 20, "t": 44, "b": 100},
        }
    }
    return _save(data_dir, "fragment_length", fig)


# ── 2.3 Dinucleotide Frequency ────────────────────────────────────────────────

def build_dinucleotide(
    rd: str,
    data_dir: str,
    frag_len: int = 167,
    source_results: dict[str, str] | None = None,
) -> dict:
    """
    Read dinucleotide.all_fragment_{pattern}.txt files.
    Column layout (bedtools nuc output):
      col index 3 = 4_usercol  (position -125 to 124)
      col index 14 = 15_user_patt_count (raw count)

    Aggregates AA/AT/TA/TT → AT-rich group, GG/GC/CG/CC → GC-rich group.
    """
    roots = [rd] + _source_roots(source_results)

    groups_map = {
        "AA/AT/TA/TT": ["AA", "AT", "TA", "TT"],
        "GG/GC/CG/CC": ["GG", "GC", "CG", "CC"],
    }

    combined = {}
    for label, patterns in groups_map.items():
        frames = []
        for pat in patterns:
            for root in roots:
                fp = os.path.join(
                    root, "2_qc", "3_dinucleotide_freq",
                    f"dinucleotide.all_fragment_{pat}.txt",
                )
                if not os.path.exists(fp):
                    continue
                try:
                    # bedtools nuc column names depend on input BED column count.
                    # pos col = *_usercol with integer values in [-200, 200]
                    # count col = last *_user_patt_count
                    header = pd.read_table(fp, nrows=5)
                    cols = list(header.columns)
                    usercols = [c for c in cols if c.endswith("_usercol")]
                    count_cols = [c for c in cols if c.endswith("_user_patt_count")]
                    if not usercols or not count_cols:
                        continue
                    # Identify position column: values are integers in [-200, 200]
                    pos_col = None
                    for uc in usercols:
                        sample_vals = pd.to_numeric(header[uc], errors="coerce").dropna()
                        if len(sample_vals) > 0:
                            vmin, vmax = sample_vals.min(), sample_vals.max()
                            if -200 <= vmin and vmax <= 200:
                                pos_col = uc
                                break
                    if pos_col is None:
                        pos_col = usercols[0]
                    count_col = count_cols[-1]
                    table = (
                        pd.read_table(fp, usecols=[pos_col, count_col])
                        .apply(pd.to_numeric, errors="coerce")
                        .dropna()
                    )
                    frames.append(table.groupby(pos_col)[count_col].sum())
                except Exception:
                    pass
        if frames:
            combined[label] = pd.concat(frames, axis=1).sum(axis=1).sort_index()

    if not combined:
        return _save(data_dir, "dinucleotide",
                     _missing_fig(
                         "Dinucleotide QC was not produced in the available results; run cftk qc -s 3 to generate it",
                         "dinucleotide",
                     ))

    # Normalise to %
    total = sum(s.sum() for s in combined.values())
    traces = []
    colors = {"AA/AT/TA/TT": "#4e79a7", "GG/GC/CG/CC": "#e15759"}

    for label, series in combined.items():
        pct = 100 * series / total if total > 0 else series
        traces.append({
            "type": "scatter", "mode": "lines", "name": label,
            "x": [int(v) for v in pct.index.tolist()],
            "y": [round(v, 5) for v in pct.values.tolist()],
            "line": {"color": colors.get(label, "#555"), "width": 2},
            "hovertemplate": f"{label}<br>Pos: %{{x}}<br>%{{y:.4f}}%<extra></extra>",
        })

    fig = {
        "data": traces,
        "layout": _layout(
            title="Dinucleotide Frequency",
            xaxis_title=f"Position relative to {frag_len}bp fragment centre (bp)",
            yaxis_title="Dinucleotide fraction (%)",
            height=400,
        ),
    }
    return _save(data_dir, "dinucleotide", fig)


# ── 3.1 PCA ──────────────────────────────────────────────────────────────────

def build_pca(rd: str, data_dir: str, group_labels: dict) -> dict:
    """
    Build PCA scatter plots per modality.

    pca_coordinates.txt format (confirmed):
      - index column (col 0 in file, unlabelled) = sample name
      - first named column = "group"
      - remaining columns  = PC1, PC2, PC3, ...

    pca_variance.txt: try multiple formats
      A) one float per line (explained variance ratio per PC, in order)
      B) tab-sep with index=PC and a numeric column
    """
    diff_base = os.path.join(rd, "3_differential")
    if not os.path.exists(diff_base):
        return _save(data_dir, "pca",
                     _missing_fig("PCA: no differential results found", "pca"))

    modalities = sorted([
        d for d in os.listdir(diff_base)
        if os.path.isdir(os.path.join(diff_base, d)) and d != "dmr"
    ])
    if not modalities:
        return _save(data_dir, "pca",
                     _missing_fig("PCA: no modality subdirectories in 3_differential/", "pca"))

    all_traces: list[dict] = []
    trace_offsets: dict[str, tuple[int, int]] = {}
    traces_per_mod: dict[str, list[dict]] = {}

    for mod in modalities:
        mod_dir    = os.path.join(diff_base, mod)
        coord_file = os.path.join(mod_dir, "pca_coordinates.txt")
        var_file   = os.path.join(mod_dir, "pca_variance.txt")

        if not os.path.exists(coord_file):
            continue

        # ── load coordinates ──────────────────────────────────────────────────
        try:
            # index_col=0 → sample names as index
            # header row: group  PC1  PC2  PC3 ...
            coords = pd.read_csv(coord_file, sep="\t", index_col=0)
            # 'group' is the first column; PC columns start from col index 1
            group_col = coords.columns[0]   # "group"
            pc_cols   = [c for c in coords.columns[1:] if c.upper().startswith("PC")]
            if len(pc_cols) < 2:
                # fallback: any column that isn't "group"
                pc_cols = [c for c in coords.columns if c != group_col]
            pc1_col = pc_cols[0]
            pc2_col = pc_cols[1]
        except Exception as e:
            print(f"[data_builder] PCA {mod}: could not parse coordinates — {e}")
            continue

        # ── load variance ─────────────────────────────────────────────────────
        var1 = var2 = None
        if os.path.exists(var_file):
            try:
                vdf = pd.read_csv(var_file, sep="\t", header=None)
                # Try format A: single column of floats
                col0 = pd.to_numeric(vdf.iloc[:, 0], errors="coerce").dropna()
                if len(col0) >= 2:
                    var1, var2 = float(col0.iloc[0]), float(col0.iloc[1])
                elif vdf.shape[1] >= 2:
                    # Format B: index=PC, value column
                    col1 = pd.to_numeric(vdf.iloc[:, 1], errors="coerce").dropna()
                    if len(col1) >= 2:
                        var1, var2 = float(col1.iloc[0]), float(col1.iloc[1])
            except Exception:
                pass

        # If variance looks like it's in fraction form (0–1) convert to %
        if var1 is not None and var1 <= 1.0:
            xlab = f"PC1 ({var1*100:.1f}%)"
            ylab = f"PC2 ({var2*100:.1f}%)" if var2 else "PC2"
        elif var1 is not None:
            xlab = f"PC1 ({var1:.1f}%)"
            ylab = f"PC2 ({var2:.1f}%)" if var2 else "PC2"
        else:
            xlab, ylab = "PC1", "PC2"

        # ── build traces — one per group ──────────────────────────────────────
        # Use the 'group' column in the file to assign colours
        # Also cross-reference with group_labels to colour-code consistently
        mod_traces: list[dict] = []

        # Build a group→colour index from group_labels (preserving order)
        group_color_idx = {g: i for i, g in enumerate(group_labels.keys())}

        # Groups present in the file
        file_groups = coords[group_col].unique().tolist() if group_col in coords.columns else []

        plotted_groups = set()

        for gi, (grp, members) in enumerate(group_labels.items()):
            # Match rows either by group_labels membership or by the file's group column
            mask_by_group_col = coords[group_col] == grp if group_col in coords.columns else pd.Series(False, index=coords.index)
            mask_by_members   = coords.index.isin(members)
            grp_coords = coords[mask_by_group_col | mask_by_members]

            if grp_coords.empty:
                continue

            x_vals = grp_coords[pc1_col].tolist()
            y_vals = grp_coords[pc2_col].tolist()
            labels = grp_coords.index.tolist()
            plotted_groups.add(grp)

            mod_traces.append({
                "type": "scatter", "mode": "markers+text",
                "name": grp,
                "x": x_vals, "y": y_vals,
                "text": labels, "textposition": "top center",
                "textfont": {"size": 9},
                "marker": {
                    "color": _GROUP_COLORS[gi % len(_GROUP_COLORS)],
                    "size": 10, "opacity": 0.85,
                    "line": {"color": "white", "width": 1},
                },
                "hovertemplate": (
                    "<b>%{text}</b><br>"
                    f"{xlab}: %{{x:.4f}}<br>{ylab}: %{{y:.4f}}"
                    "<extra></extra>"
                ),
                "_xlab": xlab,
                "_ylab": ylab,
            })

        # Any groups in the file not covered by group_labels
        for fg in file_groups:
            if fg in plotted_groups:
                continue
            grp_coords = coords[coords[group_col] == fg]
            if grp_coords.empty:
                continue
            gi = len(group_labels) + len(plotted_groups)
            mod_traces.append({
                "type": "scatter", "mode": "markers+text",
                "name": fg,
                "x": grp_coords[pc1_col].tolist(),
                "y": grp_coords[pc2_col].tolist(),
                "text": grp_coords.index.tolist(),
                "textposition": "top center",
                "textfont": {"size": 9},
                "marker": {
                    "color": _GROUP_COLORS[gi % len(_GROUP_COLORS)],
                    "size": 10, "opacity": 0.85,
                    "line": {"color": "white", "width": 1},
                },
                "hovertemplate": (
                    "<b>%{text}</b><br>"
                    f"{xlab}: %{{x:.4f}}<br>{ylab}: %{{y:.4f}}"
                    "<extra></extra>"
                ),
                "_xlab": xlab,
                "_ylab": ylab,
            })
            plotted_groups.add(fg)

        if not mod_traces:
            continue

        traces_per_mod[mod] = mod_traces
        trace_offsets[mod] = (len(all_traces), len(all_traces) + len(mod_traces))
        all_traces.extend(mod_traces)

    if not all_traces:
        return _save(data_dir, "pca",
                     _missing_fig("PCA: pca_coordinates.txt not found in any modality", "pca"))

    n_total = len(all_traces)

    # Build per-modality menu options (consumed by the HTML <select> renderer)
    options: list[dict] = []
    for mod, traces in traces_per_mod.items():
        vis = [False] * n_total
        s, e = trace_offsets[mod]
        for i in range(s, e):
            vis[i] = True
        first = traces[0]
        xlab  = first.get("_xlab", "PC1")
        ylab  = first.get("_ylab", "PC2")
        options.append({
            "label": mod.upper(),
            "visible": vis,
            "relayout": {"xaxis.title.text": xlab, "yaxis.title.text": ylab},
        })

    # Remove internal helper keys before serialisation
    for t in all_traces:
        t.pop("_xlab", None)
        t.pop("_ylab", None)

    # Initial visibility: first modality
    init_vis = [False] * n_total
    if traces_per_mod:
        first_mod = list(traces_per_mod.keys())[0]
        s, e = trace_offsets[first_mod]
        for i in range(s, e):
            init_vis[i] = True
    for i, t in enumerate(all_traces):
        t["visible"] = init_vis[i]

    # Axis labels for initial modality
    init_xlab = options[0]["relayout"]["xaxis.title.text"] if options else "PC1"
    init_ylab = options[0]["relayout"]["yaxis.title.text"] if options else "PC2"

    fig = {
        "data": all_traces,
        "layout": {
            **_layout(title="PCA", xaxis_title=init_xlab, yaxis_title=init_ylab, height=460),
        },
        # consumed by report_generator._plotly_select_render (HTML <select>)
        "_modality_menu": {"options": options},
    }
    return _save(data_dir, "pca", fig)


# ── 3.1b Correlation Heatmap (retained; no longer built by build_all) ─────────

def _find_matrix(rd: str, matrix_path=None):
    """Locate cpg_matrix.tsv across known result layouts."""
    candidates = []
    if matrix_path:
        candidates.append(matrix_path)
    candidates += [
        os.path.join(rd, "1_process", "5_merged_matrix", "cpg_matrix.tsv"),
        os.path.join(rd, "cpg_matrix", "cpg_matrix.tsv"),
        os.path.join(rd, "5_merged_matrix", "cpg_matrix.tsv"),
    ]
    for p in candidates:
        if p and os.path.exists(p):
            return p
    return None


def build_correlation(rd: str, data_dir: str, group_labels: dict,
                      matrix_path=None) -> dict:
    """Compute full Pearson correlation matrix from cpg_matrix.tsv.

    Retained for reference / manual use. Section 2.5 (Sample Correlation) was
    removed from the report, so build_all() no longer calls this.
    """
    matrix_path = _find_matrix(rd, matrix_path)

    if not matrix_path:
        return _save(data_dir, "correlation",
                     _missing_fig(
                         "Correlation: cpg_matrix.tsv not found. "
                         "Searched: results/1_process/5_merged_matrix/cpg_matrix.tsv",
                         "correlation"))

    print(f"[data_builder] computing Pearson correlation from {matrix_path} ...")

    try:
        # Read full matrix — for 10 samples × 5M CpGs ~400MB in RAM, acceptable
        df = pd.read_csv(matrix_path, sep="\t", index_col=0, low_memory=False)
        # Drop non-numeric index columns if any
        df = df.apply(pd.to_numeric, errors="coerce")
        # Drop rows where all values are NA
        df = df.dropna(how="all")

        # Fill NA with column mean for correlation computation
        df_filled = df.fillna(df.mean())

        corr = df_filled.corr(method="pearson")
        samples = corr.columns.tolist()

        print(f"[data_builder] correlation matrix: {len(samples)}×{len(samples)} samples")

    except Exception as e:
        return _save(data_dir, "correlation",
                     _missing_fig(f"Correlation: failed to compute ({e})", "correlation"))

    # Build colour labels (group membership)
    sample_to_group = {s: g for g, members in group_labels.items() for s in members}
    group_list = list(group_labels.keys())

    # Reorder samples: group_labels order
    ordered_samples = [s for g in group_list for s in group_labels[g] if s in samples]
    # Add any samples not in group_labels
    ordered_samples += [s for s in samples if s not in ordered_samples]

    corr_ordered = corr.loc[ordered_samples, ordered_samples]
    z_vals = corr_ordered.values.tolist()

    # Annotations: show correlation value in each cell
    text_grid = [[f"{v:.3f}" for v in row] for row in z_vals]

    # Custom colorscale: blue → white → red
    colorscale = [
        [0.0,  "#2166ac"],
        [0.25, "#92c5de"],
        [0.5,  "#f7f7f7"],
        [0.75, "#f4a582"],
        [1.0,  "#b2182b"],
    ]

    # Group boundary lines (shapes)
    shapes = []
    offset = 0
    for g in group_list:
        members_in_corr = [s for s in group_labels.get(g, []) if s in ordered_samples]
        n = len(members_in_corr)
        if n > 0 and offset > 0:
            shapes.append({
                "type": "line",
                "x0": offset - 0.5, "x1": offset - 0.5,
                "y0": -0.5, "y1": len(ordered_samples) - 0.5,
                "line": {"color": "#333", "width": 1.5},
            })
            shapes.append({
                "type": "line",
                "y0": offset - 0.5, "y1": offset - 0.5,
                "x0": -0.5, "x1": len(ordered_samples) - 0.5,
                "line": {"color": "#333", "width": 1.5},
            })
        offset += n

    n = len(ordered_samples)
    fig_h = max(350, n * 45 + 80)

    fig = {
        "data": [{
            "type": "heatmap",
            "z": z_vals,
            "x": ordered_samples,
            "y": ordered_samples,
            "text": text_grid,
            "texttemplate": "%{text}",
            "textfont": {"size": 9},
            "colorscale": colorscale,
            "zmin": -1, "zmax": 1,
            "colorbar": {"title": "Pearson r", "thickness": 14, "len": 0.8},
            "hovertemplate": (
                "Sample A: %{x}<br>Sample B: %{y}<br>r = %{z:.4f}<extra></extra>"
            ),
        }],
        "layout": {
            **_layout(
                title="Sample Correlation (Pearson, all CpG sites)",
                height=fig_h,
            ),
            "xaxis": {
                "tickangle": -40, "tickfont": {"size": 10},
                "automargin": True, "side": "bottom",
            },
            "yaxis": {
                "tickfont": {"size": 10}, "automargin": True,
                "autorange": "reversed",
            },
            "shapes": shapes,
        },
    }
    return _save(data_dir, "correlation", fig)


# ── 3.4 DMR Volcano ───────────────────────────────────────────────────────────

def build_dmr_volcano(rd: str, data_dir: str,
                      group_a: str = "", group_b: str = "",
                      q_thr: float = 0.05) -> dict:
    """
    Read dmr_annotated.bed and build an interactive volcano plot.
    Expected columns (from metilene + annotation):
      - mean_diff  (or similar: diff, delta_meth, mean_CpG_diff)
      - q          (q-value, from metilene)
      - gene / gene_name / nearest_gene (annotation)
    """
    bed_path = os.path.join(rd, "3_differential", "dmr", "dmr_annotated.bed")
    if not os.path.exists(bed_path):
        return _save(data_dir, "dmr_volcano",
                     _missing_fig("DMR: dmr_annotated.bed not found", "dmr_volcano"))

    try:
        df = pd.read_csv(bed_path, sep="\t")
    except Exception as e:
        return _save(data_dir, "dmr_volcano",
                     _missing_fig(f"DMR: could not read annotated BED ({e})", "dmr_volcano"))

    # Detect relevant columns
    diff_col = next(
        (c for c in df.columns
         if c.lower() in ("mean_diff", "diff", "delta_meth", "mean_cpg_diff", "mean1_mean2",
                          "cpg_diff", "methylation_diff", "meandiff")),
        None
    )
    q_col = next(
        (c for c in df.columns
         if c.lower() in ("q", "q_value", "qvalue", "fdr", "adj.p", "q_value")),
        None
    )
    # Explicit priority: only use clean gene symbol columns
    gene_col = None
    for _candidate in ["symbol", "annot.symbol"]:
        if _candidate in df.columns:
            gene_col = _candidate
            break
    print(f"[dmr_volcano] gene_col={gene_col}, diff_col={diff_col}, q_col={q_col}")
    print(f"[dmr_volcano] columns: {list(df.columns)}")

    if diff_col is None or q_col is None:
        return _save(data_dir, "dmr_volcano",
                     _missing_fig(
                         f"DMR: cannot find diff/q columns. Available: {list(df.columns[:10])}",
                         "dmr_volcano"))

    df = df.dropna(subset=[diff_col, q_col])
    df[q_col] = pd.to_numeric(df[q_col], errors="coerce")
    df[diff_col] = pd.to_numeric(df[diff_col], errors="coerce")
    df = df.dropna(subset=[diff_col, q_col])
    df["neg_log10_q"] = -np.log10(df[q_col].clip(lower=1e-300))

    sig   = df[df[q_col] < q_thr]
    unsig = df[df[q_col] >= q_thr]

    def _gene_labels(rows):
        """Return gene symbol labels, filtering out NA/nan/empty values."""
        if gene_col and gene_col in rows.columns:
            def _clean(v):
                s = str(v).strip()
                if s.lower() in ("na", "nan", "none", "", ".", "-"):
                    return ""
                return s
            return [_clean(v) for v in rows[gene_col].tolist()]
        return [""] * len(rows)

    def _top10_labels(sub_df, n=10):
        """Label the top N UNIQUE genes (by significance) in this direction.

        The annotated BED repeats the same gene/region across transcripts and
        CpGs, so we dedupe by gene symbol first — for each gene keep only its
        most significant point — then label the top N genes. Every other point
        (including duplicate occurrences of a labelled gene) gets an empty
        label so each gene name is drawn exactly once.
        """
        labels = _gene_labels(sub_df)
        if gene_col is None:
            return [""] * len(sub_df)
        idx_to_label = {idx: g for idx, g in zip(sub_df.index, labels)}
        # rows that actually carry a gene symbol
        valid_idx = [i for i in sub_df.index if idx_to_label[i]]
        if not valid_idx:
            return [""] * len(sub_df)
        valid = sub_df.loc[valid_idx].copy()
        valid["_gene"] = [idx_to_label[i] for i in valid.index]
        # representative (most significant) row index for each gene
        rep_idx = valid.groupby("_gene")["neg_log10_q"].idxmax()
        # rank genes by their best significance, take top N
        gene_best = valid.groupby("_gene")["neg_log10_q"].max().sort_values(ascending=False)
        top_genes = list(gene_best.head(n).index)
        label_idx = {rep_idx[g] for g in top_genes}
        return [idx_to_label[idx] if idx in label_idx else ""
                for idx in sub_df.index]

    traces = []
    # Non-significant — downsample to max 2000 points for render performance
    if len(unsig) > 0:
        if len(unsig) > 2000:
            unsig = unsig.sample(2000, random_state=42)
        traces.append({
            "type": "scatter", "mode": "markers",
            "name": f"Not significant (q≥{q_thr})",
            "x": unsig[diff_col].tolist(),
            "y": unsig["neg_log10_q"].tolist(),
            "text": _gene_labels(unsig),
            "marker": {"color": "#bdc3c7", "size": 4, "opacity": 0.35},
            "hovertemplate": "Δmeth: %{x:.2f}<br>-log₁₀(q): %{y:.2f}<br>Gene: %{text}<extra></extra>",
        })

    # Significant hypomethylated (diff < 0)
    hypo = sig[sig[diff_col] < 0].copy()
    if len(hypo) > 0:
        top10_labels = _top10_labels(hypo)
        # Boolean mask aligned to hypo's rows (not index values)
        top10_bool = pd.Series([bool(l) for l in top10_labels], index=hypo.index)
        hypo_top  = hypo.loc[top10_bool]
        hypo_rest = hypo.loc[~top10_bool]
        if len(hypo_rest) > 0:
            traces.append({
                "type": "scatter", "mode": "markers",
                "name": f"Hypo ({group_b or 'B'}<{group_a or 'A'}, q<{q_thr})",
                "x": hypo_rest[diff_col].tolist(), "y": hypo_rest["neg_log10_q"].tolist(),
                "text": _gene_labels(hypo_rest),
                "marker": {"color": "#4e79a7", "size": 7, "opacity": 0.8},
                "showlegend": True,
                "hovertemplate": "Δmeth: %{x:.2f}<br>-log₁₀(q): %{y:.2f}<br>Gene: %{text}<extra></extra>",
            })
        if len(hypo_top) > 0:
            traces.append({
                "type": "scatter", "mode": "markers+text",
                "name": "Hypo top10",
                "x": hypo_top[diff_col].tolist(), "y": hypo_top["neg_log10_q"].tolist(),
                "text": [l for l in top10_labels if l],
                "textposition": "top center", "textfont": {"size": 9, "color": "#4e79a7"},
                "marker": {"color": "#4e79a7", "size": 9, "opacity": 1,
                           "line": {"color": "white", "width": 1}},
                "showlegend": False,
                "hovertemplate": "Δmeth: %{x:.2f}<br>-log₁₀(q): %{y:.2f}<br>Gene: %{text}<extra></extra>",
            })

    # Significant hypermethylated (diff > 0)
    hyper = sig[sig[diff_col] > 0].copy()
    if len(hyper) > 0:
        top10_labels = _top10_labels(hyper)
        top10_bool   = pd.Series([bool(l) for l in top10_labels], index=hyper.index)
        hyper_top  = hyper.loc[top10_bool]
        hyper_rest = hyper.loc[~top10_bool]
        if len(hyper_rest) > 0:
            traces.append({
                "type": "scatter", "mode": "markers",
                "name": f"Hyper ({group_b or 'B'}>{group_a or 'A'}, q<{q_thr})",
                "x": hyper_rest[diff_col].tolist(), "y": hyper_rest["neg_log10_q"].tolist(),
                "text": _gene_labels(hyper_rest),
                "marker": {"color": "#e15759", "size": 7, "opacity": 0.8},
                "showlegend": True,
                "hovertemplate": "Δmeth: %{x:.2f}<br>-log₁₀(q): %{y:.2f}<br>Gene: %{text}<extra></extra>",
            })
        if len(hyper_top) > 0:
            traces.append({
                "type": "scatter", "mode": "markers+text",
                "name": "Hyper top10",
                "x": hyper_top[diff_col].tolist(), "y": hyper_top["neg_log10_q"].tolist(),
                "text": [l for l in top10_labels if l],
                "textposition": "top center", "textfont": {"size": 9, "color": "#e15759"},
                "marker": {"color": "#e15759", "size": 9, "opacity": 1,
                           "line": {"color": "white", "width": 1}},
                "showlegend": False,
                "hovertemplate": "Δmeth: %{x:.2f}<br>-log₁₀(q): %{y:.2f}<br>Gene: %{text}<extra></extra>",
            })

    q_line = -np.log10(q_thr)
    x_min, x_max = float(df[diff_col].min()), float(df[diff_col].max())
    x_pad = (x_max - x_min) * 0.05

    fig = {
        "data": traces,
        "layout": {
            **_layout(
                title=f"DMR Volcano ({group_a} vs {group_b})",
                xaxis_title="Mean methylation difference (Δβ)",
                yaxis_title="-log₁₀(q-value)",
                height=520,
            ),
            "shapes": [
                # Horizontal q threshold line — black dashed
                {"type": "line",
                 "x0": x_min - x_pad, "x1": x_max + x_pad,
                 "y0": q_line, "y1": q_line, "xref": "x", "yref": "y",
                 "line": {"color": "#333", "width": 1.2, "dash": "dash"}},
                # Vertical zero line
                {"type": "line",
                 "x0": 0, "x1": 0, "y0": 0, "y1": 1, "xref": "x", "yref": "paper",
                 "line": {"color": "#aaa", "width": 1}},
            ],
            "annotations": [{
                "x": x_max + x_pad, "y": q_line, "xref": "x", "yref": "y",
                "text": f"q={q_thr}", "showarrow": False,
                "font": {"size": 10, "color": "#333"},
                "xanchor": "right", "yanchor": "bottom",
            }],
        },
    }
    return _save(data_dir, "dmr_volcano", fig)


# ── 5.x MESA charts ───────────────────────────────────────────────────────────

def build_mesa_charts(rd: str, data_dir: str) -> dict[str, dict]:
    """
    Build ROC, prediction heatmap, and Spearman correlation from
    loocv_predictions.tsv.

    Expected columns: sample, true_label, {modality}_prob, [combined_prob]
    """
    pred_tsv = os.path.join(rd, "5_mesa", "loocv_predictions.tsv")
    results = {}

    if not os.path.exists(pred_tsv):
        msg = "MESA: loocv_predictions.tsv not found — run cftk mesa --loocv"
        for key in ("mesa_roc", "mesa_heatmap", "mesa_spearman"):
            results[key] = _save(data_dir, key, _missing_fig(msg, key))
        return results

    try:
        pred_df = pd.read_csv(pred_tsv, sep="\t", index_col=0)
    except Exception as e:
        msg = f"MESA: cannot read loocv_predictions.tsv ({e})"
        for key in ("mesa_roc", "mesa_heatmap", "mesa_spearman"):
            results[key] = _save(data_dir, key, _missing_fig(msg, key))
        return results

    # Detect true label column and score columns
    # true_col: any column named *true* or *label* or *y_true*
    true_col = next((c for c in pred_df.columns
                     if c.lower() in ("y_true", "true_label", "label", "true", "class")), None)
    if true_col is None:
        true_col = next((c for c in pred_df.columns
                         if "true" in c.lower() or "label" in c.lower()), None)

    # prob_cols: any numeric column that isn't the true label
    # Accepts *_prob, *_score, or bare modality names (cpg, wps, occupancy, Multimodal, ...)
    all_numeric = [c for c in pred_df.columns
                   if c != true_col and pd.to_numeric(pred_df[c], errors="coerce").notna().sum() > 0]
    # Prefer columns with *prob* or *score* in name; fall back to all numeric
    prob_cols = [c for c in all_numeric if "prob" in c.lower() or "score" in c.lower()]
    if not prob_cols:
        prob_cols = all_numeric  # use all numeric non-label columns

    if true_col is None or not prob_cols:
        msg = f"MESA: could not identify true label + score columns. Got: {list(pred_df.columns)}"
        for key in ("mesa_roc", "mesa_heatmap", "mesa_spearman"):
            results[key] = _save(data_dir, key, _missing_fig(msg, key))
        return results

    # Ensure true label is binary 0/1
    y_true_raw = pd.to_numeric(pred_df[true_col], errors="coerce").fillna(0).values
    classes = np.unique(y_true_raw)
    if len(classes) == 2:
        y_true = (y_true_raw == classes[1]).astype(int)
    else:
        y_true = y_true_raw.astype(int)



    # ── ROC curves ────────────────────────────────────────────────────────────
    roc_traces = []
    try:
        from sklearn.metrics import roc_curve, auc as sklearn_auc
        has_sklearn = True
    except ImportError:
        has_sklearn = False

    for i, col in enumerate(prob_cols):
        y_score = pred_df[col].values
        label = col.replace("_prob", "").replace("_score", "")
        if has_sklearn:
            try:
                fpr, tpr, _ = roc_curve(y_true, y_score)
                auc_val = sklearn_auc(fpr, tpr)
                roc_traces.append({
                    "type": "scatter", "mode": "lines",
                    "name": f"{label} (AUC={auc_val:.3f})",
                    "x": fpr.tolist(), "y": tpr.tolist(),
                    "line": {"color": _COLORS[i % len(_COLORS)], "width": 2},
                    "hovertemplate": f"{label}<br>FPR: %{{x:.3f}}<br>TPR: %{{y:.3f}}<extra></extra>",
                })
            except Exception:
                pass
        else:
            # Simple sort-based approximation
            order = np.argsort(-y_score)
            tpr = np.cumsum(y_true[order]) / y_true.sum()
            fpr = np.cumsum(1 - y_true[order]) / (1 - y_true).sum()
            roc_traces.append({
                "type": "scatter", "mode": "lines",
                "name": label,
                "x": fpr.tolist(), "y": tpr.tolist(),
                "line": {"color": _COLORS[i % len(_COLORS)], "width": 2},
            })

    # Diagonal reference
    roc_traces.append({
        "type": "scatter", "mode": "lines",
        "name": "Random (AUC=0.500)",
        "x": [0, 1], "y": [0, 1],
        "line": {"color": "#aaa", "width": 1, "dash": "dash"},
        "showlegend": True,
    })

    roc_fig = {
        "data": roc_traces,
        "layout": {
            **_layout(
                title="MESA LOOCV ROC Curves",
                xaxis_title="False Positive Rate",
                yaxis_title="True Positive Rate",
                height=420,
            ),
            "xaxis": {"title": "False Positive Rate", "range": [-0.03, 1.03],
                      "automargin": True, "tickfont": {"size": 10}},
            "yaxis": {"title": "True Positive Rate", "range": [-0.03, 1.03],
                      "automargin": True},
        },
    }
    results["mesa_roc"] = _save(data_dir, "mesa_roc", roc_fig)

    # ── Prediction heatmap ────────────────────────────────────────────────────
    # Sort samples by true label then by prob of first modality
    sort_key = true_col if true_col in pred_df.columns else pred_df.columns[0]
    pred_sorted = pred_df.sort_values([true_col, prob_cols[0]], ascending=[True, False])

    z_vals = pred_sorted[prob_cols].values.T.tolist()   # modalities × samples
    x_labs = pred_sorted.index.tolist()
    y_labs = [c.replace("_prob", "").replace("_score", "") for c in prob_cols]
    text_grid = [[f"{v:.3f}" for v in row] for row in pred_sorted[prob_cols].values.T.tolist()]

    hmap_fig = {
        "data": [{
            "type": "heatmap",
            "z": z_vals,
            "x": x_labs, "y": y_labs,
            "text": text_grid, "texttemplate": "%{text}",
            "textfont": {"size": 9},
            "colorscale": [
                [0.0, "#2166ac"], [0.25, "#92c5de"], [0.5, "#f7f7f7"],
                [0.75, "#f4a582"], [1.0, "#b2182b"]
            ],
            "zmin": 0, "zmax": 1,
            "colorbar": {"title": "Prob", "thickness": 14, "len": 0.8},
            "hovertemplate": "Sample: %{x}<br>Modality: %{y}<br>Prob: %{z:.4f}<extra></extra>",
        }],
        "layout": {
            **_layout(
                title="MESA LOOCV Prediction Probabilities",
                height=max(250, len(prob_cols) * 50 + 120),
            ),
            "xaxis": {"tickangle": -40, "tickfont": {"size": 9}, "automargin": True},
            "yaxis": {"tickfont": {"size": 10}, "automargin": True},
        },
    }
    results["mesa_heatmap"] = _save(data_dir, "mesa_heatmap", hmap_fig)

    # ── Spearman correlation ──────────────────────────────────────────────────
    try:
        sp_corr = pred_df[prob_cols].corr(method="spearman")
        sp_labels = [c.replace("_prob", "").replace("_score", "") for c in prob_cols]
        sp_vals = sp_corr.values.tolist()
        sp_text = [[f"{v:.3f}" for v in row] for row in sp_corr.values.tolist()]

        # Lower-triangle only: set upper triangle to None
        n_mod = len(sp_labels)
        sp_lower = []
        for ri in range(n_mod):
            row = []
            for ci in range(n_mod):
                if ci > ri:
                    row.append(None)   # upper triangle → blank
                else:
                    row.append(round(float(sp_corr.iloc[ri, ci]), 4))
            sp_lower.append(row)
        sp_text_lower = []
        for ri in range(n_mod):
            row = []
            for ci in range(n_mod):
                if ci > ri:
                    row.append("")
                else:
                    row.append(f"{sp_corr.iloc[ri, ci]:.2f}")
            sp_text_lower.append(row)

        sp_fig = {
            "data": [{
                "type": "heatmap",
                "z": sp_lower, "x": sp_labels, "y": sp_labels,
                "text": sp_text_lower, "texttemplate": "%{text}",
                "textfont": {"size": 13, "color": "white"},
                "colorscale": [
                    [0.0, "#ffffcc"], [0.25, "#fd8d3c"],
                    [0.6, "#bd0026"], [1.0, "#67000d"]
                ],
                "zmin": 0, "zmax": 1,
                "colorbar": {"title": "Spearman ρ", "thickness": 14,
                             "tickvals": [0, 0.2, 0.4, 0.6, 0.8, 1.0]},
                "hovertemplate": "%{y} vs %{x}<br>ρ = %{z:.4f}<extra></extra>",
                "showscale": True,
                "xgap": 3, "ygap": 3,
            }],
            "layout": {
                **_layout(
                    title="MESA Modality Spearman Correlation",
                    height=max(300, n_mod * 70 + 100),
                ),
                "xaxis": {"tickangle": -30, "tickfont": {"size": 11},
                          "automargin": True, "side": "bottom"},
                "yaxis": {"tickfont": {"size": 11}, "automargin": True,
                          "autorange": "reversed"},
                "plot_bgcolor": "#ffffff",
                "margin": {"l": 100, "r": 80, "t": 50, "b": 100},
            },
        }
    except Exception as e:
        sp_fig = _missing_fig(f"Spearman: {e}", "mesa_spearman")

    results["mesa_spearman"] = _save(data_dir, "mesa_spearman", sp_fig)
    return results




# ── 3.2 Violin ────────────────────────────────────────────────────────────────

def build_violin(rd: str, data_dir: str, group_labels: dict,
                 matrix_path=None, top_n: int = 2000) -> dict:
    """
    Build interactive violin/box plot of CpG beta values.
    Uses top N most variable CpG sites from cpg_matrix.tsv per modality.
    Saves data/violin.txt (all modalities in one file via dropdown).
    """
    matrix_path = _find_matrix(rd, matrix_path)
    if not matrix_path:
        return _save(data_dir, "violin",
                     _missing_fig("Violin: cpg_matrix.tsv not found", "violin"))

    try:
        print(f"[data_builder] violin: reading {matrix_path} ...")
        df = pd.read_csv(matrix_path, sep="\t", index_col=0, low_memory=False)
        df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        # Top N most variable CpGs
        var = df.var(axis=1).fillna(0)
        top_idx = var.nlargest(top_n).index
        df_top = df.loc[top_idx]
        print(f"[data_builder] violin: using top {len(df_top)} variable CpGs")
    except Exception as e:
        return _save(data_dir, "violin",
                     _missing_fig(f"Violin: failed to read matrix ({e})", "violin"))

    sample_to_group = {s: g for g, members in group_labels.items() for s in members}
    group_list = list(group_labels.keys())

    # One violin per group (aggregate all samples in group together)
    traces = []
    import random as _random
    _random.seed(42)
    for gi, (grp, members) in enumerate(group_labels.items()):
        cols = [c for c in members if c in df_top.columns]
        if not cols:
            continue
        # Pool all samples in this group
        all_vals = pd.concat([df_top[c].dropna() for c in cols]).tolist()
        # Subsample max 8000 points for performance
        if len(all_vals) > 8000:
            all_vals = _random.sample(all_vals, 8000)
        traces.append({
            "type": "violin",
            "name": grp,
            "y": all_vals,
            "x0": grp,
            "line": {"color": _GROUP_COLORS[gi % len(_GROUP_COLORS)]},
            "fillcolor": _GROUP_COLORS[gi % len(_GROUP_COLORS)],
            "opacity": 0.75,
            "meanline": {"visible": True, "color": "white", "width": 2},
            "box": {"visible": True, "fillcolor": "white", "width": 0.15},
            "points": False,
            "showlegend": True,
            "hovertemplate": f"<b>{grp}</b><br>β: %{{y:.4f}}<extra></extra>",
        })

    if not traces:
        return _save(data_dir, "violin",
                     _missing_fig("Violin: no sample data found in matrix", "violin"))

    fig = {
        "data": traces,
        "layout": {
            **_layout(
                title=f"CpG β-value Distribution (top {top_n} variable CpGs)",
                yaxis_title="β value",
                height=500,
            ),
            "violinmode": "overlay",
            "xaxis": {"title": "", "automargin": True, "tickfont": {"size": 12}},
            "yaxis": {"range": [-0.05, 1.05], "tickfont": {"size": 10},
                      "automargin": True},
            "margin": {"l": 60, "r": 20, "t": 50, "b": 60},
            "legend": {"orientation": "h", "y": -0.12},
        },
    }
    return _save(data_dir, "violin", fig)


# ── 3.3 DMC Heatmap ──────────────────────────────────────────────────────────

def build_dmc_heatmap(rd: str, data_dir: str, group_labels: dict,
                      matrix_path=None, top_n: int = 500) -> dict:
    """
    Build interactive heatmap of top N most variable CpG sites.
    Rows = CpG sites (sorted by variance desc), cols = samples (grouped).
    """
    matrix_path = _find_matrix(rd, matrix_path)
    if not matrix_path:
        return _save(data_dir, "dmc_heatmap",
                     _missing_fig("Heatmap: cpg_matrix.tsv not found", "dmc_heatmap"))

    try:
        print(f"[data_builder] heatmap: reading {matrix_path} ...")
        df = pd.read_csv(matrix_path, sep="\t", index_col=0, low_memory=False)
        df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all")
        var = df.var(axis=1).fillna(0)
        top_idx = var.nlargest(top_n).index
        df_top = df.loc[top_idx]
        print(f"[data_builder] heatmap: {len(df_top)} CpGs × {len(df_top.columns)} samples")
    except Exception as e:
        return _save(data_dir, "dmc_heatmap",
                     _missing_fig(f"Heatmap: failed to read matrix ({e})", "dmc_heatmap"))

    # Column order: by group
    group_list = list(group_labels.keys())
    ordered_cols = [s for g in group_list for s in group_labels[g] if s in df_top.columns]
    ordered_cols += [c for c in df_top.columns if c not in ordered_cols]
    df_ordered = df_top[ordered_cols]

    # Check for data completeness
    pct_nan = df_ordered.isna().sum().sum() / df_ordered.size
    if pct_nan > 0.95:
        return _save(data_dir, "dmc_heatmap",
                     _missing_fig(
                         f"Heatmap: >95% of values are NaN after numeric conversion. "
                         f"Check cpg_matrix.tsv format (expected tab-sep with numeric β values).",
                         "dmc_heatmap"))

    # Fill NaN with 0 for display (NaN renders as blank in Plotly heatmap)
    df_display = df_ordered.fillna(0)

    # Row labels: shorten if long
    row_labels = [str(idx)[:40] for idx in df_display.index.tolist()]

    z_vals = df_display.values.tolist()

    # Group boundary vertical lines
    shapes = []
    offset = 0
    for g in group_list:
        n_g = sum(1 for s in group_labels.get(g, []) if s in ordered_cols)
        if n_g > 0 and offset > 0:
            shapes.append({
                "type": "line",
                "x0": offset - 0.5, "x1": offset - 0.5,
                "y0": -0.5, "y1": len(row_labels) - 0.5,
                "line": {"color": "#333", "width": 1.5},
            })
        offset += n_g

    n_rows = len(row_labels)
    fig_h  = max(400, min(n_rows * 4 + 150, 800))  # cap at 800px

    # Hide row tick labels if too many rows
    show_row_ticks = n_rows <= 80

    fig = {
        "data": [{
            "type": "heatmap",
            "z": z_vals,
            "x": ordered_cols,
            "y": row_labels if show_row_ticks else [""] * n_rows,
            "colorscale": [
                [0.0, "#2166ac"],
                [0.5, "#f7f7f7"],
                [1.0, "#b2182b"],
            ],
            "zmin": 0, "zmax": 1,
            "colorbar": {"title": "β", "thickness": 14, "len": 0.8},
            "hovertemplate": (
                "CpG: %{y}<br>Sample: %{x}<br>β = %{z:.4f}<extra></extra>"
            ),
        }],
        "layout": {
            **_layout(
                title=f"Top {top_n} Variable CpG Sites",
                height=fig_h,
            ),
            "xaxis": {"tickangle": -40, "tickfont": {"size": 10},
                      "automargin": True, "side": "bottom"},
            "yaxis": {"tickfont": {"size": 8}, "automargin": True,
                      "autorange": "reversed",
                      "showticklabels": show_row_ticks},
            "margin": {"l": 80, "r": 80, "t": 50, "b": 80},
            "shapes": shapes,
        },
    }
    return _save(data_dir, "dmc_heatmap", fig)


# ── 3.3 DMR significant table ─────────────────────────────────────────────────

def build_dmr_table(rd: str, data_dir: str, q_thr: float = 0.05,
                   top_n: int = 50) -> dict:
    """
    Build a sortable table of the top significant, gene-annotated DMRs from
    dmr_annotated.bed.

    Columns: Gene (annot.symbol), Chr, Start, End, Δβ (mean_diff), q-value,
             -log10(q), n_CpGs, Direction.

    The annotated BED repeats the same gene/region across transcripts and CpGs,
    so this:
      1. keeps only rows that carry a gene symbol,
      2. de-duplicates by gene symbol (keeps each gene's most significant row),
      3. takes the top N (default 50) genes by q-value.
    Stored as JSON {columns, rows, ...} for the sortable HTML table.
    """
    bed_path = os.path.join(rd, "3_differential", "dmr", "dmr_annotated.bed")
    if not os.path.exists(bed_path):
        return _save(data_dir, "dmr_table",
                     {"_missing": True, "_msg": "dmr_annotated.bed not found"})

    try:
        df = pd.read_csv(bed_path, sep="\t")
    except Exception as e:
        return _save(data_dir, "dmr_table",
                     {"_missing": True, "_msg": f"could not read DMR bed ({e})"})

    # Column resolution (match real header)
    def _col(*cands):
        for c in cands:
            if c in df.columns:
                return c
        return None

    chr_col   = _col("annot.seqnames", "seqnames", "chrom", "chr")
    start_col = _col("start", "annot.start")
    end_col   = _col("end", "annot.end")
    q_col     = _col("q_value", "q", "qvalue", "fdr")
    diff_col  = _col("mean_diff", "diff", "delta_meth")
    ncpg_col  = _col("n_CpGs", "n_cpgs", "nCpGs", "num_cpgs")
    gene_col  = _col("annot.symbol", "symbol")

    if q_col is None or diff_col is None:
        return _save(data_dir, "dmr_table",
                     {"_missing": True,
                      "_msg": f"missing q/diff columns; have {list(df.columns)[:8]}"})
    if gene_col is None:
        return _save(data_dir, "dmr_table",
                     {"_missing": True,
                      "_msg": "no gene-symbol column (annot.symbol) in dmr_annotated.bed"})

    def _clean_gene(v):
        s = str(v).strip()
        return "" if s.lower() in ("na", "nan", "none", "", ".", "-") else s

    df = df.copy()
    df[q_col]    = pd.to_numeric(df[q_col], errors="coerce")
    df[diff_col] = pd.to_numeric(df[diff_col], errors="coerce")
    df = df.dropna(subset=[q_col, diff_col])
    df["_gene"] = df[gene_col].map(_clean_gene)

    # significant + gene-annotated only
    sig = df[(df[q_col] < q_thr) & (df["_gene"] != "")]
    n_total_sig_genes = sig["_gene"].nunique()
    if sig.empty:
        return _save(data_dir, "dmr_table",
                     {"_missing": True,
                      "_msg": f"no gene-annotated DMRs below q < {q_thr}"})

    # de-duplicate by gene symbol: keep each gene's most significant (min q) row
    sig = sig.sort_values(q_col, ascending=True)
    sig = sig.drop_duplicates(subset="_gene", keep="first")
    # top N genes by significance
    sig = sig.head(top_n)

    columns = ["Gene", "Chr", "Start", "End", "Δβ", "q-value",
               "-log10(q)", "n_CpGs", "Direction"]
    rows = []
    for _, r in sig.iterrows():
        gene  = r["_gene"]
        chrom = str(r.get(chr_col, "")) if chr_col else ""
        start = int(r.get(start_col, 0)) if start_col and not pd.isna(r.get(start_col)) else ""
        end   = int(r.get(end_col, 0))   if end_col   and not pd.isna(r.get(end_col))   else ""
        diff  = float(r[diff_col])
        qv    = float(r[q_col])
        neglq = round(-np.log10(max(qv, 1e-300)), 3)
        ncpg  = int(r.get(ncpg_col, 0)) if ncpg_col and not pd.isna(r.get(ncpg_col)) else ""
        direction = "Hyper" if diff > 0 else "Hypo"
        rows.append([
            gene, chrom, start, end,
            round(diff, 3), f"{qv:.2e}", neglq, ncpg, direction
        ])

    fig = {"columns": columns, "rows": rows,
           "n_total_sig": int(n_total_sig_genes),
           "shown": len(rows), "q_thr": q_thr}
    return _save(data_dir, "dmr_table", fig)


# ── 4.2 End Motif (per-group charts) ──────────────────────────────────────────

def build_end_motif(rd: str, data_dir: str, group_labels: dict,
                    top_n: int = 20) -> dict:
    """
    Build SEPARATE interactive end-motif bar charts, one per group.
    Saves data/end_motif_{group}.txt and returns {f"end_motif_{group}": fig}.
    report_generator._sec_fragmentomics() requests these per-group keys.

    Input files: 4_fragmentomics/end_motif/{sample}_4mer.tsv
      Format: tab-separated, col0 = motif (e.g. AAAA), col1 = frequency.

    For each group, produces one chart with a dropdown switching between:
      - the group mean (default)
      - each sample within that group
    The motif order (x axis) is fixed to the group-mean top N for comparability.

    Missing data is rendered as a per-group placeholder figure so the report
    shows a clean "⚠" box instead of a raw "Chart data not available" message.
    """
    em_dir    = os.path.join(rd, "4_fragmentomics", "end_motif")
    tsv_files = sorted(glob.glob(os.path.join(em_dir, "*_4mer.tsv")))
    results: dict = {}

    def _miss_all(msg):
        for grp in group_labels:
            key = f"end_motif_{grp}"
            results[key] = _save(data_dir, key, _missing_fig(msg, key))
        return results

    if not tsv_files:
        return _miss_all("End motif: no *_4mer.tsv files found — run cftk frag --end-motif")

    # Load all samples → dict[sample] = pd.Series(index=motif, value=freq)
    sample_freq: dict[str, pd.Series] = {}
    for fp in tsv_files:
        name = os.path.basename(fp).replace("_4mer.tsv", "")
        try:
            t = pd.read_csv(fp, sep="\t", header=None, names=["motif", "freq"]).dropna()
            sample_freq[name] = pd.Series(
                t["freq"].values, index=t["motif"].astype(str).values)
        except Exception as e:
            print(f"[end_motif] could not read {fp}: {e}")

    if not sample_freq:
        return _miss_all("End motif: could not parse any 4mer TSV")

    all_df = pd.DataFrame(sample_freq)            # rows=motif, cols=sample

    for gi, (grp, members) in enumerate(group_labels.items()):
        key  = f"end_motif_{grp}"
        cols = [c for c in members if c in all_df.columns]
        if not cols:
            results[key] = _save(
                data_dir, key,
                _missing_fig(f"End motif: no sample data for group '{grp}'", key))
            continue

        gmean = all_df[cols].mean(axis=1)
        # Fixed motif order = this group's mean top N
        top_motifs = gmean.sort_values(ascending=False).head(top_n).index.tolist()
        base_color = _GROUP_COLORS[gi % len(_GROUP_COLORS)]

        traces = []
        button_specs = []

        # Group mean (default visible)
        traces.append({
            "type": "bar", "name": f"{grp} (mean)",
            "x": top_motifs,
            "y": [round(float(gmean.get(m, 0.0)), 6) for m in top_motifs],
            "marker": {"color": base_color},
            "visible": True,
            "hovertemplate": f"{grp} mean<br>%{{x}}<br>freq: %{{y:.5f}}<extra></extra>",
        })
        button_specs.append((f"{grp} (mean)", 0))

        # Each sample in this group
        for sname in cols:
            s = sample_freq[sname]
            traces.append({
                "type": "bar", "name": sname,
                "x": top_motifs,
                "y": [round(float(s.get(m, 0.0)), 6) for m in top_motifs],
                "marker": {"color": base_color, "opacity": 0.85},
                "visible": False,
                "hovertemplate": f"{sname}<br>%{{x}}<br>freq: %{{y:.5f}}<extra></extra>",
            })
            button_specs.append((sname, len(traces) - 1))

        n_total = len(traces)
        options = []
        for label, tidx in button_specs:
            vis = [False] * n_total
            vis[tidx] = True
            options.append({
                "label": label,
                "visible": vis,
                "relayout": {"title.text": f"End Motif ({grp}) — {label} (top {top_n})"},
            })

        fig = {
            "data": traces,
            "layout": {
                **_layout(
                    title=f"End Motif ({grp}) — {grp} (mean) (top {top_n})",
                    xaxis_title="4-mer end motif",
                    yaxis_title="Frequency",
                    height=440,
                ),
                "xaxis": {"tickangle": -55, "tickfont": {"size": 10},
                          "automargin": True, "categoryorder": "array",
                          "categoryarray": top_motifs},
                "yaxis": {"automargin": True},
                "margin": {"l": 60, "r": 20, "t": 60, "b": 90},
                "showlegend": False,
            },
            # consumed by report_generator._plotly_select_render (HTML <select>)
            "_sample_menu": {"options": options},
        }
        results[key] = _save(data_dir, key, fig)

    return results


# ── 4.2b End Motif — grouped box plot (top 20, samples as points) ─────────────

def build_end_motif_box(rd: str, data_dir: str, group_labels: dict,
                        top_n: int = 20) -> dict:
    """
    Build ONE grouped box plot comparing groups across the top-N end motifs.

    x-axis  : top N 4-mer motifs ranked by OVERALL mean across all samples
    boxes   : one box per group at each motif (boxmode='group'), summarising
              the per-sample frequency distribution within that group
    points  : every individual sample is overlaid as a jittered point

    Saved under key 'end_motif_box'. Complements the per-group bar charts.
    """
    em_dir    = os.path.join(rd, "4_fragmentomics", "end_motif")
    tsv_files = sorted(glob.glob(os.path.join(em_dir, "*_4mer.tsv")))
    if not tsv_files:
        return _save(data_dir, "end_motif_box",
                     _missing_fig("End motif: no *_4mer.tsv files found — run cftk frag --end-motif",
                                  "end_motif_box"))

    sample_freq: dict[str, pd.Series] = {}
    for fp in tsv_files:
        name = os.path.basename(fp).replace("_4mer.tsv", "")
        try:
            t = pd.read_csv(fp, sep="\t", header=None, names=["motif", "freq"]).dropna()
            sample_freq[name] = pd.Series(
                t["freq"].values, index=t["motif"].astype(str).values)
        except Exception as e:
            print(f"[end_motif_box] could not read {fp}: {e}")

    if not sample_freq:
        return _save(data_dir, "end_motif_box",
                     _missing_fig("End motif: could not parse any 4mer TSV", "end_motif_box"))

    all_df = pd.DataFrame(sample_freq)               # rows=motif, cols=sample
    overall_mean = all_df.mean(axis=1)
    top_motifs = overall_mean.sort_values(ascending=False).head(top_n).index.tolist()

    traces = []
    for gi, (grp, members) in enumerate(group_labels.items()):
        cols = [c for c in members if c in all_df.columns]
        if not cols:
            continue
        xs, ys = [], []
        for m in top_motifs:
            for s in cols:
                v = all_df.at[m, s] if m in all_df.index else None
                if v is None or (isinstance(v, float) and np.isnan(v)):
                    continue
                xs.append(m)
                ys.append(round(float(v), 6))
        color = _GROUP_COLORS[gi % len(_GROUP_COLORS)]
        traces.append({
            "type": "box",
            "name": grp,
            "x": xs,
            "y": ys,
            "boxpoints": "all",     # show every sample as a point
            "jitter": 0.5,
            "pointpos": 0,
            "marker": {"color": color, "size": 4, "opacity": 0.7},
            "line": {"color": color},
            "fillcolor": color,
            "opacity": 0.55,
            "whiskerwidth": 0.5,
            "hovertemplate": f"{grp}<br>%{{x}}<br>freq: %{{y:.5f}}<extra></extra>",
        })

    if not traces:
        return _save(data_dir, "end_motif_box",
                     _missing_fig("End motif: no sample data in any group", "end_motif_box"))

    fig = {
        "data": traces,
        "layout": {
            **_layout(
                title=f"End Motif — group comparison (top {top_n} by overall mean)",
                xaxis_title="4-mer end motif",
                yaxis_title="Frequency",
                height=480,
            ),
            "boxmode": "group",
            "xaxis": {"tickangle": -55, "tickfont": {"size": 10},
                      "automargin": True, "categoryorder": "array",
                      "categoryarray": top_motifs},
            "yaxis": {"automargin": True},
            "margin": {"l": 60, "r": 20, "t": 54, "b": 90},
            "legend": {"orientation": "h", "y": -0.22, "x": 0},
            "showlegend": True,
        },
    }
    return _save(data_dir, "end_motif_box", fig)


# ── Main builder ──────────────────────────────────────────────────────────────

def build_all(rd: str, data_dir: str, group_labels: dict,
              group_a: str = "", group_b: str = "",
              q_thr: float = 0.05, frag_len: int = 167,
              matrix_path: str | None = None,
              source_results: dict[str, str] | None = None) -> dict[str, dict]:
    """
    Build all interactive chart data files. Returns dict key→figure.
    Missing data is handled gracefully (produces a placeholder figure).

    Section 2.5 (Sample Correlation) was removed from the report, so the
    correlation figure is no longer built (build_correlation remains available
    for manual use).
    """
    print(f"[data_builder] building interactive chart data → {data_dir}")
    results = {}

    # 2.2 Fragment Length
    results["fragment_length"] = build_fragment_length(
        rd, data_dir, group_labels, source_results=source_results
    )

    # 2.3 Dinucleotide
    results["dinucleotide"] = build_dinucleotide(
        rd, data_dir, frag_len=frag_len, source_results=source_results
    )

    # 3.1 PCA
    results["pca"] = build_pca(rd, data_dir, group_labels)

    # 3.2 Violin
    results["violin"] = build_violin(rd, data_dir, group_labels, matrix_path)

    # 3.3 DMC Heatmap
    results["dmc_heatmap"] = build_dmc_heatmap(rd, data_dir, group_labels, matrix_path)

    # 3.4 DMR Volcano
    results["dmr_volcano"] = build_dmr_volcano(rd, data_dir, group_a, group_b, q_thr)

    # 3.3 DMR significant table (sortable)
    results["dmr_table"] = build_dmr_table(rd, data_dir, q_thr=q_thr, top_n=100)

    # 4.2 End Motif — per-group charts (keys: end_motif_{group})
    results.update(build_end_motif(rd, data_dir, group_labels))

    # 4.2b End Motif — grouped box plot (key: end_motif_box)
    results["end_motif_box"] = build_end_motif_box(rd, data_dir, group_labels)

    # 5.x MESA
    mesa = build_mesa_charts(rd, data_dir)
    results.update(mesa)

    print(f"[data_builder] done. {len(results)} figures built.")
    return results
