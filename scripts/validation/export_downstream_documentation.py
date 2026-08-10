#!/usr/bin/env python3
"""Export public-safe documentation evidence from a completed 5+5 run.

This maintainer utility reads real CFTK outputs, replaces sample identifiers
with deterministic aliases, and writes only aggregate figures and metadata.
It never copies source tables or the private HTML report.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve


DPI = 180
CHUNK_ROWS = 50_000
CONTROL_COLOR = "#0072B2"
CASE_COLOR = "#D55E00"
ACCENT_COLOR = "#009E73"
MUTED_COLOR = "#6B7280"
LIGHT_COLOR = "#D9E2E8"
PANEL_SCOPE_NOTE = (
    "Targeted-panel mode restricts WPS, occupancy, and DELFI to reads and "
    "intervals overlapping the configured target BED; these outputs are not "
    "genome-wide measurements."
)
OUTPUT_NAMES = {
    "differential": "validation_10sample_differential.png",
    "fragmentomics": "validation_10sample_fragmentomics.png",
    "mesa": "validation_10sample_mesa.png",
    "report": "validation_10sample_report_preview.png",
    "summary": "validation_10sample_downstream_summary.json",
}


def _require_file(path: Path) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Required nonempty artifact not found: {path}")
    return path


def _style_axis(ax, grid_axis: str = "both") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=grid_axis, color=LIGHT_COLOR, linewidth=0.7, alpha=0.75)
    ax.set_axisbelow(True)


def _save(fig, path: Path) -> None:
    fig.savefig(
        path,
        dpi=DPI,
        facecolor="white",
        bbox_inches="tight",
        metadata={"Software": "CFTK public documentation exporter"},
    )
    plt.close(fig)


def _read_cohort(project_root: Path) -> dict:
    sample_sheet = pd.read_csv(_require_file(project_root / "samples.tsv"), sep="\t", dtype=str)
    required = {"sample", "group", "role"}
    missing = required - set(sample_sheet.columns)
    if missing:
        raise ValueError(f"samples.tsv is missing columns: {sorted(missing)}")

    sample_sheet = sample_sheet.loc[:, ["sample", "group", "role"]].copy()
    sample_sheet["role"] = sample_sheet["role"].str.lower()
    if sample_sheet["sample"].duplicated().any():
        raise ValueError("samples.tsv contains duplicate sample identifiers")

    controls = sample_sheet[sample_sheet["role"] == "control"]
    cases = sample_sheet[sample_sheet["role"] == "case"]
    if len(controls) != 5 or len(cases) != 5 or len(sample_sheet) != 10:
        raise ValueError("Documentation export requires exactly five controls and five cases")
    if set(controls["group"]) != {"Control"} or set(cases["group"]) != {"sALS"}:
        raise ValueError("Documentation export expects Control and sALS sample groups")

    ordered_samples = sample_sheet["sample"].tolist()
    aliases = {
        sample: f"Control_{index}"
        for index, sample in enumerate(controls["sample"], start=1)
    }
    aliases.update({
        sample: f"sALS_{index}"
        for index, sample in enumerate(cases["sample"], start=1)
    })
    groups = dict(zip(sample_sheet["sample"], sample_sheet["group"]))
    return {
        "ordered_samples": ordered_samples,
        "aliases": aliases,
        "groups": groups,
        "control_samples": controls["sample"].tolist(),
        "case_samples": cases["sample"].tolist(),
    }


def _read_pca(project_root: Path, modality: str, cohort: dict) -> tuple[pd.DataFrame, dict]:
    base = project_root / "results" / "3_differential" / modality
    coordinates = pd.read_csv(_require_file(base / "pca_coordinates.txt"), sep="\t")
    sample_column = coordinates.columns[0]
    coordinates = coordinates.rename(columns={sample_column: "sample"})
    required = {"sample", "group", "PC1", "PC2"}
    if not required.issubset(coordinates.columns):
        raise ValueError(f"Unexpected PCA schema for {modality}: {list(coordinates.columns)}")
    if set(coordinates["sample"]) != set(cohort["ordered_samples"]):
        raise ValueError(f"PCA sample set does not match samples.tsv for {modality}")
    expected_groups = coordinates["sample"].map(cohort["groups"])
    if not coordinates["group"].equals(expected_groups):
        raise ValueError(f"PCA groups do not match samples.tsv for {modality}")
    coordinates["alias"] = coordinates["sample"].map(cohort["aliases"])
    coordinates["order"] = coordinates["sample"].map(
        {sample: index for index, sample in enumerate(cohort["ordered_samples"])}
    )
    coordinates = coordinates.sort_values("order")

    variance = pd.read_csv(_require_file(base / "pca_variance.txt"), sep="\t")
    variance = variance.set_index("PC")["variance_explained_pct"]
    return coordinates, {"PC1": float(variance["PC1"]), "PC2": float(variance["PC2"])}


def _read_dmr(project_root: Path) -> pd.DataFrame:
    columns = [
        "chrom", "start", "end", "q_value", "mean_diff", "n_cpg",
        "p_mwu", "p_2dks", "mean_control", "mean_sals",
    ]
    dmr = pd.read_csv(
        _require_file(project_root / "results" / "3_differential" / "dmr" / "dmr_raw.bed"),
        sep="\t",
        header=None,
        names=columns,
        usecols=[3, 4],
    )
    dmr["q_value"] = pd.to_numeric(dmr["q_value"], errors="coerce")
    dmr["mean_diff"] = pd.to_numeric(dmr["mean_diff"], errors="coerce")
    dmr = dmr.dropna().reset_index(drop=True)
    if dmr.empty or not dmr["q_value"].between(0, 1).all():
        raise ValueError("dmr_raw.bed has no valid rows or contains q values outside [0, 1]")
    return dmr


def _plot_pca(ax, coordinates: pd.DataFrame, variance: dict, title: str) -> None:
    for group, color, marker in (
        ("Control", CONTROL_COLOR, "s"),
        ("sALS", CASE_COLOR, "o"),
    ):
        subset = coordinates[coordinates["group"] == group]
        ax.scatter(
            subset["PC1"], subset["PC2"], s=42, color=color, marker=marker,
            edgecolor="white", linewidth=0.5, label=f"{group} (n={len(subset)})",
        )
    ax.set(
        xlabel=f"PC1 ({variance['PC1']:.2f}%)",
        ylabel=f"PC2 ({variance['PC2']:.2f}%)",
        title=title,
    )
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax)


def _plot_differential(project_root: Path, output: Path, cohort: dict) -> dict:
    pca_data = {
        modality: _read_pca(project_root, modality, cohort)
        for modality in ("cpg", "occupancy", "wps")
    }
    dmr = _read_dmr(project_root)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.7), constrained_layout=True)
    for ax, modality, label in zip(
        axes.ravel()[:3],
        ("cpg", "occupancy", "wps"),
        ("CpG methylation PCA", "Panel occupancy PCA", "Panel WPS PCA"),
    ):
        coordinates, variance = pca_data[modality]
        _plot_pca(ax, coordinates, variance, label)

    ax = axes[1, 1]
    q_values = dmr["q_value"].clip(lower=np.finfo(float).tiny)
    significant = dmr["q_value"] < 0.05
    colors = np.where(
        significant & (dmr["mean_diff"] > 0), CASE_COLOR,
        np.where(significant & (dmr["mean_diff"] < 0), CONTROL_COLOR, "#BCC5CC"),
    )
    ax.scatter(dmr["mean_diff"], -np.log10(q_values), s=8, c=colors, alpha=0.65)
    ax.axhline(-np.log10(0.05), color=MUTED_COLOR, linestyle="--", linewidth=0.8)
    ax.axvline(0, color=MUTED_COLOR, linewidth=0.7)
    ax.set(
        xlabel="Mean difference (sALS - Control)",
        ylabel="-log10(q value)",
        title=f"Raw DMR output ({significant.sum():,} rows with q < 0.05)",
    )
    _style_axis(ax)
    fig.suptitle(
        "Observed CFTK differential outputs: five controls and five sALS samples",
        fontsize=14,
    )
    _save(fig, output)
    return {
        "modalities": ["cpg", "occupancy", "wps"],
        "pca_sample_count": 10,
        "dmr_rows": int(len(dmr)),
        "dmr_rows_q_lt_0_05": int(significant.sum()),
    }


def _stream_matrix_means(path: Path, samples: Iterable[str]) -> tuple[dict, int]:
    sample_list = list(samples)
    sums = {sample: 0.0 for sample in sample_list}
    counts = {sample: 0 for sample in sample_list}
    row_count = 0
    for chunk in pd.read_csv(_require_file(path), sep="\t", chunksize=CHUNK_ROWS):
        missing = set(sample_list) - set(chunk.columns)
        if missing:
            raise ValueError(f"Matrix {path.name} is missing samples: {sorted(missing)}")
        row_count += len(chunk)
        for sample in sample_list:
            values = pd.to_numeric(chunk[sample], errors="coerce").to_numpy(dtype=float)
            valid = np.isfinite(values)
            sums[sample] += float(values[valid].sum())
            counts[sample] += int(valid.sum())
    means = {
        sample: sums[sample] / counts[sample] if counts[sample] else float("nan")
        for sample in sample_list
    }
    if not all(np.isfinite(value) for value in means.values()):
        raise ValueError(f"Matrix {path.name} has no numeric values for one or more samples")
    return means, row_count


def _read_delfi_chromosome_means(project_root: Path, cohort: dict) -> tuple[pd.DataFrame, int]:
    delfi_dir = project_root / "results" / "4_fragmentomics" / "delfi"
    rows = []
    total_values = 0
    chromosomes = [f"chr{index}" for index in range(1, 23)]
    for sample in cohort["ordered_samples"]:
        path = _require_file(delfi_dir / f"{sample}_delfi.tsv")
        header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
        contig_column = "#contig" if "#contig" in header else "contig"
        ratio_column = "ratio_corrected" if "ratio_corrected" in header else "ratio"
        sums = {chromosome: 0.0 for chromosome in chromosomes}
        counts = {chromosome: 0 for chromosome in chromosomes}
        for chunk in pd.read_csv(
            path,
            sep="\t",
            usecols=[contig_column, ratio_column],
            chunksize=CHUNK_ROWS,
        ):
            chunk[ratio_column] = pd.to_numeric(chunk[ratio_column], errors="coerce")
            chunk = chunk.dropna(subset=[ratio_column])
            grouped = chunk.groupby(contig_column, observed=True)[ratio_column].agg(["sum", "count"])
            for chromosome, values in grouped.iterrows():
                if chromosome in sums:
                    sums[chromosome] += float(values["sum"])
                    counts[chromosome] += int(values["count"])
        for chromosome in chromosomes:
            if counts[chromosome]:
                rows.append({
                    "sample": sample,
                    "group": cohort["groups"][sample],
                    "chromosome": chromosome,
                    "ratio": sums[chromosome] / counts[chromosome],
                })
                total_values += counts[chromosome]
    return pd.DataFrame(rows), total_values


def _read_end_motifs(project_root: Path, cohort: dict) -> pd.DataFrame:
    motif_dir = project_root / "results" / "4_fragmentomics" / "end_motif"
    frames = []
    for sample in cohort["ordered_samples"]:
        path = _require_file(motif_dir / f"{sample}_4mer.tsv")
        frame = pd.read_csv(path, sep="\t", header=None, names=["motif", "frequency"])
        frame["frequency"] = pd.to_numeric(frame["frequency"], errors="coerce")
        frame["group"] = cohort["groups"][sample]
        frames.append(frame.dropna(subset=["frequency"]))
    return pd.concat(frames, ignore_index=True)


def _read_scope(project_root: Path) -> dict:
    path = _require_file(
        project_root / "results" / "4_fragmentomics" / "delfi" / "fragmentomics_scope.json"
    )
    payload = json.loads(path.read_text())
    resolved = payload.get("resolved_scope", payload)
    mode = resolved.get("mode")
    bins_count = resolved.get("bins_count")
    if mode != "panel" or not isinstance(bins_count, int) or bins_count <= 0:
        raise ValueError("Documentation export requires a positive-bin panel scope")
    return {
        "mode": mode,
        "bins_count": bins_count,
        "target_sha256": resolved.get("target_sha256"),
        "note": PANEL_SCOPE_NOTE,
    }


def _plot_sample_means(ax, means: dict, cohort: dict, title: str, ylabel: str) -> None:
    ordered = cohort["ordered_samples"]
    x = np.arange(len(ordered))
    values = np.array([means[sample] for sample in ordered], dtype=float)
    colors = [CONTROL_COLOR if cohort["groups"][sample] == "Control" else CASE_COLOR for sample in ordered]
    ax.scatter(x, values, c=colors, s=38, edgecolor="white", linewidth=0.5)
    for samples, color in (
        (cohort["control_samples"], CONTROL_COLOR),
        (cohort["case_samples"], CASE_COLOR),
    ):
        indexes = [ordered.index(sample) for sample in samples]
        group_values = [means[sample] for sample in samples]
        ax.hlines(np.mean(group_values), min(indexes) - 0.35, max(indexes) + 0.35, color=color, linewidth=2)
    ax.set(
        xticks=x,
        xticklabels=[cohort["aliases"][sample] for sample in ordered],
        ylabel=ylabel,
        title=title,
    )
    ax.tick_params(axis="x", rotation=55, labelsize=7)
    _style_axis(ax, grid_axis="y")


def _plot_fragmentomics(project_root: Path, output: Path, cohort: dict) -> dict:
    fragment_root = project_root / "results" / "4_fragmentomics"
    occupancy_means, occupancy_rows = _stream_matrix_means(
        fragment_root / "occupancy" / "occupancy_matrix.tsv", cohort["ordered_samples"]
    )
    wps_means, wps_rows = _stream_matrix_means(
        fragment_root / "wps" / "wps_matrix.tsv", cohort["ordered_samples"]
    )
    delfi, delfi_values = _read_delfi_chromosome_means(project_root, cohort)
    motifs = _read_end_motifs(project_root, cohort)
    scope = _read_scope(project_root)
    cleavage_files = list((fragment_root / "cleavage").glob("*.bw")) if (fragment_root / "cleavage").is_dir() else []

    fig, axes = plt.subplots(2, 3, figsize=(12.1, 7.6), constrained_layout=True)
    _plot_sample_means(
        axes[0, 0], occupancy_means, cohort,
        f"Occupancy matrix ({occupancy_rows:,} panel regions)", "Mean occupancy",
    )
    _plot_sample_means(
        axes[0, 1], wps_means, cohort,
        f"WPS matrix ({wps_rows:,} panel regions)", "Mean WPS",
    )

    ax = axes[0, 2]
    chromosome_order = [f"chr{index}" for index in range(1, 23)]
    for group, color in (("Control", CONTROL_COLOR), ("sALS", CASE_COLOR)):
        group_means = (
            delfi[delfi["group"] == group]
            .groupby("chromosome", observed=True)["ratio"]
            .mean()
            .reindex(chromosome_order)
        )
        ax.plot(range(1, 23), group_means, color=color, marker="o", markersize=3, label=group)
    ax.set(
        xticks=[1, 5, 10, 15, 20, 22],
        xlabel="Chromosome",
        ylabel="Mean corrected short/long ratio",
        title="DELFI-style panel-overlap summary",
    )
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax)

    ax = axes[1, 0]
    motif_group = motifs.groupby(["motif", "group"], observed=True)["frequency"].mean().unstack()
    top = motif_group.mean(axis=1).nlargest(8).index
    motif_group = motif_group.loc[top]
    positions = np.arange(len(motif_group))
    width = 0.38
    ax.bar(positions - width / 2, 100 * motif_group["Control"], width, color=CONTROL_COLOR, label="Control")
    ax.bar(positions + width / 2, 100 * motif_group["sALS"], width, color=CASE_COLOR, label="sALS")
    ax.set(
        xticks=positions,
        xticklabels=motif_group.index,
        ylabel="Mean frequency (%)",
        title="Observed top 4-mer end motifs",
    )
    ax.tick_params(axis="x", rotation=45)
    ax.legend(frameon=False, fontsize=8)
    _style_axis(ax, grid_axis="y")

    ax = axes[1, 1]
    ax.axis("off")
    if cleavage_files:
        heading = "Cleavage outputs available"
        detail = f"{len(cleavage_files)} bigWig files found\nNot rendered by this public exporter"
        color = ACCENT_COLOR
    else:
        heading = "Cleavage not run"
        detail = "No cleavage profile was produced\nfor this technical example"
        color = MUTED_COLOR
    ax.text(0.5, 0.62, heading, ha="center", va="center", fontsize=13, fontweight="bold", color=color)
    ax.text(0.5, 0.38, detail, ha="center", va="center", fontsize=10, color=MUTED_COLOR, linespacing=1.5)

    ax = axes[1, 2]
    ax.axis("off")
    ax.text(0.03, 0.86, "Assay scope", fontsize=13, fontweight="bold", color="#1F2933")
    scope_lines = [
        f"Mode: {scope['mode']}",
        f"Clipped bins: {scope['bins_count']:,}" if scope["bins_count"] is not None else "Clipped bins: unavailable",
        "WPS, occupancy, and DELFI are",
        "panel-overlap summaries, not",
        "genome-wide measurements.",
    ]
    ax.text(0.03, 0.72, "\n".join(scope_lines), va="top", fontsize=10, color=MUTED_COLOR, linespacing=1.5)

    fig.suptitle(
        "Observed fragmentomics outputs: five controls and five sALS samples",
        fontsize=14,
    )
    _save(fig, output)
    return {
        "scope": scope,
        "occupancy_matrix_rows": int(occupancy_rows),
        "wps_matrix_rows": int(wps_rows),
        "delfi_numeric_values": int(delfi_values),
        "delfi_samples": 10,
        "end_motif_samples": 10,
        "cleavage": "available" if cleavage_files else "not_run",
    }


def _read_mesa(project_root: Path, cohort: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    mesa_root = project_root / "results" / "5_mesa"
    performance = pd.read_csv(_require_file(mesa_root / "modality_performance.tsv"), sep="\t", index_col=0)
    if "best_roc_auc_mean" not in performance.columns:
        raise ValueError("modality_performance.tsv is missing best_roc_auc_mean")
    performance["best_roc_auc_mean"] = pd.to_numeric(performance["best_roc_auc_mean"], errors="raise")
    if not {"cpg", "occupancy", "wps"}.issubset(performance.index):
        raise ValueError("modality_performance.tsv is missing a required modality")

    predictions = pd.read_csv(_require_file(mesa_root / "loocv_predictions.tsv"), sep="\t")
    required = {"sample_id", "y_true", "cpg", "occupancy", "wps", "Multimodal"}
    if not required.issubset(predictions.columns):
        raise ValueError(f"Unexpected LOOCV schema: {list(predictions.columns)}")
    if set(predictions["sample_id"]) != set(cohort["ordered_samples"]):
        raise ValueError("LOOCV sample set does not match samples.tsv")
    predictions["y_true"] = pd.to_numeric(predictions["y_true"], errors="raise")
    expected_labels = predictions["sample_id"].map(
        {sample: 0 for sample in cohort["control_samples"]}
        | {sample: 1 for sample in cohort["case_samples"]}
    )
    if not predictions["y_true"].equals(expected_labels):
        raise ValueError("LOOCV y_true labels do not match control/case roles")
    for modality in ("cpg", "occupancy", "wps", "Multimodal"):
        predictions[modality] = pd.to_numeric(predictions[modality], errors="raise")
        if not np.isfinite(predictions[modality]).all():
            raise ValueError(f"LOOCV predictions contain non-finite {modality} values")
    predictions["alias"] = predictions["sample_id"].map(cohort["aliases"])
    predictions["group"] = predictions["sample_id"].map(cohort["groups"])
    predictions["order"] = predictions["sample_id"].map(
        {sample: index for index, sample in enumerate(cohort["ordered_samples"])}
    )
    return performance, predictions.sort_values("order")


def _plot_mesa(project_root: Path, output: Path, cohort: dict) -> dict:
    performance, predictions = _read_mesa(project_root, cohort)
    modalities = ["cpg", "occupancy", "wps", "Multimodal"]
    colors = ["#8E6C8A", CONTROL_COLOR, ACCENT_COLOR, CASE_COLOR]
    aucs = {
        modality: float(roc_auc_score(predictions["y_true"], predictions[modality]))
        for modality in modalities
    }

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.4), constrained_layout=True)
    ax = axes[0, 0]
    screening = performance["best_roc_auc_mean"].reindex(["cpg", "occupancy", "wps"])
    ax.bar(screening.index, screening.values, color=colors[:3])
    ax.set(ylim=(0, 1.05), ylabel="Mean ROC AUC", title="MESA modality screening (internal CV)")
    for index, value in enumerate(screening.values):
        ax.text(index, value + 0.025, f"{value:.3f}", ha="center", fontsize=8)
    _style_axis(ax, grid_axis="y")

    ax = axes[0, 1]
    for modality, color in zip(modalities, colors):
        fpr, tpr, _ = roc_curve(predictions["y_true"], predictions[modality])
        width = 2.0 if modality == "Multimodal" else 1.4
        ax.plot(fpr, tpr, color=color, linewidth=width, label=f"{modality}: {aucs[modality]:.2f}")
    ax.plot([0, 1], [0, 1], color=MUTED_COLOR, linestyle="--", linewidth=0.8)
    ax.set(xlabel="False-positive rate", ylabel="True-positive rate", title="Observed LOOCV ROC (n=10)")
    ax.legend(frameon=False, fontsize=7, title="AUC")
    _style_axis(ax)

    ax = axes[1, 0]
    x = np.arange(len(predictions))
    width = 0.17
    for offset, modality, color in zip((-1.5, -0.5, 0.5, 1.5), modalities, colors):
        ax.scatter(x + offset * width, predictions[modality], s=27, color=color, label=modality)
    ax.set(
        xticks=x,
        xticklabels=predictions["alias"],
        ylim=(-0.02, 1.02),
        ylabel="Predicted probability",
        title="Aliased per-sample LOOCV predictions",
    )
    ax.tick_params(axis="x", rotation=55, labelsize=7)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    _style_axis(ax, grid_axis="y")

    ax = axes[1, 1]
    correlation = predictions[modalities].corr(method="spearman")
    image = ax.imshow(correlation, cmap="YlGnBu", vmin=-1, vmax=1)
    ax.set(
        xticks=np.arange(len(modalities)),
        yticks=np.arange(len(modalities)),
        xticklabels=modalities,
        yticklabels=modalities,
        title="Spearman correlation of LOOCV scores",
    )
    ax.tick_params(axis="x", rotation=35)
    for row in range(len(modalities)):
        for column in range(len(modalities)):
            ax.text(column, row, f"{correlation.iloc[row, column]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Spearman r")

    fig.suptitle(
        "Observed MESA technical outputs - not biological or clinical validation",
        fontsize=14,
    )
    _save(fig, output)
    return {
        "sample_count": 10,
        "screening_best_auc": {name: float(value) for name, value in screening.items()},
        "loocv_auc": aucs,
        "interpretation": "Technical workflow output only; not model performance validation.",
    }


def _plot_report(project_root: Path, output: Path, summary: dict) -> dict:
    report_path = _require_file(project_root / "results" / "report" / "report.html")
    report_text = report_path.read_text(errors="replace")
    section_checks = {
        "Processing": "Processing" in report_text,
        "cfDNA QC": "cfDNA QC" in report_text,
        "Differential / DMR": "Differential" in report_text and "DMR" in report_text,
        "Fragmentomics": "DELFI" in report_text and "WPS" in report_text,
        "MESA": "MESA" in report_text,
    }

    fig = plt.figure(figsize=(11.2, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[0.78, 1.22])
    header = fig.add_subplot(grid[0, :])
    header.axis("off")
    header.text(0.02, 0.82, "CFTK whole-workflow report", fontsize=20, fontweight="bold", color="#1F2933")
    header.text(
        0.02,
        0.57,
        "Sanitized static overview from the completed five-control/five-sALS technical run",
        fontsize=11,
        color=MUTED_COLOR,
    )
    header.text(0.02, 0.20, "10", fontsize=24, fontweight="bold", color=CONTROL_COLOR)
    header.text(0.075, 0.22, "samples", fontsize=10, color=MUTED_COLOR)
    header.text(0.25, 0.20, "5 + 5", fontsize=24, fontweight="bold", color=CASE_COLOR)
    header.text(0.34, 0.22, "Control + sALS", fontsize=10, color=MUTED_COLOR)
    header.text(0.57, 0.20, "panel", fontsize=24, fontweight="bold", color=ACCENT_COLOR)
    header.text(0.67, 0.22, "fragmentomics scope", fontsize=10, color=MUTED_COLOR)

    ax = fig.add_subplot(grid[1, 0])
    ax.axis("off")
    ax.text(0.02, 0.95, "Report sections discovered", va="top", fontsize=12, fontweight="bold")
    y = 0.80
    for label, present in section_checks.items():
        marker = "present" if present else "not found"
        color = ACCENT_COLOR if present else CASE_COLOR
        ax.text(0.03, y, marker, color=color, fontsize=9, fontweight="bold")
        ax.text(0.28, y, label, color="#1F2933", fontsize=9)
        y -= 0.14

    ax = fig.add_subplot(grid[1, 1])
    ax.axis("off")
    differential = summary["differential"]
    ax.text(0.02, 0.95, "Comparative outputs", va="top", fontsize=12, fontweight="bold")
    lines = [
        "CpG PCA: 10 samples",
        "Occupancy PCA: 10 samples",
        "WPS PCA: 10 samples",
        f"DMR output: {differential['dmr_rows']:,} rows",
        f"DMR q < 0.05: {differential['dmr_rows_q_lt_0_05']:,} rows",
    ]
    ax.text(0.03, 0.80, "\n".join(lines), va="top", fontsize=9.5, color="#1F2933", linespacing=1.7)

    ax = fig.add_subplot(grid[1, 2])
    ax.axis("off")
    fragmentomics = summary["fragmentomics"]
    ax.text(0.02, 0.95, "Downstream outputs", va="top", fontsize=12, fontweight="bold")
    downstream_lines = [
        f"Occupancy: {fragmentomics['occupancy_matrix_rows']:,} rows",
        f"WPS: {fragmentomics['wps_matrix_rows']:,} rows",
        "DELFI-style: 10 samples",
        "End motif: 10 samples",
        f"Cleavage: {fragmentomics['cleavage'].replace('_', ' ')}",
        "MESA LOOCV: 10 samples",
    ]
    ax.text(0.03, 0.80, "\n".join(downstream_lines), va="top", fontsize=9.5, color="#1F2933", linespacing=1.7)
    fig.text(
        0.5,
        0.015,
        "Technical workflow evidence only; this cohort is not biological or clinical validation.",
        ha="center",
        fontsize=9,
        color=MUTED_COLOR,
    )
    _save(fig, output)
    return {
        "html_present": True,
        "html_size_bytes": int(report_path.stat().st_size),
        "sections_discovered": section_checks,
    }


def export_documentation_evidence(project_root: Path, output_dir: Path) -> dict:
    project_root = project_root.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cohort = _read_cohort(project_root)
    summary = {
        "schema_version": 1,
        "cohort": {
            "groups": {"Control": 5, "sALS": 5},
            "sample_aliases": [cohort["aliases"][sample] for sample in cohort["ordered_samples"]],
            "purpose": "Technical workflow example; not biological or clinical validation.",
        },
    }
    summary["differential"] = _plot_differential(
        project_root, output_dir / OUTPUT_NAMES["differential"], cohort
    )
    summary["fragmentomics"] = _plot_fragmentomics(
        project_root, output_dir / OUTPUT_NAMES["fragmentomics"], cohort
    )
    summary["mesa"] = _plot_mesa(project_root, output_dir / OUTPUT_NAMES["mesa"], cohort)
    summary["report"] = _plot_report(project_root, output_dir / OUTPUT_NAMES["report"], summary)
    (output_dir / OUTPUT_NAMES["summary"]).write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    return summary


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export sanitized documentation figures from a completed CFTK 5+5 run."
    )
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    summary = export_documentation_evidence(args.project_root, args.output_dir)
    print(
        f"Exported {len(OUTPUT_NAMES) - 1} figures and one sanitized summary "
        f"for {sum(summary['cohort']['groups'].values())} samples."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
