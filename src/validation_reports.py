"""Create private tables and sanitized figures from a CFTK run manifest."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


EVIDENCE_FILENAMES = (
    "workflow_artifact_inventory.tsv",
    "workflow_stage_evidence.tsv",
    "workflow_command_evidence.tsv",
    "workflow_stage_evidence.png",
    "workflow_resource_plan.png",
    "workflow_qc_overview.png",
    "workflow_validation_summary.json",
)


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _artifact_rows(manifest):
    rows = []
    for stage in manifest.get("stages", []):
        for artifact in stage.get("expected", []):
            path = Path(artifact["path"])
            exists = path.exists()
            nonempty_required = bool(artifact.get("nonempty", True))
            nonempty = exists and (not path.is_file() or path.stat().st_size > 0)
            rows.append({
                "stage": stage["id"],
                "stage_status": stage.get("status", "unknown"),
                "role": artifact.get("role", "output"),
                "description": artifact.get("description", ""),
                "path": str(path),
                "required": bool(artifact.get("required", True)),
                "nonempty_required": nonempty_required,
                "exists": exists,
                "nonempty": nonempty,
                "valid": exists and (not nonempty_required or nonempty),
            })
    return rows


def _stage_rows(manifest, artifacts):
    resources = {
        item["stage"]: item
        for item in manifest.get("resource_plan", {}).get("stages", [])
    }
    rows = []
    for stage in manifest.get("stages", []):
        stage_artifacts = [row for row in artifacts if row["stage"] == stage["id"]]
        required = [row for row in stage_artifacts if row["required"]]
        missing = [row for row in required if not row["valid"]]
        resource = resources.get(stage["id"], {})
        rows.append({
            "stage": stage["id"],
            "name": stage.get("name", ""),
            "status": stage.get("status", "unknown"),
            "command": stage.get("command", ""),
            "required_outputs": sum(row["role"] != "figure" for row in required),
            "required_figures": sum(row["role"] == "figure" for row in required),
            "missing_required": len(missing),
            "total_core_budget": resource.get("total_core_budget", ""),
            "concurrent_samples": resource.get("concurrent_samples", ""),
            "threads_per_sample": resource.get("threads_per_sample", ""),
            "estimated_peak_threads": resource.get("estimated_peak_threads", ""),
        })
    return rows


def _command_rows(manifest_path):
    """Read the per-attempt command mirror without altering the raw ledger."""
    ledger = Path(manifest_path).parent / "commands.jsonl"
    rows = []
    if not ledger.is_file():
        return rows, ledger
    for line_number, line in enumerate(ledger.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        record = json.loads(line)
        rows.append({
            "line": line_number,
            "event": record.get("event", ""),
            "command_id": record.get("command_id", ""),
            "label": record.get("label", ""),
            "command": record.get("command", ""),
            "returncode": "" if record.get("returncode") is None else record["returncode"],
            "timestamp": record.get("timestamp", ""),
        })
    return rows, ledger


def _command_summary(rows, ledger):
    starts = {row["command_id"] for row in rows if row["event"] == "start"}
    finishes = {row["command_id"] for row in rows if row["event"] == "finish"}
    nonzero = [
        row for row in rows
        if row["event"] == "finish" and row["returncode"] not in ("", 0)
    ]
    return {
        "path": str(ledger),
        "exists": ledger.is_file(),
        "records": len(rows),
        "starts": len(starts),
        "finishes": len(finishes),
        "unfinished": len(starts - finishes),
        "nonzero_finishes": len(nonzero),
    }


def _fragmentomics_scope(manifest):
    """Include executed stage-local scope details in machine-readable evidence."""

    scope = dict(manifest.get("fragmentomics_scope") or {})
    for stage in manifest.get("stages", []):
        for artifact in stage.get("expected", []):
            if Path(artifact.get("path", "")).name != "fragmentomics_scope.json":
                continue
            path = Path(artifact["path"])
            if not path.is_file():
                continue
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError, json.JSONDecodeError):
                continue
            resolved = payload.get("resolved_scope", payload)
            if isinstance(resolved, dict):
                scope.update({key: value for key, value in resolved.items() if value is not None})
    return scope


def _write_stage_figure(rows, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [row["stage"] for row in rows]
    outputs = [row["required_outputs"] for row in rows]
    figures = [row["required_figures"] for row in rows]
    colors = {
        "complete": "#2a9d8f", "resumed": "#2a9d8f", "adopted": "#2a9d8f",
        "planned": "#457b9d", "skipped": "#6c757d", "failed": "#e76f51",
        "pending": "#e9c46a", "running": "#e9c46a",
    }
    y = list(range(len(rows)))
    fig, axis = plt.subplots(figsize=(9, max(4.5, len(rows) * 0.62)))
    axis.barh(y, outputs, label="Outputs/reports", color="#457b9d")
    axis.barh(y, figures, left=outputs, label="Figures", color="#e9c46a")
    for index, row in enumerate(rows):
        axis.text(
            outputs[index] + figures[index] + 0.25,
            index,
            f"{row['status']} | missing {row['missing_required']}",
            va="center",
            color=colors.get(row["status"], "#333333"),
            fontsize=8,
        )
    axis.set_yticks(y, labels)
    if rows:
        axis.invert_yaxis()
    else:
        axis.text(
            0.5, 0.5, "No stage plan was recorded.", ha="center", va="center",
            transform=axis.transAxes,
        )
    axis.set_xlabel("Required artifacts")
    axis.set_title("CFTK workflow stage evidence")
    axis.grid(axis="x", color="#dddddd", linewidth=0.6)
    if rows:
        axis.legend(
            frameon=False,
            loc="lower right",
            bbox_to_anchor=(1.0, 1.02),
            ncol=2,
        )
    # Fixed margins keep sparse/partial manifests renderable without a
    # constrained-layout warning being promoted to a reporting failure.
    fig.subplots_adjust(left=0.2, right=0.98, top=0.82, bottom=0.14)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _write_resource_figure(manifest, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [
        row for row in manifest.get("resource_plan", {}).get("stages", [])
        if row.get("applicable", True)
    ]
    if not rows:
        return False
    labels = [row["stage"] for row in rows]
    budgets = [row["total_core_budget"] for row in rows]
    peaks = [row["estimated_peak_threads"] for row in rows]
    x = list(range(len(rows)))
    width = 0.36
    fig, axis = plt.subplots(figsize=(9, 4.8))
    axis.bar([value - width / 2 for value in x], budgets, width, label="Total budget", color="#457b9d")
    axis.bar([value + width / 2 for value in x], peaks, width, label="Estimated peak", color="#2a9d8f")
    axis.set_xticks(x, labels, rotation=35, ha="right")
    axis.set_ylabel("CPU threads")
    axis.set_title("CFTK recorded CPU plan")
    axis.grid(axis="y", color="#dddddd", linewidth=0.6)
    axis.legend(frameon=False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return True


def _find_qc_summary(manifest):
    for stage in manifest.get("stages", []):
        if stage.get("id") != "qc.0":
            continue
        for artifact in stage.get("expected", []):
            path = Path(artifact["path"])
            if path.name == "qc_summary.tsv" and path.is_file():
                return path
    return None


def _qc_display_labels(frame):
    """Use group/order labels so rendered QC figures omit sample identifiers."""
    groups = (
        frame["group"].fillna("Sample").astype(str).tolist()
        if "group" in frame
        else ["Sample"] * len(frame)
    )
    group_order = list(dict.fromkeys(groups))
    palette = ["#457b9d", "#e76f51", "#2a9d8f", "#e9c46a"]
    group_colors = {
        group: palette[index % len(palette)] for index, group in enumerate(group_order)
    }
    counters = Counter()
    labels = []
    for group in groups:
        counters[group] += 1
        labels.append(f"{group} {counters[group]}")
    return labels, [group_colors[group] for group in groups]


def _write_qc_figure(qc_path, path):
    if qc_path is None:
        return False
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    frame = pd.read_csv(qc_path, sep="\t")
    if frame.empty:
        return False
    metrics = [
        ("flagstat_mapped_pct", "Mapped reads (%)"),
        ("markdup_dup_pct", "Duplicates (%)"),
        ("cpg_mean_depth", "Mean CpG depth"),
        ("cpg_covered_sites", "Covered CpG sites"),
        ("cpg_global_meth_pct", "Global CpG methylation (%)"),
        ("median_frag_len", "Median fragment length (bp)"),
    ]
    labels, colors = _qc_display_labels(frame)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7.5), constrained_layout=True)
    for axis, (column, title) in zip(axes.flat, metrics):
        if column not in frame or frame[column].isna().all():
            axis.text(0.5, 0.5, "Not recorded", ha="center", va="center", transform=axis.transAxes)
        else:
            axis.bar(range(len(frame)), frame[column], color=colors)
        axis.set_xticks(range(len(frame)), labels, rotation=45, ha="right", fontsize=7)
        axis.set_title(title, fontsize=9)
        axis.grid(axis="y", color="#dddddd", linewidth=0.6)
    fig.suptitle("CFTK cohort QC overview (sample identifiers removed)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return True


def summarize(manifest_path, output_dir):
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    artifacts = _artifact_rows(manifest)
    stages = _stage_rows(manifest, artifacts)
    commands, command_ledger = _command_rows(manifest_path)
    command_summary = _command_summary(commands, command_ledger)
    output_dir = Path(output_dir)
    _write_tsv(
        output_dir / "workflow_artifact_inventory.tsv",
        artifacts,
        [
            "stage", "stage_status", "role", "description", "path", "required",
            "nonempty_required", "exists", "nonempty", "valid",
        ],
    )
    _write_tsv(
        output_dir / "workflow_stage_evidence.tsv",
        stages,
        [
            "stage", "name", "status", "command", "required_outputs",
            "required_figures", "missing_required", "total_core_budget",
            "concurrent_samples", "threads_per_sample", "estimated_peak_threads",
        ],
    )
    _write_tsv(
        output_dir / "workflow_command_evidence.tsv",
        commands,
        ["line", "event", "command_id", "label", "command", "returncode", "timestamp"],
    )
    _write_stage_figure(stages, output_dir / "workflow_stage_evidence.png")
    resource_figure = _write_resource_figure(
        manifest, output_dir / "workflow_resource_plan.png"
    )
    qc_figure = _write_qc_figure(
        _find_qc_summary(manifest), output_dir / "workflow_qc_overview.png"
    )
    missing = [
        row for row in artifacts
        if row["required"] and not row["valid"]
    ]
    summary = {
        "generated_at": _utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "manifest": str(Path(manifest_path)),
        "run_id": manifest.get("run_id"),
        "run_status": manifest.get("status"),
        "stage_status_counts": dict(Counter(row["status"] for row in stages)),
        "required_artifacts": sum(row["required"] for row in artifacts),
        "missing_required_artifacts": len(missing),
        "resource_figure_written": resource_figure,
        "qc_figure_written": qc_figure,
        "command_ledger": command_summary,
        "note": "Tables may contain private paths; review figures before public documentation use.",
    }
    if manifest.get("fragmentomics_scope") or any(
        Path(artifact.get("path", "")).name == "fragmentomics_scope.json"
        for stage in manifest.get("stages", [])
        for artifact in stage.get("expected", [])
    ):
        summary["fragmentomics_scope"] = _fragmentomics_scope(manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary["files"] = [
        filename for filename in EVIDENCE_FILENAMES
        if filename == "workflow_validation_summary.json" or (output_dir / filename).is_file()
    ]
    (output_dir / "workflow_validation_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    summarize(args.manifest, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
