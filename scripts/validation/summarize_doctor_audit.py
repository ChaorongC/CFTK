#!/usr/bin/env python3
"""Turn a CFTK doctor JSON report into beginner-readable audit artifacts.

The input report and sample sheet may contain local paths. The generated
machine-readable files are intended for a private validation directory; only
the aggregate figure and narrative should be copied into public documentation.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_samples(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or not {"sample", "group", "role"}.issubset(rows[0]):
        raise ValueError("sample sheet must contain sample, group, and role columns")
    return {row["sample"]: row for row in rows}


def _sample_check_rows(report: dict, samples: dict[str, dict[str, str]]):
    rows = []
    by_sample: dict[str, dict[str, dict]] = defaultdict(dict)
    for check in report.get("checks", []):
        check_id = str(check.get("id", ""))
        prefix = "input.bam."
        if not check_id.startswith(prefix):
            continue
        remainder = check_id[len(prefix):]
        if "." not in remainder:
            continue
        sample, check_name = remainder.split(".", 1)
        by_sample[sample][check_name] = check

    for sample, metadata in samples.items():
        checks = by_sample.get(sample, {})
        statuses = {name: item.get("status", "MISSING") for name, item in checks.items()}
        if "FAIL" in statuses.values() or "MISSING" in statuses.values():
            overall = "FAIL"
        elif "WARN" in statuses.values():
            overall = "WARN"
        else:
            overall = "PASS"
        for check_name in ("index", "dictionary", "read_group", "duplicates", "sorting"):
            item = checks.get(check_name, {})
            rows.append(
                {
                    "sample": sample,
                    "group": metadata.get("group", ""),
                    "role": metadata.get("role", ""),
                    "overall_status": overall,
                    "check": check_name,
                    "status": item.get("status", "MISSING"),
                    "summary": item.get("summary", "No doctor check was reported."),
                    "remedy": item.get("remedy", ""),
                }
            )
    return rows


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _write_figure(rows: list[dict], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    checks = ["index", "dictionary", "read_group", "duplicates", "sorting"]
    groups = sorted({row["group"] for row in rows})
    statuses = ("PASS", "WARN", "FAIL", "MISSING")
    colors = {"PASS": "#2a9d8f", "WARN": "#e9c46a", "FAIL": "#e76f51", "MISSING": "#6c757d"}
    hatches = {"PASS": "", "WARN": "//", "FAIL": "xx", "MISSING": ".."}
    fig, axes = plt.subplots(1, len(checks), figsize=(13, 4.5), sharey=True)
    for axis, check in zip(axes, checks):
        for index, group in enumerate(groups):
            counts = Counter(
                row["status"]
                for row in rows
                if row["check"] == check and row["group"] == group
            )
            bottom = 0
            for status in statuses:
                value = counts.get(status, 0)
                if value:
                    axis.bar(
                        index,
                        value,
                        bottom=bottom,
                        color=colors[status],
                        edgecolor="#333333",
                        linewidth=0.4,
                        hatch=hatches[status],
                        width=0.72,
                    )
                    bottom += value
        axis.set_title(check.replace("_", " ").title(), fontsize=9)
        axis.set_xticks(range(len(groups)), groups, rotation=35, ha="right", fontsize=8)
        axis.grid(axis="y", color="#dddddd", linewidth=0.6)
    axes[0].set_ylabel("Samples")
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=colors[s], hatch=hatches[s], edgecolor="#333333")
        for s in statuses
    ]
    fig.suptitle("CFTK doctor readiness audit: historical ALS BAM cohort", y=0.98)
    fig.legend(handles, statuses, loc="upper center", bbox_to_anchor=(0.5, 0.91), ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.82))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def summarize(doctor_json: Path, sample_sheet: Path, output_dir: Path) -> dict:
    report = json.loads(doctor_json.read_text(encoding="utf-8"))
    samples = _read_samples(sample_sheet)
    rows = _sample_check_rows(report, samples)
    if len(samples) == 0 or len(rows) != len(samples) * 5:
        raise ValueError("doctor report does not contain the expected per-sample checks")

    fieldnames = [
        "sample", "group", "role", "overall_status", "check", "status",
        "summary", "remedy",
    ]
    _write_tsv(output_dir / "cohort_readiness_checks.tsv", rows, fieldnames)

    sample_rows = []
    for sample in samples:
        sample_rows.append({
            "sample": sample,
            "group": next(row["group"] for row in rows if row["sample"] == sample),
            "role": next(row["role"] for row in rows if row["sample"] == sample),
            "overall_status": next(row["overall_status"] for row in rows if row["sample"] == sample),
            **{
                check: next(row["status"] for row in rows if row["sample"] == sample and row["check"] == check)
                for check in ("index", "dictionary", "read_group", "duplicates", "sorting")
            },
        })
    _write_tsv(
        output_dir / "cohort_readiness.tsv",
        sample_rows,
        ["sample", "group", "role", "overall_status", "index", "dictionary", "read_group", "duplicates", "sorting"],
    )

    summary = {
        "generated_at": _utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "doctor_report": str(doctor_json),
        "sample_sheet": str(sample_sheet),
        "doctor_status": report.get("status"),
        "doctor_exit_code": report.get("exit_code"),
        "sample_count": len(samples),
        "group_counts": dict(Counter(row["group"] for row in sample_rows)),
        "overall_status_counts": dict(Counter(row["overall_status"] for row in sample_rows)),
        "check_status_counts": {
            check: dict(Counter(row[check] for row in sample_rows))
            for check in ("index", "dictionary", "read_group", "duplicates", "sorting")
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "cohort_readiness_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_figure(rows, output_dir / "cohort_readiness.png")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("doctor_json", type=Path)
    parser.add_argument("sample_sheet", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    summarize(args.doctor_json, args.sample_sheet, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
