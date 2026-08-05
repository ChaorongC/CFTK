#!/usr/bin/env python3
"""Compare duplicate marking and downstream metrics from two BAM outputs.

The comparison is intentionally read-only. It treats the input alignment BAM
as the common universe, compares primary alignment duplicate flags by
``(query_name, mate)`` keys, and optionally compares Picard HsMetrics and
MethylDackel CpG bedGraph outputs. Thresholds are accepted as data, not baked
into the raw measurements, so a report cannot silently turn an arbitrary
tolerance into a scientific claim.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pysam


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _mate_key(read) -> tuple[str, int]:
    if read.is_read1:
        mate = 1
    elif read.is_read2:
        mate = 2
    else:
        mate = 0
    return read.query_name, mate


def summarize_bam(path: Path) -> dict:
    counts = Counter()
    duplicate_keys = set()
    duplicate_names = set()
    with pysam.AlignmentFile(path, "rb") as bam:
        references = tuple(zip(bam.references, bam.lengths))
        sort_order = bam.header.get("HD", {}).get("SO", "")
        for read in bam.fetch(until_eof=True):
            counts["records"] += 1
            if read.is_unmapped:
                counts["unmapped"] += 1
            else:
                counts["mapped"] += 1
            if read.is_secondary:
                counts["secondary"] += 1
                continue
            if read.is_supplementary:
                counts["supplementary"] += 1
                continue
            counts["primary"] += 1
            key = _mate_key(read)
            if read.is_paired:
                counts["paired"] += 1
            if read.is_duplicate:
                counts["duplicate_primary"] += 1
                duplicate_keys.add(key)
                duplicate_names.add(read.query_name)
            if read.is_proper_pair:
                counts["proper_pair"] += 1
        header_signature = hashlib.sha256(
            "\n".join(f"{name}\t{length}" for name, length in references).encode()
        ).hexdigest()
    counts["duplicate_primary_fraction"] = (
        counts["duplicate_primary"] / counts["primary"] if counts["primary"] else None
    )
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "references": len(references),
        "header_reference_signature": header_signature,
        "sort_order": sort_order,
        "counts": dict(counts),
        "duplicate_keys": duplicate_keys,
        "duplicate_names": duplicate_names,
    }


def _parse_picard_metrics(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header_index = next(
        (index for index, line in enumerate(lines) if line.startswith("BAIT_SET\t")),
        None,
    )
    if header_index is None or header_index + 1 >= len(lines):
        raise ValueError(f"Picard HsMetrics header not found in {path}")
    headers = lines[header_index].split("\t")
    values = lines[header_index + 1].split("\t")
    row = dict(zip(headers, values))
    numeric = {}
    for key, value in row.items():
        try:
            numeric[key] = float(value)
        except (TypeError, ValueError):
            numeric[key] = value
    return numeric


def _parse_sambamba_metrics(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"found\s+(\d[\d,]*)\s+duplicates", text)
    return {
        "reported_duplicates": int(match.group(1).replace(",", "")) if match else None,
        "metrics_path": str(path),
    }


def _parse_bedgraph(path: Path) -> dict[tuple[str, int, int], tuple[float, float, float]]:
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"Expected six BEDGRAPH columns in {path}: {line!r}")
            key = (fields[0], int(fields[1]), int(fields[2]))
            values[key] = (float(fields[3]), float(fields[4]), float(fields[5]))
    return values


def _compare_bedgraphs(sambamba_path: Path | None, picard_path: Path | None) -> dict | None:
    if not sambamba_path or not picard_path:
        return None
    sambamba = _parse_bedgraph(sambamba_path)
    picard = _parse_bedgraph(picard_path)
    shared = sorted(set(sambamba) & set(picard))
    if not shared:
        return {"sambamba_loci": len(sambamba), "picard_loci": len(picard), "shared_loci": 0}
    deltas = [abs(sambamba[key][0] - picard[key][0]) for key in shared]
    sambamba_methylated = sum(sambamba[key][1] for key in shared)
    sambamba_total = sambamba_methylated + sum(sambamba[key][2] for key in shared)
    picard_methylated = sum(picard[key][1] for key in shared)
    picard_total = picard_methylated + sum(picard[key][2] for key in shared)
    return {
        "sambamba_loci": len(sambamba),
        "picard_loci": len(picard),
        "shared_loci": len(shared),
        "shared_fraction_of_union": len(shared) / len(set(sambamba) | set(picard)),
        "mean_abs_methylation_percent_delta": sum(deltas) / len(deltas),
        "max_abs_methylation_percent_delta": max(deltas),
        "sambamba_weighted_methylation_percent": 100 * sambamba_methylated / sambamba_total if sambamba_total else None,
        "picard_weighted_methylation_percent": 100 * picard_methylated / picard_total if picard_total else None,
        "weighted_methylation_percent_absolute_delta": abs(
            100 * sambamba_methylated / sambamba_total - 100 * picard_methylated / picard_total
        ) if sambamba_total and picard_total else None,
    }


def compare_sample(
    *,
    sample: str,
    group: str,
    input_bam: Path,
    sambamba_bam: Path,
    picard_bam: Path,
    sambamba_metrics: Path | None = None,
    picard_hs_metrics: Path | None = None,
    sambamba_hs_metrics: Path | None = None,
    sambamba_bedgraph: Path | None = None,
    picard_bedgraph: Path | None = None,
) -> dict:
    source_raw = summarize_bam(input_bam)
    sambamba_raw = summarize_bam(sambamba_bam)
    picard_raw = summarize_bam(picard_bam)
    s_keys = sambamba_raw["duplicate_keys"]
    p_keys = picard_raw["duplicate_keys"]
    s_names = sambamba_raw["duplicate_names"]
    p_names = picard_raw["duplicate_names"]
    source = {key: value for key, value in source_raw.items() if not key.startswith("duplicate_")}
    sambamba = {key: value for key, value in sambamba_raw.items() if not key.startswith("duplicate_")}
    picard = {key: value for key, value in picard_raw.items() if not key.startswith("duplicate_")}
    union = s_keys | p_keys
    intersection = s_keys & p_keys
    total_primary = source["counts"]["primary"]
    hs = None
    if sambamba_hs_metrics and picard_hs_metrics:
        s_hs = _parse_picard_metrics(sambamba_hs_metrics)
        p_hs = _parse_picard_metrics(picard_hs_metrics)
        hs = {
            "sambamba": {key: s_hs.get(key) for key in ("MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_EXC_DUPE")},
            "picard": {key: p_hs.get(key) for key in ("MEAN_TARGET_COVERAGE", "PCT_TARGET_BASES_10X", "PCT_TARGET_BASES_20X", "PCT_TARGET_BASES_30X", "PCT_EXC_DUPE")},
        }
        hs["absolute_deltas"] = {
            key: abs((hs["sambamba"].get(key) or 0) - (hs["picard"].get(key) or 0))
            for key in hs["sambamba"]
        }
        hs["relative_deltas"] = {
            key: abs((hs["sambamba"].get(key) or 0) - (hs["picard"].get(key) or 0))
            / abs(hs["sambamba"].get(key))
            if isinstance(hs["sambamba"].get(key), (int, float)) and hs["sambamba"].get(key) != 0
            else None
            for key in hs["sambamba"]
        }
    return {
        "sample": sample,
        "group": group,
        "input": source,
        "sambamba": sambamba,
        "picard": picard,
        "structural_checks": {
            "record_counts_match_input": (
                source["counts"]["records"] == sambamba["counts"]["records"] == picard["counts"]["records"]
            ),
            "primary_counts_match_input": (
                source["counts"]["primary"] == sambamba["counts"]["primary"] == picard["counts"]["primary"]
            ),
            "reference_dictionaries_match": (
                source["header_reference_signature"]
                == sambamba["header_reference_signature"]
                == picard["header_reference_signature"]
            ),
            "outputs_coordinate_sorted": (
                sambamba["sort_order"] == picard["sort_order"] == "coordinate"
            ),
        },
        "duplicate_agreement": {
            "sambamba_duplicate_keys": len(s_keys),
            "picard_duplicate_keys": len(p_keys),
            "intersection_keys": len(intersection),
            "union_keys": len(union),
            "jaccard": len(intersection) / len(union) if union else 1.0,
            "symmetric_difference_keys": len(s_keys ^ p_keys),
            "primary_read_classification_agreement": 1 - len(s_keys ^ p_keys) / total_primary if total_primary else None,
            "duplicate_fraction_absolute_delta": abs(
                sambamba["counts"]["duplicate_primary_fraction"] - picard["counts"]["duplicate_primary_fraction"]
            ),
            "duplicate_name_jaccard": len(s_names & p_names) / len(s_names | p_names)
            if s_names | p_names else 1.0,
        },
        "sambamba_reported_metrics": _parse_sambamba_metrics(sambamba_metrics) if sambamba_metrics else None,
        "hs_metrics": hs,
        "methylation": _compare_bedgraphs(sambamba_bedgraph, picard_bedgraph),
    }


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_comparison_table(rows: list[dict], path: Path) -> None:
    fields = [
        "sample", "group", "structural_checks_pass", "input_primary_reads", "sambamba_duplicate_fraction",
        "picard_duplicate_fraction", "duplicate_fraction_absolute_delta",
        "primary_read_classification_agreement", "duplicate_key_jaccard",
        "shared_cpg_fraction", "mean_abs_methylation_percent_delta",
        "weighted_methylation_percent_absolute_delta",
        "mean_target_coverage_absolute_delta", "mean_target_coverage_relative_delta",
    ]
    output = []
    for row in rows:
        agreement = row["duplicate_agreement"]
        hs = row.get("hs_metrics") or {}
        hs_delta = (hs.get("absolute_deltas") or {}).get("MEAN_TARGET_COVERAGE")
        hs_relative_delta = (hs.get("relative_deltas") or {}).get("MEAN_TARGET_COVERAGE")
        methylation = row.get("methylation") or {}
        output.append({
            "sample": row["sample"],
            "group": row["group"],
            "structural_checks_pass": all(row["structural_checks"].values()),
            "input_primary_reads": row["input"]["counts"]["primary"],
            "sambamba_duplicate_fraction": row["sambamba"]["counts"]["duplicate_primary_fraction"],
            "picard_duplicate_fraction": row["picard"]["counts"]["duplicate_primary_fraction"],
            "duplicate_fraction_absolute_delta": agreement["duplicate_fraction_absolute_delta"],
            "primary_read_classification_agreement": agreement["primary_read_classification_agreement"],
            "duplicate_key_jaccard": agreement["jaccard"],
            "shared_cpg_fraction": methylation.get("shared_fraction_of_union"),
            "mean_abs_methylation_percent_delta": methylation.get("mean_abs_methylation_percent_delta"),
            "weighted_methylation_percent_absolute_delta": methylation.get("weighted_methylation_percent_absolute_delta"),
            "mean_target_coverage_absolute_delta": hs_delta,
            "mean_target_coverage_relative_delta": hs_relative_delta,
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(output)


def _write_comparison_figure(rows: list[dict], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    group_counts = Counter(row["group"] for row in rows)
    seen_groups = Counter()
    labels = []
    for row in rows:
        group = row["group"] or "sample"
        seen_groups[group] += 1
        labels.append(group if group_counts[group] == 1 else f"{group} {seen_groups[group]}")
    x = list(range(len(rows)))
    width = 0.36
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    duplicate_axis = axes[0, 0]
    agreement_axis = axes[0, 1]
    coverage_axis = axes[1, 0]
    methylation_axis = axes[1, 1]
    duplicate_axis.bar(
        [value - width / 2 for value in x],
        [row["sambamba"]["counts"]["duplicate_primary_fraction"] * 100 for row in rows],
        width,
        label="Sambamba",
        color="#457b9d",
    )
    duplicate_axis.bar(
        [value + width / 2 for value in x],
        [row["picard"]["counts"]["duplicate_primary_fraction"] * 100 for row in rows],
        width,
        label="Picard",
        color="#e76f51",
    )
    duplicate_axis.set_ylabel("Primary reads flagged duplicate (%)")
    duplicate_axis.set_title("Duplicate fraction")
    duplicate_axis.legend(frameon=False)
    duplicate_axis.grid(axis="y", color="#dddddd", linewidth=0.6)

    agreement_percent = [
        row["duplicate_agreement"]["primary_read_classification_agreement"] * 100
        for row in rows
    ]
    agreement_bars = agreement_axis.bar(
        x,
        agreement_percent,
        color="#2a9d8f",
        width=0.62,
    )
    lower = 0 if min(agreement_percent) < 95 else max(95, min(agreement_percent) - 0.25)
    agreement_axis.set_ylim(lower, 100.05)
    agreement_axis.set_ylabel("Agreement (%)")
    agreement_axis.set_title("Primary-read duplicate classification agreement")
    agreement_axis.grid(axis="y", color="#dddddd", linewidth=0.6)
    agreement_axis.bar_label(agreement_bars, fmt="%.3f%%", padding=3, fontsize=8)

    if all(row.get("hs_metrics") for row in rows):
        coverage_axis.bar(
            [value - width / 2 for value in x],
            [row["hs_metrics"]["sambamba"]["MEAN_TARGET_COVERAGE"] for row in rows],
            width,
            label="Sambamba",
            color="#457b9d",
        )
        coverage_axis.bar(
            [value + width / 2 for value in x],
            [row["hs_metrics"]["picard"]["MEAN_TARGET_COVERAGE"] for row in rows],
            width,
            label="Picard",
            color="#e76f51",
        )
    else:
        coverage_axis.text(0.5, 0.5, "HsMetrics not generated", ha="center", va="center", transform=coverage_axis.transAxes)
    coverage_axis.set_ylabel("Mean target coverage")
    coverage_axis.set_title("Twist covered-target depth")
    coverage_axis.grid(axis="y", color="#dddddd", linewidth=0.6)

    if all(row.get("methylation") for row in rows):
        methylation_axis.bar(
            [value - width / 2 for value in x],
            [row["methylation"]["sambamba_weighted_methylation_percent"] for row in rows],
            width,
            label="Sambamba",
            color="#457b9d",
        )
        methylation_axis.bar(
            [value + width / 2 for value in x],
            [row["methylation"]["picard_weighted_methylation_percent"] for row in rows],
            width,
            label="Picard",
            color="#e76f51",
        )
    else:
        methylation_axis.text(0.5, 0.5, "Methylation output not generated", ha="center", va="center", transform=methylation_axis.transAxes)
    methylation_axis.set_ylabel("Weighted CpG methylation (%)")
    methylation_axis.set_title("Shared CpG methylation")
    methylation_axis.grid(axis="y", color="#dddddd", linewidth=0.6)

    for axis in axes.flat:
        axis.set_xticks(x, labels, rotation=35, ha="right", fontsize=8)
    fig.suptitle("Sambamba versus Picard duplicate-marking comparison")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path, help="TSV with sample/group/input_bam/sambamba_bam/picard_bam")
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--metrics-dir", type=Path)
    parser.add_argument("--methylation-dir", type=Path)
    args = parser.parse_args()
    rows = []
    for row in _read_manifest(args.manifest):
        sample = row["sample"]
        metric_dir = args.metrics_dir / sample if args.metrics_dir else None
        meth_dir = args.methylation_dir / sample if args.methylation_dir else None
        rows.append(compare_sample(
            sample=sample,
            group=row.get("group", ""),
            input_bam=Path(row["input_bam"]),
            sambamba_bam=Path(row["sambamba_bam"]),
            picard_bam=Path(row["picard_bam"]),
            sambamba_metrics=Path(row["sambamba_metrics"]) if row.get("sambamba_metrics") else None,
            picard_hs_metrics=metric_dir / "picard.hs_metrics.txt" if metric_dir else None,
            sambamba_hs_metrics=metric_dir / "sambamba.hs_metrics.txt" if metric_dir else None,
            sambamba_bedgraph=meth_dir / "sambamba_CpG.bedGraph" if meth_dir else None,
            picard_bedgraph=meth_dir / "picard_CpG.bedGraph" if meth_dir else None,
        ))
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "generated_at": _utc_now(),
        "python": sys.version,
        "platform": platform.platform(),
        "manifest": str(args.manifest),
        "samples": rows,
        "acceptance_thresholds": None,
        "interpretation": "Raw measurements only; no equivalence gate was applied.",
    }
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    _write_comparison_table(rows, args.output_json.with_name("duplicate_marking_comparison.tsv"))
    _write_comparison_figure(rows, args.output_json.with_name("duplicate_marking_comparison.png"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
