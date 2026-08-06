import json
from pathlib import Path

import pysam
import pytest

from scripts.validation.compare_duplicate_marking import (
    _write_comparison_figure,
    _write_comparison_table,
    compare_sample,
)
from scripts.validation.summarize_doctor_audit import summarize
from scripts.validation.summarize_workflow_run import summarize as summarize_workflow
from validation_reports import _qc_display_labels


def _doctor_check(sample, name, status, summary=None):
    return {
        "id": f"input.bam.{sample}.{name}",
        "status": status,
        "summary": summary or f"{name}: {status}",
    }


def test_doctor_audit_writes_tables_summary_and_figure(tmp_path):
    sample_sheet = tmp_path / "samples.tsv"
    sample_sheet.write_text(
        "sample\tgroup\trole\ncontrol_1\tControl\tcontrol\ncase_1\tsALS\tcase\n"
    )
    checks = []
    for sample in ("control_1", "case_1"):
        checks.extend(
            _doctor_check(sample, name, "PASS")
            for name in ("index", "dictionary", "read_group", "duplicates", "sorting")
        )
    checks[1] = _doctor_check("control_1", "dictionary", "FAIL")
    checks[7] = _doctor_check("case_1", "read_group", "WARN")
    doctor_json = tmp_path / "doctor.json"
    doctor_json.write_text(json.dumps({"status": "FAIL", "exit_code": 1, "checks": checks}))

    output_dir = tmp_path / "report"
    result = summarize(doctor_json, sample_sheet, output_dir)

    assert result["sample_count"] == 2
    assert result["overall_status_counts"] == {"FAIL": 1, "WARN": 1}
    assert (output_dir / "cohort_readiness.tsv").is_file()
    assert (output_dir / "cohort_readiness_checks.tsv").is_file()
    assert (output_dir / "cohort_readiness.png").stat().st_size > 0


def _write_bam(path: Path, duplicate_names=()):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, name in enumerate(("read_a", "read_b", "read_c", "read_d")):
            for mate, flag, start, mate_start in (
                (1, 99, index * 100 + 10, index * 100 + 60),
                (2, 147, index * 100 + 60, index * 100 + 10),
            ):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = "A" * 50
                read.flag = flag | (1024 if name in duplicate_names else 0)
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigar = ((0, 50),)
                read.next_reference_id = 0
                read.next_reference_start = mate_start
                read.template_length = 100 if mate == 1 else -100
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                output.write(read)


def _write_hs_metrics(path: Path, mean_coverage: float):
    path.write_text(
        "## METRICS CLASS\tpicard.analysis.directed.HsMetrics\n"
        "BAIT_SET\tMEAN_TARGET_COVERAGE\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_30X\tPCT_EXC_DUPE\n"
        f"targets\t{mean_coverage}\t0.9\t0.8\t0.7\t0.1\n"
    )


def test_duplicate_comparison_uses_primary_read_keys_and_downstream_metrics(tmp_path):
    input_bam = tmp_path / "input.bam"
    sambamba_bam = tmp_path / "sambamba.bam"
    picard_bam = tmp_path / "picard.bam"
    _write_bam(input_bam)
    _write_bam(sambamba_bam, {"read_a"})
    _write_bam(picard_bam, {"read_a", "read_b"})

    sambamba_metrics = tmp_path / "sambamba.txt"
    sambamba_metrics.write_text("found 2 duplicates\n")
    sambamba_hs = tmp_path / "sambamba.hs_metrics.txt"
    picard_hs = tmp_path / "picard.hs_metrics.txt"
    _write_hs_metrics(sambamba_hs, 10.0)
    _write_hs_metrics(picard_hs, 9.9)
    sambamba_bed = tmp_path / "sambamba_CpG.bedGraph"
    picard_bed = tmp_path / "picard_CpG.bedGraph"
    sambamba_bed.write_text("track type=bedGraph\nchr1\t10\t12\t50\t5\t5\nchr1\t20\t22\t60\t6\t4\n")
    picard_bed.write_text("track type=bedGraph\nchr1\t10\t12\t50\t5\t5\nchr1\t20\t22\t50\t5\t5\n")

    result = compare_sample(
        sample="sample_1",
        group="Control",
        input_bam=input_bam,
        sambamba_bam=sambamba_bam,
        picard_bam=picard_bam,
        sambamba_metrics=sambamba_metrics,
        sambamba_hs_metrics=sambamba_hs,
        picard_hs_metrics=picard_hs,
        sambamba_bedgraph=sambamba_bed,
        picard_bedgraph=picard_bed,
    )

    agreement = result["duplicate_agreement"]
    assert all(result["structural_checks"].values())
    assert agreement["jaccard"] == 0.5
    assert agreement["primary_read_classification_agreement"] == 0.75
    assert result["hs_metrics"]["absolute_deltas"]["MEAN_TARGET_COVERAGE"] == pytest.approx(0.1)
    assert result["methylation"]["shared_fraction_of_union"] == 1.0
    assert result["methylation"]["mean_abs_methylation_percent_delta"] == 5.0
    assert result["methylation"]["weighted_methylation_percent_absolute_delta"] == pytest.approx(5.0)

    _write_comparison_table([result], tmp_path / "comparison.tsv")
    _write_comparison_figure([result], tmp_path / "comparison.png")
    assert (tmp_path / "comparison.tsv").is_file()
    assert (tmp_path / "comparison.png").stat().st_size > 0


def test_workflow_summary_records_artifacts_resources_commands_and_qc_figure(tmp_path):
    output = tmp_path / "results" / "qc_summary.tsv"
    output.parent.mkdir(parents=True)
    output.write_text(
        "sample\tgroup\tflagstat_mapped_pct\tmarkdup_dup_pct\n"
        "private_control\tControl\t95\t1.2\n"
    )
    manifest_dir = tmp_path / "run"
    manifest_dir.mkdir()
    (manifest_dir / "commands.jsonl").write_text(
        json.dumps({
            "event": "start", "command_id": "c1", "label": "trim [private_control]",
            "command": "trim_galore --cores 4", "timestamp": "now",
        }) + "\n" + json.dumps({
            "event": "finish", "command_id": "c1", "label": "trim [private_control]",
            "command": "trim_galore --cores 4", "returncode": 0, "timestamp": "later",
        }) + "\n"
    )
    manifest = {
        "run_id": "run-1",
        "status": "complete",
        "resource_plan": {
            "stages": [{
                "stage": "qc.0", "applicable": True, "total_core_budget": 8,
                "concurrent_samples": 1, "threads_per_sample": 1,
                "estimated_peak_threads": 1,
            }]
        },
        "stages": [{
            "id": "qc.0", "name": "QC", "status": "complete",
            "command": "cftk qc -s 0", "expected": [
                {"path": str(output), "role": "report", "required": True,
                 "description": "QC table"},
                {"path": str(tmp_path / "checkpoint.done"), "role": "output",
                 "required": True, "nonempty": False, "description": "Checkpoint"},
                {"path": str(tmp_path / "missing.png"), "role": "figure",
                 "required": False, "description": "Optional plot"},
            ]
        }],
    }
    (tmp_path / "checkpoint.done").touch()
    manifest_path = manifest_dir / "run.json"
    manifest_path.write_text(json.dumps(manifest))

    summary = summarize_workflow(manifest_path, tmp_path / "evidence")

    assert summary["required_artifacts"] == 2
    assert summary["missing_required_artifacts"] == 0
    assert summary["command_ledger"]["starts"] == 1
    assert summary["command_ledger"]["finishes"] == 1
    assert summary["command_ledger"]["unfinished"] == 0
    assert (tmp_path / "evidence/workflow_command_evidence.tsv").is_file()
    assert (tmp_path / "evidence/workflow_stage_evidence.png").stat().st_size > 0
    assert (tmp_path / "evidence/workflow_resource_plan.png").stat().st_size > 0
    assert (tmp_path / "evidence/workflow_qc_overview.png").stat().st_size > 0
    assert "workflow_validation_summary.json" in summary["files"]


def test_workflow_summary_records_fragmentomics_scope(tmp_path):
    manifest_dir = tmp_path / "run"
    manifest_dir.mkdir()
    scope_path = tmp_path / "results/4_fragmentomics/delfi/fragmentomics_scope.json"
    scope_path.parent.mkdir(parents=True)
    scope_path.write_text(json.dumps({
        "resolved_scope": {
            "mode": "panel",
            "target_sha256": "target-hash",
            "bins_count": 7,
            "note": "panel-only DELFI",
        }
    }))
    manifest = {
        "run_id": "run-scope",
        "status": "complete",
        "fragmentomics_scope": {"mode": "panel"},
        "stages": [{
            "id": "fragmentomics.delfi",
            "status": "complete",
            "expected": [{
                "path": str(scope_path),
                "role": "metadata",
                "required": True,
                "description": "scope",
            }],
        }],
    }
    manifest_path = manifest_dir / "run.json"
    manifest_path.write_text(json.dumps(manifest))

    summary = summarize_workflow(manifest_path, tmp_path / "evidence")

    assert summary["fragmentomics_scope"]["target_sha256"] == "target-hash"
    assert summary["fragmentomics_scope"]["bins_count"] == 7


def test_qc_figure_labels_omit_sample_identifiers():
    import pandas as pd

    frame = pd.DataFrame({
        "sample": ["private_control_001", "private_case_002"],
        "group": ["Control", "sALS"],
    })

    labels, _ = _qc_display_labels(frame)

    assert labels == ["Control 1", "sALS 1"]
    assert all("private_" not in label for label in labels)
