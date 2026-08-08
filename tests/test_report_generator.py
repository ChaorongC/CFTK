import json
import os
from pathlib import Path


def _write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload))


def _stage(stage_id, name, status, expected, **extra):
    return {
        "id": stage_id,
        "name": name,
        "status": status,
        "expected": expected,
        **extra,
    }


def test_workflow_summary_aggregates_trusted_manifests(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from report import report_generator

    results = tmp_path / "results"
    provenance = results / "provenance"
    core_output = results / "1_process" / "alignment.bam"
    wps_figure = results / "4_fragmentomics" / "wps" / "sample.wps_profile.png"
    delfi_output = results / "4_fragmentomics" / "delfi" / "sample_delfi.tsv"
    for artifact in (core_output, wps_figure, delfi_output):
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text("artifact\n")

    core_manifest = provenance / "runs" / "core" / "run.json"
    _write_json(core_manifest, {
        "run_id": "core-run",
        "stages": [_stage(
            "process.alignment", "Alignment", "complete",
            [{"path": str(core_output), "role": "output", "required": True}],
        )],
    })
    _write_json(provenance / "latest-run.json", {"manifest": str(core_manifest)})

    prior = provenance / "analysis-runs" / "prior" / "run.json"
    _write_json(prior, {
        "run_id": "prior-run",
        "stages": [_stage(
            "fragmentomics.wps", "Windowed protection score", "resumed",
            [{"path": str(wps_figure), "role": "figure", "required": True}],
            fragmentomics_scope={
                "note": "Targeted-panel mode is not genome-wide.",
            },
        )],
    })
    newer = provenance / "analysis-runs" / "newer" / "run.json"
    _write_json(newer, {
        "run_id": "newer-run",
        "stages": [
            _stage("fragmentomics.wps", "Windowed protection score", "failed", []),
            _stage(
                "fragmentomics.delfi", "DELFI-style fragmentomics", "adopted",
                [{"path": str(delfi_output), "role": "output", "required": True}],
            ),
        ],
    })
    os.utime(prior, (1, 1))
    os.utime(newer, (2, 2))

    section = report_generator._sec_workflow_summary(
        str(results), {"Control": ["sample"]}
    )

    assert "Workflow Summary" in section
    assert "Core processing / QC" in section
    assert "Alignment" in section
    assert "Windowed protection score" in section
    assert "DELFI-style fragmentomics" in section
    assert "resumed" in section
    assert "adopted" in section
    assert "Targeted-panel mode is not genome-wide." in section
    assert "prior-run" in section
    assert "newer-run" in section
    assert ">0</td>" in section


def test_fragmentomics_section_includes_occupancy(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from report import report_generator

    section = report_generator._sec_fragmentomics(
        str(tmp_path), {"Control": ["sample"]}
    )

    assert 'id="part4_1"' in section
    assert "4.1 Nucleosome Occupancy" in section
    assert "4.5 WPS" in section


def test_source_results_and_qc_scores_are_discovered_from_markdup_bams(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from report import report_generator

    local_results = tmp_path / "downstream" / "results"
    source_results = tmp_path / "processed" / "results"
    bam = source_results / "1_process" / "3_markdup" / "sample_b.markdup.bam"
    bam.parent.mkdir(parents=True)
    bam.write_text("placeholder\n")
    local_qc = local_results / "2_qc" / "qc_scores.tsv"
    source_qc = source_results / "2_qc" / "qc_scores.tsv"
    local_qc.parent.mkdir(parents=True)
    source_qc.parent.mkdir(parents=True)
    local_qc.write_text("sample\tqc_status\nsample_a\tPASS\n")
    source_qc.write_text("sample\tqc_status\nsample_b\tWARN\n")

    cfg = {
        "samples": {
            "Control": [{"name": "sample_a", "input_type": "bam", "bam": "outside.bam"}],
            "Case": [{"name": "sample_b", "input_type": "bam", "bam": str(bam)}],
        },
    }
    discovered = report_generator._source_results_from_config(cfg, str(local_results))

    assert discovered == {"sample_b": str(source_results)}
    scores = report_generator._read_qc_scores(
        str(local_results), {"Control": ["sample_a"], "Case": ["sample_b"]}, discovered
    )
    assert scores["sample"].tolist() == ["sample_a", "sample_b"]
    assert scores["qc_status"].tolist() == ["PASS", "WARN"]
