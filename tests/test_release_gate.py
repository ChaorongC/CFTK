import csv
import hashlib
import json
from pathlib import Path

import pytest

from scripts.validation.check_release_gate import (
    evaluate_release_gate,
    write_release_evidence,
)


CLEAN_ID = "20260101T010101.000000Z-11111111"
RESUME_ID = "20260101T020202.000000Z-22222222"
OTHER_ID = "20260101T030303.000000Z-33333333"
IDENTITY = {
    "config_sha256": "a" * 64,
    "lock_sha256": "b" * 64,
    "options_sha256": "c" * 64,
}


def _software_identity():
    value = {
        "software_identity_schema_version": 1,
        "name": "cftk",
        "version": "1.0.0",
        "revision": "d" * 40,
        "source": "release",
        "dirty": False,
        "source_sha256": "e" * 64,
    }
    value["identity_sha256"] = hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return value


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_artifact_tables(run_dir: Path, artifacts: list[dict]) -> None:
    for filename, rows in (
        ("expected-outputs.tsv", [value for value in artifacts if value["role"] != "figure"]),
        ("figures.tsv", [value for value in artifacts if value["role"] == "figure"]),
    ):
        with (run_dir / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=("stage", "description", "path", "required"),
                delimiter="\t",
            )
            writer.writeheader()
            for value in rows:
                writer.writerow({
                    "stage": "process.1",
                    "description": value["description"],
                    "path": value["path"],
                    "required": value["required"],
                })


def _write_attempt(
    root: Path,
    run_id: str,
    stage_status: str,
    previous_run_id: str,
    artifacts: list[dict],
    *,
    schema: int,
) -> Path:
    run_dir = root / "provenance" / "runs" / run_id
    run_dir.mkdir(parents=True)
    resource_plan = {
        "resource_plan_version": 1,
        "sample_count": 1,
        "stages": [{"stage": "process.1", "applicable": True}],
    }
    manifest = {
        "run_schema_version": schema,
        "run_id": run_id,
        "run_dir": str(run_dir.resolve()),
        "project_root": str(root.resolve()),
        "project_identity": dict(IDENTITY),
        "previous_run_id": previous_run_id,
        "status": "complete",
        "started_at": "2026-01-01T00:00:00+00:00",
        "finished_at": "2026-01-01T00:01:00+00:00",
        "stages": [{
            "id": "process.1",
            "name": "private fixture stage",
            "applicable": True,
            "status": stage_status,
            "expected": artifacts,
            "outputs": [],
            "figures": [],
        }],
    }
    if schema >= 2:
        manifest["resource_plan"] = resource_plan
        _write_json(run_dir / "resource-plan.json", resource_plan)
    if schema >= 3:
        evidence_dir = run_dir / "evidence"
        evidence_dir.mkdir()
        evidence_summary = {
            "run_id": run_id,
            "run_status": "complete",
            "required_artifacts": len(artifacts),
            "missing_required_artifacts": 0,
        }
        _write_json(evidence_dir / "workflow_validation_summary.json", evidence_summary)
        manifest["evidence"] = {
            "status": "complete",
            "directory": str(evidence_dir.resolve()),
            "summary": str((evidence_dir / "workflow_validation_summary.json").resolve()),
            "files": ["workflow_validation_summary.json"],
            "required_artifacts": len(artifacts),
            "missing_required_artifacts": 0,
        }
    if schema >= 4:
        manifest["software_identity"] = _software_identity()

    _write_json(run_dir / "run.json", manifest)
    doctor_checks = [{
        "id": "input.fixture",
        "status": "PASS",
        "summary": "private fixture is ready",
        "details": {"path": str(root.resolve())},
    }]
    if schema >= 4:
        doctor_checks.append({
            "id": "runtime.cftk",
            "status": "PASS",
            "summary": "CFTK 1.0.0 (dddddddddddd, release)",
            "details": {"software_identity": manifest["software_identity"]},
        })
    _write_json(run_dir / "doctor-before.json", {
        "status": "PASS",
        "exit_code": 0,
        "summary": {"pass": len(doctor_checks), "warn": 0, "fail": 0},
        "checks": doctor_checks,
    })
    _write_json(run_dir / "tool-versions.json", {
        "fixture_tool": {
            "status": "PASS",
            "summary": "fixture tool 1.0",
            "path": str((root / "private-bin" / "fixture-tool").resolve()),
        }
    })
    _write_artifact_tables(run_dir, artifacts)
    (run_dir / "run-summary.html").write_text(
        "<html><body>private fixture report</body></html>\n", encoding="utf-8"
    )

    event_name = "stage_completed" if stage_status == "complete" else "stage_resumed"
    events = [
        {"event": "run_started", "run_id": run_id, "timestamp": "start"},
    ]
    if stage_status == "complete":
        events.append({"event": "stage_started", "stage": "process.1", "timestamp": "start"})
    events.extend([
        {"event": event_name, "stage": "process.1", "timestamp": "finish"},
        {"event": "run_completed", "run_id": run_id, "timestamp": "finish"},
    ])
    (run_dir / "events.jsonl").write_text(
        "".join(json.dumps(value) + "\n" for value in events), encoding="utf-8"
    )

    if stage_status == "complete":
        private_name = "patient" + "_042"
        command_id = "d" * 32
        common = {
            "command_id": command_id,
            "command": f"fixture-tool --input {private_name}",
            "cwd": str(root.resolve()),
            "label": f"fixture [{private_name}]",
            "run_id": run_id,
            "shell": False,
        }
        commands = [
            {**common, "event": "start", "timestamp": "start"},
            {**common, "event": "finish", "timestamp": "finish", "returncode": 0},
        ]
        command_text = "".join(json.dumps(value) + "\n" for value in commands)
    else:
        command_text = ""
    (run_dir / "commands.jsonl").write_text(command_text, encoding="utf-8")
    return run_dir / "run.json"


def _make_pair(tmp_path: Path, *, schema: int = 1):
    private_name = "patient" + "_042"
    artifact_root = tmp_path / "private-results" / private_name
    artifact_root.mkdir(parents=True)
    report = artifact_root / "result.tsv"
    checkpoint = artifact_root / "checkpoint.done"
    figure = artifact_root / "result.png"
    report.write_text("metric\tvalue\nreads\t10\n", encoding="utf-8")
    checkpoint.touch()
    figure.write_bytes(b"not-a-real-png-but-nonempty")
    artifacts = [
        {
            "path": str(report.resolve()),
            "description": f"report for {private_name}",
            "role": "report",
            "required": True,
            "nonempty": True,
        },
        {
            "path": str(checkpoint.resolve()),
            "description": f"checkpoint for {private_name}",
            "role": "output",
            "required": True,
            "nonempty": False,
        },
        {
            "path": str(figure.resolve()),
            "description": f"figure for {private_name}",
            "role": "figure",
            "required": True,
            "nonempty": True,
        },
    ]
    clean = _write_attempt(
        tmp_path, CLEAN_ID, "complete", OTHER_ID, artifacts, schema=schema
    )
    resume = _write_attempt(
        tmp_path, RESUME_ID, "resumed", CLEAN_ID, artifacts, schema=schema
    )
    return clean, resume, artifacts


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _rewrite(path: Path, value: dict) -> None:
    _write_json(path, value)


def test_release_gate_accepts_clean_and_zero_command_resume_and_sanitizes_output(tmp_path):
    clean, resume, _ = _make_pair(tmp_path)

    report = evaluate_release_gate(clean, resume)
    first_json, first_text = write_release_evidence(report, tmp_path / "gate-one")
    second_json, second_text = write_release_evidence(report, tmp_path / "gate-two")

    assert report["status"] == "PASS"
    assert report["clean_run"]["commands"]["unique_commands"] == 1
    assert report["resume_run"]["commands"]["unique_commands"] == 0
    assert report["comparison"]["same_stage_and_artifact_contract"] is True
    assert first_json.read_bytes() == second_json.read_bytes()
    assert first_text.read_bytes() == second_text.read_bytes()
    exported = first_json.read_text(encoding="utf-8") + first_text.read_text(encoding="utf-8")
    private_name = "patient" + "_042"
    assert str(tmp_path) not in exported
    assert private_name not in exported
    assert "fixture-tool --input" not in exported
    assert ".png" not in exported


def test_release_gate_accepts_current_schema_integrated_evidence(tmp_path):
    clean, resume, _ = _make_pair(tmp_path, schema=3)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "PASS"
    assert report["clean_run"]["records"]["resource_plan_sidecar"] is True
    assert report["clean_run"]["integrated_evidence"]["status"] == "complete"
    assert report["resume_run"]["integrated_evidence"]["status"] == "complete"


def test_release_gate_accepts_schema_v4_with_clean_revision_bound_identity(tmp_path):
    clean, resume, _ = _make_pair(tmp_path, schema=4)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "PASS"
    assert report["clean_run"]["software_identity"]["revision"] == "d" * 40
    assert report["comparison"]["same_software_identity"] is True


@pytest.mark.parametrize("failure", ["dirty", "doctor_mismatch", "resume_mismatch"])
def test_release_gate_rejects_schema_v4_software_identity_failures(tmp_path, failure):
    clean, resume, _ = _make_pair(tmp_path, schema=4)
    if failure == "dirty":
        manifest = _load(clean)
        manifest["software_identity"]["dirty"] = True
        _rewrite(clean, manifest)
    elif failure == "doctor_mismatch":
        doctor = _load(clean.parent / "doctor-before.json")
        doctor["checks"][1]["details"]["software_identity"]["identity_sha256"] = "f" * 64
        _rewrite(clean.parent / "doctor-before.json", doctor)
    else:
        manifest = _load(resume)
        manifest["software_identity"] = _software_identity()
        manifest["software_identity"]["revision"] = "a" * 40
        manifest["software_identity"]["identity_sha256"] = hashlib.sha256(
            json.dumps(
                {key: value for key, value in manifest["software_identity"].items()
                 if key != "identity_sha256"},
                sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()
        _rewrite(resume, manifest)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"


@pytest.mark.parametrize("failure", ["failed", "unsupported", "malformed"])
def test_release_gate_rejects_invalid_clean_manifest(tmp_path, failure):
    clean, resume, _ = _make_pair(tmp_path)
    if failure == "malformed":
        clean.write_text("{not json", encoding="utf-8")
    else:
        manifest = _load(clean)
        if failure == "failed":
            manifest["status"] = "failed"
        else:
            manifest["run_schema_version"] = 99
        _rewrite(clean, manifest)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    assert next(value for value in report["checks"] if value["id"] == "clean.manifest")["status"] == "FAIL"


@pytest.mark.parametrize("failure", ["identity", "predecessor", "clean_status", "resume_status"])
def test_release_gate_rejects_pair_or_stage_mismatch(tmp_path, failure):
    clean, resume, _ = _make_pair(tmp_path)
    target = resume
    manifest = _load(target)
    if failure == "identity":
        manifest["project_identity"]["config_sha256"] = "e" * 64
    elif failure == "predecessor":
        manifest["previous_run_id"] = OTHER_ID
    elif failure == "clean_status":
        target = clean
        manifest = _load(target)
        manifest["stages"][0]["status"] = "adopted"
    else:
        manifest["stages"][0]["status"] = "complete"
    _rewrite(target, manifest)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"


@pytest.mark.parametrize("failure", ["nonzero", "unmatched", "malformed", "resume_command"])
def test_release_gate_rejects_invalid_command_provenance(tmp_path, failure):
    clean, resume, _ = _make_pair(tmp_path)
    clean_ledger = clean.parent / "commands.jsonl"
    resume_ledger = resume.parent / "commands.jsonl"
    if failure == "malformed":
        clean_ledger.write_text("{not json\n", encoding="utf-8")
    else:
        records = [json.loads(line) for line in clean_ledger.read_text(encoding="utf-8").splitlines()]
        if failure == "nonzero":
            records[1]["returncode"] = 2
        elif failure == "unmatched":
            records = records[:1]
        else:
            resume_ledger.write_text(
                "".join(json.dumps(value) + "\n" for value in records), encoding="utf-8"
            )
        if failure != "resume_command":
            clean_ledger.write_text(
                "".join(json.dumps(value) + "\n" for value in records), encoding="utf-8"
            )

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"


@pytest.mark.parametrize("failure", ["missing", "empty"])
def test_release_gate_rejects_invalid_required_artifact(tmp_path, failure):
    clean, resume, artifacts = _make_pair(tmp_path)
    artifact = Path(artifacts[0]["path"])
    if failure == "missing":
        artifact.unlink()
    else:
        artifact.write_text("", encoding="utf-8")

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    assert next(value for value in report["checks"] if value["id"] == "clean.artifacts")["status"] == "FAIL"


def test_release_gate_rejects_missing_run_record(tmp_path):
    clean, resume, _ = _make_pair(tmp_path)
    (clean.parent / "tool-versions.json").unlink()

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    failed = {value["id"] for value in report["checks"] if value["status"] == "FAIL"}
    assert "clean.records" in failed
    assert "clean.tools" in failed


def test_release_gate_rejects_artifact_table_mismatch(tmp_path):
    clean, resume, _ = _make_pair(tmp_path)
    with (resume.parent / "figures.tsv").open("a", encoding="utf-8") as handle:
        handle.write("process.1\textra\t/private/extra\tTrue\n")

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    check = next(value for value in report["checks"] if value["id"] == "resume.artifact_tables")
    assert check["status"] == "FAIL"


def test_release_gate_rejects_changed_resume_artifact_contract(tmp_path):
    clean, resume, _ = _make_pair(tmp_path)
    manifest = _load(resume)
    manifest["stages"][0]["expected"][0]["description"] = "changed private contract"
    _rewrite(resume, manifest)
    _write_artifact_tables(resume.parent, manifest["stages"][0]["expected"])

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    assert report["resume_run"]["artifacts"]["missing_required"] == 0
    assert report["resume_run"]["artifact_tables"]["output_rows"] == 2
    assert report["comparison"]["same_stage_and_artifact_contract"] is False


def test_release_gate_rejects_failed_doctor_and_tool_probe(tmp_path):
    clean, resume, _ = _make_pair(tmp_path)
    doctor_path = clean.parent / "doctor-before.json"
    doctor = _load(doctor_path)
    doctor.update({"status": "FAIL", "exit_code": 1, "summary": {"pass": 0, "warn": 0, "fail": 1}})
    doctor["checks"][0]["status"] = "FAIL"
    _rewrite(doctor_path, doctor)
    tools_path = clean.parent / "tool-versions.json"
    tools = _load(tools_path)
    tools["fixture_tool"]["status"] = "FAIL"
    _rewrite(tools_path, tools)

    report = evaluate_release_gate(clean, resume)

    assert report["status"] == "FAIL"
    failed = {value["id"] for value in report["checks"] if value["status"] == "FAIL"}
    assert {"clean.doctor", "clean.tools"} <= failed
