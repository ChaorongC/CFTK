#!/usr/bin/env python3
"""Validate a clean CFTK core run and its zero-command resume attempt.

This maintainer-only gate reads completed run records. It never launches CFTK,
external tools, or scheduler jobs, and its outputs omit paths, commands, sample
identifiers, filenames derived from samples, and patient-level values.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import uuid


GATE_SCHEMA_VERSION = 1
SUPPORTED_RUN_SCHEMA_VERSIONS = (1, 2, 3)
REQUIRED_IDENTITY_KEYS = (
    "config_sha256",
    "lock_sha256",
    "options_sha256",
)
BASE_RUN_FILES = (
    "run.json",
    "commands.jsonl",
    "doctor-before.json",
    "events.jsonl",
    "expected-outputs.tsv",
    "figures.tsv",
    "run-summary.html",
    "tool-versions.json",
)
CORE_STAGE_PATTERN = re.compile(r"^(?:process\.[1-4]|qc\.[0-3])$")
RUN_ID_PATTERN = re.compile(r"^\d{8}T\d{6}\.\d{6}Z-[0-9a-f]{8}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
COMMAND_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")
ARTIFACT_ROLES = {"figure", "output", "report"}
KNOWN_EVENTS = {
    "evidence_completed",
    "evidence_failed",
    "run_completed",
    "run_failed",
    "run_planned",
    "run_started",
    "stage_adopted",
    "stage_completed",
    "stage_failed",
    "stage_outputs_quarantined",
    "stage_resumed",
    "stage_skipped",
    "stage_started",
}


class GateValidationError(RuntimeError):
    """A release record does not satisfy a public-safe acceptance contract."""


class _Checks:
    def __init__(self):
        self.records = []

    def add(self, check_id, passed, success, failure, **details):
        record = {
            "id": check_id,
            "status": "PASS" if passed else "FAIL",
            "summary": success if passed else failure,
        }
        if details:
            record["details"] = details
        self.records.append(record)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _contract_sha256(artifacts: list[dict]) -> str:
    payload = json.dumps(
        sorted(
            artifacts,
            key=lambda value: (
                value["stage"], value["path"], value["role"],
                value["description"],
            ),
        ),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _read_json_object(path: Path, record_name: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise GateValidationError(f"{record_name} is missing") from exc
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise GateValidationError(f"{record_name} is unreadable or malformed") from exc
    if not isinstance(value, dict):
        raise GateValidationError(f"{record_name} must contain a JSON object")
    return value


def _read_jsonl(path: Path, record_name: str, *, allow_empty: bool = False) -> list[dict]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except FileNotFoundError as exc:
        raise GateValidationError(f"{record_name} is missing") from exc
    except (OSError, UnicodeError) as exc:
        raise GateValidationError(f"{record_name} is unreadable") from exc
    records = []
    for line in lines:
        if not line.strip():
            continue
        try:
            value = json.loads(line)
        except json.JSONDecodeError as exc:
            raise GateValidationError(f"{record_name} contains malformed JSON") from exc
        if not isinstance(value, dict):
            raise GateValidationError(f"{record_name} contains a non-object record")
        records.append(value)
    if not records and not allow_empty:
        raise GateValidationError(f"{record_name} contains no records")
    return records


def _load_manifest(manifest_path: Path) -> tuple[dict, Path]:
    manifest_path = Path(manifest_path).expanduser().resolve()
    manifest = _read_json_object(manifest_path, "run manifest")
    run_dir = manifest_path.parent

    schema = manifest.get("run_schema_version")
    if schema not in SUPPORTED_RUN_SCHEMA_VERSIONS:
        raise GateValidationError("run manifest uses an unsupported core schema")
    run_id = manifest.get("run_id")
    if not isinstance(run_id, str) or not RUN_ID_PATTERN.fullmatch(run_id):
        raise GateValidationError("run manifest has an invalid run identifier")
    if run_dir.name != run_id:
        raise GateValidationError("run manifest is not stored in its run directory")
    raw_run_dir = manifest.get("run_dir")
    if not isinstance(raw_run_dir, str) or Path(raw_run_dir).expanduser().resolve() != run_dir:
        raise GateValidationError("run manifest directory identity is inconsistent")
    if manifest.get("status") != "complete":
        raise GateValidationError("run manifest is not complete")
    if not isinstance(manifest.get("started_at"), str) or not manifest["started_at"]:
        raise GateValidationError("run manifest has no start timestamp")
    if not isinstance(manifest.get("finished_at"), str) or not manifest["finished_at"]:
        raise GateValidationError("run manifest has no finish timestamp")

    identity = manifest.get("project_identity")
    if not isinstance(identity, dict) or set(identity) != set(REQUIRED_IDENTITY_KEYS):
        raise GateValidationError("run manifest has an invalid project identity")
    if any(
        not isinstance(identity[key], str) or not SHA256_PATTERN.fullmatch(identity[key])
        for key in REQUIRED_IDENTITY_KEYS
    ):
        raise GateValidationError("run manifest has a non-SHA-256 project identity")
    if not isinstance(manifest.get("stages"), list) or not manifest["stages"]:
        raise GateValidationError("run manifest has no stage records")
    return manifest, run_dir


def _validate_run_files(manifest: dict, run_dir: Path) -> dict:
    required = list(BASE_RUN_FILES)
    if manifest["run_schema_version"] >= 2:
        required.append("resource-plan.json")
    missing = [name for name in required if not (run_dir / name).is_file()]
    nonempty = [name for name in required if name != "commands.jsonl"]
    empty = [
        name for name in nonempty
        if (run_dir / name).is_file() and (run_dir / name).stat().st_size == 0
    ]
    if missing or empty:
        raise GateValidationError(
            "required run-level records are missing or empty"
        )
    if manifest["run_schema_version"] >= 2:
        resource_plan = _read_json_object(run_dir / "resource-plan.json", "resource plan")
        if resource_plan != manifest.get("resource_plan"):
            raise GateValidationError("resource plan sidecar does not match the manifest")
    return {
        "required_files": len(required),
        "present_files": len(required) - len(missing),
        "resource_plan_sidecar": manifest["run_schema_version"] >= 2,
    }


def _normalise_artifact(stage_id: str, raw: dict) -> dict:
    if not isinstance(raw, dict):
        raise GateValidationError("stage artifact contract contains a non-object entry")
    path = raw.get("path")
    description = raw.get("description")
    role = raw.get("role", "output")
    required = raw.get("required", True)
    nonempty = raw.get("nonempty", True)
    if not isinstance(path, str) or not path or not Path(path).is_absolute():
        raise GateValidationError("stage artifact contract has an invalid path")
    if not isinstance(description, str) or not description:
        raise GateValidationError("stage artifact contract has an invalid description")
    if role not in ARTIFACT_ROLES:
        raise GateValidationError("stage artifact contract has an invalid role")
    if not isinstance(required, bool) or not isinstance(nonempty, bool):
        raise GateValidationError("stage artifact contract has invalid boolean fields")
    return {
        "stage": stage_id,
        "path": path,
        "description": description,
        "role": role,
        "required": required,
        "nonempty": nonempty,
    }


def _validate_stages(manifest: dict, expected_status: str) -> tuple[dict, dict]:
    stage_map = {}
    artifacts = []
    public_stages = []
    for stage in manifest["stages"]:
        if not isinstance(stage, dict):
            raise GateValidationError("run manifest contains a non-object stage")
        stage_id = stage.get("id")
        if not isinstance(stage_id, str) or not CORE_STAGE_PATTERN.fullmatch(stage_id):
            raise GateValidationError("run manifest contains a non-core stage identifier")
        if stage_id in stage_map:
            raise GateValidationError("run manifest contains duplicate stage identifiers")
        applicable = stage.get("applicable")
        if not isinstance(applicable, bool):
            raise GateValidationError("run manifest stage applicability is invalid")
        status = stage.get("status")
        required_status = expected_status if applicable else "skipped"
        if status != required_status:
            raise GateValidationError("run manifest stage status violates the gate contract")
        expected = stage.get("expected")
        if not isinstance(expected, list):
            raise GateValidationError("run manifest stage has no artifact contract")
        if applicable and not expected:
            raise GateValidationError("an applicable stage has an empty artifact contract")
        if not applicable and expected:
            raise GateValidationError("an inapplicable stage has unexpected artifacts")
        stage_artifacts = [_normalise_artifact(stage_id, value) for value in expected]
        artifact_paths = [value["path"] for value in stage_artifacts]
        if len(artifact_paths) != len(set(artifact_paths)):
            raise GateValidationError("a stage artifact contract contains duplicate paths")
        artifacts.extend(stage_artifacts)
        stage_map[stage_id] = {
            "applicable": applicable,
            "status": status,
            "artifacts": stage_artifacts,
        }
        role_counts = Counter(value["role"] for value in stage_artifacts)
        public_stages.append({
            "id": stage_id,
            "applicable": applicable,
            "status": status,
            "artifacts": len(stage_artifacts),
            "required_artifacts": sum(value["required"] for value in stage_artifacts),
            "nonempty_required_artifacts": sum(
                value["required"] and value["nonempty"] for value in stage_artifacts
            ),
            "role_counts": {
                role: role_counts.get(role, 0) for role in sorted(ARTIFACT_ROLES)
            },
        })
    if not any(value["applicable"] for value in stage_map.values()):
        raise GateValidationError("run manifest has no applicable core stages")
    all_paths = [value["path"] for value in artifacts]
    if len(all_paths) != len(set(all_paths)):
        raise GateValidationError("artifact paths are duplicated across core stages")
    public = {
        "applicable": sum(value["applicable"] for value in stage_map.values()),
        "skipped": sum(not value["applicable"] for value in stage_map.values()),
        "records": public_stages,
    }
    return {"map": stage_map, "artifacts": artifacts}, public


def _validate_artifacts(artifacts: list[dict]) -> dict:
    required = [value for value in artifacts if value["required"]]
    missing = 0
    empty = 0
    for artifact in required:
        path = Path(artifact["path"])
        try:
            if not path.exists():
                missing += 1
            elif artifact["nonempty"] and path.is_file() and path.stat().st_size == 0:
                empty += 1
        except OSError:
            missing += 1
    if missing or empty:
        raise GateValidationError("required artifact contracts are missing or empty")
    role_counts = Counter(value["role"] for value in artifacts)
    return {
        "expected": len(artifacts),
        "required": len(required),
        "valid_required": len(required) - missing - empty,
        "missing_required": missing,
        "empty_required": empty,
        "contract_sha256": _contract_sha256(artifacts),
        "role_counts": {
            role: role_counts.get(role, 0) for role in sorted(ARTIFACT_ROLES)
        },
    }


def _parse_table_bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise GateValidationError("artifact table contains an invalid required value")


def _read_artifact_table(path: Path, record_name: str) -> list[tuple]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != ["stage", "description", "path", "required"]:
                raise GateValidationError(f"{record_name} has an invalid header")
            rows = []
            for row in reader:
                if None in row or any(row[key] is None for key in reader.fieldnames):
                    raise GateValidationError(f"{record_name} contains a malformed row")
                rows.append((
                    row["stage"], row["description"], row["path"],
                    _parse_table_bool(row["required"]),
                ))
    except FileNotFoundError as exc:
        raise GateValidationError(f"{record_name} is missing") from exc
    except (OSError, UnicodeError, csv.Error) as exc:
        raise GateValidationError(f"{record_name} is unreadable or malformed") from exc
    return rows


def _validate_artifact_tables(run_dir: Path, artifacts: list[dict]) -> dict:
    outputs = _read_artifact_table(run_dir / "expected-outputs.tsv", "output table")
    figures = _read_artifact_table(run_dir / "figures.tsv", "figure table")
    expected_outputs = [
        (value["stage"], value["description"], value["path"], value["required"])
        for value in artifacts if value["role"] != "figure"
    ]
    expected_figures = [
        (value["stage"], value["description"], value["path"], value["required"])
        for value in artifacts if value["role"] == "figure"
    ]
    if Counter(outputs) != Counter(expected_outputs) or Counter(figures) != Counter(expected_figures):
        raise GateValidationError("artifact tables do not match the manifest contract")
    return {"output_rows": len(outputs), "figure_rows": len(figures)}


def _validate_doctor(run_dir: Path) -> tuple[dict, dict]:
    doctor = _read_json_object(run_dir / "doctor-before.json", "doctor report")
    checks = doctor.get("checks")
    summary = doctor.get("summary")
    if not isinstance(checks, list) or not checks or not isinstance(summary, dict):
        raise GateValidationError("doctor report has an invalid result structure")
    statuses = []
    for check in checks:
        if not isinstance(check, dict) or check.get("status") not in {"PASS", "WARN", "FAIL"}:
            raise GateValidationError("doctor report contains an invalid check")
        statuses.append(check["status"])
    counts = {
        "pass": statuses.count("PASS"),
        "warn": statuses.count("WARN"),
        "fail": statuses.count("FAIL"),
    }
    if any(summary.get(key) != value for key, value in counts.items()):
        raise GateValidationError("doctor summary does not match its checks")
    expected_status = "FAIL" if counts["fail"] else ("WARN" if counts["warn"] else "PASS")
    if doctor.get("exit_code") != 0 or counts["fail"] != 0 or doctor.get("status") != expected_status:
        raise GateValidationError("doctor report contains a required failure")
    public = {"status": doctor["status"], "exit_code": 0, **counts}
    return doctor, public


def _validate_tools(run_dir: Path) -> tuple[dict, dict]:
    tools = _read_json_object(run_dir / "tool-versions.json", "tool-version report")
    if not tools:
        raise GateValidationError("tool-version report is empty")
    statuses = []
    for value in tools.values():
        if not isinstance(value, dict) or value.get("status") not in {"PASS", "WARN", "FAIL"}:
            raise GateValidationError("tool-version report contains an invalid entry")
        statuses.append(value["status"])
    counts = {
        "pass": statuses.count("PASS"),
        "warn": statuses.count("WARN"),
        "fail": statuses.count("FAIL"),
    }
    if counts["fail"]:
        raise GateValidationError("tool-version report contains a failed tool probe")
    return tools, {"tools": len(tools), **counts}


def _validate_commands(run_dir: Path, manifest: dict, *, require_zero: bool) -> dict:
    ledger = run_dir / "commands.jsonl"
    records = _read_jsonl(ledger, "command ledger", allow_empty=require_zero)
    if require_zero:
        if records:
            raise GateValidationError("resume attempt executed one or more recorded commands")
        return {
            "starts": 0,
            "finishes": 0,
            "unique_commands": 0,
            "nonzero_finishes": 0,
            "unmatched_commands": 0,
            "ledger_sha256": _sha256(ledger),
        }

    by_id = {}
    for index, record in enumerate(records):
        event = record.get("event")
        command_id = record.get("command_id")
        if event not in {"start", "finish"}:
            raise GateValidationError("command ledger contains an invalid event")
        if not isinstance(command_id, str) or not COMMAND_ID_PATTERN.fullmatch(command_id):
            raise GateValidationError("command ledger contains an invalid command identifier")
        if record.get("run_id") != manifest["run_id"]:
            raise GateValidationError("command ledger run identity is inconsistent")
        if not isinstance(record.get("command"), str) or not record["command"]:
            raise GateValidationError("command ledger contains an invalid command")
        if not isinstance(record.get("cwd"), str) or not Path(record["cwd"]).is_absolute():
            raise GateValidationError("command ledger contains an invalid working directory")
        if not isinstance(record.get("label"), str) or not isinstance(record.get("shell"), bool):
            raise GateValidationError("command ledger contains invalid execution metadata")
        if not isinstance(record.get("timestamp"), str) or not record["timestamp"]:
            raise GateValidationError("command ledger contains an invalid timestamp")
        slot = by_id.setdefault(command_id, {})
        if event in slot:
            raise GateValidationError("command ledger contains a duplicate event")
        slot[event] = (index, record)

    for pair in by_id.values():
        if set(pair) != {"start", "finish"}:
            raise GateValidationError("command ledger contains an unmatched command")
        start_index, start = pair["start"]
        finish_index, finish = pair["finish"]
        if finish_index <= start_index:
            raise GateValidationError("command ledger contains an out-of-order command pair")
        for key in ("command", "cwd", "label", "run_id", "shell"):
            if start.get(key) != finish.get(key):
                raise GateValidationError("command ledger pair metadata is inconsistent")
        returncode = finish.get("returncode")
        if isinstance(returncode, bool) or not isinstance(returncode, int) or returncode != 0:
            raise GateValidationError("command ledger contains a nonzero or invalid finish")

    starts = sum(record.get("event") == "start" for record in records)
    finishes = sum(record.get("event") == "finish" for record in records)
    if not by_id or starts != finishes or starts != len(by_id):
        raise GateValidationError("command ledger does not contain paired commands")
    return {
        "starts": starts,
        "finishes": finishes,
        "unique_commands": len(by_id),
        "nonzero_finishes": 0,
        "unmatched_commands": 0,
        "ledger_sha256": _sha256(ledger),
    }


def _validate_events(
    run_dir: Path,
    manifest: dict,
    stage_internal: dict,
    expected_stage_event: str,
) -> dict:
    records = _read_jsonl(run_dir / "events.jsonl", "event ledger")
    event_counts = Counter()
    stage_events = {}
    known_stages = set(stage_internal["map"])
    for record in records:
        event = record.get("event")
        if event not in KNOWN_EVENTS:
            raise GateValidationError("event ledger contains an invalid event")
        if not isinstance(record.get("timestamp"), str) or not record["timestamp"]:
            raise GateValidationError("event ledger contains an invalid timestamp")
        if "run_id" in record and record["run_id"] != manifest["run_id"]:
            raise GateValidationError("event ledger run identity is inconsistent")
        stage_id = record.get("stage")
        if stage_id is not None:
            if stage_id not in known_stages:
                raise GateValidationError("event ledger contains an unknown stage")
            stage_events.setdefault(event, []).append(stage_id)
        event_counts[event] += 1

    if event_counts["run_started"] != 1 or event_counts["run_completed"] != 1:
        raise GateValidationError("event ledger lacks a unique start/completion pair")
    if event_counts["run_failed"] or event_counts["stage_failed"]:
        raise GateValidationError("event ledger records a workflow failure")
    applicable = {
        stage_id for stage_id, value in stage_internal["map"].items()
        if value["applicable"]
    }
    skipped = known_stages - applicable
    if Counter(stage_events.get(expected_stage_event, [])) != Counter(applicable):
        raise GateValidationError("event ledger does not match applicable stage outcomes")
    if Counter(stage_events.get("stage_skipped", [])) != Counter(skipped):
        raise GateValidationError("event ledger does not match skipped stage outcomes")
    if expected_stage_event == "stage_completed":
        if Counter(stage_events.get("stage_started", [])) != Counter(applicable):
            raise GateValidationError("clean event ledger lacks stage execution records")
        if event_counts["stage_resumed"] or event_counts["stage_adopted"]:
            raise GateValidationError("clean event ledger records reused stages")
    else:
        if event_counts["stage_started"] or event_counts["stage_completed"] or event_counts["stage_adopted"]:
            raise GateValidationError("resume event ledger records stage execution")
    return {
        "records": len(records),
        "run_started": event_counts["run_started"],
        "run_completed": event_counts["run_completed"],
        "stage_outcomes": event_counts[expected_stage_event],
        "failures": event_counts["run_failed"] + event_counts["stage_failed"],
    }


def _validate_evidence(manifest: dict, run_dir: Path, artifact_summary: dict) -> dict:
    evidence = manifest.get("evidence")
    required = manifest["run_schema_version"] >= 3
    if evidence is None and not required:
        return {"required": False, "available": False, "status": "not_required"}
    if not isinstance(evidence, dict) or evidence.get("status") != "complete":
        raise GateValidationError("integrated evidence is absent or incomplete")
    expected_dir = run_dir / "evidence"
    directory = evidence.get("directory")
    summary_path = evidence.get("summary")
    if (
        not isinstance(directory, str)
        or Path(directory).expanduser().resolve() != expected_dir
        or not isinstance(summary_path, str)
        or Path(summary_path).expanduser().resolve() != expected_dir / "workflow_validation_summary.json"
    ):
        raise GateValidationError("integrated evidence location is inconsistent")
    files = evidence.get("files")
    if not isinstance(files, list) or "workflow_validation_summary.json" not in files:
        raise GateValidationError("integrated evidence file inventory is incomplete")
    if len(files) != len(set(files)) or any(
        not isinstance(name, str) or not name or Path(name).name != name for name in files
    ):
        raise GateValidationError("integrated evidence file inventory is invalid")
    if any(
        not (expected_dir / name).is_file() or (expected_dir / name).stat().st_size == 0
        for name in files
    ):
        raise GateValidationError("integrated evidence contains missing or empty files")
    summary = _read_json_object(
        expected_dir / "workflow_validation_summary.json", "integrated evidence summary"
    )
    if (
        summary.get("run_id") != manifest["run_id"]
        or summary.get("run_status") != "complete"
        or summary.get("missing_required_artifacts") != 0
        or summary.get("required_artifacts") != artifact_summary["required"]
        or evidence.get("missing_required_artifacts") != 0
        or evidence.get("required_artifacts") != artifact_summary["required"]
    ):
        raise GateValidationError("integrated evidence summary is inconsistent")
    return {
        "required": required,
        "available": True,
        "status": "complete",
        "files": len(files),
        "missing_required_artifacts": 0,
    }


def _run_component(checks: _Checks, check_id: str, success: str, function, failure: str):
    try:
        result = function()
    except GateValidationError:
        checks.add(check_id, False, success, failure)
        return None
    except (OSError, UnicodeError, ValueError, TypeError):
        checks.add(check_id, False, success, failure)
        return None
    checks.add(check_id, True, success, failure)
    return result


def _inspect_attempt(manifest_path: Path, role: str, expected_status: str, checks: _Checks):
    manifest_result = _run_component(
        checks,
        f"{role}.manifest",
        "Completed core run manifest is structurally valid.",
        lambda: _load_manifest(manifest_path),
        "Core run manifest is missing, malformed, unsupported, or incomplete.",
    )
    if manifest_result is None:
        return None, {"available": False}
    manifest, run_dir = manifest_result
    public = {
        "available": True,
        "run_id": manifest["run_id"],
        "run_schema_version": manifest["run_schema_version"],
        "status": manifest["status"],
        "project_identity": {
            key: manifest["project_identity"][key] for key in REQUIRED_IDENTITY_KEYS
        },
    }

    public["records"] = _run_component(
        checks,
        f"{role}.records",
        "Required run-level records are present and consistent.",
        lambda: _validate_run_files(manifest, run_dir),
        "Required run-level records are missing, empty, malformed, or inconsistent.",
    ) or {"validated": False}

    stage_result = _run_component(
        checks,
        f"{role}.stages",
        f"Applicable core stages are all {expected_status}.",
        lambda: _validate_stages(manifest, expected_status),
        f"Core stage identities, applicability, statuses, or contracts are not all {expected_status}.",
    )
    stage_internal = stage_result[0] if stage_result else None
    public["stages"] = stage_result[1] if stage_result else {"validated": False}

    doctor_result = _run_component(
        checks,
        f"{role}.doctor",
        "Doctor completed with exit code 0 and no failed checks.",
        lambda: _validate_doctor(run_dir),
        "Doctor report is missing, malformed, inconsistent, or contains a failed check.",
    )
    public["doctor"] = doctor_result[1] if doctor_result else {"validated": False}

    tools_result = _run_component(
        checks,
        f"{role}.tools",
        "Recorded tool probes contain no failures.",
        lambda: _validate_tools(run_dir),
        "Tool-version report is missing, malformed, empty, or contains a failed probe.",
    )
    public["tools"] = tools_result[1] if tools_result else {"validated": False}

    public["commands"] = _run_component(
        checks,
        f"{role}.commands",
        "Command ledger satisfies the execution contract.",
        lambda: _validate_commands(run_dir, manifest, require_zero=role == "resume"),
        "Command ledger is malformed, incomplete, nonzero, or violates zero-command resume.",
    ) or {"validated": False}

    artifact_summary = None
    if stage_internal is not None:
        artifact_summary = _run_component(
            checks,
            f"{role}.artifacts",
            "All required artifacts satisfy their existing file contracts.",
            lambda: _validate_artifacts(stage_internal["artifacts"]),
            "One or more required artifacts are missing, empty, or unreadable.",
        )
        public["artifacts"] = artifact_summary or {"validated": False}
        public["artifact_tables"] = _run_component(
            checks,
            f"{role}.artifact_tables",
            "Artifact tables exactly match the manifest contracts.",
            lambda: _validate_artifact_tables(run_dir, stage_internal["artifacts"]),
            "Artifact tables are missing, malformed, or inconsistent with the manifest.",
        ) or {"validated": False}
        public["events"] = _run_component(
            checks,
            f"{role}.events",
            "Event ledger matches the terminal stage outcomes.",
            lambda: _validate_events(
                run_dir,
                manifest,
                stage_internal,
                "stage_completed" if role == "clean" else "stage_resumed",
            ),
            "Event ledger is missing, malformed, failed, or inconsistent with stage outcomes.",
        ) or {"validated": False}
    else:
        for suffix, label in (
            ("artifacts", "Artifact validation requires a valid stage contract."),
            ("artifact_tables", "Artifact-table validation requires a valid stage contract."),
            ("events", "Event validation requires a valid stage contract."),
        ):
            checks.add(f"{role}.{suffix}", False, "", label)
            public[suffix] = {"validated": False}

    if artifact_summary is not None:
        public["integrated_evidence"] = _run_component(
            checks,
            f"{role}.integrated_evidence",
            "Integrated evidence satisfies its schema-version contract.",
            lambda: _validate_evidence(manifest, run_dir, artifact_summary),
            "Integrated evidence is missing, incomplete, or inconsistent for this schema.",
        ) or {"validated": False}
    else:
        checks.add(
            f"{role}.integrated_evidence",
            False,
            "",
            "Integrated-evidence validation requires valid artifact contracts.",
        )
        public["integrated_evidence"] = {"validated": False}

    internal = {
        "manifest": manifest,
        "run_dir": run_dir,
        "stages": stage_internal,
        "doctor": doctor_result[0] if doctor_result else None,
        "tools": tools_result[0] if tools_result else None,
    }
    return internal, public


def _compare_attempts(clean: dict, resume: dict) -> tuple[bool, dict]:
    manifest_clean = clean["manifest"]
    manifest_resume = resume["manifest"]
    same_identity = manifest_clean["project_identity"] == manifest_resume["project_identity"]
    immediate_resume = manifest_resume.get("previous_run_id") == manifest_clean["run_id"]
    same_schema = manifest_clean["run_schema_version"] == manifest_resume["run_schema_version"]
    clean_stages = clean["stages"]["map"] if clean["stages"] else {}
    resume_stages = resume["stages"]["map"] if resume["stages"] else {}
    stage_contract_match = bool(clean_stages) and {
        stage_id: {
            "applicable": value["applicable"],
            "artifacts": value["artifacts"],
        }
        for stage_id, value in clean_stages.items()
    } == {
        stage_id: {
            "applicable": value["applicable"],
            "artifacts": value["artifacts"],
        }
        for stage_id, value in resume_stages.items()
    }
    same_tools = clean["tools"] is not None and clean["tools"] == resume["tools"]
    result = {
        "same_project_identity": same_identity,
        "resume_points_to_clean_run": immediate_resume,
        "same_run_schema": same_schema,
        "same_stage_and_artifact_contract": stage_contract_match,
        "same_tool_version_records": same_tools,
    }
    return all(result.values()), result


def evaluate_release_gate(clean_run: Path, resume_run: Path) -> dict:
    """Return privacy-safe acceptance evidence for two explicit run manifests."""
    checks = _Checks()
    clean_internal, clean_public = _inspect_attempt(
        clean_run, "clean", "complete", checks
    )
    resume_internal, resume_public = _inspect_attempt(
        resume_run, "resume", "resumed", checks
    )
    if clean_internal is not None and resume_internal is not None:
        comparison_passed, comparison = _compare_attempts(clean_internal, resume_internal)
        checks.add(
            "pair.identity_and_resume",
            comparison_passed,
            "Attempts share one identity and exact stage/tool contracts; resume is the immediate successor.",
            "Attempts differ in identity, schema, stage/tool contracts, or predecessor linkage.",
        )
    else:
        comparison = {"validated": False}
        checks.add(
            "pair.identity_and_resume",
            False,
            "",
            "Pair comparison requires two valid completed run manifests.",
        )
    status = "PASS" if all(value["status"] == "PASS" for value in checks.records) else "FAIL"
    return {
        "release_gate_schema_version": GATE_SCHEMA_VERSION,
        "status": status,
        "clean_run": clean_public,
        "resume_run": resume_public,
        "comparison": comparison,
        "checks": checks.records,
        "scope": "CFTK core workflow execution, provenance, artifacts, and resume behavior",
        "threshold_policy": "No biological or numerical acceptance thresholds were applied.",
        "privacy": (
            "Source paths, commands, sample-derived filenames, sample identifiers, and "
            "patient-level values are intentionally omitted."
        ),
    }


def _render_text(report: dict) -> str:
    lines = [
        f"CFTK lab release acceptance: {report['status']}",
        f"Gate schema: {report['release_gate_schema_version']}",
    ]
    for label, key in (("Clean run", "clean_run"), ("Resume run", "resume_run")):
        value = report[key]
        if value.get("available"):
            lines.append(
                f"{label}: {value['run_id']} (core schema {value['run_schema_version']}, "
                f"status {value['status']})"
            )
        else:
            lines.append(f"{label}: unavailable or invalid")
    lines.extend(["", "Checks:"])
    lines.extend(
        f"[{record['status']}] {record['id']}: {record['summary']}"
        for record in report["checks"]
    )
    lines.extend([
        "",
        report["threshold_policy"],
        report["privacy"],
    ])
    return "\n".join(lines) + "\n"


def _atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(content, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_release_evidence(report: dict, output_dir: Path) -> tuple[Path, Path]:
    """Write deterministic machine-readable and human-readable gate evidence."""
    output_dir = Path(output_dir).expanduser().resolve()
    json_path = output_dir / "release_acceptance.json"
    text_path = output_dir / "release_acceptance.txt"
    _atomic_write(json_path, json.dumps(report, indent=2, sort_keys=True) + "\n")
    _atomic_write(text_path, _render_text(report))
    return json_path, text_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean-run", required=True, type=Path, help="Completed clean run.json")
    parser.add_argument("--resume-run", required=True, type=Path, help="Completed resume run.json")
    parser.add_argument(
        "--output-dir", required=True, type=Path,
        help="Private directory for release_acceptance.json and .txt",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = evaluate_release_gate(args.clean_run, args.resume_run)
    json_path, text_path = write_release_evidence(report, args.output_dir)
    print(f"release acceptance {report['status']}: {json_path}")
    print(f"human summary: {text_path}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
