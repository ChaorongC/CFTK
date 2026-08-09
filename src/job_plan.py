"""Scheduler-neutral per-sample task-plan generation."""

from __future__ import annotations

import json
import os
import shlex
from datetime import datetime, timezone
from pathlib import Path

from init import get_all_samples, get_bam, get_work_paths, load_config


_STAGES = ("occupancy", "wps", "delfi", "end_motif", "cleavage")
_PROCESS_STAGES = tuple(f"process.{step}" for step in range(1, 5))
_CORE_STAGES = _PROCESS_STAGES + ("qc.2",)
_QC_SAMPLE_STAGES = ("qc.2",)
_FINALIZER_MARKER_SCHEMA = 1


def _timestamp():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def _plan_id(plan):
    return plan.get("plan_id") or Path(plan.get("plan_path", "")).parent.name


def finalizer_marker_path(paths, plan_id, stage):
    """Return the provenance marker path for a generated plan finalizer."""
    if not plan_id:
        return None
    return (
        Path(paths["provenance"])
        / "job-finalizers"
        / str(plan_id)
        / f"{stage}.json"
    )


def write_finalizer_marker(paths, plan_id, stage, *, task_count=None):
    """Record successful cohort finalization without changing scientific outputs."""
    marker = finalizer_marker_path(paths, plan_id, stage)
    if marker is None:
        return None
    marker.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "marker_schema_version": _FINALIZER_MARKER_SCHEMA,
        "plan_id": str(plan_id),
        "stage": stage,
        "task_count": task_count,
        "completed_at": datetime.now(timezone.utc).isoformat(),
    }
    temporary = marker.with_name(f".{marker.name}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, marker)
    return marker


def _command(config_path, stage, *, sample=None, finalize=False, scope=None):
    command = ["cftk", "--config", str(config_path), "frag", f"--{stage.replace('_', '-')}"]
    if sample is not None:
        command.extend(("--sample", str(sample), "--parallel", "1", "--no-finalize"))
    if finalize:
        command.append("--finalize")
    if scope is not None:
        command.extend(("--fragmentomics-scope", str(scope)))
    return command


def _core_command(config_path, stage, *, sample=None, finalize=False, plan_id=None):
    kind, step = stage.split(".", 1)
    command = ["cftk", "--config", str(config_path), kind, "-s", step]
    command.extend(("--parallel", "1"))
    if sample is not None:
        command.extend(("--sample", str(sample), "--no-finalize"))
    if finalize:
        command.append("--finalize")
        if plan_id:
            command.extend(("--job-plan-id", str(plan_id)))
    return command


def _write_script(path, command):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\nexec "
        + shlex.join([str(value) for value in command])
        + "\n",
        encoding="utf-8",
    )
    path.chmod(0o700)


def _write_slurm_helpers(root, stages):
    """Render optional Slurm helpers without submitting or requiring Slurm."""

    root = Path(root)
    slurm = root / "slurm"
    runner = slurm / "run-array-task.sh"
    runner.parent.mkdir(parents=True, exist_ok=True)
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "task_dir=${1:?task directory is required}\n"
        "printf -v prefix '%06d-' \"${SLURM_ARRAY_TASK_ID:?}\"\n"
        "shopt -s nullglob\n"
        "matches=(\"$task_dir\"/\"$prefix\"*.sh)\n"
        "[ ${#matches[@]} -eq 1 ] || { echo \"No unique task for index ${SLURM_ARRAY_TASK_ID}\" >&2; exit 2; }\n"
        "exec \"${matches[0]}\"\n",
        encoding="utf-8",
    )
    runner.chmod(0o700)

    submit = slurm / "submit.sh"
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "# Set CFTK_SBATCH_ARGS for site-specific resources, e.g. '--cpus-per-task=8 --mem=48G'.",
        "# This file is generated only; CFTK never submits jobs itself.",
        f"root={shlex.quote(str(root))}",
        "sbatch_args=${CFTK_SBATCH_ARGS:-}",
    ]
    has_previous_finalizer = False
    for item in stages:
        if isinstance(item, dict):
            stage = item["stage"]
            count = item["count"]
            chain = bool(item.get("chain", False))
        else:
            stage, count = item
            chain = False
        dependency = (
            'dependency="--dependency=afterok:${previous_finalizer}"'
            if chain and has_previous_finalizer
            else 'dependency=""'
        )
        lines.extend([
            dependency,
            f"job_id=$(sbatch --parsable $sbatch_args $dependency --array=0-{count - 1} \"$root/slurm/run-array-task.sh\" \"$root/tasks/{stage}\")",
            "job_id=${job_id%%;*}",
            f"final_id=$(sbatch --parsable $sbatch_args --dependency=afterok:${{job_id}} \"$root/finalize/{stage}.sh\")",
            "final_id=${final_id%%;*}",
            f"printf '%s sample array: %s; finalizer: %s\\n' {shlex.quote(stage)} \"$job_id\" \"$final_id\"",
            "previous_finalizer=$final_id",
        ])
        has_previous_finalizer = True
    submit.write_text("\n".join(lines) + "\n", encoding="utf-8")
    submit.chmod(0o700)
    return submit


def _normalize_core_stages(workflow, requested):
    defaults = {
        "core": _CORE_STAGES,
        "process": _PROCESS_STAGES,
        "qc": _QC_SAMPLE_STAGES,
    }
    if workflow not in defaults:
        raise SystemExit(
            "[plan] ERROR: core per-sample execution requires --workflow core, "
            "process, or qc"
        )
    requested = list(requested or defaults[workflow])
    normalized = []
    for value in requested:
        token = str(value).strip().lower().replace("_", ".")
        if workflow == "process" and token.isdigit():
            token = f"process.{token}"
        elif workflow == "qc" and token.isdigit():
            token = f"qc.{token}"
        if workflow == "qc" and token.startswith("qc.") and token not in _QC_SAMPLE_STAGES:
            raise SystemExit(
                f"[plan] ERROR: {token} is cohort-level and cannot be split into "
                "one-sample jobs; run it once after sample tasks finish."
            )
        if token not in _CORE_STAGES:
            raise SystemExit(
                f"[plan] ERROR: unsupported {workflow} per-sample stage {value!r}. "
                f"Use process.1-process.4 or qc.2."
            )
        if workflow == "process" and not token.startswith("process."):
            raise SystemExit(f"[plan] ERROR: {token} is not a process stage")
        if workflow == "qc" and token not in _QC_SAMPLE_STAGES:
            raise SystemExit(
                f"[plan] ERROR: {token} is cohort-level and cannot be split into "
                "one-sample jobs; run it once after sample tasks finish."
            )
        if token not in normalized:
            normalized.append(token)
    order = {stage: index for index, stage in enumerate(_CORE_STAGES)}
    return sorted(normalized, key=order.__getitem__)


def _core_stage_samples(stage, samples):
    if stage in {"process.1", "process.2"}:
        return [sample for sample in samples if sample.get("input_type") == "fastq"]
    return list(samples)


def write_core_job_plan(args):
    """Write process/QC sample tasks and success-gated cohort finalizers."""

    workflow = getattr(args, "workflow", None) or "core"
    stages = _normalize_core_stages(workflow, getattr(args, "stages", None))
    config_path = Path(args.config).expanduser().resolve()
    cfg = load_config(str(config_path))
    samples = get_all_samples(cfg)
    if not samples:
        raise SystemExit("[plan] ERROR: the project has no samples")
    paths = get_work_paths(cfg)
    plan_id = f"{workflow}-{_timestamp()}"
    root = Path(paths["provenance"]).resolve() / "job-plans" / plan_id
    root.mkdir(parents=True, exist_ok=False)

    tasks = []
    finalizers = []
    skipped_stages = []
    previous_finalizer = None
    slurm_stages = []
    for stage in stages:
        stage_samples = _core_stage_samples(stage, samples)
        if not stage_samples:
            skipped_stages.append(stage)
            continue
        stage_task_ids = []
        for index, sample in enumerate(stage_samples):
            command = _core_command(config_path, stage, sample=sample["name"])
            script = root / "tasks" / stage / f"{index:06d}-{sample['name']}.sh"
            _write_script(script, command)
            task_id = f"{stage}:{sample['name']}"
            dependencies = [previous_finalizer] if previous_finalizer else []
            tasks.append({
                "id": task_id,
                "stage": stage,
                "sample": sample["name"],
                "argv": command,
                "script": str(script),
                "depends_on": dependencies,
            })
            stage_task_ids.append(task_id)

        finalize_command = _core_command(
            config_path, stage, finalize=True, plan_id=plan_id
        )
        finalize_script = root / "finalize" / f"{stage}.sh"
        finalize_script.parent.mkdir(parents=True, exist_ok=True)
        finalize_script.write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            + "exec " + shlex.join(finalize_command) + "\n",
            encoding="utf-8",
        )
        finalize_script.chmod(0o700)
        finalizer_id = f"{stage}:finalize"
        finalizers.append({
            "id": finalizer_id,
            "stage": stage,
            "argv": finalize_command,
            "script": str(finalize_script),
            "depends_on": stage_task_ids,
            "marker": str(finalizer_marker_path(paths, plan_id, stage)),
        })
        previous_finalizer = finalizer_id
        slurm_stages.append({
            "stage": stage,
            "count": len(stage_samples),
            "chain": bool(tasks and len(slurm_stages)),
        })

    plan = {
        "job_plan_schema_version": 2,
        "workflow": workflow,
        "execution_mode": "per-sample",
        "plan_id": plan_id,
        "config": str(config_path),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stages": stages,
        "skipped_stages": skipped_stages,
        "sample_count": len(samples),
        "tasks": tasks,
        "finalizers": finalizers,
        "submission": {
            "automatic": False,
            "note": "CFTK writes task scripts but never submits scheduler jobs.",
        },
    }
    plan_path = root / "job-plan.json"
    plan["plan_path"] = str(plan_path)
    plan["task_count"] = len(tasks)
    plan["finalizer_count"] = len(finalizers)
    if getattr(args, "slurm", False):
        plan["slurm_submit_script"] = str(_write_slurm_helpers(root, slurm_stages))
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return plan


def _present(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def _artifact_groups(stage, sample, cfg, paths):
    """Describe the established per-sample artifacts used by status reporting."""
    name = sample["name"]
    if stage == "process.1":
        ext = "fq.gz" if sample.get("r1", "").endswith(".gz") else "fq"
        trimming = paths["trimming"]
        return [
            {"label": f"{name} trimmed R1", "any_of": [
                os.path.join(trimming, f"{name}_R1.{ext}"),
                os.path.join(trimming, f"{name}_val_1.{ext}"),
            ]},
            {"label": f"{name} trimmed R2", "any_of": [
                os.path.join(trimming, f"{name}_R2.{ext}"),
                os.path.join(trimming, f"{name}_val_2.{ext}"),
            ]},
        ]
    if stage == "process.2":
        bam = os.path.join(paths["alignment"], f"{name}.bam")
        return [{"label": f"{name} alignment", "paths": [
            bam, f"{bam}.bai", f"{bam}.flagstat", f"{bam}.stats",
        ]}]
    if stage == "process.3":
        bam = os.path.join(paths["markdup"], f"{name}.markdup.bam")
        return [
            {"label": f"{name} duplicate-marked BAM", "paths": [bam, f"{bam}.bai"]},
            {"label": f"{name} Picard metrics", "mode": "picard", "root": paths["markdup"], "sample": name},
        ]
    if stage == "process.4":
        prefix = os.path.join(paths["methylation"], name)
        return [{"label": f"{name} CpG calls", "paths": [f"{prefix}_CpG.bedGraph"]}]
    if stage == "qc.2":
        bam = get_bam(sample, paths)
        stem = os.path.splitext(os.path.basename(bam))[0]
        output = paths["qc"]
        return [{"label": f"{name} fragment length", "paths": [
            os.path.join(output, "2_fragment_length", f"fragment_length.{stem}.raw.csv"),
            os.path.join(output, "2_fragment_length", f"fragment_length.{stem}.hist.png"),
        ]}]
    if stage.startswith("fragmentomics."):
        kind = stage.split(".", 1)[1]
        frag = cfg.get("analysis", {}).get("frag", {})
        kmer = frag.get("end_motif", {}).get("params", {}).get("kmer", 4)
        output_specs = {
            "occupancy": (paths["occ_out"], f"{name}.occupancy.tsv"),
            "wps": (paths["wps_out"], f"{name}.wps.tsv"),
            "delfi": (paths["delfi_out"], f"{name}_delfi.tsv"),
            "end_motif": (paths["end_motif_out"], f"{name}_{kmer}mer.tsv"),
            "cleavage": (paths["cleavage_out"], f"{name}_cleavage.bw"),
        }
        if kind in output_specs:
            output_dir, filename = output_specs[kind]
            return [{"label": f"{name} {kind}", "paths": [os.path.join(output_dir, filename)]}]
    return []


def _group_complete(group):
    if group.get("mode") == "picard":
        root = Path(group["root"]) / "picard_metrics"
        name = group["sample"]
        for hs in root.rglob(f"{name}.hs_metrics.txt"):
            if all(_present(hs.parent / filename) for filename in (
                f"{name}.per_target_coverage.txt",
            )) and (hs.parent / f"{name}.multiple_metrics.done").is_file():
                return True
        return False
    if group.get("any_of"):
        return any(_present(path) for path in group["any_of"])
    return all(_present(path) for path in group.get("paths", ()))


def _status_task(task, samples, cfg, paths):
    sample = samples.get(task.get("sample"))
    if sample is None:
        return {
            **task,
            "status": "unknown",
            "missing": [f"sample {task.get('sample')!r} is not in the plan config"],
        }
    groups = _artifact_groups(task.get("stage", ""), sample, cfg, paths)
    if not groups:
        return {**task, "status": "unknown", "missing": ["artifact contract unavailable"]}
    missing = [group["label"] for group in groups if not _group_complete(group)]
    return {**task, "status": "complete" if not missing else "missing", "missing": missing}


def _read_marker(path, plan_id, stage):
    if not path or not Path(path).is_file():
        return False
    try:
        marker = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return (
        marker.get("marker_schema_version") == _FINALIZER_MARKER_SCHEMA
        and marker.get("plan_id") == plan_id
        and marker.get("stage") == stage
    )


def summarize_job_plan(plan_path):
    """Return observed artifact/finalizer state for a generated job plan."""
    plan_path = Path(plan_path).expanduser().resolve()
    try:
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"[status] ERROR: could not read job plan {plan_path}: {exc}") from exc
    config_path = Path(plan["config"]).expanduser().resolve()
    try:
        cfg = load_config(
            str(config_path), acquire_references=False,
            verify_profile=False, validate_profile=False,
        )
    except (KeyError, OSError, ValueError, SystemExit) as exc:
        raise SystemExit(f"[status] ERROR: could not resolve plan config {config_path}: {exc}") from exc
    paths = get_work_paths(cfg)
    samples = {sample["name"]: sample for sample in get_all_samples(cfg)}
    plan_id = _plan_id({**plan, "plan_path": str(plan_path)})
    tasks = [_status_task(task, samples, cfg, paths) for task in plan.get("tasks", [])]
    finalizers = []
    for finalizer in plan.get("finalizers", []):
        stage = finalizer.get("stage", "")
        stage_tasks = [task for task in tasks if task.get("stage") == stage]
        marker = finalizer.get("marker") or str(finalizer_marker_path(paths, plan_id, stage))
        artifacts_complete = stage_tasks and all(
            task["status"] == "complete" for task in stage_tasks
        )
        if _read_marker(marker, plan_id, stage) and artifacts_complete:
            state = "complete"
        elif _read_marker(marker, plan_id, stage):
            state = "stale"
        elif artifacts_complete:
            state = "ready"
        else:
            state = "pending"
        finalizers.append({
            **finalizer,
            "status": state,
            "marker": marker,
            "task_count": len(stage_tasks),
            "complete_tasks": sum(task["status"] == "complete" for task in stage_tasks),
            "missing_tasks": [task["id"] for task in stage_tasks if task["status"] != "complete"],
        })
    if not finalizers:
        overall = "no_finalizers"
    elif all(finalizer["status"] == "complete" for finalizer in finalizers):
        overall = "complete"
    elif any(finalizer["status"] == "stale" for finalizer in finalizers):
        overall = "requires_attention"
    elif any(finalizer["status"] == "ready" for finalizer in finalizers):
        overall = "ready_to_finalize"
    else:
        overall = "in_progress"
    return {
        "status_schema_version": 1,
        "status": overall,
        "plan": str(plan_path),
        "plan_id": plan_id,
        "workflow": plan.get("workflow", "unknown"),
        "config": str(config_path),
        "tasks": tasks,
        "finalizers": finalizers,
        "skipped_stages": list(plan.get("skipped_stages", [])),
        "task_count": len(tasks),
        "complete_task_count": sum(task["status"] == "complete" for task in tasks),
        "note": (
            "Artifact state only; consult Slurm/PBS/SGE for queue, exit, or "
            "scheduler failure state."
        ),
    }


def find_job_plan(config_path, workflow=None):
    """Find the newest generated job plan for the configured project."""
    try:
        cfg = load_config(
            str(Path(config_path).expanduser().resolve()), acquire_references=False,
            verify_profile=False, validate_profile=False,
        )
    except (KeyError, OSError, ValueError, SystemExit) as exc:
        raise SystemExit(f"[status] ERROR: could not resolve project config: {exc}") from exc
    root = Path(get_work_paths(cfg)["provenance"]) / "job-plans"
    candidates = [path for path in root.glob("*/job-plan.json") if path.is_file()]
    if workflow:
        candidates = [path for path in candidates if path.parent.name.startswith(f"{workflow}-")]
    if not candidates:
        qualifier = f" for workflow {workflow!r}" if workflow else ""
        raise SystemExit(
            "[status] ERROR: no generated per-sample job plan found"
            f"{qualifier}. The normal beginner path is 'cftk run --parallel N'; "
            "status applies to advanced job plans."
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def render_job_plan_status(status):
    """Render concise human-readable status without claiming scheduler state."""
    lines = [
        f"CFTK advanced job-plan status: {status['status'].upper()}",
        f"Plan: {status['plan']}",
        f"Workflow: {status['workflow']}",
        f"Tasks: {status['complete_task_count']}/{status['task_count']} artifacts complete",
    ]
    for stage in status["skipped_stages"]:
        lines.append(f"  {stage}: skipped by input type")
    for finalizer in status["finalizers"]:
        lines.append(
            f"  {finalizer['stage']}: {finalizer['complete_tasks']}/{finalizer['task_count']} "
            f"tasks complete; finalizer {finalizer['status']}"
        )
        if finalizer["missing_tasks"]:
            lines.append(f"    waiting: {', '.join(finalizer['missing_tasks'][:3])}")
    lines.append(f"Note: {status['note']}")
    return "\n".join(lines)


def write_job_plan(args):
    workflow = getattr(args, "workflow", None)
    stages = tuple(getattr(args, "stages", ()) or ())
    if workflow in {"core", "process", "qc"}:
        return write_core_job_plan(args)
    if workflow is None and stages and any(str(stage).startswith(("process", "qc")) for stage in stages):
        return write_core_job_plan(args)
    return write_fragmentomics_job_plan(args)


def write_fragmentomics_job_plan(args):
    """Write independent sample tasks and after-all finalizers for selected stages."""

    requested = tuple(dict.fromkeys(getattr(args, "stages", ()) or ()))
    invalid = [stage for stage in requested if stage not in _STAGES]
    if not requested or invalid:
        choices = ", ".join(_STAGES)
        raise SystemExit(f"[job-plan] ERROR: choose one or more fragmentomics stages: {choices}")

    config_path = Path(args.config).expanduser().resolve()
    cfg = load_config(str(config_path))
    samples = get_all_samples(cfg)
    if not samples:
        raise SystemExit("[job-plan] ERROR: the project has no samples")
    paths = get_work_paths(cfg)
    plan_id = f"fragmentomics-{_timestamp()}"
    root = Path(paths["provenance"]).resolve() / "job-plans" / plan_id
    root.mkdir(parents=True, exist_ok=False)

    scope = getattr(args, "fragmentomics_scope", None)
    tasks = []
    finalizers = []
    for stage in requested:
        for index, sample in enumerate(samples):
            command = _command(config_path, stage, sample=sample["name"], scope=scope)
            script = root / "tasks" / stage / f"{index:06d}-{sample['name']}.sh"
            _write_script(script, command)
            tasks.append({
                "id": f"fragmentomics.{stage}:{sample['name']}",
                "stage": f"fragmentomics.{stage}",
                "sample": sample["name"],
                "argv": command,
                "script": str(script),
                "depends_on": [],
            })
        finalize_command = _command(config_path, stage, finalize=True, scope=scope)
        adoption_command = [
            "cftk", "--config", str(config_path), "analyze", "--stage", stage,
            "--adopt-existing", "--job-plan-id", plan_id,
        ]
        if scope is not None:
            adoption_command.extend(("--fragmentomics-scope", str(scope)))
        finalize_script = root / "finalize" / f"{stage}.sh"
        finalize_script.parent.mkdir(parents=True, exist_ok=True)
        finalize_script.write_text(
            "#!/usr/bin/env bash\nset -euo pipefail\n"
            + shlex.join(finalize_command) + "\n"
            + "exec " + shlex.join(adoption_command) + "\n",
            encoding="utf-8",
        )
        finalize_script.chmod(0o700)
        finalizers.append({
            "id": f"fragmentomics.{stage}:finalize",
            "stage": f"fragmentomics.{stage}",
            "argv": [finalize_command, adoption_command],
            "script": str(finalize_script),
            "depends_on": [task["id"] for task in tasks if task["stage"] == f"fragmentomics.{stage}"],
            "marker": str(finalizer_marker_path(paths, plan_id, f"fragmentomics.{stage}")),
        })

    plan = {
        "job_plan_schema_version": 1,
        "workflow": "fragmentomics",
        "execution_mode": "per-sample",
        "plan_id": plan_id,
        "config": str(config_path),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "stages": list(requested),
        "sample_count": len(samples),
        "tasks": tasks,
        "finalizers": finalizers,
        "submission": {
            "automatic": False,
            "note": "CFTK writes task scripts but never submits scheduler jobs.",
        },
    }
    plan_path = root / "job-plan.json"
    plan["plan_path"] = str(plan_path)
    plan["task_count"] = len(tasks)
    plan["finalizer_count"] = len(finalizers)
    if getattr(args, "slurm", False):
        plan["slurm_submit_script"] = str(_write_slurm_helpers(
            root, [(stage, len(samples)) for stage in requested]
        ))
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return plan
