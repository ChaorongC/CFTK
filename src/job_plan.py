"""Scheduler-neutral per-sample fragmentomics task-plan generation."""

from __future__ import annotations

import json
import os
import shlex
from datetime import datetime, timezone
from pathlib import Path

from init import get_all_samples, get_work_paths, load_config


_STAGES = ("occupancy", "wps", "delfi", "end_motif", "cleavage")


def _timestamp():
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def _command(config_path, stage, *, sample=None, finalize=False, scope=None):
    command = ["cftk", "--config", str(config_path), "frag", f"--{stage.replace('_', '-')}"]
    if sample is not None:
        command.extend(("--sample", str(sample), "--parallel", "1", "--no-finalize"))
    if finalize:
        command.append("--finalize")
    if scope is not None:
        command.extend(("--fragmentomics-scope", str(scope)))
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
    for stage, count in stages:
        lines.extend([
            f"job_id=$(sbatch --parsable $sbatch_args --array=0-{count - 1} \"$root/slurm/run-array-task.sh\" \"$root/tasks/{stage}\")",
            "job_id=${job_id%%;*}",
            f"final_id=$(sbatch --parsable $sbatch_args --dependency=afterok:${{job_id}} \"$root/finalize/{stage}.sh\")",
            "final_id=${final_id%%;*}",
            f"printf '%s sample array: %s; finalizer: %s\\n' {shlex.quote(stage)} \"$job_id\" \"$final_id\"",
        ])
    submit.write_text("\n".join(lines) + "\n", encoding="utf-8")
    submit.chmod(0o700)
    return submit


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
    root = Path(paths["provenance"]).resolve() / "job-plans" / f"fragmentomics-{_timestamp()}"
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
            "--adopt-existing",
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
        })

    plan = {
        "job_plan_schema_version": 1,
        "workflow": "fragmentomics",
        "execution_mode": "per-sample",
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
