"""Fail-fast beginner workflow with validated resume and run provenance."""

import csv
import hashlib
import html
import json
import os
import shlex
import shutil
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

from init import get_all_samples, get_sequence_dictionary_path, get_work_paths, load_config
from resource_planning import detect_scheduler_allocation, plan_parallelism
from util import configure_command_log, disp


RUN_SCHEMA_VERSION = 3
DEFAULT_TOOLS = {
    "step1_trimming": "trim_galore",
    "step2_alignment": "bwameth",
    "step3_markdup": "sambamba",
    "step4_methylation": "methyldackel",
}


class RunContractError(RuntimeError):
    """A beginner-run safety or artifact contract was not satisfied."""


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _new_run_id():
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_json(path, value):
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _append_event(path, event, **details):
    record = {"event": event, "timestamp": _utc_now(), **details}
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")


def _validate_project_contract(raw, cfg, samples, lock_path):
    if raw.get("schema_version") != 2:
        raise RunContractError(
            "cftk run requires a schema-v2 project. Use the expert process/QC "
            "commands for legacy configurations."
        )
    if not Path(lock_path).is_file():
        raise RunContractError(
            f"cftk run requires cftk.lock.json beside the project config: {lock_path}. "
            "Rerun 'cftk init' to validate and lock the project."
        )
    input_types = {sample.get("input_type") for sample in samples}
    if input_types == {"fastq", "bam"}:
        raise RunContractError(
            "cftk run does not mix FASTQ and BAM samples in one beginner run. "
            "Use homogeneous inputs or the expert step commands."
        )
    if len(input_types) != 1 or next(iter(input_types), None) not in {"fastq", "bam"}:
        raise RunContractError("cftk run requires homogeneous FASTQ or BAM inputs.")
    configured = {
        key: cfg.get("process", {}).get(key, {}).get("tool")
        for key in DEFAULT_TOOLS
    }
    mismatches = {
        key: {"expected": DEFAULT_TOOLS[key], "configured": configured[key]}
        for key in DEFAULT_TOOLS
        if configured[key] != DEFAULT_TOOLS[key]
    }
    if mismatches:
        detail = ", ".join(
            f"{key}={value['configured']!r} (expected {value['expected']!r})"
            for key, value in mismatches.items()
        )
        raise RunContractError(
            "cftk run supports only the validated toolchain in schema v2: " + detail
        )
    return next(iter(input_types))


def _load_context(args):
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        raise RunContractError(f"Project config not found: {config_path}")
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RunContractError(f"Could not read project config {config_path}: {exc}") from exc
    lock_path = config_path.with_name("cftk.lock.json")
    if raw.get("schema_version") != 2:
        _validate_project_contract(raw, {}, [], lock_path)
    cfg = load_config(str(config_path))
    samples = get_all_samples(cfg)
    input_type = _validate_project_contract(raw, cfg, samples, lock_path)
    paths = get_work_paths(cfg)
    return {
        "config_path": config_path,
        "lock_path": lock_path,
        "raw": raw,
        "cfg": cfg,
        "samples": samples,
        "paths": paths,
        "identity": {
            "config_sha256": _sha256(config_path),
            "lock_sha256": _sha256(lock_path),
        },
        "input_type": input_type,
    }


def _build_stage_plan(input_type, include_dinucleotide=False):
    stages = [
        {"id": "process.1", "name": "Adapter trimming and FastQC", "kind": "process", "step": 1,
         "applicable": input_type == "fastq"},
        {"id": "process.2", "name": "Bisulfite alignment", "kind": "process", "step": 2,
         "applicable": input_type == "fastq"},
        {"id": "process.3", "name": "Duplicate marking and Picard metrics", "kind": "process", "step": 3,
         "applicable": True},
        {"id": "process.4", "name": "CpG methylation and matrix merge", "kind": "process", "step": 4,
         "applicable": True},
        {"id": "qc.2", "name": "Fragment-length QC", "kind": "qc", "step": 2,
         "applicable": True},
        {"id": "qc.0", "name": "QC metrics and scores", "kind": "qc", "step": 0,
         "applicable": True},
        {"id": "qc.1", "name": "Methylation-distribution QC", "kind": "qc", "step": 1,
         "applicable": True},
    ]
    if include_dinucleotide:
        stages.append({
            "id": "qc.3", "name": "Dinucleotide-frequency QC", "kind": "qc", "step": 3,
            "applicable": True,
        })
    return stages


def _build_resource_plan(context, stages, args):
    process = context["cfg"]["process"]
    sample_count = len(context["samples"])
    requested_parallel = args.parallel or process.get("parallel_samples", 1)
    scheduler = detect_scheduler_allocation()
    records = []
    step_keys = {
        1: "step1_trimming",
        2: "step2_alignment",
        3: "step3_markdup",
        4: "step4_methylation",
    }
    try:
        for stage in stages:
            if stage["kind"] == "process":
                total = process[step_keys[stage["step"]]]["params"].get("cores", 20)
                record = {
                    "stage": stage["id"],
                    "model": "parallel samples with per-sample tool threads",
                    **plan_parallelism(total, requested_parallel, sample_count),
                }
            else:
                total = process["step4_methylation"]["params"].get("cores", 20)
                sample_plan = plan_parallelism(total, requested_parallel, sample_count)
                if stage["id"] == "qc.2":
                    record = {
                        "stage": stage["id"],
                        "model": "parallel samples with per-sample tool threads",
                        **sample_plan,
                    }
                elif stage["id"] == "qc.0":
                    record = {
                        "stage": stage["id"],
                        "model": "sample parsing threads",
                        **sample_plan,
                        "concurrent_samples": min(total, max(1, sample_count)),
                        "threads_per_sample": 1,
                        "estimated_peak_threads": min(total, max(1, sample_count)),
                    }
                elif stage["id"] == "qc.3":
                    record = {
                        "stage": stage["id"],
                        "model": "separate sample-extraction and pattern-worker phases",
                        **sample_plan,
                        "threads_per_sample": 1,
                        "pattern_workers": min(total, 8),
                        "estimated_peak_threads": max(
                            sample_plan["concurrent_samples"], min(total, 8)
                        ),
                    }
                else:
                    record = {
                        "stage": stage["id"],
                        "model": "single Python plotting process",
                        **sample_plan,
                        "concurrent_samples": 1,
                        "threads_per_sample": 1,
                        "estimated_peak_threads": 1,
                    }
            record["applicable"] = stage["applicable"]
            records.append(record)
    except (KeyError, ValueError) as exc:
        raise RunContractError(f"Invalid CPU resource plan: {exc}") from exc

    maximum = max(
        (record["total_core_budget"] for record in records if record["applicable"]),
        default=0,
    )
    allocated = scheduler.get("allocated_cores")
    return {
        "resource_plan_version": 1,
        "sample_count": sample_count,
        "requested_parallel_samples": requested_parallel,
        "maximum_total_core_budget": maximum,
        "scheduler": scheduler,
        "scheduler_compatible": allocated is None or maximum <= allocated,
        "stages": records,
    }


def _spec(path, description, *, role="output", nonempty=True, required=True):
    return {
        "path": str(Path(path).resolve()),
        "description": description,
        "role": role,
        "required": required,
        "nonempty": nonempty,
    }


def _picard_metrics_dir(context, args):
    cfg = context["cfg"]
    paths = context["paths"]
    target = Path(args.target_bed or cfg["reference_data"].get("target_bed", "")).resolve()
    sequence_dict = Path(get_sequence_dictionary_path(cfg["reference_data"]["genome_fa"]))
    stem = target.name[:-4] if target.name.endswith(".bed") else target.name
    digest = hashlib.sha256()
    for input_path in (target, sequence_dict):
        with input_path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return Path(paths["markdup"]) / "picard_metrics" / f"{stem}.{digest.hexdigest()[:12]}"


def _artifact_specs(context, stage, args):
    paths = context["paths"]
    samples = context["samples"]
    stage_id = stage["id"]
    specs = []
    if stage_id == "process.1":
        for sample in samples:
            ext = "fq.gz" if sample["r1"].endswith(".gz") else "fq"
            name = sample["name"]
            for mate in ("R1", "R2"):
                specs.extend([
                    _spec(Path(paths["trimming"]) / f"{name}_{mate}.{ext}", f"trimmed {mate} for {name}"),
                    _spec(Path(paths["trimming"]) / f"{name}_{mate}_trimming_report.txt", f"trimming report for {name} {mate}", role="report"),
                    _spec(Path(paths["trimming"]) / f"{name}_{mate}_fastqc.html", f"FastQC HTML for {name} {mate}", role="report"),
                    _spec(Path(paths["trimming"]) / f"{name}_{mate}_fastqc.zip", f"FastQC data for {name} {mate}"),
                ])
        specs.append(_spec(
            Path(paths["trimming"]) / "multiqc/multiqc_report.html",
            "trimming MultiQC report", role="report",
        ))
    elif stage_id == "process.2":
        for sample in samples:
            bam = Path(paths["alignment"]) / f"{sample['name']}.bam"
            specs.extend([
                _spec(bam, f"aligned BAM for {sample['name']}"),
                _spec(f"{bam}.bai", f"aligned BAM index for {sample['name']}"),
                _spec(f"{bam}.flagstat", f"flagstat for {sample['name']}", role="report"),
                _spec(f"{bam}.stats", f"samtools stats for {sample['name']}", role="report"),
            ])
        specs.append(_spec(
            Path(paths["alignment"]) / "multiqc/multiqc_report.html",
            "alignment MultiQC report", role="report",
        ))
    elif stage_id == "process.3":
        metrics = _picard_metrics_dir(context, args)
        target_stem = metrics.name.rsplit(".", 1)[0]
        specs.append(_spec(metrics / f"{target_stem}.interval_list", "Picard target interval list"))
        for sample in samples:
            name = sample["name"]
            bam = Path(paths["markdup"]) / f"{name}.markdup.bam"
            prefix = metrics / f"{name}.multiple_metrics"
            specs.extend([
                _spec(bam, f"duplicate-marked BAM for {name}"),
                _spec(f"{bam}.bai", f"duplicate-marked BAM index for {name}"),
                _spec(Path(paths["markdup"]) / f"{name}.markdup_metrics.txt", f"Sambamba duplicate metrics for {name}", role="report"),
                _spec(metrics / f"{name}.hs_metrics.txt", f"Picard target metrics for {name}", role="report"),
                _spec(metrics / f"{name}.per_target_coverage.txt", f"per-target coverage for {name}", role="report"),
                _spec(f"{prefix}.alignment_summary_metrics", f"alignment summary for {name}", role="report"),
                _spec(f"{prefix}.insert_size_metrics", f"insert-size metrics for {name}", role="report"),
                _spec(f"{prefix}.gc_bias.summary_metrics", f"GC-bias summary for {name}", role="report"),
                _spec(f"{prefix}.done", f"Picard multiple-metrics checkpoint for {name}", nonempty=False),
                _spec(f"{prefix}.insert_size_histogram.pdf", f"insert-size histogram for {name}", role="figure"),
                _spec(f"{prefix}.gc_bias.pdf", f"GC-bias plot for {name}", role="figure"),
                _spec(f"{prefix}.read_length_histogram.pdf", f"read-length histogram for {name}", role="figure"),
            ])
    elif stage_id == "process.4":
        for sample in samples:
            prefix = Path(paths["methylation"]) / sample["name"]
            specs.extend([
                _spec(f"{prefix}_mbias.txt", f"M-bias table for {sample['name']}", role="report"),
                _spec(f"{prefix}_mbias_OT_OB.temp", f"OT/OB bounds for {sample['name']}", role="report"),
                _spec(f"{prefix}_CpG.bedGraph", f"CpG calls for {sample['name']}"),
                _spec(f"{prefix}_OT.svg", f"OT M-bias plot for {sample['name']}", role="figure"),
                _spec(f"{prefix}_OB.svg", f"OB M-bias plot for {sample['name']}", role="figure"),
            ])
        specs.append(_spec(
            Path(paths["cpg_matrix"]) / "cpg_matrix.tsv", "merged CpG matrix"
        ))
    elif stage_id == "qc.2":
        for sample in samples:
            stem = f"{sample['name']}.markdup"
            prefix = Path(paths["qc"]) / "2_fragment_length" / f"fragment_length.{stem}"
            specs.extend([
                _spec(f"{prefix}.raw.csv", f"raw fragment lengths for {sample['name']}"),
                _spec(f"{prefix}.hist.png", f"per-sample fragment histogram for {sample['name']}", role="figure"),
            ])
        qc2 = Path(paths["qc"]) / "2_fragment_length/fragment_length"
        specs.extend([
            _spec(f"{qc2}.png", "combined fragment-length distribution", role="figure"),
            _spec(f"{qc2}.pdf", "combined fragment-length distribution PDF", role="figure"),
        ])
    elif stage_id == "qc.0":
        specs.extend([
            _spec(Path(paths["qc"]) / "qc_summary.tsv", "QC metrics table", role="report"),
            _spec(Path(paths["qc"]) / "qc_scores.tsv", "QC scores table", role="report"),
        ])
    elif stage_id == "qc.1":
        base = Path(paths["qc"]) / "1_methylation_distribution/methylation_distribution"
        specs.extend([
            _spec(f"{base}.png", "methylation-distribution plot", role="figure"),
            _spec(f"{base}.pdf", "methylation-distribution PDF", role="figure"),
        ])
    elif stage_id == "qc.3":
        base = Path(paths["qc"]) / "3_dinucleotide_freq/dinucleotide_freq"
        specs.extend([
            _spec(f"{base}.png", "dinucleotide-frequency plot", role="figure"),
            _spec(f"{base}.pdf", "dinucleotide-frequency PDF", role="figure"),
        ])
    return specs


def _validate_artifacts(specs):
    errors = []
    for spec in specs:
        if not spec.get("required", True):
            continue
        path = Path(spec["path"])
        if not path.exists():
            errors.append(f"missing {spec['description']}: {path}")
        elif spec.get("nonempty", True) and path.is_file() and path.stat().st_size == 0:
            errors.append(f"empty {spec['description']}: {path}")
    return errors


def _can_resume(previous, identity, stage_id, specs):
    if not previous or previous.get("project_identity") != identity:
        return False
    stage = next((item for item in previous.get("stages", []) if item.get("id") == stage_id), None)
    return bool(
        stage
        and stage.get("status") in {"complete", "resumed", "adopted"}
        and not _validate_artifacts(specs)
    )


def _quarantine_paths(paths, results_root, quarantine_root):
    results_root = Path(results_root).resolve()
    quarantine_root = Path(quarantine_root).resolve()
    moved = []
    for value in sorted({Path(path).resolve() for path in paths}, key=str):
        if not value.exists():
            continue
        try:
            relative = value.relative_to(results_root)
        except ValueError as exc:
            raise RunContractError(f"Refusing to quarantine output outside results: {value}") from exc
        destination = quarantine_root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            destination = destination.with_name(f"{destination.name}.{uuid.uuid4().hex[:8]}")
        shutil.move(str(value), str(destination))
        moved.append({"source": str(value), "destination": str(destination)})
    return moved


def _run_preflight(context, args):
    from doctor import run_doctor

    steps = [3, 4] if context["input_type"] == "bam" else [1, 2, 3, 4]
    report = run_doctor(SimpleNamespace(
        config=str(context["config_path"]),
        step=steps,
        parallel=args.parallel,
        target_bed=args.target_bed,
        skip_picard_metrics=False,
    ))
    for check in report["checks"]:
        if check["id"] == "tool.multiqc" and check["status"] == "WARN":
            check["status"] = "FAIL"
            check["summary"] += " MultiQC is required by cftk run."
    path = shutil.which("bamPEFragmentSize")
    if not path:
        report["checks"].append({
            "id": "tool.bamPEFragmentSize",
            "status": "FAIL",
            "summary": "bamPEFragmentSize is not available on PATH.",
            "remedy": "Install or update the pinned CFTK environment with deepTools.",
        })
    else:
        try:
            completed = subprocess.run(
                [path, "--version"], capture_output=True, text=True, timeout=60
            )
            output = (completed.stdout or completed.stderr).strip().splitlines()
            summary = output[0][:300] if output else "version output unavailable"
            status = "PASS" if completed.returncode == 0 else "FAIL"
        except (OSError, subprocess.TimeoutExpired) as exc:
            status, summary = "FAIL", str(exc)
        report["checks"].append({
            "id": "tool.bamPEFragmentSize", "status": status,
            "summary": f"bamPEFragmentSize: {summary}", "details": {"path": path},
        })
    counts = {
        key.lower(): sum(check["status"] == key for check in report["checks"])
        for key in ("PASS", "WARN", "FAIL")
    }
    report["summary"] = counts
    report["status"] = "FAIL" if counts["fail"] else ("WARN" if counts["warn"] else "PASS")
    report["exit_code"] = 1 if counts["fail"] else 0
    return report


def _tool_versions(report):
    return {
        check["id"].split(".", 1)[1]: {
            "status": check["status"],
            "summary": check["summary"],
            **({"path": check["details"]["path"]} if check.get("details", {}).get("path") else {}),
        }
        for check in report.get("checks", [])
        if check.get("id", "").startswith("tool.")
    }


def _execute_stage(context, stage, args):
    import cftk

    if stage["kind"] == "process":
        args.step = [stage["step"]]
        args.skip_picard_metrics = False
        cftk._cmd_process(args)
    else:
        args.step = [stage["step"]]
        args.force = False
        args.title = None
        cftk._cmd_qc(args)


def _stage_command(context, stage, args):
    command = ["cftk", "--config", str(context["config_path"])]
    if stage["kind"] == "process":
        command.extend(["process", "-s", str(stage["step"])])
        if args.target_bed and stage["step"] == 3:
            command.extend(["--target-bed", str(Path(args.target_bed).expanduser().resolve())])
    else:
        command.extend(["qc", "-s", str(stage["step"])])
    if args.parallel:
        command.extend(["--parallel", str(args.parallel)])
    return shlex.join(command)


def _record_artifacts(stage_record, specs):
    values = [{
        "path": spec["path"],
        "description": spec["description"],
        "role": spec["role"],
    } for spec in specs]
    stage_record["outputs"] = [value for value in values if value["role"] != "figure"]
    stage_record["figures"] = [value for value in values if value["role"] == "figure"]


def _write_artifact_tables(manifest, run_dir):
    outputs = []
    figures = []
    for stage in manifest["stages"]:
        for value in stage.get("expected", []):
            row = {
                "stage": stage["id"], "description": value["description"],
                "path": value["path"], "required": value.get("required", True),
            }
            (figures if value["role"] == "figure" else outputs).append(row)
    for filename, rows in (("expected-outputs.tsv", outputs), ("figures.tsv", figures)):
        path = Path(run_dir) / filename
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=("stage", "description", "path", "required"), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)


def _relative_link(path, base_dir):
    try:
        return os.path.relpath(Path(path), Path(base_dir))
    except ValueError:
        return str(path)


def _summary_fragmentomics_scope(manifest):
    """Merge executed stage-local scope metadata into the run-level summary."""

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


def _write_summary_html(manifest, path, project_root=None):
    link_base = Path(path).parent
    rows = []
    artifact_sections = []
    for stage in manifest.get("stages", []):
        rows.append(
            "<tr>"
            f"<td><code>{html.escape(stage['id'])}</code></td>"
            f"<td>{html.escape(stage.get('name', ''))}</td>"
            f"<td class=\"status {html.escape(stage.get('status', 'pending'))}\">"
            f"{html.escape(stage.get('status', 'pending'))}</td>"
            "</tr>"
        )
        outputs = stage.get("outputs", [])
        figures = stage.get("figures", [])
        if not outputs and not figures and stage.get("status") == "planned":
            outputs = [
                value for value in stage.get("expected", [])
                if value.get("role") != "figure"
            ]
            figures = [
                value for value in stage.get("expected", [])
                if value.get("role") == "figure"
            ]
        links = []
        thumbnails = []
        for value in outputs + figures:
            href = _relative_link(value["path"], link_base)
            links.append(
                f'<li><a href="{html.escape(href, quote=True)}">'
                f"{html.escape(value['description'])}</a></li>"
            )
        for value in figures:
            if Path(value["path"]).suffix.lower() not in {".png", ".svg", ".jpg", ".jpeg"}:
                continue
            href = _relative_link(value["path"], link_base)
            thumbnails.append(
                '<figure>'
                f'<a href="{html.escape(href, quote=True)}">'
                f'<img loading="lazy" src="{html.escape(href, quote=True)}" '
                f'alt="{html.escape(value["description"], quote=True)}"></a>'
                f'<figcaption>{html.escape(value["description"])}</figcaption>'
                '</figure>'
            )
        if links:
            artifact_sections.append(
                f"<h3>{html.escape(stage.get('name', stage['id']))}</h3>"
                f"<ul>{''.join(links)}</ul>"
                f"<div class=\"figures\">{''.join(thumbnails)}</div>"
            )
    record_links = []
    for filename, description in (
        ("run.json", "Machine-readable run manifest"),
        ("resource-plan.json", "Resolved CPU resource plan"),
        ("events.jsonl", "Stage event ledger"),
        ("doctor-before.json", "Preflight doctor report"),
        ("tool-versions.json", "Tool versions"),
        ("commands.jsonl", "Exact command ledger"),
        ("expected-outputs.tsv", "Expected outputs"),
        ("figures.tsv", "Expected figures"),
    ):
        if (link_base / filename).exists():
            record_links.append(
                f'<li><a href="{html.escape(filename, quote=True)}">'
                f"{html.escape(description)}</a></li>"
            )
    evidence = manifest.get("evidence") or {}
    evidence_files = []
    evidence_thumbnails = []
    evidence_dir = Path(evidence["directory"]) if evidence.get("directory") else None
    if evidence_dir:
        for filename in evidence.get("files", []):
            evidence_path = evidence_dir / filename
            if not evidence_path.is_file():
                continue
            href = _relative_link(evidence_path, link_base)
            evidence_files.append(
                f'<li><a href="{html.escape(href, quote=True)}">'
                f"{html.escape(filename)}</a></li>"
            )
            if evidence_path.suffix.lower() in {".png", ".svg", ".jpg", ".jpeg"}:
                evidence_thumbnails.append(
                    "<figure>"
                    f'<a href="{html.escape(href, quote=True)}">'
                    f'<img loading="lazy" src="{html.escape(href, quote=True)}" '
                    f'alt="{html.escape(filename, quote=True)}"></a>'
                    f"<figcaption>{html.escape(filename)}</figcaption>"
                    "</figure>"
                )
    evidence_error = evidence.get("error")
    evidence_block = ""
    if evidence:
        evidence_block = (
            "<h2>Generated evidence</h2>"
            f"<p><strong>Status:</strong> {html.escape(evidence.get('status', 'unknown'))}</p>"
            + (f"<pre>{html.escape(str(evidence_error))}</pre>" if evidence_error else "")
            + (f"<ul>{''.join(evidence_files)}</ul>" if evidence_files else "")
            + (f'<div class="figures">{"".join(evidence_thumbnails)}</div>' if evidence_thumbnails else "")
        )
    resource_rows = []
    for resource in manifest.get("resource_plan", {}).get("stages", []):
        if not resource.get("applicable", True):
            continue
        resource_rows.append(
            "<tr>"
            f"<td><code>{html.escape(resource['stage'])}</code></td>"
            f"<td>{resource['total_core_budget']}</td>"
            f"<td>{resource['concurrent_samples']}</td>"
            f"<td>{resource['threads_per_sample']}</td>"
            f"<td>{resource['estimated_peak_threads']}</td>"
            f"<td>{html.escape(resource['model'])}</td>"
            "</tr>"
        )
    resource_table = (
        "<h2>CPU resource plan</h2>"
        "<table><thead><tr><th>Stage</th><th>Total budget</th>"
        "<th>Concurrent samples</th><th>Threads/sample</th>"
        "<th>Estimated peak</th><th>Execution model</th></tr></thead><tbody>"
        f"{''.join(resource_rows)}</tbody></table>"
        if resource_rows else ""
    )
    error = manifest.get("error")
    error_block = (
        f'<h2>Terminal error</h2><pre>{html.escape(error)}</pre>' if error else ""
    )
    scope = _summary_fragmentomics_scope(manifest)
    scope_block = ""
    if isinstance(scope, dict) and scope:
        scope_rows = []
        for label, key in (
            ("Mode", "mode"),
            ("Requested", "requested"),
            ("Assay", "assay"),
            ("Target BED", "target_bed"),
            ("Target BED SHA-256", "target_sha256"),
            ("WPS/occupancy regions", "region_count"),
            ("DELFI bins", "bins_count"),
            ("Derived scope directory", "scope_root"),
        ):
            value = scope.get(key)
            if value in (None, "", [], {}):
                continue
            scope_rows.append(
                "<tr>"
                f"<th>{html.escape(label)}</th>"
                f"<td><code>{html.escape(str(value))}</code></td>"
                "</tr>"
            )
        scope_block = (
            "<h2>Fragmentomics scope</h2>"
            f"<table><tbody>{''.join(scope_rows)}</tbody></table>"
            f"<p><strong>Interpretation:</strong> {html.escape(str(scope.get('note', '')))}</p>"
        )
    downstream = manifest.get("downstream") or {}
    downstream_block = ""
    if downstream:
        links = []
        for key, label in (("summary", "Downstream HTML summary"), ("manifest", "Downstream manifest")):
            value = downstream.get(key)
            if value:
                href = _relative_link(value, link_base)
                links.append(
                    f'<li><a href="{html.escape(href, quote=True)}">'
                    f"{html.escape(label)}</a></li>"
                )
        downstream_block = (
            "<h2>Downstream workflow</h2>"
            f"<p><strong>Preset:</strong> {html.escape(str(downstream.get('preset', '')))}<br>"
            f"<strong>Status:</strong> {html.escape(str(downstream.get('status', 'unknown')))}</p>"
            f"<ul>{''.join(links)}</ul>"
        )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>CFTK run {html.escape(manifest['run_id'])}</title>
<style>
body{{font-family:Arial,sans-serif;max-width:1100px;margin:32px auto;padding:0 20px;color:#202428}}
h1{{font-size:28px}} table{{border-collapse:collapse;width:100%}} th,td{{border-bottom:1px solid #d8dde2;padding:9px;text-align:left}}
.status{{font-weight:700}} .complete,.resumed,.adopted{{color:#167346}} .complete_with_reporting_error{{color:#9a6700}} .failed{{color:#b42318}} .planned,.pending{{color:#6b7280}}
code{{font-family:ui-monospace,monospace}} a{{color:#075985}}
pre{{white-space:pre-wrap;background:#f6f8fa;border:1px solid #d8dde2;padding:12px}}
.figures{{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:16px}}
figure{{margin:0}} img{{display:block;max-width:100%;height:auto;border:1px solid #d8dde2}}
figcaption{{font-size:14px;margin-top:6px;color:#4b5563}}
</style></head><body>
<h1>CFTK run summary</h1>
<p><strong>Run:</strong> <code>{html.escape(manifest['run_id'])}</code><br>
<strong>Status:</strong> {html.escape(manifest.get('status', 'unknown'))}<br>
<strong>Started:</strong> {html.escape(manifest.get('started_at', ''))}<br>
<strong>Finished:</strong> {html.escape(manifest.get('finished_at') or 'in progress')}</p>
{error_block}
{scope_block}
{downstream_block}
{resource_table}
<h2>Stages</h2><table><thead><tr><th>ID</th><th>Stage</th><th>Status</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
<h2>Run records</h2><ul>{''.join(record_links)}</ul>
{evidence_block}
<h2>Outputs and figures</h2>{''.join(artifact_sections) or '<p>No completed artifacts recorded.</p>'}
</body></html>"""
    _atomic_write_text(path, document)


def _load_previous(provenance, identity=None):
    provenance = Path(provenance).resolve()
    latest = Path(provenance) / "latest-run.json"
    candidates = []
    if latest.is_file():
        try:
            pointer = json.loads(latest.read_text(encoding="utf-8"))
            candidates.append(Path(pointer["manifest"]))
        except (OSError, KeyError, json.JSONDecodeError):
            pass
    candidates.extend(sorted((provenance / "runs").glob("*/run.json"), reverse=True))
    seen = set()
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            candidate.relative_to(provenance / "runs")
            manifest = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if identity is None or manifest.get("project_identity") == identity:
            return manifest
    return None


def _save_attempt(manifest, run_dir, provenance):
    manifest_path = Path(run_dir) / "run.json"
    _atomic_write_json(manifest_path, manifest)
    _atomic_write_json(Path(run_dir) / "resource-plan.json", manifest.get(
        "resource_plan", {"status": "NOT_PLANNED"}
    ))
    _write_summary_html(manifest, Path(run_dir) / "run-summary.html", manifest["project_root"])
    _atomic_write_json(Path(provenance) / "latest-run.json", {
        "run_id": manifest["run_id"],
        "manifest": str(manifest_path.resolve()),
        "summary": str((Path(run_dir) / "run-summary.html").resolve()),
        "status": manifest["status"],
    })


def _generate_evidence(manifest, run_dir):
    """Generate the installed evidence bundle for one terminal attempt."""
    from validation_reports import summarize

    run_dir = Path(run_dir).resolve()
    evidence_dir = run_dir / "evidence"
    summary = summarize(run_dir / "run.json", evidence_dir)
    files = list(summary.get("files", []))
    return {
        "status": "complete",
        "directory": str(evidence_dir),
        "summary": str(evidence_dir / "workflow_validation_summary.json"),
        "files": files,
        "required_artifacts": summary.get("required_artifacts", 0),
        "missing_required_artifacts": summary.get("missing_required_artifacts", 0),
    }


def _finalize_evidence(manifest, run_dir, provenance, events_path):
    """Record evidence success/failure without rerunning workflow stages."""
    run_dir = Path(run_dir).resolve()
    evidence_dir = run_dir / "evidence"
    manifest["evidence"] = {
        "status": "running",
        "directory": str(evidence_dir),
        "summary": str(evidence_dir / "workflow_validation_summary.json"),
        "files": [],
    }
    _save_attempt(manifest, run_dir, provenance)
    try:
        evidence = _generate_evidence(manifest, run_dir)
    except Exception as exc:
        manifest["evidence"] = {
            "status": "failed",
            "directory": str(evidence_dir),
            "summary": str(evidence_dir / "workflow_validation_summary.json"),
            "files": sorted(path.name for path in evidence_dir.glob("*") if path.is_file())
            if evidence_dir.is_dir() else [],
            "error": f"{type(exc).__name__}: {exc}",
        }
        if manifest.get("status") == "complete":
            manifest["status"] = "complete_with_reporting_error"
        _append_event(
            events_path, "evidence_failed", run_id=manifest["run_id"],
            error=manifest["evidence"]["error"],
        )
        _save_attempt(manifest, run_dir, provenance)
        return False
    manifest["evidence"] = evidence
    _append_event(
        events_path, "evidence_completed", run_id=manifest["run_id"],
        files=evidence.get("files", []),
    )
    _save_attempt(manifest, run_dir, provenance)
    return True


def run(args):
    try:
        context = _load_context(args)
    except RunContractError as exc:
        raise SystemExit(f"[run] ERROR: {exc}") from exc

    options = {
        "parallel": args.parallel,
        "target_bed": str(Path(args.target_bed).expanduser().resolve()) if args.target_bed else None,
        "qc_dinucleotide": bool(args.qc_dinucleotide),
    }
    identity = dict(context["identity"])
    identity["options_sha256"] = hashlib.sha256(
        json.dumps(options, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    context["identity"] = identity
    provenance = Path(context["paths"]["provenance"]).resolve()
    previous = _load_previous(provenance, identity)
    run_id = _new_run_id()
    run_dir = provenance / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    events_path = run_dir / "events.jsonl"
    configure_command_log(
        provenance / "commands.jsonl", run_id=run_id,
        mirror_paths=[run_dir / "commands.jsonl"],
    )
    _atomic_write_text(run_dir / "commands.jsonl", "")
    manifest = {
        "run_schema_version": RUN_SCHEMA_VERSION,
        "run_id": run_id,
        "run_dir": str(run_dir),
        "project_root": str(Path(context["paths"]["work_dir"]).resolve()),
        "config": str(context["config_path"]),
        "lock": str(context["lock_path"]),
        "project_identity": identity,
        "options": options,
        "resource_plan": {"status": "NOT_PLANNED"},
        "previous_run_id": previous.get("run_id") if previous else None,
        "status": "running",
        "started_at": _utc_now(),
        "finished_at": None,
        "stages": [],
    }
    not_run = {
        "report_version": 1, "status": "NOT_RUN", "exit_code": 0,
        "summary": {"pass": 0, "warn": 0, "fail": 0}, "checks": [],
        "note": "Preflight was not reached.",
    }
    _atomic_write_json(run_dir / "doctor-before.json", not_run)
    _atomic_write_json(run_dir / "tool-versions.json", {})
    _write_artifact_tables(manifest, run_dir)
    _append_event(events_path, "run_started", run_id=run_id)
    _save_attempt(manifest, run_dir, provenance)

    try:
        plan = _build_stage_plan(context["input_type"], args.qc_dinucleotide)
        manifest["resource_plan"] = _build_resource_plan(context, plan, args)
        for stage in plan:
            specs = _artifact_specs(context, stage, args) if stage["applicable"] else []
            stage_resources = next(
                value for value in manifest["resource_plan"]["stages"]
                if value["stage"] == stage["id"]
            )
            manifest["stages"].append({
                **stage,
                "command": _stage_command(context, stage, args),
                "resources": stage_resources,
                "status": "pending" if stage["applicable"] else "skipped",
                "started_at": None,
                "finished_at": None,
                "expected": specs,
                "outputs": [],
                "figures": [],
                "quarantined": [],
                "error": None,
            })
        _write_artifact_tables(manifest, run_dir)
        _save_attempt(manifest, run_dir, provenance)
    except BaseException as exc:
        manifest["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        manifest["finished_at"] = _utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _append_event(events_path, "run_failed", run_id=run_id, error=manifest["error"])
        _finalize_evidence(manifest, run_dir, provenance, events_path)
        if isinstance(exc, KeyboardInterrupt):
            raise SystemExit(130) from exc
        raise SystemExit(1) from exc

    if args.dry_run:
        preflight = {
            "report_version": 1, "status": "NOT_RUN", "exit_code": 0,
            "summary": {"pass": 0, "warn": 0, "fail": 0}, "checks": [],
            "note": "External tool probes are not executed during dry-run.",
        }
        _atomic_write_json(run_dir / "doctor-before.json", preflight)
        _atomic_write_json(run_dir / "tool-versions.json", {})
        for stage in manifest["stages"]:
            if stage["applicable"]:
                stage["status"] = "planned"
        manifest["status"] = "planned"
        manifest["finished_at"] = _utc_now()
        _append_event(events_path, "run_planned", run_id=run_id)
        if not _finalize_evidence(manifest, run_dir, provenance, events_path):
            disp(f"[run] ERROR: evidence reporting failed; inspect {run_dir / 'run.json'}")
            raise SystemExit(1)
        disp(f"[run] dry-run plan written: {run_dir / 'run-summary.html'}")
        return manifest

    try:
        preflight = _run_preflight(context, args)
        _atomic_write_json(run_dir / "doctor-before.json", preflight)
        _atomic_write_json(run_dir / "tool-versions.json", _tool_versions(preflight))
        if preflight["exit_code"] != 0:
            raise RunContractError(
                f"preflight doctor reported {preflight['summary']['fail']} required failure(s); "
                f"inspect {run_dir / 'doctor-before.json'}"
            )

        for stage_record in manifest["stages"]:
            if not stage_record["applicable"]:
                _append_event(events_path, "stage_skipped", stage=stage_record["id"], reason="input_type")
                continue
            specs = stage_record["expected"]
            if _can_resume(previous, identity, stage_record["id"], specs):
                stage_record["status"] = "resumed"
                stage_record["finished_at"] = _utc_now()
                _record_artifacts(stage_record, specs)
                _append_event(events_path, "stage_resumed", stage=stage_record["id"])
                _save_attempt(manifest, run_dir, provenance)
                continue

            issues = _validate_artifacts(specs)
            existing = [Path(spec["path"]) for spec in specs if Path(spec["path"]).exists()]
            if args.adopt_existing and not issues:
                stage_record["status"] = "adopted"
                stage_record["finished_at"] = _utc_now()
                _record_artifacts(stage_record, specs)
                _append_event(events_path, "stage_adopted", stage=stage_record["id"])
                _save_attempt(manifest, run_dir, provenance)
                continue
            previous_stage = next(
                (item for item in (previous or {}).get("stages", []) if item.get("id") == stage_record["id"]),
                None,
            )
            retryable_previous = (
                previous
                and previous.get("project_identity") == identity
                and previous_stage
                and previous_stage.get("status") in {
                    "running", "failed", "interrupted", "complete", "resumed", "adopted",
                }
            )
            if existing and (retryable_previous or args.adopt_existing):
                quarantine_root = provenance / "quarantine" / run_id
                stage_record["quarantined"] = _quarantine_paths(
                    existing, context["paths"]["results"], quarantine_root
                )
                _append_event(
                    events_path, "stage_outputs_quarantined", stage=stage_record["id"],
                    count=len(stage_record["quarantined"]),
                )
            elif existing:
                raise RunContractError(
                    f"{stage_record['id']} has outputs without a trusted matching run manifest. "
                    "Review them and rerun with --adopt-existing to validate/adopt or quarantine them."
                )

            stage_record["status"] = "running"
            stage_record["started_at"] = _utc_now()
            _append_event(events_path, "stage_started", stage=stage_record["id"], command=stage_record["command"])
            _save_attempt(manifest, run_dir, provenance)
            try:
                _execute_stage(context, stage_record, args)
                issues = _validate_artifacts(specs)
                if issues:
                    raise RunContractError(
                        f"{stage_record['id']} did not satisfy its artifact contract: "
                        + "; ".join(issues)
                    )
            except BaseException as exc:
                stage_record["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
                stage_record["finished_at"] = _utc_now()
                stage_record["error"] = f"{type(exc).__name__}: {exc}"
                _append_event(
                    events_path, "stage_failed", stage=stage_record["id"], error=stage_record["error"]
                )
                raise
            stage_record["status"] = "complete"
            stage_record["finished_at"] = _utc_now()
            _record_artifacts(stage_record, specs)
            _append_event(events_path, "stage_completed", stage=stage_record["id"])
            _save_attempt(manifest, run_dir, provenance)
    except BaseException as exc:
        manifest["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        manifest["finished_at"] = _utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _append_event(events_path, "run_failed", run_id=run_id, error=manifest["error"])
        _finalize_evidence(manifest, run_dir, provenance, events_path)
        if isinstance(exc, KeyboardInterrupt):
            raise SystemExit(130) from exc
        raise SystemExit(1) from exc

    manifest["status"] = "complete"
    manifest["finished_at"] = _utc_now()
    _append_event(events_path, "run_completed", run_id=run_id)
    if not _finalize_evidence(manifest, run_dir, provenance, events_path):
        disp(
            "[run] ERROR: analysis stages completed, but evidence reporting failed; "
            f"inspect {run_dir / 'run.json'}"
        )
        raise SystemExit(2)
    disp(f"[run] complete: {run_dir / 'run-summary.html'}")
    return manifest
