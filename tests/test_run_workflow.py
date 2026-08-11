import json
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend("src")
    for variable in (
        "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "PBS_NP", "NSLOTS",
    ):
        monkeypatch.delenv(variable, raising=False)
    import cftk
    import run_workflow

    return cftk, run_workflow


def _default_cfg(samples):
    return {
        "samples": {
            "Control": [s for s in samples if s["name"].startswith("control")],
            "Case": [s for s in samples if s["name"].startswith("case")],
        },
        "process": {
            "step1_trimming": {"tool": "trim_galore", "params": {}},
            "step2_alignment": {"tool": "bwameth", "params": {}},
            "step3_markdup": {"tool": "sambamba", "params": {}},
            "step4_methylation": {"tool": "methyldackel", "params": {}},
        },
    }


def _args(config="cftk_init.json", **overrides):
    values = {
        "config": config,
        "parallel": None,
        "target_bed": None,
        "dry_run": False,
        "adopt_existing": False,
        "qc_dinucleotide": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def _context(tmp_path, run_workflow, input_type="fastq"):
    config = tmp_path / "cftk_init.json"
    config.write_text('{"schema_version": 2}\n')
    lock = tmp_path / "cftk.lock.json"
    lock.write_text('{"lock_version": 1}\n')
    samples = [
        {"name": "control_1", "input_type": input_type},
        {"name": "case_1", "input_type": input_type},
    ]
    results = tmp_path / "results"
    paths = {
        "work_dir": str(tmp_path),
        "results": str(results),
        "provenance": str(results / "provenance"),
    }
    return {
        "config_path": config,
        "lock_path": lock,
        "raw": {"schema_version": 2},
        "cfg": _default_cfg(samples),
        "samples": samples,
        "paths": paths,
        "identity": {
            "config_sha256": "config-hash",
            "lock_sha256": "lock-hash",
        },
        "input_type": input_type,
        "run_workflow": run_workflow,
    }


def test_parser_registers_zero_argument_beginner_run(modules):
    cftk, _ = modules

    args = cftk.build_parser().parse_args(["run"])

    assert args.parallel is None
    assert args.target_bed is None
    assert not args.dry_run
    assert not args.adopt_existing
    assert not args.qc_dinucleotide
    assert args.downstream is None
    assert args.fragmentomics_scope is None


def test_beginner_run_accepts_explicit_downstream_preset(modules):
    cftk, _ = modules

    args = cftk.build_parser().parse_args(["run", "--downstream", "fragmentomics"])

    assert args.downstream == "fragmentomics"

    scoped = cftk.build_parser().parse_args([
        "run", "--downstream", "auto", "--fragmentomics-scope", "panel"
    ])
    assert scoped.fragmentomics_scope == "panel"

    differential = cftk.build_parser().parse_args([
        "run", "--downstream", "differential", "--modality", "cpg", "wps",
    ])
    assert differential.differential_modalities == ["cpg", "wps"]


def test_linked_downstream_manifest_is_rendered_in_core_summary(modules, tmp_path):
    _, run_workflow = modules
    manifest = {
        "run_id": "core-1",
        "status": "complete",
        "started_at": "now",
        "finished_at": "now",
        "stages": [],
        "downstream": {
            "preset": "auto",
            "status": "complete",
            "run_id": "analysis-1",
            "manifest": str(tmp_path / "analysis-runs/analysis-1/run.json"),
            "summary": str(tmp_path / "analysis-runs/analysis-1/run-summary.html"),
        },
    }
    output = tmp_path / "runs/core-1/run-summary.html"

    run_workflow._write_summary_html(manifest, output, tmp_path)

    rendered = output.read_text(encoding="utf-8")
    assert "Downstream workflow" in rendered
    assert "analysis-runs/analysis-1/run-summary.html" in rendered


def test_run_downstream_calls_existing_analysis_runner(modules, monkeypatch, tmp_path):
    cftk, run_workflow = modules
    core_dir = tmp_path / "results/provenance/runs/core-1"
    core_dir.mkdir(parents=True)
    core = {"run_id": "core-1", "run_dir": str(core_dir), "status": "complete"}
    analysis = {
        "run_id": "analysis-1",
        "run_dir": str(tmp_path / "results/provenance/analysis-runs/analysis-1"),
        "status": "complete",
    }
    seen = {}
    monkeypatch.setattr(run_workflow, "run", lambda args: core)
    monkeypatch.setattr(run_workflow, "_save_attempt", lambda *args: None)
    monkeypatch.setattr(run_workflow, "_append_event", lambda *args, **kwargs: None)

    import analysis_workflow
    def fake_analysis(args):
        seen["args"] = args
        return analysis

    monkeypatch.setattr(analysis_workflow, "run", fake_analysis)

    args = SimpleNamespace(
        config=str(tmp_path / "cftk_init.json"),
        parallel=3,
        target_bed=None,
        dry_run=True,
        adopt_existing=False,
        qc_dinucleotide=False,
        downstream="auto",
        differential_modalities=["cpg"],
    )
    result = cftk._cmd_run(args)

    assert result["downstream"]["preset"] == "auto"
    assert result["downstream"]["run_id"] == "analysis-1"
    assert seen["args"].preset == "auto"
    assert seen["args"].dry_run is True
    assert seen["args"].differential_modalities == ["cpg"]


def test_run_downstream_failure_keeps_core_manifest_link(modules, monkeypatch, tmp_path):
    cftk, run_workflow = modules
    core_dir = tmp_path / "results/provenance/runs/core-1"
    core_dir.mkdir(parents=True)
    core = {"run_id": "core-1", "run_dir": str(core_dir), "status": "complete"}
    monkeypatch.setattr(run_workflow, "run", lambda args: core)
    monkeypatch.setattr(run_workflow, "_save_attempt", lambda *args: None)
    monkeypatch.setattr(run_workflow, "_append_event", lambda *args, **kwargs: None)

    import analysis_workflow
    monkeypatch.setattr(analysis_workflow, "run", lambda args: (_ for _ in ()).throw(SystemExit(1)))

    args = SimpleNamespace(
        config=str(tmp_path / "cftk_init.json"),
        parallel=None,
        target_bed=None,
        dry_run=False,
        adopt_existing=False,
        qc_dinucleotide=False,
        downstream="auto",
        fragmentomics_scope=None,
    )
    with pytest.raises(SystemExit):
        cftk._cmd_run(args)

    assert core["status"] == "complete"
    assert core["downstream"]["status"] == "failed"


def test_beginner_contract_requires_schema_v2_lock_and_default_tools(
    modules, tmp_path
):
    _, run_workflow = modules
    samples = [
        {"name": "control_1", "input_type": "fastq"},
        {"name": "case_1", "input_type": "fastq"},
    ]
    cfg = _default_cfg(samples)
    lock = tmp_path / "cftk.lock.json"
    lock.write_text("{}\n")

    assert run_workflow._validate_project_contract(
        {"schema_version": 2}, cfg, samples, lock
    ) == "fastq"

    with pytest.raises(run_workflow.RunContractError, match="schema-v2"):
        run_workflow._validate_project_contract({}, cfg, samples, lock)

    lock.unlink()
    with pytest.raises(run_workflow.RunContractError, match="cftk.lock.json"):
        run_workflow._validate_project_contract(
            {"schema_version": 2}, cfg, samples, lock
        )

    lock.write_text("{}\n")
    cfg["process"]["step2_alignment"]["tool"] = "bismark"
    with pytest.raises(run_workflow.RunContractError, match="validated toolchain"):
        run_workflow._validate_project_contract(
            {"schema_version": 2}, cfg, samples, lock
        )


def test_beginner_contract_rejects_mixed_input_types(modules, tmp_path):
    _, run_workflow = modules
    samples = [
        {"name": "control_1", "input_type": "fastq"},
        {"name": "case_1", "input_type": "bam"},
    ]
    lock = tmp_path / "cftk.lock.json"
    lock.write_text("{}\n")

    with pytest.raises(run_workflow.RunContractError, match="mix FASTQ and BAM"):
        run_workflow._validate_project_contract(
            {"schema_version": 2}, _default_cfg(samples), samples, lock
        )


def test_stage_plan_has_rigorous_default_order_and_optional_dinucleotide(modules):
    _, run_workflow = modules

    fastq = run_workflow._build_stage_plan("fastq", include_dinucleotide=False)
    assert [stage["id"] for stage in fastq] == [
        "process.1", "process.2", "process.3", "process.4",
        "qc.2", "qc.0", "qc.1",
    ]
    assert all(stage["applicable"] for stage in fastq)

    bam = run_workflow._build_stage_plan("bam", include_dinucleotide=True)
    assert not bam[0]["applicable"]
    assert not bam[1]["applicable"]
    assert [stage["id"] for stage in bam][-1] == "qc.3"


def test_required_artifacts_must_exist_and_be_nonempty(modules, tmp_path):
    _, run_workflow = modules
    output = tmp_path / "output.tsv"
    specs = [{
        "path": str(output),
        "role": "output",
        "required": True,
        "nonempty": True,
        "description": "test output",
    }]

    assert "missing" in run_workflow._validate_artifacts(specs)[0]
    output.touch()
    assert "empty" in run_workflow._validate_artifacts(specs)[0]
    output.write_text("value\n")
    assert run_workflow._validate_artifacts(specs) == []


def test_manifest_resume_requires_matching_identity_and_valid_artifacts(
    modules, tmp_path
):
    _, run_workflow = modules
    output = tmp_path / "done.tsv"
    output.write_text("done\n")
    identity = {"config_sha256": "a", "lock_sha256": "b"}
    software = {"name": "cftk", "identity_sha256": "software-one"}
    previous = {
        "project_identity": identity,
        "software_identity": software,
        "stages": [{"id": "process.1", "status": "complete"}],
    }
    specs = [{
        "path": str(output), "role": "output", "required": True,
        "nonempty": True, "description": "done",
    }]

    assert run_workflow._can_resume(previous, identity, software, "process.1", specs)
    assert not run_workflow._can_resume(
        previous, {**identity, "config_sha256": "changed"}, software, "process.1", specs
    )
    assert not run_workflow._can_resume(
        previous, identity, {**software, "identity_sha256": "software-two"},
        "process.1", specs
    )
    output.unlink()
    assert not run_workflow._can_resume(previous, identity, software, "process.1", specs)


def test_matching_manifest_is_found_when_latest_attempt_is_only_a_dry_run(
    modules, tmp_path
):
    _, run_workflow = modules
    provenance = tmp_path / "provenance"
    identity = {"config_sha256": "a", "lock_sha256": "b"}
    complete_path = provenance / "runs/20260101T000000Z-a/run.json"
    planned_path = provenance / "runs/20260102T000000Z-b/run.json"
    complete_path.parent.mkdir(parents=True)
    planned_path.parent.mkdir(parents=True)
    complete_path.write_text(json.dumps({
        "run_id": "complete", "project_identity": identity,
        "software_identity": {"identity_sha256": "software-one"},
        "status": "complete",
    }))
    planned_path.write_text(json.dumps({
        "run_id": "planned", "project_identity": {**identity, "options_sha256": "dry"},
        "software_identity": {"identity_sha256": "software-one"},
        "status": "planned",
    }))
    (provenance / "latest-run.json").write_text(json.dumps({
        "manifest": str(planned_path),
    }))

    assert run_workflow._load_previous(
        provenance, identity, {"identity_sha256": "software-one"}
    )["run_id"] == "complete"


def test_quarantine_preserves_relative_paths(modules, tmp_path):
    _, run_workflow = modules
    results = tmp_path / "results"
    partial = results / "1_process/2_alignment/sample.bam"
    partial.parent.mkdir(parents=True)
    partial.write_text("partial\n")

    moved = run_workflow._quarantine_paths(
        [partial], results, results / "provenance/quarantine/run-1"
    )

    expected = (
        results / "provenance/quarantine/run-1/1_process/2_alignment/sample.bam"
    )
    assert not partial.exists()
    assert expected.read_text() == "partial\n"
    assert moved == [{"source": str(partial), "destination": str(expected)}]


def test_dry_run_writes_immutable_plan_without_executing_stages(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_run_preflight", lambda context, args: {
            "report_version": 1,
            "status": "NOT_RUN",
            "exit_code": 0,
            "summary": {"pass": 0, "warn": 0, "fail": 0},
            "checks": [],
        },
    )
    monkeypatch.setattr(
        run_workflow, "_execute_stage",
        lambda *args, **kwargs: pytest.fail("dry-run executed a stage"),
    )
    monkeypatch.setattr(run_workflow, "_artifact_specs", lambda *args: [])

    manifest = run_workflow.run(_args(str(context["config_path"]), dry_run=True))

    assert manifest["status"] == "planned"
    assert all(stage["status"] == "planned" for stage in manifest["stages"])
    run_dir = Path(manifest["run_dir"])
    assert json.loads((run_dir / "run.json").read_text())["run_id"] == manifest["run_id"]
    assert (run_dir / "resource-plan.json").is_file()
    assert (run_dir / "expected-outputs.tsv").is_file()
    assert (run_dir / "figures.tsv").is_file()
    assert (run_dir / "run-summary.html").is_file()
    assert manifest["evidence"]["status"] == "complete"
    evidence = run_dir / "evidence"
    for filename in (
        "workflow_artifact_inventory.tsv",
        "workflow_stage_evidence.tsv",
        "workflow_command_evidence.tsv",
        "workflow_stage_evidence.png",
        "workflow_resource_plan.png",
        "workflow_validation_summary.json",
    ):
        assert (evidence / filename).is_file()
    resources = json.loads((run_dir / "resource-plan.json").read_text())
    assert resources["maximum_total_core_budget"] == 20
    assert resources["stages"][0]["threads_per_sample"] == 20
    summary_html = (run_dir / "run-summary.html").read_text()
    assert "CPU resource plan" in summary_html
    assert "Generated evidence" in summary_html
    assert "evidence/workflow_stage_evidence.png" in summary_html


def test_stage_failure_stops_downstream_and_records_terminal_status(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    calls = []
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_run_preflight", lambda context, args: {
            "report_version": 1,
            "status": "PASS",
            "exit_code": 0,
            "summary": {"pass": 1, "warn": 0, "fail": 0},
            "checks": [],
        },
    )

    def artifacts(context, stage, args):
        path = tmp_path / f"{stage['id']}.txt"
        return [{
            "path": str(path), "role": "output", "required": True,
            "nonempty": True, "description": stage["id"],
        }]

    def execute(context, stage, args):
        calls.append(stage["id"])
        if stage["id"] == "process.2":
            raise RuntimeError("alignment failed")
        (tmp_path / f"{stage['id']}.txt").write_text("done\n")

    monkeypatch.setattr(run_workflow, "_artifact_specs", artifacts)
    monkeypatch.setattr(run_workflow, "_execute_stage", execute)

    with pytest.raises(SystemExit) as exc:
        run_workflow.run(_args(str(context["config_path"])))

    assert exc.value.code == 1
    assert calls == ["process.1", "process.2"]
    latest = json.loads(
        (Path(context["paths"]["provenance"]) / "latest-run.json").read_text()
    )
    manifest = json.loads(Path(latest["manifest"]).read_text())
    assert manifest["status"] == "failed"
    assert manifest["stages"][1]["status"] == "failed"
    assert manifest["stages"][2]["status"] == "pending"
    assert manifest["evidence"]["status"] == "complete"
    assert (Path(manifest["run_dir"]) / "evidence/workflow_stage_evidence.png").is_file()


def test_planning_failure_still_records_a_terminal_attempt(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_artifact_specs",
        lambda *args: (_ for _ in ()).throw(FileNotFoundError("missing target BED")),
    )

    with pytest.raises(SystemExit) as exc:
        run_workflow.run(_args(str(context["config_path"]), dry_run=True))

    assert exc.value.code == 1
    latest = json.loads(
        (Path(context["paths"]["provenance"]) / "latest-run.json").read_text()
    )
    manifest = json.loads(Path(latest["manifest"]).read_text())
    assert manifest["status"] == "failed"
    assert manifest["stages"] == []
    assert "missing target BED" in manifest["error"]
    run_dir = Path(manifest["run_dir"])
    assert (run_dir / "doctor-before.json").is_file()
    assert (run_dir / "tool-versions.json").is_file()
    assert (run_dir / "commands.jsonl").is_file()
    assert (run_dir / "expected-outputs.tsv").is_file()
    assert (run_dir / "figures.tsv").is_file()
    summary = (run_dir / "run-summary.html").read_text()
    assert "Terminal error" in summary
    assert "missing target BED" in summary
    assert "Machine-readable run manifest" in summary
    assert manifest["evidence"]["status"] == "complete"
    assert (run_dir / "evidence/workflow_stage_evidence.png").is_file()


def test_interruption_is_recorded_and_returns_shell_interrupt_status(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(run_workflow, "_run_preflight", lambda *args: {
        "report_version": 1, "status": "PASS", "exit_code": 0,
        "summary": {"pass": 1, "warn": 0, "fail": 0}, "checks": [],
    })
    monkeypatch.setattr(run_workflow, "_artifact_specs", lambda *args: [])
    monkeypatch.setattr(
        run_workflow, "_execute_stage",
        lambda *args: (_ for _ in ()).throw(KeyboardInterrupt()),
    )

    with pytest.raises(SystemExit) as exc:
        run_workflow.run(_args(str(context["config_path"])))

    assert exc.value.code == 130
    latest = json.loads(
        (Path(context["paths"]["provenance"]) / "latest-run.json").read_text()
    )
    manifest = json.loads(Path(latest["manifest"]).read_text())
    assert manifest["status"] == "interrupted"
    assert manifest["stages"][0]["status"] == "interrupted"
    assert manifest["evidence"]["status"] == "complete"


def test_complete_run_is_automatically_resumed_from_valid_manifest(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    calls = []
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_run_preflight", lambda context, args: {
            "report_version": 1,
            "status": "PASS",
            "exit_code": 0,
            "summary": {"pass": 1, "warn": 0, "fail": 0},
            "checks": [],
        },
    )

    def artifacts(context, stage, args):
        return [{
            "path": str(tmp_path / f"{stage['id']}.txt"),
            "role": "output", "required": True, "nonempty": True,
            "description": stage["id"],
        }]

    def execute(context, stage, args):
        calls.append(stage["id"])
        (tmp_path / f"{stage['id']}.txt").write_text("done\n")

    monkeypatch.setattr(run_workflow, "_artifact_specs", artifacts)
    monkeypatch.setattr(run_workflow, "_execute_stage", execute)

    first = run_workflow.run(_args(str(context["config_path"])))
    assert first["status"] == "complete"
    assert [stage["status"] for stage in first["stages"]] == ["complete"] * 7
    assert first["evidence"]["status"] == "complete"

    calls.clear()
    second = run_workflow.run(_args(str(context["config_path"])))
    assert second["status"] == "complete"
    assert calls == []
    assert [stage["status"] for stage in second["stages"]] == ["resumed"] * 7
    assert second["previous_run_id"] == first["run_id"]
    assert second["evidence"]["status"] == "complete"


def test_reporting_failure_is_distinct_and_resume_rebuilds_only_evidence(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    calls = []
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(run_workflow, "_run_preflight", lambda *args: {
        "report_version": 1, "status": "PASS", "exit_code": 0,
        "summary": {"pass": 1, "warn": 0, "fail": 0}, "checks": [],
    })

    def artifacts(context, stage, args):
        return [{
            "path": str(tmp_path / f"{stage['id']}.txt"),
            "role": "output", "required": True, "nonempty": True,
            "description": stage["id"],
        }]

    def execute(context, stage, args):
        calls.append(stage["id"])
        Path(artifacts(context, stage, args)[0]["path"]).write_text("done\n")

    monkeypatch.setattr(run_workflow, "_artifact_specs", artifacts)
    monkeypatch.setattr(run_workflow, "_execute_stage", execute)
    original_generate = run_workflow._generate_evidence
    monkeypatch.setattr(
        run_workflow, "_generate_evidence",
        lambda *args: (_ for _ in ()).throw(RuntimeError("report renderer failed")),
    )

    with pytest.raises(SystemExit) as exc:
        run_workflow.run(_args(str(context["config_path"])))

    assert exc.value.code == 2
    assert calls == [stage["id"] for stage in run_workflow._build_stage_plan("fastq")]
    latest = json.loads(
        (Path(context["paths"]["provenance"]) / "latest-run.json").read_text()
    )
    failed_report = json.loads(Path(latest["manifest"]).read_text())
    assert failed_report["status"] == "complete_with_reporting_error"
    assert failed_report["evidence"]["status"] == "failed"
    assert "report renderer failed" in failed_report["evidence"]["error"]

    monkeypatch.setattr(run_workflow, "_generate_evidence", original_generate)
    calls.clear()
    resumed = run_workflow.run(_args(str(context["config_path"])))

    assert resumed["status"] == "complete"
    assert resumed["evidence"]["status"] == "complete"
    assert calls == []
    assert [stage["status"] for stage in resumed["stages"]] == ["resumed"] * 7


def test_damaged_trusted_checkpoint_is_quarantined_and_rebuilt(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    calls = []
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(run_workflow, "_run_preflight", lambda *args: {
        "report_version": 1, "status": "PASS", "exit_code": 0,
        "summary": {"pass": 1, "warn": 0, "fail": 0}, "checks": [],
    })

    def artifacts(context, stage, args):
        return [{
            "path": str(Path(context["paths"]["results"]) / f"{stage['id']}.txt"),
            "role": "output", "required": True, "nonempty": True,
            "description": stage["id"],
        }]

    def execute(context, stage, args):
        calls.append(stage["id"])
        Path(artifacts(context, stage, args)[0]["path"]).parent.mkdir(parents=True, exist_ok=True)
        Path(artifacts(context, stage, args)[0]["path"]).write_text("done\n")

    monkeypatch.setattr(run_workflow, "_artifact_specs", artifacts)
    monkeypatch.setattr(run_workflow, "_execute_stage", execute)

    first = run_workflow.run(_args(str(context["config_path"])))
    damaged = Path(artifacts(context, first["stages"][0], _args())[0]["path"])
    damaged.write_text("")
    calls.clear()

    second = run_workflow.run(_args(str(context["config_path"])))

    assert calls == ["process.1"]
    assert second["stages"][0]["status"] == "complete"
    assert [stage["status"] for stage in second["stages"][1:]] == ["resumed"] * 6
    moved = second["stages"][0]["quarantined"]
    assert len(moved) == 1
    assert Path(moved[0]["destination"]).is_file()


def test_preflight_failure_blocks_every_stage(modules, tmp_path, monkeypatch):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_run_preflight", lambda context, args: {
            "report_version": 1,
            "status": "FAIL",
            "exit_code": 1,
            "summary": {"pass": 0, "warn": 0, "fail": 1},
            "checks": [{
                "id": "tool.missing", "status": "FAIL",
                "summary": "required tool is missing",
            }],
        },
    )
    monkeypatch.setattr(run_workflow, "_artifact_specs", lambda *args: [])
    monkeypatch.setattr(
        run_workflow, "_execute_stage",
        lambda *args: pytest.fail("preflight failure started a stage"),
    )

    with pytest.raises(SystemExit) as exc:
        run_workflow.run(_args(str(context["config_path"])))

    assert exc.value.code == 1
    latest = json.loads(
        (Path(context["paths"]["provenance"]) / "latest-run.json").read_text()
    )
    manifest = json.loads(Path(latest["manifest"]).read_text())
    assert manifest["status"] == "failed"
    assert all(stage["status"] == "pending" for stage in manifest["stages"])


def test_complete_legacy_outputs_require_explicit_adoption(
    modules, tmp_path, monkeypatch
):
    _, run_workflow = modules
    context = _context(tmp_path, run_workflow)
    monkeypatch.setattr(run_workflow, "_load_context", lambda args: context)
    monkeypatch.setattr(
        run_workflow, "_run_preflight", lambda context, args: {
            "report_version": 1, "status": "PASS", "exit_code": 0,
            "summary": {"pass": 1, "warn": 0, "fail": 0}, "checks": [],
        },
    )

    def artifacts(context, stage, args):
        path = tmp_path / f"{stage['id']}.txt"
        path.write_text("legacy\n") if not path.exists() else None
        return [{
            "path": str(path), "role": "output", "required": True,
            "nonempty": True, "description": stage["id"],
        }]

    monkeypatch.setattr(run_workflow, "_artifact_specs", artifacts)
    monkeypatch.setattr(
        run_workflow, "_execute_stage",
        lambda *args: pytest.fail("complete adopted outputs were re-executed"),
    )

    with pytest.raises(SystemExit):
        run_workflow.run(_args(str(context["config_path"])))

    adopted = run_workflow.run(
        _args(str(context["config_path"]), adopt_existing=True)
    )
    assert adopted["status"] == "complete"
    assert [stage["status"] for stage in adopted["stages"]] == ["adopted"] * 7


def test_summary_html_links_recorded_outputs_and_figures(modules, tmp_path):
    _, run_workflow = modules
    output = tmp_path / "results/output.tsv"
    figure = tmp_path / "results/figure.png"
    evidence_dir = tmp_path / "evidence"
    evidence_dir.mkdir()
    (evidence_dir / "workflow_stage_evidence.png").write_bytes(b"png")
    manifest = {
        "run_id": "run-1",
        "status": "complete",
        "started_at": "2026-08-04T00:00:00+00:00",
        "finished_at": "2026-08-04T00:01:00+00:00",
        "evidence": {
            "status": "complete", "directory": str(evidence_dir),
            "files": ["workflow_stage_evidence.png"],
        },
        "stages": [{
            "id": "qc.1", "name": "Methylation QC", "status": "complete",
            "outputs": [{"path": str(output), "description": "QC table"}],
            "figures": [{"path": str(figure), "description": "QC figure"}],
        }],
    }
    html_path = tmp_path / "run-summary.html"

    run_workflow._write_summary_html(manifest, html_path, tmp_path)

    html = html_path.read_text()
    assert "run-1" in html
    assert "Methylation QC" in html
    assert "results/output.tsv" in html
    assert "results/figure.png" in html
    assert '<img loading="lazy" src="results/figure.png"' in html
    assert "Generated evidence" in html
    assert 'src="evidence/workflow_stage_evidence.png"' in html


def test_summary_links_are_relative_to_nested_attempt_directory(modules, tmp_path):
    _, run_workflow = modules
    output = tmp_path / "results/2_qc/qc_summary.tsv"
    manifest = {
        "run_id": "run-1", "status": "complete", "started_at": "now",
        "finished_at": "now", "stages": [{
            "id": "qc.0", "name": "QC", "status": "complete",
            "outputs": [{"path": str(output), "description": "QC summary"}],
            "figures": [],
        }],
    }
    html_path = tmp_path / "results/provenance/runs/run-1/run-summary.html"

    run_workflow._write_summary_html(manifest, html_path, tmp_path)

    assert 'href="../../../2_qc/qc_summary.tsv"' in html_path.read_text()


def test_summary_html_expands_fragmentomics_scope_details(modules, tmp_path):
    _, run_workflow = modules
    scope_path = tmp_path / "results/4_fragmentomics/wps/fragmentomics_scope.json"
    scope_path.parent.mkdir(parents=True)
    scope_path.write_text(
        json.dumps({
            "schema_version": 1,
            "stage": "wps",
            "resolved_scope": {
                "mode": "panel",
                "target_bed": "/refs/twist_targets.bed",
                "target_sha256": "abc123",
                "region_count": 12,
                "note": "panel-only WPS",
            },
        }) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "run_id": "run-scope",
        "status": "complete",
        "started_at": "now",
        "finished_at": "now",
        "fragmentomics_scope": {"mode": "panel", "note": "panel-only WPS"},
        "stages": [{
            "id": "fragmentomics.wps", "name": "WPS", "status": "complete",
            "expected": [{"path": str(scope_path), "description": "scope metadata", "role": "metadata"}],
            "outputs": [], "figures": [],
        }],
    }
    html_path = tmp_path / "run-summary.html"

    run_workflow._write_summary_html(manifest, html_path, tmp_path)

    html = html_path.read_text(encoding="utf-8")
    assert "Target BED" in html
    assert "/refs/twist_targets.bed" in html
    assert "abc123" in html
    assert "12" in html
    assert "Interpretation" in html


def test_environment_and_package_register_beginner_run_dependencies():
    root = Path(__file__).resolve().parents[1]

    assert "deeptools=3.5.5" in (root / "environment.yml").read_text()
    pyproject = (root / "pyproject.toml").read_text()
    required_dependencies = pyproject.split("[project.optional-dependencies]", 1)[0]
    assert '"run_workflow"' in pyproject
    assert '"validation_reports"' in pyproject
    assert (root / "src/validation_reports.py").is_file()
    for dependency in (
        "numpy", "pandas", "scipy", "matplotlib", "seaborn", "joblib", "pysam",
    ):
        assert f'"{dependency}' in required_dependencies
    for dependency in (
        "xgboost", "scikit-learn", "adjusttext", "statsmodels",
        "pybigwig", "bx-python", "finaletoolkit", "mesa-cfdna",
    ):
        assert dependency not in required_dependencies.lower()
    assert "analysis = [" in pyproject
    assert '"mesa-cfdna>=0.7.1"' in pyproject
    assert "fragmentomics = [" in pyproject
    assert "web = [" in pyproject


def test_missing_optional_dependency_reports_install_extra(monkeypatch):
    monkeypatch.syspath_prepend("src")
    import cftk

    def missing_sklearn(_args):
        exc = ModuleNotFoundError("No module named 'sklearn'")
        exc.name = "sklearn"
        raise exc

    args = SimpleNamespace(command="diff", func=missing_sklearn)
    parser = SimpleNamespace(
        parse_args=lambda: args,
        error=lambda message: (_ for _ in ()).throw(SystemExit(message)),
    )
    monkeypatch.setattr(cftk, "build_parser", lambda: parser)

    with pytest.raises(SystemExit, match=r"cftk\[analysis\]"):
        cftk.main()
