import json
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import analysis_workflow
    import cftk

    return analysis_workflow, cftk


def _context(tmp_path):
    results = tmp_path / "results"
    config = tmp_path / "cftk_init.json"
    lock = tmp_path / "cftk.lock.json"
    config.write_text('{"schema_version": 2}\n')
    lock.write_text("{}\n")
    cfg = {
        "schema_version": 2,
        "project_name": "analysis-test",
        "output_dir": str(tmp_path),
        "comparison": None,
        "group_roles": {"Control": "control"},
        "samples": {"Control": [{"name": "control", "input_type": "bam", "bam": "control.bam"}]},
        "reference_data": {},
        "process": {
            "parallel_samples": 1,
            "step3_markdup": {"params": {"cores": 2}},
            "step4_methylation": {"params": {"cores": 2}},
        },
        "analysis": {
            "diff": {"params": {"modalities": ["cpg"]}},
            "dmr": {"params": {"cores": 2}},
            "mesa": {"params": {"modalities": ["cpg"]}},
            "frag": {},
        },
    }
    return {
        "config_path": config,
        "lock_path": lock,
        "raw": {"schema_version": 2},
        "cfg": cfg,
        "samples": [{"name": "control", "input_type": "bam", "bam": "control.bam", "group": "Control"}],
        "paths": {
            "work_dir": str(tmp_path),
            "results": str(results),
            "provenance": str(results / "provenance"),
            "differential": str(results / "3_differential"),
            "fragmentomics": str(results / "4_fragmentomics"),
            "occ_out": str(results / "4_fragmentomics/occupancy"),
            "wps_out": str(results / "4_fragmentomics/wps"),
            "delfi_out": str(results / "4_fragmentomics/delfi"),
            "end_motif_out": str(results / "4_fragmentomics/end_motif"),
            "cleavage_out": str(results / "4_fragmentomics/cleavage"),
            "mesa": str(results / "5_mesa"),
            "report": str(results / "report"),
            "methylation": str(results / "1_process/4_methylation"),
            "cpg_matrix": str(results / "1_process/5_merged_matrix"),
            "markdup": str(results / "1_process/3_markdup"),
            "alignment": str(results / "1_process/2_alignment"),
        },
        "identity": {"config_sha256": "config", "lock_sha256": "lock"},
    }


def _args(context, **overrides):
    values = {
        "config": str(context["config_path"]),
        "preset": "report",
        "stages": None,
        "parallel": None,
        "dry_run": False,
        "adopt_existing": False,
        "json": False,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def _doctor_pass(*args, **kwargs):
    return {
        "report_version": 1,
        "status": "PASS",
        "exit_code": 0,
        "summary": {"pass": 1, "warn": 0, "fail": 0},
        "checks": [],
    }


def test_presets_are_role_aware_and_do_not_infer_group_meaning(modules):
    analysis_workflow, _ = modules
    one_group = {
        "samples": {"Healthy": [{"name": "h1"}]},
        "group_roles": {"Healthy": "control"},
        "comparison": None,
    }
    two_group = {
        "samples": {"Healthy": [{"name": "h1"}], "Disease": [{"name": "d1"}]},
        "group_roles": {"Healthy": "control", "Disease": "case"},
        "comparison": "Healthy_vs_Disease",
        "analysis": {"diff": {"params": {"modalities": ["cpg", "occupancy", "wps"]}}},
    }

    assert analysis_workflow.resolve_stages("auto", one_group) == (
        "fragmentomics.occupancy", "fragmentomics.wps", "analysis.report",
    )
    assert analysis_workflow.resolve_stages("auto", two_group) == (
        "fragmentomics.occupancy",
        "fragmentomics.wps",
        "analysis.diff",
        "analysis.report",
    )
    assert analysis_workflow.comparison_role_errors(two_group) == []
    assert analysis_workflow.comparison_role_errors(one_group)
    assert "analysis.mesa" in analysis_workflow.resolve_stages("all", two_group)


def test_artifact_contract_handles_single_sample_descriptive_outputs(modules, tmp_path):
    analysis_workflow, _ = modules
    context = _context(tmp_path)

    occupancy = analysis_workflow._artifact_specs(context, "fragmentomics.occupancy")
    paths = {Path(item["path"]).name for item in occupancy}

    assert "control.occupancy.tsv" in paths
    assert "control_occupancy.png" in paths
    assert "occupancy_matrix.tsv" not in paths


def test_planned_fragment_matrix_satisfies_later_analysis_preflight(modules, tmp_path):
    analysis_workflow, _ = modules
    context = _context(tmp_path)
    context["cfg"]["samples"] = {
        "Control": [{"name": "control", "input_type": "bam", "bam": "control.bam"}],
        "Case": [{"name": "case", "input_type": "bam", "bam": "case.bam"}],
    }
    context["cfg"]["group_roles"] = {"Control": "control", "Case": "case"}
    context["cfg"]["comparison"] = "Control_vs_Case"
    context["cfg"]["analysis"]["diff"]["params"]["modalities"] = ["occupancy"]
    context["samples"] = [
        {"name": "control", "input_type": "bam", "bam": "control.bam", "group": "Control"},
        {"name": "case", "input_type": "bam", "bam": "case.bam", "group": "Case"},
    ]

    occupancy_outputs = {
        spec["path"]
        for spec in analysis_workflow._artifact_specs(
            context, "fragmentomics.occupancy"
        )
    }

    assert analysis_workflow._required_inputs(context, "analysis.diff")
    assert not analysis_workflow._required_inputs(
        context,
        "analysis.diff",
        planned_outputs=occupancy_outputs,
    )


def test_dry_run_writes_analysis_manifest_plan_and_evidence(modules, tmp_path, monkeypatch):
    analysis_workflow, _ = modules
    context = _context(tmp_path)
    monkeypatch.setattr(analysis_workflow, "_load_context", lambda args: context)
    import doctor
    monkeypatch.setattr(doctor, "run_doctor", _doctor_pass)

    manifest = analysis_workflow.run(_args(context, dry_run=True))

    assert manifest["status"] == "planned"
    run_dir = Path(manifest["run_dir"])
    assert (run_dir / "run.json").is_file()
    assert (run_dir / "analysis-plan.json").is_file()
    assert (run_dir / "doctor-before.json").is_file()
    assert (run_dir / "evidence/workflow_stage_evidence.png").is_file()
    assert json.loads((run_dir / "run.json").read_text())["workflow"] == "downstream-analysis"


def test_report_stage_resumes_only_after_valid_artifact_contract(modules, tmp_path, monkeypatch):
    analysis_workflow, _ = modules
    context = _context(tmp_path)
    monkeypatch.setattr(analysis_workflow, "_load_context", lambda args: context)
    import doctor
    monkeypatch.setattr(doctor, "run_doctor", _doctor_pass)
    calls = []

    def execute(context, stage_id, args):
        calls.append(stage_id)
        output = Path(context["paths"]["report"]) / "report.html"
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("<html>report</html>\n")

    monkeypatch.setattr(analysis_workflow, "_execute_stage", execute)
    first = analysis_workflow.run(_args(context))
    second = analysis_workflow.run(_args(context))

    assert calls == ["analysis.report"]
    assert first["status"] == "complete"
    assert second["status"] == "complete"
    assert second["stages"][0]["status"] == "resumed"
    assert (Path(second["run_dir"]) / "evidence/workflow_artifact_inventory.tsv").is_file()


def test_adopted_outputs_resume_without_requiring_adoption_again(modules, tmp_path, monkeypatch):
    analysis_workflow, _ = modules
    context = _context(tmp_path)
    monkeypatch.setattr(analysis_workflow, "_load_context", lambda args: context)
    import doctor
    monkeypatch.setattr(doctor, "run_doctor", _doctor_pass)

    output = Path(context["paths"]["report"]) / "report.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("<html>expert report</html>\n")
    monkeypatch.setattr(
        analysis_workflow,
        "_execute_stage",
        lambda *args: pytest.fail("adopted output was re-executed"),
    )

    adopted = analysis_workflow.run(_args(context, adopt_existing=True))
    resumed = analysis_workflow.run(_args(context, adopt_existing=False))

    assert adopted["options"]["adopt_existing"] is True
    assert adopted["stages"][0]["status"] == "adopted"
    assert resumed["stages"][0]["status"] == "resumed"


def test_mesa_requires_explicit_roles_and_fragmentomics_accepts_one_group(modules, tmp_path, monkeypatch):
    _, cftk = modules
    with pytest.raises(SystemExit, match="explicit control/case roles"):
        cftk._make_label(
            {
                "comparison": "Control_vs_Disease",
                "samples": {"Control": [{"name": "control"}], "Disease": [{"name": "case"}]},
            },
            {"mesa": str(tmp_path / "mesa")},
        )

    cfg = {
        "comparison": None,
        "samples": {"Control": [{"name": "control", "input_type": "bam"}]},
        "reference_data": {},
        "process": {"parallel_samples": 1, "step3_markdup": {"params": {"cores": 1}}},
        "analysis": {"frag": {}},
    }
    paths = {
        "occ_out": str(tmp_path / "occupancy"), "wps_out": str(tmp_path / "wps"),
        "delfi_out": str(tmp_path / "delfi"), "end_motif_out": str(tmp_path / "motif"),
        "cleavage_out": str(tmp_path / "cleavage"),
    }
    monkeypatch.setattr(cftk, "_load", lambda args: (cfg, paths))
    import init
    import analysis.occupancy
    import visualization.visualization
    monkeypatch.setattr(init, "get_bam", lambda sample, paths: f"{sample['name']}.bam")
    captured = []
    monkeypatch.setattr(analysis.occupancy, "run_occupancy", lambda args: captured.append(args.group_labels))
    monkeypatch.setattr(visualization.visualization, "plot_fragmentomics", lambda args, mode: None)

    cftk._cmd_frag(SimpleNamespace(
        config="unused.json", parallel=None, occupancy=True, wps=False,
        delfi=False, end_motif=False, cleavage=False,
    ))

    assert captured == [{"Control": ["control"]}]


def test_parser_exposes_planning_and_analysis_commands(modules):
    _, cftk = modules
    parser = cftk.build_parser()

    plan_args = parser.parse_args(["plan", "--preset", "comparative"])
    analyze_args = parser.parse_args(["analyze", "--stage", "diff", "report"])
    doctor_args = parser.parse_args(["doctor", "--analysis-preset", "descriptive"])

    assert plan_args.preset == "comparative"
    assert analyze_args.stages == ["diff", "report"]
    assert doctor_args.analysis_preset == "descriptive"


def test_differential_stage_receives_the_planned_cpu_budget(modules, tmp_path):
    analysis_workflow, _ = modules
    context = _context(tmp_path)

    args = _args(context, preset="differential", parallel=2)
    resource = analysis_workflow._resource_plan(
        context, ("analysis.diff",), args
    )["stages"][0]
    stage_args = analysis_workflow._stage_args(context, "analysis.diff", args)

    assert resource["threads_per_sample"] == 2
    assert stage_args.cores == resource["threads_per_sample"]


def test_dmr_r_package_probe_is_actionable(modules, monkeypatch):
    analysis_workflow, _ = modules
    import doctor

    monkeypatch.setattr(
        analysis_workflow.shutil, "which", lambda name: "/env/bin/Rscript"
    )
    monkeypatch.setattr(
        analysis_workflow.subprocess,
        "run",
        lambda *args, **kwargs: SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="Missing R packages: annotatr",
        ),
    )
    checks = doctor._Checks()

    analysis_workflow._check_dmr_r_packages(checks)

    assert checks.items == [{
        "id": "analysis.dmr.r_packages",
        "status": "FAIL",
        "summary": "DMR R annotation packages are unavailable: Missing R packages: annotatr",
        "remedy": "Install annotatr, TxDb.Hsapiens.UCSC.hg38.knownGene, and GenomicRanges in the active R environment.",
    }]


def test_doctor_runs_analysis_only_checks_for_selected_stage(modules, tmp_path, monkeypatch):
    _, _ = modules
    import doctor

    config = tmp_path / "cftk_init.json"
    config.write_text("{}\n")
    cfg = _context(tmp_path)["cfg"]
    monkeypatch.setattr(doctor, "_load_raw_config", lambda path, checks: {"schema_version": 1})
    monkeypatch.setattr(doctor, "load_config", lambda *args, **kwargs: cfg)
    monkeypatch.setattr(doctor, "_check_profile", lambda checks, raw, path: None)
    monkeypatch.setattr(doctor, "_check_lock", lambda checks, raw, path, profile: None)

    report = doctor.run_doctor(SimpleNamespace(
        config=str(config), step=[], target_bed=None, skip_picard_metrics=True,
        parallel=None, analysis_only=True, analysis_stages=["analysis.report"],
        analysis_preset=None,
    ))

    assert report["exit_code"] == 0
    assert report["steps"] == []
    assert report["analysis_stages"] == ["analysis.report"]
    assert any(item["id"] == "analysis.input" and item["status"] == "PASS" for item in report["checks"])


def test_doctor_checks_fragmentomics_inputs_and_references(modules, tmp_path, monkeypatch):
    _, _ = modules
    import doctor

    config = tmp_path / "cftk_init.json"
    config.write_text("{}\n")
    cfg = _context(tmp_path)["cfg"]
    monkeypatch.setattr(doctor, "_load_raw_config", lambda path, checks: {"schema_version": 1})
    monkeypatch.setattr(doctor, "load_config", lambda *args, **kwargs: cfg)
    monkeypatch.setattr(doctor, "_check_profile", lambda checks, raw, path: None)
    monkeypatch.setattr(doctor, "_check_lock", lambda checks, raw, path, profile: None)

    report = doctor.run_doctor(SimpleNamespace(
        config=str(config), step=[], target_bed=None, skip_picard_metrics=True,
        parallel=None, analysis_only=True, analysis_stages=["fragmentomics.wps"],
        analysis_preset=None,
    ))

    assert report["exit_code"] == 1
    assert any(item["id"] == "analysis.reference" for item in report["checks"])
    assert any(item["id"] == "analysis.input" for item in report["checks"])


def test_plan_resolves_a_real_locked_schema_v2_project(modules, tmp_path, capsys):
    analysis_workflow, _ = modules
    from test_reference_profiles import make_profile
    import init

    reference_root = tmp_path / "references"
    make_profile(reference_root, checksums=True)
    bam = tmp_path / "control.bam"
    bam.touch()
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n"
        "control\tControl\tcontrol\tbam\t\t\tcontrol.bam\n"
    )
    config = tmp_path / "cftk_init.json"
    config.write_text(json.dumps({
        "schema_version": 2,
        "project_name": "real-plan",
        "output_dir": ".",
        "assay": "twist_human_methylome",
        "genome": "hg38",
        "samples": "samples.tsv",
        "reference_mode": "local",
        "reference_root": str(reference_root),
        "reference_profile": {"id": "twist_human_methylome_hg38", "version": "1.0.0"},
    }, indent=2) + "\n")
    init.write_lockfile(config)

    payload = analysis_workflow.plan(SimpleNamespace(
        config=str(config), preset="report", stages=None, parallel=None, json=False,
    ))

    assert payload["status"] == "ready"
    plan_path = Path(capsys.readouterr().out.split("Plan: ", 1)[1].splitlines()[0])
    assert json.loads(plan_path.read_text())["resolved_stages"] == ["analysis.report"]

    manifest = analysis_workflow.run(SimpleNamespace(
        config=str(config), preset="report", stages=None, parallel=None,
        dry_run=True, adopt_existing=False, json=False,
    ))

    assert manifest["status"] == "planned"
    assert (Path(manifest["run_dir"]) / "evidence/workflow_stage_evidence.png").is_file()

    complete = analysis_workflow.run(SimpleNamespace(
        config=str(config), preset="report", stages=None, parallel=None,
        dry_run=False, adopt_existing=False, json=False,
    ))

    assert complete["status"] == "complete"
    assert (tmp_path / "results/report/report.html").is_file()
