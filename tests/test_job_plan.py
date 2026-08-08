import json
from pathlib import Path
from types import SimpleNamespace


def test_job_plan_writes_one_fragmentomics_task_per_sample_and_stage(
    monkeypatch, tmp_path
):
    monkeypatch.syspath_prepend("src")
    import job_plan

    config = tmp_path / "cftk_init.json"
    config.write_text("{}\n")
    monkeypatch.setattr(job_plan, "load_config", lambda _path: {"samples": {}})
    monkeypatch.setattr(
        job_plan,
        "get_all_samples",
        lambda _cfg: [{"name": "control"}, {"name": "case"}],
    )
    monkeypatch.setattr(
        job_plan,
        "get_work_paths",
        lambda _cfg: {"provenance": str(tmp_path / "results/provenance")},
    )

    plan = job_plan.write_fragmentomics_job_plan(SimpleNamespace(
        config=str(config), stages=["delfi", "end_motif"],
        fragmentomics_scope="panel", slurm=True,
    ))

    assert plan["execution_mode"] == "per-sample"
    assert plan["task_count"] == 4
    assert plan["finalizer_count"] == 2
    assert Path(plan["plan_path"]).is_file()
    assert json.loads(Path(plan["plan_path"]).read_text())["plan_path"] == plan["plan_path"]
    assert Path(plan["slurm_submit_script"]).is_file()
    assert "--array=0-1" in Path(plan["slurm_submit_script"]).read_text()
    task = next(item for item in plan["tasks"] if item["sample"] == "control")
    script = Path(task["script"])
    assert script.is_file()
    assert "--sample control --parallel 1 --no-finalize" in script.read_text()
    finalizer = Path(plan["finalizers"][0]["script"]).read_text()
    assert "frag --delfi --finalize" in finalizer
    assert "analyze --stage delfi --adopt-existing" in finalizer


def test_frag_parser_accepts_generated_parallel_option(monkeypatch):
    monkeypatch.syspath_prepend("src")
    import cftk

    args = cftk.build_parser().parse_args([
        "frag", "--end-motif", "--sample", "control", "--parallel", "1",
        "--no-finalize",
    ])

    assert args.samples == ["control"]
    assert args.parallel == 1
    assert args.no_finalize is True


def test_fragmentomics_sample_selection_runs_only_the_named_sample(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("src")
    import cftk
    import init
    from analysis import delfi
    from visualization import visualization

    cfg = {
        "samples": {
            "Control": [{"name": "control", "input_type": "bam"}],
            "Case": [{"name": "case", "input_type": "bam"}],
        },
        "process": {"parallel_samples": 2, "step3_markdup": {"params": {"cores": 4}}},
        "reference_data": {},
        "analysis": {"frag": {}},
    }
    paths = {
        "occ_out": str(tmp_path / "occupancy"),
        "wps_out": str(tmp_path / "wps"),
        "delfi_out": str(tmp_path / "delfi"),
        "end_motif_out": str(tmp_path / "end_motif"),
        "cleavage_out": str(tmp_path / "cleavage"),
        "fragmentomics": str(tmp_path / "fragmentomics"),
    }
    monkeypatch.setattr(cftk, "_load", lambda _args: (cfg, paths))
    monkeypatch.setattr(
        init, "get_bam", lambda sample, _paths: f"{sample['name']}.markdup.bam"
    )
    calls = []
    monkeypatch.setattr(delfi, "run_delfi", lambda args: calls.append(list(args.infile)))
    monkeypatch.setattr(visualization, "plot_fragmentomics", lambda *args, **kwargs: None)

    cftk._cmd_frag(SimpleNamespace(
        config="unused.json", parallel=None, samples=["case"],
        no_finalize=True, finalize=False, occupancy=False, wps=False,
        delfi=True, end_motif=False, cleavage=False, fragmentomics_scope=None,
    ))

    assert calls == [["case.markdup.bam"]]


def test_single_sample_fragmentomics_uses_scheduler_task_budget(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("src")
    import cftk
    import init
    from analysis import delfi
    from visualization import visualization
    import resource_planning

    cfg = {
        "reference_data": {},
        "samples": {"Control": [{"name": "control", "input_type": "bam"}]},
        "process": {"parallel_samples": 4, "step3_markdup": {"params": {"cores": 20}}},
        "analysis": {"frag": {}},
    }
    paths = {key: str(tmp_path / key) for key in (
        "occ_out", "wps_out", "delfi_out", "end_motif_out", "cleavage_out",
    )}
    paths["fragmentomics"] = str(tmp_path / "fragmentomics")
    monkeypatch.setattr(cftk, "_load", lambda _args: (cfg, paths))
    monkeypatch.setattr(init, "get_bam", lambda sample, _paths: f"{sample['name']}.bam")
    monkeypatch.setattr(
        resource_planning,
        "detect_scheduler_allocation",
        lambda: {"allocated_cores": 7, "valid": True, "scheduler": "slurm"},
    )
    observed = []
    monkeypatch.setattr(delfi, "run_delfi", lambda args: observed.append((args.cores, args.parallel)))
    monkeypatch.setattr(visualization, "plot_fragmentomics", lambda *args, **kwargs: None)

    cftk._cmd_frag(SimpleNamespace(
        config="unused.json", parallel=1, samples=["control"],
        no_finalize=True, finalize=False, occupancy=False, wps=False,
        delfi=True, end_motif=False, cleavage=False, fragmentomics_scope=None,
    ))

    assert observed == [(7, 1)]


def test_end_motif_passes_planned_workers_to_finaletoolkit(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("src")
    from analysis import end_motif

    commands = []
    monkeypatch.setattr(end_motif, "_run", lambda command, label: commands.append(command))
    args = SimpleNamespace(
        infile=[str(tmp_path / "sample.markdup.bam")],
        end_motif_out=str(tmp_path / "end_motif"),
        kmer=4, mapq=30, min_frag=100, max_frag=220, em_extra="",
        cores=7, parallel=1, genome2bit="hg38.2bit",
    )

    end_motif.run_end_motif(args)

    assert len(commands) == 1
    assert "-q 30 -w 7" in commands[0]
