from types import SimpleNamespace

import pytest


def test_parallel_plan_divides_one_total_budget(monkeypatch):
    monkeypatch.syspath_prepend("src")
    from resource_planning import plan_parallelism

    assert plan_parallelism(20, 2, 40) == {
        "total_core_budget": 20,
        "requested_parallel_samples": 2,
        "concurrent_samples": 2,
        "threads_per_sample": 10,
        "estimated_peak_threads": 20,
    }
    assert plan_parallelism(20, 5, 2)["threads_per_sample"] == 10


def test_parallel_plan_rejects_more_samples_than_cores(monkeypatch):
    monkeypatch.syspath_prepend("src")
    from resource_planning import plan_parallelism

    with pytest.raises(ValueError, match="cannot exceed"):
        plan_parallelism(4, 5, 10)


def test_scheduler_allocation_is_detected_and_enforced(monkeypatch):
    monkeypatch.syspath_prepend("src")
    from resource_planning import (
        detect_scheduler_allocation,
        ensure_scheduler_capacity,
    )

    allocation = detect_scheduler_allocation({"SLURM_CPUS_PER_TASK": "8"})

    assert allocation["allocated_cores"] == 8
    assert allocation["scheduler"] == "slurm"
    with pytest.raises(ValueError, match="scheduler allocation"):
        ensure_scheduler_capacity(20, allocation)


def test_qc_step2_receives_per_sample_share_of_total_budget(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("src")
    import cftk
    import init
    from analysis import qc
    from visualization import visualization

    cfg = {
        "comparison": "Control_vs_Case",
        "samples": {
            "Control": [{"name": "control", "group": "Control"}],
            "Case": [{"name": "case", "group": "Case"}],
        },
        "process": {
            "parallel_samples": 2,
            "step4_methylation": {"params": {"cores": 20}},
        },
        "analysis": {"qc": {"params": {}}},
        "reference_data": {"genome_fa": str(tmp_path / "reference.fa")},
    }
    paths = {
        "qc": str(tmp_path / "qc"),
        "cpg_matrix": str(tmp_path / "matrix"),
    }
    calls = []
    monkeypatch.setattr(cftk, "_load", lambda args: (cfg, paths))
    monkeypatch.setattr(init, "get_bam", lambda sample, paths: f"{sample['name']}.bam")
    monkeypatch.setattr(
        qc,
        "run_qc",
        lambda args: calls.append((args.step, args.cores, args.parallel)),
    )
    monkeypatch.setattr(visualization, "plot_qc", lambda args: None)

    cftk._cmd_qc(SimpleNamespace(
        config="unused.json", step=[2, 0], parallel=None, force=False, title=None
    ))

    assert calls == [(2, 10, 2), (0, 20, 2)]


def test_qc_accepts_one_group_processing_project(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend("src")
    import cftk
    import init
    from analysis import qc
    from visualization import visualization

    cfg = {
        "comparison": None,
        "samples": {
            "Control": [{"name": "control", "group": "Control"}],
        },
        "process": {
            "parallel_samples": 1,
            "step4_methylation": {"params": {"cores": 20}},
        },
        "analysis": {"qc": {"params": {}}},
        "reference_data": {"genome_fa": str(tmp_path / "reference.fa")},
    }
    paths = {
        "qc": str(tmp_path / "qc"),
        "cpg_matrix": str(tmp_path / "matrix"),
    }
    calls = []
    monkeypatch.setattr(cftk, "_load", lambda args: (cfg, paths))
    monkeypatch.setattr(init, "get_bam", lambda sample, paths: f"{sample['name']}.bam")
    monkeypatch.setattr(
        qc,
        "run_qc",
        lambda args: calls.append((args.step, args.cores, args.parallel, args.group_labels)),
    )
    monkeypatch.setattr(visualization, "plot_qc", lambda args: None)

    cftk._cmd_qc(SimpleNamespace(
        config="unused.json", step=[2], parallel=None, force=False, title=None
    ))

    assert calls == [(2, 20, 1, {"Control": ["control"]})]
