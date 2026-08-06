import json
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def scope_module(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    import analysis.assay_scope as assay_scope

    return assay_scope


def _write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _context(tmp_path, *, assay="twist_human_methylome", scope="auto"):
    target = _write(tmp_path / "target.bed", "chr1\t10\t30\ttarget\n")
    regions = _write(
        tmp_path / "tss_pas.bed",
        "chr1\t0\t15\tregion_a\nchr1\t25\t40\tregion_b\n",
    )
    bins = _write(
        tmp_path / "bins.bed",
        "chr1\t0\t20\tbin_a\nchr1\t20\t40\tbin_b\n",
    )
    bam = _write(tmp_path / "s1.bam", "placeholder\n")
    _write(str(bam) + ".bai", "index\n")
    cfg = {
        "assay": assay,
        "reference_profile": {"id": "twist_human_methylome_hg38"},
        "reference_data": {
            "target_bed": str(target),
            "tss_pas_bed": str(regions),
            "bins": str(bins),
        },
        "analysis": {"frag": {"scope": scope}},
    }
    paths = {"fragmentomics": str(tmp_path / "results" / "4_fragmentomics")}
    return cfg, paths, [{"name": "s1"}], [str(bam)]


def test_scope_resolution_is_panel_for_twist_and_genome_for_legacy(scope_module, tmp_path):
    cfg, _paths, _samples, _bams = _context(tmp_path)

    assert scope_module.resolve_scope(cfg)["mode"] == "panel"
    assert scope_module.resolve_scope(cfg, "genome")["mode"] == "genome"

    legacy = dict(cfg)
    legacy.pop("assay")
    legacy.pop("reference_profile")
    assert scope_module.resolve_scope(legacy)["mode"] == "genome"


def test_clip_bed_preserves_source_order_and_clips_to_panel(scope_module, tmp_path):
    source = _write(
        tmp_path / "regions.bed",
        "chr1\t0\t15\tregion_a\nchr1\t25\t40\tregion_b\nchr2\t0\t10\tother\n",
    )
    target = _write(
        tmp_path / "target.bed",
        "chr1\t10\t30\ttarget\nchr1\t28\t35\toverlap\n",
    )
    destination = tmp_path / "clipped.bed"

    count = scope_module._clip_bed(source, target, destination)

    assert count == 2
    assert destination.read_text(encoding="utf-8").splitlines() == [
        "chr1\t10\t15\tregion_a__panel_0",
        "chr1\t25\t35\tregion_b__panel_0",
    ]


def test_prepare_scope_derives_regions_bins_and_panel_bam(scope_module, tmp_path, monkeypatch):
    cfg, paths, samples, bams = _context(tmp_path)
    calls = []

    def fake_recorded_run(command, *, label="", **kwargs):
        calls.append((list(command), label))
        if command[1] == "view":
            output = Path(command[command.index("-o") + 1])
            output.write_bytes(b"panel bam")
        elif command[1] == "index":
            Path(str(command[-1]) + ".bai").write_bytes(b"index")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(scope_module, "recorded_run", fake_recorded_run)

    prepared = scope_module.prepare_scope(
        cfg,
        paths,
        samples,
        bams,
        kinds={"occupancy", "wps", "delfi"},
        cores=4,
    )

    assert prepared["info"]["mode"] == "panel"
    assert prepared["info"]["region_count"] == 2
    assert prepared["info"]["bins_count"] == 2
    assert Path(prepared["region_bed"]).is_file()
    assert Path(prepared["bins"]).is_file()
    assert Path(prepared["bam_paths"][0]).is_file()
    assert len(calls) == 2
    assert calls[0][0][0:2] == ["samtools", "view"]
    assert calls[1][0][0:2] == ["samtools", "index"]
    assert Path(prepared["info"]["scope_json"]).is_file()

    scope_module.prepare_scope(
        cfg, paths, samples, bams, kinds={"occupancy"}, cores=4
    )
    assert len(calls) == 2


def test_prepare_scope_rejects_empty_panel_overlap(scope_module, tmp_path):
    cfg, paths, samples, bams = _context(tmp_path)
    _write(tmp_path / "tss_pas.bed", "chr2\t0\t10\tno_overlap\n")

    with pytest.raises(scope_module.ScopeError, match="no overlap"):
        scope_module.prepare_scope(
            cfg, paths, samples, bams, kinds={"wps"}, cores=1
        )


def test_scope_artifacts_are_not_required_for_explicit_genome_mode(scope_module, tmp_path):
    cfg, paths, samples, _bams = _context(tmp_path)

    assert scope_module.scope_artifact_paths(
        cfg, paths, [sample["name"] for sample in samples], "wps", requested="genome"
    ) == []


def test_downstream_plan_records_panel_scope_and_intermediate_contract(
    scope_module, tmp_path, monkeypatch
):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    import analysis_workflow

    cfg, paths, samples, _bams = _context(tmp_path)
    scope = scope_module.describe_scope(
        cfg,
        paths,
        [sample["name"] for sample in samples],
        requested="panel",
        bam_paths=_bams,
    )
    context = {
        "cfg": cfg,
        "paths": {
            **paths,
            "provenance": str(tmp_path / "results" / "provenance"),
            "work_dir": str(tmp_path),
            "differential": str(tmp_path / "results" / "3_differential"),
            "occ_out": str(tmp_path / "results" / "4_fragmentomics" / "occupancy"),
            "wps_out": str(tmp_path / "results" / "4_fragmentomics" / "wps"),
            "delfi_out": str(tmp_path / "results" / "4_fragmentomics" / "delfi"),
            "end_motif_out": str(tmp_path / "results" / "4_fragmentomics" / "end_motif"),
            "cleavage_out": str(tmp_path / "results" / "4_fragmentomics" / "cleavage"),
            "mesa": str(tmp_path / "results" / "5_mesa"),
            "report": str(tmp_path / "results" / "report"),
            "methylation": str(tmp_path / "results" / "1_process" / "4_methylation"),
            "markdup": str(tmp_path / "results" / "1_process" / "3_markdup"),
            "alignment": str(tmp_path / "results" / "1_process" / "2_alignment"),
            "cpg_matrix": str(tmp_path / "results" / "1_process" / "5_merged_matrix"),
        },
        "samples": samples,
        "config_path": tmp_path / "cftk_init.json",
        "lock_path": tmp_path / "cftk.lock.json",
        "identity": {"config_sha256": "config", "lock_sha256": "lock"},
        "fragmentomics_scope_request": "panel",
        "fragmentomics_scope": scope,
    }
    args = SimpleNamespace(
        preset="fragmentomics", stages=["wps"], parallel=None,
        fragmentomics_scope="panel",
    )

    payload = analysis_workflow.build_plan(context, args)
    expected = {Path(item["path"]).name for item in payload["stages"][0]["expected"]}

    assert payload["fragmentomics_scope"]["mode"] == "panel"
    assert "target_overlap_regions.bed" in expected
    assert "--fragmentomics-scope panel" in payload["stages"][0]["command"]


def test_downstream_plan_dry_run_persists_panel_scope(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    import analysis.assay_scope as assay_scope
    import analysis_workflow

    cfg, paths, samples, bams = _context(tmp_path)
    paths = {
        **paths,
        "provenance": str(tmp_path / "results" / "provenance"),
        "work_dir": str(tmp_path),
        "differential": str(tmp_path / "results" / "3_differential"),
        "occ_out": str(tmp_path / "results" / "4_fragmentomics" / "occupancy"),
        "wps_out": str(tmp_path / "results" / "4_fragmentomics" / "wps"),
        "delfi_out": str(tmp_path / "results" / "4_fragmentomics" / "delfi"),
        "end_motif_out": str(tmp_path / "results" / "4_fragmentomics" / "end_motif"),
        "cleavage_out": str(tmp_path / "results" / "4_fragmentomics" / "cleavage"),
        "mesa": str(tmp_path / "results" / "5_mesa"),
        "report": str(tmp_path / "results" / "report"),
        "methylation": str(tmp_path / "results" / "1_process" / "4_methylation"),
        "markdup": str(tmp_path / "results" / "1_process" / "3_markdup"),
        "alignment": str(tmp_path / "results" / "1_process" / "2_alignment"),
        "cpg_matrix": str(tmp_path / "results" / "1_process" / "5_merged_matrix"),
    }
    config_path = tmp_path / "cftk_init.json"
    lock_path = tmp_path / "cftk.lock.json"
    config_path.write_text("{}\n")
    lock_path.write_text("{}\n")
    context = {
        "cfg": cfg,
        "paths": paths,
        "samples": samples,
        "config_path": config_path,
        "lock_path": lock_path,
        "identity": {"config_sha256": "config", "lock_sha256": "lock"},
        "fragmentomics_scope_request": "panel",
        "fragmentomics_scope": assay_scope.describe_scope(
            cfg, paths, [sample["name"] for sample in samples],
            requested="panel", bam_paths=bams,
        ),
    }
    monkeypatch.setattr(analysis_workflow, "_load_context", lambda args: context)
    import doctor
    monkeypatch.setattr(
        doctor,
        "run_doctor",
        lambda args: {"status": "PASS", "exit_code": 0, "checks": []},
    )

    payload = analysis_workflow.plan(SimpleNamespace(
        config=str(config_path), preset="fragmentomics", stages=["wps"],
        parallel=None, fragmentomics_scope="panel", json=False,
    ))

    latest = Path(paths["provenance"]) / "latest-analysis-plan.json"
    assert payload["status"] == "ready"
    assert json.loads(latest.read_text(encoding="utf-8"))["status"] == "ready"
