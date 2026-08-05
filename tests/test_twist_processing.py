import json
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def process_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import process

    return process


def _sample(name="sample_1"):
    return {"name": name, "input_type": "fastq"}


def _paths(tmp_path):
    return {
        "alignment": str(tmp_path / "alignment"),
        "markdup": str(tmp_path / "markdup"),
        "methylation": str(tmp_path / "methylation"),
    }


def test_bwameth_alignment_generates_sample_read_group(
    process_module, monkeypatch, tmp_path
):
    commands = []
    monkeypatch.setattr(
        process_module, "run_command", lambda cmd, label="": commands.append(cmd)
    )

    process_module._step2_align(
        _sample(), "reads_R1.fq.gz", "reads_R2.fq.gz",
        {"tool": "bwameth", "params": {}},
        {"genome_fa": "hg38.fa"}, _paths(tmp_path), 4,
    )

    assert "--read-group" in commands[0]
    assert "@RG\\tID:sample_1\\tSM:sample_1\\tLB:sample_1\\tPL:ILLUMINA" in commands[0]


def test_trim_galore_basename_reports_are_renamed_for_qc_parser(
    process_module, monkeypatch, tmp_path
):
    trimming = tmp_path / "trimming"
    paths = {**_paths(tmp_path), "trimming": str(trimming)}
    sample = {
        "name": "sample_1",
        "input_type": "fastq",
        "r1": "sample_1_R1.fastq.gz",
        "r2": "sample_1_R2.fastq.gz",
    }

    def fake_run(command, label=""):
        trimming.mkdir(parents=True, exist_ok=True)
        for filename in (
            "sample_1_val_1.fq.gz",
            "sample_1_val_2.fq.gz",
            "sample_1_val_1_fastqc.html",
            "sample_1_val_1_fastqc.zip",
            "sample_1_val_2_fastqc.html",
            "sample_1_val_2_fastqc.zip",
            "sample_1_R1.fastq.gz_trimming_report.txt",
            "sample_1_R2.fastq.gz_trimming_report.txt",
        ):
            (trimming / filename).write_text("report\n")

    monkeypatch.setattr(process_module, "run_command", fake_run)

    process_module._step1_trim(
        sample,
        {"tool": "trim_galore", "params": {}},
        {},
        paths,
        4,
    )

    assert (trimming / "sample_1_R1_trimming_report.txt").is_file()
    assert (trimming / "sample_1_R2_trimming_report.txt").is_file()
    assert not (trimming / "sample_1_R1.fastq.gz_trimming_report.txt").exists()
    assert not (trimming / "sample_1_R2.fastq.gz_trimming_report.txt").exists()


def test_trim_checkpoint_fix_preserves_canonical_reports(
    process_module, tmp_path
):
    trimming = tmp_path / "trimming"
    trimming.mkdir()
    sample = {
        "name": "sample_1",
        "input_type": "fastq",
        "r1": "sample_1_R1.fastq.gz",
        "r2": "sample_1_R2.fastq.gz",
    }
    for mate in ("R1", "R2"):
        (trimming / f"sample_1_{mate}_trimming_report.txt").write_text(
            f"Input filename: sample_1_{mate}\n"
        )

    process_module._apply_trim_rename_fix(sample, str(trimming), "fq.gz")

    for mate in ("R1", "R2"):
        report = trimming / f"sample_1_{mate}_trimming_report.txt"
        assert report.read_text() == f"Input filename: sample_1_{mate}\n"


@pytest.mark.parametrize("params, expected_depth", [({}, 10), ({"min_depth": 17}, 17)])
def test_methyldackel_merges_cpg_context_without_chh_or_chg(
    process_module, monkeypatch, tmp_path, params, expected_depth
):
    commands = []
    mbias = SimpleNamespace(
        returncode=0,
        stdout="Strand\tRead\tPosition\tMethylated\tUnmethylated\n",
        stderr="Suggested inclusion options: --OT 3,0,0,3 --OB 0,3,3,0\n",
    )
    monkeypatch.setattr("subprocess.run", lambda *args, **kwargs: mbias)
    monkeypatch.setattr(
        process_module, "run_command", lambda cmd, label="": commands.append(cmd)
    )

    process_module._step4_methylation(
        _sample(), "sample.markdup.bam",
        {"tool": "methyldackel", "params": params},
        {"genome_fa": "hg38.fa"}, _paths(tmp_path), 4,
    )

    assert len(commands) == 1
    command = commands[0]
    assert "--mergeContext" in command
    assert f"--minDepth {expected_depth}" in command
    assert "--OT 3,0,0,3 --OB 0,3,3,0" in command
    assert "--CHH" not in command
    assert "--CHG" not in command
    assert "--noCpG" not in command


def test_methyldackel_stops_when_ot_ob_cannot_be_parsed(
    process_module, monkeypatch, tmp_path
):
    commands = []
    mbias = SimpleNamespace(
        returncode=0,
        stdout="Strand\tRead\tPosition\tMethylated\tUnmethylated\n",
        stderr="No suggested bounds were emitted\n",
    )
    monkeypatch.setattr("subprocess.run", lambda *args, **kwargs: mbias)
    monkeypatch.setattr(
        process_module, "run_command", lambda cmd, label="": commands.append(cmd)
    )

    with pytest.raises(SystemExit, match="could not parse OT/OB"):
        process_module._step4_methylation(
            _sample(), "sample.markdup.bam",
            {"tool": "methyldackel", "params": {}},
            {"genome_fa": "hg38.fa"}, _paths(tmp_path), 4,
        )

    assert commands == []


def test_picard_interval_and_metrics_commands_checkpoint_independently(
    process_module, monkeypatch, tmp_path
):
    commands = []
    paths = _paths(tmp_path)
    Path(paths["markdup"]).mkdir()
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nA\n")
    reference.with_suffix(".dict").write_text("@HD\tVN:1.6\n")
    target_bed = tmp_path / "targets.bed"
    target_bed.write_text("chr1\t0\t1\n")
    monkeypatch.setattr(
        process_module, "run_command", lambda cmd, label="": commands.append(cmd)
    )

    interval_list = process_module._prepare_picard_interval_list(
        str(target_bed), str(reference), paths
    )
    process_module._run_picard_metrics(
        _sample(), "sample.markdup.bam", str(reference), interval_list, paths
    )

    assert "picard -Xmx8g BedToIntervalList" in commands[0]
    assert f"SD={reference.with_suffix('.dict')}" in commands[0]
    assert "picard -Xmx8g CollectHsMetrics" in commands[1]
    assert f"BAIT_INTERVALS={interval_list}" in commands[1]
    assert f"TARGET_INTERVALS={interval_list}" in commands[1]
    assert "PER_TARGET_COVERAGE=" in commands[1]
    assert "MINIMUM_MAPPING_QUALITY=20" in commands[1]
    assert "COVERAGE_CAP=1000" in commands[1]
    assert "NEAR_DISTANCE=500" in commands[1]
    assert "picard -Xmx8g CollectMultipleMetrics" in commands[2]
    assert "PROGRAM=CollectGcBiasMetrics" in commands[2]
    assert "PROGRAM=CollectInsertSizeMetrics" in commands[2]
    assert "PROGRAM=CollectAlignmentSummaryMetrics" in commands[2]


def test_picard_java_memory_is_tunable_and_shell_safe(process_module):
    assert process_module._picard_command("CollectHsMetrics", "12G") == (
        "picard -Xmx12g CollectHsMetrics"
    )
    with pytest.raises(SystemExit, match="picard_java_memory"):
        process_module._picard_command("CollectHsMetrics", "8g;touch bad")


def test_existing_markdup_bam_can_receive_missing_picard_metrics(
    process_module, monkeypatch, tmp_path
):
    paths = _paths(tmp_path)
    Path(paths["markdup"]).mkdir()
    bam = Path(paths["markdup"]) / "sample_1.markdup.bam"
    bam.touch()
    monkeypatch.setattr(
        process_module,
        "_step3_markdup",
        lambda *args, **kwargs: pytest.fail("markdup should remain checkpointed"),
    )
    seen = []
    monkeypatch.setattr(
        process_module,
        "_run_picard_metrics",
        lambda sample, bam_in, ref, intervals, work_paths, memory: seen.append(
            (bam_in, memory)
        ),
    )

    result = process_module._run_step3_with_metrics(
        _sample(),
        {"step3_markdup": {"tool": "sambamba", "params": {}}},
        {"genome_fa": "hg38.fa"}, paths, 4, "targets.interval_list",
    )

    assert result == str(bam)
    assert seen == [(str(bam), "8g")]


def test_target_bed_resolution_prefers_cli_override(process_module, tmp_path):
    target_bed = tmp_path / "custom.bed"
    target_bed.touch()
    assert process_module._resolve_target_bed(str(target_bed)) == str(target_bed.resolve())


def test_schema_v2_profile_target_is_used_without_cli_override(
    process_module, monkeypatch, tmp_path
):
    profile_target = str((tmp_path / "profile-targets.bed").resolve())
    reference = str((tmp_path / "hg38.fa").resolve())
    process_config = {
        "parallel_samples": 1,
        "step1_trimming": {"params": {"cores": 1}},
        "step2_alignment": {"params": {"cores": 1}},
        "step3_markdup": {"params": {"cores": 1}},
        "step4_methylation": {"params": {"cores": 1}},
    }
    config = {
        "reference_data": {"genome_fa": reference, "target_bed": profile_target},
        "process": process_config,
    }
    seen = []
    monkeypatch.setattr(process_module, "load_config", lambda path: config)
    monkeypatch.setattr(process_module, "get_work_paths", lambda cfg: _paths(tmp_path))
    monkeypatch.setattr(process_module, "get_all_samples", lambda cfg: [])
    monkeypatch.setattr(
        process_module, "_resolve_target_bed", lambda path: seen.append(path) or path
    )
    monkeypatch.setattr(
        process_module,
        "_prepare_picard_interval_list",
        lambda bed, ref, paths, memory: "targets.interval_list",
    )
    monkeypatch.setattr(process_module, "_run_multiqc", lambda *args, **kwargs: None)

    process_module.process(
        SimpleNamespace(
            step=[3], parallel=None, target_bed=None, skip_picard_metrics=False
        ),
        config_path="schema-v2.json",
    )

    assert seen == [profile_target]


def test_interval_checkpoint_changes_with_target_content(
    process_module, monkeypatch, tmp_path
):
    paths = _paths(tmp_path)
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nA\n")
    reference.with_suffix(".dict").write_text("@HD\tVN:1.6\n")
    target_bed = tmp_path / "targets.bed"
    target_bed.write_text("chr1\t0\t1\n")
    monkeypatch.setattr(process_module, "run_command", lambda *args, **kwargs: None)

    first = process_module._prepare_picard_interval_list(
        str(target_bed), str(reference), paths
    )
    target_bed.write_text("chr1\t0\t2\n")
    second = process_module._prepare_picard_interval_list(
        str(target_bed), str(reference), paths
    )

    assert first != second


def test_missing_target_bed_has_actionable_error(process_module, monkeypatch, tmp_path):
    monkeypatch.setattr(process_module, "DEFAULT_TARGET_BED", tmp_path / "missing.bed")
    with pytest.raises(SystemExit, match="--target-bed.*--skip-picard-metrics"):
        process_module._resolve_target_bed(None)


@pytest.mark.parametrize(
    "config_path", [REPO_ROOT / "cftk_init.json", REPO_ROOT / "docs/_static/cftk_init.json"]
)
def test_distributed_configs_expose_min_depth_without_clip_parameters(config_path):
    config = json.loads(config_path.read_text())
    methylation = config["process"]["step4_methylation"]["params"]
    qc_params = config["analysis"]["qc"]["params"]

    assert methylation["min_depth"] == 10
    assert "clip_r1" not in qc_params
    assert "clip_r2" not in qc_params
    assert "clip" not in config["process"]["step1_trimming"]["params"]["extra_args"]


def test_chh_conversion_is_not_part_of_default_qc_collection(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    from analysis import qc_scorer

    assert "bisulfite_conversion_rate" not in {
        rule.col for rule in qc_scorer.RULES if rule.weight > 0
    }


def test_cpg_matrix_uses_one_based_cytosine_for_merged_context(
    process_module, monkeypatch, tmp_path
):
    paths = {"cpg_matrix": str(tmp_path / "matrix")}
    bedgraphs = [str(tmp_path / "sample_1_CpG.bedGraph")]
    Path(bedgraphs[0]).touch()
    monkeypatch.setattr("shutil.which", lambda executable: "/usr/bin/bedtools")

    def fake_run(command, shell):
        tmp_output = Path(paths["cpg_matrix"]) / "cpg_matrix.tsv.tmp"
        tmp_output.write_text(
            "chrom\tstart\tend\tsample_1\n"
            "chr1\t25114\t25116\t83\n"
            "chr1\t29336\t29337\t50\n"
        )
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    output = process_module._merge_cpg(
        bedgraphs, [_sample()], paths
    )

    rows = Path(output).read_text().splitlines()
    assert rows[1].startswith("chr1_25115\t")
    assert rows[2].startswith("chr1_29337\t")
