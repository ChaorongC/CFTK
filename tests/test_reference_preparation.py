from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def init_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import init

    return init


def test_reference_preparation_falls_back_to_bwameth_py(
    init_module, monkeypatch, tmp_path
):
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nA\n")
    calls = []

    def fake_run(args, **kwargs):
        calls.append(args)
        if args[:2] == ["bwameth.py", "index"]:
            converted = Path(f"{reference}.bwameth.c2t")
            for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa"):
                Path(f"{converted}{suffix}").touch()
        elif args[:2] == ["samtools", "faidx"]:
            Path(f"{reference}.fai").touch()
        elif args[:2] == ["picard", "CreateSequenceDictionary"]:
            reference.with_suffix(".dict").touch()
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("subprocess.run", fake_run)
    init_module.prepare_reference(str(reference))

    assert calls[0] == ["bwameth", "index", str(reference)]
    assert calls[1] == ["bwameth.py", "index", str(reference)]
    assert ["samtools", "faidx", str(reference)] in calls
    assert [
        "picard", "CreateSequenceDictionary",
        f"R={reference}", f"O={reference.with_suffix('.dict')}",
    ] in calls


def test_reference_preparation_stops_when_both_bwameth_commands_fail(
    init_module, monkeypatch, tmp_path
):
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nA\n")
    monkeypatch.setattr(
        "subprocess.run", lambda *args, **kwargs: SimpleNamespace(returncode=1)
    )

    with pytest.raises(SystemExit, match="bwameth.*bwameth.py"):
        init_module.prepare_reference(str(reference))


def test_reference_preparation_rejects_fai_chromosome_mismatch(
    init_module, tmp_path
):
    reference = tmp_path / "hg38.fa"
    reference.write_text(">chr1\nAAAA\n")
    converted = Path(f"{reference}.bwameth.c2t")
    for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{converted}{suffix}").touch()
    Path(f"{reference}.fai").write_text("chr1\t5\t0\t0\t0\n")
    reference.with_suffix(".dict").touch()
    chrom_sizes = tmp_path / "hg38.chrom.sizes"
    chrom_sizes.write_text("chr1\t4\n")

    with pytest.raises(SystemExit, match="FASTA index.*chromosome sizes.*chr1"):
        init_module.prepare_reference(
            str(reference), chrom_sizes=str(chrom_sizes)
        )


@pytest.mark.parametrize("min_depth", [0, -1, 10.5, "10", True])
def test_config_rejects_invalid_min_depth(init_module, min_depth):
    config = {
        "project_name": "project",
        "output_dir": "/tmp/project",
        "comparison": "Control_vs_Disease",
        "samples": {
            "Control": [{
                "name": "control", "input_type": "bam", "bam": "control.bam"
            }],
            "Disease": [{
                "name": "disease", "input_type": "bam", "bam": "disease.bam"
            }],
        },
        "reference_data": {
            "genome_fa": "hg38.fa",
            "genome_2bit": "hg38.2bit",
            "chrom_sizes": "hg38.chrom.sizes",
        },
        "process": {
            "step1_trimming": {},
            "step2_alignment": {},
            "step3_markdup": {},
            "step4_methylation": {"params": {"min_depth": min_depth}},
        },
        "analysis": {},
    }

    assert any("min_depth" in error for error in init_module._validate(config))
