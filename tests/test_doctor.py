import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from test_reference_profiles import make_profile


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def modules(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import cftk
    import doctor
    import init

    return cftk, doctor, init


def _write_fastq(path, read_number, name="pair"):
    path.write_text(
        f"@{name}/{read_number}\n"
        "ACGT\n"
        "+\n"
        "IIII\n"
    )


def _make_project(tmp_path, init_module, *, reference_mode="local"):
    reads = tmp_path / "reads"
    reads.mkdir()
    for sample in ("control", "case"):
        _write_fastq(reads / f"{sample}_R1.fq", 1, sample)
        _write_fastq(reads / f"{sample}_R2.fq", 2, sample)
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n"
        "control\tHealthy\tcontrol\tfastq\treads/control_R1.fq\treads/control_R2.fq\t\n"
        "case\tDisease\tcase\tfastq\treads/case_R1.fq\treads/case_R2.fq\t\n"
    )
    reference_root = tmp_path / "references"
    profile_dir = make_profile(reference_root, checksums=True)
    fasta = profile_dir / "genome/hg38.fa"
    (Path(f"{fasta}.fai")).write_text("chr1\t4\t6\t4\t5\n")
    (profile_dir / "genome/hg38.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n"
    )
    converted = Path(f"{fasta}.bwameth.c2t")
    for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{converted}{suffix}").touch()

    config = {
        "schema_version": 2,
        "project_name": "doctor-test",
        "output_dir": ".",
        "assay": "twist_human_methylome",
        "genome": "hg38",
        "samples": "samples.tsv",
        "reference_mode": "local",
        "reference_root": str(reference_root),
        "reference_profile": {
            "id": "twist_human_methylome_hg38",
            "version": "1.0.0",
        },
    }
    config_path = tmp_path / "cftk_init.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    init_module.write_lockfile(config_path)
    if reference_mode == "managed":
        config["reference_mode"] = "managed"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        lock_path = tmp_path / "cftk.lock.json"
        lock = json.loads(lock_path.read_text())
        lock["project_config_sha256"] = hashlib.sha256(
            config_path.read_bytes()
        ).hexdigest()
        lock_path.write_text(json.dumps(lock, indent=2) + "\n")
    return config_path


def _fake_tools(tmp_path, monkeypatch, names):
    bindir = tmp_path / "bin"
    bindir.mkdir()
    for name in names:
        path = bindir / name
        path.write_text("#!/bin/sh\necho 'test version 1.0'\n")
        path.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir))


def _args(config, **overrides):
    values = {
        "config": str(config),
        "step": [1, 2, 3, 4],
        "target_bed": None,
        "skip_picard_metrics": False,
        "json": True,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_parser_registers_doctor_with_all_steps_by_default(modules):
    cftk, _, _ = modules
    args = cftk.build_parser().parse_args(["doctor"])

    assert args.step == [1, 2, 3, 4]
    assert args.target_bed is None
    assert not args.skip_picard_metrics
    assert not args.json
    assert '"doctor"' in (REPO_ROOT / "pyproject.toml").read_text()


def test_doctor_uses_picard_compatible_version_probe(modules):
    _, doctor, _ = modules
    cfg = {
        "samples": {},
        "process": {
            "step3_markdup": {"tool": "sambamba"},
        },
    }

    requirement = doctor._tool_requirements(cfg, {3}, False)["picard"]

    assert requirement == (
        "picard",
        ("MarkDuplicates", "--version"),
        (0, 1),
        False,
    )
    assert doctor._TOOL_PROBE_TIMEOUT == 60


def test_doctor_does_not_acquire_an_installed_managed_profile(
    modules, tmp_path, monkeypatch
):
    _, doctor, init_module = modules
    config = _make_project(tmp_path, init_module, reference_mode="managed")
    _fake_tools(
        tmp_path,
        monkeypatch,
        [
            "trim_galore", "fastqc", "bwameth.py", "bwa", "sambamba",
            "samtools", "picard", "MethylDackel", "bedtools", "multiqc",
        ],
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("doctor attempted managed acquisition")

    monkeypatch.setattr("reference_profiles._acquire_managed_profile", forbidden)
    report = doctor.run_doctor(_args(config))

    assert report["status"] == "PASS"
    assert report["exit_code"] == 0


def test_doctor_json_output_is_clean_and_missing_tools_fail(
    modules, tmp_path, monkeypatch, capsys
):
    cftk, _, init_module = modules
    config = _make_project(tmp_path, init_module)
    monkeypatch.setenv("PATH", str(tmp_path / "empty-bin"))
    args = _args(config)

    with pytest.raises(SystemExit) as exc:
        cftk._cmd_doctor(args)

    payload = json.loads(capsys.readouterr().out)
    assert exc.value.code == 1
    assert payload["status"] == "FAIL"
    assert any(
        check["id"] == "tool.trim_galore" and check["status"] == "FAIL"
        for check in payload["checks"]
    )


def test_doctor_rejects_mismatched_fastq_mates(modules, tmp_path, monkeypatch):
    _, doctor, init_module = modules
    config = _make_project(tmp_path, init_module)
    _write_fastq(tmp_path / "reads/control_R2.fq", 2, "different")
    _fake_tools(tmp_path, monkeypatch, ["trim_galore", "fastqc", "multiqc"])

    report = doctor.run_doctor(_args(config, step=[1]))

    assert report["status"] == "FAIL"
    assert any(
        check["id"] == "input.fastq.control" and "mate" in check["summary"]
        for check in report["checks"]
    )


def test_doctor_step4_rejects_bam_reference_dictionary_mismatch(
    modules, tmp_path, monkeypatch
):
    pysam = pytest.importorskip("pysam")
    _, doctor, init_module = modules
    config_path = _make_project(tmp_path, init_module)
    raw = json.loads(config_path.read_text())
    rows = ["sample\tgroup\trole\tinput_type\tr1\tr2\tbam"]
    for sample, group, role in (
        ("control", "Healthy", "control"),
        ("case", "Disease", "case"),
    ):
        bam = tmp_path / f"{sample}.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr2", "LN": 4}],
            "RG": [{"ID": sample, "SM": sample, "LB": sample, "PL": "ILLUMINA"}],
            "PG": [{"ID": "markdup", "PN": "sambamba", "CL": "sambamba markdup"}],
        }
        with pysam.AlignmentFile(bam, "wb", header=header):
            pass
        pysam.index(str(bam))
        rows.append(f"{sample}\t{group}\t{role}\tbam\t\t\t{bam.name}")
    (tmp_path / "samples.tsv").write_text("\n".join(rows) + "\n")
    init_module.write_lockfile(config_path)
    profile_dir = tmp_path / "references/twist_human_methylome_hg38/1.0.0"
    (profile_dir / "genome/hg38.dict").unlink()
    for path in profile_dir.glob("genome/hg38.fa.bwameth.c2t*"):
        path.unlink()
    _fake_tools(tmp_path, monkeypatch, ["samtools", "MethylDackel", "bedtools"])

    report = doctor.run_doctor(_args(config_path, step=[4]))

    assert report["status"] == "FAIL"
    assert any(
        check["id"] == "input.bam.control.dictionary"
        and check["status"] == "FAIL"
        for check in report["checks"]
    )
    assert not any(
        check["id"] in {"reference.dict", "reference.bwameth"}
        for check in report["checks"]
    )


def test_lock_hash_mismatch_is_a_readiness_failure(modules, tmp_path, monkeypatch):
    _, doctor, init_module = modules
    config = _make_project(tmp_path, init_module)
    lock_path = tmp_path / "cftk.lock.json"
    lock = json.loads(lock_path.read_text())
    lock["sample_sheet"]["sha256"] = hashlib.sha256(b"wrong").hexdigest()
    lock_path.write_text(json.dumps(lock))
    _fake_tools(tmp_path, monkeypatch, ["trim_galore", "fastqc", "multiqc"])

    report = doctor.run_doctor(_args(config, step=[1]))

    assert any(
        check["id"] == "project.lock" and check["status"] == "FAIL"
        for check in report["checks"]
    )
