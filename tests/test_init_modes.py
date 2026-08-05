import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from test_reference_profiles import make_profile


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def init_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import init

    return init


def make_sample_sheet(project_dir):
    reads = project_dir / "reads"
    reads.mkdir()
    for name in ("control_R1.fq.gz", "control_R2.fq.gz", "case_R1.fq.gz", "case_R2.fq.gz"):
        (reads / name).touch()
    sheet = project_dir / "samples.tsv"
    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n"
        "control\tHealthy\tcontrol\tfastq\treads/control_R1.fq.gz\treads/control_R2.fq.gz\t\n"
        "case\tDisease\tcase\tfastq\treads/case_R1.fq.gz\treads/case_R2.fq.gz\t\n"
    )
    return sheet


def compact_config(project_dir, reference_root):
    return {
        "schema_version": 2,
        "project_name": "study",
        "output_dir": ".",
        "assay": "twist_human_methylome",
        "genome": "hg38",
        "samples": "samples.tsv",
        "reference_root": str(reference_root),
        "reference_profile": {
            "id": "twist_human_methylome_hg38",
            "version": "1.0.0",
        },
    }


def init_args(config_path, **overrides):
    defaults = {
        "config": str(config_path),
        "non_interactive": True,
        "sample_sheet": None,
        "reference_root": None,
        "reference_mode": None,
        "profile": None,
        "profile_version": None,
        "project_name": None,
        "output_dir": None,
        "assay": "twist_human_methylome",
        "genome": "hg38",
        "skip_reference_prep": True,
        "ref_index": False,
        "ref_dict": False,
    }
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_sample_sheet_preserves_order_roles_and_resolves_paths(init_module, tmp_path):
    sheet = make_sample_sheet(tmp_path)
    parsed = init_module.load_sample_sheet(str(sheet))

    assert list(parsed["samples"]) == ["Healthy", "Disease"]
    assert parsed["group_roles"] == {"Healthy": "control", "Disease": "case"}
    assert parsed["samples"]["Healthy"][0]["r1"] == str(
        (tmp_path / "reads/control_R1.fq.gz").resolve()
    )


@pytest.mark.parametrize(
    "rows, message",
    [
        (
            "a\tA\tcontrol\tfastq\ta_R1.fq.gz\t\t\n",
            "requires both r1 and r2",
        ),
        (
            "a\tA\tcontrol\tbam\t\t\ta.bam\n"
            "b\tB\tcontrol\tbam\t\t\tb.bam\n"
            "c\tC\tcase\tbam\t\t\tc.bam\n",
            "exactly one group for role 'control'",
        ),
    ],
)
def test_sample_sheet_rejects_invalid_contract(init_module, tmp_path, rows, message):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n" + rows)
    with pytest.raises(SystemExit, match=message):
        init_module.load_sample_sheet(str(sheet), require_files=False)


def test_sample_sheet_rejects_unexpected_columns(init_module, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\tcondition\n"
        "a\tA\tcontrol\tbam\t\t\ta.bam\thealthy\n"
    )

    with pytest.raises(SystemExit, match="unexpected columns.*condition"):
        init_module.load_sample_sheet(str(sheet), require_files=False)


def test_sample_sheet_rejects_conflicting_or_reused_inputs(init_module, tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n"
        "a\tA\tcontrol\tfastq\ta_R1.fq.gz\ta_R2.fq.gz\ta.bam\n"
        "b\tB\tcase\tfastq\tb_R1.fq.gz\tb_R2.fq.gz\t\n"
    )
    with pytest.raises(SystemExit, match="must leave bam empty"):
        init_module.load_sample_sheet(str(sheet), require_files=False)

    sheet.write_text(
        "sample\tgroup\trole\tinput_type\tr1\tr2\tbam\n"
        "a\tA\tcontrol\tbam\t\t\tshared.bam\n"
        "b\tB\tcase\tbam\t\t\tshared.bam\n"
    )
    with pytest.raises(SystemExit, match="reused by multiple samples"):
        init_module.load_sample_sheet(str(sheet), require_files=False)


def test_schema_v2_resolves_to_legacy_contract(init_module, tmp_path):
    make_sample_sheet(tmp_path)
    make_profile(tmp_path / "references", checksums=True)
    config_path = tmp_path / "cftk_init.json"
    config_path.write_text(json.dumps(compact_config(tmp_path, tmp_path / "references")))

    config = init_module.load_config(str(config_path))

    assert config["comparison"] == "Healthy_vs_Disease"
    assert config["group_roles"] == {"Healthy": "control", "Disease": "case"}
    assert config["process"]["step4_methylation"]["params"]["min_depth"] == 10
    assert (
        config["process"]["step3_markdup"]["params"]["picard_java_memory"]
        == "8g"
    )
    assert config["process"]["step3_markdup"]["tool"] == "sambamba"
    assert config["reference_data"]["target_bed"].endswith("covered_targets.bed")
    assert config["output_dir"] == str(tmp_path.resolve())


def test_schema_v2_exposes_advanced_duplicate_marking_tool(
    init_module, tmp_path
):
    make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root)
    raw = compact_config(tmp_path, reference_root)
    raw["process"] = {"duplicate_marking_tool": "Picard"}

    config = init_module.resolve_schema_v2(raw, tmp_path / "cftk_init.json")

    assert config["process"]["step3_markdup"]["tool"] == "picard"


@pytest.mark.parametrize("value", ["unknown", 1, True])
def test_schema_v2_rejects_invalid_duplicate_marking_tool(
    init_module, tmp_path, value
):
    make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root)
    raw = compact_config(tmp_path, reference_root)
    raw["process"] = {"duplicate_marking_tool": value}

    with pytest.raises(SystemExit, match="duplicate_marking_tool"):
        init_module.resolve_schema_v2(raw, tmp_path / "cftk_init.json")


def test_environment_reference_root_overrides_config_hint(
    init_module, tmp_path, monkeypatch
):
    make_sample_sheet(tmp_path)
    environment_root = tmp_path / "portable-references"
    make_profile(environment_root)
    raw = compact_config(tmp_path, tmp_path / "other-machine-references")
    config_path = tmp_path / "cftk_init.json"
    config_path.write_text(json.dumps(raw))
    monkeypatch.setenv("CFTK_REFERENCE_ROOT", str(environment_root))

    config = init_module.load_config(str(config_path))

    assert config["reference_data"]["genome_fa"].startswith(
        str(environment_root.resolve())
    )


@pytest.mark.parametrize(
    "key, value",
    [
        ("cores", 0),
        ("cores", True),
        ("parallel_samples", -1),
        ("parallel_samples", "2"),
        ("min_depth", 2.5),
    ],
)
def test_schema_v2_rejects_invalid_positive_integer_process_settings(
    init_module, tmp_path, key, value
):
    make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root)
    raw = compact_config(tmp_path, reference_root)
    raw["process"] = {key: value}

    with pytest.raises(SystemExit, match=key):
        init_module.resolve_schema_v2(raw, tmp_path / "cftk_init.json")


@pytest.mark.parametrize("value", [0, "0g", "8", "eight", "8g;touch bad", True])
def test_schema_v2_rejects_invalid_picard_java_memory(
    init_module, tmp_path, value
):
    make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root)
    raw = compact_config(tmp_path, reference_root)
    raw["process"] = {"picard_java_memory": value}

    with pytest.raises(SystemExit, match="picard_java_memory"):
        init_module.resolve_schema_v2(raw, tmp_path / "cftk_init.json")


def test_noninteractive_init_creates_compact_config_and_portable_lock(
    init_module, tmp_path
):
    sheet = make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root, checksums=True)
    config_path = tmp_path / "cftk_init.json"

    init_module.init(init_args(
        config_path,
        sample_sheet=str(sheet),
        reference_root=str(reference_root),
        reference_mode="local",
        profile="twist_human_methylome_hg38",
        profile_version="1.0.0",
    ))

    compact = json.loads(config_path.read_text())
    lock = json.loads((tmp_path / "cftk.lock.json").read_text())
    assert compact["schema_version"] == 2
    assert compact["samples"] == "samples.tsv"
    assert lock["reference_profile"]["id"] == "twist_human_methylome_hg38"
    assert "reference_root" not in lock
    assert len(lock["sample_sheet"]["sha256"]) == 64
    assert len(lock["reference_profile"]["components"]["genome_fa"]) == 64


def test_noninteractive_init_requires_explicit_inputs(init_module, tmp_path):
    with pytest.raises(SystemExit, match="--sample-sheet.*--reference-root"):
        init_module.init(init_args(
            tmp_path / "cftk_init.json", reference_mode="local"
        ))


def test_interactive_missing_config_accepts_defaults(
    init_module, tmp_path, monkeypatch
):
    sheet = make_sample_sheet(tmp_path)
    reference_root = tmp_path / "references"
    make_profile(reference_root, checksums=True)
    config_path = tmp_path / "cftk_init.json"
    answers = iter(["", "", str(sheet), str(reference_root)])
    monkeypatch.setattr(init_module.sys, "stdin", SimpleNamespace(isatty=lambda: True))
    monkeypatch.setattr("builtins.input", lambda prompt: next(answers))

    init_module.init(init_args(
        config_path, non_interactive=False, reference_mode="local"
    ))

    compact = json.loads(config_path.read_text())
    assert compact["project_name"] == tmp_path.name
    assert compact["samples"] == "samples.tsv"
    assert compact["reference_profile"]["version"] == "1.0.0"
    assert (tmp_path / "cftk.lock.json").is_file()


def test_existing_config_is_not_overwritten(init_module, tmp_path, monkeypatch):
    config_path = tmp_path / "cftk_init.json"
    original = {"sentinel": "keep"}
    config_path.write_text(json.dumps(original))
    monkeypatch.setattr(init_module, "load_config", lambda path: {
        "project_name": "study",
        "output_dir": str(tmp_path),
        "comparison": "A_vs_B",
        "samples": {"A": [{"name": "a"}], "B": [{"name": "b"}]},
        "reference_data": {"genome_fa": "unused.fa"},
    })

    init_module.init(init_args(config_path))
    assert json.loads(config_path.read_text()) == original


def test_discovery_writes_template_without_guessing_roles(init_module, tmp_path):
    (tmp_path / "alpha_R1.fastq.gz").touch()
    (tmp_path / "alpha_R2.fastq.gz").touch()
    rows = init_module.discover_sample_inputs(tmp_path)
    sheet = init_module.write_sample_sheet_template(rows, tmp_path / "samples.tsv")

    text = Path(sheet).read_text()
    assert "alpha\t\t\tfastq" in text
    assert "control" not in text.lower()
    assert "case" not in text.lower()


def test_discovery_rejects_multilane_fastq(init_module, tmp_path):
    for lane in ("L001", "L002"):
        (tmp_path / f"alpha_{lane}_R1_001.fastq.gz").touch()
        (tmp_path / f"alpha_{lane}_R2_001.fastq.gz").touch()

    with pytest.raises(SystemExit, match="multi-lane or ambiguous"):
        init_module.discover_sample_inputs(tmp_path)


def test_valid_legacy_config_remains_unchanged(init_module, tmp_path):
    legacy = {
        "project_name": "legacy",
        "output_dir": str(tmp_path),
        "comparison": "Control_vs_Disease",
        "samples": {
            "Control": [{"name": "control", "input_type": "bam", "bam": "control.bam"}],
            "Disease": [{"name": "case", "input_type": "bam", "bam": "case.bam"}],
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
            "step4_methylation": {"params": {"min_depth": 10}},
        },
        "analysis": {},
    }
    config_path = tmp_path / "legacy.json"
    config_path.write_text(json.dumps(legacy))

    assert init_module.load_config(str(config_path)) == legacy


def test_explicit_roles_control_mesa_labels(monkeypatch, tmp_path):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import cftk

    cfg = {
        "comparison": "Alpha_vs_Beta",
        "group_roles": {"Alpha": "control", "Beta": "case"},
        "samples": {
            "Alpha": [{"name": "alpha_1"}],
            "Beta": [{"name": "beta_1"}],
        },
    }
    label_path = cftk._make_label(cfg, {"mesa": str(tmp_path / "mesa")})

    assert Path(label_path).read_text().splitlines() == ["alpha_1\t0", "beta_1\t1"]
