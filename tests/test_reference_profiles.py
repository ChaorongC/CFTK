import hashlib
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def profiles_module(monkeypatch):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "src"))
    import reference_profiles

    return reference_profiles


def make_profile(root, version="1.0.0", checksums=False):
    profile_dir = root / "twist_human_methylome_hg38" / version
    (profile_dir / "genome").mkdir(parents=True)
    (profile_dir / "assay").mkdir()
    genome = profile_dir / "genome/hg38.fa"
    genome.write_text(">chr1\nACGT\n")
    (profile_dir / "genome/hg38.2bit").write_bytes(b"2bit")
    (profile_dir / "genome/hg38.chrom.sizes").write_text("chr1\t4\n")
    (profile_dir / "assay/covered_targets.bed").write_text("chr1\t0\t4\n")
    components = {
        "genome_fa": {"path": "genome/hg38.fa"},
        "genome_2bit": {"path": "genome/hg38.2bit"},
        "chrom_sizes": {"path": "genome/hg38.chrom.sizes"},
        "target_bed": {"path": "assay/covered_targets.bed"},
    }
    if checksums:
        for value in components.values():
            path = profile_dir / value["path"]
            value["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "manifest_version": 1,
        "profile_id": "twist_human_methylome_hg38",
        "version": version,
        "assay": "twist_human_methylome",
        "genome": "hg38",
        "components": components,
    }
    (profile_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return profile_dir


def test_local_profile_resolves_components_and_verifies_checksums(
    profiles_module, tmp_path
):
    profile_dir = make_profile(tmp_path, checksums=True)
    profile = profiles_module.load_reference_profile(
        tmp_path, "twist_human_methylome_hg38", "1.0.0",
        verify_checksums=True,
    )

    assert profile["profile_dir"] == str(profile_dir.resolve())
    assert profile["components"]["genome_fa"] == str(
        (profile_dir / "genome/hg38.fa").resolve()
    )
    assert profile["components"]["target_bed"].endswith("covered_targets.bed")
    assert len(profile["component_hashes"]["genome_fa"]) == 64


def test_profile_rejects_checksum_mismatch(profiles_module, tmp_path):
    profile_dir = make_profile(tmp_path, checksums=True)
    (profile_dir / "genome/hg38.fa").write_text(">chr1\nAAAA\n")

    with pytest.raises(SystemExit, match="checksum mismatch.*genome_fa"):
        profiles_module.load_reference_profile(
            tmp_path, "twist_human_methylome_hg38", "1.0.0",
            verify_checksums=True,
        )


def test_profile_defers_component_hashing_when_verification_is_disabled(
    profiles_module, tmp_path
):
    profile_dir = make_profile(tmp_path, checksums=True)
    (profile_dir / "genome/hg38.fa").write_text(">chr1\nAAAA\n")

    profile = profiles_module.load_reference_profile(
        tmp_path, "twist_human_methylome_hg38", "1.0.0",
        verify_checksums=False,
    )

    assert len(profile["component_hashes"]["genome_fa"]) == 64


def test_profile_rejects_path_outside_profile(profiles_module, tmp_path):
    profile_dir = make_profile(tmp_path)
    manifest_path = profile_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["components"]["genome_fa"]["path"] = "../../../../outside.fa"
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(SystemExit, match="escapes profile directory"):
        profiles_module.load_reference_profile(
            tmp_path, "twist_human_methylome_hg38", "1.0.0"
        )


def test_profile_version_must_be_unambiguous(profiles_module, tmp_path):
    make_profile(tmp_path, version="1.0.0")
    make_profile(tmp_path, version="2.0.0")

    with pytest.raises(SystemExit, match="multiple versions"):
        profiles_module.load_reference_profile(
            tmp_path, "twist_human_methylome_hg38"
        )


def test_managed_mode_is_disabled_without_immutable_registry(
    profiles_module, tmp_path
):
    with pytest.raises(SystemExit, match="immutable.*managed.*available"):
        profiles_module.acquire_reference_profile(
            mode="managed",
            reference_root=tmp_path,
            profile_id="twist_human_methylome_hg38",
            version="1.0.0",
        )


@pytest.mark.parametrize(
    "bed_text, message",
    [
        ("chr2\t0\t4\n", "target BED.*unknown contig.*chr2"),
        ("chr1\t0\t5\n", "target BED.*exceeds.*chr1"),
        ("chr1\t3\t2\n", "target BED.*invalid interval"),
    ],
)
def test_profile_rejects_incompatible_target_bed(
    profiles_module, tmp_path, bed_text, message
):
    profile_dir = make_profile(tmp_path)
    (profile_dir / "assay/covered_targets.bed").write_text(bed_text)

    with pytest.raises(SystemExit, match=message):
        profiles_module.load_reference_profile(
            tmp_path, "twist_human_methylome_hg38", "1.0.0",
            validate_compatibility=True,
        )
