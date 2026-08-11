import io
import hashlib
import json
import tarfile
import zipfile

import pytest

from scripts.validation import check_distribution


def _required_members(prefix=""):
    return {
        f"{prefix}{name}": b"public fixture\n"
        for name in check_distribution.REQUIRED_PATHS
    }


def _write_wheel(path, members):
    with zipfile.ZipFile(path, "w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)


def _write_sdist(path, members):
    with tarfile.open(path, "w:gz") as archive:
        for name, content in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))


def _release_identity(revision="a" * 40):
    value = {
        "software_identity_schema_version": 1,
        "name": "cftk",
        "version": "1.0.0",
        "revision": revision,
        "source": "release",
        "dirty": False,
        "source_sha256": "b" * 64,
    }
    value["identity_sha256"] = hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return value


def test_validate_archive_accepts_required_wheel_content(tmp_path):
    wheel = tmp_path / "cftk-1.0.0-py3-none-any.whl"
    _write_wheel(wheel, _required_members())

    check_distribution.validate_archive(wheel)


def test_validate_archive_accepts_normal_sdist_metadata(tmp_path):
    sdist = tmp_path / "cftk-1.0.0.tar.gz"
    members = _required_members("cftk-1.0.0/src/")
    for name in check_distribution.SDIST_REQUIRED_PATHS:
        members[f"cftk-1.0.0/{name}"] = b"public fixture\n"
    members["cftk-1.0.0/src/cftk.egg-info/PKG-INFO"] = b"metadata\n"
    _write_sdist(sdist, members)

    check_distribution.validate_archive(sdist)


@pytest.mark.parametrize(
    "name",
    [
        "data/reference.npy",
        "src/test.ipynb",
        "analysis/__pycache__/module.pyc",
        "build/lib/cftk.py",
    ],
)
def test_validate_archive_rejects_forbidden_content(tmp_path, name):
    wheel = tmp_path / "cftk-1.0.0-py3-none-any.whl"
    members = _required_members()
    members[name] = b"private fixture\n"
    _write_wheel(wheel, members)

    with pytest.raises(check_distribution.ValidationError, match="forbidden archive paths"):
        check_distribution.validate_archive(wheel)


@pytest.mark.parametrize(
    ("name", "content"),
    [
        ("report/example.txt", b"input=/" + b"dfs99/research/project/file.bam\n"),
        ("report/example.txt", b"sample=" + b"patient" + b"_042\n"),
    ],
)
def test_validate_archive_rejects_generic_private_patterns(tmp_path, name, content):
    wheel = tmp_path / "cftk-1.0.0-py3-none-any.whl"
    members = _required_members()
    members[name] = content
    _write_wheel(wheel, members)

    with pytest.raises(check_distribution.ValidationError, match="possible private content"):
        check_distribution.validate_archive(wheel)


def test_privacy_scan_ignores_embedded_base64_data():
    encoded = (
        b"<svg><image href=\"data:image/jpeg;base64,AAAA/"
        + b"dfs99BBBB=\"/></svg>"
    )

    assert check_distribution._privacy_violations({"diagram.svg": encoded}) == []


def test_validate_release_identity_requires_matching_clean_archives(tmp_path):
    identity = _release_identity()
    wheel = tmp_path / "cftk-1.0.0-py3-none-any.whl"
    sdist = tmp_path / "cftk-1.0.0.tar.gz"
    wheel_members = _required_members()
    wheel_members["cftk_provenance/build.json"] = (
        json.dumps(identity, sort_keys=True).encode()
    )
    _write_wheel(wheel, wheel_members)
    sdist_members = _required_members("cftk-1.0.0/src/")
    for name in check_distribution.SDIST_REQUIRED_PATHS:
        sdist_members[f"cftk-1.0.0/{name}"] = b"public fixture\n"
    sdist_members["cftk-1.0.0/src/cftk_provenance/build.json"] = (
        json.dumps(identity, sort_keys=True).encode()
    )
    _write_sdist(sdist, sdist_members)

    result = check_distribution.validate_release_archives(
        [wheel, sdist], version="1.0.0", revision="a" * 40,
        source_sha256="b" * 64,
    )

    assert result == identity


def test_validate_release_identity_rejects_archive_mismatch(tmp_path):
    first = _release_identity("a" * 40)
    second = _release_identity("c" * 40)
    wheel = tmp_path / "cftk-1.0.0-py3-none-any.whl"
    members = _required_members()
    members["cftk_provenance/build.json"] = json.dumps(first).encode()
    _write_wheel(wheel, members)
    other = tmp_path / "cftk-1.0.0-other.whl"
    other_members = _required_members()
    other_members["cftk_provenance/build.json"] = json.dumps(second).encode()
    _write_wheel(other, other_members)

    with pytest.raises(check_distribution.ValidationError, match="different build identities"):
        check_distribution.validate_release_archives([wheel, other])
