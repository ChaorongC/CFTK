import io
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
