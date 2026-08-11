import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

from scripts.release import prepare_build


def _source_root(tmp_path):
    root = tmp_path / "cftk"
    package = root / "src/cftk_provenance"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("# fixture\n", encoding="utf-8")
    (root / "pyproject.toml").write_text(
        '[project]\nname = "cftk"\nversion = "1.0.0"\n', encoding="utf-8"
    )
    return root


def test_prepare_build_requires_exact_version_tag(tmp_path):
    root = _source_root(tmp_path)

    with pytest.raises(SystemExit, match="does not match project version"):
        prepare_build.prepare(
            root, revision="a" * 40, tag="v1.0.1",
            output=Path("src/cftk_provenance/build.json"),
        )


def test_prepare_build_writes_privacy_safe_identity_inside_source_root(
    tmp_path, monkeypatch
):
    root = _source_root(tmp_path)
    identity = {
        "software_identity_schema_version": 1,
        "name": "cftk",
        "version": "1.0.0",
        "revision": "a" * 40,
        "source": "release",
        "dirty": False,
        "source_sha256": "b" * 64,
        "identity_sha256": "c" * 64,
    }
    fake = SimpleNamespace(
        release_identity=lambda version, revision, source_sha256: dict(identity),
        source_tree_sha256=lambda source_root: "b" * 64,
    )
    monkeypatch.setitem(sys.modules, "cftk_provenance", fake)

    result = prepare_build.prepare(
        root, revision="a" * 40, tag="v1.0.0",
        output=Path("src/cftk_provenance/build.json"),
    )

    output = root / "src/cftk_provenance/build.json"
    assert result == identity
    assert json.loads(output.read_text(encoding="utf-8")) == identity
    assert str(tmp_path) not in output.read_text(encoding="utf-8")


def test_prepare_build_rejects_output_outside_source_root(tmp_path, monkeypatch):
    root = _source_root(tmp_path)
    fake = SimpleNamespace(
        release_identity=lambda *args: {"source": "release"},
        source_tree_sha256=lambda source_root: "b" * 64,
    )
    monkeypatch.setitem(sys.modules, "cftk_provenance", fake)

    with pytest.raises(SystemExit, match="must remain inside"):
        prepare_build.prepare(
            root, revision="a" * 40, tag="v1.0.0",
            output=tmp_path / "outside.json",
        )
