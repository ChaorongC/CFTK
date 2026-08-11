import hashlib
import json
from pathlib import Path

import pytest

import cftk_provenance as provenance


def _identity(**overrides):
    value = {
        "software_identity_schema_version": 1,
        "name": "cftk",
        "version": "1.0.0",
        "revision": "a" * 40,
        "source": "release",
        "dirty": False,
        "source_sha256": "b" * 64,
    }
    value.update(overrides)
    unsigned = {key: item for key, item in value.items() if key != "identity_sha256"}
    value["identity_sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return value


def test_hash_named_bytes_is_order_independent_and_excludes_build_identity():
    files = {
        "pyproject.toml": b"version = '1.0.0'\n",
        "src/cftk.py": b"print('cftk')\n",
        "src/cftk_provenance/build.json": b"changing build metadata\n",
    }

    first = provenance.hash_named_bytes(files)
    second = provenance.hash_named_bytes(dict(reversed(list(files.items()))))

    assert first == second
    assert first == provenance.hash_named_bytes({key: value for key, value in files.items() if "build.json" not in key})


def test_validate_identity_rejects_malformed_or_tampered_values():
    with pytest.raises(provenance.SoftwareIdentityError, match="fields are incomplete"):
        provenance.validate_software_identity({})

    value = _identity()
    value["source_sha256"] = "c" * 64
    with pytest.raises(provenance.SoftwareIdentityError, match="digest is inconsistent"):
        provenance.validate_software_identity(value)

    unknown_revision = _identity(revision="unknown", source="installed", dirty=None)
    provenance.validate_software_identity(unknown_revision)
    with pytest.raises(provenance.SoftwareIdentityError, match="revision-bound"):
        provenance.validate_software_identity(unknown_revision, require_revision=True)


def test_release_identity_is_prepared_and_deterministic():
    first = provenance.release_identity("1.0.0", "a" * 40, "b" * 64)
    second = provenance.release_identity("1.0.0", "a" * 40, "b" * 64)

    assert first == second
    assert first["source"] == "release"
    assert first["dirty"] is False
    assert provenance.validate_software_identity(first, require_revision=True) == first


def test_source_tree_identity_uses_git_revision_and_dirty_state(monkeypatch, tmp_path):
    monkeypatch.setattr(provenance, "_git_root", lambda start: tmp_path)
    monkeypatch.setattr(provenance, "source_tree_revision", lambda root: "a" * 40)
    monkeypatch.setattr(provenance, "source_tree_dirty", lambda root: True)
    monkeypatch.setattr(provenance, "source_tree_sha256", lambda root: "b" * 64)
    monkeypatch.setattr(provenance, "_read_embedded_identity", lambda: {"source": "unprepared"})
    monkeypatch.setattr(provenance, "_package_version", lambda: "1.0.0")
    provenance.get_software_identity.cache_clear()

    value = provenance.get_software_identity()

    assert value["revision"] == "a" * 40
    assert value["source"] == "git"
    assert value["dirty"] is True
    provenance.get_software_identity.cache_clear()


def test_source_tree_identity_uses_its_own_project_version(monkeypatch, tmp_path):
    (tmp_path / "pyproject.toml").write_text(
        "[project]\nversion = \"2.4.6\"\n", encoding="utf-8"
    )
    monkeypatch.setattr(provenance, "_git_root", lambda start: tmp_path)
    monkeypatch.setattr(provenance, "source_tree_revision", lambda root: "a" * 40)
    monkeypatch.setattr(provenance, "source_tree_dirty", lambda root: False)
    monkeypatch.setattr(provenance, "source_tree_sha256", lambda root: "b" * 64)
    monkeypatch.setattr(provenance, "_read_embedded_identity", lambda: {"source": "unprepared"})
    monkeypatch.setattr(provenance, "_package_version", lambda: "0.0.1")
    provenance.get_software_identity.cache_clear()

    value = provenance.get_software_identity()

    assert value["version"] == "2.4.6"
    provenance.get_software_identity.cache_clear()


def test_installed_fallback_identity_has_no_path_and_unknown_revision(monkeypatch, tmp_path):
    monkeypatch.setattr(provenance, "_git_root", lambda start: None)
    monkeypatch.setattr(provenance, "_read_embedded_identity", lambda: {"source": "unprepared"})
    monkeypatch.setattr(provenance, "_package_version", lambda: "1.0.0")
    monkeypatch.setattr(provenance, "_installed_payload_sha256", lambda: "b" * 64)
    provenance.get_software_identity.cache_clear()

    value = provenance.get_software_identity()

    assert value["source"] == "installed"
    assert value["revision"] == "unknown"
    assert value["dirty"] is None
    assert not any(str(tmp_path) in str(item) for item in value.values())
    provenance.get_software_identity.cache_clear()


def test_prepared_release_identity_is_loaded_from_embedded_build(monkeypatch):
    prepared = _identity()
    monkeypatch.setattr(provenance, "_read_embedded_identity", lambda: prepared)
    monkeypatch.setattr(provenance, "_package_version", lambda: "1.0.0")
    monkeypatch.setattr(provenance, "_git_root", lambda start: None)
    provenance.get_software_identity.cache_clear()

    assert provenance.get_software_identity() == prepared
    provenance.get_software_identity.cache_clear()
