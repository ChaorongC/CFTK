"""Resolve a privacy-safe identity for the running CFTK software payload."""

from __future__ import annotations

from functools import lru_cache
import hashlib
from importlib import metadata, resources
import json
from pathlib import Path, PurePosixPath
import re
import subprocess


SOFTWARE_IDENTITY_SCHEMA_VERSION = 1
REVISION_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
SOURCE_SUFFIXES = {".html", ".json", ".py", ".r", ".toml"}
EMBEDDED_IDENTITY_NAME = "build.json"
VERSION_PATTERN = re.compile(r"(?m)^version\s*=\s*[\"']([^\"']+)[\"']\s*$")


class SoftwareIdentityError(RuntimeError):
    """The installed or source-tree software identity cannot be trusted."""


def _eligible_source_path(path: PurePosixPath) -> bool:
    if path.as_posix() == "src/cftk_provenance/build.json":
        return False
    if any(
        part in {"__pycache__", "build"}
        or part.endswith((".dist-info", ".egg-info"))
        for part in path.parts
    ):
        return False
    return path.suffix.lower() in SOURCE_SUFFIXES


def hash_named_bytes(files: dict[str, bytes]) -> str:
    """Hash named source payloads independent of filesystem traversal order."""
    digest = hashlib.sha256()
    eligible = [
        (PurePosixPath(name), content)
        for name, content in files.items()
        if _eligible_source_path(PurePosixPath(name))
    ]
    if not eligible:
        raise SoftwareIdentityError("no release-relevant source files were found")
    for path, content in sorted(eligible, key=lambda item: item[0].as_posix()):
        encoded = path.as_posix().encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


def _git_output(root: Path, *arguments: str) -> bytes:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), *arguments],
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise SoftwareIdentityError("could not inspect the CFTK Git checkout") from exc
    return result.stdout


def _git_root(start: Path) -> Path | None:
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate.resolve()
    return None


def source_tree_files(root: Path) -> dict[str, bytes]:
    """Return tracked and untracked release-relevant files from a Git checkout."""
    root = Path(root).resolve()
    raw = _git_output(
        root,
        "ls-files",
        "-z",
        "--cached",
        "--others",
        "--exclude-standard",
        "--",
        "src",
        "pyproject.toml",
    )
    files = {}
    for value in raw.split(b"\0"):
        if not value:
            continue
        name = value.decode("utf-8")
        relative = PurePosixPath(name)
        if not _eligible_source_path(relative):
            continue
        path = root / Path(*relative.parts)
        if path.is_file():
            files[relative.as_posix()] = path.read_bytes()
    return files


def source_tree_sha256(root: Path) -> str:
    return hash_named_bytes(source_tree_files(root))


def source_tree_dirty(root: Path) -> bool:
    output = _git_output(
        Path(root).resolve(),
        "status",
        "--porcelain",
        "--untracked-files=normal",
        "--",
        "src",
        "pyproject.toml",
    )
    return bool(output.strip())


def source_tree_revision(root: Path) -> str:
    revision = _git_output(Path(root).resolve(), "rev-parse", "HEAD").decode().strip()
    if not REVISION_PATTERN.fullmatch(revision):
        raise SoftwareIdentityError("Git returned an invalid CFTK revision")
    return revision


def _package_version() -> str:
    try:
        return metadata.version("cftk")
    except metadata.PackageNotFoundError:
        embedded = _read_embedded_identity()
        version = embedded.get("version")
        return version if isinstance(version, str) and version else "development"


def _source_tree_version(root: Path) -> str:
    try:
        text = (Path(root) / "pyproject.toml").read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return _package_version()
    match = VERSION_PATTERN.search(text)
    return match.group(1) if match else _package_version()


def _read_embedded_identity() -> dict:
    try:
        text = resources.files("cftk_provenance").joinpath(
            EMBEDDED_IDENTITY_NAME
        ).read_text(encoding="utf-8")
        value = json.loads(text)
    except (FileNotFoundError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise SoftwareIdentityError("embedded CFTK build identity is unreadable") from exc
    if not isinstance(value, dict):
        raise SoftwareIdentityError("embedded CFTK build identity is not an object")
    return value


def _installed_payload_sha256() -> str:
    try:
        distribution = metadata.distribution("cftk")
    except metadata.PackageNotFoundError:
        source_root = Path(__file__).resolve().parents[1]
        files = {
            f"src/{path.relative_to(source_root).as_posix()}": path.read_bytes()
            for path in source_root.rglob("*")
            if path.is_file()
            and _eligible_source_path(
                PurePosixPath("src") / PurePosixPath(path.relative_to(source_root).as_posix())
            )
        }
        pyproject = source_root.parent / "pyproject.toml"
        if pyproject.is_file():
            files["pyproject.toml"] = pyproject.read_bytes()
        return hash_named_bytes(files)

    files = {}
    for entry in distribution.files or ():
        relative = PurePosixPath(str(entry))
        if not _eligible_source_path(relative):
            continue
        path = Path(distribution.locate_file(entry))
        if path.is_file():
            files[relative.as_posix()] = path.read_bytes()
    return hash_named_bytes(files)


def _identity_sha256(identity: dict) -> str:
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def validate_software_identity(identity: dict, *, require_revision: bool = False) -> dict:
    """Validate and return a normalized software identity without local paths."""
    if not isinstance(identity, dict):
        raise SoftwareIdentityError("software identity must be an object")
    required = {
        "software_identity_schema_version",
        "name",
        "version",
        "revision",
        "source",
        "dirty",
        "source_sha256",
        "identity_sha256",
    }
    if set(identity) != required:
        raise SoftwareIdentityError("software identity fields are incomplete")
    if identity["software_identity_schema_version"] != SOFTWARE_IDENTITY_SCHEMA_VERSION:
        raise SoftwareIdentityError("software identity schema is unsupported")
    if identity["name"] != "cftk":
        raise SoftwareIdentityError("software identity package name is invalid")
    if not isinstance(identity["version"], str) or not identity["version"]:
        raise SoftwareIdentityError("software identity version is invalid")
    revision = identity["revision"]
    if revision != "unknown" and (
        not isinstance(revision, str) or not REVISION_PATTERN.fullmatch(revision)
    ):
        raise SoftwareIdentityError("software identity revision is invalid")
    if require_revision and revision == "unknown":
        raise SoftwareIdentityError("software identity is not revision-bound")
    if identity["source"] not in {"git", "installed", "release", "source-tree"}:
        raise SoftwareIdentityError("software identity source is invalid")
    if identity["dirty"] not in {True, False, None}:
        raise SoftwareIdentityError("software identity dirty state is invalid")
    if not isinstance(identity["source_sha256"], str) or not SHA256_PATTERN.fullmatch(
        identity["source_sha256"]
    ):
        raise SoftwareIdentityError("software source SHA-256 is invalid")
    unsigned = {key: value for key, value in identity.items() if key != "identity_sha256"}
    if identity["identity_sha256"] != _identity_sha256(unsigned):
        raise SoftwareIdentityError("software identity digest is inconsistent")
    return dict(identity)


def release_identity(version: str, revision: str, source_sha256: str) -> dict:
    unsigned = {
        "software_identity_schema_version": SOFTWARE_IDENTITY_SCHEMA_VERSION,
        "name": "cftk",
        "version": version,
        "revision": revision,
        "source": "release",
        "dirty": False,
        "source_sha256": source_sha256,
    }
    return {**unsigned, "identity_sha256": _identity_sha256(unsigned)}


@lru_cache(maxsize=1)
def get_software_identity() -> dict:
    """Resolve the exact running build without exposing its filesystem path."""
    embedded = _read_embedded_identity()
    if embedded.get("source") == "release":
        version = _package_version()
        identity = release_identity(
            embedded.get("version"),
            embedded.get("revision"),
            embedded.get("source_sha256"),
        )
        if identity["version"] != version:
            raise SoftwareIdentityError(
                "installed CFTK version does not match its embedded release identity"
            )
        return validate_software_identity(identity, require_revision=True)

    root = _git_root(Path(__file__).resolve())
    if root is not None:
        version = _source_tree_version(root)
        unsigned = {
            "software_identity_schema_version": SOFTWARE_IDENTITY_SCHEMA_VERSION,
            "name": "cftk",
            "version": version,
            "revision": source_tree_revision(root),
            "source": "git",
            "dirty": source_tree_dirty(root),
            "source_sha256": source_tree_sha256(root),
        }
    else:
        version = _package_version()
        unsigned = {
            "software_identity_schema_version": SOFTWARE_IDENTITY_SCHEMA_VERSION,
            "name": "cftk",
            "version": version,
            "revision": "unknown",
            "source": "installed" if version != "development" else "source-tree",
            "dirty": None,
            "source_sha256": _installed_payload_sha256(),
        }
    return validate_software_identity({
        **unsigned,
        "identity_sha256": _identity_sha256(unsigned),
    })


def format_software_identity(identity: dict, *, short_revision: bool = True) -> str:
    value = validate_software_identity(identity)
    revision = value["revision"]
    if short_revision and revision != "unknown":
        revision = revision[:12]
    dirty = ", dirty" if value["dirty"] else ""
    return f"CFTK {value['version']} ({revision}, {value['source']}{dirty})"
