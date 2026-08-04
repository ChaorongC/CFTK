"""Versioned reference-profile discovery, validation, and acquisition."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
import threading
from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

try:
    import fcntl
except ImportError:  # pragma: no cover - Unix HPC systems provide fcntl.
    fcntl = None


DEFAULT_PROFILE_ID = "twist_human_methylome_hg38"
DEFAULT_ASSAY = "twist_human_methylome"
DEFAULT_GENOME = "hg38"
REQUIRED_COMPONENTS = {
    "genome_fa",
    "genome_2bit",
    "chrom_sizes",
    "target_bed",
}
_SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_SHA256 = re.compile(r"^[0-9a-fA-F]{64}$")
_LOCKS = {}
_LOCKS_GUARD = threading.Lock()
_DOWNLOAD_CHUNK_SIZE = 1024 * 1024
_DOWNLOAD_ATTEMPTS = 2


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_identifier(value, label):
    value = str(value or "")
    if not value or not _SAFE_ID.fullmatch(value):
        sys.exit(
            f"[references] ERROR: invalid {label} {value!r}; use only letters, "
            "digits, dot, underscore, or hyphen."
        )
    return value


def _registry_error(message):
    sys.exit(f"[references] ERROR: invalid managed registry: {message}")


def _read_registry(registry=None):
    if registry is None:
        override = os.environ.get("CFTK_REFERENCE_REGISTRY")
        if override:
            registry = Path(override).expanduser()
        else:
            try:
                text = resources.files("cftk_registry").joinpath("registry.json").read_text()
                registry = json.loads(text)
            except (OSError, json.JSONDecodeError, ModuleNotFoundError) as exc:
                sys.exit(f"[references] ERROR: could not load packaged registry: {exc}")
    if isinstance(registry, (str, os.PathLike)):
        registry_path = Path(registry).expanduser().resolve()
        try:
            registry = json.loads(registry_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            sys.exit(f"[references] ERROR: could not load registry {registry_path}: {exc}")
    if not isinstance(registry, dict):
        _registry_error("the registry must be a JSON object.")
    return registry


def _validate_metadata(value, component, field):
    if not isinstance(value, dict):
        _registry_error(f"component '{component}' requires {field} metadata.")
    for key in ("name", "url"):
        if not isinstance(value.get(key), str) or not value[key].strip():
            _registry_error(
                f"component '{component}' {field}.{key} must be a nonempty string."
            )
    parsed = urlsplit(value["url"])
    if parsed.scheme != "https" or not parsed.netloc:
        _registry_error(
            f"component '{component}' {field}.url must be an HTTPS URL."
        )


def _validate_artifact_url(url, component):
    parsed = urlsplit(url)
    if (
        parsed.scheme != "https"
        or not parsed.netloc
        or parsed.username is not None
        or parsed.password is not None
        or parsed.query
        or parsed.fragment
    ):
        _registry_error(
            f"component '{component}' artifact requires an immutable HTTPS URL."
        )


def _validate_size(value, component, field):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        _registry_error(f"component '{component}' {field} must be a positive integer.")


def _validate_hash(value, component, field):
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        _registry_error(f"component '{component}' {field} must be a SHA-256 value.")


def _validate_component_path(value, component):
    if not isinstance(value, str) or not value:
        _registry_error(f"component '{component}' path must be a relative path.")
    path = Path(value)
    if path.is_absolute() or any(part in ("", ".", "..") for part in path.parts):
        _registry_error(f"component '{component}' path must be a safe relative path.")


def validate_reference_registry(registry=None):
    """Validate and return an immutable managed-reference registry."""
    registry = _read_registry(registry)
    if registry.get("registry_version") != 1:
        _registry_error("registry_version must be 1.")
    profiles = registry.get("profiles")
    if not isinstance(profiles, dict):
        _registry_error("profiles must be an object.")

    for profile_id, versions in profiles.items():
        _safe_identifier(profile_id, "registry profile id")
        if not isinstance(versions, dict) or not versions:
            _registry_error(f"profile '{profile_id}' must contain at least one version.")
        for version, profile in versions.items():
            _safe_identifier(version, "registry profile version")
            if not isinstance(profile, dict):
                _registry_error(f"profile '{profile_id}' version '{version}' must be an object.")
            expected = {"profile_id": profile_id, "version": version}
            for key, value in expected.items():
                if profile.get(key) != value:
                    _registry_error(
                        f"profile '{profile_id}' version '{version}' {key} must be {value!r}."
                    )
            for key in ("assay", "genome"):
                if not isinstance(profile.get(key), str) or not profile[key]:
                    _registry_error(
                        f"profile '{profile_id}' version '{version}' requires {key}."
                    )
            components = profile.get("components")
            if not isinstance(components, dict):
                _registry_error(
                    f"profile '{profile_id}' version '{version}' components must be an object."
                )
            missing = sorted(REQUIRED_COMPONENTS - components.keys())
            if missing:
                _registry_error(
                    f"profile '{profile_id}' version '{version}' is missing components: {missing}."
                )
            seen_paths = set()
            for name, component in components.items():
                _safe_identifier(name, "registry component name")
                if not isinstance(component, dict):
                    _registry_error(f"component '{name}' must be an object.")
                _validate_component_path(component.get("path"), name)
                if component["path"] in seen_paths:
                    _registry_error(f"component '{name}' reuses path {component['path']!r}.")
                seen_paths.add(component["path"])
                _validate_size(component.get("size"), name, "size")
                _validate_hash(component.get("sha256"), name, "sha256")
                _validate_metadata(component.get("license"), name, "license")
                _validate_metadata(component.get("source"), name, "source")
                artifact = component.get("artifact")
                if not isinstance(artifact, dict):
                    _registry_error(f"component '{name}' requires artifact metadata.")
                urls = artifact.get("urls")
                if (
                    artifact.get("immutable") is not True
                    or not isinstance(urls, list)
                    or not urls
                    or any(not isinstance(url, str) for url in urls)
                ):
                    _registry_error(
                        f"component '{name}' artifact requires immutable HTTPS URLs."
                    )
                for url in urls:
                    _validate_artifact_url(url, name)
                _validate_size(artifact.get("size"), name, "artifact.size")
                _validate_hash(artifact.get("sha256"), name, "artifact.sha256")
                compression = artifact.get("compression", "none")
                if compression not in ("none", "gzip"):
                    _registry_error(
                        f"component '{name}' artifact.compression must be 'none' or 'gzip'."
                    )
                transform = artifact.get("transform", "none")
                if transform not in ("none", "fai_to_chrom_sizes"):
                    _registry_error(
                        f"component '{name}' artifact.transform must be 'none' "
                        "or 'fai_to_chrom_sizes'."
                    )
                if compression == "none" and transform == "none" and (
                    artifact["size"] != component["size"]
                    or artifact["sha256"].lower() != component["sha256"].lower()
                ):
                    _registry_error(
                        f"component '{name}' uncompressed artifact and installed file metadata differ."
                    )
    return registry


def _registry_entry(registry, profile_id, version):
    profiles = registry["profiles"]
    versions = profiles.get(profile_id)
    if not versions:
        sys.exit(
            f"[references] ERROR: managed registry does not contain profile '{profile_id}'."
        )
    if version is None:
        available = sorted(versions)
        if len(available) != 1:
            sys.exit(
                f"[references] ERROR: managed profile '{profile_id}' has versions "
                f"{available}; specify --profile-version."
            )
        version = available[0]
    version = _safe_identifier(version, "profile version")
    if version not in versions:
        sys.exit(
            f"[references] ERROR: managed registry does not contain profile "
            f"'{profile_id}' version '{version}'."
        )
    return version, versions[version]


def managed_profile_available(profile_id=DEFAULT_PROFILE_ID, *, registry=None):
    """Return whether a validated registry contains the requested profile."""
    registry = validate_reference_registry(registry)
    return profile_id in registry["profiles"]


def _canonical_sha256(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


@contextmanager
def _install_lock(lock_path):
    lock_key = str(lock_path.resolve())
    with _LOCKS_GUARD:
        thread_lock = _LOCKS.setdefault(lock_key, threading.Lock())
    with thread_lock:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o600)
        try:
            if fcntl is not None:
                fcntl.flock(descriptor, fcntl.LOCK_EX)
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            os.close(descriptor)


def _download_artifact(component_name, spec, destination, opener):
    artifact = spec["artifact"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    last_error = None
    for url in artifact["urls"]:
        for _ in range(_DOWNLOAD_ATTEMPTS):
            part = destination.with_suffix(destination.suffix + ".part")
            part.unlink(missing_ok=True)
            try:
                request = Request(url, headers={"User-Agent": "CFTK reference downloader/1"})
                with opener(request, timeout=60) as response, part.open("wb") as handle:
                    final_url = response.geturl() if hasattr(response, "geturl") else url
                    _validate_artifact_url(final_url, component_name)
                    digest = hashlib.sha256()
                    size = 0
                    while True:
                        block = response.read(_DOWNLOAD_CHUNK_SIZE)
                        if not block:
                            break
                        size += len(block)
                        if size > artifact["size"]:
                            sys.exit(
                                f"[references] ERROR: download size mismatch for component "
                                f"'{component_name}': received more than {artifact['size']} bytes."
                            )
                        digest.update(block)
                        handle.write(block)
                if size != artifact["size"]:
                    sys.exit(
                        f"[references] ERROR: download size mismatch for component "
                        f"'{component_name}': expected {artifact['size']}, got {size}."
                    )
                if digest.hexdigest() != artifact["sha256"].lower():
                    sys.exit(
                        f"[references] ERROR: download checksum mismatch for component "
                        f"'{component_name}'."
                    )
                os.replace(part, destination)
                return
            except (HTTPError, URLError, OSError) as exc:
                last_error = exc
                part.unlink(missing_ok=True)
    sys.exit(
        f"[references] ERROR: could not download component '{component_name}': "
        f"{last_error}"
    )


def _materialize_component(component_name, spec, artifact_path, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    part = output_path.with_suffix(output_path.suffix + ".part")
    part.unlink(missing_ok=True)
    try:
        if spec["artifact"].get("compression", "none") == "gzip":
            source = gzip.open(artifact_path, "rb")
        else:
            source = artifact_path.open("rb")
        digest = hashlib.sha256()
        size = 0
        with source, part.open("wb") as destination:
            transform = spec["artifact"].get("transform", "none")
            if transform == "fai_to_chrom_sizes":
                blocks = _project_fai_to_chrom_sizes(source, component_name)
            else:
                blocks = iter(lambda: source.read(_DOWNLOAD_CHUNK_SIZE), b"")
            for block in blocks:
                size += len(block)
                if size > spec["size"]:
                    sys.exit(
                        f"[references] ERROR: installed size mismatch for component "
                        f"'{component_name}'."
                    )
                digest.update(block)
                destination.write(block)
        if size != spec["size"]:
            sys.exit(
                f"[references] ERROR: installed size mismatch for component "
                f"'{component_name}': expected {spec['size']}, got {size}."
            )
        if digest.hexdigest() != spec["sha256"].lower():
            sys.exit(
                f"[references] ERROR: installed checksum mismatch for component "
                f"'{component_name}'."
            )
        os.replace(part, output_path)
    except (gzip.BadGzipFile, OSError) as exc:
        part.unlink(missing_ok=True)
        sys.exit(
            f"[references] ERROR: could not materialize component '{component_name}': {exc}"
        )


def _project_fai_to_chrom_sizes(source, component_name):
    for line_number, raw_line in enumerate(source, start=1):
        fields = raw_line.rstrip(b"\r\n").split(b"\t")
        if len(fields) < 2 or not fields[0] or not fields[1]:
            sys.exit(
                f"[references] ERROR: invalid FAI line {line_number} for "
                f"component '{component_name}'."
            )
        try:
            sequence_length = int(fields[1])
        except ValueError:
            sys.exit(
                f"[references] ERROR: invalid FAI sequence length on line "
                f"{line_number} for component '{component_name}'."
            )
        if sequence_length < 1:
            sys.exit(
                f"[references] ERROR: nonpositive FAI sequence length on line "
                f"{line_number} for component '{component_name}'."
            )
        yield fields[0] + b"\t" + str(sequence_length).encode("ascii") + b"\n"


def _managed_manifest(entry, registry_entry_sha256):
    components = {}
    for name, spec in entry["components"].items():
        components[name] = {
            "path": spec["path"],
            "size": spec["size"],
            "sha256": spec["sha256"].lower(),
            "artifact": spec["artifact"],
            "license": spec["license"],
            "source": spec["source"],
        }
    return {
        "manifest_version": 1,
        "profile_id": entry["profile_id"],
        "version": entry["version"],
        "assay": entry["assay"],
        "genome": entry["genome"],
        "acquisition": {
            "mode": "managed",
            "registry_version": 1,
            "registry_entry_sha256": registry_entry_sha256,
        },
        "components": components,
    }


def _acquire_managed_profile(
    reference_root,
    profile_id,
    version,
    *,
    registry=None,
    opener=None,
    verify_checksums=False,
    validate_compatibility=False,
):
    registry = validate_reference_registry(registry)
    profile_id = _safe_identifier(profile_id, "profile id")
    version, entry = _registry_entry(registry, profile_id, version)
    entry_hash = _canonical_sha256(entry)
    root = Path(reference_root).expanduser().resolve()
    profile_root = root / profile_id
    profile_dir = profile_root / version
    lock_path = profile_root / f".{version}.install.lock"
    root.mkdir(parents=True, exist_ok=True)
    opener = opener or urlopen

    with _install_lock(lock_path):
        if profile_dir.exists():
            try:
                installed = load_reference_profile(
                    root,
                    profile_id,
                    version,
                    verify_checksums=verify_checksums,
                    validate_compatibility=validate_compatibility,
                )
            except SystemExit as exc:
                sys.exit(
                    f"[references] ERROR: installed immutable profile '{profile_id}' "
                    f"version '{version}' is invalid; refusing replacement. {exc}"
                )
            if installed.get("acquisition", {}).get("registry_entry_sha256") != entry_hash:
                sys.exit(
                    f"[references] ERROR: installed immutable profile '{profile_id}' "
                    f"version '{version}' differs from the registry; refusing replacement."
                )
            return installed

        stage_root = Path(tempfile.mkdtemp(prefix=".cftk-install-", dir=root))
        stage_profile = stage_root / profile_id / version
        downloads = stage_root / ".downloads"
        try:
            stage_profile.mkdir(parents=True)
            downloads.mkdir()
            for name, spec in entry["components"].items():
                artifact_path = downloads / name
                _download_artifact(name, spec, artifact_path, opener)
                _materialize_component(
                    name, spec, artifact_path, stage_profile / spec["path"]
                )
            manifest = _managed_manifest(entry, entry_hash)
            (stage_profile / "manifest.json").write_text(
                json.dumps(manifest, indent=2) + "\n"
            )
            load_reference_profile(
                stage_root,
                profile_id,
                version,
                verify_checksums=True,
                validate_compatibility=True,
            )
            profile_root.mkdir(parents=True, exist_ok=True)
            os.replace(stage_profile, profile_dir)
        finally:
            shutil.rmtree(stage_root, ignore_errors=True)

        return load_reference_profile(
            root,
            profile_id,
            version,
            verify_checksums=verify_checksums,
            validate_compatibility=validate_compatibility,
        )


def _resolve_version(reference_root, profile_id, version):
    profile_root = reference_root / profile_id
    if version:
        return _safe_identifier(version, "profile version")
    versions = sorted(
        path.name
        for path in profile_root.iterdir()
        if path.is_dir() and (path / "manifest.json").is_file()
    ) if profile_root.is_dir() else []
    if not versions:
        sys.exit(
            f"[references] ERROR: no versions found for profile '{profile_id}' "
            f"under {profile_root}."
        )
    if len(versions) > 1:
        sys.exit(
            f"[references] ERROR: multiple versions found for profile "
            f"'{profile_id}': {versions}. Specify a version explicitly."
        )
    return versions[0]


def load_chrom_sizes(path):
    """Read a two-column chromosome-size file with strict validation."""
    chrom_sizes = {}
    try:
        with open(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or not fields[0]:
                    sys.exit(
                        "[references] ERROR: invalid chromosome sizes line "
                        f"{line_number} in {path}."
                    )
                try:
                    length = int(fields[1])
                except ValueError:
                    sys.exit(
                        "[references] ERROR: chromosome length must be an integer "
                        f"on line {line_number} in {path}."
                    )
                if length < 1:
                    sys.exit(
                        "[references] ERROR: chromosome length must be positive "
                        f"on line {line_number} in {path}."
                    )
                if fields[0] in chrom_sizes:
                    sys.exit(
                        "[references] ERROR: duplicate chromosome "
                        f"{fields[0]!r} in {path}."
                    )
                chrom_sizes[fields[0]] = length
    except OSError as exc:
        sys.exit(f"[references] ERROR: could not read chromosome sizes {path}: {exc}")
    if not chrom_sizes:
        sys.exit(f"[references] ERROR: chromosome sizes file is empty: {path}")
    return chrom_sizes


def validate_target_bed(target_bed, chrom_sizes_path):
    """Validate BED coordinates against the selected genome chromosome sizes."""
    chrom_sizes = load_chrom_sizes(chrom_sizes_path)
    records = 0
    try:
        with open(target_bed) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if (
                    not stripped
                    or stripped.startswith("#")
                    or stripped.startswith("track ")
                    or stripped.startswith("browser ")
                ):
                    continue
                fields = stripped.split("\t")
                if len(fields) < 3:
                    sys.exit(
                        f"[references] ERROR: target BED line {line_number} "
                        "requires at least three tab-separated columns."
                    )
                chrom = fields[0]
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    sys.exit(
                        f"[references] ERROR: target BED line {line_number} has "
                        "non-integer coordinates."
                    )
                if start < 0 or end <= start:
                    sys.exit(
                        f"[references] ERROR: target BED has invalid interval on "
                        f"line {line_number}: {chrom}:{start}-{end}."
                    )
                if chrom not in chrom_sizes:
                    sys.exit(
                        f"[references] ERROR: target BED uses unknown contig "
                        f"{chrom!r} on line {line_number}."
                    )
                if end > chrom_sizes[chrom]:
                    sys.exit(
                        f"[references] ERROR: target BED interval exceeds chromosome "
                        f"length for {chrom} on line {line_number}."
                    )
                records += 1
    except OSError as exc:
        sys.exit(f"[references] ERROR: could not read target BED {target_bed}: {exc}")
    if records == 0:
        sys.exit(f"[references] ERROR: target BED contains no intervals: {target_bed}")


def load_reference_profile(
    reference_root,
    profile_id=DEFAULT_PROFILE_ID,
    version=None,
    *,
    verify_checksums=False,
    validate_compatibility=False,
):
    """Load one local profile from ``root/profile_id/version/manifest.json``."""
    root = Path(reference_root).expanduser().resolve()
    profile_id = _safe_identifier(profile_id, "profile id")
    version = _resolve_version(root, profile_id, version)
    profile_dir = (root / profile_id / version).resolve()
    manifest_path = profile_dir / "manifest.json"
    if not manifest_path.is_file():
        sys.exit(f"[references] ERROR: profile manifest not found: {manifest_path}")

    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        sys.exit(f"[references] ERROR: invalid manifest {manifest_path}: {exc}")

    expected = {
        "manifest_version": 1,
        "profile_id": profile_id,
        "version": version,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            sys.exit(
                f"[references] ERROR: manifest {key} must be {value!r}, got "
                f"{manifest.get(key)!r}."
            )
    for key in ("assay", "genome"):
        if not manifest.get(key):
            sys.exit(f"[references] ERROR: manifest field '{key}' is required.")

    component_specs = manifest.get("components")
    if not isinstance(component_specs, dict):
        sys.exit("[references] ERROR: manifest 'components' must be an object.")
    missing = sorted(REQUIRED_COMPONENTS - component_specs.keys())
    if missing:
        sys.exit(
            f"[references] ERROR: profile is missing required components: {missing}."
        )

    components = {}
    component_hashes = {}
    for name, spec in component_specs.items():
        if isinstance(spec, str):
            spec = {"path": spec}
        if not isinstance(spec, dict) or not spec.get("path"):
            sys.exit(
                f"[references] ERROR: component '{name}' requires a relative path."
            )
        component_path = (profile_dir / spec["path"]).resolve()
        try:
            component_path.relative_to(profile_dir)
        except ValueError:
            sys.exit(
                f"[references] ERROR: component '{name}' escapes profile directory."
            )
        if not component_path.is_file():
            sys.exit(
                f"[references] ERROR: component '{name}' not found: {component_path}"
            )
        expected_hash = spec.get("sha256")
        if expected_hash is not None and (
            not isinstance(expected_hash, str) or not _SHA256.fullmatch(expected_hash)
        ):
            sys.exit(
                f"[references] ERROR: component '{name}' has an invalid sha256 value."
            )
        actual_hash = None
        if verify_checksums:
            actual_hash = sha256_file(component_path)
        if verify_checksums and expected_hash and actual_hash != expected_hash.lower():
            sys.exit(
                f"[references] ERROR: checksum mismatch for component '{name}'."
            )
        components[name] = str(component_path)
        component_hashes[name] = actual_hash or expected_hash

    if validate_compatibility:
        validate_target_bed(components["target_bed"], components["chrom_sizes"])

    return {
        "profile_id": profile_id,
        "version": version,
        "assay": manifest["assay"],
        "genome": manifest["genome"],
        "profile_dir": str(profile_dir),
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "components": components,
        "component_hashes": component_hashes,
        "acquisition": manifest.get("acquisition", {"mode": "local"}),
    }


def acquire_reference_profile(
    *,
    mode,
    reference_root,
    profile_id=DEFAULT_PROFILE_ID,
    version=None,
    verify_checksums=False,
    validate_compatibility=False,
    registry=None,
    opener=None,
):
    """Resolve a local profile or atomically install one from a pinned registry."""
    if mode == "managed":
        return _acquire_managed_profile(
            reference_root,
            profile_id,
            version,
            registry=registry,
            opener=opener,
            verify_checksums=verify_checksums,
            validate_compatibility=validate_compatibility,
        )
    if mode != "local":
        sys.exit("[references] ERROR: reference mode must be 'local' or 'managed'.")
    return load_reference_profile(
        reference_root,
        profile_id,
        version,
        verify_checksums=verify_checksums,
        validate_compatibility=validate_compatibility,
    )
