"""Versioned reference-profile discovery and validation for CFTK."""

from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path


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
    }


def acquire_reference_profile(
    *,
    mode,
    reference_root,
    profile_id=DEFAULT_PROFILE_ID,
    version=None,
    verify_checksums=False,
    validate_compatibility=False,
):
    """Resolve local profiles; reject managed mode until a pinned registry exists."""
    if mode == "managed":
        sys.exit(
            "[references] ERROR: no immutable managed reference profile is "
            "available yet. Use --reference-mode local with --reference-root."
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
