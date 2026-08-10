#!/usr/bin/env python3
"""Validate CFTK release archives and scan tracked public text for leaks."""

from __future__ import annotations

import argparse
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
import zipfile


REQUIRED_PATHS = {
    "analysis/__init__.py",
    "analysis/dmr_annotation.r",
    "analysis_workflow.py",
    "cftk.py",
    "cftk_registry/__init__.py",
    "cftk_registry/registry.json",
    "doctor.py",
    "init.py",
    "job_plan.py",
    "process.py",
    "reference_profiles.py",
    "report/__init__.py",
    "report/report_generator.py",
    "report/report_template.html",
    "report/software_list.json",
    "resource_planning.py",
    "run_workflow.py",
    "util.py",
    "validation_reports.py",
    "visualization/__init__.py",
}
SDIST_REQUIRED_PATHS = {
    "scripts/__init__.py",
    "scripts/validation/__init__.py",
    "scripts/validation/check_distribution.py",
}

FORBIDDEN_SUFFIXES = {
    ".ipynb",
    ".npy",
    ".npz",
    ".pickle",
    ".pkl",
    ".pyc",
}
FORBIDDEN_PARTS = {"__pycache__", ".ipynb_checkpoints", "build", "data"}

# These patterns describe categories, not identifiers from any private cohort.
PRIVATE_PATTERNS = {
    "absolute shared-filesystem path": re.compile(
        r"/(?:dfs\d+|data/homezvol\d+)(?:/[^\s\"'<>]*)?", re.IGNORECASE
    ),
    "sample-like identifier": re.compile(
        r"\b(?:control|case|patient|subject|sample|sals|als)"
        r"[_-]?[a-z]*\d{2,}\b",
        re.IGNORECASE,
    ),
}
DATA_URI_PATTERN = re.compile(
    r"data:[^;,\"'\s]+;base64,[a-z0-9+/=\r\n]+", re.IGNORECASE
)


class ValidationError(RuntimeError):
    """Raised when a release artifact violates the public package contract."""


def _archive_members(path: Path) -> dict[str, bytes]:
    if path.suffix == ".whl":
        with zipfile.ZipFile(path) as archive:
            return {
                member.filename: archive.read(member)
                for member in archive.infolist()
                if not member.is_dir()
            }
    if path.name.endswith((".tar.gz", ".tar.bz2", ".tar.xz")):
        with tarfile.open(path, "r:*") as archive:
            members = {}
            for member in archive.getmembers():
                if not member.isfile():
                    continue
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValidationError(f"cannot read archive member: {member.name}")
                members[member.name] = stream.read()
            return members
    raise ValidationError(f"unsupported distribution archive: {path}")


def _normalise_members(path: Path, raw: dict[str, bytes]) -> dict[str, bytes]:
    parsed = [PurePosixPath(name) for name in raw]
    unsafe = [
        str(name) for name in parsed if name.is_absolute() or ".." in name.parts
    ]
    if unsafe:
        raise ValidationError(f"unsafe archive paths in {path.name}: {unsafe}")
    if not parsed:
        return {}

    strip_root = path.suffix != ".whl"
    root = parsed[0].parts[0] if strip_root else None
    if strip_root and any(not name.parts or name.parts[0] != root for name in parsed):
        raise ValidationError(f"sdist members do not share one root directory: {path.name}")

    normalised = {}
    for original, member in zip(raw, parsed):
        parts = member.parts[1:] if strip_root else member.parts
        name = PurePosixPath(*parts).as_posix()
        if not name or name in normalised:
            raise ValidationError(
                f"duplicate or empty archive path in {path.name}: {name!r}"
            )
        normalised[name] = raw[original]
    return normalised


def _package_relative_paths(path: Path, members: dict[str, bytes]) -> set[str]:
    if path.suffix == ".whl":
        return set(members)
    return {
        name.removeprefix("src/")
        for name in members
        if name.startswith("src/")
    }


def _forbidden_paths(path: Path, members: dict[str, bytes]) -> list[str]:
    violations = []
    for name in members:
        member = PurePosixPath(name)
        if member.suffix.lower() in FORBIDDEN_SUFFIXES:
            violations.append(name)
            continue
        if any(part in FORBIDDEN_PARTS for part in member.parts):
            violations.append(name)
            continue
        egg_info = [part for part in member.parts if part.endswith(".egg-info")]
        normal_sdist_metadata = (
            path.suffix != ".whl"
            and egg_info
            and member.parts[:2] == ("src", "cftk.egg-info")
        )
        if egg_info and not normal_sdist_metadata:
            violations.append(name)
    return violations


def _privacy_violations(files: dict[str, bytes]) -> list[str]:
    violations = []
    for name, content in files.items():
        try:
            text = content.decode("utf-8")
        except UnicodeDecodeError:
            continue
        # Encoded image bytes can coincidentally contain path-like substrings.
        text = DATA_URI_PATTERN.sub("data:embedded-content", text)
        for label, pattern in PRIVATE_PATTERNS.items():
            match = pattern.search(text)
            if match:
                line = text.count("\n", 0, match.start()) + 1
                violations.append(f"{name}:{line}: {label}")
    return violations


def validate_archive(path: Path) -> None:
    """Validate one wheel or sdist, raising ValidationError on any violation."""
    path = Path(path)
    if not path.is_file():
        raise ValidationError(f"distribution archive not found: {path}")
    members = _normalise_members(path, _archive_members(path))
    package_paths = _package_relative_paths(path, members)
    missing = sorted(REQUIRED_PATHS - package_paths)
    if path.suffix != ".whl":
        missing.extend(sorted(SDIST_REQUIRED_PATHS - set(members)))
    forbidden = sorted(_forbidden_paths(path, members))
    private = sorted(_privacy_violations(members))
    problems = []
    if missing:
        problems.append(f"missing required package paths: {missing}")
    if forbidden:
        problems.append(f"forbidden archive paths: {forbidden}")
    if private:
        problems.append(f"possible private content: {private}")
    if problems:
        raise ValidationError(f"{path.name}: " + "; ".join(problems))


def tracked_text_files(root: Path) -> dict[str, bytes]:
    """Read Git-tracked files; undecodable binary files are ignored by the scan."""
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=True,
        capture_output=True,
    )
    files = {}
    for raw_name in result.stdout.split(b"\0"):
        if not raw_name:
            continue
        name = raw_name.decode("utf-8")
        path = root / name
        if path.is_file():
            files[name] = path.read_bytes()
    return files


def validate_source_tree(root: Path) -> None:
    root = Path(root).resolve()
    private = sorted(_privacy_violations(tracked_text_files(root)))
    if private:
        raise ValidationError("possible private content in tracked files: " + str(private))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "archives", nargs="*", type=Path, help="Wheel or sdist archives"
    )
    parser.add_argument(
        "--source-tree",
        type=Path,
        help="Also scan Git-tracked files below this repository root",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.archives and args.source_tree is None:
        raise SystemExit("provide at least one archive or --source-tree")
    try:
        for archive in args.archives:
            validate_archive(archive)
            print(f"validated distribution: {archive}")
        if args.source_tree is not None:
            validate_source_tree(args.source_tree)
            print(f"validated tracked-source privacy: {args.source_tree}")
    except (
        OSError,
        subprocess.CalledProcessError,
        tarfile.TarError,
        zipfile.BadZipFile,
        ValidationError,
    ) as exc:
        raise SystemExit(f"distribution validation failed: {exc}") from exc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
