#!/usr/bin/env python3
"""Prepare the embedded, privacy-safe identity for a CFTK release build."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys


REVISION_PATTERN = re.compile(r"^[0-9a-f]{40}$")
VERSION_PATTERN = re.compile(r"(?m)^version\s*=\s*[\"']([^\"']+)[\"']\s*$")


def _repo_root(path: Path) -> Path:
    path = Path(path).expanduser().resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / "pyproject.toml").is_file() and (
            candidate / "src/cftk_provenance/__init__.py"
        ).is_file():
            return candidate
    raise SystemExit("could not locate the CFTK source root")


def _version(root: Path) -> str:
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = VERSION_PATTERN.search(text)
    if not match:
        raise SystemExit("pyproject.toml has no explicit project version")
    return match.group(1)


def prepare(root: Path, *, revision: str, tag: str | None, output: Path) -> dict:
    root = _repo_root(root)
    if not REVISION_PATTERN.fullmatch(revision):
        raise SystemExit("release revision must be a 40-character lowercase Git SHA")
    version = _version(root)
    expected_tag = f"v{version}"
    if tag is not None and tag != expected_tag:
        raise SystemExit(f"release tag {tag!r} does not match project version {expected_tag!r}")

    sys.path.insert(0, str(root / "src"))
    from cftk_provenance import release_identity, source_tree_sha256

    identity = release_identity(version, revision, source_tree_sha256(root))
    output = Path(output)
    if not output.is_absolute():
        output = root / output
    output = output.resolve()
    try:
        output.relative_to(root)
    except ValueError as exc:
        raise SystemExit("build identity output must remain inside the source root") from exc
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return identity


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--revision", required=True)
    parser.add_argument("--tag")
    parser.add_argument(
        "--output", type=Path, default=Path("src/cftk_provenance/build.json")
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    identity = prepare(
        args.root, revision=args.revision, tag=args.tag, output=args.output
    )
    print(json.dumps(identity, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
