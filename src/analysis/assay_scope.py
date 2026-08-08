"""Assay-aware input scoping for targeted fragmentomics analyses.

The default Twist Human Methylome assay is a capture panel, not a genome-wide
cfDNA assay. This module derives deterministic panel-overlap regions and
panel-read BAMs for WPS, occupancy, and DELFI while leaving the underlying
analysis implementations unchanged.
"""

from __future__ import annotations

import bisect
import contextlib
import hashlib
import json
import os
import tempfile
from pathlib import Path

from util import recorded_run


SCOPE_CHOICES = ("auto", "panel", "genome")
TARGETED_ASSAYS = {"twist_human_methylome"}
SCOPE_OUTPUT_KEYS = {
    "occupancy": "occ_out",
    "wps": "wps_out",
    "delfi": "delfi_out",
}
PANEL_NOTE = (
    "Targeted-panel mode restricts WPS, occupancy, and DELFI to reads and "
    "intervals overlapping the configured target BED; these outputs are not "
    "genome-wide measurements."
)


class ScopeError(RuntimeError):
    """A requested fragmentomics scope cannot be prepared safely."""


def _reference(cfg):
    return cfg.get("reference_data", {}) if isinstance(cfg, dict) else {}


def _configured_scope(cfg, requested=None):
    if requested is not None:
        return requested
    return (
        cfg.get("analysis", {})
        .get("frag", {})
        .get("scope", "auto")
    )


def _is_default_targeted_assay(cfg):
    assay_value = cfg.get("assay")
    profile = cfg.get("reference_profile", {})
    if isinstance(profile, str):
        profile_id = profile.lower()
    else:
        profile_id = str(profile.get("id", "")).lower()
    if assay_value is None and not profile_id:
        # Legacy configs do not carry assay metadata. Preserve their existing
        # genome-wide fragmentomics behavior unless they opt into panel mode.
        return False
    assay = str(assay_value or "").lower()
    return assay in TARGETED_ASSAYS or profile_id.startswith("twist_human_methylome")


def _absolute(value):
    return str(Path(value).expanduser().resolve()) if value else ""


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _fingerprint(path):
    path = Path(path)
    if not path.is_file():
        return f"missing:{path.resolve()}"
    return _sha256(path)


def _bam_signature(path):
    """Use cheap input identity for large BAMs without hashing them on every plan."""

    if not path:
        return {"path": "", "size": None, "mtime_ns": None}
    path = Path(path).expanduser().resolve()
    try:
        stat = path.stat()
    except OSError:
        return {"path": str(path), "size": None, "mtime_ns": None}
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def resolve_scope(cfg, requested=None):
    """Resolve ``auto``, ``panel``, or ``genome`` for fragmentomics.

    ``auto`` selects panel mode only for the default Twist targeted assay. A
    custom targeted profile should set the schema-v2
    ``fragmentomics_scope`` field (or the resolved ``analysis.frag.scope``)
    to ``panel``; ``genome`` is available as an explicit expert override.
    """

    configured = _configured_scope(cfg, requested)
    if configured not in SCOPE_CHOICES:
        raise ScopeError(
            f"fragmentomics scope must be one of {', '.join(SCOPE_CHOICES)}, "
            f"got {configured!r}"
        )
    reference = _reference(cfg)
    target_bed = _absolute(reference.get("target_bed"))
    targeted_assay = _is_default_targeted_assay(cfg)
    if configured == "auto":
        if targeted_assay:
            if not target_bed:
                raise ScopeError(
                    "the default targeted assay requires reference_data.target_bed "
                    "for automatic fragmentomics scoping"
                )
            mode = "panel"
            reason = "default targeted assay"
        else:
            mode = "genome"
            reason = "non-targeted assay or custom profile"
    else:
        mode = configured
        reason = "explicit user override"

    if mode == "panel":
        if not target_bed:
            raise ScopeError(
                "panel fragmentomics scope requires reference_data.target_bed"
            )
        if not Path(target_bed).is_file():
            raise ScopeError(f"target BED does not exist: {target_bed}")

    return {
        "requested": configured,
        "mode": mode,
        "reason": reason,
        "assay": cfg.get("assay", "twist_human_methylome"),
        "targeted_assay": targeted_assay,
        "target_bed": target_bed if mode == "panel" else None,
        "note": PANEL_NOTE if mode == "panel" else "Genome-wide fragmentomics scope.",
    }


def scope_paths(cfg, paths, sample_names=(), bam_paths=()):
    """Return deterministic intermediate paths for a panel scope."""

    info = resolve_scope(cfg, "panel")
    reference = _reference(cfg)
    source_records = []
    for label in ("target_bed", "tss_pas_bed", "bins"):
        value = _absolute(reference.get(label))
        source_records.append({
            "label": label,
            "path": value,
            "sha256": _fingerprint(value) if value else "missing",
        })
    sample_records = [
        {"sample": str(name), **_bam_signature(path)}
        for name, path in zip(sample_names, bam_paths)
    ]
    key = hashlib.sha256(
        json.dumps(
            {"references": source_records, "samples": sample_records},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()[:16]
    root = Path(paths["fragmentomics"]) / "_scope" / f"panel-{key}"
    bam_dir = root / "bams"
    names = [str(name) for name in sample_names]
    return {
        "key": key,
        "root": root,
        "scope_json": root / "scope.json",
        "regions_bed": root / "target_overlap_regions.bed",
        "bins_bed": root / "target_overlap_bins.bed",
        "bam_dir": bam_dir,
        # Keep the public sample stem stable when fragmentomics wrappers derive
        # their output name from the scoped BAM basename.
        "bams": {name: bam_dir / f"{name}.markdup.bam" for name in names},
        "source_records": source_records,
        "sample_records": sample_records,
        "target_bed": Path(info["target_bed"]),
    }


def describe_scope(cfg, paths=None, sample_names=(), requested=None, bam_paths=()):
    """Return scope metadata suitable for plans, manifests, and doctor."""

    info = resolve_scope(cfg, requested)
    if info["mode"] != "panel" or paths is None:
        return info
    derived = scope_paths(cfg, paths, sample_names, bam_paths)
    info = dict(info)
    info.update({
        "scope_key": derived["key"],
        "scope_root": str(derived["root"]),
        "scope_json": str(derived["scope_json"]),
        "regions_bed": str(derived["regions_bed"]),
        "bins_bed": str(derived["bins_bed"]),
        "target_sha256": _fingerprint(info["target_bed"]),
        "source_records": derived["source_records"],
        "sample_inputs": derived["sample_records"],
        "sample_bams": {name: str(path) for name, path in derived["bams"].items()},
    })
    return info


def _read_bed(path):
    rows = []
    path = Path(path)
    try:
        handle = path.open(encoding="utf-8")
    except OSError as exc:
        raise ScopeError(f"could not read BED {path}: {exc}") from exc
    with handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.lstrip().startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ScopeError(f"BED {path}:{line_number} has fewer than 3 columns")
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise ScopeError(
                    f"BED {path}:{line_number} has non-integer coordinates"
                ) from exc
            if start < 0 or end <= start:
                raise ScopeError(
                    f"BED {path}:{line_number} has invalid interval {start}-{end}"
                )
            name = fields[3] if len(fields) > 3 and fields[3] else f"{fields[0]}:{start}-{end}"
            rows.append((fields[0], start, end, name))
    if not rows:
        raise ScopeError(f"BED is empty: {path}")
    return rows


def _merged_intervals(rows):
    merged = []
    for chrom, start, end, _name in sorted(rows, key=lambda row: (row[0], row[1], row[2])):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def _clip_bed(source, target, destination):
    source_rows = _read_bed(source)
    target_rows = _merged_intervals(_read_bed(target))
    target_by_chrom = {}
    for chrom, start, end in target_rows:
        target_by_chrom.setdefault(chrom, []).append((start, end))
    starts_by_chrom = {
        chrom: [start for start, _end in values]
        for chrom, values in target_by_chrom.items()
    }

    clipped = []
    for chrom, start, end, name in source_rows:
        intervals = target_by_chrom.get(chrom, ())
        starts = starts_by_chrom.get(chrom, ())
        index = max(0, bisect.bisect_left(starts, start) - 1)
        piece = 0
        while index < len(intervals) and intervals[index][0] < end:
            target_start, target_end = intervals[index]
            overlap_start = max(start, target_start)
            overlap_end = min(end, target_end)
            if overlap_end > overlap_start:
                clipped.append((
                    chrom,
                    overlap_start,
                    overlap_end,
                    f"{name}__panel_{piece}",
                ))
                piece += 1
            index += 1

    _atomic_write_text(
        destination,
        "".join(f"{chrom}\t{start}\t{end}\t{name}\n" for chrom, start, end, name in clipped),
    )
    return len(clipped)


def _nonempty(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def _atomic_write_text(destination, text):
    """Atomically replace a small shared scope metadata or BED artifact."""

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(temporary, destination)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


@contextlib.contextmanager
def _artifact_lock(path):
    """Serialize creation of one reusable scope artifact across job tasks."""

    import fcntl

    path = Path(path)
    lock = Path(f"{path}.lock")
    lock.parent.mkdir(parents=True, exist_ok=True)
    with lock.open("w", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _bam_index(path):
    path = Path(path)
    candidates = (Path(str(path) + ".bai"), path.with_suffix(".bai"))
    return next((candidate for candidate in candidates if candidate.is_file()), candidates[0])


def _panel_bam(bam, target_bed, output, cores, sample):
    bam = Path(bam).expanduser().resolve()
    output = Path(output)
    index = Path(str(output) + ".bai")
    with _artifact_lock(output):
        if _nonempty(output) and _nonempty(_bam_index(output)):
            return output
        if not bam.is_file():
            raise ScopeError(f"BAM for {sample} does not exist: {bam}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.unlink(missing_ok=True)
        index.unlink(missing_ok=True)
        threads = max(1, int(cores or 1))
        view_command = [
            "samtools", "view", "-@", str(threads), "-bh", "-L",
            str(target_bed), "-o", str(output), str(bam),
        ]
        try:
            result = recorded_run(view_command, label=f"fragmentomics scope [{sample}] view")
        except OSError as exc:
            raise ScopeError(
                f"samtools is required to create the panel-read BAM for {sample}: {exc}"
            ) from exc
        if result.returncode != 0:
            output.unlink(missing_ok=True)
            raise ScopeError(
                f"samtools could not create the panel-read BAM for {sample} "
                f"(exit {result.returncode})"
            )
        index_command = ["samtools", "index", "-@", str(threads), str(output)]
        try:
            result = recorded_run(index_command, label=f"fragmentomics scope [{sample}] index")
        except OSError as exc:
            output.unlink(missing_ok=True)
            raise ScopeError(
                f"samtools is required to index the panel-read BAM for {sample}: {exc}"
            ) from exc
        if result.returncode != 0:
            output.unlink(missing_ok=True)
            index.unlink(missing_ok=True)
            raise ScopeError(
                f"samtools could not index the panel-read BAM for {sample} "
                f"(exit {result.returncode})"
            )
        if not _nonempty(output) or not _nonempty(_bam_index(output)):
            raise ScopeError(f"samtools produced an incomplete panel-read BAM for {sample}")
    return output


def prepare_scope(
    cfg,
    paths,
    samples,
    bam_paths,
    *,
    requested=None,
    cores=1,
    kinds=(),
    materialize_samples=None,
):
    """Materialize the inputs required by selected targeted fragmentomics stages."""

    info = resolve_scope(cfg, requested)
    kinds = set(kinds)
    original_bams = [str(Path(path).expanduser().resolve()) for path in bam_paths]
    names = [sample["name"] for sample in samples]
    requested_names = set(materialize_samples or names)
    unknown = requested_names.difference(names)
    if unknown:
        raise ScopeError(
            "requested scope materialization for unknown sample(s): "
            + ", ".join(sorted(unknown))
        )
    selected_pairs = [
        (sample, bam)
        for sample, bam in zip(samples, original_bams)
        if sample["name"] in requested_names
    ]
    if info["mode"] == "genome":
        reference = _reference(cfg)
        return {
            "info": info,
            "bam_paths": [bam for _sample, bam in selected_pairs],
            "region_bed": _absolute(reference.get("tss_pas_bed")),
            "bins": _absolute(reference.get("bins")),
        }

    derived = scope_paths(cfg, paths, names, original_bams)
    root = derived["root"]
    root.mkdir(parents=True, exist_ok=True)
    reference = _reference(cfg)
    region_count = None
    bins_count = None
    if kinds.intersection({"occupancy", "wps"}):
        source = _absolute(reference.get("tss_pas_bed"))
        if not source or not Path(source).is_file():
            raise ScopeError(
                "targeted WPS/occupancy requires reference_data.tss_pas_bed"
            )
        region_count = _clip_bed(source, derived["target_bed"], derived["regions_bed"])
        if region_count == 0:
            raise ScopeError(
                "the Twist target BED has no overlap with tss_pas_bed; "
                "provide a panel-compatible region reference or use an explicit "
                "genome scope"
            )
    if "delfi" in kinds:
        source = _absolute(reference.get("bins"))
        if not source or not Path(source).is_file():
            raise ScopeError("targeted DELFI requires reference_data.bins")
        bins_count = _clip_bed(source, derived["target_bed"], derived["bins_bed"])
        if bins_count == 0:
            raise ScopeError(
                "the target BED has no overlap with the DELFI bins reference; "
                "provide panel-compatible bins or use an explicit genome scope"
            )

    scoped_bams = []
    for sample, bam in selected_pairs:
        scoped_bams.append(
            str(_panel_bam(
                bam,
                derived["target_bed"],
                derived["bams"][sample["name"]],
                cores,
                sample["name"],
            ))
        )

    info = describe_scope(
        cfg, paths, names, requested="panel", bam_paths=original_bams
    )
    info.update({
        "region_count": region_count,
        "bins_count": bins_count,
        "prepared_kinds": sorted(kinds),
    })
    metadata = {
        **info,
        "sample_bams": {
            sample["name"]: str(derived["bams"][sample["name"]])
            for sample in samples
        },
    }
    _atomic_write_text(
        derived["scope_json"], json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    return {
        "info": metadata,
        "bam_paths": scoped_bams,
        "region_bed": str(derived["regions_bed"]) if region_count is not None else _absolute(reference.get("tss_pas_bed")),
        "bins": str(derived["bins_bed"]) if bins_count is not None else _absolute(reference.get("bins")),
    }


def scope_metadata_path(paths, kind):
    """Return the user-facing scope sidecar path for one scoped stage."""

    try:
        output_dir = paths[SCOPE_OUTPUT_KEYS[kind]]
    except (KeyError, TypeError) as exc:
        raise ScopeError(f"no fragmentomics output directory is configured for {kind}") from exc
    return Path(output_dir) / "fragmentomics_scope.json"


def write_scope_metadata(paths, kind, scope):
    """Write a compact, stage-local scope record beside fragmentomics outputs."""

    if kind not in SCOPE_OUTPUT_KEYS:
        raise ScopeError(f"scope metadata is not supported for stage {kind!r}")
    destination = scope_metadata_path(paths, kind)
    resolved = dict(scope or {})
    payload = {
        "schema_version": 1,
        "stage": kind,
        "mode": resolved.get("mode"),
        "requested": resolved.get("requested"),
        "reason": resolved.get("reason"),
        "assay": resolved.get("assay"),
        "targeted_assay": resolved.get("targeted_assay"),
        "target_bed": resolved.get("target_bed"),
        "target_sha256": resolved.get("target_sha256"),
        "scope_root": resolved.get("scope_root"),
        "canonical_scope_json": resolved.get("scope_json"),
        "regions_bed": resolved.get("regions_bed"),
        "bins_bed": resolved.get("bins_bed"),
        "region_count": resolved.get("region_count"),
        "bins_count": resolved.get("bins_count"),
        "sample_inputs": resolved.get("sample_inputs"),
        "sample_bams": resolved.get("sample_bams"),
        "note": resolved.get("note"),
        "resolved_scope": resolved,
    }
    _atomic_write_text(
        destination, json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    return destination


def scope_artifact_paths(cfg, paths, sample_names, kind, requested=None, bam_paths=()):
    """Return intermediate artifact paths required by one stage."""

    info = resolve_scope(cfg, requested)
    if info["mode"] != "panel":
        return []
    derived = scope_paths(cfg, paths, sample_names, bam_paths)
    artifacts = [derived["scope_json"]]
    if kind in {"occupancy", "wps"}:
        artifacts.append(derived["regions_bed"])
    if kind == "delfi":
        artifacts.append(derived["bins_bed"])
    artifacts.extend(derived["bams"].values())
    artifacts.extend(Path(str(path) + ".bai") for path in derived["bams"].values())
    return artifacts
