"""Read-only project readiness diagnostics for the core CFTK process."""

from __future__ import annotations

import contextlib
import gzip
import importlib.metadata
import io
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from init import (
    SCHEMA_VERSION,
    get_all_samples,
    get_bam,
    get_sequence_dictionary_path,
    get_work_paths,
    load_config,
    resolve_reference_profile,
    validate_fasta_index,
)
from reference_profiles import sha256_file, validate_target_bed
from resource_planning import (
    detect_scheduler_allocation,
    ensure_scheduler_capacity,
    plan_parallelism,
)


_STEP_KEYS = {
    1: "step1_trimming",
    2: "step2_alignment",
    3: "step3_markdup",
    4: "step4_methylation",
}
_BWA_SUFFIXES = ("", ".amb", ".ann", ".bwt", ".pac", ".sa")
_TOOL_PROBE_TIMEOUT = 60


class _Checks:
    def __init__(self):
        self.items = []

    def add(self, check_id, status, summary, *, remedy=None, details=None):
        item = {"id": check_id, "status": status, "summary": str(summary)}
        if remedy:
            item["remedy"] = str(remedy)
        if details is not None:
            item["details"] = details
        self.items.append(item)

    def pass_(self, check_id, summary, **kwargs):
        self.add(check_id, "PASS", summary, **kwargs)

    def warn(self, check_id, summary, **kwargs):
        self.add(check_id, "WARN", summary, **kwargs)

    def fail(self, check_id, summary, **kwargs):
        self.add(check_id, "FAIL", summary, **kwargs)


def _exception_message(exc, captured):
    if isinstance(exc, SystemExit) and isinstance(exc.code, str):
        return exc.code
    text = captured.strip()
    if text:
        return text.splitlines()[-1]
    if isinstance(exc, SystemExit):
        return f"validation exited with status {exc.code}"
    return str(exc)


def _capture_call(function, *args, **kwargs):
    stream = io.StringIO()
    try:
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
            value = function(*args, **kwargs)
        return value, None
    except (Exception, SystemExit) as exc:
        return None, _exception_message(exc, stream.getvalue())


def _load_raw_config(path, checks):
    if not path.is_file():
        checks.fail(
            "project.config",
            f"Config not found: {path}",
            remedy="Run 'cftk init' in the project directory.",
        )
        return None
    try:
        raw = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        checks.fail(
            "project.config", f"Config is not readable JSON: {exc}",
            remedy="Correct the config JSON before running the process.",
        )
        return None
    if not isinstance(raw, dict):
        checks.fail("project.config", "Config must contain a JSON object.")
        return None
    checks.pass_("project.config", f"Readable config: {path}")
    return raw


def _check_runtime(checks):
    version = ".".join(str(value) for value in sys.version_info[:3])
    if sys.version_info < (3, 9):
        checks.fail("runtime.python", f"Python {version} is unsupported; CFTK requires 3.9+.")
    else:
        checks.pass_("runtime.python", f"Python {version}")
    try:
        version = importlib.metadata.version("cftk")
        checks.pass_("runtime.cftk", f"CFTK {version}")
    except importlib.metadata.PackageNotFoundError:
        checks.pass_("runtime.cftk", "CFTK source checkout is importable")


def _tool_requirements(cfg, steps, skip_picard_metrics):
    samples = get_all_samples(cfg) if cfg else []
    has_fastq = any(sample.get("input_type", "").lower() == "fastq" for sample in samples)
    requirements = {}

    def require(key, executable, version_args=("--version",), allowed=(0,)):
        requirements[key] = (executable, version_args, allowed, False)

    def require_picard():
        # Picard 3.0.0 prints its version but returns 1 for this probe.
        require("picard", "picard", ("MarkDuplicates", "--version"), (0, 1))

    if 1 in steps and has_fastq:
        tool = cfg["process"][_STEP_KEYS[1]].get("tool", "trim_galore").lower()
        if tool == "trim_galore":
            require("trim_galore", "trim_galore")
            require("fastqc", "fastqc")
        elif tool == "fastp":
            require("fastp", "fastp")
        else:
            requirements[f"unsupported_step1_{tool}"] = (tool, (), (), False)

    if 2 in steps and has_fastq:
        tool = cfg["process"][_STEP_KEYS[2]].get("tool", "bwameth").lower()
        if tool == "bwameth":
            require("bwameth", "bwameth.py")
            require("bwa", "bwa", (), (0, 1))
            require("sambamba", "sambamba")
            require("samtools", "samtools")
        elif tool == "bismark":
            require("bismark", "bismark")
        else:
            requirements[f"unsupported_step2_{tool}"] = (tool, (), (), False)

    if 3 in steps:
        tool = cfg["process"][_STEP_KEYS[3]].get("tool", "sambamba").lower()
        if tool == "sambamba":
            require("sambamba", "sambamba")
        elif tool == "picard":
            require_picard()
        elif tool == "samblaster":
            require("samblaster", "samblaster")
        else:
            requirements[f"unsupported_step3_{tool}"] = (tool, (), (), False)
        require("samtools", "samtools")
        if not skip_picard_metrics:
            require_picard()

    if 4 in steps:
        tool = cfg["process"][_STEP_KEYS[4]].get("tool", "methyldackel").lower()
        if tool == "methyldackel":
            require("methyldackel", "MethylDackel")
        elif tool == "bismark_extractor":
            require("bismark_extractor", "bismark_methylation_extractor")
        else:
            requirements[f"unsupported_step4_{tool}"] = (tool, (), (), False)
        if len(samples) > 1:
            require("bedtools", "bedtools")

    if steps.intersection({3, 4}) and any(
        sample.get("input_type", "").lower() == "bam" for sample in samples
    ):
        require("samtools", "samtools")

    if steps.intersection({1, 2}):
        requirements["multiqc"] = ("multiqc", ("--version",), (0,), True)
    return requirements


def _check_tools(checks, cfg, steps, skip_picard_metrics):
    for key, (executable, version_args, allowed, optional) in sorted(
        _tool_requirements(cfg, steps, skip_picard_metrics).items()
    ):
        check_id = f"tool.{key}"
        path = shutil.which(executable)
        if not path:
            method = checks.warn if optional else checks.fail
            method(
                check_id,
                f"{executable} is not available on PATH.",
                remedy=(
                    "Install the pinned processing environment or activate the "
                    "environment used for CFTK."
                ),
            )
            continue
        if not allowed:
            checks.fail(check_id, f"Configured process tool is unsupported: {executable}")
            continue
        try:
            completed = subprocess.run(
                [path, *version_args],
                capture_output=True,
                text=True,
                timeout=_TOOL_PROBE_TIMEOUT,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            method = checks.warn if optional else checks.fail
            method(check_id, f"Could not run {executable}: {exc}")
            continue
        output = (completed.stdout or completed.stderr).strip().splitlines()
        version = output[0][:300] if output else "version output unavailable"
        if completed.returncode not in allowed:
            method = checks.warn if optional else checks.fail
            method(
                check_id,
                f"{executable} version probe failed with status {completed.returncode}: {version}",
            )
        else:
            checks.pass_(check_id, f"{executable}: {version}", details={"path": path})


def _read_fai(path):
    records = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                records.append((fields[0], int(fields[1])))
    return records


def _read_dict(path):
    records = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith("@SQ\t"):
                continue
            values = {}
            for field in line.rstrip("\n").split("\t")[1:]:
                if ":" in field:
                    key, value = field.split(":", 1)
                    values[key] = value
            if "SN" in values and "LN" in values:
                records.append((values["SN"], int(values["LN"])))
    return records


def _check_reference_companions(checks, cfg, steps, skip_picard_metrics):
    ref = cfg.get("reference_data", {})
    fasta_value = ref.get("genome_fa")
    chrom_sizes = ref.get("chrom_sizes")
    if not steps.intersection({2, 3, 4}):
        return None
    if not fasta_value:
        checks.fail("reference.fasta", "No reference FASTA is configured.")
        return None
    fasta = Path(fasta_value)
    if not fasta.is_file():
        checks.fail("reference.fasta", f"Reference FASTA not found: {fasta}")
        return None
    checks.pass_("reference.fasta", f"Reference FASTA exists: {fasta}")

    fai = Path(f"{fasta}.fai")
    if not fai.is_file():
        checks.fail(
            "reference.fai", f"FASTA index not found: {fai}",
            remedy="Run 'cftk init' in an environment with samtools.",
        )
        reference_records = None
    else:
        _, error = _capture_call(validate_fasta_index, fai, chrom_sizes)
        if error:
            checks.fail("reference.fai", error)
            reference_records = None
        else:
            reference_records = _read_fai(fai)
            checks.pass_("reference.fai", f"FASTA index matches chromosome sizes ({len(reference_records)} contigs)")

    markdup_tool = cfg["process"][_STEP_KEYS[3]].get("tool", "sambamba").lower()
    needs_dict = 3 in steps and (markdup_tool == "picard" or not skip_picard_metrics)
    if needs_dict:
        sequence_dict = Path(get_sequence_dictionary_path(str(fasta)))
        if not sequence_dict.is_file():
            checks.fail(
                "reference.dict", f"Picard sequence dictionary not found: {sequence_dict}",
                remedy="Run 'cftk init' in an environment with Picard.",
            )
        else:
            try:
                dict_records = _read_dict(sequence_dict)
            except (OSError, ValueError) as exc:
                checks.fail("reference.dict", f"Invalid Picard sequence dictionary: {exc}")
            else:
                if not dict_records or (reference_records and dict_records != reference_records):
                    checks.fail("reference.dict", "Picard sequence dictionary does not exactly match the FASTA index.")
                else:
                    checks.pass_("reference.dict", f"Picard sequence dictionary has {len(dict_records)} ordered contigs")

    samples = get_all_samples(cfg)
    alignment_tool = cfg["process"][_STEP_KEYS[2]].get("tool", "bwameth").lower()
    needs_bwameth = (
        2 in steps
        and alignment_tool == "bwameth"
        and any(sample.get("input_type", "").lower() == "fastq" for sample in samples)
    )
    if needs_bwameth:
        converted = f"{fasta}.bwameth.c2t"
        missing = [
            f"{converted}{suffix}"
            for suffix in _BWA_SUFFIXES
            if not Path(f"{converted}{suffix}").is_file()
        ]
        if missing:
            checks.fail(
                "reference.bwameth",
                f"bwa-meth index is incomplete ({len(missing)} files missing).",
                remedy="Run 'cftk init' in an environment with bwa-meth and BWA.",
                details={"missing": missing},
            )
        else:
            checks.pass_("reference.bwameth", "bwa-meth converted reference and BWA indexes are complete")
    return reference_records


def _check_profile(checks, raw, config_path):
    if raw.get("schema_version") != SCHEMA_VERSION:
        checks.warn(
            "reference.profile",
            "Legacy config has no versioned profile manifest; individual reference paths are checked only.",
            remedy="Reinitialize as schema-v2 to pin reference component hashes.",
        )
        return None
    profile, error = _capture_call(
        resolve_reference_profile,
        raw,
        config_path,
        acquire_references=False,
        verify_checksums=True,
        validate_compatibility=True,
    )
    if error:
        checks.fail(
            "reference.profile", error,
            remedy="Restore the installed profile or rerun 'cftk init'; doctor never downloads or repairs it.",
        )
        return None
    expected_assay = raw.get("assay", "twist_human_methylome")
    expected_genome = raw.get("genome", "hg38")
    if profile["assay"] != expected_assay or profile["genome"] != expected_genome:
        checks.fail("reference.profile", "Reference profile assay/genome does not match the project config.")
        return profile
    checks.pass_(
        "reference.profile",
        f"Profile {profile['profile_id']} {profile['version']} passed full hashes and compatibility checks",
    )
    return profile


def _check_target_bed(checks, cfg, override):
    if override:
        target_bed = Path(override).expanduser().resolve()
        source = "command-line override"
    elif cfg.get("reference_data", {}).get("target_bed"):
        target_bed = Path(cfg["reference_data"]["target_bed"])
        source = "reference profile/config"
    else:
        from process import DEFAULT_TARGET_BED

        target_bed = Path(DEFAULT_TARGET_BED)
        source = "source-checkout fallback"
    _, error = _capture_call(
        validate_target_bed,
        target_bed,
        cfg["reference_data"].get("chrom_sizes"),
    )
    if error:
        checks.fail(
            "reference.target_bed",
            error,
            remedy=(
                "Use a schema-v2 assay profile or pass --target-bed with a BED "
                "compatible with the configured reference."
            ),
        )
    else:
        checks.pass_(
            "reference.target_bed",
            f"Picard covered-target BED is compatible ({source}): {target_bed}",
        )


def _check_lock(checks, raw, config_path, profile):
    lock_path = config_path.parent / "cftk.lock.json"
    if raw.get("schema_version") != SCHEMA_VERSION:
        checks.warn("project.lock", "Legacy config does not use cftk.lock.json.")
        return
    if not lock_path.is_file():
        checks.fail(
            "project.lock", f"Project lock not found: {lock_path}",
            remedy="Run 'cftk init' to validate the project and recreate its lock.",
        )
        return
    try:
        lock = json.loads(lock_path.read_text())
        sample_sheet = Path(raw["samples"]).expanduser()
        if not sample_sheet.is_absolute():
            sample_sheet = config_path.parent / sample_sheet
        problems = []
        if lock.get("lock_version") != 1:
            problems.append("lock version")
        if lock.get("schema_version") != SCHEMA_VERSION:
            problems.append("schema version")
        if lock.get("project_config_sha256") != sha256_file(config_path):
            problems.append("project config hash")
        if lock.get("sample_sheet", {}).get("sha256") != sha256_file(sample_sheet):
            problems.append("sample-sheet hash")
        locked_profile = lock.get("reference_profile", {})
        if profile:
            expected = {
                "id": profile["profile_id"],
                "version": profile["version"],
                "assay": profile["assay"],
                "genome": profile["genome"],
                "manifest_sha256": profile["manifest_sha256"],
            }
            problems.extend(
                f"reference {key}"
                for key, value in expected.items()
                if locked_profile.get(key) != value
            )
            if locked_profile.get("components") != profile["component_hashes"]:
                problems.append("reference component hashes")
            expected_acquisition = {
                "mode": profile.get("acquisition", {}).get("mode", "local")
            }
            registry_hash = profile.get("acquisition", {}).get(
                "registry_entry_sha256"
            )
            if registry_hash:
                expected_acquisition["registry_entry_sha256"] = registry_hash
            if locked_profile.get("acquisition") != expected_acquisition:
                problems.append("reference acquisition provenance")
    except (KeyError, OSError, json.JSONDecodeError) as exc:
        checks.fail("project.lock", f"Project lock is invalid: {exc}")
        return
    if problems:
        checks.fail(
            "project.lock",
            "Project lock is inconsistent: " + ", ".join(problems),
            remedy="Review intentional edits, then rerun 'cftk init' to recreate the lock.",
        )
    else:
        checks.pass_("project.lock", "Config, sample sheet, manifest, and component hashes match the project lock")


def _fastq_record(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt") as handle:
        lines = [handle.readline().rstrip("\r\n") for _ in range(4)]
    if not lines[0].startswith("@") or not lines[2].startswith("+"):
        raise ValueError("first record is not four-line FASTQ")
    if not lines[1] or len(lines[1]) != len(lines[3]):
        raise ValueError("first record has unequal sequence and quality lengths")
    token = lines[0][1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def _check_fastq_pair(checks, sample, r1, r2):
    check_id = f"input.fastq.{sample['name']}"
    paths = [Path(r1), Path(r2)]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        checks.fail(check_id, f"FASTQ input is missing: {missing}")
        return
    try:
        names = [_fastq_record(path) for path in paths]
    except (OSError, EOFError, ValueError) as exc:
        checks.fail(check_id, f"FASTQ first-record validation failed: {exc}")
        return
    if names[0] != names[1]:
        checks.fail(check_id, f"FASTQ mate identifiers disagree: {names[0]!r} vs {names[1]!r}")
        return
    checks.pass_(check_id, f"Readable paired FASTQ with matching first record: {names[0]}")


def _check_bam(checks, sample, bam_value, reference_records):
    name = sample["name"]
    prefix = f"input.bam.{name}"
    bam = Path(bam_value) if bam_value else None
    if not bam or not bam.is_file():
        checks.fail(prefix, f"BAM input not found: {bam_value or '<not resolved>'}")
        return
    samtools = shutil.which("samtools")
    if samtools:
        try:
            quickcheck = subprocess.run(
                [samtools, "quickcheck", "-v", str(bam)],
                capture_output=True,
                text=True,
                timeout=15,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            checks.fail(prefix, f"samtools quickcheck could not run: {exc}")
            return
        if quickcheck.returncode:
            checks.fail(prefix, f"samtools quickcheck failed: {(quickcheck.stdout or quickcheck.stderr).strip()}")
            return
    try:
        import pysam
        with pysam.AlignmentFile(bam, "rb") as alignment:
            header = alignment.header.to_dict()
            observed = [(item["SN"], int(item["LN"])) for item in header.get("SQ", [])]
    except (ImportError, OSError, ValueError) as exc:
        checks.fail(prefix, f"BAM could not be opened: {exc}")
        return
    checks.pass_(prefix, f"BAM opens successfully: {bam}")

    if header.get("HD", {}).get("SO") != "coordinate":
        checks.fail(f"{prefix}.sorting", "BAM header does not declare coordinate sorting.")
    else:
        checks.pass_(f"{prefix}.sorting", "BAM header declares coordinate sorting")

    indexes = [
        Path(f"{bam}.bai"),
        bam.with_suffix(".bai"),
        Path(f"{bam}.csi"),
        bam.with_suffix(".csi"),
    ]
    index = next((path for path in indexes if path.is_file()), None)
    if index is None:
        checks.fail(f"{prefix}.index", "BAM index is missing.")
    elif index.stat().st_mtime < bam.stat().st_mtime:
        checks.fail(
            f"{prefix}.index", f"BAM index is older than the BAM: {index}",
            remedy="Regenerate the BAM index with the same samtools environment.",
        )
    else:
        checks.pass_(f"{prefix}.index", f"BAM index exists and is not stale: {index}")

    if reference_records is None or observed != reference_records:
        details = {
            "bam_contigs": len(observed),
            "reference_contigs": len(reference_records or []),
            "first_bam_contig": list(observed[0]) if observed else None,
            "first_reference_contig": list(reference_records[0]) if reference_records else None,
        }
        checks.fail(
            f"{prefix}.dictionary",
            "BAM sequence dictionary does not exactly match the configured reference order and lengths.",
            remedy="Use the exact reference profile that was used for alignment.",
            details=details,
        )
    else:
        checks.pass_(f"{prefix}.dictionary", f"BAM dictionary exactly matches {len(observed)} reference contigs")

    read_groups = header.get("RG", [])
    if not read_groups:
        checks.fail(f"{prefix}.read_group", "BAM header has no read-group metadata.")
    else:
        incomplete = [rg.get("ID", "<unknown>") for rg in read_groups if not all(rg.get(key) for key in ("ID", "SM", "LB", "PL"))]
        if incomplete:
            checks.warn(f"{prefix}.read_group", f"Read groups lack one or more ID/SM/LB/PL fields: {incomplete}")
        else:
            checks.pass_(f"{prefix}.read_group", f"{len(read_groups)} read group(s) include ID/SM/LB/PL")

    provenance = " ".join(
        str(value)
        for record in header.get("PG", [])
        for value in (record.get("ID", ""), record.get("PN", ""), record.get("CL", ""))
    ).lower()
    if "markdup" in provenance or "markduplicates" in provenance:
        checks.pass_(f"{prefix}.duplicates", "BAM header records duplicate-marking provenance")
    else:
        checks.warn(
            f"{prefix}.duplicates",
            "BAM header does not establish duplicate-marking provenance.",
            remedy="Confirm whether duplicate marking was completed before methylation extraction.",
        )


def _trimmed_pair(sample, paths):
    ext = "fq.gz" if sample["r1"].endswith(".gz") else "fq"
    name = sample["name"]
    candidates = []
    for left, right in (("R1", "R2"), ("val_1", "val_2")):
        candidates.append(
            (
                Path(paths["trimming"]) / f"{name}_{left}.{ext}",
                Path(paths["trimming"]) / f"{name}_{right}.{ext}",
            )
        )
    return next((pair for pair in candidates if pair[0].is_file() and pair[1].is_file()), candidates[0])


def _check_inputs(checks, cfg, steps, reference_records):
    paths = get_work_paths(cfg)
    for sample in get_all_samples(cfg):
        input_type = sample.get("input_type", "").lower()
        if input_type == "bam":
            if steps.intersection({3, 4}):
                _check_bam(checks, sample, sample.get("bam"), reference_records)
            continue
        if input_type != "fastq":
            checks.fail(f"input.{sample.get('name', 'unknown')}", f"Unsupported input type: {input_type!r}")
            continue
        if 1 in steps:
            _check_fastq_pair(checks, sample, sample.get("r1"), sample.get("r2"))
        elif 2 in steps:
            _check_fastq_pair(checks, sample, *_trimmed_pair(sample, paths))
        if 3 in steps and 2 not in steps:
            _check_bam(checks, sample, get_bam(sample, paths), reference_records)
        elif 4 in steps and not steps.intersection({2, 3}):
            _check_bam(checks, sample, get_bam(sample, paths), reference_records)


def _check_output(checks, cfg):
    output = Path(cfg["output_dir"]).expanduser().resolve()
    probe_dir = output
    while not probe_dir.exists() and probe_dir != probe_dir.parent:
        probe_dir = probe_dir.parent
    try:
        with tempfile.NamedTemporaryFile(prefix=".cftk-doctor-", dir=probe_dir):
            pass
    except OSError as exc:
        checks.fail("output.writable", f"Output location is not writable via {probe_dir}: {exc}")
        return
    checks.pass_("output.writable", f"Output location is writable via {probe_dir}")
    usage = shutil.disk_usage(probe_dir)
    free_gib = usage.free / (1024 ** 3)
    if free_gib < 10:
        checks.warn("output.capacity", f"Only {free_gib:.1f} GiB is free at {probe_dir}; core processing may exhaust it.")
    else:
        checks.pass_("output.capacity", f"{free_gib:.1f} GiB free at {probe_dir}")


def _check_resources(checks, cfg, steps, parallel_override=None):
    process = cfg.get("process", {})
    requested_parallel = parallel_override or process.get("parallel_samples", 1)
    sample_count = len(get_all_samples(cfg))
    scheduler = detect_scheduler_allocation()
    plans = []
    try:
        for step in sorted(steps):
            total = process[_STEP_KEYS[step]].get("params", {}).get("cores", 20)
            ensure_scheduler_capacity(total, scheduler)
            plans.append({
                "step": step,
                **plan_parallelism(total, requested_parallel, sample_count),
            })
    except (KeyError, ValueError) as exc:
        checks.fail(
            "resource.cpu_budget",
            f"Invalid CPU resource plan: {exc}",
            remedy=(
                "Set process.cores to the total allocated CPUs and keep "
                "parallel_samples at or below that budget. Under Slurm, request "
                "at least the same --cpus-per-task value."
            ),
            details={"scheduler": scheduler, "plans": plans},
        )
        return
    allocation = scheduler.get("allocated_cores")
    allocation_text = (
        f"; scheduler allocation {allocation} via {scheduler['variable']}"
        if allocation is not None else "; no scheduler allocation detected"
    )
    checks.pass_(
        "resource.cpu_budget",
        f"CPU budget is valid for {requested_parallel} parallel sample(s)"
        f"{allocation_text}",
        details={"scheduler": scheduler, "plans": plans},
    )


def run_doctor(args):
    """Run all requested diagnostics without acquiring or repairing resources."""
    checks = _Checks()
    steps = set(getattr(args, "step", None) or [1, 2, 3, 4])
    if getattr(args, "analysis_only", False):
        steps = set()
    analysis_stages = list(getattr(args, "analysis_stages", None) or [])
    config_path = Path(args.config).expanduser().resolve()
    _check_runtime(checks)
    raw = _load_raw_config(config_path, checks)
    cfg = None
    profile = None
    reference_records = None
    if raw is not None:
        cfg, error = _capture_call(
            load_config,
            config_path,
            acquire_references=False,
            verify_profile=False,
            validate_profile=False,
        )
        if error:
            checks.fail("project.validation", error)
        else:
            checks.pass_("project.validation", f"Configuration resolves {len(get_all_samples(cfg))} samples")
        profile = _check_profile(checks, raw, config_path)
        _check_lock(checks, raw, config_path, profile)
    if cfg is not None:
        if getattr(args, "analysis_preset", None) or analysis_stages:
            try:
                from analysis_workflow import resolve_stages
                analysis_stages = list(resolve_stages(
                    getattr(args, "analysis_preset", None) or "auto",
                    cfg,
                    analysis_stages or None,
                ))
            except (RuntimeError, ValueError) as exc:
                checks.fail("analysis.plan", str(exc))
                analysis_stages = []
        reference_records = _check_reference_companions(
            checks,
            cfg,
            steps,
            getattr(args, "skip_picard_metrics", False),
        )
        if 3 in steps and not getattr(args, "skip_picard_metrics", False):
            _check_target_bed(checks, cfg, getattr(args, "target_bed", None))
        _check_tools(checks, cfg, steps, getattr(args, "skip_picard_metrics", False))
        _check_inputs(checks, cfg, steps, reference_records)
        _check_output(checks, cfg)
        _check_resources(checks, cfg, steps, getattr(args, "parallel", None))
        if analysis_stages:
            from analysis_workflow import analysis_doctor_checks
            analysis_doctor_checks(
                checks,
                cfg,
                analysis_stages,
                parallel_override=getattr(args, "parallel", None),
                fragmentomics_scope=getattr(args, "fragmentomics_scope", None),
            )
    failures = sum(item["status"] == "FAIL" for item in checks.items)
    warnings = sum(item["status"] == "WARN" for item in checks.items)
    status = "FAIL" if failures else ("WARN" if warnings else "PASS")
    return {
        "report_version": 1,
        "status": status,
        "exit_code": 1 if failures else 0,
        "config": str(config_path),
        "steps": sorted(steps),
        "analysis_preset": getattr(args, "analysis_preset", None),
        "analysis_stages": analysis_stages,
        "fragmentomics_scope": getattr(args, "fragmentomics_scope", None),
        "summary": {
            "pass": sum(item["status"] == "PASS" for item in checks.items),
            "warn": warnings,
            "fail": failures,
        },
        "checks": checks.items,
    }


def render_human(report):
    lines = [
        f"CFTK doctor: {report['status']} (steps {','.join(map(str, report['steps']))})",
        f"Config: {report['config']}",
        "",
    ]
    if report.get("analysis_stages"):
        lines.insert(1, "Analysis: " + ",".join(report["analysis_stages"]))
    for check in report["checks"]:
        lines.append(f"[{check['status']}] {check['id']}: {check['summary']}")
        if check.get("remedy"):
            lines.append(f"       Remedy: {check['remedy']}")
    summary = report["summary"]
    lines.extend(["", f"Summary: {summary['pass']} PASS, {summary['warn']} WARN, {summary['fail']} FAIL"])
    return "\n".join(lines)
