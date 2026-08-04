"""Validate configuration, prepare the genome reference, and resolve paths."""

import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from reference_profiles import (
    DEFAULT_ASSAY,
    DEFAULT_GENOME,
    DEFAULT_PROFILE_ID,
    acquire_reference_profile,
    load_chrom_sizes,
    load_reference_profile,
    managed_profile_available,
    sha256_file,
)
from util import disp

REQUIRED_TOP  = ["project_name", "output_dir", "comparison", "samples",
                  "reference_data", "process", "analysis"]
REQUIRED_REF  = ["genome_fa", "genome_2bit", "chrom_sizes"]
PROCESS_STEPS = ["step1_trimming", "step2_alignment",
                 "step3_markdup", "step4_methylation"]
SCHEMA_VERSION = 2
SAMPLE_SHEET_COLUMNS = ("sample", "group", "role", "input_type", "r1", "r2", "bam")

_SAFE_CHARS = set("abcdefghijklmnopqrstuvwxyz"
                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                  "0123456789-_.")


def load_config(
    config_path="./cftk_init.json",
    *,
    acquire_references=True,
    verify_profile=False,
    validate_profile=False,
):
    if not os.path.exists(config_path):
        sys.exit(
            f"[cftk] ERROR: config not found: {config_path}\n"
            "Run 'cftk init' in the project directory to create it."
        )
    with open(config_path) as f:
        raw = json.load(f)
    raw = _strip_comments(raw)
    cfg = (
        resolve_schema_v2(
            raw,
            config_path,
            acquire_references=acquire_references,
            verify_profile=verify_profile,
            validate_profile=validate_profile,
        )
        if raw.get("schema_version") == 2
        else raw
    )
    errors = _validate(cfg)
    if errors:
        disp("ERROR: invalid cftk_init.json:")
        for e in errors:
            disp(f"  x {e}")
        sys.exit(1)
    disp(f"Config loaded: {config_path}")
    return cfg


def _resolve_relative(path, base_dir):
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = Path(base_dir) / candidate
    return str(candidate.resolve())


def load_sample_sheet(sample_sheet, *, require_files=True):
    """Load the strict schema-v2 sample sheet and preserve its row order."""
    sheet_path = Path(sample_sheet).expanduser().resolve()
    if not sheet_path.is_file():
        sys.exit(f"[init] ERROR: sample sheet not found: {sheet_path}")
    try:
        with sheet_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = tuple(reader.fieldnames or ())
            missing_columns = [c for c in SAMPLE_SHEET_COLUMNS if c not in columns]
            if missing_columns:
                sys.exit(
                    f"[init] ERROR: sample sheet is missing columns: {missing_columns}."
                )
            unexpected_columns = [c for c in columns if c not in SAMPLE_SHEET_COLUMNS]
            if unexpected_columns:
                sys.exit(
                    "[init] ERROR: sample sheet has unexpected columns: "
                    f"{unexpected_columns}."
                )
            rows = list(reader)
    except OSError as exc:
        sys.exit(f"[init] ERROR: could not read sample sheet: {exc}")

    if not rows:
        sys.exit("[init] ERROR: sample sheet contains no samples.")

    samples = {}
    group_roles = {}
    role_groups = {"control": set(), "case": set()}
    seen_names = set()
    seen_input_paths = set()
    for line_number, row in enumerate(rows, start=2):
        name = (row.get("sample") or "").strip()
        group = (row.get("group") or "").strip()
        role = (row.get("role") or "").strip().lower()
        input_type = (row.get("input_type") or "").strip().lower()
        if not name or not group:
            sys.exit(
                f"[init] ERROR: sample sheet line {line_number} requires sample and group."
            )
        if any(char not in _SAFE_CHARS for char in name):
            sys.exit(f"[init] ERROR: invalid sample name {name!r} on line {line_number}.")
        if "_vs_" in group or any(char not in _SAFE_CHARS for char in group):
            sys.exit(f"[init] ERROR: invalid group name {group!r} on line {line_number}.")
        if name in seen_names:
            sys.exit(f"[init] ERROR: duplicate sample name {name!r} in sample sheet.")
        seen_names.add(name)
        if role not in role_groups:
            sys.exit(
                f"[init] ERROR: sample {name!r} role must be 'control' or 'case'."
            )
        if group in group_roles and group_roles[group] != role:
            sys.exit(f"[init] ERROR: group {group!r} has conflicting biological roles.")
        group_roles[group] = role
        role_groups[role].add(group)

        entry = {"name": name, "input_type": input_type}
        if input_type == "fastq":
            r1 = (row.get("r1") or "").strip()
            r2 = (row.get("r2") or "").strip()
            bam = (row.get("bam") or "").strip()
            if not r1 or not r2:
                sys.exit(f"[init] ERROR: fastq sample {name!r} requires both r1 and r2.")
            if bam:
                sys.exit(
                    f"[init] ERROR: fastq sample {name!r} must leave bam empty."
                )
            entry["r1"] = _resolve_relative(r1, sheet_path.parent)
            entry["r2"] = _resolve_relative(r2, sheet_path.parent)
            input_paths = (entry["r1"], entry["r2"])
        elif input_type == "bam":
            bam = (row.get("bam") or "").strip()
            if not bam:
                sys.exit(f"[init] ERROR: bam sample {name!r} requires bam.")
            if (row.get("r1") or "").strip() or (row.get("r2") or "").strip():
                sys.exit(
                    f"[init] ERROR: bam sample {name!r} must leave r1 and r2 empty."
                )
            entry["bam"] = _resolve_relative(bam, sheet_path.parent)
            input_paths = (entry["bam"],)
        else:
            sys.exit(
                f"[init] ERROR: sample {name!r} input_type must be fastq or bam."
            )
        if require_files:
            missing = [path for path in input_paths if not os.path.isfile(path)]
            if missing:
                sys.exit(
                    f"[init] ERROR: input file(s) not found for sample {name!r}: {missing}"
                )
        reused = [path for path in input_paths if path in seen_input_paths]
        if reused:
            sys.exit(
                f"[init] ERROR: input file(s) are reused by multiple samples: {reused}"
            )
        if len(set(input_paths)) != len(input_paths):
            sys.exit(f"[init] ERROR: sample {name!r} reuses the same input file.")
        seen_input_paths.update(input_paths)
        samples.setdefault(group, []).append(entry)

    for role, groups in role_groups.items():
        if len(groups) != 1:
            sys.exit(
                f"[init] ERROR: sample sheet must define exactly one group for "
                f"role '{role}', got {sorted(groups)}."
            )
    control_group = next(iter(role_groups["control"]))
    case_group = next(iter(role_groups["case"]))
    if len(samples) != 2:
        sys.exit("[init] ERROR: the first schema-v2 release supports exactly two groups.")
    return {
        "samples": samples,
        "group_roles": group_roles,
        "control_group": control_group,
        "case_group": case_group,
        "sample_sheet": str(sheet_path),
    }


def _default_analysis(samples, control_group, case_group):
    return {
        "qc": {"params": {"fragment": 167, "step_size": 2000}},
        "power": {"params": {
            "sample_size": 100,
            "effect_size": 0.1,
            "depth": [10, 20, 50],
            "ratio": 1.0,
            "plot_threshold": 0.8,
            "step_size": 10000,
        }},
        "diff": {"params": {
            "modalities": ["cpg", "occupancy", "wps"],
            "colors": ["#4575b4", "#d73027"],
            "top_n_heatmap": 500,
        }},
        "dmr": {
            "tool": "metilene",
            "params": {"cores": 20, "q_thr": 0.05, "top_n": 20, "extra_args": ""},
            "samples": {
                control_group: [s["name"] for s in samples[control_group]],
                case_group: [s["name"] for s in samples[case_group]],
            },
        },
        "frag": {
            "occupancy": {"tool": "danpos", "params": {"extra_args": "--paired 1 -u 0 -c 1000000"}},
            "wps": {"params": {"wps_window": 120, "wps_step": 10, "min_frag": 100, "max_frag": 220}},
            "delfi": {"tool": "finaletoolkit", "params": {"mapq": 30, "window": 20, "extra_args": ""}},
            "end_motif": {"tool": "finaletoolkit", "params": {"kmer": 4, "min_frag": 100, "max_frag": 220, "mapq": 30, "extra_args": ""}},
            "cleavage": {"tool": "finaletoolkit", "params": {"min_frag": 100, "max_frag": 220, "mapq": 30, "window": 20, "upstream": 1500, "downstream": 1500, "extra_args": ""}},
        },
        "mesa": {"params": {
            "modalities": ["cpg", "occupancy", "wps"],
            "clf": [1, 2, 3],
            "feature_size": 100,
            "subset": 0.1,
            "repeat": 3,
        }},
    }


def resolve_reference_profile(
    raw,
    config_path,
    *,
    acquire_references=True,
    verify_checksums=False,
    validate_compatibility=False,
):
    """Resolve a schema-v2 profile, optionally without managed acquisition."""
    config_dir = Path(config_path).expanduser().resolve().parent
    reference_root_value = (
        os.environ.get("CFTK_REFERENCE_ROOT")
        or raw.get("reference_root")
        or str(Path.home() / ".cache" / "cftk" / "references")
    )
    reference_root = _resolve_relative(reference_root_value, config_dir)
    profile_spec = raw.get("reference_profile", DEFAULT_PROFILE_ID)
    if isinstance(profile_spec, str):
        profile_id, profile_version = profile_spec, None
    elif isinstance(profile_spec, dict):
        profile_id = profile_spec.get("id", DEFAULT_PROFILE_ID)
        profile_version = profile_spec.get("version")
    else:
        sys.exit("[init] ERROR: reference_profile must be a string or object.")
    if acquire_references:
        profile = acquire_reference_profile(
            mode=raw.get("reference_mode", "local"),
            reference_root=reference_root,
            profile_id=profile_id,
            version=profile_version,
            verify_checksums=verify_checksums,
            validate_compatibility=validate_compatibility,
        )
    else:
        profile = load_reference_profile(
            reference_root,
            profile_id,
            profile_version,
            verify_checksums=verify_checksums,
            validate_compatibility=validate_compatibility,
        )
    return profile


def resolve_schema_v2(
    raw,
    config_path,
    *,
    acquire_references=True,
    verify_profile=False,
    validate_profile=False,
):
    """Expand a compact schema-v2 config into the established nested contract."""
    config_dir = Path(config_path).expanduser().resolve().parent
    required = ("project_name", "samples")
    missing = [key for key in required if not raw.get(key)]
    if missing:
        sys.exit(f"[init] ERROR: schema-v2 config is missing fields: {missing}.")
    assay = raw.get("assay", DEFAULT_ASSAY)
    genome = raw.get("genome", DEFAULT_GENOME)
    sample_sheet = _resolve_relative(raw["samples"], config_dir)
    sample_data = load_sample_sheet(sample_sheet)
    profile = resolve_reference_profile(
        raw,
        config_path,
        acquire_references=acquire_references,
        verify_checksums=verify_profile,
        validate_compatibility=validate_profile,
    )
    if profile["assay"] != assay or profile["genome"] != genome:
        sys.exit(
            "[init] ERROR: reference profile assay/genome does not match the "
            "project configuration."
        )

    component_keys = (
        "genome_fa", "genome_2bit", "chrom_sizes", "target_bed",
        "tss_pas_bed", "ctcf_bed", "blacklist", "gap", "bins", "cpg_std",
    )
    reference_data = {
        key: profile["components"][key]
        for key in component_keys
        if key in profile["components"]
    }
    reference_data["genome_build"] = genome

    process_settings = raw.get("process", {})
    if not isinstance(process_settings, dict):
        sys.exit("[init] ERROR: schema-v2 process must be an object.")
    cores = process_settings.get("cores", 20)
    parallel = process_settings.get("parallel_samples", 1)
    min_depth = process_settings.get("min_depth", 10)
    for name, value in (
        ("cores", cores),
        ("parallel_samples", parallel),
        ("min_depth", min_depth),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            sys.exit(
                f"[init] ERROR: schema-v2 process.{name} must be a positive integer."
            )
    process = {
        "parallel_samples": parallel,
        "step1_trimming": {"tool": "trim_galore", "params": {"cores": cores, "extra_args": ""}},
        "step2_alignment": {"tool": "bwameth", "params": {"cores": cores, "extra_args": ""}},
        "step3_markdup": {"tool": "sambamba", "params": {"cores": cores, "extra_args": ""}},
        "step4_methylation": {"tool": "methyldackel", "params": {"cores": cores, "min_depth": min_depth, "extra_args": ""}},
    }
    output_dir = _resolve_relative(raw.get("output_dir", "."), config_dir)
    control_group = sample_data["control_group"]
    case_group = sample_data["case_group"]
    cfg = {
        "schema_version": SCHEMA_VERSION,
        "project_name": raw["project_name"],
        "output_dir": output_dir,
        "comparison": f"{control_group}_vs_{case_group}",
        "control_group": control_group,
        "case_group": case_group,
        "group_roles": sample_data["group_roles"],
        "samples": sample_data["samples"],
        "assay": assay,
        "reference_profile": {"id": profile["profile_id"], "version": profile["version"]},
        "reference_data": reference_data,
        "process": process,
        "analysis": _default_analysis(sample_data["samples"], control_group, case_group),
    }
    return cfg


def _strip_comments(obj):
    if isinstance(obj, dict):
        return {k: _strip_comments(v) for k, v in obj.items()
                if not k.startswith("_comment")}
    if isinstance(obj, list):
        return [_strip_comments(v) for v in obj]
    return obj


def _validate(cfg):
    errors = []

    for k in REQUIRED_TOP:
        if k not in cfg:
            errors.append(f"Missing required top-level key: '{k}'")

    if not cfg.get("output_dir", "").strip():
        errors.append("'output_dir' must not be empty.")

    comp = cfg.get("comparison", "")
    if "_vs_" not in comp:
        errors.append("'comparison' must be formatted as 'GroupA_vs_GroupB'.")
    else:
        ga, gb = comp.split("_vs_", 1)
        samples = cfg.get("samples", {})
        if ga not in samples:
            errors.append(f"comparison group_a '{ga}' not found in samples.")
        if gb not in samples:
            errors.append(f"comparison group_b '{gb}' not found in samples.")

    samples = cfg.get("samples", {})
    if len(samples) != 2:
        errors.append(f"'samples' must define exactly 2 groups, got {len(samples)}.")

    for grp, members in samples.items():
        if not isinstance(members, list) or len(members) == 0:
            errors.append(f"Group '{grp}' must be a non-empty list.")
            continue
        for s in members:
            name  = s.get("name", "")
            if not name:
                errors.append(f"A sample in group '{grp}' is missing 'name'.")
            else:
                # sample name: alphanumeric, hyphen, underscore, dot only
                bad = [c for c in name if c not in _SAFE_CHARS]
                if bad:
                    errors.append(
                        f"Sample '{name}' in group '{grp}': name contains "
                        f"invalid characters {bad}. "
                        "Only letters, digits, - _ . are allowed."
                    )
            itype = s.get("input_type", "").lower()
            if itype not in ("fastq", "bam"):
                errors.append(
                    f"Sample '{name}': 'input_type' must be 'fastq' or 'bam'."
                )
            if itype == "fastq" and (not s.get("r1") or not s.get("r2")):
                errors.append(
                    f"Sample '{name}': fastq input requires 'r1' and 'r2'."
                )
            if itype == "bam" and not s.get("bam"):
                errors.append(f"Sample '{name}': bam input requires 'bam'.")

    ref = cfg.get("reference_data", {})
    for k in REQUIRED_REF:
        if not ref.get(k):
            errors.append(f"reference_data.{k} is required.")

    proc = cfg.get("process", {})
    for step in PROCESS_STEPS:
        if step not in proc:
            errors.append(f"process.{step} is missing.")

    methylation_params = proc.get("step4_methylation", {}).get("params", {})
    min_depth = methylation_params.get("min_depth", 10)
    if (
        isinstance(min_depth, bool)
        or not isinstance(min_depth, int)
        or min_depth < 1
    ):
        errors.append(
            "process.step4_methylation.params.min_depth must be a positive integer."
        )

    return errors


def get_all_samples(cfg):
    """Flat list of all samples with 'group' injected."""
    out = []
    for grp, members in cfg["samples"].items():
        for s in members:
            entry = dict(s)
            entry["group"] = grp
            out.append(entry)
    return out


def get_group_names(cfg):
    """Return (group_a, group_b) from comparison field."""
    ga, gb = cfg["comparison"].split("_vs_", 1)
    return ga, gb


def get_work_paths(cfg):
    """
    All standard subdirectory paths derived from output_dir.
    Directories are NOT created here — each command creates its own at runtime.
    """
    wd = cfg["output_dir"]
    r  = os.path.join(wd, "results")
    p  = os.path.join(r, "1_process")
    f  = os.path.join(r, "4_fragmentomics")
    return {
        "work_dir":      wd,
        "results":       r,
        "power":         os.path.join(r, "0_power"),
        "process":       p,
        "qc":            os.path.join(r, "2_qc"),
        "differential":  os.path.join(r, "3_differential"),
        "fragmentomics": f,
        "mesa":          os.path.join(r, "5_mesa"),
        "report":        os.path.join(r, "report"),
        # process step output dirs
        "trimming":      os.path.join(p, "1_trimming"),
        "alignment":     os.path.join(p, "2_alignment"),
        "markdup":       os.path.join(p, "3_markdup"),
        "methylation":   os.path.join(p, "4_methylation"),
        "cpg_matrix":    os.path.join(p, "5_merged_matrix"),
        # fragmentomics sub-analysis output dirs
        "occ_out":       os.path.join(f, "occupancy"),
        "wps_out":       os.path.join(f, "wps"),
        "delfi_out":     os.path.join(f, "delfi"),
        "end_motif_out": os.path.join(f, "end_motif"),
        "cleavage_out":  os.path.join(f, "cleavage"),
    }


def get_bam(sample, paths):
    """Resolve the best available BAM (markdup > alignment > direct bam)."""
    name = sample["name"]
    for c in [
        os.path.join(paths["markdup"],   f"{name}.markdup.bam"),
        os.path.join(paths["alignment"], f"{name}.bam"),
        sample.get("bam", ""),
    ]:
        if c and os.path.exists(c):
            return c
    return sample.get("bam", "")


def get_matrix_path(paths, modality):
    """Canonical location of the merged feature matrix for each modality."""
    mapping = {
        "cpg":       os.path.join(paths["cpg_matrix"], "cpg_matrix.tsv"),
        "occupancy": os.path.join(paths["occ_out"],    "occupancy_matrix.tsv"),
        "wps":       os.path.join(paths["wps_out"],    "wps_matrix.tsv"),
        "delfi":     os.path.join(paths["delfi_out"],  "delfi_matrix.tsv"),
    }
    if modality in mapping:
        return mapping[modality]
    return os.path.join(paths["fragmentomics"], f"{modality}_matrix.tsv")


_FASTQ_RE = re.compile(
    r"^(?P<sample>.+?)(?:_L(?P<lane>[0-9]{3}))?_R(?P<read>[12])"
    r"(?:_[0-9]{3})?\.(?:fastq|fq)(?:\.gz)?$",
    re.IGNORECASE,
)


def discover_sample_inputs(project_dir):
    """Discover only unambiguous single-pair FASTQs or one BAM per sample."""
    project_dir = Path(project_dir).expanduser().resolve()
    search_dirs = [project_dir]
    for name in ("data", "fastq", "fastqs", "bam", "bams", "raw_data"):
        candidate = project_dir / name
        if candidate.is_dir():
            search_dirs.append(candidate)

    fastqs = {}
    bams = {}
    for directory in search_dirs:
        for path in sorted(directory.iterdir()):
            if not path.is_file():
                continue
            match = _FASTQ_RE.match(path.name)
            if match:
                sample = match.group("sample")
                read = match.group("read")
                fastqs.setdefault(sample, {"1": [], "2": [], "lanes": set()})
                fastqs[sample][read].append(path.resolve())
                if match.group("lane"):
                    fastqs[sample]["lanes"].add(match.group("lane"))
            elif path.suffix.lower() == ".bam":
                bams.setdefault(path.stem, []).append(path.resolve())

    rows = []
    for sample in sorted(set(fastqs) | set(bams)):
        if sample in fastqs and sample in bams:
            sys.exit(
                f"[init] ERROR: sample {sample!r} has both FASTQ and BAM inputs."
            )
        if sample in fastqs:
            info = fastqs[sample]
            if len(info["1"]) != 1 or len(info["2"]) != 1 or len(info["lanes"]) > 1:
                sys.exit(
                    f"[init] ERROR: multi-lane or ambiguous FASTQ inputs for "
                    f"sample {sample!r}; create samples.tsv explicitly."
                )
            rows.append({
                "sample": sample,
                "group": "",
                "role": "",
                "input_type": "fastq",
                "r1": str(info["1"][0]),
                "r2": str(info["2"][0]),
                "bam": "",
            })
        else:
            if len(bams[sample]) != 1:
                sys.exit(
                    f"[init] ERROR: ambiguous BAM inputs for sample {sample!r}."
                )
            rows.append({
                "sample": sample,
                "group": "",
                "role": "",
                "input_type": "bam",
                "r1": "",
                "r2": "",
                "bam": str(bams[sample][0]),
            })
    if not rows:
        sys.exit(
            "[init] ERROR: no unambiguous FASTQ pairs or BAM files were found. "
            "Create samples.tsv or pass --sample-sheet PATH."
        )
    return rows


def _portable_path(path, base_dir):
    path = Path(path).resolve()
    base_dir = Path(base_dir).resolve()
    try:
        return str(path.relative_to(base_dir))
    except ValueError:
        return str(path)


def write_sample_sheet_template(rows, out_path):
    """Write discovered paths with blank biological fields for explicit editing."""
    out_path = Path(out_path).expanduser().resolve()
    if out_path.exists():
        sys.exit(f"[init] ERROR: refusing to overwrite existing sample sheet: {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_SHEET_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            portable = dict(row)
            for key in ("r1", "r2", "bam"):
                if portable.get(key):
                    portable[key] = _portable_path(portable[key], out_path.parent)
            writer.writerow(portable)
    disp(f"Sample-sheet template created: {out_path}")
    return str(out_path)


def _atomic_write_json(path, payload):
    path = Path(path)
    temp_path = path.with_name(f".{path.name}.tmp")
    with temp_path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    os.replace(temp_path, path)


def _prompt(label, default=None):
    suffix = f" [{default}]" if default is not None else ""
    value = input(f"{label}{suffix}: ").strip()
    return value or default


def _create_project_config(args):
    config_path = Path(args.config).expanduser().resolve()
    if config_path.exists():
        sys.exit(f"[init] ERROR: refusing to overwrite existing config: {config_path}")
    config_path.parent.mkdir(parents=True, exist_ok=True)
    interactive = sys.stdin.isatty() and not getattr(args, "non_interactive", False)

    sample_sheet = getattr(args, "sample_sheet", None)
    reference_root = getattr(args, "reference_root", None)
    profile_id = getattr(args, "profile", None) or DEFAULT_PROFILE_ID
    reference_mode = getattr(args, "reference_mode", None)
    if reference_mode is None:
        reference_mode = (
            "managed" if managed_profile_available(profile_id) else "local"
        )
    if not interactive and (
        not sample_sheet or (reference_mode == "local" and not reference_root)
    ):
        sys.exit(
            "[init] ERROR: noninteractive creation requires --sample-sheet PATH; "
            "local reference mode also requires --reference-root PATH."
        )

    project_name = getattr(args, "project_name", None) or config_path.parent.name
    output_dir = getattr(args, "output_dir", None) or "."
    profile_version = getattr(args, "profile_version", None)

    if not reference_root:
        reference_root = os.environ.get(
            "CFTK_REFERENCE_ROOT",
            str(Path.home() / ".cache" / "cftk" / "references"),
        )

    if interactive:
        project_name = _prompt("Project name", project_name)
        output_dir = _prompt("Output directory", output_dir)
        default_sheet = sample_sheet or str(config_path.parent / "samples.tsv")
        sample_sheet = _prompt("Sample sheet", default_sheet)
        if reference_mode == "local" and getattr(args, "reference_root", None) is None:
            reference_root = _prompt("Reference root", reference_root)

    sample_sheet_path = Path(sample_sheet).expanduser()
    if not sample_sheet_path.is_absolute():
        sample_sheet_path = config_path.parent / sample_sheet_path
    if not sample_sheet_path.is_file():
        if not interactive:
            sys.exit(f"[init] ERROR: sample sheet not found: {sample_sheet_path}")
        rows = discover_sample_inputs(config_path.parent)
        write_sample_sheet_template(rows, sample_sheet_path)
        sys.exit(
            f"[init] Edit group and role columns in {sample_sheet_path}, then rerun cftk init."
        )
    load_sample_sheet(sample_sheet_path)

    profile = acquire_reference_profile(
        mode=reference_mode,
        reference_root=reference_root,
        profile_id=profile_id,
        version=profile_version,
        verify_checksums=True,
        validate_compatibility=True,
    )
    assay = getattr(args, "assay", None) or DEFAULT_ASSAY
    genome = getattr(args, "genome", None) or DEFAULT_GENOME
    if profile["assay"] != assay or profile["genome"] != genome:
        sys.exit(
            "[init] ERROR: reference profile assay/genome does not match the "
            "requested project configuration."
        )
    compact = {
        "schema_version": SCHEMA_VERSION,
        "project_name": project_name,
        "output_dir": output_dir,
        "assay": assay,
        "genome": genome,
        "samples": _portable_path(sample_sheet_path, config_path.parent),
        "reference_mode": reference_mode,
        "reference_root": str(Path(reference_root).expanduser().resolve()),
        "reference_profile": {
            "id": profile["profile_id"],
            "version": profile["version"],
        },
        "process": {"cores": 20, "parallel_samples": 1, "min_depth": 10},
    }
    _atomic_write_json(config_path, compact)
    disp(f"Project config created: {config_path}")
    write_lockfile(config_path, profile=profile)
    return str(config_path)


def write_lockfile(config_path, *, profile=None):
    """Write portable hashes for a schema-v2 project and its local profile."""
    config_path = Path(config_path).expanduser().resolve()
    raw = _strip_comments(json.loads(config_path.read_text()))
    if raw.get("schema_version") != SCHEMA_VERSION:
        return None
    config_dir = config_path.parent
    sample_sheet = Path(_resolve_relative(raw["samples"], config_dir))
    reference_root = os.environ.get("CFTK_REFERENCE_ROOT") or raw.get("reference_root")
    if not reference_root:
        reference_root = str(Path.home() / ".cache" / "cftk" / "references")
    reference_root = _resolve_relative(reference_root, config_dir)
    profile_spec = raw.get("reference_profile", DEFAULT_PROFILE_ID)
    if isinstance(profile_spec, str):
        profile_id, version = profile_spec, None
    else:
        profile_id = profile_spec.get("id", DEFAULT_PROFILE_ID)
        version = profile_spec.get("version")
    if profile is None:
        profile = acquire_reference_profile(
            mode=raw.get("reference_mode", "local"),
            reference_root=reference_root,
            profile_id=profile_id,
            version=version,
            verify_checksums=True,
            validate_compatibility=True,
        )
    lock = {
        "lock_version": 1,
        "schema_version": SCHEMA_VERSION,
        "project_config_sha256": sha256_file(config_path),
        "sample_sheet": {
            "path": _portable_path(sample_sheet, config_dir),
            "sha256": sha256_file(sample_sheet),
        },
        "reference_profile": {
            "id": profile["profile_id"],
            "version": profile["version"],
            "assay": profile["assay"],
            "genome": profile["genome"],
            "manifest_sha256": profile["manifest_sha256"],
            "components": profile["component_hashes"],
            "acquisition": {
                "mode": profile.get("acquisition", {}).get(
                    "mode", raw.get("reference_mode", "local")
                ),
                **(
                    {
                        "registry_entry_sha256": profile["acquisition"][
                            "registry_entry_sha256"
                        ]
                    }
                    if profile.get("acquisition", {}).get("registry_entry_sha256")
                    else {}
                ),
            },
        },
    }
    lock_path = config_dir / "cftk.lock.json"
    _atomic_write_json(lock_path, lock)
    disp(f"Project lock written: {lock_path}")
    return str(lock_path)


def get_sequence_dictionary_path(reference_fa):
    """Return Picard's conventional companion dictionary path for a FASTA."""
    stem, suffix = os.path.splitext(reference_fa)
    if suffix.lower() in {".fa", ".fasta", ".fna"}:
        return f"{stem}.dict"
    return f"{reference_fa}.dict"


def _run_reference_command(command, label):
    disp(f"CMD [{label}]: {shlex.join(command)}")
    try:
        return subprocess.run(command).returncode
    except FileNotFoundError:
        return 127


def _bwameth_index_exists(reference_fa):
    converted = f"{reference_fa}.bwameth.c2t"
    return all(
        os.path.exists(f"{converted}{suffix}")
        for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa")
    )


def validate_fasta_index(fai_path, chrom_sizes_path):
    """Require FASTA-index contigs and lengths to match chromosome sizes."""
    expected = load_chrom_sizes(chrom_sizes_path)
    observed = {}
    try:
        with open(fai_path) as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 2 or not fields[0]:
                    sys.exit(
                        f"[init] ERROR: invalid FASTA index line {line_number}: {fai_path}"
                    )
                try:
                    length = int(fields[1])
                except ValueError:
                    sys.exit(
                        f"[init] ERROR: invalid FASTA index length on line "
                        f"{line_number}: {fai_path}"
                    )
                if fields[0] in observed:
                    sys.exit(
                        f"[init] ERROR: duplicate FASTA index contig {fields[0]!r}."
                    )
                observed[fields[0]] = length
    except OSError as exc:
        sys.exit(f"[init] ERROR: could not read FASTA index {fai_path}: {exc}")

    all_contigs = list(expected) + [name for name in observed if name not in expected]
    mismatched = [
        name for name in all_contigs if observed.get(name) != expected.get(name)
    ]
    if mismatched:
        name = mismatched[0]
        sys.exit(
            "[init] ERROR: FASTA index does not match chromosome sizes for "
            f"{name}: fai={observed.get(name)!r}, chrom_sizes={expected.get(name)!r}."
        )


def prepare_reference(
    reference_fa, prepare_index=True, prepare_dict=True, chrom_sizes=None
):
    """Prepare bwa-meth, samtools, and Picard indexes for a reference FASTA."""
    reference_fa = os.path.abspath(reference_fa)
    if not os.path.isfile(reference_fa):
        sys.exit(f"[init] ERROR: reference FASTA not found: {reference_fa}")

    if prepare_index:
        if _bwameth_index_exists(reference_fa):
            disp(f"bwa-meth index already present: {reference_fa}.bwameth.c2t")
        else:
            for executable in ("bwameth", "bwameth.py"):
                returncode = _run_reference_command(
                    [executable, "index", reference_fa], f"{executable} index"
                )
                if returncode == 0 and _bwameth_index_exists(reference_fa):
                    break
            else:
                sys.exit(
                    "[init] ERROR: bwa-meth indexing failed or did not create "
                    "the expected index files with either 'bwameth index' or "
                    "'bwameth.py index'."
                )

    fai = f"{reference_fa}.fai"
    if prepare_index:
        if os.path.exists(fai):
            disp(f"FASTA index already present: {fai}")
        elif _run_reference_command(
            ["samtools", "faidx", reference_fa], "samtools faidx"
        ):
            sys.exit("[init] ERROR: samtools faidx failed.")
        if not os.path.isfile(fai):
            sys.exit(f"[init] ERROR: samtools did not create FASTA index: {fai}")
        if chrom_sizes:
            validate_fasta_index(fai, chrom_sizes)

    sequence_dict = get_sequence_dictionary_path(reference_fa)
    if prepare_dict:
        if os.path.exists(sequence_dict):
            disp(f"Sequence dictionary already present: {sequence_dict}")
        elif _run_reference_command(
            [
                "picard", "CreateSequenceDictionary",
                f"R={reference_fa}", f"O={sequence_dict}",
            ],
            "picard CreateSequenceDictionary",
        ):
            sys.exit("[init] ERROR: Picard CreateSequenceDictionary failed.")
        if not os.path.isfile(sequence_dict):
            sys.exit(
                "[init] ERROR: Picard did not create sequence dictionary: "
                f"{sequence_dict}"
            )

    return {
        "reference_fa": reference_fa,
        "bwameth_reference": f"{reference_fa}.bwameth.c2t",
        "fai": fai,
        "sequence_dictionary": sequence_dict,
    }


def init(args):
    """Validate config, summarize the project, and prepare the genome FASTA."""
    created = not os.path.exists(args.config)
    if created:
        _create_project_config(args)
    cfg  = load_config(args.config)
    if cfg.get("schema_version") == SCHEMA_VERSION and not created:
        write_lockfile(args.config)
    ga, gb = get_group_names(cfg)

    disp(f"Project    : {cfg['project_name']}")
    disp(f"Output dir : {cfg['output_dir']}/results/")
    disp(f"Comparison : {ga} (label=0) vs {gb} (label=1)")
    for grp, members in cfg["samples"].items():
        names = [s["name"] for s in members]
        disp(f"  {grp} ({len(names)}): {', '.join(names)}")

    if getattr(args, "skip_reference_prep", False) and (
        getattr(args, "ref_index", False) or getattr(args, "ref_dict", False)
    ):
        sys.exit(
            "[init] ERROR: --skip-reference-prep cannot be combined with "
            "--ref-index or --ref-dict."
        )

    ref = cfg["reference_data"]["genome_fa"]
    if getattr(args, "skip_reference_prep", False):
        disp("Reference preparation skipped by user.")
    else:
        legacy_index = getattr(args, "ref_index", False)
        legacy_dict = getattr(args, "ref_dict", False)
        if legacy_index or legacy_dict:
            prepare_reference(
                ref,
                prepare_index=legacy_index,
                prepare_dict=legacy_dict,
                chrom_sizes=cfg["reference_data"].get("chrom_sizes"),
            )
        else:
            prepare_reference(
                ref, chrom_sizes=cfg["reference_data"].get("chrom_sizes")
            )

    disp("Init complete. Run: cftk process -s 1 2 3 4")
