"""Role-aware, resumable orchestration for CFTK downstream analyses.

The scientific implementations remain in ``analysis`` and ``visualization``.
This module owns only workflow planning, preflight checks, artifact contracts,
provenance, evidence, and safe resume behavior.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

from init import get_all_samples, get_bam, get_matrix_path, get_work_paths, load_config
from resource_planning import detect_scheduler_allocation, plan_parallelism
from util import configure_command_log, disp
from analysis.assay_scope import (
    ScopeError,
    describe_scope,
    resolve_scope,
    scope_artifact_paths,
    scope_metadata_path,
)

import run_workflow as _core


ANALYSIS_RUN_SCHEMA_VERSION = 1
PRESET_NAMES = (
    "auto",
    "descriptive",
    "differential",
    "dmr",
    "fragmentomics",
    "mesa",
    "comparative",
    "all",
    "report",
)
STAGE_CHOICES = (
    "analysis.diff",
    "analysis.dmr",
    "fragmentomics.occupancy",
    "fragmentomics.wps",
    "fragmentomics.delfi",
    "fragmentomics.end_motif",
    "fragmentomics.cleavage",
    "analysis.mesa",
    "analysis.report",
)

_FRAGMENT_STAGES = (
    "fragmentomics.occupancy",
    "fragmentomics.wps",
    "fragmentomics.delfi",
    "fragmentomics.end_motif",
    "fragmentomics.cleavage",
)
_COMPARATIVE_STAGES = {"analysis.diff", "analysis.dmr", "analysis.mesa"}
_FEATURE_MATRIX_STAGES = {
    "occupancy": "fragmentomics.occupancy",
    "wps": "fragmentomics.wps",
}
_DMR_R_PACKAGES = (
    "annotatr",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "GenomicRanges",
    "org.Hs.eg.db",
)
_STAGE_ALIASES = {
    "diff": ("analysis.diff",),
    "differential": ("analysis.diff",),
    "dmr": ("analysis.dmr",),
    "occupancy": ("fragmentomics.occupancy",),
    "wps": ("fragmentomics.wps",),
    "delfi": ("fragmentomics.delfi",),
    "end_motif": ("fragmentomics.end_motif",),
    "end-motif": ("fragmentomics.end_motif",),
    "cleavage": ("fragmentomics.cleavage",),
    "frag": _FRAGMENT_STAGES,
    "fragmentomics": _FRAGMENT_STAGES,
    "mesa": ("analysis.mesa",),
    "report": ("analysis.report",),
}
_PRESETS = {
    "descriptive": _FRAGMENT_STAGES[:2] + ("analysis.report",),
    "differential": ("analysis.diff",),
    "dmr": ("analysis.dmr",),
    "fragmentomics": _FRAGMENT_STAGES,
    "mesa": ("analysis.mesa",),
    "comparative": ("analysis.diff", "analysis.dmr", "analysis.mesa", "analysis.report"),
    "all": (*_FRAGMENT_STAGES, "analysis.diff", "analysis.dmr", "analysis.mesa", "analysis.report"),
    "report": ("analysis.report",),
}
_STAGE_META = {
    "analysis.diff": {
        "name": "Differential methylation, PCA, and visualizations",
        "kind": "diff",
        "comparative": True,
    },
    "analysis.dmr": {
        "name": "DMR calling and annotation",
        "kind": "dmr",
        "comparative": True,
    },
    "fragmentomics.occupancy": {
        "name": "Nucleosome occupancy",
        "kind": "occupancy",
        "comparative": False,
    },
    "fragmentomics.wps": {
        "name": "Windowed protection score",
        "kind": "wps",
        "comparative": False,
    },
    "fragmentomics.delfi": {
        "name": "DELFI-style fragmentomics",
        "kind": "delfi",
        "comparative": False,
    },
    "fragmentomics.end_motif": {
        "name": "Fragment end motifs",
        "kind": "end_motif",
        "comparative": False,
    },
    "fragmentomics.cleavage": {
        "name": "Cleavage profiles",
        "kind": "cleavage",
        "comparative": False,
    },
    "analysis.mesa": {
        "name": "MESA modality modeling and LOOCV",
        "kind": "mesa",
        "comparative": True,
    },
    "analysis.report": {
        "name": "Self-contained HTML report",
        "kind": "report",
        "comparative": False,
    },
}


class AnalysisContractError(RuntimeError):
    """A downstream workflow contract was not satisfied."""


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _new_id(prefix="analysis"):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"{prefix}-{stamp}-{os.urandom(4).hex()}"


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _role_info(cfg):
    """Return role-to-group mappings without guessing from group names."""
    groups = list(cfg.get("samples", {}).keys())
    roles = cfg.get("group_roles", {})
    control = [group for group in groups if roles.get(group) == "control"]
    case = [group for group in groups if roles.get(group) == "case"]
    return {
        "groups": groups,
        "group_roles": {group: roles.get(group) for group in groups},
        "control_groups": control,
        "case_groups": case,
        "control_group": control[0] if len(control) == 1 else None,
        "case_group": case[0] if len(case) == 1 else None,
    }


def comparison_role_errors(cfg):
    """Validate the explicit control/case contract for comparative stages."""
    info = _role_info(cfg)
    errors = []
    if len(info["groups"]) != 2:
        errors.append(
            "comparative analyses require exactly two groups with one explicit "
            "control role and one explicit case role"
        )
    if len(info["control_groups"]) != 1:
        errors.append(
            f"expected one control group in the sample sheet, found {info['control_groups'] or 'none'}"
        )
    if len(info["case_groups"]) != 1:
        errors.append(
            f"expected one case group in the sample sheet, found {info['case_groups'] or 'none'}"
        )
    expected = (
        f"{info['control_group']}_vs_{info['case_group']}"
        if info["control_group"] and info["case_group"]
        else None
    )
    if expected and cfg.get("comparison") != expected:
        errors.append(
            f"comparison must be {expected!r} from the sample-sheet roles, "
            f"not {cfg.get('comparison')!r}"
        )
    return errors


def _stage_dependencies(stage_id, cfg, selected_stages=()):
    """Return selected-stage prerequisites without changing analysis methods."""
    kind = _STAGE_META[stage_id]["kind"]
    if kind in {"diff", "mesa"} and cfg is not None:
        modalities = _config_params(
            cfg, "analysis", kind, "params", "modalities", default=["cpg"]
        ) or ["cpg"]
        dependencies = []
        for modality in modalities:
            producer = _FEATURE_MATRIX_STAGES.get(modality)
            if producer and producer not in dependencies:
                dependencies.append(producer)
        return tuple(dependencies)
    if stage_id == "analysis.report":
        return tuple(stage for stage in selected_stages if stage != stage_id)
    return ()


def _expand_and_order_stages(stages, cfg):
    """Add available feature producers and order selected stages by dependency."""
    expanded = list(stages)
    for stage_id in tuple(expanded):
        for dependency in _stage_dependencies(stage_id, cfg):
            if dependency not in expanded:
                expanded.append(dependency)

    ordered = []
    visiting = set()

    def visit(stage_id):
        if stage_id in ordered:
            return
        if stage_id in visiting:
            raise AnalysisContractError(f"cyclic analysis dependency at {stage_id}")
        visiting.add(stage_id)
        for dependency in _stage_dependencies(stage_id, cfg, expanded):
            if dependency in expanded:
                visit(dependency)
        visiting.remove(stage_id)
        ordered.append(stage_id)

    for stage_id in expanded:
        visit(stage_id)
    return tuple(ordered)


def resolve_stages(preset="auto", cfg=None, explicit=None):
    """Resolve a preset or explicit stage aliases into a stable stage list."""
    if explicit:
        stages = []
        for value in explicit:
            if value in STAGE_CHOICES:
                expanded = (value,)
            elif value in _STAGE_ALIASES:
                expanded = _STAGE_ALIASES[value]
            else:
                raise AnalysisContractError(
                    f"unknown analysis stage {value!r}; choose from {', '.join(STAGE_CHOICES)}"
                )
            for stage in expanded:
                if stage not in stages:
                    stages.append(stage)
        return _expand_and_order_stages(stages, cfg)

    preset = preset or "auto"
    if preset not in PRESET_NAMES:
        raise AnalysisContractError(
            f"unknown analysis preset {preset!r}; choose from {', '.join(PRESET_NAMES)}"
        )
    if preset == "auto":
        if cfg is not None and len(cfg.get("samples", {})) == 2:
            # The automatic path is deliberately bounded. DMR and MESA remain
            # explicit because they invoke external R or expensive modeling.
            stages = ("analysis.diff", "fragmentomics.occupancy", "fragmentomics.wps", "analysis.report")
        else:
            stages = _PRESETS["descriptive"]
        return _expand_and_order_stages(stages, cfg)
    return _expand_and_order_stages(_PRESETS[preset], cfg)


def _config_params(cfg, *keys, default=None):
    value = cfg
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def _load_context(args):
    config_path = Path(args.config).expanduser().resolve()
    if not config_path.is_file():
        raise AnalysisContractError(f"project config not found: {config_path}")
    try:
        raw = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AnalysisContractError(f"could not read project config {config_path}: {exc}") from exc
    if raw.get("schema_version") != 2:
        raise AnalysisContractError(
            "cftk analyze requires a schema-v2 project with explicit sample-sheet roles; "
            "use the expert downstream commands for legacy configurations"
        )
    lock_path = config_path.with_name("cftk.lock.json")
    if not lock_path.is_file():
        raise AnalysisContractError(
            f"project lock not found beside the config: {lock_path}; rerun 'cftk init'"
        )
    try:
        cfg = load_config(
            str(config_path),
            acquire_references=False,
            verify_profile=False,
            validate_profile=False,
        )
    except (SystemExit, OSError, KeyError, ValueError) as exc:
        raise AnalysisContractError(f"project configuration could not be resolved: {exc}") from exc
    paths = get_work_paths(cfg)
    scope_request = getattr(args, "fragmentomics_scope", None)
    return {
        "config_path": config_path,
        "lock_path": lock_path,
        "raw": raw,
        "cfg": cfg,
        "samples": get_all_samples(cfg),
        "paths": paths,
        "fragmentomics_scope_request": scope_request,
        "fragmentomics_scope": None,
        "identity": {
            "config_sha256": _sha256(config_path),
            "lock_sha256": _sha256(lock_path),
        },
    }


def _populate_scope_context(context, stages):
    """Resolve targeted scope only when a scoped fragmentomics stage is selected."""

    scoped_kinds = {"occupancy", "wps", "delfi"}
    if not any(_STAGE_META[stage]["kind"] in scoped_kinds for stage in stages):
        context["fragmentomics_scope"] = None
        return None
    context["fragmentomics_scope"] = describe_scope(
        context["cfg"],
        context["paths"],
        [sample["name"] for sample in context["samples"]],
        requested=context.get("fragmentomics_scope_request"),
        bam_paths=[get_bam(sample, context["paths"]) for sample in context["samples"]],
    )
    return context["fragmentomics_scope"]


def _sample_names(context):
    return [sample["name"] for sample in context["samples"]]


def _fragment_params(context, stage_kind):
    cfg = context["cfg"]
    frag = _config_params(cfg, "analysis", "frag", default={}) or {}
    if stage_kind == "end_motif":
        return {"kmer": _config_params(frag, "end_motif", "params", "kmer", default=4)}
    return {}


def _spec(
    path,
    description,
    *,
    role="output",
    required=True,
    nonempty=True,
    owned=True,
):
    spec = _core._spec(
        path,
        description,
        role=role,
        required=required,
        nonempty=nonempty,
    )
    spec["owned"] = owned
    return spec


def _artifact_specs(context, stage_id):
    paths = context["paths"]
    samples = _sample_names(context)
    kind = _STAGE_META[stage_id]["kind"]
    specs = []
    if kind == "diff":
        modalities = _config_params(
            context["cfg"], "analysis", "diff", "params", "modalities", default=["cpg"]
        ) or ["cpg"]
        for modality in modalities:
            base = Path(paths["differential"]) / modality
            specs.extend([
                _spec(base / "pca_coordinates.txt", f"PCA coordinates for {modality}", role="report"),
                _spec(base / "pca_variance.txt", f"PCA variance for {modality}", role="report"),
                _spec(base / "differential_result.tsv", f"differential result for {modality}", role="report"),
                _spec(base / "pca.png", f"PCA figure for {modality}", role="figure"),
                _spec(base / "pca.pdf", f"PCA PDF for {modality}", role="figure"),
                _spec(base / "violin.png", f"violin figure for {modality}", role="figure"),
                _spec(base / "violin.pdf", f"violin PDF for {modality}", role="figure"),
                _spec(base / "heatmap.png", f"heatmap figure for {modality}", role="figure"),
                _spec(base / "heatmap.pdf", f"heatmap PDF for {modality}", role="figure"),
            ])
    elif kind == "dmr":
        base = Path(paths["differential"]) / "dmr"
        specs.extend([
            _spec(base / "metilene_input.bedGraph", "metilene input table", role="report"),
            _spec(base / "dmr_raw.bed", "raw DMR calls", role="report"),
            _spec(base / "dmr_annotated.bed", "annotated DMR calls", role="report"),
            _spec(base / "dmr_volcano.png", "DMR volcano figure", role="figure"),
            _spec(base / "dmr_volcano.pdf", "DMR volcano PDF", role="figure"),
        ])
    elif kind == "occupancy":
        base = Path(paths["occ_out"])
        for sample in samples:
            specs.extend([
                _spec(base / f"{sample}.occupancy.tsv", f"occupancy table for {sample}", role="report"),
                _spec(base / f"{sample}.bw", f"occupancy bigWig for {sample}"),
                _spec(base / f"{sample}_occupancy.png", f"occupancy figure for {sample}", role="figure"),
                _spec(base / f"{sample}_occupancy.pdf", f"occupancy PDF for {sample}", role="figure"),
            ])
        if len(samples) > 1:
            specs.append(_spec(base / "occupancy_matrix.tsv", "occupancy feature matrix", role="report"))
    elif kind == "wps":
        base = Path(paths["wps_out"])
        for sample in samples:
            specs.extend([
                _spec(base / f"{sample}.wps.tsv", f"WPS table for {sample}", role="report"),
                _spec(base / f"{sample}.wps_profile.png", f"WPS figure for {sample}", role="figure"),
                _spec(base / f"{sample}.wps_profile.pdf", f"WPS PDF for {sample}", role="figure"),
            ])
        if len(samples) > 1:
            specs.append(_spec(base / "wps_matrix.tsv", "WPS feature matrix", role="report"))
    elif kind == "delfi":
        base = Path(paths["delfi_out"])
        for sample in samples:
            specs.extend([
                _spec(base / f"{sample}_delfi.tsv", f"DELFI table for {sample}", role="report"),
                _spec(base / f"{sample}_delfi_genome.png", f"DELFI figure for {sample}", role="figure"),
                _spec(base / f"{sample}_delfi_genome.pdf", f"DELFI PDF for {sample}", role="figure"),
            ])
    elif kind == "end_motif":
        base = Path(paths["end_motif_out"])
        kmer = _fragment_params(context, kind)["kmer"]
        for sample in samples:
            specs.extend([
                _spec(base / f"{sample}_{kmer}mer.tsv", f"{kmer}-mer table for {sample}", role="report"),
                _spec(base / f"{sample}_{kmer}mer_top20.png", f"end-motif figure for {sample}", role="figure"),
                _spec(base / f"{sample}_{kmer}mer_top20.pdf", f"end-motif PDF for {sample}", role="figure"),
            ])
    elif kind == "cleavage":
        base = Path(paths["cleavage_out"])
        for sample in samples:
            specs.append(_spec(base / f"{sample}_cleavage.bw", f"cleavage bigWig for {sample}", role="report"))
        # Cleavage plotting is optional when pyBigWig or a compatible site BED
        # is unavailable; the primary bigWig contract remains required.
        if context["cfg"].get("reference_data", {}).get("ctcf_bed"):
            specs.extend([
                _spec(base / "cleavage.png", "cleavage profile figure", role="figure", required=False),
                _spec(base / "cleavage.pdf", "cleavage profile PDF", role="figure", required=False),
            ])
    elif kind == "mesa":
        base = Path(paths["mesa"])
        specs.extend([
            _spec(base / "modality_performance.tsv", "MESA modality performance", role="report"),
            _spec(base / "MESA_model.pkl", "trained MESA model", role="output"),
            _spec(base / "loocv_predictions.tsv", "MESA LOOCV predictions", role="report"),
            _spec(base / "mesa_roc.png", "MESA ROC figure", role="figure"),
            _spec(base / "mesa_roc.pdf", "MESA ROC PDF", role="figure"),
            _spec(base / "mesa_heatmap.png", "MESA probability heatmap", role="figure"),
            _spec(base / "mesa_heatmap.pdf", "MESA probability heatmap PDF", role="figure"),
            _spec(base / "mesa_spearman.png", "MESA Spearman figure", role="figure"),
            _spec(base / "mesa_spearman.pdf", "MESA Spearman PDF", role="figure"),
        ])
    elif kind == "report":
        specs.append(_spec(Path(paths["report"]) / "report.html", "self-contained HTML report", role="report"))
    if kind in {"occupancy", "wps", "delfi"}:
        specs.append(
            _spec(
                scope_metadata_path(paths, kind),
                f"resolved fragmentomics scope metadata for {kind}",
                role="metadata",
            )
        )
        for path in scope_artifact_paths(
            context["cfg"],
            paths,
            samples,
            kind,
            requested=context.get("fragmentomics_scope_request"),
            bam_paths=[get_bam(sample, paths) for sample in context["samples"]],
        ):
            path = Path(path)
            role = "input" if path.suffix == ".bed" else "output"
            specs.append(
                _spec(
                    path,
                    f"targeted fragmentomics scope artifact for {kind}",
                    role=role,
                    owned=False,
                )
            )
    return specs


def _stage_requirements(context, stage_id):
    cfg = context["cfg"]
    paths = context["paths"]
    kind = _STAGE_META[stage_id]["kind"]
    requirements = {
        "python": set(),
        "executables": set(),
        "references": [],
        "inputs": [],
    }
    if kind == "diff":
        requirements["python"].update(("scipy", "sklearn", "matplotlib"))
        modalities = _config_params(cfg, "analysis", "diff", "params", "modalities", default=["cpg"]) or ["cpg"]
        requirements["inputs"].extend(get_matrix_path(paths, modality) for modality in modalities)
    elif kind == "dmr":
        requirements["executables"].update(("bedtools", "metilene", "Rscript"))
        requirements["references"].append(
            "R packages annotatr, TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicRanges, and org.Hs.eg.db"
        )
        for sample in context["samples"]:
            requirements["inputs"].append(Path(paths["methylation"]) / f"{sample['name']}_CpG.bedGraph")
    elif stage_id in _FRAGMENT_STAGES:
        requirements["inputs"].extend(get_bam(sample, paths) for sample in context["samples"])
        scope = {"mode": "genome"}
        if kind in {"occupancy", "wps", "delfi"}:
            try:
                scope = resolve_scope(
                    cfg, context.get("fragmentomics_scope_request")
                )
            except ScopeError as exc:
                requirements["references"].append(str(exc))
                scope = {"mode": "invalid"}
        if kind == "occupancy":
            requirements["executables"].update(("python", "wigToBigWig", "bigWigAverageOverBed"))
            requirements["references"].extend(("chrom_sizes", "tss_pas_bed"))
            tool = _config_params(cfg, "analysis", "frag", "occupancy", "tool", default="danpos")
            requirements["executables"].add(str(tool))
        elif kind == "wps":
            requirements["python"].update(("bx", "pysam"))
            requirements["references"].append("tss_pas_bed")
        elif kind == "delfi":
            requirements["python"].add("finaletoolkit")
            requirements["references"].extend(("chrom_sizes", "genome_2bit", "blacklist", "gap", "bins"))
        elif kind == "end_motif":
            requirements["python"].add("finaletoolkit")
            requirements["references"].append("genome_2bit")
        elif kind == "cleavage":
            requirements["python"].update(("finaletoolkit", "pyBigWig"))
            requirements["references"].extend(("chrom_sizes", "ctcf_bed"))
        if kind in {"occupancy", "wps", "delfi"} and scope.get("mode") == "panel":
            requirements["executables"].add("samtools")
            requirements["references"].append("target_bed")
    elif kind == "mesa":
        requirements["python"].update(("mesa", "sklearn", "scipy"))
        modalities = _config_params(cfg, "analysis", "mesa", "params", "modalities", default=["cpg"]) or ["cpg"]
        requirements["inputs"].extend(get_matrix_path(paths, modality) for modality in modalities)
    elif kind == "report":
        requirements["python"].update(("pandas", "numpy"))
        requirements["inputs"].append(paths["results"])
    return {
        "python": sorted(requirements["python"]),
        "executables": sorted(requirements["executables"]),
        "references": list(requirements["references"]),
        "inputs": [str(value) for value in requirements["inputs"]],
    }


def _resource_plan(context, stages, args):
    cfg = context["cfg"]
    samples = len(context["samples"])
    process = cfg.get("process", {})
    requested_parallel = getattr(args, "parallel", None) or process.get("parallel_samples", 1)
    scheduler = detect_scheduler_allocation()
    records = []
    for stage_id in stages:
        kind = _STAGE_META[stage_id]["kind"]
        if kind == "dmr":
            total = _config_params(cfg, "analysis", "dmr", "params", "cores", default=20)
        elif kind == "mesa":
            total = _config_params(cfg, "process", "step4_methylation", "params", "cores", default=20)
        elif kind == "report":
            total = 1
        else:
            total = _config_params(cfg, "process", "step3_markdup", "params", "cores", default=20)
        effective_parallel = requested_parallel if kind not in {"diff", "dmr", "mesa", "report"} else 1
        try:
            planned = plan_parallelism(total, effective_parallel, samples)
        except ValueError as exc:
            raise AnalysisContractError(f"invalid resource plan for {stage_id}: {exc}") from exc
        if kind in {"diff", "dmr", "mesa", "report"}:
            planned.update({
                "concurrent_samples": 1,
                "threads_per_sample": max(1, min(total, planned["threads_per_sample"])),
                "estimated_peak_threads": min(total, max(1, planned["threads_per_sample"])),
            })
        records.append({"stage": stage_id, "model": "stage-level downstream execution", **planned})
    maximum = max((row["total_core_budget"] for row in records), default=0)
    allocated = scheduler.get("allocated_cores")
    return {
        "resource_plan_version": 1,
        "sample_count": samples,
        "requested_parallel_samples": requested_parallel,
        "maximum_total_core_budget": maximum,
        "scheduler": scheduler,
        "scheduler_compatible": allocated is None or maximum <= allocated,
        "stages": records,
    }


def _stage_command(context, stage_id, args):
    command = ["cftk", "--config", str(context["config_path"])]
    kind = _STAGE_META[stage_id]["kind"]
    if kind == "diff":
        command.append("diff")
    elif kind == "dmr":
        command.append("dmr")
    elif kind == "mesa":
        command.extend(("mesa", "--performance", "--mesa-model", "--loocv"))
    elif kind == "report":
        command.append("report")
    else:
        command.extend(("frag", f"--{kind.replace('_', '-')}"))
    if (
        kind in {"occupancy", "wps", "delfi"}
        and getattr(args, "fragmentomics_scope", None) is not None
    ):
        command.extend(("--fragmentomics-scope", str(args.fragmentomics_scope)))
    if getattr(args, "parallel", None):
        command.extend(("--parallel", str(args.parallel)))
    return " ".join(__import__("shlex").quote(str(value)) for value in command)


def _validate_reference_fields(context, stage_id):
    cfg = context["cfg"]
    ref = cfg.get("reference_data", {})
    kind = _STAGE_META[stage_id]["kind"]
    fields = {
        "occupancy": ("chrom_sizes", "tss_pas_bed"),
        "wps": ("tss_pas_bed",),
        "delfi": ("chrom_sizes", "genome_2bit", "blacklist", "gap", "bins"),
        "end_motif": ("genome_2bit",),
        "cleavage": ("chrom_sizes", "ctcf_bed"),
    }.get(kind, ())
    if kind in {"occupancy", "wps", "delfi"}:
        try:
            scope = resolve_scope(
                cfg, context.get("fragmentomics_scope_request")
            )
        except ScopeError as exc:
            return [str(exc)]
        if scope["mode"] == "panel":
            fields = ("target_bed",) + tuple(fields)
    errors = []
    for field in fields:
        value = ref.get(field)
        if not value:
            errors.append(f"reference_data.{field} is not configured")
        elif not Path(value).is_file():
            errors.append(f"reference_data.{field} does not exist: {value}")
    return errors


def _required_inputs(context, stage_id, *, planned_outputs=()):
    paths = context["paths"]
    cfg = context["cfg"]
    kind = _STAGE_META[stage_id]["kind"]
    errors = []
    planned_paths = {Path(path).resolve() for path in planned_outputs}

    def missing_or_empty(path, description):
        path = Path(path)
        if path.resolve() in planned_paths:
            return
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"{description} is missing or empty: {path}")

    if stage_id in _FRAGMENT_STAGES:
        for sample in context["samples"]:
            bam = Path(get_bam(sample, paths))
            if not bam.is_file():
                errors.append(f"BAM for {sample['name']} is missing: {bam}")
            elif not Path(f"{bam}.bai").is_file() and not bam.with_suffix(".bai").is_file():
                errors.append(f"BAM index for {sample['name']} is missing: {bam}.bai")
    elif kind == "diff":
        modalities = _config_params(cfg, "analysis", "diff", "params", "modalities", default=["cpg"]) or ["cpg"]
        for modality in modalities:
            path = Path(get_matrix_path(paths, modality))
            missing_or_empty(path, f"matrix for modality {modality!r}")
    elif kind == "dmr":
        for sample in context["samples"]:
            path = Path(paths["methylation"]) / f"{sample['name']}_CpG.bedGraph"
            if not path.is_file() or path.stat().st_size == 0:
                errors.append(f"CpG bedGraph for {sample['name']} is missing or empty: {path}")
    elif kind == "mesa":
        modalities = _config_params(cfg, "analysis", "mesa", "params", "modalities", default=["cpg"]) or ["cpg"]
        for modality in modalities:
            path = Path(get_matrix_path(paths, modality))
            missing_or_empty(path, f"matrix for MESA modality {modality!r}")
    return errors


def _check_dmr_r_packages(checks):
    """Verify the packages imported by the bundled DMR annotation script."""
    rscript = shutil.which("Rscript")
    if not rscript:
        # The executable check already supplies the actionable failure.
        return
    package_values = ", ".join(json.dumps(package) for package in _DMR_R_PACKAGES)
    expression = (
        f"packages <- c({package_values}); "
        "missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; "
        "if (length(missing)) { message('Missing R packages: ', paste(missing, collapse = ', ')); quit(status = 1) }"
    )
    try:
        completed = subprocess.run(
            [rscript, "--vanilla", "-e", expression],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        checks.fail(
            "analysis.dmr.r_packages",
            f"Could not check DMR R annotation packages: {exc}",
            remedy="Activate the R environment containing CFTK's DMR annotation packages.",
        )
        return
    if completed.returncode == 0:
        checks.pass_(
            "analysis.dmr.r_packages",
            "Required DMR R annotation packages are available",
            details={"path": rscript, "packages": list(_DMR_R_PACKAGES)},
        )
        return
    output = (completed.stderr or completed.stdout or "package probe failed").strip()
    checks.fail(
        "analysis.dmr.r_packages",
        f"DMR R annotation packages are unavailable: {output[:300]}",
        remedy="Install annotatr, TxDb.Hsapiens.UCSC.hg38.knownGene, GenomicRanges, and org.Hs.eg.db in the active R environment.",
    )


def analysis_doctor_checks(
    checks,
    cfg,
    stages,
    *,
    parallel_override=None,
    fragmentomics_scope=None,
    adopt_existing=False,
):
    """Add read-only, workflow-specific checks to the core doctor object."""
    role_errors = []
    if any(stage in _COMPARATIVE_STAGES for stage in stages):
        role_errors = comparison_role_errors(cfg)
    if role_errors:
        checks.fail(
            "analysis.roles",
            "; ".join(role_errors),
            remedy="Edit the sample-sheet role column and rerun 'cftk init'.",
        )
    else:
        info = _role_info(cfg)
        checks.pass_(
            "analysis.roles",
            f"Explicit roles resolve {len(info['groups'])} group(s): "
            f"{info['control_group'] or 'no control/case comparison'}",
        )

    context = {
        "cfg": cfg,
        "paths": get_work_paths(cfg),
        "samples": get_all_samples(cfg),
        "fragmentomics_scope_request": fragmentomics_scope,
    }
    if any(stage in {"fragmentomics.occupancy", "fragmentomics.wps", "fragmentomics.delfi"} for stage in stages):
        try:
            scope = describe_scope(
                cfg,
                context["paths"],
                [sample["name"] for sample in context["samples"]],
                requested=fragmentomics_scope,
                bam_paths=[get_bam(sample, context["paths"]) for sample in context["samples"]],
            )
        except ScopeError as exc:
            checks.fail(
                "analysis.fragmentomics_scope",
                str(exc),
                remedy="Repair the target BED/profile or pass --fragmentomics-scope genome only for validated whole-genome inputs.",
            )
        else:
            checks.pass_(
                "analysis.fragmentomics_scope",
                scope["note"],
                details=scope,
            )
    planned_outputs = {
        spec["path"]
        for stage in stages
        for spec in _artifact_specs(context, stage)
        if spec.get("required", True)
    }
    python_modules = set()
    executables = set()
    reference_errors = []
    input_errors = []
    for stage in stages:
        req = _stage_requirements(context, stage)
        python_modules.update(req["python"])
        executables.update(req["executables"])
        reference_errors.extend(_validate_reference_fields(context, stage))
        input_errors.extend(
            _required_inputs(context, stage, planned_outputs=planned_outputs)
        )

    for module in sorted(python_modules):
        if importlib.util.find_spec(module) is None:
            checks.fail(
                f"analysis.python.{module}",
                f"Python module {module!r} is not importable",
                remedy=f"Install the matching optional dependency extra before running this preset.",
            )
        else:
            checks.pass_(f"analysis.python.{module}", f"Python module {module!r} is importable")
    for executable in sorted(executables):
        if shutil.which(executable) is None:
            checks.fail(
                f"analysis.tool.{executable}",
                f"Required executable {executable!r} was not found in PATH",
                remedy="Activate the CFTK workflow environment or install the tool listed in the installation guide.",
            )
        else:
            checks.pass_(f"analysis.tool.{executable}", f"Executable {executable!r} is available")
    if "analysis.dmr" in stages:
        _check_dmr_r_packages(checks)
    if reference_errors:
        for error in sorted(set(reference_errors)):
            checks.fail("analysis.reference", error, remedy="Install or repair the selected reference profile with 'cftk init'.")
    elif any(stage in _FRAGMENT_STAGES for stage in stages):
        checks.pass_("analysis.reference", "Fragmentomics reference components are present")
    if input_errors:
        for error in sorted(set(input_errors)):
            checks.fail("analysis.input", error, remedy="Complete core processing/QC or correct the project inputs before analysis.")
    elif stages:
        checks.pass_("analysis.input", "All selected analysis inputs are present and nonempty")

    try:
        resources = _resource_plan(
            context, stages, SimpleNamespace(parallel=parallel_override)
        )
        maximum = resources["maximum_total_core_budget"]
        scheduler = resources["scheduler"]
        allocated = scheduler.get("allocated_cores")
        complete_for_adoption = adopt_existing and all(
            not _core._validate_artifacts(_artifact_specs(context, stage))
            for stage in stages
        )
        if complete_for_adoption:
            checks.pass_(
                "analysis.resource",
                "All selected stages are complete; adoption-only validation does not require compute allocation",
                details=resources,
            )
        elif allocated is not None and maximum > allocated:
            checks.fail(
                "analysis.resource",
                f"analysis CPU budget {maximum} exceeds scheduler allocation {allocated}",
                remedy="Request enough CPUs or reduce process.cores in the project config.",
            )
        else:
            checks.pass_(
                "analysis.resource",
                "Analysis resource plan fits the configured CPU budget",
                details=resources,
            )
    except (AnalysisContractError, TypeError, ValueError) as exc:
        checks.fail("analysis.resource", f"Invalid analysis resource plan: {exc}")


def build_plan(context, args, *, doctor_report=None):
    stages = resolve_stages(getattr(args, "preset", "auto"), context["cfg"], getattr(args, "stages", None))
    role_errors = []
    if any(stage in _COMPARATIVE_STAGES for stage in stages):
        role_errors = comparison_role_errors(context["cfg"])
    resources = _resource_plan(context, stages, args)
    records = []
    for stage in stages:
        records.append({
            "id": stage,
            "name": _STAGE_META[stage]["name"],
            "kind": _STAGE_META[stage]["kind"],
            "comparative": _STAGE_META[stage]["comparative"],
            "depends_on": list(_stage_dependencies(stage, context["cfg"], stages)),
            "command": _stage_command(context, stage, args),
            "requirements": _stage_requirements(context, stage),
            "resources": next(row for row in resources["stages"] if row["stage"] == stage),
            "expected": _artifact_specs(context, stage),
            "status": "blocked" if role_errors and _STAGE_META[stage]["comparative"] else "planned",
            "fragmentomics_scope": (
                context.get("fragmentomics_scope")
                if _STAGE_META[stage]["kind"] in {"occupancy", "wps", "delfi"}
                else None
            ),
        })
    issues = [f"roles: {error}" for error in role_errors]
    if not resources["scheduler_compatible"]:
        issues.append("resource plan exceeds the detected scheduler allocation")
    if doctor_report:
        issues.extend(
            f"{check['id']}: {check['summary']}"
            for check in doctor_report.get("checks", [])
            if check.get("status") == "FAIL"
        )
    return {
        "plan_schema_version": ANALYSIS_RUN_SCHEMA_VERSION,
        "plan_id": _new_id("plan"),
        "workflow": "downstream-analysis",
        "preset": getattr(args, "preset", "auto") or "auto",
        "resolved_stages": list(stages),
        "config": str(context["config_path"]),
        "lock": str(context["lock_path"]),
        "project_identity": context["identity"],
        "roles": _role_info(context["cfg"]),
        "fragmentomics_scope": context.get("fragmentomics_scope"),
        "resource_plan": resources,
        "stages": records,
        "doctor": doctor_report,
        "status": "blocked" if issues else "ready",
        "issues": issues,
        "created_at": _utc_now(),
    }


def _save_plan(plan, context):
    root = Path(context["paths"]["provenance"]) / "analysis-plans" / plan["plan_id"]
    root.mkdir(parents=True, exist_ok=False)
    plan_path = root / "analysis-plan.json"
    _core._atomic_write_json(plan_path, plan)
    _core._atomic_write_json(Path(context["paths"]["provenance"]) / "latest-analysis-plan.json", {
        "plan_id": plan["plan_id"],
        "plan": str(plan_path.resolve()),
        "status": plan["status"],
    })
    return plan_path


def plan(args):
    try:
        context = _load_context(args)
        stages = resolve_stages(getattr(args, "preset", "auto"), context["cfg"], getattr(args, "stages", None))
        _populate_scope_context(context, stages)
        from doctor import run_doctor
        report_args = SimpleNamespace(
            config=str(context["config_path"]), step=[], target_bed=None,
            skip_picard_metrics=True, parallel=getattr(args, "parallel", None),
            analysis_stages=list(stages), analysis_only=True,
            fragmentomics_scope=getattr(args, "fragmentomics_scope", None),
        )
        doctor_report = run_doctor(report_args)
        plan_payload = build_plan(context, args, doctor_report=doctor_report)
        plan_path = _save_plan(plan_payload, context)
    except (AnalysisContractError, ScopeError, OSError, ValueError) as exc:
        raise SystemExit(f"[plan] ERROR: {exc}") from exc
    if getattr(args, "json", False):
        print(json.dumps(plan_payload, indent=2, sort_keys=True))
    else:
        print(f"CFTK analysis plan: {plan_payload['status'].upper()}")
        print(f"Preset: {plan_payload['preset']} -> {', '.join(plan_payload['resolved_stages'])}")
        print(f"Plan: {plan_path}")
        if plan_payload["issues"]:
            print("Issues:")
            for issue in plan_payload["issues"]:
                print(f"  - {issue}")
    if plan_payload["status"] != "ready":
        raise SystemExit(1)
    return plan_payload


def _load_previous(provenance, identity):
    provenance = Path(provenance).resolve()
    candidates = []
    latest = provenance / "latest-analysis.json"
    if latest.is_file():
        try:
            candidates.append(Path(json.loads(latest.read_text())["manifest"]))
        except (OSError, KeyError, json.JSONDecodeError):
            pass
    candidates.extend(sorted((provenance / "analysis-runs").glob("*/run.json"), reverse=True))
    seen = set()
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            manifest = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        if manifest.get("project_identity") == identity:
            return manifest
    return None


def _load_previous_stage(provenance, identity, stage_id):
    """Find a trusted completed stage across compatible workflow selections."""

    provenance = Path(provenance).resolve()
    current_identity = dict(identity)
    current_identity.pop("options_sha256", None)
    candidates = []
    latest = provenance / "latest-analysis.json"
    if latest.is_file():
        try:
            candidates.append(Path(json.loads(latest.read_text())["manifest"]))
        except (OSError, KeyError, json.JSONDecodeError):
            pass
    candidates.extend(sorted((provenance / "analysis-runs").glob("*/run.json"), reverse=True))
    seen = set()
    trusted = {"complete", "resumed", "adopted"}
    for candidate in candidates:
        candidate = candidate.expanduser().resolve()
        if candidate in seen:
            continue
        seen.add(candidate)
        try:
            manifest = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        manifest_identity = dict(manifest.get("project_identity", {}))
        manifest_identity.pop("options_sha256", None)
        if manifest_identity != current_identity:
            continue
        stage = next(
            (item for item in manifest.get("stages", []) if item.get("id") == stage_id),
            None,
        )
        if stage and stage.get("status") in trusted:
            return manifest
    return None


def _save_attempt(manifest, run_dir, provenance):
    manifest_path = Path(run_dir) / "run.json"
    _core._atomic_write_json(manifest_path, manifest)
    _core._atomic_write_json(Path(run_dir) / "resource-plan.json", manifest.get("resource_plan", {}))
    _core._write_summary_html(manifest, Path(run_dir) / "run-summary.html", manifest["project_root"])
    _core._atomic_write_json(Path(provenance) / "latest-analysis.json", {
        "run_id": manifest["run_id"],
        "manifest": str(manifest_path.resolve()),
        "summary": str((Path(run_dir) / "run-summary.html").resolve()),
        "status": manifest["status"],
    })


def _finalize_evidence(manifest, run_dir, provenance, events_path):
    evidence_dir = Path(run_dir) / "evidence"
    manifest["evidence"] = {
        "status": "running",
        "directory": str(evidence_dir),
        "summary": str(evidence_dir / "workflow_validation_summary.json"),
        "files": [],
    }
    _save_attempt(manifest, run_dir, provenance)
    try:
        summary = _core._generate_evidence(manifest, run_dir)
    except Exception as exc:
        manifest["evidence"] = {
            "status": "failed",
            "directory": str(evidence_dir),
            "summary": str(evidence_dir / "workflow_validation_summary.json"),
            "files": sorted(path.name for path in evidence_dir.glob("*") if path.is_file())
            if evidence_dir.is_dir() else [],
            "error": f"{type(exc).__name__}: {exc}",
        }
        if manifest.get("status") == "complete":
            manifest["status"] = "complete_with_reporting_error"
        _core._append_event(events_path, "evidence_failed", run_id=manifest["run_id"], error=manifest["evidence"]["error"])
        _save_attempt(manifest, run_dir, provenance)
        return False
    manifest["evidence"] = summary
    _core._append_event(events_path, "evidence_completed", run_id=manifest["run_id"], files=summary.get("files", []))
    _save_attempt(manifest, run_dir, provenance)
    return True


def _stage_args(context, stage_id, args):
    common = {
        "config": str(context["config_path"]),
        "parallel": getattr(args, "parallel", None),
        "fragmentomics_scope": context.get("fragmentomics_scope_request"),
        "target_bed": None,
        "skip_picard_metrics": False,
        "force": False,
    }
    kind = _STAGE_META[stage_id]["kind"]
    if kind == "diff":
        # The legacy differential command defaults Joblib to all visible CPUs.
        # Bound it to the same stage budget recorded in this workflow attempt.
        resource = _resource_plan(context, (stage_id,), args)["stages"][0]
        return SimpleNamespace(
            **common,
            modality=None,
            cores=resource["threads_per_sample"],
        )
    if kind == "dmr":
        return SimpleNamespace(**common)
    if kind == "mesa":
        return SimpleNamespace(**common, performance=True, mesa_model=True, loocv=True, perf_tsv=None)
    if kind == "report":
        return SimpleNamespace(**common)
    flags = {fragment_kind: False for fragment_kind in ("occupancy", "wps", "delfi", "end_motif", "cleavage")}
    flags[kind] = True
    return SimpleNamespace(**common, **flags)


def _execute_stage(context, stage_id, args):
    import cftk

    stage_args = _stage_args(context, stage_id, args)
    kind = _STAGE_META[stage_id]["kind"]
    if kind == "diff":
        cftk._cmd_diff(stage_args)
    elif kind == "dmr":
        cftk._cmd_dmr(stage_args)
    elif kind == "mesa":
        cftk._cmd_mesa(stage_args)
    elif kind == "report":
        cftk._cmd_report(stage_args)
    else:
        cftk._cmd_frag(stage_args)


def run(args):
    try:
        context = _load_context(args)
        stages = resolve_stages(getattr(args, "preset", "auto"), context["cfg"], getattr(args, "stages", None))
        _populate_scope_context(context, stages)
        resources = _resource_plan(context, stages, args)
    except (AnalysisContractError, ScopeError, OSError, ValueError) as exc:
        raise SystemExit(f"[analyze] ERROR: {exc}") from exc

    identity_options = {
        "preset": getattr(args, "preset", "auto") or "auto",
        "stages": list(stages),
        "parallel": getattr(args, "parallel", None),
        "fragmentomics_scope": getattr(args, "fragmentomics_scope", None),
    }
    options = {
        **identity_options,
        "adopt_existing": bool(getattr(args, "adopt_existing", False)),
    }
    identity = dict(context["identity"])
    identity["options_sha256"] = hashlib.sha256(
        json.dumps(identity_options, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    context["identity"] = identity
    provenance = Path(context["paths"]["provenance"]).resolve()
    previous = _load_previous(provenance, identity)
    run_id = _new_id("analysis")
    run_dir = provenance / "analysis-runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    events_path = run_dir / "events.jsonl"
    configure_command_log(
        provenance / "commands.jsonl", run_id=run_id, mirror_paths=[run_dir / "commands.jsonl"]
    )
    (run_dir / "commands.jsonl").write_text("")
    manifest = {
        "run_schema_version": ANALYSIS_RUN_SCHEMA_VERSION,
        "run_id": run_id,
        "run_dir": str(run_dir),
        "workflow": "downstream-analysis",
        "project_root": str(Path(context["paths"]["work_dir"]).resolve()),
        "config": str(context["config_path"]),
        "lock": str(context["lock_path"]),
        "project_identity": identity,
        "options": options,
        "roles": _role_info(context["cfg"]),
        "fragmentomics_scope": context.get("fragmentomics_scope"),
        "resource_plan": resources,
        "previous_run_id": previous.get("run_id") if previous else None,
        "status": "running",
        "started_at": _utc_now(),
        "finished_at": None,
        "stages": [],
    }
    _core._append_event(events_path, "run_started", run_id=run_id)
    try:
        for stage_id in stages:
            resource = next(row for row in resources["stages"] if row["stage"] == stage_id)
            manifest["stages"].append({
                "id": stage_id,
                "name": _STAGE_META[stage_id]["name"],
                "kind": _STAGE_META[stage_id]["kind"],
                "comparative": _STAGE_META[stage_id]["comparative"],
                "depends_on": list(_stage_dependencies(stage_id, context["cfg"], stages)),
                "applicable": True,
                "command": _stage_command(context, stage_id, args),
                "resources": resource,
                "fragmentomics_scope": (
                    context.get("fragmentomics_scope")
                    if _STAGE_META[stage_id]["kind"] in {"occupancy", "wps", "delfi"}
                    else None
                ),
                "status": "pending",
                "started_at": None,
                "finished_at": None,
                "expected": _artifact_specs(context, stage_id),
                "requirements": _stage_requirements(context, stage_id),
                "outputs": [],
                "figures": [],
                "quarantined": [],
                "error": None,
            })
        _core._write_artifact_tables(manifest, run_dir)
        _save_attempt(manifest, run_dir, provenance)
    except BaseException as exc:
        manifest["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        manifest["finished_at"] = _utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _core._append_event(events_path, "run_failed", run_id=run_id, error=manifest["error"])
        _finalize_evidence(manifest, run_dir, provenance, events_path)
        raise SystemExit(130 if isinstance(exc, KeyboardInterrupt) else 1) from exc

    from doctor import run_doctor
    doctor_args = SimpleNamespace(
        config=str(context["config_path"]), step=[], target_bed=None,
        skip_picard_metrics=True, parallel=getattr(args, "parallel", None),
        analysis_stages=list(stages), analysis_only=True,
        fragmentomics_scope=getattr(args, "fragmentomics_scope", None),
        adopt_existing=bool(getattr(args, "adopt_existing", False)),
    )
    preflight = run_doctor(doctor_args)
    _core._atomic_write_json(run_dir / "doctor-before.json", preflight)
    _core._atomic_write_json(run_dir / "analysis-plan.json", build_plan(context, args, doctor_report=preflight))
    if getattr(args, "dry_run", False):
        for stage in manifest["stages"]:
            stage["status"] = "planned"
        manifest["status"] = "planned" if preflight["exit_code"] == 0 else "blocked"
        manifest["finished_at"] = _utc_now()
        _core._append_event(events_path, "run_planned", run_id=run_id, preflight_status=preflight["status"])
        _save_attempt(manifest, run_dir, provenance)
        if not _finalize_evidence(manifest, run_dir, provenance, events_path):
            raise SystemExit(1)
        disp(f"[analyze] plan written: {run_dir / 'run-summary.html'}")
        if preflight["exit_code"] != 0:
            raise SystemExit(1)
        return manifest
    if preflight["exit_code"] != 0:
        manifest["status"] = "failed"
        manifest["finished_at"] = _utc_now()
        manifest["error"] = f"analysis preflight reported {preflight['summary']['fail']} required failure(s)"
        _core._append_event(events_path, "run_failed", run_id=run_id, error=manifest["error"])
        _save_attempt(manifest, run_dir, provenance)
        _finalize_evidence(manifest, run_dir, provenance, events_path)
        raise SystemExit(1)

    try:
        for stage_record in manifest["stages"]:
            stage_id = stage_record["id"]
            specs = stage_record["expected"]
            if _core._can_resume(previous, identity, stage_id, specs):
                stage_record["status"] = "resumed"
                stage_record["finished_at"] = _utc_now()
                _core._record_artifacts(stage_record, specs)
                _core._append_event(events_path, "stage_resumed", stage=stage_id)
                _save_attempt(manifest, run_dir, provenance)
                continue
            issues = _core._validate_artifacts(specs)
            existing = [
                Path(spec["path"])
                for spec in specs
                if spec.get("owned", True) and Path(spec["path"]).exists()
            ]
            compatible_stage = _load_previous_stage(provenance, identity, stage_id)
            if compatible_stage and not issues:
                stage_record["status"] = "resumed"
                stage_record["finished_at"] = _utc_now()
                stage_record["reused_from_run_id"] = compatible_stage.get("run_id")
                _core._record_artifacts(stage_record, specs)
                _core._append_event(
                    events_path,
                    "stage_reused",
                    stage=stage_id,
                    source_run_id=compatible_stage.get("run_id"),
                )
                _save_attempt(manifest, run_dir, provenance)
                continue
            if getattr(args, "adopt_existing", False) and not issues:
                stage_record["status"] = "adopted"
                stage_record["finished_at"] = _utc_now()
                _core._record_artifacts(stage_record, specs)
                _core._append_event(events_path, "stage_adopted", stage=stage_id)
                _save_attempt(manifest, run_dir, provenance)
                continue
            previous_stage = next((item for item in (previous or {}).get("stages", []) if item.get("id") == stage_id), None)
            retryable = bool(previous and previous.get("project_identity") == identity and previous_stage and previous_stage.get("status") in {"running", "failed", "interrupted", "complete", "resumed", "adopted"})
            if existing and (retryable or getattr(args, "adopt_existing", False)):
                stage_record["quarantined"] = _core._quarantine_paths(existing, context["paths"]["results"], provenance / "quarantine" / run_id)
                _core._append_event(events_path, "stage_outputs_quarantined", stage=stage_id, count=len(stage_record["quarantined"]))
            elif existing:
                raise AnalysisContractError(
                    f"{stage_id} has outputs without a trusted matching analysis manifest; "
                    "rerun with --adopt-existing to validate them or review them first"
                )
            stage_record["status"] = "running"
            stage_record["started_at"] = _utc_now()
            _core._append_event(events_path, "stage_started", stage=stage_id, command=stage_record["command"])
            _save_attempt(manifest, run_dir, provenance)
            try:
                _execute_stage(context, stage_id, args)
                issues = _core._validate_artifacts(specs)
                if issues:
                    raise AnalysisContractError(f"{stage_id} did not satisfy its artifact contract: {'; '.join(issues)}")
            except BaseException as exc:
                stage_record["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
                stage_record["finished_at"] = _utc_now()
                stage_record["error"] = f"{type(exc).__name__}: {exc}"
                _core._append_event(events_path, "stage_failed", stage=stage_id, error=stage_record["error"])
                raise
            stage_record["status"] = "complete"
            stage_record["finished_at"] = _utc_now()
            _core._record_artifacts(stage_record, specs)
            _core._append_event(events_path, "stage_completed", stage=stage_id)
            _save_attempt(manifest, run_dir, provenance)
    except BaseException as exc:
        manifest["status"] = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        manifest["finished_at"] = _utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        _core._append_event(events_path, "run_failed", run_id=run_id, error=manifest["error"])
        _save_attempt(manifest, run_dir, provenance)
        _finalize_evidence(manifest, run_dir, provenance, events_path)
        if isinstance(exc, AnalysisContractError):
            disp(f"[analyze] ERROR: {exc}")
        raise SystemExit(130 if isinstance(exc, KeyboardInterrupt) else 1) from exc

    manifest["status"] = "complete"
    manifest["finished_at"] = _utc_now()
    _core._append_event(events_path, "run_completed", run_id=run_id)
    if not _finalize_evidence(manifest, run_dir, provenance, events_path):
        raise SystemExit(2)
    plan_id = getattr(args, "job_plan_id", None)
    if plan_id:
        from job_plan import write_finalizer_marker
        for stage_id in stages:
            write_finalizer_marker(
                context["paths"], plan_id, stage_id,
                task_count=len(context["samples"]),
            )
    disp(f"[analyze] complete: {run_dir / 'run-summary.html'}")
    return manifest
