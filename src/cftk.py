#!/usr/bin/env python3
"""cftk — cfDNA multimodal epigenetic analysis toolkit."""

import argparse
import os
import sys

_OPTIONAL_IMPORT_EXTRAS = {
    "adjustText": "analysis",
    "mesa": "analysis",
    "sklearn": "analysis",
    "statsmodels": "analysis",
    "bx": "fragmentomics",
    "finaletoolkit": "fragmentomics",
    "pyBigWig": "fragmentomics",
}

_SRC = os.path.dirname(os.path.abspath(__file__))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


def _load(args):
    from init import load_config, get_work_paths
    from util import configure_command_log
    cfg   = load_config(args.config)
    if "comparison" not in cfg:
        cg = cfg.get("control_group", "")
        ag = cfg.get("case_group",    "")
        if cg and ag:
            cfg["comparison"] = f"{cg}_vs_{ag}"
        else:
            # Descriptive workflows (QC, fragmentomics, reporting, and merge)
            # are valid for a one-group project. Comparative handlers call
            # get_group_names themselves and retain their focused error.
            cfg["comparison"] = None
    paths = get_work_paths(cfg)
    configure_command_log(os.path.join(paths["provenance"], "commands.jsonl"))
    return cfg, paths


def _p(cfg, *keys, default=None):
    obj = cfg
    for k in keys:
        if not isinstance(obj, dict) or k not in obj:
            return default
        obj = obj[k]
    return obj


# ── Sub-command handlers ──────────────────────────────────────────────────────

def _cmd_init(args):
    from init import init
    init(args)


def _cmd_doctor(args):
    import json
    from doctor import render_human, run_doctor

    report = run_doctor(args)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(render_human(report))
    raise SystemExit(report["exit_code"])


def _cmd_run(args):
    from run_workflow import run
    return run(args)


def _cmd_plan(args):
    if getattr(args, "execution", "local") == "per-sample":
        return _write_job_plan(args)
    if getattr(args, "slurm", False):
        raise SystemExit("[plan] ERROR: --slurm requires --execution per-sample")
    from analysis_workflow import plan
    return plan(args)


def _cmd_analyze(args):
    import json
    from analysis_workflow import run

    manifest = run(args)
    if getattr(args, "json", False):
        print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def _write_job_plan(args):
    from job_plan import write_fragmentomics_job_plan

    plan = write_fragmentomics_job_plan(args)
    if getattr(args, "json", False):
        import json
        print(json.dumps(plan, indent=2, sort_keys=True))
    else:
        print(f"CFTK per-sample job plan: {plan['plan_path']}")
        print(f"Tasks: {plan['task_count']} sample task(s), {plan['finalizer_count']} finalizer(s)")
        if plan.get("slurm_submit_script"):
            print(f"Optional Slurm submission script: {plan['slurm_submit_script']}")
        else:
            print("Task scripts are scheduler-neutral; submit them with your site scheduler.")
    return plan


def _cmd_job_plan(args):
    print(
        "[job-plan] compatibility alias; use 'cftk plan --execution per-sample' for new workflows.",
        file=sys.stderr,
    )
    return _write_job_plan(args)


def _cmd_process(args):
    from process import process
    process(args, config_path=args.config)


def _cmd_qc(args):
    from init import get_all_samples, get_bam, get_matrix_path
    from analysis.qc import run_qc
    from resource_planning import (
        detect_scheduler_allocation,
        ensure_scheduler_capacity,
        plan_parallelism,
    )
    from visualization.visualization import plot_qc

    cfg, paths = _load(args)
    qc_p = _p(cfg, "analysis", "qc", "params", default={})

    all_samples = get_all_samples(cfg)

    # M2: pass all_samples and paths for step 0 (qc_parser)
    args.all_samples  = all_samples
    args.paths        = paths
    args.infile       = [get_bam(s, paths) for s in all_samples]
    args.output_dir   = paths["qc"]
    args.matrices_dir = paths["cpg_matrix"]
    args.ref_fa       = _p(cfg, "reference_data", "genome_fa", default="")
    args.fragment     = qc_p.get("fragment",  167)
    args.step_size    = qc_p.get("step_size", 2000)
    total_cores = _p(
        cfg, "process", "step4_methylation", "params", "cores", default=1
    )
    requested_parallel = getattr(args, "parallel", None) or \
                         _p(cfg, "process", "parallel_samples", default=1)
    try:
        ensure_scheduler_capacity(total_cores, detect_scheduler_allocation())
        resource_plan = plan_parallelism(
            total_cores, requested_parallel, len(all_samples)
        )
    except ValueError as exc:
        raise SystemExit(f"[qc] ERROR: invalid CPU resource plan: {exc}") from exc
    args.total_cores = total_cores
    args.parallel = resource_plan["concurrent_samples"]
    # QC labels describe available sample groups; unlike comparative analyses,
    # QC remains valid for a processing-only project with one declared group.
    args.group_labels = {
        group: [sample["name"] for sample in members]
        for group, members in cfg.get("samples", {}).items()
    }
    os.makedirs(paths["qc"], exist_ok=True)

    steps = args.step if isinstance(args.step, list) else [args.step]
    for step in steps:
        args.step = step
        args.cores = (
            resource_plan["threads_per_sample"] if step == 2 else total_cores
        )
        run_qc(args)
        if step > 0:   # step 0 has no visualization
            plot_qc(args)


def _cmd_power(args):
    from analysis.power_analysis import run_power
    from visualization.visualization import plot_power

    cfg, paths = _load(args)
    pw = _p(cfg, "analysis", "power", "params", default={})

    args.sample_size    = getattr(args, "sample_size", None) or pw.get("sample_size", 100)
    args.effect_size    = getattr(args, "effect_size", None) or pw.get("effect_size", 0.1)
    args.depth          = pw.get("depth", [10, 20, 50])
    args.ratio          = pw.get("ratio", 1.0)
    args.plot_threshold = pw.get("plot_threshold", 0.8)
    args.step_size      = pw.get("step_size", 10000)
    args.cpg_std        = _p(cfg, "reference_data", "cpg_std", default="")
    args.output_dir     = paths["power"]
    os.makedirs(paths["power"], exist_ok=True)

    run_power(args)
    plot_power(args)


def _cmd_diff(args):
    from init import get_group_names, get_matrix_path
    from analysis.differential import run_differential
    from analysis.pca_analysis import run_pca
    from visualization.visualization import plot_differential
    from util import disp

    cfg, paths = _load(args)
    diff_p = _p(cfg, "analysis", "diff", "params", default={})
    ga, gb = get_group_names(cfg)

    args.group_labels = {
        ga: [s["name"] for s in cfg["samples"].get(ga, [])],
        gb: [s["name"] for s in cfg["samples"].get(gb, [])],
    }
    args.colors     = diff_p.get("colors", None)
    args.top_n      = diff_p.get("top_n_heatmap", 500)
    args.output_dir = paths["differential"]

    modalities = (
        [args.modality] if getattr(args, "modality", None)
        else diff_p.get("modalities", ["cpg"])
    )

    for mod in modalities:
        matrix = get_matrix_path(paths, mod)
        if not os.path.exists(matrix):
            disp(f"WARNING: matrix not found for '{mod}': {matrix} — skipping.")
            continue
        mod_out = os.path.join(paths["differential"], mod)
        os.makedirs(mod_out, exist_ok=True)
        args.infile       = matrix
        args.modality     = mod
        args.feature_name = mod
        run_pca(args)
        run_differential(args)
        plot_differential(args)


def _cmd_dmr(args):
    from init import get_group_names
    from analysis.dmr import run_dmr
    from visualization.visualization import plot_dmr

    cfg, paths  = _load(args)
    dmr_p       = _p(cfg, "analysis", "dmr", "params", default={})
    ga, gb      = get_group_names(cfg)
    dmr_samples = _p(cfg, "analysis", "dmr", "samples", default={})

    args.group_a        = ga
    args.group_b        = gb
    args.q_thr          = dmr_p.get("q_thr", 0.05)
    args.top_n          = dmr_p.get("top_n", 20)
    args.threads        = dmr_p.get("cores", 20)
    args.dmr_extra_args = dmr_p.get("extra_args", "")
    args.metilene_tool  = _p(cfg, "analysis", "dmr", "tool", default="metilene")
    args.output_dir     = os.path.join(paths["differential"], "dmr")
    os.makedirs(args.output_dir, exist_ok=True)

    args.bedgraph_a = _resolve_bedgraphs(cfg, paths, ga, dmr_samples.get(ga))
    args.bedgraph_b = _resolve_bedgraphs(cfg, paths, gb, dmr_samples.get(gb))

    run_dmr(args)
    plot_dmr(args)


def _resolve_bedgraphs(cfg, paths, group_name, selected_names=None):
    all_samples = cfg["samples"].get(group_name, [])
    if selected_names:
        valid = {s["name"] for s in all_samples}
        bad   = [n for n in selected_names if n not in valid]
        if bad:
            sys.exit(
                f"[dmr] ERROR: sample(s) {bad} not found in group '{group_name}'. "
                f"Available: {sorted(valid)}"
            )
        use = [s for s in all_samples if s["name"] in selected_names]
    else:
        use = all_samples
    return [
        os.path.join(paths["methylation"], f"{s['name']}_CpG.bedGraph")
        for s in use
    ]


def _cmd_frag(args):
    from init import get_all_samples, get_bam
    from resource_planning import (
        detect_scheduler_allocation,
        ensure_scheduler_capacity,
        plan_parallelism,
    )
    from analysis.assay_scope import ScopeError, prepare_scope, write_scope_metadata
    from visualization.visualization import plot_fragmentomics

    if getattr(args, "finalize", False) and getattr(args, "no_finalize", False):
        raise SystemExit("[frag] ERROR: --finalize and --no-finalize cannot be used together")
    if getattr(args, "finalize", False) and getattr(args, "samples", None):
        raise SystemExit("[frag] ERROR: --finalize requires the complete cohort; omit --sample")

    cfg, paths = _load(args)
    ref        = cfg["reference_data"]
    frag_cfg   = _p(cfg, "analysis", "frag", default={})

    all_samples = get_all_samples(cfg)
    requested_samples = list(getattr(args, "samples", None) or [])
    known_samples = {sample["name"] for sample in all_samples}
    unknown_samples = [name for name in requested_samples if name not in known_samples]
    if unknown_samples:
        raise SystemExit(
            "[frag] ERROR: unknown --sample value(s): "
            f"{', '.join(unknown_samples)}. Available: {', '.join(sorted(known_samples))}"
        )
    selected_samples = [
        sample for sample in all_samples
        if not requested_samples or sample["name"] in set(requested_samples)
    ]
    if not selected_samples:
        raise SystemExit("[frag] ERROR: no samples selected")
    original_infile = [get_bam(s, paths) for s in all_samples]
    selected_infile = [get_bam(s, paths) for s in selected_samples]
    args.infile = list(selected_infile)
    total_cores = _p(
        cfg, "process", "step3_markdup", "params", "cores", default=20
    )
    requested_parallel = getattr(args, "parallel", None) or _p(
        cfg, "process", "parallel_samples", default=1
    )
    try:
        allocation = detect_scheduler_allocation()
        if getattr(args, "finalize", False):
            # Finalization validates/plots existing files and performs no
            # sample computation, so it needs no cohort CPU budget.
            total_cores = 1
            requested_parallel = 1
        # A per-sample scheduler task owns only its requested CPU allocation;
        # do not apply the cohort-wide budget to that isolated task.
        if requested_samples and len(selected_samples) == 1:
            allocated = allocation.get("allocated_cores")
            if allocated is not None:
                total_cores = min(total_cores, allocated)
        ensure_scheduler_capacity(total_cores, allocation)
        resource_plan = plan_parallelism(
            total_cores, requested_parallel, len(selected_samples)
        )
    except ValueError as exc:
        raise SystemExit(f"[frag] ERROR: invalid CPU resource plan: {exc}") from exc
    args.cores = resource_plan["threads_per_sample"]
    args.parallel = resource_plan["concurrent_samples"]

    args.chrom_sizes = ref.get("chrom_sizes", "")
    args.genome2bit  = ref.get("genome_2bit", "")
    args.blacklist   = ref.get("blacklist", "")
    args.gap         = ref.get("gap", "")
    args.bins        = ref.get("bins", "")
    args.region      = ref.get("tss_pas_bed", "")
    args.bed         = ref.get("ctcf_bed", "")

    def _pf(sub, key, default=None):
        return _p(frag_cfg, sub, "params", key, default=default)

    args.occ_out       = paths["occ_out"]
    args.wps_out       = paths["wps_out"]
    args.delfi_out     = paths["delfi_out"]
    args.end_motif_out = paths["end_motif_out"]
    args.cleavage_out  = paths["cleavage_out"]

    args.danpos       = _p(frag_cfg, "occupancy", "tool", default="danpos")
    args.danpos_extra = _pf("occupancy", "extra_args", default="--paired 1 -u 0 -c 1000000")

    args.wps_window = _pf("wps", "wps_window", default=120)
    args.wps_step   = _pf("wps", "wps_step",   default=10)
    args.min_frag   = _pf("end_motif", "min_frag", default=100)
    args.max_frag   = _pf("end_motif", "max_frag", default=220)

    args.delfi_mapq   = _pf("delfi", "mapq",       default=30)
    args.delfi_window = _pf("delfi", "window",     default=20)
    args.delfi_extra  = _pf("delfi", "extra_args", default="")

    args.kmer     = _pf("end_motif", "kmer",       default=4)
    args.mapq     = _pf("end_motif", "mapq",       default=30)
    args.em_extra = _pf("end_motif", "extra_args", default="")

    args.window     = _pf("cleavage", "window",     default=20)
    args.upstream   = _pf("cleavage", "upstream",   default=1500)
    args.downstream = _pf("cleavage", "downstream", default=1500)
    args.cl_mapq    = _pf("cleavage", "mapq",       default=30)
    args.cl_extra   = _pf("cleavage", "extra_args", default="")

    requested_scope = getattr(args, "fragmentomics_scope", None) or _p(
        frag_cfg, "scope", default="auto"
    )

    run_all = not any([
        getattr(args, "occupancy",  False),
        getattr(args, "wps",        False),
        getattr(args, "delfi",      False),
        getattr(args, "end_motif",  False),
        getattr(args, "cleavage",   False),
    ])

    args.group_labels = {
        group: [s["name"] for s in members]
        for group, members in cfg.get("samples", {}).items()
    }

    selected_kinds = {
        kind
        for kind, enabled in (
            ("occupancy", run_all or getattr(args, "occupancy", False)),
            ("wps", run_all or getattr(args, "wps", False)),
            ("delfi", run_all or getattr(args, "delfi", False)),
        )
        if enabled
    }
    try:
        scope = prepare_scope(
            cfg,
            paths,
            all_samples,
            original_infile,
            requested=requested_scope,
            cores=args.cores,
            kinds=selected_kinds,
            materialize_samples=[sample["name"] for sample in selected_samples],
        ) if selected_kinds else {
            "info": {"requested": requested_scope, "mode": "not_applied"},
            "bam_paths": selected_infile,
            "region_bed": args.region,
            "bins": args.bins,
        }
    except ScopeError as exc:
        raise SystemExit(f"[frag] ERROR: {exc}") from exc
    args.fragmentomics_scope = requested_scope
    args.scope_info = scope["info"]
    if selected_kinds:
        print(
            f"[frag] scope={scope['info'].get('mode', 'unknown')}: "
            f"{scope['info'].get('note', '')}",
            file=sys.stderr,
        )

    def _finalize_stage(kind):
        expected = {
            "occupancy": (
                paths["occ_out"], ".occupancy.tsv", "occupancy",
            ),
            "wps": (
                paths["wps_out"], ".wps.tsv", "wps",
            ),
            "delfi": (
                paths["delfi_out"], "_delfi.tsv", "delfi",
            ),
            "end_motif": (
                paths["end_motif_out"], f"_{args.kmer}mer.tsv", "end-motif",
            ),
            "cleavage": (
                paths["cleavage_out"], "_cleavage.bw", "cleavage",
            ),
        }[kind]
        out_dir, suffix, label = expected
        files = [
            os.path.join(out_dir, f"{sample['name']}{suffix}")
            for sample in all_samples
        ]
        missing = [path for path in files if not os.path.isfile(path) or os.path.getsize(path) == 0]
        if missing:
            raise SystemExit(
                f"[frag] ERROR: cannot finalize {label}; {len(missing)} per-sample output(s) are missing, "
                f"including {missing[0]}"
            )
        if kind == "occupancy" and len(files) > 1:
            from analysis.occupancy import _merge_occupancy
            _merge_occupancy(
                list(zip(files, [sample["name"] for sample in all_samples])), out_dir
            )
        elif kind == "wps" and len(files) > 1:
            from analysis.wps import _merge_wps
            _merge_wps(
                list(zip(files, [sample["name"] for sample in all_samples])), out_dir
            )

    def _run_stage(kind, function, output_dir, plot_mode):
        saved = (args.infile, args.region, args.bins)
        if kind in selected_kinds and scope["info"].get("mode") == "panel":
            args.infile = scope["bam_paths"]
            if kind in {"occupancy", "wps"}:
                args.region = scope["region_bed"]
            if kind == "delfi":
                args.bins = scope["bins"]
        try:
            os.makedirs(output_dir, exist_ok=True)
            if kind in selected_kinds and not getattr(args, "no_finalize", False):
                write_scope_metadata(paths, kind, scope["info"])
            if getattr(args, "finalize", False):
                _finalize_stage(kind)
                plot_fragmentomics(args, mode=plot_mode)
            else:
                function(args)
                if not getattr(args, "no_finalize", False):
                    plot_fragmentomics(args, mode=plot_mode)
        finally:
            args.infile, args.region, args.bins = saved

    if run_all or getattr(args, "occupancy", False):
        from analysis.occupancy import run_occupancy
        _run_stage("occupancy", run_occupancy, paths["occ_out"], "occupancy")
    if run_all or getattr(args, "wps", False):
        from analysis.wps import run_wps
        _run_stage("wps", run_wps, paths["wps_out"], "wps")
    if run_all or getattr(args, "delfi", False):
        from analysis.delfi import run_delfi
        _run_stage("delfi", run_delfi, paths["delfi_out"], "delfi")
    if run_all or getattr(args, "end_motif", False):
        from analysis.end_motif import run_end_motif
        _run_stage("end_motif", run_end_motif, paths["end_motif_out"], "end_motif")
    if run_all or getattr(args, "cleavage", False):
        from analysis.cleavage import run_cleavage
        _run_stage("cleavage", run_cleavage, paths["cleavage_out"], "cleavage")


def _cmd_mesa(args):
    from init import get_group_names, get_matrix_path
    from analysis.mesa import run_modality_performance, run_mesa_model, run_mesa_loocv
    from visualization.visualization import plot_mesa
    from util import disp

    disp("[mesa] starting _cmd_mesa")
    try:
        import sys as _sys
        if not _sys.stdin.isatty():
            _sys.stdin = open(os.devnull, "r")
    except Exception:
        pass
    cfg, paths = _load(args)
    disp("[mesa] config loaded")
    mesa_p = _p(cfg, "analysis", "mesa", "params", default={})
    ga, gb = get_group_names(cfg)
    disp(f"[mesa] groups: {ga} vs {gb}")

    args.output_dir = paths["mesa"]
    args.clf        = mesa_p.get("clf",          [1, 2, 3])
    args.size       = mesa_p.get("feature_size", 100)
    args.subset     = mesa_p.get("subset",       0.1)
    args.repeat     = mesa_p.get("repeat",       3)
    args.cores      = _p(cfg, "process", "step4_methylation", "params", "cores", default=-1)
    os.makedirs(paths["mesa"], exist_ok=True)

    args.modality = mesa_p.get("modalities", ["cpg"])
    args.infile   = [get_matrix_path(paths, m) for m in args.modality]
    args.label    = _make_label(cfg, paths)

    performance = None
    if getattr(args, "performance", False):
        performance = run_modality_performance(args)
    if getattr(args, "mesa_model", False):
        run_mesa_model(args, performance=performance)
    if getattr(args, "loocv", False):
        run_mesa_loocv(args, performance=performance)
        plot_mesa(args)


def _make_label(cfg, paths):
    import pandas as pd
    from init import get_group_names
    from util import disp
    ga, gb = get_group_names(cfg)

    group_roles = cfg.get("group_roles", {})
    if not group_roles:
        sys.exit(
            "[label] ERROR: MESA requires explicit control/case roles in the "
            "schema-v2 sample sheet; group-name heuristics are not supported."
        )
    if group_roles.get(ga) != "control" or group_roles.get(gb) != "case":
        sys.exit(
            "[label] ERROR: comparison order must be control_vs_case when "
            "explicit group_roles are provided."
        )
    label_a, label_b = 0, 1

    rows = [(s["name"], label_a) for s in cfg["samples"].get(ga, [])] + \
           [(s["name"], label_b) for s in cfg["samples"].get(gb, [])]
    label_path = os.path.join(paths["mesa"], "label.tsv")
    os.makedirs(paths["mesa"], exist_ok=True)
    pd.DataFrame(rows).to_csv(label_path, sep="\t", header=False, index=False)
    disp(f"[label] {ga}={label_a}, {gb}={label_b} → {label_path}")
    return label_path


def _cmd_report(args):
    from report.report_generator import generate_report

    cfg, paths = _load(args)
    os.makedirs(paths["report"], exist_ok=True)

    args.results_dir  = paths["results"]
    args.output_dir   = paths["report"]
    args.project_name = cfg.get("project_name", "cftk_project")
    args.groups       = list(cfg.get("samples", {}).keys())
    args.resolved_config = cfg
    if not hasattr(args, "config"):
        args.config = "./cftk_init.json"

    generate_report(args)


def _cmd_merge(args):
    from analysis.merge import run_merge
    from util import disp

    cfg, paths = _load(args)
    modalities = getattr(args, "modality", None) or list(cfg.get("merge", {}).keys())
    if not modalities:
        disp("[merge] ERROR: specify --modality or add a 'merge' block in config.")
        sys.exit(1)
    for mod in modalities:
        run_merge(mod, cfg, paths)


def _cmd_vis(args):
    from init import get_all_samples, get_matrix_path
    from visualization.visualization import (
        plot_qc, plot_differential, plot_dmr,
        plot_fragmentomics, plot_mesa, plot_power,
    )
    from util import disp

    cfg, paths = _load(args)
    comparison = cfg.get("comparison")
    if isinstance(comparison, str) and "_vs_" in comparison:
        ga, gb = comparison.split("_vs_", 1)
    else:
        ga = gb = None
    diff_p     = _p(cfg, "analysis", "diff",  "params", default={})
    dmr_p      = _p(cfg, "analysis", "dmr",   "params", default={})
    frag_cfg   = _p(cfg, "analysis", "frag",  default={})

    modes = args.mode if args.mode else ["all"]
    if "all" in modes:
        modes = ["power", "qc", "diff", "dmr", "frag", "mesa"]

    group_labels = {
        group: [s["name"] for s in members]
        for group, members in cfg.get("samples", {}).items()
    }

    def _pf(sub, key, default=None):
        return _p(frag_cfg, sub, "params", key, default=default)

    if "power" in modes:
        disp("[vis] power")
        args.output_dir = paths["power"]
        plot_power(args)

    if "qc" in modes:
        qc_p = _p(cfg, "analysis", "qc", "params", default={})
        args.output_dir   = paths["qc"]
        args.matrices_dir = paths["cpg_matrix"]
        args.group_labels = group_labels
        qc_step = getattr(args, "step", None)
        if isinstance(qc_step, int):
            qc_steps = [qc_step]
        elif isinstance(qc_step, list) and len(qc_step) == 1:
            qc_steps = [qc_step[0]]
        else:
            qc_steps = [1, 2, 3]   # vis: skip step 0 (no plot output)
        for step in qc_steps:
            disp(f"[vis] qc step {step}")
            args.step = step
            plot_qc(args)

    if "diff" in modes:
        if not ga or not gb:
            sys.exit("[vis] ERROR: differential visualization requires explicit control/case groups.")
        modalities = diff_p.get("modalities", ["cpg"])
        for mod in modalities:
            matrix = get_matrix_path(paths, mod)
            if not os.path.exists(matrix):
                disp(f"[vis] diff: matrix not found for '{mod}', skipping.")
                continue
            disp(f"[vis] diff — {mod}")
            args.output_dir   = paths["differential"]
            args.infile       = matrix
            args.modality     = mod
            args.feature_name = mod
            args.group_labels = group_labels
            args.colors       = diff_p.get("colors", None)
            args.top_n        = diff_p.get("top_n_heatmap", 500)
            plot_differential(args)

    if "dmr" in modes:
        if not ga or not gb:
            sys.exit("[vis] ERROR: DMR visualization requires explicit control/case groups.")
        dmr_out = os.path.join(paths["differential"], "dmr")
        ann_bed = os.path.join(dmr_out, "dmr_annotated.bed")
        if not os.path.exists(ann_bed):
            disp("[vis] dmr: dmr_annotated.bed not found, skipping.")
        else:
            disp("[vis] dmr")
            args.output_dir = dmr_out
            args.group_a    = ga
            args.group_b    = gb
            args.q_thr      = dmr_p.get("q_thr", 0.05)
            args.top_n      = dmr_p.get("top_n", 20)
            plot_dmr(args)

    if "frag" in modes:
        ref = cfg["reference_data"]
        args.occ_out       = paths["occ_out"]
        args.wps_out       = paths["wps_out"]
        args.delfi_out     = paths["delfi_out"]
        args.end_motif_out = paths["end_motif_out"]
        args.cleavage_out  = paths["cleavage_out"]
        args.region        = ref.get("tss_pas_bed", "")
        args.bed           = ref.get("ctcf_bed", "")
        args.upstream      = _pf("cleavage", "upstream",   default=1500)
        args.downstream    = _pf("cleavage", "downstream", default=1500)
        args.group_labels  = group_labels
        for mode in ["occupancy", "wps", "delfi", "end_motif", "cleavage"]:
            disp(f"[vis] frag — {mode}")
            plot_fragmentomics(args, mode=mode)

    if "mesa" in modes:
        if not ga or not gb:
            sys.exit("[vis] ERROR: MESA visualization requires explicit control/case groups.")
        pred_tsv = os.path.join(paths["mesa"], "loocv_predictions.tsv")
        if not os.path.exists(pred_tsv):
            disp("[vis] mesa: loocv_predictions.tsv not found, skipping.")
        else:
            disp("[vis] mesa")
            args.output_dir = paths["mesa"]
            plot_mesa(args)


def _cmd_run_all(args):
    from util import disp

    args.step        = [1, 2, 3, 4]
    args.performance = True
    args.mesa_model  = True
    args.loocv       = True
    args.perf_tsv    = getattr(args, "perf_tsv", None)

    def _vis(mode):
        _sa(args, "mode", [mode])
        _sa(args, "step", None)
        _cmd_vis(args)

    import glob as _glob
    cfg, paths = _load(args)
    all_samples = [s["name"] for g in cfg["samples"].values() for s in g]

    def _done_process():
        # Re-enter process checkpoints so new or target-specific Picard metrics
        # are backfilled even when BAM and CpG outputs predate this integration.
        if not getattr(args, "skip_picard_metrics", False):
            return False
        cpg = os.path.join(paths["cpg_matrix"], "cpg_matrix.tsv")
        if not os.path.exists(cpg):
            return False
        return all(
            os.path.exists(os.path.join(paths["markdup"], f"{n}.markdup.bam"))
            for n in all_samples
        )

    def _done_qc(step):
        if step == 0:
            summary = os.path.join(paths["qc"], "qc_summary.tsv")
            scores  = os.path.join(paths["qc"], "qc_scores.tsv")
            return os.path.exists(summary) and os.path.exists(scores)
        step_dirs = {
            1: "1_methylation_distribution",
            2: "2_fragment_length",
            3: "3_dinucleotide_freq",
        }
        d = os.path.join(paths["qc"], step_dirs[step])
        return os.path.isdir(d) and bool(_glob.glob(os.path.join(d, "*.png")))

    def _done_diff():
        d = paths["differential"]
        return os.path.isdir(d) and bool(
            _glob.glob(os.path.join(d, "**", "*.tsv"), recursive=True))

    def _done_dmr():
        d = os.path.join(paths["differential"], "dmr")
        return os.path.isdir(d) and bool(_glob.glob(os.path.join(d, "*.bed*")))

    def _done_frag():
        d = paths["fragmentomics"]
        return os.path.isdir(d) and bool(
            _glob.glob(os.path.join(d, "**", "*.tsv"), recursive=True))

    def _done_mesa():
        return os.path.exists(os.path.join(paths["mesa"], "loocv_predictions.tsv"))

    pipeline = [
        ("process [1-4]",
         _cmd_process,
         None,
         _done_process),

        ("qc [methylation]",
         lambda a: (_sa(a, "step", [1]), _cmd_qc(a)),
         lambda a: (_sa(a, "step", [1]), _cmd_vis(a)),
         lambda: _done_qc(1)),

        # step 2 must run BEFORE step 0 so fragment CSVs exist for median_frag_len
        ("qc [fragment]",
         lambda a: (_sa(a, "step", [2]), _cmd_qc(a)),
         lambda a: (_sa(a, "step", [2]), _cmd_vis(a)),
         lambda: _done_qc(2)),

        ("qc [dinucleotide]",
         lambda a: (_sa(a, "step", [3]), _cmd_qc(a)),
         lambda a: (_sa(a, "step", [3]), _cmd_vis(a)),
         lambda: _done_qc(3)),

        # M2: step 0 — runs AFTER step 2 so fragment CSVs are present for median_frag_len
        ("qc [summary]",
         lambda a: (_sa(a, "step", [0]), _cmd_qc(a)),
         None,
         lambda: _done_qc(0)),

        ("diff",
         _cmd_diff,
         lambda a: _vis("diff"),
         _done_diff),

        ("dmr",
         _cmd_dmr,
         lambda a: _vis("dmr"),
         _done_dmr),

        ("frag",
         _cmd_frag,
         lambda a: _vis("frag"),
         _done_frag),

        ("mesa",
         _cmd_mesa,
         lambda a: _vis("mesa"),
         _done_mesa),

        ("report",
         _cmd_report,
         None,
         None),
    ]

    for step_name, ana_fn, vis_fn, check_fn in pipeline:
        if check_fn and check_fn():
            disp(f"[run-all] ── {step_name} — already done, skipping ──")
            continue

        disp(f"[run-all] ── {step_name} ──")
        ana_ok = True
        try:
            ana_fn(args)
        except Exception as e:
            import traceback
            disp(f"[run-all] WARNING: {step_name} failed: {e}")
            disp("[run-all] full traceback:")
            traceback.print_exc()
            disp("[run-all] continuing...")
            ana_ok = False

        if ana_ok and vis_fn:
            disp(f"[run-all] ── vis [{step_name}] ──")
            try:
                vis_fn(args)
            except Exception as e:
                import traceback
                disp(f"[run-all] WARNING: vis [{step_name}] failed: {e}")
                traceback.print_exc()

    disp("[run-all] pipeline complete.")


def _sa(obj, k, v):
    setattr(obj, k, v)


# ── Argument parser ───────────────────────────────────────────────────────────

def build_parser():
    parser = argparse.ArgumentParser(
        prog="cftk",
        description="cfDNA multimodal epigenetic analysis toolkit",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--config", default="./cftk_init.json", metavar="PATH",
        help="Path to cftk_init.json  (default: ./cftk_init.json)"
    )
    sub = parser.add_subparsers(dest="mode", metavar="<command>")

    # init
    p = sub.add_parser("init",
        help="Create or validate a project and prepare its reference genome.")
    p.add_argument(
        "--non-interactive", action="store_true",
        help="Create a missing config without prompts; requires explicit inputs.",
    )
    p.add_argument("--sample-sheet", default=None, metavar="PATH")
    p.add_argument("--reference-root", default=None, metavar="PATH")
    p.add_argument(
        "--reference-mode", choices=("local", "managed"), default=None,
        help="Reference mode for a new project (default: managed).",
    )
    p.add_argument("--profile", default=None, metavar="ID")
    p.add_argument("--profile-version", default=None, metavar="VERSION")
    p.add_argument("--project-name", default=None, metavar="NAME")
    p.add_argument("--output-dir", default=None, metavar="PATH")
    p.add_argument("--assay", default="twist_human_methylome")
    p.add_argument("--genome", default="hg38")
    p.add_argument(
        "--skip-reference-prep", action="store_true",
        help="Validate the config without building bwa-meth, FASTA, or Picard indexes.",
    )
    p.add_argument("--ref-index", dest="ref_index", action="store_true",
                   help=argparse.SUPPRESS)
    p.add_argument("--ref-dict",  dest="ref_dict",  action="store_true",
                   help=argparse.SUPPRESS)
    p.set_defaults(func=_cmd_init)

    # doctor
    p = sub.add_parser(
        "doctor",
        help="Check process readiness without downloading, repairing, or running data.",
    )
    p.add_argument(
        "-s", "--step", dest="step", type=int, nargs="+",
        choices=range(1, 5), default=[1, 2, 3, 4],
        metavar="{1,2,3,4}",
        help="Process steps to check (default: 1 2 3 4).",
    )
    p.add_argument(
        "--target-bed", default=None, metavar="PATH",
        help="Covered-target BED override to validate for Picard metrics.",
    )
    p.add_argument(
        "--skip-picard-metrics", action="store_true",
        help="Do not require Picard target/alignment metrics readiness.",
    )
    p.add_argument(
        "--json", action="store_true",
        help="Write a machine-readable JSON report to stdout.",
    )
    p.add_argument(
        "--parallel", type=int, default=None, metavar="N",
        help="Validate a parallel-sample override against the total CPU budget.",
    )
    p.add_argument(
        "--analysis-preset", choices=("auto", "descriptive", "differential", "dmr",
                                       "fragmentomics", "mesa", "comparative", "all", "report"),
        default=None,
        help="Also check a downstream analysis preset (read-only).",
    )
    p.add_argument(
        "--analysis-stage", dest="analysis_stages", nargs="+", default=None,
        help="Also check explicit downstream stage IDs or aliases.",
    )
    p.add_argument(
        "--fragmentomics-scope", choices=("auto", "panel", "genome"), default=None,
        help="Scope WPS/occupancy/DELFI reads and regions (default: assay-aware auto).",
    )
    p.set_defaults(func=_cmd_doctor)

    # beginner run
    p = sub.add_parser(
        "run",
        help="Run the validated schema-v2 core workflow with fail-fast resume.",
    )
    p.add_argument("--parallel", type=int, default=None, metavar="N")
    p.add_argument("--target-bed", default=None, metavar="PATH")
    p.add_argument(
        "--dry-run", action="store_true",
        help="Write the stage/output plan without probing tools or processing data.",
    )
    p.add_argument(
        "--adopt-existing", action="store_true",
        help="Explicitly validate and adopt complete pre-manifest outputs; quarantine partial ones.",
    )
    p.add_argument(
        "--qc-dinucleotide", action="store_true",
        help="Also run the expensive dinucleotide-frequency QC stage.",
    )
    p.set_defaults(func=_cmd_run)

    # downstream planning and analysis
    p = sub.add_parser(
        "plan",
        help="Plan downstream analyses, dependencies, resources, and outputs without running them.",
    )
    p.add_argument(
        "--preset",
        choices=("auto", "descriptive", "differential", "dmr", "fragmentomics",
                 "mesa", "comparative", "all", "report"),
        default="auto",
    )
    p.add_argument("--stage", dest="stages", nargs="+", default=None, metavar="STAGE")
    p.add_argument("--parallel", type=int, default=None, metavar="N")
    p.add_argument(
        "--execution", choices=("local", "per-sample"), default="local",
        help="Execution plan: local is read-only; per-sample writes task and finalizer scripts.",
    )
    p.add_argument(
        "--fragmentomics-scope", choices=("auto", "panel", "genome"), default=None,
        help="Scope WPS/occupancy/DELFI (default: panel for Twist, genome otherwise).",
    )
    p.add_argument(
        "--slurm", action="store_true",
        help="With --execution per-sample, also write an optional Slurm-array helper without submitting it.",
    )
    p.add_argument("--json", action="store_true", help="Write only the machine-readable plan to stdout.")
    p.set_defaults(func=_cmd_plan)

    p = sub.add_parser(
        "analyze",
        help="Run role-aware downstream analyses with preflight, checkpoints, provenance, and evidence.",
    )
    p.add_argument(
        "--preset",
        choices=("auto", "descriptive", "differential", "dmr", "fragmentomics",
                 "mesa", "comparative", "all", "report"),
        default="auto",
    )
    p.add_argument("--stage", dest="stages", nargs="+", default=None, metavar="STAGE")
    p.add_argument("--parallel", type=int, default=None, metavar="N")
    p.add_argument(
        "--fragmentomics-scope", choices=("auto", "panel", "genome"), default=None,
        help="Scope WPS/occupancy/DELFI (default: panel for Twist, genome otherwise).",
    )
    p.add_argument("--dry-run", action="store_true", help="Write the downstream plan and evidence without executing stages.")
    p.add_argument("--adopt-existing", action="store_true", help="Validate and adopt complete outputs from expert commands.")
    p.add_argument("--json", action="store_true", help="Write the final manifest as JSON after completion.")
    p.set_defaults(func=_cmd_analyze)

    p = sub.add_parser(
        "job-plan",
        help="Compatibility alias for 'plan --execution per-sample'.",
    )
    p.add_argument(
        "--stage", dest="stages", nargs="+", required=True,
        choices=("occupancy", "wps", "delfi", "end_motif", "cleavage"),
        metavar="STAGE",
    )
    p.add_argument(
        "--fragmentomics-scope", choices=("auto", "panel", "genome"), default=None,
        help="Scope WPS, occupancy, and DELFI as for cftk frag.",
    )
    p.add_argument(
        "--slurm", action="store_true",
        help="Also render an optional Slurm-array submission helper; it is never submitted automatically.",
    )
    p.add_argument("--json", action="store_true", help="Write the job plan JSON to stdout.")
    p.set_defaults(func=_cmd_job_plan, execution="per-sample")

    # process
    p = sub.add_parser("process",
        help="Part 1: Raw data processing (steps 1-4).\n"
             "  1 = adapter trimming\n"
             "  2 = bisulfite alignment\n"
             "  3 = mark duplicates\n"
             "  4 = CpG methylation calling + auto cpg_matrix merge")
    p.add_argument("-s", "--step", dest="step", type=int, nargs="+",
                   required=True, choices=range(1, 5), metavar="{1,2,3,4}")
    p.add_argument("--parallel", type=int, default=None, metavar="N")
    p.add_argument(
        "--target-bed", default=None, metavar="PATH",
        help="Covered-target BED override for Picard metrics.",
    )
    p.add_argument(
        "--skip-picard-metrics", action="store_true",
        help="Skip CollectHsMetrics and CollectMultipleMetrics after markdup.",
    )
    p.set_defaults(func=_cmd_process)

    # qc — M2: step range extended to 0-3
    p = sub.add_parser("qc",
        help="Part 2: QC analysis.\n"
             "  0 = parse process outputs → qc_summary.tsv + qc_scores.tsv\n"
             "  1 = methylation distribution (needs cpg_matrix)\n"
             "  2 = fragment length distribution\n"
             "  3 = dinucleotide frequency")
    p.add_argument("-s", "--step", dest="step", type=int, nargs="+",
                   required=True, choices=range(0, 4), metavar="{0,1,2,3}",
                   help="One or more QC steps. e.g. -s 0 1 2 3")
    p.add_argument("--title",    default=None)
    p.add_argument("--parallel", type=int, default=None, metavar="N")
    p.add_argument("--force",    action="store_true",
                   help="Re-run even if output files already exist")
    p.set_defaults(func=_cmd_qc)

    # power
    p = sub.add_parser("power", help="Part 2: Statistical power analysis.")
    p.add_argument("-s", "--sample-size", dest="sample_size", type=int, default=None)
    p.add_argument("-e", "--effect-size", dest="effect_size", type=float, default=None)
    p.set_defaults(func=_cmd_power)

    # diff
    p = sub.add_parser("diff",
        help="Part 2: Differential analysis — PCA / violin / heatmap.")
    p.add_argument("--modality", default=None)
    p.set_defaults(func=_cmd_diff)

    # dmr
    p = sub.add_parser("dmr",
        help="Part 2: DMR analysis — prepare + metilene + annotation + volcano.")
    p.set_defaults(func=_cmd_dmr)

    # frag
    p = sub.add_parser("frag",
        help="Part 2: Fragmentomics (all five if no flag given).")
    p.add_argument("--occupancy", action="store_true")
    p.add_argument("--wps",       action="store_true")
    p.add_argument("--delfi",     action="store_true")
    p.add_argument("--end-motif", dest="end_motif", action="store_true")
    p.add_argument("--cleavage",  action="store_true")
    p.add_argument(
        "--sample", dest="samples", action="append", default=None, metavar="NAME",
        help="Run only one named sample; repeat only for portable local batches.",
    )
    p.add_argument(
        "--parallel", type=int, default=None, metavar="N",
        help="Run up to N selected samples concurrently (per-sample plans use 1).",
    )
    p.add_argument(
        "--no-finalize", action="store_true",
        help="Write only selected per-sample artifacts; do not merge or plot cohort outputs.",
    )
    p.add_argument(
        "--finalize", action="store_true",
        help="Validate selected per-sample artifacts, then create cohort matrices and figures without recomputing them.",
    )
    p.add_argument(
        "--fragmentomics-scope", choices=("auto", "panel", "genome"), default=None,
        help="Scope WPS/occupancy/DELFI (default: panel for Twist, genome otherwise).",
    )
    p.set_defaults(func=_cmd_frag)

    # mesa
    p = sub.add_parser("mesa", help="Part 2: MESA multimodal modeling + LOOCV.")
    p.add_argument("--modality",   nargs="+", default=None)
    p.add_argument("--infile",     nargs="+", default=None)
    p.add_argument("--label",      default=None)
    p.add_argument("--perf-tsv",   dest="perf_tsv", default=None)
    p.add_argument("-p", "--performance", dest="performance", action="store_true")
    p.add_argument("--mesa-model", dest="mesa_model", action="store_true")
    p.add_argument("--loocv",      dest="loocv",      action="store_true")
    p.set_defaults(func=_cmd_mesa)

    # merge
    p = sub.add_parser("merge",
        help="Build feature matrix from user-specified files.")
    p.add_argument("--modality", nargs="+", default=None)
    p.set_defaults(func=_cmd_merge)

    # vis
    p = sub.add_parser("vis",
        help="Re-generate visualizations from existing results.")
    p.add_argument("--mode", nargs="+", default=None,
                   choices=["power", "qc", "diff", "dmr", "frag", "mesa", "all"],
                   metavar="MODE")
    p.set_defaults(func=_cmd_vis)

    # report
    p = sub.add_parser("report",
        help="Generate self-contained HTML report.")
    p.set_defaults(func=_cmd_report)

    # run-all
    p = sub.add_parser("run-all",
        help="Run full pipeline end-to-end.")
    p.add_argument("--parallel",   type=int, default=None)
    p.add_argument("--target-bed", default=None, metavar="PATH")
    p.add_argument("--skip-picard-metrics", action="store_true")
    p.add_argument("-p", "--performance", dest="performance", action="store_true")
    p.add_argument("--mesa-model", dest="mesa_model", action="store_true")
    p.add_argument("--loocv",      dest="loocv",      action="store_true")
    p.set_defaults(func=_cmd_run_all)

    return parser


def main():
    parser = build_parser()
    args   = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(0)
    try:
        args.func(args)
    except ModuleNotFoundError as exc:
        missing = (exc.name or "").split(".", 1)[0]
        extra = _OPTIONAL_IMPORT_EXTRAS.get(missing)
        if not extra:
            raise
        parser.error(
            f"command '{getattr(args, 'command', '')}' requires the optional "
            f"'{extra}' dependencies; install with: python -m pip install "
            f"'cftk[{extra}]'"
        )


if __name__ == "__main__":
    main()
