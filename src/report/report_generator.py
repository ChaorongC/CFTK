"""
report_generator.py — self-contained HTML report with interactive Plotly.js charts.

Part 1 — Data Processing:
  1.1 Trimming      — Filtered Reads (bar) + Trimmed Sequence Lengths 3' (line)
                      → MultiQC trimming HTML (cutadapt plots)
  1.2 Trimmed QC    — FastQC modules: Sequence Counts, Quality Histograms,
                      Per-Seq Quality, GC, N Content, Length Dist, Duplication,
                      Status Checks
                      → MultiQC trimming HTML (fqc_trimmed plots)
  1.3 Alignment     — Alignment Rates, Deduplication
                      → MultiQC alignment HTML (samtools plots)
  1.4 M-bias        — CpG OT/OB R1/R2 per-position methylation bias
                      → MultiQC methylation HTML (methyldackel/bismark mbias)
  1.6 Sequencing QC Summary — Interactive sortable per-sample metrics table
                      → qc_scores.tsv

Part 2 — QC Analysis:
  2.1 Methylation Distribution
  2.2 Fragment Length Distribution
  2.3 Dinucleotide Frequency
  2.4 PCA

NOTE: Section 2.5 (Sample Correlation) has been removed.
"""

import base64
import datetime
import glob
import html
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from report.mqc_extractor import (
    PLOTLY_SCRIPT,
    merge_mqc_data,
    mqc_bar_chart,
    mqc_heatmap,
    mqc_line_chart,
    mqc_mbias_chart,
    markdup_dedup_chart,
    mbias_tsv_chart,
    mbias_legacy_chart,
)


# ── static image helpers ──────────────────────────────────────────────────────

def _b64(png_path, max_width=2400, jpeg_quality=95):
    try:
        from PIL import Image
        import io
        img = Image.open(png_path)
        if img.mode in ("RGBA", "P", "LA"):
            bg = Image.new("RGB", img.size, (255, 255, 255))
            if img.mode == "P":
                img = img.convert("RGBA")
            alpha = img.split()[-1] if img.mode in ("RGBA", "LA") else None
            bg.paste(img, mask=alpha)
            img = bg
        elif img.mode != "RGB":
            img = img.convert("RGB")
        if img.width > max_width:
            ratio = max_width / img.width
            img   = img.resize((max_width, int(img.height * ratio)), Image.LANCZOS)
        buf = io.BytesIO()
        img.save(buf, format="JPEG", quality=jpeg_quality, optimize=True)
        return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()
    except Exception:
        with open(png_path, "rb") as f:
            return "data:image/png;base64," + base64.b64encode(f.read()).decode()


def _img_tag(png_path, alt="", style="width:100%;", zoomable=True):
    if not png_path or not os.path.exists(png_path) or os.path.isdir(png_path):
        return f'<div class="fig-placeholder"><span class="hint">{alt} (not found)</span></div>'
    src = _b64(png_path)
    img = f'<img src="{src}" alt="{alt}" style="{style}"'
    if zoomable:
        img += ' onclick="openLightbox(this.src)" title="Click to zoom"'
    img += ">"
    return img


def _b64_svg(svg_path):
    if not svg_path or not os.path.exists(svg_path):
        return ""
    with open(svg_path, "rb") as f:
        return "data:image/svg+xml;base64," + base64.b64encode(f.read()).decode()


def _find(results_dir, *parts, pattern="*.png"):
    files = sorted(glob.glob(os.path.join(results_dir, *parts, pattern)))
    return files[0] if files else None


def _find_all(results_dir, *parts, pattern="*.png"):
    return sorted(glob.glob(os.path.join(results_dir, *parts, pattern)))


def _sample_names(groups: dict) -> list[str]:
    return [sample for members in groups.values() for sample in members]


def _source_results_from_config(cfg: dict | None, results_dir: str) -> dict[str, str]:
    """Infer prior CFTK result roots from canonical marked-BAM input paths.

    A downstream-only project may use marked BAMs produced by another CFTK
    project. This deliberately recognizes only the canonical CFTK hierarchy;
    arbitrary BAM inputs do not cause external paths to be scanned.
    """
    if not cfg:
        return {}
    local_root = Path(results_dir).expanduser().resolve()
    source_results: dict[str, str] = {}
    for members in cfg.get("samples", {}).values():
        for sample in members:
            if sample.get("input_type") != "bam" or not sample.get("bam"):
                continue
            bam = Path(sample["bam"]).expanduser().resolve()
            if bam.parent.name != "3_markdup" or bam.parent.parent.name != "1_process":
                continue
            source_root = bam.parent.parent.parent
            if source_root.name != "results" or not source_root.is_dir() or source_root == local_root:
                continue
            source_results[str(sample["name"])] = str(source_root)
    return source_results


def _source_roots(rd: str, source_results: dict[str, str] | None) -> list[str]:
    """Return current results followed by unique inferred source roots."""
    roots = [rd]
    for root in (source_results or {}).values():
        root = str(root)
        if root and root not in roots:
            roots.append(root)
    return roots


def _read_qc_scores(
    rd: str,
    groups: dict,
    source_results: dict[str, str] | None = None,
) -> pd.DataFrame | None:
    """Read one QC-score row per configured sample, preferring current results."""
    samples = _sample_names(groups)
    tables: dict[str, pd.DataFrame] = {}

    def table_at(root: str) -> pd.DataFrame | None:
        if root not in tables:
            path = os.path.join(root, "2_qc", "qc_scores.tsv")
            try:
                tables[root] = pd.read_csv(path, sep="\t") if os.path.exists(path) else pd.DataFrame()
            except Exception:
                tables[root] = pd.DataFrame()
        table = tables[root]
        return table if not table.empty and "sample" in table.columns else None

    rows = []
    for sample in samples:
        candidates = [rd]
        source_root = (source_results or {}).get(sample)
        if source_root and source_root not in candidates:
            candidates.append(source_root)
        for root in candidates:
            table = table_at(root)
            if table is None:
                continue
            matches = table[table["sample"].astype(str) == sample]
            if not matches.empty:
                rows.append(matches.iloc[0].to_dict())
                break
    return pd.DataFrame(rows) if rows else None


def _as_qc_dataframe(scores: str | pd.DataFrame | None) -> pd.DataFrame | None:
    if isinstance(scores, pd.DataFrame):
        return scores.copy()
    if not scores or not os.path.exists(scores):
        return None
    try:
        return pd.read_csv(scores, sep="\t")
    except Exception:
        return None


def _uid(prefix=""):
    import random, string
    safe = re.sub(r"[^a-zA-Z0-9]", "_", prefix)
    return safe + "".join(random.choices(string.ascii_lowercase, k=8))


def _dropdown_gallery(title, images, uid=None):
    if not images:
        return f"<p><em>No {title} results found.</em></p>"
    uid = uid or _uid("dd_")
    uid = re.sub(r"[^a-zA-Z0-9_]", "_", uid)
    js_parts = [f"var v=this.value;"]
    for i in range(len(images)):
        js_parts.append(
            f"document.getElementById('{uid}_p{i}').style.display="
            f"(v==='{i}')?'block':'none';"
        )
    onchange_js = "".join(js_parts)
    opts   = "".join(f'<option value="{i}">{lbl}</option>' for i, (lbl, _) in enumerate(images))
    panels = "".join(
        f'<div id="{uid}_p{i}" style="display:{"block" if i==0 else "none"};">'
        f'<div class="fig-card"><div class="fig-img">{_img_tag(p, lbl)}</div>'
        f'<div class="fig-body"><div class="fig-label">{lbl}</div></div></div></div>'
        for i, (lbl, p) in enumerate(images)
    )
    return (
        f'<div style="margin-bottom:20px;">'
        f'<div style="display:flex;align-items:center;gap:10px;margin-bottom:10px;">'
        f'<label style="font-size:14.5px;font-weight:500;color:#000000;'
        f'font-family:var(--mono);">{title}</label>'
        f'<select onchange="{onchange_js}" '
        f'style="font-size:14.5px;padding:4px 10px;border:1px solid var(--rule);'
        f'border-radius:4px;background:white;color:#000000;cursor:pointer;">'
        f'{opts}</select></div>{panels}</div>'
    )


def _missing(msg: str) -> str:
    return (
        f'<div style="padding:14px 16px;background:var(--bg2);'
        f'border:1px solid var(--rule);border-radius:8px;'
        f'font-size:15.5px;color:#000000;">{msg}</div>'
    )


def _note(msg: str, kind: str = "info") -> str:
    """Inline note box. kind: info | warn | bisulfite"""
    colors = {
        "info":       ("rgba(11,110,153,.07)",  "#0b6e99", "ℹ"),
        "warn":       ("rgba(161,98,7,.08)",    "#7a4a05", "⚠"),
        "bisulfite":  ("rgba(15,118,110,.07)",  "#0f766e", "🧬"),
    }
    bg, fg, icon = colors.get(kind, colors["info"])
    return (
        f'<div style="padding:10px 14px;margin-bottom:12px;'
        f'background:{bg};border-left:3px solid {fg};'
        f'border-radius:0 6px 6px 0;font-size:14.5px;'
        f'color:#3a4a5a;line-height:1.6;">'
        f'<strong style="color:{fg};">{icon} Note:</strong> {msg}</div>'
    )


def _read_json(path):
    try:
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, ValueError, json.JSONDecodeError):
        return None


def _manifest_from_pointer(provenance, pointer_name, run_directory):
    """Load the newest manifest named by a provenance pointer or run directory."""
    pointer = _read_json(os.path.join(provenance, pointer_name))
    if isinstance(pointer, dict) and pointer.get("manifest"):
        manifest_path = pointer["manifest"]
        if not os.path.isabs(manifest_path):
            manifest_path = os.path.join(provenance, manifest_path)
        manifest = _read_json(manifest_path)
        if isinstance(manifest, dict):
            return manifest

    candidates = glob.glob(os.path.join(provenance, run_directory, "*", "run.json"))
    if not candidates:
        return None
    return _read_json(max(candidates, key=os.path.getmtime))


def _analysis_stage_records(provenance):
    """Collect the newest trusted record for every downstream workflow stage."""
    candidates = glob.glob(os.path.join(provenance, "analysis-runs", "*", "run.json"))
    records = {}
    trusted = {"complete", "resumed", "adopted"}
    for path in sorted(candidates, key=os.path.getmtime, reverse=True):
        manifest = _read_json(path)
        if not isinstance(manifest, dict):
            continue
        for stage in manifest.get("stages", []):
            stage_id = stage.get("id")
            if not stage_id or stage_id in records or stage.get("status") not in trusted:
                continue
            records[stage_id] = {"stage": stage, "run_id": manifest.get("run_id", "")}
    return records


def _artifact_counts(stage):
    """Return reported output/figure counts and currently missing required artifacts."""
    expected = stage.get("expected", [])
    outputs = stage.get("outputs")
    figures = stage.get("figures")
    if outputs is None:
        outputs = [item for item in expected if item.get("role") != "figure"]
    if figures is None:
        figures = [item for item in expected if item.get("role") == "figure"]

    missing = 0
    for item in expected:
        if not item.get("required", True):
            continue
        path = item.get("path")
        if not path or not os.path.exists(path):
            missing += 1
        elif item.get("nonempty", True) and os.path.isfile(path) and os.path.getsize(path) == 0:
            missing += 1
    return len(outputs), len(figures), missing


def _workflow_status(status):
    colors = {
        "complete": ("#e7f5ed", "#167346"),
        "resumed": ("#e3f0f6", "#0b6e99"),
        "adopted": ("#f5edff", "#6b3fa0"),
        "failed": ("#fce8e6", "#b42318"),
        "interrupted": ("#fff4e5", "#9a6700"),
    }
    background, color = colors.get(str(status).lower(), ("#eceff3", "#46586b"))
    label = html.escape(str(status or "unknown"))
    return (
        f'<span style="display:inline-block;padding:2px 8px;border-radius:10px;'
        f'background:{background};color:{color};font-size:12.5px;font-weight:700;">{label}</span>'
    )


def _sec_workflow_summary(rd, groups):
    """Render a manifest-aware overview of the complete CFTK workflow."""
    provenance = os.path.join(rd, "provenance")
    core = _manifest_from_pointer(provenance, "latest-run.json", "runs")
    analysis_records = _analysis_stage_records(provenance)
    rows = []

    if isinstance(core, dict):
        for stage in core.get("stages", []):
            outputs, figures, missing = _artifact_counts(stage)
            rows.append((
                "Core processing / QC", stage.get("name", stage.get("id", "")),
                stage.get("status", "unknown"), outputs, figures, missing,
                core.get("run_id", ""),
            ))
    for stage_id, record in analysis_records.items():
        stage = record["stage"]
        outputs, figures, missing = _artifact_counts(stage)
        rows.append((
            "Downstream", stage.get("name", stage_id), stage.get("status", "unknown"),
            outputs, figures, missing, record["run_id"],
        ))

    if not rows:
        content = _missing(
            "No CFTK run manifest was found. The report still discovers results from standard "
            "output directories, but it cannot display manifest-backed stage status or provenance."
        )
    else:
        table_rows = "".join(
            "<tr style=\"border-top:1px solid var(--rule);\">"
            f"<td style=\"padding:8px 10px;\">{html.escape(workflow)}</td>"
            f"<td style=\"padding:8px 10px;color:var(--ink);font-weight:600;\">{html.escape(name)}</td>"
            f"<td style=\"padding:8px 10px;text-align:center;\">{_workflow_status(status)}</td>"
            f"<td style=\"padding:8px 10px;text-align:center;\">{outputs}</td>"
            f"<td style=\"padding:8px 10px;text-align:center;\">{figures}</td>"
            f"<td style=\"padding:8px 10px;text-align:center;\">{missing}</td>"
            f"<td style=\"padding:8px 10px;font-family:var(--mono);font-size:12px;\">{html.escape(run_id)}</td>"
            "</tr>"
            for workflow, name, status, outputs, figures, missing, run_id in rows
        )
        content = (
            '<div style="overflow-x:auto;border:1px solid var(--rule);border-radius:8px;">'
            '<table style="width:100%;border-collapse:collapse;font-size:14px;min-width:880px;">'
            '<thead><tr>'
            '<th>Workflow</th><th>Stage</th><th>Status</th><th>Outputs</th>'
            '<th>Figures</th><th>Missing required</th><th>Source run</th>'
            '</tr></thead><tbody>' + table_rows + '</tbody></table></div>'
        )

    scope_note = ""
    for record in analysis_records.values():
        scope = record["stage"].get("fragmentomics_scope")
        if isinstance(scope, dict) and scope.get("note"):
            scope_note = _note(html.escape(str(scope["note"])))
            break
    sample_count = sum(len(members) for members in groups.values())
    return f"""
    <section class="section" id="workflow_summary">
      <div class="section-header">
        <span class="section-num">RUN</span>
        <h2 class="section-title">Workflow Summary</h2>
        <span class="section-tag">{sample_count} samples</span>
      </div>
      <p style="font-size:14.5px;color:#000000;margin-bottom:14px;">
        Automatically discovered CFTK core-processing, QC, and downstream-stage records. Counts are
        checked against the current files so this page distinguishes completed, resumed, adopted, and
        missing artifacts.
      </p>
      {scope_note}
      {content}
    </section>"""


# ── collapsible section helper ────────────────────────────────────────────────

_CC = [0]


def _coll(title: str, content: str, open_: bool = True) -> str:
    _CC[0] += 1
    uid  = f"coll_{_CC[0]}"
    disp = "block" if open_ else "none"
    return (
        f'<div style="margin-bottom:16px;border:1px solid var(--rule);'
        f'border-radius:8px;overflow:hidden;">'
        f'<div onclick="var e=document.getElementById(\'{uid}\');'
        f'e.style.display=e.style.display===\'none\'?\'block\':\'none\';" '
        f'style="padding:10px 16px;background:var(--bg2);cursor:pointer;'
        f'display:flex;align-items:center;justify-content:space-between;'
        f'font-size:15.5px;font-weight:600;color:var(--ink);user-select:none;">'
        f'<span>{title}</span>'
        f'<span style="font-size:19px;color:#000000;">⌄</span></div>'
        f'<div id="{uid}" style="display:{disp};padding:16px;">{content}</div></div>'
    )


# ── QC Summary table (1.6) ────────────────────────────────────────────────────

_STATUS_BG = {"PASS": "rgba(46,204,113,.15)", "WARN": "rgba(243,156,18,.15)",
              "FAIL": "rgba(231,76,60,.15)",   "NA":   "rgba(189,195,199,.15)"}
_STATUS_FG = {"PASS": "#1a7a42", "WARN": "#9a6200", "FAIL": "#a93226", "NA": "#7f8c8d"}
_STATUS_BD = {"PASS": "rgba(46,204,113,.5)",  "WARN": "rgba(243,156,18,.5)",
              "FAIL": "rgba(231,76,60,.5)",    "NA":   "rgba(189,195,199,.5)"}


def _badge(status: str) -> str:
    st = str(status).upper() if status else "NA"
    return (
        f'<span style="display:inline-block;padding:2px 9px;border-radius:10px;'
        f'font-size:13px;font-weight:600;'
        f'background:{_STATUS_BG.get(st, _STATUS_BG["NA"])};'
        f'color:{_STATUS_FG.get(st, _STATUS_FG["NA"])};'
        f'border:1px solid {_STATUS_BD.get(st, _STATUS_BD["NA"])};">{st}</span>'
    )


def _score_bar(score: float, status: str, width: int = 80) -> str:
    pct   = min(max(float(score), 0), 100)
    color = _STATUS_FG.get(str(status).upper(), _STATUS_FG["NA"])
    bar_w = int(pct / 100 * width)
    return (
        f'<svg width="{width+38}" height="16" style="vertical-align:middle;">'
        f'<rect x="0" y="3" width="{width}" height="10" rx="5" fill="rgba(0,0,0,.08)"/>'
        f'<rect x="0" y="3" width="{bar_w}" height="10" rx="5" '
        f'fill="{color}" opacity=".75"/>'
        f'<text x="{width+4}" y="12" font-size="10" fill="{color}" '
        f'font-family="monospace" font-weight="600">{pct:.0f}</text>'
        f'</svg>'
    )


def _fmt(val, fmt=".1f", unit="%") -> str:
    try:
        v = float(val)
        return "—" if np.isnan(v) else f"{v:{fmt}}{unit}"
    except Exception:
        return str(val) if val else "—"


def _cell(val, status: str, fmt=".1f", unit="%") -> str:
    st  = str(status).upper() if status else "NA"
    txt = _fmt(val, fmt, unit)
    try:
        sv = float(val) if not np.isnan(float(val)) else -999
    except Exception:
        sv = -999
    return (
        f'<td style="text-align:center;padding:6px 8px;'
        f'background:{_STATUS_BG.get(st, _STATUS_BG["NA"])};" data-sort="{sv}">'
        f'<span style="color:{_STATUS_FG.get(st, _STATUS_FG["NA"])};'
        f'font-weight:500;font-size:14.5px;">{txt}</span></td>'
    )


def _qc_table(scores: str | pd.DataFrame | None) -> str:
    df = _as_qc_dataframe(scores)
    if df is None or df.empty:
        return _missing("QC scores not found — run <code>cftk qc -s 0</code> first")
    try:
        from analysis.qc_scorer import RULES
        active_rules = [r for r in RULES if r.weight > 0]
    except Exception as e:
        return _missing(f"Could not load QC scores: {e}")

    th_style = (
        "padding:8px 10px;background:var(--header-bg);color:var(--header-ink);"
        "font-size:13px;font-weight:700;letter-spacing:.02em;cursor:pointer;white-space:nowrap;"
        "border-bottom:2px solid var(--rule);"
    )
    ths = (
        f'<th style="{th_style}text-align:left;min-width:130px;" onclick="qcSort(0)">Sample</th>'
        f'<th style="{th_style}text-align:center;" onclick="qcSort(1)">Status</th>'
    )
    rule_col_map, col_idx = [], 2
    for rule in active_rules:
        val_col = f"{rule.col}_value"
        st_col  = f"{rule.col}_status"
        if val_col not in df.columns:
            continue
        rule_col_map.append((rule, col_idx, val_col, st_col))
        ths += (
            f'<th style="{th_style}text-align:center;" '
            f'title="{rule.note}" onclick="qcSort({col_idx})">{rule.label}</th>'
        )
        col_idx += 1

    tbody = ""
    for _, row in df.iterrows():
        sample = str(row.get("sample", ""))
        status = str(row.get("qc_status", "NA")).upper()
        score  = row.get("qc_score", 0)
        group  = str(row.get("group", ""))
        sk     = re.sub(r"[^a-zA-Z0-9_-]", "_", sample)
        rec    = str(row.get("recommendation", ""))

        td_s = (
            f'<td style="padding:6px 12px;font-family:var(--mono);font-size:14.5px;'
            f'color:var(--ink);font-weight:500;white-space:nowrap;" data-sort="{sample}">'
            f'{sample}</td>'
        )
        td_st = (
            f'<td style="text-align:center;padding:6px 8px;" data-sort="{status}">'
            f'{_badge(status)}</td>'
        )
        metric_tds = "".join(
            _cell(row.get(val_col, float("nan")), str(row.get(st_col, "NA")),
                  rule.fmt, rule.unit)
            for rule, _, val_col, st_col in rule_col_map
        )
        tbody += (
            f'<tr onclick="qcToggleRec(\'{sk}\')" '
            f'style="cursor:pointer;border-top:1px solid var(--rule);" '
            f'onmouseover="this.style.background=\'var(--bg2)\'" '
            f'onmouseout="this.style.background=\'\'">'
            f'{td_s}{td_st}{metric_tds}</tr>'
            f'<tr id="qcrec_row_{sk}" style="display:none;">'
            f'<td colspan="{col_idx}" style="padding:0 12px 10px;">'
            f'<div style="padding:10px 14px;margin-top:4px;'
            f'border-left:4px solid {_STATUS_FG.get(status,"#aaa")};'
            f'background:{_STATUS_BG.get(status,"#f9f9f9")};'
            f'border-radius:0 6px 6px 0;font-size:14.5px;line-height:1.7;">'
            f'{rec}</div></td></tr>'
        )

    n_pass = (df.get("qc_status", pd.Series()) == "PASS").sum()
    n_warn = (df.get("qc_status", pd.Series()) == "WARN").sum()
    n_fail = (df.get("qc_status", pd.Series()) == "FAIL").sum()
    btn = "padding:3px 12px;border-radius:10px;cursor:pointer;font-size:14.5px;border:1px solid;"
    filter_bar = (
        f'<div style="display:flex;align-items:center;gap:8px;margin-bottom:12px;flex-wrap:wrap;">'
        f'<span style="font-size:14.5px;color:#000000;font-family:var(--mono);">Filter:</span>'
        f'<button class="qc-fbtn" data-status="" onclick="qcFilter(\'\')" '
        f'style="{btn}background:var(--bg2);border-color:var(--rule);color:var(--ink);font-weight:700;">'
        f'All ({len(df)})</button>'
        f'<button class="qc-fbtn" data-status="PASS" onclick="qcFilter(\'PASS\')" '
        f'style="{btn}background:{_STATUS_BG["PASS"]};border-color:{_STATUS_BD["PASS"]};color:{_STATUS_FG["PASS"]};">'
        f'✓ PASS ({n_pass})</button>'
        f'<button class="qc-fbtn" data-status="WARN" onclick="qcFilter(\'WARN\')" '
        f'style="{btn}background:{_STATUS_BG["WARN"]};border-color:{_STATUS_BD["WARN"]};color:{_STATUS_FG["WARN"]};">'
        f'⚠ WARN ({n_warn})</button>'
        f'<button class="qc-fbtn" data-status="FAIL" onclick="qcFilter(\'FAIL\')" '
        f'style="{btn}background:{_STATUS_BG["FAIL"]};border-color:{_STATUS_BD["FAIL"]};color:{_STATUS_FG["FAIL"]};">'
        f'✕ FAIL ({n_fail})</button>'
        f'<span style="margin-left:auto;font-size:13px;color:#000000;">'
        f'Click row for recommendation · Click header to sort</span></div>'
    )
    js = """<script>
var _qcSD={};
function qcSort(col){
  var tbl=document.getElementById('qc_stats_tbl');
  var rows=Array.from(tbl.tBodies[0].rows).filter(function(r){return !r.id.startsWith('qcrec_row_');});
  var dir=(_qcSD[col]==='asc')?'desc':'asc'; _qcSD[col]=dir;
  rows.sort(function(a,b){
    var av=a.cells[col]?a.cells[col].getAttribute('data-sort'):'';
    var bv=b.cells[col]?b.cells[col].getAttribute('data-sort'):'';
    var an=parseFloat(av),bn=parseFloat(bv);
    if(!isNaN(an)&&!isNaN(bn)) return dir==='asc'?an-bn:bn-an;
    return dir==='asc'?av.localeCompare(bv):bv.localeCompare(av);
  });
  var tb=tbl.tBodies[0];
  rows.forEach(function(r){
    var k=r.cells[0].getAttribute('data-sort').replace(/[^a-zA-Z0-9_-]/g,'_');
    var rr=document.getElementById('qcrec_row_'+k);
    tb.appendChild(r); if(rr) tb.appendChild(rr);
  });
}
function qcToggleRec(s){var r=document.getElementById('qcrec_row_'+s);if(r)r.style.display=r.style.display===''?'none':'';}
function qcFilter(st){
  var rows=document.querySelectorAll('#qc_stats_tbl tbody tr:not([id^="qcrec_row_"])');
  rows.forEach(function(r){
    var s=r.cells[1]?r.cells[1].getAttribute('data-sort'):'';
    var k=r.cells[0].getAttribute('data-sort').replace(/[^a-zA-Z0-9_-]/g,'_');
    var rr=document.getElementById('qcrec_row_'+k);
    var show=!st||s===st; r.style.display=show?'':'none';
    if(rr)rr.style.display='none';
  });
  document.querySelectorAll('.qc-fbtn').forEach(function(b){
    b.style.fontWeight=b.getAttribute('data-status')===st?'700':'400';
    b.style.opacity=b.getAttribute('data-status')===st?'1':'0.65';
  });
}
</script>"""
    # Always-visible horizontal scrollbar pinned at the bottom of the table.
    # The wrapper's native scrollbar is hidden; a synced slim slider drives it.
    bottom_scrollbar = (
        '<style>#qc_tbl_wrap::-webkit-scrollbar{height:0}'
        '#qc_sb_wrap::-webkit-scrollbar{height:12px}'
        '#qc_sb_wrap::-webkit-scrollbar-thumb{background:#c4ccd4;border-radius:6px}'
        '#qc_sb_wrap::-webkit-scrollbar-track{background:var(--bg2);border-radius:6px}'
        '</style>'
        '<div id="qc_sb_wrap" style="overflow-x:auto;overflow-y:hidden;'
        'height:14px;margin-top:4px;border:1px solid var(--rule);border-radius:7px;'
        'background:var(--bg2);">'
        '<div id="qc_sb_inner" style="height:1px;width:1px;"></div></div>'
        '<script>(function(){'
        'var w=document.getElementById("qc_tbl_wrap");'
        'var sb=document.getElementById("qc_sb_wrap");'
        'var inner=document.getElementById("qc_sb_inner");'
        'if(!w||!sb||!inner)return;'
        'function sync(){inner.style.width=w.scrollWidth+"px";'
        'sb.style.display=(w.scrollWidth>w.clientWidth)?"block":"none";}'
        'sync();window.addEventListener("resize",sync);setTimeout(sync,200);'
        'var lock=false;'
        'sb.addEventListener("scroll",function(){if(lock)return;lock=true;w.scrollLeft=sb.scrollLeft;lock=false;});'
        'w.addEventListener("scroll",function(){if(lock)return;lock=true;sb.scrollLeft=w.scrollLeft;lock=false;});'
        '})();</scr'+'ipt>'
    )
    return (
        filter_bar
        + '<div id="qc_tbl_wrap" style="overflow-x:auto;overflow-y:hidden;'
        'scrollbar-width:none;-ms-overflow-style:none;'
        'border:1px solid var(--rule);border-radius:8px;">'
        f'<table id="qc_stats_tbl" style="width:100%;border-collapse:collapse;font-size:15.5px;">'
        f'<thead><tr>{ths}</tr></thead><tbody>{tbody}</tbody></table></div>'
        + bottom_scrollbar
        + js
    )


# ── remaining original helpers ────────────────────────────────────────────────

def _read_perf_table(tsv_path):
    if not tsv_path or not os.path.exists(tsv_path):
        return "<p><em>No modality performance data.</em></p>"
    try:
        import ast, re as _re
        df = pd.read_csv(tsv_path, sep="\t", index_col=0)
        rows = []
        for idx, row in df.iterrows():
            best_auc     = row.get("best_roc_auc_mean", "—")
            best_clf_raw = row.get("best_classifier, idx", "—")
            m            = _re.search(r"'([^']+)'", str(best_clf_raw))
            clf_name     = (m.group(1) if m else str(best_clf_raw)
                           ).replace("Classifier", "").replace("Regression", "Reg")
            mean_col = [c for c in df.columns if "mean" in c and "best" not in c]
            if mean_col:
                try:
                    vals     = ast.literal_eval(str(row[mean_col[0]]))
                    clf_aucs = " / ".join(f"{v:.3f}" for v in vals)
                except Exception:
                    clf_aucs = str(row[mean_col[0]])
            else:
                clf_aucs = "—"
            rows.append({"Modality": idx, "Clf AUCs (RFC/LR/SVC)": clf_aucs,
                         "Best Classifier": clf_name,
                         "Best AUC": f"{best_auc:.4f}" if isinstance(best_auc, float) else best_auc})

        cols = ["Modality", "Clf AUCs (RFC/LR/SVC)", "Best Classifier", "Best AUC"]
        th_p = ("padding:9px 14px;background:var(--header-bg);color:var(--header-ink);"
                "font-size:13px;font-weight:700;letter-spacing:.02em;text-align:left;"
                "white-space:nowrap;border-bottom:2px solid var(--rule);")
        head = "".join(f'<th style="{th_p}">{c}</th>' for c in cols)
        body = ""
        for r in rows:
            tds = ""
            for c in cols:
                val = r.get(c, "—")
                weight = "font-weight:600;" if c == "Modality" else ""
                tds += (f'<td style="padding:7px 14px;color:var(--ink);font-size:13.5px;'
                        f'white-space:nowrap;border-top:1px solid var(--rule);{weight}">{val}</td>')
            body += f"<tr>{tds}</tr>"
        return (
            '<div style="overflow-x:auto;border:1px solid var(--rule);border-radius:8px;">'
            '<table style="width:100%;border-collapse:collapse;font-size:13.5px;">'
            f'<thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'
        )
    except Exception as e:
        return f"<p><em>Could not load table: {e}</em></p>"


# ── section builders ──────────────────────────────────────────────────────────

# refined, low-saturation group palette (light bg, strong text) — replaces the
# previous pastel full-row fills. Matches the template accent variables.
_GRP_CHIP = [
    ("#e3f0f6", "#0b6e99", "rgba(11,110,153,.28)"),   # cyan / accent
    ("#d7efec", "#0f766e", "rgba(15,118,110,.28)"),    # teal
    ("#ece9f8", "#5b4baf", "rgba(91,75,175,.28)"),     # purple
    ("#fbf0d9", "#a16207", "rgba(161,98,7,.28)"),      # amber
    ("#fdebe0", "#c2410c", "rgba(194,65,12,.28)"),     # coral
    ("#e7eef5", "#2f5a8a", "rgba(47,90,138,.28)"),     # slate-blue
]


def _grp_chip(grp: str, gi: int) -> str:
    bg, fg, bd = _GRP_CHIP[gi % len(_GRP_CHIP)]
    return (
        f'<span style="display:inline-block;padding:2px 10px;border-radius:11px;'
        f'font-size:13px;font-weight:600;letter-spacing:.01em;'
        f'background:{bg};color:{fg};border:1px solid {bd};white-space:nowrap;">'
        f'{grp}</span>'
    )


def _stat_card(label: str, value: str, accent: str = "var(--accent)") -> str:
    return (
        f'<div style="flex:0 0 auto;min-width:104px;padding:11px 16px;'
        f'background:#fff;border:1px solid var(--rule);border-radius:8px;'
        f'box-shadow:var(--shadow-sm);">'
        f'<div style="font-size:12px;letter-spacing:.08em;text-transform:uppercase;'
        f'color:#000000;font-weight:600;margin-bottom:5px;">{label}</div>'
        f'<div style="font-size:26.5px;font-weight:700;color:{accent};'
        f'font-family:var(--mono);line-height:1;">{value}</div></div>'
    )


def _stat_heat_cell(val, vmin, vmax, status=None) -> str:
    """A numeric cell with a MultiQC-style proportional fill bar behind the value."""
    try:
        v = float(val)
        if np.isnan(v):
            raise ValueError
        txt = f"{v:g}"
    except Exception:
        return ('<td style="padding:6px 10px;text-align:right;color:#000000;'
                'font-family:var(--mono);font-size:14.5px;">—</td>')
    if vmax > vmin:
        pct = max(0.0, min(1.0, (v - vmin) / (vmax - vmin))) * 100
    else:
        pct = 0.0
    # tint by status if provided, else neutral accent
    st = str(status).upper() if status else ""
    fill = {
        "PASS": "rgba(38,140,80,.16)", "WARN": "rgba(176,118,12,.18)",
        "FAIL": "rgba(176,52,42,.16)",
    }.get(st, "rgba(11,110,153,.12)")
    return (
        f'<td style="padding:6px 10px;text-align:right;position:relative;'
        f'white-space:nowrap;">'
        f'<div style="position:absolute;left:0;top:3px;bottom:3px;width:{pct:.1f}%;'
        f'background:{fill};border-radius:0 3px 3px 0;"></div>'
        f'<span style="position:relative;font-family:var(--mono);font-size:14.5px;'
        f'color:var(--ink);">{txt}</span></td>'
    )


def _fragment_peaks(rd, samples, source_results=None):
    """Return {sample: peak_bp} read from 2_qc/2_fragment_length/*.raw.csv."""
    out = {}
    name_to_path = {}
    for root in _source_roots(rd, source_results):
        frag_dir = os.path.join(root, "2_qc", "2_fragment_length")
        for fp in sorted(glob.glob(os.path.join(frag_dir, "fragment_length.*.raw.csv"))):
            stem = os.path.basename(fp)
            name = (stem.replace("fragment_length.", "")
                        .replace(".raw.csv", "")
                        .replace(".markdup", ""))
            name_to_path.setdefault(name, fp)
    for s in samples:
        fp = name_to_path.get(s)
        if not fp:
            for nm, p in name_to_path.items():
                if nm == s or s in nm or nm in s:
                    fp = p
                    break
        if not fp:
            continue
        try:
            t = pd.read_table(fp, skiprows=1).iloc[:, :2]
            t.columns = ["Size", "Occurrences"]
            t = t[pd.to_numeric(t["Size"], errors="coerce") < 500].dropna()
            if t.empty:
                continue
            out[s] = int(t.loc[t["Occurrences"].astype(float).idxmax(), "Size"])
        except Exception:
            pass
    return out


def _sec_sample_statistics(rd, groups, source_results=None, qc_scores=None):
    """
    MultiQC-style per-sample statistics table. Pulls every metric available in
    2_qc/qc_scores.tsv and splits it into two column groups:
      • Sequencing & Methylation QC  — all scorer RULES (excluding beta_* density)
      • cfDNA QC                     — beta-density metrics (M-score, balance, ratio)
    Falls back to a clean grouped sample list if qc_scores.tsv is absent.
    """
    n_total = sum(len(v) for v in groups.values())
    grp_index = {g: i for i, g in enumerate(groups.keys())}
    sample_to_group = {s: g for g, m in groups.items() for s in m}

    # fragment-length peak (bp) per sample, computed from the QC raw files
    all_samples = _sample_names(groups)
    frag_peaks = _fragment_peaks(rd, all_samples, source_results=source_results)
    has_peak = bool(frag_peaks)

    df = _as_qc_dataframe(qc_scores)

    # ── fallback: no QC scores → clean grouped list ──────────────────────────
    if df is None or df.empty:
        rows = ""
        for g, members in groups.items():
            for s in members:
                peak = frag_peaks.get(s)
                peak_td = (f'<td style="padding:7px 14px;text-align:center;color:var(--ink);'
                           f'font-size:15px;">{(str(peak)+" bp") if peak else "—"}</td>'
                           ) if has_peak else ""
                rows += (
                    f'<tr style="border-top:1px solid var(--rule);">'
                    f'<td style="padding:7px 14px;color:var(--ink);font-size:15px;'
                    f'white-space:nowrap;">{s}</td>'
                    f'<td style="padding:7px 14px;">{_grp_chip(g, grp_index[g])}</td>'
                    f'{peak_td}</tr>'
                )
        fb_th = ("padding:9px 14px;background:var(--header-bg);color:var(--header-ink);"
                 "font-size:13px;font-weight:700;letter-spacing:.02em;border-bottom:2px solid var(--rule);"
                 "position:sticky;top:0;")
        peak_th = (f'<th style="{fb_th}text-align:center;">Fragment length peak</th>'
                   if has_peak else "")
        return f"""
    <section class="section" id="part_overview">
      <div class="section-header">
        <span class="section-num">STATS</span>
        <h2 class="section-title">Sample Statistics</h2>
        <span class="section-tag">{n_total} samples · {len(groups)} groups</span>
      </div>
      {_missing("Per-sample QC metrics not found (<code>2_qc/qc_scores.tsv</code>). "
                "Run <code>cftk qc -s 0</code> to populate sequencing &amp; cfDNA QC statistics.")}
      <div style="max-height:382px;overflow:auto;border:1px solid var(--rule);border-radius:8px;margin-top:8px;">
        <table style="width:100%;border-collapse:collapse;font-size:15.5px;">
          <thead><tr>
            <th style="{fb_th}text-align:left;">Sample</th>
            <th style="{fb_th}text-align:left;">Group</th>
            {peak_th}
          </tr></thead><tbody>{rows}</tbody></table>
      </div>
    </section>"""

    # ── QC scores present → full statistics table ─────────────────────────────
    # resolve sequencing-QC rule columns (exclude beta_* → those are cfDNA QC,
    # and exclude the metrics the user asked to drop)
    _SEQ_EXCLUDE = {"properly paired", "sequencing error rate",
                    "global cpg methylation", "cpg covered sites"}

    def _norm(s):
        return " ".join(str(s).lower().split())

    seq_cols = []   # (label, value_col, status_col, fmt, unit)
    try:
        from analysis.qc_scorer import RULES
        for rule in RULES:
            if getattr(rule, "weight", 0) <= 0:
                continue
            if str(rule.col).startswith("beta_"):
                continue
            if _norm(rule.label) in _SEQ_EXCLUDE:
                continue
            vcol = f"{rule.col}_value"
            if vcol in df.columns:
                seq_cols.append((rule.label, vcol, f"{rule.col}_status",
                                 getattr(rule, "fmt", ".1f"), getattr(rule, "unit", "")))
    except Exception:
        # RULES unavailable — infer from *_value columns (excluding beta_)
        for c in df.columns:
            if c.endswith("_value") and not c.startswith("beta_"):
                base = c[:-6]
                label = base.replace("_", " ").title()
                if _norm(label) in _SEQ_EXCLUDE:
                    continue
                seq_cols.append((label, c, f"{base}_status", ".2f", ""))

    # cfDNA-QC columns: fragment-length peak (value) + β M-score (PASS badge only)
    has_beta = ("beta_M_score_status" in df.columns
                or "beta_M_score_value" in df.columns or "beta_M_score" in df.columns)
    beta_vcol = ("beta_M_score_value" if "beta_M_score_value" in df.columns
                 else ("beta_M_score" if "beta_M_score" in df.columns else None))
    span_cf = (1 if has_peak else 0) + (1 if has_beta else 0)

    def _beta_badge(row):
        st = row.get("beta_M_score_status")
        if st is None or (isinstance(st, float) and np.isnan(st)) or str(st).strip() == "":
            if beta_vcol is not None:
                try:
                    st = "PASS" if float(row.get(beta_vcol)) >= 2 else "WARN"
                except Exception:
                    st = "NA"
            else:
                st = "NA"
        return _badge(str(st).upper())

    # header styles — sticky, NOT all-caps
    th = ("padding:9px 12px;background:var(--header-bg);color:var(--header-ink);"
          "font-size:14px;font-weight:700;letter-spacing:.01em;white-space:nowrap;"
          "border-bottom:2px solid var(--rule);position:sticky;top:28px;z-index:2;")
    grp_th = ("padding:0 12px;height:28px;background:var(--header-bg);"
              "color:var(--header-accent);font-size:15px;letter-spacing:.02em;"
              "font-weight:700;text-align:center;border-bottom:1px solid var(--rule);"
              "position:sticky;top:0;z-index:3;")

    # ── grouped header (two sticky rows) ──────────────────────────────────────
    top_row = (
        f'<th style="{grp_th}text-align:center;" colspan="2">Sample</th>'
        + (f'<th style="{grp_th}border-left:2px solid var(--rule);" colspan="{len(seq_cols)}">'
           f'Sequencing &amp; methylation QC</th>' if seq_cols else "")
        + (f'<th style="{grp_th}border-left:2px solid var(--rule);" colspan="{span_cf}">'
           f'cfDNA QC</th>' if span_cf else "")
    )
    sub_row = (
        f'<th style="{th}text-align:center;">Sample</th>'
        f'<th style="{th}text-align:center;">Group</th>'
    )
    for j, (label, _, _, _, _) in enumerate(seq_cols):
        bl = "border-left:2px solid var(--rule);" if j == 0 else ""
        sub_row += f'<th style="{th}text-align:center;{bl}">{label}</th>'
    cf_first = True
    if has_peak:
        sub_row += (f'<th style="{th}text-align:center;border-left:2px solid var(--rule);">'
                    f'Fragment length peak</th>')
        cf_first = False
    if has_beta:
        bl = "border-left:2px solid var(--rule);" if cf_first else ""
        sub_row += f'<th style="{th}text-align:center;{bl}">β M-score</th>'

    # ── body — centered values, no bars ───────────────────────────────────────
    def _num_cell(val, bl=""):
        try:
            v = float(val)
            if np.isnan(v):
                raise ValueError
            txt = f"{v:g}"
        except Exception:
            txt = "—"
        return (f'<td style="padding:6px 12px;text-align:center;color:var(--ink);'
                f'font-size:15px;white-space:nowrap;{bl}">{txt}</td>')

    tbody = ""
    for _, row in df.iterrows():
        sample = str(row.get("sample", ""))
        grp    = str(row.get("group", "")) or sample_to_group.get(sample, "")
        gi     = grp_index.get(grp, 0)

        td_sample = (f'<td style="padding:6px 12px;color:var(--ink);font-size:15px;'
                     f'font-weight:600;text-align:center;white-space:nowrap;">{sample}</td>')
        td_grp    = f'<td style="padding:6px 12px;text-align:center;">{_grp_chip(grp, gi) if grp else "—"}</td>'

        seq_tds = ""
        for j, (_, vc, _, _, _) in enumerate(seq_cols):
            bl = "border-left:2px solid var(--rule);" if j == 0 else ""
            seq_tds += _num_cell(row.get(vc), bl)

        cf_tds = ""
        cf_first = True
        if has_peak:
            peak = frag_peaks.get(sample)
            txt = f"{peak} bp" if peak else "—"
            cf_tds += (f'<td style="padding:6px 12px;text-align:center;color:var(--ink);'
                       f'font-size:15px;white-space:nowrap;border-left:2px solid var(--rule);">'
                       f'{txt}</td>')
            cf_first = False
        if has_beta:
            bl = "border-left:2px solid var(--rule);" if cf_first else ""
            cf_tds += f'<td style="padding:6px 12px;text-align:center;{bl}">{_beta_badge(row)}</td>'

        tbody += (
            f'<tr style="border-top:1px solid var(--rule);" '
            f'onmouseover="this.style.background=\'var(--bg2)\'" '
            f'onmouseout="this.style.background=\'\'">'
            f'{td_sample}{td_grp}{seq_tds}{cf_tds}</tr>'
        )

    hint = ""

    return f"""
    <section class="section" id="part_overview">
      <div class="section-header">
        <span class="section-num">STATS</span>
        <h2 class="section-title">Sample Statistics</h2>
        <span class="section-tag">{n_total} samples · {len(groups)} groups</span>
      </div>
      <div style="max-height:382px;overflow:auto;border:1px solid var(--rule);border-radius:8px;
           box-shadow:var(--shadow-sm);">
        <table style="width:100%;border-collapse:collapse;font-size:15.5px;">
          <thead>
            <tr>{top_row}</tr>
            <tr>{sub_row}</tr>
          </thead>
          <tbody>{tbody}</tbody>
        </table>
      </div>
      {hint}
    </section>"""


def _sec_power(rd):
    img = _img_tag(_find(rd, "0_power", pattern="power_cumulative.png"), "Power curve")
    return f"""
    <section class="section" id="part0">
      <div class="section-header">
        <span class="section-num">PART 00</span>
        <h2 class="section-title">Power Analysis</h2>
        <span class="section-tag">power</span>
      </div>
      <div class="fig-card">
        <div class="fig-img tall">{img}</div>
        <div class="fig-body">
          <div class="fig-label">Figure 0A · Power vs CpG proportion</div>
        </div>
      </div>
    </section>"""


def _mqc_key(mqc_data: dict, *candidates: str) -> str:
    """Return the first candidate key present in mqc_data, or the last candidate as default."""
    for k in candidates:
        if k in mqc_data:
            return k
        suffixed = sorted(
            key for key in mqc_data
            if key.startswith(f"{k}-") and key[len(k) + 1:].isdigit()
        )
        if suffixed:
            return suffixed[0]
    return candidates[-1]   # last = used in _missing() error message


def _sec_process(
    rd,
    groups=None,
    all_sample_names=None,
    source_results=None,
    qc_scores=None,
):
    """
    Part 1 — Data Processing.
    All interactive charts sourced from MultiQC HTML files via mqc_extractor.
    """
    groups           = groups or {}
    all_sample_names = all_sample_names or [s for m in groups.values() for s in m]

    result_roots = _source_roots(rd, source_results)

    def _multiqc_paths(step):
        return [
            os.path.join(root, "1_process", step, "multiqc", "multiqc_report.html")
            for root in result_roots
            if os.path.exists(os.path.join(root, "1_process", step, "multiqc", "multiqc_report.html"))
        ]

    trim_data = merge_mqc_data(_multiqc_paths("1_trimming"))
    align_data = merge_mqc_data(_multiqc_paths("2_alignment"))
    meth_data = merge_mqc_data(_multiqc_paths("4_methylation"))

    def _no_mqc(step_dir: str) -> str:
        return _missing(
            f"MultiQC report not found at <code>{step_dir}</code>. "
            f"Run <code>cftk process -s {step_dir.split('_')[0][-1]}</code> "
            f"and ensure MultiQC is installed."
        )

    # ── 1.1 Trimming ──────────────────────────────────────────────────────────
    # cutadapt plot key for trimmed lengths may vary: _plot_3 (3' adapter) or just _plot
    filtered_chart = mqc_bar_chart(
        _mqc_key(trim_data,
            "cutadapt_filtered_reads_plot"),
        trim_data,
        title="Filtered Reads", yaxis="Reads (M)")

    # cutadapt_trimmed_sequences_plot_3: datasets[0]=Counts, datasets[1]=Obs/Exp (confirmed)
    lengths_chart = mqc_line_chart(
        "cutadapt_trimmed_sequences_plot_3", trim_data,
        title="Trimmed Sequence Lengths (3')",
        xaxis="Length trimmed (bp)", yaxis="Count",
        dataset_idx=0)

    # ── 1.2 Trimmed QC ────────────────────────────────────────────────────────
    # Keys confirmed from real trimming MultiQC HTML (fastqc_* prefix, no fqc_trimmed_):
    #   fastqc_sequence_counts_plot              bar, orient=h, cats=[Unique, Duplicate]
    #   fastqc_per_base_sequence_quality_plot    line, datasets=[1]
    #   fastqc_per_sequence_quality_scores_plot  line, datasets=[1]
    #   fastqc_per_sequence_gc_content_plot      line, datasets=[Percentages, Counts]
    #   fastqc_per_base_n_content_plot           line, datasets=[1]
    #   fastqc_sequence_length_distribution_plot line, datasets=[1]
    #   fastqc_sequence_duplication_levels_plot  line, datasets=[1]
    #   fastqc_adapter_content_plot              line, datasets=[1]
    #   fastqc-status-check-heatmap              heatmap, xcats=modules, ycats=samples

    fqc_counts = mqc_bar_chart(
        "fastqc_sequence_counts_plot", trim_data,
        title="Sequence Counts", yaxis="Reads (M)")

    fqc_qual = mqc_line_chart(
        "fastqc_per_base_sequence_quality_plot", trim_data,
        title="Per Base Sequence Quality",
        xaxis="Position (bp)", yaxis="Phred Score")

    fqc_seq_qual = mqc_line_chart(
        "fastqc_per_sequence_quality_scores_plot", trim_data,
        title="Per Sequence Quality Scores",
        xaxis="Phred Score", yaxis="Count")

    # GC content: dataset[0]=Percentages, dataset[1]=Counts
    fqc_gc = mqc_line_chart(
        "fastqc_per_sequence_gc_content_plot", trim_data,
        title="Per Sequence GC Content",
        xaxis="GC (%)", yaxis="Percentage of reads",
        dataset_idx=0, pct_yaxis=True)

    fqc_n = mqc_line_chart(
        "fastqc_per_base_n_content_plot", trim_data,
        title="Per Base N Content",
        xaxis="Position (bp)", yaxis="N%")

    fqc_len = mqc_line_chart(
        "fastqc_sequence_length_distribution_plot", trim_data,
        title="Sequence Length Distribution",
        xaxis="Length (bp)", yaxis="Count")

    fqc_dup = mqc_line_chart(
        "fastqc_sequence_duplication_levels_plot", trim_data,
        title="Sequence Duplication Levels",
        xaxis="Duplication level", yaxis="% of reads")

    fqc_adapter = mqc_line_chart(
        "fastqc_adapter_content_plot", trim_data,
        title="Adapter Content",
        xaxis="Position (bp)", yaxis="% Adapter",
        pct_yaxis=True)

    fqc_status = mqc_heatmap(
        "fastqc-status-check-heatmap", trim_data,
        title="Sequencing Status Checks",
        exclude_modules=["Per Base Sequence Content"])

    # ── 1.3 Alignment ─────────────────────────────────────────────────────────
    # bwameth → samtools flagstat/stats → MultiQC samtools module
    # Possible plot keys depending on MultiQC version: samtools-flagstat-dp,
    # samtools_alignment_plot, etc.  Try common keys.
    # Confirmed from real MultiQC HTML: key is 'samtools_alignment_plot' (horizontal bar)
    # samtools-flagstat-dp and samtools-stats-dp are violin/table format, not bar charts
    align_key = next(
        (resolved for candidate in ["samtools_alignment_plot",
                                    "samtools-flagstat-dp",
                                    "bwameth_alignment"]
         for resolved in [_mqc_key(align_data, candidate)]
         if resolved in align_data),
        None,
    )
    if align_key:
        align_chart = mqc_bar_chart(
            align_key, align_data,
            title="Alignment Rates (bwameth)", yaxis="Reads (M)")
    else:
        align_chart = _missing(
            "Alignment plot not found in MultiQC report. "
            "bwameth alignment statistics are captured via "
            "<code>samtools flagstat</code> and <code>samtools stats</code>. "
            "Ensure MultiQC is run after <code>cftk process -s 2</code>."
        )

    # Deduplication — MultiQC has no sambamba markdup module.
    # Parse .markdup_metrics.txt from 3_markdup/ + flagstat from 2_alignment/.
    markdup_dirs = [os.path.join(root, "1_process", "3_markdup") for root in result_roots]
    flagstat_dirs = [os.path.join(root, "1_process", "2_alignment") for root in result_roots]
    dup_chart    = markdup_dedup_chart(
        markdup_dirs,
        samples=all_sample_names,
        flagstat_dir=flagstat_dirs,
        title="Deduplication (sambamba markdup)",
    )

    note_bwameth = (
        f'<div style="margin-bottom:12px;padding:10px 14px;'
        f'background:rgba(30,95,160,.06);border-left:3px solid var(--accent);'
        f'border-radius:0 6px 6px 0;font-size:14.5px;color:#000000;">'
        f'<strong>Note:</strong> bwameth does not produce Bismark-style strand '
        f'alignment or cytosine methylation context plots. Alignment rates are '
        f'derived from <code>samtools flagstat</code>; deduplication from '
        f'<code>sambamba markdup</code>.</div>'
    )

    # ── 1.4 M-bias ────────────────────────────────────────────────────────────
    # Priority order:
    #   1. MultiQC methyldackel/bismark module (if supported by MultiQC version)
    #   2. Per-position TSV files saved by qc_parser.parse_mbias_txt()
    #      (populated when MethylDackel --txt is supported)
    #   3. Legacy fallback: show OT/OB coords table from stderr-only _mbias.txt

    meth_dir = os.path.join(rd, "1_process", "4_methylation")
    mbias_tsv_dirs = [
        os.path.join(root, "1_process", "4_methylation", "mbias_data")
        for root in result_roots
    ]
    if not any(os.path.isdir(path) for path in mbias_tsv_dirs):
        mbias_tsv_dirs = [os.path.join(root, "2_qc", "mbias_data") for root in result_roots]
    mbias_tsv_paths = [
        os.path.join(source_results[sample], "1_process", "4_methylation", "mbias_data",
                     f"{sample}_mbias.tsv")
        for sample in all_sample_names
        if source_results and sample in source_results
    ]

    mbias_key = next(
        (resolved for candidate in ["methyldackel_mbias", "bismark_mbias", "mbias"]
         for resolved in [_mqc_key(meth_data, candidate)]
         if resolved in meth_data),
        None,
    )

    def _has_tsv_mbias(path):
        """
        Return True if _mbias.txt contains per-position TSV data from --txt flag.
        MethylDackel --txt outputs lines like:
            Strand  Read  Position  nMethylated  nUnmethylated
            OT      1     1         18243         18091
        Check for: "Strand" header OR numeric data lines (OT/OB + digits).
        Also handles files where first line(s) may be blank or comment.
        """
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                for _ in range(10):   # check first 10 non-empty lines
                    line = f.readline()
                    if not line:
                        break
                    s = line.strip()
                    if not s:
                        continue
                    # Header line
                    if s.startswith("Strand") or s.startswith("strand"):
                        return True
                    # Data line: starts with OT / OB / CTOT / CTOB followed by digits
                    import re as _re2
                    if _re2.match(r"^(OT|OB|CTOT|CTOB)\s+\d", s):
                        return True
                    # If first content line is stderr ("Suggested" / "["), not TSV
                    if s.startswith("Suggested") or s.startswith("["):
                        return False
        except OSError:
            pass
        return False

    def _mbias_txt_paths():
        """Yield (sample_name, path) for all readable _mbias.txt files."""
        for s in all_sample_names:
            p = os.path.join(meth_dir, f"{s}_mbias.txt")
            if os.path.exists(p) and os.path.getsize(p) > 0:
                yield s, p

    if mbias_key:
        # Path 1: MultiQC has parsed the mbias data
        mbias_chart = mqc_mbias_chart(mbias_key, meth_data,
                                       title="M-bias (MethylDackel)")

    elif any(glob.glob(os.path.join(path, "*_mbias.tsv")) for path in mbias_tsv_dirs) or any(
        os.path.exists(path) for path in mbias_tsv_paths
    ):
        # Path 2a: qc_parser already saved intermediate TSV files
        mbias_chart = mbias_tsv_chart(
            mbias_tsv_dirs, title="M-bias (MethylDackel)", tsv_paths=mbias_tsv_paths
        )

    elif any(_has_tsv_mbias(p) for _, p in _mbias_txt_paths()):
        # Path 2b: _mbias.txt files have --txt TSV data.
        # Parse directly → save to mbias_tsv_dir → render.
        # Success = TSV file written with at least one data row (not NaN check).
        mbias_tsv_dir = mbias_tsv_dirs[0]
        os.makedirs(mbias_tsv_dir, exist_ok=True)
        import sys as _sys
        _sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        tsv_written = 0
        try:
            from analysis.qc_parser import parse_mbias_txt
            for s, p in _mbias_txt_paths():
                if not _has_tsv_mbias(p):
                    continue
                parse_mbias_txt(p, save_tsv_dir=mbias_tsv_dir)
                tsv_out = os.path.join(mbias_tsv_dir, f"{s}_mbias.tsv")
                if os.path.exists(tsv_out) and os.path.getsize(tsv_out) > 50:
                    tsv_written += 1
        except Exception as _e:
            print(f"[report] mbias parse error: {_e}")
        if tsv_written > 0:
            mbias_chart = mbias_tsv_chart(mbias_tsv_dir, title="M-bias (MethylDackel)")
        else:
            mbias_chart = mbias_legacy_chart(meth_dir, samples=all_sample_names)

    else:
        # Path 3: legacy format — only OT/OB coords in stderr, no TSV data
        mbias_chart = mbias_legacy_chart(meth_dir, samples=all_sample_names)

    # ── 1.6 Sequencing QC Summary ─────────────────────────────────────────────────
    qc_table_html = _qc_table(qc_scores)

    return f"""
    <section class="section" id="part1">
      <div class="section-header">
        <span class="section-num">PART 01</span>
        <h2 class="section-title">Data Processing</h2>
        <span class="section-tag">process</span>
      </div>

      <h3 class="subsec-title" id="part1_1">1.1 Trimming</h3>
      <p style="font-size:15.5px;color:#000000;margin-bottom:14px;">
        Adapter trimming quality metrics from Trim Galore / Cutadapt.
      </p>
      {_coll("Filtered Reads", filtered_chart, open_=True)}
      {_coll("Trimmed Sequence Lengths (3')",
             _note(
        "Hover over any line to see the sample name.",
        kind="info"
    ) + lengths_chart,
             open_=True)}

      <h3 class="subsec-title" id="part1_2">1.2 Trimmed QC</h3>
      <p style="font-size:15.5px;color:#000000;margin-bottom:14px;">
        FastQC quality metrics on trimmed reads.
      </p>
      {_coll("Sequence Counts",              fqc_counts,   open_=True)}
      {_coll("Sequence Quality Histograms",  fqc_qual,     open_=False)}
      {_coll("Per Sequence Quality Scores",  fqc_seq_qual, open_=False)}
      {_coll("Per Sequence GC Content",      fqc_gc,       open_=False)}
      {_coll("Per Base N Content",           fqc_n,        open_=False)}
      {_coll("Sequence Length Distribution", fqc_len,      open_=False)}
      {_coll("Sequence Duplication Levels",  fqc_dup,      open_=False)}
      {_coll("Adapter Content",              fqc_adapter,  open_=False)}
      {_coll("Sequencing Status Checks", fqc_status, open_=True)}

      <h3 class="subsec-title" id="part1_3">1.3 Alignment</h3>
      <p style="font-size:15.5px;color:#000000;margin-bottom:14px;">
        Bisulfite sequencing alignment statistics.
      </p>
      {_coll("Alignment &amp; Deduplication",
             note_bwameth + align_chart + '<div style="margin-top:16px;">' + dup_chart + '</div>',
             open_=True)}

      <h3 class="subsec-title" id="part1_4">1.4 M-bias</h3>
      <p style="font-size:15.5px;color:#000000;margin-bottom:14px;">
        Per-position methylation bias from MethylDackel mbias.
        Dataset buttons switch between OT/OB strand and R1/R2.
        All samples shown simultaneously.
      </p>
      {_coll("M-bias Plot", mbias_chart, open_=True)}

      <h3 class="subsec-title" id="part1_6">1.5 Sequencing QC Summary</h3>
      <p style="font-size:15.5px;color:#000000;margin-bottom:14px;">
        Per-sample quality metrics. Click any row to view the recommendation.
        Column headers are sortable.
      </p>
      {_coll("Sequencing QC Metrics Table", qc_table_html, open_=True)}

    </section>"""


def _plotly_select_render(key: str, label_text: str, menu_field: str,
                          height: int = 460, square: bool = False) -> str:
    """
    Render a Plotly chart from window.__CFTK_DATA__[key] with a DELFI-style HTML
    <select> control above it (instead of a native Plotly dropdown). The fig must
    carry `menu_field` = {"options": [{"label", "visible":[...], "relayout":{...}}]}.
    Selecting an option restyles trace visibility and applies any relayout updates.
    """
    div_id = f"plotly_{key}"
    sel_id = f"sel_{key}"
    tag_end = "</scr" + "ipt>"
    container = "max-width:680px;margin:0 auto;" if square else "width:100%;"
    chart_h = "aspect-ratio:1/1;" if square else f"min-height:{height}px;"
    controls = (
        '<div style="display:flex;align-items:center;gap:10px;margin-bottom:10px;">'
        f'<label style="font-size:14.5px;font-weight:500;color:#000000;'
        f'font-family:var(--mono);">{label_text}</label>'
        f'<select id="{sel_id}" style="font-size:14.5px;padding:4px 10px;'
        f'border:1px solid var(--rule);border-radius:4px;background:white;'
        f'color:#000000;cursor:pointer;"></select>'
        '</div>'
    )
    chart = (f'<div style="{container}"><div id="{div_id}" '
             f'style="width:100%;{chart_h}"></div></div>')
    js = (
        '<script>(function(){'
        f'var fig=(window.__CFTK_DATA__||{{}})["{key}"];'
        f'var el=document.getElementById("{div_id}");var sel=document.getElementById("{sel_id}");'
        f'if(!fig){{el.innerHTML="<p style=\'color:#999;padding:20px;\'>Chart data not available for key: {key}</p>";if(sel)sel.style.display="none";return;}}'
        'if(fig._missing){el.innerHTML="<div style=\'padding:14px 16px;background:var(--bg2);border:1px solid var(--rule);border-radius:8px;font-size:15.5px;color:#000000;\'>\u26a0 "+(fig.layout&&fig.layout.title&&fig.layout.title.text||"Data not available")+"</div>";if(sel)sel.style.display="none";return;}'
        f'var menu=fig["{menu_field}"];var opts=(menu&&menu.options)||[];'
        'if(!opts.length){if(sel)sel.style.display="none";}'
        'opts.forEach(function(o,i){var op=document.createElement("option");op.value=i;op.textContent=o.label;sel.appendChild(op);});'
        'var layout=Object.assign({},fig.layout);'
        + ('layout.autosize=true;layout.height=undefined;layout.width=undefined;' if square else '')
        + f'setTimeout(function(){{'
        + f'Plotly.newPlot("{div_id}",fig.data,layout,{{responsive:true,displayModeBar:true,modeBarButtonsToRemove:["lasso2d","select2d"]}});'
        + 'if(sel&&opts.length){sel.addEventListener("change",function(){'
        + 'var o=opts[this.value];if(!o)return;'
        + 'var idx=o.visible.map(function(_,i){return i;});'
        + f'Plotly.restyle("{div_id}",{{visible:o.visible}},idx);'
        + f'if(o.relayout)Plotly.relayout("{div_id}",o.relayout);'
        + '});}'
        + f'}},0);'
        + '})();'
        + tag_end
    )
    return controls + chart + js


def _plotly_from_data(key: str, height: int = 400, square: bool = False) -> str:
    """Render a Plotly chart from window.__CFTK_DATA__[key]."""
    div_id = f"plotly_{key}"
    container_style = (
        "max-width:640px;aspect-ratio:1/1;margin:0 auto;"
        if square else "width:100%;"
    )
    return (
        f'<div style="{container_style}">'
        f'<div id="{div_id}" style="width:100%;height:100%;'
        + (f'min-height:{height}px;' if not square else '')
        + f'"></div></div>'
        f'<script>'
        f'(function(){{'
        f'var fig=(window.__CFTK_DATA__||{{}})["{key}"];'
        f'  var el=document.getElementById("{div_id}");'
        f'  if(!fig){{el.innerHTML="<p style=\'color:#999;padding:20px;\'>Chart data not available for key: {key}</p>";return;}}'
        f'  if(fig._missing){{'
        f'    el.parentElement.style.aspectRatio="auto";'
        f'    el.innerHTML="<div style=\'padding:14px 16px;background:var(--bg2);border:1px solid var(--rule);border-radius:8px;font-size:15.5px;color:#000000;\'>⚠ "+(fig.layout&&fig.layout.title&&fig.layout.title.text||"Data not available")+"</div>";'
        f'    return;}}'
        f'  if({str(square).lower()}){{'
        f'    fig.layout=fig.layout||{{}};fig.layout.autosize=true;'
        f'    fig.layout.height=undefined;fig.layout.width=undefined;}}'
        f'  Plotly.newPlot("{div_id}",fig.data,fig.layout,'
        f'  {{responsive:true,displayModeBar:true,modeBarButtonsToRemove:["lasso2d","select2d"]}});'
        f'}})();'
        f'</script>'
    )


def _frag_render() -> str:
    """Fragment length chart — fills its parent box. Per-sample traces are kept
    but hidden (and removed from the legend); clicking a sample row in the side
    table toggles that sample's line on the chart via window.fragToggleSample."""
    div_id = "plotly_fragment_length"
    tag_end = "</scr" + "ipt>"
    return (
        '<div style="width:100%;height:100%;">'
        f'<div id="{div_id}" style="width:100%;height:100%;min-height:380px;"></div>'
        '</div>'
        '<script>(function(){'
        'var fig=(window.__CFTK_DATA__||[])["fragment_length"];'
        + f'var el=document.getElementById("{div_id}");'
        + 'if(!fig){el.innerHTML="<p style=\'color:#999;padding:20px;\'>Fragment length data not available</p>";return;}'
        + 'if(fig._missing){el.innerHTML="<div style=\'padding:14px;background:var(--bg2);border:1px solid var(--rule);border-radius:8px;font-size:15.5px;color:#000000;\'>\u26a0 "+(fig.layout&&fig.layout.title&&fig.layout.title.text||"Data not available")+"</div>";return;}'
        # copy traces; hide per-sample traces + drop them from the legend; map sample name -> trace index
        + 'var data=(fig.data||[]).map(function(t){return Object.assign({},t);});'
        + 'var sidx={};'
        + 'data.forEach(function(t,i){if(t.legendgroup==="samples"){t.visible=false;t.showlegend=false;var nm=(t.name||"").split(" (")[0];sidx[nm]=i;}});'
        + 'var layout=Object.assign({},fig.layout);delete layout._trace_peaks;'
        + 'layout.autosize=true;layout.height=undefined;layout.width=undefined;'
        + f'setTimeout(function(){{'
        + f'Plotly.newPlot("{div_id}",data,layout,{{responsive:true,displayModeBar:true,modeBarButtonsToRemove:["lasso2d","select2d"]}});'
        + f'window.__fragGraph=document.getElementById("{div_id}");window.__fragSampleIdx=sidx;'
        + 'window.fragToggleSample=function(name,rowEl){'
        + 'var gd=window.__fragGraph;var idx=window.__fragSampleIdx?window.__fragSampleIdx[name]:null;'
        + 'if(!gd||idx==null)return;'
        + 'var on=!(gd.data[idx].visible===true);'
        + 'Plotly.restyle(gd,{visible:on?true:false},[idx]);'
        + 'if(rowEl){rowEl.setAttribute("data-active",on?"1":"0");rowEl.style.background=on?"rgba(11,110,153,.10)":"";}'
        + '};'
        + f'}},0);'
        + '})();'
        + tag_end
    )



def _beta_density_table(scores: str | pd.DataFrame | None, table_height: int = 420):
    """
    Render a compact per-sample table of the β M-score QC metric.
    Returns a (note_html, table_html) tuple so callers can place the
    explanatory note and the table independently (e.g. side-by-side layout).
    Read from a QC-score table or dataframe.
    """
    df = _as_qc_dataframe(scores)
    if df is None or df.empty:
        return ("", "")

    cols_needed = ["beta_M_score"]
    val_cols  = [f"{c}_value"  for c in cols_needed]
    stat_cols = [f"{c}_status" for c in cols_needed]

    # Check at least one beta column present
    any_present = any(v in df.columns or c in df.columns
                      for c, v in zip(cols_needed, val_cols))
    if not any_present:
        return ("", _missing(
            "Beta-density QC metrics not yet computed — "
            "run <code>cftk qc -s 0 --force</code> to regenerate QC scores."
        ))

    labels = {
        "beta_M_score":            ("β M-score",            ".2f", "",
                                    "min(left_peak, right_peak) / mid_density. PASS ≥ 2"),
    }

    th_s = (
        "padding:8px 12px;background:var(--header-bg);color:var(--header-ink);"
        "font-size:14px;font-weight:700;letter-spacing:.01em;white-space:nowrap;"
        "border-bottom:2px solid var(--rule);position:sticky;top:0;z-index:2;"
    )

    # Header
    ths = f'<th style="{th_s}text-align:center;">Sample</th>'
    active_cols = []
    for col, (label, fmt, unit, note) in labels.items():
        val_col  = f"{col}_value"
        stat_col = f"{col}_status"
        # Accept both with and without _value suffix
        actual_val  = val_col  if val_col  in df.columns else (col if col  in df.columns else None)
        actual_stat = stat_col if stat_col in df.columns else None
        if actual_val is None:
            continue
        active_cols.append((col, label, fmt, unit, actual_val, actual_stat))
        ths += (
            f'<th style="{th_s}text-align:center;" title="{note}">' + label + '</th>'
        )

    if not active_cols:
        return ("", _missing("Beta-density QC columns not found in qc_scores.tsv."))

    tbody = ""
    for _, row in df.iterrows():
        sample = str(row.get("sample", ""))
        group  = str(row.get("group", ""))
        td_s = (
            f'<td style="padding:6px 12px;font-family:var(--mono);font-size:14.5px;' +
            f'white-space:nowrap;text-align:center;" data-sort="{sample}">{sample}</td>'
        )
        metric_tds = ""
        for col, label, fmt, unit, actual_val, actual_stat in active_cols:
            stat = str(row.get(actual_stat, "NA")).upper() if actual_stat else "NA"
            # Per request: show FAIL as WARN (yellow) for beta-density metrics
            if stat == "FAIL":
                stat = "WARN"
            # Badge-only cell (no numeric value)
            metric_tds += (
                f'<td style="text-align:center;padding:6px 8px;'
                f'background:{_STATUS_BG.get(stat, _STATUS_BG["NA"])};">'
                + _badge(stat) + '</td>'
            )
        tbody += (
            f'<tr style="border-top:1px solid var(--rule);' +
            f'onmouseover="this.style.background=\'var(--bg2)\'"' +
            f'onmouseout="this.style.background=\'\'">'
            + td_s + metric_tds + '</tr>'
        )

    # Single note with three-point explanation of the metrics
    metric_note = _note(
        "Human cfDNA shows a bimodal distribution (peaks near 0 and 1). "
        "The <strong>β M-score</strong> assesses whether the CpG methylation "
        "distribution is healthily bimodal. "
        "β M-score = min(left_peak, right_peak) / mid_density "
        "(left_peak: max density at β≤0.15; right_peak: max density at β≥0.85; "
        "mid_density: median density at 0.35≤β≤0.65). "
        "<strong>PASS ≥ 2</strong> — a strong bimodal shape relative to the mid region.",
        kind="info"
    )

    note_html = metric_note
    table_html = (
        f'<div style="height:{table_height}px;overflow:auto;'
        'border:1px solid var(--rule);border-radius:8px;">'
        '<table style="width:100%;border-collapse:collapse;font-size:15.5px;">'
        f'<thead><tr>{ths}</tr></thead><tbody>{tbody}</tbody></table>'
        '</div>'
    )
    return (note_html, table_html)


def _fragment_peak_table(rd, groups, table_height: int = 420, source_results=None):
    """
    Per-sample fragment-length peak table (Sample | Peak site), styled and
    behaving exactly like the β M-score table (sticky centered header, hover
    rows, fixed-height scroll). Returns table HTML, or a _missing message.
    Peak sites come from _fragment_peaks (2_qc/2_fragment_length/*.raw.csv).
    """
    samples = [s for m in groups.values() for s in m]
    peaks = _fragment_peaks(rd, samples, source_results=source_results)
    if not peaks:
        return _missing(
            "Fragment-length peaks not found — run <code>cftk qc -s 2</code> "
            "to populate per-sample fragment length data."
        )

    th_s = (
        "padding:8px 12px;background:var(--header-bg);color:var(--header-ink);"
        "font-size:14px;font-weight:700;letter-spacing:.01em;white-space:nowrap;"
        "border-bottom:2px solid var(--rule);position:sticky;top:0;z-index:2;"
    )
    ths = (
        f'<th style="{th_s}text-align:center;">Sample</th>'
        f'<th style="{th_s}text-align:center;">Peak site</th>'
    )

    tbody = ""
    for g, members in groups.items():
        for s in members:
            p = peaks.get(s)
            val = f"{p} bp" if p else "—"
            td_s = (
                f'<td style="padding:6px 12px;font-family:var(--mono);font-size:14.5px;'
                f'white-space:nowrap;text-align:center;" data-sort="{s}">{s}</td>'
            )
            td_p = (
                f'<td style="text-align:center;padding:6px 8px;font-size:15.5px;'
                f'color:var(--ink);">{val}</td>'
            )
            # clickable row: toggles this sample's line on the chart; active rows
            # keep an accent tint that survives hover in/out.
            tbody += (
                f'<tr data-active="0" style="border-top:1px solid var(--rule);cursor:pointer;"'
                f' onclick="window.fragToggleSample&&window.fragToggleSample(\'{s}\',this)"'
                f' onmouseover="if(this.getAttribute(\'data-active\')!==\'1\')this.style.background=\'var(--bg2)\'"'
                f' onmouseout="if(this.getAttribute(\'data-active\')!==\'1\')this.style.background=\'\'">'
                f'{td_s}{td_p}</tr>'
            )

    return (
        f'<div style="height:{table_height}px;overflow:auto;'
        'border:1px solid var(--rule);border-radius:8px;">'
        '<table style="width:100%;border-collapse:collapse;font-size:15.5px;">'
        f'<thead><tr>{ths}</tr></thead><tbody>{tbody}</tbody></table>'
        '</div>'
    )


def _sec_qc(
    rd: str,
    groups: dict | None = None,
    source_results: dict[str, str] | None = None,
    qc_scores: pd.DataFrame | None = None,
) -> str:
    groups = groups or {}

    # 2.1 Methylation Distribution. Current cohort-level output wins; otherwise
    # retain one selectable source image per configured BAM-backed sample.
    local_meth = _find(rd, "2_qc", "1_methylation_distribution", pattern="*.png")
    source_meth = []
    if not local_meth:
        for sample in _sample_names(groups):
            root = (source_results or {}).get(sample)
            image = (_find(root, "2_qc", "1_methylation_distribution", pattern="*.png")
                     if root else None)
            if image:
                source_meth.append((sample, image))

    meth_img = _img_tag(
        local_meth or (source_meth[0][1] if len(source_meth) == 1 else None),
        "Methylation distribution",
        style="max-height:100%;max-width:100%;width:auto;object-fit:contain;display:block;margin:auto;")
    meth_gallery = (
        _dropdown_gallery("Sample", source_meth, uid="methylation_distribution")
        if len(source_meth) > 1 else ""
    )

    # 2.2 Fragment Length — interactive
    frag_chart = _frag_render()

    # 2.3 Dinucleotide — interactive
    dinuc_chart = _plotly_from_data("dinucleotide", height=420)

    # 2.4 PCA — interactive
    pca_desc = (
        '<p style="font-size:14.5px;color:#000000;margin-bottom:8px;">'
        'Select modality from dropdown. Hover points for sample names.</p>'
    )
    pca_block = _coll("PCA plot",
                      pca_desc + _plotly_select_render("pca", "Modality:",
                                                       "_modality_menu",
                                                       height=460, square=True),
                      open_=True)

    # Pre-compute collapsible blocks to avoid f-string nesting issues
    _BETA_H = 420
    beta_note, beta_table = _beta_density_table(qc_scores, table_height=_BETA_H)

    meth_body = (meth_gallery if meth_gallery else
        f'<div class="fig-card" style="height:{_BETA_H}px;display:flex;'
        f'flex-direction:column;margin:0;">'
        f'<div style="flex:1;min-height:0;display:flex;align-items:center;'
        f'justify-content:center;padding:14px;">{meth_img}</div>'
        f'</div>')
    meth_left = f'<div style="flex:1 1 58%;min-width:320px;">{meth_body}</div>'
    meth_right = (
        f'<div style="flex:1 1 38%;min-width:280px;">{beta_table}</div>'
        if beta_table else ""
    )
    meth_row = (
        f'<div style="display:flex;gap:18px;align-items:stretch;flex-wrap:wrap;">'
        f'{meth_left}{meth_right}</div>'
    )
    meth_block = _coll("β-value density plot", meth_row + beta_note, open_=True)

    frag_table = _fragment_peak_table(
        rd, groups, table_height=_BETA_H, source_results=source_results
    )
    frag_left = (
        f'<div style="flex:1 1 58%;min-width:320px;">'
        f'<div class="fig-card" style="height:{_BETA_H}px;display:flex;'
        f'flex-direction:column;margin:0;">'
        f'<div style="flex:1;min-height:0;">{frag_chart}</div>'
        f'</div></div>'
    )
    frag_right = (
        f'<div style="flex:1 1 38%;min-width:280px;">{frag_table}</div>'
        if frag_table else ""
    )
    frag_row = (
        f'<div style="display:flex;gap:18px;align-items:stretch;flex-wrap:wrap;">'
        f'{frag_left}{frag_right}</div>'
    )
    frag_desc = (
        '<p style="font-size:14.5px;color:#000000;margin-top:14px;margin-bottom:8px;">'
        'Group mean traces are shown by default (x-axis: 50–250 bp). '
        '<strong>Click a sample in the table</strong> to show its fragment-length line '
        'on the chart (click again to hide); active samples are highlighted. '
        'Hover the chart for exact values.</p>'
    )
    frag_block = _coll("Fragment length plot", frag_row + frag_desc, open_=True)

    dinuc_desc = (
        '<p style="font-size:14.5px;color:#000000;margin-bottom:8px;">'
        '10-bp periodicity of AT- and GC-rich dinucleotides around fragment centres '
        'reflects nucleosome positioning.</p>'
    )
    dinuc_block = _coll("Dinucleotide frequency plot", dinuc_desc + dinuc_chart, open_=True)

    return f"""
    <section class="section" id="part2">
      <div class="section-header">
        <span class="section-num">PART 02</span>
        <h2 class="section-title">cfDNA QC Analysis</h2>
        <span class="section-tag">qc</span>
      </div>

      <h3 class="subsec-title" id="part2_1">2.1 Methylation Distribution</h3>
      {meth_block}

      <h3 class="subsec-title" id="part2_2">2.2 Fragment Length Distribution</h3>
      {frag_block}

      <h3 class="subsec-title" id="part2_3">2.3 Dinucleotide Frequency</h3>
      {dinuc_block}

      <h3 class="subsec-title" id="part2_4">2.4 PCA</h3>
      {pca_block}
    </section>"""

def _dmr_sortable_table() -> str:
    """
    Render a client-side sortable table of top DMRs from
    window.__CFTK_DATA__["dmr_table"]. Columns sortable by clicking headers.
    The table body scrolls vertically inside a fixed ~10-row window with a
    sticky header; only gene-annotated, de-duplicated DMRs are shown.
    """
    div_id = "dmr_table_mount"
    return (
        f'<div id="{div_id}" style="margin-top:14px;"></div>'
        '<script>(function(){'
        'var d=(window.__CFTK_DATA__||{})["dmr_table"];'
        f'var el=document.getElementById("{div_id}");'
        'if(!d){el.innerHTML="<p style=\'color:#999;font-size:15.5px;\'>DMR table data not available</p>";return;}'
        'if(d._missing){el.innerHTML="<div style=\'padding:12px;background:var(--bg2);border:1px solid var(--rule);border-radius:8px;font-size:15.5px;color:#000000;\'>\u26a0 "+(d._msg||"No significant DMRs")+"</div>";return;}'
        'var cols=d.columns,rows=d.rows.slice();'
        'var sortState={col:5,asc:true};'  # default sort by q-value asc
        'function render(){'
        'var h="<div style=\'font-size:14.5px;color:#000000;margin-bottom:8px;\'>Showing top "+d.shown+" of "+d.n_total_sig+" significant DMRs (q < "+d.q_thr+"). Click a column header to sort.</div>";'
        'h+="<div style=\'max-height:362px;overflow:auto;border:1px solid var(--rule);border-radius:8px;\' id=\'dmr_scroll\'>";'
        'h+="<table style=\'width:100%;border-collapse:collapse;font-size:14.5px;min-width:760px;\'>";'
        'h+="<thead><tr>";'
        'cols.forEach(function(c,i){'
        'var arrow=sortState.col===i?(sortState.asc?" \u25b2":" \u25bc"):"";'
        'h+="<th data-ci=\'"+i+"\' style=\'position:sticky;top:0;z-index:2;padding:8px 10px;background:var(--header-bg);color:var(--header-ink);font-family:var(--sans);font-size:13px;font-weight:700;letter-spacing:.02em;white-space:nowrap;text-align:center;border-bottom:2px solid var(--rule);cursor:pointer;user-select:none;\'>"+c+arrow+"</th>";'
        '});'
        'h+="</tr></thead><tbody>";'
        'rows.forEach(function(r){'
        'h+="<tr style=\'border-top:1px solid var(--rule);\'>";'
        'r.forEach(function(v,ci){'
        'var align="center";'
        'var extra="";'
        'if(ci===8){var col=v==="Hyper"?"#a93226":"#1e5fa0";extra="color:"+col+";font-weight:600;";}'
        'if(ci===0)extra+="font-family:var(--mono);";'
        'h+="<td style=\'padding:5px 10px;text-align:"+align+";white-space:nowrap;"+extra+"\'>"+(v===""?"\u2014":v)+"</td>";'
        '});'
        'h+="</tr>";'
        '});'
        'h+="</tbody></table></div>";'
        'el.innerHTML=h;'
        'cols.forEach(function(c,i){'
        'el.querySelector("th[data-ci=\'"+i+"\']").addEventListener("click",function(){'
        'if(sortState.col===i){sortState.asc=!sortState.asc;}else{sortState.col=i;sortState.asc=true;}'
        'rows.sort(function(a,b){'
        'var x=a[i],y=b[i];'
        'var nx=parseFloat(String(x).replace(/[^0-9.eE+-]/g,"")),ny=parseFloat(String(y).replace(/[^0-9.eE+-]/g,""));'
        'var bothNum=!isNaN(nx)&&!isNaN(ny)&&String(x).trim()!==""&&String(y).trim()!=="";'
        'var cmp=bothNum?(nx-ny):String(x).localeCompare(String(y));'
        'return sortState.asc?cmp:-cmp;'
        '});'
        'render();'
        '});'
        '});'
        '}'
        'render();'
        '})();</scr'+'ipt>'
    )


def _sec_software() -> str:
    """Render the software/tools list section at the end of the report."""
    json_path = os.path.join(os.path.dirname(__file__), "software_list.json")
    try:
        with open(json_path) as f:
            data = json.load(f)
    except Exception as e:
        return (
            '<section class="section" id="part6">'
            '<div class="section-header"><span class="section-num">PART 06</span>'
            '<h2 class="section-title">Software</h2></div>'
            f'<p><em>Could not load software list: {e}</em></p></section>'
        )

    blocks = ""
    for cat in data.get("categories", []):
        stage = cat.get("stage", "")
        rows = ""
        for t in cat.get("tools", []):
            name = t.get("name", "")
            ver  = t.get("version", "")
            purp = t.get("purpose", "")
            url  = t.get("url", "")
            name_html = (f'<a href="{url}" target="_blank" '
                         f'style="color:var(--accent,#1e5fa0);text-decoration:none;">{name}</a>'
                         if url else name)
            rows += (
                '<tr style="border-top:1px solid var(--rule);">'
                f'<td style="padding:6px 12px;font-weight:600;white-space:nowrap;">{name_html}</td>'
                f'<td style="padding:6px 12px;font-family:var(--mono);font-size:14.5px;'
                f'color:#000000;white-space:nowrap;">{ver}</td>'
                f'<td style="padding:6px 12px;font-size:14.5px;color:#000000;">{purp}</td>'
                '</tr>'
            )
        blocks += (
            f'<h3 class="subsec-title">{stage}</h3>'
            '<div style="overflow-x:auto;border:1px solid var(--rule);'
            'border-radius:8px;margin-bottom:18px;">'
            '<table style="width:100%;border-collapse:collapse;font-size:15.5px;">'
            '<thead><tr>'
            '<th style="padding:7px 12px;background:var(--header-bg);color:var(--header-ink);'
            'font-size:13px;font-weight:700;letter-spacing:.02em;'
            'text-align:left;border-bottom:2px solid var(--rule);">Tool</th>'
            '<th style="padding:7px 12px;background:var(--header-bg);color:var(--header-ink);'
            'font-size:13px;font-weight:700;letter-spacing:.02em;'
            'text-align:left;border-bottom:2px solid var(--rule);">Version</th>'
            '<th style="padding:7px 12px;background:var(--header-bg);color:var(--header-ink);'
            'font-size:13px;font-weight:700;letter-spacing:.02em;'
            'text-align:left;border-bottom:2px solid var(--rule);">Purpose</th>'
            f'</tr></thead><tbody>{rows}</tbody></table></div>'
        )

    return (
        '<section class="section" id="part6">'
        '<div class="section-header">'
        '<span class="section-num">PART 06</span>'
        '<h2 class="section-title">Software &amp; Tools</h2>'
        '<span class="section-tag">env</span>'
        '</div>'
        '<p style="font-size:14.5px;color:#000000;margin-bottom:16px;">'
        'Bioinformatics tools used in this cfDNA methylation analysis pipeline, '
        'grouped by processing stage.</p>'
        f'{blocks}</section>'
    )


def _sec_differential(rd):

    diff_base  = os.path.join(rd, "3_differential")
    modalities = sorted([
        d for d in os.listdir(diff_base)
        if os.path.isdir(os.path.join(diff_base, d)) and d != "dmr"
    ]) if os.path.exists(diff_base) else []

    def _dd(stem, uid_prefix):
        items = []
        for mod in modalities:
            p = _find(rd, "3_differential", mod, pattern=f"{stem}.png")
            if p:
                items.append((mod.upper(), p))
        return _dropdown_gallery("Modality:", items, uid=f"{uid_prefix}_{stem}")

    # Pre-compute all blocks outside f-string to avoid double-brace escaping
    violin_block = _coll("Violin plots (static)", _dd("violin", "diff"), open_=True)
    heatmap_block = _coll("DMC heatmaps (static)", _dd("heatmap", "diff"), open_=True)
    dmr_chart = _plotly_from_data("dmr_volcano", height=460, square=True)
    dmr_desc = (
        '<p style="font-size:14.5px;color:#000000;margin-bottom:8px;">'
        'Each point is one DMR. Hover for gene name and exact values. '
        'Dashed line = q-value threshold.</p>'
    )
    dmr_table = _dmr_sortable_table()
    dmr_block = _coll("DMR Volcano", dmr_desc + dmr_chart, open_=True)
    dmr_table_block = _coll("Significant DMRs", dmr_table, open_=True)

    return f"""
    <section class="section" id="part3">
      <div class="section-header">
        <span class="section-num">PART 03</span>
        <h2 class="section-title">Differential Analysis</h2>
        <span class="section-tag">diff</span>
      </div>

      <h3 class="subsec-title" id="part3_1">3.1 Violin</h3>
      {violin_block}

      <h3 class="subsec-title" id="part3_2">3.2 Heatmap</h3>
      {heatmap_block}

      <h3 class="subsec-title" id="part3_3">3.3 DMR Analysis</h3>
      {dmr_block}
      {dmr_table_block}
    </section>"""

def _sec_fragmentomics(rd, groups):

    all_samples = {s for members in groups.values() for s in members}

    def _sample_dd(analysis_dir, file_pattern, strip_suffix, uid_prefix,
                   group_mean_pat=None):
        all_pngs     = sorted(glob.glob(os.path.join(analysis_dir, file_pattern)))
        name_to_path = {}
        for p in all_pngs:
            stem  = os.path.splitext(os.path.basename(p))[0]
            clean = stem.replace(".markdup", "")
            for sfx in [strip_suffix, strip_suffix.replace("*", "")]:
                clean = clean.replace(sfx, "").strip("_").strip(".")
            for s in all_samples:
                if clean == s or stem.startswith(s):
                    name_to_path[s] = p
                    break
        blocks = ""
        for grp, members in groups.items():
            items = []
            # Group mean as first dropdown option, if available
            if group_mean_pat:
                gm = sorted(glob.glob(os.path.join(
                    analysis_dir, group_mean_pat.replace("{grp}", grp))))
                gm = [p for p in gm if os.path.isfile(p)]
                if gm:
                    items.append((f"{grp} (mean)", gm[0]))
            items += [(s, name_to_path[s]) for s in members if s in name_to_path]
            if items:
                blocks += (
                    f'<div style="margin-bottom:4px;font-size:14.5px;font-weight:500;'
                    f'color:#000000;font-family:var(--mono)">{grp}</div>'
                )
                blocks += _dropdown_gallery("Sample:", items, uid=f"{uid_prefix}_{grp}")
        return blocks

    def _grp_compare(analysis_dir, grp_pat, cmp_pat):
        html = ""
        for grp in (groups if grp_pat else []):
            imgs = sorted(glob.glob(os.path.join(analysis_dir, grp_pat.replace("{grp}", grp))))
            if imgs:
                html += (
                    f'<div class="fig-card" style="margin-bottom:16px;">'
                    f'<div class="fig-img tall">{_img_tag(imgs[0], f"{grp} mean")}</div>'
                    f'<div class="fig-body"><div class="fig-label">{grp} — Group Mean</div>'
                    f'</div></div>'
                )
        if cmp_pat:
            cmp_imgs = [p for p in sorted(glob.glob(os.path.join(analysis_dir, cmp_pat)))
                        if os.path.isfile(p)]
            if cmp_imgs:
                html += (
                    f'<div class="fig-card" style="margin-bottom:16px;">'
                    f'<div class="fig-img tall">{_img_tag(cmp_imgs[0], "Comparison")}</div>'
                    f'<div class="fig-body"><div class="fig-label">Group Comparison</div>'
                    f'</div></div>'
                )
        return html

    occupancy_dir = os.path.join(rd, "4_fragmentomics", "occupancy")
    delfi_dir = os.path.join(rd, "4_fragmentomics", "delfi")
    em_dir    = os.path.join(rd, "4_fragmentomics", "end_motif")
    cl_dir    = os.path.join(rd, "4_fragmentomics", "cleavage")
    wps_dir   = os.path.join(rd, "4_fragmentomics", "wps")

    # End motif: one interactive chart per group (data keys: end_motif_{grp})
    em_charts = ""
    for grp in groups:
        em_charts += (
            f'<div style="margin-bottom:6px;font-size:14.5px;font-weight:500;'
            f'color:#000000;font-family:var(--mono)">{grp}</div>'
        )
        em_charts += _plotly_select_render(f"end_motif_{grp}", "Sample:",
                                            "_sample_menu", height=440)
    if not em_charts:
        em_charts = "<p><em>No end-motif results.</em></p>"

    return f"""
    <section class="section" id="part4">
      <div class="section-header">
        <span class="section-num">PART 04</span>
        <h2 class="section-title">Fragmentomics</h2>
        <span class="section-tag">frag</span>
      </div>
      <h3 class="subsec-title" id="part4_1">4.1 Nucleosome Occupancy</h3>
      {_sample_dd(occupancy_dir,"*_occupancy.png","_occupancy","occupancy") or "<p><em>No per-sample occupancy results.</em></p>"}
      <h3 class="subsec-title" id="part4_2">4.2 DELFI</h3>
      {_sample_dd(delfi_dir,"*_delfi_genome.png","_delfi_genome","delfi",group_mean_pat="delfi_{grp}.png") or "<p><em>No per-sample DELFI results.</em></p>"}
      {_grp_compare(delfi_dir,None,"delfi_comparison.png") or ""}
      <h3 class="subsec-title" id="part4_3">4.3 End Motif</h3>
      <p style="font-size:14.5px;color:#000000;margin-bottom:8px;">
        Top 20 4-mer end motifs, shown separately per group. Use each chart's
        dropdown to switch between the group mean and individual samples.
      </p>
      {em_charts}
      <div style="margin-top:18px;font-size:14.5px;color:#000000;margin-bottom:8px;">
        <strong style="color:#000000;">Group comparison:</strong>
        Box plots of the top 20 motifs (by overall mean). Each box summarises the
        group distribution; individual samples are overlaid as points.
      </div>
      {_plotly_from_data("end_motif_box", height=480)}
      <h3 class="subsec-title" id="part4_4">4.4 Cleavage</h3>
      {_grp_compare(cl_dir,"cleavage_{grp}_samples.png","cleavage_comparison.png") or "<p><em>No cleavage results.</em></p>"}
      <h3 class="subsec-title" id="part4_5">4.5 WPS</h3>
      {_sample_dd(wps_dir,"*.wps_profile.png",".wps_profile","wps") or "<p><em>No WPS results.</em></p>"}
    </section>"""


def _sec_mesa(rd):
    perf_tsv = os.path.join(rd, "5_mesa", "modality_performance.tsv")
    roc_chart  = _plotly_from_data("mesa_roc",      height=430, square=True)
    hmap_chart = _plotly_from_data("mesa_heatmap",  height=300)
    sp_chart   = _plotly_from_data("mesa_spearman", height=320, square=True)
    return f"""
    <section class="section" id="part5">
      <div class="section-header">
        <span class="section-num">PART 05</span>
        <h2 class="section-title">MESA Multimodal Modeling</h2>
        <span class="section-tag">mesa</span>
      </div>
      <h3 class="subsec-title" id="part5_1">5.1 Modality Performance</h3>
      {_read_perf_table(perf_tsv)}
      <h3 class="subsec-title" id="part5_2">5.2 ROC Curve</h3>
      {_coll("ROC curves", f'''
        <p style="font-size:14.5px;color:#000000;margin-bottom:8px;">
          LOOCV ROC curves per modality. AUC shown in legend. Hover for FPR/TPR values.
        </p>
        {roc_chart}''', open_=True)}
      <h3 class="subsec-title" id="part5_3">5.3 Prediction Heatmap</h3>
      {_coll("LOOCV prediction probabilities", f'''
        <p style="font-size:14.5px;color:#000000;margin-bottom:8px;">
          Samples sorted by true label then probability. Hover cells for exact values.
        </p>
        {hmap_chart}''', open_=True)}
      <h3 class="subsec-title" id="part5_4">5.4 Spearman Correlation</h3>
      {_coll("Modality Spearman correlation", f'''
        <p style="font-size:14.5px;color:#000000;margin-bottom:8px;">
          Pairwise Spearman ρ of LOOCV prediction scores across modalities.
        </p>
        {sp_chart}''', open_=True)}
    </section>"""


def _build_sidebar(rd, groups, include_power=False):
    power_link = (
        '<a href="#part0" class="nav-link nav-top"><span class="dot"></span>Power Analysis</a>'
        if include_power else ""
    )
    return f"""
    <div class="nav-section">
      <div class="nav-label">Content</div>
      <a href="#workflow_summary" class="nav-link nav-top"><span class="dot"></span>Workflow Summary</a>
      {power_link}
      <a href="#part_overview" class="nav-link nav-top"><span class="dot"></span>Sample Statistics</a>
      <a href="#part1"  class="nav-link nav-top"></span>1 Processing</a>
      <a href="#part1_1" class="nav-link nav-sub"></span>1.1 Trimming</a>
      <a href="#part1_2" class="nav-link nav-sub"></span>1.2 Trimmed QC</a>
      <a href="#part1_3" class="nav-link nav-sub"></span>1.3 Alignment</a>
      <a href="#part1_4" class="nav-link nav-sub"></span>1.4 M-bias</a>
      <a href="#part1_6" class="nav-link nav-sub"></span>1.5 Sequencing QC</a>
      <a href="#part2"  class="nav-link nav-top"></span>2 cfDNA QC Analysis</a>
      <a href="#part2_1" class="nav-link nav-sub"></span>2.1 Methylation</a>
      <a href="#part2_2" class="nav-link nav-sub"></span>2.2 Fragment Length</a>
      <a href="#part2_3" class="nav-link nav-sub"></span>2.3 Dinucleotide</a>
      <a href="#part2_4" class="nav-link nav-sub"></span>2.4 PCA</a>
      <a href="#part3"  class="nav-link nav-top"></span>3 Differential</a>
      <a href="#part3_1" class="nav-link nav-sub"></span>3.1 Violin</a>
      <a href="#part3_2" class="nav-link nav-sub"></span>3.2 Heatmap</a>
      <a href="#part3_3" class="nav-link nav-sub"></span>3.3 DMR</a>
      <a href="#part4"  class="nav-link nav-top"></span>4 Fragmentomics</a>
      <a href="#part4_1" class="nav-link nav-sub"></span>4.1 Occupancy</a>
      <a href="#part4_2" class="nav-link nav-sub"></span>4.2 DELFI</a>
      <a href="#part4_3" class="nav-link nav-sub"></span>4.3 End Motif</a>
      <a href="#part4_4" class="nav-link nav-sub"></span>4.4 Cleavage</a>
      <a href="#part4_5" class="nav-link nav-sub"></span>4.5 WPS</a>
      <a href="#part5"  class="nav-link nav-top"></span>5 MESA</a>
      <a href="#part5_1" class="nav-link nav-sub"></span>5.1 Performance</a>
      <a href="#part5_2" class="nav-link nav-sub"></span>5.2 ROC</a>
      <a href="#part5_3" class="nav-link nav-sub"></span>5.3 LOOCV Heatmap</a>
      <a href="#part5_4" class="nav-link nav-sub"></span>5.4 LOOCV Spearman</a>
      <a href="#part6"  class="nav-link nav-top"></span>6 Software &amp; Tools</a>
    </div>"""


# ── Main entry ────────────────────────────────────────────────────────────────

def generate_report(args):
    rd       = getattr(args, "results_dir", "results")
    out_dir  = getattr(args, "output_dir",  "results/report")
    project  = getattr(args, "project_name", "cftk_project")
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    cfg_path = getattr(args, "config", None)
    groups   = {}
    group_a  = ""
    group_b  = ""
    q_thr    = 0.05
    frag_len = 167
    matrix_path = None

    cfg = getattr(args, "resolved_config", None)
    if cfg is None and cfg_path and os.path.exists(cfg_path):
        with open(cfg_path) as f:
            cfg = json.load(f)
    if cfg is not None:
        for grp, members in cfg.get("samples", {}).items():
            groups[grp] = [s["name"] for s in members]
        cmp = cfg.get("comparison", "") or ""
        if "_vs_" in cmp:
            group_a, group_b = cmp.split("_vs_", 1)
        elif cfg.get("control_group") and cfg.get("case_group"):
            group_a = cfg["control_group"]
            group_b = cfg["case_group"]
        q_thr    = cfg.get("analysis", {}).get("dmr", {}).get("params", {}).get("q_thr", 0.05)
        frag_len = cfg.get("analysis", {}).get("qc", {}).get("params", {}).get("fragment", 167)
        matrix_path = os.path.join(rd, "1_process", "5_merged_matrix", "cpg_matrix.tsv")
    else:
        for g in getattr(args, "groups", []):
            groups[g] = []

    all_sample_names = [s for members in groups.values() for s in members]
    source_results = _source_results_from_config(cfg, rd)
    qc_scores = _read_qc_scores(rd, groups, source_results=source_results)
    os.makedirs(out_dir, exist_ok=True)

    # ── Build interactive chart data ──────────────────────────────────────────
    data_dir = os.path.join(out_dir, "data")
    try:
        from report.data_builder import build_all
        chart_data = build_all(
            rd=rd,
            data_dir=data_dir,
            group_labels=groups,
            group_a=group_a,
            group_b=group_b,
            q_thr=q_thr,
            frag_len=frag_len,
            matrix_path=matrix_path,
            source_results=source_results,
        )
    except Exception as e:
        import traceback
        print(f"[report] WARNING: data_builder failed: {e}")
        traceback.print_exc()
        chart_data = {}

    # Embed chart data as per-key <script> blocks.
    # Each key uses Base64-encoded JSON to avoid any quote/backslash/</script> issues.
    # The JS decodes on first access (lazy getter) so parsing is deferred.
    import base64 as _b64
    per_key_scripts = []
    for k, fig in chart_data.items():
        fig_b64 = _b64.b64encode(
            json.dumps(fig, separators=(",", ":")).encode("utf-8")
        ).decode("ascii")
        # Build JS using string concatenation — avoids f-string / brace conflicts
        js = (
            "<script>(function(){"
            "window.__CFTK_DATA__=window.__CFTK_DATA__||{};"
            "var b=\"" + fig_b64 + "\";"
            "Object.defineProperty(window.__CFTK_DATA__,\"" + k + "\",{"
            "get:function(){"
            "var v=JSON.parse(atob(b));"
            "Object.defineProperty(window.__CFTK_DATA__,\"" + k + "\"," 
            "{value:v,configurable:true,writable:true});"
            "return v;},"
            "configurable:true"
            "});})();</scr"+"ipt>"
        )
        # Fix self-closing script tag back to normal
        # tag already correct above
        per_key_scripts.append(js)
    data_js = "\n".join(per_key_scripts)

    # ── Build HTML ────────────────────────────────────────────────────────────
    include_power = bool(glob.glob(os.path.join(rd, "0_power", "*")))
    sidebar  = _build_sidebar(rd, groups, include_power=include_power)
    sections = (
        _sec_workflow_summary(rd, groups)
        + (_sec_power(rd) if include_power else "")
        + _sec_sample_statistics(
            rd, groups, source_results=source_results, qc_scores=qc_scores
        )
        + _sec_process(
            rd,
            groups=groups,
            all_sample_names=all_sample_names,
            source_results=source_results,
            qc_scores=qc_scores,
        )
        + _sec_qc(
            rd, groups=groups, source_results=source_results, qc_scores=qc_scores
        )
        + _sec_differential(rd)
        + _sec_fragmentomics(rd, groups)
        + _sec_mesa(rd)
        + _sec_software()
    )

    tmpl_path = os.path.join(os.path.dirname(__file__), "report_template.html")
    with open(tmpl_path) as f:
        template = f.read()

    # Inject Plotly.js + data block before </head>
    head_inject = f"{PLOTLY_SCRIPT}\n{data_js}"
    if PLOTLY_SCRIPT not in template:
        template = template.replace("</head>", f"{head_inject}\n</head>", 1)
    else:
        # Plotly already present; inject only data block before </head>
        template = template.replace("</head>", f"{data_js}\n</head>", 1)

    replacements = {
        "<!-- {SIDEBAR} -->":      sidebar,
        "<!-- {SECTIONS} -->":     sections,
        "<!-- {PROJECT_NAME} -->": project,
        "<!-- {GROUPS} -->":       " · ".join(groups.keys()),
        "<!-- {DATE} -->":         date_str,
        "<!-- {REPORT_DATE} -->":  date_str,
    }
    html = template
    for k, v in replacements.items():
        html = html.replace(k, v)

    out_path = os.path.join(out_dir, "report.html")
    with open(out_path, "w") as f:
        f.write(html)
    print(f"[report] saved → {out_path}")
    return out_path
