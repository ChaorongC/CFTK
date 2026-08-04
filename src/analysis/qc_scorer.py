"""
qc_scorer.py — Per-sample QC scoring for cfDNA methylation sequencing.

Reads qc_summary.tsv produced by qc_parser.collect_qc_metrics(),
applies per-metric thresholds, computes a weighted composite score (0–100),
assigns PASS / WARN / FAIL status, and generates human-readable
recommendations.

Usage (standalone):
    from analysis.qc_scorer import score_samples
    scores_df = score_samples(summary_df, out_tsv="2_qc/qc_scores.tsv")
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
from typing import NamedTuple

# ── Threshold definitions ──────────────────────────────────────────────────────
#
# Each MetricRule defines:
#   col      : column name in qc_summary.tsv
#   label    : human-readable name
#   weight   : contribution to composite score (weights sum to ~100)
#   pass_fn  : callable(value) -> bool   — PASS condition
#   warn_fn  : callable(value) -> bool   — WARN condition (else FAIL)
#   higher_is_better : used only for display / trend arrows
#   fmt      : format string for display
#   note     : short description shown in report


class MetricRule(NamedTuple):
    col: str
    label: str
    weight: float
    pass_lo: float  # lower bound for PASS (inclusive)
    pass_hi: float  # upper bound for PASS (inclusive)
    warn_lo: float  # lower bound for WARN (inclusive)
    warn_hi: float  # upper bound for WARN (inclusive)
    higher_is_better: bool = True
    fmt: str = ".1f"
    unit: str = "%"
    note: str = ""


# ── Rules (N9) ────────────────────────────────────────────────────────────────
#
# Status logic:
#   PASS : pass_lo <= value <= pass_hi
#   WARN : warn_lo <= value <= warn_hi  (and not PASS)
#   FAIL : everything else
#
# For metrics where lower is better (dup%, error rate),
# pass_lo/warn_lo are 0 and thresholds are upper bounds.

RULES: list[MetricRule] = [
    # ── Sequencing quality ───────────────────────────────────────────────────
    MetricRule(
        col="flagstat_mapped_pct",
        label="Mapping rate",
        weight=20,
        pass_lo=80, pass_hi=100,
        warn_lo=60, warn_hi=100,
        higher_is_better=True,
        fmt=".1f", unit="%",
        note="Fraction of reads mapped to reference. <60% indicates serious problems.",
    ),
    MetricRule(
        col="flagstat_properly_paired_pct",
        label="Properly paired",
        weight=8,
        pass_lo=75, pass_hi=100,
        warn_lo=55, warn_hi=100,
        higher_is_better=True,
        fmt=".1f", unit="%",
        note="Fraction of read pairs mapping concordantly.",
    ),
    MetricRule(
        col="markdup_dup_pct",
        label="Duplication rate",
        weight=12,
        pass_lo=0,  pass_hi=30,
        warn_lo=0,  warn_hi=50,
        higher_is_better=False,
        fmt=".1f", unit="%",
        note="PCR/optical duplicate fraction. >50% severely reduces usable depth.",
    ),
    MetricRule(
        col="stats_error_rate",
        label="Sequencing error rate",
        weight=5,
        pass_lo=0,  pass_hi=1.0,
        warn_lo=0,  warn_hi=2.0,
        higher_is_better=False,
        fmt=".3f", unit="%",
        note="Base substitution error rate from samtools stats.",
    ),
    MetricRule(
        col="adapter_reads_pct",
        label="Adapter reads",
        weight=0,   # excluded from QC table display
        pass_lo=0,  pass_hi=30,
        warn_lo=0,  warn_hi=50,
        higher_is_better=False,
        fmt=".1f", unit="%",
        note="Fraction of reads containing adapter sequence before trimming.",
    ),
    MetricRule(
        col="fqc_pct_modules_failed",
        label="FastQC modules failed",
        weight=0,   # excluded from QC table display
        pass_lo=0,  pass_hi=10,
        warn_lo=0,  warn_hi=30,
        higher_is_better=False,
        fmt=".0f", unit="%",
        note="Fraction of FastQC QC modules that failed.",
    ),
    # ── Methylation quality ──────────────────────────────────────────────────
    MetricRule(
        col="cpg_global_meth_pct",
        label="Global CpG methylation",
        weight=10,
        pass_lo=40, pass_hi=90,
        warn_lo=25, warn_hi=95,
        higher_is_better=True,   # context-dependent; middle range is normal
        fmt=".1f", unit="%",
        note="Genome-wide weighted mean CpG methylation. "
             "Human cfDNA typically 60–85%.",
    ),
    MetricRule(
        col="cpg_covered_sites",
        label="CpG covered sites",
        weight=8,
        pass_lo=500_000, pass_hi=float("inf"),
        warn_lo=100_000, warn_hi=float("inf"),
        higher_is_better=True,
        fmt=".0f", unit="",
        note="Number of CpG sites with ≥min_depth coverage.",
    ),
    MetricRule(
        col="cpg_mean_depth",
        label="Mean CpG depth",
        weight=5,
        pass_lo=10, pass_hi=float("inf"),
        warn_lo=5,  warn_hi=float("inf"),
        higher_is_better=True,
        fmt=".1f", unit="×",
        note="Mean read depth over covered CpG sites.",
    ),
    # ── cfDNA fragment features ──────────────────────────────────────────────
    MetricRule(
        col="median_frag_len",
        label="Median fragment length",
        weight=0,   # excluded from QC table display
        pass_lo=130, pass_hi=200,
        warn_lo=100, warn_hi=250,
        higher_is_better=True,
        fmt=".0f", unit=" bp",
        note="Median insert size. cfDNA typically peaks at ~167 bp (mononucleosome).",
    ),
    MetricRule(
        col="short_frag_ratio",
        label="Short fragment ratio (<120 bp)",
        weight=0,   # informational only — no penalty, cfDNA naturally has short frags
        pass_lo=0,  pass_hi=100,
        warn_lo=0,  warn_hi=100,
        higher_is_better=False,
        fmt=".1f", unit="%",
        note="Fraction of fragments <120 bp. Elevated in highly degraded samples.",
    ),
    MetricRule(
        col="mono_nuc_ratio",
        label="Mononucleosome ratio (140–180 bp)",
        weight=0,   # informational only
        pass_lo=0,  pass_hi=100,
        warn_lo=0,  warn_hi=100,
        higher_is_better=True,
        fmt=".1f", unit="%",
        note="Fraction of fragments in mononucleosome range (140–180 bp).",
    ),
    # ── Beta-density QC (computed from cpg_matrix.tsv) ──────────────────────
    MetricRule(
        col="beta_M_score",
        label="β M-score",
        weight=0,   # shown only in beta-density table, not main QC table
        pass_lo=2.0,  pass_hi=float("inf"),
        warn_lo=2.0,  warn_hi=float("inf"),   # no WARN — binary PASS/FAIL
        higher_is_better=True,
        fmt=".2f", unit="",
        note="min(left_peak, right_peak) / mid_density. "
             "PASS ≥ 2: bimodal distribution with strong peaks at 0 and 1.",
    ),
    MetricRule(
        col="beta_peak_balance",
        label="β Peak Balance",
        weight=0,   # shown only in beta-density table, not main QC table
        pass_lo=0.30, pass_hi=float("inf"),
        warn_lo=0.30, warn_hi=float("inf"),   # no WARN — binary PASS/FAIL
        higher_is_better=True,
        fmt=".3f", unit="",
        note="min/max of left and right peaks. "
             "PASS ≥ 0.30: balanced methylation/unmethylation peaks.",
    ),
    MetricRule(
        col="beta_boundary_mid_ratio",
        label="β Boundary/Mid Ratio",
        weight=0,   # shown only in beta-density table, not main QC table
        pass_lo=2.0,  pass_hi=float("inf"),
        warn_lo=2.0,  warn_hi=float("inf"),   # no WARN — binary PASS/FAIL
        higher_is_better=True,
        fmt=".2f", unit="",
        note="frac(β≤0.15 or β≥0.85) / frac(0.35≤β≤0.65). "
             "PASS ≥ 2: most CpGs are fully methylated or unmethylated.",
    ),
]

# Normalise weights so active rules (weight > 0) sum to 100
_ACTIVE_RULES = [r for r in RULES if r.weight > 0]
_WEIGHT_SUM   = sum(r.weight for r in _ACTIVE_RULES)


# ── Status helper ─────────────────────────────────────────────────────────────

def _get_status(value: float, rule: MetricRule) -> str:
    """Return 'PASS', 'WARN', or 'FAIL' for a single metric value."""
    if np.isnan(value):
        return "NA"
    if rule.pass_lo <= value <= rule.pass_hi:
        return "PASS"
    if rule.warn_lo <= value <= rule.warn_hi:
        return "WARN"
    return "FAIL"


def _score_metric(value: float, rule: MetricRule) -> float:
    """Return a 0–1 score contribution for one metric."""
    status = _get_status(value, rule)
    if status == "NA":
        return 0.5   # missing data — neutral contribution
    if status == "PASS":
        return 1.0
    if status == "WARN":
        return 0.5
    return 0.0       # FAIL


# ── Main scoring function ─────────────────────────────────────────────────────

def score_samples(
    summary: pd.DataFrame,
    out_tsv: str | None = None,
) -> pd.DataFrame:
    """
    Compute per-sample QC scores and PASS/WARN/FAIL status.

    Parameters
    ----------
    summary : DataFrame from qc_parser.collect_qc_metrics()
    out_tsv : optional path to write qc_scores.tsv

    Returns
    -------
    DataFrame with columns:
        sample, group, qc_score, qc_status,
        recommendation,
        {metric}_status  for every rule,
        {metric}_value   (raw value, repeated for convenience)
    """
    rows = []

    for _, row in summary.iterrows():
        result: dict = {
            "sample": row.get("sample", ""),
            "group":  row.get("group",  ""),
        }

        # Per-metric status and weighted score
        weighted_score = 0.0
        fail_metrics:  list[str] = []
        warn_metrics:  list[str] = []
        na_metrics:    list[str] = []

        for rule in RULES:
            raw = row.get(rule.col, np.nan)
            try:
                value = float(raw)
            except (TypeError, ValueError):
                value = np.nan

            status = _get_status(value, rule)
            result[f"{rule.col}_status"] = status
            result[f"{rule.col}_value"]  = value

            if rule.weight > 0:
                contrib = _score_metric(value, rule) * rule.weight / _WEIGHT_SUM
                weighted_score += contrib

            if status == "FAIL" and rule.weight > 0:
                fail_metrics.append(rule.label)
            elif status == "WARN" and rule.weight > 0:
                warn_metrics.append(rule.label)
            elif status == "NA" and rule.weight > 0:
                na_metrics.append(rule.label)

        # Composite score 0–100
        qc_score = round(weighted_score * 100, 1)
        result["qc_score"] = qc_score

        # Overall status: any FAIL → FAIL; any WARN → WARN; else PASS
        # Critical FAILs (weight >= 15) immediately force FAIL regardless of score
        critical_fail = any(
            _get_status(
                float(row.get(r.col, np.nan)) if not pd.isna(row.get(r.col, np.nan)) else np.nan,
                r
            ) == "FAIL"
            for r in _ACTIVE_RULES if r.weight >= 15
        )

        if fail_metrics or critical_fail:
            qc_status = "FAIL"
        elif warn_metrics:
            qc_status = "WARN"
        else:
            qc_status = "PASS"

        result["qc_status"] = qc_status

        # Human-readable recommendation
        recommendation = _build_recommendation(
            qc_status, qc_score, fail_metrics, warn_metrics, na_metrics
        )
        result["recommendation"] = recommendation

        rows.append(result)

    scores_df = pd.DataFrame(rows)

    if out_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
        scores_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"[qc_scorer] qc_scores.tsv → {out_tsv}  ({len(scores_df)} samples)")

    return scores_df


# ── Recommendation text ───────────────────────────────────────────────────────

def _build_recommendation(
    status: str,
    score: float,
    fail_metrics: list[str],
    warn_metrics: list[str],
    na_metrics: list[str],
) -> str:
    if status == "PASS":
        return f"Sample passes all QC thresholds (score {score}/100). Suitable for downstream analysis."

    parts: list[str] = []

    if status == "FAIL":
        parts.append(f"⚠ RECOMMEND EXCLUSION (score {score}/100).")
        if fail_metrics:
            parts.append(f"Failed metrics: {', '.join(fail_metrics)}.")
        if warn_metrics:
            parts.append(f"Also flagged: {', '.join(warn_metrics)}.")
        parts.append(
            "Review raw data quality and consider re-sequencing or excluding from analysis."
        )
    else:  # WARN
        parts.append(f"Sample passes with warnings (score {score}/100).")
        if warn_metrics:
            parts.append(f"Flagged metrics: {', '.join(warn_metrics)}.")
        parts.append(
            "Proceed with caution; verify these metrics do not affect downstream results."
        )

    if na_metrics:
        parts.append(f"Missing data for: {', '.join(na_metrics)}.")

    return " ".join(parts)


# ── Summary table for report ──────────────────────────────────────────────────

def build_display_table(scores_df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a compact display DataFrame suitable for the HTML report table.
    Columns: sample, group, qc_score, qc_status,
             + one column per active rule (value + status suffix).
    """
    cols = ["sample", "group", "qc_score", "qc_status"]
    for rule in RULES:
        val_col = f"{rule.col}_value"
        if val_col in scores_df.columns:
            cols.append(val_col)
    existing = [c for c in cols if c in scores_df.columns]
    return scores_df[existing].copy()


def get_rules() -> list[MetricRule]:
    """Expose rules list for use in report generation."""
    return RULES
