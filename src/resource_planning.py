"""Shared CPU-budget planning for CFTK processing and QC workflows."""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Optional


SCHEDULER_CPU_VARIABLES = (
    ("slurm", "SLURM_CPUS_PER_TASK"),
    ("slurm", "SLURM_CPUS_ON_NODE"),
    ("pbs", "PBS_NP"),
    ("sge", "NSLOTS"),
)


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def plan_parallelism(total_cores, parallel_samples, sample_count=None):
    """Resolve nested sample/thread parallelism under one total CPU budget."""
    total_cores = _positive_integer(total_cores, "total_cores")
    parallel_samples = _positive_integer(parallel_samples, "parallel_samples")
    if parallel_samples > total_cores:
        raise ValueError(
            f"parallel_samples ({parallel_samples}) cannot exceed the total "
            f"core budget ({total_cores})"
        )
    if sample_count is not None:
        if isinstance(sample_count, bool) or not isinstance(sample_count, int) or sample_count < 0:
            raise ValueError("sample_count must be a non-negative integer")
        concurrent_samples = min(parallel_samples, max(1, sample_count))
    else:
        concurrent_samples = parallel_samples
    threads_per_sample = max(1, total_cores // concurrent_samples)
    return {
        "total_core_budget": total_cores,
        "requested_parallel_samples": parallel_samples,
        "concurrent_samples": concurrent_samples,
        "threads_per_sample": threads_per_sample,
        "estimated_peak_threads": concurrent_samples * threads_per_sample,
    }


def detect_scheduler_allocation(
    environ: Optional[Mapping[str, str]] = None,
):
    """Return the first recognized scheduler CPU allocation, if present."""
    environ = os.environ if environ is None else environ
    for scheduler, variable in SCHEDULER_CPU_VARIABLES:
        raw = environ.get(variable)
        if raw is None:
            continue
        try:
            cores = int(raw)
        except (TypeError, ValueError):
            return {
                "scheduler": scheduler,
                "variable": variable,
                "raw_value": str(raw),
                "allocated_cores": None,
                "valid": False,
            }
        return {
            "scheduler": scheduler,
            "variable": variable,
            "raw_value": str(raw),
            "allocated_cores": cores,
            "valid": cores > 0,
        }
    return {
        "scheduler": None,
        "variable": None,
        "raw_value": None,
        "allocated_cores": None,
        "valid": True,
    }


def ensure_scheduler_capacity(total_cores, allocation=None):
    """Reject a known scheduler allocation smaller than the configured budget."""
    total_cores = _positive_integer(total_cores, "total_cores")
    allocation = detect_scheduler_allocation() if allocation is None else allocation
    if not allocation.get("valid", False):
        raise ValueError(
            f"{allocation.get('variable')} is not a positive integer: "
            f"{allocation.get('raw_value')!r}"
        )
    allocated = allocation.get("allocated_cores")
    if allocated is not None and total_cores > allocated:
        raise ValueError(
            f"total core budget ({total_cores}) exceeds the scheduler allocation "
            f"({allocated} from {allocation.get('variable')})"
        )
    return allocation
