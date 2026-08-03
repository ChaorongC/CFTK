"""Pure helpers for the Streamlit model-power calculator."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence


PRECISION_MODES = {
    "Fast": {
        "n_templates": 5,
        "simulations_per_template": 10,
    },
    "Standard": {
        "n_templates": 10,
        "simulations_per_template": 20,
    },
}

RATIO_LABELS = {
    "1:3": 1 / 3,
    "1:2": 1 / 2,
    "1:1": 1.0,
    "2:1": 2.0,
    "3:1": 3.0,
}

SD_STAT_LABELS = {
    "Central estimate": "mean",
    "Lower SD estimate": "CI_l",
    "Upper SD estimate": "CI_u",
}

EFFECT_DIRECTION_LABELS = {
    "Balanced hyper/hypomethylation": "balanced",
    "Hypermethylation only": "positive",
    "Hypomethylation only": "negative",
    "Random directions": "random",
}

MAX_TOTAL_SAMPLE_SIZE = 2000
MAX_PUBLIC_WORK_UNITS = 500_000

SUMMARY_COLUMNS = (
    "sample_size",
    "mean_depth",
    "mean_cv_auc",
    "target_attainment_probability",
    "power",
    "power_ci_low",
    "power_ci_high",
    "mean_sensitivity_at_specificity",
    "mean_feature_recall",
    "mean_feature_precision",
    "mean_selection_jaccard",
    "n_templates",
    "simulations_per_template",
    "total_simulations",
)

PROBABILITY_COLUMNS = (
    "mean_cv_auc",
    "target_attainment_probability",
    "power",
    "power_ci_low",
    "power_ci_high",
    "mean_sensitivity_at_specificity",
    "mean_feature_recall",
    "mean_feature_precision",
    "mean_selection_jaccard",
)


class AppValidationError(ValueError):
    """Validation error intended to be displayed as a concise app message."""


@dataclass(frozen=True)
class ParsedSampleSizes:
    values: tuple[int, ...]
    warnings: tuple[str, ...] = ()


def parse_sample_size_text(text: str, *, max_points: int = 8) -> ParsedSampleSizes:
    """Parse comma-separated total sample sizes."""
    if not text or not text.strip():
        raise AppValidationError("Enter at least one total sample size.")

    raw_parts = [part.strip() for part in text.split(",")]
    values = []
    for part in raw_parts:
        if not part:
            continue
        try:
            value = int(part)
        except ValueError as exc:
            raise AppValidationError(
                f"Invalid sample size {part!r}; use positive integers."
            ) from exc
        if value <= 0:
            raise AppValidationError("Sample sizes must be positive integers.")
        values.append(value)

    return _normalize_sample_sizes(values, max_points=max_points)


def parse_sample_size_range(
    *,
    start: int,
    stop: int,
    step: int,
    max_points: int = 8,
) -> ParsedSampleSizes:
    """Parse start/stop/step total sample-size settings."""
    if start <= 0 or stop <= 0 or step <= 0:
        raise AppValidationError("Sample-size start, stop, and step must be positive.")
    if stop < start:
        raise AppValidationError("Sample-size stop must be greater than or equal to start.")
    values = list(range(int(start), int(stop) + 1, int(step)))
    return _normalize_sample_sizes(values, max_points=max_points)


def _normalize_sample_sizes(
    values: Sequence[int],
    *,
    max_points: int,
) -> ParsedSampleSizes:
    if not values:
        raise AppValidationError("Enter at least one total sample size.")
    unique = tuple(sorted(set(int(value) for value in values)))
    warnings = []
    if len(unique) != len(values):
        warnings.append("Duplicate sample sizes were removed.")
    if len(unique) > max_points:
        raise AppValidationError(
            f"Select at most {max_points} sample-size points per request."
        )
    if max(unique) > MAX_TOTAL_SAMPLE_SIZE:
        raise AppValidationError(
            "The public calculator supports total sample sizes up to "
            f"{MAX_TOTAL_SAMPLE_SIZE:,}. Reduce the largest sample size."
        )
    return ParsedSampleSizes(values=unique, warnings=tuple(warnings))


def case_control_ratio(label: str, *, custom_ratio: float | None = None) -> float:
    """Convert a user-facing case:control label to n_cases / n_controls."""
    if label == "Custom":
        if custom_ratio is None:
            raise AppValidationError("Enter a custom case-to-control ratio.")
        ratio = float(custom_ratio)
    else:
        try:
            ratio = float(RATIO_LABELS[label])
        except KeyError as exc:
            raise AppValidationError(f"Unsupported case-to-control ratio: {label}.") from exc
    if ratio <= 0:
        raise AppValidationError("Case-to-control ratio must be positive.")
    return ratio


def class_counts(total_sample_size: int, ratio: float) -> tuple[int, int]:
    """Return rounded cases and controls using the simulation engine convention."""
    if total_sample_size <= 0:
        raise AppValidationError("Sample sizes must be positive integers.")
    if ratio <= 0:
        raise AppValidationError("Case-to-control ratio must be positive.")
    n_cases = int(round(total_sample_size * ratio / (1.0 + ratio)))
    n_controls = int(total_sample_size) - n_cases
    return n_cases, n_controls


def validate_class_counts(
    *,
    sample_sizes: Sequence[int],
    ratio: float,
    cv_folds: int,
) -> None:
    """Require every requested design to have enough cases and controls per fold."""
    failures = []
    for sample_size in sample_sizes:
        n_cases, n_controls = class_counts(int(sample_size), ratio)
        if n_cases < cv_folds or n_controls < cv_folds:
            failures.append((int(sample_size), n_cases, n_controls))
    if failures:
        sample_size, n_cases, n_controls = failures[0]
        raise AppValidationError(
            "Each selected sample size must include at least "
            f"{cv_folds} cases and {cv_folds} controls. "
            f"Sample size {sample_size} gives {n_cases} cases and "
            f"{n_controls} controls."
        )


def available_depths_from_manifest(reference_dir: str | Path) -> tuple[int | float, ...]:
    """Read available mean-depth values from the model-power manifest."""
    from analysis.model_power_reference import load_manifest

    try:
        manifest = load_manifest(reference_dir)
        labels = manifest["depth_labels"]
    except (OSError, TypeError, ValueError, KeyError) as exc:
        raise AppValidationError(str(exc)) from exc
    return tuple(_coerce_depth_label(value) for value in labels)


def _coerce_depth_label(value: object) -> int | float:
    numeric = float(value)
    return int(numeric) if numeric.is_integer() else numeric


def validate_depth_selection(depths: Sequence[int | float]) -> tuple[int | float, ...]:
    """Validate selected depth values."""
    selected = tuple(depths)
    if not selected:
        raise AppValidationError("Select at least one mean sequencing depth.")
    if len(selected) > 4:
        raise AppValidationError("Select at most four mean sequencing depths per run.")
    return selected


def validate_biomarker_inputs(
    *,
    n_features: int,
    n_signal_cpgs: int,
    top_k: int,
    meth_diff: float,
    effect_sd: float,
    within_block_rho: float,
    sd_stat: str,
    effect_direction: str = "balanced",
) -> None:
    """Validate biomarker and template-defining settings."""
    if not 1 <= int(n_features) <= 2000:
        raise AppValidationError("n_features must be between 1 and 2000.")
    if not 1 <= int(n_signal_cpgs) <= int(n_features):
        raise AppValidationError("n_signal_cpgs must satisfy 1 <= n_signal_cpgs <= n_features.")
    if not 1 <= int(top_k) <= int(n_features):
        raise AppValidationError("top_k must satisfy 1 <= top_k <= n_features.")
    if not 0.01 <= float(meth_diff) <= 0.15:
        raise AppValidationError("Mean absolute methylation difference must be between 0.01 and 0.15.")
    if not 0.0 <= float(effect_sd) <= 0.05:
        raise AppValidationError("Effect-size SD must be between 0.0 and 0.05.")
    if effect_direction not in {"positive", "negative", "balanced", "random"}:
        raise AppValidationError("Unsupported effect direction.")
    if not 0.0 <= float(within_block_rho) <= 0.8:
        raise AppValidationError("Within-block correlation must be between 0.0 and 0.8.")
    if sd_stat not in {"mean", "CI_l", "CI_u"}:
        raise AppValidationError("Unsupported SD uncertainty scenario.")


def validate_analysis_inputs(
    *,
    target_auc: float,
    min_observed_fraction: float,
    specificity_target: float,
    cv_folds: int,
) -> None:
    """Validate study-analysis settings that do not define templates."""
    if not 0.60 <= float(target_auc) <= 0.90:
        raise AppValidationError("Target cross-validated AUC must be between 0.60 and 0.90.")
    if int(cv_folds) not in {3, 5, 10}:
        raise AppValidationError("Cross-validation folds must be 3, 5, or 10.")
    if not 0.30 <= float(min_observed_fraction) <= 0.90:
        raise AppValidationError("Minimum observed fraction must be between 0.30 and 0.90.")
    if not 0.80 <= float(specificity_target) <= 0.99:
        raise AppValidationError("Specificity operating point must be between 0.80 and 0.99.")


def precision_settings(mode: str) -> dict[str, int]:
    """Return named precision settings."""
    try:
        return dict(PRECISION_MODES[mode])
    except KeyError as exc:
        raise AppValidationError(f"Invalid precision mode: {mode}.") from exc


def top_k_choices(n_features: int) -> tuple[int, ...]:
    """Return public top-k options bounded by n_features."""
    return tuple(value for value in (5, 10, 20, 50, 100) if value <= int(n_features))


def top_k_training_warning(
    *,
    top_k: int,
    sample_sizes: Sequence[int],
    ratio: float,
    cv_folds: int,
) -> str | None:
    """Warn when top_k is large relative to the smallest outer-training fold."""
    if not sample_sizes:
        return None
    smallest = min(int(value) for value in sample_sizes)
    n_cases, n_controls = class_counts(smallest, ratio)
    smallest_train_fold = min(n_cases, n_controls) * (cv_folds - 1) / cv_folds
    if top_k > smallest_train_fold:
        return (
            "Selected CpGs per fold is large relative to the smallest training "
            "fold; estimates may be unstable."
        )
    return None


def workload_warning(
    *,
    sample_sizes: Sequence[int],
    depths: Sequence[int | float],
    n_templates: int,
    simulations_per_template: int,
    default_sample_sizes: Sequence[int] = (50, 100, 150, 200),
    default_depths: int = 2,
    default_templates: int | None = None,
    default_simulations_per_template: int | None = None,
) -> str | None:
    """Enforce the public workload limit and flag unusually expanded designs.

    The warning baseline uses the selected precision mode's template and
    simulation counts by default. This prevents normal Fast or Standard runs
    from being mislabeled as unusually large after null calibration was added.
    """
    sample_sizes = tuple(int(value) for value in sample_sizes)
    depths = tuple(depths)
    requested = (
        sum(sample_sizes)
        * len(depths)
        * int(n_templates)
        * int(simulations_per_template)
    )

    baseline_templates = (
        int(n_templates)
        if default_templates is None
        else int(default_templates)
    )
    baseline_simulations = (
        int(simulations_per_template)
        if default_simulations_per_template is None
        else int(default_simulations_per_template)
    )
    default = (
        sum(int(value) for value in default_sample_sizes)
        * int(default_depths)
        * baseline_templates
        * baseline_simulations
    )

    if requested > MAX_PUBLIC_WORK_UNITS:
        raise AppValidationError(
            "This configuration exceeds the public calculator's computational "
            "limit. Reduce the sample-size grid, selected depths, or precision mode."
        )
    if requested > 2 * default:
        return (
            "This configuration is substantially larger than the selected "
            "precision mode's default workload and may run slowly on "
            "Streamlit Community Cloud."
        )
    return None


def build_power_kwargs(
    *,
    ratio: float,
    cv_folds: int,
    top_k: int,
    min_observed_fraction: float,
    target_auc: float,
    specificity_target: float,
) -> dict[str, object]:
    """Build fixed public-app kwargs for the discovery power runner."""
    return {
        "ratio": float(ratio),
        "models": ("logreg",),
        "cv_folds": int(cv_folds),
        "cv_repeats": 1,
        "top_k": int(top_k),
        "min_observed_fraction": float(min_observed_fraction),
        "target_auc": float(target_auc),
        "specificity_target": float(specificity_target),
        "paired_depths": True,
        "logreg_params": {
            "C": 1.0,
            "solver": "liblinear",
            "max_iter": 2000,
            "class_weight": None,
        },
    }


def build_template_kwargs(
    *,
    depths: Sequence[int | float],
    n_features: int,
    n_signal_cpgs: int,
    meth_diff: float,
    effect_sd: float,
    effect_direction: str,
    sd_stat: str,
    within_block_rho: float,
    block_size: int = 20,
) -> dict[str, object]:
    """Build template-defining kwargs for prepare_template_ensemble()."""
    return {
        "depth": tuple(depths),
        "n_features": int(n_features),
        "n_signal_cpgs": int(n_signal_cpgs),
        "meth_diff": float(meth_diff),
        "effect_sd": float(effect_sd),
        "effect_direction": effect_direction,
        "sd_stat": sd_stat,
        "within_block_rho": float(within_block_rho),
        "block_size": int(block_size),
    }


def result_metadata(
    *,
    user_inputs: Mapping[str, object],
    precision_mode: str,
    template_seed: int,
    grid_seed: int,
    ci_method: str,
    n_bootstrap: int,
) -> dict[str, object]:
    """Create reproducibility metadata for export."""
    precision = precision_settings(precision_mode)
    return {
        "user_inputs": dict(user_inputs),
        "precision_mode": precision_mode,
        "precision_settings": precision,
        "seeds": {
            "template_seed": int(template_seed),
            "grid_seed": int(grid_seed),
        },
        "fixed_settings": {
            "models": ("logreg",),
            "cv_repeats": 1,
            "block_size": 20,
            "n_jobs": 2,
            "paired_depths": True,
            "logistic_regression": {
                "C": 1.0,
                "solver": "liblinear",
                "max_iter": 2000,
                "class_weight": None,
            },
        },
        "ci_method": ci_method,
        "n_bootstrap": int(n_bootstrap) if ci_method == "hierarchical_bootstrap" else None,
    }


def display_table_columns(columns: Sequence[str], *, include_ci: bool) -> list[str]:
    """Return summary-table columns, hiding empty CI columns when disabled."""
    available = set(columns)
    result = []
    for column in SUMMARY_COLUMNS:
        if column not in available:
            continue
        if not include_ci and column in {"power_ci_low", "power_ci_high"}:
            continue
        result.append(column)
    return result
