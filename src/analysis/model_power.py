from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence
import time, os, warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.special import betaincinv, ndtr
from scipy.stats import mannwhitneyu, nbinom, truncnorm
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.proportion import proportion_confint
import matplotlib.pyplot as plt

from pathlib import Path


# Directory containing src/.
_MODULE_DIR = Path(__file__).resolve().parent.parent

# project/data/
_DEFAULT_REFERENCE_DIR = _MODULE_DIR.parent / "data"

cpg_std_summary: pd.DataFrame | None = None
cpg_mean: pd.Series | None = None


def load_default_model_power_reference(
    reference_dir: str | Path | None = None,
    *,
    depths: Sequence[float | int | str] | None = None,
    sd_stats: Sequence[str] = ("mean", "CI_l", "CI_u"),
    include_index: bool = True,
    mmap_mode: str | None = "r",
) -> Any:
    """
    Load the default model-power reference arrays and cache compatibility globals.

    This function intentionally performs all reference-data I/O explicitly, so
    importing ``analysis.model_power`` does not load large reference files.
    CFTK distributions contain the loader and simulation code, not the arrays.
    """
    from analysis.model_power_reference import (
        load_model_power_reference,
        resolve_model_power_reference_dir,
    )

    global cpg_mean, cpg_std_summary

    resolved_reference_dir = resolve_model_power_reference_dir(
        reference_dir,
        repository_dir=_DEFAULT_REFERENCE_DIR,
    )
    loaded = load_model_power_reference(
        resolved_reference_dir,
        depths=depths,
        sd_stats=sd_stats,
        include_index=include_index,
        mmap_mode=mmap_mode,
    )
    cpg_mean = loaded.cpg_mean
    cpg_std_summary = loaded.cpg_std_summary
    return loaded


def get_default_model_power_reference(
    reference_dir: str | Path | None = None,
    *,
    depths: Sequence[float | int | str] | None = None,
    sd_stats: Sequence[str] = ("mean", "CI_l", "CI_u"),
    include_index: bool = True,
    mmap_mode: str | None = "r",
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Return ``(cpg_std_summary, cpg_mean)`` for legacy notebook-style usage.

    The first call loads and caches the reference arrays. Pass ``depths`` or
    ``sd_stats`` to explicitly load a subset.
    """
    global cpg_mean, cpg_std_summary

    if (
        cpg_mean is None
        or cpg_std_summary is None
        or depths is not None
        or reference_dir is not None
    ):
        load_default_model_power_reference(
            reference_dir=reference_dir,
            depths=depths,
            sd_stats=sd_stats,
            include_index=include_index,
            mmap_mode=mmap_mode,
        )

    if cpg_std_summary is None or cpg_mean is None:
        raise RuntimeError("Model-power reference data could not be loaded.")

    return cpg_std_summary, cpg_mean


def get_nbinom_on_mean(mean_depth: float) -> tuple[float, float, float, float, float]:
    """Return the empirical CFTK negative-binomial depth parameters."""
    mean_depth = float(mean_depth)
    if mean_depth <= 0:
        raise ValueError("mean_depth must be positive.")

    k = 1.3409309326892476 * mean_depth + 3.526927718332715
    n = (mean_depth / k - 0.6291928395185754) / 0.02390549241316919
    p = n / (n + mean_depth)
    missing_pct = mean_depth * n * (-0.0007582750922856895) + 4.56905302776041

    return mean_depth, k, n, p, missing_pct / 100.0


# Backward-compatible alias for the old misspelled function name.
get_nbimon_on_mean = get_nbinom_on_mean


def _depth_label(value: float) -> str:
    value_float = float(value)
    return str(int(value_float)) if value_float.is_integer() else str(value)


def _missing_rate_at_depth(mean_depth: float) -> float:
    """Reuse the existing CFTK depth-dependent missingness model."""
    if mean_depth >= 15:
        _, _, r, p, missing_pct = get_nbinom_on_mean(mean_depth)
    else:
        _, _, r, p, missing_pct = get_nbinom_on_mean(15)
        missing_pct += (
            nbinom.cdf(15 - mean_depth, r, p)
            - nbinom.cdf(0, r, p)
        )
    return float(np.clip(missing_pct, 0.0, 1.0))


def _beta_parameters(
    mean: np.ndarray,
    variance: np.ndarray,
    *,
    beta_eps: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert beta-distribution means and variances to shape parameters."""
    mean = np.asarray(mean, dtype=np.float64)
    variance = np.asarray(variance, dtype=np.float64)

    maximum_variance = mean * (1.0 - mean)
    valid = (
        np.isfinite(mean)
        & np.isfinite(variance)
        & (mean > beta_eps)
        & (mean < 1.0 - beta_eps)
        & (variance > 0)
        & (variance < maximum_variance)
    )

    if not valid.all():
        raise ValueError(
            f"{int((~valid).sum())} mean/variance combinations cannot define "
            "beta distributions."
        )

    concentration = maximum_variance / variance - 1.0
    alpha = mean * concentration
    beta = (1.0 - mean) * concentration

    return alpha, beta


@dataclass(frozen=True)
class ModelPowerTemplate:
    """Cached biological and depth-specific parameters for fast simulation."""

    depths: tuple[float, ...]
    cpg_names: np.ndarray
    baseline_mean: np.ndarray
    raw_effect: np.ndarray
    is_signal: np.ndarray
    sd_by_depth: dict[float, np.ndarray]
    alpha_control: dict[float, np.ndarray]
    beta_control: dict[float, np.ndarray]
    alpha_case: dict[float, np.ndarray]
    beta_case: dict[float, np.ndarray]
    standardized_effect: dict[float, np.ndarray]
    missing_rate: dict[float, float]
    within_block_rho: float
    block_size: int
    beta_eps: float
    metadata: dict[str, Any]


def check_beta_feasibility(
    cpg_std_summary: pd.DataFrame,
    cpg_mean: pd.Series | pd.DataFrame,
    *,
    depth: Sequence[float] = (5, 10, 20, 30),
    sd_stat: str = "mean",
    beta_eps: float = 1e-10,
) -> pd.DataFrame:
    """Report the beta-feasible CpG fraction at each depth and jointly."""
    depths = tuple(depth)
    if not depths:
        raise ValueError("At least one depth is required.")

    if isinstance(cpg_mean, pd.DataFrame):
        if cpg_mean.shape[1] != 1:
            raise ValueError("cpg_mean DataFrame must have exactly one column.")
        mean_series = cpg_mean.iloc[:, 0].copy()
    elif isinstance(cpg_mean, pd.Series):
        mean_series = cpg_mean.copy()
    else:
        raise TypeError("cpg_mean must be a Series or one-column DataFrame.")

    mean_series = pd.to_numeric(mean_series, errors="coerce").rename("mean")
    mean_series.index = mean_series.index.astype(str)

    masks = []
    rows = []

    for mean_depth in depths:
        column = f"{_depth_label(mean_depth)}_{sd_stat}"
        if column not in cpg_std_summary.columns:
            raise KeyError(f"Missing SD column: {column}")

        sd = pd.to_numeric(cpg_std_summary[column], errors="coerce").rename("sd")
        sd.index = sd.index.astype(str)

        temp = pd.concat([mean_series, sd], axis=1, join="inner")
        variance = temp["sd"] ** 2
        max_variance = temp["mean"] * (1.0 - temp["mean"])

        mask = (
            temp["mean"].between(beta_eps, 1.0 - beta_eps, inclusive="neither")
            & np.isfinite(temp["sd"])
            & (temp["sd"] > 0)
            & (variance < max_variance)
        )
        mask.name = mean_depth
        masks.append(mask)

        rows.append(
            {
                "depth": mean_depth,
                "n_aligned_cpgs": len(temp),
                "n_feasible": int(mask.sum()),
                "feasible_fraction": float(mask.mean()),
            }
        )

    feasibility = pd.concat(masks, axis=1, join="inner")
    common = feasibility.all(axis=1)

    result = pd.DataFrame(rows).set_index("depth")
    result["n_common_feasible"] = int(common.sum())
    result["common_feasible_fraction"] = float(common.mean())
    return result


def prepare_model_power_template(
    cpg_std_summary: pd.DataFrame,
    cpg_mean: pd.Series | pd.DataFrame,
    *,
    depth: Sequence[float] = (5, 10, 20, 30),
    n_features: int | None = 500,
    step: int | None = None,
    n_signal_cpgs: int = 20,
    meth_diff: float = 0.05,
    effect_sd: float = 0.015,
    effect_direction: str = "balanced",
    sd_stat: str = "mean",
    within_block_rho: float = 0.0,
    block_size: int = 20,
    beta_eps: float = 1e-10,
    random_state: int | None = 0,
) -> ModelPowerTemplate:
    """
    Prepare and cache all parameters required for repeated power simulations.

    The selected CpGs, signal identities, and CpG-specific effects are fixed
    in the returned template and reused across Monte Carlo replicates.
    """
    rng = np.random.default_rng(random_state)
    depths = tuple(depth)

    if not depths:
        raise ValueError("At least one depth is required.")
    if len(set(depths)) != len(depths):
        raise ValueError("depth contains duplicated values.")
    if n_features is not None and n_features < 1:
        raise ValueError("n_features must be positive or None.")
    if n_features is None and (step is None or step < 1):
        raise ValueError("A positive step is required when n_features is None.")
    if n_signal_cpgs < 0:
        raise ValueError("n_signal_cpgs must be non-negative.")
    if meth_diff < 0 or effect_sd < 0:
        raise ValueError("meth_diff and effect_sd must be non-negative.")
    if n_signal_cpgs > 0 and meth_diff == 0 and effect_sd == 0:
        raise ValueError("Signal CpGs cannot have zero effect magnitude.")
    if effect_direction not in {"positive", "negative", "balanced", "random"}:
        raise ValueError(
            "effect_direction must be 'positive', 'negative', 'balanced', or 'random'."
        )
    if sd_stat not in {"mean", "CI_l", "CI_u"}:
        raise ValueError("sd_stat must be 'mean', 'CI_l', or 'CI_u'.")
    if not 0 <= within_block_rho < 1:
        raise ValueError("within_block_rho must be in [0, 1).")
    if block_size < 1:
        raise ValueError("block_size must be positive.")
    if not 0 < beta_eps < 0.5:
        raise ValueError("beta_eps must be between 0 and 0.5.")

    sd_columns = {
        mean_depth: f"{_depth_label(mean_depth)}_{sd_stat}"
        for mean_depth in depths
    }
    missing_columns = [
        column
        for column in sd_columns.values()
        if column not in cpg_std_summary.columns
    ]
    if missing_columns:
        raise KeyError("Missing SD columns: " + ", ".join(missing_columns))

    if isinstance(cpg_mean, pd.DataFrame):
        if cpg_mean.shape[1] != 1:
            raise ValueError("cpg_mean DataFrame must have exactly one column.")
        mean_series = cpg_mean.iloc[:, 0].copy()
    elif isinstance(cpg_mean, pd.Series):
        mean_series = cpg_mean.copy()
    else:
        raise TypeError("cpg_mean must be a Series or one-column DataFrame.")

    mean_series = pd.to_numeric(
        mean_series, errors="coerce"
    ).rename("baseline_mean")
    sd_table = cpg_std_summary[
        list(sd_columns.values())
    ].apply(pd.to_numeric, errors="coerce")

    mean_series.index = mean_series.index.astype(str)
    sd_table.index = sd_table.index.astype(str)

    if mean_series.index.has_duplicates:
        raise ValueError("cpg_mean contains duplicated CpG indices.")
    if sd_table.index.has_duplicates:
        raise ValueError("cpg_std_summary contains duplicated CpG indices.")

    feature_table = sd_table.join(mean_series, how="inner")
    n_aligned = len(feature_table)

    baseline_all = feature_table["baseline_mean"].to_numpy(dtype=np.float64)
    sd_matrix_all = feature_table[
        list(sd_columns.values())
    ].to_numpy(dtype=np.float64)
    variance_matrix_all = sd_matrix_all**2

    max_control_variance = (
        baseline_all[:, None] * (1.0 - baseline_all[:, None])
    )

    common_control_feasible = (
        np.isfinite(baseline_all)
        & (baseline_all > beta_eps)
        & (baseline_all < 1.0 - beta_eps)
        & np.isfinite(sd_matrix_all).all(axis=1)
        & (sd_matrix_all > 0).all(axis=1)
        & (variance_matrix_all < max_control_variance).all(axis=1)
    )

    feature_table = feature_table.loc[common_control_feasible]
    n_common_feasible = len(feature_table)

    if n_common_feasible == 0:
        raise ValueError("No CpGs are beta-feasible at every requested depth.")

    if n_features is not None:
        if n_features > n_common_feasible:
            raise ValueError(
                f"Requested {n_features} features, but only "
                f"{n_common_feasible} common beta-feasible CpGs are available."
            )
        selected_positions = rng.choice(
            n_common_feasible, size=n_features, replace=False
        )
        selected_positions.sort()
        feature_table = feature_table.iloc[selected_positions]
    else:
        feature_table = feature_table.iloc[::step]

    n_selected = len(feature_table)
    if not 0 <= n_signal_cpgs <= n_selected:
        raise ValueError(
            f"n_signal_cpgs must be between 0 and {n_selected}."
        )

    cpg_names = feature_table.index.to_numpy(dtype=str)
    baseline = feature_table["baseline_mean"].to_numpy(dtype=np.float64)
    selected_sd_matrix = feature_table[
        list(sd_columns.values())
    ].to_numpy(dtype=np.float64)
    selected_variance_matrix = selected_sd_matrix**2

    raw_effect = np.zeros(n_selected, dtype=np.float64)

    if n_signal_cpgs > 0:
        if effect_direction == "positive":
            signs = np.ones(n_signal_cpgs, dtype=np.float64)
        elif effect_direction == "negative":
            signs = -np.ones(n_signal_cpgs, dtype=np.float64)
        elif effect_direction == "balanced":
            n_positive = (n_signal_cpgs + 1) // 2
            signs = np.concatenate(
                [
                    np.ones(n_positive, dtype=np.float64),
                    -np.ones(
                        n_signal_cpgs - n_positive,
                        dtype=np.float64,
                    ),
                ]
            )
            rng.shuffle(signs)
        else:
            signs = rng.choice(
                [-1.0, 1.0], size=n_signal_cpgs
            ).astype(np.float64)

        if effect_sd == 0:
            magnitudes = np.full(
                n_signal_cpgs, meth_diff, dtype=np.float64
            )
        else:
            lower = (beta_eps - meth_diff) / effect_sd
            upper = (1.0 - beta_eps - meth_diff) / effect_sd
            magnitudes = truncnorm.rvs(
                lower,
                upper,
                loc=meth_diff,
                scale=effect_sd,
                size=n_signal_cpgs,
                random_state=rng,
            ).astype(np.float64)

        assignment_order = np.argsort(magnitudes)[::-1]
        used = np.zeros(n_selected, dtype=bool)

        for effect_index in assignment_order:
            signed_effect = signs[effect_index] * magnitudes[effect_index]
            proposed_case_mean = baseline + signed_effect
            proposed_case_max_variance = (
                proposed_case_mean[:, None]
                * (1.0 - proposed_case_mean[:, None])
            )

            feasible = (
                (~used)
                & (proposed_case_mean > beta_eps)
                & (proposed_case_mean < 1.0 - beta_eps)
                & (
                    selected_variance_matrix
                    < proposed_case_max_variance
                ).all(axis=1)
            )

            candidates = np.flatnonzero(feasible)
            if len(candidates) == 0:
                direction_name = "positive" if signed_effect > 0 else "negative"
                raise ValueError(
                    "Unable to assign all signal effects. No unused CpG can "
                    f"support a {direction_name} effect of magnitude "
                    f"{abs(signed_effect):.4f} at every requested depth. "
                    "Reduce meth_diff/effect_sd or n_signal_cpgs, or increase "
                    "n_features."
                )

            position = int(rng.choice(candidates))
            used[position] = True
            raw_effect[position] = signed_effect

    case_mean = baseline + raw_effect
    is_signal = raw_effect != 0

    sd_by_depth: dict[float, np.ndarray] = {}
    alpha_control: dict[float, np.ndarray] = {}
    beta_control: dict[float, np.ndarray] = {}
    alpha_case: dict[float, np.ndarray] = {}
    beta_case: dict[float, np.ndarray] = {}
    standardized_effect: dict[float, np.ndarray] = {}
    missing_rate: dict[float, float] = {}

    for depth_index, mean_depth in enumerate(depths):
        sd = selected_sd_matrix[:, depth_index]
        variance = sd**2

        a0, b0 = _beta_parameters(
            baseline, variance, beta_eps=beta_eps
        )
        a1, b1 = _beta_parameters(
            case_mean, variance, beta_eps=beta_eps
        )

        sd_by_depth[mean_depth] = sd
        alpha_control[mean_depth] = a0
        beta_control[mean_depth] = b0
        alpha_case[mean_depth] = a1
        beta_case[mean_depth] = b1
        standardized_effect[mean_depth] = raw_effect / sd
        missing_rate[mean_depth] = _missing_rate_at_depth(mean_depth)

    signal_effects = np.abs(raw_effect[is_signal])

    return ModelPowerTemplate(
        depths=depths,
        cpg_names=cpg_names,
        baseline_mean=baseline,
        raw_effect=raw_effect,
        is_signal=is_signal,
        sd_by_depth=sd_by_depth,
        alpha_control=alpha_control,
        beta_control=beta_control,
        alpha_case=alpha_case,
        beta_case=beta_case,
        standardized_effect=standardized_effect,
        missing_rate=missing_rate,
        within_block_rho=within_block_rho,
        block_size=block_size,
        beta_eps=beta_eps,
        metadata={
            "n_aligned_panel_cpgs": n_aligned,
            "n_common_beta_feasible_cpgs": n_common_feasible,
            "common_beta_feasible_fraction": (
                n_common_feasible / n_aligned if n_aligned else np.nan
            ),
            "n_features": n_selected,
            "n_signal_cpgs": n_signal_cpgs,
            "requested_mean_abs_effect": meth_diff,
            "requested_effect_sd": effect_sd,
            "realized_mean_abs_effect": (
                float(signal_effects.mean()) if len(signal_effects) else 0.0
            ),
            "realized_effect_sd": (
                float(signal_effects.std(ddof=1))
                if len(signal_effects) > 1
                else 0.0
            ),
            "effect_direction": effect_direction,
            "sd_stat": sd_stat,
            "random_state": random_state,
        },
    )


def _make_labels(
    total_sample_size: int,
    ratio: float,
    rng: np.random.Generator,
) -> np.ndarray:
    if total_sample_size < 4:
        raise ValueError("total_sample_size must be at least 4.")
    if ratio <= 0:
        raise ValueError("ratio must be positive.")

    n_cases = int(round(total_sample_size * ratio / (1.0 + ratio)))
    n_controls = total_sample_size - n_cases

    if min(n_cases, n_controls) < 2:
        raise ValueError("At least two cases and two controls are required.")

    y = np.concatenate(
        [
            np.zeros(n_controls, dtype=np.int8),
            np.ones(n_cases, dtype=np.int8),
        ]
    )
    rng.shuffle(y)
    return y


def _draw_copula_uniform(
    n_samples: int,
    n_features: int,
    *,
    rho: float,
    block_size: int,
    rng: np.random.Generator,
    eps: float,
) -> np.ndarray:
    if rho == 0:
        return np.clip(
            rng.random((n_samples, n_features)),
            eps,
            1.0 - eps,
        )

    latent_z = np.empty((n_samples, n_features), dtype=np.float64)

    for start in range(0, n_features, block_size):
        stop = min(start + block_size, n_features)
        width = stop - start

        shared = rng.normal(size=(n_samples, 1))
        independent = rng.normal(size=(n_samples, width))

        latent_z[:, start:stop] = (
            np.sqrt(rho) * shared
            + np.sqrt(1.0 - rho) * independent
        )

    return np.clip(ndtr(latent_z), eps, 1.0 - eps)


def simulate_from_model_power_template(
    template: ModelPowerTemplate,
    *,
    total_sample_size: int,
    ratio: float = 1.0,
    paired_depths: bool = True,
    dtype: np.dtype = np.float32,
    random_state: int | None = None,
) -> tuple[dict[float, np.ndarray], np.ndarray]:
    """
    Generate one cohort for all depths from a cached template.

    paired_depths=True reuses copula quantiles and missingness uniforms across
    depths, reducing Monte Carlo noise in depth comparisons.
    """
    rng = np.random.default_rng(random_state)
    y = _make_labels(total_sample_size, ratio, rng)

    n_samples = len(y)
    n_features = len(template.cpg_names)
    control_mask = y == 0
    case_mask = y == 1

    shared_uniform = None
    shared_missing = None

    if paired_depths:
        shared_uniform = _draw_copula_uniform(
            n_samples,
            n_features,
            rho=template.within_block_rho,
            block_size=template.block_size,
            rng=rng,
            eps=max(template.beta_eps, 1e-7),
        )
        shared_missing = rng.random((n_samples, n_features))

    X_by_depth: dict[float, np.ndarray] = {}

    for mean_depth in template.depths:
        # Direct NumPy beta sampling is faster when correlation and paired
        # depth quantiles are not requested.
        if not paired_depths and template.within_block_rho == 0:
            X = np.empty((n_samples, n_features), dtype=np.float64)
            X[control_mask] = rng.beta(
                template.alpha_control[mean_depth],
                template.beta_control[mean_depth],
                size=(int(control_mask.sum()), n_features),
            )
            X[case_mask] = rng.beta(
                template.alpha_case[mean_depth],
                template.beta_case[mean_depth],
                size=(int(case_mask.sum()), n_features),
            )
            missing_uniform = rng.random((n_samples, n_features))
        else:
            copula_uniform = (
                shared_uniform
                if paired_depths
                else _draw_copula_uniform(
                    n_samples,
                    n_features,
                    rho=template.within_block_rho,
                    block_size=template.block_size,
                    rng=rng,
                    eps=max(template.beta_eps, 1e-7),
                )
            )

            X = np.empty((n_samples, n_features), dtype=np.float64)
            X[control_mask] = betaincinv(
                template.alpha_control[mean_depth][None, :],
                template.beta_control[mean_depth][None, :],
                copula_uniform[control_mask],
            )
            X[case_mask] = betaincinv(
                template.alpha_case[mean_depth][None, :],
                template.beta_case[mean_depth][None, :],
                copula_uniform[case_mask],
            )
            missing_uniform = (
                shared_missing
                if paired_depths
                else rng.random((n_samples, n_features))
            )

        if not np.isfinite(X).all():
            n_bad = int((~np.isfinite(X)).sum())
            raise RuntimeError(
                f"{n_bad} non-finite beta values were generated at depth "
                f"{mean_depth}."
            )

        # Exact 0 and 1 are valid beta values; clipping only protects against
        # numerical excursions outside the valid range.
        X = np.clip(X, 0.0, 1.0).astype(dtype, copy=False)
        X[missing_uniform < template.missing_rate[mean_depth]] = np.nan
        X_by_depth[mean_depth] = X

    return X_by_depth, y


def simulate_train_validation_from_template(
    template: ModelPowerTemplate,
    *,
    n_train: int,
    n_validation: int,
    ratio: float = 1.0,
    paired_depths: bool = True,
    dtype: np.dtype = np.float32,
    random_state: int | None = None,
) -> dict[str, Any]:
    """Generate independent training and validation cohorts."""
    seed_sequence = np.random.SeedSequence(random_state)
    train_child, validation_child = seed_sequence.spawn(2)

    train_seed = int(train_child.generate_state(1, dtype=np.uint32)[0])
    validation_seed = int(
        validation_child.generate_state(1, dtype=np.uint32)[0]
    )

    X_train, y_train = simulate_from_model_power_template(
        template,
        total_sample_size=n_train,
        ratio=ratio,
        paired_depths=paired_depths,
        dtype=dtype,
        random_state=train_seed,
    )
    X_validation, y_validation = simulate_from_model_power_template(
        template,
        total_sample_size=n_validation,
        ratio=ratio,
        paired_depths=paired_depths,
        dtype=dtype,
        random_state=validation_seed,
    )

    return {
        "X_train_by_depth": X_train,
        "y_train": y_train,
        "X_validation_by_depth": X_validation,
        "y_validation": y_validation,
    }


def _prepare_features_for_model(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_validation: np.ndarray,
    is_signal: np.ndarray,
    *,
    min_observed_fraction: float,
    top_k: int | None,
) -> dict[str, Any]:
    """Training-only filtering, median imputation, and feature selection."""
    original_indices = np.arange(X_train.shape[1])

    observed_control = np.mean(~np.isnan(X_train[y_train == 0]), axis=0)
    observed_case = np.mean(~np.isnan(X_train[y_train == 1]), axis=0)

    keep_observed = (
        (observed_control >= min_observed_fraction)
        & (observed_case >= min_observed_fraction)
    )

    X_train = X_train[:, keep_observed]
    X_validation = X_validation[:, keep_observed]
    feature_indices = original_indices[keep_observed]

    if X_train.shape[1] == 0:
        raise RuntimeError("No CpGs pass the training observation filter.")

    train_medians = np.nanmedian(X_train, axis=0)
    keep_median = np.isfinite(train_medians)

    X_train = X_train[:, keep_median]
    X_validation = X_validation[:, keep_median]
    train_medians = train_medians[keep_median]
    feature_indices = feature_indices[keep_median]

    if X_train.shape[1] == 0:
        raise RuntimeError("No CpGs remain after median-imputation checks.")

    X_train = np.where(
        np.isnan(X_train), train_medians[None, :], X_train
    )
    X_validation = np.where(
        np.isnan(X_validation), train_medians[None, :], X_validation
    )

    training_variance = np.var(X_train, axis=0)
    keep_variable = np.isfinite(training_variance) & (training_variance > 0)

    X_train = X_train[:, keep_variable]
    X_validation = X_validation[:, keep_variable]
    feature_indices = feature_indices[keep_variable]

    if X_train.shape[1] == 0:
        raise RuntimeError("No variable CpGs remain after preprocessing.")

    n_before_selection = X_train.shape[1]

    if top_k is not None and top_k < n_before_selection:
        scores, _ = f_classif(X_train, y_train)
        scores = np.nan_to_num(
            scores,
            nan=-np.inf,
            posinf=np.finfo(np.float64).max,
            neginf=-np.inf,
        )
        selected_local = np.argsort(scores)[::-1][:top_k]
        X_train = X_train[:, selected_local]
        X_validation = X_validation[:, selected_local]
        feature_indices = feature_indices[selected_local]

    selected_signal = is_signal[feature_indices]
    n_selected_signal = int(selected_signal.sum())
    total_signal = int(is_signal.sum())

    return {
        "X_train": X_train,
        "X_validation": X_validation,
        "selected_feature_indices": feature_indices,
        "n_features_before_selection": n_before_selection,
        "n_features_selected": len(feature_indices),
        "n_selected_signal": n_selected_signal,
        "feature_recall": (
            n_selected_signal / total_signal if total_signal else np.nan
        ),
        "feature_precision": (
            n_selected_signal / len(feature_indices)
            if len(feature_indices)
            else np.nan
        ),
    }


def _auc_pvalue(y_true: np.ndarray, score: np.ndarray) -> float:
    """One-sided independent-validation rank test for AUC > 0.5."""
    return float(
        mannwhitneyu(
            score[y_true == 1],
            score[y_true == 0],
            alternative="greater",
            method="asymptotic",
        ).pvalue
    )


def _sensitivity_at_specificity(
    y_true: np.ndarray,
    score: np.ndarray,
    target_specificity: float,
) -> float:
    fpr, tpr, _ = roc_curve(y_true, score)
    allowed = fpr <= (1.0 - target_specificity + 1e-12)
    return float(np.max(tpr[allowed])) if np.any(allowed) else 0.0


def _fit_and_score_model(
    model_name: str,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_validation: np.ndarray,
    y_validation: np.ndarray,
    *,
    model_seed: int,
    logreg_params: Mapping[str, Any],
    rf_params: Mapping[str, Any],
    target_auc: float,
    alpha: float,
    specificity_target: float,
) -> dict[str, Any]:
    if model_name == "logreg":
        scaler = StandardScaler()
        X_train_model = scaler.fit_transform(X_train)
        X_validation_model = scaler.transform(X_validation)

        params = dict(logreg_params)
        params.setdefault("random_state", model_seed)
        estimator = LogisticRegression(**params)

    elif model_name == "rf":
        X_train_model = X_train
        X_validation_model = X_validation

        params = dict(rf_params)
        params.setdefault("random_state", model_seed)
        estimator = RandomForestClassifier(**params)

    else:
        raise ValueError(f"Unknown model: {model_name}")

    estimator.fit(X_train_model, y_train)
    score = estimator.predict_proba(X_validation_model)[:, 1]

    auc = float(roc_auc_score(y_validation, score))
    pvalue = _auc_pvalue(y_validation, score)
    sensitivity = _sensitivity_at_specificity(
        y_validation, score, specificity_target
    )
    brier = float(brier_score_loss(y_validation, score))

    detection_success = auc > 0.5 and pvalue < alpha
    target_success = auc >= target_auc

    return {
        "auc": auc,
        "auc_pvalue": pvalue,
        "sensitivity_at_specificity": sensitivity,
        "brier_score": brier,
        "detection_success": detection_success,
        "target_success": target_success,
        "joint_success": detection_success and target_success,
    }


def _run_one_power_replicate(
    replicate_index: int,
    replicate_seed: int,
    template: ModelPowerTemplate,
    *,
    n_train: int,
    n_validation: int,
    ratio: float,
    models: tuple[str, ...],
    top_k: int | None,
    min_observed_fraction: float,
    target_auc: float,
    alpha: float,
    specificity_target: float,
    paired_depths: bool,
    logreg_params: Mapping[str, Any],
    rf_params: Mapping[str, Any],
) -> list[dict[str, Any]]:
    simulated = simulate_train_validation_from_template(
        template,
        n_train=n_train,
        n_validation=n_validation,
        ratio=ratio,
        paired_depths=paired_depths,
        random_state=replicate_seed,
    )

    y_train = simulated["y_train"]
    y_validation = simulated["y_validation"]
    rows: list[dict[str, Any]] = []

    for depth_index, mean_depth in enumerate(template.depths):
        X_train = simulated["X_train_by_depth"][mean_depth]
        X_validation = simulated["X_validation_by_depth"][mean_depth]

        prepared = _prepare_features_for_model(
            X_train,
            y_train,
            X_validation,
            template.is_signal,
            min_observed_fraction=min_observed_fraction,
            top_k=top_k,
        )

        for model_index, model_name in enumerate(models):
            model_seed = int(
                np.random.SeedSequence(
                    [replicate_seed, depth_index, model_index]
                ).generate_state(1, dtype=np.uint32)[0]
            )

            metrics = _fit_and_score_model(
                model_name,
                prepared["X_train"],
                y_train,
                prepared["X_validation"],
                y_validation,
                model_seed=model_seed,
                logreg_params=logreg_params,
                rf_params=rf_params,
                target_auc=target_auc,
                alpha=alpha,
                specificity_target=specificity_target,
            )

            rows.append(
                {
                    "simulation": replicate_index,
                    "simulation_seed": replicate_seed,
                    "depth": mean_depth,
                    "model": model_name,
                    "n_train": n_train,
                    "n_validation": n_validation,
                    "n_train_cases": int(np.sum(y_train == 1)),
                    "n_train_controls": int(np.sum(y_train == 0)),
                    "n_validation_cases": int(
                        np.sum(y_validation == 1)
                    ),
                    "n_validation_controls": int(
                        np.sum(y_validation == 0)
                    ),
                    "n_features_before_filter": len(
                        template.cpg_names
                    ),
                    "n_features_before_selection": prepared[
                        "n_features_before_selection"
                    ],
                    "n_features_selected": prepared[
                        "n_features_selected"
                    ],
                    "n_selected_signal": prepared[
                        "n_selected_signal"
                    ],
                    "feature_recall": prepared["feature_recall"],
                    "feature_precision": prepared[
                        "feature_precision"
                    ],
                    "specificity_target": specificity_target,
                    **metrics,
                }
            )

    return rows


def _summarize_power_results(
    replicate_results: pd.DataFrame,
    *,
    confidence: float,
) -> pd.DataFrame:
    lower_q = (1.0 - confidence) / 2.0
    upper_q = 1.0 - lower_q
    confidence_alpha = 1.0 - confidence
    rows = []

    for (model_name, mean_depth), group in replicate_results.groupby(
        ["model", "depth"], sort=False
    ):
        row = {
            "model": model_name,
            "depth": mean_depth,
            "n_simulations": len(group),
            "mean_auc": group["auc"].mean(),
            "sd_auc": group["auc"].std(ddof=1),
            "median_auc": group["auc"].median(),
            "auc_interval_low": group["auc"].quantile(lower_q),
            "auc_interval_high": group["auc"].quantile(upper_q),
            "mean_sensitivity_at_specificity": group[
                "sensitivity_at_specificity"
            ].mean(),
            "mean_brier_score": group["brier_score"].mean(),
            "mean_feature_recall": group["feature_recall"].mean(),
            "mean_feature_precision": group[
                "feature_precision"
            ].mean(),
        }

        for success_column in (
            "detection_success",
            "target_success",
            "joint_success",
        ):
            successes = int(group[success_column].sum())
            total = len(group)
            low, high = proportion_confint(
                successes,
                total,
                alpha=confidence_alpha,
                method="wilson",
            )
            power_name = success_column.replace("_success", "_power")
            row[power_name] = successes / total
            row[f"{power_name}_ci_low"] = low
            row[f"{power_name}_ci_high"] = high

        rows.append(row)

    return (
        pd.DataFrame(rows)
        .sort_values(["model", "depth"])
        .reset_index(drop=True)
    )


def _adaptive_stop_reached(
    replicate_results: pd.DataFrame,
    *,
    success_column: str,
    target_ci_width: float,
    confidence: float,
    minimum_simulations: int,
) -> bool:
    confidence_alpha = 1.0 - confidence

    for _, group in replicate_results.groupby(["model", "depth"]):
        total = len(group)
        if total < minimum_simulations:
            return False

        successes = int(group[success_column].sum())
        low, high = proportion_confint(
            successes,
            total,
            alpha=confidence_alpha,
            method="wilson",
        )
        if (high - low) > target_ci_width:
            return False

    return True


def _paired_model_comparison(
    replicate_results: pd.DataFrame,
) -> pd.DataFrame:
    available = set(replicate_results["model"].unique())
    if not {"logreg", "rf"}.issubset(available):
        return pd.DataFrame()

    paired = replicate_results.pivot(
        index=["simulation", "depth"],
        columns="model",
        values="auc",
    ).dropna(subset=["logreg", "rf"])

    paired["rf_minus_logreg_auc"] = paired["rf"] - paired["logreg"]

    return (
        paired.reset_index()
        .groupby("depth")
        .agg(
            n_simulations=("rf_minus_logreg_auc", "size"),
            mean_rf_minus_logreg_auc=(
                "rf_minus_logreg_auc",
                "mean",
            ),
            median_rf_minus_logreg_auc=(
                "rf_minus_logreg_auc",
                "median",
            ),
            q025_rf_minus_logreg_auc=(
                "rf_minus_logreg_auc",
                lambda x: x.quantile(0.025),
            ),
            q975_rf_minus_logreg_auc=(
                "rf_minus_logreg_auc",
                lambda x: x.quantile(0.975),
            ),
            probability_rf_better=(
                "rf_minus_logreg_auc",
                lambda x: float(np.mean(x > 0)),
            ),
            probability_rf_improves_003=(
                "rf_minus_logreg_auc",
                lambda x: float(np.mean(x >= 0.03)),
            ),
            probability_rf_improves_005=(
                "rf_minus_logreg_auc",
                lambda x: float(np.mean(x >= 0.05)),
            ),
        )
        .reset_index()
    )


def run_fast_model_power_analysis(
    template: ModelPowerTemplate,
    *,
    n_train: int = 100,
    n_validation: int = 100,
    ratio: float = 1.0,
    models: Sequence[str] = ("logreg",),
    top_k: int | None = 50,
    min_observed_fraction: float = 0.5,
    target_auc: float = 0.75,
    alpha: float = 0.05,
    specificity_target: float = 0.90,
    initial_simulations: int = 50,
    batch_size: int = 25,
    max_simulations: int = 200,
    target_ci_width: float | None = 0.15,
    stop_metric: str = "joint_success",
    monte_carlo_ci: float = 0.95,
    paired_depths: bool = True,
    logreg_params: Mapping[str, Any] | None = None,
    rf_params: Mapping[str, Any] | None = None,
    n_jobs: int = 1,
    parallel_prefer: str = "threads",
    verbose: int = 0,
    random_state: int = 0,
) -> dict[str, Any]:
    """
    Run adaptive, web-oriented model-level power analysis.

    The primary power result is joint_power:
        P(validation AUC >= target_auc and one-sided AUC p-value < alpha).
    """
    models = tuple(models)
    allowed_models = {"logreg", "rf"}

    if not models:
        raise ValueError("At least one model is required.")
    if set(models) - allowed_models:
        raise ValueError(
            "models may contain only 'logreg' and 'rf'."
        )
    if n_train < 4 or n_validation < 4:
        raise ValueError("n_train and n_validation must each be at least 4.")
    if ratio <= 0:
        raise ValueError("ratio must be positive.")
    if top_k is not None and top_k < 1:
        raise ValueError("top_k must be positive or None.")
    if not 0 < min_observed_fraction <= 1:
        raise ValueError("min_observed_fraction must be in (0, 1].")
    if not 0.5 <= target_auc <= 1:
        raise ValueError("target_auc must be between 0.5 and 1.")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1.")
    if not 0 < specificity_target < 1:
        raise ValueError("specificity_target must be between 0 and 1.")
    if initial_simulations < 1 or batch_size < 1:
        raise ValueError("initial_simulations and batch_size must be positive.")
    if max_simulations < initial_simulations:
        raise ValueError(
            "max_simulations must be >= initial_simulations."
        )
    if target_ci_width is not None and not 0 < target_ci_width < 1:
        raise ValueError("target_ci_width must be in (0, 1) or None.")
    if stop_metric not in {
        "detection_success",
        "target_success",
        "joint_success",
    }:
        raise ValueError(
            "stop_metric must be detection_success, target_success, "
            "or joint_success."
        )
    if not 0 < monte_carlo_ci < 1:
        raise ValueError("monte_carlo_ci must be between 0 and 1.")
    if parallel_prefer not in {"threads", "processes"}:
        raise ValueError(
            "parallel_prefer must be 'threads' or 'processes'."
        )

    logreg_defaults = {
        "C": 1.0,
        "solver": "liblinear",
        "max_iter": 2000,
        "class_weight": None,
    }
    logreg_defaults.update(dict(logreg_params or {}))

    rf_defaults = {
        "n_estimators": 100,
        "max_features": "sqrt",
        "min_samples_leaf": 5,
        "class_weight": None,
        "n_jobs": 1,
    }
    rf_defaults.update(dict(rf_params or {}))

    seed_sequence = np.random.SeedSequence(random_state)
    child_sequences = seed_sequence.spawn(max_simulations)
    replicate_seeds = [
        int(child.generate_state(1, dtype=np.uint32)[0])
        for child in child_sequences
    ]

    all_rows: list[dict[str, Any]] = []
    n_completed = 0
    stopped_early = False
    start_time = time.perf_counter()

    while n_completed < max_simulations:
        requested_batch = (
            initial_simulations if n_completed == 0 else batch_size
        )
        current_batch = min(
            requested_batch, max_simulations - n_completed
        )

        indices = range(n_completed, n_completed + current_batch)

        batch_results = Parallel(
            n_jobs=n_jobs,
            prefer=parallel_prefer,
            verbose=verbose,
        )(
            delayed(_run_one_power_replicate)(
                replicate_index=index,
                replicate_seed=replicate_seeds[index],
                template=template,
                n_train=n_train,
                n_validation=n_validation,
                ratio=ratio,
                models=models,
                top_k=top_k,
                min_observed_fraction=min_observed_fraction,
                target_auc=target_auc,
                alpha=alpha,
                specificity_target=specificity_target,
                paired_depths=paired_depths,
                logreg_params=logreg_defaults,
                rf_params=rf_defaults,
            )
            for index in indices
        )

        all_rows.extend(
            row
            for replicate_rows in batch_results
            for row in replicate_rows
        )
        n_completed += current_batch

        if target_ci_width is not None:
            current_results = pd.DataFrame(all_rows)
            if _adaptive_stop_reached(
                current_results,
                success_column=stop_metric,
                target_ci_width=target_ci_width,
                confidence=monte_carlo_ci,
                minimum_simulations=initial_simulations,
            ):
                stopped_early = n_completed < max_simulations
                break

    elapsed_seconds = time.perf_counter() - start_time
    replicate_results = pd.DataFrame(all_rows)

    power_summary = _summarize_power_results(
        replicate_results,
        confidence=monte_carlo_ci,
    )

    return {
        "replicate_results": replicate_results,
        "power_summary": power_summary,
        "model_comparison": _paired_model_comparison(
            replicate_results
        ),
        "analysis_metadata": {
            "n_simulations_used": n_completed,
            "stopped_early": stopped_early,
            "elapsed_seconds": elapsed_seconds,
            "target_ci_width": target_ci_width,
            "stop_metric": stop_metric,
            "paired_depths": paired_depths,
            "n_jobs": n_jobs,
            "parallel_prefer": parallel_prefer,
        },
    }



def prepare_template_ensemble(
    cpg_std_summary,
    cpg_mean,
    *,
    n_templates=10,
    template_kwargs=None,
    random_state=1,
):
    """
    Prepare multiple independent CpG/signal templates once.

    Reuse these same templates across all sample sizes to obtain paired,
    comparable power curves.
    """
    template_kwargs = dict(template_kwargs or {})

    forbidden = {
        "cpg_std_summary",
        "cpg_mean",
        "random_state",
    }
    conflicts = forbidden.intersection(template_kwargs)

    if conflicts:
        raise ValueError(
            "Remove these keys from template_kwargs: "
            + ", ".join(sorted(conflicts))
        )

    seed_sequence = np.random.SeedSequence(random_state)
    child_sequences = seed_sequence.spawn(n_templates)

    templates = []

    for template_id, child in enumerate(child_sequences):
        template_seed = int(
            child.generate_state(1, dtype=np.uint32)[0]
        )

        template = prepare_model_power_template(
            cpg_std_summary=cpg_std_summary,
            cpg_mean=cpg_mean,
            random_state=template_seed,
            **template_kwargs,
        )

        templates.append(
            {
                "template_id": template_id,
                "template_seed": template_seed,
                "template": template,
            }
        )

    return templates

def run_power_sample_size_grid(
    templates,
    sample_sizes,
    *,
    train_fraction=0.5,
    simulations_per_template=20,
    power_kwargs=None,
    ci_method="none",
    ci_metrics=("probability_of_success",),
    confidence=0.95,
    n_bootstrap=500,
    n_jobs=8,
    random_state=1000,
):
    """
    Run model-level power analysis across total study sample sizes.

    The same pre-generated CpG templates are reused for every sample size.
    This prevents the power curves from being confounded by changes in the
    sampled CpGs or signal architecture.

    Parameters
    ----------
    templates : sequence
        Output from prepare_template_ensemble(). Each element must contain:

            {
                "template_id": ...,
                "template_seed": ...,
                "template": ModelPowerTemplate,
            }

    sample_sizes : sequence of int
        Total study sample sizes. Each total is divided into training and
        independent validation cohorts according to `train_fraction`.

    train_fraction : float, default=0.5
        Fraction of the total sample size allocated to model training.

    simulations_per_template : int, default=20
        Exact number of independently simulated train/validation cohorts for
        every template × sample-size combination.

    power_kwargs : dict or None
        Arguments passed to run_fast_model_power_analysis(), excluding:

            template
            n_train
            n_validation
            initial_simulations
            max_simulations
            target_ci_width
            n_jobs
            random_state

    ci_method : {"none", "pooled_wilson", "hierarchical_bootstrap"}
        Method used to calculate uncertainty intervals.

        - "none":
            No CI calculation. Fastest option.

        - "pooled_wilson":
            Fast Wilson interval treating cohort simulations as independent.
            This does not fully account for shared CpG templates.

        - "hierarchical_bootstrap":
            Resamples templates and then cohort simulations within templates.
            Preferred for final analyses.

    ci_metrics : sequence of str
        Power quantities for which CIs are calculated. Allowed values:

            "probability_of_success"
            "detection_power"
            "target_attainment_probability"

        For speed, the default calculates a CI only for the primary endpoint,
        probability_of_success.

    confidence : float, default=0.95
        Confidence level.

    n_bootstrap : int, default=500
        Number of hierarchical bootstrap replicates. Used only when
        ci_method="hierarchical_bootstrap".

    n_jobs : int, default=8
        Number of workers passed to run_fast_model_power_analysis().

    random_state : int, default=1000
        Master random seed.

    Returns
    -------
    dict
        replicate_results :
            One row per template × sample size × simulation × depth × model.

        template_summary :
            Power estimates within each CpG template.

        power_curve :
            Power estimates averaged equally across templates, ready for
            plotting against sample size.

        run_metadata :
            Analysis configuration.
    """
    from statsmodels.stats.proportion import proportion_confint

    # Assumes this function is available in the current namespace or imported
    # from model_power.
    # from model_power import run_fast_model_power_analysis

    power_kwargs = dict(power_kwargs or {})
    sample_sizes = [int(value) for value in sample_sizes]
    ci_metrics = tuple(ci_metrics)

    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    if not templates:
        raise ValueError("templates must contain at least one template.")

    if not sample_sizes:
        raise ValueError("sample_sizes must contain at least one value.")

    if any(value < 8 for value in sample_sizes):
        raise ValueError(
            "Every total sample size must be at least 8 so that both "
            "training and validation cohorts can contain at least 4 subjects."
        )

    if len(set(sample_sizes)) != len(sample_sizes):
        raise ValueError("sample_sizes contains duplicated values.")

    if not 0 < train_fraction < 1:
        raise ValueError("train_fraction must be between 0 and 1.")

    if simulations_per_template < 1:
        raise ValueError("simulations_per_template must be positive.")

    if ci_method not in {
        "none",
        "pooled_wilson",
        "hierarchical_bootstrap",
    }:
        raise ValueError(
            "ci_method must be 'none', 'pooled_wilson', "
            "or 'hierarchical_bootstrap'."
        )

    allowed_ci_metrics = {
        "probability_of_success",
        "detection_power",
        "target_attainment_probability",
    }

    unknown_ci_metrics = set(ci_metrics) - allowed_ci_metrics

    if unknown_ci_metrics:
        raise ValueError(
            "Unknown ci_metrics: "
            + ", ".join(sorted(unknown_ci_metrics))
        )

    if not 0 < confidence < 1:
        raise ValueError("confidence must be between 0 and 1.")

    if (
        ci_method == "hierarchical_bootstrap"
        and n_bootstrap < 100
    ):
        raise ValueError(
            "Use at least 100 bootstrap replicates for "
            "hierarchical_bootstrap."
        )

    forbidden_power_kwargs = {
        "template",
        "n_train",
        "n_validation",
        "initial_simulations",
        "max_simulations",
        "target_ci_width",
        "n_jobs",
        "random_state",
    }

    conflicts = forbidden_power_kwargs.intersection(power_kwargs)

    if conflicts:
        raise ValueError(
            "Remove these keys from power_kwargs because they are "
            "controlled by run_power_sample_size_grid: "
            + ", ".join(sorted(conflicts))
        )

    required_template_keys = {
        "template_id",
        "template_seed",
        "template",
    }

    for index, record in enumerate(templates):
        missing_keys = required_template_keys - set(record)

        if missing_keys:
            raise ValueError(
                f"Template record {index} is missing keys: "
                + ", ".join(sorted(missing_keys))
            )

    # ------------------------------------------------------------------
    # 2. CI helper
    # ------------------------------------------------------------------
    def _calculate_binary_ci(
        group,
        *,
        success_column,
        seed,
    ):
        if ci_method == "none":
            return np.nan, np.nan

        values = group[
            success_column
        ].to_numpy(dtype=bool)

        if ci_method == "pooled_wilson":
            low, high = proportion_confint(
                count=int(values.sum()),
                nobs=len(values),
                alpha=1.0 - confidence,
                method="wilson",
            )

            return float(low), float(high)

        # Hierarchical bootstrap:
        # template sampling followed by cohort sampling within template.
        rng = np.random.default_rng(seed)

        values_by_template = {
            template_id: template_group[
                success_column
            ].to_numpy(dtype=float)
            for template_id, template_group
            in group.groupby("template_id")
        }

        template_ids = np.asarray(
            list(values_by_template)
        )

        if len(template_ids) < 2:
            return np.nan, np.nan

        bootstrap_estimates = np.empty(
            n_bootstrap,
            dtype=float,
        )

        for bootstrap_index in range(n_bootstrap):
            sampled_template_ids = rng.choice(
                template_ids,
                size=len(template_ids),
                replace=True,
            )

            sampled_template_means = np.empty(
                len(sampled_template_ids),
                dtype=float,
            )

            for position, template_id in enumerate(
                sampled_template_ids
            ):
                template_values = values_by_template[
                    template_id
                ]

                sampled_cohorts = rng.choice(
                    template_values,
                    size=len(template_values),
                    replace=True,
                )

                sampled_template_means[position] = (
                    sampled_cohorts.mean()
                )

            # Equal weighting across CpG templates.
            bootstrap_estimates[bootstrap_index] = (
                sampled_template_means.mean()
            )

        tail_probability = (1.0 - confidence) / 2.0

        return (
            float(
                np.quantile(
                    bootstrap_estimates,
                    tail_probability,
                )
            ),
            float(
                np.quantile(
                    bootstrap_estimates,
                    1.0 - tail_probability,
                )
            ),
        )

    # ------------------------------------------------------------------
    # 3. Run exact Monte Carlo simulations
    # ------------------------------------------------------------------
    all_replicates = []

    for sample_size_index, total_sample_size in enumerate(
        sample_sizes
    ):
        n_train = int(
            round(total_sample_size * train_fraction)
        )
        n_validation = total_sample_size - n_train

        if min(n_train, n_validation) < 4:
            raise ValueError(
                f"sample_size={total_sample_size} gives "
                f"n_train={n_train} and "
                f"n_validation={n_validation}. "
                "Both must be at least 4."
            )

        for template_record in templates:
            template_id = template_record["template_id"]
            template_seed = template_record["template_seed"]
            template = template_record["template"]

            run_seed = int(
                np.random.SeedSequence(
                    [
                        random_state,
                        sample_size_index,
                        int(template_id),
                    ]
                ).generate_state(
                    1,
                    dtype=np.uint32,
                )[0]
            )

            result = run_fast_model_power_analysis(
                template,
                n_train=n_train,
                n_validation=n_validation,

                # Disable adaptive stopping so every template has
                # exactly equal Monte Carlo weight.
                initial_simulations=simulations_per_template,
                max_simulations=simulations_per_template,
                target_ci_width=None,

                n_jobs=n_jobs,
                random_state=run_seed,
                **power_kwargs,
            )

            replicate = result[
                "replicate_results"
            ].copy()

            replicate.insert(
                0,
                "template_id",
                template_id,
            )
            replicate.insert(
                1,
                "template_seed",
                template_seed,
            )

            replicate["sample_size"] = total_sample_size
            replicate["n_train"] = n_train
            replicate["n_validation"] = n_validation

            # User-facing alias for plotting.
            replicate["mean_depth"] = replicate["depth"]

            all_replicates.append(replicate)

    replicate_results = pd.concat(
        all_replicates,
        ignore_index=True,
    )

    # ------------------------------------------------------------------
    # 4. Calculate template-specific estimates
    # ------------------------------------------------------------------
    template_summary = (
        replicate_results
        .groupby(
            [
                "template_id",
                "sample_size",
                "n_train",
                "n_validation",
                "model",
                "mean_depth",
            ],
            as_index=False,
        )
        .agg(
            n_simulations=("joint_success", "size"),
            probability_of_success=("joint_success", "mean"),
            detection_power=("detection_success", "mean"),
            target_attainment_probability=(
                "target_success",
                "mean",
            ),
            mean_auc=("auc", "mean"),
            sd_auc=("auc", "std"),
            mean_feature_recall=("feature_recall", "mean"),
            mean_feature_precision=(
                "feature_precision",
                "mean",
            ),
        )
    )

    # ------------------------------------------------------------------
    # 5. Average equally across CpG templates
    # ------------------------------------------------------------------
    overall_rows = []

    grouped_replicates = replicate_results.groupby(
        [
            "sample_size",
            "n_train",
            "n_validation",
            "model",
            "mean_depth",
        ],
        sort=True,
    )

    for group_key, replicate_group in grouped_replicates:
        (
            total_sample_size,
            n_train,
            n_validation,
            model_name,
            mean_depth,
        ) = group_key

        template_group = template_summary.loc[
            (template_summary["sample_size"] == total_sample_size)
            & (template_summary["n_train"] == n_train)
            & (
                template_summary["n_validation"]
                == n_validation
            )
            & (template_summary["model"] == model_name)
            & (
                template_summary["mean_depth"]
                == mean_depth
            )
        ]

        # Point estimates are averaged equally across templates.
        probability_of_success = float(
            template_group[
                "probability_of_success"
            ].mean()
        )

        detection_power = float(
            template_group[
                "detection_power"
            ].mean()
        )

        target_attainment_probability = float(
            template_group[
                "target_attainment_probability"
            ].mean()
        )

        mean_auc = float(
            template_group["mean_auc"].mean()
        )

        ci_results = {
            "probability_of_success": (np.nan, np.nan),
            "detection_power": (np.nan, np.nan),
            "target_attainment_probability": (
                np.nan,
                np.nan,
            ),
        }

        source_columns = {
            "probability_of_success": "joint_success",
            "detection_power": "detection_success",
            "target_attainment_probability": (
                "target_success"
            ),
        }

        if ci_method != "none":
            for metric_index, metric_name in enumerate(
                ci_metrics
            ):
                ci_seed = int(
                    np.random.SeedSequence(
                        [
                            random_state,
                            int(total_sample_size),
                            int(float(mean_depth) * 100),
                            metric_index,
                            0 if model_name == "logreg" else 1,
                        ]
                    ).generate_state(
                        1,
                        dtype=np.uint32,
                    )[0]
                )

                ci_results[metric_name] = (
                    _calculate_binary_ci(
                        replicate_group,
                        success_column=source_columns[
                            metric_name
                        ],
                        seed=ci_seed,
                    )
                )

        overall_rows.append(
            {
                "sample_size": total_sample_size,
                "n_train": n_train,
                "n_validation": n_validation,
                "model": model_name,

                # Keep both names for compatibility.
                "depth": mean_depth,
                "mean_depth": mean_depth,

                "n_templates": int(
                    template_group[
                        "template_id"
                    ].nunique()
                ),
                "simulations_per_template": (
                    simulations_per_template
                ),
                "total_simulations": len(
                    replicate_group
                ),

                "mean_auc": mean_auc,
                "sd_auc_between_templates": float(
                    template_group[
                        "mean_auc"
                    ].std(ddof=1)
                )
                if len(template_group) > 1
                else np.nan,

                "probability_of_success": (
                    probability_of_success
                ),
                # Convenient plotting alias.
                "power": probability_of_success,
                "probability_of_success_ci_low": (
                    ci_results[
                        "probability_of_success"
                    ][0]
                ),
                "probability_of_success_ci_high": (
                    ci_results[
                        "probability_of_success"
                    ][1]
                ),

                "detection_power": detection_power,
                "detection_power_ci_low": (
                    ci_results["detection_power"][0]
                ),
                "detection_power_ci_high": (
                    ci_results["detection_power"][1]
                ),

                "target_attainment_probability": (
                    target_attainment_probability
                ),
                "target_attainment_probability_ci_low": (
                    ci_results[
                        "target_attainment_probability"
                    ][0]
                ),
                "target_attainment_probability_ci_high": (
                    ci_results[
                        "target_attainment_probability"
                    ][1]
                ),

                "between_template_power_sd": float(
                    template_group[
                        "probability_of_success"
                    ].std(ddof=1)
                )
                if len(template_group) > 1
                else np.nan,

                "template_power_q025": float(
                    template_group[
                        "probability_of_success"
                    ].quantile(0.025)
                ),
                "template_power_q975": float(
                    template_group[
                        "probability_of_success"
                    ].quantile(0.975)
                ),

                "mean_feature_recall": float(
                    template_group[
                        "mean_feature_recall"
                    ].mean()
                ),
                "mean_feature_precision": float(
                    template_group[
                        "mean_feature_precision"
                    ].mean()
                ),
                "ci_method": ci_method,
            }
        )

    power_curve = (
        pd.DataFrame(overall_rows)
        .sort_values(
            ["model", "mean_depth", "sample_size"]
        )
        .reset_index(drop=True)
    )

    return {
        "replicate_results": replicate_results,
        "template_summary": template_summary,
        "power_curve": power_curve,
        "run_metadata": {
            "sample_sizes": sample_sizes,
            "train_fraction": train_fraction,
            "n_templates": len(templates),
            "simulations_per_template": (
                simulations_per_template
            ),
            "total_cohort_simulations": (
                len(templates)
                * len(sample_sizes)
                * simulations_per_template
            ),
            "ci_method": ci_method,
            "ci_metrics": list(ci_metrics),
            "confidence": (
                confidence
                if ci_method != "none"
                else None
            ),
            "n_bootstrap": (
                n_bootstrap
                if ci_method
                == "hierarchical_bootstrap"
                else None
            ),
        },
    }


def plot_power_by_sample_size(
    power_curve,
    *,
    model="logreg",
    power_metric="probability_of_success",
    show_ci=False,
    target_power=0.80,
    sample_size_column="sample_size",
    depth_column="mean_depth",
    ci_low_column=None,
    ci_high_column=None,
    title=None,
    ax=None,
):
    """
    Plot model-level power against total sample size, with one line per
    sequencing depth.

    Parameters
    ----------
    power_curve : pandas.DataFrame
        Output result["power_curve"] from run_power_sample_size_grid().

    model : str, default="logreg"
        Model to display.

    power_metric : str, default="probability_of_success"
        Y-axis quantity. Common choices:

            "probability_of_success"
            "detection_power"
            "target_attainment_probability"
            "power"

    show_ci : bool, default=False
        Draw CI ribbons when finite CI values are available.

    target_power : float or None, default=0.80
        Horizontal target-power reference line. Set to None to omit.

    sample_size_column : str, default="sample_size"
        X-axis column.

    depth_column : str, default="mean_depth"
        Column used to define separate depth curves.

    ci_low_column, ci_high_column : str or None
        Optional explicit CI-column names. If omitted, they are inferred from
        `power_metric`.

    title : str or None
        Plot title.

    ax : matplotlib.axes.Axes or None
        Existing axis. A new figure and axis are created when omitted.

    Returns
    -------
    fig, ax
    """
    import matplotlib.pyplot as plt
    import numpy as np

    required_columns = {
        "model",
        sample_size_column,
        depth_column,
        power_metric,
    }

    missing_columns = (
        required_columns - set(power_curve.columns)
    )

    if missing_columns:
        raise KeyError(
            "power_curve is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    plot_data = (
        power_curve.loc[
            power_curve["model"] == model
        ]
        .copy()
        .sort_values(
            [depth_column, sample_size_column]
        )
    )

    if plot_data.empty:
        raise ValueError(
            f"No power-curve rows were found for model={model!r}."
        )

    # Infer CI columns when possible.
    ci_column_map = {
        "power": (
            "probability_of_success_ci_low",
            "probability_of_success_ci_high",
        ),
        "probability_of_success": (
            "probability_of_success_ci_low",
            "probability_of_success_ci_high",
        ),
        "detection_power": (
            "detection_power_ci_low",
            "detection_power_ci_high",
        ),
        "target_attainment_probability": (
            "target_attainment_probability_ci_low",
            "target_attainment_probability_ci_high",
        ),
    }

    if ci_low_column is None or ci_high_column is None:
        inferred_columns = ci_column_map.get(
            power_metric
        )

        if inferred_columns is not None:
            ci_low_column = (
                ci_low_column or inferred_columns[0]
            )
            ci_high_column = (
                ci_high_column or inferred_columns[1]
            )

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure

    ci_was_drawn = False

    for mean_depth, group in plot_data.groupby(
        depth_column,
        sort=True,
    ):
        group = group.sort_values(
            sample_size_column
        )

        x = group[
            sample_size_column
        ].to_numpy(dtype=float)

        y = group[
            power_metric
        ].to_numpy(dtype=float)

        line = ax.plot(
            x,
            y,
            marker="o",
            label=f"{mean_depth}×",
        )[0]

        if show_ci:
            ci_columns_available = (
                ci_low_column is not None
                and ci_high_column is not None
                and ci_low_column in group.columns
                and ci_high_column in group.columns
            )

            if ci_columns_available:
                lower = group[
                    ci_low_column
                ].to_numpy(dtype=float)

                upper = group[
                    ci_high_column
                ].to_numpy(dtype=float)

                valid = (
                    np.isfinite(x)
                    & np.isfinite(lower)
                    & np.isfinite(upper)
                )

                if valid.any():
                    ax.fill_between(
                        x[valid],
                        lower[valid],
                        upper[valid],
                        alpha=0.15,
                        color=line.get_color(),
                    )
                    ci_was_drawn = True

    if show_ci and not ci_was_drawn:
        warnings.warn(
            "show_ci=True, but no finite confidence intervals were "
            "available. Run run_power_sample_size_grid() with "
            "ci_method='pooled_wilson' or "
            "ci_method='hierarchical_bootstrap'.",
            RuntimeWarning,
        )

    if target_power is not None:
        ax.axhline(
            float(target_power),
            linestyle="--",
            linewidth=1,
            label=f"Target = {target_power:.0%}",
        )

    y_label_map = {
        "power": "Probability of success",
        "probability_of_success": (
            "Probability of success"
        ),
        "detection_power": "Detection power",
        "target_attainment_probability": (
            "Target-attainment probability"
        ),
    }

    ax.set_xlabel("Total sample size")
    ax.set_ylabel(
        y_label_map.get(
            power_metric,
            power_metric.replace("_", " ").title(),
        )
    )
    ax.set_ylim(0, 1.02)

    if title is None:
        title = (
            f"Model-level power by sample size: {model}"
        )

    ax.set_title(title)
    ax.legend(
        title="Mean depth",
        frameon=False,
    )

    return fig, ax
