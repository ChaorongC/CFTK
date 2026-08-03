import ast
import json
from pathlib import Path

import pytest


def test_sample_size_parsing_and_limits():
    from apps.model_power_app_helpers import (
        AppValidationError,
        parse_sample_size_range,
        parse_sample_size_text,
    )

    parsed = parse_sample_size_text("100, 50, 100, 200")
    assert parsed.values == (50, 100, 200)
    assert parsed.warnings == ("Duplicate sample sizes were removed.",)
    assert parse_sample_size_range(start=50, stop=100, step=25).values == (50, 75, 100)

    for text in ("", "10, abc", "0, 20", "1,2,3,4,5,6,7,8,9", "2001"):
        with pytest.raises(AppValidationError):
            parse_sample_size_text(text)


def test_ratio_and_class_count_validation():
    from apps.model_power_app_helpers import (
        AppValidationError,
        case_control_ratio,
        validate_class_counts,
    )

    assert case_control_ratio("1:3") == pytest.approx(1 / 3)
    assert case_control_ratio("1:1") == pytest.approx(1.0)
    assert case_control_ratio("3:1") == pytest.approx(3.0)
    assert case_control_ratio("Custom", custom_ratio=1.25) == pytest.approx(1.25)
    validate_class_counts(sample_sizes=(30,), ratio=1.0, cv_folds=5)
    with pytest.raises(AppValidationError):
        validate_class_counts(sample_sizes=(12,), ratio=1.0, cv_folds=10)


def test_manifest_depth_loading(tmp_path):
    from apps.model_power_app_helpers import available_depths_from_manifest

    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "row_count": 1,
                "depth_labels": ["5", "10", "30"],
                "stats": ["mean", "CI_l", "CI_u"],
                "index_file": "cpg_index.npz",
                "mean_file": "cpg_mean.float32.npy",
                "std_files": {"5": {}, "10": {}, "30": {}},
            }
        )
    )
    assert available_depths_from_manifest(tmp_path) == (5, 10, 30)


def test_input_validation():
    from apps.model_power_app_helpers import (
        AppValidationError,
        validate_analysis_inputs,
        validate_biomarker_inputs,
    )

    validate_biomarker_inputs(
        n_features=100,
        n_signal_cpgs=20,
        top_k=10,
        meth_diff=0.05,
        effect_sd=0.015,
        within_block_rho=0.2,
        sd_stat="mean",
    )
    validate_analysis_inputs(
        target_auc=0.75,
        min_observed_fraction=0.5,
        specificity_target=0.9,
        cv_folds=5,
    )

    with pytest.raises(AppValidationError):
        validate_biomarker_inputs(
            n_features=100,
            n_signal_cpgs=101,
            top_k=10,
            meth_diff=0.05,
            effect_sd=0.015,
            within_block_rho=0.2,
            sd_stat="mean",
        )
    with pytest.raises(AppValidationError):
        validate_biomarker_inputs(
            n_features=100,
            n_signal_cpgs=20,
            top_k=101,
            meth_diff=0.05,
            effect_sd=0.015,
            within_block_rho=0.2,
            sd_stat="mean",
        )
    with pytest.raises(AppValidationError):
        validate_analysis_inputs(
            target_auc=0.95,
            min_observed_fraction=0.5,
            specificity_target=0.9,
            cv_folds=5,
        )


def test_precision_modes():
    from apps.model_power_app_helpers import precision_settings

    assert precision_settings("Fast") == {
        "n_templates": 5,
        "simulations_per_template": 10,
    }
    assert precision_settings("Standard") == {
        "n_templates": 10,
        "simulations_per_template": 20,
    }


def test_workload_warning_uses_selected_precision_baseline():
    from apps.model_power_app_helpers import AppValidationError, workload_warning

    assert workload_warning(
        sample_sizes=(50, 100, 150, 200),
        depths=(10, 30),
        n_templates=5,
        simulations_per_template=30,
    ) is None
    assert workload_warning(
        sample_sizes=(50, 100, 150, 200),
        depths=(10, 30),
        n_templates=10,
        simulations_per_template=50,
    ) is None
    assert workload_warning(
        sample_sizes=(100, 200, 300, 400),
        depths=(5, 10, 20, 30),
        n_templates=5,
        simulations_per_template=20,
    ) is not None

    with pytest.raises(AppValidationError, match="computational limit"):
        workload_warning(
            sample_sizes=(500, 1000, 1500, 2000),
            depths=(10, 20, 30, 40),
            n_templates=10,
            simulations_per_template=20,
        )


def test_app_uses_null_calibrated_engine_and_no_explicit_penalty():
    runtime_source = Path("apps/model_power_app_runtime.py").read_text()
    runtime_tree = ast.parse(runtime_source)
    assert any(
        isinstance(node, ast.ImportFrom)
        and node.module == "analysis.model_power_operating_characteristics"
        and any(alias.name == "run_power_sample_size_grid" for alias in node.names)
        for node in ast.walk(runtime_tree)
    )

    page_source = Path("apps/model_power_app_page.py").read_text()
    assert "Probability of success (conservative)" in page_source
    assert "Target-attainment probability" in page_source
    assert "Detection power" in page_source
    assert "train_fraction" not in runtime_source + page_source

    for path in (
        Path("apps/model_power_calculator.py"),
        Path("apps/model_power_app_helpers.py"),
        Path("apps/model_power_app_runtime.py"),
        Path("apps/model_power_app_page.py"),
    ):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.keyword):
                assert node.arg != "penalty"


def test_app_reference_dir_stays_repository_local():
    runtime_tree = ast.parse(Path("apps/model_power_app_runtime.py").read_text())
    assignments = {
        target.id: node.value
        for node in runtime_tree.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name)
    }

    reference_expr = assignments["REFERENCE_DIR"]
    assert isinstance(reference_expr, ast.BinOp)
    assert isinstance(reference_expr.op, ast.Div)
    assert isinstance(reference_expr.left, ast.Name)
    assert reference_expr.left.id == "ROOT"
    assert isinstance(reference_expr.right, ast.Constant)
    assert reference_expr.right.value == "data"
