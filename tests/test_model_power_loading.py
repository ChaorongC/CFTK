import importlib
import ast
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


class ModelPowerSourcePathTests(unittest.TestCase):
    def test_default_reference_dir_source_points_to_repo_data(self):
        source_path = Path(__file__).resolve().parents[1] / "src" / "analysis" / "model_power.py"
        tree = ast.parse(source_path.read_text())

        assignments = [
            node for node in tree.body
            if isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "_DEFAULT_REFERENCE_DIR"
                for target in node.targets
            )
        ]

        self.assertEqual(len(assignments), 1)
        expression = ast.unparse(assignments[0].value)
        constants = [
            node.value
            for node in ast.walk(assignments[0].value)
            if isinstance(node, ast.Constant)
        ]
        self.assertNotIn("model_power", expression)
        self.assertIn("data", constants)


class ModelPowerLoadingTests(unittest.TestCase):
    def setUp(self):
        sys.modules.pop("analysis.model_power", None)
        self._original_statsmodels_modules = {
            name: sys.modules.get(name)
            for name in (
                "statsmodels",
                "statsmodels.stats",
                "statsmodels.stats.proportion",
            )
        }
        statsmodels = types.ModuleType("statsmodels")
        stats = types.ModuleType("statsmodels.stats")
        proportion = types.ModuleType("statsmodels.stats.proportion")

        def proportion_confint(count, nobs, alpha=0.05, method="wilson"):
            estimate = count / nobs
            return max(0.0, estimate - 0.1), min(1.0, estimate + 0.1)

        proportion.proportion_confint = proportion_confint
        stats.proportion = proportion
        statsmodels.stats = stats
        sys.modules["statsmodels"] = statsmodels
        sys.modules["statsmodels.stats"] = stats
        sys.modules["statsmodels.stats.proportion"] = proportion

    def tearDown(self):
        sys.modules.pop("analysis.model_power", None)
        for name, module in self._original_statsmodels_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module

    def _reference_dir(self, root: Path) -> Path:
        import pandas as pd

        from analysis.model_power_reference import convert_pickles_to_array_reference

        index = pd.Index(["chr1_10", "chr1_20", "chr2_30"], dtype=object)
        mean = pd.Series([0.10, 0.20, 0.30], index=index)
        std = pd.DataFrame(
            {
                "5_mean": [0.010, 0.020, 0.030],
                "5_CI_l": [0.009, 0.019, 0.029],
                "5_CI_u": [0.011, 0.021, 0.031],
                "10_mean": [0.015, 0.025, 0.035],
                "10_CI_l": [0.014, 0.024, 0.034],
                "10_CI_u": [0.016, 0.026, 0.036],
            },
            index=index,
        )
        mean_path = root / "mean.pkl"
        std_path = root / "std.pkl"
        out_dir = root / "arrays"
        mean.to_pickle(mean_path)
        std.to_pickle(std_path)
        convert_pickles_to_array_reference(mean_path, std_path, out_dir)
        return out_dir

    def test_import_does_not_read_reference_pickles(self):
        with mock.patch(
            "pandas.read_pickle",
            side_effect=AssertionError("read_pickle should not run at import"),
        ):
            module = importlib.import_module("analysis.model_power")

        self.assertTrue(hasattr(module, "load_default_model_power_reference"))

    def test_default_reference_dir_is_repo_data(self):
        module = importlib.import_module("analysis.model_power")
        expected = Path(module.__file__).resolve().parents[2] / "data"

        self.assertEqual(module._DEFAULT_REFERENCE_DIR, expected)

    def test_explicit_load_default_reference_sets_compatibility_globals(self):
        import numpy as np

        with tempfile.TemporaryDirectory() as tmp:
            reference_dir = self._reference_dir(Path(tmp))
            module = importlib.import_module("analysis.model_power")

            loaded = module.load_default_model_power_reference(
                reference_dir=reference_dir,
                depths=[10],
                sd_stats=["mean"],
                mmap_mode=None,
            )

            self.assertIs(module.cpg_mean, loaded.cpg_mean)
            self.assertIs(module.cpg_std_summary, loaded.cpg_std_summary)
            self.assertEqual(list(module.cpg_std_summary.columns), ["10_mean"])
            self.assertEqual(module.cpg_mean.dtype, np.float32)

    def test_environment_reference_dir_is_used(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_dir = self._reference_dir(Path(tmp))
            module = importlib.import_module("analysis.model_power")

            with mock.patch.dict(
                os.environ,
                {"CFTK_MODEL_POWER_DATA": str(reference_dir)},
                clear=False,
            ):
                loaded = module.load_default_model_power_reference(
                    depths=[5],
                    sd_stats=["mean"],
                    mmap_mode=None,
                )

            self.assertEqual(list(loaded.cpg_std_summary.columns), ["5_mean"])

    def test_explicit_reference_dir_overrides_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference_dir = self._reference_dir(Path(tmp))
            module = importlib.import_module("analysis.model_power")

            with mock.patch.dict(
                os.environ,
                {"CFTK_MODEL_POWER_DATA": str(Path(tmp) / "missing")},
                clear=False,
            ):
                loaded = module.load_default_model_power_reference(
                    reference_dir=reference_dir,
                    depths=[10],
                    sd_stats=["mean"],
                    mmap_mode=None,
                )

            self.assertEqual(list(loaded.cpg_std_summary.columns), ["10_mean"])

    def test_missing_reference_explains_that_arrays_are_external(self):
        with tempfile.TemporaryDirectory() as tmp:
            module = importlib.import_module("analysis.model_power")
            module._DEFAULT_REFERENCE_DIR = Path(tmp) / "missing"

            with mock.patch.dict(os.environ, {}, clear=True):
                with self.assertRaisesRegex(
                    FileNotFoundError,
                    "does not bundle its reference arrays",
                ):
                    module.load_default_model_power_reference(depths=[10])


if __name__ == "__main__":
    unittest.main()
