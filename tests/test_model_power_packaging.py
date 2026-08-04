from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_package_includes_model_power_code_but_not_reference_data():
    pyproject = (ROOT / "pyproject.toml").read_text()
    manifest = (ROOT / "MANIFEST.in").read_text()

    assert (
        'include = ["analysis*", "visualization*", "report*", '
        '"cftk_registry*"]'
    ) in pyproject
    assert 'cftk_registry = ["registry.json"]' in pyproject
    assert "include-package-data = false" in pyproject
    assert "prune data" in manifest
    assert "global-exclude *.npy *.npz *.pkl" in manifest
    assert (ROOT / "src" / "analysis" / "model_power.py").is_file()
    assert (ROOT / "src" / "cftk_registry" / "registry.json").is_file()
    assert not (ROOT / "src" / "data").exists()
