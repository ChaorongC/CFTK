import numpy as np
import pandas as pd
import pytest


def test_heatmap_accepts_duplicate_interval_features(tmp_path, monkeypatch):
    pytest.importorskip("seaborn")
    monkeypatch.syspath_prepend("src")
    from visualization.plot_differential import plot_heatmap

    matrix = pd.DataFrame(
        {
            "control_1": [1.0, 2.0, np.nan, 4.0],
            "case_1": [2.0, 3.0, 5.0, 6.0],
        },
        index=[
            "chr1:10-20",
            "chr1:10-20",
            "chr1:30-40",
            "chr1:50-60",
        ],
    )
    matrix_path = tmp_path / "wps.tsv"
    matrix.to_csv(matrix_path, sep="\t")

    plot_heatmap(
        str(matrix_path),
        {"control": ["control_1"], "case": ["case_1"]},
        str(tmp_path / "heatmap.png"),
        str(tmp_path / "heatmap.pdf"),
        feature_name="wps",
        top_n=None,
    )

    assert (tmp_path / "heatmap.png").is_file()
    assert (tmp_path / "heatmap.pdf").is_file()
