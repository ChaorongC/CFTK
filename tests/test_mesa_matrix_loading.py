import pandas as pd


def test_mesa_loader_preserves_duplicate_features_with_unique_positions(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend("src")
    from analysis.mesa import _load_matrices

    matrix = pd.DataFrame(
        {"sample_a": [1.0, 2.0, 3.0], "sample_b": [4.0, 5.0, 6.0]},
        index=["chr1:10-20", "chr1:10-20", "chr1:30-40"],
    )
    path = tmp_path / "wps.tsv"
    matrix.to_csv(path, sep="\t")

    loaded = _load_matrices([str(path)])[0]

    assert loaded.shape == (2, 3)
    assert loaded.columns.is_unique
    assert loaded.columns.tolist() == [0, 1, 2]
    assert loaded.iloc[0].tolist() == [1.0, 2.0, 3.0]
