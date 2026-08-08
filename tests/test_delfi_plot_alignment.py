from pathlib import Path


def _write_delfi(path, rows):
    path.write_text(
        "#contig\tstart\tstop\tarm\tratio\tratio_corrected\n"
        + "".join(
            f"{chrom}\t{start}\t{stop}\t{chrom}\t{ratio}\t{corrected}\n"
            for chrom, start, stop, ratio, corrected in rows
        )
    )


def test_delfi_cohort_alignment_joins_genomic_bins_not_row_positions(
    monkeypatch, tmp_path
):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from visualization.plot_fragmentomics import _align_delfi_ratios

    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    _write_delfi(first, [
        ("chr1", 0, 100, 1.0, 1.1),
        ("chr1", 100, 200, 2.0, 2.1),
        ("chr2", 0, 100, 3.0, 3.1),
    ])
    _write_delfi(second, [
        ("chr1", 0, 100, 4.0, 4.1),
        ("chr2", 0, 100, 5.0, 5.1),
    ])

    aligned, value_columns, chrom_order, chrom_positions = _align_delfi_ratios(
        [first, second]
    )

    assert list(aligned[["contig", "start", "stop"]].itertuples(index=False, name=None)) == [
        ("chr1", 0, 100),
        ("chr2", 0, 100),
    ]
    assert aligned[value_columns].to_numpy().tolist() == [[1.1, 4.1], [3.1, 5.1]]
    assert chrom_order[:2] == ["chr1", "chr2"]
    assert chrom_positions == {"chr1": 0, "chr2": 1}


def test_delfi_group_and_comparison_plots_accept_different_retained_bins(
    monkeypatch, tmp_path
):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from visualization.plot_fragmentomics import (
        plot_delfi_comparison,
        plot_delfi_group,
    )

    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    _write_delfi(first, [
        ("chr1", 0, 100, 1.0, 1.1),
        ("chr1", 100, 200, 2.0, 2.1),
        ("chr2", 0, 100, 3.0, 3.1),
    ])
    _write_delfi(second, [
        ("chr1", 0, 100, 4.0, 4.1),
        ("chr2", 0, 100, 5.0, 5.1),
    ])

    group_png = tmp_path / "group.png"
    group_pdf = tmp_path / "group.pdf"
    comparison_png = tmp_path / "comparison.png"
    comparison_pdf = tmp_path / "comparison.pdf"
    plot_delfi_group([first, second], group_png, group_pdf)
    plot_delfi_comparison(
        {"Control": [first], "Case": [second]}, comparison_png, comparison_pdf
    )

    for output in (group_png, group_pdf, comparison_png, comparison_pdf):
        assert output.is_file() and output.stat().st_size > 0
