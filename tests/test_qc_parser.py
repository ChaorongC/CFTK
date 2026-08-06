from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "contents",
    [
        (
            "#bamPEFragmentSize\n"
            "Size\tOccurrences\tSample\n"
            "100\t1\tsample.bam\n"
            "167\t2\tsample.bam\n"
            "320\t1\tsample.bam\n"
        ),
        (
            "# legacy headerless output\n"
            "100\t1\n"
            "167\t2\n"
            "320\t1\n"
        ),
    ],
)
def test_fragment_parser_supports_deeptools_and_legacy_formats(
    contents, monkeypatch, tmp_path
):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1] / "src"))
    from analysis.qc_parser import parse_fragment_csv

    fragment_csv = tmp_path / "fragment_length.sample.raw.csv"
    fragment_csv.write_text(contents)

    metrics = parse_fragment_csv(str(fragment_csv))

    assert metrics == {
        "median_frag_len": 167,
        "mononuc_peak_bp": 167,
        "short_frag_ratio": 25.0,
        "mono_nuc_ratio": 50.0,
        "di_nuc_ratio": 25.0,
    }
