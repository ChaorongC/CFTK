import base64
import gzip
import json
from pathlib import Path


def _source_root() -> Path:
    return Path(__file__).resolve().parents[1] / "src"


def test_load_mqc_data_decodes_lzstring_payload(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(_source_root()))
    from report.mqc_extractor import load_mqc_data

    html = tmp_path / "multiqc_report.html"
    html.write_text(
        '<script id="mqc_compressed_plotdata" type="text/plain">'
        "N4IgZghgzgLgjgYwPpQKZwK6oHYNUhAew2xiiQAcAbQmEALmAF8mg==="
        "</script>"
    )

    assert load_mqc_data(str(html)) == {"fastqc_sequence_counts_plot": {}}


def test_load_mqc_data_skips_empty_placeholder_payload(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(_source_root()))
    from report.mqc_extractor import load_mqc_data

    html = tmp_path / "multiqc_report.html"
    html.write_text(
        '<script id="mqc_compressed_plotdata" type="text/plain">N4XyA===</script>'
        '<script id="mqc_compressed_plotdata" type="text/plain">'
        "N4IgZghgzgLgjgYwPpQKZwK6oHYNUhAew2xiiQAcAbQmEALmAF8mg==="
        "</script>"
    )

    assert load_mqc_data(str(html)) == {"fastqc_sequence_counts_plot": {}}


def test_load_mqc_data_keeps_legacy_gzip_compatibility(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(_source_root()))
    from report.mqc_extractor import load_mqc_data

    encoded = base64.b64encode(gzip.compress(json.dumps({"legacy": {}}).encode())).decode()
    html = tmp_path / "multiqc_report.html"
    html.write_text(
        f'<script id="mqc_compressed_plotdata" type="text/plain">{encoded}</script>'
    )

    assert load_mqc_data(str(html)) == {"legacy": {}}


def test_merge_mqc_dicts_retains_bar_line_and_heatmap_samples(monkeypatch):
    monkeypatch.syspath_prepend(str(_source_root()))
    from report.mqc_extractor import merge_mqc_dicts

    def plots(sample, value):
        return {
            "bar": {"datasets": [{"samples": [sample], "cats": [
                {"name": "Reads", "data": [value], "data_pct": [100]},
            ]}]},
            "line": {"datasets": [{"lines": [
                {"name": sample, "pairs": [[1, value]]},
            ]}]},
            "heatmap": {"datasets": [{
                "xcats": ["Basic Statistics"],
                "ycats": [sample],
                "rows": [[1]],
            }]},
        }

    merged = merge_mqc_dicts([plots("sample_a", 10), plots("sample_b", 20)])

    bar = merged["bar"]["datasets"][0]
    assert bar["samples"] == ["sample_a", "sample_b"]
    assert bar["cats"][0]["data"] == [10, 20]
    assert len(merged["line"]["datasets"][0]["lines"]) == 2
    heatmap = merged["heatmap"]["datasets"][0]
    assert heatmap["ycats"] == ["sample_a", "sample_b"]
    assert heatmap["rows"] == [[1], [1]]


def test_normalise_legacy_scalar_line_data(monkeypatch):
    monkeypatch.syspath_prepend(str(_source_root()))
    from report.mqc_extractor import _normalise_mqc_payload

    payload = {
        "duplication": {
            "plot_type": "xy_line",
            "datasets": [[{"name": "sample", "data": [0.2, 0.8]}]],
            "config": {"categories": [1, 2]},
        }
    }

    normalized = _normalise_mqc_payload(payload)
    assert normalized["duplication"]["datasets"][0]["lines"][0]["pairs"] == [
        [1, 0.2], [2, 0.8]
    ]
