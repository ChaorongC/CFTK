"""
mqc_extractor.py — Extract MultiQC HTML compressed plot data and convert to
                   Plotly.js interactive charts for the cftk report.

MultiQC stores all plot data as Base64-compressed JSON in a DOM element:
  <script id="mqc_compressed_plotdata" type="application/json">
    N4IgZghg...base64...
  </script>

Current MultiQC releases use LZString Base64 compression. Older reports used
gzip+Base64, which remains supported for compatibility.

This module:
  1. Loads and decompresses that data (with file-level caching).
  2. Converts each MultiQC plot format to a Plotly.js figure dict.
  3. Renders the figure as an inline HTML div with embedded JSON.

Supported MultiQC plot types → Plotly equivalents:
  bar plot   → plotly bar (stacked), with Counts/Percentages toggle
  x/y line   → plotly scatter (lines), with dataset selector buttons
  heatmap    → plotly heatmap, with Pass/Warn/Fail colour mapping

Usage:
  from report.mqc_extractor import load_mqc_data, mqc_bar_chart, mqc_line_chart
  from report.mqc_extractor import mqc_heatmap, mqc_mbias_chart

  data = load_mqc_data("/path/to/multiqc_report.html")
  html = mqc_bar_chart("cutadapt_filtered_reads_plot", data,
                        title="Filtered Reads", yaxis="Reads")
"""

from __future__ import annotations

import base64
import gzip
import json
import os
import re

# ── file-level cache: html_path → decompressed dict ──────────────────────────
_MQC_CACHE: dict[str, dict] = {}

_LZSTRING_BASE64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
_LZSTRING_INDEX = {char: index for index, char in enumerate(_LZSTRING_BASE64)}


def _lzstring_decompress_from_base64(value: str) -> str | None:
    """Decode ``LZString.decompressFromBase64`` using only the stdlib."""
    value = "".join(value.split())
    if not value:
        return ""
    if any(char not in _LZSTRING_INDEX for char in value):
        return None

    value_len = len(value)
    data_value = _LZSTRING_INDEX[value[0]]
    data_position = 32
    data_index = 1

    def read_bits(width: int) -> int | None:
        nonlocal data_value, data_position, data_index
        bits = 0
        power = 1
        limit = 1 << width
        while power != limit:
            bit = data_value & data_position
            data_position >>= 1
            if data_position == 0:
                data_position = 32
                if data_index >= value_len:
                    return None
                data_value = _LZSTRING_INDEX[value[data_index]]
                data_index += 1
            if bit:
                bits |= power
            power <<= 1
        return bits

    dictionary: dict[int, str] = {0: "0", 1: "1", 2: "2"}
    next_code = read_bits(2)
    if next_code is None:
        return None
    if next_code == 0:
        char_code = read_bits(8)
    elif next_code == 1:
        char_code = read_bits(16)
    elif next_code == 2:
        return ""
    else:
        return None
    if char_code is None:
        return None

    dictionary[3] = chr(char_code)
    w = dictionary[3]
    result = [w]
    dict_size = 4
    enlarge_in = 4
    num_bits = 3

    while True:
        code = read_bits(num_bits)
        if code is None:
            return None
        if code == 0:
            char_code = read_bits(8)
            if char_code is None:
                return None
            dictionary[dict_size] = chr(char_code)
            code = dict_size
            dict_size += 1
            enlarge_in -= 1
        elif code == 1:
            char_code = read_bits(16)
            if char_code is None:
                return None
            dictionary[dict_size] = chr(char_code)
            code = dict_size
            dict_size += 1
            enlarge_in -= 1
        elif code == 2:
            return "".join(result)

        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1

        if code in dictionary:
            entry = dictionary[code]
        elif code == dict_size:
            entry = w + w[0]
        else:
            return None

        result.append(entry)
        dictionary[dict_size] = w + entry[0]
        dict_size += 1
        enlarge_in -= 1
        w = entry

        if enlarge_in == 0:
            enlarge_in = 1 << num_bits
            num_bits += 1


def _decode_mqc_payload(encoded: str) -> dict | None:
    """Return a decoded MultiQC payload, accepting current and legacy formats."""
    errors = []
    try:
        decoded = _lzstring_decompress_from_base64(encoded)
        if decoded is not None:
            payload = json.loads(decoded)
            if isinstance(payload, dict):
                return payload
            errors.append("LZString payload is not a JSON object")
    except Exception as exc:
        errors.append(f"LZString: {exc}")

    try:
        raw = base64.b64decode(encoded)
        payload = json.loads(gzip.decompress(raw).decode("utf-8"))
        if isinstance(payload, dict):
            return payload
        errors.append("gzip payload is not a JSON object")
    except Exception as exc:
        errors.append(f"gzip: {exc}")

    raise ValueError("; ".join(errors) or "unsupported MultiQC payload")


def _normalise_mqc_payload(data: dict) -> dict:
    """Normalise legacy MultiQC plot records to this module's renderer schema."""
    for payload in data.values():
        if not isinstance(payload, dict):
            continue
        plot_type = payload.get("plot_type")
        datasets = payload.get("datasets")

        # MultiQC 1.10 writes bar and line datasets as lists of lists, with
        # sample labels parallel to the bar datasets and point pairs in ``data``.
        if plot_type == "bar_graph" and isinstance(datasets, list) and datasets and isinstance(datasets[0], list):
            sample_sets = payload.get("samples", [])
            if sample_sets and not isinstance(sample_sets[0], (list, tuple)):
                sample_sets = [sample_sets] * len(datasets)
            payload["datasets"] = [
                {
                    "samples": sample_sets[index] if index < len(sample_sets) else [],
                    "cats": categories,
                }
                for index, categories in enumerate(datasets)
            ]
        elif plot_type == "xy_line" and isinstance(datasets, list) and datasets and isinstance(datasets[0], list):
            categories = payload.get("config", {}).get("categories", [])

            def normalised_pairs(line: dict) -> list:
                pairs = line.get("pairs", line.get("data", []))
                if not pairs or isinstance(pairs[0], (list, tuple)):
                    return pairs
                x_values = categories if len(categories) == len(pairs) else list(range(1, len(pairs) + 1))
                return [[x, y] for x, y in zip(x_values, pairs)]

            payload["datasets"] = [
                {
                    "label": f"Dataset {index + 1}",
                    "lines": [
                        {
                            **line,
                            "pairs": normalised_pairs(line),
                        }
                        for line in lines
                    ],
                }
                for index, lines in enumerate(datasets)
            ]
        elif plot_type == "heatmap" and "datasets" not in payload:
            xcats = payload.get("xcats", [])
            ycats = payload.get("ycats", [])
            rows = [[0] * len(xcats) for _ in ycats]
            for point in payload.get("data", []):
                if not isinstance(point, list) or len(point) < 3:
                    continue
                x, y, value = point[:3]
                if isinstance(x, int) and isinstance(y, int) and 0 <= y < len(rows) and 0 <= x < len(xcats):
                    rows[y][x] = value
            payload["datasets"] = [{"xcats": xcats, "ycats": ycats, "rows": rows}]
    return data


def load_mqc_data(html_path: str) -> dict:
    """
    Load and decompress the MultiQC plot data from a multiqc_report.html file.
    Returns an empty dict if the file doesn't exist or cannot be parsed.
    Results are cached by file path.
    """
    if not html_path or not os.path.exists(html_path):
        return {}

    if html_path in _MQC_CACHE:
        return _MQC_CACHE[html_path]

    try:
        with open(html_path, encoding="utf-8", errors="replace") as f:
            content = f.read()

        matches = list(re.finditer(
            r'id=["\']mqc_compressed_plotdata["\'][^>]*>(.*?)</(?:script|div)',
            content, re.DOTALL,
        ))
        if not matches:
            _MQC_CACHE[html_path] = {}
            return {}

        decoded_empty = None
        errors = []
        for match in matches:
            try:
                data = _normalise_mqc_payload(_decode_mqc_payload(match.group(1).strip()))
                if data:
                    _MQC_CACHE[html_path] = data
                    return data
                decoded_empty = data
            except Exception as exc:
                errors.append(str(exc))
        if decoded_empty is not None:
            _MQC_CACHE[html_path] = decoded_empty
            return decoded_empty
        raise ValueError("; ".join(errors) or "no usable MultiQC payload")

    except Exception as e:
        print(f"[mqc_extractor] WARNING: could not load {html_path}: {e}")
        _MQC_CACHE[html_path] = {}
        return {}


def _merge_mqc_dataset(target: dict, source: dict) -> None:
    """Merge the bar, line, or heatmap shapes consumed by this report."""
    if "lines" in source:
        target.setdefault("lines", []).extend(source.get("lines", []))
        return

    if "ycats" in source and "rows" in source:
        target.setdefault("ycats", []).extend(source.get("ycats", []))
        target.setdefault("rows", []).extend(source.get("rows", []))
        return

    if "samples" not in source or "cats" not in source:
        return

    target.setdefault("samples", []).extend(source.get("samples", []))
    target_cats = target.setdefault("cats", [])
    by_name = {str(cat.get("name", index)): cat for index, cat in enumerate(target_cats)}
    for index, source_cat in enumerate(source.get("cats", [])):
        key = str(source_cat.get("name", index))
        target_cat = by_name.get(key)
        if target_cat is None:
            target_cat = {
                name: value[:] if isinstance(value, list) else value
                for name, value in source_cat.items()
            }
            target_cats.append(target_cat)
            by_name[key] = target_cat
            continue
        for field in ("data", "data_pct"):
            if field in source_cat:
                target_cat.setdefault(field, []).extend(source_cat[field])


def merge_mqc_dicts(mqc_dicts: list[dict]) -> dict:
    """Merge MultiQC plot dictionaries while retaining the first plot metadata."""
    merged: dict = {}
    for mqc_data in mqc_dicts:
        for key, payload in mqc_data.items():
            if not isinstance(payload, dict):
                continue
            if key not in merged:
                merged[key] = json.loads(json.dumps(payload))
                continue

            target_sets = merged[key].get("datasets", [])
            source_sets = payload.get("datasets", [])
            for index, source_set in enumerate(source_sets):
                if index >= len(target_sets):
                    target_sets.append(json.loads(json.dumps(source_set)))
                elif isinstance(source_set, dict) and isinstance(target_sets[index], dict):
                    _merge_mqc_dataset(target_sets[index], source_set)
    return merged


def merge_mqc_data(html_paths: list[str]) -> dict:
    """Load and merge the readable MultiQC reports in ``html_paths``."""
    return merge_mqc_dicts([load_mqc_data(path) for path in html_paths])


# ── Plotly CDN (injected once by report_generator into <head>) ────────────────
PLOTLY_CDN = "https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.1/plotly.min.js"
PLOTLY_SCRIPT = f'<script src="{PLOTLY_CDN}"></script>'

# ── colour palette matching MultiQC defaults ──────────────────────────────────
_MQC_COLORS = [
    "#7cb5ec", "#434348", "#90ed7d", "#f7a35c",
    "#8085e9", "#f15c80", "#e4d354", "#2b908f",
    "#f45b5b", "#91e8e1",
]

# ── status heatmap colour mapping (MultiQC values → colours) ─────────────────
_STATUS_COLORSCALE = [
    [0.00, "rgba(255,255,255,0)"],   # 0.0 → N/A (transparent)
    [0.24, "rgba(255,255,255,0)"],
    [0.25, "#d9534f"],               # 0.25 → Fail (red)
    [0.49, "#d9534f"],
    [0.50, "#f0ad4e"],               # 0.5  → Warn (amber)
    [0.99, "#f0ad4e"],
    [1.00, "#5cb85c"],               # 1.0  → Pass (green)
]

_STATUS_TEXT = {1.0: "Pass", 0.5: "Warn", 0.25: "Fail", 0.0: "N/A"}


# ── uid generator ─────────────────────────────────────────────────────────────

_UID_COUNTER = [0]


def _uid(prefix: str = "mqc") -> str:
    _UID_COUNTER[0] += 1
    safe = re.sub(r"[^a-zA-Z0-9_]", "_", prefix)
    return f"{safe}_{_UID_COUNTER[0]}"


# ── base layout ──────────────────────────────────────────────────────────────

def _layout(title: str = "", xaxis_title: str = "", yaxis_title: str = "",
            barmode: str = "stack", height: int = 400) -> dict:
    return {
        "title":     {"text": title, "font": {"size": 13}},
        "height":    height,
        "margin":    {"l": 60, "r": 20, "t": 45, "b": 90},
        "barmode":   barmode,
        "legend":    {"orientation": "h", "y": -0.28, "x": 0},
        "xaxis":     {"title": xaxis_title, "automargin": True,
                      "tickangle": -40, "tickfont": {"size": 10}},
        "yaxis":     {"title": yaxis_title, "automargin": True},
        "plot_bgcolor":  "#ffffff",
        "paper_bgcolor": "#ffffff",
        "font":      {"family": "Times New Roman, Times, serif", "size": 11},
        "hoverlabel": {"bgcolor": "white", "font": {"size": 11}},
    }


# ── render helper ─────────────────────────────────────────────────────────────

def _render(div_id: str, fig: dict, height: int = 400) -> str:
    fig_json = json.dumps(fig)
    return (
        f'<div id="{div_id}" style="width:100%;height:{height}px;"></div>'
        f"<script>"
        f"(function(){{"
        f'  var fig={fig_json};'
        f'  var plt=Plotly.newPlot("{div_id}",fig.data,fig.layout,'
        f'  {{responsive:true,displayModeBar:true,'
        f'  modeBarButtonsToRemove:["lasso2d","select2d"]}});'
        # Disable legend click (toggle trace visibility) and doubleclick (isolate)
        f'  document.getElementById("{div_id}").on("plotly_legendclick",function(){{return false;}});'
        f'  document.getElementById("{div_id}").on("plotly_legenddoubleclick",function(){{return false;}});'
        f"}})();"
        f"</script>"
    )


def _missing(msg: str) -> str:
    return (
        f'<div style="padding:14px 16px;background:var(--bg2,#f8f9fa);'
        f'border:1px solid var(--rule,#dde3ea);border-radius:8px;'
        f'font-size:13px;color:var(--ink3,#7a8a9a);">{msg}</div>'
    )


def _resolve_plot_key(plot_key: str, mqc_data: dict) -> str:
    """Resolve the ``-1`` suffix used by newer MultiQC embedded plot keys."""
    if plot_key in mqc_data:
        return plot_key
    matches = sorted(
        key for key in mqc_data
        if key.startswith(f"{plot_key}-") and key[len(plot_key) + 1:].isdigit()
    )
    return matches[0] if matches else plot_key


# ── mqc_bar_chart ─────────────────────────────────────────────────────────────

def mqc_bar_chart(
    plot_key: str,
    mqc_data: dict,
    title: str = "",
    yaxis: str = "",
    dataset_idx: int = 0,
    height: int = 420,
) -> str:
    """
    Convert a MultiQC bar plot to a Plotly stacked bar chart.

    MultiQC bar format:
      datasets[i].cats = [{name, color, data:[per-sample counts],
                                          data_pct:[per-sample pct]}, ...]
      datasets[i].samples = [sample_name, ...]

    Renders with Counts/Percentages toggle buttons if data_pct is present.
    """
    plot_key = _resolve_plot_key(plot_key, mqc_data)
    if plot_key not in mqc_data:
        return _missing(f"{title}: data not found in MultiQC report "
                        f"(key: <code>{plot_key}</code>)")

    pdata = mqc_data[plot_key]
    datasets = pdata.get("datasets", [])
    if not datasets or dataset_idx >= len(datasets):
        return _missing(f"{title}: no datasets in MultiQC plot")

    ds       = datasets[dataset_idx]
    cats     = ds.get("cats", [])
    samples  = ds.get("samples", [])

    if not cats or not samples:
        return _missing(f"{title}: empty data")

    # Determine if pct switch is available
    has_pct = all("data_pct" in cat for cat in cats)
    # Determine orientation from trace_params (most reliable) or axis_controlled_by_switches
    trace_params  = ds.get("trace_params", {})
    orientation   = trace_params.get("orientation", "")
    if orientation:
        is_horizontal = (orientation == "h")
    else:
        # Fallback: axis_controlled_by_switches=['xaxis'] → horizontal (Counts on x-axis)
        pct_axis      = pdata.get("axis_controlled_by_switches", ["yaxis"])[0]
        is_horizontal = (pct_axis == "xaxis")

    traces_count = []
    traces_pct   = []

    for i, cat in enumerate(cats):
        color = cat.get("color", _MQC_COLORS[i % len(_MQC_COLORS)])
        # Normalise rgba(...) colors from MultiQC
        name  = cat.get("name", f"Cat {i}")
        count_vals = cat.get("data", [])
        pct_vals   = cat.get("data_pct", []) if has_pct else []

        common = {
            "type":        "bar",
            "name":        name,
            "marker":      {"color": color},
            "hovertemplate": f"{name}<br>%{{x}}: %{{y:,.0f}}<extra></extra>"
            if not is_horizontal else
            f"{name}<br>%{{y}}: %{{x:,.0f}}<extra></extra>",
        }
        if is_horizontal:
            traces_count.append({**common, "x": count_vals, "y": samples,
                                  "orientation": "h"})
            if has_pct:
                traces_pct.append({**common, "x": pct_vals, "y": samples,
                                   "orientation": "h"})
        else:
            traces_count.append({**common, "x": samples, "y": count_vals})
            if has_pct:
                traces_pct.append({**common, "x": samples, "y": pct_vals})

    div_id = _uid(f"bar_{plot_key[:20]}")
    layout = _layout(title=title, yaxis_title=yaxis, barmode="stack", height=height)
    if is_horizontal:
        layout["xaxis"]["title"] = yaxis
        layout["yaxis"]["title"] = ""
        layout["yaxis"]["automargin"] = True
        layout["margin"]["b"] = 40
        layout["margin"]["l"] = 130

    fig = {"data": traces_count, "layout": layout}

    if has_pct:
        # Add Counts/Percentages toggle buttons via updatemenus
        pct_axis_update = pdata.get("pct_axis_update", {})
        tick_suffix = pct_axis_update.get("ticksuffix", "%")

        fig["layout"]["updatemenus"] = [{
            "type": "buttons", "direction": "left",
            "x": 0, "y": 1.12, "xanchor": "left",
            "bgcolor": "#f0f0f0", "bordercolor": "#ccc",
            "font": {"size": 11},
            "showactive": True,
            "buttons": [
                {"label": "Counts",
                 "method": "update",
                 "args": [
                     {"x" if is_horizontal else "y": [t["x" if is_horizontal else "y"]
                                                      for t in traces_count]},
                     {f"{'x' if is_horizontal else 'y'}axis": {"ticksuffix": "",
                                                                "title": yaxis,
                                                                "automargin": True}},
                 ]},
                {"label": "Percentages",
                 "method": "update",
                 "args": [
                     {"x" if is_horizontal else "y": [t["x" if is_horizontal else "y"]
                                                      for t in traces_pct]},
                     {f"{'x' if is_horizontal else 'y'}axis": {"ticksuffix": tick_suffix,
                                                                "title": "%",
                                                                "automargin": True}},
                 ]},
            ],
        }]

    return _render(div_id, fig, height=height)


# ── mqc_line_chart ────────────────────────────────────────────────────────────

def mqc_line_chart(
    plot_key: str,
    mqc_data: dict,
    title: str = "",
    xaxis: str = "",
    yaxis: str = "",
    dataset_idx: int = 0,
    height: int = 420,
    pct_yaxis: bool = False,
) -> str:
    """
    Convert a MultiQC x/y line plot to Plotly scatter lines.

    MultiQC line format:
      datasets[i].lines = [{name, pairs:[[x,y],...], color, width, dash}, ...]

    If multiple datasets exist (e.g. Counts vs Obs/Exp), renders dataset
    selector buttons.
    """
    plot_key = _resolve_plot_key(plot_key, mqc_data)
    if plot_key not in mqc_data:
        return _missing(f"{title}: data not found (key: <code>{plot_key}</code>)")

    pdata    = mqc_data[plot_key]
    datasets = pdata.get("datasets", [])
    if not datasets:
        return _missing(f"{title}: no datasets")

    # Build traces for all datasets, mark first dataset visible
    all_traces: list[dict] = []
    dataset_labels: list[str] = []

    for di, ds in enumerate(datasets):
        label = ds.get("label", f"Dataset {di}")
        dataset_labels.append(label)
        lines = ds.get("lines", [])
        for i, line in enumerate(lines):
            name   = line.get("name", f"Sample {i}")
            pairs  = line.get("pairs", [])
            color  = line.get("color", _MQC_COLORS[i % len(_MQC_COLORS)])
            width  = line.get("width", 1.5)
            dash   = line.get("dash", "solid")
            xs     = [p[0] for p in pairs]
            ys     = [p[1] for p in pairs]
            all_traces.append({
                "type":      "scatter",
                "mode":      "lines",
                "name":      name,
                "x":         xs,
                "y":         ys,
                "visible":   (di == dataset_idx),
                "line":      {"color": color, "width": width, "dash": dash},
                "hovertemplate": f"{name}<br>x: %{{x}}<br>y: %{{y:.3f}}<extra></extra>",
            })

    div_id = _uid(f"line_{plot_key[:20]}")
    layout_kwargs = dict(title=title, xaxis_title=xaxis, yaxis_title=yaxis,
                         barmode="overlay", height=height)
    layout = _layout(**layout_kwargs)

    if pct_yaxis:
        layout["yaxis"]["ticksuffix"] = "%"

    # Hide legend for line charts — with 20+ sample traces it becomes unreadable.
    # Hover tooltip already shows the sample name on mouseover.
    layout["showlegend"] = False

    fig = {"data": all_traces, "layout": layout}

    # If multiple datasets: add dropdown/buttons to switch
    if len(datasets) > 1:
        n_per_ds = len(datasets[0].get("lines", []))
        buttons  = []
        for di, label in enumerate(dataset_labels):
            vis = [False] * len(all_traces)
            for j in range(n_per_ds):
                idx = di * n_per_ds + j
                if idx < len(vis):
                    vis[idx] = True
            buttons.append({
                "label":  label,
                "method": "update",
                "args":   [{"visible": vis}],
            })
        fig["layout"]["updatemenus"] = [{
            "type": "buttons", "direction": "left",
            "x": 0, "y": 1.12, "xanchor": "left",
            "bgcolor": "#f0f0f0", "bordercolor": "#ccc",
            "font": {"size": 11}, "showactive": True,
            "buttons": buttons,
        }]

    return _render(div_id, fig, height=height)


# ── mqc_heatmap ───────────────────────────────────────────────────────────────

def mqc_heatmap(
    plot_key: str,
    mqc_data: dict,
    title: str = "",
    height: int | None = None,
    exclude_modules: list | None = None,
) -> str:
    """
    Convert a MultiQC heatmap (status check) to Plotly heatmap.

    MultiQC heatmap format:
      datasets[0].xcats = [module_name, ...]   (x axis = modules)
      datasets[0].ycats = [sample_name, ...]   (y axis = samples)
      datasets[0].rows  = [[val, ...], ...]    (values: 1=Pass, 0.5=Warn, 0.25=Fail)

    exclude_modules: list of module names (case-insensitive) to drop from x axis.
    """
    plot_key = _resolve_plot_key(plot_key, mqc_data)
    if plot_key not in mqc_data:
        return _missing(f"{title}: data not found (key: <code>{plot_key}</code>)")

    pdata    = mqc_data[plot_key]
    datasets = pdata.get("datasets", [])
    if not datasets:
        return _missing(f"{title}: no datasets")

    ds      = datasets[0]
    xcats   = ds.get("xcats", [])
    ycats   = ds.get("ycats", [])
    rows    = ds.get("rows",  [])

    if not xcats or not ycats or not rows:
        return _missing(f"{title}: empty heatmap data")

    # Filter out excluded modules from x axis (and corresponding columns in rows)
    if exclude_modules:
        excl_lower = {m.lower().strip() for m in exclude_modules}
        keep_idx = [i for i, x in enumerate(xcats)
                    if str(x).lower().strip() not in excl_lower]
        xcats = [xcats[i] for i in keep_idx]
        rows  = [[row[i] for i in keep_idx if i < len(row)] for row in rows]

    # Build text grid for cell annotations
    text_grid = [
        [_STATUS_TEXT.get(v, "N/A") for v in row]
        for row in rows
    ]

    n_samples = len(ycats)
    h = height or max(300, n_samples * 32 + 120)

    div_id = _uid(f"hm_{plot_key[:20]}")
    fig = {
        "data": [{
            "type":          "heatmap",
            "z":             rows,
            "x":             xcats,
            "y":             ycats,
            "text":          text_grid,
            "texttemplate":  "%{text}",
            "textfont":      {"size": 10},
            "colorscale":    _STATUS_COLORSCALE,
            "zmin":          0,
            "zmax":          1,
            "showscale":     False,
            "hovertemplate": "Sample: %{y}<br>Module: %{x}<br>Status: %{text}<extra></extra>",
        }],
        "layout": {
            **_layout(title=title, height=h),
            "xaxis": {"tickangle": -35, "tickfont": {"size": 10},
                      "automargin": True, "side": "bottom"},
            "yaxis": {"tickfont": {"size": 10}, "automargin": True},
            "margin": {"l": 160, "r": 20, "t": 45, "b": 120},
        },
    }
    return _render(div_id, fig, height=h)


# ── mqc_mbias_chart ───────────────────────────────────────────────────────────

def mqc_mbias_chart(
    plot_key: str,
    mqc_data: dict,
    title: str = "M-bias",
    height: int = 450,
) -> str:
    """
    Convert a MultiQC M-bias line plot to Plotly.
    All samples shown simultaneously; dataset buttons switch context
    (e.g. CpG OT R1 / CpG OT R2 / CpG OB R1 / CpG OB R2).

    MultiQC M-bias format (bismark_mbias or methyldackel_mbias):
      datasets[0]: label='CpG R1' / 'CpG OT R1' etc.
        lines[i]: {name: sample_name, pairs: [[pos, pct], ...]}

    The pct_axis_update field specifies Y-axis tick suffix (%).
    """
    plot_key = _resolve_plot_key(plot_key, mqc_data)
    if plot_key not in mqc_data:
        return _missing(
            f"M-bias interactive data not found "
            f"(key: <code>{plot_key}</code>). "
            f"Run <code>cftk process -s 4</code> to generate mbias data, "
            f"then ensure MultiQC is run on the methylation directory."
        )

    pdata    = mqc_data[plot_key]
    datasets = pdata.get("datasets", [])
    if not datasets:
        return _missing("M-bias: no datasets in MultiQC plot")

    all_traces: list[dict]     = []
    dataset_labels: list[str] = []
    traces_per_ds: list[int]  = []

    for di, ds in enumerate(datasets):
        label = ds.get("label", f"Context {di}")
        dataset_labels.append(label)
        lines = ds.get("lines", [])
        traces_per_ds.append(len(lines))

        for i, line in enumerate(lines):
            name  = line.get("name", f"Sample {i}")
            pairs = line.get("pairs", [])
            color = line.get("color", _MQC_COLORS[i % len(_MQC_COLORS)])
            xs    = [p[0] for p in pairs]
            ys    = [p[1] for p in pairs]
            all_traces.append({
                "type":    "scatter",
                "mode":    "lines",
                "name":    name,
                "x":       xs,
                "y":       ys,
                "visible": (di == 0),
                "line":    {"color": color, "width": 1.5},
                "hovertemplate": (
                    f"{name}<br>Position: %{{x}} bp"
                    f"<br>Methylation: %{{y:.1f}}%<extra></extra>"
                ),
            })

    # Buttons to switch between contexts (datasets)
    buttons = []
    offset = 0
    for di, (label, n) in enumerate(zip(dataset_labels, traces_per_ds)):
        vis = [False] * len(all_traces)
        for j in range(n):
            if offset + j < len(vis):
                vis[offset + j] = True
        buttons.append({
            "label":  label,
            "method": "update",
            "args":   [{"visible": vis}],
        })
        offset += n

    pct_update = pdata.get("pct_axis_update", {})

    div_id = _uid("mbias")
    fig = {
        "data": all_traces,
        "layout": {
            **_layout(title=title,
                      xaxis_title="Position in Read (bp)",
                      yaxis_title="% Methylation",
                      height=height),
            "yaxis": {
                "title":      "% Methylation",
                "ticksuffix": pct_update.get("ticksuffix", "%"),
                "automargin": True,
            },
            "showlegend": False,  # sample names visible on hover; legend too cluttered
            "updatemenus": [{
                "type":      "buttons",
                "direction": "left",
                "x": 0, "y": 1.13, "xanchor": "left",
                "bgcolor":     "#f0f0f0",
                "bordercolor": "#ccc",
                "font":        {"size": 11},
                "showactive":  True,
                "buttons":     buttons,
            }],
            "annotations": [{
                "text":      "Context:",
                "showarrow": False,
                "x": 0, "y": 1.18,
                "xref": "paper", "yref": "paper",
                "xanchor": "left",
                "font": {"size": 11},
            }],
        },
    }
    return _render(div_id, fig, height=height)


# ── mqc_general_stats_table ───────────────────────────────────────────────────

def mqc_general_stats_table(html_path: str) -> dict:
    """
    Extract the general_stats_table violin/metrics data from a MultiQC HTML.
    Returns a dict: {sample_name: {metric_key: value, ...}, ...}
    Useful for extracting per-sample flagstat, FastQC, and cutadapt summary values.
    """
    mqc_data = load_mqc_data(html_path)
    if "general_stats_table" not in mqc_data:
        return {}

    pdata    = mqc_data["general_stats_table"]
    datasets = pdata.get("datasets", [])
    if not datasets:
        return {}

    ds      = datasets[0]
    metrics = ds.get("metrics", [])
    hbm     = ds.get("header_by_metric", {})

    # The violin plot stores per-sample data differently from bar/line
    # Fall back to reading multiqc_data/ TSV files for tabular data
    return {"metrics": metrics, "headers": hbm}


# ── markdup_dedup_chart ────────────────────────────────────────────────────────
# MultiQC has no built-in sambamba markdup module.
# We parse .markdup_metrics.txt files directly and render a Plotly bar chart.
#
# sambamba markdup stderr format (captured to .markdup_metrics.txt via P1 patch):
#   [...] found 1234567 duplicates in 9876543 reads (12.50% duplicate reads)
#   OR
#   [...] 1234567 / 9876543 = 12.50% duplicates
#
# We show: per-sample % duplicate and % unique as a stacked horizontal bar.

def markdup_dedup_chart(
    markdup_dir: str | list[str],
    samples: list[str],
    flagstat_dir: str | list[str] = "",
    title: str = "Deduplication (sambamba markdup)",
    height: int = 420,
) -> str:
    """
    Parse sambamba markdup metrics and render per-sample deduplication chart.

    sambamba markdup output format (captured to .markdup_metrics.txt):
        sorted N end pairs
        found M duplicates
        total time elapsed: ...

    Duplication % = found_duplicates / total_reads * 100
    total_reads sourced from {name}.bam.flagstat (in flagstat_dir / markdup_dir).

    markdup_dir:  one or more 3_markdup/ directories containing metrics files.
    flagstat_dir: one or more 2_alignment/ directories containing flagstat
                  files. Falls back to markdup_dir if not provided.
    samples:      list of sample names (ordering of bars).
    """
    import glob, re as _re, os as _os

    dup_pcts:    list[float] = []
    uniq_pcts:   list[float] = []
    valid_names: list[str]   = []
    missing:     list[str]   = []

    def _dirs(value):
        if not value:
            return []
        if isinstance(value, (list, tuple)):
            return [str(path) for path in value if path]
        return [str(value)]

    markdup_dirs = _dirs(markdup_dir)
    flagstat_dirs = _dirs(flagstat_dir) or markdup_dirs

    def _parse_flagstat_total(name):
        """Return total read count from {name}.bam.flagstat, or None."""
        # Also search markdup_dir as fallback (some pipelines put flagstat there)
        search_dirs = list(dict.fromkeys(flagstat_dirs + markdup_dirs))
        for d in search_dirs:
            for fname in [f"{name}.bam.flagstat",
                          f"{name}.markdup.bam.flagstat",
                          f"{name}.flagstat"]:
                p = _os.path.join(d, fname)
                if _os.path.exists(p):
                    try:
                        with open(p, encoding="utf-8", errors="replace") as f:
                            for line in f:
                                # "N + 0 in total" is always first line
                                m = _re.match(r"(\d+)", line.strip())
                                if m:
                                    return int(m.group(1))
                    except Exception:
                        pass
        return None

    for name in samples:
        mpath = next(
            (_os.path.join(directory, f"{name}.markdup_metrics.txt")
             for directory in markdup_dirs
             if _os.path.exists(_os.path.join(directory, f"{name}.markdup_metrics.txt"))),
            None,
        )
        if not mpath:
            missing.append(name)
            continue

        dup_count  = None
        pair_count = None   # "sorted N end pairs"
        try:
            with open(mpath, encoding="utf-8", errors="replace") as f:
                text = f.read()

            # Primary: "found N duplicates"
            m = _re.search(r"found\s+(\d[\d,]*)\s+duplicates", text)
            if m:
                dup_count = int(m.group(1).replace(",", ""))

            # Paired reads from "sorted N end pairs"
            m2 = _re.search(r"sorted\s+(\d[\d,]*)\s+end\s+pairs", text)
            if m2:
                pair_count = int(m2.group(1).replace(",", "")) * 2  # pairs → reads

        except Exception:
            pass

        if dup_count is None:
            missing.append(name)
            continue

        # Total reads: prefer flagstat (most accurate), fall back to pair_count
        total = _parse_flagstat_total(name)
        used_pairs_fallback = False
        if not total and pair_count:
            total = pair_count
            used_pairs_fallback = True
        if not total:
            missing.append(name)
            continue

        pct = round(dup_count / total * 100, 2)
        dup_pcts.append(pct)
        uniq_pcts.append(round(100.0 - pct, 2))
        valid_names.append(name)

    if not valid_names:
        files_found = [
            path for directory in markdup_dirs
            for path in glob.glob(_os.path.join(directory, "*.markdup_metrics.txt"))
        ]
        hint = (
            f"Found {len(files_found)} .markdup_metrics.txt file(s) but could not extract "
            f"duplicate count or total reads.<br>"
            f"Expected: <em>found N duplicates</em> line in sambamba output, "
            f"and a matching <code>.bam.flagstat</code> file."
            if files_found else
            f"No .markdup_metrics.txt files found. "
            f"Ensure the P1 patch is applied (sambamba stderr → file) and "
            f"<code>cftk process -s 3</code> has been re-run."
        )
        return _missing(hint)

    # Build both counts and pct traces
    # Counts: need raw unique/dup read counts
    # We have dup_pcts and uniq_pcts; reconstruct counts from flagstat totals
    total_reads = []
    for name in valid_names:
        t = _parse_flagstat_total(name)
        total_reads.append(t or 0)

    uniq_counts = [round((u / 100) * t) for u, t in zip(uniq_pcts, total_reads)]
    dup_counts  = [round((d / 100) * t) for d, t in zip(dup_pcts,  total_reads)]

    def _M(vals):
        return [round(v / 1e6, 3) for v in vals]

    traces_pct = [
        {"type": "bar", "name": "Unique",
         "x": uniq_pcts, "y": valid_names, "orientation": "h",
         "marker": {"color": "rgba(30,95,160,0.82)"},
         "hovertemplate": "%{y}<br>Unique: %{x:.2f}%<extra></extra>"},
        {"type": "bar", "name": "Duplicates",
         "x": dup_pcts, "y": valid_names, "orientation": "h",
         "marker": {"color": "rgba(160,160,160,0.75)"},
         "hovertemplate": "%{y}<br>Duplicates: %{x:.2f}%<extra></extra>"},
    ]
    traces_cnt = [
        {"type": "bar", "name": "Unique",
         "x": _M(uniq_counts), "y": valid_names, "orientation": "h",
         "marker": {"color": "rgba(30,95,160,0.82)"},
         "hovertemplate": "%{y}<br>Unique: %{x:.3f}M<extra></extra>"},
        {"type": "bar", "name": "Duplicates",
         "x": _M(dup_counts), "y": valid_names, "orientation": "h",
         "marker": {"color": "rgba(160,160,160,0.75)"},
         "hovertemplate": "%{y}<br>Duplicates: %{x:.3f}M<extra></extra>"},
    ]

    n = len(valid_names)
    h = height if n <= 12 else max(height, n * 28 + 120)

    layout = {
        **_layout(title=title, height=h),
        "barmode":  "stack",
        "xaxis":    {"title": "%", "ticksuffix": "%", "range": [0, 100],
                     "automargin": True},
        "yaxis":    {"automargin": True, "tickfont": {"size": 11}},
        "margin":   {"l": 140, "r": 20, "t": 60, "b": 50},
        "legend":   {"orientation": "h", "y": -0.12},
        "showlegend": True,
        "updatemenus": [{
            "type": "buttons", "direction": "left",
            "x": 0, "y": 1.12, "xanchor": "left",
            "bgcolor": "#f0f0f0", "bordercolor": "#ccc",
            "font": {"size": 11}, "showactive": True,
            "buttons": [
                {"label": "Percentages", "method": "update",
                 "args": [
                     {"x": [t["x"] for t in traces_pct]},
                     {"xaxis": {"title": "%", "ticksuffix": "%",
                                "range": [0, 100], "automargin": True}},
                 ]},
                {"label": "Counts", "method": "update",
                 "args": [
                     {"x": [t["x"] for t in traces_cnt]},
                     {"xaxis": {"title": "Reads (M)", "ticksuffix": "",
                                "range": None, "automargin": True}},
                 ]},
            ],
        }],
    }

    warning = ""
    if missing:
        warning = _missing(
            f"Could not parse {len(missing)} sample(s): "
            f"{', '.join(f'<code>{s}</code>' for s in missing[:5])}"
            + (" …" if len(missing) > 5 else "")
        )

    div_id = _uid("dedup")
    return warning + _render(div_id, {"data": traces_pct, "layout": layout}, height=h)


# ── mbias_tsv_chart ────────────────────────────────────────────────────────────
# When MethylDackel --txt is NOT supported or the old format was used,
# _mbias.txt contains only stderr (OT/OB coords, no per-position data).
# In that case MultiQC also has no M-bias plot data.
#
# Fallback: qc_parser.parse_mbias_txt() saves per-position TSV to
#   {qc_dir}/mbias_data/{name}_mbias.tsv
# Format: sample  strand  read  pos  pct_meth
#
# This function reads those TSVs and renders the same interactive chart
# as mqc_mbias_chart(), with OT/OB × R1/R2 context buttons.
#
# If --txt WAS supported but MultiQC just doesn't have a methyldackel module,
# this is the primary rendering path.

def mbias_tsv_chart(
    mbias_tsv_dir: str | list[str],
    title: str = "M-bias (MethylDackel)",
    height: int = 450,
    tsv_paths: list[str] | None = None,
) -> str:
    """
    Read per-position mbias TSV files saved by qc_parser.parse_mbias_txt()
    and render an interactive Plotly M-bias chart.

    mbias_tsv_dir: one or more mbias_data/ directories.
    tsv_paths: explicit per-sample TSV paths, used when processing artifacts
               are discovered from BAM-backed source projects.
    """
    import glob as _glob, os as _os
    import pandas as _pd

    if isinstance(mbias_tsv_dir, (list, tuple)):
        mbias_dirs = [str(path) for path in mbias_tsv_dir if path]
    else:
        mbias_dirs = [str(mbias_tsv_dir)] if mbias_tsv_dir else []
    tsv_files = [
        path for directory in mbias_dirs
        for path in sorted(_glob.glob(_os.path.join(directory, "*_mbias.tsv")))
    ]
    tsv_files.extend(path for path in (tsv_paths or []) if _os.path.isfile(path))
    tsv_files = list(dict.fromkeys(tsv_files))
    if not tsv_files:
        return _missing(
            "No M-bias TSV files found in the available processing results. "
            f"This directory is populated by <code>cftk qc -s 1</code> "
            f"(which calls <code>parse_mbias_txt()</code>). "
            f"If <code>_mbias.txt</code> only contains "
            f"<em>Suggested inclusion options</em>, the "
            f"<code>MethylDackel mbias --txt</code> flag is not supported by "
            f"your MethylDackel version — see below for the legacy fallback."
        )

    # Load and concatenate all sample TSVs
    dfs = []
    for path in tsv_files:
        try:
            df = _pd.read_csv(path, sep="\t")
            dfs.append(df)
        except Exception:
            continue

    if not dfs:
        return _missing(f"Could not read M-bias TSV files from <code>{mbias_tsv_dir}</code>")

    df = _pd.concat(dfs, ignore_index=True)
    # Expected columns: sample, strand, read, pos, pct_meth

    samples = sorted(df["sample"].unique())
    # Contexts: OT R1, OT R2, OB R1, OB R2 (and optionally CTOT/CTOB)
    contexts = []
    for strand in ["OT", "OB", "CTOT", "CTOB"]:
        for read in sorted(df["read"].unique()):
            sub = df[(df["strand"] == strand) & (df["read"] == read)]
            if not sub.empty:
                contexts.append((f"{strand} R{read}", strand, read))

    if not contexts:
        return _missing("M-bias TSV files are empty or have unexpected format.")

    all_traces: list[dict]   = []
    traces_per_ctx: list[int] = []

    for ci, (label, strand, read) in enumerate(contexts):
        n = 0
        for i, name in enumerate(samples):
            sub = df[(df["sample"] == name) &
                     (df["strand"] == strand) &
                     (df["read"] == read)].sort_values("pos")
            if sub.empty:
                continue
            color = _MQC_COLORS[i % len(_MQC_COLORS)]
            all_traces.append({
                "type":    "scatter",
                "mode":    "lines",
                "name":    name,
                "x":       sub["pos"].tolist(),
                "y":       sub["pct_meth"].tolist(),
                "visible": (ci == 0),
                "line":    {"color": color, "width": 1.5},
                "hovertemplate": (
                    f"{name}<br>Position: %{{x}} bp"
                    f"<br>Methylation: %{{y:.2f}}%<extra></extra>"
                ),
            })
            n += 1
        traces_per_ctx.append(n)

    # Context selector buttons
    buttons = []
    offset  = 0
    for ci, (label, strand, read) in enumerate(contexts):
        n   = traces_per_ctx[ci]
        vis = [False] * len(all_traces)
        for j in range(n):
            if offset + j < len(vis):
                vis[offset + j] = True
        buttons.append({
            "label":  label,
            "method": "update",
            "args":   [{"visible": vis}],
        })
        offset += n

    div_id = _uid("mbias_tsv")
    fig = {
        "data": all_traces,
        "layout": {
            **_layout(title=title,
                      xaxis_title="Position in Read (bp)",
                      yaxis_title="% Methylation",
                      height=height),
            "yaxis":       {"title": "% Methylation", "ticksuffix": "%",
                            "automargin": True},
            "showlegend":  False,
            "updatemenus": [{
                "type":       "buttons",
                "direction":  "left",
                "x": 0, "y": 1.13, "xanchor": "left",
                "bgcolor":    "#f0f0f0", "bordercolor": "#ccc",
                "font":       {"size": 11}, "showactive": True,
                "buttons":    buttons,
            }],
            "annotations": [{
                "text":      "Context:",
                "showarrow": False,
                "x": 0, "y": 1.18,
                "xref": "paper", "yref": "paper",
                "xanchor": "left",
                "font":    {"size": 11},
            }],
        },
    }
    return _render(div_id, fig, height=height)


# ── mbias_legacy_chart ─────────────────────────────────────────────────────────
# When MethylDackel version does NOT support --txt flag, _mbias.txt only contains
# stderr lines like: "Suggested inclusion options: --OT 0,134,0,117 --OB 20,0,24,135"
# In this case there is no per-position data available from any source.
# Render an informative fallback explaining the version issue and showing the
# OT/OB coordinates that were extracted for the extract step.

def mbias_legacy_chart(
    meth_dir: str,
    samples: list[str],
) -> str:
    """
    Render an informative table of OT/OB coordinates from legacy _mbias.txt files
    (those containing only 'Suggested inclusion options:' stderr output).
    """
    import os as _os, re as _re

    rows_html = ""
    n_found   = 0

    for name in samples:
        path = _os.path.join(meth_dir, f"{name}_mbias.txt")
        if not _os.path.exists(path):
            continue
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                text = f.read()
        except Exception:
            continue

        m = _re.search(r"--OT\s+(\S+)\s+--OB\s+(\S+)", text)
        if not m:
            continue

        ot, ob = m.group(1), m.group(2)
        rows_html += (
            f"<tr>"
            f"<td style='padding:6px 12px;font-family:monospace;font-size:12px;"
            f"border-top:1px solid #e5e7eb;'>{name}</td>"
            f"<td style='padding:6px 12px;font-family:monospace;font-size:12px;"
            f"border-top:1px solid #e5e7eb;color:#1e5fa0;'>{ot}</td>"
            f"<td style='padding:6px 12px;font-family:monospace;font-size:12px;"
            f"border-top:1px solid #e5e7eb;color:#2e6b47;'>{ob}</td>"
            f"</tr>"
        )
        n_found += 1

    if not n_found:
        return _missing("No M-bias OT/OB coordinates found in _mbias.txt files.")

    return (
        f'<div style="padding:12px 16px;margin-bottom:14px;'
        f'background:rgba(243,156,18,.07);border-left:3px solid #9a6200;'
        f'border-radius:0 6px 6px 0;font-size:12px;color:var(--ink2);line-height:1.7;">'
        f'<strong style="color:#9a6200;">⚠ Legacy MethylDackel format detected</strong><br>'
        f'Your <code>_mbias.txt</code> files contain only '
        f'<em>Suggested inclusion options</em> (stderr), not per-position TSV data. '
        f'This means <code>MethylDackel mbias --txt</code> is not supported by your '
        f'installed version. Per-position M-bias plots are unavailable.<br>'
        f'The OT/OB inclusion coordinates below were correctly extracted and used '
        f'for the CpG extraction step.</div>'
        f'<div style="overflow-x:auto;border:1px solid #e5e7eb;border-radius:8px;">'
        f'<table style="width:100%;border-collapse:collapse;font-size:13px;">'
        f'<thead><tr>'
        f'<th style="padding:8px 12px;background:#f8f9fa;text-align:left;'
        f'font-size:11px;font-family:monospace;text-transform:uppercase;">Sample</th>'
        f'<th style="padding:8px 12px;background:#f8f9fa;text-align:left;'
        f'font-size:11px;font-family:monospace;text-transform:uppercase;">--OT</th>'
        f'<th style="padding:8px 12px;background:#f8f9fa;text-align:left;'
        f'font-size:11px;font-family:monospace;text-transform:uppercase;">--OB</th>'
        f'</tr></thead>'
        f'<tbody>{rows_html}</tbody>'
        f'</table></div>'
    )
