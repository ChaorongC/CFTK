from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


FORMAT_VERSION = 1
DEFAULT_STATS = ("mean", "CI_l", "CI_u")
REFERENCE_DATA_ENV = "CFTK_MODEL_POWER_DATA"
_REQUIRED_MANIFEST_FIELDS = {
    "depth_labels",
    "format_version",
    "index_file",
    "mean_file",
    "row_count",
    "stats",
    "std_files",
}
_STD_COLUMN_RE = re.compile(
    r"^(?P<depth>[0-9]+(?:\.[0-9]+)?)_(?P<stat>mean|CI_l|CI_u)$"
)
_CPG_INDEX_RE = re.compile(r"^chr(?P<chrom>[^_]+)_(?P<pos>[0-9]+)$")


@dataclass(frozen=True)
class ModelPowerReference:
    """Model-power reference data reconstructed from array files."""

    cpg_mean: pd.Series
    cpg_std_summary: pd.DataFrame
    manifest: dict


def resolve_model_power_reference_dir(
    reference_dir: str | Path | None = None,
    *,
    repository_dir: str | Path | None = None,
) -> Path:
    """Resolve an external model-power reference directory.

    Resolution order is an explicit argument, ``CFTK_MODEL_POWER_DATA``, then
    the supplied repository checkout fallback. Reference arrays are not part
    of the installed CFTK distribution.
    """
    if reference_dir is not None:
        resolved = Path(reference_dir).expanduser()
    else:
        configured = os.environ.get(REFERENCE_DATA_ENV, "").strip()
        if configured:
            resolved = Path(configured).expanduser()
        elif repository_dir is not None:
            resolved = Path(repository_dir).expanduser()
        else:
            resolved = Path(__file__).resolve().parents[2] / "data"

    manifest_path = resolved / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Model-power reference data not found at {resolved}. CFTK installs "
            "the calculator code but does not bundle its reference arrays. "
            "Pass reference_dir='/path/to/CFTK/data', set "
            f"{REFERENCE_DATA_ENV}, or run from a repository checkout that "
            "contains data/manifest.json."
        )
    return resolved


def _depth_label(value: float | int | str) -> str:
    value_float = float(value)
    return str(int(value_float)) if value_float.is_integer() else str(value)


def _depth_manifest_value(label: str) -> int | float:
    value = float(label)
    return int(value) if value.is_integer() else value


def _normalize_depths(depths: Sequence[float | int | str]) -> list[str]:
    labels = [_depth_label(depth) for depth in depths]
    if len(set(labels)) != len(labels):
        raise ValueError("depths contains duplicated values.")
    return labels


def _file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_std_columns(columns: Sequence[object]) -> tuple[list[str], dict[str, dict[str, str]]]:
    parsed: dict[str, dict[str, str]] = {}
    for column in columns:
        column_label = str(column)
        match = _STD_COLUMN_RE.match(column_label)
        if match is None:
            raise ValueError(
                f"Unexpected std column {column_label!r}; expected '<depth>_<stat>'."
            )
        depth = _depth_label(match.group("depth"))
        stat = match.group("stat")
        parsed.setdefault(depth, {})[stat] = column_label

    depths = sorted(parsed, key=lambda label: float(label))
    for depth in depths:
        missing = [stat for stat in DEFAULT_STATS if stat not in parsed[depth]]
        if missing:
            raise ValueError(
                f"Depth {depth} is missing std columns for: {', '.join(missing)}."
            )
    return depths, parsed


def _split_cpg_index(index: pd.Index) -> tuple[np.ndarray, np.ndarray]:
    chrom = np.empty(len(index), dtype=object)
    pos = np.empty(len(index), dtype=np.uint32)

    for row_idx, value in enumerate(index.astype(str)):
        match = _CPG_INDEX_RE.match(value)
        if match is None:
            raise ValueError(
                f"CpG index value {value!r} does not match expected chr<chrom>_<pos>."
            )
        position = int(match.group("pos"))
        if position > np.iinfo(np.uint32).max:
            raise ValueError(f"CpG index position is too large for uint32: {value!r}.")
        chrom[row_idx] = match.group("chrom")
        pos[row_idx] = position

    return chrom.astype(str), pos


def _join_cpg_index(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    chrom = chrom.astype(str)
    pos = pos.astype(str)
    return np.char.add(np.char.add("chr", chrom), np.char.add("_", pos))


def _read_mean(path: Path) -> pd.Series:
    mean = pd.read_pickle(path)
    if isinstance(mean, pd.DataFrame):
        if mean.shape[1] != 1:
            raise ValueError("Mean pickle DataFrame must have exactly one column.")
        mean = mean.iloc[:, 0]
    if not isinstance(mean, pd.Series):
        raise TypeError("Mean pickle must contain a Series or one-column DataFrame.")
    if not mean.index.is_unique:
        raise ValueError("Mean index contains duplicated CpG IDs.")
    return pd.to_numeric(mean, errors="raise")


def _read_std(path: Path) -> pd.DataFrame:
    std = pd.read_pickle(path)
    if not isinstance(std, pd.DataFrame):
        raise TypeError("STD pickle must contain a DataFrame.")
    if not std.index.is_unique:
        raise ValueError("STD index contains duplicated CpG IDs.")
    return std.apply(pd.to_numeric, errors="raise")


def convert_pickles_to_array_reference(
    mean_pickle: str | Path,
    std_pickle: str | Path,
    output_dir: str | Path,
    *,
    dtype: str | np.dtype = "float32",
    overwrite: bool = False,
) -> dict:
    """
    Convert existing model-power pickles to shared-index NumPy arrays.

    The output layout is:
      cpg_index.npz
      cpg_mean.<dtype>.npy
      std_by_depth/depth_<depth>_mean.<dtype>.npy
      std_by_depth/depth_<depth>_CI.<dtype>.npy
      manifest.json
    """
    mean_pickle = Path(mean_pickle)
    std_pickle = Path(std_pickle)
    output_dir = Path(output_dir)
    array_dtype = np.dtype(dtype)

    if array_dtype.kind not in {"f"}:
        raise ValueError("dtype must be a floating-point dtype.")

    manifest_path = output_dir / "manifest.json"
    if output_dir.exists() and manifest_path.exists() and not overwrite:
        raise FileExistsError(
            f"{output_dir} already contains manifest.json; pass overwrite=True."
        )

    mean = _read_mean(mean_pickle)
    std = _read_std(std_pickle)

    mean_index = mean.index.astype(str)
    std_index = std.index.astype(str)
    if not mean_index.equals(std_index):
        raise ValueError(
            "Mean and STD pickles must have identical row order and CpG IDs."
        )

    depths, parsed_columns = _parse_std_columns(std.columns)

    output_dir.mkdir(parents=True, exist_ok=True)
    std_dir = output_dir / "std_by_depth"
    std_dir.mkdir(parents=True, exist_ok=True)

    cpg_chrom, cpg_pos = _split_cpg_index(mean_index)
    mean_array = mean.to_numpy(dtype=array_dtype, copy=True)

    np.savez_compressed(output_dir / "cpg_index.npz", chrom=cpg_chrom, pos=cpg_pos)
    np.save(output_dir / f"cpg_mean.{array_dtype.name}.npy", mean_array)

    std_files = {}
    for depth in depths:
        mean_column = parsed_columns[depth]["mean"]
        ci_columns = [parsed_columns[depth]["CI_l"], parsed_columns[depth]["CI_u"]]

        mean_file = f"depth_{depth}_mean.{array_dtype.name}.npy"
        ci_file = f"depth_{depth}_CI.{array_dtype.name}.npy"

        mean_matrix = std.loc[:, mean_column].to_numpy(dtype=array_dtype, copy=True)
        ci_matrix = std.loc[:, ci_columns].to_numpy(dtype=array_dtype, copy=True)

        np.save(std_dir / mean_file, mean_matrix)
        np.save(std_dir / ci_file, ci_matrix)

        std_files[depth] = {
            "mean": f"std_by_depth/{mean_file}",
            "CI": f"std_by_depth/{ci_file}",
        }

    manifest = {
        "format_version": FORMAT_VERSION,
        "row_count": int(len(cpg_chrom)),
        "depths": [_depth_manifest_value(depth) for depth in depths],
        "depth_labels": depths,
        "stats": list(DEFAULT_STATS),
        "dtype": array_dtype.name,
        "index_file": "cpg_index.npz",
        "index_format": "chrom_pos",
        "mean_file": f"cpg_mean.{array_dtype.name}.npy",
        "ci_columns": ["CI_l", "CI_u"],
        "std_files": std_files,
        "source_files": {
            "mean_pickle": str(mean_pickle),
            "std_pickle": str(std_pickle),
        },
        "source_sha256": {
            "mean_pickle": _file_sha256(mean_pickle),
            "std_pickle": _file_sha256(std_pickle),
        },
    }

    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def load_manifest(reference_dir: str | Path) -> dict:
    reference_dir = Path(reference_dir)
    manifest_path = reference_dir / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Model-power reference manifest not found: {manifest_path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Model-power reference manifest is not valid JSON: {manifest_path}"
        ) from exc

    if not isinstance(manifest, dict):
        raise ValueError(
            f"Model-power reference manifest must contain a JSON object: "
            f"{manifest_path}"
        )

    if manifest.get("format_version") != FORMAT_VERSION:
        raise ValueError(
            f"Unsupported model-power reference format_version: "
            f"{manifest.get('format_version')!r}."
        )

    missing_fields = sorted(_REQUIRED_MANIFEST_FIELDS - set(manifest))
    if missing_fields:
        raise ValueError(
            f"Model-power reference manifest is missing required fields: "
            f"{', '.join(missing_fields)}."
        )

    try:
        row_count = int(manifest["row_count"])
    except (TypeError, ValueError) as exc:
        raise ValueError("Model-power manifest row_count must be an integer.") from exc
    if row_count < 1:
        raise ValueError("Model-power manifest row_count must be positive.")

    depth_labels = [str(value) for value in manifest["depth_labels"]]
    if not depth_labels or len(set(depth_labels)) != len(depth_labels):
        raise ValueError(
            "Model-power manifest depth_labels must contain unique values."
        )
    if not isinstance(manifest["stats"], list) or not manifest["stats"]:
        raise ValueError("Model-power manifest stats must be a non-empty list.")
    if not isinstance(manifest["std_files"], dict):
        raise ValueError("Model-power manifest std_files must be an object.")
    missing_depth_files = [
        depth for depth in depth_labels if depth not in manifest["std_files"]
    ]
    if missing_depth_files:
        raise ValueError(
            "Model-power manifest std_files is missing depths: "
            + ", ".join(missing_depth_files)
        )
    return manifest


def _validate_reference_request(
    manifest: dict,
    *,
    depths: Sequence[float | int | str] | None,
    sd_stats: Sequence[str],
    include_index: bool,
) -> tuple[list[str], list[str], list[str]]:
    available_depths = [str(depth) for depth in manifest["depth_labels"]]
    requested_depths = (
        available_depths if depths is None else _normalize_depths(depths)
    )
    missing_depths = sorted(set(requested_depths) - set(available_depths), key=float)
    if missing_depths:
        raise ValueError(
            "Depths not available in model-power reference: "
            + ", ".join(missing_depths)
        )

    requested_stats = list(sd_stats)
    missing_stats = sorted(set(requested_stats) - set(manifest["stats"]))
    if missing_stats:
        raise ValueError(
            "STD stats not available in model-power reference: "
            + ", ".join(missing_stats)
        )

    required_files = [str(manifest["mean_file"])]
    if include_index:
        required_files.append(str(manifest["index_file"]))
    for depth in requested_depths:
        depth_files = manifest["std_files"][depth]
        if "mean" in requested_stats:
            try:
                required_files.append(str(depth_files["mean"]))
            except (KeyError, TypeError) as exc:
                raise ValueError(
                    f"Model-power manifest has no mean array for depth {depth}."
                ) from exc
        if {"CI_l", "CI_u"}.intersection(requested_stats):
            try:
                required_files.append(str(depth_files["CI"]))
            except (KeyError, TypeError) as exc:
                raise ValueError(
                    f"Model-power manifest has no CI array for depth {depth}."
                ) from exc
    return requested_depths, requested_stats, required_files


def validate_model_power_reference(
    reference_dir: str | Path,
    *,
    depths: Sequence[float | int | str] | None = None,
    sd_stats: Sequence[str] = DEFAULT_STATS,
    include_index: bool = True,
) -> dict:
    """Validate the manifest and files required for a reference-data request."""
    reference_dir = Path(reference_dir)
    manifest = load_manifest(reference_dir)
    _, _, required_files = _validate_reference_request(
        manifest,
        depths=depths,
        sd_stats=sd_stats,
        include_index=include_index,
    )
    missing = [
        reference_dir / relative_path
        for relative_path in dict.fromkeys(required_files)
        if not (reference_dir / relative_path).is_file()
    ]
    if missing:
        missing_text = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(
            "Model-power reference data is incomplete. Missing required files:\n"
            f"{missing_text}"
        )
    return manifest


def load_model_power_reference(
    reference_dir: str | Path,
    *,
    depths: Sequence[float | int | str] | None = None,
    sd_stats: Sequence[str] = DEFAULT_STATS,
    include_index: bool = True,
    mmap_mode: str | None = "r",
) -> ModelPowerReference:
    """Load selected model-power reference depths from the array layout."""
    reference_dir = Path(reference_dir)
    manifest = validate_model_power_reference(
        reference_dir,
        depths=depths,
        sd_stats=sd_stats,
        include_index=include_index,
    )
    requested_depths, requested_stats, _ = _validate_reference_request(
        manifest,
        depths=depths,
        sd_stats=sd_stats,
        include_index=include_index,
    )

    row_count = int(manifest["row_count"])
    mean_array = np.load(
        reference_dir / manifest["mean_file"],
        mmap_mode=mmap_mode,
        allow_pickle=False,
    )
    if mean_array.shape != (row_count,):
        raise ValueError("Mean array shape does not match manifest row_count.")

    output_columns = {}
    for depth in requested_depths:
        depth_files = manifest["std_files"][depth]

        if "mean" in requested_stats:
            depth_mean = np.load(
                reference_dir / depth_files["mean"],
                mmap_mode=mmap_mode,
                allow_pickle=False,
            )
            if depth_mean.shape != (row_count,):
                raise ValueError(
                    f"Mean STD array for depth {depth} has shape "
                    f"{depth_mean.shape}, expected {(row_count,)}."
                )
            output_columns[f"{depth}_mean"] = depth_mean

        ci_stats = [stat for stat in requested_stats if stat in {"CI_l", "CI_u"}]
        if ci_stats:
            ci_matrix = np.load(
                reference_dir / depth_files["CI"],
                mmap_mode=mmap_mode,
                allow_pickle=False,
            )
            expected_shape = (row_count, 2)
            if ci_matrix.shape != expected_shape:
                raise ValueError(
                    f"CI STD array for depth {depth} has shape {ci_matrix.shape}, "
                    f"expected {expected_shape}."
                )
            ci_positions = {"CI_l": 0, "CI_u": 1}
            for stat in ci_stats:
                output_columns[f"{depth}_{stat}"] = ci_matrix[:, ci_positions[stat]]

    if include_index:
        index_path = reference_dir / manifest["index_file"]
        cpg_index_data = np.load(index_path, allow_pickle=False)
        if isinstance(cpg_index_data, np.lib.npyio.NpzFile):
            with cpg_index_data:
                if {"chrom", "pos"}.issubset(cpg_index_data.files):
                    cpg_index = _join_cpg_index(
                        cpg_index_data["chrom"],
                        cpg_index_data["pos"],
                    )
                else:
                    cpg_index = cpg_index_data["cpg_index"]
        else:
            cpg_index = cpg_index_data
        if len(cpg_index) != row_count:
            raise ValueError("CpG index row count does not match manifest.")
        index = pd.Index(cpg_index.astype(str), dtype=object)
    else:
        index = pd.RangeIndex(row_count)

    cpg_mean = pd.Series(mean_array, index=index, name="baseline_mean")
    cpg_std_summary = pd.DataFrame(output_columns, index=index)

    return ModelPowerReference(
        cpg_mean=cpg_mean,
        cpg_std_summary=cpg_std_summary,
        manifest=manifest,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert CFTK model-power reference pickles to NumPy arrays."
    )
    parser.add_argument("--mean-pickle", required=True, type=Path)
    parser.add_argument("--std-pickle", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--dtype", default="float32", choices=["float32", "float64"])
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    manifest = convert_pickles_to_array_reference(
        args.mean_pickle,
        args.std_pickle,
        args.output_dir,
        dtype=args.dtype,
        overwrite=args.overwrite,
    )
    print(
        f"Wrote model-power reference arrays to {args.output_dir} "
        f"({manifest['row_count']} CpGs, {len(manifest['depths'])} depths)."
    )


if __name__ == "__main__":
    main()
