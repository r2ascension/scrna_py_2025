#!/usr/bin/env python3
from __future__ import annotations

import copy
import importlib.util
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np
import pandas as pd


RUN_ROOT_DEFAULT = Path("/home/h2048/output/program_full_parallel_methods_20260507")
HELPER_PATH = Path("/home/h2048/script/py/cnmf_helper_20260419_v1_1.py")
LEVEL_CONFIG = {
    "l2": {
        "subdir": "cnmf_by_celltype_l2",
        "status_file": "cnmf_l2_status.json",
        "meta_file": "cell_metadata_for_l2_cnmf.csv",
        "celltype_status_keys": ("celltype_l2", "cell_type_L2"),
    },
    "l3": {
        "subdir": "cnmf_by_celltype",
        "status_file": "cnmf_l3_status.json",
        "meta_file": "cell_metadata_for_l3_cnmf.csv",
        "celltype_status_keys": ("celltype_l3", "cell_type_L3"),
        "parent_status_keys": ("celltype_l2", "cell_type_L2"),
    },
}
DEFAULT_GROUP_COLUMNS = ["tissue", "sample", "dataset", "condition"]
DEFAULT_ROE_THRESHOLDS = (0.5, 0.8, 1.2, 2.0)


def env_flag(name: str, default: bool = False) -> bool:
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "y", "on"}


def split_env_csv(name: str, default: Optional[Sequence[str]] = None) -> List[str]:
    raw = os.getenv(name, "").strip()
    if not raw:
        return list(default or [])
    return [part.strip() for part in raw.split(",") if part.strip()]


def natural_sort_key(value: object) -> List[object]:
    parts = re.split(r"(\d+)", str(value))
    key: List[object] = []
    for part in parts:
        if not part:
            continue
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part.lower()))
    return key


def sort_gep_labels(labels: Iterable[object]) -> List[str]:
    def key(label: object) -> Tuple[str, int, str]:
        text = str(label)
        if "__GEP_" in text:
            prefix, suffix = text.rsplit("__GEP_", 1)
            if suffix.isdigit():
                return (prefix.lower(), int(suffix), text)
        if text.startswith("GEP_") and text[4:].isdigit():
            return ("", int(text[4:]), text)
        return (text.lower(), math.inf, text)

    return [str(label) for label in sorted(labels, key=key)]


def sanitize_group_series(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    invalid = out.isna() | out.isin(["", "NA", "N/A", "None", "nan", "<NA>"])
    return out.mask(invalid)


def roe_to_symbol(value: float, thresholds: Sequence[float] = DEFAULT_ROE_THRESHOLDS) -> str:
    if value is None or pd.isna(value):
        return "-"
    if np.isinf(value):
        return "+++" if value > 0 else "---"
    if value >= thresholds[3]:
        return "+++"
    if value >= thresholds[2]:
        return "+"
    if value > thresholds[1]:
        return "+/-"
    if value >= thresholds[0]:
        return "-"
    return "---"


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def save_table(df: pd.DataFrame, path: Path) -> None:
    ensure_parent(path)
    df.to_csv(path, sep="\t")


def load_helper_module():
    if not HELPER_PATH.exists():
        raise FileNotFoundError(f"Missing helper module: {HELPER_PATH}")
    spec = importlib.util.spec_from_file_location("cnmf_helper_20260419_v1_1", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"Failed to load helper module from {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


HELPER = load_helper_module()
load_usage_matrix = HELPER.load_usage_matrix


def discover_k_values(source_dir: Path, run_name: str) -> List[int]:
    ks = set()
    gene_dir = source_dir / "gep_gene_tables"
    for path in gene_dir.glob("gep_gene_scores_k*.tsv"):
        match = re.search(r"k(\d+)", path.name)
        if match:
            ks.add(int(match.group(1)))
    if ks:
        return sorted(ks)

    run_dir = source_dir / "cnmf_output" / run_name
    for path in run_dir.glob(f"{run_name}.gene_spectra_score.k_*.txt"):
        match = re.search(r"\.k_(\d+)\.", path.name)
        if match:
            ks.add(int(match.group(1)))
    return sorted(ks)


def choose_status_value(status: Dict[str, object], keys: Sequence[str]) -> Optional[str]:
    for key in keys:
        value = status.get(key)
        if value is not None and str(value).strip():
            return str(value)
    return None


def discover_sources(
    run_root: Path,
    level: str,
    lineages: Optional[Sequence[str]] = None,
    celltypes: Optional[Sequence[str]] = None,
) -> List[Dict[str, object]]:
    cfg = LEVEL_CONFIG[level]
    lineage_filter = set(lineages or [])
    celltype_filter = set(celltypes or [])
    discovered: List[Dict[str, object]] = []

    for lineage_dir in sorted(run_root.iterdir() if run_root.exists() else [], key=lambda p: natural_sort_key(p.name)):
        if not lineage_dir.is_dir():
            continue
        lineage = lineage_dir.name
        if lineage_filter and lineage not in lineage_filter:
            continue

        base_dir = lineage_dir / cfg["subdir"]
        if not base_dir.exists():
            continue

        candidate_dirs = list(base_dir.glob("*/cnmf_full")) + list(base_dir.glob("cnmf_*/cnmf_full"))
        dedup: Dict[Tuple[str, str], Dict[str, object]] = {}

        for source_dir in sorted(candidate_dirs, key=lambda p: natural_sort_key(str(p))):
            status_path = source_dir / cfg["status_file"]
            meta_path = source_dir / cfg["meta_file"]
            if not status_path.exists() or not meta_path.exists():
                continue
            try:
                status = json.loads(status_path.read_text())
            except Exception:
                continue
            if str(status.get("status", "")).lower() != "ok":
                continue

            celltype = choose_status_value(status, cfg["celltype_status_keys"]) or source_dir.parent.name.removeprefix("cnmf_")
            if celltype_filter and celltype not in celltype_filter:
                continue
            run_name = status.get("run_name")
            if not run_name:
                run_dirs = [p.name for p in (source_dir / "cnmf_output").iterdir() if p.is_dir()] if (source_dir / "cnmf_output").exists() else []
                if len(run_dirs) == 1:
                    run_name = run_dirs[0]
            if not run_name:
                continue

            source_label = celltype
            if level == "l3":
                parent = choose_status_value(status, cfg.get("parent_status_keys", ()))
                if parent and parent != celltype:
                    source_label = f"{parent}>{celltype}"

            key = (lineage, celltype)
            is_canonical = source_dir.parent.name == celltype
            record = {
                "level": level,
                "lineage": lineage,
                "celltype": celltype,
                "source_label": source_label,
                "source_dir": source_dir,
                "status_path": status_path,
                "meta_path": meta_path,
                "status": status,
                "run_name": str(run_name),
                "cnmf_run_dir": source_dir / "cnmf_output" / str(run_name),
                "is_canonical": is_canonical,
            }
            current = dedup.get(key)
            if current is None or (record["is_canonical"] and not current["is_canonical"]):
                dedup[key] = record

        discovered.extend(dedup.values())

    return sorted(discovered, key=lambda r: (str(r["lineage"]), str(r["source_label"])))


def load_metadata(meta_path: Path) -> pd.DataFrame:
    meta = pd.read_csv(meta_path, sep=None, engine="python")
    if "cell" not in meta.columns:
        raise ValueError(f"Metadata missing required 'cell' column: {meta_path}")
    meta["cell"] = meta["cell"].astype(str)
    meta = meta.drop_duplicates(subset=["cell"], keep="first")
    return meta


def build_assignments(
    record: Dict[str, object],
    metadata: pd.DataFrame,
    usage_df: pd.DataFrame,
    k: int,
) -> Optional[pd.DataFrame]:
    if usage_df is None or usage_df.empty:
        return None

    usage = usage_df.copy()
    usage.index = usage.index.astype(str)
    usage = usage.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    meta_indexed = metadata.set_index("cell", drop=False)
    common = meta_indexed.index.intersection(usage.index)
    if len(common) == 0:
        return None

    meta_indexed = meta_indexed.loc[common].copy()
    usage = usage.loc[common].copy()
    total_usage = usage.sum(axis=1)
    positive = total_usage > 0

    dominant_gep = pd.Series(pd.NA, index=usage.index, dtype="object")
    dominant_usage = pd.Series(np.nan, index=usage.index, dtype="float64")
    if positive.any():
        dominant_gep.loc[positive] = usage.loc[positive].idxmax(axis=1).astype(str)
        dominant_usage.loc[positive] = usage.loc[positive].max(axis=1).astype(float)

    usage_norm = usage.div(total_usage.where(positive, np.nan), axis=0)
    usage_norm = usage_norm.clip(lower=1e-12)
    gep_entropy = -(usage_norm * np.log2(usage_norm)).sum(axis=1)
    gep_entropy.loc[~positive] = np.nan

    assignments = meta_indexed.reset_index(drop=True).copy()
    assignments["k"] = int(k)
    assignments["level"] = str(record["level"])
    assignments["lineage"] = str(record["lineage"])
    assignments["source_celltype"] = str(record["celltype"])
    assignments["source_label"] = str(record["source_label"])
    assignments["dominant_gep"] = dominant_gep.to_numpy()
    assignments["dominant_usage"] = dominant_usage.to_numpy()
    assignments["total_usage"] = total_usage.to_numpy(dtype=float)
    assignments["gep_entropy"] = gep_entropy.to_numpy(dtype=float)
    assignments["gep_celltype"] = pd.array(
        [
            f"{record['source_label']}__{gep}" if isinstance(gep, str) and gep else pd.NA
            for gep in assignments["dominant_gep"].tolist()
        ],
        dtype="string",
    )
    return assignments


def build_dominant_group_tables(
    assignments: pd.DataFrame,
    row_key: str,
    group_col: str,
    min_row_sum: int,
) -> Optional[Dict[str, pd.DataFrame]]:
    if row_key not in assignments.columns or group_col not in assignments.columns:
        return None

    work = assignments[[row_key, group_col]].copy()
    work = work[work[row_key].notna()].copy()
    work[group_col] = sanitize_group_series(work[group_col])
    work = work.dropna(subset=[group_col])
    if work.empty or work[group_col].nunique() < 2:
        return None

    observed = pd.crosstab(work[row_key].astype(str), work[group_col].astype(str))
    observed = observed.loc[observed.sum(axis=1) >= min_row_sum].copy()
    if observed.empty or observed.shape[1] < 2:
        return None

    observed = observed.loc[sort_gep_labels(observed.index), sorted(observed.columns, key=natural_sort_key)]
    fraction = observed.div(observed.sum(axis=0), axis=1).fillna(0.0)
    expected = pd.DataFrame(
        np.outer(observed.sum(axis=1), observed.sum(axis=0)) / float(observed.values.sum()),
        index=observed.index,
        columns=observed.columns,
    )
    roe = observed / expected
    symbol = roe.applymap(roe_to_symbol)
    summary = pd.DataFrame(
        {
            row_key: observed.index,
            f"n_{group_col}_nonzero": (observed > 0).sum(axis=1).to_numpy(dtype=int),
            "n_cells": observed.sum(axis=1).to_numpy(dtype=int),
            f"dominant_{group_col}": observed.idxmax(axis=1).to_numpy(),
            f"dominant_{group_col}_roe": roe.max(axis=1).to_numpy(dtype=float),
            f"dominant_{group_col}_fraction": fraction.max(axis=1).to_numpy(dtype=float),
        }
    )
    return {
        "observed": observed,
        "fraction": fraction,
        "expected": expected,
        "roe": roe,
        "symbol": symbol,
        "summary": summary,
    }


def build_mean_usage_by_group(
    assignments: pd.DataFrame,
    usage_df: pd.DataFrame,
    group_col: str,
    row_prefix: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    if group_col not in assignments.columns:
        return None

    work = assignments[["cell", group_col]].copy()
    work[group_col] = sanitize_group_series(work[group_col])
    work = work.dropna(subset=[group_col])
    if work.empty or work[group_col].nunique() < 2:
        return None

    usage = usage_df.copy()
    usage.index = usage.index.astype(str)
    usage = usage.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    joined = work.drop_duplicates(subset=["cell"]).set_index("cell").join(usage, how="inner")
    if joined.empty:
        return None

    grouped = joined.groupby(group_col)[list(usage.columns)].mean().T
    grouped = grouped.loc[sort_gep_labels(grouped.index), sorted(grouped.columns, key=natural_sort_key)]
    if row_prefix:
        grouped.index = [f"{row_prefix}__{idx}" for idx in grouped.index]
    return grouped


def maybe_copy_cmap(name: str):
    cmap = plt.get_cmap(name)
    try:
        cmap = cmap.copy()
    except AttributeError:
        cmap = copy.copy(cmap)
    cmap.set_bad("#f0f0f0")
    return cmap


def plot_numeric_heatmap(
    df: pd.DataFrame,
    output_prefix: Path,
    title: str,
    cmap_name: str,
    cbar_label: str,
    annotate: bool = True,
    annotation_fmt: str = "{:.2f}",
) -> None:
    if df is None or df.empty:
        return

    values = df.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)
    n_rows, n_cols = df.shape
    fig_width = max(6.0, 0.75 * n_cols + 2.5)
    fig_height = max(4.0, 0.38 * n_rows + 2.0)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(masked, aspect="auto", cmap=maybe_copy_cmap(cmap_name), interpolation="nearest")
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(df.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(df.index)
    ax.set_title(title)
    ax.set_xlabel(df.columns.name or "group")
    ax.set_ylabel(df.index.name or "GEP")

    if annotate and n_rows * n_cols <= 180:
        for i in range(n_rows):
            for j in range(n_cols):
                if np.isfinite(values[i, j]):
                    ax.text(j, i, annotation_fmt.format(values[i, j]), ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(image, ax=ax, shrink=0.85)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    ensure_parent(output_prefix.with_suffix(".png"))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_roe_heatmap(
    roe_df: pd.DataFrame,
    symbol_df: pd.DataFrame,
    output_prefix: Path,
    title: str,
) -> None:
    if roe_df is None or roe_df.empty:
        return

    display = np.log2(roe_df.clip(lower=1e-6))
    values = display.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)
    bound = max(0.5, float(np.nanmax(np.abs(values)))) if np.isfinite(values).any() else 1.0
    norm = mcolors.TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)

    n_rows, n_cols = roe_df.shape
    fig_width = max(6.0, 0.75 * n_cols + 2.5)
    fig_height = max(4.0, 0.38 * n_rows + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    image = ax.imshow(masked, aspect="auto", cmap=maybe_copy_cmap("RdBu_r"), norm=norm, interpolation="nearest")
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(roe_df.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(roe_df.index)
    ax.set_title(title)
    ax.set_xlabel(roe_df.columns.name or "group")
    ax.set_ylabel(roe_df.index.name or "GEP")

    if n_rows * n_cols <= 220:
        for i in range(n_rows):
            for j in range(n_cols):
                ax.text(j, i, str(symbol_df.iat[i, j]), ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(image, ax=ax, shrink=0.85)
    cbar.set_label("log2(Ro/e)")
    fig.tight_layout()
    ensure_parent(output_prefix.with_suffix(".png"))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def export_group_bundle(
    output_dir: Path,
    prefix: str,
    group_col: str,
    dominant_tables: Optional[Dict[str, pd.DataFrame]],
    mean_usage: Optional[pd.DataFrame],
    max_plot_columns: int,
) -> Dict[str, object]:
    manifest: Dict[str, object] = {
        "group_col": group_col,
        "tables": {},
        "plots": {},
        "skipped": False,
        "reason": None,
    }

    if dominant_tables is None and mean_usage is None:
        manifest["skipped"] = True
        manifest["reason"] = "not_enough_group_variation"
        return manifest

    if dominant_tables is not None:
        for key in ["observed", "fraction", "expected", "roe", "symbol", "summary"]:
            path = output_dir / f"{prefix}_{key}.tsv"
            save_table(dominant_tables[key], path)
            manifest["tables"][key] = str(path)

        if dominant_tables["observed"].shape[1] <= max_plot_columns:
            frac_prefix = output_dir / f"{prefix}_fraction_heatmap"
            plot_numeric_heatmap(
                dominant_tables["fraction"],
                frac_prefix,
                title=f"{prefix.replace('_', ' ')} fraction",
                cmap_name="viridis",
                cbar_label="fraction within group",
            )
            manifest["plots"]["fraction_heatmap"] = str(frac_prefix.with_suffix(".png"))

            roe_prefix = output_dir / f"{prefix}_roe_heatmap"
            plot_roe_heatmap(
                dominant_tables["roe"],
                dominant_tables["symbol"],
                roe_prefix,
                title=f"{prefix.replace('_', ' ')} Ro/e",
            )
            manifest["plots"]["roe_heatmap"] = str(roe_prefix.with_suffix(".png"))
        else:
            manifest["plots"]["fraction_heatmap"] = None
            manifest["plots"]["roe_heatmap"] = None
            manifest["reason"] = f"too_many_{group_col}_levels_for_plot"

    if mean_usage is not None:
        mean_path = output_dir / f"{prefix}_mean_usage.tsv"
        save_table(mean_usage, mean_path)
        manifest["tables"]["mean_usage"] = str(mean_path)
        if mean_usage.shape[1] <= max_plot_columns:
            mean_prefix = output_dir / f"{prefix}_mean_usage_heatmap"
            plot_numeric_heatmap(
                mean_usage,
                mean_prefix,
                title=f"{prefix.replace('_', ' ')} mean usage",
                cmap_name="magma",
                cbar_label="mean usage",
            )
            manifest["plots"]["mean_usage_heatmap"] = str(mean_prefix.with_suffix(".png"))
        else:
            manifest["plots"]["mean_usage_heatmap"] = None
            if manifest["reason"] is None:
                manifest["reason"] = f"too_many_{group_col}_levels_for_mean_usage_plot"

    return manifest


def process_source(
    record: Dict[str, object],
    group_columns: Sequence[str],
    min_row_sum: int,
    source_plot_max_cols: int,
) -> Tuple[List[Dict[str, object]], Dict[int, Dict[str, object]]]:
    metadata = load_metadata(Path(record["meta_path"]))
    ks = discover_k_values(Path(record["source_dir"]), str(record["run_name"]))
    source_rows: List[Dict[str, object]] = []
    aggregate_payloads: Dict[int, Dict[str, object]] = {}

    for k in ks:
        usage_df = load_usage_matrix(Path(record["cnmf_run_dir"]), str(record["run_name"]), int(k))
        assignments = build_assignments(record, metadata, usage_df, int(k)) if usage_df is not None else None
        if assignments is None or assignments.empty:
            continue

        source_out_dir = Path(record["source_dir"]) / "gep_celltype_derived" / f"k{k}"
        source_out_dir.mkdir(parents=True, exist_ok=True)
        assignments_path = source_out_dir / f"cell_gep_assignments_k{k}.tsv"
        save_table(assignments, assignments_path)

        group_manifests: Dict[str, object] = {}
        prefixed_mean_usage_by_group: Dict[str, pd.DataFrame] = {}
        for group_col in group_columns:
            dominant_tables = build_dominant_group_tables(assignments, "dominant_gep", group_col, min_row_sum)
            mean_usage = build_mean_usage_by_group(assignments, usage_df, group_col)
            if mean_usage is not None:
                prefixed = mean_usage.copy()
                prefixed.index = [f"{record['source_label']}__{idx}" for idx in prefixed.index]
                prefixed_mean_usage_by_group[group_col] = prefixed
            bundle_prefix = f"dominant_gep_by_{group_col}"
            group_manifests[group_col] = export_group_bundle(
                source_out_dir,
                bundle_prefix,
                group_col,
                dominant_tables,
                mean_usage,
                max_plot_columns=source_plot_max_cols,
            )

        manifest = {
            "level": record["level"],
            "lineage": record["lineage"],
            "source_celltype": record["celltype"],
            "source_label": record["source_label"],
            "run_name": record["run_name"],
            "k": int(k),
            "n_cells": int(len(assignments)),
            "n_geps": int(usage_df.shape[1]),
            "assignments_path": str(assignments_path),
            "groups": group_manifests,
        }
        manifest_path = source_out_dir / "manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))

        source_rows.append(
            {
                "level": record["level"],
                "lineage": record["lineage"],
                "source_celltype": record["celltype"],
                "source_label": record["source_label"],
                "k": int(k),
                "n_cells": int(len(assignments)),
                "assignments_path": str(assignments_path),
                "manifest_path": str(manifest_path),
            }
        )
        aggregate_payloads[int(k)] = {
            "assignments": assignments,
            "mean_usage_by_group": prefixed_mean_usage_by_group,
            "source_label": record["source_label"],
        }

    return source_rows, aggregate_payloads


def write_aggregate_outputs(
    run_root: Path,
    level: str,
    k: int,
    payloads: Sequence[Dict[str, object]],
    group_columns: Sequence[str],
    min_row_sum: int,
    aggregate_plot_max_cols: int,
) -> Optional[Path]:
    if not payloads:
        return None

    out_dir = run_root / "cnmf_gep_celltype_derived" / level / f"k{k}"
    out_dir.mkdir(parents=True, exist_ok=True)
    combined_assignments = pd.concat([payload["assignments"] for payload in payloads], ignore_index=True)
    combined_path = out_dir / f"{level}_gep_celltype_assignments_k{k}.tsv"
    save_table(combined_assignments, combined_path)

    group_manifests: Dict[str, object] = {}
    for group_col in group_columns:
        dominant_tables = build_dominant_group_tables(combined_assignments, "gep_celltype", group_col, min_row_sum)
        mean_frames = [
            payload["mean_usage_by_group"].get(group_col)
            for payload in payloads
            if payload["mean_usage_by_group"].get(group_col) is not None
        ]
        mean_usage = None
        if mean_frames:
            mean_usage = pd.concat(mean_frames, axis=0, sort=False)
            mean_usage = mean_usage.loc[sort_gep_labels(mean_usage.index), sorted(mean_usage.columns, key=natural_sort_key)]
        group_manifests[group_col] = export_group_bundle(
            out_dir,
            f"gep_celltype_by_{group_col}",
            group_col,
            dominant_tables,
            mean_usage,
            max_plot_columns=aggregate_plot_max_cols,
        )

    manifest = {
        "level": level,
        "k": int(k),
        "n_sources": int(len(payloads)),
        "n_cells": int(len(combined_assignments)),
        "assignments_path": str(combined_path),
        "groups": group_manifests,
        "source_labels": [str(payload["source_label"]) for payload in payloads],
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    return manifest_path


def main() -> int:
    run_root = Path(os.getenv("CNMF_GEP_PLOT_RUN_ROOT", str(RUN_ROOT_DEFAULT)))
    levels = split_env_csv("CNMF_GEP_PLOT_LEVELS", default=["l2"])
    lineages = split_env_csv("CNMF_GEP_PLOT_LINEAGES")
    celltypes = split_env_csv("CNMF_GEP_PLOT_CELLTYPES")
    group_columns = split_env_csv("CNMF_GEP_PLOT_GROUP_COLUMNS", default=DEFAULT_GROUP_COLUMNS)
    min_row_sum = int(os.getenv("CNMF_GEP_ROE_MIN_ROW_SUM", "10"))
    source_plot_max_cols = int(os.getenv("CNMF_GEP_SOURCE_MAX_COLUMNS", "18"))
    aggregate_plot_max_cols = int(os.getenv("CNMF_GEP_AGG_MAX_COLUMNS", "30"))

    summary_rows: List[Dict[str, object]] = []
    aggregate_summary_rows: List[Dict[str, object]] = []

    for level in levels:
        if level not in LEVEL_CONFIG:
            print(f"[WARN] Unsupported level: {level}")
            continue

        sources = discover_sources(run_root, level, lineages=lineages or None, celltypes=celltypes or None)
        if not sources:
            print(f"[WARN] No completed sources discovered for level={level}")
            continue

        aggregate_payloads_by_k: Dict[int, List[Dict[str, object]]] = {}
        for record in sources:
            source_rows, payloads = process_source(
                record,
                group_columns=group_columns,
                min_row_sum=min_row_sum,
                source_plot_max_cols=source_plot_max_cols,
            )
            summary_rows.extend(source_rows)
            for k, payload in payloads.items():
                aggregate_payloads_by_k.setdefault(int(k), []).append(payload)

        for k, payloads in sorted(aggregate_payloads_by_k.items()):
            manifest_path = write_aggregate_outputs(
                run_root,
                level=level,
                k=int(k),
                payloads=payloads,
                group_columns=group_columns,
                min_row_sum=min_row_sum,
                aggregate_plot_max_cols=aggregate_plot_max_cols,
            )
            if manifest_path is not None:
                aggregate_summary_rows.append(
                    {
                        "level": level,
                        "k": int(k),
                        "n_sources": len(payloads),
                        "manifest_path": str(manifest_path),
                    }
                )

    derived_root = run_root / "cnmf_gep_celltype_derived"
    derived_root.mkdir(parents=True, exist_ok=True)
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(derived_root / "source_manifest.tsv", sep="\t", index=False)
    if aggregate_summary_rows:
        pd.DataFrame(aggregate_summary_rows).to_csv(derived_root / "aggregate_manifest.tsv", sep="\t", index=False)

    print(json.dumps({
        "run_root": str(run_root),
        "levels": levels,
        "n_source_rows": len(summary_rows),
        "n_aggregate_rows": len(aggregate_summary_rows),
        "source_manifest": str(derived_root / "source_manifest.tsv"),
        "aggregate_manifest": str(derived_root / "aggregate_manifest.tsv"),
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
