#!/usr/bin/env python3
"""Multi-route annotation consistency workflow around scExtract, CellTypist, ScType, and mLLMCelltype.

Notes
-----
The current upstream `yxwucq/scExtract` repository on GitHub main provides the
built-in literature-guided `auto_extract.auto_extract()` function, but does *not*
currently ship `add_celltypist_annotation` or `add_sctype_annotation` helpers on
that branch. Therefore this orchestrator integrates upstream `scExtract` at the
function level where available, and adds companion local function-based branches
for CellTypist, ScType-style annotation, plus a pluggable `mLLMCelltype` adapter
branch that can consume a precomputed h5ad or a configured Python callable.
"""

from __future__ import annotations

import argparse
import gc
import importlib
import inspect
import json
import math
import os
import re
import sys
import traceback
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover - handled at runtime
    matplotlib = None
    plt = None

try:
    import yaml
except ImportError:  # pragma: no cover - handled at runtime
    yaml = None

from scextract_annotation_bridge_20260524 import apply_bridge

try:
    ad.settings.allow_write_nullable_strings = True
except Exception:
    pass


WORKSPACE_ROOT = Path("/home/h2048")
DEFAULT_ENV_PYTHON = WORKSPACE_ROOT / "miniconda3" / "envs" / "scextract_env" / "bin" / "python"

CELLTYPIST_LINEAGE_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("stromal_endothelial", (r"stromal[_-]?endothelial", r"(?:^|[/_-])endothelial(?:[/_-]|$)")),
    ("stromal_fibroblast", (r"stromal[_-]?fibro(?:blast)?", r"(?:^|[/_-])fibro(?:blast)?(?:[/_-]|$)")),
    ("stromal_smc", (r"stromal[_-]?smc", r"smooth[_-]?muscle", r"pericyte", r"(?:^|[/_-])smc(?:[/_-]|$)")),
    ("epithelial", (r"(?:^|[/_-])epithelial(?:[/_-]|$)",)),
    ("myeloid", (r"(?:^|[/_-])myeloid(?:[/_-]|$)", r"monocyte", r"macrophage", r"dendritic")),
    ("tnk", (r"(?:^|[/_-])tnk(?:[/_-]|$)", r"t[_-]?nk", r"(?:^|[/_-])nk(?:[/_-]|$)", r"tcell", r"t_cell")),
    ("bcell", (r"(?:^|[/_-])bcell(?:[/_-]|$)", r"b[_-]?cell", r"b[_-]?plasma", r"(?:^|[/_-])plasma(?:[/_-]|$)")),
]


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def stringify(value: Any) -> str:
    if value is None:
        return ""
    return str(value)


def close_backed_adata(adata: ad.AnnData | None) -> None:
    if adata is None:
        return
    file_handle = getattr(adata, "file", None)
    if file_handle is None:
        return
    try:
        file_handle.close()
    except Exception:
        pass


def should_materialize_preflight_copy(cfg: dict[str, Any]) -> bool:
    return bool(deep_get(cfg, "run", "materialize_preflight_copy", default=False))


def should_write_branch_h5ad(cfg: dict[str, Any]) -> bool:
    return bool(deep_get(cfg, "run", "write_branch_h5ad", default=False))


def write_annotation_table(df: pd.DataFrame, path: Path) -> Path:
    out = df.copy()
    out.index = out.index.astype(str)
    out.to_csv(path, sep="\t", compression="gzip", index_label="obs_name")
    return path


def read_annotation_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def write_tsv(df: pd.DataFrame, path: Path) -> Path:
    df.to_csv(path, sep="\t", index=False)
    return path


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if yaml is not None:
        cfg = yaml.safe_load(text)
        return cfg or {}
    try:
        return json.loads(text)
    except Exception as exc:  # pragma: no cover - runtime guard
        raise RuntimeError("PyYAML is not installed and config is not valid JSON.") from exc


def deep_get(obj: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = obj
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def now_stamp() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def resolve_run_dir(cfg: dict[str, Any], cli_output_dir: str | None) -> Path:
    if cli_output_dir:
        return ensure_dir(Path(cli_output_dir))
    output_root = Path(deep_get(cfg, "run", "output_root", default=str(WORKSPACE_ROOT / "data" / "py" / datetime.now().strftime("%Y%m%d"))))
    prefix = deep_get(cfg, "run", "run_name_prefix", default="scextract_three_way")
    return ensure_dir(output_root / f"{prefix}_{now_stamp()}")


def resolve_scextract_python(cfg: dict[str, Any]) -> str:
    configured = deep_get(cfg, "scextract", "python_executable", default="auto")
    if configured and str(configured).lower() != "auto":
        return str(configured)
    if DEFAULT_ENV_PYTHON.exists():
        return str(DEFAULT_ENV_PYTHON)
    return sys.executable


def infer_lineage_name_from_source_path(source_h5ad_path: Path | str | None) -> str | None:
    if source_h5ad_path is None:
        return None
    path = Path(source_h5ad_path)
    searchable = "/".join(part.lower() for part in path.parts[-6:])
    for lineage_name, patterns in CELLTYPIST_LINEAGE_PATTERNS:
        if any(re.search(pattern, searchable) for pattern in patterns):
            return lineage_name
    return None


def resolve_celltypist_model_path(cfg: dict[str, Any], source_h5ad_path: Path | str | None) -> tuple[Path, dict[str, Any]]:
    auto_select = bool(deep_get(cfg, "celltypist", "auto_select_model", default=False))
    configured_lineage = stringify(
        deep_get(
            cfg,
            "celltypist",
            "lineage_name",
            default=deep_get(cfg, "run", "lineage_name", default=""),
        )
    ).strip()
    configured_lineage = configured_lineage if configured_lineage.lower() not in {"", "auto"} else ""

    configured_model_path = stringify(deep_get(cfg, "celltypist", "model_path", default="")).strip()
    fallback_model_path = stringify(
        deep_get(cfg, "celltypist", "default_model_path", default=configured_model_path)
    ).strip()
    lineage_model_map = deep_get(cfg, "celltypist", "lineage_model_map", default={}) or {}
    if not isinstance(lineage_model_map, dict):
        raise TypeError("celltypist.lineage_model_map must be a mapping of lineage -> model path")

    resolved_lineage = configured_lineage or (infer_lineage_name_from_source_path(source_h5ad_path) if auto_select else None)
    resolved_model_raw = ""
    selection_source = "configured_model_path"

    if auto_select and resolved_lineage and resolved_lineage in lineage_model_map:
        resolved_model_raw = stringify(lineage_model_map[resolved_lineage]).strip()
        selection_source = "configured_lineage_name" if configured_lineage else "inferred_from_input_path"
        selection_source = f"{selection_source}_mapped"

    if not resolved_model_raw:
        resolved_model_raw = configured_model_path or fallback_model_path
        if auto_select:
            if configured_lineage:
                selection_source = "configured_lineage_name_fallback"
            elif resolved_lineage:
                selection_source = "inferred_from_input_path_fallback"
            else:
                selection_source = "default_model_path_fallback"

    if not resolved_model_raw:
        raise ValueError("CellTypist model path is empty after resolution")

    return Path(resolved_model_raw), {
        "auto_select_model": auto_select,
        "requested_lineage_name": configured_lineage or None,
        "resolved_lineage_name": resolved_lineage,
        "selection_source": selection_source,
        "source_h5ad_path": str(source_h5ad_path) if source_h5ad_path is not None else None,
        "default_model_path": fallback_model_path or None,
        "configured_model_path": configured_model_path or None,
    }


def load_workspace_dotenv(dotenv_path: Path = WORKSPACE_ROOT / ".env") -> bool:
    if not dotenv_path.exists():
        return False
    try:
        from dotenv import load_dotenv
    except Exception:
        return False
    load_dotenv(dotenv_path=dotenv_path, override=False)
    return True


def normalize_obs_strings(adata: ad.AnnData) -> None:
    for col in adata.obs.columns:
        dtype_name = getattr(adata.obs[col].dtype, "name", str(adata.obs[col].dtype))
        if dtype_name in {"object", "category", "string"}:
            adata.obs[col] = adata.obs[col].astype("string")


def is_count_like_matrix(matrix, max_sample: int = 2000) -> bool:
    if sp.issparse(matrix):
        data = np.asarray(matrix.data)
    else:
        data = np.asarray(matrix).ravel()
    if data.size == 0:
        return False
    if data.size > max_sample:
        rng = np.random.default_rng(42)
        idx = rng.choice(data.size, size=max_sample, replace=False)
        data = data[idx]
    data = data[np.isfinite(data)]
    if data.size == 0:
        return False
    return bool(np.all(data >= 0) and np.allclose(data, np.round(data), atol=1e-6))


def ensure_counts_layer(adata: ad.AnnData, counts_layer: str) -> dict[str, Any]:
    if counts_layer in adata.layers:
        return {"counts_layer": counts_layer, "created": False, "source": counts_layer}
    if is_count_like_matrix(adata.X):
        adata.layers[counts_layer] = adata.X.copy()
        return {"counts_layer": counts_layer, "created": True, "source": "X"}
    return {"counts_layer": counts_layer, "created": False, "source": None, "warning": "Missing counts layer and X is not count-like"}


def normalize_gene_names(adata: ad.AnnData) -> None:
    if adata.var_names.duplicated().sum() > 0:
        adata.var_names_make_unique()
    symbol_cols = ["gene_symbol", "gene_symbols", "symbol", "feature_name", "gene_name"]
    for col in symbol_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str).str.strip()
            non_empty = (vals != "").mean() if len(vals) else 0.0
            if non_empty > 0.5:
                adata.var["symbol_base"] = vals.values
                return
    adata.var["symbol_base"] = adata.var_names.astype(str).values


def preserve_full_raw_if_missing(adata: ad.AnnData, counts_layer: str) -> None:
    if adata.raw is not None and adata.raw.n_vars >= adata.n_vars:
        if "symbol_base" in adata.var.columns and "symbol_base" not in adata.raw.var.columns:
            common = adata.raw.var_names.intersection(adata.var_names)
            if len(common) > 0:
                adata.raw.var.loc[common, "symbol_base"] = adata.var.loc[common, "symbol_base"].values
        return
    if counts_layer not in adata.layers:
        raise KeyError(f"Required counts layer not found: {counts_layer}")
    adata.raw = ad.AnnData(X=adata.layers[counts_layer], obs=adata.obs.copy(), var=adata.var.copy())


def build_celltypist_input(adata: ad.AnnData, model_features: pd.Index | None, counts_layer: str) -> ad.AnnData:
    import scanpy as sc

    if adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var.copy()
    else:
        if counts_layer not in adata.layers:
            raise KeyError(f"Required counts layer not found: {counts_layer}")
        X = adata.layers[counts_layer]
        var = adata.var.copy()

    if "symbol_base" in var.columns:
        base_symbols = var["symbol_base"].astype(str).values
    elif "symbol_base" in adata.var.columns:
        mapping = adata.var["symbol_base"].astype(str).to_dict()
        base_symbols = np.array([mapping.get(g, str(g)) for g in var.index], dtype=object)
    else:
        base_symbols = var.index.astype(str).values

    base_index = pd.Index(base_symbols)
    if model_features is not None:
        keep = base_index.isin(model_features) & (~base_index.duplicated(keep="first"))
    else:
        keep = ~base_index.duplicated(keep="first")

    X_ct = X[:, keep].copy()
    var_ct = pd.DataFrame(index=base_index[keep])
    obs_ct = pd.DataFrame(index=adata.obs_names)
    ad_ct = ad.AnnData(X=X_ct, obs=obs_ct, var=var_ct)
    sc.pp.normalize_total(ad_ct, target_sum=1e4)
    sc.pp.log1p(ad_ct)
    return ad_ct


def extract_celltypist_confidence(pred) -> np.ndarray:
    if hasattr(pred, "probability_matrix") and pred.probability_matrix is not None:
        try:
            return np.asarray(pred.probability_matrix.max(axis=1)).ravel()
        except Exception:
            pass
    pred_df = getattr(pred, "predicted_labels", None)
    if isinstance(pred_df, pd.DataFrame):
        for col in ("conf_score", "confidence", "confidence_score", "prob"):
            if col in pred_df.columns:
                return pred_df[col].to_numpy()
    return np.ones(shape=(len(pred.predicted_labels),), dtype=float)


def normalize_annotation_values(series: pd.Series, unknown_label: str) -> pd.Series:
    values = series.astype("string").fillna(unknown_label).astype(str).str.strip()
    return values.replace({"": unknown_label, "nan": unknown_label, "None": unknown_label, "<NA>": unknown_label})


def pretty_annotation_label(column_name: str) -> str:
    mapping = {
        "cell_type_scextract": "scExtract",
        "cell_type_celltypist_scextract": "CellTypist",
        "cell_type_sctype_scextract": "ScType",
        "cell_type_mllmcelltype_scextract": "mLLMCelltype",
    }
    if column_name in mapping:
        return mapping[column_name]
    return column_name.replace("cell_type_", "").replace("_scextract", "")


def shorten_text(value: str, limit: int = 28) -> str:
    text = stringify(value).strip()
    if len(text) <= limit:
        return text
    return f"{text[: limit - 1]}…"


def extract_obsm_matrix(h5ad_path: Path, key: str) -> np.ndarray | None:
    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        try:
            keys = list(backed.obsm.keys())
        except Exception:
            keys = []
        if key not in keys:
            return None
        matrix = np.asarray(backed.obsm[key])
        if matrix.ndim != 2 or matrix.shape[1] < 2:
            return None
        return matrix[:, :2]
    finally:
        close_backed_adata(backed)


def save_matplotlib_figure(fig, figure_dir: Path, stem: str, dpi: int, write_pdf: bool) -> dict[str, str | None]:
    png_path = figure_dir / f"{stem}.png"
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    pdf_path = figure_dir / f"{stem}.pdf" if write_pdf else None
    if pdf_path is not None:
        fig.savefig(pdf_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return {
        "png_path": str(png_path),
        "pdf_path": str(pdf_path) if pdf_path is not None else None,
    }


def plot_label_count_panels(
    count_tables: dict[str, pd.DataFrame],
    available_columns: list[str],
    figure_dir: Path,
    top_n: int,
    dpi: int,
    write_pdf: bool,
) -> dict[str, Any]:
    if plt is None or not available_columns:
        return {
            "figure_type": "label_count_panels",
            "status": "skipped",
            "note": "matplotlib unavailable or no available annotation columns",
            "png_path": None,
            "pdf_path": None,
        }

    n_panels = len(available_columns)
    ncols = min(2, n_panels)
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax in axes[n_panels:]:
        ax.axis("off")

    for ax, column in zip(axes, available_columns):
        counts = count_tables.get(column, pd.DataFrame(columns=["label", "count"]))
        counts = counts.head(top_n).copy()
        if counts.empty:
            ax.text(0.5, 0.5, "No labels", ha="center", va="center")
            ax.axis("off")
            continue
        counts["display_label"] = counts["label"].map(lambda x: shorten_text(x, limit=32))
        counts = counts.sort_values("count", ascending=True)
        ax.barh(counts["display_label"], counts["count"], color="#4C78A8")
        ax.set_title(f"{pretty_annotation_label(column)} top labels")
        ax.set_xlabel("Cells")
        ax.set_ylabel("Label")
        for y, value in enumerate(counts["count"]):
            ax.text(float(value), y, f" {int(value)}", va="center", fontsize=8)

    fig.suptitle("Top label counts by annotation method", fontsize=14)
    fig.tight_layout()
    saved = save_matplotlib_figure(fig, figure_dir, "label_count_panels", dpi=dpi, write_pdf=write_pdf)
    return {
        "figure_type": "label_count_panels",
        "status": "ok",
        "note": f"top_n={top_n}",
        **saved,
    }


def plot_pairwise_agreement_heatmap(
    pairwise: pd.DataFrame,
    columns: list[str],
    available_columns: list[str],
    figure_dir: Path,
    dpi: int,
    write_pdf: bool,
) -> dict[str, Any]:
    if plt is None or not columns:
        return {
            "figure_type": "pairwise_agreement_heatmap",
            "status": "skipped",
            "note": "matplotlib unavailable or no annotation columns configured",
            "png_path": None,
            "pdf_path": None,
        }

    matrix = pd.DataFrame(np.nan, index=columns, columns=columns, dtype=float)
    for column in columns:
        if column in available_columns:
            matrix.loc[column, column] = 1.0
    for row in pairwise.to_dict(orient="records"):
        left = row["left"]
        right = row["right"]
        value = row["exact_match_rate"]
        matrix.loc[left, right] = value
        matrix.loc[right, left] = value

    display = [pretty_annotation_label(col) for col in columns]
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    fig, ax = plt.subplots(figsize=(1.8 * len(columns) + 2.5, 1.5 * len(columns) + 2.0))
    im = ax.imshow(masked, cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(columns)))
    ax.set_yticks(range(len(columns)))
    ax.set_xticklabels(display, rotation=30, ha="right")
    ax.set_yticklabels(display)
    ax.set_title("Pairwise annotation agreement")
    for i in range(len(columns)):
        for j in range(len(columns)):
            value = matrix.iat[i, j]
            if pd.notna(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Exact match rate")
    fig.tight_layout()
    saved = save_matplotlib_figure(fig, figure_dir, "pairwise_agreement_heatmap", dpi=dpi, write_pdf=write_pdf)
    return {
        "figure_type": "pairwise_agreement_heatmap",
        "status": "ok",
        "note": f"n_methods={len(columns)}",
        **saved,
    }


def plot_truth_benchmark_heatmap(
    truth_df: pd.DataFrame,
    figure_dir: Path,
    dpi: int,
    write_pdf: bool,
) -> dict[str, Any]:
    if plt is None:
        return {
            "figure_type": "truth_benchmark_heatmap",
            "status": "skipped",
            "note": "matplotlib unavailable",
            "png_path": None,
            "pdf_path": None,
        }
    if truth_df.empty:
        return {
            "figure_type": "truth_benchmark_heatmap",
            "status": "skipped",
            "note": "truth benchmark unavailable",
            "png_path": None,
            "pdf_path": None,
        }

    metric_cols = [col for col in ("exact_match_rate", "ari", "nmi") if col in truth_df.columns]
    if not metric_cols:
        return {
            "figure_type": "truth_benchmark_heatmap",
            "status": "skipped",
            "note": "no benchmark metrics present",
            "png_path": None,
            "pdf_path": None,
        }

    matrix = truth_df.set_index("prediction_column")[metric_cols].copy()
    matrix.index = [pretty_annotation_label(x) for x in matrix.index]
    masked = np.ma.masked_invalid(matrix.to_numpy(dtype=float))
    fig, ax = plt.subplots(figsize=(1.8 * len(metric_cols) + 2.5, 1.1 * len(matrix.index) + 2.0))
    im = ax.imshow(masked, cmap="magma", vmin=0.0, vmax=1.0)
    ax.set_xticks(range(len(metric_cols)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(metric_cols)
    ax.set_yticklabels(matrix.index)
    ax.set_title("Benchmark vs truth labels")
    for i in range(len(matrix.index)):
        for j in range(len(metric_cols)):
            value = matrix.iat[i, j]
            if pd.notna(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Score")
    fig.tight_layout()
    saved = save_matplotlib_figure(fig, figure_dir, "truth_benchmark_heatmap", dpi=dpi, write_pdf=write_pdf)
    return {
        "figure_type": "truth_benchmark_heatmap",
        "status": "ok",
        "note": f"metrics={','.join(metric_cols)}",
        **saved,
    }


def plot_annotation_umap_panels(
    obs: pd.DataFrame,
    available_columns: list[str],
    merged_h5ad: Path,
    figure_dir: Path,
    cfg: dict[str, Any],
    dpi: int,
    write_pdf: bool,
) -> dict[str, Any]:
    if plt is None:
        return {
            "figure_type": "annotation_umap_panels",
            "status": "skipped",
            "note": "matplotlib unavailable",
            "png_path": None,
            "pdf_path": None,
        }
    if not available_columns:
        return {
            "figure_type": "annotation_umap_panels",
            "status": "skipped",
            "note": "no available annotation columns",
            "png_path": None,
            "pdf_path": None,
        }

    umap_key = stringify(deep_get(cfg, "summary", "umap_key", default="X_umap")).strip() or "X_umap"
    if not bool(deep_get(cfg, "summary", "generate_umap_if_present", default=True)):
        return {
            "figure_type": "annotation_umap_panels",
            "status": "skipped",
            "note": "disabled by config",
            "png_path": None,
            "pdf_path": None,
        }

    coords = extract_obsm_matrix(merged_h5ad, umap_key)
    if coords is None:
        return {
            "figure_type": "annotation_umap_panels",
            "status": "skipped",
            "note": f"{umap_key} not present in merged h5ad",
            "png_path": None,
            "pdf_path": None,
        }

    max_legend = int(deep_get(cfg, "summary", "max_legend_labels", default=12))
    point_size = float(deep_get(cfg, "summary", "umap_point_size", default=4.0))
    alpha = float(deep_get(cfg, "summary", "umap_alpha", default=0.8))
    n_panels = len(available_columns)
    ncols = min(2, n_panels)
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.5 * ncols, 5.8 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax in axes[n_panels:]:
        ax.axis("off")

    for ax, column in zip(axes, available_columns):
        series = obs[column].astype("string").fillna("Unknown").astype(str)
        top_labels = list(series.value_counts().head(max_legend).index)
        collapsed = series.where(series.isin(top_labels), other="Other")
        categories = list(pd.Index(top_labels + (["Other"] if "Other" in collapsed.values else [])))
        if "Unknown" in collapsed.values and "Unknown" not in categories:
            categories.append("Unknown")
        categories = [cat for cat in categories if cat in set(collapsed.values)]
        cmap = plt.get_cmap("tab20", max(len(categories), 1))
        for idx, category in enumerate(categories):
            mask = collapsed.values == category
            if not np.any(mask):
                continue
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=point_size,
                alpha=alpha,
                c=[cmap(idx)],
                label=shorten_text(category, limit=28),
                linewidths=0,
                rasterized=False,
            )
        ax.set_title(pretty_annotation_label(column))
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.legend(loc="best", fontsize=7, markerscale=2, frameon=False)

    fig.suptitle("Annotation distributions on existing UMAP", fontsize=14)
    fig.tight_layout()
    saved = save_matplotlib_figure(fig, figure_dir, "annotation_umap_panels", dpi=dpi, write_pdf=write_pdf)
    return {
        "figure_type": "annotation_umap_panels",
        "status": "ok",
        "note": f"umap_key={umap_key}",
        **saved,
    }


def generate_summary_visualizations(
    obs: pd.DataFrame,
    count_tables: dict[str, pd.DataFrame],
    columns: list[str],
    available_columns: list[str],
    pairwise: pd.DataFrame,
    truth_df: pd.DataFrame,
    merged_h5ad: Path,
    summary_dir: Path,
    cfg: dict[str, Any],
) -> tuple[list[dict[str, Any]], Path]:
    figure_dir = ensure_dir(summary_dir / "figures")
    top_n = int(deep_get(cfg, "summary", "top_n_labels_plot", default=12))
    dpi = int(deep_get(cfg, "summary", "figure_dpi", default=150))
    write_pdf = bool(deep_get(cfg, "summary", "write_pdf", default=True))

    manifest_rows = [
        plot_label_count_panels(count_tables, available_columns, figure_dir, top_n=top_n, dpi=dpi, write_pdf=write_pdf),
        plot_pairwise_agreement_heatmap(pairwise, columns, available_columns, figure_dir, dpi=dpi, write_pdf=write_pdf),
        plot_truth_benchmark_heatmap(truth_df, figure_dir, dpi=dpi, write_pdf=write_pdf),
        plot_annotation_umap_panels(obs, available_columns, merged_h5ad, figure_dir, cfg, dpi=dpi, write_pdf=write_pdf),
    ]

    manifest_path = write_tsv(pd.DataFrame(manifest_rows), summary_dir / "figure_manifest.tsv")
    return manifest_rows, manifest_path


def series_numeric_like(series: pd.Series) -> bool:
    sample = series.astype(str).replace({"nan": None, "<NA>": None}).dropna().head(20)
    if sample.empty:
        return False
    return bool(sample.map(lambda x: re.fullmatch(r"-?\d+(\.\d+)?", x) is not None).all())


def infer_scextract_label_column(adata: ad.AnnData) -> str:
    return infer_scextract_label_column_from_obs(adata.obs)


def infer_scextract_label_column_from_obs(obs: pd.DataFrame) -> str:
    candidates = []
    for col in ("leiden", "louvain"):
        if col in obs.columns and not series_numeric_like(obs[col]):
            candidates.append(col)
    if candidates:
        return candidates[0]
    for col in obs.columns:
        if col.endswith("_rough") and not series_numeric_like(obs[col]):
            return col
    raise RuntimeError("Unable to infer scExtract annotation label column from output AnnData")


def extract_obs_columns_from_h5ad(h5ad_path: Path, columns: list[str] | None = None) -> pd.DataFrame:
    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        obs = backed.obs.copy()
    finally:
        close_backed_adata(backed)
    if columns is None:
        return obs
    keep = [col for col in columns if col in obs.columns]
    if not keep:
        return pd.DataFrame(index=obs.index.copy())
    return obs.loc[:, keep].copy()


def load_branch_annotation_obs(
    artifact: dict[str, str | None] | None,
    fallback_columns: list[str] | None = None,
) -> pd.DataFrame | None:
    if not artifact:
        return None

    annotations_tsv = artifact.get("annotations_tsv")
    if annotations_tsv:
        table_path = Path(annotations_tsv)
        if table_path.exists():
            return read_annotation_table(table_path)

    h5ad_path = artifact.get("h5ad")
    if h5ad_path:
        resolved = Path(h5ad_path)
        if resolved.exists():
            return extract_obs_columns_from_h5ad(resolved, columns=fallback_columns)

    return None


def save_status(branch_dir: Path, branch_name: str, status: str, message: str = "", **extra) -> dict[str, Any]:
    payload = {"branch": branch_name, "status": status, "message": message, **extra}
    write_json(payload, branch_dir / f"{branch_name}_status.json")
    return payload


def ensure_batch_key(adata: ad.AnnData, batch_key: str, auto_create: bool) -> dict[str, Any]:
    if batch_key in adata.obs.columns:
        return {"batch_key": batch_key, "created": False}
    if auto_create:
        adata.obs[batch_key] = pd.Categorical(["batch_0"] * adata.n_obs)
        return {"batch_key": batch_key, "created": True}
    return {"batch_key": batch_key, "created": False, "warning": f"Missing batch key: {batch_key}"}


def ensure_cluster_key(adata: ad.AnnData, cluster_key: str, auto_compute: bool) -> dict[str, Any]:
    if cluster_key in adata.obs.columns:
        return {"cluster_key": cluster_key, "created": False}
    if not auto_compute:
        raise KeyError(f"Missing required cluster key: {cluster_key}")

    import scanpy as sc

    work_key = "X_pca"
    if work_key not in adata.obsm:
        n_comps = max(2, min(50, adata.n_obs - 1, adata.n_vars - 1))
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=42)
    n_neighbors = max(2, min(15, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=42)
    if cluster_key == "louvain":
        sc.tl.louvain(adata, key_added=cluster_key, random_state=42)
    else:
        sc.tl.leiden(adata, key_added=cluster_key, random_state=42)
    return {"cluster_key": cluster_key, "created": True}


def resolve_cluster_key(adata: ad.AnnData, cfg: dict[str, Any]) -> tuple[str, str]:
    requested = stringify(deep_get(cfg, "data_contract", "cluster_key", default="leiden")).strip() or "leiden"
    if requested in adata.obs.columns:
        return requested, "requested"

    candidates = deep_get(
        cfg,
        "data_contract",
        "cluster_key_candidates",
        default=["leiden", "seurat_clusters", "louvain", "choir_l3", "choir_cluster", "cluster", "cell_type_L3", "cell_type_L2"],
    )
    if not isinstance(candidates, list):
        candidates = [requested]
    for candidate in candidates:
        candidate_name = stringify(candidate).strip()
        if candidate_name and candidate_name in adata.obs.columns:
            return candidate_name, "fallback_existing"
    return requested, "to_be_computed"


def prepare_preflight_adata(input_h5ad: Path, cfg: dict[str, Any], runtime_dir: Path) -> tuple[Path, dict[str, Any]]:
    materialize_reason = "config_requested"
    if not should_materialize_preflight_copy(cfg):
        counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
        batch_key = deep_get(cfg, "data_contract", "batch_key", default="Batch")
        auto_create_batch = bool(deep_get(cfg, "data_contract", "auto_create_batch_key", default=True))
        backed = ad.read_h5ad(input_h5ad, backed="r")
        try:
            cluster_key, cluster_key_source = resolve_cluster_key(backed, cfg)
            counts_source = counts_layer if counts_layer in backed.layers else None
            batch_present = batch_key in backed.obs.columns
            if batch_present or not auto_create_batch:
                return input_h5ad, {
                    "counts": {
                        "counts_layer": counts_layer,
                        "created": False,
                        "source": counts_source,
                        "warning": None if counts_source else f"Missing counts layer: {counts_layer}",
                    },
                    "batch": {
                        "batch_key": batch_key,
                        "created": False,
                        "warning": None if batch_present else f"Missing batch key: {batch_key}",
                    },
                    "raw": {
                        "raw_ready": backed.raw is not None,
                        "raw_n_vars": int(backed.raw.n_vars) if backed.raw is not None else 0,
                    },
                    "cluster": {
                        "cluster_key": cluster_key,
                        "created": False,
                        "prepared": cluster_key in backed.obs.columns,
                        "cluster_key_source": cluster_key_source,
                    },
                    "n_obs": int(backed.n_obs),
                    "n_vars": int(backed.n_vars),
                    "preflight_path": str(input_h5ad),
                    "materialized_copy": False,
                }
            materialize_reason = f"missing_batch_key:{batch_key}"
        finally:
            close_backed_adata(backed)

    adata = ad.read_h5ad(input_h5ad)
    normalize_obs_strings(adata)
    normalize_gene_names(adata)

    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    batch_key = deep_get(cfg, "data_contract", "batch_key", default="Batch")
    cluster_key = deep_get(cfg, "data_contract", "cluster_key", default="leiden")

    counts_info = ensure_counts_layer(adata, counts_layer)
    batch_info = ensure_batch_key(adata, batch_key, bool(deep_get(cfg, "data_contract", "auto_create_batch_key", default=True)))

    preserve_info: dict[str, Any]
    try:
        preserve_full_raw_if_missing(adata, counts_layer)
        preserve_info = {"raw_ready": True, "raw_n_vars": int(adata.raw.n_vars) if adata.raw is not None else 0}
    except Exception as exc:
        preserve_info = {"raw_ready": False, "warning": str(exc)}

    cluster_info = {"cluster_key": cluster_key, "created": False, "prepared": cluster_key in adata.obs.columns}
    preflight_path = runtime_dir / "preflight_input.h5ad"
    adata.write_h5ad(preflight_path, compression="gzip")
    return preflight_path, {
        "counts": counts_info,
        "batch": batch_info,
        "raw": preserve_info,
        "cluster": cluster_info,
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "preflight_path": str(preflight_path),
        "materialized_copy": True,
        "materialize_reason": materialize_reason,
    }


def ensure_scextract_repo_ready(repo_dir: Path) -> None:
    main_py = repo_dir / "main.py"
    config_py = repo_dir / "auto_extract" / "config.py"
    if not main_py.exists():
        raise FileNotFoundError(f"scExtract main.py not found: {main_py}")
    if not config_py.exists():
        raise FileNotFoundError(f"scExtract auto_extract/config.py not found: {config_py}")


def load_scextract_auto_extract(repo_dir: Path):
    repo_dir_str = str(repo_dir)
    if repo_dir_str not in sys.path:
        sys.path.insert(0, repo_dir_str)
    module = importlib.import_module("auto_extract.auto_extract")
    auto_extract_func = getattr(module, "auto_extract", None)
    if auto_extract_func is None:
        raise AttributeError("Could not resolve auto_extract.auto_extract from scExtract repository")
    return auto_extract_func


def load_configured_callable(module_path: str, callable_name: str, repo_dir: Path | None = None):
    if repo_dir is not None:
        repo_dir = Path(repo_dir)
        if not repo_dir.exists():
            raise FileNotFoundError(f"Configured repo_dir does not exist: {repo_dir}")
        repo_dir_str = str(repo_dir)
        if repo_dir_str not in sys.path:
            sys.path.insert(0, repo_dir_str)
    module = importlib.import_module(module_path)
    target = module
    for attr in callable_name.split("."):
        target = getattr(target, attr)
    return target


def call_with_supported_kwargs(func, candidate_kwargs: dict[str, Any]):
    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        return func(**candidate_kwargs)
    accepts_var_kw = any(param.kind == inspect.Parameter.VAR_KEYWORD for param in sig.parameters.values())
    supported_kwargs = candidate_kwargs if accepts_var_kw else {k: v for k, v in candidate_kwargs.items() if k in sig.parameters}
    return func(**supported_kwargs)


def resolve_callable_result(result: Any, output_path: Path) -> dict[str, Any]:
    if isinstance(result, ad.AnnData):
        return {"kind": "adata", "adata": result}
    if isinstance(result, (str, Path)):
        result_path = Path(result)
        if not result_path.exists():
            raise FileNotFoundError(f"Callable returned output path that does not exist: {result_path}")
        kind = "h5ad" if result_path.suffix == ".h5ad" else "table"
        return {"kind": kind, "path": result_path}
    if isinstance(result, dict):
        annotation_table = stringify(result.get("annotation_table_path", result.get("obs_table_path", ""))).strip()
        output_h5ad = stringify(result.get("output_h5ad", "")).strip()
        if annotation_table:
            table_path = Path(annotation_table)
            if not table_path.exists():
                raise FileNotFoundError(f"Callable returned annotation_table_path that does not exist: {table_path}")
            return {
                "kind": "table",
                "path": table_path,
                "label_column": result.get("label_column"),
                "confidence_column": result.get("confidence_column"),
                "score_column": result.get("score_column"),
            }
        if output_h5ad:
            h5ad_path = Path(output_h5ad)
            if not h5ad_path.exists():
                raise FileNotFoundError(f"Callable returned output_h5ad that does not exist: {h5ad_path}")
            return {
                "kind": "h5ad",
                "path": h5ad_path,
                "label_column": result.get("label_column"),
                "confidence_column": result.get("confidence_column"),
                "score_column": result.get("score_column"),
            }
    if result is None:
        if output_path.exists():
            return {"kind": "h5ad", "path": output_path}
        raise FileNotFoundError(f"Callable returned None and expected output was not created: {output_path}")
    raise TypeError(f"Unsupported callable return type: {type(result)!r}")


def resolve_bridge_confidence_settings(cfg: dict[str, Any], selected_label_source: str) -> tuple[str | None, float | None]:
    if selected_label_source == "celltypist_annotation":
        return "celltypist_confidence_scextract", deep_get(cfg, "celltypist", "min_confidence", default=None)
    if selected_label_source == "mllmcelltype_annotation":
        return "mllmcelltype_confidence_scextract", deep_get(cfg, "mllmcelltype", "min_confidence", default=None)
    return None, None


def run_scextract_branch(preflight_h5ad: Path, pdf_path: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], dict[str, str | None] | None]:
    branch_name = "scExtract_builtin"
    output_name = deep_get(cfg, "scextract", "output_name", default="scextract_processed.h5ad")
    repo_dir = Path(deep_get(cfg, "scextract", "repo_dir", default=str(WORKSPACE_ROOT / "tools" / "scExtract")))
    ensure_scextract_repo_ready(repo_dir)

    if not os.getenv(deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY")):
        return save_status(branch_dir, branch_name, "error", "Missing DEEPSEEK_API_KEY in environment"), None

    deepseek_env_key = stringify(deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY")).strip() or "DEEPSEEK_API_KEY"
    resolved_api_key = os.getenv(deepseek_env_key)
    if resolved_api_key:
        os.environ["SCEXTRACT_API_KEY"] = resolved_api_key
    os.environ["SCEXTRACT_API_TYPE"] = stringify(deep_get(cfg, "deepseek", "type", default="openai")).strip() or "openai"
    os.environ["SCEXTRACT_API_BASE_URL"] = stringify(deep_get(cfg, "deepseek", "base_url", default="https://api.deepseek.com/v1")).strip() or "https://api.deepseek.com/v1"
    os.environ["SCEXTRACT_MODEL"] = stringify(deep_get(cfg, "deepseek", "model", default="deepseek-chat")).strip() or "deepseek-chat"
    os.environ["SCEXTRACT_TOOL_MODEL"] = stringify(deep_get(cfg, "deepseek", "tool_model", default=deep_get(cfg, "deepseek", "model", default="deepseek-chat"))).strip() or os.environ["SCEXTRACT_MODEL"]

    output_path = branch_dir / output_name
    function_log = branch_dir / "scExtract_builtin_function_call.json"
    write_json(
        {
            "function": "auto_extract.auto_extract",
            "repo_dir": str(repo_dir),
            "adata_path": str(preflight_h5ad),
            "pdf_path": str(pdf_path),
            "output_dir": str(branch_dir),
            "output_name": output_name,
        },
        function_log,
    )

    try:
        auto_extract_func = load_scextract_auto_extract(repo_dir)
        auto_extract_func(
            adata_path=str(preflight_h5ad),
            pdf_path=str(pdf_path),
            output_dir=str(branch_dir),
            output_name=output_name,
        )
    except Exception as exc:
        tb_path = branch_dir / "scExtract_builtin_traceback.txt"
        tb_path.write_text(traceback.format_exc(), encoding="utf-8")
        return save_status(
            branch_dir,
            branch_name,
            "error",
            f"scExtract function call failed: {exc}",
            traceback_file=str(tb_path),
            function="auto_extract.auto_extract",
            function_log=str(function_log),
        ), None

    log_path = branch_dir / "auto_extract.log"
    if not output_path.exists():
        return save_status(branch_dir, branch_name, "error", "Expected scExtract output h5ad not found", log_path=str(log_path), function_log=str(function_log)), None

    backed_out = ad.read_h5ad(output_path, backed="r")
    try:
        label_col = infer_scextract_label_column(backed_out)
        obs_cols = [label_col, "Certainty", "Source", "Tissue"]
        annotations = backed_out.obs.loc[:, [col for col in obs_cols if col in backed_out.obs.columns]].copy()
        annotation_path = write_annotation_table(annotations, branch_dir / "scExtract_annotations.tsv.gz")
        n_obs = int(backed_out.n_obs)
        n_vars = int(backed_out.n_vars)
    finally:
        close_backed_adata(backed_out)

    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"Built-in scExtract function completed with label column {label_col}",
        output_h5ad=str(output_path),
        annotation_table=str(annotation_path),
        label_column=label_col,
        function="auto_extract.auto_extract",
        function_log=str(function_log),
        log_path=str(log_path),
        n_obs=n_obs,
        n_vars=n_vars,
    ), {"annotations_tsv": str(annotation_path), "h5ad": str(output_path)}


def run_celltypist_branch(preflight_h5ad: Path, cfg: dict[str, Any], branch_dir: Path, source_h5ad_path: Path | None = None) -> tuple[dict[str, Any], dict[str, str | None] | None]:
    branch_name = "celltypist"
    if not bool(deep_get(cfg, "celltypist", "enabled", default=True)):
        return save_status(branch_dir, branch_name, "skipped", "CellTypist branch disabled in config"), None

    try:
        import celltypist
        from celltypist import models
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"CellTypist import failed: {exc}"), None

    adata = ad.read_h5ad(preflight_h5ad)
    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    counts_info = ensure_counts_layer(adata, counts_layer)
    normalize_gene_names(adata)
    preserve_full_raw_if_missing(adata, counts_layer)

    model_path, model_selection = resolve_celltypist_model_path(cfg, source_h5ad_path)
    if not model_path.exists():
        return save_status(
            branch_dir,
            branch_name,
            "error",
            f"CellTypist model not found: {model_path}",
            **model_selection,
        ), None

    os.environ.setdefault("CELLTYPIST_FOLDER", str(model_path.parent))
    ct_model = models.Model.load(str(model_path))

    model_features = None
    for attr in ("features", "genes", "feature_names"):
        if hasattr(ct_model, attr):
            try:
                model_features = pd.Index([str(x) for x in getattr(ct_model, attr)])
                break
            except Exception:
                model_features = None

    ad_ct = build_celltypist_input(adata, model_features=model_features, counts_layer=counts_layer)
    majority_voting = bool(deep_get(cfg, "celltypist", "majority_voting", default=True))
    pred = celltypist.annotate(ad_ct, model=ct_model, majority_voting=majority_voting)
    pred_df = pred.predicted_labels.copy()
    pred_df.index = ad_ct.obs_names

    obs_out = pd.DataFrame(index=adata.obs_names)
    obs_out["celltypist_annotation"] = pred_df["predicted_labels"].astype(str).reindex(adata.obs_names).values
    if majority_voting and "majority_voting" in pred_df.columns:
        obs_out["celltypist_majority_voting"] = pred_df["majority_voting"].astype(str).reindex(adata.obs_names).values
    obs_out["celltypist_confidence_scextract"] = extract_celltypist_confidence(pred)

    annotation_path = write_annotation_table(obs_out, branch_dir / "celltypist_annotations.tsv.gz")
    output_path: Path | None = None
    if should_write_branch_h5ad(cfg):
        adata.obs["celltypist_annotation"] = obs_out["celltypist_annotation"].values
        if "celltypist_majority_voting" in obs_out.columns:
            adata.obs["celltypist_majority_voting"] = obs_out["celltypist_majority_voting"].values
        adata.obs["celltypist_confidence_scextract"] = pd.to_numeric(obs_out["celltypist_confidence_scextract"], errors="coerce").values
        adata.uns["celltypist_model_path"] = str(model_path)
        adata.uns["celltypist_model_selection"] = model_selection
        adata.uns["celltypist_counts_layer_source"] = counts_info.get("source")
        output_path = branch_dir / "celltypist_annotation.h5ad"
        adata.write_h5ad(output_path, compression="gzip")

    mean_confidence = float(np.nanmean(pd.to_numeric(obs_out["celltypist_confidence_scextract"], errors="coerce")))
    n_unique_labels = int(pd.Series(obs_out["celltypist_annotation"].astype(str)).nunique())
    del pred_df, pred, ad_ct, ct_model
    gc.collect()
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"CellTypist completed with model {model_path.name}",
        output_h5ad=str(output_path) if output_path is not None else None,
        annotation_table=str(annotation_path),
        model_path=str(model_path),
        model_name=model_path.name,
        n_unique_labels=n_unique_labels,
        mean_confidence=mean_confidence,
        **model_selection,
    ), {"annotations_tsv": str(annotation_path), "h5ad": str(output_path) if output_path is not None else None}


def parse_marker_string(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    tokens = re.split(r"[,;/\s]+", text)
    return [tok.strip().upper() for tok in tokens if tok and tok.strip()]


def build_gene_index(adata: ad.AnnData) -> tuple[dict[str, int], Any]:
    if adata.raw is not None:
        var = adata.raw.var
        matrix = adata.raw.X
    else:
        var = adata.var
        matrix = adata.X
    if "symbol_base" in var.columns:
        genes = var["symbol_base"].astype(str).values
    else:
        genes = var.index.astype(str).values
    mapping: dict[str, int] = {}
    for idx, gene in enumerate(genes):
        key = gene.strip().upper()
        if key and key not in mapping:
            mapping[key] = idx
    return mapping, matrix


def run_sctype_branch(preflight_h5ad: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], dict[str, str | None] | None]:
    branch_name = "sctype"
    if not bool(deep_get(cfg, "sctype", "enabled", default=True)):
        return save_status(branch_dir, branch_name, "skipped", "ScType branch disabled in config"), None

    try:
        import scanpy as sc
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"scanpy import failed: {exc}"), None

    adata = ad.read_h5ad(preflight_h5ad)
    counts_layer = deep_get(cfg, "data_contract", "counts_layer", default="counts")
    cluster_key, cluster_key_source = resolve_cluster_key(adata, cfg)
    db_path = Path(deep_get(cfg, "sctype", "db_path", default=""))
    tissue = stringify(deep_get(cfg, "sctype", "tissue", default="Lung"))
    unknown_min_score = float(deep_get(cfg, "sctype", "unknown_min_score", default=0.0))

    if not db_path.exists():
        return save_status(branch_dir, branch_name, "error", f"ScType DB not found: {db_path}"), None

    normalize_gene_names(adata)
    counts_info = ensure_counts_layer(adata, counts_layer)
    preserve_full_raw_if_missing(adata, counts_layer)
    try:
        cluster_info = ensure_cluster_key(adata, cluster_key, bool(deep_get(cfg, "data_contract", "auto_compute_cluster_key", default=True)))
    except Exception as exc:
        return save_status(branch_dir, branch_name, "error", f"Could not prepare cluster key {cluster_key}: {exc}"), None

    db = pd.read_excel(db_path)
    required_cols = {"tissueType", "cellName", "geneSymbolmore1", "geneSymbolmore2"}
    if not required_cols.issubset(db.columns):
        return save_status(branch_dir, branch_name, "error", f"Unexpected ScType DB columns: {list(db.columns)}"), None

    if tissue and tissue.lower() not in {"all", "auto"}:
        filtered = db[db["tissueType"].astype(str).str.contains(tissue, case=False, na=False)].copy()
        if filtered.empty:
            filtered = db.copy()
    else:
        filtered = db.copy()

    gene_index, matrix = build_gene_index(adata)
    cluster_labels = adata.obs[cluster_key].astype(str).to_numpy()
    unique_clusters = pd.Index(pd.unique(cluster_labels))
    cluster_means: dict[str, np.ndarray] = {}

    if sp.issparse(matrix):
        matrix = matrix.tocsr()

    for cluster in unique_clusters:
        mask = cluster_labels == cluster
        if sp.issparse(matrix):
            mean_vec = np.asarray(matrix[mask].mean(axis=0)).ravel()
        else:
            mean_vec = np.asarray(matrix[mask]).mean(axis=0)
        cluster_means[str(cluster)] = np.asarray(mean_vec, dtype=float)

    score_rows: list[dict[str, Any]] = []
    for row in filtered.itertuples(index=False):
        pos_markers = parse_marker_string(getattr(row, "geneSymbolmore1", None))
        neg_markers = parse_marker_string(getattr(row, "geneSymbolmore2", None))
        pos_idx = [gene_index[g] for g in pos_markers if g in gene_index]
        neg_idx = [gene_index[g] for g in neg_markers if g in gene_index]
        if not pos_idx:
            continue
        for cluster, mean_vec in cluster_means.items():
            pos_score = float(np.mean(mean_vec[pos_idx])) if pos_idx else 0.0
            neg_score = float(np.mean(mean_vec[neg_idx])) if neg_idx else 0.0
            score_rows.append(
                {
                    "cluster": cluster,
                    "cell_name": getattr(row, "cellName", None),
                    "short_name": getattr(row, "shortName", None),
                    "tissueType": getattr(row, "tissueType", None),
                    "score": pos_score - neg_score,
                    "positive_markers_present": len(pos_idx),
                    "negative_markers_present": len(neg_idx),
                }
            )

    if not score_rows:
        return save_status(branch_dir, branch_name, "error", "No ScType marker genes overlapped the dataset"), None

    score_df = pd.DataFrame(score_rows)
    best_df = score_df.sort_values(["cluster", "score"], ascending=[True, False]).groupby("cluster", as_index=False).first()
    best_df["assigned_label"] = best_df["cell_name"].astype(str)
    best_df.loc[best_df["score"] <= unknown_min_score, "assigned_label"] = deep_get(cfg, "data_contract", "unknown_label", default="Unknown")
    cluster_to_label = dict(zip(best_df["cluster"].astype(str), best_df["assigned_label"].astype(str)))
    cluster_to_score = dict(zip(best_df["cluster"].astype(str), best_df["score"].astype(float)))

    obs_out = pd.DataFrame(index=adata.obs_names)
    obs_out["sctype_annotation"] = pd.Series(cluster_labels, index=adata.obs_names).map(cluster_to_label).fillna(deep_get(cfg, "data_contract", "unknown_label", default="Unknown")).astype(str).values
    obs_out["sctype_score"] = pd.Series(cluster_labels, index=adata.obs_names).map(cluster_to_score).astype(float).values

    score_df.to_csv(branch_dir / "sctype_cluster_scores.tsv", sep="\t", index=False)
    best_df.to_csv(branch_dir / "sctype_cluster_best.tsv", sep="\t", index=False)

    annotation_path = write_annotation_table(obs_out, branch_dir / "sctype_annotations.tsv.gz")
    output_path: Path | None = None
    if should_write_branch_h5ad(cfg):
        adata.obs["sctype_annotation"] = obs_out["sctype_annotation"].values
        adata.obs["sctype_score"] = pd.to_numeric(obs_out["sctype_score"], errors="coerce").values
        adata.uns["sctype_db_path"] = str(db_path)
        adata.uns["sctype_tissue_requested"] = tissue
        adata.uns["sctype_cluster_key"] = cluster_key
        adata.uns["sctype_cluster_key_source"] = cluster_key_source
        adata.uns["sctype_counts_layer_source"] = counts_info.get("source")
        output_path = branch_dir / "sctype_annotation.h5ad"
        adata.write_h5ad(output_path, compression="gzip")

    n_unique_labels = int(pd.Series(obs_out["sctype_annotation"].astype(str)).nunique())
    del score_df, best_df
    gc.collect()
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"ScType completed using {cluster_key} and tissue={tissue}",
        output_h5ad=str(output_path) if output_path is not None else None,
        annotation_table=str(annotation_path),
        db_path=str(db_path),
        cluster_key=cluster_key,
        cluster_key_source=cluster_key_source,
        cluster_key_created=bool(cluster_info.get("created", False)),
        n_unique_labels=n_unique_labels,
    ), {"annotations_tsv": str(annotation_path), "h5ad": str(output_path) if output_path is not None else None}


def run_mllmcelltype_branch(preflight_h5ad: Path, cfg: dict[str, Any], branch_dir: Path) -> tuple[dict[str, Any], dict[str, str | None] | None]:
    branch_name = "mllmcelltype"
    if not bool(deep_get(cfg, "mllmcelltype", "enabled", default=True)):
        return save_status(branch_dir, branch_name, "skipped", "mLLMCelltype branch disabled in config"), None

    adapter_mode = stringify(deep_get(cfg, "mllmcelltype", "adapter_mode", default="auto")).strip().lower() or "auto"
    precomputed_h5ad = stringify(deep_get(cfg, "mllmcelltype", "precomputed_h5ad", default="")).strip()
    repo_dir_raw = stringify(deep_get(cfg, "mllmcelltype", "repo_dir", default="")).strip()
    module_path = stringify(deep_get(cfg, "mllmcelltype", "module_path", default="")).strip()
    callable_name = stringify(deep_get(cfg, "mllmcelltype", "callable_name", default="")).strip()
    label_column = stringify(deep_get(cfg, "mllmcelltype", "label_column", default="mllmcelltype_annotation")).strip() or "mllmcelltype_annotation"
    confidence_column = stringify(deep_get(cfg, "mllmcelltype", "confidence_column", default="")).strip() or None
    score_column = stringify(deep_get(cfg, "mllmcelltype", "score_column", default="")).strip() or None
    output_name = stringify(deep_get(cfg, "mllmcelltype", "output_name", default="mllmcelltype_annotation.h5ad")).strip() or "mllmcelltype_annotation.h5ad"
    output_path = branch_dir / output_name
    annotation_path = branch_dir / "mllmcelltype_annotations.tsv.gz"
    unknown_label = deep_get(cfg, "data_contract", "unknown_label", default="Unknown")
    write_branch_h5ad = should_write_branch_h5ad(cfg)

    repo_dir = Path(repo_dir_raw) if repo_dir_raw else None
    requested_kwargs = deep_get(cfg, "mllmcelltype", "function_kwargs", default={}) or {}
    if not isinstance(requested_kwargs, dict):
        return save_status(branch_dir, branch_name, "error", "mLLMCelltype function_kwargs must be a mapping"), None

    mode_candidates: list[str]
    if adapter_mode in {"precomputed", "precomputed_h5ad"}:
        mode_candidates = ["precomputed_h5ad"]
    elif adapter_mode in {"callable", "function"}:
        mode_candidates = ["callable"]
    else:
        mode_candidates = ["precomputed_h5ad", "callable"]

    result_adata: ad.AnnData | None = None
    result_payload: dict[str, Any] | None = None
    source_mode_used: str | None = None

    try:
        for candidate in mode_candidates:
            if candidate == "precomputed_h5ad" and precomputed_h5ad:
                precomputed_path = Path(precomputed_h5ad)
                if not precomputed_path.exists():
                    return save_status(
                        branch_dir,
                        branch_name,
                        "error",
                        f"Configured precomputed_h5ad does not exist: {precomputed_path}",
                        adapter_mode=adapter_mode,
                    ), None
                result_payload = {"kind": "h5ad", "path": precomputed_path}
                source_mode_used = "precomputed_h5ad"
                break

            if candidate == "callable" and module_path and callable_name:
                callable_func = load_configured_callable(module_path, callable_name, repo_dir=repo_dir)
                candidate_kwargs = {
                    "adata_path": str(preflight_h5ad),
                    "input_h5ad": str(preflight_h5ad),
                    "output_h5ad": str(output_path),
                    "output_path": str(output_path),
                    "branch_dir": str(branch_dir),
                    "config": deep_get(cfg, "mllmcelltype", default={}) or {},
                    "workflow_config": cfg,
                    "unknown_label": unknown_label,
                    "write_h5ad": write_branch_h5ad,
                    "annotation_table_path": str(annotation_path),
                }
                candidate_kwargs.update(requested_kwargs)
                result = call_with_supported_kwargs(callable_func, candidate_kwargs)
                result_payload = resolve_callable_result(result, output_path)
                source_mode_used = "callable"
                break
    except Exception as exc:
        tb_path = branch_dir / "mllmcelltype_traceback.txt"
        tb_path.write_text(traceback.format_exc(), encoding="utf-8")
        return save_status(
            branch_dir,
            branch_name,
            "error",
            f"mLLMCelltype adapter failed: {exc}",
            adapter_mode=adapter_mode,
            traceback_file=str(tb_path),
            repo_dir=str(repo_dir) if repo_dir is not None else None,
            module_path=module_path or None,
            callable_name=callable_name or None,
        ), None

    if result_payload is None:
        return save_status(
            branch_dir,
            branch_name,
            "skipped",
            "mLLMCelltype branch has no runnable source; set precomputed_h5ad or module_path+callable_name",
            adapter_mode=adapter_mode,
            repo_dir=str(repo_dir) if repo_dir is not None else None,
            module_path=module_path or None,
            callable_name=callable_name or None,
            precomputed_h5ad=precomputed_h5ad or None,
        ), None

    base_backed = ad.read_h5ad(preflight_h5ad, backed="r")
    try:
        base_obs_names = pd.Index(base_backed.obs_names.astype(str))
    finally:
        close_backed_adata(base_backed)

    result_kind = stringify(result_payload.get("kind", "")).strip().lower()
    label_col_for_read = stringify(result_payload.get("label_column", label_column)).strip() or label_column
    confidence_col_for_read = stringify(result_payload.get("confidence_column", confidence_column)).strip() or confidence_column
    score_col_for_read = stringify(result_payload.get("score_column", score_column)).strip() or score_column

    if result_kind == "adata":
        result_adata = result_payload.get("adata")
        if not isinstance(result_adata, ad.AnnData):
            return save_status(branch_dir, branch_name, "error", "Callable result kind=adata but payload had no AnnData"), None
        result_obs = result_adata.obs.copy()
    elif result_kind == "h5ad":
        result_path = Path(result_payload["path"])
        result_obs = extract_obs_columns_from_h5ad(result_path, columns=[label_col_for_read, confidence_col_for_read, score_col_for_read])
    elif result_kind == "table":
        result_path = Path(result_payload["path"])
        result_obs = read_annotation_table(result_path)
    else:
        return save_status(branch_dir, branch_name, "error", f"Unsupported resolved mLLMCelltype result kind: {result_kind}"), None

    if label_col_for_read not in result_obs.columns:
        return save_status(
            branch_dir,
            branch_name,
            "error",
            f"Configured mLLMCelltype label column not found: {label_col_for_read}",
            adapter_mode_used=source_mode_used,
            available_obs_columns=list(result_obs.columns),
        ), None

    overlap_cells = int(len(base_obs_names.intersection(pd.Index(result_obs.index.astype(str)))))
    if overlap_cells == 0:
        return save_status(
            branch_dir,
            branch_name,
            "error",
            "mLLMCelltype result has zero overlapping cells with workflow input",
            adapter_mode_used=source_mode_used,
            label_column=label_col_for_read,
        ), None

    normalized_labels = normalize_annotation_values(result_obs[label_col_for_read], unknown_label)
    normalized_labels = normalized_labels.reindex(base_obs_names).fillna(unknown_label)
    obs_out = pd.DataFrame(index=base_obs_names)
    obs_out["mllmcelltype_annotation"] = normalized_labels.astype(str).values

    if confidence_col_for_read and confidence_col_for_read in result_obs.columns:
        confidence = pd.to_numeric(result_obs[confidence_col_for_read], errors="coerce").reindex(base_obs_names)
        obs_out["mllmcelltype_confidence_scextract"] = confidence.values

    if score_col_for_read and score_col_for_read in result_obs.columns:
        score = pd.to_numeric(result_obs[score_col_for_read], errors="coerce").reindex(base_obs_names)
        obs_out["mllmcelltype_score"] = score.values

    annotation_path = write_annotation_table(obs_out, annotation_path)

    output_h5ad_for_status: str | None = None
    if write_branch_h5ad and result_kind != "h5ad":
        base_adata = ad.read_h5ad(preflight_h5ad)
        base_adata.obs["mllmcelltype_annotation"] = obs_out["mllmcelltype_annotation"].values
        if "mllmcelltype_confidence_scextract" in obs_out.columns:
            base_adata.obs["mllmcelltype_confidence_scextract"] = pd.to_numeric(obs_out["mllmcelltype_confidence_scextract"], errors="coerce").values
        if "mllmcelltype_score" in obs_out.columns:
            base_adata.obs["mllmcelltype_score"] = pd.to_numeric(obs_out["mllmcelltype_score"], errors="coerce").values
        base_adata.uns["mllmcelltype_source_mode"] = source_mode_used
        base_adata.uns["mllmcelltype_label_column_source"] = label_col_for_read
        if confidence_col_for_read:
            base_adata.uns["mllmcelltype_confidence_column_source"] = confidence_col_for_read
        if score_col_for_read:
            base_adata.uns["mllmcelltype_score_column_source"] = score_col_for_read
        if repo_dir is not None:
            base_adata.uns["mllmcelltype_repo_dir"] = str(repo_dir)
        if module_path:
            base_adata.uns["mllmcelltype_module_path"] = module_path
        if callable_name:
            base_adata.uns["mllmcelltype_callable_name"] = callable_name
        base_adata.write_h5ad(output_path, compression="gzip")
        output_h5ad_for_status = str(output_path)
    elif result_kind == "h5ad":
        output_h5ad_for_status = str(result_payload["path"])

    missing_after_reindex = int((normalized_labels.astype(str) == unknown_label).sum())
    return save_status(
        branch_dir,
        branch_name,
        "ok",
        f"mLLMCelltype completed via {source_mode_used}",
        output_h5ad=output_h5ad_for_status,
        annotation_table=str(annotation_path),
        adapter_mode=adapter_mode,
        adapter_mode_used=source_mode_used,
        label_column=label_col_for_read,
        confidence_column=confidence_col_for_read,
        score_column=score_col_for_read,
        overlap_cells=overlap_cells,
        n_unique_labels=int(pd.Series(obs_out["mllmcelltype_annotation"].astype(str)).nunique()),
        n_unknown_after_reindex=missing_after_reindex,
    ), {"annotations_tsv": str(annotation_path), "h5ad": output_h5ad_for_status}


def merge_branch_outputs(
    preflight_h5ad: Path,
    branch_outputs: dict[str, dict[str, str | None] | None],
    cfg: dict[str, Any],
    merged_dir: Path,
) -> tuple[Path, dict[str, Any]]:
    merged = ad.read_h5ad(preflight_h5ad)
    aliases = deep_get(cfg, "aliases", default={})
    merge_summary: dict[str, Any] = {}

    sc_obs = load_branch_annotation_obs(branch_outputs.get("scExtract_builtin"), fallback_columns=["leiden", "louvain", "Certainty", "Source", "Tissue"])
    if sc_obs is not None and not sc_obs.empty:
        label_col = infer_scextract_label_column_from_obs(sc_obs)
        alias_col = aliases.get("scExtract", "cell_type_scextract")
        merged.obs[alias_col] = sc_obs[label_col].astype(str).reindex(merged.obs_names).values
        for col in ("Certainty", "Source", "Tissue"):
            if col in sc_obs.columns:
                merged.obs[f"scextract_{col.lower()}"] = sc_obs[col].astype(str).reindex(merged.obs_names).values
        merge_summary["scExtract"] = {"label_column": label_col, "alias_column": alias_col}

    ct_obs = load_branch_annotation_obs(branch_outputs.get("celltypist"), fallback_columns=["celltypist_annotation", "celltypist_majority_voting", "celltypist_confidence_scextract"])
    if ct_obs is not None and not ct_obs.empty:
        alias_col = aliases.get("celltypist_annotation", "cell_type_celltypist_scextract")
        merged.obs[alias_col] = ct_obs["celltypist_annotation"].astype(str).reindex(merged.obs_names).values
        if "celltypist_majority_voting" in ct_obs.columns:
            merged.obs["celltypist_majority_voting_scextract"] = ct_obs["celltypist_majority_voting"].astype(str).reindex(merged.obs_names).values
        if "celltypist_confidence_scextract" in ct_obs.columns:
            merged.obs["celltypist_confidence_scextract"] = pd.to_numeric(ct_obs["celltypist_confidence_scextract"], errors="coerce").reindex(merged.obs_names).values
        merge_summary["celltypist"] = {"alias_column": alias_col}

    sctype_obs = load_branch_annotation_obs(branch_outputs.get("sctype"), fallback_columns=["sctype_annotation", "sctype_score"])
    if sctype_obs is not None and not sctype_obs.empty:
        alias_col = aliases.get("sctype_annotation", "cell_type_sctype_scextract")
        merged.obs[alias_col] = sctype_obs["sctype_annotation"].astype(str).reindex(merged.obs_names).values
        if "sctype_score" in sctype_obs.columns:
            merged.obs["sctype_score"] = pd.to_numeric(sctype_obs["sctype_score"], errors="coerce").reindex(merged.obs_names).values
        merge_summary["sctype"] = {"alias_column": alias_col}

    mllm_obs = load_branch_annotation_obs(branch_outputs.get("mllmcelltype"), fallback_columns=["mllmcelltype_annotation", "mllmcelltype_confidence_scextract", "mllmcelltype_score"])
    if mllm_obs is not None and not mllm_obs.empty:
        alias_col = aliases.get("mllmcelltype_annotation", "cell_type_mllmcelltype_scextract")
        merged.obs[alias_col] = mllm_obs["mllmcelltype_annotation"].astype(str).reindex(merged.obs_names).values
        if "mllmcelltype_confidence_scextract" in mllm_obs.columns:
            merged.obs["mllmcelltype_confidence_scextract"] = pd.to_numeric(mllm_obs["mllmcelltype_confidence_scextract"], errors="coerce").reindex(merged.obs_names).values
        if "mllmcelltype_score" in mllm_obs.columns:
            merged.obs["mllmcelltype_score"] = pd.to_numeric(mllm_obs["mllmcelltype_score"], errors="coerce").reindex(merged.obs_names).values
        merge_summary["mllmcelltype"] = {"alias_column": alias_col}

    selected_label_source = deep_get(cfg, "run", "selected_label_source", default="compare_only")
    confidence_column, min_confidence = resolve_bridge_confidence_settings(cfg, selected_label_source)

    bridge_summary = apply_bridge(
        merged,
        selected_label_source=selected_label_source,
        unknown_label=deep_get(cfg, "data_contract", "unknown_label", default="Unknown"),
        min_cells_per_type=int(deep_get(cfg, "data_contract", "min_cells_per_type", default=10)),
        confidence_column=confidence_column,
        min_confidence=min_confidence,
    )
    merge_summary["bridge"] = bridge_summary

    normalize_obs_strings(merged)

    output_path = merged_dir / "scextract_three_way_merged.h5ad"
    merged.write_h5ad(output_path, compression="gzip")
    return output_path, merge_summary


def compute_pairwise_agreement(adata: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for left, right in combinations(columns, 2):
        s1 = adata[left].astype("string") if left in adata.columns else pd.Series(pd.NA, index=adata.index)
        s2 = adata[right].astype("string") if right in adata.columns else pd.Series(pd.NA, index=adata.index)
        valid = s1.notna() & s2.notna() & (s1.astype(str) != "") & (s2.astype(str) != "")
        if valid.any():
            agreement = float((s1[valid].astype(str) == s2[valid].astype(str)).mean())
            compared = int(valid.sum())
        else:
            agreement = np.nan
            compared = 0
        rows.append(
            {
                "left": left,
                "right": right,
                "left_available": left in adata.columns,
                "right_available": right in adata.columns,
                "n_compared": compared,
                "exact_match_rate": agreement,
            }
        )
    return pd.DataFrame(rows)


def compute_truth_benchmark(adata: pd.DataFrame, columns: list[str], true_group_key: str) -> pd.DataFrame:
    if true_group_key not in adata.columns:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    truth = adata[true_group_key].astype(str)
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    except Exception:
        adjusted_rand_score = None
        normalized_mutual_info_score = None

    for col in columns:
        if col not in adata.columns:
            continue
        pred = adata[col].astype(str)
        valid = truth.notna() & pred.notna() & (truth != "") & (pred != "")
        if not valid.any():
            continue
        truth_valid = truth[valid].astype(str)
        pred_valid = pred[valid].astype(str)
        exact = float((truth_valid.str.strip().str.lower() == pred_valid.str.strip().str.lower()).mean())
        row: dict[str, Any] = {"prediction_column": col, "n_compared": int(valid.sum()), "exact_match_rate": exact}
        if adjusted_rand_score is not None and normalized_mutual_info_score is not None:
            row["ari"] = float(adjusted_rand_score(truth_valid, pred_valid))
            row["nmi"] = float(normalized_mutual_info_score(truth_valid, pred_valid))
        rows.append(row)
    return pd.DataFrame(rows)


def build_summary_outputs(merged_h5ad: Path, cfg: dict[str, Any], summary_dir: Path) -> dict[str, Any]:
    backed = ad.read_h5ad(merged_h5ad, backed="r")
    try:
        obs = backed.obs.copy()
        n_obs = int(backed.n_obs)
    finally:
        close_backed_adata(backed)

    columns = [
        deep_get(cfg, "aliases", "scExtract", default="cell_type_scextract"),
        deep_get(cfg, "aliases", "celltypist_annotation", default="cell_type_celltypist_scextract"),
        deep_get(cfg, "aliases", "sctype_annotation", default="cell_type_sctype_scextract"),
        deep_get(cfg, "aliases", "mllmcelltype_annotation", default="cell_type_mllmcelltype_scextract"),
    ]
    available_columns = [col for col in columns if col in obs.columns]

    branch_rows = []
    count_tables: dict[str, pd.DataFrame] = {}
    for col in columns:
        if col in obs.columns:
            series = obs[col].astype("string")
            missing = int(series.isna().sum() + (series.astype(str) == "").sum())
            branch_rows.append(
                {
                    "column": col,
                    "available": True,
                    "n_cells": n_obs,
                    "n_missing": missing,
                    "n_unique_labels": int(series.dropna().astype(str).replace("", pd.NA).dropna().nunique()),
                }
            )
            counts = series.astype(str).value_counts().rename_axis("label").reset_index(name="count")
        else:
            branch_rows.append(
                {
                    "column": col,
                    "available": False,
                    "n_cells": n_obs,
                    "n_missing": n_obs,
                    "n_unique_labels": 0,
                }
            )
            counts = pd.DataFrame(columns=["label", "count"])
        count_tables[col] = counts.copy()
        write_tsv(counts, summary_dir / f"{col}_counts.tsv")

    pairwise = compute_pairwise_agreement(obs, columns)
    write_tsv(pairwise, summary_dir / "pairwise_agreement.tsv")

    true_group_key = deep_get(cfg, "data_contract", "true_group_key", default=None)
    truth_df = compute_truth_benchmark(obs, columns, true_group_key) if true_group_key else pd.DataFrame()
    if not truth_df.empty:
        write_tsv(truth_df, summary_dir / "truth_benchmark.tsv")

    branch_df = pd.DataFrame(branch_rows)
    write_tsv(branch_df, summary_dir / "branch_summary.tsv")
    figure_rows, figure_manifest_path = generate_summary_visualizations(
        obs=obs,
        count_tables=count_tables,
        columns=columns,
        available_columns=available_columns,
        pairwise=pairwise,
        truth_df=truth_df,
        merged_h5ad=merged_h5ad,
        summary_dir=summary_dir,
        cfg=cfg,
    )
    summary = {
        "columns": columns,
        "available_columns": available_columns,
        "branch_summary_tsv": str(summary_dir / "branch_summary.tsv"),
        "pairwise_agreement_tsv": str(summary_dir / "pairwise_agreement.tsv"),
        "truth_benchmark_tsv": str(summary_dir / "truth_benchmark.tsv") if not truth_df.empty else None,
        "figure_manifest_tsv": str(figure_manifest_path),
        "figures_dir": str(summary_dir / "figures"),
        "n_figures_ok": int(sum(1 for row in figure_rows if row.get("status") == "ok")),
    }
    write_json(summary, summary_dir / "summary_manifest.json")
    return summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run annotation consistency workflow across scExtract/CellTypist/ScType/mLLMCelltype.")
    parser.add_argument("--adata", required=True, help="Input AnnData h5ad path")
    parser.add_argument("--pdf", required=True, help="Source article PDF/txt path")
    parser.add_argument("--config-yaml", required=True, help="Workflow YAML config")
    parser.add_argument("--output-dir", default=None, help="Explicit run directory (optional)")
    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    config_path = Path(args.config_yaml)
    cfg = load_config(config_path)
    dotenv_loaded = load_workspace_dotenv()
    run_dir = resolve_run_dir(cfg, args.output_dir)
    runtime_dir = ensure_dir(run_dir / "runtime")
    summary_dir = ensure_dir(run_dir / "summary")
    merged_dir = ensure_dir(run_dir / "merged")

    branch_dirs = {
        "scExtract_builtin": ensure_dir(run_dir / "branches" / "scExtract_builtin"),
        "celltypist": ensure_dir(run_dir / "branches" / "celltypist"),
        "sctype": ensure_dir(run_dir / "branches" / "sctype"),
        "mllmcelltype": ensure_dir(run_dir / "branches" / "mllmcelltype"),
    }

    manifest: dict[str, Any] = {
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "input_h5ad": str(Path(args.adata).resolve()),
        "input_pdf": str(Path(args.pdf).resolve()),
        "config_yaml": str(config_path.resolve()),
        "run_dir": str(run_dir.resolve()),
        "dotenv_loaded": dotenv_loaded,
        "dotenv_path": str((WORKSPACE_ROOT / ".env").resolve()),
        "scextract_python": resolve_scextract_python(cfg),
        "deepseek_env_key": deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY"),
        "deepseek_env_present": bool(os.getenv(deep_get(cfg, "deepseek", "env_key", default="DEEPSEEK_API_KEY"))),
        "configured_sources": [
            "scExtract",
            "celltypist_annotation",
            "sctype_annotation",
            "mllmcelltype_annotation",
        ],
        "mllmcelltype_enabled": bool(deep_get(cfg, "mllmcelltype", "enabled", default=True)),
    }
    write_json(manifest, run_dir / "run_manifest.json")

    preflight_h5ad, preflight_summary = prepare_preflight_adata(Path(args.adata), cfg, runtime_dir)
    manifest["preflight"] = preflight_summary
    write_json(manifest, run_dir / "run_manifest.json")

    statuses: list[dict[str, Any]] = []
    branch_outputs: dict[str, dict[str, str | None] | None] = {}
    branch_runners = [
        ("scExtract_builtin", lambda: run_scextract_branch(preflight_h5ad, Path(args.pdf), cfg, branch_dirs["scExtract_builtin"])),
        ("celltypist", lambda: run_celltypist_branch(preflight_h5ad, cfg, branch_dirs["celltypist"], source_h5ad_path=Path(args.adata))),
        ("sctype", lambda: run_sctype_branch(preflight_h5ad, cfg, branch_dirs["sctype"])),
        ("mllmcelltype", lambda: run_mllmcelltype_branch(preflight_h5ad, cfg, branch_dirs["mllmcelltype"])),
    ]

    for branch_name, runner in branch_runners:
        try:
            status, output_path = runner()
        except Exception as exc:  # pragma: no cover - runtime guard
            tb_path = branch_dirs[branch_name] / f"{branch_name}_traceback.txt"
            tb_path.write_text(traceback.format_exc(), encoding="utf-8")
            status = save_status(branch_dirs[branch_name], branch_name, "error", str(exc), traceback_file=str(tb_path))
            output_path = None
        statuses.append(status)
        branch_outputs[branch_name] = output_path
        gc.collect()

    merged_h5ad, merge_summary = merge_branch_outputs(preflight_h5ad, branch_outputs, cfg, merged_dir)
    summary_manifest = build_summary_outputs(merged_h5ad, cfg, summary_dir)

    manifest.update(
        {
            "finished_at": datetime.now().isoformat(timespec="seconds"),
            "statuses": statuses,
            "branch_outputs": branch_outputs,
            "merged_h5ad": str(merged_h5ad),
            "merge_summary": merge_summary,
            "summary_manifest": summary_manifest,
        }
    )
    write_json(manifest, run_dir / "run_manifest.json")

    pd.DataFrame(statuses).to_csv(run_dir / "branch_status.tsv", sep="\t", index=False)
    print(json.dumps({
        "run_dir": str(run_dir),
        "merged_h5ad": str(merged_h5ad),
        "status_tsv": str(run_dir / "branch_status.tsv"),
        "summary_manifest": str(summary_dir / "summary_manifest.json"),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
