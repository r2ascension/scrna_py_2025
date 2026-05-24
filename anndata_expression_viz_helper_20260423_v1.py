#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Generic AnnData expression + visualization helper.

This module extracts the reusable parts of lineage-specific notebook analysis:
- choosing expression source from `.raw` or `.X`
- checking marker availability
- safely extracting dense/sparse gene vectors
- building grouped mean/pct expression matrices
- organizing output folders and writing simple artifacts
- lightweight wrappers for Scanpy dotplots and seaborn heatmaps
- building generic review/dotplot tables for cluster-vs-marker comparisons

The module is intentionally import-safe: heavy plotting dependencies are imported
lazily inside rendering functions.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import sparse


def _is_categorical_series(series: pd.Series) -> bool:
    return isinstance(series.dtype, pd.CategoricalDtype)

__all__ = [
    "deduplicate_preserve_order",
    "normalize_string_series",
    "safe_obs_column",
    "ensure_categorical_with_order",
    "resolve_expression_source",
    "gene_available",
    "extract_gene_vector",
    "get_available_markers",
    "build_mean_expression_matrix",
    "build_pct_expression_matrix",
    "build_review_dotplot_table",
    "prepare_output_dirs",
    "write_tsv",
    "write_json",
    "save_current_figure",
    "compute_dynamic_figsize",
    "render_scanpy_dotplot",
    "render_seaborn_heatmap",
]


def deduplicate_preserve_order(items: Iterable[Any]) -> list[Any]:
    seen: set[Any] = set()
    ordered: list[Any] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def normalize_string_series(series: pd.Series) -> pd.Series:
    normalized = pd.Series(series, copy=True)
    return normalized.astype("string").str.strip()


def safe_obs_column(adata: Any, column: str) -> pd.Series:
    if column not in adata.obs.columns:
        raise KeyError(f"Column not found in adata.obs: {column}")
    return adata.obs[column]


def ensure_categorical_with_order(series: pd.Series, categories: Sequence[str] | None = None) -> pd.Series:
    if categories is None:
        return pd.Series(pd.Categorical(series), index=series.index, name=series.name)
    return pd.Series(
        pd.Categorical(series, categories=list(categories), ordered=True),
        index=series.index,
        name=series.name,
    )


def resolve_expression_source(adata: Any, prefer_raw: bool = True) -> dict[str, Any]:
    has_raw = getattr(adata, "raw", None) is not None and getattr(adata.raw, "var_names", None) is not None
    if prefer_raw and has_raw:
        gene_universe = pd.Index(adata.raw.var_names)
        return {
            "use_raw": True,
            "source_name": ".raw",
            "gene_universe": gene_universe,
            "n_genes": int(len(gene_universe)),
        }

    gene_universe = pd.Index(adata.var_names)
    return {
        "use_raw": False,
        "source_name": ".X",
        "gene_universe": gene_universe,
        "n_genes": int(len(gene_universe)),
    }


def gene_available(adata: Any, gene: str) -> bool:
    source = resolve_expression_source(adata)
    if gene in source["gene_universe"]:
        return True
    return gene in pd.Index(adata.var_names)


def _to_1d_array(values: Any) -> np.ndarray:
    if sparse.issparse(values):
        return values.toarray().ravel()
    return np.asarray(values).ravel()


def extract_gene_vector(adata: Any, gene: str, use_raw: bool | None = None) -> np.ndarray | None:
    raw = getattr(adata, "raw", None)
    raw_has_gene = raw is not None and gene in pd.Index(raw.var_names)
    x_has_gene = gene in pd.Index(adata.var_names)

    if use_raw is True and raw_has_gene:
        return _to_1d_array(raw[:, gene].X)
    if use_raw is False and x_has_gene:
        return _to_1d_array(adata[:, gene].X)

    if raw_has_gene:
        return _to_1d_array(raw[:, gene].X)
    if x_has_gene:
        return _to_1d_array(adata[:, gene].X)
    return None


def get_available_markers(
    marker_dict: Mapping[str, Sequence[str]],
    var_names: Sequence[str] | pd.Index,
) -> tuple[dict[str, list[str]], list[str]]:
    universe = set(map(str, var_names))
    available: dict[str, list[str]] = {}
    missing_all: list[str] = []

    for label, markers in marker_dict.items():
        keep = [marker for marker in markers if str(marker) in universe]
        missing = [marker for marker in markers if str(marker) not in universe]
        if keep:
            available[str(label)] = keep
        missing_all.extend(map(str, missing))

    return available, deduplicate_preserve_order(missing_all)


def _resolve_group_order(series: pd.Series, categories: Sequence[str] | None = None) -> list[str]:
    if categories is not None:
        return [str(cat) for cat in categories]
    if _is_categorical_series(series):
        return [str(cat) for cat in series.cat.categories]
    values = normalize_string_series(series).dropna().tolist()
    return sorted(set(values))


def build_mean_expression_matrix(
    adata: Any,
    genes: Sequence[str],
    groupby: str,
    use_raw: bool | None = None,
    categories: Sequence[str] | None = None,
) -> pd.DataFrame:
    group_series = safe_obs_column(adata, groupby)
    ordered_groups = _resolve_group_order(group_series, categories=categories)
    valid_genes = [gene for gene in deduplicate_preserve_order(genes) if gene_available(adata, gene)]
    if not valid_genes:
        raise ValueError("No requested genes are available in the AnnData object")

    group_clean = normalize_string_series(group_series)
    matrix: list[list[float]] = []
    row_names: list[str] = []

    for group in ordered_groups:
        mask = (group_clean == str(group)).to_numpy()
        if mask.sum() == 0:
            continue
        row = []
        for gene in valid_genes:
            values = extract_gene_vector(adata, gene, use_raw=use_raw)
            if values is None:
                row.append(np.nan)
                continue
            row.append(float(np.asarray(values)[mask].mean()))
        matrix.append(row)
        row_names.append(str(group))

    return pd.DataFrame(matrix, index=row_names, columns=valid_genes)


def build_pct_expression_matrix(
    adata: Any,
    genes: Sequence[str],
    groupby: str,
    use_raw: bool | None = None,
    categories: Sequence[str] | None = None,
) -> pd.DataFrame:
    group_series = safe_obs_column(adata, groupby)
    ordered_groups = _resolve_group_order(group_series, categories=categories)
    valid_genes = [gene for gene in deduplicate_preserve_order(genes) if gene_available(adata, gene)]
    if not valid_genes:
        raise ValueError("No requested genes are available in the AnnData object")

    group_clean = normalize_string_series(group_series)
    matrix: list[list[float]] = []
    row_names: list[str] = []

    for group in ordered_groups:
        mask = (group_clean == str(group)).to_numpy()
        if mask.sum() == 0:
            continue
        row = []
        for gene in valid_genes:
            values = extract_gene_vector(adata, gene, use_raw=use_raw)
            if values is None:
                row.append(0.0)
                continue
            selected = np.asarray(values)[mask]
            row.append(float((selected > 0).mean()))
        matrix.append(row)
        row_names.append(str(group))

    return pd.DataFrame(matrix, index=row_names, columns=valid_genes)


def _minmax_scale(series: pd.Series) -> pd.Series:
    series = series.astype(float)
    finite = series[np.isfinite(series)]
    if finite.empty:
        return pd.Series(np.zeros(len(series)), index=series.index)
    smin = finite.min()
    smax = finite.max()
    if smax == smin:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (series - smin) / (smax - smin)


def build_review_dotplot_table(
    adata: Any,
    review_clusters: Sequence[str],
    panel_specs: Sequence[Mapping[str, Any]],
    groupby_col: str,
    use_raw: bool | None = None,
) -> pd.DataFrame:
    review_clusters = [str(cluster) for cluster in review_clusters]
    group_clean = normalize_string_series(safe_obs_column(adata, groupby_col))
    masks = {cluster: (group_clean == cluster).to_numpy() for cluster in review_clusters}

    rows: list[dict[str, Any]] = []
    plot_label_order: list[str] = []

    for spec in panel_specs:
        spec_cluster = str(spec.get("cluster", spec.get("spec_cluster", "")))
        side = str(spec.get("side", "panel"))
        label = str(spec.get("label", spec.get(f"{side}_label", "")))
        note = str(spec.get("note", ""))
        genes = [gene for gene in spec.get("genes", []) if gene_available(adata, gene)]
        for gene in genes:
            plot_label = f"{spec_cluster}|{side}|{gene}"
            plot_label_order.append(plot_label)
            values = extract_gene_vector(adata, gene, use_raw=use_raw)
            if values is None:
                continue
            values = np.asarray(values)
            for cluster in review_clusters:
                mask = masks[cluster]
                if mask.sum() == 0:
                    mean_expr = np.nan
                    pct_expr = 0.0
                else:
                    selected = values[mask]
                    mean_expr = float(np.log1p(selected).mean())
                    pct_expr = float((selected > 0).mean())
                rows.append(
                    {
                        "review_cluster": cluster,
                        "spec_cluster": spec_cluster,
                        "side": side,
                        "label": label,
                        "note": note,
                        "gene": gene,
                        "plot_label": plot_label,
                        "mean_expr": mean_expr,
                        "pct_expr": pct_expr,
                    }
                )

    plot_df = pd.DataFrame(rows)
    if plot_df.empty:
        return plot_df

    plot_df["mean_expr_scaled"] = plot_df.groupby("plot_label")["mean_expr"].transform(_minmax_scale)
    x_map = {label: idx for idx, label in enumerate(deduplicate_preserve_order(plot_label_order))}
    y_map = {cluster: idx for idx, cluster in enumerate(review_clusters)}
    plot_df["x"] = plot_df["plot_label"].map(x_map)
    plot_df["y"] = plot_df["review_cluster"].map(y_map)
    return plot_df


def prepare_output_dirs(
    base_dir: str | Path,
    subdirs: Sequence[str] = ("figures", "tables", "artifacts"),
) -> dict[str, Path]:
    base = Path(base_dir)
    base.mkdir(parents=True, exist_ok=True)
    result = {"base": base}
    for name in subdirs:
        path = base / str(name)
        path.mkdir(parents=True, exist_ok=True)
        result[str(name)] = path
    return result


def write_tsv(df: pd.DataFrame, path: str | Path, index: bool = True) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=index)
    return out


def write_json(payload: Mapping[str, Any], path: str | Path) -> Path:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
    return out


def save_current_figure(
    path: str | Path,
    dpi: int = 300,
    close: bool = True,
    figure: Any | None = None,
) -> Path:
    import matplotlib.pyplot as plt

    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    if figure is None:
        plt.savefig(out, dpi=dpi, bbox_inches="tight")
    else:
        figure.savefig(out, dpi=dpi, bbox_inches="tight")
    if close:
        if figure is None:
            plt.close()
        else:
            plt.close(figure)
    return out


def compute_dynamic_figsize(n_rows: int, n_cols: int, kind: str = "heatmap") -> tuple[float, float]:
    if kind == "dotplot":
        return max(12.0, n_cols * 0.3), max(6.0, n_rows * 0.4)
    if kind == "review":
        return max(18.0, n_cols * 0.42), max(6.0, n_rows * 0.72)
    return max(12.0, n_cols * 0.4), max(10.0, n_rows * 0.15)


def render_scanpy_dotplot(
    adata: Any,
    var_names: Sequence[str],
    groupby: str,
    output_path: str | Path | None = None,
    use_raw: bool | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 300,
    close: bool = True,
    **kwargs: Any,
) -> Any:
    import matplotlib.pyplot as plt
    import scanpy as sc

    dotplot = sc.pl.dotplot(
        adata,
        var_names=list(var_names),
        groupby=groupby,
        use_raw=use_raw,
        show=False,
        return_fig=True,
        figsize=figsize,
        **kwargs,
    )
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        if hasattr(dotplot, "savefig"):
            dotplot.savefig(str(out), dpi=dpi, bbox_inches="tight")
        else:
            plt.savefig(out, dpi=dpi, bbox_inches="tight")
        if close:
            plt.close("all")
    return dotplot


def render_seaborn_heatmap(
    dataframe: pd.DataFrame,
    output_path: str | Path | None = None,
    title: str | None = None,
    cmap: str = "RdYlBu_r",
    center: float | None = 0.0,
    robust: bool = True,
    cbar_label: str | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 300,
    close: bool = True,
    **kwargs: Any,
) -> Any:
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=figsize)
    heatmap = sns.heatmap(
        dataframe,
        cmap=cmap,
        center=center,
        robust=robust,
        yticklabels=True,
        xticklabels=True,
        cbar_kws={"label": cbar_label} if cbar_label else None,
        linewidths=0.1,
        linecolor="lightgray",
        ax=ax,
        **kwargs,
    )
    if title:
        ax.set_title(title, fontsize=14, pad=20)
    ax.set_xlabel(dataframe.columns.name or "", fontsize=12)
    ax.set_ylabel(dataframe.index.name or "", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=8)
    plt.tight_layout()
    if output_path is not None:
        save_current_figure(output_path, dpi=dpi, close=close, figure=fig)
    return heatmap


if __name__ == "__main__":
    print("[INFO] anndata_expression_viz_helper_20260423_v1.py loaded successfully")
