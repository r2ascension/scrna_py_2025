#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generic lineage merge + scHPL comparison pipeline.

Design goal
-----------
Port the logic of `bcell_merge_schpl_20260330_v3_0.ipynb` into a reusable
script that works for lineages where the most recent artifacts are stored as a
single merged h5ad containing both reference and query cells (`data_source`).

Unified methodology
-------------------
This file implements the **global lineage-wise** variant of the shared
scHPL/treeArches workflow used in this directory.

- scHPL is treated as a **secondary hierarchical validator** after the base
    mapping model (typically scANVI / scArches), not as a replacement for it.
- Training and prediction are run on a shared integrated latent space
    (`cfg.latent_key`), while UMAP is visualization-only.
- Query cells are written back with the standard `schpl_*` columns:
    `schpl_pred_raw`, `schpl_pred`, `schpl_prob`, `schpl_rejected`, and
    `schpl_reject_type`.
- `Rejected` means the cell is not stably absorbed by the current reference
    hierarchy and should enter follow-up analysis; it is **not** an automatic
    claim of a new cell type.

Expected input
--------------
A merged h5ad that already contains:
- reference cells + query cells in one object
- `data_source` (or configured source key)
- one latent representation shared by ref/query (default: `X_scANVI`)
- a common UMAP basis (default: `X_umap`)
- reference truth labels on ref cells (configured key)
- query predictions + confidence on query cells (configured keys)

Outputs
-------
- merged reference + query h5ad with new `schpl_*` columns on query cells
- trained scHPL tree pickle
- config JSON
- multiple QC figures and CSV summaries

For the branch-wise stromal variant, see the dedicated
`stromal_schpl_treeArches_*.ipynb` notebooks and their follow-up script.
"""

from __future__ import annotations

import gc
import json
import pickle
import re
import time
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pandas.api.types import CategoricalDtype
from scipy.sparse import issparse, csr_matrix

from scHPL import predict as schpl_predict
from scHPL import train as schpl_train
from scHPL import utils as schpl_utils

warnings.filterwarnings("ignore")


@dataclass
class LineageMergeSchPLConfig:
    lineage_name: str
    lineage_slug: str
    version: str
    input_merged_h5ad: str
    output_dir: str
    reference_label_key: str
    query_pred_key: str
    query_final_key: str
    query_conf_key: str
    marker_genes: list[str] = field(default_factory=list)
    source_key: str = "data_source"
    reference_source_value: str = "reference"
    query_source_value: str = "query"
    batch_key: str = "sample"
    tissue_key: str = "tissue"
    latent_key: str = "X_scANVI"
    umap_key: str = "X_umap"
    marker_obsm_key: str = "marker_expr"
    marker_layer_preferred: str = "log1p"
    unlabeled: str = "Unknown"
    schpl_classifier: str = "knn"
    schpl_dimred: bool = True
    schpl_use_re: bool = True
    schpl_n_neighbors: int = 50
    schpl_fn: float = 0.5
    schpl_rej_threshold: float = 0.5
    schpl_qc_n_neighbors: int = 30
    schpl_qc_leiden_resolution: float = 0.5
    schpl_novel_rej_threshold: float = 0.50
    schpl_novel_min_cells: int = 50
    random_seed: int = 42
    dpi: int = 300
    figure_format: str = "pdf"
    source_palette: dict[str, str] = field(
        default_factory=lambda: {"reference": "#1f77b4", "query": "#ff7f0e"}
    )
    confidence_cmap: str = "viridis"
    expression_cmap: str = "Reds"


def build_global_categories(
    adata: sc.AnnData,
    ref_label_key: str,
    qry_label_key: str,
    source_key: str,
    ref_value: str,
    qry_value: str,
) -> list[str]:
    m_ref = adata.obs[source_key].astype(str).eq(ref_value)
    m_qry = adata.obs[source_key].astype(str).eq(qry_value)
    ref = adata.obs.loc[m_ref, ref_label_key].astype("string").dropna()
    qry = adata.obs.loc[m_qry, qry_label_key].astype("string").dropna()
    ref_order = pd.unique(ref)
    qry_extra = [x for x in pd.unique(qry) if x not in set(ref_order)]
    return list(ref_order) + sorted(qry_extra)


def make_palette(categories: Sequence[str]) -> dict[str, object]:
    n = len(categories)
    if n <= 20:
        colors = list(plt.get_cmap("tab20").colors)[:n]
    elif n <= 102:
        colors = sc.pl.palettes.default_102[:n]
    else:
        colors = [plt.cm.hsv(i / max(n, 1)) for i in range(n)]
        print(f"[WARN] {n} categories (>102), using HSV colormap")
    return dict(zip(categories, colors))


def attach_marker_layer(
    adata: sc.AnnData,
    genes: Sequence[str],
    layer_out: str = "marker_expr",
    prefer_layer: str = "log1p",
    allow_fallback: bool = True,
) -> tuple[sc.AnnData, list[str]]:
    if not genes:
        print("  [INFO] No marker genes configured; marker panel disabled")
        adata.uns[f"{layer_out}_genes"] = []
        adata.uns[f"{layer_out}_source"] = None
        adata.uns[f"{layer_out}_disabled_reason"] = "no_marker_genes_configured"
        return adata, []

    genes_avail = [g for g in genes if g in adata.var_names]
    if not genes_avail:
        print(
            f"  [WARN] None of the configured marker genes were found in var_names; "
            f"marker panel disabled ({len(genes)} requested)"
        )
        adata.uns[f"{layer_out}_genes"] = []
        adata.uns[f"{layer_out}_source"] = None
        adata.uns[f"{layer_out}_disabled_reason"] = "marker_genes_missing_from_var_names"
        return adata, []

    if prefer_layer in (adata.layers or {}):
        X = adata[:, genes_avail].layers[prefer_layer]
        source = f"layer:{prefer_layer}"
    elif allow_fallback and adata.raw is not None and all(g in adata.raw.var_names for g in genes_avail):
        X = adata.raw[:, genes_avail].X
        source = "raw"
    elif allow_fallback:
        X = adata[:, genes_avail].X
        source = "X"
    else:
        print(
            f"  [WARN] Required comparable layer '{prefer_layer}' not found and fallback disabled; "
            "marker panel disabled"
        )
        adata.uns[f"{layer_out}_genes"] = []
        adata.uns[f"{layer_out}_source"] = None
        adata.uns[f"{layer_out}_disabled_reason"] = (
            f"missing_preferred_layer:{prefer_layer}"
        )
        return adata, []

    adata.obsm[layer_out] = X.toarray() if issparse(X) else np.asarray(X)
    adata.uns[f"{layer_out}_genes"] = genes_avail
    adata.uns[f"{layer_out}_source"] = source
    print(f"  Stored {len(genes_avail)} markers in .obsm['{layer_out}'] (from {source})")
    return adata, genes_avail


def marker_sources_are_comparable(ref_source: str | None, qry_source: str | None) -> bool:
    return (
        ref_source is not None
        and qry_source is not None
        and ref_source == qry_source
        and (str(ref_source).startswith("layer:") or ref_source == "X")
    )


def make_safe_newick_label_map(labels: Iterable[str]) -> tuple[dict[str, str], dict[str, str]]:
    labels = [str(x) for x in pd.unique(pd.Series(labels, dtype="string").dropna())]
    label_to_token: dict[str, str] = {}
    token_to_label: dict[str, str] = {}
    used = {"root", "root2"}

    for i, label in enumerate(labels):
        base = re.sub(r"[^0-9A-Za-z_]", "_", label).strip("_") or f"label_{i}"
        token = base
        suffix = 1
        while token in used:
            token = f"{base}_{suffix}"
            suffix += 1
        used.add(token)
        label_to_token[label] = token
        token_to_label[token] = label

    return label_to_token, token_to_label


def _normalize_schpl_token(value: object) -> str:
    return re.sub(r"\s+", "", str(value).strip()).lower()


_SCHPL_REJECTION_TOKENS = {
    "dist": {"rejection(dist)", "rejected(dist)", "reject(dist)"},
    "RE": {
        "rejection(re)",
        "rejected(re)",
        "reject(re)",
        "rejection(reconstructionerror)",
        "rejected(reconstructionerror)",
    },
    "prob": {
        "rejection(prob)",
        "rejected(prob)",
        "reject(prob)",
        "root",
        "root2",
    },
}


def parse_schpl_predictions(
    y_pred_raw: Sequence[object],
    known_labels: Iterable[str],
    token_to_label: dict[str, str] | None = None,
) -> dict[str, np.ndarray]:
    token_to_label = token_to_label or {}
    known_labels = {str(x) for x in known_labels}

    raw_tokens = np.asarray(y_pred_raw, dtype=str)
    raw_norm = np.array([str(x).strip() for x in raw_tokens], dtype=object)
    raw_normalized = np.array([_normalize_schpl_token(x) for x in raw_norm], dtype=object)
    raw_mapped = np.array([token_to_label.get(x, x) for x in raw_norm], dtype=object)

    mask_rej_dist = np.isin(raw_normalized, list(_SCHPL_REJECTION_TOKENS["dist"]))
    mask_rej_re = np.isin(raw_normalized, list(_SCHPL_REJECTION_TOKENS["RE"]))
    mask_rej_prob = np.isin(raw_normalized, list(_SCHPL_REJECTION_TOKENS["prob"]))

    mask_rejected = mask_rej_dist | mask_rej_re | mask_rej_prob
    mask_accepted = np.array([x in known_labels for x in raw_mapped], dtype=bool) & ~mask_rejected
    mask_unknown = ~(mask_rejected | mask_accepted)

    if mask_unknown.any():
        unknown_vals = sorted(pd.unique(raw_norm[mask_unknown]).tolist())
        raise AssertionError(
            "Unparsed scHPL outputs encountered. "
            f"Inspect unique raw tokens and extend the explicit mapping if needed: {unknown_vals[:10]}"
        )

    pred_clean = raw_mapped.astype(object).copy()
    pred_clean[mask_rejected] = "Rejected"
    rej_type = np.where(
        mask_rej_dist,
        "dist",
        np.where(mask_rej_re, "RE", np.where(mask_rej_prob, "prob", "accepted")),
    )

    return {
        "raw_tokens": raw_tokens,
        "raw_mapped": raw_mapped,
        "pred_clean": pred_clean,
        "mask_rej_dist": mask_rej_dist,
        "mask_rej_re": mask_rej_re,
        "mask_rej_prob": mask_rej_prob,
        "mask_rejected": mask_rejected,
        "mask_accepted": mask_accepted,
        "rej_type": rej_type,
    }


def extract_schpl_confidence(
    y_prob_raw: object,
    raw_tokens: Sequence[object],
    raw_mapped: Sequence[object],
    token_to_label: dict[str, str] | None = None,
) -> np.ndarray:
    n_cells = len(raw_tokens)
    token_to_label = token_to_label or {}

    if y_prob_raw is None:
        return np.full(n_cells, np.nan, dtype=np.float32)

    prob_columns: list[str] | None = None
    if isinstance(y_prob_raw, pd.DataFrame):
        prob_values = y_prob_raw.to_numpy(dtype=np.float32, copy=False)
        prob_columns = [token_to_label.get(str(col), str(col)) for col in y_prob_raw.columns]
    else:
        prob_values = np.asarray(y_prob_raw, dtype=np.float32)

    if prob_values.ndim == 1:
        if prob_values.shape[0] != n_cells:
            raise ValueError(
                f"1D scHPL probability output length mismatch: {prob_values.shape[0]} != {n_cells}"
            )
        return prob_values.astype(np.float32, copy=False)

    if prob_values.ndim != 2:
        raise ValueError(f"Unsupported scHPL probability output ndim={prob_values.ndim}")

    if prob_values.shape[0] != n_cells:
        raise ValueError(
            f"2D scHPL probability output row mismatch: {prob_values.shape[0]} != {n_cells}"
        )

    if prob_values.shape[1] == 0:
        warnings.warn("scHPL returned an empty probability matrix; storing NaN confidence")
        print("[WARN] scHPL returned an empty probability matrix; storing NaN confidence")
        return np.full(n_cells, np.nan, dtype=np.float32)

    if prob_values.shape[1] == 1:
        return prob_values[:, 0].astype(np.float32, copy=False)

    row_max = prob_values.max(axis=1).astype(np.float32, copy=False)
    if not prob_columns:
        warnings.warn(
            "scHPL returned a 2D probability matrix without column labels; using row-wise max as cell-level confidence"
        )
        print(
            "[WARN] scHPL returned a 2D probability matrix without column labels; "
            "using row-wise max as schpl_prob"
        )
        return row_max

    col_index = {str(label): idx for idx, label in enumerate(prob_columns)}
    aligned = np.full(n_cells, np.nan, dtype=np.float32)
    matched = 0
    for i, label in enumerate(np.asarray(raw_mapped, dtype=object)):
        idx = col_index.get(str(label))
        if idx is not None:
            aligned[i] = prob_values[i, idx]
            matched += 1

    if matched == 0:
        warnings.warn(
            "Could not align predicted labels to probability matrix columns; using row-wise max as cell-level confidence"
        )
        print(
            "[WARN] Could not align predicted labels to probability matrix columns; "
            "using row-wise max as schpl_prob"
        )
        return row_max

    if matched < n_cells:
        warnings.warn(
            f"Aligned scHPL probabilities for {matched}/{n_cells} cells; filling unmatched rows with row-wise max"
        )
        print(
            f"[WARN] Aligned scHPL probabilities for {matched}/{n_cells} cells; "
            "filling unmatched rows with row-wise max"
        )
        aligned[np.isnan(aligned)] = row_max[np.isnan(aligned)]

    return aligned


def remove_unused_categories(adata: sc.AnnData, col: str) -> None:
    if col in adata.obs.columns:
        s = adata.obs[col]
        if isinstance(s.dtype, CategoricalDtype):
            adata.obs[col] = s.cat.remove_unused_categories()


def save_rasterized_figure(fig: plt.Figure, path: Path, dpi: int = 300) -> None:
    for ax in fig.axes:
        for coll in ax.collections:
            try:
                coll.set_rasterized(True)
            except Exception:
                pass
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def sanitize_object_columns(df: pd.DataFrame, df_name: str) -> None:
    bool_fixed: list[str] = []
    str_fixed: list[str] = []
    for col in df.columns:
        s = df[col]
        dtype_str = str(s.dtype)
        if dtype_str.startswith("string"):
            df[col] = s.astype("string")
            str_fixed.append(col)
            continue
        if s.dtype != object:
            continue
        non_null = s.dropna()
        if len(non_null) == 0:
            df[col] = pd.Series(pd.array([pd.NA] * len(s), dtype="string"), index=s.index)
            str_fixed.append(col)
            continue
        py_types = set(non_null.map(lambda x: type(x).__name__))
        if py_types <= {"bool"}:
            df[col] = s.astype("boolean")
            bool_fixed.append(col)
        elif py_types <= {"str"}:
            df[col] = s.astype("string")
            str_fixed.append(col)
        else:
            df[col] = s.map(lambda x: pd.NA if pd.isna(x) else str(x)).astype("string")
            str_fixed.append(col)
    if bool_fixed:
        print(f"  [sanitize] {df_name} object(bool+NA) -> pandas BooleanDtype: {bool_fixed}")
    if str_fixed:
        print(f"  [sanitize] {df_name} string/object -> pandas StringDtype: {str_fixed}")


def make_scanpy_plot_safe_series(series: pd.Series) -> pd.Series:
    """
    Work around scanpy 1.11.x + pandas 1.5.x categorical plotting issues.

    Some categorical color columns hit `Categorical.map(..., na_action=...)`,
    which is only supported by newer pandas. Converting plotting-only copies to
    plain object/string values avoids that code path without changing the saved
    artifacts written earlier in the pipeline.
    """
    if isinstance(series.dtype, CategoricalDtype):
        s = series.astype("string")
        return s.where(s.notna(), np.nan).astype(object)
    dtype_str = str(series.dtype)
    if dtype_str.startswith("string"):
        return series.where(series.notna(), np.nan).astype(object)
    return series


def _check_keys(adata: sc.AnnData, required_obs: Sequence[str], required_obsm: Sequence[str], adata_name: str) -> None:
    missing_obs = [k for k in required_obs if k not in adata.obs.columns]
    missing_obsm = [k for k in required_obsm if k not in adata.obsm]
    if missing_obs or missing_obsm:
        raise KeyError(
            f"[{adata_name}] Missing keys. obs={missing_obs} obsm={missing_obsm} | "
            f"available obs sample={list(adata.obs.columns)[:20]} obsm={list(adata.obsm.keys())}"
        )


def load_reference_and_query_from_merged(
    cfg: LineageMergeSchPLConfig,
) -> tuple[sc.AnnData, sc.AnnData, sc.AnnData]:
    path = Path(cfg.input_merged_h5ad)
    if not path.exists():
        raise FileNotFoundError(f"Merged input not found: {path}")

    print("=" * 80)
    print(f"LOAD MERGED INPUT — {cfg.lineage_name}")
    print("=" * 80)
    print(f"Input  : {path}")
    adata_in = sc.read_h5ad(path)
    adata_in.var_names_make_unique()
    adata_in.obs_names_make_unique()

    if cfg.source_key not in adata_in.obs.columns:
        raise KeyError(f"Merged input missing source key '{cfg.source_key}'")
    if cfg.batch_key not in adata_in.obs.columns:
        placeholder_batch = "unknown_batch"
        print(
            f"[WARN] Missing batch key '{cfg.batch_key}' in merged input; "
            f"creating placeholder column='{placeholder_batch}'"
        )
        adata_in.obs[cfg.batch_key] = placeholder_batch
    _check_keys(
        adata_in,
        [cfg.source_key, cfg.reference_label_key, cfg.query_pred_key, cfg.query_conf_key],
        [cfg.latent_key, cfg.umap_key],
        "merged_input",
    )

    src = adata_in.obs[cfg.source_key].astype(str)
    ref_mask = src.eq(cfg.reference_source_value)
    qry_mask = src.eq(cfg.query_source_value)
    if int(ref_mask.sum()) == 0 or int(qry_mask.sum()) == 0:
        raise ValueError(
            f"Could not split merged input by {cfg.source_key}. "
            f"reference={int(ref_mask.sum())}, query={int(qry_mask.sum())}"
        )

    adata_ref = adata_in[ref_mask].copy()
    adata_query = adata_in[qry_mask].copy()

    if cfg.query_final_key not in adata_query.obs.columns and cfg.query_pred_key in adata_query.obs.columns:
        adata_query.obs[cfg.query_final_key] = adata_query.obs[cfg.query_pred_key].astype(object)
    if cfg.query_pred_key not in adata_query.obs.columns and cfg.query_final_key in adata_query.obs.columns:
        adata_query.obs[cfg.query_pred_key] = adata_query.obs[cfg.query_final_key].astype(object)

    _check_keys(adata_ref, [cfg.reference_label_key], [cfg.latent_key, cfg.umap_key], "reference")
    _check_keys(
        adata_query,
        [cfg.query_pred_key, cfg.query_final_key, cfg.query_conf_key],
        [cfg.latent_key, cfg.umap_key],
        "query",
    )

    adata_ref.obsm["X_umap"] = np.asarray(adata_ref.obsm[cfg.umap_key]).copy()
    adata_query.obsm["X_umap"] = np.asarray(adata_query.obsm[cfg.umap_key]).copy()

    print(f"Merged shape    : {adata_in.shape}")
    print(f"Reference cells : {adata_ref.n_obs:,}")
    print(f"Query cells     : {adata_query.n_obs:,}")
    print(f"Latent key      : {cfg.latent_key}")
    print(f"UMAP key        : {cfg.umap_key}")

    return adata_in, adata_ref, adata_query


def run_pipeline(cfg: LineageMergeSchPLConfig) -> None:
    np.random.seed(cfg.random_seed)
    sc.settings.set_figure_params(
        dpi=cfg.dpi,
        facecolor="white",
        format=cfg.figure_format,
        vector_friendly=False,
    )
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    pipeline_start = time.time()
    output_dir = Path(cfg.output_dir)
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"{cfg.lineage_name} Merge + Visualization + scHPL")
    print("=" * 80)
    print(f"scanpy : {sc.__version__}")
    print(f"pandas : {pd.__version__}")
    print(f"numpy  : {np.__version__}")
    print(f"Input  : {cfg.input_merged_h5ad}")
    print(f"Output : {cfg.output_dir}")

    adata_in, adata_ref, adata_query = load_reference_and_query_from_merged(cfg)

    query_pred_key = cfg.query_pred_key
    query_final_key = cfg.query_final_key
    if query_final_key == query_pred_key:
        query_final_key = f"{query_pred_key}_final"
        if query_final_key not in adata_query.obs.columns:
            adata_query.obs[query_final_key] = adata_query.obs[query_pred_key].astype(object)

    print("=" * 80)
    print("ATTACH MARKER EXPRESSION")
    print("=" * 80)
    adata_ref, markers_ref = attach_marker_layer(
        adata_ref,
        cfg.marker_genes,
        layer_out=cfg.marker_obsm_key,
        prefer_layer=cfg.marker_layer_preferred,
        allow_fallback=True,
    )
    adata_query, markers_qry = attach_marker_layer(
        adata_query,
        cfg.marker_genes,
        layer_out=cfg.marker_obsm_key,
        prefer_layer=cfg.marker_layer_preferred,
        allow_fallback=True,
    )
    marker_source_ref = adata_ref.uns.get(f"{cfg.marker_obsm_key}_source")
    marker_source_qry = adata_query.uns.get(f"{cfg.marker_obsm_key}_source")
    marker_expr_comparable = marker_sources_are_comparable(marker_source_ref, marker_source_qry)
    markers_common = sorted(set(markers_ref) & set(markers_qry))
    print(f"Common markers : {markers_common}")
    print(f"Comparable     : {marker_expr_comparable} (ref={marker_source_ref}, query={marker_source_qry})")

    print("=" * 80)
    print("scHPL TRAIN (reference)")
    print("=" * 80)
    X_ref = np.asarray(adata_ref.obsm[cfg.latent_key], dtype=np.float32)
    y_ref_original = adata_ref.obs[cfg.reference_label_key].astype(str).values
    known_mask = y_ref_original != cfg.unlabeled
    n_unknown_ref = int((~known_mask).sum())
    if n_unknown_ref > 0:
        print(f"Dropping {n_unknown_ref} Unknown cells from reference training")
        X_ref = X_ref[known_mask]
        y_ref_original = y_ref_original[known_mask]

    class_counts = pd.Series(y_ref_original).value_counts()
    small_classes = class_counts[class_counts < cfg.schpl_n_neighbors]
    if len(small_classes) > 0:
        print(f"[INFO] {len(small_classes)} classes < {cfg.schpl_n_neighbors}; dynamic_neighbors=True will adjust")

    label_to_token, token_to_label = make_safe_newick_label_map(y_ref_original)
    y_ref = np.array([label_to_token[x] for x in y_ref_original], dtype=object)
    unique_tokens = np.array([label_to_token[x] for x in np.unique(y_ref_original)], dtype=object)
    tree_newick = f"({','.join(unique_tokens)})root;"
    tree_init = schpl_utils.create_tree(tree_newick)
    t0 = time.time()
    tree_trained = schpl_train.train_tree(
        data=X_ref,
        labels=y_ref,
        tree=tree_init,
        classifier=cfg.schpl_classifier,
        dimred=cfg.schpl_dimred,
        useRE=cfg.schpl_use_re,
        FN=cfg.schpl_fn,
        n_neighbors=cfg.schpl_n_neighbors,
        dynamic_neighbors=True,
    )
    print(f"[OK] Training done: {time.time() - t0:.1f} s")
    tree_path = output_dir / f"schpl_tree_{cfg.lineage_slug}_{cfg.version}.pkl"
    with open(tree_path, "wb") as f:
        pickle.dump(tree_trained, f)
    print(f"[OK] Tree saved: {tree_path}")

    print("=" * 80)
    print("scHPL PREDICT (query)")
    print("=" * 80)
    X_qry = np.asarray(adata_query.obsm[cfg.latent_key], dtype=np.float32)
    t0 = time.time()
    y_pred_raw, y_prob_raw = schpl_predict.predict_labels(
        X_qry,
        tree=tree_trained,
        threshold=cfg.schpl_rej_threshold,
    )
    elapsed_pred = time.time() - t0
    print(f"[OK] Prediction done: {elapsed_pred:.1f} s")

    parsed = parse_schpl_predictions(
        y_pred_raw,
        known_labels=np.unique(y_ref_original),
        token_to_label=token_to_label,
    )
    y_pred_arr = np.asarray(parsed["raw_mapped"], dtype=object)
    y_pred_clean = np.asarray(parsed["pred_clean"], dtype=object)
    mask_rej_dist = parsed["mask_rej_dist"]
    mask_rej_re = parsed["mask_rej_re"]
    mask_rej_prob = parsed["mask_rej_prob"]
    mask_rejected = parsed["mask_rejected"]
    mask_accepted = parsed["mask_accepted"]
    rej_type = parsed["rej_type"]

    col_schpl_raw = "schpl_pred_raw"
    col_schpl_pred = "schpl_pred"
    col_schpl_prob = "schpl_prob"
    col_schpl_rejected = "schpl_rejected"
    col_schpl_rej_type = "schpl_reject_type"

    adata_query.obs[col_schpl_raw] = y_pred_arr
    adata_query.obs[col_schpl_pred] = y_pred_clean
    adata_query.obs[col_schpl_rejected] = mask_rejected
    adata_query.obs[col_schpl_rej_type] = rej_type
    y_prob_arr = extract_schpl_confidence(
        y_prob_raw,
        raw_tokens=parsed["raw_tokens"],
        raw_mapped=parsed["raw_mapped"],
        token_to_label=token_to_label,
    )
    adata_query.obs[col_schpl_prob] = y_prob_arr

    n_total = adata_query.n_obs
    n_rejected = int(mask_rejected.sum())
    n_accepted = int(mask_accepted.sum())
    rej_rate = n_rejected / max(n_total, 1) * 100
    print(f"Accepted : {n_accepted:,} ({n_accepted/max(n_total,1)*100:.1f}%)")
    print(f"Rejected : {n_rejected:,} ({rej_rate:.1f}%)")

    print("=" * 80)
    print("NOVEL CELL TYPE CANDIDATE DETECTION")
    print("=" * 80)
    if "neighbors_schpl" not in adata_query.uns:
        sc.pp.neighbors(
            adata_query,
            use_rep=cfg.latent_key,
            n_neighbors=cfg.schpl_qc_n_neighbors,
            random_state=cfg.random_seed,
            key_added="neighbors_schpl",
        )
    sc.tl.leiden(
        adata_query,
        resolution=cfg.schpl_qc_leiden_resolution,
        key_added="leiden_schpl_qc",
        neighbors_key="neighbors_schpl",
        random_state=cfg.random_seed,
    )
    cluster_rej = (
        adata_query.obs.groupby("leiden_schpl_qc")[col_schpl_rejected]
        .agg(n_rejected="sum", n_total="count", rej_rate="mean")
        .sort_values("rej_rate", ascending=False)
    )
    cluster_rej["rej_rate_pct"] = (cluster_rej["rej_rate"] * 100).round(1)
    novel_clusters = cluster_rej[
        (cluster_rej["rej_rate"] > cfg.schpl_novel_rej_threshold)
        & (cluster_rej["n_rejected"] >= cfg.schpl_novel_min_cells)
    ].index.tolist()
    adata_query.obs["schpl_novel_candidate"] = (
        adata_query.obs["leiden_schpl_qc"].isin(novel_clusters)
        & adata_query.obs[col_schpl_rejected]
    )
    cluster_csv = output_dir / f"cluster_rejection_summary_{cfg.version}.csv"
    cluster_rej.to_csv(cluster_csv)
    print(f"Novel clusters : {novel_clusters}")
    print(f"[OK] {cluster_csv}")

    print("=" * 80)
    print("COMPARE BASE MODEL vs scHPL")
    print("=" * 80)
    df_accepted = adata_query.obs.loc[
        ~adata_query.obs[col_schpl_rejected],
        [query_pred_key, query_final_key, cfg.query_conf_key, col_schpl_pred],
    ].copy()
    agree_pred = df_accepted[query_pred_key].astype(str) == df_accepted[col_schpl_pred].astype(str)
    agree_final = df_accepted[query_final_key].astype(str) == df_accepted[col_schpl_pred].astype(str)
    agree_rate_pred = float(agree_pred.mean() * 100) if len(df_accepted) else float("nan")
    agree_rate_final = float(agree_final.mean() * 100) if len(df_accepted) else float("nan")
    print(f"Agreement (pred  vs scHPL): {agree_rate_pred:.1f}%")
    print(f"Agreement (final vs scHPL): {agree_rate_final:.1f}%")

    print("=" * 80)
    print("GENE ALIGNMENT + MERGE")
    print("=" * 80)
    ref_genes = list(adata_ref.var_names)
    common_genes = [g for g in ref_genes if g in adata_query.var_names]
    overlap_pct = len(common_genes) / max(len(ref_genes), 1) * 100
    print(f"Reference genes : {len(ref_genes):,}")
    print(f"Common genes    : {len(common_genes):,} ({overlap_pct:.1f}%)")
    adata_ref_m = adata_ref[:, common_genes].copy()
    adata_qry_m = adata_query[:, common_genes].copy()

    keep_obsm = {"X_umap", cfg.latent_key}
    if marker_expr_comparable:
        keep_obsm.add(cfg.marker_obsm_key)
    for adata in [adata_ref_m, adata_qry_m]:
        for key in list(adata.obsm.keys()):
            if key not in keep_obsm:
                del adata.obsm[key]
        adata.uns = {}
        if issparse(adata.X):
            adata.X = csr_matrix(adata.X)

    adata_ref_m.obs_names = [f"ref::{x}" for x in adata_ref_m.obs_names]
    adata_qry_m.obs_names = [f"qry::{x}" for x in adata_qry_m.obs_names]

    schpl_new_cols = [
        col_schpl_raw,
        col_schpl_pred,
        col_schpl_rej_type,
        "leiden_schpl_qc",
        "schpl_novel_candidate",
    ]
    for col in [query_pred_key, query_final_key] + schpl_new_cols:
        if col not in adata_ref_m.obs.columns:
            adata_ref_m.obs[col] = pd.NA
    for col in [cfg.query_conf_key, col_schpl_prob]:
        if col not in adata_ref_m.obs.columns:
            adata_ref_m.obs[col] = np.nan
    for col in [col_schpl_rejected, "schpl_novel_candidate"]:
        if col not in adata_ref_m.obs.columns:
            adata_ref_m.obs[col] = False
    if cfg.reference_label_key not in adata_qry_m.obs.columns:
        adata_qry_m.obs[cfg.reference_label_key] = pd.NA

    adata_merged = sc.concat(
        {cfg.reference_source_value: adata_ref_m, cfg.query_source_value: adata_qry_m},
        axis=0,
        join="inner",
        merge="unique",
        label=cfg.source_key,
    )
    expected_obsm = ["X_umap", cfg.latent_key] + ([cfg.marker_obsm_key] if marker_expr_comparable else [])
    for key in expected_obsm:
        assert key in adata_merged.obsm, f"{key} lost during merge"
        assert adata_merged.obsm[key].shape[0] == adata_merged.n_obs
    if marker_expr_comparable:
        adata_merged.uns[f"{cfg.marker_obsm_key}_genes"] = markers_common
    else:
        adata_merged.uns["marker_expr_disabled_reason"] = (
            f"Merged marker panels disabled because ref source={marker_source_ref}, query source={marker_source_qry}"
        )

    print("=" * 80)
    print("GLOBAL CATEGORIES + VIZ COLUMNS")
    print("=" * 80)
    global_cats = build_global_categories(
        adata_merged,
        cfg.reference_label_key,
        query_final_key,
        cfg.source_key,
        cfg.reference_source_value,
        cfg.query_source_value,
    )
    ct_palette = make_palette(global_cats)
    global_cats_ext = global_cats + (["Rejected"] if "Rejected" not in global_cats else [])
    ct_palette_ext = make_palette(global_cats_ext)
    ct_palette_ext["Rejected"] = "#BDBDBD"

    ref_mask = adata_merged.obs[cfg.source_key].astype(str).eq(cfg.reference_source_value)
    qry_mask = adata_merged.obs[cfg.source_key].astype(str).eq(cfg.query_source_value)

    def _make_label_col(src_col: str, mask: pd.Series, cats: Sequence[str]) -> pd.Categorical:
        s = pd.Series(pd.NA, index=adata_merged.obs_names, dtype="string")
        if src_col in adata_merged.obs.columns:
            s.loc[mask] = adata_merged.obs.loc[mask, src_col].astype("string")
        return pd.Categorical(s, categories=cats)

    viz_ref_only = "viz__ref_label_only"
    viz_qry_final = "viz__qry_label_final_only"
    viz_qry_pred = "viz__qry_label_pred_only"
    viz_qry_conf = "viz__qry_conf_only"
    viz_schpl_pred = "viz__schpl_pred_qry_only"

    adata_merged.obs[viz_qry_final] = _make_label_col(query_final_key, qry_mask, global_cats)
    adata_merged.obs[viz_qry_pred] = _make_label_col(query_pred_key, qry_mask, global_cats)
    adata_merged.obs[viz_ref_only] = _make_label_col(cfg.reference_label_key, ref_mask, global_cats)
    conf_arr = np.full(adata_merged.n_obs, np.nan)
    conf_arr[qry_mask.values] = pd.to_numeric(
        adata_merged.obs.loc[qry_mask, cfg.query_conf_key], errors="coerce"
    ).values
    adata_merged.obs[viz_qry_conf] = conf_arr
    adata_merged.obs[viz_schpl_pred] = _make_label_col(col_schpl_pred, qry_mask, global_cats_ext)

    print("=" * 80)
    print("SAVE MERGED DATA + CONFIG")
    print("=" * 80)
    sanitize_object_columns(adata_merged.obs, "obs")
    sanitize_object_columns(adata_merged.var, "var")
    merged_output = output_dir / f"{cfg.lineage_slug}_reference_plus_query_schpl_{cfg.version}.h5ad"
    try:
        adata_merged.write_h5ad(merged_output, compression="gzip")
    except RuntimeError as e:
        err_msg = str(e)
        if "allow_write_nullable_strings" not in err_msg:
            raise
        if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
            print(
                "[WARN] anndata refused nullable string write; enabling "
                "anndata.settings.allow_write_nullable_strings=True and retrying"
            )
            ad.settings.allow_write_nullable_strings = True
            adata_merged.write_h5ad(merged_output, compression="gzip")
        else:
            raise
    print(f"[OK] {merged_output}")

    viz_config = {
        **asdict(cfg),
        "timestamp": datetime.now().isoformat(),
        "n_reference_cells": int(ref_mask.sum()),
        "n_query_cells": int(qry_mask.sum()),
        "query_pred_key_resolved": query_pred_key,
        "query_final_key_resolved": query_final_key,
        "marker_expr": {
            "preferred_layer": cfg.marker_layer_preferred,
            "reference_source": marker_source_ref,
            "query_source": marker_source_qry,
            "comparable_for_merged_plots": bool(marker_expr_comparable),
        },
        "schpl": {
            "n_ref_cells": int(X_ref.shape[0]),
            "n_ref_classes": int(len(np.unique(y_ref_original))),
            "n_rejected": n_rejected,
            "rejection_rate": round(rej_rate, 3),
            "qc_n_neighbors": int(cfg.schpl_qc_n_neighbors),
            "qc_leiden_resolution": float(cfg.schpl_qc_leiden_resolution),
            "novel_rej_threshold": float(cfg.schpl_novel_rej_threshold),
            "novel_min_cells": int(cfg.schpl_novel_min_cells),
            "novel_clusters": novel_clusters,
            "agreement_pred": round(float(agree_rate_pred), 3),
            "agreement_final": round(float(agree_rate_final), 3),
        },
        "global_categories": global_cats,
        "marker_genes_found": markers_common,
    }
    cfg_path = output_dir / f"config_{cfg.version}.json"
    with open(cfg_path, "w") as f:
        json.dump(viz_config, f, indent=2)
    print(f"[OK] {cfg_path}")

    print("=" * 80)
    print("FIGURES")
    print("=" * 80)
    for plot_col in [
        cfg.source_key,
        viz_ref_only,
        viz_qry_final,
        viz_qry_pred,
        viz_schpl_pred,
        col_schpl_rej_type,
    ]:
        if plot_col in adata_merged.obs.columns:
            adata_merged.obs[plot_col] = make_scanpy_plot_safe_series(adata_merged.obs[plot_col])

    marker_genes_list = adata_merged.uns.get(f"{cfg.marker_obsm_key}_genes", []) if marker_expr_comparable else []
    umap_coords = adata_merged.obsm["X_umap"]

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    sc.pl.umap(
        adata_merged,
        color=cfg.source_key,
        ax=axes[0, 0],
        show=False,
        title="Data Source",
        palette=cfg.source_palette,
        frameon=False,
        s=20,
    )
    sc.pl.umap(
        adata_merged,
        color=viz_ref_only,
        ax=axes[0, 1],
        show=False,
        title="Reference Labels",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    sc.pl.umap(
        adata_merged,
        color=viz_qry_final,
        ax=axes[0, 2],
        show=False,
        title="Base Model Final Label (Query)",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    sc.pl.umap(
        adata_merged,
        color=viz_schpl_pred,
        ax=axes[1, 0],
        show=False,
        title="scHPL Prediction (Rejected=grey)",
        legend_loc="right margin",
        palette=ct_palette_ext,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    sc.pl.umap(
        adata_merged,
        color=viz_qry_conf,
        ax=axes[1, 1],
        show=False,
        title="Base Model Confidence (Query)",
        cmap=cfg.confidence_cmap,
        vmin=0,
        vmax=1,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    if marker_expr_comparable and marker_genes_list:
        marker_gene = marker_genes_list[0]
        idx = marker_genes_list.index(marker_gene)
        expr = adata_merged.obsm[cfg.marker_obsm_key][:, idx]
        sc_obj = axes[1, 2].scatter(
            umap_coords[:, 0],
            umap_coords[:, 1],
            c=expr,
            cmap=cfg.expression_cmap,
            s=20,
            rasterized=True,
        )
        axes[1, 2].set_title(f"{marker_gene} Expression", fontsize=14)
        axes[1, 2].axis("off")
        plt.colorbar(sc_obj, ax=axes[1, 2], fraction=0.046, pad=0.04)
    else:
        axes[1, 2].text(
            0.5,
            0.5,
            "Marker panel disabled",
            ha="center",
            va="center",
            transform=axes[1, 2].transAxes,
        )
        axes[1, 2].axis("off")
    plt.suptitle(f"{cfg.lineage_name} Merge + scHPL ({cfg.version})", fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_rasterized_figure(fig, figures_dir / f"merged_overview_6panel_{cfg.version}.{cfg.figure_format}", dpi=cfg.dpi)

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    sc.pl.umap(adata_merged, color=viz_qry_final, ax=axes[0, 0], show=False,
               title="Base Model Final Label", legend_loc="right margin",
               palette=ct_palette, frameon=False, s=20, na_color="lightgray")
    sc.pl.umap(adata_merged, color=viz_schpl_pred, ax=axes[0, 1], show=False,
               title="scHPL Prediction", legend_loc="right margin",
               palette=ct_palette_ext, frameon=False, s=20, na_color="lightgray")
    sc.pl.umap(adata_merged, color=col_schpl_rejected, ax=axes[0, 2], show=False,
               title="scHPL Rejected Cells", frameon=False, s=20)
    sc.pl.umap(adata_merged, color=col_schpl_rej_type, ax=axes[1, 0], show=False,
               title="Rejection Type", legend_loc="right margin",
               legend_fontsize=8, frameon=False, s=20)
    sc.pl.umap(adata_merged, color=col_schpl_prob, ax=axes[1, 1], show=False,
               title="scHPL Posterior Probability", cmap=cfg.confidence_cmap,
               vmin=0, vmax=1, frameon=False, s=20)
    sc.pl.umap(adata_merged, color="schpl_novel_candidate", ax=axes[1, 2], show=False,
               title="Novel Candidate Cells", frameon=False, s=20)
    plt.suptitle(f"scHPL Overview — {cfg.lineage_name}", fontsize=13, fontweight="bold")
    plt.tight_layout()
    save_rasterized_figure(fig, figures_dir / f"schpl_overview_6panel_{cfg.version}.{cfg.figure_format}", dpi=cfg.dpi)

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    qry_m = adata_merged.obs[cfg.source_key].astype(str).eq(cfg.query_source_value)
    prob_all = adata_merged.obs.loc[qry_m, col_schpl_prob]
    rej_flag = adata_merged.obs.loc[qry_m, col_schpl_rejected]
    prob_acc = prob_all[~rej_flag].dropna()
    prob_rej = prob_all[rej_flag].dropna()
    axes[0].hist(prob_acc, bins=50, alpha=0.7, color="steelblue", label="Accepted", density=True)
    axes[0].hist(prob_rej, bins=50, alpha=0.7, color="tomato", label="Rejected", density=True)
    axes[0].axvline(cfg.schpl_rej_threshold, color="red", linestyle="--",
                    label=f"rej_threshold={cfg.schpl_rej_threshold}")
    axes[0].legend(); axes[0].grid(alpha=0.3); axes[0].set_title("scHPL Probability")
    rej_type_counts = adata_merged.obs.loc[qry_m, col_schpl_rej_type].value_counts()
    colors_pie = ["steelblue" if t == "accepted" else "tomato" for t in rej_type_counts.index]
    axes[1].pie(rej_type_counts.values, labels=rej_type_counts.index,
                autopct="%1.1f%%", colors=colors_pie, startangle=90)
    axes[1].set_title("Rejection Type Distribution")
    conf_acc = adata_merged.obs.loc[qry_m & ~adata_merged.obs[col_schpl_rejected], cfg.query_conf_key]
    conf_rej = adata_merged.obs.loc[qry_m & adata_merged.obs[col_schpl_rejected], cfg.query_conf_key]
    axes[2].hist(pd.to_numeric(conf_acc, errors="coerce").dropna(), bins=50,
                 alpha=0.7, color="steelblue", label="scHPL Accepted", density=True)
    axes[2].hist(pd.to_numeric(conf_rej, errors="coerce").dropna(), bins=50,
                 alpha=0.7, color="tomato", label="scHPL Rejected", density=True)
    axes[2].axvline(0.5, color="black", linestyle="--", label="base-model thr=0.5")
    axes[2].legend(); axes[2].grid(alpha=0.3); axes[2].set_title("Base Model Confidence")
    plt.tight_layout()
    save_rasterized_figure(fig, figures_dir / f"schpl_rejection_analysis_{cfg.version}.{cfg.figure_format}", dpi=cfg.dpi)

    print("=" * 80)
    print("QC SUMMARY TABLE")
    print("=" * 80)
    df_qc = pd.DataFrame(
        {
            "cell_type": adata_merged.obs.loc[qry_m, viz_qry_final].astype("string"),
            "confidence": pd.to_numeric(adata_merged.obs.loc[qry_m, viz_qry_conf], errors="coerce"),
            "schpl_pred": adata_merged.obs.loc[qry_m, col_schpl_pred].astype(str),
            "schpl_rej": adata_merged.obs.loc[qry_m, col_schpl_rejected],
        }
    ).dropna(subset=["cell_type", "confidence"])
    summary = df_qc.groupby("cell_type")["confidence"].agg(
        Count="count", Mean="mean", Median="median", Std="std", Min="min", Max="max"
    ).round(3)
    summary["Pct"] = (summary["Count"] / max(summary["Count"].sum(), 1) * 100).round(2)
    summary["Rej_n"] = df_qc.groupby("cell_type")["schpl_rej"].sum().astype(int)
    summary["Rej_rate"] = (summary["Rej_n"] / summary["Count"] * 100).round(1)
    summary = summary.sort_values("Count", ascending=False)
    qc_csv = output_dir / f"qc_summary_by_celltype_{cfg.version}.csv"
    summary.to_csv(qc_csv)
    print(summary.to_string())
    print(f"[OK] {qc_csv}")

    top_cts = summary.head(10).index
    df_top = df_qc[df_qc["cell_type"].isin(top_cts)]
    if len(df_top) > 0:
        ct_order = df_top.groupby("cell_type")["confidence"].median().sort_values(ascending=False).index
        fig, axes = plt.subplots(1, 2, figsize=(18, 6))
        sns.violinplot(data=df_top, x="cell_type", y="confidence", order=ct_order, ax=axes[0], palette="Set2")
        axes[0].axhline(0.5, color="red", linestyle="--", linewidth=2)
        axes[0].tick_params(axis="x", rotation=45)
        axes[0].grid(axis="y", alpha=0.3)
        axes[0].set_title("Base Model Confidence by Cell Type (Violin)")
        sns.boxplot(data=df_top, x="cell_type", y="confidence", order=ct_order, ax=axes[1], palette="Set2")
        axes[1].axhline(0.5, color="red", linestyle="--", linewidth=2)
        axes[1].tick_params(axis="x", rotation=45)
        axes[1].grid(axis="y", alpha=0.3)
        axes[1].set_title("Base Model Confidence by Cell Type (Box)")
        plt.tight_layout()
        save_rasterized_figure(fig, figures_dir / f"confidence_by_celltype_{cfg.version}.{cfg.figure_format}", dpi=cfg.dpi)

    elapsed = (time.time() - pipeline_start) / 60
    label_csv = output_dir / f"label_comparison_{cfg.version}.csv"
    label_export_cols = [
        cfg.source_key,
        cfg.batch_key,
        cfg.reference_label_key,
        query_pred_key,
        query_final_key,
        cfg.query_conf_key,
        col_schpl_raw,
        col_schpl_pred,
        col_schpl_prob,
        col_schpl_rejected,
        col_schpl_rej_type,
        "leiden_schpl_qc",
        "schpl_novel_candidate",
    ]
    label_export_cols = [col for col in label_export_cols if col in adata_merged.obs.columns]
    adata_merged.obs[label_export_cols].to_csv(label_csv)
    print(f"[OK] {label_csv}")

    print("\n" + "=" * 80)
    print(f"FINAL SUMMARY — {cfg.lineage_name} Merge + scHPL ({cfg.version})")
    print("=" * 80)
    print(f"Runtime      : {elapsed:.1f} min")
    print(f"Merged shape : {adata_merged.shape}")
    print(f"Reference    : {int(ref_mask.sum()):,}")
    print(f"Query        : {int(qry_mask.sum()):,}")
    print(f"Accepted     : {n_accepted:,} ({n_accepted / max(n_total, 1) * 100:.1f}%)")
    print(f"Rejected     : {n_rejected:,} ({rej_rate:.1f}%)")
    print(f"Novel cand   : {int(adata_merged.obs['schpl_novel_candidate'].sum()):,}")
    print(f"Outputs      : {merged_output}")
    print(f"Tree         : {tree_path}")
    print(f"Config       : {cfg_path}")
    print(f"Figures dir  : {figures_dir}/ ({len(list(figures_dir.glob('*')))} files)")
    print("=" * 80)

    del adata_in, adata_ref, adata_query, adata_ref_m, adata_qry_m, adata_merged
    gc.collect()
