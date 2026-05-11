#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
B-cell full rerun for query UMAP same-space pairing (2026-04-19).

This script rebuilds the B-cell ref→query→merge chain in a fresh 0419 output tree.
It intentionally reuses the validated 2026-04-17 reference retraining logic, then
maps query cells with the saved UMAP operator so query coordinates are guaranteed to
live in the same space as the fresh reference embedding.

Pipeline stages
---------------
1. Reference rerun (reuse `bcell_scvi_scanvi_ref_c13c25drop_20260417.py` logic)
2. Query scArches mapping into the fresh reference model + UMAP operator transform
3. Full-gene merge with HVG `.X`, `layers['counts']`, `layers['log1p']`, and `.raw`
4. Re-run B-cell schpl wrapper on the fresh merged object
5. Re-run schpl reject follow-up on the fresh schpl output

Design notes
------------
- Outputs are written to new 0419 directories for safe comparison.
- The reference stage creates compatibility aliases (`reference_with_L2_umap.h5ad`,
  `scanvi_existing_model`, `umap_operator.joblib`) so downstream tooling keeps the
  old production contract while using the new reference.
- Query predictions come from the fresh L3 scANVI model, then collapse to L2 for the
  downstream B-cell schpl workflow.
"""

from __future__ import annotations

import argparse
import gc
import importlib.util
import json
import os
import pickle
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse, vstack as sp_vstack

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path("/home/h2048")
SCRIPT_DIR = PROJECT_ROOT / "script" / "py"
REF_MODULE_PATH = SCRIPT_DIR / "bcell_scvi_scanvi_ref_c13c25drop_20260417.py"
SCHPL_SCRIPT = SCRIPT_DIR / "bcell_merge_schpl_20260412_v1.py"
FOLLOWUP_SCRIPT = SCRIPT_DIR / "bcell_schpl_reject_followup_20260412_v1.py"
DEFAULT_SCHPL_PYTHON = PROJECT_ROOT / "miniconda3" / "envs" / "scarches_stable" / "bin" / "python"
SCHPL_PYTHON_ENVVAR = "B_CELL_SCHPL_PYTHON"

OUTPUT_ROOT = PROJECT_ROOT / "data" / "py" / "0419" / "bcell_query_umap_full_rerun_20260419"
REFERENCE_DIR = OUTPUT_ROOT / "reference"
QUERY_DIR = OUTPUT_ROOT / "query_mapping"
MERGE_DIR = OUTPUT_ROOT / "merged"
MERGE_FIG_DIR = MERGE_DIR / "figures"
SUMMARY_JSON = OUTPUT_ROOT / "pipeline_summary.json"

QUERY_INPUT_H5AD = PROJECT_ROOT / "data" / "py" / "0127" / "scarches_mapping_FIXED_v1_2" / "subsets" / "b_cells.h5ad"
REFERENCE_ALIAS_H5AD = REFERENCE_DIR / "reference_with_L2_umap.h5ad"
REFERENCE_UMAP_ALIAS = REFERENCE_DIR / "umap_operator.joblib"
REFERENCE_SCANVI_ALIAS = REFERENCE_DIR / "scanvi_existing_model"
QUERY_OUTPUT_H5AD = QUERY_DIR / "query_mapped_L2_20260419.h5ad"
MERGED_OUTPUT_H5AD = MERGE_DIR / "reference_plus_query_merged_fullgenes_L2_20260419.h5ad"
MERGED_OVERVIEW_FIG = MERGE_FIG_DIR / "merged_umap_overview_20260419.pdf"
MERGED_SUMMARY_JSON = MERGE_DIR / "merge_summary_20260419.json"

SCHPL_VERSION = "v1_3_scanvi_rerun_20260419"
SCHPL_OUTPUT_DIR = PROJECT_ROOT / "data" / "py" / "0419" / f"bcell_merge_schpl_{SCHPL_VERSION}"
SCHPL_OUTPUT_H5AD = SCHPL_OUTPUT_DIR / f"bcell_reference_plus_query_schpl_{SCHPL_VERSION}.h5ad"
SCHPL_CONFIG_JSON = SCHPL_OUTPUT_DIR / f"config_{SCHPL_VERSION}.json"
FOLLOWUP_OUTPUT_DIR = PROJECT_ROOT / "data" / "py" / "0419" / f"bcell_schpl_reject_followup_{SCHPL_VERSION}"

QUERY_MIN_CELLS = 3
QUERY_CONFIDENCE_THRESHOLD = 0.5
QUERY_MAX_EPOCHS = 200
QUERY_BATCH_SIZE = 256
QUERY_LEARNING_RATE = 5e-4
QUERY_WEIGHT_DECAY = 0.0
RANDOM_SEED = 42
DPI = 300
FIGURE_FORMAT = "pdf"

sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor="white", format=FIGURE_FORMAT)


@dataclass
class ReferenceArtifacts:
    output_dir: Path
    h5ad: Path
    compat_h5ad: Path
    compat_umap: Path
    compat_model_dir: Path
    training_config: Path
    hvg_file: Path


@dataclass
class QueryArtifacts:
    output_h5ad: Path
    n_cells: int
    n_l2_labels: int


@dataclass
class MergeArtifacts:
    output_h5ad: Path
    n_cells: int
    n_reference: int
    n_query: int
    n_hvgs: int
    n_full_genes: int


def print_banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def remove_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def resolve_schpl_python_executable() -> Path:
    override = os.environ.get(SCHPL_PYTHON_ENVVAR, "").strip()
    candidate = Path(override).expanduser() if override else DEFAULT_SCHPL_PYTHON
    candidate = candidate.resolve()
    if not candidate.exists():
        raise FileNotFoundError(
            "schpl/follow-up python executable not found: "
            f"{candidate}\n"
            f"Set {SCHPL_PYTHON_ENVVAR} to a valid interpreter if needed."
        )
    return candidate


def link_or_copy(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        remove_path(dst)
    try:
        os.symlink(src, dst)
    except OSError:
        if src.is_dir():
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)


def load_module(module_path: Path, module_name: str) -> Any:
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # type: ignore[union-attr]
    return module


def _to_csr(X: Any) -> csr_matrix:
    if issparse(X):
        return csr_matrix(X)
    return csr_matrix(np.asarray(X))


def filter_genes_min_cells(X_csr: csr_matrix, genes: pd.Index, min_cells: int) -> tuple[np.ndarray, pd.Index, np.ndarray]:
    nnz = np.asarray(X_csr.getnnz(axis=0)).ravel()
    keep = nnz >= int(min_cells)
    kept_genes = genes[keep]
    return keep, kept_genes, nnz


def align_full_to_reference_genes(X_qry_csr: csr_matrix, qry_genes: pd.Index, ref_genes: pd.Index) -> tuple[csr_matrix, dict[str, float | int]]:
    ref_index = pd.Index(ref_genes)
    ref_pos = ref_index.get_indexer(qry_genes)
    keep_mask = ref_pos >= 0
    n_shared = int(keep_mask.sum())
    if n_shared == 0:
        raise ValueError("No shared genes between query and reference full gene list")

    X_shared = X_qry_csr[:, np.where(keep_mask)[0]]
    ref_pos_shared = ref_pos[keep_mask]

    order = np.argsort(ref_pos_shared)
    X_shared = X_shared[:, order]
    ref_pos_sorted = ref_pos_shared[order]

    X_coo = X_shared.tocoo()
    new_col = ref_pos_sorted[X_coo.col]
    X_aligned = csr_matrix(
        (X_coo.data, (X_coo.row, new_col)),
        shape=(X_qry_csr.shape[0], len(ref_genes)),
    )
    stats = {
        "n_qry_genes": int(len(qry_genes)),
        "n_ref_genes": int(len(ref_genes)),
        "n_shared": n_shared,
        "shared_frac_vs_ref": n_shared / max(1, len(ref_genes)),
        "shared_frac_vs_qry": n_shared / max(1, len(qry_genes)),
    }
    return X_aligned, stats


def compute_log1p_layer(counts: csr_matrix) -> csr_matrix:
    temp = AnnData(X=counts.copy())
    sc.pp.normalize_total(temp, target_sum=1e4)
    sc.pp.log1p(temp)
    return _to_csr(temp.X)


def _as_python_scalar(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    serializable = json.loads(json.dumps(payload, default=lambda x: _as_python_scalar(x)))
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(serializable, handle, indent=2, ensure_ascii=False)


def prepare_reference_aliases(ref_module: Any) -> ReferenceArtifacts:
    compat_h5ad = REFERENCE_ALIAS_H5AD
    compat_model_dir = REFERENCE_SCANVI_ALIAS
    compat_umap = REFERENCE_UMAP_ALIAS

    link_or_copy(ref_module.OUTPUT_H5AD, compat_h5ad)
    link_or_copy(ref_module.OUTPUT_SCANVI_DIR, compat_model_dir)
    link_or_copy(ref_module.OUTPUT_UMAP, compat_umap)

    ref_adata = sc.read_h5ad(ref_module.OUTPUT_H5AD, backed="r")
    try:
        tissue_series = ref_adata.obs[ref_module.TISSUE_KEY]
        batch_series = ref_adata.obs[ref_module.BATCH_KEY]
        if hasattr(tissue_series.dtype, "categories"):
            tissue_categories = list(map(str, tissue_series.dtype.categories))
        else:
            tissue_categories = sorted(map(str, pd.Index(tissue_series.astype(str).unique())))
        if hasattr(batch_series.dtype, "categories"):
            batch_categories = list(map(str, batch_series.dtype.categories))
        else:
            batch_categories = sorted(map(str, pd.Index(batch_series.astype(str).unique())))
    finally:
        if getattr(ref_adata, "file", None) is not None:
            ref_adata.file.close()

    config_path = ref_module.OUTPUT_CONFIG
    with open(config_path, "r", encoding="utf-8") as handle:
        config = json.load(handle)
    config.update(
        {
            "reference_with_umap_h5ad": str(compat_h5ad),
            "scanvi_existing_model_dir": str(compat_model_dir),
            "umap_operator_joblib": str(compat_umap),
            "tissue_categories": tissue_categories,
            "batch_categories": batch_categories,
            "query_mapping_confidence_threshold": QUERY_CONFIDENCE_THRESHOLD,
            "query_mapping_input_h5ad": str(QUERY_INPUT_H5AD),
            "compat_alias_created_20260419": True,
        }
    )
    write_json(config_path, config)

    return ReferenceArtifacts(
        output_dir=ref_module.OUTPUT_DIR,
        h5ad=ref_module.OUTPUT_H5AD,
        compat_h5ad=compat_h5ad,
        compat_umap=compat_umap,
        compat_model_dir=compat_model_dir,
        training_config=config_path,
        hvg_file=ref_module.OUTPUT_HVG,
    )


def run_reference_rerun() -> tuple[Any, ReferenceArtifacts]:
    print_banner("Stage 1/5: fresh reference rerun (reuse 2026-04-17 logic)")
    ensure_dir(REFERENCE_DIR)
    ref_module = load_module(REF_MODULE_PATH, "bcell_ref_rerun_20260417_dynamic")

    ref_module.OUTPUT_DIR = REFERENCE_DIR
    ref_module.OUTPUT_H5AD = REFERENCE_DIR / "bcell_reference_c22_c13_c25drop_scanvi_L3_ref_20260419.h5ad"
    ref_module.OUTPUT_SCVI_DIR = REFERENCE_DIR / "scvi_bcell_ref_c22_c13_c25drop_l3_20260419"
    ref_module.OUTPUT_SCANVI_DIR = REFERENCE_DIR / "scanvi_bcell_L3_ref_c22_c13_c25drop_20260419"
    ref_module.OUTPUT_UMAP = REFERENCE_DIR / "umap_operator_scanvi.joblib"
    ref_module.OUTPUT_REMOVED = REFERENCE_DIR / "removed_choir_clusters_c22_c13_c25_cells.csv"
    ref_module.OUTPUT_TARGETS = REFERENCE_DIR / "target_choir_cluster_union_c22_c13_c25.csv"
    ref_module.OUTPUT_CONFIG = REFERENCE_DIR / "training_config.json"
    ref_module.OUTPUT_HVG = REFERENCE_DIR / "hvg_genes.txt"

    ref_module.main()
    artifacts = prepare_reference_aliases(ref_module)
    print(f"[OK] compatibility aliases prepared in {artifacts.output_dir}")
    return ref_module, artifacts


def make_query_work_adata(query_original: AnnData, ref_module: Any, ref_cfg: dict[str, Any]) -> AnnData:
    adata = query_original.copy()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    ref_module.prepare_covariates(adata)

    tissue_categories = list(map(str, ref_cfg.get("tissue_categories", [])))
    if not tissue_categories:
        tissue_categories = ["unknown_tissue"]
    if "unknown_tissue" not in tissue_categories:
        tissue_categories.append("unknown_tissue")

    tissue_values = adata.obs[ref_module.TISSUE_KEY].astype(str).fillna("unknown_tissue")
    tissue_values = tissue_values.where(tissue_values.isin(set(tissue_categories)), "unknown_tissue")
    adata.obs[ref_module.TISSUE_KEY] = pd.Categorical(tissue_values, categories=tissue_categories)

    if ref_module.BATCH_KEY not in adata.obs.columns:
        adata.obs[ref_module.BATCH_KEY] = "query_batch"
    adata.obs[ref_module.BATCH_KEY] = adata.obs[ref_module.BATCH_KEY].astype(str)

    adata.obs[ref_module.SCANVI_LABELS_KEY] = pd.Categorical(
        [ref_module.UNLABELED_CATEGORY] * adata.n_obs,
        categories=[ref_module.UNLABELED_CATEGORY],
    )
    return adata


def restore_query_metadata(
    query_original: AnnData,
    obs_names: pd.Index,
    latent: np.ndarray,
    umap_coords: np.ndarray,
    l3_pred: pd.Series,
    l2_pred: pd.Series,
    confidence: np.ndarray,
    ref_module: Any,
) -> AnnData:
    adata = query_original[obs_names].copy()
    adata.obs_names_make_unique()

    if "X_umap" in adata.obsm and "X_umap_pre_20260419_rerun" not in adata.obsm:
        adata.obsm["X_umap_pre_20260419_rerun"] = np.asarray(adata.obsm["X_umap"]).copy()

    final_labels = np.where(confidence >= QUERY_CONFIDENCE_THRESHOLD, l2_pred.astype(str).to_numpy(), ref_module.UNLABELED_CATEGORY)

    adata.obs[ref_module.L3_PRED_KEY] = pd.Categorical(l3_pred.astype(str))
    adata.obs[ref_module.LEGACY_L3_OUTPUT_KEY] = pd.Categorical(l3_pred.astype(str))
    adata.obs["L2_scanvi_pred"] = pd.Categorical(l2_pred.astype(str))
    adata.obs["Cell_Type_L2_pred"] = pd.Categorical(l2_pred.astype(str))
    adata.obs["Cell_Type_L2_final"] = pd.Categorical(final_labels.astype(str))
    adata.obs[ref_module.L3_CONFIDENCE_KEY] = pd.to_numeric(confidence, errors="coerce")
    adata.obs["mapping_confidence"] = pd.to_numeric(confidence, errors="coerce")
    adata.obs["L2_scanvi_confidence"] = pd.to_numeric(confidence, errors="coerce")

    adata.obsm["X_scanvi"] = np.asarray(latent, dtype=np.float32)
    adata.obsm["X_scANVI_L2"] = np.asarray(latent, dtype=np.float32)
    adata.obsm["X_scANVI_L3"] = np.asarray(latent, dtype=np.float32)
    adata.obsm["X_umap_mapped"] = np.asarray(umap_coords, dtype=np.float32)
    adata.obsm["X_umap"] = np.asarray(umap_coords, dtype=np.float32)

    adata.uns["query_mapping_reference"] = {
        "reference_output_dir": str(REFERENCE_DIR),
        "reference_h5ad": str(REFERENCE_ALIAS_H5AD),
        "reference_scanvi_model": str(REFERENCE_SCANVI_ALIAS),
        "umap_operator": str(REFERENCE_UMAP_ALIAS),
        "confidence_threshold": QUERY_CONFIDENCE_THRESHOLD,
    }
    return adata


def run_query_mapping(ref_module: Any, ref_artifacts: ReferenceArtifacts) -> QueryArtifacts:
    print_banner("Stage 2/5: query mapping into fresh reference + same-space UMAP transform")
    ensure_dir(QUERY_DIR)

    with open(ref_artifacts.training_config, "r", encoding="utf-8") as handle:
        ref_cfg = json.load(handle)

    accelerator, devices, device_name = ref_module.resolve_training_device()
    print(f"[INFO] query mapping device: {accelerator}{f' ({device_name})' if device_name else ''}")

    query_original = sc.read_h5ad(QUERY_INPUT_H5AD)
    query_original.var_names_make_unique()
    query_original.obs_names_make_unique()
    query_work = make_query_work_adata(query_original, ref_module, ref_cfg)

    scvi.model.SCANVI.prepare_query_anndata(query_work, str(ref_artifacts.compat_model_dir))
    scanvi_query = scvi.model.SCANVI.load_query_data(query_work, str(ref_artifacts.compat_model_dir))

    train_kwargs: dict[str, Any] = {
        "max_epochs": QUERY_MAX_EPOCHS,
        "batch_size": QUERY_BATCH_SIZE,
        "early_stopping": True,
        "early_stopping_patience": 30,
        "plan_kwargs": {"lr": QUERY_LEARNING_RATE, "weight_decay": QUERY_WEIGHT_DECAY},
        "accelerator": accelerator,
        "devices": 1,
    }
    if accelerator == "gpu" and devices is not None:
        train_kwargs["devices"] = devices

    scanvi_query.train(**train_kwargs)

    latent = np.asarray(scanvi_query.get_latent_representation(query_work), dtype=np.float32)
    l3_pred = pd.Series(scanvi_query.predict(query_work), index=query_work.obs_names, dtype="object")
    l3_pred_guard, n_guard = ref_module.collapse_non_nose_gcb_predictions_to_memory(l3_pred, query_work.obs[ref_module.TISSUE_KEY])
    l2_pred = pd.Series(ref_module.derive_l2_from_l3_predictions(l3_pred_guard.astype(str)), index=query_work.obs_names, dtype="object")
    soft = scanvi_query.predict(query_work, soft=True)
    if hasattr(soft, "to_numpy"):
        soft = soft.to_numpy()
    soft = np.asarray(soft, dtype=np.float32)
    confidence = np.asarray(soft.max(axis=1), dtype=np.float32)

    umap_op = joblib.load(ref_artifacts.compat_umap)
    umap_coords = np.asarray(umap_op.transform(latent), dtype=np.float32)

    mapped = restore_query_metadata(
        query_original=query_original,
        obs_names=query_work.obs_names,
        latent=latent,
        umap_coords=umap_coords,
        l3_pred=l3_pred_guard,
        l2_pred=l2_pred,
        confidence=confidence,
        ref_module=ref_module,
    )
    ref_module.sanitize_anndata_for_legacy_h5ad(mapped)
    mapped.write_h5ad(QUERY_OUTPUT_H5AD, compression="gzip")

    summary = {
        "input_query_h5ad": str(QUERY_INPUT_H5AD),
        "output_query_h5ad": str(QUERY_OUTPUT_H5AD),
        "n_cells": int(mapped.n_obs),
        "n_l2_labels": int(mapped.obs["Cell_Type_L2_final"].astype(str).nunique()),
        "median_confidence": float(np.nanmedian(confidence)),
        "n_non_nose_gc_b_prediction_guardrail": int(n_guard),
    }
    write_json(QUERY_DIR / "query_mapping_summary_20260419.json", summary)
    print(f"[OK] query mapped output saved: {QUERY_OUTPUT_H5AD}")
    return QueryArtifacts(output_h5ad=QUERY_OUTPUT_H5AD, n_cells=int(mapped.n_obs), n_l2_labels=int(mapped.obs["Cell_Type_L2_final"].astype(str).nunique()))


def plot_merge_overview(adata_all: AnnData) -> None:
    ensure_dir(MERGE_FIG_DIR)
    ref_mask = adata_all.obs["data_source"].astype(str).eq("reference")
    qry_mask = adata_all.obs["data_source"].astype(str).eq("query")

    fig, axes = plt.subplots(2, 3, figsize=(20, 13))
    sc.pl.umap(
        adata_all,
        color="data_source",
        ax=axes[0, 0],
        show=False,
        title="Data Source",
        s=5,
        palette={"reference": "#1f77b4", "query": "#ff7f0e"},
    )

    if "Cell_Type_L2" in adata_all.obs.columns:
        adata_all.obs["_ref_l2_only"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
        adata_all.obs.loc[ref_mask, "_ref_l2_only"] = adata_all.obs.loc[ref_mask, "Cell_Type_L2"].astype(str).values
        adata_all.obs["_ref_l2_only"] = adata_all.obs["_ref_l2_only"].astype("category")
        sc.pl.umap(adata_all, color="_ref_l2_only", ax=axes[0, 1], show=False, title="Reference L2 (ref only)", legend_loc="right margin", s=5)
        adata_all.obs.drop(columns=["_ref_l2_only"], inplace=True)

    if "Cell_Type_L2_final" in adata_all.obs.columns:
        adata_all.obs["_qry_l2_only"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
        adata_all.obs.loc[qry_mask, "_qry_l2_only"] = adata_all.obs.loc[qry_mask, "Cell_Type_L2_final"].astype(str).values
        adata_all.obs["_qry_l2_only"] = adata_all.obs["_qry_l2_only"].astype("category")
        sc.pl.umap(adata_all, color="_qry_l2_only", ax=axes[0, 2], show=False, title="Query L2 final (query only)", legend_loc="right margin", s=5)
        adata_all.obs.drop(columns=["_qry_l2_only"], inplace=True)

    if "mapping_confidence" in adata_all.obs.columns:
        adata_all.obs["_qry_conf_only"] = np.nan
        adata_all.obs.loc[qry_mask, "_qry_conf_only"] = pd.to_numeric(adata_all.obs.loc[qry_mask, "mapping_confidence"], errors="coerce").values
        sc.pl.umap(adata_all, color="_qry_conf_only", ax=axes[1, 0], show=False, title="Query confidence", cmap="viridis", vmin=0, vmax=1, s=5)
        adata_all.obs.drop(columns=["_qry_conf_only"], inplace=True)

    sc.pl.umap(adata_all, color="sample", ax=axes[1, 1], show=False, title="Sample", s=5)

    marker = "CD19"
    use_raw = marker not in adata_all.var_names and adata_all.raw is not None and marker in adata_all.raw.var_names
    if marker in adata_all.var_names or use_raw:
        sc.pl.umap(adata_all, color=marker, ax=axes[1, 2], show=False, title=marker, cmap="Reds", s=5, use_raw=use_raw)

    plt.tight_layout()
    plt.savefig(MERGED_OVERVIEW_FIG, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def run_merge(ref_module: Any, ref_artifacts: ReferenceArtifacts, query_artifacts: QueryArtifacts) -> MergeArtifacts:
    print_banner("Stage 3/5: full-gene merge with same-space query UMAP")
    ensure_dir(MERGE_DIR)
    ensure_dir(MERGE_FIG_DIR)

    adata_ref = sc.read_h5ad(ref_artifacts.compat_h5ad)
    adata_qry = sc.read_h5ad(query_artifacts.output_h5ad)
    adata_ref.var_names_make_unique()
    adata_ref.obs_names_make_unique()
    adata_qry.var_names_make_unique()
    adata_qry.obs_names_make_unique()

    if "X_umap_mapped" in adata_qry.obsm:
        adata_qry.obsm["X_umap"] = np.asarray(adata_qry.obsm["X_umap_mapped"]).copy()

    if adata_ref.raw is None:
        raise ValueError("Fresh reference h5ad is missing .raw; cannot build full-gene merge")

    ref_full_genes = pd.Index(adata_ref.raw.var_names)
    ref_full_var = adata_ref.raw.var.copy()
    ref_full_X = _to_csr(adata_ref.raw.X)
    ref_hvg_genes = pd.Index(adata_ref.var_names)

    if "counts" in adata_qry.layers:
        qry_full_X = _to_csr(adata_qry.layers["counts"])
    else:
        qry_full_X = _to_csr(adata_qry.X)
    qry_gene_axis = pd.Index(adata_qry.var_names)

    keep_mask, qry_genes_kept, _ = filter_genes_min_cells(qry_full_X, qry_gene_axis, QUERY_MIN_CELLS)
    qry_full_X_filt = qry_full_X[:, keep_mask]
    qry_gene_axis_filt = qry_genes_kept
    qry_full_X_aligned, align_stats = align_full_to_reference_genes(qry_full_X_filt, qry_gene_axis_filt, ref_full_genes)

    hvgs_for_merge = [g for g in ref_hvg_genes.tolist() if g in set(qry_gene_axis_filt.tolist())]
    if not hvgs_for_merge:
        raise ValueError("No overlapping HVGs found between fresh reference and mapped query")

    adata_ref_hvg = adata_ref[:, hvgs_for_merge].copy()
    adata_qry_hvg = adata_qry[:, hvgs_for_merge].copy()

    if "Cell_Type_L2" not in adata_ref_hvg.obs.columns:
        if "Cell_Type_L2_input" in adata_ref_hvg.obs.columns:
            adata_ref_hvg.obs["Cell_Type_L2"] = adata_ref_hvg.obs["Cell_Type_L2_input"].astype(object)
        elif "L2_scanvi_pred" in adata_ref_hvg.obs.columns:
            adata_ref_hvg.obs["Cell_Type_L2"] = adata_ref_hvg.obs["L2_scanvi_pred"].astype(object)
        else:
            adata_ref_hvg.obs["Cell_Type_L2"] = pd.NA

    if "Cell_Type_L2_pred" not in adata_qry_hvg.obs.columns:
        if "L2_scanvi_pred" in adata_qry_hvg.obs.columns:
            adata_qry_hvg.obs["Cell_Type_L2_pred"] = adata_qry_hvg.obs["L2_scanvi_pred"].astype(object)
        else:
            adata_qry_hvg.obs["Cell_Type_L2_pred"] = pd.NA

    if "Cell_Type_L2_final" not in adata_qry_hvg.obs.columns:
        if "Cell_Type_L2_pred" in adata_qry_hvg.obs.columns:
            adata_qry_hvg.obs["Cell_Type_L2_final"] = adata_qry_hvg.obs["Cell_Type_L2_pred"].astype(object)
        else:
            adata_qry_hvg.obs["Cell_Type_L2_final"] = pd.NA

    if "mapping_confidence" not in adata_qry_hvg.obs.columns:
        if "L2_scanvi_confidence" in adata_qry_hvg.obs.columns:
            adata_qry_hvg.obs["mapping_confidence"] = pd.to_numeric(adata_qry_hvg.obs["L2_scanvi_confidence"], errors="coerce")
        elif "L3_scanvi_confidence" in adata_qry_hvg.obs.columns:
            adata_qry_hvg.obs["mapping_confidence"] = pd.to_numeric(adata_qry_hvg.obs["L3_scanvi_confidence"], errors="coerce")
        else:
            adata_qry_hvg.obs["mapping_confidence"] = np.nan

    for col in ["Cell_Type_L2", "Cell_Type_L2_pred", "Cell_Type_L2_final"]:
        if col not in adata_ref_hvg.obs.columns:
            adata_ref_hvg.obs[col] = pd.NA
        if col not in adata_qry_hvg.obs.columns:
            adata_qry_hvg.obs[col] = pd.NA
    if "mapping_confidence" not in adata_ref_hvg.obs.columns:
        adata_ref_hvg.obs["mapping_confidence"] = np.nan
    if "mapping_confidence" not in adata_qry_hvg.obs.columns:
        adata_qry_hvg.obs["mapping_confidence"] = np.nan

    critical_obsm = {"X_umap", "X_scANVI_L2", "X_scanvi", "X_umap_mapped"}
    for adata in (adata_ref_hvg, adata_qry_hvg):
        for key in list(adata.obsm.keys()):
            if key not in critical_obsm:
                del adata.obsm[key]
        for key in list(adata.obsp.keys()):
            del adata.obsp[key]
        adata.uns = {}
        adata.obs_names_make_unique()

    adata_ref_hvg.obs_names = pd.Index([f"ref::{x}" for x in adata_ref_hvg.obs_names])
    adata_qry_hvg.obs_names = pd.Index([f"qry::{x}" for x in adata_qry_hvg.obs_names])

    adata_all = sc.concat(
        {"reference": adata_ref_hvg, "query": adata_qry_hvg},
        axis=0,
        join="inner",
        merge="unique",
        label="data_source",
    )

    if "X_umap" not in adata_all.obsm or adata_all.obsm["X_umap"].shape[0] != adata_all.n_obs:
        raise RuntimeError("X_umap was lost or misaligned during reference-query concat")

    X_full_merged = sp_vstack([ref_full_X, qry_full_X_aligned], format="csr")
    adata_all.raw = AnnData(X=X_full_merged, obs=adata_all.obs.copy(), var=ref_full_var.copy())

    hvg_pos = ref_full_genes.get_indexer(pd.Index(hvgs_for_merge))
    if (hvg_pos < 0).any():
        raise RuntimeError("Internal merge error: some HVGs are absent from reference full gene universe")
    counts_hvg = X_full_merged[:, hvg_pos]
    adata_all.layers["counts"] = csr_matrix(counts_hvg)
    adata_all.layers["log1p"] = compute_log1p_layer(csr_matrix(counts_hvg))

    adata_all.obs["Cell_Type_L2_unified"] = pd.NA
    ref_mask = adata_all.obs["data_source"].astype(str).eq("reference")
    qry_mask = adata_all.obs["data_source"].astype(str).eq("query")
    if "Cell_Type_L2" in adata_all.obs.columns:
        adata_all.obs.loc[ref_mask, "Cell_Type_L2_unified"] = adata_all.obs.loc[ref_mask, "Cell_Type_L2"].astype(str).values
    if "Cell_Type_L2_final" in adata_all.obs.columns:
        adata_all.obs.loc[qry_mask, "Cell_Type_L2_unified"] = adata_all.obs.loc[qry_mask, "Cell_Type_L2_final"].astype(str).values
    adata_all.obs["Cell_Type_L2_unified"] = adata_all.obs["Cell_Type_L2_unified"].astype("category")

    ref_module.sanitize_anndata_for_legacy_h5ad(adata_all)
    adata_all.write_h5ad(MERGED_OUTPUT_H5AD, compression="gzip")
    plot_merge_overview(adata_all)

    merge_summary = {
        "reference_h5ad": str(ref_artifacts.compat_h5ad),
        "query_h5ad": str(query_artifacts.output_h5ad),
        "merged_h5ad": str(MERGED_OUTPUT_H5AD),
        "n_cells": int(adata_all.n_obs),
        "n_reference": int(ref_mask.sum()),
        "n_query": int(qry_mask.sum()),
        "n_hvgs": int(adata_all.n_vars),
        "n_full_genes": int(adata_all.raw.n_vars) if adata_all.raw is not None else None,
        "query_alignment_stats": align_stats,
        "merged_overview_figure": str(MERGED_OVERVIEW_FIG),
    }
    write_json(MERGED_SUMMARY_JSON, merge_summary)
    print(f"[OK] merged h5ad saved: {MERGED_OUTPUT_H5AD}")
    return MergeArtifacts(
        output_h5ad=MERGED_OUTPUT_H5AD,
        n_cells=int(adata_all.n_obs),
        n_reference=int(ref_mask.sum()),
        n_query=int(qry_mask.sum()),
        n_hvgs=int(adata_all.n_vars),
        n_full_genes=int(adata_all.raw.n_vars) if adata_all.raw is not None else 0,
    )


def run_python_subprocess(args: list[str]) -> None:
    print(f"[run] {' '.join(args)}")
    subprocess.run(args, check=True)


def run_schpl_and_followup(merge_artifacts: MergeArtifacts) -> None:
    schpl_python = resolve_schpl_python_executable()
    print(f"[INFO] schpl/follow-up python: {schpl_python}")

    print_banner("Stage 4/5: rerun B-cell schpl wrapper")
    run_python_subprocess(
        [
            str(schpl_python),
            str(SCHPL_SCRIPT),
            "--input-h5ad",
            str(merge_artifacts.output_h5ad),
            "--version",
            SCHPL_VERSION,
            "--output-dir",
            str(SCHPL_OUTPUT_DIR),
            "--allow-overwrite",
        ]
    )

    print_banner("Stage 5/5: rerun schpl reject follow-up")
    run_python_subprocess(
        [
            str(schpl_python),
            str(FOLLOWUP_SCRIPT),
            "--input-h5ad",
            str(SCHPL_OUTPUT_H5AD),
            "--config-json",
            str(SCHPL_CONFIG_JSON),
            "--output-dir",
            str(FOLLOWUP_OUTPUT_DIR),
        ]
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="B-cell full rerun for same-space query UMAP pairing")
    parser.add_argument("--skip-reference", action="store_true")
    parser.add_argument("--skip-query", action="store_true")
    parser.add_argument("--skip-merge", action="store_true")
    parser.add_argument("--skip-schpl", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_dir(OUTPUT_ROOT)
    ensure_dir(REFERENCE_DIR)
    ensure_dir(QUERY_DIR)
    ensure_dir(MERGE_DIR)

    ref_module: Any | None = None
    ref_artifacts: ReferenceArtifacts | None = None
    query_artifacts: QueryArtifacts | None = None
    merge_artifacts: MergeArtifacts | None = None

    if args.skip_reference:
        ref_module = load_module(REF_MODULE_PATH, "bcell_ref_rerun_20260417_existing")
        ref_artifacts = ReferenceArtifacts(
            output_dir=REFERENCE_DIR,
            h5ad=REFERENCE_DIR / "bcell_reference_c22_c13_c25drop_scanvi_L3_ref_20260419.h5ad",
            compat_h5ad=REFERENCE_ALIAS_H5AD,
            compat_umap=REFERENCE_UMAP_ALIAS,
            compat_model_dir=REFERENCE_SCANVI_ALIAS,
            training_config=REFERENCE_DIR / "training_config.json",
            hvg_file=REFERENCE_DIR / "hvg_genes.txt",
        )
    else:
        ref_module, ref_artifacts = run_reference_rerun()

    assert ref_module is not None and ref_artifacts is not None

    if args.skip_query:
        query_artifacts = QueryArtifacts(output_h5ad=QUERY_OUTPUT_H5AD, n_cells=0, n_l2_labels=0)
    else:
        query_artifacts = run_query_mapping(ref_module, ref_artifacts)

    if args.skip_merge:
        merge_artifacts = MergeArtifacts(output_h5ad=MERGED_OUTPUT_H5AD, n_cells=0, n_reference=0, n_query=0, n_hvgs=0, n_full_genes=0)
    else:
        merge_artifacts = run_merge(ref_module, ref_artifacts, query_artifacts)

    if not args.skip_schpl:
        run_schpl_and_followup(merge_artifacts)

    summary = {
        "reference_dir": str(REFERENCE_DIR),
        "reference_h5ad": str(ref_artifacts.compat_h5ad),
        "reference_model_dir": str(ref_artifacts.compat_model_dir),
        "query_h5ad": str(query_artifacts.output_h5ad),
        "merged_h5ad": str(merge_artifacts.output_h5ad),
        "schpl_output_dir": str(SCHPL_OUTPUT_DIR),
        "schpl_output_h5ad": str(SCHPL_OUTPUT_H5AD),
        "schpl_config_json": str(SCHPL_CONFIG_JSON),
        "followup_output_dir": str(FOLLOWUP_OUTPUT_DIR),
    }
    write_json(SUMMARY_JSON, summary)

    print_banner("DONE: B-cell full rerun for query UMAP pairing")
    print(f"Reference dir : {REFERENCE_DIR}")
    print(f"Query output  : {query_artifacts.output_h5ad}")
    print(f"Merged output : {merge_artifacts.output_h5ad}")
    print(f"SCHPL output  : {SCHPL_OUTPUT_H5AD}")
    print(f"Follow-up dir : {FOLLOWUP_OUTPUT_DIR}")


if __name__ == "__main__":
    main()
