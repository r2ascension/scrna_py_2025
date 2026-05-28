#!/usr/bin/env python3
"""Shared helpers for normal-airway anatomical-site ML workflows.

The first implementation tranche focuses on:
1. robust metadata contract resolution from heterogeneous h5ad objects,
2. healthy-only site standardization,
3. sample-level feature extraction for composition / latent summaries /
   pseudobulk PCA, and
4. small deterministic toy and real-subset smoke helpers.
"""

from __future__ import annotations

import json
import warnings
from collections import OrderedDict
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

warnings.filterwarnings("ignore", category=FutureWarning)
try:
    ad.settings.allow_write_nullable_strings = True
except Exception:
    pass

WORKSPACE_ROOT = Path("/home/h2048")
DATE_TAG = "20260527"
DEFAULT_OUTPUT_ROOT = WORKSPACE_ROOT / "data" / "py" / DATE_TAG
DEFAULT_LOG_ROOT = WORKSPACE_ROOT / "logs" / DATE_TAG
DEFAULT_COUNTS_LAYER = "counts"
DEFAULT_UNKNOWN_LABEL = "Unknown"

DEFAULT_SITE_MAP = OrderedDict(
    {
        "nose": "nasal",
        "nasal": "nasal",
        "sinus": "sinus",
        "respiratory airway": "bronchus",
        "airway": "bronchus",
        "bronchus": "bronchus",
        "bronchial": "bronchus",
        "trachea": "bronchus",
        "large airway": "bronchus",
        "lung parenchyma": "lung_parenchyma",
        "parenchyma": "lung_parenchyma",
        "alveolar": "lung_parenchyma",
        "lung": "lung_parenchyma",
    }
)
SITE_AXIS_MAP = {
    "nasal": "upper",
    "sinus": "upper",
    "bronchus": "lower",
    "lung_parenchyma": "lower",
}
HEALTHY_VALUES = {"healthy", "control", "normal", "donor", "non-diseased", "nondiseased"}


def default_config() -> dict[str, Any]:
    return {
        "run": {
            "output_root": str(DEFAULT_OUTPUT_ROOT),
            "run_name_prefix": "normal_airway_ml_mvp",
            "random_seed": 20260527,
        },
        "data_contract": {
            "input_h5ad": str(WORKSPACE_ROOT / "data" / "coredata0524" / "coredata0524_allcells_scanvi_l1_20260524.h5ad"),
            "counts_layer": DEFAULT_COUNTS_LAYER,
            "unknown_label": DEFAULT_UNKNOWN_LABEL,
            "sample_key_candidates": ["sample", "Sample", "orig.ident"],
            "dataset_key_candidates": ["dataset", "study", "batch"],
            "study_key_candidates": ["study", "dataset", "batch"],
            "batch_key_candidates": ["batch", "dataset", "study"],
            "tissue_key_candidates": ["tissue", "anatomical_site", "site", "tissue_level_2"],
            "condition_key_candidates": ["condition", "disease_level_2", "group", "disease", "status"],
            "cell_type_key_candidates": ["cell_type_L2", "cell_type_L3", "cell_type_L1", "Annotation_2"],
            "latent_key_candidates": ["X_scvi", "X_scanvi"],
            "healthy_values": sorted(HEALTHY_VALUES),
            "site_map": dict(DEFAULT_SITE_MAP),
            "site_axis_map": dict(SITE_AXIS_MAP),
        },
        "feature_export": {
            "cell_type_key": "cell_type_L2",
            "min_cells_per_sample_celltype": 10,
            "min_samples_per_celltype": 3,
            "pseudobulk_top_genes": 200,
            "pseudobulk_n_pcs": 5,
            "latent_max_dims": 16,
            "epsilon": 1e-6,
        },
        "train": {
            "label_column": "site_label",
            "dataset_group_column": "dataset",
            "n_splits": 4,
            "logistic_max_iter": 2000,
            "random_forest_estimators": 300,
            "random_seed": 20260527,
            "permutation_repeats": 10,
        },
        "smoke": {
            "real_input_h5ad": str(
                WORKSPACE_ROOT
                / "data"
                / "R"
                / "0508"
                / "bcell_tissue_comparison_v2_6_8_c22_c13_c25_c14drop_20260508"
                / "bcell_tissue_comparison_final_fullgene.h5ad"
            ),
            "toy_healthy_samples_per_site": 3,
            "toy_disease_samples_per_site": 1,
            "toy_cells_per_sample": 36,
            "toy_n_genes": 120,
            "toy_latent_dim": 16,
            "real_samples_per_site": 2,
            "real_max_cells_per_sample": 180,
            "real_min_cells_per_sample": 30,
        },
    }


# ---------------------------------------------------------------------------
# Generic utilities
# ---------------------------------------------------------------------------


def deep_get(obj: dict[str, Any], *keys: str, default: Any = None) -> Any:
    cur: Any = obj
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def merge_nested_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge_nested_dict(out[key], value)
        else:
            out[key] = value
    return out


def load_config(config_path: Path | str | None = None) -> dict[str, Any]:
    cfg = default_config()
    if config_path is None:
        return cfg
    path = Path(config_path)
    text = path.read_text(encoding="utf-8")
    if yaml is not None:
        override = yaml.safe_load(text) or {}
    else:  # pragma: no cover
        override = json.loads(text)
    return merge_nested_dict(cfg, override)


def timestamp_slug() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def ensure_dir(path: Path | str) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_json(payload: dict[str, Any], path: Path | str) -> Path:
    target = Path(path)
    target.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    return target


def normalize_string_series(values: Iterable[Any] | pd.Series, missing_value: str = "") -> pd.Series:
    series = values if isinstance(values, pd.Series) else pd.Series(list(values))
    out = series.astype("string")
    out = out.fillna(missing_value)
    out = out.astype(str).str.strip()
    out = out.replace({"nan": missing_value, "None": missing_value, "<NA>": missing_value})
    return out


def infer_first_existing(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    colset = set(columns)
    for candidate in candidates:
        if candidate in colset:
            return candidate
    return None


def safe_mode(series: pd.Series, fallback: str = DEFAULT_UNKNOWN_LABEL) -> str:
    clean = normalize_string_series(series, missing_value="")
    clean = clean[clean != ""]
    if clean.empty:
        return fallback
    counts = clean.value_counts(sort=True)
    return str(counts.index[0])


def summarize_top_counts(series: pd.Series, top_n: int = 8) -> dict[str, int]:
    clean = normalize_string_series(series, missing_value="")
    clean = clean[clean != ""]
    return {str(k): int(v) for k, v in clean.value_counts().head(top_n).items()}


# ---------------------------------------------------------------------------
# Metadata contract resolution
# ---------------------------------------------------------------------------


def standardize_site_series(
    values: pd.Series,
    site_map: dict[str, str] | None = None,
    unknown_value: str = DEFAULT_UNKNOWN_LABEL,
) -> pd.Series:
    mapping = {str(k).strip().lower(): str(v).strip() for k, v in (site_map or DEFAULT_SITE_MAP).items()}
    raw = normalize_string_series(values, missing_value=unknown_value)
    normalized = raw.str.lower().map(mapping)
    out = normalized.fillna(unknown_value).astype(str)
    return out


def build_healthy_mask(
    obs: pd.DataFrame,
    condition_key_candidates: Iterable[str],
    healthy_values: Iterable[str] | None = None,
) -> tuple[pd.Series, str | None]:
    healthy_vocab = {str(v).strip().lower() for v in (healthy_values or HEALTHY_VALUES)}
    condition_col = infer_first_existing(obs.columns, condition_key_candidates)
    if condition_col is None:
        return pd.Series([False] * len(obs), index=obs.index, dtype=bool), None
    condition = normalize_string_series(obs[condition_col], missing_value="")
    mask = condition.str.lower().isin(healthy_vocab)
    return mask.astype(bool), condition_col


def prepare_obs_contract(
    obs: pd.DataFrame,
    available_obsm_keys: Iterable[str] | None,
    cfg: dict[str, Any],
    cell_type_key_override: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    contract_cfg = deep_get(cfg, "data_contract", default={}) or {}
    sample_col = infer_first_existing(obs.columns, deep_get(contract_cfg, "sample_key_candidates", default=["sample"]))
    dataset_col = infer_first_existing(obs.columns, deep_get(contract_cfg, "dataset_key_candidates", default=["dataset", "study"]))
    study_col = infer_first_existing(obs.columns, deep_get(contract_cfg, "study_key_candidates", default=["study", "dataset"]))
    batch_col = infer_first_existing(obs.columns, deep_get(contract_cfg, "batch_key_candidates", default=["batch", "dataset"]))
    tissue_col = infer_first_existing(obs.columns, deep_get(contract_cfg, "tissue_key_candidates", default=["tissue"]))
    cell_type_col = cell_type_key_override or infer_first_existing(
        obs.columns, deep_get(contract_cfg, "cell_type_key_candidates", default=["cell_type_L2", "cell_type_L3"])
    )
    healthy_mask, condition_col = build_healthy_mask(
        obs,
        deep_get(contract_cfg, "condition_key_candidates", default=["condition", "disease_level_2"]),
        healthy_values=deep_get(contract_cfg, "healthy_values", default=sorted(HEALTHY_VALUES)),
    )
    site_label = (
        standardize_site_series(
            obs[tissue_col] if tissue_col else pd.Series([DEFAULT_UNKNOWN_LABEL] * len(obs), index=obs.index),
            site_map=deep_get(contract_cfg, "site_map", default=dict(DEFAULT_SITE_MAP)),
            unknown_value=deep_get(contract_cfg, "unknown_label", default=DEFAULT_UNKNOWN_LABEL),
        )
        if tissue_col
        else pd.Series([DEFAULT_UNKNOWN_LABEL] * len(obs), index=obs.index, dtype="string")
    )
    site_axis_map = {str(k): str(v) for k, v in deep_get(contract_cfg, "site_axis_map", default=dict(SITE_AXIS_MAP)).items()}
    site_axis = site_label.map(site_axis_map).fillna(deep_get(contract_cfg, "unknown_label", default=DEFAULT_UNKNOWN_LABEL))

    prepared = obs.copy()
    prepared["sample_resolved"] = normalize_string_series(
        obs[sample_col] if sample_col else pd.Series([f"sample_{i:05d}" for i in range(len(obs))], index=obs.index),
        missing_value="sample_missing",
    )
    prepared["dataset_resolved"] = normalize_string_series(
        obs[dataset_col] if dataset_col else prepared["sample_resolved"],
        missing_value="dataset_missing",
    )
    prepared["study_resolved"] = normalize_string_series(
        obs[study_col] if study_col else prepared["dataset_resolved"],
        missing_value="study_missing",
    )
    prepared["batch_resolved"] = normalize_string_series(
        obs[batch_col] if batch_col else prepared["dataset_resolved"],
        missing_value="batch_missing",
    )
    prepared["condition_resolved"] = normalize_string_series(
        obs[condition_col] if condition_col else pd.Series([DEFAULT_UNKNOWN_LABEL] * len(obs), index=obs.index),
        missing_value=DEFAULT_UNKNOWN_LABEL,
    )
    prepared["raw_tissue_label"] = normalize_string_series(
        obs[tissue_col] if tissue_col else pd.Series([DEFAULT_UNKNOWN_LABEL] * len(obs), index=obs.index),
        missing_value=DEFAULT_UNKNOWN_LABEL,
    )
    prepared["site_label"] = site_label.astype(str).values
    prepared["site_axis"] = site_axis.astype(str).values
    prepared["is_healthy"] = healthy_mask.astype(bool).values
    prepared["cell_type_resolved"] = normalize_string_series(
        obs[cell_type_col] if cell_type_col else pd.Series([DEFAULT_UNKNOWN_LABEL] * len(obs), index=obs.index),
        missing_value=DEFAULT_UNKNOWN_LABEL,
    )

    latent_key = infer_first_existing(list(available_obsm_keys or []), deep_get(contract_cfg, "latent_key_candidates", default=["X_scvi", "X_scanvi"]))
    analysis_mask = prepared["is_healthy"].astype(bool) & (prepared["site_label"].astype(str) != deep_get(contract_cfg, "unknown_label", default=DEFAULT_UNKNOWN_LABEL))

    contract = {
        "sample_col": sample_col,
        "dataset_col": dataset_col,
        "study_col": study_col,
        "batch_col": batch_col,
        "tissue_col": tissue_col,
        "condition_col": condition_col,
        "cell_type_col": cell_type_col,
        "latent_key": latent_key,
        "analysis_mask": analysis_mask,
        "n_analysis_cells": int(analysis_mask.sum()),
    }
    return prepared, contract


# ---------------------------------------------------------------------------
# AnnData access helpers
# ---------------------------------------------------------------------------


def ensure_counts_layer(adata: ad.AnnData, counts_layer: str = DEFAULT_COUNTS_LAYER) -> dict[str, Any]:
    if counts_layer in adata.layers:
        return {"counts_layer": counts_layer, "created": False, "source": counts_layer}
    if adata.raw is not None and adata.raw.n_vars >= adata.n_vars:
        raw_index = pd.Index(adata.raw.var_names.astype(str))
        current_index = pd.Index(adata.var_names.astype(str))
        if current_index.isin(raw_index).all():
            pos = raw_index.get_indexer(current_index)
            adata.layers[counts_layer] = adata.raw.X[:, pos].copy()
            return {"counts_layer": counts_layer, "created": True, "source": "raw"}
    adata.layers[counts_layer] = adata.X.copy()
    return {"counts_layer": counts_layer, "created": True, "source": "X"}


def load_backed_obs_summary(h5ad_path: Path | str, top_n: int = 8) -> dict[str, Any]:
    path = Path(h5ad_path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        obs = backed.obs.copy()
        summary = {
            "h5ad_path": str(path),
            "shape": [int(backed.n_obs), int(backed.n_vars)],
            "obs_ncols": int(obs.shape[1]),
            "obs_columns": list(obs.columns),
            "var_ncols": int(backed.var.shape[1]),
            "layers": list(backed.layers.keys()),
            "obsm": list(backed.obsm.keys()),
        }
        return summary
    finally:
        try:
            backed.file.close()
        except Exception:
            pass


def load_subset_by_obs_names(h5ad_path: Path | str, obs_names: Iterable[str]) -> ad.AnnData:
    path = Path(h5ad_path)
    requested = pd.Index(pd.Series(list(obs_names), dtype="string").astype(str))
    backed = ad.read_h5ad(path, backed="r")
    try:
        index = pd.Index(backed.obs_names.astype(str))
        pos = index.get_indexer(requested)
        pos = pos[pos >= 0]
        subset = backed[pos]
        if hasattr(subset, "to_memory"):
            return subset.to_memory()
        return subset.copy()
    finally:
        try:
            backed.file.close()
        except Exception:
            pass


def load_real_smoke_subset(
    h5ad_path: Path | str,
    cfg: dict[str, Any],
    seed: int,
) -> tuple[ad.AnnData, dict[str, Any]]:
    path = Path(h5ad_path)
    backed = ad.read_h5ad(path, backed="r")
    try:
        obs = backed.obs.copy()
        prepared, contract = prepare_obs_contract(obs, list(backed.obsm.keys()), cfg)
        sample_col = "sample_resolved"
        site_col = "site_label"
        analysis = prepared.loc[contract["analysis_mask"]].copy()
        min_cells = int(deep_get(cfg, "smoke", "real_min_cells_per_sample", default=30))
        max_samples_per_site = int(deep_get(cfg, "smoke", "real_samples_per_site", default=2))
        counts = (
            analysis.groupby([site_col, sample_col], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        thresholds = []
        for candidate in [min_cells, max(1, min_cells // 2), 10, 5, 1]:
            if candidate not in thresholds:
                thresholds.append(candidate)

        selected_samples: list[str] = []
        selected_threshold: int | None = None
        selected_site_counts: dict[str, int] = {}
        for threshold in thresholds:
            filtered = counts[counts["n_cells"] >= threshold].copy()
            candidate_samples: list[str] = []
            candidate_site_counts: dict[str, int] = {}
            for site in ["nasal", "sinus", "bronchus", "lung_parenchyma"]:
                site_rows = filtered[filtered[site_col] == site].sort_values(["n_cells", sample_col], ascending=[False, True])
                picked = site_rows[sample_col].head(max_samples_per_site).astype(str).tolist()
                if picked:
                    candidate_samples.extend(picked)
                    candidate_site_counts[site] = len(picked)
            if len(candidate_site_counts) >= 2:
                selected_samples = sorted(set(candidate_samples))
                selected_threshold = threshold
                selected_site_counts = candidate_site_counts
                break

        if not selected_samples:
            filtered = counts.sort_values([site_col, "n_cells", sample_col], ascending=[True, False, True]).copy()
            candidate_site_counts = {}
            for site in ["nasal", "sinus", "bronchus", "lung_parenchyma"]:
                site_rows = filtered[filtered[site_col] == site]
                picked = site_rows[sample_col].head(max_samples_per_site).astype(str).tolist()
                if picked:
                    selected_samples.extend(picked)
                    candidate_site_counts[site] = len(picked)
            selected_samples = sorted(set(selected_samples))
            selected_threshold = 0
            selected_site_counts = candidate_site_counts
        if not selected_samples:
            raise ValueError("No real smoke samples passed healthy/site/min-cell filters")
        keep_obs = analysis[analysis[sample_col].astype(str).isin(selected_samples)].copy()
        if keep_obs[site_col].nunique() < 2:
            raise ValueError(
                "Real smoke subset resolved to fewer than 2 anatomical sites; "
                f"selected sites={sorted(keep_obs[site_col].astype(str).unique().tolist())}"
            )
        subset = backed[keep_obs.index]
        adata = subset.to_memory() if hasattr(subset, "to_memory") else subset.copy()
    finally:
        try:
            backed.file.close()
        except Exception:
            pass

    max_cells_per_sample = int(deep_get(cfg, "smoke", "real_max_cells_per_sample", default=180))
    rng = np.random.default_rng(seed)
    prepared_sub, contract_sub = prepare_obs_contract(adata.obs.copy(), list(adata.obsm.keys()), cfg)
    selected_idx: list[str] = []
    for sample_id, sample_obs in prepared_sub.loc[contract_sub["analysis_mask"]].groupby("sample_resolved", observed=True):
        idx = sample_obs.index.to_numpy()
        if len(idx) > max_cells_per_sample:
            idx = rng.choice(idx, size=max_cells_per_sample, replace=False)
        selected_idx.extend(idx.tolist())
    adata = adata[selected_idx].copy()
    return adata, {
        "source_h5ad": str(path),
        "selected_samples": selected_samples,
        "n_selected_samples": len(selected_samples),
        "selected_min_cells_threshold": selected_threshold,
        "selected_site_sample_counts": selected_site_counts,
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "site_counts": prepared_sub.loc[selected_idx, "site_label"].value_counts().to_dict(),
    }


# ---------------------------------------------------------------------------
# Toy data
# ---------------------------------------------------------------------------


def make_toy_airway_adata(cfg: dict[str, Any], seed: int | None = None) -> ad.AnnData:
    smoke_cfg = deep_get(cfg, "smoke", default={}) or {}
    rng = np.random.default_rng(int(seed if seed is not None else deep_get(cfg, "run", "random_seed", default=20260527)))
    sites = ["nasal", "sinus", "bronchus", "lung_parenchyma"]
    raw_tissue = {
        "nasal": "nose",
        "sinus": "sinus",
        "bronchus": "respiratory airway",
        "lung_parenchyma": "lung parenchyma",
    }
    celltypes = ["Secretory_Lineage", "Alveolar_Macrophage", "Fibroblast"]
    n_healthy = int(smoke_cfg.get("toy_healthy_samples_per_site", 3))
    n_disease = int(smoke_cfg.get("toy_disease_samples_per_site", 1))
    cells_per_sample = int(smoke_cfg.get("toy_cells_per_sample", 36))
    n_genes = int(smoke_cfg.get("toy_n_genes", 120))
    latent_dim = int(smoke_cfg.get("toy_latent_dim", 16))

    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]

    def _allocate_marker_blocks(names: list[str], start: int, block_size: int) -> tuple[dict[str, np.ndarray], int]:
        blocks: dict[str, np.ndarray] = {}
        cursor = start
        for name in names:
            end = min(cursor + block_size, n_genes)
            blocks[name] = np.arange(cursor, end)
            cursor = end
        return blocks, cursor

    total_groups = len(sites) + len(celltypes) + 1
    block_size = max(4, n_genes // max(total_groups + 1, 1))
    site_marker_blocks, cursor = _allocate_marker_blocks(sites, 0, block_size)
    celltype_marker_blocks, cursor = _allocate_marker_blocks(celltypes, cursor, block_size)
    disease_block = np.arange(cursor, min(cursor + block_size, n_genes))

    obs_rows: list[dict[str, Any]] = []
    matrices: list[np.ndarray] = []
    latent_rows: list[np.ndarray] = []

    site_vectors = {
        "nasal": np.array([2.0, 0.0, 0.0, 0.0]),
        "sinus": np.array([0.0, 2.0, 0.0, 0.0]),
        "bronchus": np.array([0.0, 0.0, 2.0, 0.0]),
        "lung_parenchyma": np.array([0.0, 0.0, 0.0, 2.0]),
    }
    celltype_vectors = {
        "Secretory_Lineage": np.array([1.2, 0.0, 0.0]),
        "Alveolar_Macrophage": np.array([0.0, 1.2, 0.0]),
        "Fibroblast": np.array([0.0, 0.0, 1.2]),
    }

    sample_counter = 0
    for dataset_id, site in enumerate(sites, start=1):
        for sample_type, n_samples, condition in (("healthy", n_healthy, "Healthy"), ("disease", n_disease, "Disease")):
            for local_idx in range(n_samples):
                sample_counter += 1
                sample_id = f"toy_{site}_{sample_type}_{local_idx+1:02d}"
                dataset = f"toy_dataset_{dataset_id}"
                batch = f"toy_batch_{dataset_id}"
                study = f"toy_study_{1 + (dataset_id % 2)}"
                for cell_idx in range(cells_per_sample):
                    celltype = celltypes[cell_idx % len(celltypes)]
                    lam = np.full(n_genes, 2.0, dtype=float)
                    lam[site_marker_blocks[site]] += 10.0
                    lam[celltype_marker_blocks[celltype]] += 8.0
                    if condition != "Healthy" and disease_block.size:
                        lam[disease_block] += 3.5
                    counts = rng.poisson(lam=lam).astype(np.float32)
                    matrices.append(counts)
                    latent = np.zeros(latent_dim, dtype=np.float32)
                    latent[:4] = site_vectors[site]
                    latent[4:7] = celltype_vectors[celltype]
                    latent[7] = 0.7 if condition == "Healthy" else -0.7
                    latent += rng.normal(0.0, 0.15, size=latent_dim).astype(np.float32)
                    latent_rows.append(latent)
                    obs_rows.append(
                        {
                            "sample": sample_id,
                            "dataset": dataset,
                            "study": study,
                            "batch": batch,
                            "tissue": raw_tissue[site],
                            "condition": condition,
                            "cell_type_L1": "Toy",
                            "cell_type_L2": celltype,
                            "cell_type_L3": f"{celltype}_sub",
                            "scanvi_confidence": float(np.clip(rng.normal(0.995, 0.002), 0.90, 1.0)),
                        }
                    )

    X = np.vstack(matrices).astype(np.float32)
    obs = pd.DataFrame(obs_rows, index=[f"toy_cell_{i:05d}" for i in range(len(obs_rows))])
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers[DEFAULT_COUNTS_LAYER] = adata.X.copy()
    adata.layers["log1p"] = np.log1p(X).astype(np.float32)
    latent_matrix = np.vstack(latent_rows).astype(np.float32)
    adata.obsm["X_scvi"] = latent_matrix
    adata.obsm["X_scanvi"] = latent_matrix.copy()
    return adata


# ---------------------------------------------------------------------------
# Feature engineering
# ---------------------------------------------------------------------------


def build_sample_metadata(prepared_obs: pd.DataFrame, analysis_mask: pd.Series) -> pd.DataFrame:
    frame = prepared_obs.loc[analysis_mask].copy()
    grouped = frame.groupby("sample_resolved", observed=True)
    sample_meta = grouped.agg(
        sample=("sample_resolved", lambda s: safe_mode(s, fallback="sample_missing")),
        dataset=("dataset_resolved", safe_mode),
        study=("study_resolved", safe_mode),
        batch=("batch_resolved", safe_mode),
        site_label=("site_label", safe_mode),
        site_axis=("site_axis", safe_mode),
        condition=("condition_resolved", safe_mode),
        n_cells=("sample_resolved", "size"),
        n_cell_types=("cell_type_resolved", "nunique"),
    )
    sample_meta.index.name = "sample"
    sample_meta = sample_meta.reset_index(drop=True).set_index("sample", drop=False)
    return sample_meta


def aggregate_composition_features(
    prepared_obs: pd.DataFrame,
    analysis_mask: pd.Series,
    epsilon: float = 1e-6,
) -> dict[str, pd.DataFrame]:
    frame = prepared_obs.loc[analysis_mask, ["sample_resolved", "cell_type_resolved"]].copy()
    counts = pd.crosstab(frame["sample_resolved"], frame["cell_type_resolved"])
    counts.index.name = "sample"
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    log_fraction = np.log(fractions + epsilon)
    clr = log_fraction.sub(log_fraction.mean(axis=1), axis=0)

    wide = pd.DataFrame(index=counts.index)
    for col in counts.columns:
        wide[f"composition_count__{col}"] = counts[col].astype(float)
        wide[f"composition_fraction__{col}"] = fractions[col].astype(float)
        wide[f"composition_clr__{col}"] = clr[col].astype(float)
    wide.index.name = "sample"

    long = counts.reset_index().melt(id_vars="sample", var_name="cell_type", value_name="n_cells")
    long["fraction"] = long.apply(lambda row: float(fractions.loc[row["sample"], row["cell_type"]]), axis=1)
    long["clr"] = long.apply(lambda row: float(clr.loc[row["sample"], row["cell_type"]]), axis=1)
    return {"counts": counts, "fractions": fractions, "clr": clr, "wide": wide, "long": long}


def aggregate_latent_summary_features(
    adata: ad.AnnData,
    prepared_obs: pd.DataFrame,
    analysis_mask: pd.Series,
    latent_key: str | None,
    max_dims: int = 16,
) -> dict[str, Any]:
    if latent_key is None or latent_key not in adata.obsm:
        return {"status": "skipped", "reason": "latent_key_missing", "wide": pd.DataFrame()}
    latent = np.asarray(adata.obsm[latent_key])
    if latent.ndim != 2 or latent.shape[0] != adata.n_obs:
        return {"status": "skipped", "reason": "invalid_latent_matrix", "wide": pd.DataFrame()}
    n_dims = int(min(max_dims, latent.shape[1]))
    frame = prepared_obs.loc[analysis_mask, ["sample_resolved", "cell_type_resolved"]].copy()
    frame["row_idx"] = np.flatnonzero(np.asarray(analysis_mask, dtype=bool))

    feature_map: dict[str, dict[str, float]] = {}
    for (sample_id, celltype), sub in frame.groupby(["sample_resolved", "cell_type_resolved"], observed=True):
        coords = latent[sub["row_idx"].to_numpy(), :n_dims]
        mean_vec = coords.mean(axis=0)
        sample_features = feature_map.setdefault(str(sample_id), {})
        for dim_idx, value in enumerate(mean_vec, start=1):
            sample_features[f"latent_mean__{celltype}__dim_{dim_idx:03d}"] = float(value)
    wide = pd.DataFrame.from_dict(feature_map, orient="index").sort_index().fillna(0.0)
    wide.index.name = "sample"
    return {"status": "ok", "latent_key": latent_key, "n_dims": n_dims, "wide": wide}


def _matrix_take_rows(matrix: Any, row_idx: np.ndarray) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix[row_idx].toarray().astype(np.float32)
    return np.asarray(matrix[row_idx], dtype=np.float32)


def aggregate_pseudobulk_pca_features(
    adata: ad.AnnData,
    prepared_obs: pd.DataFrame,
    analysis_mask: pd.Series,
    counts_layer: str,
    min_cells_per_group: int = 10,
    min_samples_per_celltype: int = 3,
    top_genes: int = 200,
    n_pcs: int = 5,
) -> dict[str, Any]:
    if counts_layer not in adata.layers:
        return {"status": "skipped", "reason": f"missing_counts_layer:{counts_layer}", "wide": pd.DataFrame()}

    counts_matrix = adata.layers[counts_layer]
    frame = prepared_obs.loc[analysis_mask, ["sample_resolved", "cell_type_resolved"]].copy()
    frame["row_idx"] = np.flatnonzero(np.asarray(analysis_mask, dtype=bool))
    group_sizes = (
        frame.groupby(["sample_resolved", "cell_type_resolved"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    valid_groups = group_sizes[group_sizes["n_cells"] >= min_cells_per_group]
    if valid_groups.empty:
        return {"status": "skipped", "reason": "no_groups_passed_min_cells", "wide": pd.DataFrame()}

    feature_map: dict[str, dict[str, float]] = {}
    skipped_celltypes: dict[str, str] = {}
    for celltype, celltype_rows in valid_groups.groupby("cell_type_resolved", observed=True):
        samples = celltype_rows["sample_resolved"].astype(str).tolist()
        if len(samples) < min_samples_per_celltype:
            skipped_celltypes[str(celltype)] = "insufficient_samples"
            continue
        sample_vectors: list[np.ndarray] = []
        sample_ids: list[str] = []
        for sample_id in samples:
            idx = frame.loc[
                (frame["sample_resolved"].astype(str) == str(sample_id))
                & (frame["cell_type_resolved"].astype(str) == str(celltype)),
                "row_idx",
            ].to_numpy(dtype=int)
            if idx.size < min_cells_per_group:
                continue
            summed = _matrix_take_rows(counts_matrix, idx).sum(axis=0)
            sample_vectors.append(np.asarray(summed, dtype=np.float32))
            sample_ids.append(str(sample_id))
        if len(sample_ids) < min_samples_per_celltype:
            skipped_celltypes[str(celltype)] = "insufficient_post_agg_samples"
            continue
        pb = np.vstack(sample_vectors)
        library = pb.sum(axis=1, keepdims=True)
        library[library == 0] = 1.0
        log_cpm = np.log1p((pb / library) * 1e4)
        gene_var = np.var(log_cpm, axis=0)
        top_idx = np.argsort(gene_var)[::-1][: min(top_genes, log_cpm.shape[1])]
        X = log_cpm[:, top_idx]
        X = StandardScaler(with_mean=True, with_std=True).fit_transform(X)
        n_components = int(min(n_pcs, max(1, X.shape[0] - 1), X.shape[1]))
        if n_components < 1:
            skipped_celltypes[str(celltype)] = "pca_rank_zero"
            continue
        coords = PCA(n_components=n_components, random_state=42).fit_transform(X)
        for row_idx, sample_id in enumerate(sample_ids):
            sample_features = feature_map.setdefault(sample_id, {})
            for pc_idx in range(n_components):
                sample_features[f"pb_pca__{celltype}__PC{pc_idx + 1}"] = float(coords[row_idx, pc_idx])
    wide = pd.DataFrame.from_dict(feature_map, orient="index").sort_index().fillna(0.0)
    wide.index.name = "sample"
    return {
        "status": "ok" if not wide.empty else "skipped",
        "skipped_celltypes": skipped_celltypes,
        "wide": wide,
    }


def join_feature_blocks(sample_meta: pd.DataFrame, blocks: list[pd.DataFrame]) -> pd.DataFrame:
    out = sample_meta.copy()
    for block in blocks:
        if block is None or block.empty:
            continue
        aligned = block.copy()
        aligned.index = aligned.index.astype(str)
        out = out.join(aligned, how="left")
    return out.fillna(0.0)


def split_feature_metadata(feature_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_cols = [col for col in ["sample", "dataset", "study", "batch", "site_label", "site_axis", "condition", "n_cells", "n_cell_types"] if col in feature_table.columns]
    meta = feature_table.loc[:, meta_cols].copy()
    features = feature_table.drop(columns=meta_cols).copy()
    return meta, features
