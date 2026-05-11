#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stromal branch-wise scArches / scANVI query mapping v1.0
========================================================

Purpose
-------
Map stromal query subsets against the *matching branch-specific* references:
  - endothelial  -> endothelial reference
  - fibroblast   -> fibroblast reference
  - smc          -> smooth muscle / pericyte reference

Why this exists
---------------
The old `stromal_scarches_query_20260315_v2_1.py` used a single global stromal
reference and produced one global `X_scanvi`. That latent is not compatible with
branch-specific scHPL trees trained from separate branch references. Feeding the
old global query latent into new branchwise scHPL greatly inflates rejection.

This script fixes that by:
  1. loading one query subset per modeled branch,
  2. computing the same covariates used by the branchwise reference training,
  3. aligning genes to the matching branch HVG space,
  4. running scArches / scANVI query mapping *per branch*,
  5. saving one query h5ad per branch plus a root manifest for downstream scHPL.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import random
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

for _k in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(_k, "8")

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from anndata import AnnData
from scipy import sparse

PIPELINE_VERSION = "v1.0-branchwise-query"
PIPELINE_DATE = "2026-04-07"
PIPELINE_TAG = "v1_0_branchwise_query"
UNLABELED = "Unknown"
RANDOM_SEED = 42

BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
QRY_LABEL_KEY = "cell_type_scarches_pred"
QRY_FINAL_KEY = "cell_type_scarches_final"
QRY_CONF_KEY = "scarches_confidence"
QRY_MARGIN_KEY = "scarches_margin"
QRY_LATENT_KEY = "X_scanvi"

QUERY_DIR = Path("/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets")
QUERY_CONFIGS: Dict[str, Dict[str, str]] = {
    "fibroblast": {
        "h5ad": str(QUERY_DIR / "fibroblast_cells.h5ad"),
        "expected_l2": "Fibroblast",
    },
    "endothelial": {
        "h5ad": str(QUERY_DIR / "endothelial_cells.h5ad"),
        "expected_l2": "Endothelial",
    },
    "smc": {
        "h5ad": str(QUERY_DIR / "smc_cells.h5ad"),
        "expected_l2": "Smooth_Muscle",
    },
}

L3_TO_L2 = {
    "Endothelia_Lymphatic": "Endothelial",
    "Endothelia_vascular_Cap_a": "Endothelial",
    "Endothelia_vascular_Cap_g": "Endothelial",
    "Endothelia_vascular_arterial_pulmonary": "Endothelial",
    "Endothelia_vascular_arterial_systemic": "Endothelial",
    "Endothelia_vascular_venous_pulmonary": "Endothelial",
    "Endothelia_vascular_venous_systemic": "Endothelial",
    "Fibro_adventitial": "Fibroblast",
    "Fibro_alveolar": "Fibroblast",
    "Fibro_myofibroblast": "Fibroblast",
    "Fibro_peribronchial": "Fibroblast",
    "Fibro_stress_activated": "Fibroblast",
    "Muscle_pericyte_pulmonary": "Pericyte",
    "Muscle_pericyte_systemic": "Pericyte",
    "Muscle_perivascular_immune_recruiting": "Smooth_Muscle",
    "Muscle_smooth_arterial_systemic": "Smooth_Muscle",
    "Muscle_smooth_pulmonary": "Smooth_Muscle",
    "Schwann_nonmyelinating": "Schwann",
}

CONTINUOUS_COVARIATES = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
STRESS_SIGNATURE_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1",
    "CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6",
    "CXXC1","DCP1A","DCP2","DCTN1","DDIT4","DDX10","DKC1","DNAJA4","DNAJB9",
    "DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1","EIF4A1","EIF4A2","EIF4A3",
    "EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A","EXOC2","EXOSC1","EXOSC10",
    "EXOSC2","EXOSC4","EXOSC5","EXOSC9","FKBP14","FUS","GEMIN4","GOSR2","H2AX",
    "HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IARS1","IFIT1","IGFBP1","IMP3",
    "KDELR3","KHSRP","KIF5B","LSM1","LSM4","MTHFD2","NFYA","NFYB","NHP2","NOLC1",
    "NOP14","NOP56","NPM1","NABP1","PAIP1","PARN","PDIA5","PDIA6","POP4","PREB",
    "PSAT1","RPS14","RRP9","SDAD1","SEC11A","SEC31A","SERP1","SHC1","MTREX",
    "SLC1A4","SLC30A5","SLC7A5","SPCS1","SPCS3","SRPRA","SRPRB","SSR1","STC2",
    "TARS1","TATDN2","TSPYL2","SKIC3","TUBB2A","VEGFA","WFS1","WIPI1","XBP1",
    "XPOT","YIF1A","YWHAZ","ZBTB17",
]
S_GENES = [
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
    "CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP",
    "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
    "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
    "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8",
]
G2M_GENES = [
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
    "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
    "HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
    "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
    "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA",
]


@dataclass
class PipelineConfig:
    reference_manifest: Path
    output_dir: Path
    confidence_threshold: float = 0.5
    max_epochs_scvi: int = 200
    max_epochs_scanvi: int = 100
    batch_size: int = 256
    use_gpu: bool = torch.cuda.is_available()
    skip_figures: bool = True
    limit_cells_per_branch: int | None = None
    branches: Tuple[str, ...] = ("endothelial", "fibroblast", "smc")


def set_seed(seed: int = RANDOM_SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    scvi.settings.seed = seed
    scvi.settings.dl_num_workers = 0
    sc.settings.verbosity = 2
    sc.settings.n_jobs = 16


def log_header(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def _ensure_category_value(series: pd.Series, value: str) -> pd.Series:
    if not pd.api.types.is_categorical_dtype(series):
        series = series.astype("category")
    if value not in series.cat.categories:
        series = series.cat.add_categories([value])
    return series


def _sanitize_object_cols(df: pd.DataFrame) -> None:
    for col in df.columns:
        s = df[col]
        if s.dtype != object:
            continue
        non_na = s.dropna()
        if len(non_na) == 0:
            df[col] = s.fillna("").astype(str)
            continue
        if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
            df[col] = s.fillna(False).astype(np.int8)
            continue
        if non_na.map(lambda x: not isinstance(x, str)).any():
            df[col] = s.map(lambda x: "" if pd.isna(x) else str(x))


def load_reference_manifest(path: Path) -> Dict[str, object]:
    if not path.exists():
        raise FileNotFoundError(f"[ERROR] Reference manifest not found: {path}")
    with open(path) as f:
        manifest = json.load(f)
    if "branches" not in manifest:
        raise KeyError("[ERROR] Reference manifest missing 'branches'")
    return manifest


def load_query_subset(branch_name: str, cfg: Dict[str, str], limit_cells: int | None = None) -> AnnData:
    query_path = Path(cfg["h5ad"])
    if not query_path.exists():
        raise FileNotFoundError(f"[ERROR] Query subset not found: {query_path}")

    ad = sc.read_h5ad(query_path)
    if not ad.var_names.is_unique:
        ad.var_names_make_unique()

    if "counts" not in ad.layers:
        if ad.raw is not None:
            ad.layers["counts"] = sparse.csr_matrix(ad.raw[:, ad.var_names].X).astype(np.float32)
        else:
            ad.layers["counts"] = sparse.csr_matrix(ad.X).astype(np.float32)
    else:
        ad.layers["counts"] = sparse.csr_matrix(ad.layers["counts"]).astype(np.float32)

    ad.X = ad.layers["counts"]

    if BATCH_KEY not in ad.obs.columns:
        candidates = [c for c in ad.obs.columns if any(k in c.lower() for k in ("sample", "batch", "donor"))]
        if candidates:
            ad.obs[BATCH_KEY] = ad.obs[candidates[0]].astype(str)
        else:
            ad.obs[BATCH_KEY] = f"query_{branch_name}"
    ad.obs[BATCH_KEY] = "qry_" + branch_name + "_" + ad.obs[BATCH_KEY].astype(str)

    if TISSUE_KEY not in ad.obs.columns:
        ad.obs[TISSUE_KEY] = f"query_{branch_name}"

    ad.obs["query_subset"] = branch_name
    ad.obs["expected_L2"] = cfg["expected_l2"]

    if limit_cells is not None and ad.n_obs > limit_cells:
        rng = np.random.default_rng(RANDOM_SEED)
        keep_idx = rng.choice(ad.n_obs, size=limit_cells, replace=False)
        ad = ad[keep_idx].copy()
        print(f"  [INFO] limited to {ad.n_obs:,} cells for branch '{branch_name}'")

    for col in [BATCH_KEY, TISSUE_KEY, "query_subset", "expected_L2"]:
        ad.obs[col] = ad.obs[col].astype(str)

    return ad


def compute_covariates(adata: AnnData, branch_name: str) -> None:
    print(f"  Computing covariates for {branch_name}...")
    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
        adata.layers["counts"] = counts
    else:
        counts = counts.tocsr()
        adata.layers["counts"] = counts

    gene_names = pd.Index(adata.var_names.astype(str))

    mt_mask = gene_names.str.upper().str.startswith("MT-") | gene_names.str.upper().str.startswith("MT.")
    total_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    if mt_mask.any():
        mt_counts = np.asarray(counts[:, np.where(mt_mask)[0]].sum(axis=1)).ravel().astype(np.float64)
        pct_counts_mt = np.zeros_like(total_counts, dtype=np.float32)
        nz = total_counts > 0
        pct_counts_mt[nz] = (mt_counts[nz] / total_counts[nz] * 100.0).astype(np.float32)
    else:
        pct_counts_mt = np.zeros(adata.n_obs, dtype=np.float32)
    adata.obs["pct_counts_mt"] = pct_counts_mt

    stress_genes_present = [g for g in STRESS_SIGNATURE_GENES if g in gene_names]
    if len(stress_genes_present) >= 10:
        adata_tmp = AnnData(X=counts.copy(), var=adata.var.copy())
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        stress_idx = gene_names.get_indexer(stress_genes_present)
        stress_expr = adata_tmp.X[:, stress_idx]
        if sparse.issparse(stress_expr):
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
        else:
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
        del stress_expr, adata_tmp
    else:
        stress_score = np.zeros(adata.n_obs, dtype=np.float32)
    adata.obs["stress_score"] = stress_score

    s_genes_present = [g for g in S_GENES if g in gene_names]
    g2m_genes_present = [g for g in G2M_GENES if g in gene_names]
    if len(s_genes_present) >= 10 and len(g2m_genes_present) >= 10:
        adata_tmp = AnnData(X=counts.copy(), obs=adata.obs[[BATCH_KEY]].copy(), var=adata.var.copy())
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        sc.tl.score_genes_cell_cycle(adata_tmp, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
        adata.obs["S_score"] = adata_tmp.obs["S_score"].to_numpy(dtype=np.float32)
        adata.obs["G2M_score"] = adata_tmp.obs["G2M_score"].to_numpy(dtype=np.float32)
        adata.obs["phase"] = adata_tmp.obs["phase"].astype(str).to_numpy()
        del adata_tmp
    else:
        adata.obs["S_score"] = np.zeros(adata.n_obs, dtype=np.float32)
        adata.obs["G2M_score"] = np.zeros(adata.n_obs, dtype=np.float32)
        adata.obs["phase"] = np.repeat("G1", adata.n_obs)

    gc.collect()


def get_reference_var_names(branch_info: Dict[str, object]) -> List[str]:
    scanvi_model_dir = Path(branch_info["scanvi_model_dir"])
    scvi_model_dir = Path(branch_info["scvi_model_dir"])
    var_file_candidates = [
        scanvi_model_dir / "SCANVI_var_names.csv",
        scvi_model_dir / "var_names.csv",
    ]
    for path in var_file_candidates:
        if path.exists():
            return pd.read_csv(path, header=None)[0].astype(str).tolist()

    ref_h5ad = Path(branch_info["reference_h5ad"])
    ad_ref = sc.read_h5ad(ref_h5ad)
    try:
        return ad_ref.var_names.astype(str).tolist()
    finally:
        del ad_ref
        gc.collect()


def align_to_reference_hvg(adata_full: AnnData, ref_var_names: List[str]) -> AnnData:
    ref_set = set(ref_var_names)
    qry_genes = set(adata_full.var_names.astype(str))
    overlap = qry_genes & ref_set
    if len(overlap) == 0:
        raise ValueError("[ERROR] Query and reference HVG spaces have 0 overlapping genes")

    gene_to_qry_col = {g: i for i, g in enumerate(adata_full.var_names.astype(str))}
    present_ref_idx = []
    present_qry_idx = []
    for j, g in enumerate(ref_var_names):
        if g in gene_to_qry_col:
            present_ref_idx.append(j)
            present_qry_idx.append(gene_to_qry_col[g])

    present_ref_idx = np.asarray(present_ref_idx, dtype=np.int32)
    present_qry_idx = np.asarray(present_qry_idx, dtype=np.int32)

    src_sub = sparse.coo_matrix(adata_full.layers["counts"][:, present_qry_idx])
    aligned_counts = sparse.csr_matrix(
        (src_sub.data.astype(np.float32), (src_sub.row, present_ref_idx[src_sub.col])),
        shape=(adata_full.n_obs, len(ref_var_names)),
    )

    var_df = pd.DataFrame(index=pd.Index(ref_var_names, dtype=str))
    var_df["in_query"] = var_df.index.isin(qry_genes)

    adata_hvg = AnnData(X=aligned_counts, obs=adata_full.obs.copy(), var=var_df)
    adata_hvg.layers["counts"] = aligned_counts.copy()
    adata_hvg.raw = AnnData(
        X=adata_full.layers["counts"],
        obs=adata_full.obs.copy(),
        var=adata_full.var.copy(),
    )
    return adata_hvg


def is_cuda_training_error(exc: Exception) -> bool:
    msg = str(exc).lower()
    return any(token in msg for token in ["cuda", "cudnn", "nccl", "device-side", "acceleratorerror"])


def cleanup_cuda_cache() -> None:
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


def train_query_scvi(adata: AnnData, ref_scvi_dir: Path, cfg: PipelineConfig, branch_out_dir: Path) -> Path:
    log_header(f"SCVI QUERY MAPPING: {branch_out_dir.name}")
    scvi.model.SCVI.prepare_query_anndata(adata, str(ref_scvi_dir))
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"]

    qry_scvi = scvi.model.SCVI.load_query_data(adata, str(ref_scvi_dir))
    model_dir = branch_out_dir / "models" / f"query_scvi_{branch_out_dir.name}_{PIPELINE_TAG}"
    model_dir.parent.mkdir(parents=True, exist_ok=True)

    accelerator = "gpu" if cfg.use_gpu and torch.cuda.is_available() else "cpu"
    t0 = time.time()
    qry_scvi.train(
        max_epochs=cfg.max_epochs_scvi,
        batch_size=cfg.batch_size,
        early_stopping=True,
        early_stopping_patience=20,
        train_size=0.9,
        accelerator=accelerator,
        devices=1,
        plan_kwargs={"weight_decay": 0.0, "lr": 5e-4},
    )
    print(f"  [OK] scVI query fine-tune done in {(time.time() - t0) / 60:.1f} min")
    qry_scvi.save(str(model_dir), overwrite=True)
    adata.obsm["X_scvi"] = qry_scvi.get_latent_representation()
    del qry_scvi
    cleanup_cuda_cache()
    return model_dir


def train_query_scanvi(
    adata: AnnData,
    ref_scanvi_dir: Path,
    cfg: PipelineConfig,
    branch_out_dir: Path,
    label_key: str,
) -> Path:
    adata.obs[label_key] = pd.Categorical(np.repeat(UNLABELED, adata.n_obs))
    adata.obs[label_key] = _ensure_category_value(adata.obs[label_key], UNLABELED)

    scvi.model.SCANVI.prepare_query_anndata(adata, str(ref_scanvi_dir))
    adata.obs[label_key] = _ensure_category_value(adata.obs[label_key], UNLABELED)

    qry_scanvi = scvi.model.SCANVI.load_query_data(adata, str(ref_scanvi_dir))
    model_dir = branch_out_dir / "models" / f"query_scanvi_{branch_out_dir.name}_{PIPELINE_TAG}"
    model_dir.parent.mkdir(parents=True, exist_ok=True)

    accelerator = "gpu" if cfg.use_gpu and torch.cuda.is_available() else "cpu"
    train_kwargs = dict(
        max_epochs=cfg.max_epochs_scanvi,
        batch_size=cfg.batch_size,
        early_stopping=True,
        early_stopping_patience=15,
        train_size=0.9,
        devices=1,
        plan_kwargs={"weight_decay": 0.0, "lr": 5e-4},
    )

    try:
        t0 = time.time()
        qry_scanvi.train(accelerator=accelerator, **train_kwargs)
        print(f"  [OK] scANVI query fine-tune done in {(time.time() - t0) / 60:.1f} min")
    except Exception as exc:
        if accelerator != "gpu" or not is_cuda_training_error(exc):
            raise
        print(f"  [WARN] scANVI GPU query training failed; retrying on CPU. Error: {exc}")
        del qry_scanvi
        cleanup_cuda_cache()
        qry_scanvi = scvi.model.SCANVI.load_query_data(adata, str(ref_scanvi_dir))
        t0 = time.time()
        qry_scanvi.train(accelerator="cpu", **train_kwargs)
        print(f"  [OK] scANVI CPU fallback done in {(time.time() - t0) / 60:.1f} min")

    qry_scanvi.save(str(model_dir), overwrite=True)
    adata.obsm[QRY_LATENT_KEY] = qry_scanvi.get_latent_representation()

    pred_labels = qry_scanvi.predict()
    if hasattr(pred_labels, "reindex"):
        pred_labels = pred_labels.reindex(adata.obs_names)
    adata.obs[QRY_LABEL_KEY] = pred_labels.astype(str)

    soft_df = qry_scanvi.predict(soft=True)
    if hasattr(soft_df, "reindex"):
        soft_df = soft_df.reindex(adata.obs_names)
    soft_arr = soft_df.to_numpy() if hasattr(soft_df, "to_numpy") else np.asarray(soft_df)
    confidence = soft_arr.max(axis=1).astype(np.float32)
    if soft_arr.shape[1] >= 2:
        sorted_soft = np.sort(soft_arr, axis=1)[:, ::-1]
        margin = (sorted_soft[:, 0] - sorted_soft[:, 1]).astype(np.float32)
    else:
        margin = confidence.copy()

    adata.obs[QRY_CONF_KEY] = confidence
    adata.obs[QRY_MARGIN_KEY] = margin
    adata.obs[QRY_FINAL_KEY] = np.where(confidence >= cfg.confidence_threshold, adata.obs[QRY_LABEL_KEY], UNLABELED)
    adata.obs["predicted_L2"] = adata.obs[QRY_LABEL_KEY].map(L3_TO_L2).fillna("Unknown")
    adata.uns["scarches_label_classes"] = soft_df.columns.tolist() if hasattr(soft_df, "columns") else []

    del qry_scanvi
    cleanup_cuda_cache()
    return model_dir


def save_branch_query_outputs(
    adata: AnnData,
    branch_name: str,
    branch_out_dir: Path,
    branch_info: Dict[str, object],
    scvi_query_dir: Path,
    scanvi_query_dir: Path,
    cfg: PipelineConfig,
) -> Dict[str, object]:
    label_summary = adata.obs[[
        BATCH_KEY,
        TISSUE_KEY,
        "query_subset",
        "expected_L2",
        "pct_counts_mt",
        "stress_score",
        "S_score",
        "G2M_score",
        "phase",
        QRY_LABEL_KEY,
        "predicted_L2",
        QRY_CONF_KEY,
        QRY_MARGIN_KEY,
        QRY_FINAL_KEY,
    ]].copy()
    label_summary.to_csv(branch_out_dir / f"{branch_name}_label_summary.csv")

    _sanitize_object_cols(adata.obs)
    _sanitize_object_cols(adata.var)
    if adata.raw is not None:
        _sanitize_object_cols(adata.raw.var)

    out_h5ad = branch_out_dir / f"adata_{branch_name}_query_{PIPELINE_TAG}.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip", compression_opts=9)

    accepted_mask = adata.obs[QRY_FINAL_KEY].astype(str) != UNLABELED
    metrics = {
        "branch": branch_name,
        "query_h5ad": str(out_h5ad),
        "query_subset": branch_name,
        "expected_l2": QUERY_CONFIGS[branch_name]["expected_l2"],
        "n_cells": int(adata.n_obs),
        "n_hvg": int(adata.n_vars),
        "n_raw_genes": int(adata.raw.n_vars if adata.raw is not None else adata.n_vars),
        "n_accepted": int(accepted_mask.sum()),
        "n_rejected": int((~accepted_mask).sum()),
        "accepted_pct": float(accepted_mask.mean() * 100.0),
        "confidence_median": float(np.median(adata.obs[QRY_CONF_KEY].to_numpy(dtype=np.float32))),
        "reference_h5ad": str(branch_info["reference_h5ad"]),
        "reference_scvi_model_dir": str(branch_info["scvi_model_dir"]),
        "reference_scanvi_model_dir": str(branch_info["scanvi_model_dir"]),
        "query_scvi_model_dir": str(scvi_query_dir),
        "query_scanvi_model_dir": str(scanvi_query_dir),
        "label_key": QRY_LABEL_KEY,
        "final_label_key": QRY_FINAL_KEY,
        "confidence_key": QRY_CONF_KEY,
        "latent_key": QRY_LATENT_KEY,
        "continuous_covariate_keys": list(CONTINUOUS_COVARIATES),
        "confidence_threshold": float(cfg.confidence_threshold),
    }
    print(json.dumps(metrics, indent=2, ensure_ascii=False))
    return metrics


def run_pipeline(cfg: PipelineConfig) -> Dict[str, object]:
    set_seed(RANDOM_SEED)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    ref_manifest = load_reference_manifest(cfg.reference_manifest)
    branch_results: Dict[str, Dict[str, object]] = {}

    for branch_name in cfg.branches:
        log_header(f"BRANCH QUERY MAPPING: {branch_name}")
        if branch_name not in QUERY_CONFIGS:
            raise KeyError(f"[ERROR] Unknown query branch: {branch_name}")
        if branch_name not in ref_manifest["branches"]:
            raise KeyError(f"[ERROR] Branch '{branch_name}' missing from reference manifest")

        branch_info = ref_manifest["branches"][branch_name]
        if not branch_info.get("scanvi_model_dir"):
            raise ValueError(f"[ERROR] Branch '{branch_name}' reference lacks scanvi_model_dir")

        adata_query = load_query_subset(
            branch_name=branch_name,
            cfg=QUERY_CONFIGS[branch_name],
            limit_cells=cfg.limit_cells_per_branch,
        )
        print(f"  Query shape: {adata_query.n_obs:,} x {adata_query.n_vars:,}")

        compute_covariates(adata_query, branch_name)
        ref_var_names = get_reference_var_names(branch_info)
        adata_hvg = align_to_reference_hvg(adata_query, ref_var_names)
        print(f"  Aligned to reference HVG: {adata_hvg.n_obs:,} x {adata_hvg.n_vars:,}")

        branch_out_dir = cfg.output_dir / branch_name
        branch_out_dir.mkdir(parents=True, exist_ok=True)

        scvi_query_dir = train_query_scvi(
            adata=adata_hvg,
            ref_scvi_dir=Path(branch_info["scvi_model_dir"]),
            cfg=cfg,
            branch_out_dir=branch_out_dir,
        )
        scanvi_query_dir = train_query_scanvi(
            adata=adata_hvg,
            ref_scanvi_dir=Path(branch_info["scanvi_model_dir"]),
            cfg=cfg,
            branch_out_dir=branch_out_dir,
            label_key=branch_info.get("label_key", "scanvi_label"),
        )

        branch_results[branch_name] = save_branch_query_outputs(
            adata=adata_hvg,
            branch_name=branch_name,
            branch_out_dir=branch_out_dir,
            branch_info=branch_info,
            scvi_query_dir=scvi_query_dir,
            scanvi_query_dir=scanvi_query_dir,
            cfg=cfg,
        )

        del adata_query, adata_hvg
        gc.collect()

    summary_df = pd.DataFrame(branch_results.values()).sort_values("branch")
    summary_df.to_csv(cfg.output_dir / "branch_query_summary.csv", index=False)

    manifest = {
        "pipeline_version": PIPELINE_VERSION,
        "pipeline_date": PIPELINE_DATE,
        "pipeline_tag": PIPELINE_TAG,
        "reference_manifest": str(cfg.reference_manifest),
        "output_dir": str(cfg.output_dir),
        "confidence_threshold": float(cfg.confidence_threshold),
        "continuous_covariate_keys": list(CONTINUOUS_COVARIATES),
        "branches": branch_results,
    }
    manifest_path = cfg.output_dir / "branch_query_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[OK] Query branch manifest saved: {manifest_path}")
    return manifest


def parse_args() -> PipelineConfig:
    parser = argparse.ArgumentParser(description="Branch-wise stromal scArches query mapping")
    parser.add_argument(
        "--reference-manifest",
        type=Path,
        default=Path("/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/branch_reference_manifest.json"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/home/h2048/data/py/0407/stromal_scarches_query_branchwise_v1_0"),
    )
    parser.add_argument("--confidence-threshold", type=float, default=0.5)
    parser.add_argument("--max-epochs-scvi", type=int, default=200)
    parser.add_argument("--max-epochs-scanvi", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--limit-cells-per-branch", type=int, default=None)
    parser.add_argument(
        "--branches",
        nargs="+",
        default=["endothelial", "fibroblast", "smc"],
        choices=["endothelial", "fibroblast", "smc"],
    )
    args = parser.parse_args()

    return PipelineConfig(
        reference_manifest=args.reference_manifest,
        output_dir=args.output_dir,
        confidence_threshold=args.confidence_threshold,
        max_epochs_scvi=args.max_epochs_scvi,
        max_epochs_scanvi=args.max_epochs_scanvi,
        batch_size=args.batch_size,
        use_gpu=(not args.cpu_only) and torch.cuda.is_available(),
        limit_cells_per_branch=args.limit_cells_per_branch,
        branches=tuple(args.branches),
    )


def main() -> None:
    cfg = parse_args()
    start = time.time()
    manifest = run_pipeline(cfg)
    elapsed = (time.time() - start) / 60
    print("\n" + "=" * 80)
    print(f"PIPELINE COMPLETE ({elapsed:.1f} min)")
    print("=" * 80)
    for branch_name, info in manifest["branches"].items():
        print(
            f"  {branch_name:<12} cells={info['n_cells']:,} "
            f"accepted={info['n_accepted']:,} rejected={info['n_rejected']:,} "
            f"query_h5ad={info['query_h5ad']}"
        )


if __name__ == "__main__":
    main()
