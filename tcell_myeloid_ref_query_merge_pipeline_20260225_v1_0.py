#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
T Cell + Myeloid Reference-Query Merge Pipeline v1.0
==============================================================================

Purpose:
--------
Merge reference and query data for T cells and Myeloid cells, then run
complete scVI + CellTypist + scANVI pipeline.

Key Design:
-----------
1. Direct merge (NOT scArches mapping) - allows discovery of new cell types
2. Reference cells keep their labels
3. Query cells marked as "Unknown" for scANVI to reclassify
4. Full gene preservation in .raw (for downstream DE analysis)
5. HVG-based training (memory efficient)

Input:
------
- Reference h5ad: with cell type annotations
- Query h5ad: raw counts, may contain new cell types

Output:
-------
- Merged h5ad with scVI/scANVI latent spaces
- CellTypist predictions
- scANVI predictions with confidence scores
- UMAP visualizations
- Saved UMAP operator for future projection

Based on: merged_scanvi_training_20260208_v2_5_5.py (HOTFIX)
Author: Claude Code
Date: 2026-02-25
Version: 1.0
==============================================================================
"""

import sys
import warnings
import json
import gc
import joblib
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix, vstack as sp_vstack
from scipy import sparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
import scvi
import celltypist
from celltypist import models
from umap import UMAP

import torch

# ===== P1-FIX: Thread Control (HPC Stability) =====
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# Silence warnings
warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None

print("=" * 80)
print("T Cell + Myeloid Ref-Query Merge Pipeline v1.0")
print("=" * 80)

# Check GPU
gpu_available = torch.cuda.is_available()
print(f"GPU available: {gpu_available}")
if gpu_available:
    print(f"GPU Device: {torch.cuda.get_device_name(0)}")

# ===== P1-FIX: scVI Stability Settings =====
scvi.settings.dl_num_workers = 0
print(f"scVI dl_num_workers: {scvi.settings.dl_num_workers}")

# ==============================================================================
# 1. CONFIGURATION - EDIT THESE PATHS
# ==============================================================================

# --- Input Files ---
REFERENCE_H5AD = "/path/to/your/reference_tcell_myeloid.h5ad"  # EDIT THIS
QUERY_H5AD = "/path/to/your/query_tcell_myeloid.h5ad"          # EDIT THIS

# --- Column Configuration ---
REF_LABEL_COLUMN = "cell_type"  # Reference annotation column
BATCH_KEY = "sample"            # Batch/sample column
TISSUE_KEY = "tissue"           # Tissue column

# --- Output Configuration ---
OUTPUT_DIR = "/home/h2048/data/py/20260225_tcell_myeloid_merged_pipeline"
OUTPUT_PREFIX = "tcell_myeloid_merged"

# --- Cell Type Selection ---
# Specify which cell types to include from reference (None = all)
INCLUDE_CELL_TYPES = None  # e.g., ["CD4_T", "CD8_T", "Monocyte", "Macrophage"]

# --- Pipeline Parameters ---
N_HVG = 4000
FORCE_MARKERS_IN_HVG = True

# scVI Parameters
SCVI_N_LATENT = 50
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_DROPOUT = 0.1
MAX_EPOCHS_SCVI = 400

# scANVI Parameters
MAX_EPOCHS_SCANVI = 200
UNLABELED_CATEGORY = "Unknown"

# Training Parameters
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 0.0

# CellTypist Parameters
CELLTYPIST_MODEL = "Human_Lung_Atlas.pkl"
CELLTYPIST_MAJORITY_VOTE = True

# --- Gene Signatures ---
TCELL_MARKERS = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B",
    "CCR7", "SELL", "LEF1", "TCF7",           # Naive
    "IL7R", "GZMK", "GZMB", "PRF1", "GNLY", "NKG7",  # Effector
    "FOXP3", "IL2RA", "CTLA4",                 # Treg
    "GATA3", "IL13", "IL5",                    # Th2
    "TBX21", "IFNG",                           # Th1
    "RORC", "IL17A", "IL17F",                  # Th17
]

MYELOID_MARKERS = [
    "LYZ", "CD14", "FCGR3A", "MS4A7",         # Monocytes
    "CD68", "CD163", "MRC1", "MARCO",         # Macrophages
    "FCGR1A", "S100A8", "S100A9", "FCER1A",   # DCs
    "CLEC9A", "XCR1", "CLEC10A",              # cDC1, cDC2
    "LILRA4", "CLEC4C",                        # pDC
]

STRESS_SIGNATURE_GENES = [
    "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "IER2"
]

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UHRF1",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "HELLS", "RFC2", "RPA2",
    "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1"
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB1", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E"
]

FORCED_MARKERS = list(set(TCELL_MARKERS + MYELOID_MARKERS))

# --- Reproducibility ---
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

def ensure_counts_layer(adata, counts_layer="counts"):
    """
    P1-FIX: Robust counts validation with float32 tolerance.
    """
    if counts_layer not in (adata.layers or {}):
        raise ValueError(
            f"CRITICAL ERROR: layers['{counts_layer}'] not found. "
            "scVI requires raw counts!"
        )

    # Sample validation
    X_counts = adata.layers[counts_layer]
    if issparse(X_counts):
        sample_data = X_counts.data[:1000]
    else:
        sample_data = X_counts.flat[:1000]

    sample = np.asarray(sample_data, dtype=np.float64)

    # Check 1: Non-negative
    if np.any(sample < 0):
        raise ValueError("counts contains negative values!")

    # Check 2: Near-integer (P1-FIX: tolerance for float32 storage)
    if not np.allclose(sample, np.round(sample), atol=1e-6):
        raise ValueError(
            "counts looks non-integer (possible normalized/log data). "
            f"Sample range: [{sample.min():.4f}, {sample.max():.4f}]"
        )

    # Ensure CSR for speed
    if issparse(X_counts) and not isinstance(X_counts, csr_matrix):
        adata.layers[counts_layer] = csr_matrix(X_counts)

    return counts_layer


def ensure_batch_tissue(adata, batch_key, tissue_key):
    """Ensure batch and tissue columns exist and are categorical."""
    if batch_key not in adata.obs.columns:
        print(f"  WARNING: {batch_key} not found, creating placeholder")
        adata.obs[batch_key] = "unknown_batch"
    adata.obs[batch_key] = adata.obs[batch_key].astype("category")

    if tissue_key not in adata.obs.columns:
        print(f"  WARNING: {tissue_key} not found, creating placeholder")
        adata.obs[tissue_key] = "unknown_tissue"
    adata.obs[tissue_key] = adata.obs[tissue_key].astype("category")
    if "unknown_tissue" not in adata.obs[tissue_key].cat.categories:
        adata.obs[tissue_key] = adata.obs[tissue_key].cat.add_categories(["unknown_tissue"])
    adata.obs[tissue_key] = adata.obs[tissue_key].fillna("unknown_tissue")


def ensure_pct_mt(adata, counts_layer="counts"):
    """Ensure MT percentage is computed."""
    for col in ["pct_counts_mt", "percent.mt", "percent_mt"]:
        if col in adata.obs.columns:
            if col != "pct_counts_mt":
                adata.obs["pct_counts_mt"] = adata.obs[col]
            return
    # Compute
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, layer=counts_layer)


def ensure_stress_score(adata, stress_genes, counts_layer="counts"):
    """Compute stress signature score."""
    if "stress_score" in adata.obs.columns:
        return

    genes = [g for g in stress_genes if g in adata.var_names]
    if len(genes) < 5:
        adata.obs["stress_score"] = 0.0
        return

    X = adata.layers[counts_layer]
    idx = adata.var_names.get_indexer(genes)
    idx = idx[idx >= 0]
    sub = X[:, idx]

    if issparse(X):
        tot = np.asarray(X.sum(axis=1)).ravel()
    else:
        tot = X.sum(axis=1)
    tot[tot == 0] = 1.0

    if issparse(sub):
        scale = sparse.diags(1e4 / tot)
        sub_norm = scale @ sub
        sub_norm = sub_norm.tocsr()
        sub_norm.data = np.log1p(sub_norm.data)
        expr = np.asarray(sub_norm.mean(axis=1)).ravel()
    else:
        sub_norm = (sub / tot[:, None]) * 1e4
        expr = np.log1p(sub_norm).mean(axis=1)

    mn, mx = float(expr.min()), float(expr.max())
    adata.obs["stress_score"] = (expr - mn) / (mx - mn) if mx > mn else 0.0


def ensure_cell_cycle_scores(adata, s_genes, g2m_genes, counts_layer="counts"):
    """Compute cell cycle scores."""
    if all(k in adata.obs.columns for k in ["S_score", "G2M_score", "phase"]):
        return

    s_in = [g for g in s_genes if g in adata.var_names]
    g_in = [g for g in g2m_genes if g in adata.var_names]

    if len(s_in) < 5 or len(g_in) < 5:
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"
        return

    # Lightweight scoring
    cc_union = list(dict.fromkeys(s_in + g_in))
    cc_idx = adata.var_names.get_indexer(cc_union)
    X_cc = adata.layers[counts_layer][:, cc_idx].copy()
    var_cc = adata.var.iloc[cc_idx].copy()

    ad_tmp = sc.AnnData(X=X_cc, var=var_cc)
    sc.pp.normalize_total(ad_tmp, target_sum=1e4)
    sc.pp.log1p(ad_tmp)
    sc.tl.score_genes_cell_cycle(ad_tmp, s_genes=s_in, g2m_genes=g_in)

    adata.obs["S_score"] = ad_tmp.obs["S_score"].values
    adata.obs["G2M_score"] = ad_tmp.obs["G2M_score"].values
    adata.obs["phase"] = ad_tmp.obs["phase"].values

    del ad_tmp
    gc.collect()


def prepare_covariates(adata):
    """Prepare all covariates for scVI training."""
    print("  -> Validating counts (robust)...")
    ensure_counts_layer(adata, "counts")

    print("  -> Batch/Tissue...")
    ensure_batch_tissue(adata, BATCH_KEY, TISSUE_KEY)

    print("  -> MT% / Stress / Cell Cycle...")
    ensure_pct_mt(adata, "counts")
    ensure_stress_score(adata, STRESS_SIGNATURE_GENES, "counts")
    ensure_cell_cycle_scores(adata, S_GENES, G2M_GENES, "counts")


def run_celltypist(adata, model_name=CELLTYPIST_MODEL, majority_vote=True):
    """Run CellTypist annotation."""
    print("\n[CellTypist] Starting annotation...")

    # Download model if needed
    models.download_models(model=model_name)
    model = models.Model.load(model_name)

    # Prepare data (log-normalized)
    adata_ct = adata.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    # Predict
    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=majority_vote,
        mode='best match'
    )

    # Add results to original adata
    adata.obs["celltypist_pred"] = predictions.predicted_labels.predicted_labels
    adata.obs["celltypist_confidence"] = predictions.probability_matrix.max(axis=1).values

    if majority_vote:
        adata.obs["celltypist_majority"] = predictions.predicted_labels.majority_voting

    print(f"  -> CellTypist predictions added")
    print(f"  -> Top predictions:")
    print(adata.obs["celltypist_pred"].value_counts().head(10))

    del adata_ct
    gc.collect()

    return predictions


# ==============================================================================
# 3. MAIN PIPELINE
# ==============================================================================

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

# --- Step 1: Load Reference and Query ---
print("\n[Step 1] Loading Reference and Query...")

print(f"  Loading reference: {REFERENCE_H5AD}")
adata_ref = sc.read_h5ad(REFERENCE_H5AD)
adata_ref.var_names_make_unique()
print(f"    Reference shape: {adata_ref.shape}")

print(f"  Loading query: {QUERY_H5AD}")
adata_qry = sc.read_h5ad(QUERY_H5AD)
adata_qry.var_names_make_unique()
print(f"    Query shape: {adata_qry.shape}")

# --- Step 2: Subset Reference if needed ---
if INCLUDE_CELL_TYPES is not None:
    print(f"\n[Step 2] Subsetting reference to: {INCLUDE_CELL_TYPES}")
    if REF_LABEL_COLUMN in adata_ref.obs.columns:
        mask = adata_ref.obs[REF_LABEL_COLUMN].isin(INCLUDE_CELL_TYPES)
        adata_ref = adata_ref[mask].copy()
        print(f"    Reference after subset: {adata_ref.shape}")
    else:
        print(f"    WARNING: {REF_LABEL_COLUMN} not found, skipping subset")

# --- Step 3: Add data_source column and prepare labels ---
print("\n[Step 3] Preparing labels...")

adata_ref.obs["data_source"] = "reference"
adata_qry.obs["data_source"] = "query"

# Reference labels
if REF_LABEL_COLUMN in adata_ref.obs.columns:
    adata_ref.obs["cell_type_original"] = adata_ref.obs[REF_LABEL_COLUMN].astype(str)
    print(f"  -> Reference labels from: {REF_LABEL_COLUMN}")
    print(f"     Categories: {adata_ref.obs['cell_type_original'].nunique()}")
else:
    print(f"  WARNING: {REF_LABEL_COLUMN} not found in reference!")
    adata_ref.obs["cell_type_original"] = "reference_cell"

# Query labels = Unknown (will be reclassified by scANVI)
adata_qry.obs["cell_type_original"] = UNLABELED_CATEGORY

# --- Step 4: Find common genes and subset ---
print("\n[Step 4] Finding common genes...")

ref_genes = set(adata_ref.var_names)
qry_genes = set(adata_qry.var_names)
common_genes = list(ref_genes.intersection(qry_genes))

print(f"  Reference genes: {len(ref_genes):,}")
print(f"  Query genes: {len(qry_genes):,}")
print(f"  Common genes: {len(common_genes):,}")

if len(common_genes) < 1000:
    raise ValueError(f"Too few common genes ({len(common_genes)}). Check gene naming!")

adata_ref = adata_ref[:, common_genes].copy()
adata_qry = adata_qry[:, common_genes].copy()

# --- Step 5: Ensure counts layers ---
print("\n[Step 5] Validating counts...")
ensure_counts_layer(adata_ref, "counts")
ensure_counts_layer(adata_qry, "counts")

# --- Step 6: Concatenate ---
print("\n[Step 6] Concatenating reference and query...")

# Make obs_names unique
adata_ref.obs_names = pd.Index([f"ref_{x}" for x in adata_ref.obs_names])
adata_qry.obs_names = pd.Index([f"qry_{x}" for x in adata_qry.obs_names])

adata_merged = sc.concat(
    {"reference": adata_ref, "query": adata_qry},
    axis=0,
    join="inner",
    merge="unique",
    label="data_source"
)

print(f"  Merged shape: {adata_merged.shape}")
print(f"  Reference cells: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
print(f"  Query cells: {(adata_merged.obs['data_source'] == 'query').sum():,}")

# Clean up
del adata_ref, adata_qry
gc.collect()

# --- Step 7: Prepare covariates ---
print("\n[Step 7] Preparing covariates...")
prepare_covariates(adata_merged)

# P1-FIX: Add symbol_base for robust marker matching
print("  -> Adding symbol_base column (P1-FIX)...")
adata_merged.var["symbol_base"] = adata_merged.var_names.str.replace(r"-\d+$", "", regex=True)

# --- Step 8: HVG Selection ---
print("\n[Step 8] Selecting HVGs (Robust)...")
hvg_method = "unknown"
try:
    print("  -> Attempting: batch-aware seurat_v3")
    sc.pp.highly_variable_genes(
        adata_merged, layer="counts", n_top_genes=N_HVG,
        batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
    )
    hvg_method = "batch_seurat_v3"
except Exception as e1:
    try:
        print(f"  -> Failed ({str(e1)[:50]}). Attempting: standard seurat_v3")
        sc.pp.highly_variable_genes(
            adata_merged, layer="counts", n_top_genes=N_HVG,
            flavor="seurat_v3", subset=False
        )
        hvg_method = "standard_seurat_v3"
    except Exception as e2:
        print(f"  -> Failed ({str(e2)[:50]}). Fallback: cell_ranger")
        sc.pp.highly_variable_genes(
            adata_merged, layer="counts", n_top_genes=N_HVG,
            flavor="cell_ranger", subset=False
        )
        hvg_method = "cell_ranger"

print(f"  -> Method used: {hvg_method}")

# P1-FIX: Use symbol_base for marker matching
if FORCE_MARKERS_IN_HVG:
    n_added = 0
    marker_set = set(FORCED_MARKERS)

    for idx, symbol_base in enumerate(adata_merged.var["symbol_base"]):
        if symbol_base in marker_set:
            real_name = adata_merged.var_names[idx]
            if not adata_merged.var.loc[real_name, "highly_variable"]:
                adata_merged.var.loc[real_name, "highly_variable"] = True
                n_added += 1
                print(f"    Added: {real_name} (base: {symbol_base})")

    print(f"  -> Forced {n_added} / {len(FORCED_MARKERS)} markers into HVG")

n_hvg_final = adata_merged.var["highly_variable"].sum()
print(f"  -> Final HVG count: {n_hvg_final}")

# Save HVG list
hvg_genes = adata_merged.var_names[adata_merged.var["highly_variable"]].tolist()
with open(output_dir / f"{OUTPUT_PREFIX}_hvg_genes.txt", "w") as f:
    f.write("\n".join(hvg_genes))

# --- Step 9: Build FULL matrix for .raw (all cells, all genes) ---
print("\n[Step 9] Building full matrix for .raw...")

# Get full counts for all cells
full_counts = adata_merged.layers["counts"]
if issparse(full_counts) and not isinstance(full_counts, csr_matrix):
    full_counts = csr_matrix(full_counts)

# Save var dataframe for .raw
raw_var = adata_merged.var.copy()

print(f"  Full matrix shape for .raw: {full_counts.shape}")

# --- Step 10: Create training subset (HVG only, single-copy counts) ---
print("\n[Step 10] Creating scVI Training Data (P1-FIX: single-copy HVG)...")

hvg_mask = adata_merged.var["highly_variable"].values
X_hvg = adata_merged.layers["counts"][:, hvg_mask]

# Ensure CSR
if issparse(X_hvg) and not isinstance(X_hvg, csr_matrix):
    X_hvg = csr_matrix(X_hvg)

# P1-FIX: Single copy, layer points to X (0 extra memory cost)
adata_train = sc.AnnData(
    X=X_hvg.copy(),  # Only ONE copy needed
    obs=adata_merged.obs.copy(),
    var=adata_merged.var.iloc[hvg_mask].copy()
)
adata_train.var_names = adata_merged.var_names[hvg_mask]
adata_train.layers["counts"] = adata_train.X  # Shared reference, no .copy()

print(f"  Training data: {adata_train.shape}")
print(f"  Layers: {list(adata_train.layers.keys())}")
print(f"  Memory: counts layer is shared reference to X (0 extra cost)")

# Prepare unified labels for scANVI
adata_train.obs["scanvi_labels"] = adata_train.obs["cell_type_original"].astype(str)
adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")

# Ensure Unknown category exists
if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
    adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])

print(f"  -> Label distribution:")
print(adata_train.obs["scanvi_labels"].value_counts())

gc.collect()

# --- Step 11: CellTypist Annotation ---
print("\n[Step 11] Running CellTypist...")
try:
    run_celltypist(adata_train, model_name=CELLTYPIST_MODEL, majority_vote=CELLTYPIST_MAJORITY_VOTE)

    # Also add to full adata
    adata_merged.obs["celltypist_pred"] = adata_train.obs["celltypist_pred"]
    adata_merged.obs["celltypist_confidence"] = adata_train.obs["celltypist_confidence"]
    if CELLTYPIST_MAJORITY_VOTE:
        adata_merged.obs["celltypist_majority"] = adata_train.obs["celltypist_majority"]
except Exception as e:
    print(f"  WARNING: CellTypist failed: {e}")

# --- Step 12: scVI Training ---
print("\n[Step 12] Training scVI...")

setup_kwargs = {
    "layer": "counts",
    "batch_key": BATCH_KEY,
    "continuous_covariate_keys": ["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
    "categorical_covariate_keys": [TISSUE_KEY]
}

scvi.model.SCVI.setup_anndata(adata_train, **setup_kwargs)

scvi_model = scvi.model.SCVI(
    adata_train,
    n_latent=SCVI_N_LATENT,
    n_layers=SCVI_N_LAYERS,
    n_hidden=SCVI_N_HIDDEN,
    dropout_rate=SCVI_DROPOUT
)

train_kwargs = {
    "max_epochs": MAX_EPOCHS_SCVI,
    "batch_size": BATCH_SIZE,
    "early_stopping": True,
    "early_stopping_patience": 30,
    "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
}
if gpu_available:
    train_kwargs["accelerator"] = "gpu"
    train_kwargs["devices"] = 1

scvi_model.train(**train_kwargs)
print("  -> scVI training complete")

# --- Step 13: scANVI Training ---
print("\n[Step 13] Training scANVI...")

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata_train,
    labels_key="scanvi_labels",
    unlabeled_category=UNLABELED_CATEGORY
)

scanvi_train_kwargs = {
    "max_epochs": MAX_EPOCHS_SCANVI,
    "batch_size": BATCH_SIZE,
    "early_stopping": True,
    "early_stopping_patience": 20,
    "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
}
if gpu_available:
    scanvi_train_kwargs["accelerator"] = "gpu"
    scanvi_train_kwargs["devices"] = 1

scanvi_model.train(**scanvi_train_kwargs)
print("  -> scANVI training complete")

# --- Step 14: Export Results (P1-FIX: Strict Validation) ---
print("\n[Step 14] Exporting Results (P1-FIX: Strict Validation)...")

# Extract latent
latent_scanvi = scanvi_model.get_latent_representation(adata_train)
print(f"  -> Latent shape: {latent_scanvi.shape}")

# Verify completeness
if latent_scanvi.shape[0] != len(adata_train):
    raise ValueError(
        f"CRITICAL: Latent dimension mismatch! "
        f"Expected {len(adata_train)}, got {latent_scanvi.shape[0]}"
    )

# Index-aligned writeback
latent_df = pd.DataFrame(
    latent_scanvi,
    index=adata_train.obs_names,
    columns=[f"scANVI_{i}" for i in range(latent_scanvi.shape[1])]
)

latent_aligned = latent_df.reindex(adata_merged.obs_names)

# P1-FIX: Hard error on missing cells (not silent fillna)
n_missing = latent_aligned.isna().any(axis=1).sum()
if n_missing > 0:
    raise ValueError(
        f"CRITICAL: {n_missing} cells missing latent representation! "
        "This should NOT happen - check data consistency."
    )

adata_merged.obsm["X_scANVI"] = latent_aligned.values

# Predictions
pred_labels = scanvi_model.predict(adata_train)
pred_df = pd.Series(pred_labels, index=adata_train.obs_names)
pred_aligned = pred_df.reindex(adata_merged.obs_names)

if pred_aligned.isna().sum() > 0:
    raise ValueError(f"CRITICAL: {pred_aligned.isna().sum()} cells missing predictions!")

adata_merged.obs["scanvi_pred"] = pred_aligned.values

# P2-ENHANCEMENT: Soft probabilities with label column names
print("  -> Extracting probabilities (P2-FIX: with label names)...")
proba = scanvi_model.predict(adata_train, soft=True).astype(np.float32)

# Get label order from model
label_order = scanvi_model.labels
proba_df = pd.DataFrame(
    proba,
    index=adata_train.obs_names,
    columns=label_order
)
proba_aligned = proba_df.reindex(adata_merged.obs_names)

if proba_aligned.isna().any().any():
    raise ValueError("CRITICAL: Missing probability values!")

adata_merged.obsm["scanvi_proba"] = proba_aligned.values
adata_merged.obs["scanvi_confidence"] = proba_aligned.values.max(axis=1)

print(f"  -> scANVI predictions:")
print(adata_merged.obs["scanvi_pred"].value_counts().head(15))

# --- Step 15: Attach .raw (FULL genes) ---
print("\n[Step 15] Attaching .raw (FULL genes)...")

from anndata import AnnData

adata_merged.raw = AnnData(
    X=full_counts,
    obs=adata_merged.obs.copy(),
    var=raw_var,
)

assert adata_merged.raw.n_obs == adata_merged.n_obs
assert adata_merged.raw.n_vars == len(raw_var)
print(f"  OK .raw attached: {adata_merged.raw.n_vars} genes")

# --- Step 16: Multiple UMAPs (scVI, scANVI, etc.) ---
print("\n[Step 16] Computing Multiple UMAPs...")

# 16a: UMAP on scVI latent
print("  -> Computing UMAP on scVI latent...")
sc.pp.neighbors(adata_merged, use_rep="X_scVI", n_neighbors=30, random_state=RANDOM_SEED)
sc.tl.umap(adata_merged, random_state=RANDOM_SEED)
adata_merged.obsm["X_umap_scVI"] = adata_merged.obsm["X_umap"].copy()
print(f"     Saved to X_umap_scVI")

# 16b: UMAP on scANVI latent (DEFAULT)
print("  -> Computing UMAP on scANVI latent (DEFAULT)...")
sc.pp.neighbors(adata_merged, use_rep="X_scANVI", n_neighbors=30, random_state=RANDOM_SEED)
sc.tl.umap(adata_merged, random_state=RANDOM_SEED)
adata_merged.obsm["X_umap_scANVI"] = adata_merged.obsm["X_umap"].copy()

# Save as DEFAULT X_umap (for sc.pl.umap default behavior)
adata_merged.obsm["X_umap"] = adata_merged.obsm["X_umap_scANVI"].copy()
print(f"     Saved to X_umap_scVI and X_umap (default)")

# 16c: Save UMAP operator for scANVI (for future projection)
umap_op_scanvi = UMAP(
    n_neighbors=30,
    n_components=2,
    min_dist=0.5,
    spread=1.0,
    metric="euclidean",
    random_state=RANDOM_SEED
)
umap_op_scanvi.fit(adata_merged.obsm["X_scANVI"])
joblib.dump(umap_op_scanvi, output_dir / f"{OUTPUT_PREFIX}_umap_scanvi_operator.joblib")
print(f"  -> UMAP operator (scANVI) saved")

# P2-ENHANCEMENT: Confidence histogram
print("\n[P2-DIAGNOSTIC] Plotting confidence distribution...")
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(adata_merged.obs["scanvi_confidence"], bins=50, edgecolor="black")
ax.set_xlabel("Prediction Confidence")
ax.set_ylabel("Cell Count")
ax.set_title("scANVI Prediction Confidence Distribution")
plt.tight_layout()
plt.savefig(output_dir / f"{OUTPUT_PREFIX}_confidence_histogram.png", dpi=150)
plt.close()
print(f"  -> Confidence histogram saved")

# --- Step 17: Save Results ---
print("\n[Step 17] Saving Results...")

# Save models
scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
scvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scvi_model", overwrite=True)

# Save config
config = {
    "version": "1.0",
    "timestamp": datetime.now().isoformat(),
    "input": {
        "reference": REFERENCE_H5AD,
        "query": QUERY_H5AD,
        "ref_label_column": REF_LABEL_COLUMN
    },
    "output_dir": str(output_dir),
    "parameters": {
        "n_hvg": int(n_hvg_final),
        "hvg_method": hvg_method,
        "scvi_n_latent": SCVI_N_LATENT,
        "scvi_n_layers": SCVI_N_LAYERS,
        "max_epochs_scvi": MAX_EPOCHS_SCVI,
        "max_epochs_scanvi": MAX_EPOCHS_SCANVI
    },
    "cell_counts": {
        "reference": int((adata_merged.obs["data_source"] == "reference").sum()),
        "query": int((adata_merged.obs["data_source"] == "query").sum()),
        "total": int(adata_merged.n_obs)
    },
    "scanvi_labels": list(label_order),
    "umap_spaces": {
        "X_umap": "DEFAULT - scANVI-based UMAP",
        "X_umap_scVI": "scVI latent space UMAP",
        "X_umap_scANVI": "scANVI latent space UMAP"
    },
    "fixes_applied": {
        "P0": ["X_umap_default", "multiple_umaps", "save_umap_operator"],
        "P1": ["robust_int_check", "single_copy_counts", "symbol_base_markers",
               "hard_error_missing_latent", "dl_num_workers_0", "thread_control"],
        "P2": ["proba_label_names", "confidence_hist"]
    }
}

with open(output_dir / f"{OUTPUT_PREFIX}_config.json", "w") as f:
    json.dump(config, f, indent=2)

# Save H5AD (with .raw)
output_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"
adata_merged.write_h5ad(output_h5ad, compression="gzip")
print(f"  -> Saved: {output_h5ad}")
print(f"     - Shape: {adata_merged.shape}")
print(f"     - .raw genes: {adata_merged.raw.n_vars}")

# Save training subset
train_h5ad = output_dir / f"{OUTPUT_PREFIX}_train_HVG.h5ad"
adata_train.write_h5ad(train_h5ad, compression="gzip")
print(f"  -> Saved: {train_h5ad}")

# Clean up
del adata_train, scvi_model
gc.collect()

# --- Step 18: Visualization ---
print("\n[Step 18] Creating visualizations...")

# Figure 1: Overview using DEFAULT X_umap (scANVI-based)
fig, axes = plt.subplots(3, 3, figsize=(21, 18))

# Row 1: Data source and original labels
sc.pl.umap(adata_merged, color="data_source", ax=axes[0, 0], show=False,
           title="Data Source (scANVI UMAP)", s=10)

# Reference original labels (ref only)
adata_merged.obs["_ref_labels_only"] = pd.Series(pd.NA, index=adata_merged.obs_names, dtype="object")
mask_ref = adata_merged.obs["data_source"] == "reference"
adata_merged.obs.loc[mask_ref, "_ref_labels_only"] = adata_merged.obs.loc[mask_ref, "cell_type_original"].astype(str).values
adata_merged.obs["_ref_labels_only"] = adata_merged.obs["_ref_labels_only"].astype("category")

sc.pl.umap(adata_merged, color="_ref_labels_only", ax=axes[0, 1], show=False,
           title="Reference Original Labels", legend_loc="right margin", s=10)

# scANVI predictions
sc.pl.umap(adata_merged, color="scanvi_pred", ax=axes[0, 2], show=False,
           title="scANVI Predictions", legend_loc="right margin", s=10)

# Row 2: CellTypist results
if "celltypist_pred" in adata_merged.obs.columns:
    sc.pl.umap(adata_merged, color="celltypist_pred", ax=axes[1, 0], show=False,
               title="CellTypist Predictions", legend_loc="right margin", s=10)
    sc.pl.umap(adata_merged, color="celltypist_confidence", ax=axes[1, 1], show=False,
               title="CellTypist Confidence", cmap="viridis", vmin=0, vmax=1, s=10)
else:
    axes[1, 0].text(0.5, 0.5, "CellTypist\nN/A", ha="center", va="center", transform=axes[1, 0].transAxes)
    axes[1, 1].text(0.5, 0.5, "CellTypist\nN/A", ha="center", va="center", transform=axes[1, 1].transAxes)

sc.pl.umap(adata_merged, color="scanvi_confidence", ax=axes[1, 2], show=False,
           title="scANVI Confidence", cmap="viridis", vmin=0, vmax=1, s=10)

# Row 3: Marker genes
for idx, marker in enumerate(["CD3E", "CD4", "CD8A"]):
    if marker in adata_merged.var_names:
        sc.pl.umap(adata_merged, color=marker, ax=axes[2, idx], show=False,
                   title=f"{marker} (T Cell)", cmap="Reds", s=10, use_raw=True)
    elif marker in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color=marker, ax=axes[2, idx], show=False,
                   title=f"{marker} (T Cell)", cmap="Reds", s=10, use_raw=True)
    else:
        axes[2, idx].text(0.5, 0.5, f"{marker}\nN/A", ha="center", va="center",
                         transform=axes[2, idx].transAxes)

plt.tight_layout()
fig_path = output_dir / f"{OUTPUT_PREFIX}_overview_scANVI.pdf"
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"  -> Saved: {fig_path}")

# Figure 2: Compare scVI vs scANVI UMAPs
fig, axes = plt.subplots(2, 3, figsize=(21, 12))

# scVI UMAP
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="data_source", ax=axes[0, 0], show=False,
           title="Data Source (scVI UMAP)", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="scanvi_pred", ax=axes[0, 1], show=False,
           title="scANVI Predictions (scVI UMAP)", legend_loc="right margin", s=10)

# Cell type original on scVI
adata_merged.obs["_ref_labels_only"] = pd.Series(pd.NA, index=adata_merged.obs_names, dtype="object")
adata_merged.obs.loc[mask_ref, "_ref_labels_only"] = adata_merged.obs.loc[mask_ref, "cell_type_original"].astype(str).values
adata_merged.obs["_ref_labels_only"] = adata_merged.obs["_ref_labels_only"].astype("category")
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="_ref_labels_only", ax=axes[0, 2], show=False,
           title="Reference Labels (scVI UMAP)", legend_loc="right margin", s=10)

# scANVI UMAP
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="data_source", ax=axes[1, 0], show=False,
           title="Data Source (scANVI UMAP)", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="scanvi_pred", ax=axes[1, 1], show=False,
           title="scANVI Predictions (scANVI UMAP)", legend_loc="right margin", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="_ref_labels_only", ax=axes[1, 2], show=False,
           title="Reference Labels (scANVI UMAP)", legend_loc="right margin", s=10)

plt.tight_layout()
fig_path2 = output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf"
plt.savefig(fig_path2, dpi=300, bbox_inches="tight")
plt.close()
print(f"  -> Saved: {fig_path2}")

# Clean up temporary columns
adata_merged.obs.drop(columns=["_ref_labels_only"], inplace=True, errors="ignore")

# --- Step 19: Summary ---
print("\n" + "=" * 80)
print("PIPELINE COMPLETE")
print("=" * 80)
print(f"\nOutput Directory: {output_dir}")
print(f"\nKey Files:")
print(f"  - Results: {OUTPUT_PREFIX}_results.h5ad")
print(f"  - Training data: {OUTPUT_PREFIX}_train_HVG.h5ad")
print(f"  - scVI model: {OUTPUT_PREFIX}_scvi_model/")
print(f"  - scANVI model: {OUTPUT_PREFIX}_scanvi_model/")
print(f"  - UMAP operator (scANVI): {OUTPUT_PREFIX}_umap_scanvi_operator.joblib")
print(f"  - Config: {OUTPUT_PREFIX}_config.json")
print(f"  - Overview plot (scANVI): {OUTPUT_PREFIX}_overview_scANVI.pdf")
print(f"  - UMAP comparison: {OUTPUT_PREFIX}_umap_comparison.pdf")
print(f"  - Confidence hist: {OUTPUT_PREFIX}_confidence_histogram.png")

print(f"\nData Summary:")
print(f"  - Total cells: {adata_merged.n_obs:,}")
print(f"  - Reference: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
print(f"  - Query: {(adata_merged.obs['data_source'] == 'query').sum():,}")
print(f"  - HVG genes: {n_hvg_final:,}")
print(f"  - Full genes in .raw: {adata_merged.raw.n_vars:,}")

print(f"\nObs columns added:")
for col in ["data_source", "cell_type_original", "scanvi_pred", "scanvi_confidence",
            "celltypist_pred", "celltypist_confidence"]:
    if col in adata_merged.obs.columns:
        print(f"  - {col}")

print(f"\nObsm keys:")
for key in adata_merged.obsm.keys():
    print(f"  - {key}: {adata_merged.obsm[key].shape}")

print(f"\nUMAP spaces available:")
print(f"  - X_umap: DEFAULT (scANVI-based)")
print(f"  - X_umap_scVI: scVI latent space")
print(f"  - X_umap_scANVI: scANVI latent space")
print(f"  -> Use sc.pl.umap(adata) for default, or sc.pl.umap(adata, basis='X_umap_scVI') for scVI")

print("\n" + "=" * 80)
print("NEXT STEPS:")
print("=" * 80)
print("1. Check overview plots for visualization")
print("2. Load results: adata = sc.read_h5ad('..._results.h5ad')")
print("")
print("UMAP visualization options:")
print("  - sc.pl.umap(adata)                    # DEFAULT (scANVI)")
print("  - sc.pl.umap(adata, basis='X_umap_scVI')   # scVI latent space")
print("  - sc.pl.umap(adata, basis='X_umap_scANVI') # scANVI latent space")
print("")
print("3. Compare scanvi_pred between reference and query")
print("4. Identify novel cell types in query")
print("5. Use adata.raw for DE analysis (full gene set)")
print("=" * 80)
