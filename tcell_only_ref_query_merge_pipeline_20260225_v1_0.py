#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
T Cell ONLY Reference-Query Merge Pipeline v1.0
==============================================================================

Purpose:
--------
Specialized pipeline for T cell analysis with focus on CD4/CD8 distinction.

Key Features:
-------------
1. T cell specific marker prioritization (CD3, CD4, CD8, etc.)
2. CD4/CD8 scoring for validation
3. Dual annotation: Coarse (CD4/CD8) + Fine (subtypes)
4. Direct merge (NOT scArches) - allows discovery of new T cell subsets
5. Reference cells keep labels, query cells marked as "Unknown"

CD4/CD8 Classification Strategy:
--------------------------------
- Primary markers: CD4 (CD4+), CD8A/CD8B (CD8+)
- Scoring: CD4_score vs CD8_score for validation
- Handles: DP (double positive), DN (double negative), NK-T

Input:
------
- Reference h5ad: T cells with annotations (coarse + fine labels preferred)
- Query h5ad: T cells raw counts, may contain new subsets

Output:
-------
- Merged h5ad with scVI/scANVI latent spaces
- CD4/CD8 scores for validation
- Coarse (CD4/CD8) and Fine (subtypes) annotations
- CellTypist predictions
- UMAP visualizations (scVI, scANVI)

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
from scipy.sparse import issparse, csr_matrix
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
print("T Cell ONLY Ref-Query Merge Pipeline v1.0")
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
REFERENCE_H5AD = "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"  # EDIT THIS
QUERY_H5AD = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/t_cells.h5ad"          # EDIT THIS

# --- Column Configuration ---
# Reference annotation columns
REF_LABEL_COARSE = "cell_type_L2"   # e.g., "CD4", "CD8", "NK-T", "DN", "DP"
REF_LABEL_COARSE = "cell_type_L3"       # e.g., "CD4_Naive", "CD8_TEM", etc.

# If you only have one label column, set both to the same:
# REF_LABEL_COARSE = "cell_type"
# REF_LABEL_FINE = "cell_type"

BATCH_KEY = "sample"      # Batch/sample column
TISSUE_KEY = "tissue"     # Tissue column

# --- Output Configuration ---
OUTPUT_DIR = "/home/h2048/data/py/20260225/tcell_only_merged_pipeline"
OUTPUT_PREFIX = "tcell_only_merged"

# --- Cell Type Selection ---
# Include only specific T cell types from reference (None = all)
INCLUDE_COARSE_TYPES = None  # e.g., ["CD4", "CD8"]

# --- Pipeline Parameters ---
N_HVG = 4000
FORCE_MARKERS_IN_HVG = True

# scVI Parameters
SCVI_N_LATENT = 100
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
CELLTYPIST_MODEL = "/home/h2048/data/source/reference/celltypist_models/Immune_All_Low.pkl"
CELLTYPIST_MAJORITY_VOTE = True

# --- T Cell Gene Signatures ---

# CD4/CD8 Score Thresholds
CD4_SCORE_THRESHOLD = 0.3
CD8_SCORE_THRESHOLD = 0.3

# Core T cell markers (pan-T)
TCELL_CORE_MARKERS = ["CD3D", "CD3E", "CD3G", "PTPRC"]

# CD4 vs CD8 markers
CD4_MARKERS = ["CD4", "IL7R", "CD40LG"]
CD8_MARKERS = ["CD8A", "CD8B"]

# Naive/Memory markers
NAIVE_MARKERS = ["CCR7", "SELL", "TCF7", "LEF1", "CD27"]
CM_MARKERS = ["CCR7", "CD27", "IL7R"]  # Central memory
EM_MARKERS = ["GZMK", "CXCR3", "CCR5"]  # Effector memory
TEMRA_MARKERS = ["GZMB", "PRF1", "GNLY", "NKG7"]  # Effector/terminally differentiated

# Treg markers
TREG_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "IKZF2"]

# Th1/Th2/Th17 markers
TH1_MARKERS = ["TBX21", "IFNG", "CXCR3"]
TH2_MARKERS = ["GATA3", "IL4", "IL5", "IL13"]
TH17_MARKERS = ["RORC", "IL17A", "IL17F", "CCR6"]

# Tissue resident
TRM_MARKERS = ["CD69", "ITGAE", "CXCR6"]

# Proliferation
PROLIF_MARKERS = ["MKI67", "TOP2A", "PCNA"]

# Stress signature
STRESS_SIGNATURE_GENES = [
    "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "IER2"
]

# Cell cycle
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

# All forced markers for HVG selection
FORCED_MARKERS = list(set(
    TCELL_CORE_MARKERS + CD4_MARKERS + CD8_MARKERS +
    NAIVE_MARKERS + CM_MARKERS + EM_MARKERS + TEMRA_MARKERS +
    TREG_MARKERS + TH1_MARKERS + TH2_MARKERS + TH17_MARKERS +
    TRM_MARKERS + PROLIF_MARKERS
))

# --- Reproducibility ---
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

def ensure_counts_layer(adata, counts_layer="counts"):
    """P1-FIX: Robust counts validation with float32 tolerance."""
    if counts_layer not in (adata.layers or {}):
        raise ValueError(
            f"CRITICAL ERROR: layers['{counts_layer}'] not found. "
            "scVI requires raw counts!"
        )

    X_counts = adata.layers[counts_layer]
    if issparse(X_counts):
        sample_data = X_counts.data[:1000]
    else:
        sample_data = X_counts.flat[:1000]

    sample = np.asarray(sample_data, dtype=np.float64)

    if np.any(sample < 0):
        raise ValueError("counts contains negative values!")

    if not np.allclose(sample, np.round(sample), atol=1e-6):
        raise ValueError(
            "counts looks non-integer (possible normalized/log data). "
            f"Sample range: [{sample.min():.4f}, {sample.max():.4f}]"
        )

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


def compute_module_score(adata, gene_list, score_name, counts_layer="counts"):
    """Compute module score using scanpy's score_genes for better accuracy."""
    genes = [g for g in gene_list if g in adata.var_names]
    if len(genes) < 3:
        adata.obs[score_name] = 0.0
        return

    # Create temporary AnnData with log1p-normalized data for scoring
    adata_tmp = sc.AnnData(X=adata.layers[counts_layer].copy(), var=adata.var.copy())
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)

    # Use scanpy's score_genes (subtracts random background genes)
    sc.tl.score_genes(adata_tmp, gene_list=genes, score_name=score_name, random_state=RANDOM_SEED)

    # Copy scores back
    adata.obs[score_name] = adata_tmp.obs[score_name].values

    # Normalize to 0-1 range for thresholding
    scores = adata.obs[score_name]
    mn, mx = float(scores.min()), float(scores.max())
    if mx > mn:
        adata.obs[f"{score_name}_norm"] = (scores - mn) / (mx - mn)
    else:
        adata.obs[f"{score_name}_norm"] = 0.0

    del adata_tmp
    gc.collect()


def compute_cd4_cd8_scores(adata):
    """Compute CD4 and CD8 module scores for validation."""
    print("  -> Computing CD4/CD8 module scores...")

    compute_module_score(adata, CD4_MARKERS, "CD4_score")
    compute_module_score(adata, CD8_MARKERS, "CD8_score")

    # Classify based on normalized scores (using configurable thresholds)
    cd4_high = adata.obs["CD4_score_norm"] > CD4_SCORE_THRESHOLD
    cd8_high = adata.obs["CD8_score_norm"] > CD8_SCORE_THRESHOLD

    classifications = []
    for i in range(len(adata)):
        cd4 = cd4_high.iloc[i]
        cd8 = cd8_high.iloc[i]
        if cd4 and not cd8:
            classifications.append("CD4_single")
        elif cd8 and not cd4:
            classifications.append("CD8_single")
        elif cd4 and cd8:
            classifications.append("DP")
        else:
            classifications.append("DN")

    adata.obs["cd4_cd8_by_score"] = classifications
    print(f"    CD4_single: {sum([c == 'CD4_single' for c in classifications])}")
    print(f"    CD8_single: {sum([c == 'CD8_single' for c in classifications])}")
    print(f"    DP: {sum([c == 'DP' for c in classifications])}")
    print(f"    DN: {sum([c == 'DN' for c in classifications])}")


def prepare_covariates(adata):
    """Prepare covariates for scVI training."""
    print("  -> Validating counts...")
    ensure_counts_layer(adata, "counts")

    print("  -> Batch/Tissue...")
    ensure_batch_tissue(adata, BATCH_KEY, TISSUE_KEY)

    print("  -> Computing signature scores...")
    # MT%
    if "pct_counts_mt" not in adata.obs.columns:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, layer="counts")

    # CD4/CD8 scores
    compute_cd4_cd8_scores(adata)

    # Stress score
    compute_module_score(adata, STRESS_SIGNATURE_GENES, "stress_score")

    # Cell cycle
    if not all(k in adata.obs.columns for k in ["S_score", "G2M_score"]):
        s_in = [g for g in S_GENES if g in adata.var_names]
        g_in = [g for g in G2M_GENES if g in adata.var_names]

        if len(s_in) >= 5 and len(g_in) >= 5:
            cc_union = list(dict.fromkeys(s_in + g_in))
            cc_idx = adata.var_names.get_indexer(cc_union)
            X_cc = adata.layers["counts"][:, cc_idx].copy()
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
        else:
            adata.obs["S_score"] = 0.0
            adata.obs["G2M_score"] = 0.0
            adata.obs["phase"] = "G1"


def run_celltypist(adata, model_name=CELLTYPIST_MODEL, majority_vote=True):
    """Run CellTypist annotation."""
    print("\n[CellTypist] Starting annotation...")

    models.download_models(model=model_name)
    model = models.Model.load(model_name)

    adata_ct = adata.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=majority_vote,
        mode='best match'
    )

    adata.obs["celltypist_pred"] = predictions.predicted_labels.predicted_labels
    adata.obs["celltypist_confidence"] = predictions.probability_matrix.max(axis=1).values

    if majority_vote:
        adata.obs["celltypist_majority"] = predictions.predicted_labels.majority_voting

    print(f"  -> CellTypist predictions added")
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
if INCLUDE_COARSE_TYPES is not None and REF_LABEL_COARSE in adata_ref.obs.columns:
    print(f"\n[Step 2] Subsetting reference to: {INCLUDE_COARSE_TYPES}")
    mask = adata_ref.obs[REF_LABEL_COARSE].isin(INCLUDE_COARSE_TYPES)
    adata_ref = adata_ref[mask].copy()
    print(f"    Reference after subset: {adata_ref.shape}")

# --- Step 3: Prepare labels ---
print("\n[Step 3] Preparing labels...")

adata_ref.obs["data_source"] = "reference"
adata_qry.obs["data_source"] = "query"

# Coarse labels (CD4/CD8 level)
if REF_LABEL_COARSE in adata_ref.obs.columns:
    adata_ref.obs["cell_type_coarse"] = adata_ref.obs[REF_LABEL_COARSE].astype(str)
    print(f"  -> Coarse labels from: {REF_LABEL_COARSE}")
    print(adata_ref.obs["cell_type_coarse"].value_counts())
else:
    print(f"  WARNING: {REF_LABEL_COARSE} not found, using 'T_cell'")
    adata_ref.obs["cell_type_coarse"] = "T_cell"

# Fine labels (subtype level)
if REF_LABEL_FINE in adata_ref.obs.columns:
    adata_ref.obs["cell_type_fine"] = adata_ref.obs[REF_LABEL_FINE].astype(str)
    print(f"  -> Fine labels from: {REF_LABEL_FINE}")
else:
    # Use coarse as fine if no separate fine labels
    adata_ref.obs["cell_type_fine"] = adata_ref.obs["cell_type_coarse"]
    print(f"  -> Using coarse labels as fine labels")

# Query labels = Unknown
adata_qry.obs["cell_type_coarse"] = UNLABELED_CATEGORY
adata_qry.obs["cell_type_fine"] = UNLABELED_CATEGORY

# --- Step 4: Find common genes ---
print("\n[Step 4] Finding common genes...")

ref_genes = set(adata_ref.var_names)
qry_genes = set(adata_qry.var_names)
common_genes = list(ref_genes.intersection(qry_genes))

print(f"  Reference: {len(ref_genes):,}, Query: {len(qry_genes):,}, Common: {len(common_genes):,}")

if len(common_genes) < 1000:
    raise ValueError(f"Too few common genes ({len(common_genes)}). Check gene naming!")

adata_ref = adata_ref[:, common_genes].copy()
adata_qry = adata_qry[:, common_genes].copy()

# --- Step 5: Validate counts ---
print("\n[Step 5] Validating counts...")
ensure_counts_layer(adata_ref, "counts")
ensure_counts_layer(adata_qry, "counts")

# --- Step 6: Concatenate ---
print("\n[Step 6] Concatenating...")

adata_ref.obs_names = pd.Index([f"ref_{x}" for x in adata_ref.obs_names])
adata_qry.obs_names = pd.Index([f"qry_{x}" for x in adata_qry.obs_names])

adata_merged = sc.concat(
    {"reference": adata_ref, "query": adata_qry},
    axis=0,
    join="inner",
    merge="unique",
    label="data_source"
)

print(f"  Merged: {adata_merged.shape}")
print(f"  Reference: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
print(f"  Query: {(adata_merged.obs['data_source'] == 'query').sum():,}")

del adata_ref, adata_qry
gc.collect()

# --- Step 7: Prepare covariates ---
print("\n[Step 7] Preparing covariates...")
prepare_covariates(adata_merged)

# P1-FIX: Add symbol_base
print("  -> Adding symbol_base column...")
adata_merged.var["symbol_base"] = adata_merged.var_names.str.replace(r"-\d+$", "", regex=True)

# --- Step 8: HVG Selection ---
print("\n[Step 8] Selecting HVGs...")
hvg_method = "unknown"
try:
    sc.pp.highly_variable_genes(
        adata_merged, layer="counts", n_top_genes=N_HVG,
        batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
    )
    hvg_method = "batch_seurat_v3"
except Exception as e1:
    print(f"  -> batch-aware failed ({str(e1)[:50]}), trying standard...")
    try:
        sc.pp.highly_variable_genes(
            adata_merged, layer="counts", n_top_genes=N_HVG,
            flavor="seurat_v3", subset=False
        )
        hvg_method = "standard_seurat_v3"
    except Exception as e2:
        print(f"  -> standard failed ({str(e2)[:50]}), fallback to cell_ranger")
        sc.pp.highly_variable_genes(
            adata_merged, layer="counts", n_top_genes=N_HVG,
            flavor="cell_ranger", subset=False
        )
        hvg_method = "cell_ranger"

print(f"  -> Method: {hvg_method}")

# Force T cell markers into HVG
if FORCE_MARKERS_IN_HVG:
    n_added = 0
    marker_set = set(FORCED_MARKERS)
    for idx, symbol_base in enumerate(adata_merged.var["symbol_base"]):
        if symbol_base in marker_set:
            real_name = adata_merged.var_names[idx]
            if not adata_merged.var.loc[real_name, "highly_variable"]:
                adata_merged.var.loc[real_name, "highly_variable"] = True
                n_added += 1
    print(f"  -> Forced {n_added}/{len(FORCED_MARKERS)} markers into HVG")

n_hvg_final = adata_merged.var["highly_variable"].sum()
print(f"  -> Final HVG count: {n_hvg_final}")

hvg_genes = adata_merged.var_names[adata_merged.var["highly_variable"]].tolist()
with open(output_dir / f"{OUTPUT_PREFIX}_hvg_genes.txt", "w") as f:
    f.write("\n".join(hvg_genes))

# --- Step 9: Build FULL matrix for .raw ---
print("\n[Step 9] Building full matrix for .raw...")
full_counts = adata_merged.layers["counts"]
if issparse(full_counts) and not isinstance(full_counts, csr_matrix):
    full_counts = csr_matrix(full_counts)
raw_var = adata_merged.var.copy()
print(f"  Full matrix shape: {full_counts.shape}")

# --- Step 10: Create training subset ---
print("\n[Step 10] Creating training subset...")

hvg_mask = adata_merged.var["highly_variable"].values
X_hvg = adata_merged.layers["counts"][:, hvg_mask]

if issparse(X_hvg) and not isinstance(X_hvg, csr_matrix):
    X_hvg = csr_matrix(X_hvg)

adata_train = sc.AnnData(
    X=X_hvg.copy(),
    obs=adata_merged.obs.copy(),
    var=adata_merged.var.iloc[hvg_mask].copy()
)
adata_train.var_names = adata_merged.var_names[hvg_mask]
adata_train.layers["counts"] = adata_train.X

print(f"  Training data: {adata_train.shape}")

# Use FINE labels for scANVI training
adata_train.obs["scanvi_labels"] = adata_train.obs["cell_type_fine"].astype(str)
adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")

if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
    adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])

print(f"  -> Label distribution:")
print(adata_train.obs["scanvi_labels"].value_counts())

gc.collect()

# --- Step 11: CellTypist ---
print("\n[Step 11] Running CellTypist...")
try:
    run_celltypist(adata_train)
    adata_merged.obs["celltypist_pred"] = adata_train.obs["celltypist_pred"]
    adata_merged.obs["celltypist_confidence"] = adata_train.obs["celltypist_confidence"]
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
print("  -> scVI complete")

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
print("  -> scANVI complete")

# --- Step 14: Export Results ---
print("\n[Step 14] Exporting Results...")

latent_scanvi = scanvi_model.get_latent_representation(adata_train)
latent_df = pd.DataFrame(
    latent_scanvi,
    index=adata_train.obs_names,
    columns=[f"scANVI_{i}" for i in range(latent_scanvi.shape[1])]
)
latent_aligned = latent_df.reindex(adata_merged.obs_names)

if latent_aligned.isna().any().any():
    raise ValueError("CRITICAL: Missing latent representation!")

adata_merged.obsm["X_scANVI"] = latent_aligned.values

pred_labels = scanvi_model.predict(adata_train)
pred_df = pd.Series(pred_labels, index=adata_train.obs_names)
pred_aligned = pred_df.reindex(adata_merged.obs_names)
adata_merged.obs["scanvi_pred"] = pred_aligned.values

proba = scanvi_model.predict(adata_train, soft=True).astype(np.float32)
label_order = scanvi_model.labels
proba_df = pd.DataFrame(proba, index=adata_train.obs_names, columns=label_order)
proba_aligned = proba_df.reindex(adata_merged.obs_names)
adata_merged.obsm["scanvi_proba"] = proba_aligned.values
adata_merged.obs["scanvi_confidence"] = proba_aligned.values.max(axis=1)

print(f"  -> scANVI predictions:")
print(adata_merged.obs["scanvi_pred"].value_counts().head(15))

# --- Step 15: Attach .raw ---
print("\n[Step 15] Attaching .raw...")
from anndata import AnnData

adata_merged.raw = AnnData(X=full_counts, obs=adata_merged.obs.copy(), var=raw_var)
print(f"  OK .raw: {adata_merged.raw.n_vars} genes")

# --- Step 16: Multiple UMAPs (scVI + scANVI) ---
print("\n[Step 16] Computing Multiple UMAPs...")

# Get scVI latent for UMAP
print("  -> Getting scVI latent representation...")
adata_merged.obsm["X_scVI"] = scvi_model.get_latent_representation(adata_train)

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
adata_merged.obsm["X_umap"] = adata_merged.obsm["X_umap_scANVI"].copy()
print(f"     Saved to X_umap_scANVI and X_umap (default)")

# Save UMAP operator for scANVI
umap_op_scanvi = UMAP(n_neighbors=30, n_components=2, min_dist=0.5, spread=1.0,
                      metric="euclidean", random_state=RANDOM_SEED)
umap_op_scanvi.fit(adata_merged.obsm["X_scANVI"])
joblib.dump(umap_op_scanvi, output_dir / f"{OUTPUT_PREFIX}_umap_scanvi_operator.joblib")
print(f"  -> UMAP operator (scANVI) saved")

# P2: Confidence histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(adata_merged.obs["scanvi_confidence"], bins=50, edgecolor="black")
ax.set_xlabel("Prediction Confidence")
ax.set_ylabel("Cell Count")
ax.set_title("scANVI Confidence Distribution")
plt.tight_layout()
plt.savefig(output_dir / f"{OUTPUT_PREFIX}_confidence_histogram.png", dpi=150)
plt.close()

# --- Step 17: Save Results ---
print("\n[Step 17] Saving Results...")

scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
scvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scvi_model", overwrite=True)

config = {
    "version": "1.0",
    "timestamp": datetime.now().isoformat(),
    "input": {"reference": REFERENCE_H5AD, "query": QUERY_H5AD},
    "n_hvg": int(n_hvg_final),
    "scanvi_labels": list(label_order),
    "tcell_specific": {
        "cd4_markers": CD4_MARKERS,
        "cd8_markers": CD8_MARKERS,
        "cd4_score_threshold": CD4_SCORE_THRESHOLD,
        "cd8_score_threshold": CD8_SCORE_THRESHOLD,
        "forced_markers_count": len(FORCED_MARKERS)
    },
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

output_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"
adata_merged.write_h5ad(output_h5ad, compression="gzip")
print(f"  -> Saved: {output_h5ad}")

train_h5ad = output_dir / f"{OUTPUT_PREFIX}_train_HVG.h5ad"
adata_train.write_h5ad(train_h5ad, compression="gzip")

# --- Step 18: T Cell Specific Visualization ---
print("\n[Step 18] Creating T cell visualizations...")

# Figure 1: Overview with CD4/CD8 validation (scANVI UMAP)
fig = plt.figure(figsize=(24, 20))
gs = fig.add_gridspec(5, 4, hspace=0.3, wspace=0.3)

# Row 1: Data source and labels
ax1 = fig.add_subplot(gs[0, 0])
sc.pl.umap(adata_merged, color="data_source", ax=ax1, show=False, title="Data Source (scANVI)", s=15)

ax2 = fig.add_subplot(gs[0, 1])
adata_merged.obs["_ref_coarse"] = pd.Series(pd.NA, index=adata_merged.obs_names, dtype="object")
mask_ref = adata_merged.obs["data_source"] == "reference"
adata_merged.obs.loc[mask_ref, "_ref_coarse"] = adata_merged.obs.loc[mask_ref, "cell_type_coarse"].astype(str).values
adata_merged.obs["_ref_coarse"] = adata_merged.obs["_ref_coarse"].astype("category")
sc.pl.umap(adata_merged, color="_ref_coarse", ax=ax2, show=False, title="Reference Coarse Labels", legend_loc="on data", s=15)

ax3 = fig.add_subplot(gs[0, 2])
sc.pl.umap(adata_merged, color="scanvi_pred", ax=ax3, show=False, title="scANVI Predictions", legend_loc="on data", s=15)

ax4 = fig.add_subplot(gs[0, 3])
sc.pl.umap(adata_merged, color="scanvi_confidence", ax=ax4, show=False, title="Confidence", cmap="viridis", vmin=0, vmax=1, s=15)

# Row 2: CD4/CD8 scores and classification
ax5 = fig.add_subplot(gs[1, 0])
sc.pl.umap(adata_merged, color="CD4_score", ax=ax5, show=False, title="CD4 Module Score", cmap="Reds", s=15)

ax6 = fig.add_subplot(gs[1, 1])
sc.pl.umap(adata_merged, color="CD8_score", ax=ax6, show=False, title="CD8 Module Score", cmap="Blues", s=15)

ax7 = fig.add_subplot(gs[1, 2])
sc.pl.umap(adata_merged, color="cd4_cd8_by_score", ax=ax7, show=False, title="CD4/CD8 by Score", legend_loc="on data", s=15)

ax8 = fig.add_subplot(gs[1, 3])
# CD4 vs CD8 score scatter
query_mask = adata_merged.obs["data_source"] == "query"
ax8.scatter(
    adata_merged.obs.loc[query_mask, "CD4_score"],
    adata_merged.obs.loc[query_mask, "CD8_score"],
    c=adata_merged.obs.loc[query_mask, "scanvi_confidence"],
    cmap="viridis", s=5, alpha=0.5
)
ax8.set_xlabel("CD4 Score")
ax8.set_ylabel("CD8 Score")
ax8.set_title(f"Query: CD4 vs CD8 Score\n(threshold={CD4_SCORE_THRESHOLD})")
ax8.axhline(y=CD8_SCORE_THRESHOLD, color='k', linestyle='--', alpha=0.3)
ax8.axvline(x=CD4_SCORE_THRESHOLD, color='k', linestyle='--', alpha=0.3)

# Row 3: Core T cell markers + CD8B
ax9 = fig.add_subplot(gs[2, 0])
if "CD3E" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="CD3E", ax=ax9, show=False, title="CD3E (pan-T)", cmap="Reds", s=15, use_raw=True)

ax10 = fig.add_subplot(gs[2, 1])
if "CD4" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="CD4", ax=ax10, show=False, title="CD4", cmap="Reds", s=15, use_raw=True)

ax11 = fig.add_subplot(gs[2, 2])
if "CD8A" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="CD8A", ax=ax11, show=False, title="CD8A", cmap="Reds", s=15, use_raw=True)

ax12 = fig.add_subplot(gs[2, 3])
if "CD8B" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="CD8B", ax=ax12, show=False, title="CD8B", cmap="Reds", s=15, use_raw=True)

# Row 4: Functional markers + Proliferation
ax13 = fig.add_subplot(gs[3, 0])
if "CCR7" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="CCR7", ax=ax13, show=False, title="CCR7 (Naive/CM)", cmap="Reds", s=15, use_raw=True)

ax14 = fig.add_subplot(gs[3, 1])
if "FOXP3" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="FOXP3", ax=ax14, show=False, title="FOXP3 (Treg)", cmap="Reds", s=15, use_raw=True)

ax15 = fig.add_subplot(gs[3, 2])
if "GZMB" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="GZMB", ax=ax15, show=False, title="GZMB (Effector)", cmap="Reds", s=15, use_raw=True)

ax16 = fig.add_subplot(gs[3, 3])
if "MKI67" in adata_merged.raw.var_names:
    sc.pl.umap(adata_merged, color="MKI67", ax=ax16, show=False, title="MKI67 (Proliferating)", cmap="Reds", s=15, use_raw=True)

# Row 5: Query-specific analysis
ax17 = fig.add_subplot(gs[4, 0])
# Compare reference vs query CD4/CD8 distribution
ref_cd4cd8 = adata_merged.obs.loc[adata_merged.obs["data_source"] == "reference", "cell_type_coarse"].value_counts()
qry_cd4cd8 = adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "cd4_cd8_by_score"].value_counts()
x = np.arange(len(ref_cd4cd8.index))
width = 0.35
ax17.bar(x - width/2, ref_cd4cd8.values, width, label="Reference (annotated)", alpha=0.8)
ax17.bar(x + width/2, [qry_cd4cd8.get(k, 0) for k in ref_cd4cd8.index], width, label="Query (by score)", alpha=0.8)
ax17.set_xticks(x)
ax17.set_xticklabels(ref_cd4cd8.index, rotation=45, ha="right")
ax17.set_ylabel("Cell Count")
ax17.set_title("CD4/CD8 Distribution")
ax17.legend()

ax18 = fig.add_subplot(gs[4, 1])
# Novel cell types in query (high confidence, not matching reference labels)
query_cells = adata_merged[adata_merged.obs["data_source"] == "query"]
high_conf = query_cells.obs["scanvi_confidence"] >= 0.7
novel_types = query_cells.obs.loc[high_conf, "scanvi_pred"].value_counts().head(10)
ax18.barh(range(len(novel_types)), novel_types.values)
ax18.set_yticks(range(len(novel_types)))
ax18.set_yticklabels(novel_types.index)
ax18.set_xlabel("Cell Count")
ax18.set_title("Top Predicted Types in Query (conf >= 0.7)")

ax19 = fig.add_subplot(gs[4, 2])
# Confidence by data source
conf_ref = adata_merged.obs.loc[adata_merged.obs["data_source"] == "reference", "scanvi_confidence"]
conf_qry = adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "scanvi_confidence"]
ax19.hist([conf_ref, conf_qry], bins=30, label=["Reference", "Query"], alpha=0.7)
ax19.set_xlabel("Confidence")
ax19.set_ylabel("Cell Count")
ax19.set_title("Confidence Distribution by Source")
ax19.legend()

ax20 = fig.add_subplot(gs[4, 3])
sc.pl.umap(adata_merged, color="celltypist_pred", ax=ax20, show=False, title="CellTypist Predictions", legend_loc="on data", s=15)

plt.savefig(output_dir / f"{OUTPUT_PREFIX}_tcell_overview.pdf", dpi=300, bbox_inches="tight")
plt.close()
print(f"  -> Saved: {OUTPUT_PREFIX}_tcell_overview.pdf")

# Figure 2: scVI vs scANVI UMAP comparison
print("  -> Creating scVI vs scANVI UMAP comparison...")
fig2, axes2 = plt.subplots(2, 3, figsize=(18, 12))

# scVI UMAP plots
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="data_source", ax=axes2[0, 0], show=False,
           title="Data Source (scVI UMAP)", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="_ref_coarse", ax=axes2[0, 1], show=False,
           title="Reference Labels (scVI)", legend_loc="on data", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scVI", color="scanvi_pred", ax=axes2[0, 2], show=False,
           title="scANVI Predictions (scVI)", legend_loc="on data", s=10)

# scANVI UMAP plots
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="data_source", ax=axes2[1, 0], show=False,
           title="Data Source (scANVI UMAP)", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="_ref_coarse", ax=axes2[1, 1], show=False,
           title="Reference Labels (scANVI)", legend_loc="on data", s=10)
sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="scanvi_pred", ax=axes2[1, 2], show=False,
           title="scANVI Predictions (scANVI)", legend_loc="on data", s=10)

plt.tight_layout()
plt.savefig(output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf", dpi=300, bbox_inches="tight")
plt.close()
print(f"  -> Saved: {OUTPUT_PREFIX}_umap_comparison.pdf")

# Clean up
adata_merged.obs.drop(columns=["_ref_coarse"], inplace=True, errors="ignore")

# --- Summary ---
print("\n" + "=" * 80)
print("T CELL PIPELINE COMPLETE")
print("=" * 80)
print(f"\nOutput: {output_dir}")
print(f"\nKey observations:")
print(f"  - Total cells: {adata_merged.n_obs:,}")
print(f"  - Reference: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
print(f"  - Query: {(adata_merged.obs['data_source'] == 'query').sum():,}")

# CD4/CD8 breakdown
print(f"\nCD4/CD8 Distribution in Query (by score):")
print(adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "cd4_cd8_by_score"].value_counts())

print(f"\nTop scANVI predictions in Query:")
print(adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "scanvi_pred"].value_counts().head(10))

print(f"\nUMAP spaces available:")
print(f"  - X_umap: DEFAULT (scANVI-based)")
print(f"  - X_umap_scVI: scVI latent space")
print(f"  - X_umap_scANVI: scANVI latent space")
print(f"  -> Use sc.pl.umap(adata) for default, or sc.pl.umap(adata, basis='X_umap_scVI') for scVI")

print("\n" + "=" * 80)
print("CD4/CD8 Validation Strategy:")
print("=" * 80)
print("1. Compare 'scanvi_pred' (from reference labels) with 'cd4_cd8_by_score'")
print(f"   - Score thresholds: CD4 >= {CD4_SCORE_THRESHOLD}, CD8 >= {CD8_SCORE_THRESHOLD}")
print("2. High-confidence query cells with mismatches may be novel subsets")
print("3. Use CD4_score vs CD8_score scatter to identify DP/DN populations")
print("4. Check core T cell markers (CD3E) to filter out non-T contaminants")
print("=" * 80)
