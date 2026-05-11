#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
Production Script: MERGED scVI+scANVI Training - HOTFIX v2.5.5
==============================================================================

Version: 2.5.5-HOTFIX (Production-Grade)
Date: 2026-02-08
Author: r2end

Critical Fixes (P0):
1. ✅ UMAP writes to default X_umap (not X_umap_scanvi)
2. ✅ No adata.raw (use layer-based visualization instead)

High-Risk Fixes (P1):
3. ✅ Robust integer validation (float32 tolerance)
4. ✅ Single-copy counts in adata_train (no duplicate)
5. ✅ Minimal .obs construction for compatibility checks
6. ✅ Symbol_base matching for FORCED_MARKERS
7. ✅ Hard error on missing latent (no silent fillna)
8. ✅ HPC-safe settings (dl_num_workers=0, thread control)

Enhancements (P2):
9. ✅ Proba columns with label names
10. ✅ Comprehensive diagnostics (confidence hist, full counts)

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
print("MERGED scVI+scANVI Training (v2.5.5-HOTFIX)")
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
# 1. CONFIGURATION
# ==============================================================================

# --- IO Paths ---
INPUT_MERGED_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3_MERGE_v2_0/reference_plus_query_merged_v2_0.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0207/merged_scanvi_L2_prod_v2_5_5_HOTFIX"

# --- Column Keys ---
DATASOURCE_KEY = "data_source"
REF_LABEL_KEY = "Cell_Type_L2"
QRY_LABEL_KEY = "Cell_Type_L2_final"
QRY_CONF_KEY = "mapping_confidence"
TRAIN_LABEL_KEY = "Cell_Type_L2_train"
UNLABELED_CATEGORY = "Unknown"

# --- Covariate Keys ---
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"

# --- Query Filter Settings ---
USE_QUERY_CONF_FILTER = True
QUERY_CONF_THR = 0.5

# --- HVG Settings ---
N_HVG = 4000
FORCE_MARKERS_IN_HVG = True

# --- Training Params ---
SCVI_N_LATENT = 50
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_DROPOUT = 0.1

MAX_EPOCHS_SCVI = 400
MAX_EPOCHS_SCANVI = 200
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 0.0

# --- Gene Signatures ---
STRESS_SIGNATURE_GENES = [
    "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "IER2"
]

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UHRF1",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS",
    "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP",
    "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45",
    "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN",
    "POLA1", "CHAF1B", "BRIP1", "E2F8"
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB1", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E",
    "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
    "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8",
    "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5",
    "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
]

FORCED_MARKERS = [
    "CD19","MS4A1","CD79A","CD79B","BANK1","BLK","FCRL2","FCRL3","FCRL5",
    "TCL1A","FCER2","IGHD","IL4R","CD27","TNFRSF13B","AIM2",
    "AICDA","BCL6","MME","LMO2","STMN1","CXCR4","RGS13",
    "SDC1","CD38","XBP1","PRDM1","JCHAIN","IGHA1","IGHA2","IGHG1","IGHG2",
    "IGHG3","IGHG4","MZB1","DERL3","SSR4","ITGAX","TBX21","CXCR3"
]

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
            "scVI requires raw counts! Please regenerate merged object."
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

# ==============================================================================
# 3. MAIN PIPELINE
# ==============================================================================

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

# --- Step 1: Load Data ---
print("\n[Step 1] Loading MERGED data...")
adata_full = sc.read_h5ad(INPUT_MERGED_H5AD)
adata_full.var_names_make_unique()
print(f"Shape: {adata_full.shape}")
print(f"Data sources: {adata_full.obs[DATASOURCE_KEY].value_counts().to_dict()}")

# P1-FIX: Add symbol_base for robust marker matching
print("  -> Adding symbol_base column (P1-FIX)...")
adata_full.var["symbol_base"] = adata_full.var_names.str.replace(r"-\d+$", "", regex=True)

# --- Step 2: Prepare Covariates ---
print("\n[Step 2] Preparing Covariates...")
prepare_covariates(adata_full)

# --- Step 3: Build Unified Labels ---
print("\n[Step 3] Building Unified Labels (Ref + High-Conf Query)...")
is_ref = adata_full.obs[DATASOURCE_KEY] == "reference"
is_qry = adata_full.obs[DATASOURCE_KEY] == "query"

unified_labels = pd.Series(
    [UNLABELED_CATEGORY] * len(adata_full),
    index=adata_full.obs_names,
    dtype="object"
)

# Fill Reference
ref_labels = adata_full.obs.loc[is_ref, REF_LABEL_KEY].values
unified_labels.loc[is_ref] = ref_labels

# Fill Query (high-confidence)
if USE_QUERY_CONF_FILTER:
    conf_scores = pd.to_numeric(adata_full.obs[QRY_CONF_KEY], errors="coerce").fillna(0.0)
    valid_qry_mask = is_qry & (conf_scores >= QUERY_CONF_THR)
    
    qry_labels = adata_full.obs.loc[valid_qry_mask, QRY_LABEL_KEY].values
    unified_labels.loc[valid_qry_mask] = qry_labels
    
    print(f"  -> Query cells retained: {valid_qry_mask.sum():,} / {is_qry.sum():,} (Conf >= {QUERY_CONF_THR})")
else:
    qry_labels = adata_full.obs.loc[is_qry, QRY_LABEL_KEY].values
    unified_labels.loc[is_qry] = qry_labels

# Clean conversion to categorical
unified_labels = unified_labels.fillna(UNLABELED_CATEGORY)
unified_labels = pd.Categorical(unified_labels)

if UNLABELED_CATEGORY not in unified_labels.categories:
    unified_labels = unified_labels.add_categories([UNLABELED_CATEGORY])

adata_full.obs[TRAIN_LABEL_KEY] = unified_labels

print(f"  -> Label Distribution (Top 15):")
label_counts = adata_full.obs[TRAIN_LABEL_KEY].value_counts()
print(label_counts.head(15))
print(f"  -> Total categories: {len(label_counts)}")
print(f"  -> Unknown cells: {(adata_full.obs[TRAIN_LABEL_KEY] == UNLABELED_CATEGORY).sum():,}")

# P2-ENHANCEMENT: Save full label counts
label_counts.to_csv(output_dir / "train_label_counts_full.csv")
print(f"  -> Full label counts saved")

# --- Step 4: HVG Selection ---
print("\n[Step 4] Selecting HVGs (Robust)...")
hvg_method = "unknown"
try:
    print("  -> Attempting: batch-aware seurat_v3")
    sc.pp.highly_variable_genes(
        adata_full, layer="counts", n_top_genes=N_HVG, 
        batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
    )
    hvg_method = "batch_seurat_v3"
except Exception as e1:
    try:
        print(f"  -> Failed ({str(e1)[:50]}). Attempting: standard seurat_v3")
        sc.pp.highly_variable_genes(
            adata_full, layer="counts", n_top_genes=N_HVG, 
            flavor="seurat_v3", subset=False
        )
        hvg_method = "standard_seurat_v3"
    except Exception as e2:
        print(f"  -> Failed ({str(e2)[:50]}). Fallback: cell_ranger")
        sc.pp.highly_variable_genes(
            adata_full, layer="counts", n_top_genes=N_HVG, 
            flavor="cell_ranger", subset=False
        )
        hvg_method = "cell_ranger"

print(f"  -> Method used: {hvg_method}")

# P1-FIX: Use symbol_base for marker matching
if FORCE_MARKERS_IN_HVG:
    n_added = 0
    marker_set = set(FORCED_MARKERS)
    
    for idx, symbol_base in enumerate(adata_full.var["symbol_base"]):
        if symbol_base in marker_set:
            real_name = adata_full.var_names[idx]
            if not adata_full.var.loc[real_name, "highly_variable"]:
                adata_full.var.loc[real_name, "highly_variable"] = True
                n_added += 1
                print(f"    Added: {real_name} (base: {symbol_base})")
    
    print(f"  -> Forced {n_added} / {len(FORCED_MARKERS)} markers into HVG")

n_hvg_final = adata_full.var["highly_variable"].sum()
print(f"  -> Final HVG count: {n_hvg_final}")

# Save HVG list
hvg_genes = adata_full.var_names[adata_full.var["highly_variable"]].tolist()
with open(output_dir / "hvg_genes.txt", "w") as f:
    f.write("\n".join(hvg_genes))

# --- Step 5: P0-FIX: NO adata.raw (Use Layer-Based Viz) ---
print("\n[Step 5] SKIPPING adata.raw (P0-FIX: Use layer-based viz)")
print("  -> Recommendation: sc.pl.umap(adata, color=['gene'], layer='counts')")
print("  -> This avoids disk/memory doubling and scale confusion")

# --- Step 6: Create Training Subset (P1-FIX: Single Copy) ---
print("\n[Step 6] Creating scVI Training Data (P1-FIX: Single-copy HVG)...")

hvg_mask = adata_full.var["highly_variable"].values
X_hvg = adata_full.layers["counts"][:, hvg_mask]

# Ensure CSR
if issparse(X_hvg) and not isinstance(X_hvg, csr_matrix):
    X_hvg = csr_matrix(X_hvg)

# P1-FIX: Single copy, layer points to X
adata_train = sc.AnnData(
    X=X_hvg.copy(),  # Only ONE copy needed
    obs=adata_full.obs.copy(),
    var=adata_full.var.iloc[hvg_mask].copy()
)
adata_train.var_names = adata_full.var_names[hvg_mask]
adata_train.layers["counts"] = adata_train.X  # ✅ No second .copy()

print(f"  -> Training data: {adata_train.shape}")
print(f"  -> Layers: {list(adata_train.layers.keys())}")
print(f"  -> Memory: counts layer is shared reference to X (0 extra cost)")

# Verify critical columns
assert TRAIN_LABEL_KEY in adata_train.obs.columns
assert UNLABELED_CATEGORY in adata_train.obs[TRAIN_LABEL_KEY].cat.categories
assert "counts" in adata_train.layers

gc.collect()

# --- Step 7: scVI Training ---
print("\n[Step 7] scVI Training...")

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

print("  -> Training scVI...")
scvi_model.train(**train_kwargs)

try:
    final_loss = scvi_model.history['reconstruction_loss_train'].iloc[-1]
    print(f"  ✅ scVI Final Loss: {final_loss:.4f}")
except:
    print("  ✅ scVI Training Complete")

# --- Step 8: scANVI Training ---
print("\n[Step 8] scANVI Training...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata_train,
    labels_key=TRAIN_LABEL_KEY,
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
print("  ✅ scANVI Training Complete")

# --- Step 9: Index-Aligned Writeback (P1-FIX: Hard Error on Missing) ---
print("\n[Step 9] Exporting Results (P1-FIX: Strict Validation)...")

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
    columns=[f"scANVI_L2_{i}" for i in range(latent_scanvi.shape[1])]
)

latent_aligned = latent_df.reindex(adata_full.obs_names)

# P1-FIX: Hard error on missing cells (not silent fillna)
n_missing = latent_aligned.isna().any(axis=1).sum()
if n_missing > 0:
    missing_ids = latent_aligned[latent_aligned.isna().any(axis=1)].index.tolist()
    with open(output_dir / "MISSING_LATENT_CELLS.txt", "w") as f:
        f.write("\n".join(missing_ids))
    
    raise ValueError(
        f"CRITICAL: {n_missing} cells missing latent representation! "
        f"IDs saved to MISSING_LATENT_CELLS.txt. "
        "This should NOT happen - check data consistency."
    )

adata_full.obsm["X_scANVI_L2"] = latent_aligned.values

# Predictions
pred_labels = scanvi_model.predict(adata_train)
pred_df = pd.Series(pred_labels, index=adata_train.obs_names)
pred_aligned = pred_df.reindex(adata_full.obs_names)

if pred_aligned.isna().sum() > 0:
    raise ValueError(
        f"CRITICAL: {pred_aligned.isna().sum()} cells missing predictions!"
    )

adata_full.obs["L2_scanvi_pred"] = pred_aligned.values

# P2-ENHANCEMENT: Soft probabilities with label column names
print("  -> Extracting probabilities (P2-FIX: with label names)...")
proba = scanvi_model.predict(adata_train, soft=True).astype(np.float32)

# Get label order from model
label_order = scanvi_model.labels  # This is the order of proba columns
proba_df = pd.DataFrame(
    proba, 
    index=adata_train.obs_names,
    columns=label_order  # P2-ENHANCEMENT: Use actual label names
)
proba_aligned = proba_df.reindex(adata_full.obs_names)

if proba_aligned.isna().any().any():
    raise ValueError("CRITICAL: Missing probability values!")

adata_full.obsm["proba_L2_scanvi"] = proba_aligned.values
adata_full.obs["L2_scanvi_confidence"] = proba_aligned.values.max(axis=1)

# Save proba column mapping
proba_col_map = {i: label for i, label in enumerate(label_order)}
with open(output_dir / "proba_column_mapping.json", "w") as f:
    json.dump(proba_col_map, f, indent=2)
print(f"  -> Proba column mapping saved")

print(f"  ✅ Results written to adata_full")

# P2-ENHANCEMENT: Confidence histogram
print("\n[P2-DIAGNOSTIC] Plotting confidence distribution...")
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(adata_full.obs["L2_scanvi_confidence"], bins=50, edgecolor="black")
ax.set_xlabel("Prediction Confidence")
ax.set_ylabel("Cell Count")
ax.set_title("scANVI L2 Prediction Confidence Distribution")
plt.tight_layout()
plt.savefig(output_dir / "confidence_histogram.png", dpi=150)
plt.close()
print(f"  -> Confidence histogram saved")

# --- Step 10: UMAP (P0-FIX: Write to Default X_umap) ---
print("\n[Step 10] Fitting UMAP (P0-FIX: Write to X_umap)...")

umap_op = UMAP(
    n_neighbors=30,
    n_components=2,
    min_dist=0.5,
    spread=1.0,
    metric="euclidean",
    random_state=RANDOM_SEED
)

print(f"  -> Fitting UMAP on {len(adata_full):,} cells...")
X_umap = umap_op.fit_transform(adata_full.obsm["X_scANVI_L2"])

# P0-FIX: Write to DEFAULT X_umap (not X_umap_scanvi)
adata_full.obsm["X_umap"] = X_umap
print(f"  ✅ UMAP saved to X_umap (Scanpy default)")

# Save operator
joblib.dump(umap_op, output_dir / "umap_operator.joblib")
print(f"  -> UMAP operator saved")

# --- Step 11: Save Outputs ---
print("\n[Step 11] Saving Results...")

# 1. Model
scanvi_model.save(output_dir / "scanvi_model", overwrite=True)

# 2. Config
config = {
    "version": "2.5.5-HOTFIX",
    "timestamp": datetime.now().isoformat(),
    "input_shape": list(adata_full.shape),
    "n_hvg": int(n_hvg_final),
    "hvg_method": hvg_method,
    "scvi_params": {
        "n_latent": SCVI_N_LATENT,
        "n_layers": SCVI_N_LAYERS,
        "n_hidden": SCVI_N_HIDDEN,
        "dropout": SCVI_DROPOUT
    },
    "covariates": {
        "continuous": setup_kwargs.get("continuous_covariate_keys"),
        "categorical": setup_kwargs.get("categorical_covariate_keys")
    },
    "training_stats": {
        "train_labels": adata_train.obs[TRAIN_LABEL_KEY].value_counts().to_dict(),
        "query_confidence_threshold": QUERY_CONF_THR if USE_QUERY_CONF_FILTER else None
    },
    "fixes_applied": {
        "P0": ["X_umap_default", "no_raw_layer_viz"],
        "P1": ["robust_int_check", "single_copy_counts", "symbol_base_markers", 
               "hard_error_missing_latent", "dl_num_workers_0", "thread_control"],
        "P2": ["proba_label_names", "full_label_counts", "confidence_hist"]
    }
}
with open(output_dir / "training_config.json", "w") as f:
    json.dump(config, f, indent=2)

# 3. H5AD (P0-FIX: No .raw, smaller file)
out_path = output_dir / "merged_scanvi_L2_prod_v2_5_5.h5ad"
adata_full.write_h5ad(out_path, compression="gzip")
print(f"  ✅ Full data saved: {out_path}")
print(f"     - Shape: {adata_full.shape}")
print(f"     - No .raw (use layer='counts' for visualization)")

# 4. Save training subset
train_path = output_dir / "merged_scanvi_L2_train_HVG.h5ad"
adata_train.write_h5ad(train_path, compression="gzip")
print(f"  ✅ Training data saved: {train_path}")

# Clean up
del adata_train, scvi_model
gc.collect()

print(f"\n{'='*80}")
print(f"[SUCCESS] Pipeline Finished (v2.5.5-HOTFIX)")
print(f"Output Directory: {OUTPUT_DIR}")
print(f"{'='*80}")
print(f"\nKey Files:")
print(f"  - Full data: merged_scanvi_L2_prod_v2_5_5.h5ad")
print(f"  - Model: scanvi_model/")
print(f"  - UMAP operator: umap_operator.joblib")
print(f"  - Config: training_config.json")
print(f"  - Diagnostics: confidence_histogram.png, train_label_counts_full.csv")
print(f"\n{'='*80}")
print(f"USAGE EXAMPLES (P0-FIX Applied):")
print(f"{'='*80}")
print(f"# Load")
print(f"adata = sc.read_h5ad('merged_scanvi_L2_prod_v2_5_5.h5ad')")
print(f"")
print(f"# Default UMAP (no basis parameter needed)")
print(f"sc.pl.umap(adata, color=['L2_scanvi_pred', 'data_source'])")
print(f"")
print(f"# Visualize any gene (use layer parameter)")
print(f"sc.pl.umap(adata, color=['CD19', 'MS4A1'], layer='counts')")
print(f"")
print(f"# Or log-transform for better viz")
print(f"adata.layers['log1p_viz'] = sc.pp.log1p(adata.layers['counts'].copy(), copy=True)")
print(f"sc.pl.umap(adata, color=['CD19'], layer='log1p_viz', vmin=0, vmax=4)")
print(f"")
print(f"# DE analysis on full genes")
print(f"# (Need to construct temporary full-gene object if needed)")
