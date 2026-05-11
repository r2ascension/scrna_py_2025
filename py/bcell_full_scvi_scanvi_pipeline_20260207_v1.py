#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
Production Script: MERGED scVI+scANVI Training (Reference + Query)
==============================================================================

Version: 2.5.3-MERGED-PROD (Fixed & Robust)
Date: 2026-02-07
Author: r2end (Refactored)

Key Features:
- STRICT counts layer validation.
- Unified Label Construction (Ref Truth + High-Conf Query Pred).
- 3-Level HVG Fallback (Batch-Seurat -> Seurat -> Cell Ranger).
- Consistent UMAP: Embedding matches the saved operator exactly.
- Full Covariates: MT%, Stress, Cell Cycle, Tissue, Batch.

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
from umap import UMAP  # pip install umap-learn

import torch

# Silence warnings
warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None 

print("=" * 80)
print("MERGED scVI+scANVI Training (v2.5.3-PROD)")
print("=" * 80)

# Check GPU
gpu_available = torch.cuda.is_available()
print(f"GPU available: {gpu_available}")
if gpu_available:
    print(f"GPU Device: {torch.cuda.get_device_name(0)}")

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# --- IO Paths ---
INPUT_MERGED_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3_MERGE_v2_0/reference_plus_query_merged_v2_0.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0207/merged_scanvi_L2_prod_v1"

# --- Column Keys ---
DATASOURCE_KEY = "data_source"         # reference / query
REF_LABEL_KEY = "Cell_Type_L2"         # Reference Ground Truth
QRY_LABEL_KEY = "Cell_Type_L2_final"   # Query Prediction
QRY_CONF_KEY = "mapping_confidence"    # Query Confidence
TRAIN_LABEL_KEY = "Cell_Type_L2_train" # Unified Label Key (Target)
UNLABELED_CATEGORY = "Unknown"

# --- Covariate Keys ---
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"

# --- Query Filter Settings ---
USE_QUERY_CONF_FILTER = True
QUERY_CONF_THR = 0.5  # Only use query labels with conf >= 0.5

# --- HVG Settings ---
N_HVG = 4000
FORCE_MARKERS_IN_HVG = True

# --- Training Params (Matches v2.5.3) ---
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
    """Strictly ensures raw counts exist. Fails if uncertain."""
    if counts_layer in (adata.layers or {}):
        adata.X = adata.layers[counts_layer]
        # Ensure CSR for speed
        if issparse(adata.X) and not isinstance(adata.X, csr_matrix):
            adata.X = adata.X.tocsr()
            adata.layers[counts_layer] = adata.X
        return counts_layer
    else:
        # Fallback check: is .X raw integers?
        if issparse(adata.X):
            data_sample = adata.X.data[:1000]
        else:
            data_sample = adata.X.flat[:1000]
            
        is_int = np.all(np.mod(data_sample, 1) == 0)
        if is_int:
            print("WARNING: layers['counts'] missing, but .X appears to be integer counts. Using .X.")
            adata.layers[counts_layer] = adata.X.copy()
            return counts_layer
        else:
            raise ValueError(
                "CRITICAL ERROR: layers['counts'] not found and .X contains non-integers (log/scaled?). "
                "scVI requires raw counts! Please regenerate merged object."
            )

def ensure_batch_tissue(adata, batch_key, tissue_key):
    # Batch
    if batch_key not in adata.obs.columns:
        adata.obs[batch_key] = "unknown_batch"
    adata.obs[batch_key] = adata.obs[batch_key].astype("category")
    
    # Tissue
    if tissue_key not in adata.obs.columns:
        adata.obs[tissue_key] = "unknown_tissue"
    adata.obs[tissue_key] = adata.obs[tissue_key].astype("category")
    if "unknown_tissue" not in adata.obs[tissue_key].cat.categories:
        adata.obs[tissue_key] = adata.obs[tissue_key].cat.add_categories(["unknown_tissue"])
    adata.obs[tissue_key] = adata.obs[tissue_key].fillna("unknown_tissue")

def ensure_pct_mt(adata, counts_layer="counts"):
    for col in ["pct_counts_mt", "percent.mt", "percent_mt"]:
        if col in adata.obs.columns:
            if col != "pct_counts_mt":
                adata.obs["pct_counts_mt"] = adata.obs[col]
            return
    # Compute
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, layer=counts_layer)

def ensure_stress_score(adata, stress_genes, counts_layer="counts"):
    if "stress_score" in adata.obs.columns: return

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
    if all(k in adata.obs.columns for k in ["S_score", "G2M_score", "phase"]): return

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

def prepare_covariates(adata):
    print("  -> Validating counts...")
    l = ensure_counts_layer(adata, "counts")
    
    print("  -> Batch/Tissue...")
    ensure_batch_tissue(adata, BATCH_KEY, TISSUE_KEY)
    
    print("  -> MT% / Stress / Cell Cycle...")
    ensure_pct_mt(adata, l)
    ensure_stress_score(adata, STRESS_SIGNATURE_GENES, l)
    ensure_cell_cycle_scores(adata, S_GENES, G2M_GENES, l)

# ==============================================================================
# 3. MAIN PIPELINE
# ==============================================================================

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

# --- Step 1: Load Data ---
print("\n[Step 1] Loading MERGED data...")
adata = sc.read_h5ad(INPUT_MERGED_H5AD)
adata.var_names_make_unique()
print(f"Shape: {adata.shape}")

# --- Step 2: Prepare Covariates ---
print("\n[Step 2] Preparing Covariates...")
prepare_covariates(adata)

# --- Step 3: Build Unified Labels (FIXED) ---
print("\n[Step 3] Merging Labels (Ref + High-Conf Query)...")
is_ref = adata.obs[DATASOURCE_KEY] == "reference"
is_qry = adata.obs[DATASOURCE_KEY] == "query"

# Initialize with correct string type
adata.obs[TRAIN_LABEL_KEY] = pd.Series(UNLABELED_CATEGORY, index=adata.obs_names).astype("string")

# 1. Fill Reference
adata.obs.loc[is_ref, TRAIN_LABEL_KEY] = adata.obs.loc[is_ref, REF_LABEL_KEY].astype("string")

# 2. Fill Query (Simplified & Safe)
conf_full = pd.to_numeric(adata.obs[QRY_CONF_KEY], errors="coerce").fillna(0.0)
valid_qry_idx = is_qry & (conf_full >= QUERY_CONF_THR)
adata.obs.loc[valid_qry_idx, TRAIN_LABEL_KEY] = adata.obs.loc[valid_qry_idx, QRY_LABEL_KEY].astype("string")

print(f"  -> Query retained: {valid_qry_idx.sum()} / {is_qry.sum()} (Thr={QUERY_CONF_THR})")

# 3. Clean Categories (NA -> Unknown)
adata.obs[TRAIN_LABEL_KEY] = adata.obs[TRAIN_LABEL_KEY].fillna(UNLABELED_CATEGORY)
adata.obs[TRAIN_LABEL_KEY] = adata.obs[TRAIN_LABEL_KEY].astype("category")
if UNLABELED_CATEGORY not in adata.obs[TRAIN_LABEL_KEY].cat.categories:
    adata.obs[TRAIN_LABEL_KEY] = adata.obs[TRAIN_LABEL_KEY].cat.add_categories([UNLABELED_CATEGORY])

print("  -> Label Counts (Top 10):")
print(adata.obs[TRAIN_LABEL_KEY].value_counts().head(10))

# --- Step 4: HVG (3-Level Fallback) ---
print("\n[Step 4] Selecting HVGs (Robust)...")
hvg_method = "unknown"
try:
    print("  -> Attempting: batch-aware seurat_v3")
    sc.pp.highly_variable_genes(
        adata, layer="counts", n_top_genes=N_HVG, batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
    )
    hvg_method = "batch_seurat_v3"
except Exception as e1:
    try:
        print(f"  -> Failed ({str(e1)[:50]}). Attempting: standard seurat_v3")
        sc.pp.highly_variable_genes(
            adata, layer="counts", n_top_genes=N_HVG, flavor="seurat_v3", subset=False
        )
        hvg_method = "standard_seurat_v3"
    except Exception as e2:
        print(f"  -> Failed ({str(e2)[:50]}). Fallback: cell_ranger")
        sc.pp.highly_variable_genes(
            adata, layer="counts", n_top_genes=N_HVG, flavor="cell_ranger", subset=False
        )
        hvg_method = "cell_ranger"

print(f"  -> Method used: {hvg_method}")

if FORCE_MARKERS_IN_HVG:
    n_added = 0
    for g in FORCED_MARKERS:
        if g in adata.var_names:
            if not adata.var.loc[g, "highly_variable"]:
                adata.var.loc[g, "highly_variable"] = True
                n_added += 1
    print(f"  -> Forced {n_added} markers into HVG list.")

# Save HVG list
hvg_genes = adata.var_names[adata.var["highly_variable"]]
with open(output_dir / "hvg_genes.txt", "w") as f:
    f.write("\n".join(hvg_genes))

adata = adata[:, adata.var["highly_variable"]].copy()
gc.collect()

# --- Step 5: scVI Setup & Assertion ---
print("\n[Step 5] scVI Setup...")

# Critical Assertions
assert "counts" in adata.layers, "Lost counts layer during HVG subset!"
assert TRAIN_LABEL_KEY in adata.obs.columns
assert UNLABELED_CATEGORY in adata.obs[TRAIN_LABEL_KEY].cat.categories

setup_kwargs = {
    "layer": "counts", # Explicitly point to the layer
    "batch_key": BATCH_KEY,
    "continuous_covariate_keys": ["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
    "categorical_covariate_keys": [TISSUE_KEY]
}
# Filter None
setup_kwargs = {k: v for k, v in setup_kwargs.items() if v is not None}

scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)

model = scvi.model.SCVI(
    adata,
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
model.train(**train_kwargs)

try:
    hist = model.history['reconstruction_loss_train'].iloc[-1]
    print(f"  -> scVI Final Loss: {hist:.4f}")
except:
    pass

# --- Step 6: scANVI Training ---
print("\n[Step 6] scANVI Training (Semi-supervised)...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model,
    adata=adata,
    labels_key=TRAIN_LABEL_KEY,
    unlabeled_category=UNLABELED_CATEGORY
)

scanvi_model.train(
    max_epochs=MAX_EPOCHS_SCANVI,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    early_stopping_patience=20,
    plan_kwargs={"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
    **{k:v for k,v in train_kwargs.items() if k in ["accelerator", "devices"]}
)
print("  -> scANVI Training Complete.")

# --- Step 7: Export Representations ---
print("\n[Step 7] Exporting Data...")
adata.obsm["X_scANVI_L2"] = scanvi_model.get_latent_representation(adata)
adata.obs["L2_scanvi_pred"] = scanvi_model.predict(adata)

# Soft proba (Float32 to save space)
proba = scanvi_model.predict(adata, soft=True).astype(np.float32)
adata.obsm["proba_L2_scanvi"] = proba
adata.obs["L2_scanvi_confidence"] = np.max(proba, axis=1)

# --- Step 8: UMAP (Consistent Operator) ---
print("\n[Step 8] Fitting & Saving UMAP Operator (Consistency Fix)...")
# We use umap-learn DIRECTLY here to ensure the embedding matches the saved operator.
# Using 'euclidean' as standard default for Latent spaces in v2.5.3 style.

umap_op = UMAP(
    n_neighbors=30,
    n_components=2,
    min_dist=0.5,
    spread=1.0,
    metric="euclidean",
    random_state=RANDOM_SEED
)

# Fit and Transform in one go
X_umap = umap_op.fit_transform(adata.obsm["X_scANVI_L2"])
adata.obsm["X_umap_scanvi"] = X_umap

# Save the fitted operator
joblib.dump(umap_op, output_dir / "umap_operator.joblib")
print("  -> UMAP operator saved. X_umap_scanvi updated.")

# --- Step 9: Save Outputs ---
print("\n[Step 9] Saving Results...")

# 1. Model
scanvi_model.save(output_dir / "scanvi_model", overwrite=True)

# 2. Config
config = {
    "timestamp": datetime.now().isoformat(),
    "n_hvg": int(adata.n_vars),
    "hvg_method": hvg_method,
    "scvi_latent": SCVI_N_LATENT,
    "covariates_continuous": setup_kwargs.get("continuous_covariate_keys"),
    "covariates_categorical": setup_kwargs.get("categorical_covariate_keys"),
    "training_stats": {
        "train_labels": adata.obs[TRAIN_LABEL_KEY].value_counts().to_dict(),
        "query_confidence_threshold": QUERY_CONF_THR
    }
}
with open(output_dir / "training_config.json", "w") as f:
    json.dump(config, f, indent=2)

# 3. H5AD
out_path = output_dir / "merged_scanvi_L2_prod.h5ad"
adata.write_h5ad(out_path, compression="gzip")

print(f"\n[SUCCESS] Pipeline Finished. Output: {OUTPUT_DIR}")