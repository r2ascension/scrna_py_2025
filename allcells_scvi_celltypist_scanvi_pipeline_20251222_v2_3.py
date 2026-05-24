#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
All Cells Analysis: scVI → CellTypist → Dual scANVI Pipeline v2.2.0

Author: Clinical-Bioinformatics Team
Date: 2025-12-22
Version: 2.2.0 (Dual scANVI: existing cell_type + CellTypist)

⭐ What's New in v2.2.0:
========================
NEW FEATURE:
- **Dual scANVI Pipeline**: Train two scANVI models
  1. scANVI-Existing: Uses existing cell_type metadata (no confidence filtering)
  2. scANVI-CellTypist: Uses CellTypist predictions (with confidence filtering)
- Compare both annotation strategies
- Each has separate neighbors graphs, UMAPs, and predictions

All v2.1.3.1 fixes included:
- HVG consistency (P0-1)
- CellTypist memory efficiency (P0-2)
- Gene name strategy with symbol_base (P0-3)
- True semi-supervised learning (P1-1)
- Graph separation (P1-2)
- NaN assertions (P1-3)
- scipy sparse indexing hotfix

Pipeline Overview:
==================
Input: adata_bbknn_annotated_corrected.h5ad
  ↓
Step 1: Data Preparation
  ↓
Step 2: scVI Integration
  ↓
Step 3: CellTypist Annotation
  ↓
Step 4A: scANVI-Existing (new!)
  - Use existing cell_type metadata
  - No confidence filtering
  - Save as scanvi_predictions_existing
  ↓
Step 4B: scANVI-CellTypist
  - Use CellTypist predictions
  - With confidence filtering
  - Save as scanvi_predictions_celltypist
  ↓
Step 5: Comparison & Visualization
  ↓
Output: Fully annotated h5ad with both scANVI results

Requirements:
=============
pip install "scvi-tools>=1.1.4" "lightning>=2.2,<3" scanpy celltypist mygene
"""

# ==============================================================================
# Step 0: Import Libraries and Configuration
# ==============================================================================

import sys
import os
from pathlib import Path
import warnings
import json
import time
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import gc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scanpy as sc
import scvi
import celltypist
from celltypist import models

import mygene

warnings.filterwarnings('ignore')

print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"celltypist: {celltypist.__version__}")
print(f"Python: {sys.version}")

import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")


# ==============================================================================
# Configuration Section
# ==============================================================================

# Input/Output
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected_filtered.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1207/allcells_scvi_analysis"
CELLTYPIST_MODEL_PATH = "/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl"

# Keys
BATCH_KEY = "dataset"
TISSUE_KEY = "tissue_sampling_method"
EXISTING_CELLTYPE_KEY = "cell_type"  # ⭐ Used for scANVI-Existing

# Model loading
LOAD_SCVI_IF_EXISTS = True
LOAD_SCANVI_EXISTING_IF_EXISTS = True   # ⭐ NEW
LOAD_SCANVI_CELLTYPIST_IF_EXISTS = True  # ⭐ NEW

# HVG
USE_HVG_FOR_SCVI = True
N_HVG_SCVI = 4000
HVG_FLAVOR = "seurat_v3"

# ⭐ Semi-supervised for CellTypist-based scANVI only
LOW_CONFIDENCE_THRESHOLD = 0.5
RARE_TYPE_THRESHOLD = 10

# scVI
SCVI_N_LATENT = 150
SCVI_N_LAYERS = 4
SCVI_DROPOUT_RATE = 0.1
SCVI_MAX_EPOCHS = 400
SCVI_LEARNING_RATE = 1e-3
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# scANVI
SCANVI_MAX_EPOCHS = 200
SCANVI_LEARNING_RATE = 1e-3
SCANVI_EARLY_STOPPING_PATIENCE = 30

# UMAP
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
UMAP_N_NEIGHBORS = 50

# Visualization
DPI = 300
FIGURE_FORMAT = 'pdf'

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
scvi.settings.seed = RANDOM_SEED

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)
(output_dir / "models").mkdir(exist_ok=True)

print(f"\nOutput directory: {output_dir}")

# Save configuration
config = {
    'version': '2.2.0',
    'input_h5ad': str(INPUT_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'existing_celltype_key': EXISTING_CELLTYPE_KEY,
    'use_hvg_for_scvi': USE_HVG_FOR_SCVI,
    'n_hvg_scvi': N_HVG_SCVI,
    'low_confidence_threshold': LOW_CONFIDENCE_THRESHOLD,
    'rare_type_threshold': RARE_TYPE_THRESHOLD,
    'dual_scanvi': True,  # ⭐ NEW
    'random_seed': RANDOM_SEED,
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'pipeline_config_v2.2.0.json', 'w') as f:
    json.dump(config, f, indent=2)

print("✓ Configuration saved")


# ==============================================================================
# Step 1: Data Loading and Preparation
# ==============================================================================

print("\n" + "="*70)
print("Step 1: Data Loading and Preparation")
print("="*70)

print(f"\nLoading data from: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"\nData loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found")

if EXISTING_CELLTYPE_KEY not in adata.obs.columns:
    raise ValueError(f"Existing cell type key '{EXISTING_CELLTYPE_KEY}' not found")

print(f"\nBatch distribution ({BATCH_KEY}):")
print(adata.obs[BATCH_KEY].value_counts())

print(f"\nExisting cell type distribution ({EXISTING_CELLTYPE_KEY}):")
print(adata.obs[EXISTING_CELLTYPE_KEY].value_counts())

# Extract raw counts
print("\n" + "-"*70)
print("Extracting raw counts")
print("-"*70)

if 'counts' in adata.layers:
    print("✓ Found raw counts in layers['counts']")
elif hasattr(adata, 'raw') and adata.raw is not None:
    adata.layers['counts'] = adata.raw.X.copy()
    print("✓ Extracted counts from .raw.X")
else:
    raise ValueError("Cannot find raw counts")

# Set adata.raw.X
print("\nSetting adata.raw.X (shared memory)...")
adata_raw = sc.AnnData(
    X=adata.layers["counts"],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
adata.raw = adata_raw
print("✓ adata.raw.X set")

# Gene name handling
print("\n" + "-"*70)
print("Gene Name Handling")
print("-"*70)

sample_gene = str(adata.var_names[0])
print(f"Sample gene: {sample_gene}")

if 'symbol' not in adata.var.columns:
    if sample_gene.startswith('ENSG'):
        print("✓ Detected ENSEMBL IDs")
        
        symbol_candidates = ['gene_symbols', 'feature_name', 'gene_name']
        existing_symbol_col = None
        for col in symbol_candidates:
            if col in adata.var.columns:
                existing_symbol_col = col
                break
        
        if existing_symbol_col:
            adata.var['symbol'] = adata.var[existing_symbol_col].astype(str)
            print(f"  Using existing column: {existing_symbol_col}")
        else:
            print("  Using mygene for conversion...")
            mg = mygene.MyGeneInfo()
            ensembl_ids = adata.var_names.tolist()
            
            query_result = mg.querymany(
                ensembl_ids,
                scopes='ensembl.gene',
                fields='symbol',
                species='human',
                returnall=True,
                verbose=False
            )
            
            ensembl_to_symbol = {}
            for item in query_result['out']:
                if 'symbol' in item:
                    ensembl_to_symbol[item['query']] = item['symbol']
            
            adata.var['symbol'] = [ensembl_to_symbol.get(e, e) for e in adata.var_names]
            print(f"  Converted: {len(ensembl_to_symbol):,}/{len(ensembl_ids):,}")
    else:
        adata.var['symbol'] = adata.var_names.astype(str)
else:
    print("✓ Symbol column exists")

adata.var['symbol_base'] = adata.var['symbol'].str.replace(r'-\d+$', '', regex=True)
print("✓ Created symbol_base column")

print("\n✓ Data preparation complete")


# ==============================================================================
# Step 2: scVI Integration
# ==============================================================================

print("\n" + "="*70)
print("Step 2: scVI Integration")
print("="*70)

scvi_model_dir = output_dir / "models" / "scvi_model"
scvi_model_exists = scvi_model_dir.exists() and (scvi_model_dir / "model.pt").exists()
hvg_file = scvi_model_dir / "hvg_genes.txt"

hvg_indices = None
actual_hvg_count = adata.n_vars
hvg_selection_method = "all_genes"
hvg_genes_list = None

if USE_HVG_FOR_SCVI:
    print(f"\n⭐ Using HVG ({N_HVG_SCVI} genes)")
    
    # Load HVG from file if exists
    if scvi_model_exists and hvg_file.exists() and LOAD_SCVI_IF_EXISTS:
        print("\nLoading HVG from file...")
        try:
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_genes_in_data = [g for g in hvg_genes_saved if g in set(adata.var_names)]
            
            if len(hvg_genes_in_data) >= len(hvg_genes_saved) * 0.95:
                hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_in_data))
                hvg_genes_list = hvg_genes_in_data
                actual_hvg_count = len(hvg_genes_list)
                hvg_selection_method = "loaded_from_file"
                print(f"✓ Loaded {actual_hvg_count:,} HVG genes")
            else:
                hvg_genes_list = None
        except Exception as e:
            print(f"⚠️  Failed to load HVG: {e}")
            hvg_genes_list = None
    
    # Compute HVG if not loaded
    if hvg_genes_list is None:
        print("\nComputing HVG...")
        try:
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                batch_key=BATCH_KEY, flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "batch-aware"
        except:
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "non-batch-aware"
        
        hvg_mask = adata.var["highly_variable"].values
        hvg_indices = np.flatnonzero(hvg_mask)
        hvg_genes_list = adata.var_names[hvg_indices].tolist()
        actual_hvg_count = len(hvg_genes_list)
        print(f"✓ Selected {actual_hvg_count:,} HVG")
    
    # Construct scVI AnnData
    adata_scvi = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[[BATCH_KEY]].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scvi.var_names = pd.Index(hvg_genes_list)
    print(f"✓ scVI AnnData: {adata_scvi.n_obs:,} × {adata_scvi.n_vars:,}")
else:
    adata_scvi = adata.copy()
    adata_scvi.X = adata.layers['counts'].copy()

# Load or train scVI
print("\n" + "-"*70)
print("scVI Model: Load or Train")
print("-"*70)

if scvi_model_exists and LOAD_SCVI_IF_EXISTS:
    print(f"✓ Loading existing scVI model...")
    try:
        scvi_model = scvi.model.SCVI.load(scvi_model_dir, adata=adata_scvi)
        print("✓ scVI model loaded")
    except Exception as e:
        print(f"⚠️  Load failed: {e}")
        scvi_model_exists = False

if not scvi_model_exists or not LOAD_SCVI_IF_EXISTS:
    print("Training new scVI model...")
    
    scvi.model.SCVI.setup_anndata(adata_scvi, layer=None, batch_key=BATCH_KEY)
    scvi_model = scvi.model.SCVI(
        adata_scvi, n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS,
        dropout_rate=SCVI_DROPOUT_RATE, gene_likelihood="nb"
    )
    
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'early_stopping': SCVI_EARLY_STOPPING,
        'early_stopping_patience': SCVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {'lr': SCVI_LEARNING_RATE},
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scvi_model.train(**train_kwargs)
    print("✓ scVI training complete")
    
    scvi_model.save(scvi_model_dir, overwrite=True)
    print(f"✓ Model saved")
    
    if USE_HVG_FOR_SCVI and hvg_genes_list:
        pd.Series(hvg_genes_list).to_csv(hvg_file, index=False, header=False)
        print("✓ HVG list saved")

# Generate latent
print("\nGenerating scVI latent...")
latent_df = pd.DataFrame(
    scvi_model.get_latent_representation(),
    index=adata_scvi.obs_names
)
latent_aligned = latent_df.reindex(adata.obs_names)
if latent_aligned.isna().any().any():
    raise RuntimeError("scVI latent reindex produced NaNs")
adata.obsm['X_scvi'] = latent_aligned.to_numpy().astype(np.float32)
print(f"✓ scVI latent: {adata.obsm['X_scvi'].shape}")

del adata_scvi
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scvi)...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS, 
                key_added='neighbors_scvi', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, 
           neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scvi', 
             neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
print(f"✓ Leiden: {adata.obs['leiden_scvi'].nunique()} clusters")

print("✓ scVI integration complete")


# ==============================================================================
# Step 3: CellTypist Annotation
# ==============================================================================

print("\n" + "="*70)
print("Step 3: CellTypist Annotation")
print("="*70)

print(f"\nLoading CellTypist model...")
celltypist_model = models.Model.load(CELLTYPIST_MODEL_PATH)

# Get model features
model_features = None
for attr in ['features', 'genes', 'var_names']:
    if hasattr(celltypist_model, attr):
        mf = getattr(celltypist_model, attr)
        if mf is not None and len(mf) > 100:
            model_features = pd.Index(mf).astype(str)
            break

if model_features is None:
    raise RuntimeError("Cannot retrieve CellTypist features")

# Match genes
print("\nMatching genes using symbol_base...")
common_genes = adata.var['symbol_base'].isin(model_features)
n_common = common_genes.sum()
print(f"  Overlapping: {n_common:,}/{len(model_features)}")

if n_common < 500:
    raise RuntimeError(f"Too few overlapping genes: {n_common}")

# Construct CellTypist AnnData (memory efficient)
print(f"\nConstructing CellTypist AnnData (only {n_common:,} genes)...")
X_subset = adata.layers["counts"][:, common_genes.values]

adata_celltypist = sc.AnnData(
    X=X_subset.copy(),
    obs=adata.obs.copy(),
    var=adata.var.loc[common_genes].copy()
)
adata_celltypist.var_names = adata.var.loc[common_genes, 'symbol_base'].values

sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
sc.pp.log1p(adata_celltypist)
print("✓ Data prepared")

# Run CellTypist
print("\nRunning CellTypist prediction...")
predictions = celltypist.annotate(adata_celltypist, model=celltypist_model, majority_voting=True)

# Extract results
predicted_labels = predictions.predicted_labels
pred_df = predicted_labels.reindex(adata.obs_names)
if pred_df.isna().any().any():
    raise RuntimeError("CellTypist reindex produced NaNs")

adata.obs['celltypist_predicted'] = pred_df.predicted_labels.astype(str).values
if 'majority_voting' in pred_df.columns:
    adata.obs['celltypist_majority_voting'] = pred_df.majority_voting.astype(str).values
else:
    adata.obs['celltypist_majority_voting'] = pred_df.predicted_labels.astype(str).values

# Confidence
conf_candidates = ['conf_score', 'confidence', 'confidence_score', 'prob']
confidence_column = next((col for col in conf_candidates if col in predicted_labels.columns), None)
if confidence_column:
    adata.obs['celltypist_conf_score'] = pred_df[confidence_column].reindex(adata.obs_names).values
else:
    adata.obs['celltypist_conf_score'] = 1.0

del adata_celltypist, X_subset
gc.collect()

print(f"\nCellTypist summary:")
print(f"  Unique types: {adata.obs['celltypist_majority_voting'].nunique()}")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")

print("✓ CellTypist complete")


# ==============================================================================
# Step 4A: scANVI-Existing (Using existing cell_type metadata)
# ==============================================================================

print("\n" + "="*70)
print("Step 4A: scANVI-Existing (Using cell_type metadata)")
print("="*70)

print(f"\n⭐ Training scANVI using '{EXISTING_CELLTYPE_KEY}' metadata")
print(f"   No confidence filtering applied")

scanvi_existing_model_dir = output_dir / "models" / "scanvi_existing_model"
scanvi_existing_exists = scanvi_existing_model_dir.exists() and (scanvi_existing_model_dir / "model.pt").exists()

# Prepare data
print("\nPreparing scANVI-Existing data...")
if USE_HVG_FOR_SCVI:
    if hvg_indices is None:
        if 'highly_variable' in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            hvg_indices = np.flatnonzero(hvg_mask)
        elif hvg_file.exists():
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_saved))
    
    adata_scanvi_existing = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[[BATCH_KEY, EXISTING_CELLTYPE_KEY]].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scanvi_existing.var_names = adata.var_names[hvg_indices]
else:
    adata_scanvi_existing = sc.AnnData(
        X=adata.layers["counts"].copy(),
        obs=adata.obs[[BATCH_KEY, EXISTING_CELLTYPE_KEY]].copy(),
        var=adata.var.copy()
    )

print(f"✓ scANVI-Existing AnnData: {adata_scanvi_existing.n_obs:,} × {adata_scanvi_existing.n_vars:,}")

# Ensure 'Unknown' category exists
labels_existing = adata_scanvi_existing.obs[EXISTING_CELLTYPE_KEY].astype('category')
if 'Unknown' not in labels_existing.cat.categories:
    labels_existing = labels_existing.cat.add_categories(['Unknown'])
    adata_scanvi_existing.obs[EXISTING_CELLTYPE_KEY] = labels_existing

# Load or train
print("\n" + "-"*70)
print("scANVI-Existing Model: Load or Train")
print("-"*70)

if scanvi_existing_exists and LOAD_SCANVI_EXISTING_IF_EXISTS:
    print("Loading existing scANVI-Existing model...")
    try:
        scanvi_existing_model = scvi.model.SCANVI.load(scanvi_existing_model_dir, adata=adata_scanvi_existing)
        print("✓ Model loaded")
    except Exception as e:
        print(f"⚠️  Load failed: {e}")
        scanvi_existing_exists = False

if not scanvi_existing_exists or not LOAD_SCANVI_EXISTING_IF_EXISTS:
    print("Training new scANVI-Existing model...")
    
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_existing, layer=None, batch_key=BATCH_KEY,
        labels_key=EXISTING_CELLTYPE_KEY, unlabeled_category="Unknown"
    )
    
    scanvi_existing_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category="Unknown",
        adata=adata_scanvi_existing, labels_key=EXISTING_CELLTYPE_KEY
    )
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'early_stopping': True,
        'early_stopping_patience': SCANVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {'lr': SCANVI_LEARNING_RATE},
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scanvi_existing_model.train(**train_kwargs)
    print("✓ Training complete")
    
    scanvi_existing_model.save(scanvi_existing_model_dir, overwrite=True)
    print("✓ Model saved")

# Generate results
print("\nGenerating scANVI-Existing results...")

# Latent
latent_scanvi_existing_df = pd.DataFrame(
    scanvi_existing_model.get_latent_representation(),
    index=adata_scanvi_existing.obs_names
)
latent_existing_aligned = latent_scanvi_existing_df.reindex(adata.obs_names)
if latent_existing_aligned.isna().any().any():
    raise RuntimeError("scANVI-Existing latent reindex NaNs")
adata.obsm['X_scanvi_existing'] = latent_existing_aligned.to_numpy().astype(np.float32)

# Predictions
predictions_existing_series = pd.Series(
    scanvi_existing_model.predict(),
    index=adata_scanvi_existing.obs_names
)
predictions_existing_aligned = predictions_existing_series.reindex(adata.obs_names)
if predictions_existing_aligned.isna().any():
    raise RuntimeError("scANVI-Existing predictions reindex NaNs")
adata.obs['scanvi_predictions_existing'] = predictions_existing_aligned.astype(str).values

# Probabilities
probs_existing = np.asarray(scanvi_existing_model.predict(soft=True))
probs_existing_df = pd.DataFrame(probs_existing, index=adata_scanvi_existing.obs_names)
probs_existing_aligned = probs_existing_df.reindex(adata.obs_names)
if probs_existing_aligned.isna().any().any():
    raise RuntimeError("scANVI-Existing probs reindex NaNs")
adata.obsm['scanvi_probabilities_existing'] = probs_existing_aligned.to_numpy().astype(np.float32)

del adata_scanvi_existing
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scanvi_existing)...")
sc.pp.neighbors(adata, use_rep='X_scanvi_existing', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_existing'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_existing',
             neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print(f"\n✓ scANVI-Existing complete")
print(f"  Cell types: {adata.obs['scanvi_predictions_existing'].nunique()}")
print(f"  Top 5:")
for ct, count in adata.obs['scanvi_predictions_existing'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")


# ==============================================================================
# Step 4B: scANVI-CellTypist (Using CellTypist with filtering)
# ==============================================================================

print("\n" + "="*70)
print("Step 4B: scANVI-CellTypist (Using CellTypist predictions)")
print("="*70)

print(f"\n⭐ Training scANVI using CellTypist predictions")
print(f"   With confidence filtering (threshold: {LOW_CONFIDENCE_THRESHOLD})")
print(f"   With rare type filtering (threshold: {RARE_TYPE_THRESHOLD} cells)")

# Create semi-supervised labels
print("\nCreating semi-supervised labels...")
labels_celltypist = adata.obs['celltypist_majority_voting'].astype(str).copy()
conf = adata.obs['celltypist_conf_score'].astype(float)

n_original = len(labels_celltypist)
n_low_conf = (conf < LOW_CONFIDENCE_THRESHOLD).sum()
labels_celltypist[conf < LOW_CONFIDENCE_THRESHOLD] = "Unknown"

value_counts = labels_celltypist.value_counts()
rare_types = value_counts[value_counts < RARE_TYPE_THRESHOLD].index
rare_types = rare_types[rare_types != "Unknown"]
n_rare = labels_celltypist.isin(rare_types).sum()
labels_celltypist[labels_celltypist.isin(rare_types)] = "Unknown"

n_unknown = (labels_celltypist == "Unknown").sum()
n_labeled = n_original - n_unknown

print(f"  Total: {n_original:,}")
print(f"  Low confidence → Unknown: {n_low_conf:,} ({n_low_conf/n_original*100:.1f}%)")
print(f"  Rare types → Unknown: {n_rare:,} ({n_rare/n_original*100:.1f}%)")
print(f"  Final labeled: {n_labeled:,} ({n_labeled/n_original*100:.1f}%)")
print(f"  Final Unknown: {n_unknown:,} ({n_unknown/n_original*100:.1f}%)")

adata.obs['scanvi_label_celltypist'] = pd.Categorical(labels_celltypist)

scanvi_celltypist_model_dir = output_dir / "models" / "scanvi_celltypist_model"
scanvi_celltypist_exists = scanvi_celltypist_model_dir.exists() and (scanvi_celltypist_model_dir / "model.pt").exists()

# Prepare data
print("\nPreparing scANVI-CellTypist data...")
if USE_HVG_FOR_SCVI:
    adata_scanvi_celltypist = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[[BATCH_KEY, 'scanvi_label_celltypist']].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scanvi_celltypist.var_names = adata.var_names[hvg_indices]
else:
    adata_scanvi_celltypist = sc.AnnData(
        X=adata.layers["counts"].copy(),
        obs=adata.obs[[BATCH_KEY, 'scanvi_label_celltypist']].copy(),
        var=adata.var.copy()
    )

print(f"✓ scANVI-CellTypist AnnData: {adata_scanvi_celltypist.n_obs:,} × {adata_scanvi_celltypist.n_vars:,}")

# Ensure 'Unknown' exists
labels_cat = adata_scanvi_celltypist.obs['scanvi_label_celltypist'].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi_celltypist.obs['scanvi_label_celltypist'] = labels_cat

# Load or train
print("\n" + "-"*70)
print("scANVI-CellTypist Model: Load or Train")
print("-"*70)

if scanvi_celltypist_exists and LOAD_SCANVI_CELLTYPIST_IF_EXISTS:
    print("Loading existing scANVI-CellTypist model...")
    try:
        scanvi_celltypist_model = scvi.model.SCANVI.load(scanvi_celltypist_model_dir, adata=adata_scanvi_celltypist)
        print("✓ Model loaded")
    except Exception as e:
        print(f"⚠️  Load failed: {e}")
        scanvi_celltypist_exists = False

if not scanvi_celltypist_exists or not LOAD_SCANVI_CELLTYPIST_IF_EXISTS:
    print("Training new scANVI-CellTypist model...")
    
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_celltypist, layer=None, batch_key=BATCH_KEY,
        labels_key='scanvi_label_celltypist', unlabeled_category="Unknown"
    )
    
    scanvi_celltypist_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category="Unknown",
        adata=adata_scanvi_celltypist, labels_key='scanvi_label_celltypist'
    )
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'early_stopping': True,
        'early_stopping_patience': SCANVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {'lr': SCANVI_LEARNING_RATE},
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scanvi_celltypist_model.train(**train_kwargs)
    print("✓ Training complete")
    
    scanvi_celltypist_model.save(scanvi_celltypist_model_dir, overwrite=True)
    print("✓ Model saved")

# Generate results
print("\nGenerating scANVI-CellTypist results...")

# Latent
latent_scanvi_celltypist_df = pd.DataFrame(
    scanvi_celltypist_model.get_latent_representation(),
    index=adata_scanvi_celltypist.obs_names
)
latent_celltypist_aligned = latent_scanvi_celltypist_df.reindex(adata.obs_names)
if latent_celltypist_aligned.isna().any().any():
    raise RuntimeError("scANVI-CellTypist latent reindex NaNs")
adata.obsm['X_scanvi_celltypist'] = latent_celltypist_aligned.to_numpy().astype(np.float32)

# Predictions
predictions_celltypist_series = pd.Series(
    scanvi_celltypist_model.predict(),
    index=adata_scanvi_celltypist.obs_names
)
predictions_celltypist_aligned = predictions_celltypist_series.reindex(adata.obs_names)
if predictions_celltypist_aligned.isna().any():
    raise RuntimeError("scANVI-CellTypist predictions reindex NaNs")
adata.obs['scanvi_predictions_celltypist'] = predictions_celltypist_aligned.astype(str).values

# Probabilities
probs_celltypist = np.asarray(scanvi_celltypist_model.predict(soft=True))
probs_celltypist_df = pd.DataFrame(probs_celltypist, index=adata_scanvi_celltypist.obs_names)
probs_celltypist_aligned = probs_celltypist_df.reindex(adata.obs_names)
if probs_celltypist_aligned.isna().any().any():
    raise RuntimeError("scANVI-CellTypist probs reindex NaNs")
adata.obsm['scanvi_probabilities_celltypist'] = probs_celltypist_aligned.to_numpy().astype(np.float32)

del adata_scanvi_celltypist
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scanvi_celltypist)...")
sc.pp.neighbors(adata, use_rep='X_scanvi_celltypist', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_celltypist'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_celltypist',
             neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print(f"\n✓ scANVI-CellTypist complete")
print(f"  Cell types: {adata.obs['scanvi_predictions_celltypist'].nunique()}")
print(f"  Top 5:")
for ct, count in adata.obs['scanvi_predictions_celltypist'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")


# ==============================================================================
# Step 5: Comparison & Final Export
# ==============================================================================

print("\n" + "="*70)
print("Step 5: Comparison & Final Export")
print("="*70)

# Comparison visualization
print("\nGenerating comparison plots...")

fig, axes = plt.subplots(2, 4, figsize=(24, 12))

# Row 1: scANVI-Existing
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_existing'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scANVI-Existing: Batch')
sc.pl.umap(adata, color=EXISTING_CELLTYPE_KEY, ax=axes[0, 1], show=False, title='scANVI-Existing: Input')
sc.pl.umap(adata, color='scanvi_predictions_existing', ax=axes[0, 2], show=False, title='scANVI-Existing: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_existing', ax=axes[0, 3], show=False, title='scANVI-Existing: Leiden')

# Row 2: scANVI-CellTypist
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0], show=False, title='scANVI-CellTypist: Batch')
sc.pl.umap(adata, color='scanvi_label_celltypist', ax=axes[1, 1], show=False, title='scANVI-CellTypist: Input')
sc.pl.umap(adata, color='scanvi_predictions_celltypist', ax=axes[1, 2], show=False, title='scANVI-CellTypist: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_celltypist', ax=axes[1, 3], show=False, title='scANVI-CellTypist: Leiden')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"dual_scanvi_comparison.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.close()

print("✓ Comparison plots saved")

# Set default UMAP to CellTypist-based
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()

# Save final dataset
output_file = output_dir / "adata_allcells_dual_scanvi_final_v2.2.0.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"✓ Saved ({file_size:.2f} GB)")

# Summary
print("\n" + "="*70)
print("DUAL scANVI ANALYSIS SUMMARY")
print("="*70)

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("Dual scANVI Analysis Summary v2.2.0")
summary_lines.append("="*70)

summary_lines.append("\n[ scANVI-Existing (cell_type metadata) ]")
summary_lines.append(f"  Input: All cells with existing labels")
summary_lines.append(f"  Filtering: None")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_existing'].nunique()}")

summary_lines.append("\n[ scANVI-CellTypist (CellTypist predictions) ]")
summary_lines.append(f"  Input: CellTypist predictions")
summary_lines.append(f"  Labeled: {n_labeled:,}/{n_original:,} ({n_labeled/n_original*100:.1f}%)")
summary_lines.append(f"  Unknown: {n_unknown:,}/{n_original:,} ({n_unknown/n_original*100:.1f}%)")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_celltypist'].nunique()}")

summary_lines.append("\n[ Key Outputs ]")
summary_lines.append("  scANVI-Existing:")
summary_lines.append(f"    - scanvi_predictions_existing")
summary_lines.append(f"    - X_scanvi_existing, X_umap_scanvi_existing")
summary_lines.append(f"    - neighbors_scanvi_existing")
summary_lines.append("  scANVI-CellTypist:")
summary_lines.append(f"    - scanvi_predictions_celltypist")
summary_lines.append(f"    - X_scanvi_celltypist, X_umap_scanvi_celltypist")
summary_lines.append(f"    - neighbors_scanvi_celltypist")

summary_lines.append("\n" + "="*70)
summary_lines.append("Analysis complete!")
summary_lines.append("="*70)

summary_text = '\n'.join(summary_lines)
print(summary_text)

summary_file = output_dir / "analysis_summary_v2.2.0.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

print(f"\n✓ Summary saved to: {summary_file}")

print("\n" + "="*70)
print("🎉 DUAL scANVI PIPELINE COMPLETE (v2.2.0)")
print("="*70)

print(f"\n📊 Final Results:")
print(f"  • Dual scANVI annotations available")
print(f"  • Compare: scanvi_predictions_existing vs scanvi_predictions_celltypist")
print(f"  • Both have separate UMAPs and neighbor graphs")
print(f"  • Use either for downstream analysis")

print("\n" + "="*70)