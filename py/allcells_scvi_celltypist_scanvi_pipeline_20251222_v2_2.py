#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
All Cells Analysis: scVI → CellTypist → scANVI Pipeline v2.1.3.1 (Hotfix)

Author: Clinical-Bioinformatics Team
Date: 2025-12-22
Version: 2.1.3.1 (Hotfix: scipy sparse indexing)

⭐ What's New in v2.1.3.1:
========================
HOTFIX:
- Fixed scipy sparse matrix indexing error with pandas Series
  (convert to numpy array: common_genes.values)

v2.1.3 Changes:
========================
CRITICAL FIXES (P0 - Must Fix):
1. **HVG Consistency** - Save/load HVG gene list for model reuse
2. **CellTypist Memory** - Only copy model features (not all genes)
3. **Gene Name Strategy** - Use symbol_base column, don't rewrite var_names

IMPORTANT FIXES (P1 - Strongly Recommended):
4. **True Semi-supervised** - Set low-confidence + rare types to Unknown
5. **Graph Separation** - Use neighbors_key for scVI/scANVI graphs
6. **NaN Assertions** - Check for missing values after reindex

OPTIMIZATIONS (P2 - Nice to Have):
7. **Explicit lr** - Pass learning rate in plan_kwargs
8. **Float32 storage** - Reduce file size for probabilities/latent
9. **Remove unused** - Clean up seaborn import

Pipeline Overview:
==================
Input: adata_bbknn_annotated_corrected.h5ad
  ↓
Step 1: Data Preparation
  - Load and verify raw counts
  - Set adata.raw.X (shared memory)
  - Handle gene names via symbol column
  ↓
Step 2: scVI Integration (HVG optimized + reusable)
  - Load HVG list if model exists, or compute and save
  - Train/load scVI model
  - Save neighbors as "neighbors_scvi"
  ↓
Step 3: CellTypist Annotation (memory efficient)
  - Only use model feature genes
  - Match via symbol_base
  ↓
Step 4: scANVI Refinement (true semi-supervised)
  - Set low-confidence + rare types to Unknown
  - Train/load scANVI model
  - Save neighbors as "neighbors_scanvi"
  ↓
Output: Fully annotated h5ad

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
matplotlib.use('Agg')  # Non-interactive backend for HPC
import matplotlib.pyplot as plt

# Single-cell analysis
import scanpy as sc
import scvi
import celltypist
from celltypist import models

# Gene name conversion
import mygene

warnings.filterwarnings('ignore')

# Check versions
print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"celltypist: {celltypist.__version__}")
print(f"Python: {sys.version}")

# Check GPU availability
import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")
    print(f"GPU count: {torch.cuda.device_count()}")


# ==============================================================================
# Configuration Section
# ==============================================================================

# Input/Output
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1207/allcells_scvi_analysis"
CELLTYPIST_MODEL_PATH = "/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl"

# Analysis parameters
BATCH_KEY = "dataset"
TISSUE_KEY = "tissue_sampling_method"
EXISTING_CELLTYPE_KEY = "cell_type"

# Model loading flags
LOAD_SCVI_IF_EXISTS = True
LOAD_SCANVI_IF_EXISTS = True

# HVG selection for scVI
USE_HVG_FOR_SCVI = True
N_HVG_SCVI = 4000
HVG_FLAVOR = "seurat_v3"

# ⭐ NEW in v2.1.3: Semi-supervised parameters
LOW_CONFIDENCE_THRESHOLD = 0.5  # CellTypist confidence threshold
RARE_TYPE_THRESHOLD = 10  # Minimum cells per type

# scVI parameters
SCVI_N_LATENT = 150
SCVI_N_LAYERS = 4
SCVI_DROPOUT_RATE = 0.1
SCVI_MAX_EPOCHS = 400
SCVI_LEARNING_RATE = 1e-3  # ⭐ Explicit lr
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# scANVI parameters
SCANVI_MAX_EPOCHS = 200
SCANVI_LEARNING_RATE = 1e-3  # ⭐ Explicit lr
SCANVI_EARLY_STOPPING_PATIENCE = 30

# UMAP parameters
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
    'version': '2.1.3.1',
    'input_h5ad': str(INPUT_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'tissue_key': TISSUE_KEY,
    'use_hvg_for_scvi': USE_HVG_FOR_SCVI,
    'n_hvg_scvi': N_HVG_SCVI,
    'low_confidence_threshold': LOW_CONFIDENCE_THRESHOLD,
    'rare_type_threshold': RARE_TYPE_THRESHOLD,
    'scvi_n_latent': SCVI_N_LATENT,
    'scvi_learning_rate': SCVI_LEARNING_RATE,
    'scvi_max_epochs': SCVI_MAX_EPOCHS,
    'scanvi_learning_rate': SCANVI_LEARNING_RATE,
    'scanvi_max_epochs': SCANVI_MAX_EPOCHS,
    'load_scvi_if_exists': LOAD_SCVI_IF_EXISTS,
    'load_scanvi_if_exists': LOAD_SCANVI_IF_EXISTS,
    'random_seed': RANDOM_SEED,
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'pipeline_config_v2.1.3.1.json', 'w') as f:
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
print(f"  Memory: ~{adata.X.data.nbytes / 1e9:.2f} GB" if issparse(adata.X) else f"  Memory: ~{adata.X.nbytes / 1e9:.2f} GB")

print(f"\n.obs columns: {list(adata.obs.columns[:10])}...")
print(f".layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f".obsm: {list(adata.obsm.keys()) if adata.obsm else 'None'}")

# Verify keys
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found in .obs")

print(f"\nBatch distribution ({BATCH_KEY}):")
print(adata.obs[BATCH_KEY].value_counts())

if TISSUE_KEY in adata.obs.columns:
    print(f"\nTissue distribution ({TISSUE_KEY}):")
    print(adata.obs[TISSUE_KEY].value_counts())

# Extract raw counts
print("\n" + "-"*70)
print("Extracting raw counts for scVI/scANVI")
print("-"*70)

if 'counts' in adata.layers:
    print("✓ Found raw counts in layers['counts']")
    
    sample_values = adata.layers['counts'][:100, :100]
    if issparse(sample_values):
        sample_values = sample_values.toarray()
    
    has_decimals = np.any(sample_values != sample_values.astype(int))
    if has_decimals:
        print("⚠️  Warning: 'counts' layer contains non-integer values")
    else:
        print("✓ Verified: counts are integer values")
    
elif hasattr(adata, 'raw') and adata.raw is not None:
    print("⚠️  'counts' not in layers, extracting from .raw")
    adata.layers['counts'] = adata.raw.X.copy()
    print("✓ Extracted counts from .raw.X")
else:
    print("⚠️  Checking if .X is raw counts...")
    sample_values = adata.X[:100, :100]
    if issparse(sample_values):
        sample_values = sample_values.toarray()
    
    has_decimals = np.any(sample_values != sample_values.astype(int))
    max_val = np.max(sample_values)
    
    if not has_decimals and max_val > 10:
        print("✓ .X appears to be raw counts")
        adata.layers['counts'] = adata.X.copy()
    else:
        raise ValueError("Cannot find raw counts in data")

# Set adata.raw.X (shared memory)
print("\n" + "-"*70)
print("Setting adata.raw.X (shared memory)")
print("-"*70)

print("Creating adata.raw with raw UMI counts...")
adata_raw = sc.AnnData(
    X=adata.layers["counts"],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
adata.raw = adata_raw

print("✓ adata.raw.X now points to raw UMI counts")
print(f"  Memory saved: ~{adata.layers['counts'].data.nbytes / 1e9:.1f} GB (shared)")

# ⭐ P0-3 FIX: Gene name handling strategy
print("\n" + "-"*70)
print("Gene Name Handling (P0-3 Fix)")
print("-"*70)

sample_gene = str(adata.var_names[0])
print(f"Sample gene: {sample_gene}")

# Check if we need symbol column
if 'symbol' not in adata.var.columns:
    if sample_gene.startswith('ENSG'):
        print("✓ Detected ENSEMBL IDs, will create symbol column")
        
        # Try to use existing annotations first
        symbol_candidates = ['gene_symbols', 'feature_name', 'gene_name']
        existing_symbol_col = None
        for col in symbol_candidates:
            if col in adata.var.columns:
                existing_symbol_col = col
                print(f"  Found existing symbol column: {col}")
                break
        
        if existing_symbol_col:
            adata.var['symbol'] = adata.var[existing_symbol_col].astype(str)
        else:
            print("  No existing symbol column, using mygene for conversion...")
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
            
            # Keep ENSEMBL if no symbol found
            adata.var['symbol'] = [ensembl_to_symbol.get(e, e) for e in adata.var_names]
            
            print(f"  Converted: {len(ensembl_to_symbol):,}/{len(ensembl_ids):,} genes")
    else:
        # Already symbols
        print("✓ Using gene symbols from var_names")
        adata.var['symbol'] = adata.var_names.astype(str)
else:
    print("✓ Symbol column already exists")

# Create symbol_base (remove -1/-2 suffixes from make_unique)
adata.var['symbol_base'] = adata.var['symbol'].str.replace(r'-\d+$', '', regex=True)
print("✓ Created symbol_base column for CellTypist matching")
print(f"  Example: {adata.var['symbol'].iloc[0]} → {adata.var['symbol_base'].iloc[0]}")

# Keep var_names unchanged (critical for model reuse)
print("⚠️  NOTE: var_names unchanged (ENSEMBL or original)")
print(f"  This ensures model reusability across runs")

print("\n✓ Data preparation complete")


# ==============================================================================
# Step 2: scVI Integration (with HVG Consistency)
# ==============================================================================

print("\n" + "="*70)
print("Step 2: scVI Integration (P0-1 Fix: HVG Consistency)")
print("="*70)

scvi_model_dir = output_dir / "models" / "scvi_model"
scvi_model_exists = scvi_model_dir.exists() and (scvi_model_dir / "model.pt").exists()
hvg_file = scvi_model_dir / "hvg_genes.txt"

hvg_indices = None
actual_hvg_count = adata.n_vars
hvg_selection_method = "all_genes"
hvg_genes_list = None

if USE_HVG_FOR_SCVI:
    print(f"\n⭐ Using HVG approach ({N_HVG_SCVI} genes)")
    
    # ⭐ P0-1 FIX: Load HVG from saved file if model exists
    if scvi_model_exists and hvg_file.exists() and LOAD_SCVI_IF_EXISTS:
        print("\n" + "-"*70)
        print("Loading HVG list from saved model (P0-1 Fix)")
        print("-"*70)
        
        try:
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            print(f"✓ Loaded {len(hvg_genes_saved):,} HVG genes from {hvg_file.name}")
            
            # Match with current adata
            hvg_genes_in_data = [g for g in hvg_genes_saved if g in set(adata.var_names)]
            
            if len(hvg_genes_in_data) < len(hvg_genes_saved) * 0.95:
                print(f"⚠️  Only {len(hvg_genes_in_data)}/{len(hvg_genes_saved)} HVG found in current data")
                print(f"   Will recompute HVG (data mismatch)")
                hvg_genes_list = None
            else:
                hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_in_data))
                # Preserve order from saved file
                hvg_genes_list = hvg_genes_in_data
                actual_hvg_count = len(hvg_genes_list)
                hvg_selection_method = "loaded_from_file"
                print(f"✓ HVG genes matched and ordered: {actual_hvg_count:,}")
        
        except Exception as e:
            print(f"⚠️  Failed to load HVG file: {e}")
            print("   Will recompute HVG")
            hvg_genes_list = None
    
    # Compute HVG if not loaded
    if hvg_genes_list is None:
        print("\n" + "-"*70)
        print("Computing HVG")
        print("-"*70)
        
        print(f"Computing HVG using {HVG_FLAVOR} method...")
        
        try:
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=N_HVG_SCVI,
                batch_key=BATCH_KEY,
                flavor=HVG_FLAVOR,
                subset=False
            )
            hvg_selection_method = "batch-aware"
            print("✓ Batch-aware HVG selection succeeded")
        except Exception as e:
            print(f"⚠️  Batch-aware HVG failed: {type(e).__name__}")
            print("   Falling back to non-batch-aware")
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=N_HVG_SCVI,
                flavor=HVG_FLAVOR,
                subset=False
            )
            hvg_selection_method = "non-batch-aware"
            print("✓ Non-batch-aware HVG selection succeeded")
        
        hvg_mask = adata.var["highly_variable"].values
        hvg_indices = np.flatnonzero(hvg_mask)
        hvg_genes_list = adata.var_names[hvg_indices].tolist()
        actual_hvg_count = len(hvg_genes_list)
        
        print(f"✓ Selected {actual_hvg_count:,} HVG genes")
    
    # Construct scVI AnnData
    print("\n" + "-"*70)
    print("Preparing scVI AnnData")
    print("-"*70)
    
    adata_scvi = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[[BATCH_KEY]].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scvi.var_names = pd.Index(hvg_genes_list)
    
    print(f"✓ scVI AnnData: {adata_scvi.n_obs:,} cells × {adata_scvi.n_vars:,} genes")
    
else:
    print("\n⚠️  Using ALL genes")
    adata_scvi = adata.copy()
    adata_scvi.X = adata.layers['counts'].copy()
    print(f"✓ scVI AnnData: {adata_scvi.n_obs:,} cells × {adata_scvi.n_vars:,} genes")

# Load or train scVI model
print("\n" + "-"*70)
print("scVI Model: Load or Train")
print("-"*70)

if scvi_model_exists and LOAD_SCVI_IF_EXISTS:
    print(f"\n✓ Found existing scVI model at: {scvi_model_dir}")
    print("  Loading model...")
    
    try:
        scvi_model = scvi.model.SCVI.load(scvi_model_dir, adata=adata_scvi)
        print("✓ scVI model loaded successfully")
        
        if hasattr(scvi_model, 'is_trained_') and scvi_model.is_trained_:
            print("✓ Model is trained")
        else:
            print("⚠️  Model not trained, will train now")
            scvi_model_exists = False
    except Exception as e:
        print(f"⚠️  Failed to load: {type(e).__name__}: {e}")
        print("  Will train new model")
        scvi_model_exists = False

if not scvi_model_exists or not LOAD_SCVI_IF_EXISTS:
    print("\nTraining new scVI model...")
    
    # Setup
    scvi.model.SCVI.setup_anndata(
        adata_scvi,
        layer=None,
        batch_key=BATCH_KEY
    )
    
    # Initialize
    scvi_model = scvi.model.SCVI(
        adata_scvi,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        dropout_rate=SCVI_DROPOUT_RATE,
        gene_likelihood="nb"
    )
    
    try:
        n_params = scvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
    print(f"✓ scVI model initialized ({n_params:,} parameters)")
    
    # Train
    print(f"\nTraining scVI (max {SCVI_MAX_EPOCHS} epochs)...")
    start_time = time.time()
    
    # ⭐ P2-1 FIX: Explicit lr in plan_kwargs
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'early_stopping': SCVI_EARLY_STOPPING,
        'early_stopping_patience': SCVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {'lr': SCVI_LEARNING_RATE},  # ⭐ Explicit learning rate
        'enable_model_summary': True,
        'enable_progress_bar': True,
    }
    
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scvi_model.train(**train_kwargs)
    
    training_time = time.time() - start_time
    print(f"\n✓ scVI training completed in {training_time/60:.1f} minutes")
    
    # Save model
    scvi_model.save(scvi_model_dir, overwrite=True)
    print(f"✓ scVI model saved to: {scvi_model_dir}")
    
    # ⭐ P0-1 FIX: Save HVG gene list
    if USE_HVG_FOR_SCVI and hvg_genes_list is not None:
        pd.Series(hvg_genes_list).to_csv(hvg_file, index=False, header=False)
        print(f"✓ HVG gene list saved to: {hvg_file.name}")
    
    # Plot training history
    print("\nPlotting training history...")
    train_history = scvi_model.history
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    train_elbo = None
    val_elbo = None
    if hasattr(train_history, '__contains__'):
        if 'elbo_train' in train_history:
            train_elbo = train_history['elbo_train']
        if 'elbo_validation' in train_history:
            val_elbo = train_history['elbo_validation']
    
    if train_elbo is not None and len(train_elbo) > 0:
        axes[0].plot(train_elbo, label='Train', linewidth=2)
    if val_elbo is not None and len(val_elbo) > 0:
        axes[0].plot(val_elbo, label='Validation', linewidth=2)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('ELBO Loss')
    axes[0].set_title('scVI Training - ELBO')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    train_recon = None
    val_recon = None
    if hasattr(train_history, '__contains__'):
        if 'reconstruction_loss_train' in train_history:
            train_recon = train_history['reconstruction_loss_train']
        if 'reconstruction_loss_validation' in train_history:
            val_recon = train_history['reconstruction_loss_validation']
    
    if train_recon is not None and len(train_recon) > 0:
        axes[1].plot(train_recon, label='Train', linewidth=2)
    if val_recon is not None and len(val_recon) > 0:
        axes[1].plot(val_recon, label='Validation', linewidth=2)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Reconstruction Loss')
    axes[1].set_title('scVI Training - Reconstruction')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"scvi_training_history.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.close()
    
    print("✓ Training history saved")

# Generate scVI latent representation
print("\n" + "-"*70)
print("Generating scVI latent representation")
print("-"*70)

print("Computing scVI latent space...")
latent_df = pd.DataFrame(
    scvi_model.get_latent_representation(),
    index=adata_scvi.obs_names
)

# ⭐ P1-3 FIX: Assert no NaNs after reindex
latent_aligned = latent_df.reindex(adata.obs_names)
if latent_aligned.isna().any().any():
    raise RuntimeError("scVI latent reindex produced NaNs - obs_names mismatch")

# ⭐ P2-3 FIX: Use float32 for storage
adata.obsm['X_scvi'] = latent_aligned.to_numpy().astype(np.float32)

print(f"✓ scVI latent space: {adata.obsm['X_scvi'].shape}")

# Clean up
del adata_scvi
gc.collect()
print("✓ Cleaned up temporary scVI data")

# ⭐ P1-2 FIX: Save neighbors with key
print("\nComputing neighborhood graph (neighbors_scvi)...")
sc.pp.neighbors(
    adata, 
    use_rep='X_scvi', 
    n_neighbors=UMAP_N_NEIGHBORS, 
    key_added='neighbors_scvi',  # ⭐ Separate key
    random_state=RANDOM_SEED
)
print("✓ Neighbor graph saved as 'neighbors_scvi'")

print("Computing UMAP from scVI latent...")
sc.tl.umap(
    adata, 
    min_dist=UMAP_MIN_DIST, 
    spread=UMAP_SPREAD, 
    neighbors_key='neighbors_scvi',  # ⭐ Use specific neighbors
    random_state=RANDOM_SEED
)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print("✓ scVI UMAP saved as 'X_umap_scvi'")

# Clustering
print("\nClustering on scVI latent...")
sc.tl.leiden(
    adata, 
    resolution=0.5, 
    key_added='leiden_scvi', 
    neighbors_key='neighbors_scvi',  # ⭐ Use specific neighbors
    random_state=RANDOM_SEED
)
print(f"✓ Identified {adata.obs['leiden_scvi'].nunique()} clusters")

# Visualize
print("\nVisualizing scVI integration...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scVI: Batch')
if TISSUE_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 1], show=False, title='scVI: Tissue')
sc.pl.umap(adata, color='leiden_scvi', ax=axes[0, 2], show=False, title='scVI: Leiden')

if EXISTING_CELLTYPE_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=EXISTING_CELLTYPE_KEY, ax=axes[1, 0], show=False, title='scVI: Existing Type')

if 'n_genes_by_counts' in adata.obs.columns:
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1, 1], show=False, title='scVI: # Genes', cmap='viridis')
if 'total_counts' in adata.obs.columns:
    sc.pl.umap(adata, color='total_counts', ax=axes[1, 2], show=False, title='scVI: Total Counts', cmap='viridis')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"scvi_umap_overview.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.close()

print("✓ scVI visualization complete")


# ==============================================================================
# Step 3: CellTypist Annotation (Memory Efficient)
# ==============================================================================

print("\n" + "="*70)
print("Step 3: CellTypist Annotation (P0-2 Fix: Memory Efficient)")
print("="*70)

print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL_PATH}")
celltypist_model = models.Model.load(CELLTYPIST_MODEL_PATH)
print(f"✓ Model loaded")

# ⭐ P0-2 FIX: Get model features and subset first
print("\n" + "-"*70)
print("Getting CellTypist model features (P0-2 Fix)")
print("-"*70)

model_features = None
for attr in ['features', 'genes', 'var_names']:
    if hasattr(celltypist_model, attr):
        mf = getattr(celltypist_model, attr)
        if mf is not None and len(mf) > 100:
            model_features = pd.Index(mf).astype(str)
            print(f"✓ Found {len(model_features)} model features via '{attr}'")
            break

if model_features is None:
    raise RuntimeError("Cannot retrieve CellTypist model features")

# Match model features with symbol_base
print("\nMatching genes using symbol_base...")
common_genes = adata.var['symbol_base'].isin(model_features)
n_common = common_genes.sum()

print(f"  Overlapping genes: {n_common:,}/{len(model_features)} model features")

if n_common < 500:
    raise RuntimeError(f"Too few overlapping genes: {n_common}")

# ⭐ P0-2 FIX: Only copy subset of counts
print("\nConstructing CellTypist AnnData (memory efficient)...")
print(f"  Only copying {n_common:,} genes (not all {adata.n_vars:,})")

# ⭐ HOTFIX: Convert pandas Series to numpy array for scipy sparse indexing
X_subset = adata.layers["counts"][:, common_genes.values]

adata_celltypist = sc.AnnData(
    X=X_subset.copy(),  # ⭐ Only subset
    obs=adata.obs.copy(),
    var=adata.var.loc[common_genes].copy()
)

# Use symbol_base as var_names for CellTypist
adata_celltypist.var_names = adata.var.loc[common_genes, 'symbol_base'].values

print("✓ CellTypist AnnData created")
print(f"  Shape: {adata_celltypist.n_obs:,} cells × {adata_celltypist.n_vars:,} genes")
print(f"  Memory saved: ~{(adata.n_vars - n_common) / adata.n_vars * 100:.1f}%")

# Normalize and log-transform
print("  Normalizing...")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)

print("  Log1p transforming...")
sc.pp.log1p(adata_celltypist)

print("✓ Data prepared for CellTypist")

# Run CellTypist
print("\nRunning CellTypist prediction...")
predictions = celltypist.annotate(
    adata_celltypist,
    model=celltypist_model,
    majority_voting=True
)

print("✓ CellTypist annotation complete")

# Extract results with index alignment
predicted_labels = predictions.predicted_labels

# ⭐ P1-3 FIX: Assert no NaNs
pred_df = predicted_labels.reindex(adata.obs_names)
if pred_df.isna().any().any():
    raise RuntimeError("CellTypist predictions reindex produced NaNs")

adata.obs['celltypist_predicted'] = pred_df.predicted_labels.astype(str).values
if 'majority_voting' in pred_df.columns:
    adata.obs['celltypist_majority_voting'] = pred_df.majority_voting.astype(str).values
else:
    adata.obs['celltypist_majority_voting'] = pred_df.predicted_labels.astype(str).values

# Confidence score
conf_candidates = ['conf_score', 'confidence', 'confidence_score', 'prob']
confidence_column = next((col for col in conf_candidates if col in predicted_labels.columns), None)

if confidence_column is not None:
    adata.obs['celltypist_conf_score'] = pred_df[confidence_column].reindex(adata.obs_names).values
else:
    adata.obs['celltypist_conf_score'] = 1.0

# Clean up
del adata_celltypist, X_subset
gc.collect()
print("✓ Cleaned up temporary CellTypist data")

# Summary
print(f"\nCellTypist annotation summary:")
print(f"  Predicted labels: {adata.obs['celltypist_predicted'].nunique()} unique types")
print(f"  After majority voting: {adata.obs['celltypist_majority_voting'].nunique()} unique types")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
print(f"  Median confidence: {adata.obs['celltypist_conf_score'].median():.3f}")

print("\nTop 10 cell types:")
print(adata.obs['celltypist_majority_voting'].value_counts().head(10))

# Visualize
print("\nVisualizing CellTypist annotations...")

fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.umap(adata, color='celltypist_predicted', ax=axes[0], show=False, 
           title='CellTypist: Predicted')
sc.pl.umap(adata, color='celltypist_majority_voting', ax=axes[1], show=False, 
           title='CellTypist: Majority Voting')
sc.pl.umap(adata, color='celltypist_conf_score', ax=axes[2], show=False, 
           title='CellTypist: Confidence', cmap='RdYlGn')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"celltypist_annotations.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.close()

print("✓ CellTypist visualization complete")


# ==============================================================================
# Step 4: scANVI Semi-supervised Refinement (True Semi-supervised)
# ==============================================================================

print("\n" + "="*70)
print("Step 4: scANVI Refinement (P1-1 Fix: True Semi-supervised)")
print("="*70)

# ⭐ P1-1 FIX: Create true semi-supervised labels
print("\n" + "-"*70)
print("Creating semi-supervised labels (P1-1 Fix)")
print("-"*70)

print(f"Strategy:")
print(f"  1. Set confidence < {LOW_CONFIDENCE_THRESHOLD} → Unknown")
print(f"  2. Set types with < {RARE_TYPE_THRESHOLD} cells → Unknown")

# Start with CellTypist labels
labels = adata.obs['celltypist_majority_voting'].astype(str).copy()
conf = adata.obs['celltypist_conf_score'].astype(float)

# Track changes
n_original = len(labels)
n_low_conf = (conf < LOW_CONFIDENCE_THRESHOLD).sum()

# Set low confidence to Unknown
labels[conf < LOW_CONFIDENCE_THRESHOLD] = "Unknown"

# Set rare types to Unknown
value_counts = labels.value_counts()
rare_types = value_counts[value_counts < RARE_TYPE_THRESHOLD].index
rare_types = rare_types[rare_types != "Unknown"]  # Don't count Unknown itself
n_rare = labels.isin(rare_types).sum()

labels[labels.isin(rare_types)] = "Unknown"

# Summary
n_unknown = (labels == "Unknown").sum()
n_labeled = n_original - n_unknown

print(f"\nSemi-supervised label statistics:")
print(f"  Total cells: {n_original:,}")
print(f"  Low confidence → Unknown: {n_low_conf:,} ({n_low_conf/n_original*100:.1f}%)")
print(f"  Rare types → Unknown: {n_rare:,} ({n_rare/n_original*100:.1f}%)")
print(f"  Final labeled: {n_labeled:,} ({n_labeled/n_original*100:.1f}%)")
print(f"  Final Unknown: {n_unknown:,} ({n_unknown/n_original*100:.1f}%)")
print(f"  Unique types (excl Unknown): {labels[labels != 'Unknown'].nunique()}")

# Save as categorical
adata.obs['scanvi_label'] = pd.Categorical(labels)
scanvi_label_key = 'scanvi_label'

print(f"\n✓ Semi-supervised labels created: '{scanvi_label_key}'")

# Prepare scANVI data
scanvi_model_dir = output_dir / "models" / "scanvi_model"
scanvi_model_exists = scanvi_model_dir.exists() and (scanvi_model_dir / "model.pt").exists()

print("\n" + "-"*70)
print("Preparing scANVI AnnData")
print("-"*70)

if USE_HVG_FOR_SCVI:
    print(f"  Using same {actual_hvg_count} HVG subset as scVI")
    
    if hvg_indices is None:
        if 'highly_variable' in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            hvg_indices = np.flatnonzero(hvg_mask)
        else:
            # Load from file
            if hvg_file.exists():
                hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
                hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_saved))
    
    adata_scanvi = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[[BATCH_KEY, scanvi_label_key]].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scanvi.var_names = adata.var_names[hvg_indices]
else:
    print("  Using full gene set")
    adata_scanvi = sc.AnnData(
        X=adata.layers["counts"].copy(),
        obs=adata.obs[[BATCH_KEY, scanvi_label_key]].copy(),
        var=adata.var.copy()
    )

print(f"✓ scANVI AnnData: {adata_scanvi.n_obs:,} cells × {adata_scanvi.n_vars:,} genes")

# Ensure 'Unknown' category exists (safety check)
labels_cat = adata_scanvi.obs[scanvi_label_key].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi.obs[scanvi_label_key] = labels_cat
    print("✓ Added 'Unknown' category (safety check)")

# Load or train scANVI
print("\n" + "-"*70)
print("scANVI Model: Load or Train")
print("-"*70)

if scanvi_model_exists and LOAD_SCANVI_IF_EXISTS:
    print(f"\n✓ Found existing scANVI model at: {scanvi_model_dir}")
    print("  Loading model...")
    
    try:
        scanvi_model = scvi.model.SCANVI.load(scanvi_model_dir, adata=adata_scanvi)
        print("✓ scANVI model loaded successfully")
        
        if hasattr(scanvi_model, 'is_trained_') and scanvi_model.is_trained_:
            print("✓ Model is trained")
        else:
            print("⚠️  Model not trained, will train now")
            scanvi_model_exists = False
    except Exception as e:
        print(f"⚠️  Failed to load: {type(e).__name__}: {e}")
        print("  Will train new model")
        scanvi_model_exists = False

if not scanvi_model_exists or not LOAD_SCANVI_IF_EXISTS:
    print("\nTraining new scANVI model...")
    
    # Setup
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi,
        layer=None,
        batch_key=BATCH_KEY,
        labels_key=scanvi_label_key,
        unlabeled_category="Unknown"
    )
    
    # Initialize from scVI
    print("\nInitializing scANVI from pre-trained scVI...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        adata=adata_scanvi,
        labels_key=scanvi_label_key
    )
    
    try:
        n_params = scanvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)
    print(f"✓ scANVI model initialized ({n_params:,} parameters)")
    
    # Train
    print(f"\nTraining scANVI (max {SCANVI_MAX_EPOCHS} epochs)...")
    start_time = time.time()
    
    # ⭐ P2-1 FIX: Explicit lr
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'early_stopping': True,
        'early_stopping_patience': SCANVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {'lr': SCANVI_LEARNING_RATE},  # ⭐ Explicit learning rate
        'enable_model_summary': True,
        'enable_progress_bar': True,
    }
    
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scanvi_model.train(**train_kwargs)
    
    training_time = time.time() - start_time
    print(f"\n✓ scANVI training completed in {training_time/60:.1f} minutes")
    
    # Save
    scanvi_model.save(scanvi_model_dir, overwrite=True)
    print(f"✓ scANVI model saved to: {scanvi_model_dir}")
    
    # Plot training history
    print("\nPlotting scANVI training history...")
    train_history = scanvi_model.history
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    train_elbo = None
    val_elbo = None
    if hasattr(train_history, '__contains__'):
        if 'elbo_train' in train_history:
            train_elbo = train_history['elbo_train']
        if 'elbo_validation' in train_history:
            val_elbo = train_history['elbo_validation']
    
    if train_elbo is not None and len(train_elbo) > 0:
        axes[0].plot(train_elbo, label='Train', linewidth=2)
    if val_elbo is not None and len(val_elbo) > 0:
        axes[0].plot(val_elbo, label='Validation', linewidth=2)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('ELBO Loss')
    axes[0].set_title('scANVI Training - ELBO')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    train_class = None
    val_class = None
    if hasattr(train_history, '__contains__'):
        if 'classification_loss_train' in train_history:
            train_class = train_history['classification_loss_train']
        if 'classification_loss_validation' in train_history:
            val_class = train_history['classification_loss_validation']
    
    if train_class is not None and len(train_class) > 0:
        axes[1].plot(train_class, label='Train', linewidth=2)
    if val_class is not None and len(val_class) > 0:
        axes[1].plot(val_class, label='Validation', linewidth=2)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Classification Loss')
    axes[1].set_title('scANVI Training - Classification')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"scanvi_training_history.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.close()
    
    print("✓ Training history saved")

# Generate scANVI results
print("\n" + "-"*70)
print("Generating scANVI predictions and embeddings")
print("-"*70)

# Latent representation
latent_scanvi_df = pd.DataFrame(
    scanvi_model.get_latent_representation(),
    index=adata_scanvi.obs_names
)

# ⭐ P1-3 + P2-3 FIX: Assert no NaNs + float32
latent_scanvi_aligned = latent_scanvi_df.reindex(adata.obs_names)
if latent_scanvi_aligned.isna().any().any():
    raise RuntimeError("scANVI latent reindex produced NaNs")

adata.obsm['X_scanvi'] = latent_scanvi_aligned.to_numpy().astype(np.float32)
print(f"✓ scANVI latent space: {adata.obsm['X_scanvi'].shape}")

# Predictions
predictions_scanvi_series = pd.Series(
    scanvi_model.predict(),
    index=adata_scanvi.obs_names
)

# ⭐ P1-3 FIX: Assert no NaNs
predictions_aligned = predictions_scanvi_series.reindex(adata.obs_names)
if predictions_aligned.isna().any():
    raise RuntimeError("scANVI predictions reindex produced NaNs")

adata.obs['scanvi_predictions'] = predictions_aligned.astype(str).values
print(f"✓ scANVI predictions: {adata.obs['scanvi_predictions'].nunique()} types")

# Probabilities
probs = np.asarray(scanvi_model.predict(soft=True))

if hasattr(scanvi_model, 'classifier_') and hasattr(scanvi_model.classifier_, 'classifier'):
    try:
        celltype_names = scanvi_model.classifier_.classifier.classes_
        adata.uns['scanvi_celltype_order'] = celltype_names.tolist()
    except:
        adata.uns['scanvi_celltype_order'] = [f'Type_{i}' for i in range(probs.shape[1])]
else:
    adata.uns['scanvi_celltype_order'] = [f'Type_{i}' for i in range(probs.shape[1])]

probs_df = pd.DataFrame(probs, index=adata_scanvi.obs_names)

# ⭐ P1-3 + P2-3 FIX: Assert no NaNs + float32
probs_aligned = probs_df.reindex(adata.obs_names)
if probs_aligned.isna().any().any():
    raise RuntimeError("scANVI probabilities reindex produced NaNs")

adata.obsm['scanvi_probabilities'] = probs_aligned.to_numpy().astype(np.float32)
print("✓ Posterior probabilities saved")

# Clean up
del adata_scanvi
gc.collect()
print("✓ Cleaned up temporary scANVI data")

# ⭐ P1-2 FIX: Separate neighbors graph
print("\nComputing neighborhood graph (neighbors_scanvi)...")
sc.pp.neighbors(
    adata, 
    use_rep='X_scanvi', 
    n_neighbors=UMAP_N_NEIGHBORS, 
    key_added='neighbors_scanvi',  # ⭐ Separate key
    random_state=RANDOM_SEED
)
print("✓ Neighbor graph saved as 'neighbors_scanvi'")

print("Computing UMAP from scANVI latent...")
sc.tl.umap(
    adata, 
    min_dist=UMAP_MIN_DIST, 
    spread=UMAP_SPREAD, 
    neighbors_key='neighbors_scanvi',  # ⭐ Use specific neighbors
    random_state=RANDOM_SEED
)
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
print("✓ scANVI UMAP saved as 'X_umap_scanvi'")

# Clustering
print("\nClustering on scANVI latent...")
sc.tl.leiden(
    adata, 
    resolution=0.5, 
    key_added='leiden_scanvi', 
    neighbors_key='neighbors_scanvi',  # ⭐ Use specific neighbors
    random_state=RANDOM_SEED
)
print(f"✓ Identified {adata.obs['leiden_scanvi'].nunique()} clusters")

# Summary
print(f"\nscANVI annotation summary:")
print(f"  Predicted types: {adata.obs['scanvi_predictions'].nunique()}")
print(f"  Types (excl Unknown): {adata.obs[adata.obs['scanvi_predictions'] != 'Unknown']['scanvi_predictions'].nunique()}")

print("\nTop 10 cell types:")
print(adata.obs['scanvi_predictions'].value_counts().head(10))

# Visualize
print("\nVisualizing scANVI results...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scANVI: Batch')
if TISSUE_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 1], show=False, title='scANVI: Tissue')
sc.pl.umap(adata, color='leiden_scanvi', ax=axes[0, 2], show=False, title='scANVI: Leiden')

sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1, 0], show=False, title='scANVI: Predictions')
sc.pl.umap(adata, color='scanvi_label', ax=axes[1, 1], show=False, title='scANVI: Input Labels')
sc.pl.umap(adata, color='celltypist_conf_score', ax=axes[1, 2], show=False, 
           title='CellTypist: Confidence', cmap='RdYlGn')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"scanvi_umap_overview.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.close()

print("✓ scANVI visualization complete")


# ==============================================================================
# Step 5: Final Data Export
# ==============================================================================

print("\n" + "="*70)
print("Step 5: Saving Final Results")
print("="*70)

# Set default UMAP to scANVI
if 'X_umap_scanvi' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()
    print("✓ Default UMAP set to scANVI")
elif 'X_umap_scvi' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()
    print("✓ Default UMAP set to scVI")

# Save final dataset
output_file = output_dir / "adata_allcells_scvi_celltypist_scanvi_final_v2.1.3.1.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"✓ Final dataset saved ({file_size:.2f} GB)")

# Generate summary
print("\n" + "="*70)
print("ANALYSIS SUMMARY v2.1.3.1")
print("="*70)

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("All Cells Analysis Summary v2.1.3.1 (Hotfix)")
summary_lines.append("="*70)
summary_lines.append(f"\nPipeline Version: 2.1.3.1 (scipy indexing hotfix)")
summary_lines.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_lines.append(f"Random Seed: {RANDOM_SEED}")

summary_lines.append("\n[ Input Data ]")
summary_lines.append(f"  File: {INPUT_H5AD}")
summary_lines.append(f"  Cells: {adata.n_obs:,}")
summary_lines.append(f"  Genes: {adata.n_vars:,}")
summary_lines.append(f"  Batches: {adata.obs[BATCH_KEY].nunique()}")

summary_lines.append("\n[ scVI Integration ]")
summary_lines.append(f"  Method: {'HVG subset' if USE_HVG_FOR_SCVI else 'All genes'}")
if USE_HVG_FOR_SCVI:
    summary_lines.append(f"  HVG genes: {actual_hvg_count:,}")
    summary_lines.append(f"  HVG method: {hvg_selection_method}")
summary_lines.append(f"  Model loaded: {scvi_model_exists and LOAD_SCVI_IF_EXISTS}")
summary_lines.append(f"  Neighbors key: neighbors_scvi")

summary_lines.append("\n[ CellTypist Annotation ]")
summary_lines.append(f"  Model: {Path(CELLTYPIST_MODEL_PATH).name}")
summary_lines.append(f"  Predicted types: {adata.obs['celltypist_majority_voting'].nunique()}")
summary_lines.append(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
summary_lines.append(f"  Memory optimization: Only {n_common:,}/{adata.n_vars:,} genes copied")

summary_lines.append("\n[ scANVI Refinement (Semi-supervised) ]")
summary_lines.append(f"  Input labeled: {n_labeled:,}/{n_original:,} ({n_labeled/n_original*100:.1f}%)")
summary_lines.append(f"  Input unlabeled: {n_unknown:,}/{n_original:,} ({n_unknown/n_original*100:.1f}%)")
summary_lines.append(f"  Low confidence threshold: {LOW_CONFIDENCE_THRESHOLD}")
summary_lines.append(f"  Rare type threshold: {RARE_TYPE_THRESHOLD} cells")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions'].nunique()}")
summary_lines.append(f"  Model loaded: {scanvi_model_exists and LOAD_SCANVI_IF_EXISTS}")
summary_lines.append(f"  Neighbors key: neighbors_scanvi")

summary_lines.append("\n  Top 10 cell types:")
for celltype, count in adata.obs['scanvi_predictions'].value_counts().head(10).items():
    pct = count / adata.n_obs * 100
    summary_lines.append(f"    {celltype}: {count:,} ({pct:.1f}%)")

summary_lines.append("\n[ P0 Fixes (Critical) ]")
summary_lines.append("  ✓ P0-1: HVG consistency - Save/load HVG gene list")
summary_lines.append("  ✓ P0-2: CellTypist memory - Only model features copied")
summary_lines.append("  ✓ P0-3: Gene names - Use symbol_base, keep var_names")

summary_lines.append("\n[ P1 Fixes (Important) ]")
summary_lines.append("  ✓ P1-1: True semi-supervised - Low-conf + rare → Unknown")
summary_lines.append("  ✓ P1-2: Graph separation - neighbors_scvi/scanvi keys")
summary_lines.append("  ✓ P1-3: NaN assertions - Check reindex results")

summary_lines.append("\n[ P2 Optimizations ]")
summary_lines.append("  ✓ P2-1: Explicit lr in plan_kwargs")
summary_lines.append("  ✓ P2-3: Float32 for latent/probs storage")
summary_lines.append("  ✓ Matplotlib Agg backend for HPC")

summary_lines.append("\n[ Output Files ]")
summary_lines.append(f"  Final data: {output_file.name} ({file_size:.2f} GB)")
summary_lines.append(f"  scVI model: {scvi_model_dir}")
summary_lines.append(f"  scANVI model: {scanvi_model_dir}")
if USE_HVG_FOR_SCVI:
    summary_lines.append(f"  HVG genes: {hvg_file.name}")

summary_lines.append("\n" + "="*70)
summary_lines.append("Analysis completed successfully!")
summary_lines.append("="*70)

summary_text = '\n'.join(summary_lines)
print(summary_text)

# Save summary
summary_file = output_dir / "analysis_summary_v2.1.3.1.txt"
with open(summary_file, 'w', encoding='utf-8') as f:
    f.write(summary_text)

print(f"\n✓ Summary saved to: {summary_file}")


# ==============================================================================
# Final Notes
# ==============================================================================

print("\n" + "="*70)
print("🎉 PIPELINE COMPLETE (v2.1.3.1 - Production + Hotfix)")
print("="*70)

print(f"\n📁 Output Directory: {output_dir}")
print(f"\n📊 Main Results:")
print(f"  • Final data: {output_file.name}")
print(f"  • Summary: {summary_file.name}")
print(f"  • Figures: figures/")
print(f"  • Models: models/ (reusable!)")

print(f"\n🔬 Key Improvements in v2.1.3:")
print(f"  • HVG Consistency: {hvg_file.name if USE_HVG_FOR_SCVI else 'N/A'}")
print(f"  • Memory Efficient: CellTypist uses {n_common:,} genes (not {adata.n_vars:,})")
print(f"  • True Semi-supervised: {n_unknown:,} unlabeled cells")
print(f"  • Separate Graphs: neighbors_scvi + neighbors_scanvi")
print(f"  • NaN Safe: All reindex operations validated")

print(f"\n💡 Next Steps:")
print(f"  1. Validate cell types with marker genes")
print(f"  2. Check Unknown cells - may need manual review")
print(f"  3. Use neighbors_scvi/scanvi for downstream analysis")
print(f"  4. Models are reusable - set LOAD_*_IF_EXISTS=True")

print("\n" + "="*70)
print()