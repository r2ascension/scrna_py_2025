#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Epithelial Cell: scVI → CellTypist → scANVI Pipeline
=====================================================

Complete 3-stage semi-supervised annotation pipeline:
1. scVI: Unsupervised batch correction
2. CellTypist: Automated annotation on scVI latent space
3. Combine CellTypist + Manual annotations
4. scANVI: Semi-supervised fine-tuning

Memory optimized for 270k+ cells with zero-copy strategy.

Input: adata_with_manual_annotations.h5ad
Output: adata_epithelial_FINAL_annotated.h5ad

Author: r2end
Date: 2024-12-09
Version: 3.0 - Complete 3-stage pipeline with scvi-tools 1.1+ compatibility
"""

import os
import warnings
warnings.filterwarnings('ignore')

# Core libraries
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import time

# scvi-tools (CRITICAL: Version 1.1+ required)
import scvi
import torch

# CellTypist
import celltypist
from celltypist import models

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output paths
INPUT_FILE = "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/annotation_results/adata_with_manual_annotations.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1210/bbknn_celltype_analysis/Epithelial/scanvi_results"
MODEL_DIR = os.path.join(OUTPUT_DIR, "models")

# Create directories
os.makedirs(OUTPUT_DIR, exist_ok=True)  
os.makedirs(MODEL_DIR, exist_ok=True)

# Model hyperparameters
N_LATENT = 100              # Latent space dimensions
N_LAYERS = 3               # Hidden layers
N_HVG = 4000               # Highly variable genes

# Training parameters
SCVI_MAX_EPOCHS = 800      # scVI epochs
SCANVI_MAX_EPOCHS = 600    # scANVI epochs
BATCH_SIZE = 2048          # Training batch size
LEARNING_RATE = 1e-3       # Learning rate
EARLY_STOPPING = True      # Early stopping
EARLY_STOPPING_PATIENCE = 50

# CellTypist model
CELLTYPIST_MODEL = '/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl'  # Can change based on cell type
CELLTYPIST_MAJORITY_VOTING = True

# Column names
BATCH_KEY = 'batch'
MANUAL_LABEL_KEY = 'cell_type_manual'
CELLTYPIST_LABEL_KEY = 'cell_type_celltypist'
SCANVI_LABEL_KEY = 'cell_type_scanvi_pred'
CLUSTER_KEY = 'leiden_bbknn'

# Random seed
RANDOM_SEED = 42

print("="*80)
print("EPITHELIAL 3-STAGE scVI-CellTypist-scANVI PIPELINE")
print("="*80)
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Workflow: scVI → CellTypist → Manual merge → scANVI")
print(f"\nModel Parameters:")
print(f"  n_latent: {N_LATENT}, n_layers: {N_LAYERS}, n_hvg: {N_HVG}")
print(f"  scVI epochs: {SCVI_MAX_EPOCHS}, scANVI epochs: {SCANVI_MAX_EPOCHS}")

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print(f"\n{'='*80}")
print("ENVIRONMENT CHECK")
print("="*80)

# Set random seeds
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

print(f"\nPackage versions:")
print(f"  - scanpy: {sc.__version__}")
print(f"  - scvi-tools: {scvi.__version__}")
print(f"  - torch: {torch.__version__}")
print(f"  - celltypist: {celltypist.__version__}")

# GPU check
gpu_available = torch.cuda.is_available()
if gpu_available:
    gpu_name = torch.cuda.get_device_name(0)
    gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"\n✓ GPU available: {gpu_name}")
    print(f"  Memory: {gpu_memory:.1f} GB")
    accelerator = 'gpu'
    devices = 'auto'
    scvi.settings.dl_num_workers = 0  # Stable for HPC
else:
    print(f"\n⚠️  GPU not available, using CPU")
    print(f"  Training will be slower (~5-8 hours vs ~40-60 min)")
    accelerator = 'cpu'
    devices = 'auto'

# ============================================================================
# LOAD DATA
# ============================================================================

print(f"\n{'='*80}")
print("DATA LOADING")
print("="*80)

print(f"\nLoading: {INPUT_FILE}")
adata = sc.read_h5ad(INPUT_FILE)

print(f"✓ Data loaded")
print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Validate required columns
required_cols = [BATCH_KEY, MANUAL_LABEL_KEY, CLUSTER_KEY]
missing_cols = [col for col in required_cols if col not in adata.obs.columns]
if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

print(f"✓ All required columns present")

# Check annotation status
n_unknown = (adata.obs[MANUAL_LABEL_KEY] == 'Unknown').sum()
n_skip = (adata.obs[MANUAL_LABEL_KEY] == 'Skip').sum()
n_labeled = adata.n_obs - n_unknown - n_skip

print(f"\nManual annotation status:")
print(f"  Total cells: {adata.n_obs:,}")
print(f"  Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  Unknown: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(f"  Skip: {n_skip:,} ({n_skip/adata.n_obs*100:.1f}%)")

# ============================================================================
# PREPARE RAW COUNTS - MEMORY OPTIMIZED (Zero-Copy Strategy)
# ============================================================================

print(f"\n{'='*80}")
print("PREPARE DATA - MEMORY OPTIMIZED")
print("="*80)

print(f"\n⚠️  CRITICAL: scVI requires raw count data (not log-normalized)")
print(f"⚠️  Using zero-copy strategy to minimize memory usage")

# Select counts layer - NO COPYING!
counts_layer = None

if 'counts' in adata.layers:
    counts_layer = 'counts'
    print(f"✓ Using adata.layers['counts'] (no copy)")
elif 'raw_counts' in adata.layers:
    counts_layer = 'raw_counts'
    print(f"✓ Using adata.layers['raw_counts'] (no copy)")
elif adata.raw is not None:
    print(f"✓ Found adata.raw - creating reference (no full copy)")
    adata.layers['counts'] = adata.raw.X
    counts_layer = 'counts'
else:
    print(f"\n❌ No raw counts found! Checking current adata.X...")
    test_values = adata.X[:1000, :100]
    if hasattr(adata.X, 'toarray'):
        test_values = test_values.toarray()
    max_val = test_values.max()
    has_decimals = np.any(test_values != test_values.astype(int))
    
    if max_val < 20 and has_decimals:
        raise ValueError(
            f"❌ Data appears log-normalized (max={max_val:.2f})\n"
            f"   scVI requires raw counts! Please reload before normalization."
        )
    else:
        print(f"  Data looks like counts (max={max_val:.0f}), using adata.X")
        counts_layer = None

# Point X to counts (reference, not copy)
if counts_layer is not None:
    adata.X = adata.layers[counts_layer]

# Basic filtering
print(f"\nBasic filtering...")
print(f"  Before: {adata.shape[1]:,} genes")
sc.pp.filter_genes(adata, min_cells=3)
print(f"  After gene filter: {adata.shape[1]:,} genes")

n_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
if n_before > adata.n_obs:
    print(f"  Removed {n_before - adata.n_obs:,} low-quality cells")

# Remove Skip cells
print(f"\nRemoving 'Skip' cells...")
n_before = adata.n_obs
adata = adata[adata.obs[MANUAL_LABEL_KEY] != 'Skip'].copy()
n_removed = n_before - adata.n_obs
print(f"  Removed {n_removed:,} 'Skip' cells")
print(f"  Remaining: {adata.n_obs:,} cells")

# ⭐ CRITICAL: Save a copy for CellTypist BEFORE HVG selection
# CellTypist works better with full genes, not HVG subset
print(f"\nCreating copy for CellTypist (full genes)...")
adata_ct = adata.copy()
print(f"✓ CellTypist copy created: {adata_ct.n_obs:,} cells × {adata_ct.n_vars:,} genes")

# ===== HVG Selection with Robust Fallback =====
print(f"\nSelecting {N_HVG} highly variable genes...")

# Robust HVG selection strategy
try:
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        flavor='seurat_v3',
        subset=False  # ⭐ 不删除基因
    )
    hvg_method = "batch-aware"
    print(f"✓ Batch-aware HVG selection succeeded")
except Exception as e:
    print(f"⚠️  Batch-aware HVG failed: {e}")
    print(f"   Reason: Likely small batches or other batch-related issue")
    print(f"   Falling back to non-batch-aware method...")
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=N_HVG,
        flavor='seurat_v3',
        subset=False
    )
    hvg_method = "non-batch-aware"
    print(f"✓ Fallback HVG selection succeeded")

# Record method used
adata.uns['hvg_method'] = hvg_method
print(f"   Method used: {hvg_method}")
print(f"   HVG selected: {adata.var['highly_variable'].sum():,}")

# ⭐ CRITICAL: Preserve full genes to .raw BEFORE subsetting
print(f"\n⭐ Preserving full gene data to .raw (shared memory)...")
adata.raw = sc.AnnData(
    X=adata.layers["counts"],  # Shared memory, no .copy()
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ Full genes preserved: {adata.raw.n_vars:,} genes")
print(f"   Memory cost: 0 GB (shared memory)")

# Subset to HVG for model training
print(f"\n⭐ Subsetting to HVG for scVI/scANVI training...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"✓ Training data ready: {adata.n_obs:,} cells × {adata.n_vars:,} HVG")

# Verify structure
print(f"\n✓ Data structure verified:")
print(f"   adata.X: {adata.n_vars:,} HVG (for training)")
print(f"   adata.raw.X: {adata.raw.n_vars:,} genes (for analysis)")
print(f"   Extra memory: 0 GB (shared)")

# Save checkpoint
checkpoint_file = os.path.join(OUTPUT_DIR, "adata_prepared.h5ad")
print(f"\nSaving prepared data: {checkpoint_file}")
adata.write_h5ad(checkpoint_file)

print(f"\n✓ Data preparation complete")
print(f"  Final shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ============================================================================
# STAGE 1: scVI PRE-TRAINING (Unsupervised Batch Correction)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 1: scVI PRE-TRAINING")
print("="*80)

print(f"\nModel configuration:")
print(f"  - Latent dimensions: {N_LATENT}")
print(f"  - Hidden layers: {N_LAYERS}")
print(f"  - Max epochs: {SCVI_MAX_EPOCHS}")
print(f"  - Batch size: {BATCH_SIZE}")
print(f"  - Device: {accelerator}")

print(f"\nEstimated time:")
print(f"  GPU: ~15-30 minutes")
print(f"  CPU: ~2-4 hours")

# Setup scVI
print(f"\nSetting up scVI model...")
scvi.model.SCVI.setup_anndata(
    adata,
    batch_key=BATCH_KEY,
    layer=None  # Use adata.X (raw counts)
)

vae = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    gene_likelihood='nb'
)

print(f"✓ scVI model created")
print(f"  Parameters: {sum(p.numel() for p in vae.module.parameters()):,}")

# Train scVI
print(f"\n{'='*80}")
print(f"Training scVI...")
print(f"{'='*80}\n")

start_time = time.time()

# CRITICAL: scvi-tools 1.1+ - lr goes in plan_kwargs, not top-level
train_kwargs = {
    'max_epochs': SCVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'early_stopping': EARLY_STOPPING,
    'accelerator': accelerator,
    'devices': devices,
    'plan_kwargs': {'lr': LEARNING_RATE},
}

if EARLY_STOPPING:
    train_kwargs['early_stopping_patience'] = EARLY_STOPPING_PATIENCE

vae.train(**train_kwargs)

elapsed = time.time() - start_time
print(f"\n{'='*80}")
print(f"✓ scVI training complete!")
print(f"  Time: {int(elapsed//60)} min {int(elapsed%60)} sec")
print(f"{'='*80}")

# Save scVI model
scvi_model_dir = os.path.join(MODEL_DIR, "scvi_model")
print(f"\nSaving scVI model: {scvi_model_dir}")
vae.save(scvi_model_dir, overwrite=True)

# Generate latent representation
print(f"Generating scVI latent representation...")
adata.obsm['X_scvi'] = vae.get_latent_representation()
print(f"✓ X_scvi saved, shape: {adata.obsm['X_scvi'].shape}")

# Plot training history (with error handling for scvi-tools 1.1+)
print(f"\nPlotting training history...")

history = vae.history

train_elbo = None
val_elbo = None

if hasattr(history, '__contains__') and 'elbo_train' in history:
    train_elbo = history['elbo_train']
if hasattr(history, '__contains__') and 'elbo_validation' in history:
    val_elbo = history['elbo_validation']

fig, ax = plt.subplots(figsize=(10, 6))

if train_elbo is not None and len(train_elbo) > 0:
    ax.plot(train_elbo, label='Training', linewidth=2)
if val_elbo is not None and len(val_elbo) > 0:
    ax.plot(val_elbo, label='Validation', linewidth=2)

ax.set_xlabel('Epoch', fontsize=12)
ax.set_ylabel('ELBO Loss', fontsize=12)
ax.set_title('scVI Training History', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '01_scvi_training_history.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Stage 1 complete!")

# ============================================================================
# STAGE 2: CellTypist ANNOTATION (Automated)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 2: CellTypist AUTOMATED ANNOTATION")
print("="*80)

print(f"\nCellTypist configuration:")
print(f"  - Model: {CELLTYPIST_MODEL}")
print(f"  - Majority voting: {CELLTYPIST_MAJORITY_VOTING}")
print(f"  - Using: X_scvi latent space")

# Load CellTypist model
print(f"\nLoading CellTypist model...")
try:
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"✓ Model loaded: {CELLTYPIST_MODEL}")
except:
    print(f"⚠️  Model not found locally, downloading...")
    models.download_models(model=CELLTYPIST_MODEL, force_update=True)
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"✓ Model downloaded and loaded")

# Run CellTypist prediction
print(f"\nRunning CellTypist prediction...")
print(f"  This may take 5-15 minutes for {adata_ct.n_obs:,} cells...")
print(f"  Using full gene matrix: {adata_ct.n_vars:,} genes")

start_time = time.time()

# CellTypist expects normalized log1p data
# Check and normalize adata_ct if needed
if 'log1p' in adata_ct.layers:
    print(f"  Using adata_ct.layers['log1p'] for CellTypist")
    adata_ct_for_pred = adata_ct.copy()
    adata_ct_for_pred.X = adata_ct.layers['log1p']
elif adata_ct.X.max() > 20:  # Looks like raw counts
    print(f"  Normalizing data for CellTypist...")
    adata_ct_for_pred = adata_ct.copy()
    sc.pp.normalize_total(adata_ct_for_pred, target_sum=1e4)
    sc.pp.log1p(adata_ct_for_pred)
else:
    print(f"  Using current adata_ct.X for CellTypist")
    adata_ct_for_pred = adata_ct.copy()

predictions = celltypist.annotate(
    adata_ct_for_pred,
    model=model,
    majority_voting=CELLTYPIST_MAJORITY_VOTING
)

elapsed = time.time() - start_time
print(f"✓ CellTypist prediction complete!")
print(f"  Time: {int(elapsed//60)} min {int(elapsed%60)} sec")

# Extract predictions and map to main adata (HVG subset)
adata.obs[CELLTYPIST_LABEL_KEY] = predictions.predicted_labels['predicted_labels'].values
if CELLTYPIST_MAJORITY_VOTING:
    adata.obs[f'{CELLTYPIST_LABEL_KEY}_majority'] = predictions.predicted_labels['majority_voting'].values
adata.obs['celltypist_conf_score'] = predictions.predicted_labels['conf_score'].values

# Clean up to save memory
del adata_ct, adata_ct_for_pred
import gc
gc.collect()

# Show CellTypist results
print(f"\nCellTypist results:")
celltypist_counts = adata.obs[CELLTYPIST_LABEL_KEY].value_counts()
print(f"  Unique cell types: {len(celltypist_counts)}")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
print(f"\n  Top 10 cell types:")
for cell_type, count in celltypist_counts.head(10).items():
    pct = count / adata.n_obs * 100
    print(f"    - {cell_type}: {count:,} ({pct:.1f}%)")

print(f"\n✓ Stage 2 complete!")

# ============================================================================
# STAGE 3: MERGE ANNOTATIONS (CellTypist + Manual)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 3: MERGE CELLTYPIST + MANUAL ANNOTATIONS")
print("="*80)

print(f"\nMerging strategy:")
print(f"  1. Use manual annotation if available (not 'Unknown')")
print(f"  2. Otherwise use CellTypist prediction")
print(f"  3. Final label → 'combined_cell_type'")

# Create combined labels
def merge_annotations(row):
    manual = row[MANUAL_LABEL_KEY]
    celltypist = row[CELLTYPIST_LABEL_KEY]
    
    if manual != 'Unknown':
        return manual  # Use manual if available
    else:
        return celltypist  # Use CellTypist for unknowns

adata.obs['combined_cell_type'] = adata.obs.apply(merge_annotations, axis=1)

# Track annotation source
def get_annotation_source(row):
    if row[MANUAL_LABEL_KEY] != 'Unknown':
        return 'Manual'
    else:
        return 'CellTypist'

adata.obs['annotation_source'] = adata.obs.apply(get_annotation_source, axis=1)

# Statistics
n_manual = (adata.obs['annotation_source'] == 'Manual').sum()
n_celltypist = (adata.obs['annotation_source'] == 'CellTypist').sum()

print(f"\nCombined annotation statistics:")
print(f"  From manual: {n_manual:,} ({n_manual/adata.n_obs*100:.1f}%)")
print(f"  From CellTypist: {n_celltypist:,} ({n_celltypist/adata.n_obs*100:.1f}%)")

# Show combined distribution
combined_counts = adata.obs['combined_cell_type'].value_counts()
print(f"\n  Combined cell types: {len(combined_counts)}")
print(f"  Top 10:")
for cell_type, count in combined_counts.head(10).items():
    pct = count / adata.n_obs * 100
    source_counts = adata.obs[adata.obs['combined_cell_type'] == cell_type]['annotation_source'].value_counts()
    print(f"    - {cell_type}: {count:,} ({pct:.1f}%)")
    print(f"      [Manual: {source_counts.get('Manual', 0)}, CellTypist: {source_counts.get('CellTypist', 0)}]")

# Save checkpoint
checkpoint_file = os.path.join(OUTPUT_DIR, "adata_combined_annotations.h5ad")
print(f"\nSaving combined annotations: {checkpoint_file}")
adata.write_h5ad(checkpoint_file)

print(f"\n✓ Stage 3 complete!")

# ============================================================================
# STAGE 4: scANVI FINE-TUNING (Semi-supervised)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 4: scANVI FINE-TUNING")
print("="*80)

print(f"\nModel configuration:")
print(f"  - Base: scVI (pre-trained)")
print(f"  - Max epochs: {SCANVI_MAX_EPOCHS}")
print(f"  - Labels: combined_cell_type")
print(f"  - Unlabeled: 'Unknown' (placeholder, all cells labeled)")

print(f"\nEstimated time:")
print(f"  GPU: ~10-20 minutes")
print(f"  CPU: ~1-2 hours")

# Prepare labels for scANVI
adata.obs['cell_type_for_scanvi'] = adata.obs['combined_cell_type'].copy()

# ⭐ CRITICAL: Ensure 'Unknown' is in categories (scvi-tools 1.1+ requirement)
# Even if no cells use it, scANVI needs it as a valid category
print(f"\nPreparing label categories for scANVI...")
labels = adata.obs['cell_type_for_scanvi'].astype('category')
if 'Unknown' not in labels.cat.categories:
    print(f"  Adding 'Unknown' as placeholder category (no cells use it)")
    labels = labels.cat.add_categories(['Unknown'])
    adata.obs['cell_type_for_scanvi'] = labels
else:
    adata.obs['cell_type_for_scanvi'] = labels

print(f"\nLabel distribution for scANVI:")
scanvi_label_counts = adata.obs['cell_type_for_scanvi'].value_counts()
print(f"  Unique types: {len(scanvi_label_counts)}")
n_unknown = (adata.obs['cell_type_for_scanvi'] == 'Unknown').sum()
print(f"  'Unknown' cells: {n_unknown} (should be 0)")
print(f"  All cells have labels: {n_unknown == 0}")

# ⭐ CRITICAL: Initialize scANVI from scVI (in-memory, NO re-setup)
# scvi-tools 1.1+: DO NOT call SCVI.setup_anndata again before from_scvi_model
# The vae object already has its AnnDataManager configured
print(f"\nInitializing scANVI from pre-trained scVI...")
print(f"  ⚠️  Using in-memory vae object (no model reload)")

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category='Unknown',
    labels_key='cell_type_for_scanvi'
)

print(f"✓ scANVI model created")
print(f"  Parameters: {sum(p.numel() for p in lvae.module.parameters()):,}")

# Train scANVI
print(f"\n{'='*80}")
print(f"Training scANVI...")
print(f"{'='*80}\n")

start_time = time.time()

# CRITICAL: scvi-tools 1.1+ - lr goes in plan_kwargs
train_kwargs = {
    'max_epochs': SCANVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'accelerator': accelerator,
    'devices': devices,
    'plan_kwargs': {'lr': LEARNING_RATE},
}

lvae.train(**train_kwargs)

elapsed = time.time() - start_time
print(f"\n{'='*80}")
print(f"✓ scANVI training complete!")
print(f"  Time: {int(elapsed//60)} min {int(elapsed%60)} sec")
print(f"{'='*80}")

# Save scANVI model
scanvi_model_dir = os.path.join(MODEL_DIR, "scanvi_model")
print(f"\nSaving scANVI model: {scanvi_model_dir}")
lvae.save(scanvi_model_dir, overwrite=True)

# Generate predictions
print(f"\nGenerating final predictions...")
adata.obs[SCANVI_LABEL_KEY] = lvae.predict()
predictions_probs = lvae.predict(soft=True)
adata.obs['scanvi_confidence'] = predictions_probs.max(axis=1)
adata.obsm['X_scanvi'] = lvae.get_latent_representation()

print(f"✓ Predictions saved")
print(f"  Mean confidence: {adata.obs['scanvi_confidence'].mean():.3f}")

# Plot training history (with error handling for scvi-tools 1.1+)
print(f"\nPlotting training history...")

history = lvae.history

train_elbo = None
val_elbo = None

if hasattr(history, '__contains__') and 'elbo_train' in history:
    train_elbo = history['elbo_train']
if hasattr(history, '__contains__') and 'elbo_validation' in history:
    val_elbo = history['elbo_validation']

fig, ax = plt.subplots(figsize=(10, 6))

if train_elbo is not None and len(train_elbo) > 0:
    ax.plot(train_elbo, label='Training', linewidth=2)
if val_elbo is not None and len(val_elbo) > 0:
    ax.plot(val_elbo, label='Validation', linewidth=2)

ax.set_xlabel('Epoch', fontsize=12)
ax.set_ylabel('ELBO Loss', fontsize=12)
ax.set_title('scANVI Training History', fontsize=14, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '02_scanvi_training_history.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"\n✓ Stage 4 complete!")

# ============================================================================
# VALIDATION
# ============================================================================

print(f"\n{'='*80}")
print("VALIDATION")
print("="*80)

# Compare annotations
print(f"\n1. ANNOTATION AGREEMENT")
print(f"{'='*60}")

# Compare: combined vs scanvi
agreements = (adata.obs['combined_cell_type'] == adata.obs[SCANVI_LABEL_KEY]).sum()
agreement_rate = agreements / adata.n_obs * 100

print(f"Combined vs scANVI agreement: {agreement_rate:.2f}%")
print(f"  Agreements: {agreements:,} / {adata.n_obs:,}")

if agreement_rate > 90:
    print(f"\n✓ Excellent: >90% agreement")
elif agreement_rate > 80:
    print(f"\n✓ Good: >80% agreement")
elif agreement_rate > 70:
    print(f"\n⚠️  Moderate: 70-80% agreement")
else:
    print(f"\n❌ Low: <70% agreement - review disagreements")

# Show top disagreements
if agreement_rate < 100:
    disagreements = adata.obs[adata.obs['combined_cell_type'] != adata.obs[SCANVI_LABEL_KEY]]
    if len(disagreements) > 0:
        print(f"\nTop 5 disagreements (combined → scANVI):")
        from collections import Counter
        disagree_pairs = list(zip(disagreements['combined_cell_type'], disagreements[SCANVI_LABEL_KEY]))
        for (combined, scanvi), count in Counter(disagree_pairs).most_common(5):
            print(f"  {combined} → {scanvi}: {count} cells")

# Confidence distribution
print(f"\n2. CONFIDENCE DISTRIBUTION")
print(f"{'='*60}")

mean_conf = adata.obs['scanvi_confidence'].mean()
median_conf = adata.obs['scanvi_confidence'].median()

adata.obs['confidence_category'] = pd.cut(
    adata.obs['scanvi_confidence'],
    bins=[0, 0.5, 0.8, 1.0],
    labels=['Low (<0.5)', 'Medium (0.5-0.8)', 'High (>0.8)']
)

conf_dist = adata.obs['confidence_category'].value_counts()

print(f"Mean confidence: {mean_conf:.3f}")
print(f"Median confidence: {median_conf:.3f}")
print(f"\nConfidence distribution:")
for cat, count in conf_dist.items():
    pct = count / adata.n_obs * 100
    print(f"  {cat}: {count:,} ({pct:.1f}%)")

# Final cell type distribution
print(f"\n3. FINAL CELL TYPE DISTRIBUTION")
print(f"{'='*60}")

final_counts = adata.obs[SCANVI_LABEL_KEY].value_counts()
print(f"\n{'Cell Type':<30} {'Count':<12} {'Percent':<10} {'Mean Conf':<12}")
print("-"*70)
for cell_type, count in final_counts.head(15).items():
    pct = count / adata.n_obs * 100
    mean_conf = adata.obs[adata.obs[SCANVI_LABEL_KEY] == cell_type]['scanvi_confidence'].mean()
    print(f"{cell_type:<30} {count:<12,} {pct:<10.2f}% {mean_conf:<12.3f}")

# ============================================================================
# UMAP VISUALIZATION
# ============================================================================

print(f"\n{'='*80}")
print("UMAP VISUALIZATION ON scANVI LATENT SPACE")
print("="*80)

print(f"\nComputing UMAP...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)
sc.tl.umap(adata)

print(f"✓ UMAP computed")

# Generate comprehensive UMAP
print(f"\nGenerating UMAP visualizations...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

# Panel 1: scANVI predictions
sc.pl.umap(adata, color=SCANVI_LABEL_KEY, ax=axes[0], show=False,
          title='scANVI Predicted Cell Types', frameon=False, s=5)

# Panel 2: Confidence
sc.pl.umap(adata, color='scanvi_confidence', ax=axes[1], show=False,
          title='Prediction Confidence', frameon=False, cmap='viridis', vmin=0, vmax=1, s=5)

# Panel 3: Batch
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[2], show=False,
          title='Batch (Check Mixing)', frameon=False, s=3)

# Panel 4: Annotation source (ensure categorical for palette compatibility)
adata.obs['annotation_source'] = adata.obs['annotation_source'].astype('category')
sc.pl.umap(adata, color='annotation_source', ax=axes[3], show=False,
          title='Annotation Source', frameon=False,
          palette={'Manual': 'blue', 'CellTypist': 'red'}, s=5)

# Panel 5: Combined annotations
sc.pl.umap(adata, color='combined_cell_type', ax=axes[4], show=False,
          title='Combined Annotations', frameon=False, s=5)

# Panel 6: Confidence categories
sc.pl.umap(adata, color='confidence_category', ax=axes[5], show=False,
          title='Confidence Categories', frameon=False,
          palette={'High (>0.8)': 'green', 'Medium (0.5-0.8)': 'orange', 'Low (<0.5)': 'red'}, s=5)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '03_umap_comprehensive.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ UMAP saved: 03_umap_comprehensive.png")

# ============================================================================
# SAVE FINAL RESULTS
# ============================================================================

print(f"\n{'='*80}")
print("SAVE FINAL RESULTS")
print("="*80)

final_file = os.path.join(OUTPUT_DIR, "adata_epithelial_FINAL_annotated.h5ad")

print(f"\nSaving final dataset...")
print(f"File: {final_file}")

# Add metadata
adata.uns['pipeline_params'] = {
    'n_latent': N_LATENT,
    'n_layers': N_LAYERS,
    'scvi_epochs': SCVI_MAX_EPOCHS,
    'scanvi_epochs': SCANVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'n_hvg': N_HVG,
    'celltypist_model': CELLTYPIST_MODEL,
    'random_seed': RANDOM_SEED
}

adata.write_h5ad(final_file)

print(f"\n✓ Final dataset saved!")
print(f"  Size: {os.path.getsize(final_file) / 1e9:.2f} GB")

# Generate README
readme = f"""# Epithelial 3-Stage scVI-CellTypist-scANVI Results

## Overview
Complete semi-supervised annotation of {adata.n_obs:,} epithelial cells.

## Workflow
1. **scVI**: Unsupervised batch correction
2. **CellTypist**: Automated annotation on scVI space
3. **Merge**: Combined CellTypist + Manual annotations
4. **scANVI**: Semi-supervised fine-tuning

## Main Output
**adata_epithelial_FINAL_annotated.h5ad**

### Key Columns
- `adata.obs['{SCANVI_LABEL_KEY}']` - ⭐ Final cell type (USE THIS)
- `adata.obs['combined_cell_type']` - Pre-scANVI combined labels
- `adata.obs['{CELLTYPIST_LABEL_KEY}']` - CellTypist predictions
- `adata.obs['{MANUAL_LABEL_KEY}']` - Manual annotations
- `adata.obs['annotation_source']` - Manual vs CellTypist
- `adata.obs['scanvi_confidence']` - Confidence (0-1)

### Key Representations
- `adata.obsm['X_scanvi']` - ⭐ Batch-corrected latent space
- `adata.obsm['X_scvi']` - scVI latent space
- `adata.obsm['X_umap']` - UMAP coordinates

## Statistics
- Total cells: {adata.n_obs:,}
- Genes (HVG): {adata.n_vars:,}
- Cell types: {len(final_counts)}
- Mean confidence: {mean_conf:.3f}
- Agreement (combined vs scANVI): {agreement_rate:.1f}%
- Manual annotations: {n_manual:,} ({n_manual/adata.n_obs*100:.1f}%)
- CellTypist filled: {n_celltypist:,} ({n_celltypist/adata.n_obs*100:.1f}%)

## Models
- scVI: `models/scvi_model/` (reusable)
- scANVI: `models/scanvi_model/` (reusable)

## Figures
1. `01_scvi_training_history.png` - scVI training
2. `02_scanvi_training_history.png` - scANVI training
3. `03_umap_comprehensive.png` - 6-panel UMAP overview

## Next Steps
Use `adata.obs['{SCANVI_LABEL_KEY}']` for:
- Differential expression
- Trajectory inference
- Cell-cell communication
- Publication figures

---
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Pipeline version: 3.0 (scVI → CellTypist → scANVI)
"""

with open(os.path.join(OUTPUT_DIR, 'README.md'), 'w') as f:
    f.write(readme)

print(f"\n✓ README.md created")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"🎉 PIPELINE COMPLETE!")
print(f"{'='*80}")

print(f"\n📁 All results saved to:")
print(f"   {OUTPUT_DIR}")

print(f"\n⭐ Main output:")
print(f"   {final_file}")

print(f"\n📊 Summary:")
print(f"   - Cells: {adata.n_obs:,}")
print(f"   - Cell types: {len(final_counts)}")
print(f"   - Mean confidence: {mean_conf:.3f}")
print(f"   - Agreement: {agreement_rate:.1f}%")

print(f"\n💡 Use adata.obs['{SCANVI_LABEL_KEY}'] for final annotations!")

print(f"\n✅ Ready for downstream analysis!")
print(f"{'='*80}")
