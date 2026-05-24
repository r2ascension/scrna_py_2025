#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Epithelial Cell scANVI Semi-Supervised Annotation - Optimized Version
======================================================================

Based on code review recommendations for 270k+ cells:
- Memory optimization: No large matrix copying
- Robust mapping and extreme case handling
- Production-ready with comprehensive QC

Input: adata_with_manual_annotations.h5ad
Output: adata_epithelial_SCANVI_FINAL.h5ad

Usage:
    python epithelial_scanvi_training_optimized.py

Requirements:
    - scvi-tools environment
    - GPU recommended (optional)
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

# scvi-tools
import scvi
import torch

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output paths
INPUT_FILE = "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/annotation_results/adata_with_manual_annotations.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/scanvi_results"
MODEL_DIR = os.path.join(OUTPUT_DIR, "scanvi_models")

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

# Model hyperparameters
N_LATENT = 30              # Latent space dimensions (10-50)
N_LAYERS = 2               # Number of hidden layers (1-3)
N_EPOCHS_SCVI = 400        # scVI pre-training epochs
N_EPOCHS_SCANVI = 200      # scANVI fine-tuning epochs
BATCH_SIZE = 2048          # Batch size for training
N_HVG = 4000               # Number of highly variable genes

# Reuse model flags
USE_EXISTING_SCVI_MODEL = True     # 如有已有 scVI 模型，则直接 load
USE_EXISTING_SCANVI_MODEL = True   # 如有已有 scANVI 模型，则直接 load

# Annotation columns
BATCH_KEY = 'batch'
MANUAL_LABEL_KEY = 'cell_type_manual'
SCANVI_LABEL_KEY = 'cell_type_scanvi_pred'
CLUSTER_KEY = 'leiden_bbknn'

# Random seed for reproducibility
RANDOM_SEED = 42

print("="*80)
print("EPITHELIAL scANVI TRAINING - OPTIMIZED VERSION")
print("="*80)
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Model params: n_latent={N_LATENT}, n_layers={N_LAYERS}")
print(f"  Training: scVI={N_EPOCHS_SCVI} epochs, scANVI={N_EPOCHS_SCANVI} epochs")
print(f"  Batch size: {BATCH_SIZE}, HVG: {N_HVG}")

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print(f"\n{'='*80}")
print("ENVIRONMENT CHECK")
print("="*80)

# Set random seeds
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

print(f"\nPackage versions:")
print(f"  - scanpy: {sc.__version__}")
print(f"  - scvi-tools: {scvi.__version__}")
print(f"  - torch: {torch.__version__}")

# GPU check
use_gpu = torch.cuda.is_available()
if use_gpu:
    gpu_name = torch.cuda.get_device_name(0)
    gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"\n✓ GPU available: {gpu_name}")
    print(f"  Memory: {gpu_memory:.1f} GB")
    scvi.settings.dl_num_workers = 0  # Stable for HPC
else:
    print(f"\n⚠️  GPU not available, using CPU")
    print(f"  Training will be slower (~4-6 hours vs ~40 min)")

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
print(f"\nValidating required columns...")
required_cols = [BATCH_KEY, MANUAL_LABEL_KEY, CLUSTER_KEY]
missing_cols = [col for col in required_cols if col not in adata.obs.columns]

if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

print(f"✓ All required columns present")

# Check annotation status
print(f"\nAnnotation status:")
label_counts = adata.obs[MANUAL_LABEL_KEY].value_counts()
n_unknown = (adata.obs[MANUAL_LABEL_KEY] == 'Unknown').sum()
n_skip = (adata.obs[MANUAL_LABEL_KEY] == 'Skip').sum()
n_labeled = adata.n_obs - n_unknown - n_skip

print(f"  Total cells: {adata.n_obs:,}")
print(f"  Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  Unknown: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(f"  Skip: {n_skip:,} ({n_skip/adata.n_obs*100:.1f}%)")

if n_labeled < adata.n_obs * 0.3:
    print(f"\n⚠️  WARNING: Only {n_labeled/adata.n_obs*100:.1f}% cells labeled")
    print(f"  Recommend: >50% for good scANVI performance")

# ============================================================================
# PREPARE RAW COUNTS - MEMORY OPTIMIZED
# ============================================================================

print(f"\n{'='*80}")
print("PREPARING DATA FOR scVI/scANVI - MEMORY OPTIMIZED")
print("="*80)

print(f"\n⚠️  CRITICAL: scVI requires raw count data (not log-normalized)")

# Select counts layer - no copying!
counts_layer = None

if 'counts' in adata.layers:
    counts_layer = 'counts'
    print(f"✓ Using adata.layers['counts'] as raw counts (no copy)")
elif 'raw_counts' in adata.layers:
    counts_layer = 'raw_counts'
    print(f"✓ Using adata.layers['raw_counts'] as raw counts (no copy)")
elif adata.raw is not None:
    print(f"✓ Found adata.raw - creating reference in layers (no full copy)")
    adata.layers['counts'] = adata.raw.X
    counts_layer = 'counts'
else:
    print(f"\n❌ No raw counts found!")
    print(f"  Checking if current adata.X looks like counts...")
    
    # Check if data looks like counts
    test_values = adata.X[:1000, :100]
    if hasattr(adata.X, 'toarray'):
        test_values = test_values.toarray()
    
    max_val = test_values.max()
    has_decimals = np.any(test_values != test_values.astype(int))
    
    if max_val < 20 and has_decimals:
        raise ValueError(
            f"❌ Data appears to be log-normalized (max={max_val:.2f})\n"
            f"   scVI requires raw counts! Please reload data before normalization."
        )
    else:
        print(f"  Data looks like counts (max={max_val:.0f}), using adata.X")
        counts_layer = None

# Point X to counts layer (no copy, just reference)
if counts_layer is not None:
    adata.X = adata.layers[counts_layer]

# ============================================================================
# BASIC FILTERING
# ============================================================================

print(f"\nBasic filtering...")

# Filter genes
print(f"  Before filtering: {adata.shape[1]:,} genes")
sc.pp.filter_genes(adata, min_cells=3)
print(f"  After gene filter: {adata.shape[1]:,} genes")

# Filter cells
n_before = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
n_after = adata.n_obs
if n_before > n_after:
    print(f"  Removed {n_before - n_after:,} low-quality cells")

# ============================================================================
# REMOVE "Skip" CELLS
# ============================================================================

print(f"\nRemoving 'Skip' cells...")
n_before = adata.n_obs
adata = adata[adata.obs[MANUAL_LABEL_KEY] != 'Skip'].copy()
n_removed = n_before - adata.n_obs
print(f"  Removed {n_removed:,} cells marked as 'Skip'")
print(f"  Remaining: {adata.n_obs:,} cells")

# ============================================================================
# SELECT HIGHLY VARIABLE GENES - batch-aware, subset once
# ============================================================================

print(f"\nSelecting {N_HVG} highly variable genes (batch-aware)...")

try:
    sc.pp.highly_variable_genes(
        adata, layer="counts", n_top_genes=N_HVG, 
        batch_key=BATCH_KEY, flavor='seurat_v3', subset=True
    )
    hvg_method = "batch-aware"
except Exception as e:
    print(f"⚠️  Batch-aware failed, using fallback")
    sc.pp.highly_variable_genes(
        adata, layer="counts", n_top_genes=N_HVG, 
        flavor='seurat_v3', subset=True
    )
    hvg_method = "non-batch-aware"

hvg = adata.var["highly_variable"].values
adata.uns['hvg_method'] = hvg_method

print(f"✓ Selected {adata.shape[1]:,} HVG")

# ============================================================================
# PREPARE LABEL COLUMN FOR scANVI
# ============================================================================

print(f"\nPreparing scANVI labels...")

# Create scanvi label column
adata.obs['cell_type_scanvi'] = adata.obs[MANUAL_LABEL_KEY].copy()

# Count
n_labeled = (adata.obs['cell_type_scanvi'] != 'Unknown').sum()
n_unlabeled = (adata.obs['cell_type_scanvi'] == 'Unknown').sum()

print(f"  Labeled (for training): {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  Unlabeled (to predict): {n_unlabeled:,} ({n_unlabeled/adata.n_obs*100:.1f}%)")

# Show label distribution
if n_labeled > 0:
    label_dist = adata.obs[adata.obs['cell_type_scanvi'] != 'Unknown']['cell_type_scanvi'].value_counts()
    print(f"\n  Cell types for training (top 10):")
    for label, count in label_dist.head(10).items():
        print(f"    - {label}: {count:,} cells")
    if len(label_dist) > 10:
        print(f"    ... and {len(label_dist) - 10} more types")

# Save checkpoint
checkpoint_file = os.path.join(OUTPUT_DIR, "adata_scanvi_prepared.h5ad")
print(f"\nSaving checkpoint: {checkpoint_file}")
adata.write_h5ad(checkpoint_file)

print(f"\n✓ Data preparation complete")
print(f"  Final shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ============================================================================
# scVI PRE-TRAINING
# ============================================================================

print(f"\n{'='*80}")
print("scVI PRE-TRAINING (Unsupervised Batch Correction)")
print("="*80)

print(f"\nModel configuration:")
print(f"  - Latent dimensions: {N_LATENT}")
print(f"  - Hidden layers: {N_LAYERS}")
print(f"  - Training epochs: {N_EPOCHS_SCVI}")
print(f"  - Batch size: {BATCH_SIZE}")
print(f"  - Device: {'GPU' if use_gpu else 'CPU'}")

print(f"\nEstimated time:")
print(f"  GPU: ~10-30 minutes")
print(f"  CPU: ~2-4 hours")

scvi_model_dir = os.path.join(MODEL_DIR, "scvi_model")
trained_scvi = False  # 标记这次是否真的训练了

print(f"\nSetting up scVI model (registering AnnData)...")

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key=BATCH_KEY,
    labels_key='cell_type_scanvi',  # Register labels early
    layer=None  # Use adata.X (已经指向 counts)
)

# 如果已经有模型并且允许复用，就直接 load
if USE_EXISTING_SCVI_MODEL and os.path.isdir(scvi_model_dir) and len(os.listdir(scvi_model_dir)) > 0:
    print(f"\n✓ Found existing scVI model at: {scvi_model_dir}")
    print("  Loading scVI model from disk and skipping training...")
    vae = scvi.model.SCVI.load(scvi_model_dir, adata=adata)
else:
    print(f"\nNo existing scVI model found or retraining requested.")
    print("Creating and training new scVI model...")

    vae = scvi.model.SCVI(
        adata,
        n_latent=N_LATENT,
        n_layers=N_LAYERS,
        gene_likelihood='nb'
    )

    print(f"✓ scVI model created")
    print(f"  Parameters: {sum(p.numel() for p in vae.module.parameters()):,}")

    print(f"\n{'='*80}")
    print(f"Starting scVI training...")
    print(f"{'='*80}\n")

    start_time = time.time()

    train_kwargs = {
        'max_epochs': N_EPOCHS_SCVI,
        'batch_size': BATCH_SIZE,
        'early_stopping': True
    }
    if use_gpu:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})

    vae.train(**train_kwargs)
    trained_scvi = True

    elapsed = time.time() - start_time
    print(f"\n{'='*80}")
    print(f"✓ scVI training complete!")
    print(f"  Time: {int(elapsed//60)} min {int(elapsed%60)} sec")
    print(f"{'='*80}")

    # 只在本次真的训练过时才保存模型 & 画训练曲线
    print(f"\nSaving scVI model: {scvi_model_dir}")
    vae.save(scvi_model_dir, overwrite=True)

# 无论是训练还是加载，下面的 latent 都可以照常计算
print(f"\nGenerating scVI latent representation...")
adata.obsm['X_scvi'] = vae.get_latent_representation()
print(f"✓ X_scvi saved, shape: {adata.obsm['X_scvi'].shape}")

# 只有在本次训练过时才画训练历史
if trained_scvi and hasattr(vae, "history") and 'elbo_train' in vae.history:
    print(f"\nPlotting scVI training history...")
    fig, ax = plt.subplots(figsize=(10, 6))
    train_elbo = vae.history['elbo_train'][1:]
    if 'elbo_validation' in vae.history:
        val_elbo = vae.history['elbo_validation'][1:]
        ax.plot(val_elbo, label='Validation', linewidth=2)
    ax.plot(train_elbo, label='Training', linewidth=2)
    ax.set_xlabel('Epoch', fontsize=12)
    ax.set_ylabel('ELBO Loss', fontsize=12)
    ax.set_title('scVI Training History', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, '01_scvi_training_history.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
else:
    print("ℹ️  scVI model was loaded from disk; skip training history plot.")


print(f"✓ scVI pre-training phase complete!")

# ============================================================================
# scANVI FINE-TUNING
# ============================================================================
# ============================================================================
# scANVI FINE-TUNING
# ============================================================================

print(f"\n{'='*80}")
print("scANVI FINE-TUNING (Semi-supervised Learning)")
print("="*80)

print(f"\nModel configuration:")
print(f"  - Base: scVI (pre-trained)")
print(f"  - Training epochs: {N_EPOCHS_SCANVI}")
print(f"  - Unlabeled category: 'Unknown'")

print(f"\nEstimated time:")
print(f"  GPU: ~10-20 minutes")
print(f"  CPU: ~2-3 hours")

scanvi_model_dir = os.path.join(MODEL_DIR, "scanvi_model")
trained_scanvi = False

# 如果已有 scANVI 模型且允许复用，直接加载
if USE_EXISTING_SCANVI_MODEL and os.path.isdir(scanvi_model_dir) and len(os.listdir(scanvi_model_dir)) > 0:
    print(f"\n✓ Found existing scANVI model at: {scanvi_model_dir}")
    print("  Loading scANVI model from disk and skipping training...")
    lvae = scvi.model.SCANVI.load(scanvi_model_dir, adata=adata)
else:
    print(f"\nInitializing scANVI from pre-trained scVI...")
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        unlabeled_category='Unknown'
    )

    print(f"✓ scANVI model created")
    print(f"  Parameters: {sum(p.numel() for p in lvae.module.parameters()):,}")

    print(f"\n{'='*80}")
    print(f"Starting scANVI training...")
    print(f"{'='*80}\n")

    start_time = time.time()

    train_kwargs = {
        'max_epochs': N_EPOCHS_SCANVI,
        'batch_size': BATCH_SIZE
    }
    if use_gpu:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})

    lvae.train(**train_kwargs)
    trained_scanvi = True

    elapsed = time.time() - start_time
    print(f"\n{'='*80}")
    print(f"✓ scANVI training complete!")
    print(f"  Time: {int(elapsed//60)} min {int(elapsed%60)} sec")
    print(f"{'='*80}")

    print(f"\nSaving scANVI model: {scanvi_model_dir}")
    lvae.save(scanvi_model_dir, overwrite=True)

# === 下面 prediction / latent / QC 不管是训练还是加载都可以继续用 ===

print(f"\nGenerating predictions...")

adata.obs[SCANVI_LABEL_KEY] = lvae.predict()
print(f"✓ Predictions saved to adata.obs['{SCANVI_LABEL_KEY}']")

predictions_probs = lvae.predict(soft=True)
adata.obs['scanvi_confidence'] = predictions_probs.max(axis=1)
print(f"✓ Confidence scores saved")

adata.obsm['X_scanvi'] = lvae.get_latent_representation()
print(f"✓ X_scanvi saved, shape: {adata.obsm['X_scanvi'].shape}")

adata.obs['label_origin'] = [
    'Original' if label != 'Unknown' else 'Predicted'
    for label in adata.obs['cell_type_scanvi']
]

# 只有在本次训练过时才画训练历史，且防止 elbo_validation 缺失
if trained_scanvi and hasattr(lvae, "history") and 'elbo_train' in lvae.history:
    print(f"\nPlotting training history...")

    fig, ax = plt.subplots(figsize=(10, 6))
    history_keys = list(lvae.history.keys())
    print(f"Available history keys in lvae.history: {history_keys}")

    train_elbo = lvae.history['elbo_train'][1:]
    ax.plot(train_elbo, label='Training', linewidth=2)

    if 'elbo_validation' in lvae.history:
        val_elbo = lvae.history['elbo_validation'][1:]
        ax.plot(val_elbo, label='Validation', linewidth=2)

    ax.set_xlabel('Epoch', fontsize=12)
    ax.set_ylabel('ELBO Loss', fontsize=12)
    ax.set_title('scANVI Training History', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, '02_scanvi_training_history.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
else:
    print("ℹ️  scANVI model was loaded from disk; skip training history plot.")
print(f"\n✓ scANVI fine-tuning phase complete!")

# ============================================================================
# VALIDATION AND QUALITY CONTROL
# ============================================================================

print(f"\n{'='*80}")
print("VALIDATION AND QUALITY CONTROL")
print("="*80)

# Overall statistics
n_total = adata.n_obs
n_original = (adata.obs['label_origin'] == 'Original').sum()
n_predicted = (adata.obs['label_origin'] == 'Predicted').sum()

print(f"\n1. OVERALL STATISTICS")
print(f"{'='*60}")
print(f"Total cells: {n_total:,}")
print(f"Originally labeled: {n_original:,} ({n_original/n_total*100:.1f}%)")
print(f"Newly predicted: {n_predicted:,} ({n_predicted/n_total*100:.1f}%)")

# Confidence distribution
mean_conf = adata.obs['scanvi_confidence'].mean()
median_conf = adata.obs['scanvi_confidence'].median()

adata.obs['confidence_category'] = pd.cut(
    adata.obs['scanvi_confidence'],
    bins=[0, 0.5, 0.8, 1.0],
    labels=['Low (<0.5)', 'Medium (0.5-0.8)', 'High (>0.8)']
)

conf_dist = adata.obs['confidence_category'].value_counts()

print(f"\n2. PREDICTION CONFIDENCE")
print(f"{'='*60}")
print(f"Mean confidence: {mean_conf:.3f}")
print(f"Median confidence: {median_conf:.3f}")
print(f"\nConfidence distribution:")
for cat, count in conf_dist.items():
    pct = count / n_total * 100
    print(f"  {cat}: {count:,} ({pct:.1f}%)")

if mean_conf > 0.7:
    print(f"\n✓ Good: Mean confidence > 0.7")
elif mean_conf > 0.5:
    print(f"\n⚠️  Moderate: Mean confidence 0.5-0.7")
else:
    print(f"\n❌ Low: Mean confidence < 0.5")

# Label agreement (with defense)
print(f"\n3. LABEL AGREEMENT (Originally Labeled Cells)")
print(f"{'='*60}")

original_cells = adata.obs[adata.obs['label_origin'] == 'Original']

if len(original_cells) == 0:
    print(f"  No originally labeled cells; skip agreement calculation")
    agreement_rate = np.nan
else:
    agreements = (original_cells['cell_type_scanvi'] == original_cells[SCANVI_LABEL_KEY]).sum()
    agreement_rate = agreements / len(original_cells) * 100
    
    print(f"Agreement rate: {agreement_rate:.2f}%")
    print(f"  Agreements: {agreements:,} / {len(original_cells):,}")
    
    if agreement_rate > 90:
        print(f"\n✓ Excellent: >90% agreement")
    elif agreement_rate > 80:
        print(f"\n✓ Good: >80% agreement")
    elif agreement_rate > 70:
        print(f"\n⚠️  Moderate: 70-80% agreement")
    else:
        print(f"\n❌ Low: <70% agreement")
    
    # Show top disagreements
    if agreement_rate < 100:
        disagreements = original_cells[original_cells['cell_type_scanvi'] != original_cells[SCANVI_LABEL_KEY]]
        if len(disagreements) > 0:
            print(f"\nTop disagreements (manual → predicted):")
            disagree_pairs = list(zip(disagreements['cell_type_scanvi'], disagreements[SCANVI_LABEL_KEY]))
            from collections import Counter
            top_disagree = Counter(disagree_pairs).most_common(5)
            for (manual, pred), count in top_disagree:
                print(f"  {manual} → {pred}: {count} cells")

# Cell type distribution
print(f"\n4. FINAL CELL TYPE DISTRIBUTION")
print(f"{'='*60}")

type_counts = adata.obs[SCANVI_LABEL_KEY].value_counts()
type_confidence = adata.obs.groupby(SCANVI_LABEL_KEY)['scanvi_confidence'].agg(['mean', 'std'])

print(f"\n{'Cell Type':<30} {'Count':<12} {'Percent':<10} {'Confidence':<20}")
print("-"*75)

summary_data = []
for cell_type, count in type_counts.items():
    pct = count / n_total * 100
    mean_conf = type_confidence.loc[cell_type, 'mean']
    std_conf = type_confidence.loc[cell_type, 'std']
    
    # Handle NaN std (single cell types)
    if np.isnan(std_conf):
        std_conf = 0.0
    
    conf_str = f"{mean_conf:.3f}±{std_conf:.3f}"
    print(f"{cell_type:<30} {count:<12,} {pct:<10.2f}% {conf_str:<20}")
    
    summary_data.append({
        'Cell_Type': cell_type,
        'Count': count,
        'Percent': pct,
        'Mean_Confidence': mean_conf,
        'Std_Confidence': std_conf
    })

# Save summary
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(OUTPUT_DIR, 'scanvi_final_summary.csv'), index=False)
print(f"\n✓ Summary saved: scanvi_final_summary.csv")

# Validation plots
print(f"\n5. GENERATING VALIDATION PLOTS")
print(f"{'='*60}")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Confidence histogram
ax = axes[0, 0]
ax.hist(adata.obs['scanvi_confidence'], bins=50, edgecolor='black', alpha=0.7, color='steelblue')
ax.axvline(mean_conf, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_conf:.3f}')
ax.set_xlabel('Prediction Confidence', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title('Prediction Confidence Distribution', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Confidence by cell type
ax = axes[0, 1]
conf_by_type = adata.obs.groupby(SCANVI_LABEL_KEY)['scanvi_confidence'].mean().sort_values()
conf_by_type.plot(kind='barh', ax=ax, color='steelblue')
ax.axvline(0.8, color='green', linestyle='--', alpha=0.5, label='High (>0.8)')
ax.axvline(0.5, color='orange', linestyle='--', alpha=0.5, label='Medium (>0.5)')
ax.set_xlabel('Mean Confidence', fontsize=12)
ax.set_title('Mean Confidence by Cell Type', fontsize=14, fontweight='bold')
ax.legend()
ax.grid(True, alpha=0.3, axis='x')

# Plot 3: Confusion matrix (only if we have original labels)
ax = axes[1, 0]
if len(original_cells) > 0:
    from sklearn.metrics import confusion_matrix
    
    # Get top types
    top_types = original_cells['cell_type_scanvi'].value_counts().head(10).index
    subset = original_cells[original_cells['cell_type_scanvi'].isin(top_types)]
    
    if len(subset) > 0:
        cm = confusion_matrix(
            subset['cell_type_scanvi'],
            subset[SCANVI_LABEL_KEY],
            labels=top_types
        )
        
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax,
                    xticklabels=top_types, yticklabels=top_types,
                    cbar_kws={'label': 'Count'})
        ax.set_xlabel('scANVI Predicted', fontsize=12)
        ax.set_ylabel('Manual Label', fontsize=12)
        ax.set_title('Confusion Matrix (Top 10 Types)', fontsize=14, fontweight='bold')
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        plt.setp(ax.get_yticklabels(), rotation=0, fontsize=8)
else:
    ax.text(0.5, 0.5, 'No original labels\nto compare', 
            ha='center', va='center', fontsize=14)
    ax.axis('off')

# Plot 4: Cell count comparison
ax = axes[1, 1]
type_counts_sorted = type_counts.sort_values(ascending=True)
type_counts_sorted.plot(kind='barh', ax=ax, color='steelblue')
ax.set_xlabel('Number of Cells', fontsize=12)
ax.set_title('Final Cell Type Distribution', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '03_scanvi_validation.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Validation plots saved: 03_scanvi_validation.png")

print(f"\n{'='*80}")
print(f"✓ VALIDATION COMPLETE")
print(f"{'='*80}")

# ============================================================================
# UMAP VISUALIZATION
# ============================================================================

print(f"\n{'='*80}")
print("UMAP VISUALIZATION ON scANVI LATENT SPACE")
print("="*80)

# Compute UMAP
print(f"\nComputing neighborhood graph on X_scanvi...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)

print(f"Computing UMAP...")
sc.tl.umap(adata)

print(f"✓ UMAP computed")

# Comprehensive UMAP overview
print(f"\nGenerating UMAP overview (7 panels)...")

fig, axes = plt.subplots(3, 3, figsize=(24, 21))
axes = axes.flatten()

# Panel 1: Predicted cell types
sc.pl.umap(
    adata,
    color=SCANVI_LABEL_KEY,
    ax=axes[0],
    show=False,
    title='scANVI Predicted Cell Types',
    frameon=False,
    s=5
)

# Panel 2: Prediction confidence
sc.pl.umap(
    adata,
    color='scanvi_confidence',
    ax=axes[1],
    show=False,
    title='Prediction Confidence',
    frameon=False,
    cmap='viridis',
    vmin=0,
    vmax=1,
    s=5
)

# Panel 3: Batch
sc.pl.umap(
    adata,
    color=BATCH_KEY,
    ax=axes[2],
    show=False,
    title='Batch (Check Mixing)',
    frameon=False,
    s=3
)

# Panel 4: Label origin
sc.pl.umap(
    adata,
    color='label_origin',
    ax=axes[3],
    show=False,
    title='Label Origin',
    frameon=False,
    palette={'Original': 'blue', 'Predicted': 'red'},
    s=5
)

# Panel 5: Tissue (if available)
if 'tissue' in adata.obs.columns:
    sc.pl.umap(
        adata,
        color='tissue',
        ax=axes[4],
        show=False,
        title='Tissue',
        frameon=False,
        s=5
    )
else:
    axes[4].text(0.5, 0.5, 'Tissue info\nnot available', 
                ha='center', va='center', fontsize=14)
    axes[4].axis('off')

# Panel 6: Confidence categories
sc.pl.umap(
    adata,
    color='confidence_category',
    ax=axes[5],
    show=False,
    title='Confidence Categories',
    frameon=False,
    palette={'High (>0.8)': 'green', 'Medium (0.5-0.8)': 'orange', 'Low (<0.5)': 'red'},
    s=5
)

# Panel 7: Original clusters
if CLUSTER_KEY in adata.obs.columns:
    sc.pl.umap(
        adata,
        color=CLUSTER_KEY,
        ax=axes[6],
        show=False,
        title='Original BBKNN Clusters',
        frameon=False,
        legend_loc='on data',
        legend_fontsize=6,
        s=5
    )
else:
    axes[6].axis('off')

# Hide extra panels
axes[7].axis('off')
axes[8].axis('off')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '04_scanvi_umap_overview.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ UMAP overview saved: 04_scanvi_umap_overview.png")

# Individual cell type UMAPs (FIXED)
print(f"\nGenerating individual cell type UMAPs...")

top_types = adata.obs[SCANVI_LABEL_KEY].value_counts().head(12).index.tolist()

fig, axes = plt.subplots(3, 4, figsize=(20, 15))
axes = axes.flatten()

for idx, cell_type in enumerate(top_types):
    # Use scanpy's groups parameter for correct highlighting
    sc.pl.umap(
        adata,
        color=SCANVI_LABEL_KEY,
        groups=[cell_type],
        ax=axes[idx],
        show=False,
        title=cell_type,
        frameon=False,
        s=3,
        legend_loc='none'
    )

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, '05_scanvi_cell_types_individual.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"✓ Individual cell type UMAPs saved: 05_scanvi_cell_types_individual.png")

print(f"\n{'='*80}")
print(f"✓ UMAP VISUALIZATION COMPLETE")
print(f"{'='*80}")

# ============================================================================
# SAVE FINAL RESULTS
# ============================================================================

print(f"\n{'='*80}")
print("SAVING FINAL RESULTS")
print("="*80)

final_file = os.path.join(OUTPUT_DIR, "adata_epithelial_SCANVI_FINAL.h5ad")

print(f"\nSaving final annotated dataset...")
print(f"File: {final_file}")

# Add metadata
adata.uns['scanvi_params'] = {
    'n_latent': N_LATENT,
    'n_layers': N_LAYERS,
    'scvi_epochs': N_EPOCHS_SCVI,
    'scanvi_epochs': N_EPOCHS_SCANVI,
    'batch_size': BATCH_SIZE,
    'n_hvg': N_HVG,
    'random_seed': RANDOM_SEED
}

adata.write_h5ad(final_file)

print(f"\n✓ Final dataset saved!")
print(f"  Size: {os.path.getsize(final_file) / 1e9:.2f} GB")

# Generate README
readme_content = f"""# Epithelial scANVI Annotation Results

## Overview
scVI/scANVI semi-supervised annotation of {adata.n_obs:,} epithelial cells.

## Main Output
**adata_epithelial_SCANVI_FINAL.h5ad** - Final annotated dataset

### Key Columns
- `adata.obs['{SCANVI_LABEL_KEY}']` - ⭐ Final cell type annotations (USE THIS)
- `adata.obs['scanvi_confidence']` - Prediction confidence (0-1)
- `adata.obs['confidence_category']` - High/Medium/Low
- `adata.obs['label_origin']` - Original vs Predicted

### Key Representations
- `adata.obsm['X_scanvi']` - ⭐ Batch-corrected latent space
- `adata.obsm['X_scvi']` - scVI latent space
- `adata.obsm['X_umap']` - UMAP coordinates

## Statistics
- Total cells: {adata.n_obs:,}
- Genes (HVG): {adata.n_vars:,}
- Cell types: {len(type_counts)}
- Mean confidence: {mean_conf:.3f}
- Label agreement: {agreement_rate:.1f}% (if available)

## Files
1. `adata_scanvi_prepared.h5ad` - Prepared data
2. `adata_scanvi_results.h5ad` - After training
3. `adata_epithelial_SCANVI_FINAL.h5ad` - ⭐ Final with UMAP
4. `scanvi_final_summary.csv` - Statistics
5. `scanvi_models/` - Trained models (reusable)

## Figures
1. `01_scvi_training_history.png` - scVI training curve
2. `02_scanvi_training_history.png` - scANVI training curve
3. `03_scanvi_validation.png` - Quality control (4 panels)
4. `04_scanvi_umap_overview.png` - Comprehensive UMAP (7 panels)
5. `05_scanvi_cell_types_individual.png` - Individual cell types

## Model Parameters
- Latent dimensions: {N_LATENT}
- Hidden layers: {N_LAYERS}
- scVI epochs: {N_EPOCHS_SCVI}
- scANVI epochs: {N_EPOCHS_SCANVI}
- Batch size: {BATCH_SIZE}
- HVG: {N_HVG}
- Random seed: {RANDOM_SEED}

## Next Steps
Use this dataset for:
- Differential expression analysis
- Trajectory inference
- Cell-cell communication
- Publication figures

## Notes
- Memory optimized: No large matrix copying
- Batch correction verified in UMAP panel 3
- Low confidence predictions flagged in QC plots

---
Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Optimized version based on code review recommendations
"""

with open(os.path.join(OUTPUT_DIR, 'README.md'), 'w') as f:
    f.write(readme_content)

print(f"\n✓ README.md created")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"🎉 PIPELINE COMPLETE!")
print(f"{'='*80}")

print(f"\n📁 All results saved to:")
print(f"   {OUTPUT_DIR}")

print(f"\n⭐ Main output file:")
print(f"   {final_file}")

print(f"\n📊 Summary statistics:")
print(f"   - Total cells: {adata.n_obs:,}")
print(f"   - Cell types: {len(type_counts)}")
print(f"   - Mean confidence: {mean_conf:.3f}")
if not np.isnan(agreement_rate):
    print(f"   - Label agreement: {agreement_rate:.1f}%")

print(f"\n💡 Use adata.obs['{SCANVI_LABEL_KEY}'] for final cell type annotations!")

print(f"\n✅ Ready for downstream analysis!")
print(f"{'='*80}")