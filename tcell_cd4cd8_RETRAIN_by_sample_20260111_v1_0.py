# ============================================================================
# T Cell CD4/CD8 scVI/scANVI Re-training by Sample (v1.0 - FIXED v2)
# ============================================================================
# Purpose: Re-train scVI and scANVI models using 'Sample' as batch_key
#          (instead of 'dataset')
# 
# Author: r2end
# Date: 2025-01-11
# Fixed: 2025-01-12 (v2)
#   - Fix 1: Robust parameter counting (QUICK_REFERENCE_MEMORY Pitfall 8)
#   - Fix 2: Replace use_gpu with accelerator/devices (QUICK_REFERENCE_MEMORY Error 0)
# 
# CRITICAL CHANGE:
# - Old: BATCH_KEY = 'dataset'
# - New: BATCH_KEY = 'Sample'
# 
# FIXES APPLIED:
# - ✓ Robust parameter counting compatible with all scvi-tools versions
# - ✓ Modern training API (accelerator/devices) for scvi-tools ≥1.1
# ============================================================================

# %% [markdown]
## Imports

# %%
import os
import sys
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Headless server compatible
import matplotlib.pyplot as plt
import time
import gc

import scvi
import torch

print("Libraries imported successfully")

# %% [markdown]
## Configuration

# %%
# ============================================================================
# INPUT/OUTPUT PATHS
# ============================================================================

# ⚠️ UPDATE THESE PATHS
INPUT_FILE = "/home/h2048/data/py/0109/tcell_cd4cd8_scanvi_v1_3/adata_tcell_cd4cd8_FINAL.h5ad"  # Your final annotated file
OUTPUT_DIR = f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/cd4cd8_tcell"

# ============================================================================
# TRAINING PARAMETERS
# ============================================================================

BATCH_KEY = 'sample'  # ⭐ KEY CHANGE: Sample instead of dataset

# scVI parameters
SCVI_N_LATENT = 75
SCVI_N_HIDDEN = 128
SCVI_N_LAYERS = 2
SCVI_DROPOUT = 0.1
SCVI_MAX_EPOCHS = 400

# scANVI parameters
SCANVI_MAX_EPOCHS = 600
SCANVI_UNLABELED_CATEGORY = 'Unknown'

# Performance
BATCH_SIZE = 2048
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 50

# Label columns
LABEL_KEY = 'cell_type_scanvi_cd4cd8_filt'  # ⚠️ Adjust if different in your file

# GPU/CPU
USE_GPU = True

RANDOM_SEED = 42

# Set seeds
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

# Create output directory
OUTPUT_DIR = Path(OUTPUT_DIR)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "models").mkdir(exist_ok=True)

print("="*80)
print("T CELL CD4/CD8 RE-TRAINING BY SAMPLE (FIXED v2)")
print("="*80)
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Batch key: {BATCH_KEY} ⭐")
print(f"  Label key: {LABEL_KEY}")
print(f"  GPU available: {torch.cuda.is_available()}")
print(f"  Random seed: {RANDOM_SEED}")

# %% [markdown]
## STAGE 1: Load Data

# %%
print(f"\n{'='*80}")
print("STAGE 1: LOAD DATA")
print("="*80)

adata = sc.read_h5ad(INPUT_FILE)

print(f"\n📊 Data loaded:")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")

# Check critical columns
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"❌ Batch key '{BATCH_KEY}' not found in adata.obs")
if LABEL_KEY not in adata.obs.columns:
    raise ValueError(f"❌ Label key '{LABEL_KEY}' not found in adata.obs")

print(f"\n✓ Batch column: {BATCH_KEY}")
print(f"   Unique batches: {adata.obs[BATCH_KEY].nunique()}")
print(f"   Batch sizes:")
batch_counts = adata.obs[BATCH_KEY].value_counts()
print(f"      Mean: {batch_counts.mean():.0f} cells")
print(f"      Median: {batch_counts.median():.0f} cells")
print(f"      Range: {batch_counts.min()}-{batch_counts.max()} cells")

print(f"\n✓ Label column: {LABEL_KEY}")
print(f"   Unique cell types: {adata.obs[LABEL_KEY].nunique()}")
label_counts = adata.obs[LABEL_KEY].value_counts()
print(f"   Top 10 types:")
for i, (label, count) in enumerate(label_counts.head(10).items(), 1):
    pct = count / adata.n_obs * 100
    print(f"      {i:2d}. {label:40s} {count:6,} ({pct:5.1f}%)")

# Check for required data
if adata.X is None or adata.X.shape[0] == 0:
    raise ValueError("❌ No expression data in adata.X")

print(f"\n✓ Expression matrix:")
print(f"   Shape: {adata.X.shape}")
print(f"   Type: {type(adata.X)}")
print(f"   Dtype: {adata.X.dtype if hasattr(adata.X, 'dtype') else 'N/A'}")

# %% [markdown]
## STAGE 2: Check/Prepare HVG Subset

# %%
print(f"\n{'='*80}")
print("STAGE 2: HVG PREPARATION")
print("="*80)

# Check if already HVG subset
if adata.raw is not None:
    print(f"\n✓ .raw detected: {adata.raw.n_vars:,} genes")
    print(f"   Current subset: {adata.n_vars:,} genes")
    
    if 'highly_variable' in adata.var.columns:
        n_hvg = adata.var['highly_variable'].sum()
        print(f"   HVG marker in .var: {n_hvg} genes marked")
        
        if n_hvg == adata.n_vars:
            print(f"   ✓ Already subsetted to HVG")
        else:
            print(f"   ⚠️ Warning: HVG count mismatch")
    else:
        print(f"   ℹ️ No 'highly_variable' column, assuming current is HVG subset")
else:
    print(f"\n⚠️ No .raw detected")
    print(f"   Will use current {adata.n_vars:,} genes for training")
    print(f"   Note: Training on all genes may require more memory")

# Memory estimation
n_cells = adata.n_obs
n_genes = adata.n_vars
approx_mem_gb = (n_cells * n_genes * 4) / 1e9  # Rough estimate
print(f"\n💾 Approximate memory requirement: {approx_mem_gb:.1f} GB")
if approx_mem_gb > 40:
    print(f"   ⚠️ Warning: High memory usage expected")

# %% [markdown]
## STAGE 3: scVI Training

# %%
print(f"\n{'='*80}")
print("STAGE 3: scVI TRAINING")
print("="*80)

# Setup scVI
print(f"\nSetting up scVI model...")
print(f"   Batch key: {BATCH_KEY}")
print(f"   n_latent: {SCVI_N_LATENT}")
print(f"   n_hidden: {SCVI_N_HIDDEN}")
print(f"   n_layers: {SCVI_N_LAYERS}")

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key=BATCH_KEY
)

# Create model
vae = scvi.model.SCVI(
    adata,
    n_latent=SCVI_N_LATENT,
    n_hidden=SCVI_N_HIDDEN,
    n_layers=SCVI_N_LAYERS,
    dropout_rate=SCVI_DROPOUT,
    gene_likelihood="nb"
)

print(f"✓ Model created")

# ⭐ FIXED: Robust parameter counting (QUICK_REFERENCE_MEMORY Pitfall 8)
try:
    n_params = vae.module.n_params
except AttributeError:
    # Fallback: PyTorch standard method
    n_params = sum(p.numel() for p in vae.module.parameters() if p.requires_grad)
print(f"   Trainable parameters: {n_params:,}")

# Train
print(f"\nTraining scVI...")
print(f"   Max epochs: {SCVI_MAX_EPOCHS}")
print(f"   Batch size: {BATCH_SIZE}")
print(f"   Early stopping: {EARLY_STOPPING}")
print(f"   Accelerator: {'gpu' if USE_GPU else 'cpu'}")

t0 = time.time()

# ⭐ FIXED: Replace use_gpu with accelerator/devices (scvi-tools ≥1.1)
train_kwargs = {
    'max_epochs': SCVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'early_stopping': EARLY_STOPPING,
    'early_stopping_patience': EARLY_STOPPING_PATIENCE,
}

if USE_GPU:
    train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
else:
    train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})

vae.train(**train_kwargs)

train_time = time.time() - t0
print(f"\n✓ Training complete")
print(f"   Time: {train_time/60:.1f} minutes")

# Save model
scvi_model_path = OUTPUT_DIR / "models" / "scvi_model_by_sample"
print(f"\nSaving scVI model...")
vae.save(scvi_model_path, overwrite=True)
print(f"✓ Model saved: {scvi_model_path}")

# Get latent representation
print(f"\nComputing latent representation...")
latent = vae.get_latent_representation()
adata.obsm['X_scvi'] = latent
print(f"✓ Latent space: {latent.shape}")

# Optional: Basic visualization
print(f"\nComputing UMAP from scVI latent...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
sc.tl.umap(adata)
print(f"✓ UMAP computed")

# Save checkpoint
checkpoint_file = OUTPUT_DIR / "adata_after_scvi.h5ad"
print(f"\nSaving checkpoint...")
adata.write_h5ad(checkpoint_file, compression='gzip')
print(f"✓ Checkpoint saved: {checkpoint_file.name}")

# Clear memory
del latent
gc.collect()

# %% [markdown]
## STAGE 4: scANVI Training

# %%
print(f"\n{'='*80}")
print("STAGE 4: scANVI TRAINING")
print("="*80)

# Check labels
print(f"\nPreparing labels...")
print(f"   Label column: {LABEL_KEY}")

# Create clean labels for scANVI (must be strings, no NaN)
labels = adata.obs[LABEL_KEY].astype(str).copy()
labels = labels.replace({'nan': SCANVI_UNLABELED_CATEGORY, 
                         'None': SCANVI_UNLABELED_CATEGORY,
                         '': SCANVI_UNLABELED_CATEGORY})

n_labeled = (labels != SCANVI_UNLABELED_CATEGORY).sum()
n_total = len(labels)
print(f"   Labeled cells: {n_labeled:,} ({n_labeled/n_total*100:.1f}%)")
print(f"   Unlabeled: {n_total - n_labeled:,}")
print(f"   Unique types: {labels.nunique()}")

# Temporary column for scANVI
adata.obs['scanvi_labels_temp'] = labels

# Setup scANVI from trained scVI
print(f"\nInitializing scANVI from scVI...")
print(f"   Unlabeled category: {SCANVI_UNLABELED_CATEGORY}")

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key='scanvi_labels_temp',
    unlabeled_category=SCANVI_UNLABELED_CATEGORY
)

print(f"✓ scANVI model initialized")

# ⭐ FIXED: Robust parameter counting
try:
    n_params = lvae.module.n_params
except AttributeError:
    n_params = sum(p.numel() for p in lvae.module.parameters() if p.requires_grad)
print(f"   Trainable parameters: {n_params:,}")

# Train
print(f"\nTraining scANVI...")
print(f"   Max epochs: {SCANVI_MAX_EPOCHS}")
print(f"   Batch size: {BATCH_SIZE}")
print(f"   Accelerator: {'gpu' if USE_GPU else 'cpu'}")

t0 = time.time()

# ⭐ FIXED: Replace use_gpu with accelerator/devices (scvi-tools ≥1.1)
train_kwargs = {
    'max_epochs': SCANVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'early_stopping': EARLY_STOPPING,
    'early_stopping_patience': EARLY_STOPPING_PATIENCE,
}

if USE_GPU:
    train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
else:
    train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})

lvae.train(**train_kwargs)

train_time = time.time() - t0
print(f"\n✓ Training complete")
print(f"   Time: {train_time/60:.1f} minutes")

# Save model
scanvi_model_path = OUTPUT_DIR / "models" / "scanvi_model_by_sample"
print(f"\nSaving scANVI model...")
lvae.save(scanvi_model_path, overwrite=True)
print(f"✓ Model saved: {scanvi_model_path}")

# Get predictions
print(f"\nGenerating predictions...")
predictions = lvae.predict()
adata.obs['cell_type_scanvi_by_sample'] = pd.Categorical(predictions)

# Get latent representation
latent_scanvi = lvae.get_latent_representation()
adata.obsm['X_scanvi'] = latent_scanvi
print(f"✓ Latent space: {latent_scanvi.shape}")

# Compute UMAP from scANVI latent
print(f"\nComputing UMAP from scANVI latent...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)
sc.tl.umap(adata)
print(f"✓ UMAP computed")

# Clean up temporary column
adata.obs.drop(columns=['scanvi_labels_temp'], inplace=True)

# Clear memory
del latent_scanvi
gc.collect()

# %% [markdown]
## STAGE 5: Save Final Results

# %%
print(f"\n{'='*80}")
print("STAGE 5: SAVE FINAL RESULTS")
print("="*80)

# Add metadata
adata.uns['retrain_info'] = {
    'version': '1.0_retrain_by_sample_FIXED_v2',
    'date': '2025-01-12',
    'batch_key': BATCH_KEY,
    'original_label_key': LABEL_KEY,
    'new_prediction_key': 'cell_type_scanvi_by_sample',
    'scvi_params': {
        'n_latent': SCVI_N_LATENT,
        'n_hidden': SCVI_N_HIDDEN,
        'n_layers': SCVI_N_LAYERS,
        'dropout': SCVI_DROPOUT,
        'max_epochs': SCVI_MAX_EPOCHS
    },
    'scanvi_params': {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'unlabeled_category': SCANVI_UNLABELED_CATEGORY
    },
    'random_seed': RANDOM_SEED,
    'fixes_applied': [
        'Robust parameter counting (QUICK_REFERENCE_MEMORY Pitfall 8)',
        'Replace use_gpu with accelerator/devices (QUICK_REFERENCE_MEMORY Error 0)'
    ]
}

# Save final h5ad
output_file = OUTPUT_DIR / "adata_tcell_RETRAINED_by_sample.h5ad"
print(f"\nSaving final h5ad...")
adata.write_h5ad(output_file, compression='gzip')
size_gb = output_file.stat().st_size / 1e9
print(f"✓ Saved: {output_file.name}")
print(f"   Size: {size_gb:.2f} GB")

# Export predictions
print(f"\nExporting predictions...")
predictions_df = adata.obs[[
    BATCH_KEY,
    LABEL_KEY,
    'cell_type_scanvi_by_sample'
]].copy()
predictions_df.to_csv(OUTPUT_DIR / "predictions_by_sample.csv")
print(f"✓ Predictions saved")

# Summary statistics
print(f"\n📊 Prediction Summary:")
pred_counts = adata.obs['cell_type_scanvi_by_sample'].value_counts()
print(f"   Unique types: {len(pred_counts)}")
print(f"   Top 10 predictions:")
for i, (ct, count) in enumerate(pred_counts.head(10).items(), 1):
    pct = count / adata.n_obs * 100
    print(f"      {i:2d}. {ct:40s} {count:6,} ({pct:5.1f}%)")

# Compare with original labels if available
if LABEL_KEY in adata.obs.columns:
    print(f"\n🔄 Comparison with Original Labels:")
    agreement = (adata.obs[LABEL_KEY] == adata.obs['cell_type_scanvi_by_sample']).sum()
    agreement_pct = agreement / adata.n_obs * 100
    print(f"   Agreement: {agreement:,} / {adata.n_obs:,} ({agreement_pct:.1f}%)")

# %% [markdown]
## STAGE 6: Generate Visualizations

# %%
print(f"\n{'='*80}")
print("STAGE 6: VISUALIZATIONS")
print("="*80)

# UMAP colored by predictions
print(f"\nGenerating UMAP plots...")

fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# Original labels
sc.pl.umap(adata, color=LABEL_KEY, ax=axes[0], show=False, title='Original Labels')

# New predictions
sc.pl.umap(adata, color='cell_type_scanvi_by_sample', ax=axes[1], show=False, 
          title='scANVI Predictions (by Sample)')

plt.tight_layout()
plt.savefig(OUTPUT_DIR / "umap_comparison_by_sample.pdf", dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: umap_comparison_by_sample.pdf")

# Batch integration check
print(f"\nGenerating batch integration plot...")
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, 
          title=f'Batch Integration ({BATCH_KEY})')
plt.savefig(OUTPUT_DIR / "umap_batch_integration.pdf", dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: umap_batch_integration.pdf")

# %% [markdown]
## Final Summary

# %%
print(f"\n{'='*80}")
print("✨ RE-TRAINING COMPLETE")
print("="*80)

print(f"\n📊 Summary:")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")
print(f"   Batches ({BATCH_KEY}): {adata.obs[BATCH_KEY].nunique()}")
print(f"   Original types ({LABEL_KEY}): {adata.obs[LABEL_KEY].nunique()}")
print(f"   Predicted types: {adata.obs['cell_type_scanvi_by_sample'].nunique()}")

print(f"\n📁 Output Files:")
print(f"   Main data: {output_file.name}")
print(f"   Predictions: predictions_by_sample.csv")
print(f"   scVI model: models/scvi_model_by_sample/")
print(f"   scANVI model: models/scanvi_model_by_sample/")
print(f"   Checkpoint: adata_after_scvi.h5ad")

print(f"\n📈 Visualizations:")
print(f"   umap_comparison_by_sample.pdf")
print(f"   umap_batch_integration.pdf")

print(f"\n🔑 New Column:")
print(f"   cell_type_scanvi_by_sample: Final predictions ⭐")

print(f"\n💡 Key Differences:")
print(f"   Old batch key: 'dataset'")
print(f"   New batch key: 'Sample' ⭐")
print(f"   This accounts for sample-level batch effects")
print(f"   instead of dataset-level batch effects")

print(f"\n🔧 Fixes Applied:")
print(f"   1. Robust parameter counting for scVI/scANVI models")
print(f"      (QUICK_REFERENCE_MEMORY Pitfall 8)")
print(f"   2. Modern training API: accelerator/devices instead of use_gpu")
print(f"      (QUICK_REFERENCE_MEMORY Error 0)")
print(f"   Compatible with scvi-tools ≥1.1 and all future versions")

print("="*80)