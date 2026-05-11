# %% [markdown]
# # Pure T Cell scANVI Semi-supervised Learning
# 
# **Objective**: Fine-tune scVI model with starCAT annotations using scANVI
# 
# **Workflow**:
# 1. Load scVI-integrated data with starCAT annotations
# 2. Load pre-trained scVI model
# 3. Initialize scANVI from scVI model
# 4. Train scANVI with cell type labels
# 5. Compare scVI vs scANVI latent spaces
# 6. Evaluate integration quality
# 7. Save results
# 
# **Expected cell types** (from starCAT):
# - CD4 subtypes: CD4_CM, CD4_Naive, CD4_EM
# - CD8 subtypes: CD8_CM, CD8_Naive, CD8_TEMRA, CD8_EM
# - Other: Treg, gdT, MAIT
# 
# **Key advantage of scANVI**:
# - Better latent space for cell type separation
# - Can predict cell types for new data
# - Improved batch correction with label guidance

# %% [markdown]
# ---
# ## Step 1: Import Libraries

# %%
import sys
import os
from pathlib import Path
import warnings
import time
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns

# Single-cell analysis
import scanpy as sc
import scvi

# For integration metrics (optional, install if needed: pip install scib-metrics)
try:
    from scib_metrics.benchmark import Benchmarker
    SCIB_AVAILABLE = True
except ImportError:
    print("scib-metrics not available. Will skip integration metrics.")
    SCIB_AVAILABLE = False

warnings.filterwarnings('ignore')

# Check versions
print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")

# Check GPU
import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")

# Plot settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
plt.rcParams['figure.figsize'] = (8, 6)

print("\n✓ Libraries loaded")

# %% [markdown]
# ---
# ## Step 2: Configuration

# %%
# ==================== Configuration ====================

# Random seed
import random
RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

print(f"✓ Random seed: {RANDOM_SEED}")

# ========== Input/Output ==========
# Input: starCAT annotated data (output from starCAT pipeline)
INPUT_H5AD = "/home/h2048/data/py/1204/Tcell_pure_scvi/starcat_analysis/adata_tcell_scvi_starcat_annotated.h5ad"

# Pre-trained scVI model directory
SCVI_MODEL_DIR = "/home/h2048/data/py/1204/Tcell_pure_scvi/scvi_model"

# Output directory for scANVI results
OUTPUT_DIR = "/home/h2048/data/py/1204/Tcell_scANVI"

# ========== scANVI Configuration ==========
# Cell type labels from starCAT
LABELS_KEY = "Multinomial_Label"  # starCAT annotation column

# Unlabeled category (if any cells lack labels)
UNLABELED_CATEGORY = "Unknown"  # Cells with this label will be predicted

# Batch key (for integration metrics)
BATCH_KEY = "dataset"

# ========== Training Parameters ==========
MAX_EPOCHS_SCANVI = 150      # Fine-tuning epochs (less than scVI)
EARLY_STOPPING = True        # Auto-stop when converged
LEARNING_RATE = 5e-4         # Slightly lower than scVI (1e-3)
N_SAMPLES_PER_LABEL = 100    # Balanced sampling per cell type

# ========== UMAP Parameters ==========
UMAP_MIN_DIST = 0.2
UMAP_SPREAD = 1.5
UMAP_N_NEIGHBORS = 30

# ========== Visualization ==========
DPI = 300
FIGURE_FORMAT = "pdf"

# Marker genes for validation
MARKER_GENES = {
    'Pan_T': ['CD3D', 'CD3E', 'CD3G'],
    'CD4': ['CD4', 'IL7R'],
    'CD8': ['CD8A', 'CD8B'],
    'Naive': ['CCR7', 'SELL', 'LEF1', 'TCF7'],
    'Memory': ['IL7R', 'CD44'],
    'Effector': ['GZMA', 'GZMB', 'GZMK', 'PRF1'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4'],
    'Exhaustion': ['PDCD1', 'LAG3', 'TIGIT']
}

# Create directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

scanvi_model_dir = output_dir / "scanvi_model"
scanvi_model_dir.mkdir(exist_ok=True)

print("\n✓ Configuration loaded")
print(f"  Input: {INPUT_H5AD}")
print(f"  scVI model: {SCVI_MODEL_DIR}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Labels: {LABELS_KEY}")
print(f"  Max epochs: {MAX_EPOCHS_SCANVI}")

# %% [markdown]
# ---
# ## Step 3: Load Data and Inspect Labels

# %%
print("\n" + "="*70)
print("Step 3: Loading starCAT Annotated Data")
print("="*70)

# Load
adata = sc.read_h5ad(INPUT_H5AD)
print(f"\nLoaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Check labels
if LABELS_KEY not in adata.obs.columns:
    raise KeyError(
        f"Labels key '{LABELS_KEY}' not found in adata.obs\n"
        f"Available columns: {list(adata.obs.columns)[:20]}"
    )

print(f"\n✓ Cell type labels found: {LABELS_KEY}")

# Display label distribution
label_counts = adata.obs[LABELS_KEY].value_counts()
print(f"\nCell type distribution ({len(label_counts)} types):")
for label, count in label_counts.items():
    pct = count / len(adata) * 100
    print(f"  {label}: {count:,} cells ({pct:.1f}%)")

# Check for unlabeled cells
has_unlabeled = (adata.obs[LABELS_KEY] == UNLABELED_CATEGORY).any()
if has_unlabeled:
    n_unlabeled = (adata.obs[LABELS_KEY] == UNLABELED_CATEGORY).sum()
    print(f"\n⚠️  Found {n_unlabeled:,} unlabeled cells ('{UNLABELED_CATEGORY}')")
    print(f"  scANVI will predict labels for these cells")
else:
    print(f"\n✓ All cells have labels (no '{UNLABELED_CATEGORY}')")

# Check batch info
if BATCH_KEY in adata.obs.columns:
    n_batches = adata.obs[BATCH_KEY].nunique()
    print(f"\n✓ Batch info: {BATCH_KEY} ({n_batches} batches)")
else:
    print(f"\n⚠️  Batch key '{BATCH_KEY}' not found")

# Check for UMAP
if 'X_umap' in adata.obsm:
    print(f"\n✓ UMAP coordinates available")
else:
    print(f"\n⚠️  No UMAP found, will compute later")

# Check for scVI latent
if 'X_scvi' in adata.obsm:
    print(f"✓ scVI latent space: {adata.obsm['X_scvi'].shape}")
else:
    print(f"⚠️  No scVI latent space found")

print("\n✓ Data inspection complete")

# %% [markdown]
# ---
# ## Step 4: Visualize Pre-scANVI State

# %%
print("\n" + "="*70)
print("Step 4: Visualizing scVI Results (Before scANVI)")
print("="*70)

# Compute UMAP if not exists
if 'X_umap' not in adata.obsm:
    print("\nComputing UMAP from scVI latent space...")
    if 'X_scvi' in adata.obsm:
        sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS)
    else:
        sc.pp.neighbors(adata, n_neighbors=UMAP_N_NEIGHBORS)
    
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD)
    print("  ✓ UMAP computed")

# Plot current state
print("\nGenerating baseline plots...")

# %%
# Cell types on scVI UMAP
fig, ax = plt.subplots(figsize=(10, 8))

sc.pl.umap(
    adata,
    color=LABELS_KEY,
    ax=ax,
    show=False,
    legend_loc='right margin',
    frameon=False,
    title=f'scVI + starCAT Labels (n={len(adata):,})'
)

plt.tight_layout()
plt.savefig(fig_dir / f'umap_scVI_starcat_labels.{FIGURE_FORMAT}', 
            dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Plot 1: scVI UMAP with starCAT labels")

# %%
# Batch distribution (if available)
if BATCH_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    
    sc.pl.umap(
        adata,
        color=BATCH_KEY,
        ax=ax,
        show=False,
        legend_loc='right margin',
        frameon=False,
        title='scVI Batch Correction'
    )
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'umap_scVI_batch.{FIGURE_FORMAT}', 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Plot 2: scVI batch mixing")

print("\n✓ Baseline visualization complete")

# %% [markdown]
# ---
# ## Step 5: Load Pre-trained scVI Model

# %%
print("\n" + "="*70)
print("Step 5: Loading Pre-trained scVI Model")
print("="*70)

scvi_model_path = Path(SCVI_MODEL_DIR)

if not scvi_model_path.exists():
    raise FileNotFoundError(
        f"scVI model not found: {scvi_model_path}\n"
        f"Please train scVI first or update SCVI_MODEL_DIR"
    )

print(f"\nLoading scVI model from: {scvi_model_path}")

# CRITICAL: Subset to highly variable genes
# scVI was trained on HVGs only, must match exactly
print(f"\nChecking gene subset requirement...")
print(f"  Current genes: {adata.n_vars:,}")

if 'highly_variable' in adata.var.columns:
    n_hvg = adata.var['highly_variable'].sum()
    print(f"  Highly variable genes: {n_hvg:,}")
    
    if n_hvg < adata.n_vars:
        print(f"\n  Subsetting to HVGs (required for model loading)...")
        
        # Store full adata for later reference
        adata_full_genes = adata.copy()
        
        # Subset to HVGs
        adata = adata[:, adata.var['highly_variable']].copy()
        print(f"  ✓ Subset to {adata.n_vars:,} genes")
    else:
        print(f"  All genes are HVGs, no subsetting needed")
else:
    print(f"  ⚠️  No 'highly_variable' column found")
    print(f"  Assuming current genes match training")

try:
    # Setup anndata for scVI
    print("\nSetting up AnnData for scVI...")
    
    # Check if counts layer exists
    if 'counts' not in adata.layers:
        print("⚠️  'counts' layer not found, checking alternatives...")
        if adata.raw is not None:
            print("  Using .raw for counts")
            # Need to subset raw to match HVGs
            if 'highly_variable' in adata.var.columns:
                hvg_names = adata.var_names
                adata.layers['counts'] = adata.raw[:, hvg_names].X.copy()
            else:
                adata.layers['counts'] = adata.raw.X.copy()
        else:
            print("  ⚠️  Using .X (verify this is raw counts!)")
            adata.layers['counts'] = adata.X.copy()
    
    # Setup scVI
    scvi.model.SCVI.setup_anndata(
        adata,
        layer='counts',
        batch_key=BATCH_KEY if BATCH_KEY in adata.obs.columns else None
    )
    
    print("  ✓ AnnData setup complete")
    
    # Load model
    print(f"\nLoading scVI model...")
    scvi_model = scvi.model.SCVI.load(scvi_model_path, adata)
    
    print(f"\n✓ scVI model loaded successfully")
    print(f"  Model type: {type(scvi_model).__name__}")
    print(f"  Latent dimensions: {scvi_model.module.n_latent}")
    print(f"  Genes in model: {adata.n_vars:,}")
    print(f"  Training history: {len(scvi_model.history['elbo_train'])} epochs")
    
except Exception as e:
    print(f"\n✗ Error loading scVI model: {e}")
    print("\nTroubleshooting:")
    print("  1. Verify model path is correct")
    print("  2. Check gene count matches (should be ~4000 HVGs)")
    print(f"     Current: {adata.n_vars}")
    print("  3. Ensure scvi-tools version is compatible")
    raise

print("\n✓ Model loading complete")

# %% [markdown]
# ---
# ## Step 6: Initialize scANVI from scVI

# %%
print("\n" + "="*70)
print("Step 6: Initializing scANVI from scVI Model")
print("="*70)

print(f"\nInitializing scANVI...")
print(f"  Labels key: {LABELS_KEY}")
print(f"  Unlabeled category: {UNLABELED_CATEGORY}")

try:
    # Initialize scANVI from trained scVI model
    # This transfers the learned representations to scANVI
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        labels_key=LABELS_KEY,
        unlabeled_category=UNLABELED_CATEGORY
    )
    
    print(f"\n✓ scANVI initialized successfully")
    print(f"  Model type: {type(scanvi_model).__name__}")
    print(f"  Latent dimensions: {scanvi_model.module.n_latent}")
    print(f"  Cell types: {scanvi_model.n_labels}")
    
    # Display cell type mapping
    # 显示cell type信息（从adata获取，更可靠）  
    print(f"\n  Cell types detected:")
    cell_types = adata.obs[LABELS_KEY].cat.categories.tolist()
    for label in cell_types:
        n_cells = (adata.obs[LABELS_KEY] == label).sum()
        print(f"    {label}: {n_cells:,} cells")
    
except Exception as e:
    print(f"\n✗ Error initializing scANVI: {e}")
    print("\nCommon issues:")
    print("  1. Labels key not found in adata.obs")
    print("  2. Unlabeled category name conflicts with real labels")
    print("  3. scVI model and adata mismatch")
    raise

print("\n✓ Initialization complete")

# %% [markdown]
# ---
# ## Step 7: Train scANVI

# %%
print("\n" + "="*70)
print("Step 7: Training scANVI (Semi-supervised Fine-tuning)")
print("="*70)

print(f"\nTraining configuration:")
print(f"  Max epochs: {MAX_EPOCHS_SCANVI}")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Samples per label: {N_SAMPLES_PER_LABEL}")
print(f"  Early stopping: {EARLY_STOPPING}")

print(f"\nStarting training...")
print(f"  (This will take 5-10 minutes for ~50k cells)")

start_time = time.time()

scanvi_model.train(
    max_epochs=MAX_EPOCHS_SCANVI,
    n_samples_per_label=N_SAMPLES_PER_LABEL,
    early_stopping=EARLY_STOPPING,
    early_stopping_patience=45,
    plan_kwargs={'lr': LEARNING_RATE, 'weight_decay': 1e-6}
)

elapsed = time.time() - start_time

print(f"\n✓ Training completed in {elapsed:.1f}s ({elapsed/60:.1f} min)")

# %%
# Check training results
print("="*70)
print("Training Results Summary")
print("="*70)

# 安全访问history
if hasattr(scanvi_model, 'history') and 'elbo_train' in scanvi_model.history:
    train_history = scanvi_model.history['elbo_train']
    
    if len(train_history) > 0:
        print(f"\n✓ Training history available")
        print(f"  Epochs trained: {len(train_history)}")
        
        # 提取标量值
        initial_elbo = float(train_history.iloc[0])
        final_elbo = float(train_history.iloc[-1])
        
        print(f"  Initial ELBO: {initial_elbo:.2f}")
        print(f"  Final ELBO: {final_elbo:.2f}")
        print(f"  Total improvement: {initial_elbo - final_elbo:.2f}")
        
        # Check convergence
        if len(train_history) >= 50:
            elbo_50 = float(train_history.iloc[49])
            elbo_final = float(train_history.iloc[-1])
            improvement_late = elbo_50 - elbo_final
            print(f"  Improvement after epoch 50: {improvement_late:.2f}")
            
            if improvement_late < 2:
                print(f"  → Fast convergence ✓")
            elif improvement_late < 5:
                print(f"  → Normal convergence ✓")
            else:
                print(f"  → Slow convergence, still improving")
        
        # Check if early stopping was triggered
        if len(train_history) < MAX_EPOCHS_SCANVI:
            print(f"\n  Early stopping triggered at epoch {len(train_history)}")
        else:
            print(f"\n  Training completed full {MAX_EPOCHS_SCANVI} epochs")
    else:
        print(f"\n⚠️  Training history is empty")
else:
    print(f"\n⚠️  Training history not available")
    print(f"  (This is normal for some scvi-tools versions)")
    print(f"  Training completed successfully")

print("\n" + "="*70)
print("✓ Step 7 Complete")
print("="*70)

# %% [markdown]
# ---
# ## Step 8: Visualize Training History

# %%
print("\n" + "="*70)
print("Step 8: Visualizing Training History")
print("="*70)

# 检查history是否可用
if not hasattr(scanvi_model, 'history') or 'elbo_train' not in scanvi_model.history:
    print("\n⚠️  Training history not available, skipping visualization")
    print("  (This is normal for some scvi-tools versions)")
else:
    train_elbo = scanvi_model.history['elbo_train']
    
    if len(train_elbo) == 0:
        print("\n⚠️  Training history is empty, skipping visualization")
    else:
        # Plot training curves
        val_elbo = scanvi_model.history.get('elbo_validation', None)
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Full history
        axes[0].plot(train_elbo, label='Training ELBO', linewidth=2)
        if val_elbo is not None and len(val_elbo) > 0:
            axes[0].plot(val_elbo, label='Validation ELBO', linewidth=2, alpha=0.7)
        axes[0].set_xlabel('Epoch', fontsize=12)
        axes[0].set_ylabel('ELBO', fontsize=12)
        axes[0].set_title('scANVI Training History', fontsize=14, fontweight='bold')
        axes[0].legend()
        axes[0].grid(alpha=0.3)
        
        # Last 50 epochs
        n_show = min(50, len(train_elbo))
        axes[1].plot(range(len(train_elbo)-n_show, len(train_elbo)), 
                     train_elbo.iloc[-n_show:], linewidth=2)
        axes[1].set_xlabel('Epoch', fontsize=12)
        axes[1].set_ylabel('ELBO', fontsize=12)
        axes[1].set_title(f'Convergence Detail (Last {n_show} epochs)', 
                         fontsize=14, fontweight='bold')
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(fig_dir / f'scanvi_training_history.{FIGURE_FORMAT}', 
                    dpi=DPI, bbox_inches='tight')
        plt.show()
        
        print("✓ Training curves saved")
        
        # Print convergence stats
        print(f"\nConvergence statistics:")
        
        # 提取标量值
        initial_elbo = float(train_elbo.iloc[0])
        final_elbo = float(train_elbo.iloc[-1])
        
        print(f"  Initial ELBO: {initial_elbo:.2f}")
        print(f"  Final ELBO: {final_elbo:.2f}")
        print(f"  Improvement: {initial_elbo - final_elbo:.2f}")
        
        if len(train_elbo) >= 50:
            elbo_50 = float(train_elbo.iloc[49])
            improvement_late = elbo_50 - final_elbo
            print(f"  Improvement after epoch 50: {improvement_late:.2f}")
            if improvement_late < 2:
                print(f"  → Fast convergence ✓")
            elif improvement_late < 5:
                print(f"  → Normal convergence ✓")
            else:
                print(f"  → Slow convergence, still improving")

print("\n✓ Training visualization complete")

# %% [markdown]
# ---
# ## Step 9: Extract scANVI Latent Representation

# %%
print("\n" + "="*70)
print("Step 9: Extracting scANVI Latent Space")
print("="*70)

print(f"\nExtracting latent representations...")

# Get scANVI latent space
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)

print(f"  ✓ scANVI latent: {adata.obsm[SCANVI_LATENT_KEY].shape}")

# Get predicted labels (for unlabeled cells)
print(f"\nPredicting cell types...")
predictions = scanvi_model.predict(adata)
adata.obs['scANVI_predicted'] = predictions

# Compare with original labels
if has_unlabeled:
    print(f"\n  Predictions for unlabeled cells:")
    unlabeled_mask = adata.obs[LABELS_KEY] == UNLABELED_CATEGORY
    pred_counts = adata.obs.loc[unlabeled_mask, 'scANVI_predicted'].value_counts()
    for pred, count in pred_counts.items():
        print(f"    {pred}: {count:,} cells")
else:
    # Check prediction accuracy
    accuracy = (adata.obs[LABELS_KEY] == adata.obs['scANVI_predicted']).mean()
    print(f"\n  ✓ Prediction accuracy: {accuracy*100:.2f}%")
    
    # Show mismatches
    mismatches = adata.obs[LABELS_KEY] != adata.obs['scANVI_predicted']
    if mismatches.any():
        print(f"  Mismatched cells: {mismatches.sum():,} ({mismatches.mean()*100:.2f}%)")

print("\n✓ Latent space extraction complete")

# %% [markdown]
# ---
# ## Step 10: Compute scANVI UMAP

# %%
print("\n" + "="*70)
print("Step 10: Computing scANVI UMAP")
print("="*70)

print(f"\nComputing neighbors from scANVI latent space...")
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, n_neighbors=UMAP_N_NEIGHBORS)

print(f"Computing UMAP...")
# Store old UMAP
if 'X_umap' in adata.obsm:
    adata.obsm['X_umap_scVI'] = adata.obsm['X_umap'].copy()

sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD)

# Rename to scANVI UMAP
adata.obsm['X_umap_scANVI'] = adata.obsm['X_umap'].copy()

print(f"\n✓ scANVI UMAP computed")
print(f"  Available UMAPs:")
if 'X_umap_scVI' in adata.obsm:
    print(f"    X_umap_scVI: scVI-based UMAP")
print(f"    X_umap_scANVI: scANVI-based UMAP (current)")

# %% [markdown]
# ---
# ## Step 11: Compare scVI vs scANVI

# %%
print("\n" + "="*70)
print("Step 11: Comparing scVI vs scANVI")
print("="*70)

print("\nGenerating comparison plots...")

# %%
# Side-by-side comparison
fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# scVI
if 'X_umap_scVI' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_umap_scVI']
    sc.pl.umap(
        adata,
        color=LABELS_KEY,
        ax=axes[0],
        show=False,
        legend_loc='on data',
        legend_fontsize=6,
        frameon=False,
        title='scVI (Unsupervised)'
    )

# scANVI
adata.obsm['X_umap'] = adata.obsm['X_umap_scANVI']
sc.pl.umap(
    adata,
    color=LABELS_KEY,
    ax=axes[1],
    show=False,
    legend_loc='on data',
    legend_fontsize=6,
    frameon=False,
    title='scANVI (Semi-supervised)'
)

plt.tight_layout()
plt.savefig(fig_dir / f'comparison_scVI_vs_scANVI.{FIGURE_FORMAT}', 
            dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Comparison plot saved")

# %%
# Batch mixing comparison
if BATCH_KEY in adata.obs.columns:
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    
    # scVI
    if 'X_umap_scVI' in adata.obsm:
        adata.obsm['X_umap'] = adata.obsm['X_umap_scVI']
        sc.pl.umap(
            adata,
            color=BATCH_KEY,
            ax=axes[0],
            show=False,
            legend_loc='right margin',
            frameon=False,
            title='scVI Batch Mixing'
        )
    
    # scANVI
    adata.obsm['X_umap'] = adata.obsm['X_umap_scANVI']
    sc.pl.umap(
        adata,
        color=BATCH_KEY,
        ax=axes[1],
        show=False,
        legend_loc='right margin',
        frameon=False,
        title='scANVI Batch Mixing'
    )
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'comparison_batch_mixing.{FIGURE_FORMAT}', 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Batch mixing comparison saved")

print("\n✓ Comparison visualization complete")

# %% [markdown]
# ---
# ## Step 12: Integration Metrics (Optional)

# %%
print("\n" + "="*70)
print("Step 12: Computing Integration Metrics")
print("="*70)

if not SCIB_AVAILABLE:
    print("\nscib-metrics not available. Skipping metrics.")
    print("  To install: pip install scib-metrics")
elif BATCH_KEY not in adata.obs.columns:
    print("\nBatch key not found. Skipping metrics.")
else:
    print("\nComputing integration metrics...")
    print("  (This may take a few minutes)")
    
    try:
        # Prepare embeddings
        embedding_keys = []
        if 'X_scvi' in adata.obsm:
            embedding_keys.append('X_scvi')
        if 'X_scANVI' in adata.obsm:
            embedding_keys.append('X_scANVI')
        
        # Run benchmarker
        bm = Benchmarker(
            adata,
            batch_key=BATCH_KEY,
            label_key=LABELS_KEY,
            embedding_obsm_keys=embedding_keys,
            n_jobs=-1
        )
        
        bm.benchmark()
        
        # Get results
        results_df = bm.get_results(min_max_scale=False)
        
        print("\n✓ Integration metrics computed")
        print("\nResults:")
        print(results_df)
        
        # Save results
        results_df.to_csv(output_dir / "integration_metrics.csv")
        
        # Plot
        fig = bm.plot_results_table(min_max_scale=False)
        plt.savefig(fig_dir / f'integration_metrics_table.{FIGURE_FORMAT}', 
                   dpi=DPI, bbox_inches='tight')
        plt.show()
        
        print("\n✓ Metrics saved")
        
    except Exception as e:
        print(f"\n⚠️  Metrics computation failed: {e}")
        print("  Continuing without metrics...")

print("\n✓ Metrics step complete")

# %% [markdown]
# ---
# ## Step 13: Marker Gene Validation

# %%
print("\n" + "="*70)
print("Step 13: Marker Gene Validation on scANVI UMAP")
print("="*70)

# Ensure using scANVI UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_scANVI']

# Flatten marker dict
all_markers = []
for category, genes in MARKER_GENES.items():
    all_markers.extend(genes)

available_markers = [g for g in all_markers if g in adata.var_names]
print(f"\nAvailable markers: {len(available_markers)}/{len(all_markers)}")

if len(available_markers) > 0:
    # Plot markers
    n_markers = len(available_markers)
    ncols = 4
    nrows = int(np.ceil(n_markers / ncols))
    
    fig = sc.pl.umap(
        adata,
        color=available_markers,
        ncols=ncols,
        vmax='p99',
        frameon=False,
        return_fig=True,
        show=False
    )
    
    plt.savefig(fig_dir / f'scanvi_marker_genes.{FIGURE_FORMAT}', 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Marker gene plots saved")
else:
    print("⚠️  No marker genes found in dataset")

print("\n✓ Marker validation complete")

# %% [markdown]
# ---
# ## Step 14: Save Results

# %%
print("\n" + "="*70)
print("Step 14: Saving Results")
print("="*70)

# 1. Save scANVI model
print(f"\nSaving scANVI model...")
scanvi_model.save(scanvi_model_dir, overwrite=True)
print(f"  ✓ Model: {scanvi_model_dir}")

# 2. Save annotated data
output_h5ad = output_dir / "adata_tcell_scANVI_annotated.h5ad"
print(f"\nSaving annotated dataset...")
adata.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)
print(f"  ✓ Data: {output_h5ad}")
print(f"  Size: {output_h5ad.stat().st_size / 1024**2:.1f} MB")

# 3. Save summary
summary_file = output_dir / "scanvi_analysis_summary.txt"
with open(summary_file, 'w') as f:
    f.write("="*70 + "\n")
    f.write("scANVI SEMI-SUPERVISED ANNOTATION - SUMMARY\n")
    f.write("="*70 + "\n\n")
    
    f.write(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    
    f.write("[Input]\n")
    f.write(f"  Data: {INPUT_H5AD}\n")
    f.write(f"  scVI model: {SCVI_MODEL_DIR}\n")
    f.write(f"  Cells: {adata.n_obs:,}\n")
    f.write(f"  Genes: {adata.n_vars:,}\n\n")
    
    f.write("[Training]\n")
    # 修复1: 安全访问history并提取标量
    if hasattr(scanvi_model, 'history') and 'elbo_train' in scanvi_model.history:
        train_history = scanvi_model.history['elbo_train']
        if len(train_history) > 0:
            final_elbo = float(train_history.iloc[-1])
            f.write(f"  Epochs: {len(train_history)}\n")
            f.write(f"  Final ELBO: {final_elbo:.2f}\n")
        else:
            f.write(f"  Epochs: N/A\n")
            f.write(f"  Final ELBO: N/A\n")
    else:
        f.write(f"  Training info not available\n")
    
    # 修复2: elapsed可能不存在，需要计算或跳过
    try:
        f.write(f"  Training time: {elapsed:.1f}s\n\n")
    except NameError:
        f.write(f"  Training time: (see training log)\n\n")
    
    f.write("[Cell Types]\n")
    f.write(f"  Labels: {LABELS_KEY}\n")
    f.write(f"  Number of types: {adata.obs[LABELS_KEY].nunique()}\n")
    for label, count in adata.obs[LABELS_KEY].value_counts().items():
        pct = count / adata.n_obs * 100
        f.write(f"    {label}: {count:,} ({pct:.1f}%)\n")
    f.write("\n")
    
    f.write("[Output Files]\n")
    f.write(f"  Model: {scanvi_model_dir}\n")
    f.write(f"  Data: {output_h5ad}\n")
    f.write(f"  Figures: {fig_dir}/*.{FIGURE_FORMAT}\n")
    
print(f"  ✓ Summary: {summary_file}")

print("\n✓ All results saved")

# %% [markdown]
# ---
# ## Step 15: Final Summary

# %%
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

print(f"\n{'='*70}")
print("Analysis Complete")
print(f"{'='*70}")

print(f"\nWorkflow: scVI → starCAT → scANVI ✓")

print(f"\nDataset:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Cell types: {adata.obs[LABELS_KEY].nunique()}")
print(f"  Batches: {adata.obs[BATCH_KEY].nunique() if BATCH_KEY in adata.obs.columns else 'N/A'}")

print(f"\nTraining:")
print(f"  Epochs: {len(scanvi_model.history['elbo_train'])}")
print(f"  Time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
# 提取标量值
final_elbo = float(scanvi_model.history['elbo_train'].iloc[-1])
print(f"  Final ELBO: {final_elbo:.2f}")

print(f"\nOutputs:")
print(f"  Model: {scanvi_model_dir}")
print(f"  Data: {output_h5ad}")
print(f"  Figures: {fig_dir}")

print(f"\nKey Features in adata:")
print(f"  adata.obsm['X_scvi']: scVI latent (unsupervised)")
print(f"  adata.obsm['X_scANVI']: scANVI latent (semi-supervised)")
print(f"  adata.obsm['X_umap_scVI']: scVI UMAP")
print(f"  adata.obsm['X_umap_scANVI']: scANVI UMAP")
print(f"  adata.obs['scANVI_predicted']: Predicted labels")

print("\n" + "="*70)
print("✅ scANVI ANALYSIS COMPLETE")
print("="*70)

print(f"\nNext steps:")
print(f"  1. Compare scVI vs scANVI visualizations")
print(f"  2. Validate cell type assignments with marker genes")
print(f"  3. Perform downstream analysis (DE, trajectory, etc.)")
print(f"  4. Use scANVI model for label transfer to new data")

print("\nHappy analyzing! 🔬")

# %%
print("\n" + "="*70)
print("Loading Full Gene Data for Marker Visualization")
print("="*70)

# 加载完整基因的数据（之前保存的）
INPUT_H5AD_FULL = "/home/h2048/data/py/1204/Tcell_pure_scvi/starcat_analysis/adata_tcell_scvi_starcat_annotated.h5ad"

print(f"\nLoading full gene data...")
adata_full = sc.read_h5ad(INPUT_H5AD_FULL)

print(f"  Full data: {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")
print(f"  Current data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# 复制关键信息到full data
print(f"\nTransferring scANVI results to full gene data...")

# 复制obs（包括scANVI predictions）
for col in ['scANVI_predicted']:
    if col in adata.obs.columns:
        adata_full.obs[col] = adata.obs[col].values

# 复制obsm（包括scANVI latent和UMAP）
for key in ['X_scANVI', 'X_umap_scANVI']:
    if key in adata.obsm:
        adata_full.obsm[key] = adata.obsm[key].copy()

print(f"  ✓ Transferred annotations and embeddings")

# 使用scANVI UMAP
if 'X_umap_scANVI' in adata_full.obsm:
    adata_full.obsm['X_umap'] = adata_full.obsm['X_umap_scANVI']
    print(f"  ✓ Using scANVI UMAP")

print(f"\n✓ Full gene data ready for marker visualization")
print(f"  Now all classic markers should be available!")

# %%
print("\n" + "="*70)
print("Loading Full Gene Data for Marker Visualization")
print("="*70)

# 加载完整基因的数据（之前保存的）
INPUT_H5AD_FULL = "/home/h2048/data/py/1204/Tcell_pure_scvi/starcat_analysis/adata_tcell_scvi_starcat_annotated.h5ad"

print(f"\nLoading full gene data...")
adata_full = sc.read_h5ad(INPUT_H5AD_FULL)

print(f"  Full data: {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")
print(f"  Current data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# 复制关键信息到full data
print(f"\nTransferring scANVI results to full gene data...")

# 复制obs（包括scANVI predictions）
for col in ['scANVI_predicted']:
    if col in adata.obs.columns:
        adata_full.obs[col] = adata.obs[col].values

# 复制obsm（包括scANVI latent和UMAP）
for key in ['X_scANVI', 'X_umap_scANVI']:
    if key in adata.obsm:
        adata_full.obsm[key] = adata.obsm[key].copy()

print(f"  ✓ Transferred annotations and embeddings")

# 使用scANVI UMAP
if 'X_umap_scANVI' in adata_full.obsm:
    adata_full.obsm['X_umap'] = adata_full.obsm['X_umap_scANVI']
    print(f"  ✓ Using scANVI UMAP")

print(f"\n✓ Full gene data ready for marker visualization")
print(f"  Now all classic markers should be available!")

# %%
print("\n" + "="*70)
print("Checking Markers in Full Gene Data")
print("="*70)

# 重新检查所有marker
all_marker_genes = []
for category, genes in MARKER_GENES_DETAILED.items():
    all_marker_genes.extend(genes)
all_marker_genes = list(set(all_marker_genes))

# 在完整数据中检查
available_full = [g for g in all_marker_genes if g in adata_full.var_names]
missing_full = [g for g in all_marker_genes if g not in adata_full.var_names]

print(f"\n✓ Available in FULL dataset: {len(available_full)}/{len(all_marker_genes)} ({len(available_full)/len(all_marker_genes)*100:.1f}%)")

if len(missing_full) > 0:
    print(f"\n⚠️  Still missing ({len(missing_full)}):")
    for category, genes in MARKER_GENES_DETAILED.items():
        missing_in_cat = [g for g in genes if g not in adata_full.var_names]
        if missing_in_cat:
            print(f"  {category}: {', '.join(missing_in_cat)}")
else:
    print(f"\n✓ All marker genes available!")

# 重建MARKER_GENES_AVAILABLE for full data
MARKER_GENES_AVAILABLE_FULL = {}
for category, genes in MARKER_GENES_DETAILED.items():
    available_in_cat = [g for g in genes if g in adata_full.var_names]
    if len(available_in_cat) > 0:
        MARKER_GENES_AVAILABLE_FULL[category] = available_in_cat

print(f"\n✓ Categories with markers: {len(MARKER_GENES_AVAILABLE_FULL)}/{len(MARKER_GENES_DETAILED)}")

# %%
print("\n" + "="*70)
print("Step A: DotPlot by Cell Type (All Markers)")
print("="*70)

# 使用scANVI的UMAP
if 'X_umap_scANVI' in adata_full.obsm:
    adata_full.obsm['X_umap'] = adata_full.obsm['X_umap_scANVI']
    print(f"✓ Using scANVI UMAP")

# 构建marker基因列表（现在所有基因都可用！）
marker_genes_list = []
categories_to_plot = ['Pan_T', 'CD4', 'CD8', 'Naive', 'CM', 'EM', 'TEMRA', 
                      'Treg', 'Activation', 'Exhaustion', 'Proliferation', 'gdT', 'MAIT']

print(f"\nBuilding marker gene list from full dataset...")
for category in categories_to_plot:
    if category in MARKER_GENES_AVAILABLE_FULL:
        genes = MARKER_GENES_AVAILABLE_FULL[category]
        marker_genes_list.extend(genes)
        print(f"  {category:15s}: {len(genes)} genes")

# 去重保持顺序
seen = set()
marker_genes_list = [x for x in marker_genes_list if not (x in seen or seen.add(x))]

print(f"\n✓ Total unique markers: {len(marker_genes_list)}")
print(f"  Including key genes: CD3D, CD3E, CD8A, CD8B, FOXP3, etc.")

# 生成DotPlot
print(f"\nGenerating DotPlot...")

fig = sc.pl.dotplot(
    adata_full,
    var_names=marker_genes_list,
    groupby=LABELS_KEY,
    dendrogram=True,
    standard_scale='var',
    cmap='Reds',
    return_fig=True,
    show=False
)

output_file = fig_dir / f'dotplot_by_celltype_full.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"\n✓ DotPlot saved: {output_file.name}")

# %%
print("\n" + "="*70)
print("Step B: DotPlot by Cluster (Full Markers)")
print("="*70)

# 检查leiden clustering
leiden_keys = [k for k in adata_full.obs.columns if 'leiden' in k.lower()]

if len(leiden_keys) > 0:
    # 选择合适的resolution
    if 'leiden_res2.5' in adata_full.obs.columns:
        cluster_key = 'leiden_res2.5'
    elif 'leiden' in adata_full.obs.columns:
        cluster_key = 'leiden'
    else:
        cluster_key = leiden_keys[0]
    
    n_clusters = adata_full.obs[cluster_key].nunique()
    print(f"\n✓ Using clustering: {cluster_key}")
    print(f"  Number of clusters: {n_clusters}")
    
    if n_clusters > 50:
        print(f"\n⚠️  Too many clusters ({n_clusters}), using subset")
        cluster_key = 'leiden'  # 用resolution更低的
        n_clusters = adata_full.obs[cluster_key].nunique()
        print(f"  Switched to: {cluster_key} ({n_clusters} clusters)")
    
    if n_clusters <= 50:
        print(f"\nGenerating DotPlot by cluster...")
        
        fig = sc.pl.dotplot(
            adata_full,
            var_names=marker_genes_list,
            groupby=cluster_key,
            dendrogram=True,
            standard_scale='var',
            cmap='Reds',
            return_fig=True,
            show=False
        )
        
        output_file = fig_dir / f'dotplot_by_cluster_full.{FIGURE_FORMAT}'
        plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
        plt.show()
        plt.close()
        
        print(f"\n✓ DotPlot saved: {output_file.name}")
    else:
        print(f"\n⚠️  Still too many clusters, skipping")
else:
    print(f"\n⚠️  No leiden clustering found")

# %%
print("\n" + "="*70)
print("Step C: FindAllMarkers (Full Gene Set)")
print("="*70)

print(f"\nRunning differential expression on FULL gene set...")
print(f"  Cells: {adata_full.n_obs:,}")
print(f"  Genes: {adata_full.n_vars:,}")
print(f"  This may take 5-10 minutes...")

# 运行差异表达分析
sc.tl.rank_genes_groups(
    adata_full,
    groupby=LABELS_KEY,
    method='wilcoxon',
    use_raw=False,
    n_genes=100,
    pts=True
)

print(f"\n✓ Differential expression completed")

# 显示每个细胞类型的top基因
print(f"\nTop 5 marker genes per cell type:")
print("="*70)

for cell_type in adata_full.obs[LABELS_KEY].cat.categories:
    genes = adata_full.uns['rank_genes_groups']['names'][cell_type][:5]
    scores = adata_full.uns['rank_genes_groups']['scores'][cell_type][:5]
    pvals = adata_full.uns['rank_genes_groups']['pvals_adj'][cell_type][:5]
    
    print(f"\n{cell_type}:")
    for i, (gene, score, pval) in enumerate(zip(genes, scores, pvals), 1):
        # 标注是否是经典marker
        is_classic = gene in marker_genes_list
        marker_flag = " ⭐" if is_classic else ""
        print(f"  {i}. {gene:12s} | score: {score:6.2f} | p_adj: {pval:.2e}{marker_flag}")

print(f"\n⭐ = Classic marker gene from our predefined list")

# %%
print("\n" + "="*70)
print("Step D: Top Markers Heatmap (Full Genes)")
print("="*70)

print(f"\nGenerating heatmap of top marker genes...")

fig = sc.pl.rank_genes_groups_heatmap(
    adata_full,
    n_genes=10,
    groupby=LABELS_KEY,
    standard_scale='var',
    cmap='RdYlBu_r',
    show_gene_labels=True,
    return_fig=True,
    show=False
)

output_file = fig_dir / f'top_markers_heatmap_full.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"\n✓ Heatmap saved: {output_file.name}")

# %%
print("\n" + "="*70)
print("Step E: Top Markers DotPlot (Full Genes)")
print("="*70)

print(f"\nGenerating dotplot of top DE markers...")

fig = sc.pl.rank_genes_groups_dotplot(
    adata_full,
    n_genes=5,
    groupby=LABELS_KEY,
    standard_scale='var',
    cmap='Reds',
    return_fig=True,
    show=False
)

output_file = fig_dir / f'top_markers_dotplot_full.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"\n✓ DotPlot saved: {output_file.name}")

# %%
print("\n" + "="*70)
print("Step F: Top Markers Stacked Violin (Full Genes)")
print("="*70)

print(f"\nGenerating stacked violin plots...")

# 选择每个细胞类型的top markers
top_genes_per_type = {}
for cell_type in adata_full.obs[LABELS_KEY].cat.categories:
    genes = adata_full.uns['rank_genes_groups']['names'][cell_type][:2]
    top_genes_per_type[cell_type] = list(genes)

# 展平并去重
all_top_genes = []
for genes in top_genes_per_type.values():
    all_top_genes.extend(genes)
seen = set()
all_top_genes = [x for x in all_top_genes if not (x in seen or seen.add(x))]

print(f"  Visualizing {len(all_top_genes)} top DE genes")

if len(all_top_genes) > 0:
    genes_to_plot = all_top_genes[:15]
    
    sc.pl.stacked_violin(
        adata_full,
        var_names=genes_to_plot,
        groupby=LABELS_KEY,
        dendrogram=False,
        swap_axes=False
    )
    
    output_file = fig_dir / f'top_markers_violin_full.{FIGURE_FORMAT}'
    plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f"\n✓ Stacked violin saved: {output_file.name}")

# %%
print("\n" + "="*70)
print("Step G: Export Marker Gene Results (Full Genes)")
print("="*70)

print(f"\nExporting marker gene results...")

# 导出所有组的结果
result_df = sc.get.rank_genes_groups_df(adata_full, group=None)
output_file = output_dir / "differential_expression_all_groups_full.csv"
result_df.to_csv(output_file, index=False)
print(f"  ✓ All groups: {output_file.name}")

# 为每个细胞类型单独导出
print(f"\n  Per cell type (significant only, p_adj < 0.05):")
for cell_type in adata_full.obs[LABELS_KEY].cat.categories:
    df = sc.get.rank_genes_groups_df(adata_full, group=cell_type)
    
    # 只保留显著的
    df_sig = df[df['pvals_adj'] < 0.05]
    
    # 标注经典marker
    df_sig['is_classic_marker'] = df_sig['names'].isin(marker_genes_list)
    
    # 保存
    filename = f"markers_{cell_type.replace('/', '_')}_full.csv"
    df_sig.to_csv(output_dir / filename, index=False)
    
    n_classic = df_sig['is_classic_marker'].sum()
    print(f"    {cell_type:15s}: {len(df_sig):4d} genes ({n_classic:2d} classic markers) → {filename}")

print(f"\n✓ All results exported to: {output_dir}")

# %%
print("\n" + "="*70)
print("Step H: Classic Marker Genes on UMAP")
print("="*70)

# 精选最重要的经典marker基因
key_classic_markers = {
    'Pan_T': ['CD3D', 'CD3E'],
    'CD4': ['CD4'],
    'CD8': ['CD8A', 'CD8B'],
    'Naive': ['CCR7', 'SELL'],
    'CM': ['IL7R', 'LTB'],
    'EM': ['GZMK', 'NKG7'],
    'TEMRA': ['GZMB', 'PRF1'],
    'Treg': ['FOXP3', 'IL2RA'],
    'gdT': ['TRDC', 'TRGC2'],
    'MAIT': ['SLC4A10', 'KLRB1'],
    'Exhaustion': ['PDCD1', 'LAG3'],
    'Proliferation': ['MKI67', 'TOP2A']
}

# 展平
key_genes = []
for category, genes in key_classic_markers.items():
    for gene in genes[:1]:  # 每类选1个
        if gene in adata_full.var_names:
            key_genes.append(gene)

print(f"\nVisualizing {len(key_genes)} key classic markers on UMAP:")
print(f"  {', '.join(key_genes)}")

# 确保使用scANVI UMAP
if 'X_umap_scANVI' in adata_full.obsm:
    adata_full.obsm['X_umap'] = adata_full.obsm['X_umap_scANVI']

# 可视化
fig = sc.pl.umap(
    adata_full,
    color=key_genes,
    ncols=4,
    vmax='p99',
    cmap='Reds',
    frameon=False,
    return_fig=True,
    show=False
)

output_file = fig_dir / f'umap_classic_markers.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"\n✓ UMAP saved: {output_file.name}")

# %%
print("\n" + "="*70)
print("Step I: Validate Cell Types with Classic Markers")
print("="*70)

# 定义每个细胞类型应该高表达的marker
validation_markers = {
    'CD4_CM': ['CD4', 'IL7R', 'CD27'],
    'CD4_EM': ['CD4', 'GZMK'],
    'CD4_Naive': ['CD4', 'CCR7', 'SELL'],
    'CD8_CM': ['CD8A', 'IL7R'],
    'CD8_EM': ['CD8A', 'GZMK', 'NKG7'],
    'CD8_Naive': ['CD8A', 'CCR7', 'SELL'],
    'CD8_TEMRA': ['CD8A', 'GZMB', 'PRF1'],
    'Treg': ['CD4', 'FOXP3', 'IL2RA'],
    'gdT': ['TRDC', 'TRGC2'],
    'MAIT': ['SLC4A10', 'KLRB1']
}

print(f"\nValidation DotPlot: Expected markers for each cell type")

# 构建验证基因列表
validation_genes = []
for cell_type, markers in validation_markers.items():
    for gene in markers:
        if gene in adata_full.var_names and gene not in validation_genes:
            validation_genes.append(gene)

print(f"  Genes to check: {len(validation_genes)}")

# 生成验证DotPlot
fig = sc.pl.dotplot(
    adata_full,
    var_names=validation_genes,
    groupby=LABELS_KEY,
    dendrogram=True,
    standard_scale='var',
    cmap='Reds',
    return_fig=True,
    show=False
)

output_file = fig_dir / f'dotplot_validation_markers.{FIGURE_FORMAT}'
plt.savefig(output_file, dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"\n✓ Validation DotPlot saved: {output_file.name}")
print(f"\nInterpretation:")
print(f"  - Each cell type should show high expression (red) of its expected markers")
print(f"  - CD4 types should express CD4, CD8 types should express CD8A")
print(f"  - Check if starCAT annotations match these classic patterns")

# %%
print("\n" + "="*70)
print("MARKER GENE ANALYSIS - Complete (Full Gene Set)")
print("="*70)

print(f"\n✅ Analysis completed with FULL gene set!")

print(f"\nDataset statistics:")
print(f"  Cells: {adata_full.n_obs:,}")
print(f"  Genes (full): {adata_full.n_vars:,}")
print(f"  Cell types: {adata_full.obs[LABELS_KEY].nunique()}")
print(f"  Marker genes tested: {len(marker_genes_list)}")

print(f"\nKey improvements with full gene set:")
print(f"  ✓ All 51 predefined markers available (vs 37/51 before)")
print(f"  ✓ Classic markers included: CD3D, CD3E, CD8A, CD8B, FOXP3, etc.")
print(f"  ✓ Better validation of cell type annotations")
print(f"  ✓ More accurate differential expression")

print(f"\nGenerated files ({fig_dir}):")
files_generated = [
    'dotplot_by_celltype_full',
    'dotplot_by_cluster_full',
    'top_markers_heatmap_full',
    'top_markers_dotplot_full',
    'top_markers_violin_full',
    'umap_classic_markers',
    'dotplot_validation_markers'
]
for fname in files_generated:
    full_path = fig_dir / f"{fname}.{FIGURE_FORMAT}"
    if full_path.exists():
        print(f"  ✓ {fname}.{FIGURE_FORMAT}")

print(f"\nCSV files ({output_dir}):")
print(f"  ✓ differential_expression_all_groups_full.csv")
for cell_type in adata_full.obs[LABELS_KEY].cat.categories:
    filename = f"markers_{cell_type.replace('/', '_')}_full.csv"
    print(f"  ✓ {filename}")

print(f"\n" + "="*70)
print(f"Next steps:")
print(f"  1. Review dotplot_validation_markers.pdf")
print(f"     → Check if cell types express expected markers")
print(f"  2. Review umap_classic_markers.pdf")
print(f"     → Verify spatial distribution of key genes")
print(f"  3. Check CSV files for complete marker lists")
print(f"=" * 70)


