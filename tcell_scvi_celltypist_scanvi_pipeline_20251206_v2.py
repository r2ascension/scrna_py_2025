# %% [markdown]
# # T/NK Cell Analysis: scVI → CellTypist → scANVI Pipeline v2.0
# 
# **Author:** Clinical-Bioinformatics Team  
# **Date:** 2025-12-06  
# **Version:** 2.0 (Critical Fixes Applied)
# 
# ## 📋 Pipeline Overview
# 
# ```
# Input: adata_T_bbknn.h5ad (~40k T/NK cells)
#   ↓
# Step 1: Data Preparation
#   - Load BBKNN preprocessed data
#   - Copy layers['counts'] → .raw
#   - Verify data structure
#   ↓
# Step 2: scVI Integration (Batch Correction)
#   - Train scVI model on raw counts
#   - Generate X_scvi latent space (75 dims)
#   - Compute X_umap_scvi for visualization ⭐ FIXED
#   ↓
# Step 3: CellTypist Annotation
#   - Auto-annotate using Immune_All_Low.pkl
#   - Simplified gene name conversion (ENSEMBL→SYMBOL) ⭐ FIXED
#   - Generate predicted_labels + majority_voting
#   ↓
# Step 4: scANVI Refinement (Semi-supervised)
#   - Use CellTypist labels as reference
#   - Train scANVI for fine-grained classification
#   - Generate X_umap_scanvi (separate from scVI) ⭐ FIXED
#   ↓
# Step 5: Comprehensive Visualization
#   - Compare batch mixing (scVI vs scANVI UMAPs) ⭐ FIXED
#   - Cell type composition analysis
#   - Confidence score evaluation
#   ↓
# Output: Fully annotated h5ad with multiple embeddings
# ```
# 
# ## 🆕 What's New in v2.0
# 
# ### Critical Fixes:
# 1. **UMAP Overwriting Fixed:** scVI and scANVI UMAPs now correctly saved as separate embeddings
# 2. **CellTypist Gene Detection:** Simplified logic, no longer relies on model internals
# 3. **CellTypist Input:** Explicitly uses `layers['counts']`, no guessing
# 4. **Reproducibility:** Added `scvi.settings.seed` for consistent results
# 5. **Gene Names:** Added `var_names_make_unique()` after conversion
# 6. **Probabilities Saved:** scANVI posterior probabilities now stored in `.obsm`
# 
# See `CHANGELOG_v2.md` for full details.
# 
# ## ⚠️ Requirements
# 
# ```bash
# # Core packages
# pip install scvi-tools scanpy celltypist mygene
# 
# # Download CellTypist model
# celltypist.models.download_models(model='Immune_All_Low.pkl')
# ```

# %% [markdown]
# ---
# ## Step 0: Import Libraries and Configuration

# %%
# Import core libraries
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
from scipy.stats import entropy
from sklearn.neighbors import NearestNeighbors

import matplotlib.pyplot as plt
import seaborn as sns

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

print("\n✓ All libraries loaded successfully")

# %%
# ==================== Configuration Parameters ====================

# Random seed for reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

# ⭐ NEW in v2.0: Set scvi-tools seed for full reproducibility
scvi.settings.seed = RANDOM_SEED
print(f"✓ Random seeds set: {RANDOM_SEED}")

# ========== Input/Output Configuration ==========
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1206/tcell_scvi_celltypist_scanvi_v2"

# CellTypist model path (预留路径，请提前下载模型)
CELLTYPIST_MODEL_PATH = "/home/h2048/data/source/reference/celltypist_models/Immune_All_Low.pkl"
# 如果模型已在默认路径，可以直接用模型名：
# CELLTYPIST_MODEL_PATH = "Immune_All_Low.pkl"

# ========== Batch Configuration ==========
BATCH_KEY = "dataset"
CONFOUNDER_KEYS = ["tissue"]  # Biological covariates to preserve

# ========== scVI Model Parameters ==========
SCVI_N_LATENT = 75
SCVI_N_LAYERS = 3
SCVI_DROPOUT = 0.1
SCVI_MAX_EPOCHS = 400
SCVI_EARLY_STOPPING = True
SCVI_LEARNING_RATE = 1e-3

# ========== CellTypist Parameters ==========
CELLTYPIST_MAJORITY_VOTING = True
CELLTYPIST_MODE = "best match"  # or "prob match"

# ========== scANVI Parameters ==========
SCANVI_MAX_EPOCHS = 300
SCANVI_EARLY_STOPPING = True
UNLABELED_CATEGORY = "Unknown"  # For low-confidence predictions
LOW_CONF_THRESHOLD = 0.5  # Mark cells below this as Unknown for scANVI

# ========== Dimensionality Reduction ==========
UMAP_N_NEIGHBORS = 50
UMAP_MIN_DIST = 0.3
UMAP_SPREAD = 1.5

# ========== Visualization ==========
DPI = 300
FIGURE_FORMAT = "pdf"

# ========== Checkpoint Configuration ==========
SAVE_CHECKPOINTS = True

# ========== Create Output Directories ==========
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

model_dir = output_dir / "models"
model_dir.mkdir(exist_ok=True)

scvi_model_dir = model_dir / "scvi_model"
scvi_model_dir.mkdir(exist_ok=True)

scanvi_model_dir = model_dir / "scanvi_model"
scanvi_model_dir.mkdir(exist_ok=True)

checkpoint_dir = output_dir / "checkpoints"
if SAVE_CHECKPOINTS:
    checkpoint_dir.mkdir(exist_ok=True)

print("\n✓ Configuration loaded")
print(f"  Input: {INPUT_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  CellTypist model: {CELLTYPIST_MODEL_PATH}")
print(f"\n✓ Directories created:")
print(f"  Main output: {output_dir}")
print(f"  Figures: {fig_dir}")
print(f"  Models: {model_dir}")
if SAVE_CHECKPOINTS:
    print(f"  Checkpoints: {checkpoint_dir}")

# %% [markdown]
# ---
# ## Step 1: Data Loading and Raw Counts Preparation

# %%
print("\n" + "="*70)
print("Step 1: Loading Data and Preparing Raw Counts")
print("="*70)

# Load data
adata = sc.read_h5ad(INPUT_H5AD)
print(f"\nLoaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Display basic info
print(f"\nData structure:")
print(f"  .X shape: {adata.X.shape}")
print(f"  .X type: {type(adata.X)}")
if issparse(adata.X):
    print(f"  .X sparse format: {adata.X.getformat()}")
print(f"  .obs columns ({len(adata.obs.columns)}): {list(adata.obs.columns[:10])}...")
print(f"  .var columns ({len(adata.var.columns)}): {list(adata.var.columns)}")
print(f"  .layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f"  .raw: {'Present' if adata.raw is not None else 'None'}")

# %%
# Critical: Check and prepare raw counts for scVI and CellTypist
print("\n" + "="*70)
print("Preparing Raw Counts")
print("="*70)

# Strategy: Copy layers['counts'] to .raw
if 'counts' not in adata.layers:
    raise KeyError(
        "layers['counts'] not found in adata.\n"
        "This pipeline requires raw counts stored in layers['counts']."
    )

print("\n✓ Found layers['counts']")

# Check if counts are truly raw (not log-transformed)
counts_sample = adata.layers['counts'][:1000, :1000] if adata.n_obs > 1000 else adata.layers['counts']
if issparse(counts_sample):
    counts_sample = counts_sample.toarray()

counts_max = np.max(counts_sample)
counts_mean = np.mean(counts_sample[counts_sample > 0]) if np.any(counts_sample > 0) else 0

print(f"\nValidating raw counts (sample statistics):")
print(f"  Max value: {counts_max:.2f}")
print(f"  Mean (non-zero): {counts_mean:.2f}")

if counts_max < 10 and counts_mean < 1:
    raise ValueError(
        f"layers['counts'] appears to be log-transformed!\n"
        f"Max={counts_max:.2f}, Mean={counts_mean:.2f}\n"
        f"scVI and downstream tools require raw counts (integers > 100 typical)."
    )

print("  ✓ Data appears to be raw counts (not log-transformed)")

# Create .raw from layers['counts']
print("\nCreating .raw from layers['counts']...")
adata.raw = sc.AnnData(
    X=adata.layers['counts'].copy(),
    var=adata.var.copy(),
    obs=adata.obs.copy()
)

print(f"✓ .raw created: {adata.raw.shape[0]:,} cells × {adata.raw.shape[1]:,} genes")
print(f"  .raw.X type: {type(adata.raw.X)}")
if issparse(adata.raw.X):
    print(f"  .raw.X format: {adata.raw.X.getformat()}")

# %%
# Verify batch key
print("\n" + "="*70)
print("Verifying Batch Information")
print("="*70)

if BATCH_KEY not in adata.obs.columns:
    raise KeyError(
        f"BATCH_KEY '{BATCH_KEY}' not found in adata.obs.\n"
        f"Available columns: {list(adata.obs.columns)}"
    )

n_batches = adata.obs[BATCH_KEY].nunique()
print(f"\nBatch key: '{BATCH_KEY}' ✓")
print(f"Number of batches: {n_batches}\n")

batch_counts = adata.obs[BATCH_KEY].value_counts()
for batch, count in batch_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  {batch}: {count:,} cells ({pct:.1f}%)")

# Check confounders
missing_confounders = [c for c in CONFOUNDER_KEYS if c not in adata.obs.columns]
if missing_confounders:
    print(f"\n⚠️  Warning: Confounders {missing_confounders} not found")
    print("   Will proceed without biological covariates")
    CONFOUNDER_KEYS = [c for c in CONFOUNDER_KEYS if c in adata.obs.columns]
else:
    print(f"\n✓ Confounders found: {CONFOUNDER_KEYS}")
    for conf in CONFOUNDER_KEYS:
        n_levels = adata.obs[conf].nunique()
        print(f"  {conf}: {n_levels} levels")

# %%
# Save checkpoint 1
if SAVE_CHECKPOINTS:
    checkpoint_file = checkpoint_dir / 'checkpoint_01_raw_prepared.h5ad'
    print(f"\nSaving checkpoint: {checkpoint_file}")
    adata.write_h5ad(checkpoint_file)
    print("✓ Checkpoint saved")

# %% [markdown]
# ---
# ## Step 2: scVI Batch Integration

# %%
print("\n" + "="*70)
print("Step 2: scVI Setup and Training")
print("="*70)

# Setup AnnData for scVI
print("\nSetting up AnnData for scVI...")

setup_kwargs = {
    'layer': 'counts',
    'batch_key': BATCH_KEY
}

if CONFOUNDER_KEYS:
    setup_kwargs['categorical_covariate_keys'] = CONFOUNDER_KEYS
    print(f"  Including confounders: {CONFOUNDER_KEYS}")

scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
print("✓ AnnData setup complete")

# Create scVI model
print(f"\nCreating scVI model...")
print(f"  Latent dimensions: {SCVI_N_LATENT}")
print(f"  Neural network layers: {SCVI_N_LAYERS}")
print(f"  Dropout rate: {SCVI_DROPOUT}")

vae = scvi.model.SCVI(
    adata,
    n_latent=SCVI_N_LATENT,
    n_layers=SCVI_N_LAYERS,
    dropout_rate=SCVI_DROPOUT,
    gene_likelihood='nb'
)
print("✓ scVI model created")

# %%
# Train scVI model
print("\n" + "="*70)
print("Training scVI Model")
print("="*70)
print(f"Max epochs: {SCVI_MAX_EPOCHS}")
print(f"Early stopping: {SCVI_EARLY_STOPPING}")
print(f"Learning rate: {SCVI_LEARNING_RATE}")
print(f"GPU available: {gpu_available}")
print("\nTraining started (estimated 10-20 minutes for 40k cells)...\n")

start_time = time.time()

vae.train(
    max_epochs=SCVI_MAX_EPOCHS,
    early_stopping=SCVI_EARLY_STOPPING,
    early_stopping_patience=15,
    plan_kwargs={'lr': SCVI_LEARNING_RATE}
)

training_time = time.time() - start_time

# Extract training metrics
train_elbo_raw = vae.history['elbo_train']
# ⭐ FIXED: More robust handling of history
if hasattr(train_elbo_raw, 'values'):
    train_elbo = list(train_elbo_raw.values)
else:
    train_elbo = list(train_elbo_raw)

n_training_epochs = len(train_elbo)
final_elbo_raw = train_elbo[-1]
if hasattr(final_elbo_raw, 'item'):
    final_elbo = float(final_elbo_raw.item())
else:
    final_elbo = float(final_elbo_raw)

print(f"\n{'='*70}")
print(f"✓ Training completed")
print(f"  Time: {training_time:.1f}s ({training_time/60:.1f} min)")
print(f"  Epochs: {n_training_epochs}")
print(f"  Final ELBO: {final_elbo:.2f}")
print("="*70)

# %%
# Plot training history
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Full history
axes[0].plot(train_elbo, linewidth=2, color='#1f77b4')
axes[0].set_xlabel('Epoch', fontsize=12)
axes[0].set_ylabel('ELBO', fontsize=12)
axes[0].set_title('scVI Training History', fontsize=14, fontweight='bold')
axes[0].grid(alpha=0.3)

# Last 50 epochs
n_last = min(50, len(train_elbo))
axes[1].plot(train_elbo[-n_last:], linewidth=2, color='#ff7f0e')
axes[1].set_xlabel(f'Epoch (last {n_last})', fontsize=12)
axes[1].set_ylabel('ELBO', fontsize=12)
axes[1].set_title('Convergence Detail', fontsize=14, fontweight='bold')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / f'scvi_training_history.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Training history plot saved")

# %%
# Extract latent representation and save model
print("\nExtracting latent representation...")
latent = vae.get_latent_representation()
adata.obsm['X_scvi'] = latent
print(f"✓ Latent representation extracted: {latent.shape}")

# Save scVI model
print(f"\nSaving scVI model to: {scvi_model_dir}")
vae.save(scvi_model_dir, overwrite=True)
print("✓ scVI model saved")

# Compute neighbors and UMAP from scVI latent space
print("\nComputing neighbors from scVI latent space...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS)
print(f"✓ Neighbors computed (k={UMAP_N_NEIGHBORS})")

print("\nComputing UMAP from scVI...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=RANDOM_SEED)
print(f"✓ UMAP computed (min_dist={UMAP_MIN_DIST}, spread={UMAP_SPREAD})")

# ⭐ CRITICAL FIX: Save scVI UMAP immediately to prevent overwriting
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print("\n✓ scVI UMAP saved as 'X_umap_scvi' (prevents overwriting)")

# %%
# Visualize scVI results
print("\nGenerating scVI batch mixing visualization...")

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0],
           title='scVI Batch Integration', show=False,
           legend_loc='right margin', frameon=False)

# Color by original annotation if available
if 'Annotation' in adata.obs.columns:
    sc.pl.umap(adata, color='Annotation', ax=axes[1],
               title='Original BBKNN Annotation', show=False,
               legend_loc='on data', legend_fontsize=7, frameon=False)
else:
    axes[1].text(0.5, 0.5, 'No annotation available',
                ha='center', va='center', transform=axes[1].transAxes,
                fontsize=14)
    axes[1].axis('off')

plt.tight_layout()
plt.savefig(fig_dir / f'umap_scvi_batch_mixing.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ scVI visualization saved")

# %%
# Save checkpoint 2
if SAVE_CHECKPOINTS:
    checkpoint_file = checkpoint_dir / 'checkpoint_02_scvi_trained.h5ad'
    print(f"\nSaving checkpoint: {checkpoint_file}")
    adata.write_h5ad(checkpoint_file)
    print("✓ Checkpoint saved")

# %% [markdown]
# ---
# ## Step 3: CellTypist Automatic Annotation
# 
# ### 📌 CellTypist Requirements (v2.0 Updated)
# 
# 1. **Log-normalized counts** (not raw counts!)
# 2. **No missing values** in `.X`
# 3. **Gene names** must match model:
#    - CellTypist human models (including Immune_All_Low) expect **HGNC gene symbols**
#    - We only need to convert if data uses ENSEMBL IDs
# 
# ### 🔄 Simplified Gene Conversion (v2.0)
# 
# ```
# 1. Check if adata.var_names are ENSEMBL IDs (start with 'ENSG')
# 2. If yes → Convert ENSEMBL → SYMBOL using mygene
# 3. If no → Assume already gene symbols, proceed
# 4. Run var_names_make_unique() to handle duplicates
# 5. Run CellTypist annotation
# ```

# %%
print("\n" + "="*70)
print("Step 3: CellTypist Annotation Preparation")
print("="*70)

# Load CellTypist model
print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL_PATH}")

try:
    celltypist_model = models.Model.load(model=CELLTYPIST_MODEL_PATH)
    print("✓ Model loaded successfully")
    
    print(f"\nModel information:")
    print(f"  Cell types: {len(celltypist_model.cell_types)}")
    print(f"  Example cell types: {list(celltypist_model.cell_types[:5])}")
    
except Exception as e:
    print(f"\n✗ ERROR loading model: {e}")
    print("\nPlease check:")
    print(f"  1. Model path: {CELLTYPIST_MODEL_PATH}")
    print("  2. Model file exists and is valid")
    print("  3. CellTypist is properly installed")
    raise

# %%
# ⭐ FIXED: Simplified gene format detection
print("\n" + "="*70)
print("Gene Name Format Detection (Simplified)")
print("="*70)

sample_adata_gene = str(adata.var_names[0])
print(f"\nData gene format:")
print(f"  First gene: {sample_adata_gene}")

# CellTypist human models expect HGNC gene symbols
# Only need to convert if data uses ENSEMBL IDs
if sample_adata_gene.startswith('ENSG'):
    print(f"  Format: ENSEMBL")
    print("\n⚠️  Need conversion: ENSEMBL → SYMBOL")
    print("   CellTypist models expect gene symbols")
    need_conversion = True
else:
    print(f"  Format: Gene Symbol (assumed)")
    print("\n✓ No conversion needed - proceeding with current gene names")
    need_conversion = False

# %%
# Gene name conversion (if needed)
if need_conversion:
    print("\n" + "="*70)
    print("Gene Name Conversion: ENSEMBL → SYMBOL")
    print("="*70)
    
    print("\nQuerying mygene database (may take 1-2 minutes)...")
    
    mg = mygene.MyGeneInfo()
    ensembl_ids = adata.var_names.tolist()
    
    query_results = mg.querymany(
        ensembl_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species='human'
    )
    
    ensembl_to_symbol = {}
    for result in query_results:
        if 'symbol' in result:
            symbol = result['symbol']
            ensembl_to_symbol[result['query']] = symbol
    
    # Update var_names
    original_n = len(adata.var_names)
    adata.var_names = [ensembl_to_symbol.get(ensembl_id, ensembl_id) for ensembl_id in adata.var_names]
    
    converted_n = sum([1 for gene in adata.var_names if not gene.startswith('ENSG')])
    
    print(f"\nConversion results:")
    print(f"  Original genes: {original_n}")
    print(f"  Converted: {converted_n} ({converted_n/original_n*100:.1f}%)")
    print(f"  Unmapped: {original_n - converted_n}")
    
    # ⭐ FIXED: Handle duplicate gene names from many-to-one mapping
    print("\nHandling duplicate gene names...")
    n_duplicates_before = adata.var_names.duplicated().sum()
    if n_duplicates_before > 0:
        print(f"  Found {n_duplicates_before} duplicate gene names")
        adata.var_names_make_unique()
        print(f"  ✓ Duplicates resolved (appended suffixes)")
    else:
        print("  No duplicates found")
    
    print("\n✓ Gene name conversion completed")
    print(f"  Updated var_names format: SYMBOL")
    print(f"  Example genes: {list(adata.var_names[:5])}")

# %%
# ⭐ FIXED: Explicit raw counts preparation (no guessing)
print("\n" + "="*70)
print("Preparing Data for CellTypist (Explicit Method)")
print("="*70)

# Create a copy and explicitly use raw counts
print("\nCreating CellTypist input from raw counts...")
adata_celltypist = adata.copy()

# Explicitly use layers['counts'] as starting point
print("  Using layers['counts'] as input")
adata_celltypist.X = adata.layers['counts'].copy()

# Apply normalization and log transformation
print("  Applying normalize_total (target_sum=1e4)")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)

print("  Applying log1p transformation")
sc.pp.log1p(adata_celltypist)

# Check for missing values
print("\nChecking for missing values...")
if issparse(adata_celltypist.X):
    has_nan = np.isnan(adata_celltypist.X.data).any()
else:
    has_nan = np.isnan(adata_celltypist.X).any()

if has_nan:
    print("  ⚠️  Found NaN values, filling with 0...")
    if issparse(adata_celltypist.X):
        adata_celltypist.X.data = np.nan_to_num(adata_celltypist.X.data, nan=0.0)
    else:
        adata_celltypist.X = np.nan_to_num(adata_celltypist.X, nan=0.0)
    print("  ✓ NaN values handled")
else:
    print("  ✓ No missing values found")

# Remove .raw to avoid interference
if adata_celltypist.raw is not None:
    print("\n  Removing .raw to avoid CellTypist issues...")
    del adata_celltypist.raw

# Verify final state
x_sample = adata_celltypist.X[:100, :100]
if issparse(x_sample):
    x_sample = x_sample.toarray()
x_max = np.max(x_sample)
x_mean = np.mean(x_sample[x_sample > 0])

print("\nFinal .X statistics:")
print(f"  Max value: {x_max:.2f}")
print(f"  Mean (non-zero): {x_mean:.2f}")
print(f"  Expected: max ~8-12, mean ~1-3 for log-normalized data")

print("\n✓ Data ready for CellTypist annotation")

# %%
# Run CellTypist annotation
print("\n" + "="*70)
print("Running CellTypist Annotation")
print("="*70)

print(f"\nModel: {CELLTYPIST_MODEL_PATH}")
print(f"Majority voting: {CELLTYPIST_MAJORITY_VOTING}")
print(f"Mode: {CELLTYPIST_MODE}")
print("\nAnnotation started (estimated 2-5 minutes)...\n")

start_time = time.time()

predictions = celltypist.annotate(
    adata_celltypist,
    model=CELLTYPIST_MODEL_PATH,
    majority_voting=CELLTYPIST_MAJORITY_VOTING,
    mode=CELLTYPIST_MODE
)

annotation_time = time.time() - start_time

print(f"\n✓ Annotation completed in {annotation_time:.1f}s ({annotation_time/60:.1f} min)")

# Convert predictions to AnnData
adata_celltypist = predictions.to_adata()

print("\nAnnotation results added to .obs:")
celltypist_cols = [col for col in adata_celltypist.obs.columns 
                   if col in ['predicted_labels', 'majority_voting', 'over_clustering', 'conf_score']]
for col in celltypist_cols:
    print(f"  - {col}")

# %%
# Display annotation results
print("\n" + "="*70)
print("CellTypist Annotation Summary")
print("="*70)

# Get total cell count for consistent percentage calculations
n_total_cells = adata_celltypist.n_obs

# Predicted labels distribution
print("\n[ Predicted Labels Distribution ]")
label_counts = adata_celltypist.obs['predicted_labels'].value_counts()
print(f"\nTotal cell types: {len(label_counts)}")
print(f"\nTop 10 cell types:")
for label, count in label_counts.head(10).items():
    pct = count / n_total_cells * 100
    print(f"  {label}: {count:,} cells ({pct:.1f}%)")

if len(label_counts) > 10:
    print(f"  ... and {len(label_counts) - 10} more cell types")

# Majority voting comparison (if enabled)
if CELLTYPIST_MAJORITY_VOTING and 'majority_voting' in adata_celltypist.obs.columns:
    print("\n[ Majority Voting vs Predicted Labels ]")
    
    # Convert to string to avoid categorical comparison issues
    n_changed = (adata_celltypist.obs['predicted_labels'].astype(str) != 
                 adata_celltypist.obs['majority_voting'].astype(str)).sum()
    pct_changed = n_changed / n_total_cells * 100
    
    print(f"  Cells changed by majority voting: {n_changed:,} ({pct_changed:.1f}%)")

# Confidence scores (if available)
if 'conf_score' in adata_celltypist.obs.columns:
    print("\n[ Confidence Scores ]")
    conf_scores = adata_celltypist.obs['conf_score']
    print(f"  Mean: {conf_scores.mean():.3f}")
    print(f"  Median: {conf_scores.median():.3f}")
    print(f"  Min: {conf_scores.min():.3f}")
    print(f"  Max: {conf_scores.max():.3f}")
    
    # Low confidence cells
    n_low_conf = (conf_scores < LOW_CONF_THRESHOLD).sum()
    pct_low_conf = n_low_conf / n_total_cells * 100
    print(f"\n  Cells with confidence < {LOW_CONF_THRESHOLD}: {n_low_conf:,} ({pct_low_conf:.1f}%)")

# %%
# Transfer CellTypist annotations to main adata
print("\n" + "="*70)
print("Transferring Annotations to Main Dataset")
print("="*70)

# Transfer annotation columns
for col in celltypist_cols:
    if col in adata_celltypist.obs.columns:
        adata.obs[col] = adata_celltypist.obs[col].values
        print(f"  ✓ Transferred: {col}")

print("\n✓ CellTypist annotations added to main dataset")

# Clean up temporary object
del adata_celltypist
print("  Removed temporary CellTypist object")

# %%
# Visualize CellTypist annotations on scVI UMAP
print("\nGenerating CellTypist annotation visualization...")

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# Predicted labels
sc.pl.umap(adata, color='predicted_labels', ax=axes[0],
           title='CellTypist Predicted Labels', show=False,
           legend_loc='right margin', frameon=False)

# Majority voting (if available)
if 'majority_voting' in adata.obs.columns:
    sc.pl.umap(adata, color='majority_voting', ax=axes[1],
               title='CellTypist Majority Voting', show=False,
               legend_loc='right margin', frameon=False)
else:
    axes[1].text(0.5, 0.5, 'Majority voting not performed',
                ha='center', va='center', transform=axes[1].transAxes,
                fontsize=14)
    axes[1].axis('off')

plt.tight_layout()
plt.savefig(fig_dir / f'umap_celltypist_annotations.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ CellTypist visualization saved")

# %%
# Save checkpoint 3
if SAVE_CHECKPOINTS:
    checkpoint_file = checkpoint_dir / 'checkpoint_03_celltypist_annotated.h5ad'
    print(f"\nSaving checkpoint: {checkpoint_file}")
    adata.write_h5ad(checkpoint_file)
    print("✓ Checkpoint saved")

# %% [markdown]
# ---
# ## Step 4: scANVI Semi-supervised Refinement
# 
# ### 📌 scANVI Strategy
# 
# scANVI is a **semi-supervised extension of scVI** that:
# 1. Uses existing cell type labels (from CellTypist) as references
# 2. Refines annotations for low-confidence cells
# 3. Predicts labels for unlabeled/uncertain cells
# 4. Generates a new latent space (`X_scanvi`) optimized for classification
# 
# ### 🎯 Label Processing
# 
# We'll use `majority_voting` (if available) or `predicted_labels` as the reference.
# Low-confidence cells (< 0.5) will be marked as "Unknown" for scANVI to refine.

# %%
print("\n" + "="*70)
print("Step 4: scANVI Semi-supervised Training")
print("="*70)

# Choose labels for scANVI
if 'majority_voting' in adata.obs.columns:
    labels_key = 'majority_voting'
    print(f"\nUsing 'majority_voting' as reference labels")
else:
    labels_key = 'predicted_labels'
    print(f"\nUsing 'predicted_labels' as reference labels")

# Create labels column for scANVI (convert to string to avoid categorical issues)
adata.obs['scanvi_labels'] = adata.obs[labels_key].astype(str).copy()

# Mark low-confidence cells as unlabeled
if 'conf_score' in adata.obs.columns:
    low_conf_mask = adata.obs['conf_score'] < LOW_CONF_THRESHOLD
    n_low_conf = low_conf_mask.sum()
    
    print(f"\nMarking low-confidence cells as '{UNLABELED_CATEGORY}':")
    print(f"  Threshold: {LOW_CONF_THRESHOLD}")
    print(f"  Cells to refine: {n_low_conf:,} ({n_low_conf/adata.n_obs*100:.1f}%)")
    
    # Now safe to assign since we converted to string above
    adata.obs.loc[low_conf_mask, 'scanvi_labels'] = UNLABELED_CATEGORY
else:
    print("\n⚠️  No confidence scores available, using all labels as-is")

# Display label distribution
label_counts = adata.obs['scanvi_labels'].value_counts()
print(f"\nLabel distribution for scANVI:")
print(f"  Total categories: {len(label_counts)}")
if UNLABELED_CATEGORY in label_counts.index:
    print(f"  {UNLABELED_CATEGORY}: {label_counts[UNLABELED_CATEGORY]:,} cells")
print(f"  Labeled cell types: {len(label_counts) - (1 if UNLABELED_CATEGORY in label_counts.index else 0)}")

# %%
# Initialize scANVI from trained scVI model
print("\n" + "="*70)
print("Initializing scANVI Model")
print("="*70)

print("\nConverting scVI model to scANVI...")

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key='scanvi_labels',
    unlabeled_category=UNLABELED_CATEGORY
)

print("✓ scANVI model initialized")
print(f"  Inherits scVI latent space: {SCVI_N_LATENT} dimensions")
print(f"  Reference labels: scanvi_labels")
print(f"  Unlabeled category: {UNLABELED_CATEGORY}")

# %%
# Train scANVI model
print("\n" + "="*70)
print("Training scANVI Model")
print("="*70)
print(f"Max epochs: {SCANVI_MAX_EPOCHS}")
print(f"Early stopping: {SCANVI_EARLY_STOPPING}")
print("\nTraining started (estimated 10-15 minutes)...\n")

start_time = time.time()

lvae.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    early_stopping=SCANVI_EARLY_STOPPING,
    early_stopping_patience=15
)

scanvi_training_time = time.time() - start_time

# Extract training metrics
scanvi_train_elbo_raw = lvae.history['elbo_train']
if hasattr(scanvi_train_elbo_raw, 'values'):
    scanvi_train_elbo = list(scanvi_train_elbo_raw.values)
else:
    scanvi_train_elbo = list(scanvi_train_elbo_raw)

n_scanvi_epochs = len(scanvi_train_elbo)
final_scanvi_elbo_raw = scanvi_train_elbo[-1]
if hasattr(final_scanvi_elbo_raw, 'item'):
    final_scanvi_elbo = float(final_scanvi_elbo_raw.item())
else:
    final_scanvi_elbo = float(final_scanvi_elbo_raw)

print(f"\n{'='*70}")
print(f"✓ scANVI training completed")
print(f"  Time: {scanvi_training_time:.1f}s ({scanvi_training_time/60:.1f} min)")
print(f"  Epochs: {n_scanvi_epochs}")
print(f"  Final ELBO: {final_scanvi_elbo:.2f}")
print("="*70)

# %%
# Plot scANVI training history
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Full history
axes[0].plot(scanvi_train_elbo, linewidth=2, color='#2ca02c')
axes[0].set_xlabel('Epoch', fontsize=12)
axes[0].set_ylabel('ELBO', fontsize=12)
axes[0].set_title('scANVI Training History', fontsize=14, fontweight='bold')
axes[0].grid(alpha=0.3)

# Last 50 epochs
n_last = min(50, len(scanvi_train_elbo))
axes[1].plot(scanvi_train_elbo[-n_last:], linewidth=2, color='#d62728')
axes[1].set_xlabel(f'Epoch (last {n_last})', fontsize=12)
axes[1].set_ylabel('ELBO', fontsize=12)
axes[1].set_title('Convergence Detail', fontsize=14, fontweight='bold')
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(fig_dir / f'scanvi_training_history.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ scANVI training history plot saved")

# %%
# Extract scANVI results
print("\n" + "="*70)
print("Extracting scANVI Results")
print("="*70)

# Get latent representation
print("\nExtracting latent representation...")
latent_scanvi = lvae.get_latent_representation()
adata.obsm['X_scanvi'] = latent_scanvi
print(f"✓ Latent representation extracted: {latent_scanvi.shape}")

# Get refined predictions
print("\nGenerating refined predictions...")
adata.obs['scanvi_predictions'] = lvae.predict()
print("✓ Predictions generated")

# Get prediction probabilities
print("\nCalculating prediction probabilities...")
predictions_df = lvae.predict(soft=True)
adata.obs['scanvi_confidence'] = predictions_df.max(axis=1).values
print("✓ Confidence scores calculated")

# ⭐ NEW in v2.0: Save full probability matrix for posterior analysis
print("\nSaving probability matrix...")
adata.obsm['scanvi_probabilities'] = predictions_df.values
adata.uns['scanvi_celltype_order'] = predictions_df.columns.tolist()
print(f"✓ Probabilities saved: {predictions_df.shape}")
print(f"  Access via: adata.obsm['scanvi_probabilities']")
print(f"  Cell type order: adata.uns['scanvi_celltype_order']")

# Save scANVI model
print(f"\nSaving scANVI model to: {scanvi_model_dir}")
lvae.save(scanvi_model_dir, overwrite=True)
print("✓ scANVI model saved")

# %%
# Compute UMAP from scANVI latent space
print("\nComputing UMAP from scANVI latent space...")

sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=UMAP_N_NEIGHBORS, key_added='scanvi')
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, neighbors_key='scanvi')

# ⭐ CRITICAL FIX: Save scANVI UMAP separately (don't overwrite scVI UMAP)
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()

print(f"✓ UMAP computed and saved as 'X_umap_scanvi'")
print("\n✓ Now we have TWO separate UMAPs:")
print("  - X_umap_scvi: From scVI (batch-corrected)")
print("  - X_umap_scanvi: From scANVI (classification-optimized)")

# %%
# scANVI prediction summary
print("\n" + "="*70)
print("scANVI Prediction Summary")
print("="*70)

# Compare original CellTypist labels with scANVI predictions
print("\n[ Label Comparison ]")

# For cells that were labeled by CellTypist
labeled_mask = adata.obs['scanvi_labels'] != UNLABELED_CATEGORY
n_labeled = labeled_mask.sum()

# For cells that scANVI had to predict
unlabeled_mask = ~labeled_mask
n_unlabeled = unlabeled_mask.sum()

print(f"\nCells with CellTypist labels: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"Cells predicted by scANVI: {n_unlabeled:,} ({n_unlabeled/adata.n_obs*100:.1f}%)")

# Agreement between CellTypist and scANVI (for labeled cells)
if n_labeled > 0:
    agreement = (adata.obs.loc[labeled_mask, 'scanvi_labels'] == 
                 adata.obs.loc[labeled_mask, 'scanvi_predictions']).sum()
    print(f"\nAgreement (labeled cells): {agreement:,}/{n_labeled:,} ({agreement/n_labeled*100:.1f}%)")

# Confidence scores
print("\n[ scANVI Confidence Scores ]")
print(f"  Mean: {adata.obs['scanvi_confidence'].mean():.3f}")
print(f"  Median: {adata.obs['scanvi_confidence'].median():.3f}")
print(f"  Min: {adata.obs['scanvi_confidence'].min():.3f}")
print(f"  Max: {adata.obs['scanvi_confidence'].max():.3f}")

# High confidence predictions
high_conf_threshold = 0.8
n_high_conf = (adata.obs['scanvi_confidence'] > high_conf_threshold).sum()
print(f"\n  High confidence (> {high_conf_threshold}): {n_high_conf:,} ({n_high_conf/adata.n_obs*100:.1f}%)")

# Final cell type distribution
print("\n[ Final Cell Type Distribution ]")
final_counts = adata.obs['scanvi_predictions'].value_counts()
print(f"\nTotal cell types: {len(final_counts)}")
print(f"\nTop 10 cell types:")
for celltype, count in final_counts.head(10).items():
    pct = count / adata.n_obs * 100
    print(f"  {celltype}: {count:,} cells ({pct:.1f}%)")

# %%
# Visualize scANVI results (using scANVI UMAP)
print("\nGenerating scANVI visualization...")

# Temporarily set UMAP to scANVI version for plotting
adata.obsm['X_umap_temp'] = adata.obsm['X_umap'].copy()
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()

fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# scANVI predictions
sc.pl.umap(adata, color='scanvi_predictions', ax=axes[0, 0],
           title='scANVI Refined Predictions', show=False,
           legend_loc='right margin', frameon=False)

# Confidence scores
sc.pl.umap(adata, color='scanvi_confidence', ax=axes[0, 1],
           title='scANVI Confidence Scores', show=False,
           cmap='viridis', frameon=False)

# Batch mixing
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0],
           title='Batch Distribution (scANVI UMAP)', show=False,
           legend_loc='right margin', frameon=False)

# Original CellTypist labels
sc.pl.umap(adata, color=labels_key, ax=axes[1, 1],
           title=f'Original {labels_key.replace("_", " ").title()}', show=False,
           legend_loc='right margin', frameon=False)

plt.tight_layout()
plt.savefig(fig_dir / f'umap_scanvi_results.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

# Restore original UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_temp'].copy()
del adata.obsm['X_umap_temp']

print("✓ scANVI visualization saved")

# %%
# Save checkpoint 4
if SAVE_CHECKPOINTS:
    checkpoint_file = checkpoint_dir / 'checkpoint_04_scanvi_refined.h5ad'
    print(f"\nSaving checkpoint: {checkpoint_file}")
    adata.write_h5ad(checkpoint_file)
    print("✓ Checkpoint saved")

# %% [markdown]
# ---
# ## Step 5: Comprehensive Comparison and Visualization
# 
# ### ⭐ v2.0 Fixed: Correct scVI vs scANVI Comparison
# 
# This step now correctly compares the TWO separate UMAPs:
# - **Row 1:** scVI UMAP (`X_umap_scvi`) - batch-corrected embedding
# - **Row 2:** scANVI UMAP (`X_umap_scanvi`) - classification-optimized embedding

# %%
print("\n" + "="*70)
print("Step 5: Comprehensive Comparison (FIXED in v2.0)")
print("="*70)

# ⭐ CRITICAL FIX: Side-by-side comparison using SEPARATE UMAPs
print("\nGenerating side-by-side comparison...")
print("Row 1: scVI UMAP (X_umap_scvi)")
print("Row 2: scANVI UMAP (X_umap_scanvi)")

fig, axes = plt.subplots(2, 3, figsize=(20, 12))

# Row 1: scVI UMAP
adata.obsm['X_umap_temp'] = adata.obsm['X_umap'].copy()
adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0],
           title='scVI: Batch', show=False, legend_loc='right margin', frameon=False)

sc.pl.umap(adata, color='predicted_labels', ax=axes[0, 1],
           title='scVI: CellTypist Labels', show=False, legend_loc='right margin', frameon=False)

if 'conf_score' in adata.obs.columns:
    sc.pl.umap(adata, color='conf_score', ax=axes[0, 2],
               title='scVI: CellTypist Confidence', show=False, cmap='viridis', frameon=False)
else:
    axes[0, 2].axis('off')

# Row 2: scANVI UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0],
           title='scANVI: Batch', show=False, legend_loc='right margin', frameon=False)

sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1, 1],
           title='scANVI: Refined Predictions', show=False, legend_loc='right margin', frameon=False)

sc.pl.umap(adata, color='scanvi_confidence', ax=axes[1, 2],
           title='scANVI: Confidence', show=False, cmap='viridis', frameon=False)

# Restore original UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_temp'].copy()
del adata.obsm['X_umap_temp']

plt.tight_layout()
plt.savefig(fig_dir / f'comparison_scvi_vs_scanvi.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Comparison plot saved (using correct UMAPs)")

# %%
# Cell type composition by batch
print("\nGenerating cell type composition analysis...")

# Create composition DataFrame
composition_df = pd.crosstab(
    adata.obs[BATCH_KEY],
    adata.obs['scanvi_predictions'],
    normalize='index'
) * 100

# Plot top 10 cell types
top_celltypes = adata.obs['scanvi_predictions'].value_counts().head(10).index
composition_top = composition_df[top_celltypes]

fig, ax = plt.subplots(figsize=(14, 8))
composition_top.T.plot(kind='bar', stacked=False, ax=ax, width=0.8)

ax.set_xlabel('Cell Type', fontsize=12)
ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Cell Type Composition by Batch (Top 10)', fontsize=14, fontweight='bold')
ax.legend(title='Batch', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(fig_dir / f'composition_by_batch.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Composition plot saved")

# %%
# Confidence score comparison: CellTypist vs scANVI
if 'conf_score' in adata.obs.columns:
    print("\nGenerating confidence score comparison...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    axes[0].hist(adata.obs['conf_score'], bins=50, alpha=0.5, label='CellTypist', color='blue')
    axes[0].hist(adata.obs['scanvi_confidence'], bins=50, alpha=0.5, label='scANVI', color='red')
    axes[0].set_xlabel('Confidence Score', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('Confidence Score Distribution', fontsize=14, fontweight='bold')
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # Scatter plot
    axes[1].scatter(adata.obs['conf_score'], adata.obs['scanvi_confidence'], 
                   alpha=0.3, s=1, c='black')
    axes[1].plot([0, 1], [0, 1], 'r--', linewidth=2, label='Identity line')
    axes[1].set_xlabel('CellTypist Confidence', fontsize=12)
    axes[1].set_ylabel('scANVI Confidence', fontsize=12)
    axes[1].set_title('Confidence Correlation', fontsize=14, fontweight='bold')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'confidence_comparison.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Confidence comparison plot saved")

# %% [markdown]
# ---
# ## Step 6: Final Data Export and Structure Documentation

# %%
print("\n" + "="*70)
print("Step 6: Saving Final Results")
print("="*70)

# Save final annotated dataset
output_file = output_dir / "adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"✓ Final dataset saved ({file_size:.2f} GB)")

# %%
# Document complete data structure
print("\n" + "="*70)
print("COMPLETE DATA STRUCTURE (v2.0)")
print("="*70)

print(f"\nAnnData object: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# .X
print("\n[ .X - Main Expression Matrix ]")
print(f"  Type: {type(adata.X)}")
if issparse(adata.X):
    print(f"  Format: {adata.X.getformat()}")
print(f"  Content: Processed expression (likely log-normalized from BBKNN)")

# .layers
print("\n[ .layers - Alternative Representations ]")
if adata.layers:
    for layer_name in adata.layers.keys():
        print(f"  {layer_name}: {adata.layers[layer_name].shape}")
        if layer_name == 'counts':
            print(f"    → Raw UMI counts (integer values)")
else:
    print("  None")

# .raw
print("\n[ .raw - Raw Counts Storage ]")
if adata.raw is not None:
    print(f"  Shape: {adata.raw.shape[0]:,} cells × {adata.raw.shape[1]:,} genes")
    print(f"  Content: Copy of layers['counts'] for marker gene visualization")
    print(f"  Type: {type(adata.raw.X)}")
else:
    print("  None")

# .obs (cell metadata)
print("\n[ .obs - Cell Metadata ]")
print(f"  Columns: {len(adata.obs.columns)}")
print("\n  Key columns (v2.0):")

key_columns = {
    BATCH_KEY: "Batch identifier",
    'predicted_labels': "CellTypist initial prediction",
    'majority_voting': "CellTypist majority voting result",
    'conf_score': "CellTypist confidence score",
    'scanvi_labels': "Labels used for scANVI training",
    'scanvi_predictions': "scANVI refined predictions (FINAL ANNOTATION)",
    'scanvi_confidence': "scANVI prediction confidence",
    'over_clustering': "CellTypist over-clustering result"
}

for col, desc in key_columns.items():
    if col in adata.obs.columns:
        dtype = adata.obs[col].dtype
        print(f"    {col} ({dtype}): {desc}")

# .var (gene metadata)
print("\n[ .var - Gene Metadata ]")
print(f"  Columns: {list(adata.var.columns)}")
print(f"  Gene names: {adata.var_names[0]} ... {adata.var_names[-1]}")

# .obsm (dimensional reductions)
print("\n[ .obsm - Dimensional Reductions (v2.0 FIXED) ]")
for key in adata.obsm.keys():
    shape = adata.obsm[key].shape
    print(f"  {key}: {shape}")
    if key == 'X_scvi':
        print(f"    → scVI latent space ({SCVI_N_LATENT} dims, batch-corrected)")
    elif key == 'X_scanvi':
        print(f"    → scANVI latent space ({SCVI_N_LATENT} dims, classification-optimized)")
    elif key == 'X_umap_scvi':  # ⭐ FIXED
        print(f"    → UMAP from scVI (2D, batch-corrected visualization)")
    elif key == 'X_umap_scanvi':  # ⭐ FIXED
        print(f"    → UMAP from scANVI (2D, classification visualization)")
    elif key == 'scanvi_probabilities':  # ⭐ NEW in v2.0
        print(f"    → scANVI posterior probabilities (for all cell types)")
    elif key == 'X_pca':
        print(f"    → PCA from original BBKNN analysis")

# .obsp (pairwise relationships)
print("\n[ .obsp - Cell-Cell Graphs ]")
if adata.obsp:
    for key in adata.obsp.keys():
        shape = adata.obsp[key].shape
        print(f"  {key}: {shape}")
else:
    print("  None")

# .uns (unstructured metadata)
print("\n[ .uns - Unstructured Metadata ]")
if adata.uns:
    print(f"  Keys: {list(adata.uns.keys())}")
    if 'scanvi_celltype_order' in adata.uns:  # ⭐ NEW in v2.0
        print(f"\n  scanvi_celltype_order: Order of cell types in scanvi_probabilities")
else:
    print("  None")

# %%
# Generate analysis summary report
print("\n" + "="*70)
print("ANALYSIS SUMMARY REPORT (v2.0)")
print("="*70)

summary = []
summary.append("T/NK Cell Analysis: scVI → CellTypist → scANVI Pipeline v2.0")
summary.append("="*70)
summary.append(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary.append("")

summary.append("[ Dataset Information ]")
summary.append(f"  Cells: {adata.shape[0]:,}")
summary.append(f"  Genes: {adata.shape[1]:,}")
summary.append(f"  Batches: {n_batches}")
summary.append("")

summary.append("[ Step 1: scVI Integration ]")
summary.append(f"  Latent dimensions: {SCVI_N_LATENT}")
summary.append(f"  Training epochs: {n_training_epochs}")
summary.append(f"  Training time: {training_time/60:.1f} min")
summary.append(f"  Final ELBO: {final_elbo:.2f}")
summary.append("")

summary.append("[ Step 2: CellTypist Annotation ]")
summary.append(f"  Model: {CELLTYPIST_MODEL_PATH}")
summary.append(f"  Cell types detected: {len(label_counts)}")
summary.append(f"  Annotation time: {annotation_time/60:.1f} min")
if 'conf_score' in adata.obs.columns:
    summary.append(f"  Mean confidence: {adata.obs['conf_score'].mean():.3f}")
summary.append("")

summary.append("[ Step 3: scANVI Refinement ]")
summary.append(f"  Training epochs: {n_scanvi_epochs}")
summary.append(f"  Training time: {scanvi_training_time/60:.1f} min")
summary.append(f"  Final ELBO: {final_scanvi_elbo:.2f}")
summary.append(f"  Mean confidence: {adata.obs['scanvi_confidence'].mean():.3f}")
summary.append("")

summary.append("[ Final Cell Type Distribution (Top 10) ]")
final_counts = adata.obs['scanvi_predictions'].value_counts().head(10)
for celltype, count in final_counts.items():
    pct = count / adata.n_obs * 100
    summary.append(f"  {celltype}: {count:,} cells ({pct:.1f}%)")
summary.append("")

summary.append("[ Output Files ]")
summary.append(f"  Main dataset: {output_file}")
summary.append(f"  scVI model: {scvi_model_dir}")
summary.append(f"  scANVI model: {scanvi_model_dir}")
summary.append(f"  Figures: {fig_dir}/*.{FIGURE_FORMAT}")
if SAVE_CHECKPOINTS:
    summary.append(f"  Checkpoints: {checkpoint_dir}/*.h5ad")
summary.append("")

summary.append("[ Key Annotations in .obs ]")
summary.append("  scanvi_predictions: FINAL refined cell type annotation")
summary.append("  scanvi_confidence: Confidence score for each prediction")
summary.append("  predicted_labels: Original CellTypist prediction")
summary.append("  conf_score: CellTypist confidence score")
summary.append("")

summary.append("[ Embeddings in .obsm (v2.0 FIXED) ]")
summary.append("  X_scvi: scVI latent space (batch-corrected)")
summary.append("  X_scanvi: scANVI latent space (classification-optimized)")
summary.append("  X_umap_scvi: UMAP from scVI (batch-corrected visualization)")
summary.append("  X_umap_scanvi: UMAP from scANVI (classification visualization)")
summary.append("  scanvi_probabilities: Full posterior probability matrix (NEW)")
summary.append("")

summary.append("[ v2.0 Changes ]")
summary.append("  ✓ Fixed UMAP overwriting issue (scVI and scANVI now separate)")
summary.append("  ✓ Simplified CellTypist gene format detection")
summary.append("  ✓ Explicit raw counts preparation for CellTypist")
summary.append("  ✓ Added scvi.settings.seed for reproducibility")
summary.append("  ✓ Added var_names_make_unique() after gene conversion")
summary.append("  ✓ Saved scANVI probability matrix for posterior analysis")
summary.append("")

summary.append("="*70)
summary.append("Analysis completed successfully")
summary.append("="*70)

summary_text = '\n'.join(summary)

# Print and save summary
print("\n" + summary_text)

summary_file = output_dir / "analysis_summary_v2.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

print(f"\n✓ Summary saved to: {summary_file}")

# %% [markdown]
# ---
# ## ✅ Analysis Complete (v2.0)
# 
# ### 🆕 What's New in v2.0
# 
# **Critical Fixes:**
# 1. **UMAP Separation:** scVI and scANVI UMAPs are now correctly saved as separate embeddings (`X_umap_scvi` and `X_umap_scanvi`)
# 2. **Gene Conversion:** Simplified logic, only converts ENSEMBL→SYMBOL when needed
# 3. **CellTypist Input:** Now explicitly uses `layers['counts']`, no guessing
# 4. **Reproducibility:** Added `scvi.settings.seed` for full reproducibility
# 5. **Gene Names:** Added `var_names_make_unique()` to handle duplicates
# 6. **Posterior Probs:** Saved scANVI probabilities in `.obsm['scanvi_probabilities']`
# 
# ### 📊 What You Have Now
# 
# 1. **Fully annotated dataset** with multiple annotation layers:
#    - `scanvi_predictions`: **Primary annotation** (refined by scANVI)
#    - `scanvi_probabilities`: Full probability matrix for each cell type
#    - `predicted_labels`: CellTypist automatic annotation
#    - `majority_voting`: CellTypist with local refinement
# 
# 2. **Multiple embeddings** for different purposes:
#    - `X_scvi`: Batch-corrected latent space (75D)
#    - `X_scanvi`: Classification-optimized latent space (75D)
#    - `X_umap_scvi`: UMAP from scVI (batch-corrected visualization)
#    - `X_umap_scanvi`: UMAP from scANVI (classification visualization)
# 
# 3. **Raw counts preserved** in `.raw` for downstream analysis
# 
# 4. **Trained models** saved for reproducibility
# 
# ### 🔬 Recommended Next Steps
# 
# 1. **Validate annotations** using known T/NK cell markers
# 2. **Perform differential expression** between cell types
# 3. **Visualize posterior probabilities:**
#    ```python
#    import seaborn as sns
#    probs = adata.obsm['scanvi_probabilities']
#    celltypes = adata.uns['scanvi_celltype_order']
#    sns.heatmap(probs[:100, :], xticklabels=celltypes)
#    ```
# 4. **Compare UMAPs:**
#    ```python
#    # scVI UMAP
#    adata.obsm['X_umap'] = adata.obsm['X_umap_scvi']
#    sc.pl.umap(adata, color='scanvi_predictions')
#    
#    # scANVI UMAP
#    adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi']
#    sc.pl.umap(adata, color='scanvi_predictions')
#    ```
# 
# ### 📖 Loading Results Later
# 
# ```python
# import scanpy as sc
# import scvi
# 
# # Load annotated data
# adata = sc.read_h5ad("adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad")
# 
# # Access different UMAPs
# adata.obsm['X_umap'] = adata.obsm['X_umap_scvi']      # scVI UMAP
# adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi']    # scANVI UMAP
# 
# # Access posterior probabilities
# probs = adata.obsm['scanvi_probabilities']
# celltypes = adata.uns['scanvi_celltype_order']
# ```
# 
# ---
# 
# **Version:** 2.0  
# **Date:** 2025-12-06  
# **Status:** Production-ready with critical fixes applied


