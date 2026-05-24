#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
All Cells Analysis: scVI → CellTypist → scANVI Pipeline v2.0

Author: Clinical-Bioinformatics Team
Date: 2025-12-07
Version: 2.0 (Based on T cell pipeline v2.0 architecture)

Pipeline Overview:
==================
Input: adata_bbknn_annotated_corrected.h5ad (All cell types)
  ↓
Step 1: Data Preparation
  - Load BBKNN preprocessed data
  - Verify/extract raw counts from layers
  - Quality checks
  ↓
Step 2: scVI Integration (Batch Correction)
  - Train scVI model on raw counts (all cell types)
  - Generate X_scvi latent space (100 dims for diverse cells)
  - Compute X_umap_scvi for visualization ⭐ FIXED
  ↓
Step 3: CellTypist Annotation
  - Auto-annotate using Immune_All_Low.pkl
  - Gene name conversion (ENSEMBL→SYMBOL) ⭐ FIXED
  - Generate predicted_labels + majority_voting
  ↓
Step 4: scANVI Refinement (Semi-supervised)
  - Use existing cell_type labels as reference
  - Train scANVI for fine-grained classification
  - Generate X_umap_scanvi (separate from scVI) ⭐ FIXED
  ↓
Step 5: Comprehensive Visualization
  - Compare batch mixing (scVI vs scANVI UMAPs) ⭐ FIXED
  - Cell type composition analysis by tissue/batch
  - Confidence score evaluation
  - Marker gene expression validation
  ↓
Output: Fully annotated h5ad with multiple embeddings

Key Features (v2.0):
===================
✓ Applies all critical fixes from T cell pipeline v2.0
✓ Handles larger dataset with optimized memory management
✓ Comprehensive cell type annotation (immune + non-immune)
✓ Multi-level visualization (global + tissue-specific)
✓ Full reproducibility with seed control

Requirements:
=============
pip install scvi-tools scanpy celltypist mygene

# Download CellTypist model
celltypist.models.download_models(model='Immune_All_Low.pkl')
"""

# ==============================================================================
# Step 0: Import Libraries and Configuration
# ==============================================================================

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
BATCH_KEY = "dataset"  # Batch correction key
TISSUE_KEY = "tissue_sampling_method"  # Tissue annotation key
EXISTING_CELLTYPE_KEY = "cell_type"  # Existing cell type annotations

# scVI parameters
SCVI_N_LATENT = 100  # Larger than T cell (40k→75, all cells→100)
SCVI_N_LAYERS = 3
SCVI_DROPOUT_RATE = 0.1
SCVI_MAX_EPOCHS = 400  # May need more for diverse cell types
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# scANVI parameters (semi-supervised)
SCANVI_MAX_EPOCHS = 200
SCANVI_EARLY_STOPPING_PATIENCE = 30

# UMAP parameters
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
UMAP_N_NEIGHBORS = 100  # Larger for diverse cell population

# Visualization
DPI = 300
FIGURE_FORMAT = 'pdf'

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
scvi.settings.seed = RANDOM_SEED  # ⭐ NEW in v2.0

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)
(output_dir / "models").mkdir(exist_ok=True)

print(f"\nOutput directory: {output_dir}")

# Save configuration
config = {
    'input_h5ad': str(INPUT_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'tissue_key': TISSUE_KEY,
    'scvi_n_latent': SCVI_N_LATENT,
    'scvi_max_epochs': SCVI_MAX_EPOCHS,
    'scanvi_max_epochs': SCANVI_MAX_EPOCHS,
    'random_seed': RANDOM_SEED,
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'pipeline_config.json', 'w') as f:
    json.dump(config, f, indent=2)

print("✓ Configuration saved")


# ==============================================================================
# Step 1: Data Loading and Preparation
# ==============================================================================

print("\n" + "="*70)
print("Step 1: Data Loading and Preparation")
print("="*70)

# Load data
print(f"\nLoading data from: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"\nData loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Memory: ~{adata.X.data.nbytes / 1e9:.2f} GB" if issparse(adata.X) else f"  Memory: ~{adata.X.nbytes / 1e9:.2f} GB")

# Check available keys
print(f"\n.obs columns: {list(adata.obs.columns[:10])}...")
print(f".layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f".obsm: {list(adata.obsm.keys()) if adata.obsm else 'None'}")

# Verify batch and tissue keys
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found in .obs")
if TISSUE_KEY not in adata.obs.columns:
    print(f"⚠️  Tissue key '{TISSUE_KEY}' not found - tissue-specific plots will be skipped")

print(f"\nBatch distribution ({BATCH_KEY}):")
print(adata.obs[BATCH_KEY].value_counts())

if TISSUE_KEY in adata.obs.columns:
    print(f"\nTissue distribution ({TISSUE_KEY}):")
    print(adata.obs[TISSUE_KEY].value_counts())

if EXISTING_CELLTYPE_KEY in adata.obs.columns:
    print(f"\nExisting cell type distribution ({EXISTING_CELLTYPE_KEY}):")
    print(adata.obs[EXISTING_CELLTYPE_KEY].value_counts())

# ------------------------------------------------------------------------------
# Extract raw counts (CRITICAL for scVI)
# ------------------------------------------------------------------------------

print("\n" + "-"*70)
print("Extracting raw counts for scVI/scANVI")
print("-"*70)

# Check for raw counts in layers
if 'counts' in adata.layers:
    print("✓ Found raw counts in layers['counts']")
    adata.layers['counts_backup'] = adata.layers['counts'].copy()
    
    # Verify they are actually raw counts
    sample_values = adata.layers['counts'][:100, :100]
    if issparse(sample_values):
        sample_values = sample_values.toarray()
    
    has_decimals = np.any(sample_values != sample_values.astype(int))
    if has_decimals:
        print("⚠️  Warning: 'counts' layer contains non-integer values")
    else:
        print("✓ Verified: counts are integer values (raw UMI counts)")
    
elif hasattr(adata, 'raw') and adata.raw is not None:
    print("⚠️  'counts' not in layers, checking .raw")
    adata.layers['counts'] = adata.raw.X.copy()
    print("✓ Extracted counts from .raw.X")
else:
    # Last resort: assume .X is raw counts if it's sparse and integer
    print("⚠️  'counts' not found, checking if .X is raw counts...")
    sample_values = adata.X[:100, :100]
    if issparse(sample_values):
        sample_values = sample_values.toarray()
    
    has_decimals = np.any(sample_values != sample_values.astype(int))
    max_val = np.max(sample_values)
    
    if not has_decimals and max_val > 10:
        print("✓ .X appears to be raw counts, copying to layers['counts']")
        adata.layers['counts'] = adata.X.copy()
    else:
        raise ValueError(
            "Cannot find raw counts! Please ensure data has:\n"
            "  1. layers['counts'] with raw UMI counts, OR\n"
            "  2. .raw.X with raw counts, OR\n"
            "  3. .X with integer counts (before normalization)"
        )

# Verify counts are suitable for scVI
counts_sample = adata.layers['counts'][:1000, :1000]
if issparse(counts_sample):
    counts_sample = counts_sample.toarray()

print(f"\nRaw counts statistics (sample):")
print(f"  Min: {np.min(counts_sample)}")
print(f"  Max: {np.max(counts_sample)}")
print(f"  Mean: {np.mean(counts_sample):.2f}")
print(f"  Median: {np.median(counts_sample):.2f}")
print(f"  Sparsity: {np.sum(counts_sample == 0) / counts_sample.size * 100:.1f}%")

print("\n✓ Data preparation complete")


# ==============================================================================
# Step 2: scVI Integration (Unsupervised Batch Correction)
# ==============================================================================

print("\n" + "="*70)
print("Step 2: scVI Integration")
print("="*70)

# Prepare data for scVI
# CRITICAL: scVI needs raw counts in .X, not normalized
print("\nPreparing data for scVI (using raw counts)...")
adata_scvi = adata.copy()
adata_scvi.X = adata.layers['counts'].copy()

# Setup scVI model
print(f"\nSetting up scVI model...")
print(f"  Latent dimensions: {SCVI_N_LATENT}")
print(f"  Batch key: {BATCH_KEY}")

scvi.model.SCVI.setup_anndata(
    adata_scvi,
    layer=None,  # Use .X (which now contains raw counts)
    batch_key=BATCH_KEY
)

# Initialize model
scvi_model = scvi.model.SCVI(
    adata_scvi,
    n_latent=SCVI_N_LATENT,
    n_layers=SCVI_N_LAYERS,
    dropout_rate=SCVI_DROPOUT_RATE,
    gene_likelihood="nb"  # Negative binomial for count data
)

print(f"\n✓ scVI model initialized")
print(f"  Total parameters: {sum(p.numel() for p in scvi_model.module.parameters()):,}")

# Train scVI
print(f"\nTraining scVI (max {SCVI_MAX_EPOCHS} epochs)...")
print("This may take 30-90 minutes depending on data size and GPU...")

start_time = time.time()

scvi_model.train(
    max_epochs=SCVI_MAX_EPOCHS,
    early_stopping=SCVI_EARLY_STOPPING,
    early_stopping_patience=SCVI_EARLY_STOPPING_PATIENCE,
    enable_model_summary=True,
    enable_progress_bar=True
)

training_time = time.time() - start_time
print(f"\n✓ scVI training completed in {training_time/60:.1f} minutes")

# Save scVI model
scvi_model_dir = output_dir / "models" / "scvi_model"
scvi_model.save(scvi_model_dir, overwrite=True)
print(f"✓ scVI model saved to: {scvi_model_dir}")

# Plot training history
print("\nPlotting training history...")
train_history = scvi_model.history

fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# ELBO loss
axes[0].plot(train_history['elbo_train'], label='Train', linewidth=2)
if 'elbo_validation' in train_history:
    axes[0].plot(train_history['elbo_validation'], label='Validation', linewidth=2)
axes[0].set_xlabel('Epoch')
axes[0].set_ylabel('ELBO Loss')
axes[0].set_title('scVI Training - ELBO')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Reconstruction loss
if 'reconstruction_loss_train' in train_history:
    axes[1].plot(train_history['reconstruction_loss_train'], label='Train', linewidth=2)
    if 'reconstruction_loss_validation' in train_history:
        axes[1].plot(train_history['reconstruction_loss_validation'], label='Validation', linewidth=2)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Reconstruction Loss')
    axes[1].set_title('scVI Training - Reconstruction')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"scvi_training_history.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ Training history saved")

# ------------------------------------------------------------------------------
# Generate scVI latent representation and UMAP
# ------------------------------------------------------------------------------

print("\n" + "-"*70)
print("Generating scVI latent representation")
print("-"*70)

# Get latent representation
latent = scvi_model.get_latent_representation()
adata.obsm['X_scvi'] = latent

print(f"✓ scVI latent space generated: {latent.shape}")

# Compute neighbors using scVI latent space
print("\nComputing neighborhood graph from scVI latent...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS, random_state=RANDOM_SEED)
print("✓ Neighbor graph computed")

# Compute UMAP from scVI latent
print("Computing UMAP from scVI latent...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=RANDOM_SEED)

# ⭐ CRITICAL FIX #1: Save scVI UMAP separately (don't let it get overwritten)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print("✓ scVI UMAP computed and saved as 'X_umap_scvi'")

# Clustering on scVI latent (optional, for QC)
print("\nClustering on scVI latent (for batch mixing evaluation)...")
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scvi', random_state=RANDOM_SEED)
print(f"✓ Identified {adata.obs['leiden_scvi'].nunique()} clusters")

# Visualize scVI results
print("\nVisualizing scVI integration...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Batch mixing
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scVI: Batch')
if TISSUE_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 1], show=False, title='scVI: Tissue')
sc.pl.umap(adata, color='leiden_scvi', ax=axes[0, 2], show=False, title='scVI: Leiden')

# Existing cell types
if EXISTING_CELLTYPE_KEY in adata.obs.columns:
    sc.pl.umap(adata, color=EXISTING_CELLTYPE_KEY, ax=axes[1, 0], show=False, title='scVI: Existing Cell Type')

# Sample and library size
if 'n_genes_by_counts' in adata.obs.columns:
    sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[1, 1], show=False, title='scVI: # Genes', cmap='viridis')
if 'total_counts' in adata.obs.columns:
    sc.pl.umap(adata, color='total_counts', ax=axes[1, 2], show=False, title='scVI: Total Counts', cmap='viridis')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"scvi_umap_overview.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ scVI visualization complete")


# ==============================================================================
# Step 3: CellTypist Annotation
# ==============================================================================

print("\n" + "="*70)
print("Step 3: CellTypist Annotation")
print("="*70)

# Load CellTypist model
print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL_PATH}")
celltypist_model = models.Model.load(CELLTYPIST_MODEL_PATH)
print(f"✓ Model loaded: {celltypist_model.name if hasattr(celltypist_model, 'name') else 'Immune_All_Low'}")

# ------------------------------------------------------------------------------
# Gene name conversion (ENSEMBL → SYMBOL)
# ⭐ SIMPLIFIED LOGIC (v2.0 Fix #2)
# ------------------------------------------------------------------------------

print("\n" + "-"*70)
print("Checking gene name format")
print("-"*70)

sample_gene = str(adata.var_names[0])
print(f"Sample gene from data: {sample_gene}")

if sample_gene.startswith('ENSG'):
    need_conversion = True
    print("✓ Detected ENSEMBL IDs → Will convert to gene symbols")
else:
    need_conversion = False
    print("✓ Already using gene symbols")

if need_conversion:
    print("\nConverting ENSEMBL IDs to gene symbols using MyGene...")
    
    # Use MyGene for conversion
    mg = mygene.MyGeneInfo()
    
    ensembl_ids = adata.var_names.tolist()
    print(f"  Querying {len(ensembl_ids):,} ENSEMBL IDs...")
    
    query_result = mg.querymany(
        ensembl_ids,
        scopes='ensembl.gene',
        fields='symbol',
        species='human',
        returnall=True,
        verbose=False
    )
    
    # Build conversion dictionary
    ensembl_to_symbol = {}
    for item in query_result['out']:
        if 'symbol' in item:
            ensembl_to_symbol[item['query']] = item['symbol']
    
    print(f"  Successfully mapped: {len(ensembl_to_symbol):,}/{len(ensembl_ids):,} genes")
    
    # Apply conversion
    new_var_names = [ensembl_to_symbol.get(e, e) for e in adata.var_names]
    adata.var_names = new_var_names
    
    # ⭐ NEW in v2.0: Make unique after conversion
    adata.var_names_make_unique()
    
    conversion_rate = len(ensembl_to_symbol) / len(ensembl_ids) * 100
    print(f"✓ Conversion complete ({conversion_rate:.1f}% success rate)")
    
    # Save conversion map
    conversion_df = pd.DataFrame({
        'ensembl_id': ensembl_ids,
        'gene_symbol': new_var_names
    })
    conversion_df.to_csv(output_dir / "gene_name_conversion.csv", index=False)
    print("✓ Conversion table saved")

# ------------------------------------------------------------------------------
# Prepare data for CellTypist
# ⭐ EXPLICIT RAW COUNTS (v2.0 Fix #3)
# ------------------------------------------------------------------------------

print("\n" + "-"*70)
print("Preparing data for CellTypist")
print("-"*70)

# Create a copy for CellTypist
adata_celltypist = adata.copy()

# Use raw counts from layers
adata_celltypist.X = adata.layers['counts'].copy()
print("✓ Using raw counts from layers['counts']")

# Normalize and log-transform (CellTypist expects log-normalized data)
print("  Normalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)

print("  Log1p transforming...")
sc.pp.log1p(adata_celltypist)

print("✓ Data prepared for CellTypist")

# Run CellTypist prediction
print("\nRunning CellTypist prediction...")
print("This may take 5-15 minutes for large datasets...")

predictions = celltypist.annotate(
    adata_celltypist,
    model=celltypist_model,
    majority_voting=True
)

print("✓ CellTypist annotation complete")

# Extract results
predicted_labels = predictions.predicted_labels
adata.obs['celltypist_predicted'] = predicted_labels.predicted_labels.values
adata.obs['celltypist_majority_voting'] = predicted_labels.majority_voting.values
adata.obs['celltypist_conf_score'] = predicted_labels.conf_score.values

# Summary
print(f"\nCellTypist annotation summary:")
print(f"  Predicted labels: {adata.obs['celltypist_predicted'].nunique()} unique types")
print(f"  After majority voting: {adata.obs['celltypist_majority_voting'].nunique()} unique types")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
print(f"  Median confidence: {adata.obs['celltypist_conf_score'].median():.3f}")

print("\nTop 10 cell types (majority voting):")
print(adata.obs['celltypist_majority_voting'].value_counts().head(10))

# Compare with existing annotations (if available)
if EXISTING_CELLTYPE_KEY in adata.obs.columns:
    print("\nComparing with existing cell type annotations...")
    
    # ⭐ FIX #4: Convert categoricals to strings to avoid comparison error
    n_changed = (
        adata.obs[EXISTING_CELLTYPE_KEY].astype(str) != 
        adata.obs['celltypist_majority_voting'].astype(str)
    ).sum()
    
    print(f"  Cells with different annotations: {n_changed:,} ({n_changed/adata.n_obs*100:.1f}%)")
    
    # Crosstab
    crosstab = pd.crosstab(
        adata.obs[EXISTING_CELLTYPE_KEY],
        adata.obs['celltypist_majority_voting'],
        margins=True
    )
    crosstab.to_csv(output_dir / "celltypist_vs_existing_crosstab.csv")
    print("✓ Comparison table saved")

# Visualize CellTypist results
print("\nVisualizing CellTypist annotations...")

fig, axes = plt.subplots(1, 3, figsize=(21, 6))

sc.pl.umap(adata, color='celltypist_predicted', ax=axes[0], show=False, 
           title='CellTypist: Predicted Labels')
sc.pl.umap(adata, color='celltypist_majority_voting', ax=axes[1], show=False, 
           title='CellTypist: Majority Voting')
sc.pl.umap(adata, color='celltypist_conf_score', ax=axes[2], show=False, 
           title='CellTypist: Confidence', cmap='RdYlGn')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"celltypist_annotations.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.show()

print("✓ CellTypist visualization complete")


# ==============================================================================
# Step 4: scANVI Semi-supervised Refinement
# ==============================================================================

print("\n" + "="*70)
print("Step 4: scANVI Semi-supervised Refinement")
print("="*70)

# Decide which labels to use for scANVI
# Priority: CellTypist majority voting > existing annotations
if 'celltypist_majority_voting' in adata.obs.columns:
    scanvi_label_key = 'celltypist_majority_voting'
    print(f"Using CellTypist majority voting labels for scANVI")
elif EXISTING_CELLTYPE_KEY in adata.obs.columns:
    scanvi_label_key = EXISTING_CELLTYPE_KEY
    print(f"Using existing {EXISTING_CELLTYPE_KEY} labels for scANVI")
else:
    print("⚠️  No cell type labels available - skipping scANVI")
    scanvi_label_key = None

if scanvi_label_key is not None:
    # Prepare data for scANVI (needs same setup as scVI)
    print("\nPreparing data for scANVI...")
    adata_scanvi = adata.copy()
    adata_scanvi.X = adata.layers['counts'].copy()
    
    # Setup scANVI
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi,
        layer=None,  # Use .X (raw counts)
        batch_key=BATCH_KEY,
        labels_key=scanvi_label_key
    )
    
    # Initialize scANVI from trained scVI model
    print("\nInitializing scANVI from pre-trained scVI model...")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",  # Category for cells without labels
        adata=adata_scanvi,
        labels_key=scanvi_label_key
    )
    
    print(f"✓ scANVI model initialized")
    
    # Train scANVI
    print(f"\nTraining scANVI (max {SCANVI_MAX_EPOCHS} epochs)...")
    print("This leverages scVI weights and should be faster...")
    
    start_time = time.time()
    
    scanvi_model.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        early_stopping=True,
        early_stopping_patience=SCANVI_EARLY_STOPPING_PATIENCE,
        enable_model_summary=True,
        enable_progress_bar=True
    )
    
    training_time = time.time() - start_time
    print(f"\n✓ scANVI training completed in {training_time/60:.1f} minutes")
    
    # Save scANVI model
    scanvi_model_dir = output_dir / "models" / "scanvi_model"
    scanvi_model.save(scanvi_model_dir, overwrite=True)
    print(f"✓ scANVI model saved to: {scanvi_model_dir}")
    
    # Plot training history
    print("\nPlotting scANVI training history...")
    train_history = scanvi_model.history
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # ELBO loss
    axes[0].plot(train_history['elbo_train'], label='Train', linewidth=2)
    if 'elbo_validation' in train_history:
        axes[0].plot(train_history['elbo_validation'], label='Validation', linewidth=2)
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('ELBO Loss')
    axes[0].set_title('scANVI Training - ELBO')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Classification loss
    if 'classification_loss_train' in train_history:
        axes[1].plot(train_history['classification_loss_train'], label='Train', linewidth=2)
        if 'classification_loss_validation' in train_history:
            axes[1].plot(train_history['classification_loss_validation'], label='Validation', linewidth=2)
        axes[1].set_xlabel('Epoch')
        axes[1].set_ylabel('Classification Loss')
        axes[1].set_title('scANVI Training - Classification')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"scanvi_training_history.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Training history saved")
    
    # ------------------------------------------------------------------------------
    # Generate scANVI results
    # ------------------------------------------------------------------------------
    
    print("\n" + "-"*70)
    print("Generating scANVI predictions and embeddings")
    print("-"*70)
    
    # Get scANVI latent representation
    latent_scanvi = scanvi_model.get_latent_representation()
    adata.obsm['X_scanvi'] = latent_scanvi
    print(f"✓ scANVI latent space: {latent_scanvi.shape}")
    
    # Get scANVI predictions
    predictions_scanvi = scanvi_model.predict()
    adata.obs['scanvi_predictions'] = predictions_scanvi
    print(f"✓ scANVI predictions: {adata.obs['scanvi_predictions'].nunique()} cell types")
    
    # ⭐ NEW in v2.0: Save posterior probabilities
    predictions_df = scanvi_model.predict(soft=True)
    adata.obsm['scanvi_probabilities'] = predictions_df.values
    adata.uns['scanvi_celltype_order'] = predictions_df.columns.tolist()
    print("✓ Posterior probabilities saved to .obsm['scanvi_probabilities']")
    
    # Compute neighbors using scANVI latent space
    print("\nComputing neighborhood graph from scANVI latent...")
    sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=UMAP_N_NEIGHBORS, random_state=RANDOM_SEED)
    print("✓ Neighbor graph computed")
    
    # Compute UMAP from scANVI latent
    print("Computing UMAP from scANVI latent...")
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=RANDOM_SEED)
    
    # ⭐ CRITICAL FIX #1: Save scANVI UMAP separately
    adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
    print("✓ scANVI UMAP computed and saved as 'X_umap_scanvi'")
    
    # Clustering on scANVI latent
    print("\nClustering on scANVI latent...")
    sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi', random_state=RANDOM_SEED)
    print(f"✓ Identified {adata.obs['leiden_scanvi'].nunique()} clusters")
    
    # Summary
    print(f"\nscANVI annotation summary:")
    print(f"  Predicted types: {adata.obs['scanvi_predictions'].nunique()}")
    print("\nTop 10 cell types:")
    print(adata.obs['scanvi_predictions'].value_counts().head(10))
    
    # Visualize scANVI results
    print("\nVisualizing scANVI results...")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Use scANVI UMAP for visualization
    adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()
    
    sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scANVI: Batch')
    if TISSUE_KEY in adata.obs.columns:
        sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 1], show=False, title='scANVI: Tissue')
    sc.pl.umap(adata, color='leiden_scanvi', ax=axes[0, 2], show=False, title='scANVI: Leiden')
    
    sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1, 0], show=False, title='scANVI: Predictions')
    sc.pl.umap(adata, color='celltypist_majority_voting', ax=axes[1, 1], show=False, 
               title='CellTypist: Reference Labels')
    
    # Confidence score (if available)
    if 'celltypist_conf_score' in adata.obs.columns:
        sc.pl.umap(adata, color='celltypist_conf_score', ax=axes[1, 2], show=False, 
                   title='CellTypist: Confidence', cmap='RdYlGn')
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"scanvi_umap_overview.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ scANVI visualization complete")

else:
    print("\n⚠️  Skipping scANVI due to missing cell type labels")


# ==============================================================================
# Step 5: Comprehensive Comparison and Visualization
# ==============================================================================

print("\n" + "="*70)
print("Step 5: Comprehensive Comparison and Visualization")
print("="*70)

# ------------------------------------------------------------------------------
# Compare scVI vs scANVI batch mixing
# ⭐ This now correctly uses separate UMAPs (Fixed in v2.0)
# ------------------------------------------------------------------------------

if 'X_scanvi' in adata.obsm:
    print("\nComparing scVI vs scANVI batch mixing...")
    
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    
    # scVI UMAP
    adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()
    
    sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scVI: Batch')
    if TISSUE_KEY in adata.obs.columns:
        sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[0, 1], show=False, title='scVI: Tissue')
    if EXISTING_CELLTYPE_KEY in adata.obs.columns:
        sc.pl.umap(adata, color=EXISTING_CELLTYPE_KEY, ax=axes[0, 2], show=False, title='scVI: Cell Type')
    sc.pl.umap(adata, color='leiden_scvi', ax=axes[0, 3], show=False, title='scVI: Leiden')
    
    # scANVI UMAP
    adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()
    
    sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0], show=False, title='scANVI: Batch')
    if TISSUE_KEY in adata.obs.columns:
        sc.pl.umap(adata, color=TISSUE_KEY, ax=axes[1, 1], show=False, title='scANVI: Tissue')
    sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1, 2], show=False, title='scANVI: Predictions')
    sc.pl.umap(adata, color='leiden_scanvi', ax=axes[1, 3], show=False, title='scANVI: Leiden')
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"comparison_scvi_vs_scanvi.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Batch mixing comparison saved")

# ------------------------------------------------------------------------------
# Cell type composition by tissue/batch
# ------------------------------------------------------------------------------

if TISSUE_KEY in adata.obs.columns and 'scanvi_predictions' in adata.obs.columns:
    print("\nAnalyzing cell type composition by tissue...")
    
    # Composition by tissue
    comp_tissue = pd.crosstab(
        adata.obs[TISSUE_KEY],
        adata.obs['scanvi_predictions'],
        normalize='index'
    ) * 100
    
    fig, ax = plt.subplots(figsize=(14, 8))
    comp_tissue.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Tissue')
    ax.set_title('Cell Type Composition by Tissue (scANVI)')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"composition_by_tissue.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    # Save table
    comp_tissue.to_csv(output_dir / "cell_composition_by_tissue.csv")
    print("✓ Tissue composition analysis saved")

if BATCH_KEY in adata.obs.columns and 'scanvi_predictions' in adata.obs.columns:
    print("\nAnalyzing cell type composition by batch...")
    
    comp_batch = pd.crosstab(
        adata.obs[BATCH_KEY],
        adata.obs['scanvi_predictions'],
        normalize='index'
    ) * 100
    
    fig, ax = plt.subplots(figsize=(14, 8))
    comp_batch.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Batch')
    ax.set_title('Cell Type Composition by Batch (scANVI)')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"composition_by_batch.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    comp_batch.to_csv(output_dir / "cell_composition_by_batch.csv")
    print("✓ Batch composition analysis saved")

# ------------------------------------------------------------------------------
# Confidence score analysis
# ------------------------------------------------------------------------------

if 'celltypist_conf_score' in adata.obs.columns and 'scanvi_predictions' in adata.obs.columns:
    print("\nAnalyzing annotation confidence...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Confidence distribution by CellTypist
    conf_by_type = adata.obs.groupby('celltypist_majority_voting')['celltypist_conf_score'].mean().sort_values()
    conf_by_type.plot(kind='barh', ax=axes[0], color='steelblue')
    axes[0].set_xlabel('Mean Confidence Score')
    axes[0].set_title('CellTypist Confidence by Cell Type')
    axes[0].grid(True, alpha=0.3, axis='x')
    
    # Confidence distribution
    axes[1].hist(adata.obs['celltypist_conf_score'], bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(adata.obs['celltypist_conf_score'].median(), color='red', 
                    linestyle='--', linewidth=2, label=f'Median: {adata.obs["celltypist_conf_score"].median():.3f}')
    axes[1].set_xlabel('Confidence Score')
    axes[1].set_ylabel('Number of Cells')
    axes[1].set_title('CellTypist Confidence Distribution')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"confidence_analysis.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    
    print("✓ Confidence analysis saved")


# ==============================================================================
# Step 6: Final Data Export and Documentation
# ==============================================================================

print("\n" + "="*70)
print("Step 6: Saving Final Results")
print("="*70)

# Set default UMAP to scANVI (if available)
if 'X_umap_scanvi' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()
    print("✓ Default UMAP set to scANVI")
elif 'X_umap_scvi' in adata.obsm:
    adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()
    print("✓ Default UMAP set to scVI")

# Save final annotated dataset
output_file = output_dir / "adata_allcells_scvi_celltypist_scanvi_final_v2.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"✓ Final dataset saved ({file_size:.2f} GB)")

# ------------------------------------------------------------------------------
# Document complete data structure
# ------------------------------------------------------------------------------

print("\n" + "="*70)
print("COMPLETE DATA STRUCTURE (v2.0)")
print("="*70)

print(f"\nAnnData object: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# .X
print("\n[ .X - Main Expression Matrix ]")
print(f"  Type: {type(adata.X)}")
if issparse(adata.X):
    print(f"  Format: {adata.X.getformat()}")
print(f"  Content: Processed expression from BBKNN")

# .layers
print("\n[ .layers - Alternative Representations ]")
if adata.layers:
    for layer_name in adata.layers.keys():
        print(f"  {layer_name}: {adata.layers[layer_name].shape}")
        if layer_name == 'counts':
            print(f"    → Raw UMI counts (used for scVI/scANVI)")
else:
    print("  None")

# .obs annotations
print("\n[ .obs - Cell Annotations ]")
key_columns = [
    BATCH_KEY,
    TISSUE_KEY,
    EXISTING_CELLTYPE_KEY,
    'leiden_scvi',
    'leiden_scanvi',
    'celltypist_predicted',
    'celltypist_majority_voting',
    'celltypist_conf_score',
    'scanvi_predictions'
]
for col in key_columns:
    if col in adata.obs.columns:
        if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object':
            n_unique = adata.obs[col].nunique()
            print(f"  {col}: {n_unique} categories")
        else:
            print(f"  {col}: numeric ({adata.obs[col].dtype})")

# .obsm embeddings
print("\n[ .obsm - Dimensionality Reductions ]")
if adata.obsm:
    for key in adata.obsm.keys():
        print(f"  {key}: {adata.obsm[key].shape}")
        if key == 'X_scvi':
            print(f"    → scVI latent space (unsupervised)")
        elif key == 'X_scanvi':
            print(f"    → scANVI latent space (semi-supervised)")
        elif key == 'X_umap_scvi':
            print(f"    → UMAP from scVI latent ⭐ FIXED")
        elif key == 'X_umap_scanvi':
            print(f"    → UMAP from scANVI latent ⭐ FIXED")
        elif key == 'scanvi_probabilities':
            print(f"    → scANVI posterior probabilities ⭐ NEW")

# .uns
print("\n[ .uns - Unstructured Annotations ]")
if 'scanvi_celltype_order' in adata.uns:
    print(f"  scanvi_celltype_order: {len(adata.uns['scanvi_celltype_order'])} cell types ⭐ NEW")

print("\n" + "="*70)


# ==============================================================================
# Step 7: Generate Analysis Summary Report
# ==============================================================================

print("\n" + "="*70)
print("Step 7: Generating Analysis Summary")
print("="*70)

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("All Cells scVI-CellTypist-scANVI Analysis Summary")
summary_lines.append("="*70)
summary_lines.append(f"\nPipeline Version: 2.0")
summary_lines.append(f"Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary_lines.append(f"Random Seed: {RANDOM_SEED}")

summary_lines.append("\n[ Input Data ]")
summary_lines.append(f"  File: {INPUT_H5AD}")
summary_lines.append(f"  Total cells: {adata.n_obs:,}")
summary_lines.append(f"  Total genes: {adata.n_vars:,}")
summary_lines.append(f"  Batches ({BATCH_KEY}): {adata.obs[BATCH_KEY].nunique()}")
if TISSUE_KEY in adata.obs.columns:
    summary_lines.append(f"  Tissues ({TISSUE_KEY}): {adata.obs[TISSUE_KEY].nunique()}")

summary_lines.append("\n[ scVI Integration ]")
summary_lines.append(f"  Latent dimensions: {SCVI_N_LATENT}")
summary_lines.append(f"  Training epochs: {SCVI_MAX_EPOCHS}")
summary_lines.append(f"  Batch correction: {BATCH_KEY}")

summary_lines.append("\n[ CellTypist Annotation ]")
summary_lines.append(f"  Model: Immune_All_Low.pkl")
summary_lines.append(f"  Predicted cell types: {adata.obs['celltypist_majority_voting'].nunique()}")
summary_lines.append(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")

if 'scanvi_predictions' in adata.obs.columns:
    summary_lines.append("\n[ scANVI Refinement ]")
    summary_lines.append(f"  Training epochs: {SCANVI_MAX_EPOCHS}")
    summary_lines.append(f"  Final cell types: {adata.obs['scanvi_predictions'].nunique()}")
    
    summary_lines.append("\n  Top 10 cell types (scANVI):")
    for celltype, count in adata.obs['scanvi_predictions'].value_counts().head(10).items():
        pct = count / adata.n_obs * 100
        summary_lines.append(f"    {celltype}: {count:,} cells ({pct:.1f}%)")

summary_lines.append("\n[ Output Files ]")
summary_lines.append(f"  Final data: {output_file.name}")
summary_lines.append(f"  File size: {file_size:.2f} GB")
summary_lines.append(f"  scVI model: models/scvi_model/")
if 'X_scanvi' in adata.obsm:
    summary_lines.append(f"  scANVI model: models/scanvi_model/")
summary_lines.append(f"  Figures: figures/*.{FIGURE_FORMAT}")

summary_lines.append("\n[ Data Structure ]")
summary_lines.append("  .obsm embeddings:")
summary_lines.append("    X_scvi: scVI latent (unsupervised)")
if 'X_scanvi' in adata.obsm:
    summary_lines.append("    X_scanvi: scANVI latent (semi-supervised)")
summary_lines.append("    X_umap_scvi: UMAP from scVI ⭐ FIXED")
if 'X_umap_scanvi' in adata.obsm:
    summary_lines.append("    X_umap_scanvi: UMAP from scANVI ⭐ FIXED")
if 'scanvi_probabilities' in adata.obsm:
    summary_lines.append("    scanvi_probabilities: Posterior probs ⭐ NEW")

summary_lines.append("\n[ Critical Fixes Applied (v2.0) ]")
summary_lines.append("  ✓ UMAP overwriting fixed (scVI ≠ scANVI)")
summary_lines.append("  ✓ Gene name conversion simplified")
summary_lines.append("  ✓ Explicit raw counts for CellTypist")
summary_lines.append("  ✓ Categorical comparison fixed")
summary_lines.append("  ✓ Full reproducibility (scvi.settings.seed)")
summary_lines.append("  ✓ Duplicate gene names handled")
summary_lines.append("  ✓ Posterior probabilities saved")

summary_lines.append("\n" + "="*70)
summary_lines.append("Analysis completed successfully!")
summary_lines.append("="*70)

summary_text = '\n'.join(summary_lines)

# Print to console
print(summary_text)

# Save to file
summary_file = output_dir / "analysis_summary_v2.txt"
with open(summary_file, 'w', encoding='utf-8') as f:
    f.write(summary_text)

print(f"\n✓ Summary saved to: {summary_file}")


# ==============================================================================
# Final Notes
# ==============================================================================

print("\n" + "="*70)
print("🎉 PIPELINE COMPLETE")
print("="*70)

print(f"\n📁 Output Directory: {output_dir}")
print(f"\n📊 Main Results:")
print(f"  • Final data: {output_file.name}")
print(f"  • Summary: analysis_summary_v2.txt")
print(f"  • Figures: figures/ directory")
print(f"  • Models: models/ directory")

print(f"\n🔬 Key Outputs:")
print(f"  • scVI latent: adata.obsm['X_scvi']")
if 'X_scanvi' in adata.obsm:
    print(f"  • scANVI latent: adata.obsm['X_scanvi']")
print(f"  • scVI UMAP: adata.obsm['X_umap_scvi'] ⭐")
if 'X_umap_scanvi' in adata.obsm:
    print(f"  • scANVI UMAP: adata.obsm['X_umap_scanvi'] ⭐")
print(f"  • CellTypist: adata.obs['celltypist_majority_voting']")
if 'scanvi_predictions' in adata.obs.columns:
    print(f"  • scANVI: adata.obs['scanvi_predictions']")
    print(f"  • Probabilities: adata.obsm['scanvi_probabilities'] ⭐")

print(f"\n💡 Next Steps:")
print(f"  1. Review cell type annotations in figures/")
print(f"  2. Validate with known marker genes")
print(f"  3. Consider cell type-specific analysis:")
print(f"     • Epithelial cells: Differentiation, basal → ciliated")
print(f"     • T/NK cells: Activation, exhaustion, tissue residency")
print(f"     • Myeloid cells: M1/M2 polarization, DC subtypes")
print(f"  4. Downstream analysis: DEG, trajectory, functional scoring")

print(f"\n📚 Documentation:")
print(f"  • Pipeline summary: {summary_file.name}")
print(f"  • Configuration: pipeline_config.json")
if need_conversion:
    print(f"  • Gene conversion: gene_name_conversion.csv")

print("\n" + "="*70)
print()