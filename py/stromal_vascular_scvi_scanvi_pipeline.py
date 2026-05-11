#!/usr/bin/env python3
"""
stromal_vascular_scvi_scanvi_pipeline.py

Stromal/Vascular Cells (Endothelial + Fibroblast + SMC) 
scVI/scANVI Complete Pipeline - IMPROVED VERSION

Purpose:
1. Batch correction using scVI (memory-optimized with HVG)
2. Cell type annotation (automated + manual)
3. Semi-supervised refinement using scANVI
4. Return annotations to main dataset for integration

Author: Clinical-Bioinformatics Team
Date: 2024-12-09
Version: v1.1 - IMPROVED with robustness patches

Key Improvements (v1.1):
- GPU auto-detection with unified random seeds
- Batch key auto-detection with fallback
- Optimized HVG selection (direct on counts layer)
- Shared-memory .raw construction (0 extra cost)
- UMAP naming collision avoidance (X_umap_scvi)
- Expanded Fibroblast markers (FRC-like, CAF, IFN-response)
- Removed misleading n_pcs parameter
"""

import os
import sys
import warnings
from pathlib import Path
from datetime import datetime
import time
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
from typing import Optional, Dict, List

# Suppress warnings
warnings.filterwarnings('ignore')

# ============================================================================
# REPRODUCIBILITY & GPU AUTO-DETECTION
# ============================================================================

# --- Unified Random Seeds ---
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED
scvi.settings.dl_num_workers = 0  # More stable DataLoader

# --- GPU Auto-Detection ---
GPU_AVAILABLE = torch.cuda.is_available()
print("\n" + "="*80)
print("GPU & REPRODUCIBILITY CHECK")
print("="*80)
print(f"Random seed: {RANDOM_SEED}")
print(f"GPU available: {GPU_AVAILABLE}")
if GPU_AVAILABLE:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")
    print(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
else:
    print("Running on CPU (this will be slow)")

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# ========== Input/Output Paths ==========
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"  # Main integrated dataset
OUTPUT_DIR = "/home/h2048/data/py/1209/stromal_vascular_scvi_analysis"


# Cell type selection key in adata.obs
CELLTYPE_KEY = "cell_type"  # Column containing major cell types
TARGET_CELLTYPES = ["Endothelial", "Fibroblast", "SMC"]  # Cells to analyze

# ========== Batch Correction Settings ==========
# Note: BATCH_KEY will be auto-detected from data
PREFERRED_BATCH_KEYS = ["sample", "dataset", "batch"]  # Priority order
N_HVG = 4000  # Highly variable genes for training (memory optimization)

# ========== scVI Model Parameters ==========
N_LATENT = 50  # Latent dimensions (30-50 for stromal cells)
N_LAYERS = 2
N_HIDDEN = 128
DROPOUT_RATE = 0.1
GENE_LIKELIHOOD = "nb"  # Negative binomial

# ========== Training Parameters ==========
N_EPOCHS_SCVI = 400  # scVI training epochs
N_EPOCHS_SCANVI = 200  # scANVI fine-tuning epochs
BATCH_SIZE = 256  # Reduce if GPU memory limited
TRAIN_SIZE = 0.9
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 45

# ========== scANVI Settings ==========
SCANVI_LABEL_KEY = "cell_subtype_manual"  # For semi-supervised learning
UNLABELED_CATEGORY = "Unknown"

# ========== Clustering Parameters ==========
CLUSTERING_RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]
N_NEIGHBORS = 30

# ========== Visualization Settings ==========
DPI = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE = 2

# ========== Computational Resources ==========
N_JOBS = 48

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print("\n" + "="*80)
print("STROMAL/VASCULAR CELLS scVI/scANVI PIPELINE - v1.1 IMPROVED")
print("="*80)
print(f"\nPipeline Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\nTarget Cell Types: {', '.join(TARGET_CELLTYPES)}")

# Create output directory structure
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

model_dir = output_dir / "scvi_models"
model_dir.mkdir(exist_ok=True)

scanvi_model_dir = output_dir / "scanvi_models"
scanvi_model_dir.mkdir(exist_ok=True)

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.n_jobs = N_JOBS
sc.settings.figdir = fig_dir
sc.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)

print(f"\n✓ Output directories created:")
print(f"   Base: {output_dir}")
print(f"   Figures: {fig_dir}")
print(f"   Models: {model_dir}")

# ============================================================================
# MARKER GENE DEFINITIONS - EXPANDED VERSION
# ============================================================================

print("\n" + "="*80)
print("MARKER GENE DEFINITIONS - EXPANDED")
print("="*80)

MARKER_GENES = {
    "Endothelial": {
        # Core endothelial markers
        "Pan_Endothelial":  ['PECAM1', 'CDH5', 'VWF', 'KDR', 'CLDN5', 'ENG', 'ESAM', 'ECSCR'],
        
        # Arterial endothelium
        "Arterial_EC":      ['GJA5', 'EFNB2', 'BMX', 'SEMA3G', 'HEY1', 'DLL4', 'GJA4'],
        
        # Venous endothelium
        "Venous_EC":        ['NR2F2', 'ACKR1', 'SELP', 'VWF', 'VCAM1'],
        
        # Capillary endothelium
        "Capillary_EC":     ['CA4', 'RGCC', 'APLN', 'PLVAP', 'KLF2', 'KLF4'],
        
        # Tip cells (angiogenic)
        "Tip_EC":           ['APLN', 'KDR', 'ESM1', 'ANGPT2'],
        
        # Lymphatic endothelium
        "Lymphatic_EC":     ['PROX1', 'PDPN', 'LYVE1', 'FLT4', 'CCL21', 'MMRN1'],
        
        # Activated endothelium
        "Activated_EC":     ['SELE', 'VCAM1', 'ICAM1', 'IL6'],
        
        # Proliferating
        "Proliferating":    ['MKI67', 'TOP2A']
    },

    "Fibroblast": {
        # ========== Core/Homeostatic Fibroblasts ==========
        "Pan_Fibroblast":   ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRA', 
                             'DPT', 'FBLN1', 'COL6A3'],
        
        # PI16+ adventitial/reticular fibroblasts (universal backbone)
        "Adventitial_PI16": ['PI16', 'DPT', 'MFAP5', 'CXCL14', 'C7', 'CFD', 
                             'APOD', 'GSN', 'PTGDS'],
        
        # COL15A1+ parenchymal fibroblasts (tissue-specific backbone)
        "Parenchymal_COL15A1": ['COL15A1', 'LAMA2', 'LAMB1', 'LAMC1', 'COL14A1', 
                                'THY1', 'FBLN1', 'FBLN2'],
        
        # ========== Immune-Interacting Fibroblasts ==========
        # FRC-like (Fibroblastic Reticular Cell-like, lymphoid tissue)
        "FRC_like":         ['CCL19', 'CCL21', 'CXCL13', 'CXCL12', 'CD74', 
                             'HLA-DRA', 'HLA-DPA1', 'HLA-DPB1', 'TNFSF11'],
        
        # Antigen-presenting fibroblasts
        "Antigen_Presenting": ['HLA-DRA', 'HLA-DRB1', 'CD74', 'CIITA', 'CTSB', 'CTSS'],
        
        # Inflammatory IL6+ fibroblasts
        "Inflammatory_IL6": ['IL6', 'IL1B', 'CXCL8', 'CCL2', 'PTGS2', 'ICAM1', 
                             'NFKBIA', 'TNFAIP3'],
        
        # IFN-response fibroblasts
        "IFN_Response":     ['ISG15', 'IFI6', 'IFIT1', 'IFIT3', 'IFITM3', 'CXCL10', 
                             'STAT1', 'GBP1', 'OAS1', 'MX1'],
        
        # ========== Matrix Remodeling & Myofibroblasts ==========
        # SPARC+/COL3A1+ matrix fibroblasts
        "Matrix_SPARC_COL3A1": ['SPARC', 'COL3A1', 'COL5A1', 'COL5A2', 'COL1A1', 
                                'COL1A2', 'COL6A3', 'MMP2', 'TIMP1', 'TIMP3'],
        
        # ACTA2+ myofibroblasts
        "Myofibroblast_ACTA2": ['ACTA2', 'TAGLN', 'MYL9', 'CNN1', 'POSTN', 
                                'ITGA11', 'THBS2', 'COL1A1'],
        
        # ========== Cancer-Associated Fibroblasts (CAF) ==========
        # LRRC15+ CAF (associated with poor prognosis)
        "CAF_LRRC15":       ['LRRC15', 'ITGA11', 'COL1A1', 'THBS2', 'FAP', 'PDPN'],
        
        # MMP1+ inflammatory CAF
        "CAF_MMP1":         ['MMP1', 'MMP3', 'MMP10', 'MMP11', 'IL11', 'OSMR'],
        
        # Cell cycle
        "Proliferating":    ['MKI67', 'TOP2A', 'PCNA']
    },

    "SMC": {
        # Vascular smooth muscle cells
        "Vascular_SMC":     ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'MYLK', 'TPM1', 
                             'TPM2', 'CARMN', 'SMTN'],
        
        # Contractile SMC
        "Contractile_SMC":  ['MYH11', 'ACTA2', 'TAGLN', 'CARMN', 'SMTN'],
        
        # Synthetic SMC
        "Synthetic_SMC":    ['FN1', 'COL1A1', 'COL3A1', 'MMP2'],
        
        # Pericytes
        "Pericyte":         ['RGS5', 'PDGFRB', 'NOTCH3', 'MCAM', 'CSPG4', 'KCNJ8', 'ABCC9'],
        
        # Proliferating
        "Proliferating":    ['MKI67', 'TOP2A']
    }
}

# Flatten all markers for validation
ALL_MARKERS = sorted({
    gene 
    for cell_type_markers in MARKER_GENES.values() 
    for marker_list in cell_type_markers.values() 
    for gene in marker_list
})

print(f"\n✓ Marker genes defined:")
print(f"   Endothelial subtypes: {len(MARKER_GENES['Endothelial'])}")
print(f"   Fibroblast subtypes: {len(MARKER_GENES['Fibroblast'])} (⭐ EXPANDED)")
print(f"   SMC subtypes: {len(MARKER_GENES['SMC'])}")
print(f"   Total unique markers: {len(ALL_MARKERS)}")
print(f"\n   ⭐ New Fibroblast markers include:")
print(f"      - PI16+/COL15A1+ dual backbones")
print(f"      - FRC-like, Antigen-presenting")
print(f"      - Inflammatory IL6+, IFN-response")
print(f"      - LRRC15+/MMP1+ CAF subtypes")

# ============================================================================
# STEP 1: DATA LOADING & SUBSET
# ============================================================================

print("\n" + "="*80)
print("STEP 1: DATA LOADING & SUBSETTING")
print("="*80)

print(f"\nLoading main dataset: {INPUT_H5AD}")
if not Path(INPUT_H5AD).exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_H5AD}")

adata_full = sc.read_h5ad(INPUT_H5AD)
print(f"✓ Full dataset loaded:")
print(f"   Cells: {adata_full.n_obs:,}")
print(f"   Genes: {adata_full.n_vars:,}")

# ========== Batch Key Auto-Detection ==========
print(f"\n⭐ Auto-detecting batch key...")
BATCH_KEY = None
for key in PREFERRED_BATCH_KEYS:
    if key in adata_full.obs.columns:
        BATCH_KEY = key
        print(f"✓ Found batch key: '{BATCH_KEY}'")
        break

if BATCH_KEY is None:
    raise ValueError(f"No batch key found in obs. Tried: {PREFERRED_BATCH_KEYS}")

# Check cell type distribution
print(f"\nCell type distribution (key: {CELLTYPE_KEY}):")
if CELLTYPE_KEY not in adata_full.obs.columns:
    raise ValueError(f"Cell type key '{CELLTYPE_KEY}' not found in adata.obs")

celltype_counts = adata_full.obs[CELLTYPE_KEY].value_counts()
print(celltype_counts)

# Subset to target cell types
print(f"\nSubsetting to target cell types: {TARGET_CELLTYPES}")
mask = adata_full.obs[CELLTYPE_KEY].isin(TARGET_CELLTYPES)
adata = adata_full[mask].copy()

print(f"✓ Subset complete:")
print(f"   Cells retained: {adata.n_obs:,} ({adata.n_obs/adata_full.n_obs*100:.1f}%)")
print(f"   Genes: {adata.n_vars:,}")

# Cell type distribution in subset
print(f"\nTarget cell distribution:")
target_counts = adata.obs[CELLTYPE_KEY].value_counts()
for ct, count in target_counts.items():
    pct = count / adata.n_obs * 100
    print(f"   {ct}: {count:,} ({pct:.1f}%)")

# Batch distribution
print(f"\nBatch distribution (key: {BATCH_KEY}):")
batch_counts = adata.obs[BATCH_KEY].value_counts()
for batch, count in batch_counts.items():
    pct = count / adata.n_obs * 100
    print(f"   {batch}: {count:,} ({pct:.1f}%)")

# Clean up
del adata_full
gc.collect()
print(f"\n✓ Memory cleaned")

# ============================================================================
# STEP 2: DATA PREPROCESSING
# ============================================================================

print("\n" + "="*80)
print("STEP 2: DATA PREPROCESSING")
print("="*80)

# Filter genes (min 3 cells)
print(f"\nFiltering genes (min_cells=3)...")
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
n_genes_after = adata.n_vars
n_removed = n_genes_before - n_genes_after
print(f"✓ Removed {n_removed:,} genes ({n_removed/n_genes_before*100:.1f}%)")
print(f"✓ Remaining: {n_genes_after:,} genes")

# ========== Ensure counts layer & build .raw (IMPROVED) ==========
print(f"\n⭐ Preserving raw counts with shared memory...")

# Ensure layers['counts'] exists
if 'counts' not in adata.layers:
    # Priority: .raw.X > .X (assume counts)
    if adata.raw is not None:
        adata.layers['counts'] = adata.raw.X
        print("✓ layers['counts'] <- .raw.X (no copy)")
    else:
        if sparse.issparse(adata.X):
            adata.layers['counts'] = adata.X
        else:
            adata.layers['counts'] = sparse.csr_matrix(adata.X)
        print("✓ layers['counts'] <- .X (assumed counts)")
else:
    print("✓ layers['counts'] already exists")

# Build .raw with shared reference (0 extra memory cost)
adata.raw = sc.AnnData(
    X=adata.layers['counts'],  # Shared reference, NOT .copy()
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ .raw set from layers['counts'] (shared reference)")
print(f"   Memory savings: ~{adata.n_obs * adata.n_vars * 8 / 1e9:.1f} GB")

# Basic QC metrics
print(f"\nCalculating QC metrics...")
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
print(f"✓ QC metrics calculated")

# ============================================================================
# STEP 3: HIGHLY VARIABLE GENES (HVG) SELECTION - IMPROVED
# ============================================================================

print("\n" + "="*80)
print("STEP 3: HVG SELECTION (Memory Optimization) - IMPROVED")
print("="*80)

print(f"\n⭐ Selecting {N_HVG} HVGs directly on counts layer...")
print(f"   Strategy: HVG for training, full genes via .raw for analysis")

# HVG selection with batch-aware strategy and fallback
# IMPROVED: Direct on layer='counts', no .X manipulation
try:
    print(f"   Attempting batch-aware HVG (batch_key={BATCH_KEY})...")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',        # ⭐ KEY: Direct on counts, no .X write
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        subset=False,          # Don't subset yet
        flavor='seurat_v3'
    )
    hvg_method = "batch-aware"
    print(f"   ✓ Batch-aware HVG successful")
    
except Exception as e:
    print(f"   ⚠ Batch-aware HVG failed: {e}")
    print(f"   Falling back to non-batch-aware...")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        subset=False,
        flavor='seurat_v3'
    )
    hvg_method = "non-batch-aware"
    print(f"   ✓ Non-batch-aware HVG successful")

n_hvg = int(adata.var['highly_variable'].sum())
print(f"\n✓ Selected {n_hvg:,} HVGs using {hvg_method} method")

# Check marker gene coverage
markers_in_data = [g for g in ALL_MARKERS if g in adata.var_names]
markers_in_hvg = [g for g in markers_in_data if adata.var.loc[g, 'highly_variable']]
print(f"\n⭐ Marker gene coverage:")
print(f"   Total markers defined: {len(ALL_MARKERS)}")
print(f"   Markers in dataset: {len(markers_in_data)} ({len(markers_in_data)/len(ALL_MARKERS)*100:.1f}%)")
print(f"   Markers in HVG: {len(markers_in_hvg)} ({len(markers_in_hvg)/len(markers_in_data)*100:.1f}%)")

# Subset to HVG for training (full genes preserved in .raw)
print(f"\n⭐ Subsetting to HVG for scVI training...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"✓ Subset complete: {adata.n_vars:,} genes")
print(f"✓ Full gene access via adata.raw: {adata.raw.n_vars:,} genes")

# Normalize for visualization (not for training)
print(f"\nNormalizing HVG data for visualization...")
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
print(f"✓ Normalization complete")

# Save metadata
adata.uns['hvg_method'] = hvg_method
adata.uns['n_hvg'] = n_hvg

# ============================================================================
# STEP 4: scVI MODEL TRAINING
# ============================================================================

print("\n" + "="*80)
print("STEP 4: scVI MODEL TRAINING")
print("="*80)

print(f"\nPreparing data for scVI...")
print(f"   Using layers['counts'] for model input")
print(f"   Batch key: {BATCH_KEY}")

# Setup scVI
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=BATCH_KEY
)
print(f"✓ scVI data setup complete")

scvi_model = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    n_hidden=N_HIDDEN,
    dropout_rate=DROPOUT_RATE,
    gene_likelihood=GENE_LIKELIHOOD
)

print(f"✓ Model created")

# ---- Robust parameter counting (兼容无 n_params 的版本) ----
try:
    # 某些版本/示例里可能定义过 n_params
    n_params = scvi_model.module.n_params
except AttributeError:
    # 通用做法：统计所有可训练参数
    n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)

print(f"   Trainable parameters: {n_params:,}")


# Train scVI model
print(f"\nTraining scVI model...")
print(f"   Max epochs: {N_EPOCHS_SCVI}")
print(f"   Batch size: {BATCH_SIZE}")
print(f"   Train size: {TRAIN_SIZE}")
print(f"   Early stopping: {EARLY_STOPPING}")
print(f"   Device: {'GPU' if GPU_AVAILABLE else 'CPU'}")

train_start = time.time()

scvi_model.train(
    max_epochs=N_EPOCHS_SCVI,
    batch_size=BATCH_SIZE,
    train_size=TRAIN_SIZE,
    early_stopping=EARLY_STOPPING,
    early_stopping_patience=EARLY_STOPPING_PATIENCE
    # ,use_gpu=GPU_AVAILABLE
)

train_elapsed = time.time() - train_start
print(f"\n✓ scVI training complete: {train_elapsed:.1f}s ({train_elapsed/60:.1f} min)")

# Save model
scvi_model_path = model_dir / "scvi_model"
scvi_model.save(scvi_model_path, overwrite=True)
print(f"✓ Model saved: {scvi_model_path}")

# Training history visualization
print(f"\nGenerating training history plot...")

try:
    train_history = scvi_model.history.get('elbo_train', None)
    val_history = scvi_model.history.get('elbo_validation', None)
    
    if train_history is not None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Handle both Series and list/array
        if hasattr(train_history, 'index'):
            ax.plot(train_history.index, train_history.values, label='Training ELBO', linewidth=2)
        else:
            ax.plot(train_history, label='Training ELBO', linewidth=2)
        
        if val_history is not None:
            if hasattr(val_history, 'index'):
                ax.plot(val_history.index, val_history.values, label='Validation ELBO', 
                       linewidth=2, linestyle='--')
            else:
                ax.plot(val_history, label='Validation ELBO', linewidth=2, linestyle='--')
        
        ax.set_xlabel('Epoch', fontsize=12)
        ax.set_ylabel('ELBO', fontsize=12)
        ax.set_title('scVI Training History', fontsize=14, fontweight='bold')
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(fig_dir / f'01_scvi_training_history.{FIGURE_FORMAT}', dpi=DPI)
        plt.close()
        print(f"✓ Training plot saved: 01_scvi_training_history.{FIGURE_FORMAT}")
    else:
        print(f"⚠ Training history not available")
        
except Exception as e:
    print(f"⚠ Training plot generation failed: {e}")

# ============================================================================
# STEP 5: scVI LATENT SPACE & BATCH CORRECTION - IMPROVED
# ============================================================================

print("\n" + "="*80)
print("STEP 5: scVI LATENT SPACE EXTRACTION - IMPROVED")
print("="*80)

print(f"\nExtracting scVI latent representation...")
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"✓ Latent space extracted: {adata.obsm['X_scvi'].shape}")

# Compute neighbors on scVI latent space
# ⭐ IMPROVED: Remove n_pcs when use_rep='X_scvi' (it's ignored anyway)
print(f"\n⭐ Computing neighborhood graph...")
print(f"   use_rep: X_scvi (latent space)")
print(f"   n_neighbors: {N_NEIGHBORS}")
print(f"   Note: n_pcs parameter removed (not applicable for latent space)")

sc.pp.neighbors(
    adata,
    use_rep='X_scvi',
    n_neighbors=N_NEIGHBORS
    # ⭐ n_pcs removed - it's meaningless with use_rep='X_scvi'
)
print(f"✓ Neighborhood graph computed")

# Compute UMAP
# ⭐ IMPROVED: Save as X_umap_scvi to avoid collision with scANVI
print(f"\n⭐ Computing UMAP (saved as X_umap_scvi)...")
sc.tl.umap(adata)
adata.obsm['X_umap_scvi'] = adata.obsm.pop('X_umap')
print(f"✓ UMAP computed: {adata.obsm['X_umap_scvi'].shape}")
print(f"   Key: X_umap_scvi (avoids collision with future X_umap_scanvi)")

# Leiden clustering at multiple resolutions
print(f"\nPerforming Leiden clustering at multiple resolutions...")
for res in CLUSTERING_RESOLUTIONS:
    print(f"   Resolution {res}...")
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_scvi_res{res}')
    n_clusters = adata.obs[f'leiden_scvi_res{res}'].nunique()
    print(f"      → {n_clusters} clusters")

# Default clustering
default_res = 0.5
adata.obs['leiden_scvi'] = adata.obs[f'leiden_scvi_res{default_res}'].copy()
print(f"\n✓ Default clustering: leiden_scvi (resolution={default_res})")

# ============================================================================
# STEP 6: BATCH CORRECTION VALIDATION
# ============================================================================

print("\n" + "="*80)
print("STEP 6: BATCH CORRECTION VALIDATION")
print("="*80)

print(f"\nGenerating batch correction validation plots...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: UMAP by batch
# ⭐ IMPROVED: Use obsm='X_umap_scvi'
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, 
           title='Batch Distribution', s=UMAP_SIZE, frameon=False,
           obsm='X_umap_scvi')

# Plot 2: UMAP by cell type
sc.pl.umap(adata, color=CELLTYPE_KEY, ax=axes[0, 1], show=False,
           title='Cell Type', s=UMAP_SIZE, frameon=False, legend_loc='right margin',
           obsm='X_umap_scvi')

# Plot 3: UMAP by leiden clustering
sc.pl.umap(adata, color='leiden_scvi', ax=axes[0, 2], show=False,
           title='Leiden Clustering', s=UMAP_SIZE, frameon=False, legend_loc='right margin',
           obsm='X_umap_scvi')

# Plot 4-6: Pan markers for each cell type
pan_markers = ['PECAM1', 'COL1A1', 'ACTA2']  # EC, Fib, SMC
for idx, marker in enumerate(pan_markers):
    if marker in adata.var_names:
        sc.pl.umap(adata, color=marker, ax=axes[1, idx], show=False,
                   title=f'{marker} Expression', s=UMAP_SIZE, frameon=False, 
                   use_raw=False, vmax='p99', obsm='X_umap_scvi')
    else:
        axes[1, idx].text(0.5, 0.5, f'{marker}\nNot in HVG', 
                          ha='center', va='center', fontsize=12)
        axes[1, idx].axis('off')

plt.tight_layout()
plt.savefig(fig_dir / f'02_scvi_batch_correction_validation.{FIGURE_FORMAT}', 
            dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✓ Validation plot saved: 02_scvi_batch_correction_validation.{FIGURE_FORMAT}")

# ============================================================================
# STEP 7: MARKER GENE VISUALIZATION (using .raw for full gene access)
# ============================================================================

print("\n" + "="*80)
print("STEP 7: MARKER GENE VISUALIZATION")
print("="*80)

print(f"\n⭐ Using .raw for full gene marker visualization...")
print(f"   Available genes: {adata.raw.n_vars:,}")

# Check marker availability in raw data
markers_in_raw = [g for g in ALL_MARKERS if g in adata.raw.var_names]
print(f"\nMarker gene availability in .raw:")
print(f"   Total markers: {len(ALL_MARKERS)}")
print(f"   Available: {len(markers_in_raw)} ({len(markers_in_raw)/len(ALL_MARKERS)*100:.1f}%)")

# Generate dotplots for each cell type
print(f"\nGenerating marker gene dotplots...")

for cell_type, marker_dict in MARKER_GENES.items():
    print(f"\n   {cell_type}...")
    
    # Collect markers for this cell type
    cell_markers = []
    for subtype, genes in marker_dict.items():
        cell_markers.extend([g for g in genes if g in adata.raw.var_names])
    
    cell_markers = sorted(set(cell_markers))
    
    if len(cell_markers) == 0:
        print(f"      ⚠ No markers available, skipping")
        continue
    
    print(f"      Markers available: {len(cell_markers)}")
    
    # Create dotplot
    fig, ax = plt.subplots(figsize=(max(12, len(cell_markers)*0.3), 8))
    
    try:
        sc.pl.dotplot(
            adata,
            var_names=cell_markers,
            groupby='leiden_scvi',
            use_raw=True,  # ⭐ Use full gene data
            ax=ax,
            show=False,
            standard_scale='var'
        )
        
        plt.tight_layout()
        filename = f'03_markers_{cell_type.lower()}_dotplot.{FIGURE_FORMAT}'
        plt.savefig(fig_dir / filename, dpi=DPI, bbox_inches='tight')
        plt.close()
        print(f"      ✓ Saved: {filename}")
        
    except Exception as e:
        print(f"      ⚠ Dotplot failed: {e}")
        plt.close()

# Generate comprehensive UMAP with key markers
print(f"\nGenerating comprehensive marker UMAP...")

key_markers_plot = [
    'PECAM1', 'CDH5',  # Endothelial
    'COL1A1', 'DCN',   # Fibroblast
    'ACTA2', 'MYH11',  # SMC
    'RGS5', 'PDGFRB',  # Pericyte
    'PI16', 'COL15A1', # ⭐ NEW: Fibroblast backbones
    'MKI67'  # Proliferation
]

available_markers = [g for g in key_markers_plot if g in adata.raw.var_names]

if len(available_markers) > 0:
    n_cols = 3
    n_rows = (len(available_markers) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    axes = axes.flatten()
    
    for idx, marker in enumerate(available_markers):
        sc.pl.umap(
            adata,
            color=marker,
            ax=axes[idx],
            show=False,
            title=marker,
            use_raw=True,  # ⭐ Use full gene data
            vmax='p99',
            frameon=False,
            s=UMAP_SIZE,
            obsm='X_umap_scvi'  # ⭐ Use correct UMAP key
        )
    
    # Hide extra subplots
    for idx in range(len(available_markers), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'04_key_markers_umap.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"✓ Key markers UMAP saved: 04_key_markers_umap.{FIGURE_FORMAT}")

# ============================================================================
# STEP 8: DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================

print("\n" + "="*80)
print("STEP 8: DIFFERENTIAL EXPRESSION ANALYSIS")
print("="*80)

print(f"\nPerforming differential expression analysis...")
print(f"   Method: Wilcoxon rank-sum test")
print(f"   Using .raw for full gene analysis")

sc.tl.rank_genes_groups(
    adata,
    groupby='leiden_scvi',
    method='wilcoxon',
    use_raw=True,  # ⭐ Use full gene data
    n_genes=50
)
print(f"✓ Differential expression complete")

# Export top markers per cluster
print(f"\nExporting top markers per cluster...")
result_df = sc.get.rank_genes_groups_df(adata, group=None)
result_df.to_csv(output_dir / "leiden_scvi_markers.csv", index=False)
print(f"✓ Markers exported: leiden_scvi_markers.csv")

# Generate marker heatmap (top 5 per cluster)
print(f"\nGenerating marker gene heatmap...")

try:
    fig = sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=5,
        use_raw=True,
        show=False,
        cmap='RdBu_r',
        figsize=(12, 10)
    )
    
    plt.savefig(fig_dir / f'05_leiden_markers_heatmap.{FIGURE_FORMAT}', 
                dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"✓ Heatmap saved: 05_leiden_markers_heatmap.{FIGURE_FORMAT}")
    
except Exception as e:
    print(f"⚠ Heatmap generation failed: {e}")

# ============================================================================
# STEP 9: MANUAL ANNOTATION PREPARATION
# ============================================================================

print("\n" + "="*80)
print("STEP 9: MANUAL ANNOTATION PREPARATION")
print("="*80)

print(f"\nPreparing annotation template...")

# Create annotation template
clusters = adata.obs['leiden_scvi'].cat.categories.tolist()
annotation_template = pd.DataFrame({
    'Cluster': clusters,
    'CellType': ['Unknown'] * len(clusters),
    'Notes': [''] * len(clusters)
})

annotation_file = output_dir / "annotation_template.csv"
annotation_template.to_csv(annotation_file, index=False)

print(f"✓ Annotation template created: {annotation_file}")
print(f"\n⭐ Next steps for manual annotation:")
print(f"1. Review marker gene dotplots and UMAPs")
print(f"2. Examine differential expression results")
print(f"3. Fill in CellType column in: {annotation_file}")
print(f"4. ⭐ NEW: Use expanded Fibroblast markers for accurate classification")
print(f"5. Save and continue to scANVI refinement")

# Generate cluster summary
print(f"\nCluster summary:")
cluster_summary = adata.obs.groupby('leiden_scvi').agg({
    CELLTYPE_KEY: lambda x: x.value_counts().to_dict(),
    BATCH_KEY: 'count'
}).rename(columns={BATCH_KEY: 'n_cells'})

cluster_summary.to_csv(output_dir / "cluster_summary.csv")
print(f"✓ Cluster summary saved: cluster_summary.csv")

for cluster in clusters:
    n_cells = (adata.obs['leiden_scvi'] == cluster).sum()
    cell_types = adata.obs[adata.obs['leiden_scvi'] == cluster][CELLTYPE_KEY].value_counts()
    print(f"\n   Cluster {cluster}: {n_cells} cells")
    for ct, count in cell_types.head(3).items():
        pct = count / n_cells * 100
        print(f"      {ct}: {count} ({pct:.1f}%)")

# ============================================================================
# STEP 10: SAVE CHECKPOINT
# ============================================================================

print("\n" + "="*80)
print("STEP 10: SAVING CHECKPOINT")
print("="*80)

checkpoint_file = output_dir / "adata_scvi_checkpoint.h5ad"

print(f"\nSaving analysis checkpoint...")
print(f"File: {checkpoint_file}")

# Add metadata
adata.uns['scvi_params'] = {
    'n_latent': N_LATENT,
    'n_layers': N_LAYERS,
    'n_hidden': N_HIDDEN,
    'epochs': N_EPOCHS_SCVI,
    'batch_size': BATCH_SIZE,
    'n_hvg': N_HVG,
    'hvg_method': hvg_method,
    'batch_key': BATCH_KEY,
    'random_seed': RANDOM_SEED,
    'gpu_used': GPU_AVAILABLE
}

adata.uns['improvements_v1_1'] = {
    'gpu_auto_detection': True,
    'batch_key_auto_detection': True,
    'hvg_on_counts_layer': True,
    'shared_memory_raw': True,
    'umap_naming_fixed': 'X_umap_scvi',
    'expanded_fibroblast_markers': True
}

adata.write_h5ad(checkpoint_file, compression='gzip')

print(f"✓ Checkpoint saved")
print(f"   Size: {checkpoint_file.stat().st_size / 1e9:.2f} GB")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print("🎉 scVI ANALYSIS COMPLETE (Stage 1) - v1.1 IMPROVED")
print("="*80)

print(f"\n📊 Analysis Summary:")
print(f"   Cells analyzed: {adata.n_obs:,}")
print(f"   Genes (HVG): {adata.n_vars:,}")
print(f"   Genes (full): {adata.raw.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")
print(f"   Batch key (auto-detected): {BATCH_KEY}")
print(f"   Clusters: {adata.obs['leiden_scvi'].nunique()}")
print(f"   GPU used: {GPU_AVAILABLE}")
print(f"   Random seed: {RANDOM_SEED}")

print(f"\n⭐ v1.1 Improvements Applied:")
print(f"   ✅ GPU auto-detection with unified random seeds")
print(f"   ✅ Batch key auto-detection (fallback strategy)")
print(f"   ✅ HVG selection directly on counts layer")
print(f"   ✅ Shared-memory .raw construction (0 extra cost)")
print(f"   ✅ UMAP naming fixed (X_umap_scvi)")
print(f"   ✅ Removed misleading n_pcs parameter")
print(f"   ✅ Expanded Fibroblast markers (+8 subtypes)")

print(f"\n📁 Output Files:")
print(f"   Checkpoint: {checkpoint_file}")
print(f"   scVI model: {scvi_model_path}")
print(f"   Markers: {output_dir / 'leiden_scvi_markers.csv'}")
print(f"   Annotation template: {annotation_file}")

print(f"\n📈 Figures:")
print(f"   Training history: 01_scvi_training_history.{FIGURE_FORMAT}")
print(f"   Batch validation: 02_scvi_batch_correction_validation.{FIGURE_FORMAT}")
print(f"   Marker dotplots: 03_markers_*.{FIGURE_FORMAT}")
print(f"   Key markers UMAP: 04_key_markers_umap.{FIGURE_FORMAT}")
print(f"   Marker heatmap: 05_leiden_markers_heatmap.{FIGURE_FORMAT}")

print(f"\n🔬 Next Steps:")
print(f"1. Review visualization results in: {fig_dir}")
print(f"2. Perform manual annotation: {annotation_file}")
print(f"3. ⭐ Use expanded Fibroblast marker dictionary for classification")
print(f"4. Run scANVI semi-supervised refinement (stromal_vascular_scanvi_refinement.py)")
print(f"5. Return annotations to main dataset")

print(f"\n💡 Memory Optimization:")
print(f"   HVG training: {N_HVG:,} genes (reduced memory by ~90%)")
print(f"   Full gene access: via .raw (shared memory, 0 cost)")
print(f"   Peak memory: ~{adata.n_obs * N_HVG * 8 / 1e9:.1f} GB (estimate)")

print(f"\n✅ Ready for manual annotation and scANVI refinement!")
print("="*80)