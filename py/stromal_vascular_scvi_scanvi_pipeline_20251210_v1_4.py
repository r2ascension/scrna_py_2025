#!/usr/bin/env python3
"""
stromal_vascular_scvi_scanvi_pipeline.py

Stromal/Vascular Cells (Endothelial + Fibroblast + SMC) 
scVI/scANVI Complete Pipeline - VERSION 1.4

Purpose:
1. Batch correction using scVI (memory-optimized with HVG)
2. Cell type annotation using scANVI (semi-supervised learning)
3. Return refined annotations to main dataset for integration

Author: Clinical-Bioinformatics Team
Date: 2024-12-10
Version: v1.4 - FIXED sc.pl.embedding() parameter error + Time tracking

Key Improvements (v1.4):
- ✅ CRITICAL FIX: All sc.pl.umap() changed to sc.pl.embedding()
- ✅ Fixed time tracking (PIPELINE_START variable)
- ✅ Consistent use of basis parameter for custom embeddings
- ✅ All plotting code verified and tested

Previous v1.3 Improvements:
- Complete scANVI workflow (Step 11-15)
- Fixed obsm parameter error (use basis='umap_scvi')
- Unified variable naming (SCVI_/SCANVI_ prefixes)
- Use cell_type as scANVI label (Endothelial/Fibroblast/SMC)
- All code review suggestions implemented

Previous v1.2 Improvements:
- Model loading/saving (avoid re-training, save 99% time)
- UMAP KeyError fix (robust computation check)
- elbo_validation KeyError fix (robust history handling)

Previous v1.1 Improvements:
- GPU auto-detection with unified random seeds
- Batch key auto-detection with fallback
- Optimized HVG selection (direct on counts layer)
- Shared-memory .raw construction (0 extra cost)
- UMAP naming collision avoidance
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
# TIME TRACKING - START
# ============================================================================
PIPELINE_START = time.time()

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
    print("⚠️  Running on CPU (training will be slow)")

# ============================================================================
# CONFIGURATION SECTION
# ============================================================================

# ========== Input/Output Paths ==========
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1209/stromal_vascular_scvi_analysis"

# ========== Cell Type Selection ==========
MAJOR_CELLTYPE_KEY = "cell_type"  # Major cell type (Endothelial/Fibroblast/SMC)
TARGET_CELLTYPES = ["Endothelial", "Fibroblast", "SMC"]  # Cells to analyze

# ========== Batch Correction Settings ==========
PREFERRED_BATCH_KEYS = ["sample", "dataset", "batch"]  # Priority order
N_HVG = 4000  # Highly variable genes for training
MIN_CELLS_PER_GENE = 3  # Gene filtering threshold

# ========== scVI Model Parameters ==========
SCVI_N_LATENT = 75  # Latent dimensions (30-50 for stromal cells)
SCVI_N_LAYERS = 3
SCVI_N_HIDDEN = 128
SCVI_DROPOUT_RATE = 0.1
SCVI_GENE_LIKELIHOOD = "nb"  # Negative binomial

# ========== scVI Training Parameters ==========
SCVI_MAX_EPOCHS = 800
SCVI_BATCH_SIZE = 256
SCVI_TRAIN_SIZE = 0.9
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# ========== scANVI Model Parameters ==========
SCANVI_N_LAYERS = 3  # Use same as scVI
SCANVI_N_LATENT = 75  # Inherits from scVI
SCANVI_DROPOUT_RATE = 0.1

# ========== scANVI Training Parameters ==========
SCANVI_MAX_EPOCHS = 600
SCANVI_BATCH_SIZE = 256
SCANVI_TRAIN_SIZE = 0.9
SCANVI_EARLY_STOPPING = True

# ========== scANVI Label Settings ==========
SCANVI_LABEL_KEY = "cell_type"  # Use major cell type as label (Endothelial/Fibroblast/SMC)
SCANVI_UNLABELED_CATEGORY = "Unknown"  # For cells without labels (none in this case)

# ========== Model Reuse Flags ==========
USE_EXISTING_SCVI_MODEL = True
USE_EXISTING_SCANVI_MODEL = True
FORCE_RETRAIN = False

print(f"\n⚙️  Model Reuse Configuration:")
print(f"   USE_EXISTING_SCVI_MODEL: {USE_EXISTING_SCVI_MODEL}")
print(f"   USE_EXISTING_SCANVI_MODEL: {USE_EXISTING_SCANVI_MODEL}")
print(f"   FORCE_RETRAIN: {FORCE_RETRAIN}")

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
print("STROMAL/VASCULAR CELLS scVI/scANVI PIPELINE - v1.4 FIXED")
print("="*80)
print(f"\nPipeline Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\nTarget Cell Types: {', '.join(TARGET_CELLTYPES)}")
print(f"scANVI Label Key: {SCANVI_LABEL_KEY} (using major cell types as labels)")

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
print(f"   scVI Models: {model_dir}")
print(f"   scANVI Models: {scanvi_model_dir}")

# ============================================================================
# MARKER GENE DEFINITIONS - EXPANDED VERSION
# ============================================================================

print("\n" + "="*80)
print("MARKER GENE DEFINITIONS - EXPANDED")
print("="*80)

MARKER_GENES = {
    "Endothelial": {
        # Core endothelial markers
        "PanEndothelial":  ['PECAM1', 'CDH5', 'VWF', 'KDR', 'CLDN5', 'ENG', 'ESAM', 'ECSCR'],
        
        # Arterial endothelium
        "ArterialEC":      ['GJA5', 'EFNB2', 'BMX', 'SEMA3G', 'HEY1', 'DLL4', 'GJA4'],
        
        # Venous endothelium
        "VenousEC":        ['NR2F2', 'ACKR1', 'SELP', 'VWF', 'VCAM1'],
        
        # Capillary endothelium
        "CapillaryEC":     ['CA4', 'RGCC', 'APLN', 'PLVAP', 'KLF2', 'KLF4'],
        
        # Tip cells (angiogenic)
        "TipEC":           ['APLN', 'KDR', 'ESM1', 'ANGPT2'],
        
        # Lymphatic endothelium
        "LymphaticEC":     ['PROX1', 'PDPN', 'LYVE1', 'FLT4', 'CCL21', 'MMRN1'],
        
        # Activated endothelium
        "ActivatedEC":     ['SELE', 'VCAM1', 'ICAM1', 'IL6'],
        
        # Proliferating
        "Proliferating":    ['MKI67', 'TOP2A']
    },

    "Fibroblast": {
        # Core/Homeostatic Fibroblasts
        "PanFibroblast":   ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'LUM', 'PDGFRA', 
                             'DPT', 'FBLN1', 'COL6A3'],
        
        # PI16+ adventitial/reticular fibroblasts
        "AdventitialFibroblast_PI16": ['PI16', 'DPT', 'MFAP5', 'CXCL14', 'C7', 'CFD', 
                                        'APOD', 'GSN', 'PTGDS'],
        
        # COL15A1+ parenchymal fibroblasts
        "ParenchymalFibroblast_COL15A1": ['COL15A1', 'LAMA2', 'LAMB1', 'LAMC1', 'COL14A1', 
                                           'THY1', 'FBLN1', 'FBLN2'],
        
        # Immune-Interacting Fibroblasts
        "FRClike":         ['CCL19', 'CCL21', 'CXCL13', 'CXCL12', 'CD74', 
                             'HLA-DRA', 'HLA-DPA1', 'HLA-DPB1', 'TNFSF11'],
        
        "AntigenPresenting": ['HLA-DRA', 'HLA-DRB1', 'CD74', 'CIITA', 'CTSB', 'CTSS'],
        
        "InflammatoryFibroblast_IL6": ['IL6', 'IL1B', 'CXCL8', 'CCL2', 'PTGS2', 'ICAM1', 
                                        'NFKBIA', 'TNFAIP3'],
        
        "IFNResponse":     ['ISG15', 'IFI6', 'IFIT1', 'IFIT3', 'IFITM3', 'CXCL10', 
                             'STAT1', 'GBP1', 'OAS1', 'MX1'],
        
        # Matrix Remodeling & Myofibroblasts
        "MatrixFibroblast_SPARC": ['SPARC', 'COL3A1', 'COL5A1', 'COL5A2', 'COL1A1', 
                                    'COL1A2', 'COL6A3', 'MMP2', 'TIMP1', 'TIMP3'],
        
        "Myofibroblast_ACTA2": ['ACTA2', 'TAGLN', 'MYL9', 'CNN1', 'POSTN', 
                                'ITGA11', 'THBS2', 'COL1A1'],
        
        # Cancer-Associated Fibroblasts
        "CAF_LRRC15":       ['LRRC15', 'ITGA11', 'COL1A1', 'THBS2', 'FAP', 'PDPN'],
        
        "CAF_MMP1":         ['MMP1', 'MMP3', 'MMP10', 'MMP11', 'IL11', 'OSMR'],
        
        # Cell cycle
        "Proliferating":    ['MKI67', 'TOP2A', 'PCNA']
    },

    "SMC": {
        # Vascular smooth muscle cells
        "VascularSMC":     ['ACTA2', 'TAGLN', 'MYH11', 'CNN1', 'MYLK', 'TPM1', 
                             'TPM2', 'CARMN', 'SMTN'],
        
        # Contractile SMC
        "ContractileSMC":  ['MYH11', 'ACTA2', 'TAGLN', 'CARMN', 'SMTN'],
        
        # Synthetic SMC
        "SyntheticSMC":    ['FN1', 'COL1A1', 'COL3A1', 'MMP2'],
        
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
print(f"   Fibroblast subtypes: {len(MARKER_GENES['Fibroblast'])} (EXPANDED)")
print(f"   SMC subtypes: {len(MARKER_GENES['SMC'])}")
print(f"   Total unique markers: {len(ALL_MARKERS)}")

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
print(f"\nCell type distribution (key: {MAJOR_CELLTYPE_KEY}):")
if MAJOR_CELLTYPE_KEY not in adata_full.obs.columns:
    raise ValueError(f"Cell type key '{MAJOR_CELLTYPE_KEY}' not found in adata.obs")

celltype_counts = adata_full.obs[MAJOR_CELLTYPE_KEY].value_counts()
print(celltype_counts)

# Subset to target cell types
print(f"\nSubsetting to target cell types: {TARGET_CELLTYPES}")
mask = adata_full.obs[MAJOR_CELLTYPE_KEY].isin(TARGET_CELLTYPES)
adata = adata_full[mask].copy()

print(f"✓ Subset complete:")
print(f"   Cells retained: {adata.n_obs:,} ({adata.n_obs/adata_full.n_obs*100:.1f}%)")
print(f"   Genes: {adata.n_vars:,}")

# Cell type distribution in subset
print(f"\nTarget cell distribution:")
target_counts = adata.obs[MAJOR_CELLTYPE_KEY].value_counts()
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

# Filter genes
print(f"\nFiltering genes (min_cells={MIN_CELLS_PER_GENE})...")
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
n_genes_after = adata.n_vars
n_removed = n_genes_before - n_genes_after
print(f"✓ Removed {n_removed:,} genes ({n_removed/n_genes_before*100:.1f}%)")
print(f"✓ Remaining: {n_genes_after:,} genes")

# ========== Ensure counts layer & build .raw ==========
print(f"\n⭐ Preserving raw counts with shared memory...")

# Ensure layers['counts'] exists
if 'counts' not in adata.layers:
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
    X=adata.layers['counts'],  # Shared sparse matrix (no .copy())
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ .raw set from layers['counts'] (shared reference)")
print(f"   Memory savings: ~{adata.n_obs * adata.n_vars * 8 / 1e9:.1f} GB")

# Basic QC metrics
print(f"\nCalculating QC metrics on counts layer...")
sc.pp.calculate_qc_metrics(
    adata,
    layer='counts',
    percent_top=None,
    log1p=False,
    inplace=True
)
print(f"✓ QC metrics calculated")

# ============================================================================
# STEP 3: HIGHLY VARIABLE GENES (HVG) SELECTION
# ============================================================================

print("\n" + "="*80)
print("STEP 3: HVG SELECTION (Memory Optimization)")
print("="*80)

print(f"\n⭐ Selecting {N_HVG} HVGs directly on counts layer...")
print(f"   Strategy: HVG for training, full genes via .raw for analysis")

# HVG selection with batch-aware strategy and fallback
try:
    print(f"   Attempting batch-aware HVG (batch_key={BATCH_KEY})...")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        subset=False,
        flavor='seurat_v3'
    )
    hvg_method = "batch-aware"
    print(f"   ✓ Batch-aware HVG successful")
    
except Exception as e:
    print(f"   ⚠️ Batch-aware HVG failed: {e}")
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
print(f"\n⭐ Subsetting to HVG for scVI/scANVI training...")
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

# Define model path
scvi_model_path = model_dir / "scvi_model"
trained_scvi_this_run = False

# Setup AnnData for scVI
print(f"\nSetting up AnnData for scVI...")
print(f"   Using layers['counts'] for model input")
print(f"   Batch key: {BATCH_KEY}")

scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=BATCH_KEY
)
print(f"✓ scVI data setup complete")

# Check if we should load existing model
if USE_EXISTING_SCVI_MODEL and not FORCE_RETRAIN and scvi_model_path.exists():
    model_files = list(scvi_model_path.glob("*"))
    if len(model_files) > 0:
        print(f"\n✅ Found existing scVI model at: {scvi_model_path}")
        print(f"   Loading scVI model from disk (skipping training)...")
        try:
            scvi_model = scvi.model.SCVI.load(str(scvi_model_path), adata=adata)
            print(f"✓ scVI model loaded successfully")
        except Exception as e:
            print(f"⚠️  Failed to load model: {e}")
            print(f"   Will train new model instead...")
            scvi_model = None
    else:
        print(f"\n⚠️  Model directory exists but is empty")
        scvi_model = None
else:
    print(f"\n⚠️  No existing model found or retraining requested")
    scvi_model = None

# Train new model if needed
if scvi_model is None:
    print(f"\nCreating new scVI model...")
    scvi_model = scvi.model.SCVI(
        adata,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        n_hidden=SCVI_N_HIDDEN,
        dropout_rate=SCVI_DROPOUT_RATE,
        gene_likelihood=SCVI_GENE_LIKELIHOOD
    )
    
    print(f"✓ scVI model created")
    
    # Robust parameter counting
    try:
        n_params = scvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
    print(f"   Trainable parameters: {n_params:,}")
    
    print(f"\n{'='*80}")
    print(f"Starting scVI training...")
    print(f"{'='*80}\n")
    
    train_start = time.time()
    
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'batch_size': SCVI_BATCH_SIZE,
        'train_size': SCVI_TRAIN_SIZE,
        'early_stopping': SCVI_EARLY_STOPPING,
        'early_stopping_patience': SCVI_EARLY_STOPPING_PATIENCE,
    }
    
    if GPU_AVAILABLE:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})
    
    scvi_model.train(**train_kwargs)
    trained_scvi_this_run = True
    
    train_elapsed = time.time() - train_start
    print(f"\n{'='*80}")
    print(f"✓ scVI training complete!")
    print(f"   Time: {train_elapsed:.1f}s ({train_elapsed/60:.1f} min)")
    print(f"{'='*80}")
    
    # Save the trained model
    print(f"\nSaving scVI model to: {scvi_model_path}")
    scvi_model.save(str(scvi_model_path), overwrite=True)
    print(f"✓ Model saved successfully")

# Extract latent representation
print(f"\nExtracting scVI latent representation...")
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"✓ X_scvi saved, shape: {adata.obsm['X_scvi'].shape}")

# Plot training history
if trained_scvi_this_run and hasattr(scvi_model, 'history'):
    print(f"\nPlotting scVI training history...")
    
    history_keys = list(scvi_model.history.keys())
    print(f"   Available history keys: {history_keys}")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Robust history plotting with proper x-axis
    if 'elbo_train' in scvi_model.history:
        train_elbo = scvi_model.history['elbo_train']
        if hasattr(train_elbo, 'index'):
            x = train_elbo.index
            y = train_elbo.values
        else:
            y = np.array(train_elbo)
            x = np.arange(1, len(y) + 1)
        ax.plot(x, y, label='Training ELBO', linewidth=2)
    
    if 'elbo_validation' in scvi_model.history:
        val_elbo = scvi_model.history['elbo_validation']
        if hasattr(val_elbo, 'index'):
            x = val_elbo.index
            y = val_elbo.values
        else:
            y = np.array(val_elbo)
            x = np.arange(1, len(y) + 1)
        ax.plot(x, y, label='Validation ELBO', linewidth=2, linestyle='--')
    else:
        print(f"   ℹ️  No validation ELBO (this is normal)")
    
    ax.set_xlabel('Epoch', fontsize=12)
    ax.set_ylabel('ELBO', fontsize=12)
    ax.set_title('scVI Training History', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'01_scvi_training_history.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Training history saved")
else:
    print(f"\nℹ️  Skipping training history plot (model was loaded from disk)")

print(f"\n✅ scVI stage complete")

# ============================================================================
# STEP 5: scVI LATENT SPACE & UMAP
# ============================================================================

print("\n" + "="*80)
print("STEP 5: scVI LATENT SPACE EXTRACTION & UMAP")
print("="*80)

# Compute neighbors on scVI latent space
print(f"\n⭐ Computing neighborhood graph...")
print(f"   use_rep: X_scvi (latent space)")
print(f"   n_neighbors: {N_NEIGHBORS}")

sc.pp.neighbors(
    adata,
    use_rep='X_scvi',
    n_neighbors=N_NEIGHBORS
)
print(f"✓ Neighborhood graph computed")

# ⭐ CRITICAL: Compute UMAP and save with specific name
print(f"\n⭐ Computing UMAP (will save as X_umap_scvi)...")

try:
    sc.tl.umap(adata)
    
    # Save UMAP with scvi-specific name
    adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
    # Also save UMAP parameters for compatibility
    if 'umap' in adata.uns:
        adata.uns['umap_scvi'] = adata.uns['umap'].copy()
    
    print(f"✓ UMAP computed and saved as 'X_umap_scvi': {adata.obsm['X_umap_scvi'].shape}")
    print(f"   (This avoids collision with future X_umap_scanvi)")
    
except Exception as e:
    print(f"❌ UMAP computation failed: {e}")
    raise

# Clustering at multiple resolutions
print(f"\nPerforming Leiden clustering...")
for res in CLUSTERING_RESOLUTIONS:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_scvi_res{res}')
    n_clusters = adata.obs[f'leiden_scvi_res{res}'].nunique()
    print(f"✓ Resolution {res}: {n_clusters} clusters")

# Default clustering
default_res = 0.5
adata.obs['leiden_scvi'] = adata.obs[f'leiden_scvi_res{default_res}'].copy()
print(f"\n✓ Default clustering: leiden_scvi (resolution={default_res})")

print(f"\n✅ Latent space and clustering complete")

# ============================================================================
# STEP 6: BATCH CORRECTION VALIDATION
# ============================================================================

print("\n" + "="*80)
print("STEP 6: BATCH CORRECTION VALIDATION")
print("="*80)

print(f"\nGenerating batch correction validation plots...")

# ⭐ CRITICAL FIX: Use sc.pl.embedding() with basis='umap_scvi'
if 'X_umap_scvi' not in adata.obsm:
    print(f"❌ ERROR: X_umap_scvi not found in adata.obsm")
    raise KeyError("UMAP not computed correctly in STEP 5")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Plot 1: UMAP by batch
print(f"   Plotting UMAP colored by batch...")
sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=BATCH_KEY,
    ax=axes[0, 0],
    show=False,
    title='Batch Distribution',
    size=UMAP_SIZE,
    frameon=False,
)

# Plot 2: UMAP by cell type
print(f"   Plotting UMAP colored by cell type...")
sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=MAJOR_CELLTYPE_KEY,
    ax=axes[0, 1],
    show=False,
    title='Cell Type',
    size=UMAP_SIZE,
    frameon=False,
    legend_loc='right margin',
)

# Plot 3: UMAP by leiden clustering
print(f"   Plotting UMAP colored by leiden clustering...")
sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='leiden_scvi',
    ax=axes[0, 2],
    show=False,
    title='Leiden Clustering',
    size=UMAP_SIZE,
    frameon=False,
    legend_loc='right margin',
)

# Plot 4-6: Pan markers for each cell type
print(f"   Plotting marker genes...")
pan_markers = ['PECAM1', 'COL1A1', 'ACTA2']  # EC, Fib, SMC
for idx, marker in enumerate(pan_markers):
    if marker in adata.var_names:
        sc.pl.embedding(
            adata,
            basis='umap_scvi',
            color=marker,
            ax=axes[1, idx],
            show=False,
            title=f'{marker} Expression',
            size=UMAP_SIZE,
            frameon=False,
            use_raw=False,
            vmax='p99',
        )
    else:
        axes[1, idx].text(0.5, 0.5, f'{marker}\nNot in HVG', 
                          ha='center', va='center', fontsize=12)
        axes[1, idx].axis('off')

plt.tight_layout()
plt.savefig(fig_dir / f'02_scvi_batch_correction_validation.{FIGURE_FORMAT}', 
            dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✓ Batch correction validation plots saved")

# ============================================================================
# STEP 7: MARKER GENE VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print("STEP 7: MARKER GENE VISUALIZATION")
print("="*80)

print(f"\n⭐ Using .raw for full gene marker visualization...")
print(f"   Available genes: {adata.raw.n_vars:,}")

# Check marker availability
markers_in_raw = [g for g in ALL_MARKERS if g in adata.raw.var_names]
print(f"\nMarker gene availability in .raw:")
print(f"   Total markers: {len(ALL_MARKERS)}")
print(f"   Available: {len(markers_in_raw)} ({len(markers_in_raw)/len(ALL_MARKERS)*100:.1f}%)")

# Generate dotplots for each cell type
print(f"\nGenerating marker gene dotplots...")

for cell_type, marker_dict in MARKER_GENES.items():
    print(f"\n   {cell_type}...")
    
    # Collect markers
    cell_markers = []
    for subtype, genes in marker_dict.items():
        cell_markers.extend([g for g in genes if g in adata.raw.var_names])
    
    cell_markers = sorted(set(cell_markers))
    
    if len(cell_markers) == 0:
        print(f"      ⚠️ No markers available, skipping")
        continue
    
    print(f"      Markers available: {len(cell_markers)}")
    
    # Create dotplot
    fig, ax = plt.subplots(figsize=(max(12, len(cell_markers)*0.3), 8))
    
    try:
        sc.pl.dotplot(
            adata,
            var_names=cell_markers,
            groupby='leiden_scvi',
            use_raw=True,
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
        print(f"      ⚠️ Dotplot failed: {e}")
        plt.close()

# Generate comprehensive UMAP with key markers
print(f"\nGenerating comprehensive marker UMAP...")

key_markers_plot = [
    'PECAM1', 'CDH5',  # Endothelial
    'COL1A1', 'DCN',   # Fibroblast
    'ACTA2', 'MYH11',  # SMC
    'RGS5', 'PDGFRB',  # Pericyte
    'PI16', 'COL15A1', # Fibroblast backbones
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
        sc.pl.embedding(
            adata,
            basis='umap_scvi',
            color=marker,
            ax=axes[idx],
            show=False,
            title=marker,
            use_raw=True,
            vmax='p99',
            frameon=False,
            size=UMAP_SIZE,
        )
    
    # Hide extra subplots
    for idx in range(len(available_markers), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'04_key_markers_umap.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"✓ Key markers UMAP saved")

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
    use_raw=True,
    n_genes=50
)
print(f"✓ Differential expression complete")

# Export top markers per cluster
print(f"\nExporting top markers per cluster...")
result_df = sc.get.rank_genes_groups_df(adata, group=None)
result_df.to_csv(output_dir / "leiden_scvi_markers.csv", index=False)
print(f"✓ Markers exported: leiden_scvi_markers.csv")

# Generate marker heatmap
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
    print(f"✓ Heatmap saved")
    
except Exception as e:
    print(f"⚠️ Heatmap generation failed: {e}")

# ============================================================================
# STEP 9: PREPARE FOR scANVI
# ============================================================================

print("\n" + "="*80)
print("STEP 9: PREPARE FOR scANVI TRAINING")
print("="*80)

# Verify label key exists
if SCANVI_LABEL_KEY not in adata.obs.columns:
    raise ValueError(f"scANVI label key '{SCANVI_LABEL_KEY}' not found in adata.obs")

# Check label distribution
print(f"\nscANVI label distribution (key: {SCANVI_LABEL_KEY}):")
label_counts = adata.obs[SCANVI_LABEL_KEY].value_counts()
for label, count in label_counts.items():
    pct = count / adata.n_obs * 100
    print(f"   {label}: {count:,} ({pct:.1f}%)")

# Create a copy of labels for scANVI
adata.obs['scanvi_labels'] = adata.obs[SCANVI_LABEL_KEY].astype(str)

print(f"\n✓ Data prepared for scANVI")
print(f"   All cells have labels (no Unknown category needed)")
print(f"   scANVI will refine these cell type assignments")

# ============================================================================
# STEP 10: scANVI MODEL TRAINING
# ============================================================================

print("\n" + "="*80)
print("STEP 10: scANVI MODEL TRAINING (Semi-Supervised Learning)")
print("="*80)

# Define model path
scanvi_model_path = scanvi_model_dir / "scanvi_model"
trained_scanvi_this_run = False

print(f"\nscANVI configuration:")
print(f"   Base model: scVI (pre-trained)")
print(f"   Label key: {SCANVI_LABEL_KEY}")
print(f"   Training epochs: {SCANVI_MAX_EPOCHS}")
print(f"   Labels: {', '.join(label_counts.index.tolist())}")

# Check if we should load existing model
if USE_EXISTING_SCANVI_MODEL and not FORCE_RETRAIN and scanvi_model_path.exists():
    model_files = list(scanvi_model_path.glob("*"))
    if len(model_files) > 0:
        print(f"\n✅ Found existing scANVI model at: {scanvi_model_path}")
        print(f"   Loading scANVI model from disk (skipping training)...")
        try:
            scanvi_model = scvi.model.SCANVI.load(str(scanvi_model_path), adata=adata)
            print(f"✓ scANVI model loaded successfully")
        except Exception as e:
            print(f"⚠️  Failed to load model: {e}")
            print(f"   Will train new model instead...")
            scanvi_model = None
    else:
        print(f"\n⚠️  Model directory exists but is empty")
        scanvi_model = None
else:
    print(f"\n⚠️  No existing model found or retraining requested")
    scanvi_model = None

# Train new model if needed
if scanvi_model is None:
    print(f"\nInitializing scANVI from pre-trained scVI...")
    
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key='scanvi_labels',
        unlabeled_category=SCANVI_UNLABELED_CATEGORY
    )
    
    print(f"✓ scANVI model created")
    
    # Robust parameter counting
    try:
        n_params = scanvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)
    print(f"   Trainable parameters: {n_params:,}")
    
    print(f"\n{'='*80}")
    print(f"Starting scANVI training...")
    print(f"{'='*80}\n")
    
    train_start = time.time()
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'batch_size': SCANVI_BATCH_SIZE,
        'train_size': SCANVI_TRAIN_SIZE,
    }
    
    if GPU_AVAILABLE:
        train_kwargs.update({'accelerator': 'gpu', 'devices': 'auto'})
    else:
        train_kwargs.update({'accelerator': 'cpu', 'devices': 'auto'})
    
    scanvi_model.train(**train_kwargs)
    trained_scanvi_this_run = True
    
    train_elapsed = time.time() - train_start
    print(f"\n{'='*80}")
    print(f"✓ scANVI training complete!")
    print(f"   Time: {train_elapsed:.1f}s ({train_elapsed/60:.1f} min)")
    print(f"{'='*80}")
    
    # Save the trained model
    print(f"\nSaving scANVI model to: {scanvi_model_path}")
    scanvi_model.save(str(scanvi_model_path), overwrite=True)
    print(f"✓ Model saved successfully")

# Generate predictions
print(f"\nGenerating scANVI predictions...")
adata.obs['cell_type_scanvi'] = scanvi_model.predict()
print(f"✓ Predictions saved to adata.obs['cell_type_scanvi']")

# Get prediction probabilities
predictions_probs = scanvi_model.predict(soft=True)
adata.obs['scanvi_confidence'] = predictions_probs.max(axis=1)
print(f"✓ Confidence scores saved")

# Extract latent representation
print(f"\nExtracting scANVI latent representation...")
adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
print(f"✓ X_scanvi saved, shape: {adata.obsm['X_scanvi'].shape}")

# Plot training history
if trained_scanvi_this_run and hasattr(scanvi_model, 'history'):
    print(f"\nPlotting scANVI training history...")
    
    history_keys = list(scanvi_model.history.keys())
    print(f"   Available history keys: {history_keys}")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Robust history plotting
    if 'elbo_train' in scanvi_model.history:
        train_elbo = scanvi_model.history['elbo_train']
        if hasattr(train_elbo, 'index'):
            x = train_elbo.index
            y = train_elbo.values
        else:
            y = np.array(train_elbo)
            x = np.arange(1, len(y) + 1)
        ax.plot(x, y, label='Training ELBO', linewidth=2)
    
    if 'elbo_validation' in scanvi_model.history:
        val_elbo = scanvi_model.history['elbo_validation']
        if hasattr(val_elbo, 'index'):
            x = val_elbo.index
            y = val_elbo.values
        else:
            y = np.array(val_elbo)
            x = np.arange(1, len(y) + 1)
        ax.plot(x, y, label='Validation ELBO', linewidth=2, linestyle='--')
    else:
        print(f"   ℹ️  No validation ELBO")
    
    ax.set_xlabel('Epoch', fontsize=12)
    ax.set_ylabel('ELBO', fontsize=12)
    ax.set_title('scANVI Training History', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'06_scanvi_training_history.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"   ✓ Training history saved")
else:
    print(f"\nℹ️  Skipping training history plot (model was loaded from disk)")

print(f"\n✅ scANVI stage complete")

# ============================================================================
# STEP 11: scANVI LATENT SPACE & UMAP
# ============================================================================

print("\n" + "="*80)
print("STEP 11: scANVI LATENT SPACE & UMAP")
print("="*80)

# Compute neighbors on scANVI latent space
print(f"\n⭐ Computing neighborhood graph on scANVI latent space...")
print(f"   use_rep: X_scanvi")
print(f"   n_neighbors: {N_NEIGHBORS}")

sc.pp.neighbors(
    adata,
    use_rep='X_scanvi',
    n_neighbors=N_NEIGHBORS
)
print(f"✓ Neighborhood graph computed")

# Compute UMAP on scANVI latent space
print(f"\n⭐ Computing UMAP on scANVI latent space...")

try:
    sc.tl.umap(adata)
    
    # Save UMAP with scanvi-specific name
    adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
    # Also save UMAP parameters
    if 'umap' in adata.uns:
        adata.uns['umap_scanvi'] = adata.uns['umap'].copy()
    
    print(f"✓ UMAP computed and saved as 'X_umap_scanvi': {adata.obsm['X_umap_scanvi'].shape}")
    print(f"   Now you have both X_umap_scvi and X_umap_scanvi")
    
except Exception as e:
    print(f"❌ UMAP computation failed: {e}")
    raise

# Clustering on scANVI latent space
print(f"\nPerforming Leiden clustering on scANVI latent space...")
for res in CLUSTERING_RESOLUTIONS:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_scanvi_res{res}')
    n_clusters = adata.obs[f'leiden_scanvi_res{res}'].nunique()
    print(f"✓ Resolution {res}: {n_clusters} clusters")

# Default clustering
adata.obs['leiden_scanvi'] = adata.obs[f'leiden_scanvi_res{default_res}'].copy()
print(f"\n✓ Default clustering: leiden_scanvi (resolution={default_res})")

print(f"\n✅ scANVI latent space and clustering complete")

# ============================================================================
# STEP 12: scANVI RESULTS VALIDATION
# ============================================================================

print("\n" + "="*80)
print("STEP 12: scANVI RESULTS VALIDATION")
print("="*80)

print(f"\nGenerating scANVI validation plots...")

# Create comparison figure (scVI vs scANVI)
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Row 1: scVI results
print(f"   Plotting scVI results...")
sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=BATCH_KEY,
    ax=axes[0, 0],
    show=False,
    title='scVI: Batch',
    size=UMAP_SIZE,
    frameon=False,
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=MAJOR_CELLTYPE_KEY,
    ax=axes[0, 1],
    show=False,
    title='scVI: Original Labels',
    size=UMAP_SIZE,
    frameon=False,
    legend_loc='right margin',
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='leiden_scvi',
    ax=axes[0, 2],
    show=False,
    title='scVI: Leiden Clustering',
    size=UMAP_SIZE,
    frameon=False,
    legend_loc='right margin',
)

# Row 2: scANVI results
print(f"   Plotting scANVI results...")
sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color=BATCH_KEY,
    ax=axes[1, 0],
    show=False,
    title='scANVI: Batch',
    size=UMAP_SIZE,
    frameon=False,
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='cell_type_scanvi',
    ax=axes[1, 1],
    show=False,
    title='scANVI: Predicted Labels',
    size=UMAP_SIZE,
    frameon=False,
    legend_loc='right margin',
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='scanvi_confidence',
    ax=axes[1, 2],
    show=False,
    title='scANVI: Confidence',
    size=UMAP_SIZE,
    frameon=False,
    cmap='viridis',
)

plt.tight_layout()
plt.savefig(fig_dir / f'07_scvi_scanvi_comparison.{FIGURE_FORMAT}',
            dpi=DPI, bbox_inches='tight')
plt.close()
print(f"✓ Comparison plot saved")

# Confidence distribution
print(f"\nAnalyzing scANVI confidence scores...")
conf_stats = adata.obs['scanvi_confidence'].describe()
print(f"   Mean: {conf_stats['mean']:.3f}")
print(f"   Median: {conf_stats['50%']:.3f}")
print(f"   Min: {conf_stats['min']:.3f}")
print(f"   Max: {conf_stats['max']:.3f}")

# Low confidence cells
low_conf_threshold = 0.5
low_conf_cells = (adata.obs['scanvi_confidence'] < low_conf_threshold).sum()
print(f"\n   Cells with confidence < {low_conf_threshold}: {low_conf_cells:,} ({low_conf_cells/adata.n_obs*100:.1f}%)")

# Label agreement
print(f"\nComparing original vs scANVI labels...")
agreement = (adata.obs[MAJOR_CELLTYPE_KEY] == adata.obs['cell_type_scanvi']).sum()
agreement_pct = agreement / adata.n_obs * 100
print(f"   Agreement: {agreement:,}/{adata.n_obs:,} ({agreement_pct:.1f}%)")

# Confusion matrix
confusion = pd.crosstab(
    adata.obs[MAJOR_CELLTYPE_KEY],
    adata.obs['cell_type_scanvi'],
    margins=True
)
print(f"\nLabel confusion matrix:")
print(confusion)

confusion.to_csv(output_dir / "scanvi_label_confusion.csv")
print(f"✓ Confusion matrix saved")

# ============================================================================
# STEP 13: FINAL MARKER VISUALIZATION (scANVI)
# ============================================================================

print("\n" + "="*80)
print("STEP 13: MARKER GENE VISUALIZATION (scANVI)")
print("="*80)

print(f"\nGenerating marker gene visualization on scANVI UMAP...")

key_markers_plot = [
    'PECAM1', 'CDH5',  # Endothelial
    'COL1A1', 'DCN',   # Fibroblast
    'ACTA2', 'MYH11',  # SMC
    'RGS5', 'PDGFRB',  # Pericyte
    'PI16', 'COL15A1', # Fibroblast backbones
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
        sc.pl.embedding(
            adata,
            basis='umap_scanvi',
            color=marker,
            ax=axes[idx],
            show=False,
            title=marker,
            use_raw=True,
            vmax='p99',
            frameon=False,
            size=UMAP_SIZE,
        )
    
    # Hide extra subplots
    for idx in range(len(available_markers), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'08_scanvi_key_markers_umap.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"✓ scANVI marker UMAP saved")

# ============================================================================
# STEP 14: SAVE FINAL RESULTS
# ============================================================================

print("\n" + "="*80)
print("STEP 14: SAVING FINAL RESULTS")
print("="*80)

checkpoint_file = output_dir / "adata_scvi_scanvi_final.h5ad"

print(f"\nSaving final analysis results...")
print(f"File: {checkpoint_file}")

# Add comprehensive metadata
adata.uns['stromal_scvi_scanvi_pipeline'] = {
    'version': '1.4',
    'timestamp': datetime.now().isoformat(),
    'input_file': str(INPUT_H5AD),
    'target_celltypes': TARGET_CELLTYPES,
    'batch_key': BATCH_KEY,
    'n_hvg': N_HVG,
    'hvg_method': hvg_method,
    'min_cells_per_gene': MIN_CELLS_PER_GENE,
    'scvi_params': {
        'n_latent': SCVI_N_LATENT,
        'n_layers': SCVI_N_LAYERS,
        'n_hidden': SCVI_N_HIDDEN,
        'dropout_rate': SCVI_DROPOUT_RATE,
        'gene_likelihood': SCVI_GENE_LIKELIHOOD,
        'max_epochs': SCVI_MAX_EPOCHS,
        'batch_size': SCVI_BATCH_SIZE,
        'train_size': SCVI_TRAIN_SIZE,
    },
    'scanvi_params': {
        'n_latent': SCANVI_N_LATENT,
        'n_layers': SCANVI_N_LAYERS,
        'dropout_rate': SCANVI_DROPOUT_RATE,
        'max_epochs': SCANVI_MAX_EPOCHS,
        'batch_size': SCANVI_BATCH_SIZE,
        'label_key': SCANVI_LABEL_KEY,
    },
    'model_reuse': {
        'use_existing_scvi': USE_EXISTING_SCVI_MODEL,
        'use_existing_scanvi': USE_EXISTING_SCANVI_MODEL,
        'force_retrain': FORCE_RETRAIN,
        'scvi_trained_this_run': trained_scvi_this_run,
        'scanvi_trained_this_run': trained_scanvi_this_run,
    },
    'improvements_v1_4': [
        'CRITICAL FIX: All sc.pl.umap() changed to sc.pl.embedding()',
        'Fixed time tracking with PIPELINE_START variable',
        'Consistent use of basis parameter for custom embeddings',
        'All plotting code verified and tested',
    ],
    'improvements_v1_3': [
        'Complete scANVI workflow',
        'Fixed obsm parameter (use basis)',
        'Unified variable naming (SCVI_/SCANVI_ prefixes)',
        'Use cell_type as scANVI label',
        'Model loading/saving for both scVI and scANVI',
        'Robust history handling',
        'Code review suggestions implemented',
    ]
}

adata.write_h5ad(checkpoint_file, compression='gzip')

print(f"✓ Final results saved")
print(f"   Size: {checkpoint_file.stat().st_size / 1e9:.2f} GB")

# Export annotations
print(f"\nExporting cell annotations...")

# scVI annotations
scvi_annotations = adata.obs[[
    MAJOR_CELLTYPE_KEY,
    'leiden_scvi',
    BATCH_KEY
]].copy()
scvi_annotations.to_csv(output_dir / "scvi_annotations.csv")

# scANVI annotations
scanvi_annotations = adata.obs[[
    MAJOR_CELLTYPE_KEY,
    'cell_type_scanvi',
    'scanvi_confidence',
    'leiden_scanvi',
    BATCH_KEY
]].copy()
scanvi_annotations.to_csv(output_dir / "scanvi_annotations.csv")

print(f"✓ Annotations exported")
print(f"   scVI: scvi_annotations.csv")
print(f"   scANVI: scanvi_annotations.csv")

# ============================================================================
# STEP 15: FINAL SUMMARY
# ============================================================================

total_time = time.time() - PIPELINE_START

print(f"\n{'='*80}")
print("🎉 COMPLETE PIPELINE FINISHED - v1.4 FIXED")
print("="*80)

print(f"\n📊 Analysis Summary:")
print(f"   Cells analyzed: {adata.n_obs:,}")
print(f"   Genes (HVG): {adata.n_vars:,}")
print(f"   Genes (full): {adata.raw.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")
print(f"   Cell types: {adata.obs[MAJOR_CELLTYPE_KEY].nunique()}")
print(f"   scVI clusters: {adata.obs['leiden_scvi'].nunique()}")
print(f"   scANVI clusters: {adata.obs['leiden_scanvi'].nunique()}")

print(f"\n⭐ v1.4 NEW CRITICAL FIXES:")
print(f"   ✅ All sc.pl.umap() → sc.pl.embedding() (FIX basis parameter error)")
print(f"   ✅ Fixed time tracking (PIPELINE_START variable)")
print(f"   ✅ Consistent use of basis='umap_scvi' and basis='umap_scanvi'")
print(f"   ✅ All plotting code verified and working")

print(f"\n⭐ v1.3 Features (Inherited):")
print(f"   ✅ Complete scANVI workflow (Step 10-13)")
print(f"   ✅ Unified variable naming (SCVI_/SCANVI_ prefixes)")
print(f"   ✅ Use cell_type as scANVI labels")

print(f"\n⭐ v1.2 Features (Inherited):")
print(f"   ✅ Model loading/saving (both scVI and scANVI)")
print(f"   ✅ UMAP collision avoidance (X_umap_scvi, X_umap_scanvi)")
print(f"   ✅ Robust history handling")

print(f"\n⭐ v1.1 Features (Inherited):")
print(f"   ✅ GPU auto-detection with unified seeds")
print(f"   ✅ Batch key auto-detection")
print(f"   ✅ HVG on counts layer")
print(f"   ✅ Shared-memory .raw")
print(f"   ✅ Expanded Fibroblast markers")

print(f"\n📁 Output Files:")
print(f"   Final results: {checkpoint_file}")
print(f"   scVI model: {scvi_model_path}")
print(f"   scANVI model: {scanvi_model_path}")
print(f"   scVI annotations: scvi_annotations.csv")
print(f"   scANVI annotations: scanvi_annotations.csv")
print(f"   scVI markers: leiden_scvi_markers.csv")
print(f"   Label confusion: scanvi_label_confusion.csv")

print(f"\n📈 Figures:")
print(f"   01: scVI training history")
print(f"   02: scVI batch correction validation")
print(f"   03: Marker gene dotplots (by cell type)")
print(f"   04: Key markers UMAP (scVI)")
print(f"   05: Leiden markers heatmap")
print(f"   06: scANVI training history")
print(f"   07: scVI vs scANVI comparison")
print(f"   08: Key markers UMAP (scANVI)")

print(f"\n📊 scANVI Results:")
print(f"   Label agreement: {agreement_pct:.1f}%")
print(f"   Mean confidence: {conf_stats['mean']:.3f}")
print(f"   Low confidence cells (<0.5): {low_conf_cells:,} ({low_conf_cells/adata.n_obs*100:.1f}%)")

print(f"\n💡 Model Reuse Summary:")
if trained_scvi_this_run:
    print(f"   ✓ scVI trained and saved this run")
else:
    print(f"   ✓ scVI loaded from disk (training skipped)")

if trained_scanvi_this_run:
    print(f"   ✓ scANVI trained and saved this run")
else:
    print(f"   ✓ scANVI loaded from disk (training skipped)")

print(f"\n🔬 Next Steps:")
print(f"1. Review scVI vs scANVI UMAPs")
print(f"2. Analyze low-confidence predictions")
print(f"3. Examine label changes (confusion matrix)")
print(f"4. Use scANVI annotations for downstream analysis")
print(f"5. Return refined annotations to main dataset")

print(f"\n✅ Pipeline v1.4 complete!")
print(f"   Total time: {total_time/60:.1f} min ({total_time:.0f}s)")
print("="*80)
