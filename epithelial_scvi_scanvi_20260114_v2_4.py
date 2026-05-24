#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Epithelial Cell scVI-scANVI Training Pipeline (v2.4-FIXED)
With Robust HVG Selection & Hierarchical Label Consolidation

New in v2.4-FIXED:
- ⭐ Robust HVG selection with try-except fallback
- ⭐ Data validation for infinity/NaN values
- ⭐ Increased MIN_CELLS_PER_BATCH to prevent LOESS failures
- ⭐ Complete integration of QUICK_REFERENCE_MEMORY best practices
- ⭐ Version-compatible parameter counting

Based on v2.3-REUSABLE with critical bug fixes

Author: r2end
Date: 2025-01-14
Version: 2.4-FIXED - Production with robust error handling
"""

import os
import gc
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from scipy import sparse
from sklearn.metrics import confusion_matrix
from pathlib import Path

warnings.filterwarnings('ignore')

# Set plotting style
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

print("="*80)
print("Epithelial scVI-scANVI Pipeline (v2.4-FIXED)")
print("With Robust HVG Selection & Model Reuse")
print("="*80)
print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"scanpy version: {sc.__version__}")
print(f"scvi-tools version: {scvi.__version__}")
print(f"PyTorch version: {torch.__version__}")
print(f"Device: {'GPU' if torch.cuda.is_available() else 'CPU'}")
print("="*80)

# ===== CONFIGURATION =====

# Input/Output paths
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0114/celltypist_epithelial"
CHECKPOINT_DIR = f"{OUTPUT_DIR}/checkpoints"

# ⭐ Model Reuse Configuration
AUTO_LOAD_MODELS = False  # Auto-detect and load existing models
FORCE_RETRAIN_SCVI = False  # Force retrain scVI even if exists
FORCE_RETRAIN_SCANVI = False  # Force retrain scANVI even if exists

# Paths for pretrained models (optional)
PRETRAINED_SCVI_PATH = None
PRETRAINED_SCANVI_PATH = None

# Create directories
for dir_path in [OUTPUT_DIR, f"{OUTPUT_DIR}/figures", f"{OUTPUT_DIR}/scvi_model", 
                 f"{OUTPUT_DIR}/scanvi_model", CHECKPOINT_DIR]:
    os.makedirs(dir_path, exist_ok=True)

# Key column names
BATCH_KEY = 'sample'
CELLTYPE_KEY = 'celltypist_pred'
UNLABELED_CATEGORY = 'Unknown'

# ===== Hierarchical Label Configuration =====
USE_HIERARCHICAL_LABELS = True
VISUALIZE_HIERARCHY = True

# Hierarchical label mapping
MAJOR_LINEAGE_MAP = {
    # Alveolar lineage
    'AT1': 'Alveolar',
    'AT1 ': 'Alveolar',
    'AT2': 'Alveolar',
    
    # Basal lineage
    'Basal': 'Basal_Lineage',
    'Suprabasal': 'Basal_Lineage',
    'SMG_Basal': 'Basal_Lineage',
    'Dividing_Basal': 'Basal_Lineage',
    
    # Ciliated lineage
    'Ciliated': 'Ciliated_Lineage',
    'Deuterosome': 'Ciliated_Lineage',
    'Deuterosomal': 'Ciliated_Lineage',
    
    # Secretory lineage
    'Secretory_Goblet': 'Secretory_Lineage',
    'Secretory_Club': 'Secretory_Lineage',
    'SMG_Mucous': 'Secretory_Lineage',
    'SMG_Serous': 'Secretory_Lineage',
    'SCGB1A1+': 'Secretory_Lineage',
    
    # Duct cells
    'SMG_Duct': 'Duct',
    
    # Rare specialized cells
    'Ionocyte_n_Brush': 'Rare_Specialized',
    'Ionocyte': 'Rare_Specialized',
    'Brush': 'Rare_Specialized'
}

# Force-include marker genes
FORCE_INCLUDE_MARKERS = [
    # Alveolar markers
    'SFTPC', 'SFTPB', 'ABCA3', 'SLC34A2', 'ETV5',  # AT2
    'AGER', 'PDPN', 'CLDN18', 'CAV1',  # AT1
    
    # Basal markers
    'TP63', 'KRT5', 'KRT14', 'KRT15',
    
    # Suprabasal markers
    'KRT13', 'KRT4',
    
    # Ciliated markers
    'FOXJ1', 'RSPH1', 'PIFO', 'DNAH5',
    
    # Deuterosome markers
    'DEUP1', 'CCNO',
    
    # Secretory markers
    'MUC5AC', 'MUC5B', 'SCGB1A1', 'SPDEF',
    'BPIFA1', 'WFDC2',
    
    # Rare cells
    'CFTR', 'FOXI1',  # Ionocyte
    'DCLK1', 'TRPM5',  # Brush
    
    # SMG markers
    'LYZ', 'LTF',  # Serous
    'KRT7', 'KRT19'  # Duct
]

# Color schemes
MAJOR_LINEAGE_COLORS = {
    'Alveolar': '#E74C3C',
    'Basal_Lineage': '#3498DB',
    'Ciliated_Lineage': '#2ECC71',
    'Secretory_Lineage': '#F39C12',
    'Duct': '#9B59B6',
    'Rare_Specialized': '#95A5A6'
}

# scVI parameters
N_HVG = 4000
N_LATENT_SCVI = 100
N_LAYERS = 3
N_HIDDEN = 256
DROPOUT_RATE = 0.1
GENE_LIKELIHOOD = "nb"
DISPERSION = "gene-batch"

# scANVI parameters
N_LATENT_SCANVI = 75

# Training parameters
SCVI_MAX_EPOCHS = 400
SCANVI_MAX_EPOCHS = 200
EARLY_STOPPING = True
BATCH_SIZE = 1024

# QC thresholds
MIN_CELLS_PER_BATCH = 30  # ⭐ INCREASED from 10 to prevent LOESS failures
LOW_CONFIDENCE_THRESHOLD = 0.3
MAX_LOW_CONFIDENCE_PCT = 80

# Options
SAVE_CHECKPOINTS = True
MARK_LOW_CONFIDENCE_AS_UNKNOWN = True

# Random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# Device
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Model path detection
scvi_model_path = PRETRAINED_SCVI_PATH or f"{OUTPUT_DIR}/scvi_model"
scanvi_model_path = PRETRAINED_SCANVI_PATH or f"{OUTPUT_DIR}/scanvi_model"

scvi_model_exists = Path(scvi_model_path).exists() and Path(f"{scvi_model_path}/model.pt").exists()
scanvi_model_exists = Path(scanvi_model_path).exists() and Path(f"{scanvi_model_path}/model.pt").exists()

# Training decision logic
TRAIN_SCVI = FORCE_RETRAIN_SCVI or (not scvi_model_exists) or (not AUTO_LOAD_MODELS)
TRAIN_SCANVI = FORCE_RETRAIN_SCANVI or (not scanvi_model_exists) or (not AUTO_LOAD_MODELS)

print("\n" + "="*80)
print("CONFIGURATION")
print("="*80)
print(f"Input file: {INPUT_H5AD}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Device: {DEVICE}")

print(f"\n⭐ Model Reuse Strategy:")
print(f"  AUTO_LOAD_MODELS: {AUTO_LOAD_MODELS}")
print(f"  FORCE_RETRAIN_SCVI: {FORCE_RETRAIN_SCVI}")
print(f"  FORCE_RETRAIN_SCANVI: {FORCE_RETRAIN_SCANVI}")

print(f"\n⭐ Model Status:")
print(f"  scVI model path: {scvi_model_path}")
print(f"  scVI exists: {scvi_model_exists}")
print(f"  Will train scVI: {TRAIN_SCVI}")
print(f"\n  scANVI model path: {scanvi_model_path}")
print(f"  scANVI exists: {scanvi_model_exists}")
print(f"  Will train scANVI: {TRAIN_SCANVI}")

print(f"\n⭐ Hierarchical Labels: {USE_HIERARCHICAL_LABELS}")
if USE_HIERARCHICAL_LABELS:
    print(f"  Fine types: {len(MAJOR_LINEAGE_MAP)}")
    print(f"  Major lineages: {len(set(MAJOR_LINEAGE_MAP.values()))}")

print(f"\nBatch key: {BATCH_KEY}")
print(f"Min cells per batch: {MIN_CELLS_PER_BATCH} (⭐ INCREASED for stability)")
print(f"Low confidence threshold: {LOW_CONFIDENCE_THRESHOLD}")

print(f"\nscVI parameters:")
print(f"  HVG: {N_HVG}, Latent: {N_LATENT_SCVI}, Batch size: {BATCH_SIZE}")
print(f"  Dispersion: {DISPERSION}")

print(f"\nscANVI parameters:")
print(f"  Latent: {N_LATENT_SCANVI}, Max epochs: {SCANVI_MAX_EPOCHS}")
print("="*80)


# ===== STEP 1: Load and Validate Data =====

print("\n" + "="*80)
print("STEP 1: Loading and Validating Data")
print("="*80)

adata = sc.read_h5ad(INPUT_H5AD)

print(f"\n📊 Data Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# Ensure unique cell IDs
if not adata.obs_names.is_unique:
    print("\n⚠️  Duplicate cell IDs found. Making unique...")
    adata.obs_names_make_unique()
    print("✓ Cell IDs are now unique")

# Check required columns
required_columns = [CELLTYPE_KEY, 'celltypist_conf']
missing_columns = [col for col in required_columns if col not in adata.obs.columns]
if missing_columns:
    raise ValueError(f"Missing required columns: {missing_columns}")

# Check batch key
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found in adata.obs")

# Check for counts layer
if 'counts' not in adata.layers:
    print("\n⚠️  No 'counts' layer found")
    if sparse.issparse(adata.X):
        if adata.X.max() > 50 and not np.any(adata.X.data % 1):
            print("   X appears to be raw counts, using it")
            adata.layers['counts'] = adata.X.copy()
        else:
            raise ValueError("X does not appear to be raw counts. Need counts layer")
    else:
        raise ValueError("counts layer required for scVI training")
else:
    print("✓ counts layer found")

# ⭐ NEW: Validate data for invalid values
print("\n1.4 Checking for invalid values in counts layer...")

# Check for infinity
if sparse.issparse(adata.layers["counts"]):
    counts_data = adata.layers["counts"].data
else:
    counts_data = adata.layers["counts"]

if np.isinf(counts_data).any():
    print("⚠️  Found infinity values in counts layer!")
    if sparse.issparse(adata.layers["counts"]):
        adata.layers["counts"].data[np.isinf(adata.layers["counts"].data)] = 0
    else:
        adata.layers["counts"][np.isinf(adata.layers["counts"])] = 0
    print("✓ Replaced infinity with 0")

# Check for NaN
if np.isnan(counts_data).any():
    print("⚠️  Found NaN values in counts layer!")
    if sparse.issparse(adata.layers["counts"]):
        adata.layers["counts"].data[np.isnan(adata.layers["counts"].data)] = 0
    else:
        adata.layers["counts"][np.isnan(adata.layers["counts"])] = 0
    print("✓ Replaced NaN with 0")

print("✓ Data validation complete")


# ===== STEP 2: Quality Control & Filtering =====

print("\n" + "="*80)
print("STEP 2: Quality Control & Filtering")
print("="*80)

# Filter small batches
batch_counts = adata.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index

if len(small_batches) > 0:
    n_small_batches = len(small_batches)
    n_cells_removed = batch_counts[small_batches].sum()
    
    print(f"\n⚠️  Removing {n_small_batches} batches with <{MIN_CELLS_PER_BATCH} cells")
    print(f"    ({n_cells_removed:,} cells will be removed)")
    
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"✓ Remaining cells: {adata.n_obs:,}")
    
    # Show remaining batch distribution
    print(f"\n📊 Batch size distribution after filtering:")
    remaining_batch_counts = adata.obs[BATCH_KEY].value_counts()
    print(f"  Min: {remaining_batch_counts.min()}")
    print(f"  Max: {remaining_batch_counts.max()}")
    print(f"  Mean: {remaining_batch_counts.mean():.1f}")
    print(f"  Median: {remaining_batch_counts.median():.1f}")

print("\n✓ QC complete")


# ===== STEP 3: HVG Selection with Robust Error Handling =====

print("\n" + "="*80)
print("STEP 3: Highly Variable Gene Selection (ROBUST)")
print("="*80)

print(f"\n3.1 Selecting {N_HVG} highly variable genes with fallback strategy...")

# ⭐ CRITICAL FIX: Robust HVG selection with try-except fallback
hvg_method = None
hvg_flavor = None

try:
    # Attempt batch-aware HVG selection with seurat_v3
    print("   Attempting batch-aware HVG (seurat_v3)...")
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        flavor="seurat_v3",
        subset=False
    )
    hvg_method = "batch-aware"
    hvg_flavor = "seurat_v3"
    print("✓ Successfully used batch-aware HVG selection (seurat_v3)")
    
except Exception as e:
    print(f"⚠️  Batch-aware seurat_v3 failed: {str(e)[:100]}")
    print("   Falling back to batch-aware seurat...")
    
    try:
        # Fallback 1: batch-aware with seurat
        sc.pp.highly_variable_genes(
            adata,
            layer="counts",
            n_top_genes=N_HVG,
            batch_key=BATCH_KEY,
            flavor="seurat",
            subset=False
        )
        hvg_method = "batch-aware"
        hvg_flavor = "seurat"
        print("✓ Successfully used batch-aware HVG selection (seurat)")
        
    except Exception as e2:
        print(f"⚠️  Batch-aware seurat also failed: {str(e2)[:100]}")
        print("   Falling back to non-batch-aware HVG...")
        
        try:
            # Fallback 2: non-batch-aware with seurat_v3
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=N_HVG,
                flavor="seurat_v3",
                subset=False
            )
            hvg_method = "non-batch-aware"
            hvg_flavor = "seurat_v3"
            print("✓ Successfully used non-batch-aware HVG selection (seurat_v3)")
            
        except Exception as e3:
            print(f"⚠️  Non-batch-aware seurat_v3 also failed: {str(e3)[:100]}")
            print("   Final fallback to non-batch-aware seurat...")
            
            # Fallback 3: non-batch-aware with seurat (most conservative)
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=N_HVG,
                flavor="seurat",
                subset=False
            )
            hvg_method = "non-batch-aware"
            hvg_flavor = "seurat"
            print("✓ Used non-batch-aware HVG selection (seurat) as final fallback")

# Store HVG method info
adata.uns['hvg_method'] = hvg_method
adata.uns['hvg_flavor'] = hvg_flavor

hvg_standard = set(adata.var_names[adata.var['highly_variable']])
print(f"\n✓ Standard HVG: {len(hvg_standard)} genes")
print(f"  Method: {hvg_method}, Flavor: {hvg_flavor}")

# Force-include critical markers
print(f"\n3.2 Force-including {len(FORCE_INCLUDE_MARKERS)} marker genes...")
markers_in_data = [g for g in FORCE_INCLUDE_MARKERS if g in adata.var_names]
markers_missing = [g for g in FORCE_INCLUDE_MARKERS if g not in adata.var_names]

print(f"  Found: {len(markers_in_data)} markers")
if len(markers_missing) > 0:
    print(f"  Missing: {len(markers_missing)} markers")
    if len(markers_missing) <= 10:
        print(f"    {markers_missing}")
    else:
        print(f"    {markers_missing[:10]} ...")

# Union HVG + markers
hvg_final = hvg_standard | set(markers_in_data)
adata.var['highly_variable_final'] = adata.var_names.isin(hvg_final)

print(f"\n✓ Final HVG: {sum(adata.var['highly_variable_final'])} genes")
print(f"  Added {len(hvg_final) - len(hvg_standard)} marker genes")

# ⭐ CRITICAL: Save full gene matrix to .raw BEFORE subsetting
print("\n3.3 Preserving full gene matrix to .raw (shared memory)...")
adata.raw = sc.AnnData(
    X=adata.layers['counts'],  # ⭐ No .copy() - shared memory!
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ adata.raw: {adata.raw.n_obs:,} cells × {adata.raw.n_vars:,} genes")
print(f"  Memory: Shared with layers['counts'] (0 extra cost)")

# Subset to HVG for model training
print("\n3.4 Subsetting to HVG for efficient model training...")
adata = adata[:, adata.var['highly_variable_final']].copy()
print(f"✓ adata: {adata.n_obs:,} cells × {adata.n_vars:,} genes (HVG only)")

print("\n✓ HVG selection complete")
print(f"\n📊 Memory Strategy:")
print(f"  Training: {adata.n_vars} HVG (fast, memory-efficient)")
print(f"  Analysis: {adata.raw.n_vars} genes (complete, via .raw)")
print(f"  Extra cost: 0 GB (shared memory)")


# ===== STEP 4: scVI Model Training or Loading =====

print("\n" + "="*80)
print("STEP 4: scVI Model Training or Loading")
print("="*80)

if TRAIN_SCVI:
    print("\n⭐ Training new scVI model...")
    
    # Setup for scVI
    print("\n4.1 Setting up scVI...")
    scvi.model.SCVI.setup_anndata(
        adata,
        layer='counts',
        batch_key=BATCH_KEY
    )
    print(f"✓ AnnData registered with scVI")
    
    # Create scVI model
    print("\n4.2 Creating scVI model...")
    scvi_model = scvi.model.SCVI(
        adata,
        n_latent=N_LATENT_SCVI,
        n_layers=N_LAYERS,
        n_hidden=N_HIDDEN,
        dropout_rate=DROPOUT_RATE,
        gene_likelihood=GENE_LIKELIHOOD,
        dispersion=DISPERSION
    )
    
    # ⭐ Robust parameter counting (version-compatible)
    try:
        n_params = scvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
    
    print(f"✓ scVI model created")
    print(f"  Parameters: {n_params:,}")
    print(f"  Latent dimensions: {N_LATENT_SCVI}")
    print(f"  Hidden layers: {N_LAYERS}")
    print(f"  Hidden units: {N_HIDDEN}")
    
    # Train scVI
    print("\n4.3 Training scVI model...")
    print(f"  Max epochs: {SCVI_MAX_EPOCHS}")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Device: {DEVICE}")
    
    train_start = datetime.now()
    
    scvi_model.train(
        max_epochs=SCVI_MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        train_size=0.9,
        plan_kwargs={'lr': 1e-3}
    )
    
    train_end = datetime.now()
    train_duration = (train_end - train_start).total_seconds() / 60
    
    print(f"\n✓ scVI training completed in {train_duration:.1f} minutes")
    
    # Save scVI model
    print(f"\n4.4 Saving scVI model to: {scvi_model_path}")
    scvi_model.save(scvi_model_path, overwrite=True)
    print(f"✓ Model saved")

else:
    print(f"\n⭐ Loading existing scVI model from: {scvi_model_path}")
    
    try:
        # ⭐ Try with map_location first (newer versions)
        try:
            scvi_model = scvi.model.SCVI.load(
                scvi_model_path,
                adata=adata,
                map_location=DEVICE
            )
        except TypeError:
            # Fallback for versions that don't support map_location
            print("   (map_location not supported, using default loading)")
            scvi_model = scvi.model.SCVI.load(
                scvi_model_path,
                adata=adata
            )
        
        print(f"✓ scVI model loaded successfully")
        
        # ⭐ Robust parameter counting
        try:
            n_params = scvi_model.module.n_params
        except AttributeError:
            n_params = sum(p.numel() for p in scvi_model.module.parameters() if p.requires_grad)
        
        print(f"  Parameters: {n_params:,}")
        
    except Exception as e:
        raise RuntimeError(f"Failed to load scVI model. Error: {e}\n"
                         f"Set FORCE_RETRAIN_SCVI=True to retrain.")

# Extract scVI latent representation
print("\n4.5 Extracting scVI latent representation...")
adata.obsm['X_scvi'] = scvi_model.get_latent_representation()
print(f"✓ Latent representation: {adata.obsm['X_scvi'].shape}")

print("\n✓ Step 4 complete")

# Clean up to save memory
if not TRAIN_SCANVI:  # Only delete if we don't need it for scANVI
    print("\n4.6 Cleaning up scVI model from memory...")
    # Note: Keep scvi_model if we need it for scANVI training
    # del scvi_model
    # gc.collect()
else:
    print("\n4.6 Keeping scVI model in memory for scANVI initialization...")


# ===== STEP 5: Dimensionality Reduction =====

print("\n" + "="*80)
print("STEP 5: Dimensionality Reduction")
print("="*80)

print("\n5.1 Computing PCA on scVI latent space...")
sc.tl.pca(adata, use_highly_variable=False)  # Use all current genes (HVG)
print(f"✓ PCA computed: {adata.obsm['X_pca'].shape}")

print("\n5.2 Computing neighbors on scVI latent...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15)
print(f"✓ Neighbors computed")

print("\n5.3 Computing UMAP...")
sc.tl.umap(adata)
print(f"✓ UMAP computed: {adata.obsm['X_umap'].shape}")

print("\n✓ Step 5 complete")


# ===== STEP 6: Hierarchical Label Preparation =====

print("\n" + "="*80)
print("STEP 6: Hierarchical Label Preparation")
print("="*80)

if USE_HIERARCHICAL_LABELS:
    print("\n6.1 Creating hierarchical labels...")
    
    # Store original fine-grained types
    adata.obs['fine_type'] = adata.obs[CELLTYPE_KEY].copy()
    
    # Map to major lineages
    adata.obs['major_lineage'] = adata.obs['fine_type'].map(MAJOR_LINEAGE_MAP)
    
    # Handle unmapped types
    unmapped_mask = adata.obs['major_lineage'].isna()
    if unmapped_mask.any():
        print(f"\n⚠️  {unmapped_mask.sum()} cells have unmapped cell types:")
        unmapped_types = adata.obs.loc[unmapped_mask, 'fine_type'].unique()
        for utype in unmapped_types[:10]:
            count = (adata.obs['fine_type'] == utype).sum()
            print(f"    {utype}: {count} cells")
        if len(unmapped_types) > 10:
            print(f"    ... and {len(unmapped_types) - 10} more")
        
        # Set unmapped to Unknown
        adata.obs.loc[unmapped_mask, 'major_lineage'] = UNLABELED_CATEGORY
        print(f"✓ Set unmapped types to '{UNLABELED_CATEGORY}'")
    
    # Summary
    print(f"\n📊 Hierarchical Label Summary:")
    print(f"  Fine types: {adata.obs['fine_type'].nunique()}")
    print(f"  Major lineages: {adata.obs['major_lineage'].nunique()}")
    
    print(f"\n📋 Fine type → Major lineage mapping:")
    mapping_df = pd.DataFrame({
        'fine_type': adata.obs['fine_type'],
        'major_lineage': adata.obs['major_lineage']
    }).drop_duplicates()
    
    for _, row in mapping_df.sort_values(['major_lineage', 'fine_type']).iterrows():
        fine_type = row['fine_type']
        lineage = row['major_lineage']
        count = ((adata.obs['fine_type'] == fine_type) & 
                 (adata.obs['major_lineage'] == lineage)).sum()
        pct = 100 * count / adata.n_obs
        print(f"  {fine_type:25s} → {lineage:20s}: {count:7,} ({pct:5.2f}%)")
    
    print(f"\n✓ Using major lineages for scANVI training")

# Mark low-confidence cells as Unknown
print(f"\n6.2 Preparing scANVI labels...")

if USE_HIERARCHICAL_LABELS:
    adata.obs['scanvi_labels'] = adata.obs['major_lineage'].copy()
else:
    adata.obs['scanvi_labels'] = adata.obs[CELLTYPE_KEY].copy()

if MARK_LOW_CONFIDENCE_AS_UNKNOWN:
    low_conf_mask = adata.obs['celltypist_conf'] < LOW_CONFIDENCE_THRESHOLD
    n_low_conf = low_conf_mask.sum()
    adata.obs.loc[low_conf_mask, 'scanvi_labels'] = UNLABELED_CATEGORY
    print(f"✓ Marked {n_low_conf:,} low-confidence cells ({100*n_low_conf/adata.n_obs:.2f}%) as '{UNLABELED_CATEGORY}'")

print(f"\n📊 scANVI Label Summary:")
print(f"  Total: {adata.n_obs:,}")
labeled_cells = (adata.obs['scanvi_labels'] != UNLABELED_CATEGORY).sum()
print(f"  Labeled: {labeled_cells:,} ({100*labeled_cells/adata.n_obs:.2f}%)")
unknown_cells = (adata.obs['scanvi_labels'] == UNLABELED_CATEGORY).sum()
print(f"  Unknown: {unknown_cells:,} ({100*unknown_cells/adata.n_obs:.2f}%)")
print(f"  Categories: {adata.obs['scanvi_labels'].nunique()}")

print(f"\n📋 Label distribution:")
for label, count in adata.obs['scanvi_labels'].value_counts().items():
    pct = 100 * count / adata.n_obs
    print(f"  {label:25s}: {count:7,} ({pct:5.2f}%)")

print("\n✓ Labels prepared for scANVI")


# ===== STEP 7: scANVI Model Training or Loading =====

print("\n" + "="*80)
print("STEP 7: scANVI Model Training or Loading")
print("="*80)

if TRAIN_SCANVI:
    print("\n⭐ Training new scANVI model...")
    
    # ⭐ Check scVI-scANVI compatibility
    print("\n7.0 Checking scVI model compatibility...")
    n_labels_expected = adata.obs['scanvi_labels'].nunique()
    print(f"   Expected label categories: {n_labels_expected}")
    
    # Check if scVI model dimensions match current data
    try:
        scvi_n_vars = scvi_model.adata_manager.adata.n_vars
        current_n_vars = adata.n_vars
        
        print(f"   scVI model genes: {scvi_n_vars}")
        print(f"   Current data genes: {current_n_vars}")
        
        if scvi_n_vars != current_n_vars:
            print(f"\n⚠️  WARNING: Gene count mismatch detected!")
            print(f"   This usually means HVG selection changed.")
            print(f"   scANVI initialization may fail.")
            print(f"\n💡 Recommendation: Set FORCE_RETRAIN_SCVI=True and rerun")
    except Exception as e:
        print(f"   Could not verify compatibility: {e}")
    
    print("\n7.1 Creating scANVI model from scVI...")
    
    try:
        # Try standard initialization (scvi-tools v1.0+)
        print("   Attempting standard initialization...")
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            adata=adata,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key='scanvi_labels'
        )
        print("✓ scANVI initialized successfully")
        
    except (ValueError, RuntimeError) as e:
        error_msg = str(e)
        print(f"⚠️  Standard init failed: {error_msg[:150]}")
        
        # Check if it's a size mismatch error (model incompatibility)
        if "size mismatch" in error_msg and ("px_r" in error_msg or "decoder" in error_msg):
            print("\n" + "="*80)
            print("❌ CRITICAL ERROR: Model Architecture Mismatch")
            print("="*80)
            print("\n🔍 Root Cause:")
            print("   The existing scVI model was trained with different data:")
            print("   - Different number of HVG genes")
            print("   - Different HVG selection method")
            print("   - Possibly different label categories")
            print("\n💡 Solution:")
            print("   Set the following parameters at the top of the script:")
            print("   FORCE_RETRAIN_SCVI = True")
            print("   FORCE_RETRAIN_SCANVI = True")
            print("\n📝 Then rerun:")
            print("   nohup python epithelial_scvi_scanvi_20260114_v2_4_FIXED.py > epi_retrain.log 2>&1 &")
            print("\n⏱️  Expected time: ~90-120 minutes (scVI + scANVI)")
            print("="*80)
            raise RuntimeError("scVI-scANVI model compatibility issue. See instructions above.")
        
        print("   Attempting explicit setup fallback...")
        
        # Fallback: explicit setup
        scvi.model.SCANVI.setup_anndata(
            adata,
            layer='counts',
            batch_key=BATCH_KEY,
            labels_key='scanvi_labels',
            unlabeled_category=UNLABELED_CATEGORY
        )
        
        try:
            scanvi_model = scvi.model.SCANVI.from_scvi_model(
                scvi_model,
                adata=adata,
                unlabeled_category=UNLABELED_CATEGORY,
                labels_key='scanvi_labels'
            )
            print("✓ scANVI initialized with fallback method")
        except RuntimeError as e2:
            # If fallback also fails with size mismatch, provide clear instructions
            if "size mismatch" in str(e2):
                print("\n" + "="*80)
                print("❌ CRITICAL ERROR: Model Architecture Mismatch (Confirmed)")
                print("="*80)
                print("\n🔍 Both initialization methods failed due to incompatible model dimensions.")
                print("\n💡 REQUIRED ACTION:")
                print("   1. Edit the script and set:")
                print("      FORCE_RETRAIN_SCVI = True")
                print("      FORCE_RETRAIN_SCANVI = True")
                print("\n   2. Rerun the script:")
                print("      nohup python epithelial_scvi_scanvi_20260114_v2_4_FIXED.py > epi_retrain.log 2>&1 &")
                print("\n⏱️  Total training time: ~90-120 minutes")
                print("="*80)
                raise
            else:
                raise
    
    # ⭐ Robust parameter counting
    try:
        n_params = scanvi_model.module.n_params
    except AttributeError:
        n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)
    
    print(f"\n✓ scANVI model created")
    print(f"  Parameters: {n_params:,}")
    print(f"  Unlabeled category: {UNLABELED_CATEGORY}")
    
    # Train scANVI
    print("\n7.2 Training scANVI model...")
    print(f"  Max epochs: {SCANVI_MAX_EPOCHS}")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Device: {DEVICE}")
    
    train_start = datetime.now()
    
    scanvi_model.train(
        max_epochs=SCANVI_MAX_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOPPING,
        train_size=0.9,
        plan_kwargs={'lr': 1e-3}
    )
    
    train_end = datetime.now()
    train_duration = (train_end - train_start).total_seconds() / 60
    
    print(f"\n✓ scANVI training completed in {train_duration:.1f} minutes")
    
    # Save scANVI model
    print(f"\n7.3 Saving scANVI model to: {scanvi_model_path}")
    scanvi_model.save(scanvi_model_path, overwrite=True)
    print(f"✓ Model saved")

else:
    print(f"\n⭐ Loading existing scANVI model from: {scanvi_model_path}")
    
    try:
        # ⭐ Try with map_location first (newer versions)
        try:
            scanvi_model = scvi.model.SCANVI.load(
                scanvi_model_path,
                adata=adata,
                map_location=DEVICE
            )
        except TypeError:
            # Fallback for versions that don't support map_location
            print("   (map_location not supported, using default loading)")
            scanvi_model = scvi.model.SCANVI.load(
                scanvi_model_path,
                adata=adata
            )
        
        print(f"✓ scANVI model loaded successfully")
        
        # ⭐ Robust parameter counting
        try:
            n_params = scanvi_model.module.n_params
        except AttributeError:
            n_params = sum(p.numel() for p in scanvi_model.module.parameters() if p.requires_grad)
        
        print(f"  Parameters: {n_params:,}")
        
    except Exception as e:
        raise RuntimeError(f"Failed to load scANVI model. Error: {e}\n"
                         f"Set FORCE_RETRAIN_SCANVI=True to retrain.")

# Extract predictions
print("\n7.4 Extracting scANVI predictions...")

adata.obsm['X_scanvi'] = scanvi_model.get_latent_representation()
print(f"✓ Latent representation: {adata.obsm['X_scanvi'].shape}")

# Get predictions
predictions = scanvi_model.predict()
adata.obs['scanvi_predictions'] = predictions
print(f"✓ Predictions extracted")

# Get prediction probabilities (confidence)
pred_probs = scanvi_model.predict(soft=True)
if isinstance(pred_probs, np.ndarray):
    adata.obs['scanvi_confidence'] = pred_probs.max(axis=1)
else:
    adata.obs['scanvi_confidence'] = np.asarray(pred_probs).max(axis=1)
print(f"✓ Confidence scores calculated")

# Compute neighbors and UMAP on scANVI latent
print("\n7.5 Computing scANVI-based neighbors and UMAP...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)
sc.tl.umap(adata)

print("\n✓ Step 7 complete")

# Clean up models from memory
print("\n7.6 Cleaning up models from memory...")
del scvi_model, scanvi_model
gc.collect()
print("✓ Memory cleanup complete")


# ===== STEP 8: Final Results & Visualization =====

print("\n" + "="*80)
print("STEP 8: Final Results & Visualization")
print("="*80)

# Generate summary statistics
print("\n8.1 Generating summary statistics...")

stats_df = pd.DataFrame({
    'Cell_Type': adata.obs['scanvi_predictions'].value_counts().index,
    'Count': adata.obs['scanvi_predictions'].value_counts().values,
    'Percentage': 100 * adata.obs['scanvi_predictions'].value_counts(normalize=True).values
})

# Add confidence statistics
conf_by_type = adata.obs.groupby('scanvi_predictions')['scanvi_confidence'].agg(['mean', 'std', 'min', 'max'])
stats_df = stats_df.merge(conf_by_type, left_on='Cell_Type', right_index=True, how='left')
stats_df.columns = ['Cell_Type', 'Count', 'Percentage', 'Mean_Confidence', 'Std_Confidence', 'Min_Confidence', 'Max_Confidence']

stats_df = stats_df.sort_values('Count', ascending=False)
stats_df.to_csv(f"{OUTPUT_DIR}/scanvi_cell_type_statistics.csv", index=False)
print("✓ Saved: scanvi_cell_type_statistics.csv")

# Create visualizations
print("\n8.2 Creating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

if USE_HIERARCHICAL_LABELS:
    sc.pl.umap(adata, color='major_lineage',
               ax=axes[0, 0],
               frameon=False,
               title='Major Lineages (scANVI Training Labels)',
               legend_loc='right margin',
               s=5,
               show=False)
    
    sc.pl.umap(adata, color='fine_type',
               ax=axes[0, 1],
               frameon=False,
               title='Fine Cell Types (Original CellTypist)',
               legend_loc='right margin',
               legend_fontsize=8,
               s=5,
               show=False)
else:
    sc.pl.umap(adata, color=CELLTYPE_KEY,
               ax=axes[0, 0],
               frameon=False,
               title='Original CellTypist Annotations',
               legend_loc='right margin',
               s=5,
               show=False)
    
    sc.pl.umap(adata, color=BATCH_KEY,
               ax=axes[0, 1],
               frameon=False,
               title='Sample Batch',
               legend_loc='right margin',
               s=5,
               show=False)

sc.pl.umap(adata, color='scanvi_predictions',
           ax=axes[1, 0],
           frameon=False,
           title='scANVI Predictions',
           legend_loc='right margin',
           s=5,
           show=False)

sc.pl.umap(adata, color='scanvi_confidence',
           ax=axes[1, 1],
           frameon=False,
           title='scANVI Prediction Confidence',
           cmap='viridis',
           vmin=0, vmax=1,
           s=5,
           show=False)

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/figures/scanvi_results_overview.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Saved: scanvi_results_overview.png")

# Save final annotated dataset
print("\n8.3 Saving final annotated dataset...")
final_file = f"{OUTPUT_DIR}/epithelial_scanvi_final.h5ad"

# Ensure data structure completeness
print("\n📊 Final data structure:")
print(f"  adata.X shape: {adata.X.shape}")
print(f"  adata.raw.X shape: {adata.raw.X.shape}")
print(f"  Layers: {list(adata.layers.keys())}")
print(f"  Obsm: {list(adata.obsm.keys())}")

adata.write_h5ad(final_file, compression='gzip')
print(f"✓ Saved: {final_file}")

print("\n" + "="*80)
print("✓ PIPELINE COMPLETE")
print("="*80)
print(f"Finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

print(f"\n🎯 Key outputs:")
print(f"  1. scVI model: {scvi_model_path}")
print(f"  2. scANVI model: {scanvi_model_path}")
print(f"  3. Final dataset: {final_file}")
print(f"  4. Statistics: {OUTPUT_DIR}/scanvi_cell_type_statistics.csv")
print(f"  5. Figures: {OUTPUT_DIR}/figures/")

print(f"\n📊 HVG Strategy Summary:")
print(f"  Method: {adata.uns['hvg_method']}")
print(f"  Flavor: {adata.uns['hvg_flavor']}")
print(f"  Training genes: {adata.n_vars}")
print(f"  Analysis genes: {adata.raw.n_vars}")

print("\n⭐ Key improvements in v2.4-FIXED:")
print("  1. Robust HVG selection with multi-level fallback")
print("  2. Data validation for infinity/NaN values")
print("  3. Increased MIN_CELLS_PER_BATCH to 30 for stability")
print("  4. Version-compatible parameter counting")
print("  5. Shared memory for .raw (0 extra cost)")
print("  6. Complete data structure preservation")

print("="*80)
