#!/usr/bin/env python3
"""
All Cells: Integrate scANVI Annotations → scVI → scANVI Pipeline v1.1
=======================================================================

Purpose:
--------
1. Extract scANVI + CellTypist annotations from cell-type-specific analyses
2. Merge annotations back to denormalized all-cells data
3. Remove Ji_Hoon_Ahn_2021 dataset
4. Filter labels (low confidence + rare cell types < 200 cells)
5. Train scVI for batch correction (memory-optimized)
6. Train scANVI with filtered labels (semi-supervised)
7. Generate comparative UMAP visualizations

Input Files:
------------
- Main data: /home/h2048/data/py/20251215/denormalize_counts_test/adata_denormalized.h5ad
- B cell: /home/h2048/data/py/1208/b_scvi_celltypist_scanvi_v1_1/adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad
- Myeloid: /home/h2048/data/py/1208/myeloid_scvi_celltypist_scanvi_v1_1/adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad
- T cell: /home/h2048/data/py/1206/tcell_scvi_celltypist_scanvi_v2/adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad
- Stromal: /home/h2048/data/py/1209/stromal_vascular_scvi_analysis/adata_scvi_scanvi_final.h5ad
- Epithelial: /home/h2048/data/py/1214/epithelial_scanvi_v3_5_1_hotfix/adata_epithelial_FINAL.h5ad

Output:
-------
- Integrated h5ad with scVI and scANVI embeddings
- UMAP visualizations
- Analysis summary

Environment:
-----------
conda activate scvi

Author: Claude Code
Date: 2025-12-16
Version: 1.1 (Memory Optimized + CellTypist)
"""

import sys
import os
from pathlib import Path
from datetime import datetime
import warnings
import json
import gc  # For memory management

import numpy as np
import pandas as pd
from scipy.sparse import issparse

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to save memory
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import scvi

warnings.filterwarnings('ignore')

# =============================================================================
# Configuration
# =============================================================================

# Input/Output paths
INPUT_MAIN = "/home/h2048/data/py/20251215/denormalize_counts_test/adata_denormalized.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%Y%m%d')}/allcells_integrated_scanvi")

# Cell-type-specific scANVI result files
SCANVI_RESULTS = {
    'B_cell': {
        'path': '/home/h2048/data/py/1208/b_scvi_celltypist_scanvi_v1_1/adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_annotation_col': 'scanvi_predictions',
        'scanvi_confidence_col': 'scanvi_confidence',
        'celltypist_cols': ['predicted_labels', 'majority_voting', 'conf_score']
    },
    'Myeloid': {
        'path': '/home/h2048/data/py/1208/myeloid_scvi_celltypist_scanvi_v1_1/adata_myeloid_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_annotation_col': 'cell_type_scanvi_filt',
        'scanvi_confidence_col': 'scanvi_confidence',
        'celltypist_cols': ['predicted_labels', 'majority_voting', 'conf_score']
    },
    'T_cell': {
        'path': '/home/h2048/data/py/1206/tcell_scvi_celltypist_scanvi_v2/adata_tcell_scvi_celltypist_scanvi_final_v2.h5ad',
        'scanvi_annotation_col': 'scanvi_predictions',
        'scanvi_confidence_col': 'scanvi_confidence',
        'celltypist_cols': ['predicted_labels', 'majority_voting', 'conf_score']
    },
    'Stromal': {
        'path': '/home/h2048/data/py/1209/stromal_vascular_scvi_analysis/adata_scvi_scanvi_final.h5ad',
        'scanvi_annotation_col': 'cell_type_scanvi_filt',
        'scanvi_confidence_col': 'scanvi_confidence',
        'celltypist_cols': []  # Stromal may not have CellTypist
    },
    'Epithelial': {
        'path': '/home/h2048/data/py/1214/epithelial_scanvi_v3_5_1_hotfix/adata_epithelial_FINAL.h5ad',
        'scanvi_annotation_col': 'cell_type_scanvi_filt',
        'scanvi_confidence_col': 'scanvi_confidence',
        'celltypist_cols': ['celltypist_predicted_labels', 'celltypist_majority_voting', 'conf_score']
    }
}

# Unified column names in output
UNIFIED_SCANVI_ANNOTATION = 'scanvi_annotation'
UNIFIED_SCANVI_CONFIDENCE = 'scanvi_confidence'
UNIFIED_CELLTYPIST_PREDICTED = 'celltypist_predicted_labels'
UNIFIED_CELLTYPIST_MAJORITY = 'celltypist_majority_voting'
UNIFIED_CELLTYPIST_CONFIDENCE = 'celltypist_conf_score'
UNIFIED_LINEAGE_COL = 'cell_lineage'

# Dataset to remove
DATASET_TO_REMOVE = 'Ji_Hoon_Ahn_2021'

# Filtering criteria for training labels
MIN_CONFIDENCE = 0.5  # Cells below this confidence → unlabeled
MIN_CELLS_PER_TYPE = 200  # Cell types with fewer cells → unlabeled
UNLABELED_CATEGORY = "Unknown"

# scVI/scANVI parameters
BATCH_KEY = 'dataset'
N_LATENT = 30
N_LAYERS = 2
SCVI_MAX_EPOCHS = 400
SCANVI_MAX_EPOCHS = 200

# UMAP parameters
UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.3
UMAP_SPREAD = 1.0

# Visualization
FIGURE_FORMAT = 'pdf'
DPI = 300

# Memory optimization
ENABLE_MEMORY_OPTIMIZATION = True
CHECKPOINT_COMPRESSION = 'gzip'

# Random seed
RANDOM_SEED = 42
scvi.settings.seed = RANDOM_SEED
np.random.seed(RANDOM_SEED)

# =============================================================================
# Setup
# =============================================================================

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
checkpoint_dir = OUTPUT_DIR / "checkpoints"
checkpoint_dir.mkdir(exist_ok=True)
fig_dir = OUTPUT_DIR / "figures"
fig_dir.mkdir(exist_ok=True)
model_dir = OUTPUT_DIR / "models"
model_dir.mkdir(exist_ok=True)

LOG_FILE = OUTPUT_DIR / f"integration_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

def log_print(msg):
    """Print and log simultaneously"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {msg}"
    print(formatted_msg)
    with open(LOG_FILE, 'a') as f:
        f.write(formatted_msg + '\n')

def memory_usage_mb():
    """Get current memory usage in MB"""
    import psutil
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024

def log_memory(label=""):
    """Log current memory usage"""
    mem_mb = memory_usage_mb()
    log_print(f"  [Memory] {label}: {mem_mb:.1f} MB ({mem_mb/1024:.2f} GB)")

# =============================================================================
# Step 1: Extract Annotations (Memory-Optimized)
# =============================================================================

log_print("="*80)
log_print("STEP 1: EXTRACT ANNOTATIONS FROM CELL-TYPE-SPECIFIC ANALYSES")
log_print("="*80)
log_memory("Initial")

annotations_dict = {}

for lineage, config in SCANVI_RESULTS.items():
    log_print(f"\n{lineage}:")
    log_print(f"  File: {config['path']}")

    # Check if file exists
    if not Path(config['path']).exists():
        log_print(f"  ⚠️  WARNING: File not found, skipping this lineage")
        continue

    try:
        # Read only obs (not full h5ad) to save memory
        log_print(f"  Loading annotations only (minimal read)...")
        adata_subset = sc.read_h5ad(config['path'], backed='r')  # Backed mode for memory efficiency

        log_print(f"  ✓ Loaded: {adata_subset.shape[0]} cells × {adata_subset.shape[1]} genes")

        # Extract scANVI annotations
        scanvi_col = config['scanvi_annotation_col']
        scanvi_conf_col = config['scanvi_confidence_col']

        if scanvi_col not in adata_subset.obs.columns:
            log_print(f"  ✗ ERROR: Column '{scanvi_col}' not found!")
            log_print(f"    Available columns: {list(adata_subset.obs.columns[:10])}")
            del adata_subset
            gc.collect()
            continue

        # Initialize DataFrame with scANVI annotations
        annotations = pd.DataFrame(index=adata_subset.obs.index)
        annotations[UNIFIED_SCANVI_ANNOTATION] = adata_subset.obs[scanvi_col].astype(str)

        # scANVI confidence
        if scanvi_conf_col in adata_subset.obs.columns:
            annotations[UNIFIED_SCANVI_CONFIDENCE] = adata_subset.obs[scanvi_conf_col]
        else:
            log_print(f"  ⚠️  WARNING: '{scanvi_conf_col}' not found, using default 1.0")
            annotations[UNIFIED_SCANVI_CONFIDENCE] = 1.0

        # Extract CellTypist annotations (if available)
        celltypist_cols = config.get('celltypist_cols', [])

        for ct_col in celltypist_cols:
            if ct_col in adata_subset.obs.columns:
                # Map to unified names
                if 'predicted_labels' in ct_col:
                    unified_name = UNIFIED_CELLTYPIST_PREDICTED
                elif 'majority_voting' in ct_col:
                    unified_name = UNIFIED_CELLTYPIST_MAJORITY
                elif 'conf_score' in ct_col:
                    unified_name = UNIFIED_CELLTYPIST_CONFIDENCE
                else:
                    unified_name = f"celltypist_{ct_col}"

                annotations[unified_name] = adata_subset.obs[ct_col]
                log_print(f"  ✓ Extracted CellTypist: {ct_col} → {unified_name}")

        # Add lineage label
        annotations[UNIFIED_LINEAGE_COL] = lineage

        # Store in dict
        annotations_dict[lineage] = annotations

        # Summary
        n_unique_types = annotations[UNIFIED_SCANVI_ANNOTATION].nunique()
        mean_conf = annotations[UNIFIED_SCANVI_CONFIDENCE].mean()
        log_print(f"  ✓ Extracted {len(annotations)} cells")
        log_print(f"    scANVI cell types: {n_unique_types}")
        log_print(f"    scANVI mean confidence: {mean_conf:.3f}")

        # Check CellTypist presence
        if UNIFIED_CELLTYPIST_PREDICTED in annotations.columns:
            n_ct_types = annotations[UNIFIED_CELLTYPIST_PREDICTED].nunique()
            log_print(f"    CellTypist types: {n_ct_types}")

        # Clean up to free memory
        del adata_subset
        gc.collect()
        log_memory(f"After {lineage}")

    except Exception as e:
        log_print(f"  ✗ ERROR: {e}")
        import traceback
        log_print(traceback.format_exc())
        gc.collect()

# Combine all annotations
log_print("\n" + "-"*80)
log_print("Combining annotations from all lineages...")
log_print("-"*80)

if len(annotations_dict) == 0:
    log_print("FATAL ERROR: No annotations were extracted!")
    sys.exit(1)

all_annotations = pd.concat(annotations_dict.values(), axis=0)
log_print(f"✓ Combined annotations: {len(all_annotations)} cells total")
log_print(f"  Unique scANVI types: {all_annotations[UNIFIED_SCANVI_ANNOTATION].nunique()}")
log_print(f"  Lineage distribution:")
for lineage, count in all_annotations[UNIFIED_LINEAGE_COL].value_counts().items():
    log_print(f"    {lineage}: {count} cells")

# Check for duplicate barcodes
duplicates = all_annotations.index.duplicated()
if duplicates.any():
    n_dup = duplicates.sum()
    log_print(f"\n⚠️  WARNING: Found {n_dup} duplicate cell barcodes!")
    log_print(f"  Keeping first occurrence only...")
    all_annotations = all_annotations[~duplicates]

# Save extracted annotations
annotations_file = OUTPUT_DIR / "extracted_annotations.csv"
all_annotations.to_csv(annotations_file)
log_print(f"\n✓ Saved extracted annotations to: {annotations_file.name}")

# Clear intermediate data
del annotations_dict
gc.collect()
log_memory("After combining annotations")

# =============================================================================
# Step 2: Load Main Data and Merge Annotations
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 2: LOAD MAIN DATA AND MERGE ANNOTATIONS")
log_print("="*80)

log_print(f"\nLoading main data: {INPUT_MAIN}")
log_memory("Before loading main data")

adata = sc.read_h5ad(INPUT_MAIN)

log_print(f"✓ Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
log_print(f"  Layers: {list(adata.layers.keys())}")
log_memory("After loading main data")

# Check for batch key
if BATCH_KEY not in adata.obs.columns:
    log_print(f"\n✗ ERROR: Batch key '{BATCH_KEY}' not found in obs!")
    log_print(f"  Available columns: {list(adata.obs.columns)}")
    sys.exit(1)

log_print(f"\nDataset distribution (before filtering):")
for dataset, count in adata.obs[BATCH_KEY].value_counts().items():
    log_print(f"  {dataset}: {count} cells")

# =============================================================================
# Step 3: Remove Ji_Hoon_Ahn_2021 Dataset
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 3: REMOVE SPECIFIC DATASET")
log_print("="*80)

cells_before = adata.shape[0]
if DATASET_TO_REMOVE in adata.obs[BATCH_KEY].values:
    n_remove = (adata.obs[BATCH_KEY] == DATASET_TO_REMOVE).sum()
    log_print(f"\nRemoving dataset: {DATASET_TO_REMOVE}")
    log_print(f"  Cells to remove: {n_remove}")

    adata = adata[adata.obs[BATCH_KEY] != DATASET_TO_REMOVE].copy()
    gc.collect()

    cells_after = adata.shape[0]
    log_print(f"✓ Removed {cells_before - cells_after} cells")
    log_print(f"  Remaining: {cells_after} cells")
    log_memory("After dataset removal")
else:
    log_print(f"\n✓ Dataset '{DATASET_TO_REMOVE}' not found, no removal needed")

log_print(f"\nDataset distribution (after filtering):")
for dataset, count in adata.obs[BATCH_KEY].value_counts().items():
    log_print(f"  {dataset}: {count} cells")

# =============================================================================
# Step 4: Merge Annotations by Barcode
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 4: MERGE ANNOTATIONS BY CELL BARCODE")
log_print("="*80)

# Merge annotations based on cell barcode (index)
log_print("\nMatching cell barcodes...")
matched_cells = adata.obs.index.isin(all_annotations.index)
n_matched = matched_cells.sum()
n_total = len(adata)

log_print(f"  Total cells in main data: {n_total}")
log_print(f"  Cells with annotations: {n_matched} ({n_matched/n_total*100:.1f}%)")
log_print(f"  Cells without annotations: {n_total - n_matched} ({(n_total-n_matched)/n_total*100:.1f}%)")

# Initialize all annotation columns with default values
adata.obs[UNIFIED_SCANVI_ANNOTATION] = UNLABELED_CATEGORY
adata.obs[UNIFIED_SCANVI_CONFIDENCE] = 0.0
adata.obs[UNIFIED_LINEAGE_COL] = 'Unknown'

# Initialize CellTypist columns (if present in annotations)
if UNIFIED_CELLTYPIST_PREDICTED in all_annotations.columns:
    adata.obs[UNIFIED_CELLTYPIST_PREDICTED] = 'Unknown'
if UNIFIED_CELLTYPIST_MAJORITY in all_annotations.columns:
    adata.obs[UNIFIED_CELLTYPIST_MAJORITY] = 'Unknown'
if UNIFIED_CELLTYPIST_CONFIDENCE in all_annotations.columns:
    adata.obs[UNIFIED_CELLTYPIST_CONFIDENCE] = 0.0

# Merge annotations using reindex (handles missing values)
matched_annotations = all_annotations.reindex(adata.obs.index)
valid_mask = ~matched_annotations[UNIFIED_SCANVI_ANNOTATION].isna()

# Merge scANVI columns
adata.obs.loc[valid_mask, UNIFIED_SCANVI_ANNOTATION] = matched_annotations.loc[valid_mask, UNIFIED_SCANVI_ANNOTATION]
adata.obs.loc[valid_mask, UNIFIED_SCANVI_CONFIDENCE] = matched_annotations.loc[valid_mask, UNIFIED_SCANVI_CONFIDENCE]
adata.obs.loc[valid_mask, UNIFIED_LINEAGE_COL] = matched_annotations.loc[valid_mask, UNIFIED_LINEAGE_COL]

# Merge CellTypist columns (if available)
for col in [UNIFIED_CELLTYPIST_PREDICTED, UNIFIED_CELLTYPIST_MAJORITY, UNIFIED_CELLTYPIST_CONFIDENCE]:
    if col in matched_annotations.columns:
        valid_col_mask = ~matched_annotations[col].isna()
        adata.obs.loc[valid_col_mask, col] = matched_annotations.loc[valid_col_mask, col]
        log_print(f"  ✓ Merged CellTypist column: {col}")

log_print(f"\n✓ Merged annotations")
log_print(f"  Unique scANVI types (including Unknown): {adata.obs[UNIFIED_SCANVI_ANNOTATION].nunique()}")
log_print(f"  Mean scANVI confidence (excluding Unknown): {adata.obs[adata.obs[UNIFIED_SCANVI_ANNOTATION] != UNLABELED_CATEGORY][UNIFIED_SCANVI_CONFIDENCE].mean():.3f}")

if UNIFIED_CELLTYPIST_PREDICTED in adata.obs.columns:
    n_ct = (adata.obs[UNIFIED_CELLTYPIST_PREDICTED] != 'Unknown').sum()
    log_print(f"  Cells with CellTypist labels: {n_ct}")

# Clean up
del all_annotations, matched_annotations
gc.collect()
log_memory("After merging annotations")

# Save checkpoint
checkpoint_file = checkpoint_dir / "checkpoint_01_annotations_merged.h5ad"
log_print(f"\nSaving checkpoint: {checkpoint_file.name}")
adata.write_h5ad(checkpoint_file, compression=CHECKPOINT_COMPRESSION)
log_print(f"✓ Checkpoint saved")
log_memory("After checkpoint save")

# =============================================================================
# Step 5: Filter Labels for Training
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 5: FILTER LABELS FOR scANVI TRAINING")
log_print("="*80)

# Create a copy of annotations for training labels
adata.obs['labels_for_scanvi'] = adata.obs[UNIFIED_SCANVI_ANNOTATION].copy()

# Count cells before filtering
labeled_before = (adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY).sum()
log_print(f"\nLabeled cells before filtering: {labeled_before}")

# Filter 1: Low confidence
log_print(f"\n[Filter 1] Removing low confidence predictions (< {MIN_CONFIDENCE})...")
low_conf_mask = (adata.obs[UNIFIED_SCANVI_CONFIDENCE] < MIN_CONFIDENCE) & \
                (adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY)
n_low_conf = low_conf_mask.sum()

if n_low_conf > 0:
    adata.obs.loc[low_conf_mask, 'labels_for_scanvi'] = UNLABELED_CATEGORY
    log_print(f"  Removed {n_low_conf} low confidence labels")
else:
    log_print(f"  No low confidence labels to remove")

# Filter 2: Rare cell types (< MIN_CELLS_PER_TYPE)
log_print(f"\n[Filter 2] Removing rare cell types (< {MIN_CELLS_PER_TYPE} cells)...")
labeled_cells = adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY
celltype_counts = adata.obs.loc[labeled_cells, 'labels_for_scanvi'].value_counts()
rare_types = celltype_counts[celltype_counts < MIN_CELLS_PER_TYPE].index.tolist()

if len(rare_types) > 0:
    log_print(f"  Found {len(rare_types)} rare cell types:")
    for ct in rare_types[:10]:  # Show first 10
        log_print(f"    - {ct}: {celltype_counts[ct]} cells")
    if len(rare_types) > 10:
        log_print(f"    ... and {len(rare_types) - 10} more")

    rare_mask = adata.obs['labels_for_scanvi'].isin(rare_types)
    n_rare = rare_mask.sum()
    adata.obs.loc[rare_mask, 'labels_for_scanvi'] = UNLABELED_CATEGORY
    log_print(f"  Removed {n_rare} cells from rare types")
else:
    log_print(f"  No rare cell types found")

# Summary after filtering
labeled_after = (adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY).sum()
unlabeled_after = (adata.obs['labels_for_scanvi'] == UNLABELED_CATEGORY).sum()

log_print(f"\n" + "="*60)
log_print("Label Filtering Summary:")
log_print("="*60)
log_print(f"  Before filtering: {labeled_before} labeled cells")
log_print(f"  After filtering:  {labeled_after} labeled cells")
log_print(f"  Removed: {labeled_before - labeled_after} labels")
log_print(f"  Unlabeled cells (for semi-supervised learning): {unlabeled_after}")
log_print(f"  Final unique cell types for training: {adata.obs[adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY]['labels_for_scanvi'].nunique()}")

# Top cell types for training
log_print(f"\nTop 10 cell types for training:")
for ct, count in adata.obs[adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY]['labels_for_scanvi'].value_counts().head(10).items():
    log_print(f"  {ct}: {count} cells")

# Save checkpoint
checkpoint_file = checkpoint_dir / "checkpoint_02_labels_filtered.h5ad"
log_print(f"\nSaving checkpoint: {checkpoint_file.name}")
adata.write_h5ad(checkpoint_file, compression=CHECKPOINT_COMPRESSION)
log_print(f"✓ Checkpoint saved")
log_memory("After label filtering")

# =============================================================================
# Step 6: Data Preparation for scVI
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 6: PREPARE DATA FOR scVI")
log_print("="*80)

# Check counts layer
if 'counts' not in adata.layers:
    log_print("\n✗ ERROR: 'counts' layer not found!")
    log_print(f"  Available layers: {list(adata.layers.keys())}")
    sys.exit(1)

log_print(f"\n✓ Using layers['counts'] for scVI training")
log_print(f"  Counts shape: {adata.layers['counts'].shape}")
log_print(f"  Counts dtype: {adata.layers['counts'].dtype}")

# Check for raw counts (should be integers after denormalization)
counts_sample = adata.layers['counts'][:100, :100]
if issparse(counts_sample):
    counts_sample = counts_sample.toarray()
is_integer = np.allclose(counts_sample, np.round(counts_sample))
log_print(f"  Integer-like values: {is_integer}")

if not is_integer:
    log_print(f"  ⚠️  WARNING: Counts may not be raw integers!")

# Store counts in .X temporarily for scVI
log_print(f"\nPreparing data for scVI...")
adata.layers['original_X'] = adata.X.copy()  # Backup original X
adata.X = adata.layers['counts'].copy()

# Filter genes (keep genes expressed in at least 10 cells)
log_print(f"\nFiltering genes...")
genes_before = adata.shape[1]
sc.pp.filter_genes(adata, min_cells=10)
genes_after = adata.shape[1]
log_print(f"  Genes before: {genes_before}")
log_print(f"  Genes after: {genes_after}")
log_print(f"  Removed: {genes_before - genes_after} genes")

log_print(f"\n✓ Final data shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
log_memory("After data preparation")

# Save checkpoint
checkpoint_file = checkpoint_dir / "checkpoint_03_data_prepared.h5ad"
log_print(f"\nSaving checkpoint: {checkpoint_file.name}")
adata.write_h5ad(checkpoint_file, compression=CHECKPOINT_COMPRESSION)
log_print(f"✓ Checkpoint saved")

# =============================================================================
# Step 7: Train scVI Model (Memory-Optimized)
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 7: TRAIN scVI MODEL (BATCH CORRECTION)")
log_print("="*80)

log_print(f"\nSetting up scVI model...")
log_print(f"  Batch key: {BATCH_KEY}")
log_print(f"  N latent: {N_LATENT}")
log_print(f"  N layers: {N_LAYERS}")
log_print(f"  Max epochs: {SCVI_MAX_EPOCHS}")

scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    batch_key=BATCH_KEY
)

vae = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    gene_likelihood='nb'
)

log_print(f"\n✓ scVI model created")
log_memory("After scVI model creation")

log_print(f"\nTraining scVI model...")
import time
scvi_start = time.time()

vae.train(
    max_epochs=SCVI_MAX_EPOCHS,
    early_stopping=True,
    early_stopping_patience=20,
    batch_size=128  # Smaller batch size for memory
)

scvi_time = time.time() - scvi_start

log_print(f"\n✓ scVI training completed")
log_print(f"  Time: {scvi_time:.1f}s ({scvi_time/60:.1f} min)")
log_memory("After scVI training")

# Save scVI model
scvi_model_dir = model_dir / "scvi_model"
vae.save(scvi_model_dir, overwrite=True)
log_print(f"\n✓ scVI model saved to: {scvi_model_dir}")

# Extract scVI latent representation
log_print(f"\nExtracting scVI latent representation...")
latent_scvi = vae.get_latent_representation()
adata.obsm['X_scvi'] = latent_scvi
log_print(f"✓ Latent representation: {latent_scvi.shape}")

del latent_scvi
gc.collect()
log_memory("After scVI latent extraction")

# Compute UMAP from scVI
log_print(f"\nComputing UMAP from scVI latent space...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS, key_added='scvi')
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, neighbors_key='scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
log_print(f"✓ UMAP computed and saved as 'X_umap_scvi'")
log_memory("After scVI UMAP")

# Save checkpoint
checkpoint_file = checkpoint_dir / "checkpoint_04_scvi_trained.h5ad"
log_print(f"\nSaving checkpoint: {checkpoint_file.name}")
adata.write_h5ad(checkpoint_file, compression=CHECKPOINT_COMPRESSION)
log_print(f"✓ Checkpoint saved")

# =============================================================================
# Step 8: Train scANVI Model (Semi-supervised)
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 8: TRAIN scANVI MODEL (SEMI-SUPERVISED)")
log_print("="*80)

log_print(f"\nSetting up scANVI model...")
log_print(f"  Using filtered labels from: 'labels_for_scanvi'")
log_print(f"  Unlabeled category: '{UNLABELED_CATEGORY}'")

n_labeled = (adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY).sum()
n_unlabeled = (adata.obs['labels_for_scanvi'] == UNLABELED_CATEGORY).sum()
log_print(f"  Labeled cells: {n_labeled} ({n_labeled/adata.n_obs*100:.1f}%)")
log_print(f"  Unlabeled cells: {n_unlabeled} ({n_unlabeled/adata.n_obs*100:.1f}%)")

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key='labels_for_scanvi',
    unlabeled_category=UNLABELED_CATEGORY
)

log_print(f"\n✓ scANVI model created from scVI")
log_memory("After scANVI model creation")

# Free up scVI model memory
del vae
gc.collect()
log_memory("After deleting scVI model")

log_print(f"\nTraining scANVI model...")
scanvi_start = time.time()

lvae.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    early_stopping=True,
    early_stopping_patience=15,
    batch_size=128  # Smaller batch size for memory
)

scanvi_time = time.time() - scanvi_start

log_print(f"\n✓ scANVI training completed")
log_print(f"  Time: {scanvi_time:.1f}s ({scanvi_time/60:.1f} min)")
log_memory("After scANVI training")

# Save scANVI model
scanvi_model_dir = model_dir / "scanvi_model"
lvae.save(scanvi_model_dir, overwrite=True)
log_print(f"\n✓ scANVI model saved to: {scanvi_model_dir}")

# Extract scANVI results
log_print(f"\nExtracting scANVI results...")

# Latent representation
latent_scanvi = lvae.get_latent_representation()
adata.obsm['X_scanvi'] = latent_scanvi
log_print(f"✓ Latent representation: {latent_scanvi.shape}")

del latent_scanvi
gc.collect()

# Predictions
adata.obs['scanvi_predictions'] = lvae.predict()
log_print(f"✓ Predictions: {adata.obs['scanvi_predictions'].nunique()} cell types")

# Confidence scores
predictions_df = lvae.predict(soft=True)
adata.obs['scanvi_pred_confidence'] = predictions_df.max(axis=1).values
log_print(f"✓ Confidence scores computed")

# Probability matrix
adata.obsm['scanvi_probabilities'] = predictions_df.values
adata.uns['scanvi_celltype_order'] = predictions_df.columns.tolist()
log_print(f"✓ Probability matrix saved: {predictions_df.shape}")

del predictions_df
gc.collect()
log_memory("After scANVI extraction")

# Compute UMAP from scANVI
log_print(f"\nComputing UMAP from scANVI latent space...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=UMAP_N_NEIGHBORS, key_added='scanvi')
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, neighbors_key='scanvi')
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
log_print(f"✓ UMAP computed and saved as 'X_umap_scanvi'")
log_memory("After scANVI UMAP")

# Save checkpoint
checkpoint_file = checkpoint_dir / "checkpoint_05_scanvi_trained.h5ad"
log_print(f"\nSaving checkpoint: {checkpoint_file.name}")
adata.write_h5ad(checkpoint_file, compression=CHECKPOINT_COMPRESSION)
log_print(f"✓ Checkpoint saved")

# Free up scANVI model
del lvae
gc.collect()
log_memory("After deleting scANVI model")

# =============================================================================
# Step 9: Visualizations (Memory-Efficient)
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 9: GENERATE VISUALIZATIONS")
log_print("="*80)

# Set plotting style
sc.set_figure_params(dpi=DPI, frameon=False, figsize=(8, 6), facecolor='white')

# Figure 1: scVI UMAP (Batch Mixing)
log_print(f"\n[Figure 1] scVI UMAP - Batch Mixing...")
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()

sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0], show=False,
           title='scVI Integration: Batch Effect Correction',
           legend_loc='right margin', frameon=False)

sc.pl.umap(adata, color=UNIFIED_LINEAGE_COL, ax=axes[1], show=False,
           title='scVI Integration: Cell Lineages',
           legend_loc='right margin', frameon=False)

plt.tight_layout()
fig.savefig(fig_dir / f'scvi_umap_batch_mixing.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close(fig)
gc.collect()
log_print(f"  ✓ Saved: scvi_umap_batch_mixing.{FIGURE_FORMAT}")

# Figure 2: scANVI UMAP (Cell Type Annotations)
log_print(f"\n[Figure 2] scANVI UMAP - Cell Type Annotations...")
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()

# Original annotations
sc.pl.umap(adata, color=UNIFIED_SCANVI_ANNOTATION, ax=axes[0, 0], show=False,
           title='Original scANVI Annotations (from subsets)',
           legend_loc='right margin', frameon=False, legend_fontsize=6)

# scANVI predictions
sc.pl.umap(adata, color='scanvi_predictions', ax=axes[0, 1], show=False,
           title='scANVI Refined Predictions',
           legend_loc='right margin', frameon=False, legend_fontsize=6)

# Confidence scores (original)
sc.pl.umap(adata, color=UNIFIED_SCANVI_CONFIDENCE, ax=axes[1, 0], show=False,
           title='Original Annotation Confidence',
           cmap='viridis', frameon=False)

# Confidence scores (scANVI predictions)
sc.pl.umap(adata, color='scanvi_pred_confidence', ax=axes[1, 1], show=False,
           title='scANVI Prediction Confidence',
           cmap='viridis', frameon=False)

plt.tight_layout()
fig.savefig(fig_dir / f'scanvi_umap_annotations.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close(fig)
gc.collect()
log_print(f"  ✓ Saved: scanvi_umap_annotations.{FIGURE_FORMAT}")

# Figure 3: Comparison (scVI vs scANVI)
log_print(f"\n[Figure 3] Comparison: scVI vs scANVI...")
fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# scVI UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_scvi'].copy()
sc.pl.umap(adata, color='scanvi_predictions', ax=axes[0], show=False,
           title='scVI UMAP (Batch Corrected)',
           legend_loc='right margin', frameon=False, legend_fontsize=6)

# scANVI UMAP
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi'].copy()
sc.pl.umap(adata, color='scanvi_predictions', ax=axes[1], show=False,
           title='scANVI UMAP (Cell Type Refined)',
           legend_loc='right margin', frameon=False, legend_fontsize=6)

plt.tight_layout()
fig.savefig(fig_dir / f'comparison_scvi_vs_scanvi.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close(fig)
gc.collect()
log_print(f"  ✓ Saved: comparison_scvi_vs_scanvi.{FIGURE_FORMAT}")

# Figure 4: Cell Type Distribution by Batch
log_print(f"\n[Figure 4] Cell Type Distribution by Batch...")
fig, ax = plt.subplots(figsize=(14, 8))

# Create crosstab
ct_batch = pd.crosstab(
    adata.obs['scanvi_predictions'],
    adata.obs[BATCH_KEY],
    normalize='columns'
) * 100

# Plot heatmap
sns.heatmap(ct_batch, cmap='viridis', annot=False, fmt='.1f',
            cbar_kws={'label': '% of cells'}, ax=ax)
ax.set_xlabel('Dataset', fontsize=12)
ax.set_ylabel('Cell Type', fontsize=12)
ax.set_title('Cell Type Distribution Across Datasets (%)', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)

plt.tight_layout()
fig.savefig(fig_dir / f'celltype_by_batch.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
plt.close(fig)
gc.collect()
log_print(f"  ✓ Saved: celltype_by_batch.{FIGURE_FORMAT}")

log_memory("After all visualizations")

# =============================================================================
# Step 10: Final Summary and Save
# =============================================================================

log_print("\n" + "="*80)
log_print("STEP 10: FINAL SUMMARY AND SAVE")
log_print("="*80)

# Restore original X
adata.X = adata.layers['original_X'].copy()
del adata.layers['original_X']
gc.collect()

# Save final integrated data
final_file = OUTPUT_DIR / "adata_allcells_scvi_scanvi_integrated_v1_1.h5ad"
log_print(f"\nSaving final integrated data: {final_file.name}")
log_memory("Before final save")

adata.write_h5ad(final_file, compression='gzip')

log_print(f"✓ Final integrated data saved")
log_print(f"  Size: {final_file.stat().st_size / 1e9:.2f} GB")
log_memory("After final save")

# Generate summary report
summary_lines = [
    "="*80,
    "ALL CELLS scVI→scANVI INTEGRATION PIPELINE - SUMMARY REPORT",
    "="*80,
    "",
    f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    f"Output directory: {OUTPUT_DIR}",
    f"Version: 1.1 (Memory Optimized + CellTypist)",
    "",
    "="*80,
    "DATA SUMMARY",
    "="*80,
    f"  Final cells: {adata.shape[0]:,}",
    f"  Final genes: {adata.shape[1]:,}",
    f"  Batches: {adata.obs[BATCH_KEY].nunique()}",
    "",
    "="*80,
    "ANNOTATION SUMMARY",
    "="*80,
    f"  Cells with original scANVI annotations: {(adata.obs[UNIFIED_SCANVI_ANNOTATION] != UNLABELED_CATEGORY).sum():,}",
    f"  Unique cell types (original scANVI): {adata.obs[adata.obs[UNIFIED_SCANVI_ANNOTATION] != UNLABELED_CATEGORY][UNIFIED_SCANVI_ANNOTATION].nunique()}",
    f"  Mean confidence (original scANVI): {adata.obs[adata.obs[UNIFIED_SCANVI_ANNOTATION] != UNLABELED_CATEGORY][UNIFIED_SCANVI_CONFIDENCE].mean():.3f}",
    "",
]

# Add CellTypist summary if available
if UNIFIED_CELLTYPIST_PREDICTED in adata.obs.columns:
    n_ct = (adata.obs[UNIFIED_CELLTYPIST_PREDICTED] != 'Unknown').sum()
    summary_lines.append(f"  Cells with CellTypist annotations: {n_ct:,}")
    if n_ct > 0:
        summary_lines.append(f"  Unique CellTypist types: {adata.obs[adata.obs[UNIFIED_CELLTYPIST_PREDICTED] != 'Unknown'][UNIFIED_CELLTYPIST_PREDICTED].nunique()}")

summary_lines.extend([
    "",
    f"  Cells used for training: {(adata.obs['labels_for_scanvi'] != UNLABELED_CATEGORY).sum():,}",
    f"  Unlabeled cells: {(adata.obs['labels_for_scanvi'] == UNLABELED_CATEGORY).sum():,}",
    "",
    f"  scANVI predicted cell types: {adata.obs['scanvi_predictions'].nunique()}",
    f"  Mean prediction confidence: {adata.obs['scanvi_pred_confidence'].mean():.3f}",
    "",
    "Top 10 cell types (scANVI predictions):",
])

for ct, count in adata.obs['scanvi_predictions'].value_counts().head(10).items():
    pct = count / adata.n_obs * 100
    summary_lines.append(f"  {ct}: {count:,} cells ({pct:.1f}%)")

summary_lines.extend([
    "",
    "="*80,
    "TRAINING SUMMARY",
    "="*80,
    f"  scVI training time: {scvi_time:.1f}s ({scvi_time/60:.1f} min)",
    f"  scANVI training time: {scanvi_time:.1f}s ({scanvi_time/60:.1f} min)",
    f"  Total training time: {(scvi_time + scanvi_time):.1f}s ({(scvi_time + scanvi_time)/60:.1f} min)",
    "",
    "="*80,
    "OUTPUT FILES",
    "="*80,
    f"  Main data: {final_file.name}",
    f"  scVI model: models/scvi_model/",
    f"  scANVI model: models/scanvi_model/",
    f"  Figures: figures/",
    f"  Checkpoints: checkpoints/",
    f"  Log file: {LOG_FILE.name}",
    "",
    "="*80,
    "KEY COLUMNS IN adata.obs",
    "="*80,
    f"  {UNIFIED_SCANVI_ANNOTATION}: Original scANVI annotations from cell-type-specific analyses",
    f"  {UNIFIED_SCANVI_CONFIDENCE}: Original scANVI annotation confidence",
    f"  {UNIFIED_LINEAGE_COL}: Cell lineage source (B_cell, Myeloid, T_cell, Stromal, Epithelial)",
])

# Add CellTypist columns if present
if UNIFIED_CELLTYPIST_PREDICTED in adata.obs.columns:
    summary_lines.append(f"  {UNIFIED_CELLTYPIST_PREDICTED}: CellTypist predicted labels")
if UNIFIED_CELLTYPIST_MAJORITY in adata.obs.columns:
    summary_lines.append(f"  {UNIFIED_CELLTYPIST_MAJORITY}: CellTypist majority voting labels")
if UNIFIED_CELLTYPIST_CONFIDENCE in adata.obs.columns:
    summary_lines.append(f"  {UNIFIED_CELLTYPIST_CONFIDENCE}: CellTypist confidence scores")

summary_lines.extend([
    f"  labels_for_scanvi: Filtered labels used for scANVI training",
    f"  scanvi_predictions: scANVI refined predictions",
    f"  scanvi_pred_confidence: scANVI prediction confidence",
    "",
    "="*80,
    "KEY EMBEDDINGS IN adata.obsm",
    "="*80,
    f"  X_scvi: scVI latent representation ({adata.obsm['X_scvi'].shape})",
    f"  X_scanvi: scANVI latent representation ({adata.obsm['X_scanvi'].shape})",
    f"  X_umap_scvi: UMAP from scVI (batch mixing)",
    f"  X_umap_scanvi: UMAP from scANVI (cell type refined)",
    f"  scanvi_probabilities: scANVI prediction probabilities ({adata.obsm['scanvi_probabilities'].shape})",
    "",
    "="*80,
    "ENVIRONMENT",
    "="*80,
    f"  conda activate scvi",
    "",
    "="*80,
    "PIPELINE COMPLETE",
    "="*80
])

summary_text = "\n".join(summary_lines)
log_print("\n" + summary_text)

# Save summary report
summary_file = OUTPUT_DIR / "analysis_summary_v1_1.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

log_print(f"\n✓ Summary report saved: {summary_file.name}")
log_print(f"\n{'='*80}")
log_print(f"ALL PROCESSING COMPLETE!")
log_print(f"{'='*80}")
log_print(f"\nOutput directory: {OUTPUT_DIR}")
log_print(f"Check {LOG_FILE.name} for full details")
log_memory("Final memory usage")
