#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Step 2: scANVI Reference Mapping via scArches (Production Version)

Purpose:
- Map new query dataset to scanvi_existing_model
- Transfer cell_type labels with confidence filtering
- Project query cells into reference UMAP space

Adaptations for Your Training Pipeline:
- batch_key = "dataset" (matches your training)
- labels_key = "cell_type" (EXISTING_CELLTYPE_KEY)
- Gene naming: symbol (verified)
- HVG alignment: reads hvg_genes.txt from reference
- .X = counts (layer=None, matches your setup)

Requirements:
- Reference model must have UMAP operator (run step1_add_umap_operator.py first)
- Query genes must be in symbol format
- Query must have raw counts in layers['counts']

Author: r2end
Date: 2025-01-13
Version: 1.0 (Production)
"""

# ==============================================================================
# Step 0: Imports and Configuration
# ==============================================================================

import sys
import os
from pathlib import Path
import warnings
import json
import gc
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import scvi

warnings.filterwarnings('ignore')

print("=" * 70)
print("scArches scANVI Reference Mapping Pipeline")
print("=" * 70)
print(f"\nscanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"Python: {sys.version}")

import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")

# ==============================================================================
# CONFIGURATION SECTION (  MODIFY THIS)
# ==============================================================================

# ===== Input/Output Paths =====
SCANVI_REF_DIR = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX/models/scanvi_existing_model"
HVG_FILE = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX/models/scvi_model/hvg_genes.txt"

#   UPDATE THESE PATHS:
QUERY_H5AD = "/home/h2048/data/source/polyp/polyp_obj_updated_20260113.h5ad"  # Your query data
OUTPUT_DIR = "/home/h2048/data/py/0126/scarches_mapping_output"                 # Output directory

# ===== Keys (Must Match Reference Training) =====
BATCH_KEY = "sample"                # From your training Line 89
LABELS_KEY = "cell_type"             # EXISTING_CELLTYPE_KEY, Line 90
UNLABELED_CATEGORY = "Unknown"

# ===== scArches Training Parameters =====
MAX_EPOCHS = 200                     # scArches standard: 100-200
WEIGHT_DECAY = 0.0                   # Must be 0 for scArches
LEARNING_RATE = 5e-4                 # Same as your training Line 98
BATCH_SIZE = 256                     # Same as your training Line 99
EARLY_STOPPING_PATIENCE = 30
CHECK_VAL_EVERY_N_EPOCH = 10

# ===== Label Transfer Confidence Threshold =====
CONFIDENCE_THRESHOLD = 0.5           # Cells below this → Unknown

# ===== UMAP Parameters (Match Training) =====
UMAP_N_NEIGHBORS = 50                # From your training Line 108
UMAP_MIN_DIST = 0.5                  # From your training Line 107
UMAP_SPREAD = 1.0                    # From your training Line 108

# ===== Quality Control (Optional Pre-filtering) =====
APPLY_QC = False                     # Set True to enable QC filtering
MIN_GENES = 200
MAX_PCT_MT = 20

# ===== Visualization =====
DPI = 300
FIGURE_FORMAT = 'pdf'

# ===== Reproducibility =====
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
scvi.settings.seed = RANDOM_SEED

# Validate paths before proceeding
if not Path(SCANVI_REF_DIR).exists():
    print(f"\n❌ ERROR: Reference model not found at {SCANVI_REF_DIR}")
    print("   Please update SCANVI_REF_DIR in the configuration section")
    sys.exit(1)

if not Path(HVG_FILE).exists():
    print(f"\n❌ ERROR: HVG file not found at {HVG_FILE}")
    print("   Please update HVG_FILE in the configuration section")
    sys.exit(1)

if not Path(QUERY_H5AD).exists():
    print(f"\n❌ ERROR: Query data not found at {QUERY_H5AD}")
    print("   Please update QUERY_H5AD in the configuration section")
    sys.exit(1)

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)

print(f"\nOutput directory: {output_dir}")

# Save configuration
config = {
    'pipeline': 'scArches_scANVI_mapping',
    'version': '1.0',
    'reference_model': str(SCANVI_REF_DIR),
    'query_h5ad': str(QUERY_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'labels_key': LABELS_KEY,
    'confidence_threshold': CONFIDENCE_THRESHOLD,
    'max_epochs': MAX_EPOCHS,
    'weight_decay': WEIGHT_DECAY,
    'random_seed': RANDOM_SEED,
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'mapping_config.json', 'w') as f:
    json.dump(config, f, indent=2)

print("✓ Configuration saved")

# ==============================================================================
# Step 1: Load and Validate Query Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Loading and Validating Query Data")
print("=" * 70)

print(f"\nLoading query from: {QUERY_H5AD}")
adata_query = sc.read_h5ad(QUERY_H5AD)

print(f"\nQuery data loaded:")
print(f"  Cells: {adata_query.n_obs:,}")
print(f"  Genes: {adata_query.n_vars:,}")
print(f"  Obs columns: {list(adata_query.obs.columns[:10])}...")

# ==============================================================================
# Automatic Data Format Detection and Conversion
# ==============================================================================

print("\n" + "-" * 70)
print("Automatic Data Format Detection")
print("-" * 70)

def detect_data_type(adata):
    """
    Detect if .X contains raw counts or normalized data
    Returns: 'counts', 'log1p', or 'unknown'
    """
    X_data = adata.X
    if issparse(X_data):
        max_val = X_data.max()
        mean_val = X_data.mean()
    else:
        max_val = X_data.max()
        mean_val = X_data.mean()
    
    # Decision logic (from QUICK_REFERENCE_MEMORY)
    if max_val > 100 and mean_val < 50:
        return 'counts', max_val, mean_val
    elif max_val < 15 and mean_val < 5:
        return 'log1p', max_val, mean_val
    else:
        return 'unknown', max_val, mean_val

data_type, max_val, mean_val = detect_data_type(adata_query)

print(f"\n.X Statistics:")
print(f"  Max value: {max_val:.2f}")
print(f"  Mean value: {mean_val:.2f}")
print(f"  Detected type: {data_type.upper()}")

# Validate that .X contains raw counts
if data_type != 'counts':
    print(f"\n❌ ERROR: .X does not appear to contain raw counts")
    print(f"   Max={max_val:.2f}, Mean={mean_val:.2f}")
    print(f"\n   For scArches mapping, .X MUST contain raw counts (not log1p normalized)")
    print(f"\n   If your data is normalized:")
    print(f"     1. Re-load from original raw counts file")
    print(f"     2. Or check if raw counts exist in .raw or .layers")
    sys.exit(1)

print(f"✓ Confirmed: .X contains raw counts")

# Check if counts layer exists
if 'counts' not in adata_query.layers:
    print(f"\n  'counts' layer not found")
    print(f"   Creating layers['counts'] from .X (shared memory)...")
    
    # Use shared memory (no copy) - from QUICK_REFERENCE_MEMORY
    adata_query.layers['counts'] = adata_query.X
    print(f"   ✓ layers['counts'] created (0 extra memory)")
else:
    print(f"\n✓ Found existing layers['counts']")
    
    # Verify counts layer matches .X
    if issparse(adata_query.layers['counts']):
        layer_max = adata_query.layers['counts'].max()
    else:
        layer_max = adata_query.layers['counts'].max()
    
    if abs(layer_max - max_val) > 1:
        print(f"\n  Warning: layers['counts'] differs from .X")
        print(f"   .X max: {max_val:.2f}")
        print(f"   layers['counts'] max: {layer_max:.2f}")
        print(f"   Using .X as authoritative source...")
        adata_query.layers['counts'] = adata_query.X

# Set .raw (shared memory, no copy) - from QUICK_REFERENCE_MEMORY
if not hasattr(adata_query, 'raw') or adata_query.raw is None:
    print(f"\n  .raw not found")
    print(f"   Creating .raw from layers['counts'] (shared memory)...")
    
    adata_query.raw = sc.AnnData(
        X=adata_query.layers['counts'],  # Shared memory, no .copy()
        obs=adata_query.obs.copy(),
        var=adata_query.var.copy()
    )
    print(f"   ✓ .raw created (0 extra memory)")
else:
    print(f"\n✓ Found existing .raw")

print(f"\n✓ Data structure validated:")
print(f"  - .X: raw counts (max={max_val:.2f})")
print(f"  - layers['counts']: available")
print(f"  - .raw: available")

# Ensure batch_key exists
if BATCH_KEY not in adata_query.obs.columns:
    print(f"\n  '{BATCH_KEY}' not found, creating default batch")
    adata_query.obs[BATCH_KEY] = "query_batch"
else:
    print(f"\n✓ Found batch key: {BATCH_KEY}")
    print(f"  Batch distribution:")
    batch_counts = adata_query.obs[BATCH_KEY].value_counts()
    if len(batch_counts) <= 20:
        print(batch_counts)
    else:
        print(f"  Total batches: {len(batch_counts)}")
        print(f"  Top 10:")
        print(batch_counts.head(10))

# ==============================================================================
# Gene Naming Validation
# ==============================================================================

print("\n" + "-" * 70)
print("Gene Naming Validation")
print("-" * 70)

# Check gene naming format
sample_genes = adata_query.var_names[:5].tolist()
print(f"\nSample gene names (first 5): {sample_genes}")

# Detect gene naming system
if all(str(g).startswith('ENSG') for g in sample_genes):
    gene_format = 'ENSEMBL'
    print(f"✓ Detected: ENSEMBL IDs")
    print(f"\n❌ ERROR: Reference uses gene symbols, but query uses ENSEMBL IDs")
    print(f"\n   You must convert ENSEMBL IDs to symbols before mapping:")
    print(f"\n   import mygene")
    print(f"   mg = mygene.MyGeneInfo()")
    print(f"   result = mg.querymany(adata.var_names.tolist(),")
    print(f"                         scopes='ensembl.gene',")
    print(f"                         fields='symbol',")
    print(f"                         species='human')")
    print(f"   # Update adata.var_names with symbols")
    sys.exit(1)
elif all(str(g).replace('-', '').replace('_', '').isalnum() and len(str(g)) < 20 for g in sample_genes):
    gene_format = 'symbol'
    print(f"✓ Detected: Gene symbols")
    print(f"  Gene count: {adata_query.n_vars:,} (full genome scale)")
else:
    gene_format = 'unknown'
    print(f"  Warning: Cannot determine gene naming format")
    print(f"  Proceeding anyway, but verify gene names match reference")

# Verify gene names are unique
if adata_query.var_names.duplicated().any():
    n_dup = adata_query.var_names.duplicated().sum()
    print(f"\n  Warning: {n_dup} duplicate gene names found")
    print(f"  Making gene names unique...")
    adata_query.var_names_make_unique()
    print(f"  ✓ Gene names are now unique")

print(f"\n✓ Gene naming validation complete")

# Optional: Basic QC (only cell-level, NO gene filtering per P1 issue)
if APPLY_QC:
    print("\n" + "-" * 70)
    print("Applying Basic Quality Control (Cell-level Only)")
    print("-" * 70)
    
    n_cells_before = adata_query.n_obs
    
    # Calculate QC metrics
    adata_query.var['mt'] = adata_query.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_query, qc_vars=['mt'], percent_top=None, inplace=True)
    
    # Filter cells
    sc.pp.filter_cells(adata_query, min_genes=MIN_GENES)
    adata_query = adata_query[adata_query.obs.pct_counts_mt < MAX_PCT_MT].copy()
    
    n_cells_after = adata_query.n_obs
    print(f"  Filtered: {n_cells_before:,} → {n_cells_after:,} cells "
          f"({(1 - n_cells_after/n_cells_before)*100:.1f}% removed)")

print("\n✓ Query data validated")

# ==============================================================================
# Step 2: Gene Alignment with Reference (HVG)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Aligning Query Genes with Reference HVG")
print("=" * 70)

# Load HVG list from reference
print(f"\nLoading HVG list from: {HVG_FILE}")
hvg_genes_ref = pd.read_csv(HVG_FILE, header=None)[0].astype(str).tolist()
print(f"  Reference HVG genes: {len(hvg_genes_ref):,}")

# Check gene overlap
query_genes_set = set(adata_query.var_names.astype(str))
hvg_in_query = [g for g in hvg_genes_ref if g in query_genes_set]
overlap_pct = len(hvg_in_query) / len(hvg_genes_ref) * 100

print(f"\nGene overlap:")
print(f"  Query genes: {len(query_genes_set):,}")
print(f"  HVG in query: {len(hvg_in_query):,}/{len(hvg_genes_ref):,} ({overlap_pct:.1f}%)")

if overlap_pct < 70:
    print(f"\n❌ ERROR: Low gene overlap ({overlap_pct:.1f}%)")
    print(f"\n   Possible causes:")
    print(f"     1. Query uses different gene naming (ENSEMBL vs symbol)")
    print(f"     2. Query is from different organism/annotation version")
    print(f"     3. Query has been heavily filtered")
    print(f"\n   Reference genes (first 5): {hvg_genes_ref[:5]}")
    print(f"   Query genes (first 5): {list(query_genes_set)[:5]}")
    sys.exit(1)

print(f"\n✓ Gene overlap acceptable ({overlap_pct:.1f}%)")

# ==============================================================================
# Step 3: Prepare Query AnnData for scArches
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Preparing Query for scArches Mapping")
print("=" * 70)

# ⭐ CRITICAL: Move counts to .X (matches reference setup: layer=None)
# No .copy() needed - scvi doesn't modify X during prepare_query_anndata
print("\nSetting query.X = counts (matching reference setup)...")
adata_query.X = adata_query.layers['counts']  # No copy, saves memory
print("✓ Query.X set to raw counts (shared memory, 0 extra cost)")

# Prepare query data (this will pad missing genes with zeros)
print("\nRunning scvi prepare_query_anndata...")
print("  This will:")
print("    - Reorder genes to match reference")
print("    - Pad missing HVG genes with zeros")
print("    - Validate data structure")

try:
    scvi.model.SCANVI.prepare_query_anndata(
        adata_query,
        SCANVI_REF_DIR
    )
    print("✓ Query data prepared for scArches")
except Exception as e:
    print(f"\n❌ ERROR preparing query data: {e}")
    sys.exit(1)

# Verify preparation
print(f"\nQuery after preparation:")
print(f"  Shape: {adata_query.n_obs:,} × {adata_query.n_vars:,}")
print(f"  Genes match reference: {adata_query.n_vars == len(hvg_genes_ref)}")

# ==============================================================================
# Step 3.5: Metadata Cleaning (CRITICAL for load_query_data) [FIXED]
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3.5: Cleaning Metadata for scArches Compatibility (FIXED)")
print("=" * 70)

print("\nCleaning query metadata to prevent data type conflicts...")
print("  Strategy: Force strict string→categorical typing for batch_key and labels_key")

# Backup original obs columns (for later restore logic)
adata_query.uns['original_obs_columns'] = list(adata_query.obs.columns)
print(f"  ✓ Original metadata backed up: {len(adata_query.obs.columns)} columns")

# Build a minimal obs with strict typing
obs_clean = pd.DataFrame(index=adata_query.obs_names)

# --- Batch key: force string, fill NA, then categorical ---
print(f"\n1. Processing batch_key ('{BATCH_KEY}'):")
if BATCH_KEY in adata_query.obs.columns:
    b = adata_query.obs[BATCH_KEY]
else:
    print(f"     '{BATCH_KEY}' not found, creating default batch")
    b = pd.Series(["query_batch"] * adata_query.n_obs, index=adata_query.obs_names)

# Force to string (removes mixed float/str issues)
b = b.astype("string")
n_na_before = b.isna().sum()
b = b.fillna("unknown_batch")
obs_clean[BATCH_KEY] = pd.Categorical(b.astype(str))

print(f"   ✓ Converted to strict categorical")
print(f"   ✓ NA values filled: {n_na_before} → 0")
print(f"   ✓ Categories: {obs_clean[BATCH_KEY].nunique()}")

# --- Labels key: ALWAYS set to Unknown for scArches query mapping ---
# (This avoids treating query as partially labeled and avoids mixed-type categories)
print(f"\n2. Processing labels_key ('{LABELS_KEY}'):")
if LABELS_KEY in adata_query.obs.columns:
    # Keep backup for later comparison
    backup_col = f"{LABELS_KEY}_query_original"
    obs_clean[backup_col] = adata_query.obs[LABELS_KEY].astype("string").astype(str)
    print(f"   ✓ Original labels backed up to '{backup_col}'")
    print(f"   ✓ Original categories: {adata_query.obs[LABELS_KEY].nunique()}")

# Set ALL cells to Unknown (standard scArches query practice)
obs_clean[LABELS_KEY] = pd.Categorical([UNLABELED_CATEGORY] * adata_query.n_obs)
print(f"   ✓ All {adata_query.n_obs:,} cells set to '{UNLABELED_CATEGORY}' (scArches standard)")

# Replace obs with cleaned version
adata_query.obs = obs_clean

# Sanity checks
print(f"\n3. Final validation:")
for col in [BATCH_KEY, LABELS_KEY]:
    s = adata_query.obs[col]
    n_na = int(pd.isna(s).sum())
    cats = list(s.cat.categories[:10]) if hasattr(s, "cat") else []
    print(f"   ✓ {col}:")
    print(f"     - dtype: {s.dtype}")
    print(f"     - NA count: {n_na}")
    print(f"     - n_categories: {len(s.cat.categories)}")
    if len(cats) <= 10:
        print(f"     - categories: {cats}")
    else:
        print(f"     - categories (first 10): {cats}")

print("\n✓ Metadata cleaned successfully (strict categorical, no mixed types)")
print("  This eliminates '<' comparison errors during load_query_data")

# ==============================================================================
# Step 4: Load Query into Reference Model (scArches Surgery)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Loading Query into Reference Model")
print("=" * 70)

print(f"\nLoading reference model from: {SCANVI_REF_DIR}")
print("  This performs 'scArches surgery':")
print("    - Clones reference encoder/decoder")
print("    - Attaches query data")
print("    - Prepares for fine-tuning")
print(f"\n  ⭐ Note: Query cells are all '{UNLABELED_CATEGORY}' (scArches best practice)")
print(f"     This avoids mixed-type issues and ensures proper label transfer")

try:
    scanvi_query = scvi.model.SCANVI.load_query_data(
        adata_query,
        SCANVI_REF_DIR
    )
    print("✓ Query loaded into reference model")
except Exception as e:
    print(f"\n❌ ERROR loading query into model: {e}")
    sys.exit(1)

# ==============================================================================
# Step 5: Fine-tune Model on Query Data (scArches Training)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: Fine-tuning Model on Query Data")
print("=" * 70)

print(f"\nTraining parameters:")
print(f"  Max epochs: {MAX_EPOCHS}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Weight decay: {WEIGHT_DECAY} (must be 0 for scArches)")
print(f"  Early stopping patience: {EARLY_STOPPING_PATIENCE}")

train_kwargs = {
    'max_epochs': MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'early_stopping': True,
    'early_stopping_patience': EARLY_STOPPING_PATIENCE,
    'check_val_every_n_epoch': CHECK_VAL_EVERY_N_EPOCH,
    'plan_kwargs': {
        'lr': LEARNING_RATE,
        'weight_decay': WEIGHT_DECAY  # MUST be 0.0 for scArches
    },
    'enable_progress_bar': True
}

if gpu_available:
    train_kwargs['accelerator'] = 'gpu'
    train_kwargs['devices'] = 'auto'
    print(f"  Using GPU: {torch.cuda.get_device_name(0)}")

print("\nStarting training...")
try:
    scanvi_query.train(**train_kwargs)
    print("\n✓ Fine-tuning completed")
except Exception as e:
    print(f"\n❌ ERROR during training: {e}")
    sys.exit(1)

# ==============================================================================
# Step 6: Extract Latent Representations and Predictions
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Extracting Results (Latent, Predictions, Confidence)")
print("=" * 70)

# Get scANVI latent
print("\n1. Extracting latent representation...")
LATENT_KEY = "X_scANVI_mapped"
try:
    adata_query.obsm[LATENT_KEY] = scanvi_query.get_latent_representation()
    print(f"   ✓ Saved: .obsm['{LATENT_KEY}'] {adata_query.obsm[LATENT_KEY].shape}")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# Get hard predictions
print("\n2. Extracting cell type predictions...")
PRED_KEY = "cell_type_mapped"
try:
    adata_query.obs[PRED_KEY] = scanvi_query.predict()
    print(f"   ✓ Saved: .obs['{PRED_KEY}']")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# Get soft probabilities and confidence
print("\n3. Calculating prediction confidence...")
try:
    probs = scanvi_query.predict(soft=True)  # Shape: (n_cells, n_cell_types)
    confidence = probs.max(axis=1)
    
    CONF_KEY = "mapping_confidence"
    adata_query.obs[CONF_KEY] = confidence
    print(f"   ✓ Saved: .obs['{CONF_KEY}']")
    
    # Calculate margin (top1 - top2)
    top2_probs = np.partition(probs, -2, axis=1)[:, -2:]
    margin = top2_probs[:, 1] - top2_probs[:, 0]
    adata_query.obs['mapping_margin'] = margin
    print(f"   ✓ Saved: .obs['mapping_margin'] (top1 - top2 probability)")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# Apply confidence threshold
print(f"\n4. Applying confidence threshold ({CONFIDENCE_THRESHOLD})...")
FINAL_KEY = "cell_type_final"
adata_query.obs[FINAL_KEY] = np.where(
    confidence < CONFIDENCE_THRESHOLD,
    UNLABELED_CATEGORY,
    adata_query.obs[PRED_KEY].astype(str)
)
print(f"   ✓ Saved: .obs['{FINAL_KEY}']")

print("\n✓ All predictions extracted")

# ==============================================================================
# Step 7: UMAP Projection into Reference Space
# ==============================================================================

print("\n" + "=" * 70)
print("Step 7: UMAP Projection into Reference Space")
print("=" * 70)

UMAP_PROJECT_KEY = "X_umap_mapped"

if hasattr(scanvi_query, 'umap_op_'):
    print("✓ Reference model has UMAP operator")
    print("  Projecting query into reference UMAP coordinates...")
    
    try:
        adata_query.obsm[UMAP_PROJECT_KEY] = scanvi_query.umap_op_.transform(
            adata_query.obsm[LATENT_KEY]
        )
        
        print(f"  ✓ Saved: .obsm['{UMAP_PROJECT_KEY}'] {adata_query.obsm[UMAP_PROJECT_KEY].shape}")
        print("  ✓ Query and reference are now in the SAME UMAP space")
    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        sys.exit(1)
else:
    print("  Warning: Reference model doesn't have UMAP operator")
    print("  Computing new UMAP on query latent...")
    print("  Note: These coordinates will differ from reference")
    print("  → Run step1_add_umap_operator.py to fix this")
    
    import umap
    
    umap_op = umap.UMAP(
        n_neighbors=UMAP_N_NEIGHBORS,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=RANDOM_SEED
    )
    
    adata_query.obsm[UMAP_PROJECT_KEY] = umap_op.fit_transform(
        adata_query.obsm[LATENT_KEY]
    )
    
    print(f"  ✓ Saved: .obsm['{UMAP_PROJECT_KEY}']")

print("\n✓ UMAP projection completed")

# ==============================================================================
# Step 8: Quality Control and Summary Statistics
# ==============================================================================

print("\n" + "=" * 70)
print("Step 8: Quality Control Summary")
print("=" * 70)

# 1. Confidence distribution
print("\n1. Prediction Confidence Distribution:")
print(adata_query.obs[CONF_KEY].describe())

low_conf_count = (adata_query.obs[CONF_KEY] < CONFIDENCE_THRESHOLD).sum()
low_conf_pct = low_conf_count / adata_query.n_obs * 100
print(f"\n   Cells with confidence < {CONFIDENCE_THRESHOLD}: {low_conf_count:,} ({low_conf_pct:.1f}%)")

# 2. Margin distribution
print("\n2. Prediction Margin (Top1 - Top2):")
print(adata_query.obs['mapping_margin'].describe())

low_margin_count = (adata_query.obs['mapping_margin'] < 0.1).sum()
low_margin_pct = low_margin_count / adata_query.n_obs * 100
print(f"\n   Cells with margin < 0.1: {low_margin_count:,} ({low_margin_pct:.1f}%)")
print(f"   (Low margin = ambiguous predictions)")

# 3. Final cell type distribution
print("\n3. Final Cell Type Distribution:")
celltype_counts = adata_query.obs[FINAL_KEY].value_counts()
print(celltype_counts.head(20))

# 4. Unknown cells
unknown_count = (adata_query.obs[FINAL_KEY] == UNLABELED_CATEGORY).sum()
unknown_pct = unknown_count / adata_query.n_obs * 100
print(f"\n4. {UNLABELED_CATEGORY} Cells:")
print(f"   Count: {unknown_count:,} ({unknown_pct:.1f}%)")

if unknown_pct > 30:
    print(f"\n     High percentage of Unknown cells ({unknown_pct:.1f}%)")
    print(f"   Possible causes:")
    print(f"     - Reference doesn't cover query cell types")
    print(f"     - Query data quality issues")
    print(f"     - Confidence threshold too strict ({CONFIDENCE_THRESHOLD})")

# 5. Batch effects (if multiple batches)
if adata_query.obs[BATCH_KEY].nunique() > 1:
    print("\n5. Cell Type Distribution by Batch:")
    batch_celltype = pd.crosstab(
        adata_query.obs[BATCH_KEY],
        adata_query.obs[FINAL_KEY],
        normalize='index'
    ) * 100
    print(batch_celltype.round(1).head())

print("\n✓ Quality control completed")

# ==============================================================================
# Step 9: Restore Original Metadata and Save
# ==============================================================================

print("\n" + "=" * 70)
print("Step 9: Restoring Original Metadata and Saving")
print("=" * 70)

# Restore original metadata
if 'original_obs_columns' in adata_query.uns:
    print("\nRestoring original metadata columns...")
    
    # Load original h5ad again to get full metadata
    print(f"  Re-loading original file: {QUERY_H5AD}")
    adata_original = sc.read_h5ad(QUERY_H5AD)
    
    # Verify cell order matches
    if not (adata_original.obs_names == adata_query.obs_names).all():
        print("  Warning: Cell order differs, aligning...")
        adata_original = adata_original[adata_query.obs_names].copy()
    
    # Merge original metadata with new predictions
    print(f"  Original metadata: {adata_original.obs.shape[1]} columns")
    print(f"  New predictions: {len([c for c in adata_query.obs.columns if c not in [BATCH_KEY, LABELS_KEY]])} columns")
    
    # Keep new prediction columns
    new_cols = [c for c in adata_query.obs.columns if c not in adata_original.obs.columns]
    
    # Merge
    for col in new_cols:
        adata_original.obs[col] = adata_query.obs[col].values
    
    # Replace obsm (latent, UMAP)
    for key in adata_query.obsm.keys():
        adata_original.obsm[key] = adata_query.obsm[key]
    
    # Use the merged version for saving
    adata_query = adata_original
    
    print(f"  ✓ Restored {adata_original.obs.shape[1]} original columns")
    print(f"  ✓ Added {len(new_cols)} new prediction columns")
    
    del adata_original
    gc.collect()
else:
    print("\nNo metadata restoration needed")

output_file = output_dir / "query_mapped_to_reference.h5ad"
print(f"\nSaving to: {output_file}")

try:
    adata_query.write_h5ad(output_file, compression='gzip', compression_opts=9)
    file_size = output_file.stat().st_size / (1024**3)
    print(f"✓ Saved ({file_size:.2f} GB)")
except Exception as e:
    print(f"❌ ERROR saving file: {e}")
    sys.exit(1)

# Memory cleanup
del scanvi_query
gc.collect()

# ==============================================================================
# Step 10: Visualization
# ==============================================================================

print("\n" + "=" * 70)
print("Step 10: Generating Visualizations")
print("=" * 70)

# Set UMAP for plotting
adata_query.obsm['X_umap'] = adata_query.obsm[UMAP_PROJECT_KEY].copy()

# Create figure
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

try:
    # Row 1: Basic
    sc.pl.umap(adata_query, color=BATCH_KEY, ax=axes[0, 0], show=False, 
               title='Batch Distribution')
    sc.pl.umap(adata_query, color=PRED_KEY, ax=axes[0, 1], show=False,
               title='Predicted Cell Types (Raw)')
    sc.pl.umap(adata_query, color=FINAL_KEY, ax=axes[0, 2], show=False,
               title=f'Final Cell Types (conf>{CONFIDENCE_THRESHOLD})')
    
    # Row 2: Quality metrics
    sc.pl.umap(adata_query, color=CONF_KEY, ax=axes[1, 0], show=False,
               title='Prediction Confidence', cmap='viridis')
    sc.pl.umap(adata_query, color='mapping_margin', ax=axes[1, 1], show=False,
               title='Prediction Margin (Top1-Top2)', cmap='viridis')
    
    # Confidence histogram
    axes[1, 2].hist(adata_query.obs[CONF_KEY], bins=50, edgecolor='black')
    axes[1, 2].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--',
                       label=f'Threshold={CONFIDENCE_THRESHOLD}')
    axes[1, 2].set_xlabel('Prediction Confidence')
    axes[1, 2].set_ylabel('Number of Cells')
    axes[1, 2].set_title('Confidence Distribution')
    axes[1, 2].legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / "figures" / f"mapping_overview.{FIGURE_FORMAT}",
                dpi=DPI, bbox_inches='tight')
    plt.savefig(output_dir / "figures" / "mapping_overview.png",
                dpi=150, bbox_inches='tight')
    print(f"✓ Saved: mapping_overview.{FIGURE_FORMAT}")
    
    plt.close()
except Exception as e:
    print(f"  Warning: Visualization failed: {e}")

# ==============================================================================
# Step 11: Export Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Step 11: Generating Summary Report")
print("=" * 70)

# Create summary
summary = {
    'Query': {
        'cells': int(adata_query.n_obs),
        'genes': int(adata_query.n_vars),
        'batches': int(adata_query.obs[BATCH_KEY].nunique())
    },
    'Gene_Overlap': {
        'reference_hvg': len(hvg_genes_ref),
        'overlap_count': len(hvg_in_query),
        'overlap_percent': float(overlap_pct)
    },
    'Training': {
        'epochs': MAX_EPOCHS,
        'batch_size': BATCH_SIZE,
        'learning_rate': LEARNING_RATE,
        'weight_decay': WEIGHT_DECAY
    },
    'Predictions': {
        'unique_cell_types': int(adata_query.obs[FINAL_KEY].nunique()),
        'median_confidence': float(adata_query.obs[CONF_KEY].median()),
        'median_margin': float(adata_query.obs['mapping_margin'].median()),
        'low_confidence_cells': int(low_conf_count),
        'low_confidence_percent': float(low_conf_pct),
        'unknown_cells': int(unknown_count),
        'unknown_percent': float(unknown_pct)
    },
    'Top_Cell_Types': dict(celltype_counts.head(10).to_dict())
}

# Save as JSON
summary_json = output_dir / "mapping_summary.json"
with open(summary_json, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"✓ Saved JSON summary: {summary_json}")

# Save as text
summary_text = []
summary_text.append("=" * 70)
summary_text.append("scArches scANVI Mapping Summary")
summary_text.append("=" * 70)
summary_text.append(f"\nTimestamp: {datetime.now().isoformat()}")
summary_text.append(f"\nReference Model: {SCANVI_REF_DIR}")
summary_text.append(f"Query Data: {QUERY_H5AD}")
summary_text.append(f"Output: {output_file}")

summary_text.append("\n\n[ Query Dataset ]")
summary_text.append(f"  Cells: {adata_query.n_obs:,}")
summary_text.append(f"  Genes: {adata_query.n_vars:,}")
summary_text.append(f"  Batches: {adata_query.obs[BATCH_KEY].nunique()}")

summary_text.append("\n[ Gene Overlap ]")
summary_text.append(f"  Reference HVG: {len(hvg_genes_ref):,}")
summary_text.append(f"  Overlap: {len(hvg_in_query):,} ({overlap_pct:.1f}%)")

summary_text.append("\n[ Mapping Quality ]")
summary_text.append(f"  Median confidence: {adata_query.obs[CONF_KEY].median():.3f}")
summary_text.append(f"  Median margin: {adata_query.obs['mapping_margin'].median():.3f}")
summary_text.append(f"  Low confidence cells: {low_conf_count:,} ({low_conf_pct:.1f}%)")
summary_text.append(f"  Unknown cells: {unknown_count:,} ({unknown_pct:.1f}%)")

summary_text.append("\n[ Cell Type Distribution ]")
for ct, count in celltype_counts.head(10).items():
    pct = count / adata_query.n_obs * 100
    summary_text.append(f"  {ct}: {count:,} ({pct:.1f}%)")

summary_text.append("\n" + "=" * 70)
summary_text.append("Mapping completed successfully")
summary_text.append("=" * 70)

summary_txt = '\n'.join(summary_text)
print("\n" + summary_txt)

summary_txt_file = output_dir / "mapping_summary.txt"
with open(summary_txt_file, 'w') as f:
    f.write(summary_txt)
print(f"\n✓ Saved text summary: {summary_txt_file}")

# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

print("\n" + "=" * 70)
print("✅ scArches scANVI MAPPING COMPLETE")
print("=" * 70)

print(f"\n📊 Key Outputs:")
print(f"  1. Mapped data: {output_file}")
print(f"  2. Visualizations: {output_dir}/figures/")
print(f"  3. Summary: {summary_json}")

print(f"\n📈 Next Steps:")
print(f"  - Review confidence distribution (low confidence = uncertain cells)")
print(f"  - Check marker gene expression for validation")
print(f"  - Compare query UMAP with reference (if you have reference h5ad)")
print(f"  - Adjust confidence threshold if needed ({CONFIDENCE_THRESHOLD})")

print("\n" + "=" * 70)
print("Thank you for using scArches scANVI mapping pipeline!")
print("=" * 70)