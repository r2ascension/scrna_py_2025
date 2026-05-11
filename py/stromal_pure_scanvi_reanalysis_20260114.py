#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pure Stromal Cell Re-analysis: scVI → CellTypist → scANVI Pipeline
===================================================================

Purpose:
- Load adata_stromal_FINAL.h5ad (contains contaminating epithelial/immune cells)
- Filter to PURE stromal cells: Endothelial, Fibroblasts, Smooth Muscle, Schwann
- Remove: Basal, CD8_TRM/EM, Ciliated, SMG_Serous (epithelial/immune contaminants)
- Re-run scVI → CellTypist → scANVI for refined annotation

Based on: stromal_vascular_scvi_celltypist_scanvi_pipeline_v3.5.1
Critical Fixes: v3.5.1 (P0/P1/P2 fixes)

Author: r2end
Date: 2025-01-15
Version: 1.0
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import time
import gc
from datetime import datetime

import scvi
import torch
import celltypist
from celltypist import models

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output
INPUT_FILE = "/home/h2048/data/py/0111/celltypist_stromal/adata_stromal_FINAL.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/stromal_pure_scanvi")

# Cell type filtering
# ⭐ KEEP these cell types (pure stromal)
KEEP_PATTERNS = [
    'Endothelia',      # All endothelial subtypes
    'Fibro',           # All fibroblast subtypes
    'Muscle',          # Smooth muscle and pericytes
    'Schwann'          # Neural cells
]

# ⭐ REMOVE these cell types (contaminants)
REMOVE_PATTERNS = [
    'Basal',           # Epithelial basal cells
    'Ciliated',        # Epithelial ciliated cells
    'CD8',             # T cells
    'CD4',             # T cells
    'TRM',             # Tissue-resident memory T cells
    'SMG',             # Salivary/mucous glands
    'Goblet',          # Epithelial goblet cells
    'Secretory',       # Epithelial secretory cells
]

# Model parameters
N_LATENT = 75
N_LAYERS = 3
N_HVG = 4000

# Training
SCVI_MAX_EPOCHS = 800
SCANVI_MAX_EPOCHS = 600
BATCH_SIZE = 2048
LEARNING_RATE = 1e-3
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 50

# CellTypist
CELLTYPIST_MODEL = "/home/h2048/data/source/reference/celltypist_models/Cells_Lung_Airway.pkl"
CELLTYPIST_MAJORITY_VOTING = True
CELLTYPIST_MIN_CONFIDENCE = 0.5

# Filtering
MIN_CELLS_PER_TYPE = 50

# Keys
BATCH_KEY = 'sample'
ORIGINAL_LABEL_KEY = 'cell_type_scanvi_filt'  # From previous analysis

# Seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

print("="*80)
print("PURE STROMAL CELL RE-ANALYSIS PIPELINE")
print("="*80)
print(f"\nGoal: Remove epithelial/immune contaminants, re-annotate pure stromal cells")
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Keep patterns: {KEEP_PATTERNS}")
print(f"  Remove patterns: {REMOVE_PATTERNS}")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def merge_rare_types(s, min_cells=50, other="Unknown"):
    """Merge rare cell types"""
    s = s.astype(str).copy()
    vc = s.value_counts()
    rare = vc[vc < min_cells].index
    
    if len(rare) > 0:
        print(f"   Merging {len(rare)} rare types (<{min_cells} cells) → {other}")
        for rt in rare[:5]:
            print(f"      {rt}: {vc[rt]} cells")
        if len(rare) > 5:
            print(f"      ... and {len(rare)-5} more")
        s.loc[s.isin(rare)] = other
    
    return s

def save_celltype_counts(counts_dict, output_dir):
    """Save cell type counts to CSV"""
    for name, series in counts_dict.items():
        df = pd.DataFrame({
            'cell_type': series.index,
            'count': series.values,
            'percentage': (series.values / series.sum() * 100).round(2)
        })
        df = df.sort_values('count', ascending=False)
        filepath = output_dir / f"{name}.csv"
        df.to_csv(filepath, index=False)
        print(f"   ✓ Saved: {filepath.name}")

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print(f"\n{'='*80}")
print("ENVIRONMENT SETUP")
print("="*80)

# GPU
GPU_AVAILABLE = torch.cuda.is_available()
print(f"\nGPU available: {GPU_AVAILABLE}")

if GPU_AVAILABLE:
    gpu_name = torch.cuda.get_device_name(0)
    gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"GPU: {gpu_name}")
    print(f"Memory: {gpu_memory:.1f} GB")
    accelerator = 'gpu'
    devices = 1
else:
    print("⚠️  Running on CPU")
    accelerator = 'cpu'
    devices = 'auto'

# Directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR = OUTPUT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)
MODEL_DIR = OUTPUT_DIR / "models"
MODEL_DIR.mkdir(exist_ok=True)

# Scanpy settings
sc.settings.verbosity = 3
sc.settings.n_jobs = 48
sc.settings.figdir = FIG_DIR
sc.set_figure_params(dpi=300, facecolor='white', format='pdf')

print(f"✓ Output directory: {OUTPUT_DIR}")

# ============================================================================
# STAGE 0: DATA LOADING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 0: DATA LOADING")
print("="*80)

print(f"\nLoading: {INPUT_FILE}")
adata = sc.read_h5ad(INPUT_FILE)

print(f"\n✓ Data loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Batches: {adata.obs[BATCH_KEY].nunique()}")

# Check annotation column
if ORIGINAL_LABEL_KEY not in adata.obs.columns:
    print(f"\n⚠️  '{ORIGINAL_LABEL_KEY}' not found, checking alternatives...")
    alt_cols = ['cell_type_scanvi_raw', 'cell_type_celltypist_filt', 
                'celltypist_majority_voting', 'predicted_labels']
    for col in alt_cols:
        if col in adata.obs.columns:
            ORIGINAL_LABEL_KEY = col
            print(f"  ✓ Using: {col}")
            break
    
    if ORIGINAL_LABEL_KEY not in adata.obs.columns:
        raise ValueError(f"Cannot find cell type annotation column")

# Show original distribution
print(f"\n{'='*80}")
print("ORIGINAL CELL TYPE DISTRIBUTION")
print("="*80)
original_counts = adata.obs[ORIGINAL_LABEL_KEY].value_counts()
print(f"\nTotal unique types: {len(original_counts)}")
print(f"\nTop 20 types:")
for ct, count in original_counts.head(20).items():
    pct = count / adata.n_obs * 100
    print(f"  {ct}: {count:,} ({pct:.2f}%)")

# ============================================================================
# STAGE 1: CELL TYPE FILTERING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 1: FILTER TO PURE STROMAL CELLS")
print("="*80)

print(f"\n⭐ Filtering strategy:")
print(f"  KEEP: {', '.join(KEEP_PATTERNS)}")
print(f"  REMOVE: {', '.join(REMOVE_PATTERNS)}")

# Create filter
cell_types = adata.obs[ORIGINAL_LABEL_KEY].astype(str)

# Step 1: Mark cells to KEEP
keep_mask = pd.Series(False, index=adata.obs_names)
for pattern in KEEP_PATTERNS:
    matches = cell_types.str.contains(pattern, case=False, na=False)
    n_match = matches.sum()
    print(f"\n  Pattern '{pattern}': {n_match:,} cells")
    keep_mask |= matches

print(f"\n  Total KEEP candidates: {keep_mask.sum():,}")

# Step 2: Mark cells to REMOVE (override KEEP)
remove_mask = pd.Series(False, index=adata.obs_names)
for pattern in REMOVE_PATTERNS:
    matches = cell_types.str.contains(pattern, case=False, na=False)
    n_match = matches.sum()
    if n_match > 0:
        print(f"  Pattern '{pattern}': {n_match:,} cells (will remove)")
        remove_mask |= matches

print(f"\n  Total REMOVE: {remove_mask.sum():,}")

# Final filter
final_keep = keep_mask & ~remove_mask
n_keep = final_keep.sum()
n_remove = (~final_keep).sum()

print(f"\n{'='*80}")
print(f"FILTERING RESULTS:")
print(f"  Before: {adata.n_obs:,} cells")
print(f"  Keep: {n_keep:,} ({n_keep/adata.n_obs*100:.1f}%)")
print(f"  Remove: {n_remove:,} ({n_remove/adata.n_obs*100:.1f}%)")
print("="*80)

if n_keep < 1000:
    print(f"\n⚠️  WARNING: Only {n_keep} cells remaining!")
    print(f"  Check if filtering patterns are too strict")

# Show what's being removed
print(f"\n📊 Cell types being REMOVED:")
removed_types = adata.obs.loc[~final_keep, ORIGINAL_LABEL_KEY].value_counts()
for ct, count in removed_types.items():
    print(f"  ✗ {ct}: {count:,} cells")

# Show what's being kept
print(f"\n📊 Cell types being KEPT:")
kept_types = adata.obs.loc[final_keep, ORIGINAL_LABEL_KEY].value_counts()
for ct, count in kept_types.items():
    pct = count / n_keep * 100
    print(f"  ✓ {ct}: {count:,} ({pct:.2f}%)")

# Save filtering statistics
filter_stats = pd.DataFrame({
    'cell_type': original_counts.index,
    'original_count': original_counts.values,
    'kept_count': [kept_types.get(ct, 0) for ct in original_counts.index],
    'removed_count': [removed_types.get(ct, 0) for ct in original_counts.index],
    'action': ['KEEP' if ct in kept_types.index else 'REMOVE' 
               for ct in original_counts.index]
})
filter_stats.to_csv(OUTPUT_DIR / "filtering_statistics.csv", index=False)
print(f"\n✓ Filtering statistics saved")

# Apply filter
print(f"\n⭐ Applying filter...")
adata = adata[final_keep, :].copy()
print(f"✓ Filtered data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Clean up metadata
gc.collect()

# ============================================================================
# STAGE 2: DATA STRUCTURE CHECK & RESTORATION
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 2: DATA STRUCTURE CHECK")
print("="*80)

print(f"\nChecking data structure...")
print(f"  Layers: {list(adata.layers.keys())}")
print(f"  .raw: {adata.raw is not None}")
if adata.raw is not None:
    print(f"  .raw genes: {adata.raw.n_vars:,}")

# Ensure counts layer exists
if 'counts' not in adata.layers:
    print(f"\n⚠️  'counts' layer missing, reconstructing...")
    if adata.raw is not None:
        adata.layers['counts'] = adata.raw.X.copy()
        print(f"  ✓ Created from .raw")
    elif 'log1p' in adata.layers:
        print(f"  ⚠️  WARNING: Only log1p available, cannot reconstruct counts")
        raise ValueError("Need raw counts for scVI training")
    else:
        raise ValueError("Cannot find counts data")

# Ensure .raw exists (will be updated after HVG selection)
if adata.raw is None or adata.raw.n_vars < adata.n_vars:
    print(f"\n⚠️  .raw missing or incomplete, preserving full genes...")
    adata.raw = sc.AnnData(
        X=adata.layers['counts'],
        obs=adata.obs.copy(),
        var=adata.var.copy()
    )
    print(f"  ✓ Created adata.raw: {adata.raw.n_vars:,} genes")

# Gene name standardization
print(f"\n⭐ Checking gene names...")
if 'symbol_base' not in adata.var.columns:
    print(f"  Creating 'symbol_base' column...")
    for col in ['gene_symbol', 'gene_symbols', 'symbol', 'feature_name']:
        if col in adata.var.columns:
            adata.var['symbol_base'] = adata.var[col].astype(str)
            print(f"  ✓ Using: {col}")
            break
    
    if 'symbol_base' not in adata.var.columns:
        adata.var['symbol_base'] = adata.var_names.astype(str)
        print(f"  ✓ Using var_names")

# Patch symbol_base into raw if needed
if 'symbol_base' not in adata.raw.var.columns:
    common = adata.raw.var_names.intersection(adata.var_names)
    if len(common) > 0:
        adata.raw.var.loc[common, 'symbol_base'] = adata.var.loc[common, 'symbol_base'].values
        print(f"  ✓ Patched symbol_base into .raw ({len(common):,} genes)")

print(f"\n✓ Data structure validated")

# ============================================================================
# STAGE 3: BATCH SIZE VALIDATION
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 3: BATCH SIZE VALIDATION")
print("="*80)

batch_counts = adata.obs[BATCH_KEY].value_counts()
print(f"\nBatch distribution:")
print(f"  Total batches: {len(batch_counts)}")
print(f"  Min size: {batch_counts.min()}")
print(f"  Max size: {batch_counts.max()}")
print(f"  Mean size: {batch_counts.mean():.0f}")

# Remove small batches (<3 cells)
small_batches = batch_counts[batch_counts < 3].index
if len(small_batches) > 0:
    print(f"\n⚠️  Found {len(small_batches)} small batches (<3 cells):")
    for batch in small_batches:
        print(f"  {batch}: {batch_counts[batch]} cells")
    
    print(f"\n  Removing small batches...")
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  ✓ After removal: {adata.n_obs:,} cells")
    
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    print(f"  ✓ Remaining batches: {len(batch_counts)}")

print(f"\n✓ Batch validation complete")

# ============================================================================
# STAGE 4: HVG SELECTION
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 4: HVG SELECTION")
print("="*80)

print(f"\n⭐ Selecting {N_HVG} highly variable genes...")

try:
    print(f"   Attempting batch-aware HVG selection...")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        subset=False,
        flavor='seurat_v3'
    )
    hvg_method = "batch-aware"
    print(f"   ✓ Success")
except Exception as e:
    print(f"   ⚠️  Batch-aware failed: {str(e)[:80]}")
    print(f"   Falling back to non-batch-aware...")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        subset=False,
        flavor='seurat_v3'
    )
    hvg_method = "non-batch-aware"

hvg_mask = adata.var['highly_variable'].to_numpy()
n_hvg_selected = int(hvg_mask.sum())
print(f"\n✓ Selected {n_hvg_selected:,} HVGs ({hvg_method})")

# Save HVG list
hvg_file = MODEL_DIR / "hvg_genes.txt"
genes_to_save = adata.var_names[hvg_mask].astype(str)
pd.Series(genes_to_save.values).to_csv(hvg_file, index=False, header=False)
print(f"✓ Gene list saved: {hvg_file}")

# ============================================================================
# STAGE 5: FORCE INCLUDE STROMAL MARKERS
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 5: FORCE INCLUDE STROMAL MARKERS")
print("="*80)

# Stromal core markers
STROMAL_MARKERS = {
    'Pan_Endothelial': ['PECAM1', 'CDH5', 'VWF', 'CLDN5'],
    'Arterial_EC': ['GJA5', 'EFNB2', 'BMX'],
    'Venous_EC': ['ACKR1', 'NR2F2', 'EPHB4'],
    'Capillary_EC': ['CA4', 'RGCC', 'PLVAP'],
    'Lymphatic_EC': ['PROX1', 'LYVE1', 'PDPN', 'FLT4'],
    'Pan_Fibroblast': ['COL1A1', 'COL1A2', 'COL3A1', 'DCN', 'VIM'],
    'Myofibroblast': ['ACTA2', 'TAGLN', 'MYH11'],
    'Pan_SMC': ['ACTA2', 'MYH11', 'TAGLN', 'MYLK', 'CNN1'],
    'Pericytes': ['RGS5', 'PDGFRB', 'CSPG4', 'MCAM'],
    'Schwann': ['S100B', 'MPZ', 'PMP22', 'SOX10']
}

ALL_MARKERS = []
for genes in STROMAL_MARKERS.values():
    ALL_MARKERS.extend(genes)
ALL_MARKERS = list(set(ALL_MARKERS))

print(f"\nTotal markers to check: {len(ALL_MARKERS)}")

# Check markers
gene_symbols = adata.var['symbol_base'].astype(str)
missing_in_data = []
missing_in_hvg = []
present_in_hvg = []

print(f"\nChecking markers:")
for category, genes in STROMAL_MARKERS.items():
    print(f"\n{category}:")
    for gene in genes:
        matches = gene_symbols.str.upper() == gene.upper()
        
        if matches.sum() == 0:
            print(f"  ✗ {gene}: NOT IN DATA")
            missing_in_data.append(gene)
        else:
            idx = matches[matches].index[0]
            is_hvg = adata.var.loc[idx, 'highly_variable']
            
            if is_hvg:
                print(f"  ✓ {gene}: in HVG")
                present_in_hvg.append(gene)
            else:
                print(f"  ⚠️  {gene}: adding to HVG")
                missing_in_hvg.append(gene)

# Add missing markers to HVG
if len(missing_in_hvg) > 0:
    print(f"\n🔧 Adding {len(missing_in_hvg)} markers to HVG...")
    
    n_added = 0
    for gene in missing_in_hvg:
        matches = gene_symbols.str.upper() == gene.upper()
        if matches.sum() > 0:
            idx = matches[matches].index[0]
            adata.var.loc[idx, 'highly_variable'] = True
            hvg_mask[adata.var_names.get_loc(idx)] = True
            n_added += 1
    
    new_hvg_count = adata.var['highly_variable'].sum()
    print(f"✅ HVG updated: {n_hvg_selected} → {new_hvg_count} (+{n_added})")
    
    # Re-save gene list
    genes_to_save = adata.var_names[hvg_mask].astype(str)
    pd.Series(genes_to_save.values).to_csv(hvg_file, index=False, header=False)

print(f"\n✓ Final HVG count: {adata.var['highly_variable'].sum()}")

# ============================================================================
# STAGE 6: CellTypist ANNOTATION
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 6: CellTypist ANNOTATION")
print("="*80)

# Load model
print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL}")
try:
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"✓ Model loaded")
except:
    print(f"  Downloading...")
    models.download_models(model=CELLTYPIST_MODEL)
    model = models.Model.load(model=CELLTYPIST_MODEL)

# Get model features
print(f"\n⭐ Validating model features...")
model_features = None
for attr in ['features', 'genes', 'var_names']:
    if hasattr(model, attr):
        mf = getattr(model, attr)
        if mf is not None and len(mf) > 1000:
            model_features = pd.Index(mf).astype(str)
            print(f"  ✓ Found {len(model_features)} features via '{attr}'")
            break

if model_features is None:
    raise ValueError("Cannot retrieve CellTypist model features")

# Build input
print(f"\n⭐ Building CellTypist input...")
sym = pd.Index(adata.raw.var['symbol_base']).astype(str)
keep_mask = sym.isin(model_features) & ~sym.duplicated(keep='first')
n_overlap = keep_mask.sum()

print(f"  Our data: {len(sym)} genes")
print(f"  Overlap: {n_overlap} genes")
print(f"  Coverage: {n_overlap/len(model_features)*100:.1f}%")

X_ct = adata.raw.X[:, keep_mask].copy()
var_ct = pd.DataFrame(index=sym[keep_mask])

adata_celltypist = sc.AnnData(
    X=X_ct,
    obs=pd.DataFrame(index=adata.obs_names),
    var=var_ct
)

print(f"✓ CellTypist input: {adata_celltypist.n_obs:,} × {adata_celltypist.n_vars:,}")

# Normalize
print(f"\n  Normalizing...")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
sc.pp.log1p(adata_celltypist)

# Predict
print(f"\n⭐ Running CellTypist prediction...")
start_time = time.time()

predictions = celltypist.annotate(
    adata_celltypist,
    model=model,
    majority_voting=CELLTYPIST_MAJORITY_VOTING
)

elapsed = time.time() - start_time
print(f"✓ Prediction complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# Extract results
pred_df = predictions.predicted_labels
pred_df = pred_df.reindex(adata.obs_names)

if CELLTYPIST_MAJORITY_VOTING and 'majority_voting' in pred_df.columns:
    use_col = 'majority_voting'
else:
    use_col = 'predicted_labels'

adata.obs['cell_type_celltypist_raw'] = pred_df[use_col].astype(str).values

# Confidence
conf_col = None
for name in ['conf_score', 'confidence', 'confidence_score', 'prob']:
    if name in pred_df.columns:
        conf_col = name
        break

if conf_col:
    adata.obs['celltypist_confidence'] = pred_df[conf_col].values
else:
    adata.obs['celltypist_confidence'] = 1.0

# Clean up
del adata_celltypist, X_ct, var_ct
gc.collect()

# Summary
celltypist_counts = adata.obs['cell_type_celltypist_raw'].value_counts()
print(f"\n✓ CellTypist results:")
print(f"  Unique types: {len(celltypist_counts)}")
print(f"  Mean confidence: {adata.obs['celltypist_confidence'].mean():.3f}")

print(f"\nTop 10 types:")
for ct, count in celltypist_counts.head(10).items():
    print(f"  {ct}: {count:,}")

# Save counts
save_celltype_counts(
    {'celltypist_counts_raw': celltypist_counts},
    OUTPUT_DIR
)

# ============================================================================
# STAGE 7: PREPARE LABELS FOR scANVI
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 7: PREPARE scANVI LABELS")
print("="*80)

print(f"\n⭐ Creating scANVI training labels...")
print(f"  Low confidence threshold: {CELLTYPIST_MIN_CONFIDENCE}")
print(f"  Rare type threshold: {MIN_CELLS_PER_TYPE} cells")

# Filter by confidence
ct_raw = adata.obs['cell_type_celltypist_raw'].astype(str)
ct_lc = ct_raw.where(
    adata.obs['celltypist_confidence'] >= CELLTYPIST_MIN_CONFIDENCE,
    'Unknown'
)
n_low_conf = (ct_lc == 'Unknown').sum()
print(f"\n  Low confidence → Unknown: {n_low_conf:,} cells")

# Filter rare types
ct_filt = merge_rare_types(ct_lc, min_cells=MIN_CELLS_PER_TYPE, other='Unknown')
adata.obs['cell_type_celltypist_filt'] = ct_filt

# Convert to categorical
labels = ct_filt.astype('category')
if 'Unknown' not in labels.cat.categories:
    labels = labels.cat.add_categories(['Unknown'])

adata.obs['labels_for_scanvi'] = labels

# Statistics
n_unknown = int((labels == 'Unknown').sum())
n_labeled = int(adata.n_obs - n_unknown)
n_unique = int(labels.nunique()) - (1 if 'Unknown' in labels.cat.categories else 0)

print(f"\n✓ scANVI labels:")
print(f"  Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  Unknown: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(f"  Unique types: {n_unique}")

# Save filtered counts
save_celltype_counts(
    {'celltypist_counts_filt': adata.obs['cell_type_celltypist_filt'].value_counts()},
    OUTPUT_DIR
)

# ============================================================================
# STAGE 8: BUILD adata_model FOR TRAINING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 8: BUILD adata_model (HVG-ONLY)")
print("="*80)

print(f"\n⭐⭐⭐ P0-FIX: Creating dedicated adata_model...")

X_hvg = adata.layers['counts'][:, hvg_mask].copy()
obs_hvg = adata.obs[[BATCH_KEY, 'labels_for_scanvi']].copy()
var_hvg = adata.var.loc[hvg_mask, []].copy()

adata_model = sc.AnnData(X=X_hvg, obs=obs_hvg, var=var_hvg)
adata_model.layers['counts'] = adata_model.X.copy()

print(f"\n✓ adata_model:")
print(f"  Shape: {adata_model.n_obs:,} × {adata_model.n_vars:,}")
print(f"  obs: {list(adata_model.obs.columns)}")

del X_hvg, obs_hvg, var_hvg
gc.collect()

# ============================================================================
# STAGE 9: scVI TRAINING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 9: scVI BATCH CORRECTION")
print("="*80)

# Setup
print(f"\n⭐ Setting up scVI...")
scvi.model.SCVI.setup_anndata(
    adata_model,
    layer='counts',
    batch_key=BATCH_KEY
)
print(f"✓ scVI setup complete")

# Train
print(f"\n⭐ Training scVI model...")
vae = scvi.model.SCVI(
    adata_model,
    n_latent=N_LATENT,
    n_layers=N_LAYERS,
    gene_likelihood="nb"
)

n_params = sum(p.numel() for p in vae.module.parameters())
print(f"  Parameters: {n_params:,}")

print(f"\n{'='*80}")
print(f"Training scVI...")
print(f"{'='*80}\n")

start_time = time.time()
train_kwargs = {
    'max_epochs': SCVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'accelerator': accelerator,
    'devices': devices,
    'plan_kwargs': {'lr': LEARNING_RATE},
}
if EARLY_STOPPING:
    train_kwargs['early_stopping'] = True
    train_kwargs['early_stopping_patience'] = EARLY_STOPPING_PATIENCE

vae.train(**train_kwargs)
elapsed = time.time() - start_time
print(f"\n✓ Training complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# Save model
scvi_model_dir = MODEL_DIR / "scvi_model"
print(f"\n⭐ Saving scVI model: {scvi_model_dir}")
vae.save(scvi_model_dir, overwrite=True)

# Get latent representation
print(f"\n⭐ Writing X_scvi to main adata (index-aligned)...")
latent_scvi = pd.DataFrame(
    vae.get_latent_representation(),
    index=adata_model.obs_names
)
adata.obsm['X_scvi'] = latent_scvi.reindex(adata.obs_names).to_numpy()
print(f"✓ X_scvi: {adata.obsm['X_scvi'].shape}")

# ============================================================================
# STAGE 10: scANVI TRAINING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 10: scANVI FINE-TUNING")
print("="*80)

print(f"\n⭐ P0-FIX: Initializing scANVI from scVI...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata_model,
    labels_key='labels_for_scanvi',
    unlabeled_category='Unknown'
)
print(f"✓ scANVI created")

# Train
print(f"\n{'='*80}")
print(f"Training scANVI...")
print(f"{'='*80}\n")

start_time = time.time()
train_kwargs = {
    'max_epochs': SCANVI_MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'train_size': 0.9,
    'accelerator': accelerator,
    'devices': devices,
    'plan_kwargs': {'lr': LEARNING_RATE},
}
if EARLY_STOPPING:
    train_kwargs['early_stopping'] = True
    train_kwargs['early_stopping_patience'] = EARLY_STOPPING_PATIENCE

lvae.train(**train_kwargs)
elapsed = time.time() - start_time
print(f"\n✓ Training complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# Save model
scanvi_model_dir = MODEL_DIR / "scanvi_model"
print(f"\n⭐ Saving scANVI model: {scanvi_model_dir}")
lvae.save(scanvi_model_dir, overwrite=True)

# Get predictions
print(f"\n⭐ Generating scANVI predictions (index-aligned)...")

pred = pd.Series(lvae.predict(), index=adata_model.obs_names)
adata.obs['cell_type_scanvi_raw'] = pred.reindex(adata.obs_names).astype(str).values

probs = np.asarray(lvae.predict(soft=True))
conf = pd.Series(probs.max(axis=1), index=adata_model.obs_names)
adata.obs['scanvi_confidence'] = conf.reindex(adata.obs_names).values

z = pd.DataFrame(lvae.get_latent_representation(), index=adata_model.obs_names)
adata.obsm['X_scanvi'] = z.reindex(adata.obs_names).to_numpy()

# Filter rare types
print(f"\n⭐ Filtering rare types in scANVI predictions...")
adata.obs['cell_type_scanvi_filt'] = merge_rare_types(
    adata.obs['cell_type_scanvi_raw'],
    min_cells=MIN_CELLS_PER_TYPE,
    other='Unknown'
).astype('category')

print(f"\n✓ scANVI predictions:")
print(f"  X_scanvi: {adata.obsm['X_scanvi'].shape}")
print(f"  Unique types (raw): {adata.obs['cell_type_scanvi_raw'].nunique()}")
print(f"  Unique types (filtered): {adata.obs['cell_type_scanvi_filt'].nunique()}")

# Save counts
save_celltype_counts(
    {
        'scanvi_counts_raw': adata.obs['cell_type_scanvi_raw'].value_counts(),
        'scanvi_counts_filt': adata.obs['cell_type_scanvi_filt'].value_counts()
    },
    OUTPUT_DIR
)

# Clean up
del adata_model
gc.collect()

# ============================================================================
# STAGE 11: VALIDATION
# ============================================================================

print(f"\n{'='*80}")
print("VALIDATION")
print("="*80)

# Confidence
mean_conf = adata.obs['scanvi_confidence'].mean()
print(f"\n1. CONFIDENCE")
print(f"  Mean: {mean_conf:.3f}")

# Distribution
final_counts = adata.obs['cell_type_scanvi_filt'].value_counts()
print(f"\n2. FINAL DISTRIBUTION")
print(f"  Cell types: {len(final_counts)}")
print(f"\nTop 15 types:")
for ct, count in final_counts.head(15).items():
    pct = count / adata.n_obs * 100
    print(f"  {ct}: {count:,} ({pct:.2f}%)")

# Agreement
ct_filt = adata.obs['cell_type_celltypist_filt'].astype('object')
scanvi_filt = adata.obs['cell_type_scanvi_filt'].astype('object')
agreement = (ct_filt == scanvi_filt).sum()
agreement_rate = agreement / adata.n_obs * 100

print(f"\n3. AGREEMENT")
print(f"  CellTypist-scANVI: {agreement_rate:.1f}%")

# ============================================================================
# STAGE 12: VISUALIZATION
# ============================================================================

print(f"\n{'='*80}")
print("VISUALIZATION")
print("="*80)

# Create log1p layer
print(f"\n⭐ Creating log1p layer...")
if 'log1p' not in adata.layers:
    from scipy.sparse import issparse
    if issparse(adata.layers['counts']):
        counts_norm = adata.layers['counts'].copy()
        temp_adata = sc.AnnData(X=counts_norm)
        sc.pp.normalize_total(temp_adata, target_sum=1e4)
        sc.pp.log1p(temp_adata)
        adata.layers['log1p'] = temp_adata.X
        del temp_adata, counts_norm
    else:
        adata.layers['log1p'] = adata.layers['counts'].copy()
        temp_adata = sc.AnnData(X=adata.layers['log1p'])
        sc.pp.normalize_total(temp_adata, target_sum=1e4)
        sc.pp.log1p(temp_adata)
        adata.layers['log1p'] = temp_adata.X
        del temp_adata
    print(f"✓ log1p layer created")

adata.X = adata.layers['log1p']

# scVI UMAP
print(f"\n⭐ Computing scVI UMAP...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=15, key_added='neighbors_scvi')
sc.tl.umap(adata, neighbors_key='neighbors_scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print(f"✓ X_umap_scvi: {adata.obsm['X_umap_scvi'].shape}")

# scANVI UMAP
print(f"\n⭐ Computing scANVI UMAP...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15, key_added='neighbors_scanvi')
sc.tl.umap(adata, neighbors_key='neighbors_scanvi')
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
print(f"✓ X_umap_scanvi: {adata.obsm['X_umap_scanvi'].shape}")

# Figure 1: scVI space
print(f"\n⭐ Generating scVI comparison plot...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

sc.pl.embedding(
    adata, basis='umap_scvi', color='cell_type_celltypist_filt',
    ax=axes[0], show=False, title='CellTypist (scVI)', size=2,
    legend_loc='right margin', frameon=False
)

sc.pl.embedding(
    adata, basis='umap_scvi', color='celltypist_confidence',
    ax=axes[1], show=False, title='CellTypist Confidence', size=2,
    cmap='viridis', frameon=False
)

sc.pl.embedding(
    adata, basis='umap_scvi', color=BATCH_KEY,
    ax=axes[2], show=False, title='Batch (scVI)', size=1, frameon=False
)

plt.tight_layout()
plt.savefig(FIG_DIR / 'pure_stromal_scvi.pdf', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: pure_stromal_scvi.pdf")

# Figure 2: scANVI space
print(f"\n⭐ Generating scANVI comparison plot...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

sc.pl.embedding(
    adata, basis='umap_scanvi', color='cell_type_scanvi_filt',
    ax=axes[0], show=False, title='scANVI Final', size=2,
    legend_loc='right margin', frameon=False
)

sc.pl.embedding(
    adata, basis='umap_scanvi', color='scanvi_confidence',
    ax=axes[1], show=False, title='scANVI Confidence', size=2,
    cmap='viridis', frameon=False
)

sc.pl.embedding(
    adata, basis='umap_scanvi', color=BATCH_KEY,
    ax=axes[2], show=False, title='Batch (scANVI)', size=1, frameon=False
)

plt.tight_layout()
plt.savefig(FIG_DIR / 'pure_stromal_scanvi.pdf', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Saved: pure_stromal_scanvi.pdf")

# Figure 3: Before/After filtering
if ORIGINAL_LABEL_KEY in adata.obs.columns:
    print(f"\n⭐ Generating before/after comparison...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    sc.pl.embedding(
        adata, basis='umap_scanvi', color=ORIGINAL_LABEL_KEY,
        ax=axes[0], show=False, title='Original Annotation (Mixed)', size=2,
        legend_loc='right margin', frameon=False
    )
    
    sc.pl.embedding(
        adata, basis='umap_scanvi', color='cell_type_scanvi_filt',
        ax=axes[1], show=False, title='New Annotation (Pure Stromal)', size=2,
        legend_loc='right margin', frameon=False
    )
    
    plt.tight_layout()
    plt.savefig(FIG_DIR / 'before_after_filtering.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved: before_after_filtering.pdf")

# ============================================================================
# STAGE 13: SAVE RESULTS
# ============================================================================

print(f"\n{'='*80}")
print("SAVE RESULTS")
print("="*80)

output_file = OUTPUT_DIR / "adata_stromal_PURE_FINAL.h5ad"

# Add metadata
adata.uns['pipeline_info'] = {
    'version': '1.0-PURE-STROMAL',
    'date': datetime.now().strftime('%Y-%m-%d'),
    'purpose': 'Pure stromal cells after filtering epithelial/immune contaminants',
    'input_file': str(INPUT_FILE),
    'cells_before_filtering': len(final_keep) + len(~final_keep),
    'cells_after_filtering': adata.n_obs,
    'removed_cell_types': removed_types.index.tolist(),
    'kept_patterns': KEEP_PATTERNS,
    'remove_patterns': REMOVE_PATTERNS,
    'batch_key': BATCH_KEY,
    'n_hvg': int(adata.var['highly_variable'].sum()),
    'hvg_method': hvg_method,
    'celltypist_model': str(CELLTYPIST_MODEL),
    'critical_fixes': [
        'P0: Dedicated adata_model for scVI/scANVI',
        'P0: Index-aligned writeback',
        'P0: predict(soft=True) ndarray handling',
        'P1: CellTypist majority voting',
        'P1: Model features validation',
        'P2: Memory optimization'
    ]
}

adata.write_h5ad(output_file, compression='gzip')
size_gb = output_file.stat().st_size / 1e9
print(f"\n✓ Saved: {output_file}")
print(f"  Size: {size_gb:.2f} GB")

# Export annotations
annotations = adata.obs[[
    ORIGINAL_LABEL_KEY,
    'cell_type_celltypist_raw',
    'cell_type_celltypist_filt',
    'celltypist_confidence',
    'cell_type_scanvi_raw',
    'cell_type_scanvi_filt',
    'scanvi_confidence',
    BATCH_KEY
]].copy()
annotations.to_csv(OUTPUT_DIR / "annotations_complete.csv")
print(f"✓ Annotations exported")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print("🎉 PURE STROMAL CELL RE-ANALYSIS COMPLETE")
print("="*80)

print(f"\n📊 Summary:")
print(f"  Original cells: {len(final_keep) + len(~final_keep):,}")
print(f"  Filtered cells: {adata.n_obs:,}")
print(f"  Removal rate: {len(~final_keep)/(len(final_keep)+len(~final_keep))*100:.1f}%")
print(f"  Genes (full): {adata.n_vars:,}")
print(f"  Genes (HVG): {adata.var['highly_variable'].sum()}")
print(f"  Batches: {adata.obs[BATCH_KEY].nunique()}")

print(f"\n⭐ Cell Types Removed:")
for ct in removed_types.index[:10]:
    print(f"  ✗ {ct}")
if len(removed_types) > 10:
    print(f"  ... and {len(removed_types)-10} more")

print(f"\n⭐ Final Cell Types (scANVI filtered):")
for ct, count in final_counts.head(15).items():
    pct = count / adata.n_obs * 100
    print(f"  ✓ {ct}: {count:,} ({pct:.2f}%)")

print(f"\n📁 Output Files:")
print(f"  h5ad: {output_file.name}")
print(f"  Annotations: annotations_complete.csv")
print(f"  Filtering stats: filtering_statistics.csv")
print(f"  Cell type counts: *_counts_*.csv")
print(f"  Models: scvi_model/, scanvi_model/")

print(f"\n📈 Figures:")
print(f"  pure_stromal_scvi.pdf - scVI UMAP")
print(f"  pure_stromal_scanvi.pdf - scANVI UMAP")
print(f"  before_after_filtering.pdf - Comparison")

print(f"\n" + "="*80)
print("NEXT STEPS:")
print("="*80)
print("1. Validate cell type annotations in UMAP plots")
print("2. Check endothelial subtype separation (arterial/venous/lymphatic)")
print("3. Examine fibroblast heterogeneity and myofibroblast markers")
print("4. Compare smooth muscle and pericyte populations")
print("5. Assess Schwann cell markers if present")
print("6. Downstream: marker gene analysis, trajectory, differential expression")
print("="*80)

print(f"\n✅ Pipeline complete! Output: {OUTPUT_DIR}")
