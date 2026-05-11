#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Epithelial Cell: scVI → CellTypist → scANVI Pipeline v3.4
==========================================================

Production-ready pipeline with all critical fixes:
- ✅ Robust gene name conversion (local-first, HPC-compatible)
- ✅ Proper HVG subset handling (explicit counts layer management)
- ✅ CellTypist feature matching with base symbols (no suffix pollution)
- ✅ Index-aligned result writing (no order assumptions)
- ✅ Pre-trained model loading with correct HVG coordination

Critical Fixes in v3.4 (from review feedback):
1. HVG subset: Explicit counts layer management post-subsetting
2. Gene conversion: Local-first with mygene fallback (HPC-safe)
3. Symbol deduplication: Base symbols for CellTypist, no -1/-2 suffixes
4. Result alignment: Index-based merge (not position-based)
5. HVG coordination: Load HVG list with pre-trained model

Input: Raw epithelial data with counts layer
Output: adata_epithelial_FINAL.h5ad + reusable models

Author: r2end
Date: 2024-12-11  
Version: 3.4 - Production fixes for HPC/large-scale deployment
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Core libraries
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import time

# scvi-tools
import scvi
import torch

# CellTypist  
import celltypist
from celltypist import models

# Optional: mygene (with fallback)
try:
    import mygene
    MYGENE_AVAILABLE = True
except ImportError:
    MYGENE_AVAILABLE = False
    print("⚠️  mygene not available, will use local gene names only")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output paths
INPUT_FILE = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1212/epithelial_scanvi_v3_4"

# ⭐ Pre-trained model loading
PRETRAINED_SCVI_MODEL = '/home/h2048/data/py/1206/bbknn_celltype_analysis/Epithelial/scanvi_results/models/scvi_model'  # Set path to reuse model, or None to train
# Example: "/path/to/previous_run/models/scvi_model"

# Model hyperparameters
N_LATENT = 100
N_LAYERS = 3
N_HVG = 4000

# Training parameters
SCVI_MAX_EPOCHS = 800
SCANVI_MAX_EPOCHS = 600
BATCH_SIZE = 2048
LEARNING_RATE = 1e-3
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 50

# CellTypist
CELLTYPIST_MODEL = '/home/h2048/data/source/reference/celltypist_models/Human_Lung_Atlas.pkl'
CELLTYPIST_MAJORITY_VOTING = True
CELLTYPIST_MIN_CONFIDENCE = 0.5

# Column names
BATCH_KEY = 'dataset'
CELLTYPIST_LABEL_KEY = 'cell_type_celltypist'
SCANVI_LABEL_KEY = 'cell_type_scanvi'

# Random seed
RANDOM_SEED = 42

print("="*80)
print("EPITHELIAL PIPELINE v3.4 (PRODUCTION)")
print("="*80)
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Pre-trained scVI: {PRETRAINED_SCVI_MODEL or 'None (train fresh)'}")

# Create directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
MODEL_DIR = Path(OUTPUT_DIR) / "models"
MODEL_DIR.mkdir(exist_ok=True)

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

print(f"\n{'='*80}")
print("ENVIRONMENT CHECK")
print("="*80)

# Set seeds
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

print(f"\nPackage versions:")
print(f"  scanpy: {sc.__version__}")
print(f"  scvi-tools: {scvi.__version__}")
print(f"  torch: {torch.__version__}")
print(f"  celltypist: {celltypist.__version__}")
print(f"  mygene: {'available' if MYGENE_AVAILABLE else 'not available'}")

# GPU check
gpu_available = torch.cuda.is_available()
if gpu_available:
    gpu_name = torch.cuda.get_device_name(0)
    print(f"\n✓ GPU: {gpu_name}")
    accelerator, devices = 'gpu', 'auto'
    scvi.settings.dl_num_workers = 0
else:
    print(f"\n⚠️  No GPU, using CPU")
    accelerator, devices = 'cpu', 'auto'

# ============================================================================
# LOAD DATA
# ============================================================================

print(f"\n{'='*80}")
print("DATA LOADING")
print("="*80)

adata = sc.read_h5ad(INPUT_FILE)
print(f"✓ Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Validate
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found")
if 'counts' not in adata.layers:
    raise ValueError("adata.layers['counts'] not found")

print(f"✓ Batch key '{BATCH_KEY}' found ({adata.obs[BATCH_KEY].nunique()} batches)")

# Basic QC
n_genes_before = adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
print(f"✓ QC: {n_genes_before:,} → {adata.n_vars:,} genes")

# ============================================================================
# GENE NAME CONVERSION (LOCAL-FIRST, HPC-SAFE)
# ============================================================================

print(f"\n{'='*80}")
print("GENE NAME CONVERSION (HPC-SAFE)")
print("="*80)

print(f"\nStrategy: Local annotation > mygene online > keep original")

# Check current gene format
sample_gene = str(adata.var_names[0])
print(f"  Current format: {sample_gene}")

# Try local annotation columns first
adata.var['symbol_base'] = None

for col in ['gene_symbol', 'gene_symbols', 'symbol', 'features', 'feature_name']:
    if col in adata.var.columns:
        print(f"  ✓ Found local column: '{col}'")
        adata.var['symbol_base'] = adata.var[col].astype(str).values
        conversion_method = f"local ({col})"
        break
else:
    # Try mygene if available and looks like ENSEMBL
    if MYGENE_AVAILABLE and sample_gene.startswith('ENSG'):
        print(f"  Attempting mygene conversion (ENSEMBL detected)...")
        try:
            mg = mygene.MyGeneInfo()
            results = mg.querymany(
                adata.var_names.tolist(),
                scopes='ensembl.gene',
                fields='symbol',
                species='human',
                as_dataframe=False
            )
            
            symbol_map = {}
            for res in results:
                if 'symbol' in res:
                    symbol_map[res['query']] = res['symbol']
            
            adata.var['symbol_base'] = [
                symbol_map.get(g, g) for g in adata.var_names
            ]
            
            n_converted = sum([s != g for s, g in zip(adata.var['symbol_base'], adata.var_names)])
            print(f"  ✓ mygene: {n_converted}/{len(adata.var_names)} converted")
            conversion_method = "mygene"
            
        except Exception as e:
            print(f"  ⚠️  mygene failed ({str(e)[:50]}), keeping original")
            adata.var['symbol_base'] = adata.var_names.astype(str)
            conversion_method = "original"
    else:
        print(f"  Using original gene names")
        adata.var['symbol_base'] = adata.var_names.astype(str)
        conversion_method = "original"

adata.uns['gene_conversion_method'] = conversion_method
print(f"\n✓ Gene conversion complete: {conversion_method}")

# ============================================================================
# HVG SELECTION (3-TIER FALLBACK)
# ============================================================================

print(f"\n{'='*80}")
print("HVG SELECTION - ROBUST STRATEGY")
print("="*80)

# Check if loading pre-trained model with HVG list
hvg_file = None
if PRETRAINED_SCVI_MODEL:
    # ⭐ FIX: HVG file should be in model's parent directory
    model_parent = Path(PRETRAINED_SCVI_MODEL).parent.parent
    hvg_file = model_parent / "hvg_genes.txt"
    if hvg_file.exists():
        print(f"\n⭐ Loading HVG list from pre-trained model:")
        print(f"  {hvg_file}")
        hvg_genes = pd.read_csv(hvg_file, header=None)[0].tolist()
        
        # Mark HVG in adata.var
        adata.var['highly_variable'] = adata.var['symbol_base'].isin(hvg_genes)
        n_hvg = adata.var['highly_variable'].sum()
        
        print(f"✓ Loaded {n_hvg} HVG from file")
        hvg_method = "loaded_from_model"
    else:
        print(f"  ⚠️  HVG file not found, will compute fresh")
        hvg_file = None

# Compute HVG if not loaded
if hvg_file is None or n_hvg == 0:
    print(f"\nSelecting {N_HVG} highly variable genes...")
    
    # Strategy 1: Batch-aware Seurat v3
    try:
        print(f"  Attempt 1: Batch-aware (seurat_v3)...")
        sc.pp.highly_variable_genes(
            adata, layer='counts', n_top_genes=N_HVG,
            batch_key=BATCH_KEY, flavor='seurat_v3', subset=False
        )
        hvg_method = "batch-aware (seurat_v3)"
        print(f"  ✓ Success")
    
    except Exception as e:
        print(f"  ✗ Failed: {str(e)[:80]}")
        
        # Strategy 2: Non-batch-aware Seurat v3
        try:
            print(f"\n  Attempt 2: Non-batch-aware (seurat_v3)...")
            sc.pp.highly_variable_genes(
                adata, layer='counts', n_top_genes=N_HVG,
                flavor='seurat_v3', subset=False
            )
            hvg_method = "non-batch-aware (seurat_v3)"
            print(f"  ✓ Success")
        
        except Exception as e2:
            print(f"  ✗ Failed: {str(e2)[:80]}")
            
            # Strategy 3: Original Seurat (most stable)
            print(f"\n  Attempt 3: Fallback (seurat)...")
            sc.pp.highly_variable_genes(
                adata, layer='counts', n_top_genes=N_HVG,
                flavor='seurat', subset=False
            )
            hvg_method = "fallback (seurat)"
            print(f"  ✓ Success")
    
    n_hvg = adata.var['highly_variable'].sum()
    
    # Save HVG list for future runs
    hvg_genes = adata.var[adata.var['highly_variable']]['symbol_base'].tolist()
    hvg_file_out = Path(OUTPUT_DIR) / "hvg_genes.txt"
    pd.Series(hvg_genes).to_csv(hvg_file_out, index=False, header=False)
    print(f"  ✓ HVG list saved: {hvg_file_out}")

adata.uns['hvg_method'] = hvg_method
print(f"\n✓ HVG selection complete")
print(f"  Method: {hvg_method}")
print(f"  HVG count: {n_hvg:,}")

# ⭐ CRITICAL: Preserve full genes to .raw BEFORE subsetting
print(f"\n⭐ Preserving full genes to .raw (shared memory)...")
adata.raw = sc.AnnData(
    X=adata.layers["counts"],  # Shared memory
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ Full genes preserved: {adata.raw.n_vars:,}")

# Subset to HVG
print(f"\n⭐ Subsetting to HVG for training...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"✓ Training data: {adata.n_obs:,} cells × {adata.n_vars:,} HVG")

# ⭐ CRITICAL FIX: Ensure counts layer exists post-subsetting
# After subsetting, adata.X is HVG counts, but layers might be lost
print(f"\n⭐ Ensuring counts layer exists post-subsetting...")
if 'counts' not in adata.layers:
    print(f"  Creating adata.layers['counts'] from adata.X")
    adata.layers['counts'] = adata.X.copy()
else:
    print(f"  adata.layers['counts'] already exists")

# Verify counts layer is correct
test_max = adata.layers['counts'][:100, :10].max()
if hasattr(test_max, 'item'):
    test_max = test_max.item()
print(f"  Counts layer check: max={test_max:.1f} (should be >10 for counts)")

# ============================================================================
# STAGE 1: scVI (Train or Load)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 1: scVI BATCH CORRECTION")
print("="*80)

vae = None
scvi_model_dir = MODEL_DIR / "scvi_model"

# Load pre-trained model if specified
if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
    print(f"\n⭐ Loading pre-trained scVI model:")
    print(f"  {PRETRAINED_SCVI_MODEL}")
    
    scvi.model.SCVI.setup_anndata(
        adata, batch_key=BATCH_KEY, layer='counts'
    )
    vae = scvi.model.SCVI.load(PRETRAINED_SCVI_MODEL, adata=adata)
    
    print(f"✓ Model loaded")
    
else:
    print(f"\n⭐ Training new scVI model...")
    
    scvi.model.SCVI.setup_anndata(
        adata, batch_key=BATCH_KEY, layer='counts'
    )
    vae = scvi.model.SCVI(
        adata, n_latent=N_LATENT, n_layers=N_LAYERS,
        gene_likelihood='nb'
    )
    
    print(f"✓ Model created")
    
    # Train
    print(f"\n{'='*80}")
    print(f"Training scVI...")
    print(f"{'='*80}\n")
    
    start_time = time.time()
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'batch_size': BATCH_SIZE,
        'train_size': 0.9,
        'early_stopping': EARLY_STOPPING,
        'accelerator': accelerator,
        'devices': devices,
        'plan_kwargs': {'lr': LEARNING_RATE},
    }
    if EARLY_STOPPING:
        train_kwargs['early_stopping_patience'] = EARLY_STOPPING_PATIENCE
    
    vae.train(**train_kwargs)
    elapsed = time.time() - start_time
    print(f"\n✓ Training complete: {int(elapsed//60)}m {int(elapsed%60)}s")
    
    # Save
    print(f"\nSaving scVI model: {scvi_model_dir}")
    vae.save(scvi_model_dir, overwrite=True)

# Generate latent space
print(f"\nGenerating scVI latent space...")
adata.obsm['X_scvi'] = vae.get_latent_representation()
print(f"✓ X_scvi saved: {adata.obsm['X_scvi'].shape}")

print(f"\n✓ Stage 1 complete")

# ============================================================================
# STAGE 2: CellTypist ANNOTATION (WITH PROPER FEATURE MATCHING)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 2: CellTypist ANNOTATION")
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

# ⭐ CRITICAL: Build CellTypist input with proper feature matching
print(f"\n⭐ Building CellTypist input with base symbol matching...")

# Get model features
model_features = pd.Index(getattr(model, 'features', []))
print(f"  Model expects {len(model_features)} features")

# Use base symbols (no -1/-2 suffixes) from .raw
sym = pd.Index(adata.raw.var['symbol_base']).astype(str)

# Find overlap: genes in model AND not duplicated in our data
keep_mask = sym.isin(model_features) & ~sym.duplicated(keep='first')
n_overlap = keep_mask.sum()

print(f"  Our data: {len(sym)} genes")
print(f"  Overlap: {n_overlap} genes")
print(f"  Coverage: {n_overlap/len(model_features)*100:.1f}%")

if n_overlap < 1000:
    print(f"\n  ⚠️  WARNING: Low overlap ({n_overlap} genes)")
    print(f"  Consider checking gene name format")

# Build CellTypist AnnData (only overlapping genes, deduplicated)
print(f"\n  Building CellTypist AnnData...")
X_ct = adata.raw.X[:, keep_mask].copy()
var_ct = pd.DataFrame(index=sym[keep_mask])

adata_celltypist = sc.AnnData(
    X=X_ct,
    obs=pd.DataFrame(index=adata.obs_names),
    var=var_ct
)

print(f"✓ CellTypist input: {adata_celltypist.n_obs:,} cells × {adata_celltypist.n_vars:,} genes")

# Normalize
print(f"\n  Normalizing...")
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
sc.pp.log1p(adata_celltypist)

# Run CellTypist
print(f"\nRunning CellTypist prediction...")
start_time = time.time()

predictions = celltypist.annotate(
    adata_celltypist,
    model=model,
    majority_voting=CELLTYPIST_MAJORITY_VOTING
)

elapsed = time.time() - start_time
print(f"✓ Prediction complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# ⭐ CRITICAL: Index-aligned result writing
print(f"\n⭐ Writing results with index alignment...")
pred_df = predictions.predicted_labels
pred_df = pred_df.reindex(adata.obs_names)  # Force alignment

adata.obs[CELLTYPIST_LABEL_KEY] = pred_df['predicted_labels'].astype(str).values

if CELLTYPIST_MAJORITY_VOTING and 'majority_voting' in pred_df.columns:
    adata.obs[f'{CELLTYPIST_LABEL_KEY}_majority'] = pred_df['majority_voting'].astype(str).values

# Confidence detection
conf_col = None
for name in ['conf_score', 'confidence', 'confidence_score', 'prob']:
    if name in pred_df.columns:
        conf_col = name
        break

if conf_col:
    adata.obs['celltypist_confidence'] = pred_df[conf_col].values
    print(f"  ✓ Confidence: {conf_col}")
else:
    adata.obs['celltypist_confidence'] = 1.0
    print(f"  ⚠️  No confidence, using 1.0")

# Clean up
del adata_celltypist, X_ct, var_ct
import gc
gc.collect()

# Summary
celltypist_counts = adata.obs[CELLTYPIST_LABEL_KEY].value_counts()
mean_conf = adata.obs['celltypist_confidence'].mean()

print(f"\n✓ CellTypist results:")
print(f"  Unique types: {len(celltypist_counts)}")
print(f"  Mean confidence: {mean_conf:.3f}")

print(f"\n✓ Stage 2 complete")

# ============================================================================
# STAGE 3: PREPARE LABELS FOR scANVI
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 3: PREPARE LABELS FOR scANVI")
print("="*80)

print(f"\nLabel strategy:")
print(f"  Source: CellTypist predictions")
print(f"  Low confidence (<{CELLTYPIST_MIN_CONFIDENCE}): mark as 'Unknown'")

# Create labels
adata.obs['labels_for_scanvi'] = adata.obs[CELLTYPIST_LABEL_KEY].copy()

# Mark low confidence
low_conf_mask = adata.obs['celltypist_confidence'] < CELLTYPIST_MIN_CONFIDENCE
n_low_conf = low_conf_mask.sum()

if n_low_conf > 0:
    adata.obs.loc[low_conf_mask, 'labels_for_scanvi'] = 'Unknown'

# Convert to categorical with 'Unknown'
labels = adata.obs['labels_for_scanvi'].astype('category')
if 'Unknown' not in labels.cat.categories:
    labels = labels.cat.add_categories(['Unknown'])
adata.obs['labels_for_scanvi'] = labels

n_unknown = (adata.obs['labels_for_scanvi'] == 'Unknown').sum()
n_labeled = adata.n_obs - n_unknown

print(f"\n✓ Labels ready:")
print(f"  Total: {adata.n_obs:,}")
print(f"  Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"  Unknown: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")

# ============================================================================
# STAGE 4: scANVI FINE-TUNING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 4: scANVI FINE-TUNING")
print("="*80)

print(f"\nInitializing scANVI from scVI...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category='Unknown',
    labels_key='labels_for_scanvi'
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

# Save
scanvi_model_dir = MODEL_DIR / "scanvi_model"
print(f"\nSaving scANVI model: {scanvi_model_dir}")
lvae.save(scanvi_model_dir, overwrite=True)

# Generate predictions
print(f"\nGenerating final predictions...")
adata.obs[SCANVI_LABEL_KEY] = lvae.predict()

# Confidence (with cleanup)
pred_probs = lvae.predict(soft=True)
adata.obs['scanvi_confidence'] = pred_probs.max(axis=1).values
del pred_probs  # Clean up large matrix
gc.collect()

adata.obsm['X_scanvi'] = lvae.get_latent_representation()

print(f"✓ Predictions saved")
print(f"\n✓ Stage 4 complete")

# ============================================================================
# VALIDATION & SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print("VALIDATION")
print("="*80)

# Confidence
overall_mean_conf = adata.obs['scanvi_confidence'].mean()
print(f"\n1. CONFIDENCE")
print(f"  Mean: {overall_mean_conf:.3f}")

# Distribution
final_counts = adata.obs[SCANVI_LABEL_KEY].value_counts()
print(f"\n2. FINAL DISTRIBUTION")
print(f"  Cell types: {len(final_counts)}")

# Agreement
agreements = (adata.obs[CELLTYPIST_LABEL_KEY] == adata.obs[SCANVI_LABEL_KEY]).sum()
agreement_rate = agreements / adata.n_obs * 100
print(f"\n3. AGREEMENT")
print(f"  CellTypist-scANVI: {agreement_rate:.1f}%")

# ============================================================================
# UMAP & SAVE
# ============================================================================

print(f"\n{'='*80}")
print("UMAP & SAVE")
print("="*80)

print(f"\nComputing UMAP...")
sc.pp.neighbors(adata, use_rep='X_scanvi', n_neighbors=15)
sc.tl.umap(adata)
print(f"✓ UMAP computed")

# Save
final_file = Path(OUTPUT_DIR) / "adata_epithelial_FINAL.h5ad"
print(f"\nSaving final dataset: {final_file}")

adata.uns['pipeline_params'] = {
    'version': '3.4_production',
    'n_latent': N_LATENT,
    'n_hvg': N_HVG,
    'gene_conversion': adata.uns['gene_conversion_method'],
    'hvg_method': hvg_method,
    'celltypist_model': CELLTYPIST_MODEL,
    'random_seed': RANDOM_SEED
}

adata.write_h5ad(final_file, compression='gzip')
file_size = final_file.stat().st_size / 1e9
print(f"✓ Saved ({file_size:.2f} GB)")

# README
readme = f"""# Epithelial Pipeline v3.4 Results

## Key Improvements (v3.4)
- ✅ HPC-safe gene conversion (local-first)
- ✅ Proper CellTypist feature matching (base symbols)
- ✅ Index-aligned result writing
- ✅ Explicit counts layer management post-HVG
- ✅ Pre-trained model HVG coordination

## Main Output
**{final_file.name}**

### Use This Column
- `adata.obs['{SCANVI_LABEL_KEY}']` - Final cell type (scANVI refined)

### Statistics
- Cells: {adata.n_obs:,}
- HVG: {adata.n_vars:,}
- Full genes (.raw): {adata.raw.n_vars:,}
- Cell types: {len(final_counts)}
- Mean confidence: {overall_mean_conf:.3f}
- Agreement: {agreement_rate:.1f}%

### Gene Conversion
- Method: {adata.uns['gene_conversion_method']}
- HVG selection: {hvg_method}

### Models (Reusable)
- scVI: `models/scvi_model/`
- scANVI: `models/scanvi_model/`
- HVG list: `hvg_genes.txt`

---
Pipeline: v3.4 Production
Generated: {pd.Timestamp.now()}
"""

with open(Path(OUTPUT_DIR) / 'README.md', 'w') as f:
    f.write(readme)

print(f"✓ README created")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"🎉 PIPELINE COMPLETE!")
print(f"{'='*80}")
print(f"\n📁 Output: {OUTPUT_DIR}")
print(f"⭐ Main file: {final_file.name}")
print(f"\n📊 Summary:")
print(f"  - Cells: {adata.n_obs:,}")
print(f"  - Cell types: {len(final_counts)}")
print(f"  - Confidence: {overall_mean_conf:.3f}")
print(f"\n💡 Use adata.obs['{SCANVI_LABEL_KEY}'] for annotations")
print(f"\n✅ Ready for downstream analysis!")
print(f"{'='*80}")
