#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
==============================================================================
scANVI Reference Mapping via scArches - FIXED v1.2.1 PRODUCTION
==============================================================================

CRITICAL FIXES in v1.2.1:
==========================
1. ✅ Calculate ALL covariates matching training EXACTLY
2. ✅ Proper covariate timing (BEFORE prepare_query_anndata)
3. ✅ Tissue_key handling (categorical type, proper validation)
4. ✅ Metadata cleaning preserves ALL covariates
5. ✅ Shared memory for .X (no unnecessary copies)
6. ✅ UMAP parameters match training (n_neighbors=30)
7. ✅ Enhanced validation at each critical step
8. ✅ Clear error messages with debugging info

Version: 1.2.1 (PRODUCTION-READY)
Author: r2end
Date: 2025-01-27

Pipeline Flow:
==============
Step 0: Configuration & Path Validation
Step 1: Load and Validate Query Data
  ├─ Data format detection (counts vs log1p)
  ├─ Gene naming validation
  └─ Structure validation (layers, raw)
Step 1.5: Calculate Required Covariates (CRITICAL)
  ├─ MT% recalculation (handle zeros)
  ├─ Stress score (normalized [0,1])
  ├─ Cell cycle scores (S, G2M)
  └─ Tissue key validation/creation
Step 2: Gene Alignment (HVG overlap check)
Step 3: Prepare Query (scVI prepare_query_anndata)
Step 3.5: Metadata Cleaning (PRESERVE COVARIATES)
Step 4: Load Query into Reference Model
Step 5: Fine-tune on Query Data
Step 6: Extract Results (latent, predictions, confidence)
Step 7: UMAP Projection
Step 8: Quality Control Summary
Step 9: Restore Metadata and Save
Step 10: Visualization
Step 11: Export Summary Report

Requirements:
=============
pip install "scvi-tools>=1.1.4" scanpy numpy pandas scipy matplotlib seaborn

Reference:
==========
Training script: allcells_scvi_celltypist_scanvi_pipeline_20251222_v2_3.py
QUICK_REFERENCE_MEMORY: v2.12
"""

# ==============================================================================
# IMPORTS
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
print("scArches scANVI Reference Mapping - FIXED v1.2.1 PRODUCTION")
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
# CONFIGURATION SECTION (⚠️ MODIFY THESE PATHS)
# ==============================================================================

print("\n" + "=" * 70)
print("CONFIGURATION")
print("=" * 70)

# ===== Input/Output Paths =====
SCANVI_REF_DIR = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX/models/scanvi_existing_model"
HVG_FILE = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX/models/scvi_model/hvg_genes.txt"

# ⚠️ UPDATE THESE PATHS:
QUERY_H5AD = "/home/h2048/data/source/polyp/polyp_obj_updated_20260113.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2"

# ===== Keys (Must Match Reference Training) =====
# Reference: Training script Lines 89-90
BATCH_KEY = "sample"                 # Categorical covariate
TISSUE_KEY = "tissue"                # Categorical covariate (CRITICAL)
LABELS_KEY = "cell_type"             # For scANVI
UNLABELED_CATEGORY = "Unknown"

# ===== Required Covariates =====
# CRITICAL: Must match training EXACTLY
# Reference: Training script Lines 91-98
REQUIRED_COVARIATES = [
    BATCH_KEY,           # Categorical
    TISSUE_KEY,          # Categorical
    'pct_counts_mt',     # Continuous
    'stress_score',      # Continuous
    'S_score',           # Continuous
    'G2M_score'          # Continuous
]

# Categorical vs Continuous separation (for validation)
CATEGORICAL_COVARIATES = [BATCH_KEY, TISSUE_KEY]
CONTINUOUS_COVARIATES = ['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']

# ===== Stress Signature Genes =====
# Reference: Training script Lines 153-164
STRESS_SIGNATURE_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1",
    "CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6",
    "CXXC1","DCP1A","DCP2","DCTN1","DDIT4","DDX10","DKC1","DNAJA4","DNAJB9",
    "DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1","EIF4A1","EIF4A2","EIF4A3",
    "EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A","EXOC2","EXOSC1","EXOSC10",
    "EXOSC2","EXOSC4","EXOSC5","EXOSC9","FKBP14","FUS","GEMIN4","GOSR2","H2AX",
    "HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IARS1","IFIT1","IGFBP1","IMP3",
    "KDELR3","KHSRP","KIF5B","LSM1","LSM4","MTHFD2","NFYA","NFYB","NHP2","NOLC1",
    "NOP14","NOP56","NPM1","NABP1","PAIP1","PARN","PDIA5","PDIA6","POP4","PREB",
    "PSAT1","RPS14","RRP9","SDAD1","SEC11A","SEC31A","SERP1","SHC1","MTREX",
    "SLC1A4","SLC30A5","SLC7A5","SPCS1","SPCS3","SRPRA","SRPRB","SSR1","STC2",
    "TARS1","TATDN2","TSPYL2","SKIC3","TUBB2A","VEGFA","WFS1","WIPI1","XBP1",
    "XPOT","YIF1A","YWHAZ","ZBTB17"
]

# ===== Cell Cycle Genes =====
# Reference: Training script Lines 167-188
S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
    'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
    'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
    'USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8'
]

G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
    'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B',
    'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2',
    'DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN',
    'LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA'
]

# ===== scArches Training Parameters =====
# Reference: Training script Lines 224-226
MAX_EPOCHS = 200
WEIGHT_DECAY = 0.0                   # CRITICAL: Must be 0 for scArches
LEARNING_RATE = 5e-4
BATCH_SIZE = 256
EARLY_STOPPING_PATIENCE = 30
CHECK_VAL_EVERY_N_EPOCH = 10

# ===== Label Transfer Confidence =====
CONFIDENCE_THRESHOLD = 0.5

# ===== UMAP Parameters =====
# Reference: Training script v2.5.2 PRODUCTION Lines 227-229
# CRITICAL: Must match training exactly for comparable UMAP projection
UMAP_N_NEIGHBORS = 30     # Standard for large datasets
UMAP_MIN_DIST = 0.5       # Increased from default 0.3
UMAP_SPREAD = 1.0         # Standard

# ===== Quality Control =====
APPLY_QC = False
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

print("\nConfiguration loaded:")
print(f"  Reference model: {SCANVI_REF_DIR}")
print(f"  Query data: {QUERY_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  Batch key: {BATCH_KEY}")
print(f"  Tissue key: {TISSUE_KEY}")
print(f"  Required covariates: {len(REQUIRED_COVARIATES)}")


# ==============================================================================
# PATH VALIDATION
# ==============================================================================

print("\n" + "=" * 70)
print("Path Validation")
print("=" * 70)

validation_failed = False

if not Path(SCANVI_REF_DIR).exists():
    print(f"\n❌ ERROR: Reference model not found")
    print(f"   Path: {SCANVI_REF_DIR}")
    validation_failed = True
else:
    if not (Path(SCANVI_REF_DIR) / "model.pt").exists():
        print(f"\n❌ ERROR: model.pt not found in reference directory")
        validation_failed = True
    else:
        print(f"✓ Reference model found")

if not Path(HVG_FILE).exists():
    print(f"\n❌ ERROR: HVG file not found")
    print(f"   Path: {HVG_FILE}")
    validation_failed = True
else:
    print(f"✓ HVG file found")

if not Path(QUERY_H5AD).exists():
    print(f"\n❌ ERROR: Query data not found")
    print(f"   Path: {QUERY_H5AD}")
    validation_failed = True
else:
    print(f"✓ Query data found")

if validation_failed:
    print(f"\n❌ Path validation failed. Please check configuration.")
    sys.exit(1)

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)
print(f"✓ Output directory created: {output_dir}")

# Save configuration
config = {
    'pipeline': 'scArches_scANVI_mapping_FIXED',
    'version': '1.2.1_PRODUCTION',
    'reference_model': str(SCANVI_REF_DIR),
    'hvg_file': str(HVG_FILE),
    'query_h5ad': str(QUERY_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'tissue_key': TISSUE_KEY,
    'labels_key': LABELS_KEY,
    'required_covariates': REQUIRED_COVARIATES,
    'categorical_covariates': CATEGORICAL_COVARIATES,
    'continuous_covariates': CONTINUOUS_COVARIATES,
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
print("Step 1: Load and Validate Query Data")
print("=" * 70)

print(f"\nLoading: {QUERY_H5AD}")
adata_query = sc.read_h5ad(QUERY_H5AD)

print(f"\nLoaded successfully:")
print(f"  Shape: {adata_query.n_obs:,} cells × {adata_query.n_vars:,} genes")
print(f"  Obs columns: {len(adata_query.obs.columns)}")
print(f"  Layers: {list(adata_query.layers.keys())}")


# ==============================================================================
# Data Format Detection
# ==============================================================================

print("\n" + "-" * 70)
print("Data Format Detection")
print("-" * 70)

def detect_data_type(adata):
    """
    Detect if .X contains raw counts or normalized data.
    Reference: QUICK_REFERENCE_MEMORY v2.12
    """
    X_data = adata.X
    if issparse(X_data):
        max_val = X_data.max()
        mean_val = X_data.mean()
    else:
        max_val = X_data.max()
        mean_val = X_data.mean()
    
    # Decision rules
    if max_val > 100 and mean_val < 50:
        return 'counts', max_val, mean_val
    elif max_val < 15 and mean_val < 5:
        return 'log1p', max_val, mean_val
    else:
        return 'unknown', max_val, mean_val

data_type, max_val, mean_val = detect_data_type(adata_query)

print(f"\n.X Statistics:")
print(f"  Max: {max_val:.2f}")
print(f"  Mean: {mean_val:.2f}")
print(f"  Detected: {data_type.upper()}")

if data_type != 'counts':
    print(f"\n❌ ERROR: .X must contain raw counts for scArches")
    print(f"   Current state: {data_type}")
    print(f"\n   Options:")
    print(f"   1. Re-load from raw counts file")
    print(f"   2. Check if counts exist in .raw or .layers")
    sys.exit(1)

print(f"✓ Confirmed: .X contains raw counts")

# Ensure layers['counts'] exists (shared memory)
if 'counts' not in adata_query.layers:
    print(f"\n  Creating layers['counts'] from .X...")
    adata_query.layers['counts'] = adata_query.X  # Shared memory
    print(f"  ✓ layers['counts'] created (0 memory cost)")
else:
    print(f"✓ layers['counts'] exists")

# Ensure .raw exists (shared memory)
if not hasattr(adata_query, 'raw') or adata_query.raw is None:
    print(f"\n  Creating .raw from layers['counts']...")
    adata_query.raw = sc.AnnData(
        X=adata_query.layers['counts'],  # Shared memory
        obs=adata_query.obs.copy(),
        var=adata_query.var.copy()
    )
    print(f"  ✓ .raw created (0 memory cost)")
else:
    print(f"✓ .raw exists")


# ==============================================================================
# Gene Naming Validation
# ==============================================================================

print("\n" + "-" * 70)
print("Gene Naming Validation")
print("-" * 70)

sample_genes = adata_query.var_names[:5].tolist()
print(f"\nSample genes: {sample_genes}")

# Check for ENSEMBL IDs
if all(str(g).startswith('ENSG') for g in sample_genes):
    print(f"\n❌ ERROR: Query uses ENSEMBL IDs")
    print(f"   Reference uses gene symbols")
    print(f"\n   Solution: Convert ENSEMBL to symbols using mygene")
    sys.exit(1)
else:
    print(f"✓ Gene symbols detected")

# Check for duplicates
if adata_query.var_names.duplicated().any():
    n_dup = adata_query.var_names.duplicated().sum()
    print(f"\n  Warning: {n_dup} duplicate gene names")
    print(f"  Making unique...")
    adata_query.var_names_make_unique()
    print(f"  ✓ Gene names unique")

print(f"\n✓ Gene validation complete")


# ==============================================================================
# Step 1.5: Calculate Required Covariates
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1.5: Calculate Required Covariates (CRITICAL)")
print("=" * 70)

print("\n⚠️ TIMING IS CRITICAL:")
print("  Covariates MUST be calculated BEFORE prepare_query_anndata")
print("  Reason: prepare_query_anndata reorders/pads genes")
print("  If we calculate after, gene indices will be wrong!")

print("\nThis step ensures query data has EXACT covariates as training:")
print(f"  Categorical: {CATEGORICAL_COVARIATES}")
print(f"  Continuous: {CONTINUOUS_COVARIATES}")

# Check existing
print("\nChecking existing covariates...")
existing = [c for c in REQUIRED_COVARIATES if c in adata_query.obs.columns]
missing = [c for c in REQUIRED_COVARIATES if c not in adata_query.obs.columns]

if existing:
    print(f"  Existing: {existing}")
if missing:
    print(f"  Missing: {missing}")

print(f"\n  Strategy: Recalculate ALL for consistency")


# --- 1. MT% Calculation ---
print("\n" + "-" * 70)
print("1. MT% Calculation")
print("-" * 70)

mt_genes = (adata_query.var_names.str.startswith('MT-') | 
            adata_query.var_names.str.startswith('Mt-') | 
            adata_query.var_names.str.startswith('mt-'))
n_mt_genes = mt_genes.sum()
print(f"  MT genes: {n_mt_genes}")

if n_mt_genes > 0:
    counts = adata_query.layers['counts']
    if issparse(counts):
        total_counts = np.array(counts.sum(axis=1)).flatten()
        mt_counts = np.array(counts[:, mt_genes].sum(axis=1)).flatten()
    else:
        total_counts = counts.sum(axis=1)
        mt_counts = counts[:, mt_genes].sum(axis=1)
    
    # Handle zero total counts
    mt_pct = np.zeros_like(total_counts, dtype=np.float32)
    nonzero_mask = total_counts > 0
    mt_pct[nonzero_mask] = (mt_counts[nonzero_mask] / total_counts[nonzero_mask]) * 100
    
    adata_query.obs['pct_counts_mt'] = mt_pct
    
    print(f"  ✓ Calculated")
    print(f"    Range: [{mt_pct.min():.2f}, {mt_pct.max():.2f}]%")
    print(f"    Mean: {mt_pct.mean():.2f}%")
    print(f"    Median: {np.median(mt_pct):.2f}%")
    print(f"    Zero-count cells: {(~nonzero_mask).sum():,}")
else:
    adata_query.obs['pct_counts_mt'] = 0.0
    print(f"  ⚠ No MT genes found, set to 0")


# --- 2. Stress Score Calculation ---
print("\n" + "-" * 70)
print("2. Stress Score Calculation")
print("-" * 70)

stress_genes_in_data = [g for g in STRESS_SIGNATURE_GENES if g in adata_query.var_names]
n_stress = len(stress_genes_in_data)
print(f"  Signature genes: {len(STRESS_SIGNATURE_GENES)}")
print(f"  Found in data: {n_stress} ({n_stress/len(STRESS_SIGNATURE_GENES)*100:.1f}%)")

if n_stress >= 10:
    # Temporary normalized data
    adata_temp = sc.AnnData(
        X=adata_query.layers['counts'].copy(),
        var=adata_query.var.copy()
    )
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    
    stress_indices = [i for i, g in enumerate(adata_query.var_names) if g in stress_genes_in_data]
    stress_expr = adata_temp.X[:, stress_indices]
    
    if issparse(stress_expr):
        stress_score_raw = np.array(stress_expr.mean(axis=1)).flatten()
    else:
        stress_score_raw = stress_expr.mean(axis=1)
    
    # Normalize to [0, 1]
    stress_min = stress_score_raw.min()
    stress_max = stress_score_raw.max()
    
    if stress_max > stress_min:
        stress_score = (stress_score_raw - stress_min) / (stress_max - stress_min)
    else:
        stress_score = np.zeros_like(stress_score_raw)
    
    adata_query.obs['stress_score'] = stress_score
    
    del adata_temp, stress_expr
    gc.collect()
    
    print(f"  ✓ Calculated")
    print(f"    Range: [{stress_score.min():.3f}, {stress_score.max():.3f}]")
    print(f"    Mean: {stress_score.mean():.3f}")
    print(f"    Median: {np.median(stress_score):.3f}")
else:
    adata_query.obs['stress_score'] = 0.0
    print(f"  ⚠ Too few genes (<10), set to 0")


# --- 3. Cell Cycle Scores ---
print("\n" + "-" * 70)
print("3. Cell Cycle Scores")
print("-" * 70)

s_genes = [g for g in S_GENES if g in adata_query.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata_query.var_names]

print(f"  S genes: {len(s_genes)}/{len(S_GENES)} ({len(s_genes)/len(S_GENES)*100:.1f}%)")
print(f"  G2M genes: {len(g2m_genes)}/{len(G2M_GENES)} ({len(g2m_genes)/len(G2M_GENES)*100:.1f}%)")

if len(s_genes) >= 10 and len(g2m_genes) >= 10:
    adata_temp = sc.AnnData(
        X=adata_query.layers['counts'].copy(),
        var=adata_query.var.copy()
    )
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    
    sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_genes, g2m_genes=g2m_genes)
    
    adata_query.obs['S_score'] = adata_temp.obs['S_score'].values
    adata_query.obs['G2M_score'] = adata_temp.obs['G2M_score'].values
    adata_query.obs['phase'] = adata_temp.obs['phase'].values
    
    del adata_temp
    gc.collect()
    
    print(f"  ✓ Calculated")
    print(f"    S_score: mean={adata_query.obs['S_score'].mean():.3f}")
    print(f"    G2M_score: mean={adata_query.obs['G2M_score'].mean():.3f}")
    print(f"    Phase distribution:")
    phase_counts = adata_query.obs['phase'].value_counts()
    for phase, count in phase_counts.items():
        print(f"      {phase}: {count:,} ({count/adata_query.n_obs*100:.1f}%)")
else:
    adata_query.obs['S_score'] = 0.0
    adata_query.obs['G2M_score'] = 0.0
    adata_query.obs['phase'] = 'G1'
    print(f"  ⚠ Insufficient genes, set to 0")


# --- 4. Batch Key Validation ---
print("\n" + "-" * 70)
print("4. Batch Key Validation")
print("-" * 70)

if BATCH_KEY not in adata_query.obs.columns:
    print(f"  '{BATCH_KEY}' not found, creating default")
    adata_query.obs[BATCH_KEY] = "query_batch"
    print(f"  ✓ Created: 'query_batch'")
else:
    print(f"  ✓ Found: {BATCH_KEY}")
    n_batches = adata_query.obs[BATCH_KEY].nunique()
    print(f"    Unique batches: {n_batches}")
    
    if n_batches <= 10:
        batch_counts = adata_query.obs[BATCH_KEY].value_counts()
        print(f"    Distribution:")
        for batch, count in batch_counts.items():
            print(f"      {batch}: {count:,}")


# --- 5. Tissue Key Validation (CRITICAL) ---
print("\n" + "-" * 70)
print("5. Tissue Key Validation (CRITICAL)")
print("-" * 70)

if TISSUE_KEY not in adata_query.obs.columns:
    print(f"  '{TISSUE_KEY}' not found, creating default")
    adata_query.obs[TISSUE_KEY] = "query_tissue"
    print(f"  ✓ Created: 'query_tissue'")
else:
    print(f"  ✓ Found: {TISSUE_KEY}")
    n_tissues = adata_query.obs[TISSUE_KEY].nunique()
    print(f"    Unique tissues: {n_tissues}")
    
    if n_tissues <= 10:
        tissue_counts = adata_query.obs[TISSUE_KEY].value_counts()
        print(f"    Distribution:")
        for tissue, count in tissue_counts.items():
            print(f"      {tissue}: {count:,}")


# --- Validation Summary ---
print("\n" + "-" * 70)
print("Covariate Validation Summary")
print("-" * 70)

print("\nCategorical covariates:")
for cov in CATEGORICAL_COVARIATES:
    if cov in adata_query.obs.columns:
        n_cats = adata_query.obs[cov].nunique()
        print(f"  ✓ {cov}: {n_cats} categories")
    else:
        print(f"  ❌ {cov}: MISSING")

print("\nContinuous covariates:")
for cov in CONTINUOUS_COVARIATES:
    if cov in adata_query.obs.columns:
        mean_val = adata_query.obs[cov].mean()
        print(f"  ✓ {cov}: mean={mean_val:.3f}")
    else:
        print(f"  ❌ {cov}: MISSING")

# Final check
missing = [c for c in REQUIRED_COVARIATES if c not in adata_query.obs.columns]
if missing:
    print(f"\n❌ ERROR: Missing required covariates: {missing}")
    sys.exit(1)

print(f"\n✓ All required covariates present")


# ==============================================================================
# Step 2: Gene Alignment with Reference
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Gene Alignment with Reference")
print("=" * 70)

print(f"\nLoading reference HVG list: {HVG_FILE}")
hvg_genes_ref = pd.read_csv(HVG_FILE, header=None)[0].astype(str).tolist()
print(f"  Reference HVG: {len(hvg_genes_ref):,}")

query_genes_set = set(adata_query.var_names.astype(str))
hvg_in_query = [g for g in hvg_genes_ref if g in query_genes_set]
overlap_pct = len(hvg_in_query) / len(hvg_genes_ref) * 100

print(f"\nGene overlap analysis:")
print(f"  Query genes: {len(query_genes_set):,}")
print(f"  Overlap: {len(hvg_in_query):,}/{len(hvg_genes_ref):,}")
print(f"  Percentage: {overlap_pct:.1f}%")

if overlap_pct < 70:
    print(f"\n❌ ERROR: Insufficient gene overlap ({overlap_pct:.1f}%)")
    print(f"\n   Minimum required: 70%")
    print(f"   Possible causes:")
    print(f"     1. Different gene naming (ENSEMBL vs symbols)")
    print(f"     2. Different organism/genome version")
    print(f"     3. Excessive gene filtering")
    print(f"\n   Sample reference genes: {hvg_genes_ref[:5]}")
    print(f"   Sample query genes: {list(query_genes_set)[:5]}")
    sys.exit(1)

print(f"\n✓ Gene overlap acceptable")


# ==============================================================================
# Step 3: Prepare Query for scArches
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Prepare Query for scArches")
print("=" * 70)

print("\n⚠️ CRITICAL STEP: prepare_query_anndata")
print("  This function will:")
print("    1. Validate query has required covariates")
print("    2. Reorder genes to match reference HVG order")
print("    3. Pad missing genes with zeros")
print("    4. Validate data structure for scArches")
print("\n  ⚠️ Why .X must be counts:")
print("    - Reference was trained with layer=None (X=counts)")
print("    - Query must match this setup exactly")

print("\nSetting .X to raw counts...")
adata_query.X = adata_query.layers['counts']  # Shared memory, 0 cost
print("✓ .X = counts (0 memory cost)")

print("\nRunning scvi.model.SCANVI.prepare_query_anndata...")
print("  (This may take a few minutes for large datasets)")

try:
    scvi.model.SCANVI.prepare_query_anndata(
        adata_query,
        SCANVI_REF_DIR
    )
    print("✓ Query prepared successfully")
except Exception as e:
    print(f"\n❌ ERROR during prepare_query_anndata:")
    print(f"   {e}")
    print(f"\n   Debugging info:")
    print(f"   - Query shape: {adata_query.shape}")
    print(f"   - Reference dir: {SCANVI_REF_DIR}")
    sys.exit(1)

print(f"\nQuery after preparation:")
print(f"  Shape: {adata_query.n_obs:,} × {adata_query.n_vars:,}")
print(f"  Genes match reference: {adata_query.n_vars == len(hvg_genes_ref)}")

# ⚠️ CRITICAL VALIDATION: Ensure covariates survived preparation
print(f"\n⚠️ CRITICAL: Validating covariates after preparation...")
missing_after_prep = [c for c in REQUIRED_COVARIATES if c not in adata_query.obs.columns]
if missing_after_prep:
    print(f"  ❌ ERROR: Covariates lost during preparation!")
    print(f"     Missing: {missing_after_prep}")
    print(f"\n     This should NOT happen if covariates were calculated before prepare.")
    print(f"     Please report this as a bug.")
    sys.exit(1)
else:
    print(f"  ✓ All {len(REQUIRED_COVARIATES)} covariates preserved")
    
    # Verify values are still correct
    for cov in CONTINUOUS_COVARIATES:
        if cov in adata_query.obs.columns:
            mean_val = adata_query.obs[cov].mean()
            if np.isnan(mean_val) or mean_val == 0:
                print(f"     ⚠️ Warning: {cov} may have been reset (mean={mean_val:.3f})")
            else:
                print(f"     ✓ {cov}: mean={mean_val:.3f}")


# ==============================================================================
# Step 3.5: Metadata Cleaning (PRESERVE COVARIATES)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3.5: Metadata Cleaning (CRITICAL)")
print("=" * 70)

print("\nStrategy:")
print("  1. Backup original obs columns")
print("  2. Build clean obs with:")
print("     - Strict categorical types (batch, tissue, labels)")
print("     - ALL continuous covariates (preserved as-is)")
print("  3. Set all labels to Unknown (scArches standard)")

# Backup
adata_query.uns['original_obs_columns'] = list(adata_query.obs.columns)
n_original_cols = len(adata_query.obs.columns)
print(f"\n✓ Backed up {n_original_cols} original columns")

# Build clean obs
obs_clean = pd.DataFrame(index=adata_query.obs_names)

# --- Batch key (categorical) ---
print(f"\n1. Processing {BATCH_KEY}:")
if BATCH_KEY in adata_query.obs.columns:
    b = adata_query.obs[BATCH_KEY]
else:
    b = pd.Series(["query_batch"] * adata_query.n_obs, index=adata_query.obs_names)

# Force to string, fill NA, convert to categorical
b_str = b.astype("string").fillna("unknown_batch").astype(str)
obs_clean[BATCH_KEY] = pd.Categorical(b_str)

print(f"   ✓ Type: category")
print(f"   ✓ Categories: {obs_clean[BATCH_KEY].nunique()}")
print(f"   ✓ NA count: 0")

# --- Tissue key (categorical) ---
print(f"\n2. Processing {TISSUE_KEY}:")
if TISSUE_KEY in adata_query.obs.columns:
    t = adata_query.obs[TISSUE_KEY]
else:
    t = pd.Series(["query_tissue"] * adata_query.n_obs, index=adata_query.obs_names)

t_str = t.astype("string").fillna("unknown_tissue").astype(str)
obs_clean[TISSUE_KEY] = pd.Categorical(t_str)

print(f"   ✓ Type: category")
print(f"   ✓ Categories: {obs_clean[TISSUE_KEY].nunique()}")
print(f"   ✓ NA count: 0")

# --- Labels key (categorical, ALL Unknown) ---
print(f"\n3. Processing {LABELS_KEY}:")

# Backup original labels if they exist
if LABELS_KEY in adata_query.obs.columns:
    backup_col = f"{LABELS_KEY}_query_original"
    obs_clean[backup_col] = adata_query.obs[LABELS_KEY].astype("string").astype(str)
    print(f"   ✓ Original labels backed up to: {backup_col}")

# Set ALL to Unknown (scArches standard for query mapping)
obs_clean[LABELS_KEY] = pd.Categorical([UNLABELED_CATEGORY] * adata_query.n_obs)
print(f"   ✓ All {adata_query.n_obs:,} cells set to: {UNLABELED_CATEGORY}")

# --- Continuous covariates (preserve as-is) ---
print(f"\n4. Processing continuous covariates:")
for cov in CONTINUOUS_COVARIATES:
    if cov in adata_query.obs.columns:
        obs_clean[cov] = adata_query.obs[cov].values
        print(f"   ✓ {cov}: preserved (mean={obs_clean[cov].mean():.3f})")
    else:
        print(f"   ❌ {cov}: MISSING")

# --- Replace obs ---
adata_query.obs = obs_clean

# --- Final validation ---
print(f"\n5. Final validation:")

validation_errors = []

for cov in CATEGORICAL_COVARIATES:
    if cov not in adata_query.obs.columns:
        validation_errors.append(f"Missing categorical: {cov}")
    else:
        s = adata_query.obs[cov]
        if not hasattr(s, 'cat'):
            validation_errors.append(f"Not categorical: {cov}")
        else:
            n_na = pd.isna(s).sum()
            if n_na > 0:
                validation_errors.append(f"Has NA values: {cov}")
            print(f"   ✓ {cov}: category ({len(s.cat.categories)} cats, 0 NA)")

for cov in CONTINUOUS_COVARIATES:
    if cov not in adata_query.obs.columns:
        validation_errors.append(f"Missing continuous: {cov}")
    else:
        s = adata_query.obs[cov]
        print(f"   ✓ {cov}: {s.dtype} (mean={s.mean():.3f})")

if validation_errors:
    print(f"\n❌ Validation errors:")
    for err in validation_errors:
        print(f"   - {err}")
    sys.exit(1)

print(f"\n✓ Metadata cleaning complete")
print(f"  - All categorical covariates: strict typing, 0 NA")
print(f"  - All continuous covariates: preserved")
print(f"  - Labels: all set to Unknown")


# ==============================================================================
# Step 4: Load Query into Reference Model
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Load Query into Reference Model")
print("=" * 70)

print(f"\nLoading reference model: {SCANVI_REF_DIR}")
print("  Performing scArches surgery:")
print("    - Clone reference architecture")
print("    - Attach query data with covariates")
print("    - Initialize for fine-tuning")

try:
    scanvi_query = scvi.model.SCANVI.load_query_data(
        adata_query,
        SCANVI_REF_DIR
    )
    print("✓ Query loaded successfully")
except Exception as e:
    print(f"\n❌ ERROR during load_query_data:")
    print(f"   {e}")
    print(f"\n   Debugging info:")
    print(f"   - Query shape: {adata_query.shape}")
    print(f"   - Obs columns: {list(adata_query.obs.columns)}")
    print(f"   - Batch key: {BATCH_KEY} (exists: {BATCH_KEY in adata_query.obs})")
    print(f"   - Tissue key: {TISSUE_KEY} (exists: {TISSUE_KEY in adata_query.obs})")
    print(f"   - Labels key: {LABELS_KEY} (exists: {LABELS_KEY in adata_query.obs})")
    
    # Check data types
    print(f"\n   Data types:")
    for col in [BATCH_KEY, TISSUE_KEY, LABELS_KEY]:
        if col in adata_query.obs.columns:
            dtype = adata_query.obs[col].dtype
            is_cat = hasattr(adata_query.obs[col], 'cat')
            print(f"   - {col}: {dtype} (categorical: {is_cat})")
    
    sys.exit(1)

print(f"\n✓ scArches surgery complete")
print(f"  Query integrated into reference architecture")


# ==============================================================================
# Step 5: Fine-tune on Query Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: Fine-tune Model on Query Data")
print("=" * 70)

print(f"\nTraining configuration:")
print(f"  Max epochs: {MAX_EPOCHS}")
print(f"  Batch size: {BATCH_SIZE}")
print(f"  Learning rate: {LEARNING_RATE}")
print(f"  Weight decay: {WEIGHT_DECAY} ⚠️ MUST be 0 for scArches")
print(f"  Early stopping patience: {EARLY_STOPPING_PATIENCE}")

train_kwargs = {
    'max_epochs': MAX_EPOCHS,
    'batch_size': BATCH_SIZE,
    'early_stopping': True,
    'early_stopping_patience': EARLY_STOPPING_PATIENCE,
    'check_val_every_n_epoch': CHECK_VAL_EVERY_N_EPOCH,
    'plan_kwargs': {
        'lr': LEARNING_RATE,
        'weight_decay': WEIGHT_DECAY
    },
    'enable_progress_bar': True
}

if gpu_available:
    train_kwargs['accelerator'] = 'gpu'
    train_kwargs['devices'] = 'auto'
    print(f"\n  Using GPU: {torch.cuda.get_device_name(0)}")
else:
    print(f"\n  Using CPU (warning: will be slow)")

print(f"\nStarting fine-tuning...")
train_start = datetime.now()

try:
    scanvi_query.train(**train_kwargs)
    train_end = datetime.now()
    train_duration = (train_end - train_start).total_seconds() / 60
    print(f"\n✓ Fine-tuning complete")
    print(f"  Duration: {train_duration:.1f} minutes")
except Exception as e:
    print(f"\n❌ ERROR during training:")
    print(f"   {e}")
    sys.exit(1)


# ==============================================================================
# Step 6: Extract Results
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Extract Results")
print("=" * 70)

# --- Latent representation ---
print("\n1. Extracting latent representation...")
LATENT_KEY = "X_scANVI_mapped"
try:
    adata_query.obsm[LATENT_KEY] = scanvi_query.get_latent_representation()
    print(f"   ✓ .obsm['{LATENT_KEY}']: {adata_query.obsm[LATENT_KEY].shape}")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# --- Hard predictions ---
print("\n2. Extracting cell type predictions...")
PRED_KEY = "cell_type_mapped"
try:
    adata_query.obs[PRED_KEY] = scanvi_query.predict()
    n_types = adata_query.obs[PRED_KEY].nunique()
    print(f"   ✓ .obs['{PRED_KEY}']: {n_types} cell types")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# --- Soft probabilities and confidence ---
print("\n3. Calculating prediction confidence...")
try:
    probs = scanvi_query.predict(soft=True)
    confidence = probs.max(axis=1)
    
    CONF_KEY = "mapping_confidence"
    adata_query.obs[CONF_KEY] = confidence
    print(f"   ✓ .obs['{CONF_KEY}']: range [{confidence.min():.3f}, {confidence.max():.3f}]")
    
    # Margin (top1 - top2)
    top2_probs = np.partition(probs, -2, axis=1)[:, -2:]
    margin = top2_probs[:, 1] - top2_probs[:, 0]
    adata_query.obs['mapping_margin'] = margin
    print(f"   ✓ .obs['mapping_margin']: mean={margin.mean():.3f}")
except Exception as e:
    print(f"   ❌ ERROR: {e}")
    sys.exit(1)

# --- Apply confidence threshold ---
print(f"\n4. Applying confidence threshold: {CONFIDENCE_THRESHOLD}")
FINAL_KEY = "cell_type_final"
adata_query.obs[FINAL_KEY] = np.where(
    confidence < CONFIDENCE_THRESHOLD,
    UNLABELED_CATEGORY,
    adata_query.obs[PRED_KEY].astype(str)
)

n_high_conf = (confidence >= CONFIDENCE_THRESHOLD).sum()
n_low_conf = (confidence < CONFIDENCE_THRESHOLD).sum()
print(f"   High confidence: {n_high_conf:,} ({n_high_conf/adata_query.n_obs*100:.1f}%)")
print(f"   Low confidence: {n_low_conf:,} ({n_low_conf/adata_query.n_obs*100:.1f}%)")

print(f"\n✓ All results extracted")


# ==============================================================================
# Step 7: UMAP Projection
# ==============================================================================

print("\n" + "=" * 70)
print("Step 7: UMAP Projection")
print("=" * 70)

UMAP_PROJECT_KEY = "X_umap_mapped"

if hasattr(scanvi_query, 'umap_op_'):
    print("✓ Reference has UMAP operator")
    print("  Projecting query into reference UMAP space...")
    
    try:
        adata_query.obsm[UMAP_PROJECT_KEY] = scanvi_query.umap_op_.transform(
            adata_query.obsm[LATENT_KEY]
        )
        print(f"  ✓ .obsm['{UMAP_PROJECT_KEY}']: {adata_query.obsm[UMAP_PROJECT_KEY].shape}")
        print(f"  ✓ Query and reference in SAME UMAP space")
    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        sys.exit(1)
else:
    print("⚠️ No UMAP operator in reference model")
    print("  Computing new UMAP on query latent...")
    
    try:
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
        print(f"  ✓ .obsm['{UMAP_PROJECT_KEY}']: {adata_query.obsm[UMAP_PROJECT_KEY].shape}")
        print(f"  ⚠️ Coordinates differ from reference")
    except Exception as e:
        print(f"  ❌ ERROR: {e}")
        sys.exit(1)

print(f"\n✓ UMAP projection complete")


# ==============================================================================
# Step 8: Quality Control Summary
# ==============================================================================

print("\n" + "=" * 70)
print("Step 8: Quality Control Summary")
print("=" * 70)

# --- Confidence distribution ---
print("\n1. Prediction Confidence:")
print(adata_query.obs[CONF_KEY].describe())

low_conf_count = (adata_query.obs[CONF_KEY] < CONFIDENCE_THRESHOLD).sum()
low_conf_pct = low_conf_count / adata_query.n_obs * 100
print(f"\n   Low confidence (<{CONFIDENCE_THRESHOLD}): {low_conf_count:,} ({low_conf_pct:.1f}%)")

# --- Margin distribution ---
print("\n2. Prediction Margin (Top1 - Top2):")
print(adata_query.obs['mapping_margin'].describe())

low_margin_count = (adata_query.obs['mapping_margin'] < 0.1).sum()
low_margin_pct = low_margin_count / adata_query.n_obs * 100
print(f"\n   Low margin (<0.1): {low_margin_count:,} ({low_margin_pct:.1f}%)")
print(f"   (Ambiguous predictions)")

# --- Cell type distribution ---
print("\n3. Final Cell Type Distribution:")
celltype_counts = adata_query.obs[FINAL_KEY].value_counts()
print(celltype_counts.head(20))

# --- Unknown cells ---
unknown_count = (adata_query.obs[FINAL_KEY] == UNLABELED_CATEGORY).sum()
unknown_pct = unknown_count / adata_query.n_obs * 100
print(f"\n4. {UNLABELED_CATEGORY} Cells:")
print(f"   Count: {unknown_count:,} ({unknown_pct:.1f}%)")

if unknown_pct > 30:
    print(f"\n   ⚠️ High Unknown percentage ({unknown_pct:.1f}%)")
    print(f"   Possible causes:")
    print(f"     - Reference doesn't cover query cell types")
    print(f"     - Query data quality issues")
    print(f"     - Threshold too strict ({CONFIDENCE_THRESHOLD})")

# --- Batch effects ---
if adata_query.obs[BATCH_KEY].nunique() > 1:
    print("\n5. Cell Type by Batch:")
    batch_celltype = pd.crosstab(
        adata_query.obs[BATCH_KEY],
        adata_query.obs[FINAL_KEY],
        normalize='index'
    ) * 100
    print(batch_celltype.round(1).head())

print(f"\n✓ Quality control complete")


# ==============================================================================
# Step 9: Restore Metadata and Save
# ==============================================================================

print("\n" + "=" * 70)
print("Step 9: Restore Metadata and Save")
print("=" * 70)

if 'original_obs_columns' in adata_query.uns:
    print("\nRestoring original metadata...")
    
    # Reload original h5ad
    print(f"  Loading original: {QUERY_H5AD}")
    adata_original = sc.read_h5ad(QUERY_H5AD)
    
    # Align cell order
    if not (adata_original.obs_names == adata_query.obs_names).all():
        print(f"  Aligning cell order...")
        adata_original = adata_original[adata_query.obs_names].copy()
    
    # Merge: add new columns to original
    new_cols = [c for c in adata_query.obs.columns if c not in adata_original.obs.columns]
    
    print(f"  Merging metadata:")
    print(f"    Original: {adata_original.obs.shape[1]} columns")
    print(f"    New: {len(new_cols)} columns")
    
    for col in new_cols:
        adata_original.obs[col] = adata_query.obs[col].values
    
    # Copy obsm (latent, UMAP)
    for key in adata_query.obsm.keys():
        adata_original.obsm[key] = adata_query.obsm[key]
    
    adata_query = adata_original
    
    print(f"  ✓ Restored {adata_original.obs.shape[1]} total columns")
    
    del adata_original
    gc.collect()

# Save
output_file = output_dir / "query_mapped_to_reference.h5ad"
print(f"\nSaving: {output_file}")

try:
    adata_query.write_h5ad(output_file, compression='gzip', compression_opts=9)
    file_size = output_file.stat().st_size / (1024**3)
    print(f"✓ Saved ({file_size:.2f} GB)")
except Exception as e:
    print(f"❌ ERROR saving: {e}")
    sys.exit(1)

# Cleanup
del scanvi_query
gc.collect()


# ==============================================================================
# Step 10: Visualization
# ==============================================================================

print("\n" + "=" * 70)
print("Step 10: Visualization")
print("=" * 70)

adata_query.obsm['X_umap'] = adata_query.obsm[UMAP_PROJECT_KEY].copy()

print("\nGenerating overview figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

try:
    # Row 1
    sc.pl.umap(adata_query, color=BATCH_KEY, ax=axes[0, 0], show=False, 
               title='Batch Distribution', legend_loc='right margin')
    sc.pl.umap(adata_query, color=PRED_KEY, ax=axes[0, 1], show=False,
               title='Predicted Cell Types', legend_loc='right margin')
    sc.pl.umap(adata_query, color=FINAL_KEY, ax=axes[0, 2], show=False,
               title=f'Final (conf>{CONFIDENCE_THRESHOLD})', legend_loc='right margin')
    
    # Row 2
    sc.pl.umap(adata_query, color=CONF_KEY, ax=axes[1, 0], show=False,
               title='Confidence', cmap='viridis', vmin=0, vmax=1)
    sc.pl.umap(adata_query, color='mapping_margin', ax=axes[1, 1], show=False,
               title='Margin (Top1-Top2)', cmap='viridis')
    
    # Confidence histogram
    axes[1, 2].hist(adata_query.obs[CONF_KEY], bins=50, edgecolor='black', alpha=0.7)
    axes[1, 2].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', linewidth=2,
                       label=f'Threshold={CONFIDENCE_THRESHOLD}')
    axes[1, 2].set_xlabel('Prediction Confidence', fontsize=12)
    axes[1, 2].set_ylabel('Number of Cells', fontsize=12)
    axes[1, 2].set_title('Confidence Distribution', fontsize=14)
    axes[1, 2].legend(fontsize=10)
    axes[1, 2].grid(alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    fig_pdf = output_dir / "figures" / f"mapping_overview.{FIGURE_FORMAT}"
    fig_png = output_dir / "figures" / "mapping_overview.png"
    
    plt.savefig(fig_pdf, dpi=DPI, bbox_inches='tight')
    plt.savefig(fig_png, dpi=150, bbox_inches='tight')
    
    print(f"✓ Saved: mapping_overview.{FIGURE_FORMAT}")
    print(f"✓ Saved: mapping_overview.png")
    
    plt.close()
except Exception as e:
    print(f"⚠️ Visualization warning: {e}")


# ==============================================================================
# Step 11: Export Summary
# ==============================================================================

print("\n" + "=" * 70)
print("Step 11: Export Summary Report")
print("=" * 70)

# Create summary dict
summary = {
    'pipeline': 'scArches_scANVI_mapping_FIXED',
    'version': '1.2.1_PRODUCTION',
    'timestamp': datetime.now().isoformat(),
    'query': {
        'cells': int(adata_query.n_obs),
        'genes': int(adata_query.n_vars),
        'batches': int(adata_query.obs[BATCH_KEY].nunique()),
        'tissues': int(adata_query.obs[TISSUE_KEY].nunique())
    },
    'gene_overlap': {
        'reference_hvg': len(hvg_genes_ref),
        'overlap_count': len(hvg_in_query),
        'overlap_percent': float(overlap_pct)
    },
    'covariates': {
        'pct_counts_mt': {
            'mean': float(adata_query.obs['pct_counts_mt'].mean()),
            'median': float(np.median(adata_query.obs['pct_counts_mt']))
        },
        'stress_score': {
            'mean': float(adata_query.obs['stress_score'].mean()),
            'median': float(np.median(adata_query.obs['stress_score']))
        },
        'S_score': {
            'mean': float(adata_query.obs['S_score'].mean())
        },
        'G2M_score': {
            'mean': float(adata_query.obs['G2M_score'].mean())
        }
    },
    'predictions': {
        'unique_cell_types': int(adata_query.obs[FINAL_KEY].nunique()),
        'confidence_median': float(adata_query.obs[CONF_KEY].median()),
        'margin_median': float(adata_query.obs['mapping_margin'].median()),
        'low_confidence_count': int(low_conf_count),
        'low_confidence_percent': float(low_conf_pct),
        'unknown_count': int(unknown_count),
        'unknown_percent': float(unknown_pct)
    },
    'top_cell_types': dict(celltype_counts.head(10).to_dict())
}

# Save JSON
summary_json = output_dir / "mapping_summary.json"
with open(summary_json, 'w') as f:
    json.dump(summary, f, indent=2)
print(f"✓ Saved: {summary_json}")

# Generate text report
report_lines = []
report_lines.append("=" * 70)
report_lines.append("scArches scANVI Mapping Summary (FIXED v1.2.1 PRODUCTION)")
report_lines.append("=" * 70)
report_lines.append(f"\nTimestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
report_lines.append(f"\n[ Configuration ]")
report_lines.append(f"  Reference model: {SCANVI_REF_DIR}")
report_lines.append(f"  Query data: {QUERY_H5AD}")
report_lines.append(f"  Output: {output_file}")

report_lines.append(f"\n[ Query Dataset ]")
report_lines.append(f"  Cells: {adata_query.n_obs:,}")
report_lines.append(f"  Genes: {adata_query.n_vars:,}")
report_lines.append(f"  Batches: {adata_query.obs[BATCH_KEY].nunique()}")
report_lines.append(f"  Tissues: {adata_query.obs[TISSUE_KEY].nunique()}")

report_lines.append(f"\n[ Gene Overlap ]")
report_lines.append(f"  Reference HVG: {len(hvg_genes_ref):,}")
report_lines.append(f"  Overlap: {len(hvg_in_query):,} ({overlap_pct:.1f}%)")

report_lines.append(f"\n[ Covariates Calculated ]")
report_lines.append(f"  pct_counts_mt: mean={adata_query.obs['pct_counts_mt'].mean():.2f}%")
report_lines.append(f"  stress_score: mean={adata_query.obs['stress_score'].mean():.3f}")
report_lines.append(f"  S_score: mean={adata_query.obs['S_score'].mean():.3f}")
report_lines.append(f"  G2M_score: mean={adata_query.obs['G2M_score'].mean():.3f}")

report_lines.append(f"\n[ Mapping Quality ]")
report_lines.append(f"  Confidence median: {adata_query.obs[CONF_KEY].median():.3f}")
report_lines.append(f"  Margin median: {adata_query.obs['mapping_margin'].median():.3f}")
report_lines.append(f"  Low confidence: {low_conf_count:,} ({low_conf_pct:.1f}%)")
report_lines.append(f"  Unknown cells: {unknown_count:,} ({unknown_pct:.1f}%)")

report_lines.append(f"\n[ Top 10 Cell Types ]")
for ct, count in celltype_counts.head(10).items():
    pct = count / adata_query.n_obs * 100
    report_lines.append(f"  {ct}: {count:,} ({pct:.1f}%)")

report_lines.append(f"\n" + "=" * 70)
report_lines.append("Mapping completed successfully")
report_lines.append("=" * 70)

report_text = '\n'.join(report_lines)

# Print to console
print("\n" + report_text)

# Save to file
summary_txt = output_dir / "mapping_summary.txt"
with open(summary_txt, 'w') as f:
    f.write(report_text)
print(f"\n✓ Saved: {summary_txt}")


# ==============================================================================
# PIPELINE COMPLETE
# ==============================================================================

print("\n" + "=" * 70)
print("✅ PIPELINE COMPLETE")
print("=" * 70)

print(f"\n📊 Key Outputs:")
print(f"  1. Mapped data: {output_file}")
print(f"  2. Figures: {output_dir / 'figures'}")
print(f"  3. Summary JSON: {summary_json}")
print(f"  4. Summary text: {summary_txt}")

print(f"\n🔬 Critical Features:")
print(f"  ✓ All covariates calculated (MT%, stress, cell cycle)")
print(f"  ✓ Tissue key validated/created")
print(f"  ✓ Metadata cleaning preserved covariates")
print(f"  ✓ Strict categorical typing (no mixed types)")
print(f"  ✓ Gene overlap validated ({overlap_pct:.1f}%)")

print(f"\n📈 Quality Metrics:")
print(f"  - Confidence median: {adata_query.obs[CONF_KEY].median():.3f}")
print(f"  - Low confidence: {low_conf_pct:.1f}%")
print(f"  - Unknown cells: {unknown_pct:.1f}%")
print(f"  - Predicted cell types: {adata_query.obs[FINAL_KEY].nunique()}")

print(f"\n📖 Next Steps:")
print(f"  1. Review confidence distribution")
print(f"  2. Validate with marker genes")
print(f"  3. Compare UMAP with reference (if available)")
print(f"  4. Adjust threshold if needed (current: {CONFIDENCE_THRESHOLD})")

print(f"\n💡 Tips:")
print(f"  - Low confidence cells: consider manual curation")
print(f"  - High Unknown %: may need different reference model")
print(f"  - Check batch effects in visualizations")

print("\n" + "=" * 70)
print("Thank you for using scArches mapping pipeline!")
print("=" * 70)
