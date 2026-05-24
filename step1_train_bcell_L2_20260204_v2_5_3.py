#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
==============================================================================
Step 1: Train B Cell L2 Reference Model - v2.5.3 PRODUCTION (ALL FIXES)
==============================================================================

CRITICAL FIXES Applied:
-----------------------
[OK] P0-1: Save UMAP operator (joblib) for same-space projection
[OK] P0-2: Output reference_with_umap_L2.h5ad for merging
[OK] P0-3: Remove .raw construction (avoid memory risk)
[OK] P0-4: Ensure Unknown category always exists
[OK] P1-1: Pre-add unknown_tissue to categories + save to config
[OK] P1-3: Three-level HVG fallback (seurat_v3 → cell_ranger)
[OK] P1-4: var_names_make_unique for consistency
[OK] P2-1: Unified naming (L2_KEY, SCANVI_MODEL_DIR)

Pipeline Flow:
--------------
1. Load reference data (with L3 labels)
2. Map L3 → L2 (aggregate to coarse labels)
3. Calculate covariates (MT%, stress, cell cycle) - OPTIMIZED
4. Prepare categorical covariates for scArches
5. Select HVG with 3-level fallback
6. Force biological markers into HVG
7. Train scVI (unsupervised integration)
8. Train scANVI (semi-supervised on L2 labels)
9. [OK] Fit and save UMAP operator for query projection
10. [OK] Save reference with UMAP for merging
11. Save model + HVG list + config

Version: 2.5.3 PRODUCTION
Date: 2025-02-04
Author: r2end
"""

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

import scanpy as sc
import scvi

# For UMAP operator saving (P0-1 fix)
import umap
import joblib

warnings.filterwarnings('ignore')

print("=" * 70)
print("B Cell L2 Reference Training - v2.5.3 PRODUCTION (ALL FIXES)")
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
# CONFIGURATION
# ==============================================================================

print("\n" + "=" * 70)
print("CONFIGURATION")
print("=" * 70)

# ===== Input/Output Paths =====
INPUT_REF_H5AD = "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3"

# ===== Keys (UNIFIED NAMING - P2-1 fix) =====
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
L3_KEY = "cell_type_expert"  # Source (fine-grained)
L2_KEY = "Cell_Type_L2"  # Target (coarse-grained)
UNLABELED_CATEGORY = "Unknown"

# ===== L3 → L2 Mapping =====
L3_TO_L2_MAP = {
    # Atypical Memory B
    "Atypical_Memory_B": "Atypical_Memory_B",
    "IGHEplus_Atypical_Memory_B": "Atypical_Memory_B",
    
    # GC B
    "GC_B_Light_Zone_Centrocyte": "GC_B",
    "GC_B_Transitional": "GC_B",
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B": "GC_B",
    
    # Memory B
    "Memory_B": "Memory_B",
    
    # Naive B
    "Naive_B": "Naive_B",
    
    # Plasma
    "Plasma_IgA": "Plasma",
    "Plasma_IgG": "Plasma",
    "Plasma": "Plasma",
    
    # Unknown
    "Unknown": "Unknown"
}

# ===== Covariate Gene Lists =====
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
    "XPOT","YIF1A","YWHAZ","ZBTB17","ARFGEF1","ATG9A","CAND1","CCAR1","CD58",
    "CUL4A","DBI","DDX17","DDX18","DDX21","DYNC1I2","EIF2B4","EIF3A","FAM120A",
    "FKBP1A","GANAB","GNB1","HNRNPA0","HNRNPD","HSPD1","ILF3","ISG15","ITGB1",
    "KRR1","MRPL11","MRPS18C","MVP","NCL","NMD3","NOP2","OGFOD1","PABPC1",
    "PHB","PHB2","POLR1B","PPP1CA","RANBP1","RBM3","RPL18A","RPL3","RPL7",
    "RPS3","RPS6","RTN4","SART1","SET","SSRP1","TCEB2","TMED10","TUBB","VDAC1",
    "VIM","VMP1"
]

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

# ===== HVG Settings =====
N_HVG = 4000
HVG_FLAVOR = "seurat_v3"

# ===== Force B Cell Markers into HVG =====
FORCE_MARKERS_IN_HVG = True
FORCED_MARKERS_BCELL = [
    # Pan B cell
    "CD19","MS4A1","CD79A","CD79B","BANK1","BLK","FCRL2","FCRL3","FCRL5",
    # Naive B
    "TCL1A","FCER2","IGHD","IL4R",
    # Memory B
    "CD27","TNFRSF13B","AIM2",
    # GC B
    "AICDA","BCL6","MME","LMO2","STMN1","CXCR4","RGS13",
    # Plasma
    "SDC1","CD38","XBP1","PRDM1","JCHAIN","IGHA1","IGHA2","IGHG1","IGHG2",
    "IGHG3","IGHG4","MZB1","DERL3","SSR4",
    # Atypical Memory B
    "ITGAX","TBX21","FCRL5","CXCR3"
]

# ===== scVI Parameters =====
SCVI_N_LATENT = 50
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_DROPOUT_RATE = 0.1
SCVI_GENE_LIKELIHOOD = "nb"

# ===== Training Parameters =====
MAX_EPOCHS_SCVI = 400
MAX_EPOCHS_SCANVI = 200
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 0.0  # scArches compatible

# ===== UMAP Parameters (P0-1 fix) =====
UMAP_N_NEIGHBORS = 30
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0

# ===== Reproducibility =====
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED

print("\nConfiguration:")
print(f"  Input: {INPUT_REF_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  L3 → L2 mapping: {len(L3_TO_L2_MAP)} entries")
print(f"  HVG: {N_HVG}")
print(f"  Forced markers: {len(FORCED_MARKERS_BCELL)}")


# ==============================================================================
# Create Output Directory
# ==============================================================================

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"\nOK Output directory: {output_dir}")


# ==============================================================================
# Step 1: Load and Validate Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Load and Validate Reference Data")
print("=" * 70)

print(f"\nLoading: {INPUT_REF_H5AD}")
adata = sc.read_h5ad(INPUT_REF_H5AD)

# P1-4 FIX: Make gene names unique
adata.var_names_make_unique()
print(f"OK Gene names made unique")

print(f"Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Check for L3 labels
if L3_KEY not in adata.obs.columns:
    raise ValueError(f"Column '{L3_KEY}' not found in reference data!")

print(f"\nL3 label distribution:")
l3_counts = adata.obs[L3_KEY].value_counts()
for label, count in l3_counts.head(10).items():
    print(f"  {label}: {count:,}")


# ==============================================================================
# Step 2: L3 → L2 Mapping
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Map L3 → L2 Labels")
print("=" * 70)

print(f"\nMapping {L3_KEY} → {L2_KEY}...")

l3_series = adata.obs[L3_KEY].astype(str)
adata.obs[L2_KEY] = l3_series.map(L3_TO_L2_MAP)

# Check unmapped
n_missing = adata.obs[L2_KEY].isna().sum()
if n_missing > 0:
    print(f"\nWARNING  Warning: {n_missing} cells have unmapped L3 labels")
    unmapped_labels = l3_series[adata.obs[L2_KEY].isna()].unique()
    print(f"Unmapped L3 labels: {list(unmapped_labels)}")
    print("Setting to 'Unknown'...")
    adata.obs[L2_KEY] = adata.obs[L2_KEY].fillna(UNLABELED_CATEGORY)

# Convert to categorical
adata.obs[L2_KEY] = adata.obs[L2_KEY].astype("category")

# P0-4 FIX: Ensure Unknown category exists
if UNLABELED_CATEGORY not in adata.obs[L2_KEY].cat.categories:
    adata.obs[L2_KEY] = adata.obs[L2_KEY].cat.add_categories([UNLABELED_CATEGORY])
    print(f"\nOK Added '{UNLABELED_CATEGORY}' category (required by scANVI)")

print(f"\nL2 label distribution:")
l2_counts = adata.obs[L2_KEY].value_counts()
for label, count in l2_counts.items():
    print(f"  {label}: {count:,} ({count/adata.n_obs*100:.1f}%)")


# ==============================================================================
# Step 3: Ensure .X is Raw Counts
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Data Format Validation")
print("=" * 70)

if 'counts' in adata.layers:
    adata.X = adata.layers['counts']
    print("Set .X = layers['counts'] (shared memory)")
else:
    adata.layers['counts'] = adata.X
    print("Created layers['counts'] from .X")


# ==============================================================================
# Step 4: Calculate Covariates
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Calculate Covariates")
print("=" * 70)

# --- 4.1: MT% ---
print("\n4.1: MT% Calculation")
print("-" * 40)

mt_col_found = False
for col_name in ['pct_counts_mt', 'percent.mt', 'percent_mt']:
    if col_name in adata.obs.columns:
        if col_name != 'pct_counts_mt':
            adata.obs['pct_counts_mt'] = adata.obs[col_name]
            print(f"OK Renamed '{col_name}' → 'pct_counts_mt'")
        else:
            print(f"OK 'pct_counts_mt' exists")
        mt_col_found = True
        break

if not mt_col_found:
    print("WARNING  No MT% column, calculating...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, layer='counts')
    print(f"OK Calculated")

print(f"  Mean: {adata.obs['pct_counts_mt'].mean():.2f}%")

# --- 4.2: Stress Score (OPTIMIZED - P1-2) ---
print("\n4.2: Stress Score (Memory-Optimized)")
print("-" * 40)

counts = adata.layers['counts']
stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
n_stress = len(stress_genes)
print(f"Stress genes: {n_stress}/{len(STRESS_SIGNATURE_GENES)}")

if n_stress >= 10:
    idx = adata.var_names.get_indexer(stress_genes)
    idx = idx[idx >= 0]
    
    sub = counts[:, idx]
    
    if issparse(counts):
        tot = np.asarray(counts.sum(axis=1)).ravel()
    else:
        tot = counts.sum(axis=1)
    tot[tot == 0] = 1.0
    
    from scipy import sparse
    scale = sparse.diags(1e4 / tot)
    sub_norm = scale @ sub
    
    if issparse(sub_norm):
        sub_norm = sub_norm.tocsr(copy=True)
        sub_norm.data = np.log1p(sub_norm.data)
        expr = np.asarray(sub_norm.mean(axis=1)).ravel()
    else:
        expr = np.log1p(sub_norm).mean(axis=1)
    
    mn, mx = expr.min(), expr.max()
    adata.obs['stress_score'] = (expr - mn) / (mx - mn) if mx > mn else 0.0
    
    print(f"OK Mean: {adata.obs['stress_score'].mean():.3f}")
    
    del sub, sub_norm, scale
    gc.collect()
else:
    adata.obs['stress_score'] = 0.0
    print("WARNING  Too few genes (<10), set to 0")

# --- 4.3: Cell Cycle ---
print("\n4.3: Cell Cycle Scores")
print("-" * 40)

s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]

print(f"S genes: {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes: {len(g2m_genes)}/{len(G2M_GENES)}")

if len(s_genes) >= 10 and len(g2m_genes) >= 10:
    adata_temp = sc.AnnData(
        X=adata.layers['counts'].copy(),
        var=adata.var.copy()
    )
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_genes, g2m_genes=g2m_genes)
    
    adata.obs['S_score'] = adata_temp.obs['S_score'].values
    adata.obs['G2M_score'] = adata_temp.obs['G2M_score'].values
    adata.obs['phase'] = adata_temp.obs['phase'].values
    
    del adata_temp
    gc.collect()
    
    print(f"OK S mean: {adata.obs['S_score'].mean():.3f}")
    print(f"OK G2M mean: {adata.obs['G2M_score'].mean():.3f}")
else:
    adata.obs['S_score'] = 0.0
    adata.obs['G2M_score'] = 0.0
    adata.obs['phase'] = 'G1'
    print("WARNING  Too few genes (<10), set to 0")

# --- 4.4: Batch and Tissue ---
print("\n4.4: Batch and Tissue Keys")
print("-" * 40)

if BATCH_KEY not in adata.obs.columns:
    adata.obs[BATCH_KEY] = "ref_batch"
    print(f"OK Created '{BATCH_KEY}'")
else:
    print(f"OK '{BATCH_KEY}' exists: {adata.obs[BATCH_KEY].nunique()} batches")

if TISSUE_KEY not in adata.obs.columns:
    adata.obs[TISSUE_KEY] = "ref_tissue"
    print(f"OK Created '{TISSUE_KEY}'")
else:
    print(f"OK '{TISSUE_KEY}' exists: {adata.obs[TISSUE_KEY].nunique()} tissues")

print("\nOK All covariates calculated")


# ==============================================================================
# Step 4.5: Prepare Categorical Covariates for scArches (P1-1 FIX)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4.5: Prepare Categorical Covariates for scArches")
print("=" * 70)

print("\nNOTE Adding 'unknown' categories for scArches compatibility...")

# Tissue: add unknown_tissue
adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].astype("category")
if "unknown_tissue" not in adata.obs[TISSUE_KEY].cat.categories:
    adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].cat.add_categories(["unknown_tissue"])
    print(f"OK Added 'unknown_tissue' to {TISSUE_KEY}")

adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].fillna("unknown_tissue")

# Batch: ensure categorical
adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].astype("category")
print(f"OK {BATCH_KEY} is categorical")

print(f"\nCategorical covariate summary:")
print(f"  {BATCH_KEY}: {adata.obs[BATCH_KEY].nunique()} categories")
print(f"  {TISSUE_KEY}: {adata.obs[TISSUE_KEY].nunique()} categories")
print(f"  {L2_KEY}: {adata.obs[L2_KEY].nunique()} categories")


# ==============================================================================
# Step 5: HVG Selection with Three-Level Fallback (P1-3 FIX)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: HVG Selection (3-Level Fallback)")
print("=" * 70)

n_batches = adata.obs[BATCH_KEY].nunique()
print(f"\nBatches: {n_batches}")

hvg_method = None

# Level 1: Batch-aware seurat_v3
if n_batches > 1:
    print("\nLevel 1: Attempting batch-aware seurat_v3...")
    try:
        sc.pp.highly_variable_genes(
            adata, layer='counts', n_top_genes=N_HVG,
            batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
        )
        hvg_method = "batch-aware seurat_v3"
        print(f"OK {hvg_method}")
    except Exception as e:
        print(f"WARNING  Failed: {str(e)[:100]}")

# Level 2: Non-batch seurat_v3
if hvg_method is None:
    print("\nLevel 2: Attempting non-batch seurat_v3...")
    try:
        sc.pp.highly_variable_genes(
            adata, layer='counts', n_top_genes=N_HVG,
            flavor="seurat_v3", subset=False
        )
        hvg_method = "non-batch seurat_v3"
        print(f"OK {hvg_method}")
    except Exception as e:
        print(f"WARNING  Failed: {str(e)[:100]}")

# Level 3: cell_ranger (always works)
if hvg_method is None:
    print("\nLevel 3: Using cell_ranger...")
    sc.pp.highly_variable_genes(
        adata, layer='counts', n_top_genes=N_HVG,
        flavor="cell_ranger", subset=False
    )
    hvg_method = "cell_ranger (final fallback)"
    print(f"OK {hvg_method}")

n_hvg_initial = adata.var['highly_variable'].sum()
print(f"\nInitial HVG: {n_hvg_initial:,}")

# Force markers
if FORCE_MARKERS_IN_HVG:
    print(f"\nForcing {len(FORCED_MARKERS_BCELL)} B cell markers...")
    
    markers_in_data = [m for m in FORCED_MARKERS_BCELL if m in adata.var_names]
    for marker in markers_in_data:
        adata.var.loc[marker, 'highly_variable'] = True
    
    n_hvg_final = adata.var['highly_variable'].sum()
    n_added = n_hvg_final - n_hvg_initial
    print(f"  Added: {n_added} markers")
    print(f"  Final HVG: {n_hvg_final:,}")
else:
    n_hvg_final = n_hvg_initial

print(f"\nHVG method: {hvg_method}")


# ==============================================================================
# Step 6: Subset to HVG for Training (P0-3 FIX: NO .raw)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Subset to HVG for Training")
print("=" * 70)

print(f"\nNOTE Memory Optimization (P0-3 fix):")
print(f"  - Training uses HVG only")
print(f"  - NO .raw construction (avoids memory risk)")
print(f"  - For marker validation: use original h5ad file")

print(f"\nSubsetting to HVG...")
adata = adata[:, adata.var['highly_variable']].copy()

print(f"OK Training shape: {adata.shape}")
print(f"  HVG genes: {adata.n_vars:,}")

# Save HVG list
hvg_file = output_dir / "hvg_genes.txt"
with open(hvg_file, 'w') as f:
    for gene in adata.var_names:
        f.write(f"{gene}\n")
print(f"OK HVG list saved: {hvg_file}")


# ==============================================================================
# Step 7: scVI Training
# ==============================================================================

print("\n" + "=" * 70)
print("Step 7: scVI Training")
print("=" * 70)

print(f"\nSetting up scVI...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer=None,
    batch_key=BATCH_KEY,
    continuous_covariate_keys=['pct_counts_mt', 'stress_score', 'S_score', 'G2M_score'],
    categorical_covariate_keys=[TISSUE_KEY]
)
print("OK scVI setup complete")

print(f"\nTraining scVI...")
print(f"  Latent: {SCVI_N_LATENT}, Layers: {SCVI_N_LAYERS}, Epochs: {MAX_EPOCHS_SCVI}")

scvi_model = scvi.model.SCVI(
    adata,
    n_latent=SCVI_N_LATENT,
    n_layers=SCVI_N_LAYERS,
    n_hidden=SCVI_N_HIDDEN,
    dropout_rate=SCVI_DROPOUT_RATE,
    gene_likelihood=SCVI_GENE_LIKELIHOOD
)

train_kwargs = {
    'max_epochs': MAX_EPOCHS_SCVI,
    'batch_size': BATCH_SIZE,
    'early_stopping': True,
    'early_stopping_patience': 30,
    'plan_kwargs': {
        'lr': LEARNING_RATE,
        'weight_decay': WEIGHT_DECAY
    }
}

if gpu_available:
    train_kwargs['accelerator'] = 'gpu'
    train_kwargs['devices'] = 'auto'

scvi_model.train(**train_kwargs)
print("OK scVI training complete")


# ==============================================================================
# Step 8: scANVI Training
# ==============================================================================

print("\n" + "=" * 70)
print("Step 8: scANVI Training")
print("=" * 70)

print(f"\nTraining scANVI on '{L2_KEY}' labels...")

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key=L2_KEY,
    unlabeled_category=UNLABELED_CATEGORY
)

train_kwargs_scanvi = {
    'max_epochs': MAX_EPOCHS_SCANVI,
    'batch_size': BATCH_SIZE,
    'early_stopping': True,
    'early_stopping_patience': 30,
    'plan_kwargs': {
        'lr': LEARNING_RATE,
        'weight_decay': WEIGHT_DECAY
    }
}

if gpu_available:
    train_kwargs_scanvi['accelerator'] = 'gpu'
    train_kwargs_scanvi['devices'] = 'auto'

scanvi_model.train(**train_kwargs_scanvi)
print("OK scANVI training complete")


# ==============================================================================
# Step 9: Fit and Save UMAP Operator (P0-1 FIX - CRITICAL)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 9: Fit and Save UMAP Operator (P0-1 FIX)")
print("=" * 70)

print("\nNOTE CRITICAL: Fitting UMAP operator for query projection...")
print("This ensures reference and query will be in THE SAME UMAP space.\n")

# 1. Get reference latent
print("Extracting reference latent...")
X_ref_latent = scanvi_model.get_latent_representation(adata)
adata.obsm['X_scANVI_L2'] = X_ref_latent
print(f"OK Reference latent: {X_ref_latent.shape}")

# 2. Fit UMAP operator
print(f"\nFitting UMAP operator...")
print(f"  Parameters: n_neighbors={UMAP_N_NEIGHBORS}, min_dist={UMAP_MIN_DIST}")

umap_operator = umap.UMAP(
    n_neighbors=UMAP_N_NEIGHBORS,
    min_dist=UMAP_MIN_DIST,
    spread=UMAP_SPREAD,
    random_state=RANDOM_SEED
).fit(X_ref_latent)

adata.obsm['X_umap'] = umap_operator.embedding_
print(f"OK UMAP fitted: {adata.obsm['X_umap'].shape}")

# 3. Save UMAP operator
umap_operator_file = output_dir / "umap_operator.joblib"
joblib.dump(umap_operator, umap_operator_file)
print(f"\nOK UMAP operator saved: {umap_operator_file}")
print(f"  This file is REQUIRED for Step 2 query projection")


# ==============================================================================
# Step 10: Save Reference with UMAP (P0-2 FIX)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 10: Save Reference with UMAP for Merge (P0-2 FIX)")
print("=" * 70)

ref_with_umap_file = output_dir / "reference_with_L2_umap.h5ad"
print(f"\nSaving reference with UMAP: {ref_with_umap_file}")

adata.write_h5ad(ref_with_umap_file, compression='gzip')
file_size = ref_with_umap_file.stat().st_size / 1024**3

print(f"OK Saved ({file_size:.2f} GB)")
print(f"  Contents:")
print(f"    - X_scANVI_L2: latent representation")
print(f"    - X_umap: UMAP coordinates (for same-space merge)")
print(f"    - {L2_KEY}: L2 cell type labels")


# ==============================================================================
# Step 11: Save scANVI Model
# ==============================================================================

print("\n" + "=" * 70)
print("Step 11: Save scANVI Model")
print("=" * 70)

scanvi_model_dir = output_dir / "scanvi_existing_model"
print(f"\nSaving scANVI model: {scanvi_model_dir}")
scanvi_model.save(scanvi_model_dir, overwrite=True)
print("OK Model saved")


# ==============================================================================
# Step 12: Save Configuration
# ==============================================================================

print("\n" + "=" * 70)
print("Step 12: Save Configuration")
print("=" * 70)

config = {
    'pipeline': 'B_Cell_L2_Reference_Training',
    'version': '2.5.3_PRODUCTION_ALL_FIXES',
    'timestamp': datetime.now().isoformat(),
    'input_ref_h5ad': str(INPUT_REF_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'l3_to_l2_mapping': L3_TO_L2_MAP,
    'batch_key': BATCH_KEY,
    'tissue_key': TISSUE_KEY,
    'l2_key': L2_KEY,
    'unlabeled_category': UNLABELED_CATEGORY,
    'n_cells': int(adata.n_obs),
    'n_hvg': int(adata.n_vars),
    'hvg_method': hvg_method,
    'tissue_categories': list(adata.obs[TISSUE_KEY].cat.categories),  # P1-1 FIX
    'batch_categories': list(adata.obs[BATCH_KEY].cat.categories),
    'scvi_n_latent': SCVI_N_LATENT,
    'umap_n_neighbors': UMAP_N_NEIGHBORS,
    'umap_min_dist': UMAP_MIN_DIST,
    'random_seed': RANDOM_SEED,
    'l2_distribution': adata.obs[L2_KEY].value_counts().to_dict()
}

config_file = output_dir / 'training_config.json'
with open(config_file, 'w') as f:
    json.dump(config, f, indent=2)
print(f"OK Configuration saved: {config_file}")


# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("[OK] TRAINING COMPLETE - v2.5.3 ALL FIXES APPLIED")
print("=" * 70)

print(f"\nINFO Dataset Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  HVG: {adata.n_vars:,}")
print(f"  L2 types: {adata.obs[L2_KEY].nunique()}")

print(f"\nTARGET L2 Distribution:")
for label, count in adata.obs[L2_KEY].value_counts().items():
    print(f"  {label}: {count:,} ({count/adata.n_obs*100:.1f}%)")

print(f"\nFILE Key Outputs:")
print(f"  Model: {scanvi_model_dir}")
print(f"  UMAP operator: {umap_operator_file}")
print(f"  Reference w/ UMAP: {ref_with_umap_file}")
print(f"  HVG list: {hvg_file}")
print(f"  Config: {config_file}")

print(f"\n✨ Critical Fixes Applied:")
print(f"  [OK] P0-1: UMAP operator saved for same-space projection")
print(f"  [OK] P0-2: Reference with X_umap for merge")
print(f"  [OK] P0-3: No .raw construction (memory safe)")
print(f"  [OK] P0-4: Unknown category guaranteed")
print(f"  [OK] P1-1: unknown_tissue pre-added + saved to config")
print(f"  [OK] P1-3: 3-level HVG fallback")
print(f"  [OK] P1-4: var_names_make_unique")

print("\n" + "=" * 70)
