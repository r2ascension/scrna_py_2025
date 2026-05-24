#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
All Cells Analysis: scVI â†’ CellTypist â†’ Dual scANVI Pipeline v2.3.1 HOTFIX

Author: Clinical-Bioinformatics Team
Date: 2025-01-15
Version: 2.3.1 HOTFIX

â­ Critical Fixes in HOTFIX:
===========================
P0-2: Corrected neighborhood purity calculation
  - Now measures agreement with cell's OWN label (not dominant neighbor label)
  - Excludes self-loops
  - Uses efficient CSR indexing

P0-4: Fixed summary variable scoping
  - Separate variables for existing vs celltypist Unknown counts
  - Accurate reporting in final summary

Pipeline Overview:
==================
Input: adata_bbknn_annotated_corrected.h5ad
  â†“
Step 0: Whitelist Filtering
  â†“
Step 1: Data Preparation + Covariates Calculation
  - MT% recalculation for zeros
  - Normalized stress_score & cell_cycle_score
  â†“
Step 2: scVI Integration (scArches-ready + proper covariates)
  â†“
Step 3: CellTypist Annotation
  â†“
Step 4A: scANVI-Existing (with CORRECTED Unknown cleaning B+C + covariates fix)
  â†“
Step 4B: scANVI-CellTypist (with covariates fix)
  â†“
Step 5: Comparison & Export

Requirements:
=============
pip install "scvi-tools>=1.1.4" "lightning>=2.2,<3" scanpy celltypist mygene h5py
"""

# ==============================================================================
# Step 0: Import Libraries and Configuration
# ==============================================================================

import sys
import os
from pathlib import Path
import warnings
import json
import time
from datetime import datetime
from collections import Counter

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import gc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scanpy as sc
import scvi
import celltypist
from celltypist import models

import mygene

warnings.filterwarnings('ignore')

# Version check for scANVI fix
from packaging import version
if version.parse(scvi.__version__) < version.parse("1.1.4"):
    raise RuntimeError(
        f"âš ï¸  scvi-tools {scvi.__version__} detected. "
        "Please upgrade to >=1.1.4 for critical scANVI fixes."
    )

print(f"scanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"celltypist: {celltypist.__version__}")
print(f"Python: {sys.version}")

import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")


# ==============================================================================
# Configuration Section
# ==============================================================================

# Input/Output
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected_filtered.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX"
CELLTYPIST_MODEL_PATH = "/home/h2048/data/source/reference/celltypist_models/Cells_Lung_Airway.pkl"

# NEW: Whitelist filtering from downstream processed h5ad
DOWNSTREAM_H5AD_DIR = "/home/h2048/data/core20260115"
DOWNSTREAM_H5AD_GLOB = "**/*.h5ad"
CELLID_SET_MODE = "union"  # "union" or "intersection"
ENABLE_WHITELIST_FILTERING = True

# Keys
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"  # NEW: will be used as categorical covariate
EXISTING_CELLTYPE_KEY = "cell_type"

# Required covariates for scVI/scANVI models
REQUIRED_COVARIATES = [BATCH_KEY, TISSUE_KEY, 'pct_counts_mt', 'stress_score', 'S_score', 'G2M_score']

# Model loading (set to False for clean retrain)
LOAD_SCVI_IF_EXISTS = False
LOAD_SCANVI_EXISTING_IF_EXISTS = False
LOAD_SCANVI_CELLTYPIST_IF_EXISTS = False

# HVG
USE_HVG_FOR_SCVI = True
N_HVG_SCVI = 6000
HVG_FLAVOR = "seurat_v3"

# NEW: Force important markers into HVG
FORCE_MARKERS_IN_HVG = True
FORCE_MARKERS_CASE_INSENSITIVE = True

# Lineage anchors (å¿…å…¥ï¼Œä¸å½±å“ä¸»ç»“æž„)
FORCED_MARKERS_LINEAGE = [
    # Epithelial
    "FXYD3","EPCAM","ELF3","IGFBP2","SERPINF1","TSPAN1","SCGB1A1",
    "AGER","SFTPC","FOXJ1","KRT5","MUC5B","KRT8",
    # Immune general / Lymphoid
    "CD53","PTPRC","CORO1A","ISG20","CCL5",
    # B cells
    "MS4A1","TNFRSF17","CD19","CD79A","SDC1",
    # T cells
    "CD40LG","TNFRSF25","CD28","CD4","CD3D","CD3E","CD2","TRBC2",
    "CD8A","CD8B","TRGC2",
    # Myeloid
    "FCER1G","C1ORF162","CLEC7A","CD1C","CD86","CD14","XCR1","HLA-DRA",
    # Stromal / Fibroblast / SMC
    "COL1A2","DCN","MFAP4","LUM","COL6A3","CFD","COL1A1","PDGFRA",
    "MXRA8","NBL1","VCAN","LEPR","MYH11","TINAGL1","PLN","DES",
    "ACTA2","CNN1","TAGLN",
    # Endothelial
    "CLDN5","ECSCR","CLEC14A","VWF","PECAM1","ACKR1","PTPRB","PDE2A",
    "PLAT","GJA5","SPARCL1","AQP1","RNASE1","MMRN1","CCL21","TFF3",
]

# State markers (å¯é€‰ï¼Œé»˜è®¤ä¸åŠ å…¥é¿å… cell cycle ä¸»å¯¼)
FORCED_MARKERS_STATE = [
    "MKI67","TOP2A","TK1","CENPW"  # Proliferation
]

# Default: only lineage (state ä½œä¸ºå¼€å…³)
INCLUDE_STATE_MARKERS = False  # Set True if you want proliferation axis
FORCED_MARKERS = FORCED_MARKERS_LINEAGE + (FORCED_MARKERS_STATE if INCLUDE_STATE_MARKERS else [])

# NEW: Stress signature genes
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

# Cell cycle genes (local, no network dependency)
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

# Semi-supervised for CellTypist-based scANVI only
LOW_CONFIDENCE_THRESHOLD = 0.5
RARE_TYPE_THRESHOLD = 10

# NEW: Unknown cleaning for existing labels (Strategy B+C)
EXISTING_LABEL_PURITY_THRESHOLD = 0.5  # Neighborhood purity
EXISTING_LABEL_MT_THRESHOLD = 20  # mt% cutoff
EXISTING_LABEL_STRESS_PERCENTILE = 95  # stress percentile cutoff

# scVI (FIXED: separate from scANVI parameters)
SCVI_N_LATENT = 150  # â­ Set to 150 for stability
SCVI_N_LAYERS = 2  # â­ Reduced from 4 for scArches compatibility
SCVI_DROPOUT_RATE = 0.2  # â­ Increased from 0.1 for regularization
SCVI_MAX_EPOCHS = 400
SCVI_LEARNING_RATE = 1e-3
SCVI_BATCH_SIZE = 256
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 45

# â­ NEW: scArches-ready parameters (CRITICAL for mapping)
SCVI_ENCODE_COVARIATES = True  # Essential for reference mapping
SCVI_USE_LAYER_NORM = "both"
SCVI_USE_BATCH_NORM = "none"

# scANVI
SCANVI_MAX_EPOCHS = 200
SCANVI_LEARNING_RATE = 5e-4
SCANVI_BATCH_SIZE = 256
SCANVI_EARLY_STOPPING_PATIENCE = 30
SCANVI_N_SAMPLES_PER_LABEL = 2000  # â­ NEW: balance class imbalance

# UMAP
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
UMAP_N_NEIGHBORS = 50

# Visualization
DPI = 300
FIGURE_FORMAT = 'pdf'

# Reproducibility
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
scvi.settings.seed = RANDOM_SEED

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)
(output_dir / "models").mkdir(exist_ok=True)

print(f"\nOutput directory: {output_dir}")

# Save configuration
config = {
    'version': '2.3.1-HOTFIX',
    'input_h5ad': str(INPUT_H5AD),
    'output_dir': str(OUTPUT_DIR),
    'batch_key': BATCH_KEY,
    'tissue_key': TISSUE_KEY,
    'existing_celltype_key': EXISTING_CELLTYPE_KEY,
    'whitelist_filtering': ENABLE_WHITELIST_FILTERING,
    'downstream_h5ad_dir': str(DOWNSTREAM_H5AD_DIR) if ENABLE_WHITELIST_FILTERING else None,
    'use_hvg_for_scvi': USE_HVG_FOR_SCVI,
    'n_hvg_scvi': N_HVG_SCVI,
    'force_markers': FORCE_MARKERS_IN_HVG,
    'n_forced_markers': len(FORCED_MARKERS),
    'scvi_encode_covariates': SCVI_ENCODE_COVARIATES,
    'scvi_n_latent': SCVI_N_LATENT,
    'scvi_n_layers': SCVI_N_LAYERS,
    'low_confidence_threshold': LOW_CONFIDENCE_THRESHOLD,
    'rare_type_threshold': RARE_TYPE_THRESHOLD,
    'existing_unknown_cleaning': True,
    'dual_scanvi': True,
    'n_samples_per_label': SCANVI_N_SAMPLES_PER_LABEL,
    'random_seed': RANDOM_SEED,
    'hotfixes': ['P0-2_corrected_purity', 'P0-4_fixed_summary_vars'],
    'timestamp': datetime.now().isoformat()
}

with open(output_dir / 'pipeline_config_v2.3.1_HOTFIX.json', 'w') as f:
    json.dump(config, f, indent=2)

print("âœ“ Configuration saved")


# ==============================================================================
# Utility Functions
# ==============================================================================

def _read_obs_index_from_h5ad(h5ad_path: str):
    """
    Read obs_names from h5ad without loading the whole AnnData.
    Memory efficient for whitelist collection.
    """
    try:
        import h5py
        with h5py.File(h5ad_path, "r") as f:
            if "obs" in f and "_index" in f["obs"]:
                arr = f["obs"]["_index"][()]
                if hasattr(arr, "dtype") and arr.dtype.kind in ("S", "O"):
                    arr = np.array([x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr])
                else:
                    arr = arr.astype(str)
                return arr.tolist()
    except Exception:
        pass

    # Fallback to scanpy
    import scanpy as sc
    ad = sc.read_h5ad(h5ad_path, backed="r")
    ids = ad.obs_names.astype(str).tolist()
    try:
        if getattr(ad, "isbacked", False) and hasattr(ad, "file"):
            ad.file.close()
    except Exception:
        pass
    return ids


def collect_cell_ids_from_folder(folder: str, pattern: str="**/*.h5ad", mode: str="union"):
    """
    Collect cell IDs from downstream processed h5ad files.
    
    mode:
      - 'union': keep cells that appear in ANY downstream h5ad
      - 'intersection': keep cells that appear in ALL downstream h5ad
    """
    folder = Path(folder)
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No h5ad found under: {folder} (pattern={pattern})")

    cell_set = None
    for i, fp in enumerate(files, 1):
        ids = _read_obs_index_from_h5ad(str(fp))
        ids = set(map(str, ids))

        if cell_set is None:
            cell_set = ids
        else:
            if mode.lower() == "intersection":
                cell_set &= ids
            else:
                cell_set |= ids

        if i % 20 == 0:
            print(f"  Processed {i}/{len(files)} files, current cells: {len(cell_set):,}")

    print(f"âœ“ Downstream h5ad files: {len(files)}")
    print(f"âœ“ Collected cell ids ({mode}): {len(cell_set):,}")
    return cell_set


def force_include_markers_in_hvg(
    adata,
    n_top_genes: int,
    forced_markers: list,
    symbol_col: str = "symbol_base",
    case_insensitive: bool = True,
    hvg_col: str = "highly_variable"
):
    """
    Force include biological markers into HVG set while maintaining target size.
    Strategy: add missing markers, then drop worst-ranked non-forced HVGs.
    """
    if hvg_col not in adata.var.columns:
        raise ValueError(f"'{hvg_col}' not found. Run sc.pp.highly_variable_genes first.")

    # Match genes
    if symbol_col in adata.var.columns:
        gene_vec = adata.var[symbol_col].astype(str)
    else:
        gene_vec = adata.var_names.astype(str)

    if case_insensitive:
        gene_vec_cmp = gene_vec.str.upper()
        forced_cmp = set([str(x).upper() for x in forced_markers])
    else:
        gene_vec_cmp = gene_vec
        forced_cmp = set([str(x) for x in forced_markers])

    forced_mask = gene_vec_cmp.isin(forced_cmp).values
    present_forced = int(forced_mask.sum())

    hvg_mask = adata.var[hvg_col].values.astype(bool)
    original_hvg_n = int(hvg_mask.sum())

    # Add forced markers
    missing_forced_mask = forced_mask & (~hvg_mask)
    add_n = int(missing_forced_mask.sum())
    hvg_mask = hvg_mask | forced_mask

    # Trim to n_top_genes if exceeded
    if int(hvg_mask.sum()) > n_top_genes:
        # Use ranking to drop worst non-forced HVGs
        if "highly_variable_rank" in adata.var.columns:
            rank = adata.var["highly_variable_rank"].fillna(1e18).values.astype(float)
            sort_key = rank
        elif "variances_norm" in adata.var.columns:
            vn = adata.var["variances_norm"].fillna(-1e18).values.astype(float)
            sort_key = -vn
        else:
            print("âš ï¸  No HVG ranking column found; keep expanded HVG set.")
            adata.var[hvg_col] = hvg_mask
            return {
                "original_hvg_n": original_hvg_n,
                "present_forced": present_forced,
                "added_forced": add_n,
                "final_hvg_n": int(hvg_mask.sum()),
                "trimmed": 0
            }

        candidates_drop = np.where(hvg_mask & (~forced_mask))[0]
        candidates_drop = candidates_drop[np.argsort(sort_key[candidates_drop])[::-1]]

        need_drop = int(hvg_mask.sum()) - n_top_genes
        drop_idx = candidates_drop[:need_drop]
        hvg_mask[drop_idx] = False

        trimmed = len(drop_idx)
    else:
        trimmed = 0

    adata.var[hvg_col] = hvg_mask

    return {
        "original_hvg_n": original_hvg_n,
        "present_forced": present_forced,
        "added_forced": add_n,
        "final_hvg_n": int(hvg_mask.sum()),
        "trimmed": trimmed
    }


def calculate_label_agreement_purity(adata, label_key, neighbors_key='neighbors'):
    """
    â­ HOTFIX P0-2: Corrected neighborhood purity calculation
    
    Calculate the fraction of neighbors that share the SAME label as each cell.
    This is the correct definition for detecting mislabeled/ambiguous cells.
    
    Previous incorrect version measured "dominant neighbor label fraction",
    which could be high even for mislabeled cells if neighbors form a pure cluster.
    
    Parameters:
    -----------
    adata : AnnData
    label_key : str
        Column in adata.obs containing cell type labels
    neighbors_key : str
        Prefix for neighbors graph in obsp (default: 'neighbors')
    
    Returns:
    --------
    np.ndarray
        Purity scores for each cell (0 to 1)
        High purity = cell's label matches most neighbors (likely correct)
        Low purity = cell's label differs from neighbors (likely mislabeled/ambiguous)
    """
    # Get connectivity matrix
    connectivities = adata.obsp[f'{neighbors_key}_connectivities'].tocsr()
    
    # Convert labels to categorical codes for fast comparison
    labels = pd.Categorical(adata.obs[label_key].astype(str))
    label_codes = labels.codes
    
    # Use CSR format's indptr and indices for efficient row iteration
    indptr = connectivities.indptr
    indices = connectivities.indices
    
    purities = np.empty(adata.n_obs, dtype=np.float32)
    
    for i in range(adata.n_obs):
        # Get neighbors for cell i
        neighbors = indices[indptr[i]:indptr[i+1]]
        
        # Exclude self (if present in connectivity matrix)
        neighbors = neighbors[neighbors != i]
        
        if neighbors.size == 0:
            purities[i] = 0.0
        else:
            # Calculate fraction of neighbors with same label as cell i
            same_label = label_codes[neighbors] == label_codes[i]
            purities[i] = np.mean(same_label)
    
    return purities


# ==============================================================================
# Step 1: Data Loading and Preparation
# ==============================================================================

print("\n" + "="*70)
print("Step 1: Data Loading and Preparation")
print("="*70)

print(f"\nLoading data from: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"\nData loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found")

if EXISTING_CELLTYPE_KEY not in adata.obs.columns:
    raise ValueError(f"Existing cell type key '{EXISTING_CELLTYPE_KEY}' not found")

print(f"\nBatch distribution ({BATCH_KEY}):")
print(adata.obs[BATCH_KEY].value_counts())

print(f"\nExisting cell type distribution ({EXISTING_CELLTYPE_KEY}):")
print(adata.obs[EXISTING_CELLTYPE_KEY].value_counts())


# ==============================================================================
# NEW: Whitelist Filtering
# ==============================================================================

if ENABLE_WHITELIST_FILTERING:
    print("\n" + "-"*70)
    print("NEW: Cleaning cells by downstream whitelist")
    print("-"*70)

    valid_cell_ids = collect_cell_ids_from_folder(
        DOWNSTREAM_H5AD_DIR,
        pattern=DOWNSTREAM_H5AD_GLOB,
        mode=CELLID_SET_MODE
    )

    mask_keep = adata.obs_names.astype(str).isin(valid_cell_ids)
    n_before = adata.n_obs
    n_keep = mask_keep.sum()
    print(f"  Keep: {n_keep:,}/{n_before:,} cells ({n_keep/n_before*100:.2f}%)")
    print(f"  Remove: {n_before - n_keep:,} cells")
    
    adata = adata[mask_keep].copy()
    print(f"âœ“ Filtered to {adata.n_obs:,} cells")


# Extract raw counts
print("\n" + "-"*70)
print("Extracting raw counts")
print("-"*70)

if 'counts' in adata.layers:
    print("âœ“ Found raw counts in layers['counts']")
elif hasattr(adata, 'raw') and adata.raw is not None:
    adata.layers['counts'] = adata.raw.X.copy()
    print("âœ“ Extracted counts from .raw.X")
else:
    raise ValueError("Cannot find raw counts")

# Set adata.raw.X (after whitelist filtering)
print("\nSetting adata.raw.X (shared memory)...")
adata_raw = sc.AnnData(
    X=adata.layers["counts"],
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
adata.raw = adata_raw
print("âœ“ adata.raw.X set")


# ==============================================================================
# PATCH 1: Normalize covariate column names & MT% recalculation
# ==============================================================================

print("\n" + "-"*70)
print("PATCH 1: Normalizing covariate column names")
print("-"*70)

# Handle mt% column name variations
mt_col_found = False
if 'percent.mt' in adata.obs.columns:
    adata.obs['pct_counts_mt'] = adata.obs['percent.mt']
    print("âœ“ Renamed 'percent.mt' â†’ 'pct_counts_mt'")
    mt_col_found = True
elif 'percent_mt' in adata.obs.columns:
    adata.obs['pct_counts_mt'] = adata.obs['percent_mt']
    print("âœ“ Renamed 'percent_mt' â†’ 'pct_counts_mt'")
    mt_col_found = True
elif 'pct_counts_mt' in adata.obs.columns:
    print("âœ“ 'pct_counts_mt' already exists")
    mt_col_found = True

# Check for zeros and recalculate if needed
if mt_col_found:
    n_zeros = (adata.obs['pct_counts_mt'] == 0).sum()
    pct_zeros = n_zeros / adata.n_obs * 100
    print(f"\nâš ï¸  MT% diagnostics:")
    print(f"  Cells with mt%=0: {n_zeros:,} ({pct_zeros:.1f}%)")
    print(f"  Mean mt%: {adata.obs['pct_counts_mt'].mean():.2f}%")
    print(f"  Median mt%: {adata.obs['pct_counts_mt'].median():.2f}%")
    
    # If >5% cells have mt%=0, recalculate for those cells
    if pct_zeros > 5:
        print(f"\nâš ï¸  {pct_zeros:.1f}% cells have mt%=0 (suspicious)")
        print("  Recalculating mt% from counts layer...")
        
        # Identify MT genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        n_mt_genes = adata.var['mt'].sum()
        print(f"  Found {n_mt_genes} MT genes")
        
        if n_mt_genes > 0:
            # Calculate for cells with mt%=0
            zero_mask = adata.obs['pct_counts_mt'] == 0
            if zero_mask.sum() > 0:
                # Use counts layer
                counts = adata.layers['counts']
                if issparse(counts):
                    total_counts = np.asarray(counts.sum(axis=1)).flatten()
                    mt_counts = np.asarray(counts[:, adata.var['mt'].values].sum(axis=1)).flatten()
                else:
                    total_counts = counts.sum(axis=1)
                    mt_counts = counts[:, adata.var['mt'].values].sum(axis=1)
                
                # Recalculate mt%
                mt_pct_recalc = np.zeros(adata.n_obs)
                nonzero_total = total_counts > 0
                mt_pct_recalc[nonzero_total] = (mt_counts[nonzero_total] / total_counts[nonzero_total]) * 100
                
                # Update only zero cells
                adata.obs.loc[zero_mask, 'pct_counts_mt'] = mt_pct_recalc[zero_mask]
                
                n_updated = (adata.obs.loc[zero_mask, 'pct_counts_mt'] > 0).sum()
                print(f"âœ“ Recalculated mt% for {n_updated:,} cells")
                print(f"  New mean mt%: {adata.obs['pct_counts_mt'].mean():.2f}%")
                print(f"  New median mt%: {adata.obs['pct_counts_mt'].median():.2f}%")
                print(f"  Remaining zeros: {(adata.obs['pct_counts_mt']==0).sum():,}")
        else:
            print("  âš ï¸  No MT genes found, cannot recalculate")
else:
    # No mt column found, calculate from scratch
    print("âš ï¸  No mt% column found, calculating from scratch...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, layer='counts')
    print("âœ“ Calculated 'pct_counts_mt'")


# ==============================================================================
# Gene Name Handling (MUST come before stress/cycle scoring)
# ==============================================================================

print("\n" + "-"*70)
print("Gene Name Handling")
print("-"*70)

sample_gene = str(adata.var_names[0])
print(f"Sample gene: {sample_gene}")

if 'symbol' not in adata.var.columns:
    if sample_gene.startswith('ENSG'):
        print("âœ“ Detected ENSEMBL IDs")
        
        symbol_candidates = ['gene_symbols', 'feature_name', 'gene_name']
        existing_symbol_col = None
        for col in symbol_candidates:
            if col in adata.var.columns:
                existing_symbol_col = col
                break
        
        if existing_symbol_col:
            adata.var['symbol'] = adata.var[existing_symbol_col].astype(str)
            print(f"  Using existing column: {existing_symbol_col}")
        else:
            print("  Using mygene for conversion...")
            mg = mygene.MyGeneInfo()
            ensembl_ids = adata.var_names.tolist()
            
            query_result = mg.querymany(
                ensembl_ids,
                scopes='ensembl.gene',
                fields='symbol',
                species='human',
                returnall=True,
                verbose=False
            )
            
            ensembl_to_symbol = {}
            for item in query_result['out']:
                if 'symbol' in item:
                    ensembl_to_symbol[item['query']] = item['symbol']
            
            adata.var['symbol'] = [ensembl_to_symbol.get(e, e) for e in adata.var_names]
            print(f"  Converted: {len(ensembl_to_symbol):,}/{len(ensembl_ids):,}")
    else:
        adata.var['symbol'] = adata.var_names.astype(str)
else:
    print("âœ“ Symbol column exists")

adata.var['symbol_base'] = adata.var['symbol'].str.replace(r'-\d+$', '', regex=True)
print("âœ“ Created symbol_base column")


# ==============================================================================
# PATCH 3: Calculate Covariates (Normalized)
# ==============================================================================

print("\n" + "-"*70)
print("PATCH 3: Calculating Covariates (Normalized)")
print("-"*70)

# 2. Stress score (è§„èŒƒåŒ–ç‰ˆæœ¬)
print("\nCalculating stress_score...")
adata_temp = sc.AnnData(
    X=adata.layers['counts'].copy(),
    var=adata.var.copy()
)
sc.pp.normalize_total(adata_temp, target_sum=1e4)
sc.pp.log1p(adata_temp)

# Match using var_names (since already confirmed as symbols)
stress_genes_in_data = [g for g in STRESS_SIGNATURE_GENES if g in adata_temp.var_names]
print(f"  Found {len(stress_genes_in_data)}/{len(STRESS_SIGNATURE_GENES)} stress genes")

if len(stress_genes_in_data) > 10:
    stress_gene_idx = [list(adata_temp.var_names).index(g) for g in stress_genes_in_data]
    if issparse(adata_temp.X):
        stress_expr = np.asarray(adata_temp.X[:, stress_gene_idx].mean(axis=1)).flatten()
    else:
        stress_expr = adata_temp.X[:, stress_gene_idx].mean(axis=1)
    adata.obs['stress_score'] = stress_expr
    print(f"  Mean stress_score: {adata.obs['stress_score'].mean():.3f}")
    print(f"  Median stress_score: {adata.obs['stress_score'].median():.3f}")
else:
    adata.obs['stress_score'] = 0.0
    print("  âš ï¸  Too few stress genes found, set stress_score=0")

del adata_temp
gc.collect()

# 3. Cell cycle scores (è§„èŒƒåŒ–ç‰ˆæœ¬ + æœ¬åœ°åŸºå› åˆ—è¡¨)
print("\nCalculating cell_cycle_scores...")

s_genes = [x for x in S_GENES if x in adata.var_names]
g2m_genes = [x for x in G2M_GENES if x in adata.var_names]

print(f"  Found S genes: {len(s_genes)}/{len(S_GENES)}")
print(f"  Found G2M genes: {len(g2m_genes)}/{len(G2M_GENES)}")

if len(s_genes) > 10 and len(g2m_genes) > 10:
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
    
    print(f"  S_score mean: {adata.obs['S_score'].mean():.3f}")
    print(f"  G2M_score mean: {adata.obs['G2M_score'].mean():.3f}")
    print(f"  Phase distribution:")
    print(adata.obs['phase'].value_counts())
else:
    adata.obs['S_score'] = 0.0
    adata.obs['G2M_score'] = 0.0
    adata.obs['phase'] = 'G1'
    print(f"  âš ï¸  Too few cell cycle genes found")
    print("  Set scores=0")

print("\nâœ“ All covariates calculated")
print("\nâœ“ Data preparation complete")


# ==============================================================================
# Step 2: scVI Integration
# ==============================================================================

print("\n" + "="*70)
print("Step 2: scVI Integration (scArches-ready)")
print("="*70)

scvi_model_dir = output_dir / "models" / "scvi_model"
scvi_model_exists = scvi_model_dir.exists() and (scvi_model_dir / "model.pt").exists()
hvg_file = scvi_model_dir / "hvg_genes.txt"

hvg_indices = None
actual_hvg_count = adata.n_vars
hvg_selection_method = "all_genes"
hvg_genes_list = None

if USE_HVG_FOR_SCVI:
    print(f"\nâ­ Using HVG ({N_HVG_SCVI} genes)")
    
    # Load HVG from file if exists
    if scvi_model_exists and hvg_file.exists() and LOAD_SCVI_IF_EXISTS:
        print("\nLoading HVG from file...")
        try:
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_genes_in_data = [g for g in hvg_genes_saved if g in set(adata.var_names)]
            
            if len(hvg_genes_in_data) >= len(hvg_genes_saved) * 0.95:
                hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_in_data))
                hvg_genes_list = hvg_genes_in_data
                actual_hvg_count = len(hvg_genes_list)
                hvg_selection_method = "loaded_from_file"
                print(f"âœ“ Loaded {actual_hvg_count:,} HVG genes")
            else:
                hvg_genes_list = None
        except Exception as e:
            print(f"âš ï¸  Failed to load HVG: {e}")
            hvg_genes_list = None
    
    # Compute HVG if not loaded
    if hvg_genes_list is None:
        print("\nComputing HVG...")
        try:
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                batch_key=BATCH_KEY, flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "batch-aware"
        except:
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "non-batch-aware"
        
        hvg_mask = adata.var["highly_variable"].values
        hvg_indices = np.flatnonzero(hvg_mask)
        hvg_genes_list = adata.var_names[hvg_indices].tolist()
        actual_hvg_count = len(hvg_genes_list)
        print(f"âœ“ Selected {actual_hvg_count:,} HVG")
        
        # NEW: Force markers into HVG
        if FORCE_MARKERS_IN_HVG:
            print("\n" + "-"*70)
            print("NEW: Forcing biological markers into HVG")
            print("-"*70)
            
            adata.var["highly_variable_original"] = adata.var["highly_variable"].values
            
            report = force_include_markers_in_hvg(
                adata,
                n_top_genes=N_HVG_SCVI,
                forced_markers=FORCED_MARKERS,
                symbol_col="symbol_base",
                case_insensitive=FORCE_MARKERS_CASE_INSENSITIVE,
                hvg_col="highly_variable"
            )
            print("âœ“ Forced markers into HVG:")
            for k, v in report.items():
                print(f"  {k}: {v}")
            
            # Refresh HVG after forcing
            hvg_mask = adata.var["highly_variable"].values
            hvg_indices = np.flatnonzero(hvg_mask)
            hvg_genes_list = adata.var_names[hvg_indices].tolist()
            actual_hvg_count = len(hvg_genes_list)
            print(f"âœ“ Final HVG count: {actual_hvg_count:,}")
    
    # Construct scVI AnnData with HVG
    print("\nConstructing scVI AnnData with HVG...")
    adata_scvi = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[REQUIRED_COVARIATES].copy(),
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scvi.var_names = pd.Index(hvg_genes_list)
    print(f"âœ“ scVI AnnData: {adata_scvi.n_obs:,} Ã— {adata_scvi.n_vars:,}")
else:
    adata_scvi = adata.copy()
    adata_scvi.X = adata.layers['counts'].copy()

# Verify covariate keys exist
missing_covariates = [k for k in REQUIRED_COVARIATES if k not in adata_scvi.obs.columns]
if missing_covariates:
    raise ValueError(f"Missing covariate columns: {missing_covariates}")

# Load or train scVI
print("\n" + "-"*70)
print("scVI Model: Load or Train (scArches-ready)")
print("-"*70)

if scvi_model_exists and LOAD_SCVI_IF_EXISTS:
    print(f"âœ“ Loading existing scVI model...")
    try:
        scvi_model = scvi.model.SCVI.load(scvi_model_dir, adata=adata_scvi)
        print("âœ“ scVI model loaded")
    except Exception as e:
        print(f"âš ï¸  Load failed: {e}")
        scvi_model_exists = False

if not scvi_model_exists or not LOAD_SCVI_IF_EXISTS:
    print("Training new scVI model (scArches-ready)...")
    
    # Setup with covariates
    scvi.model.SCVI.setup_anndata(
        adata_scvi, 
        layer=None, 
        batch_key=BATCH_KEY,
        categorical_covariate_keys=[TISSUE_KEY],
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES if k not in [BATCH_KEY, TISSUE_KEY]]
    )
    
    # Initialize model with scArches-ready parameters
    scvi_model = scvi.model.SCVI(
        adata_scvi, 
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        dropout_rate=SCVI_DROPOUT_RATE,
        gene_likelihood="nb",
        encode_covariates=SCVI_ENCODE_COVARIATES,
        use_layer_norm=SCVI_USE_LAYER_NORM,
        use_batch_norm=SCVI_USE_BATCH_NORM
    )
    
    print(f"âœ“ Model initialized:")
    print(f"  n_latent: {SCVI_N_LATENT}")
    print(f"  n_layers: {SCVI_N_LAYERS}")
    print(f"  encode_covariates: {SCVI_ENCODE_COVARIATES}")
    print(f"  use_layer_norm: {SCVI_USE_LAYER_NORM}")
    
    # âœ… FIXED: Use SCVI_* parameters
    train_kwargs = {
        'max_epochs': SCVI_MAX_EPOCHS,
        'batch_size': SCVI_BATCH_SIZE,
        'early_stopping': SCVI_EARLY_STOPPING,
        'early_stopping_patience': SCVI_EARLY_STOPPING_PATIENCE,
        'plan_kwargs': {
            'lr': SCVI_LEARNING_RATE
        },
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scvi_model.train(**train_kwargs)
    print("âœ“ scVI training complete")
    
    scvi_model.save(scvi_model_dir, overwrite=True)
    print(f"âœ“ Model saved")
    
    if USE_HVG_FOR_SCVI and hvg_genes_list:
        pd.Series(hvg_genes_list).to_csv(hvg_file, index=False, header=False)
        print("âœ“ HVG list saved")

# Generate latent
print("\nGenerating scVI latent...")
latent_df = pd.DataFrame(
    scvi_model.get_latent_representation(),
    index=adata_scvi.obs_names
)
latent_aligned = latent_df.reindex(adata.obs_names)
if latent_aligned.isna().any().any():
    raise RuntimeError("scVI latent reindex produced NaNs")
adata.obsm['X_scvi'] = latent_aligned.to_numpy().astype(np.float32)
print(f"âœ“ scVI latent: {adata.obsm['X_scvi'].shape}")

del adata_scvi
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scvi)...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS, 
                key_added='neighbors_scvi', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, 
           neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scvi', 
             neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
print(f"âœ“ Leiden: {adata.obs['leiden_scvi'].nunique()} clusters")

print("âœ“ scVI integration complete")


# ==============================================================================
# Step 3: CellTypist Annotation
# ==============================================================================

print("\n" + "="*70)
print("Step 3: CellTypist Annotation")
print("="*70)

print(f"\nLoading CellTypist model...")
celltypist_model = models.Model.load(CELLTYPIST_MODEL_PATH)

# Get model features
model_features = None
for attr in ['features', 'genes', 'var_names']:
    if hasattr(celltypist_model, attr):
        mf = getattr(celltypist_model, attr)
        if mf is not None and len(mf) > 100:
            model_features = pd.Index(mf).astype(str)
            break

if model_features is None:
    raise RuntimeError("Cannot retrieve CellTypist features")

# Match genes
print("\nMatching genes using symbol_base...")
common_genes = adata.var['symbol_base'].isin(model_features)
n_common = common_genes.sum()
print(f"  Overlapping: {n_common:,}/{len(model_features)}")

if n_common < 500:
    raise RuntimeError(f"Too few overlapping genes: {n_common}")

# Construct CellTypist AnnData (memory efficient)
print(f"\nConstructing CellTypist AnnData (only {n_common:,} genes)...")
X_subset = adata.layers["counts"][:, common_genes.values]

adata_celltypist = sc.AnnData(
    X=X_subset.copy(),
    obs=adata.obs.copy(),
    var=adata.var.loc[common_genes].copy()
)
adata_celltypist.var_names = adata.var.loc[common_genes, 'symbol_base'].values

sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
sc.pp.log1p(adata_celltypist)
print("âœ“ Data prepared")

# Run CellTypist
print("\nRunning CellTypist prediction...")
predictions = celltypist.annotate(adata_celltypist, model=celltypist_model, majority_voting=True)

# Extract results
predicted_labels = predictions.predicted_labels
pred_df = predicted_labels.reindex(adata.obs_names)
if pred_df.isna().any().any():
    raise RuntimeError("CellTypist reindex produced NaNs")

adata.obs['celltypist_predicted'] = pred_df.predicted_labels.astype(str).values
if 'majority_voting' in pred_df.columns:
    adata.obs['celltypist_majority_voting'] = pred_df.majority_voting.astype(str).values
else:
    adata.obs['celltypist_majority_voting'] = pred_df.predicted_labels.astype(str).values

# Confidence
conf_candidates = ['conf_score', 'confidence', 'confidence_score', 'prob']
confidence_column = next((col for col in conf_candidates if col in predicted_labels.columns), None)
if confidence_column:
    adata.obs['celltypist_conf_score'] = pred_df[confidence_column].reindex(adata.obs_names).values
else:
    adata.obs['celltypist_conf_score'] = 1.0

del adata_celltypist, X_subset
gc.collect()

print(f"\nCellTypist summary:")
print(f"  Unique types: {adata.obs['celltypist_majority_voting'].nunique()}")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")

print("âœ“ CellTypist complete")


# ==============================================================================
# Step 4A: scANVI-Existing (HOTFIX P0-2: Corrected Unknown cleaning)
# ==============================================================================

print("\n" + "="*70)
print("Step 4A: scANVI-Existing (Using cell_type metadata)")
print("="*70)

print(f"\nâ­ Training scANVI using '{EXISTING_CELLTYPE_KEY}' metadata")
print(f"   With Unknown cleaning (Strategy B+C: CORRECTED purity + quality)")
print(f"   With n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")
print(f"   With covariates (FIXED)")

scanvi_existing_model_dir = output_dir / "models" / "scanvi_existing_model"
scanvi_existing_exists = scanvi_existing_model_dir.exists() and (scanvi_existing_model_dir / "model.pt").exists()

# HOTFIX P0-2: Clean existing labels with CORRECTED Unknown strategy
print("\nCleaning existing labels (Strategy B+C - HOTFIX)...")
labels_existing = adata.obs[EXISTING_CELLTYPE_KEY].astype(str).copy()
n_original_existing = len(labels_existing)  # â­ P0-4 FIX: unique variable name

# Strategy B: CORRECTED Neighborhood purity
print("  Computing neighborhood purity (CORRECTED)...")
print("    (measuring agreement with cell's OWN label, not dominant neighbor)")

purity_array = calculate_label_agreement_purity(adata, EXISTING_CELLTYPE_KEY, 'neighbors_scvi')
adata.obs['label_purity_existing'] = purity_array

low_purity_mask = purity_array < EXISTING_LABEL_PURITY_THRESHOLD
n_low_purity_existing = low_purity_mask.sum()
labels_existing[low_purity_mask] = "Unknown"

print(f"    Threshold: {EXISTING_LABEL_PURITY_THRESHOLD}")
print(f"    Mean purity: {purity_array.mean():.3f}")
print(f"    Median purity: {np.median(purity_array):.3f}")
print(f"    Low purity â†’ Unknown: {n_low_purity_existing:,} ({n_low_purity_existing/n_original_existing*100:.1f}%)")

# Strategy C: Quality filtering
quality_mask = (
    (adata.obs['pct_counts_mt'] > EXISTING_LABEL_MT_THRESHOLD) |
    (adata.obs['stress_score'] > np.percentile(adata.obs['stress_score'], EXISTING_LABEL_STRESS_PERCENTILE))
)
n_low_quality_existing = quality_mask.sum()
labels_existing[quality_mask] = "Unknown"

print(f"    MT% threshold: {EXISTING_LABEL_MT_THRESHOLD}%")
print(f"    Stress percentile: {EXISTING_LABEL_STRESS_PERCENTILE}th")
print(f"    Low quality â†’ Unknown: {n_low_quality_existing:,} ({n_low_quality_existing/n_original_existing*100:.1f}%)")

# Summary
n_unknown_existing = (labels_existing == "Unknown").sum()  # â­ P0-4 FIX: unique variable name
n_labeled_existing = n_original_existing - n_unknown_existing

print(f"\n  Summary:")
print(f"    Total: {n_original_existing:,}")
print(f"    Final labeled: {n_labeled_existing:,} ({n_labeled_existing/n_original_existing*100:.1f}%)")
print(f"    Final Unknown: {n_unknown_existing:,} ({n_unknown_existing/n_original_existing*100:.1f}%)")

adata.obs['scanvi_label_existing_cleaned'] = pd.Categorical(labels_existing)

# Prepare data (PATCH 2B: Include covariates)
print("\nPreparing scANVI-Existing data...")

if USE_HVG_FOR_SCVI:
    if hvg_indices is None:
        if 'highly_variable' in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            hvg_indices = np.flatnonzero(hvg_mask)
        elif hvg_file.exists():
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_saved))
    
    adata_scanvi_existing = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_existing_cleaned']].copy(),  # â­ FIXED
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scanvi_existing.var_names = adata.var_names[hvg_indices]
else:
    adata_scanvi_existing = sc.AnnData(
        X=adata.layers["counts"].copy(),
        obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_existing_cleaned']].copy(),  # â­ FIXED
        var=adata.var.copy()
    )

print(f"âœ“ scANVI-Existing AnnData: {adata_scanvi_existing.n_obs:,} Ã— {adata_scanvi_existing.n_vars:,}")

# Ensure 'Unknown' category exists
labels_cat = adata_scanvi_existing.obs['scanvi_label_existing_cleaned'].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi_existing.obs['scanvi_label_existing_cleaned'] = labels_cat

# Load or train
print("\n" + "-"*70)
print("scANVI-Existing Model: Load or Train")
print("-"*70)

if scanvi_existing_exists and LOAD_SCANVI_EXISTING_IF_EXISTS:
    print("Loading existing scANVI-Existing model...")
    try:
        scanvi_existing_model = scvi.model.SCANVI.load(scanvi_existing_model_dir, adata=adata_scanvi_existing)
        print("âœ“ Model loaded")
    except Exception as e:
        print(f"âš ï¸  Load failed: {e}")
        scanvi_existing_exists = False

if not scanvi_existing_exists or not LOAD_SCANVI_EXISTING_IF_EXISTS:
    print("Training new scANVI-Existing model...")
    
    # PATCH 2C: Setup with covariates (ä¸Ž scVI ä¸€è‡´)
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_existing,
        layer=None,
        batch_key=BATCH_KEY,
        labels_key='scanvi_label_existing_cleaned',
        unlabeled_category="Unknown",
        categorical_covariate_keys=[TISSUE_KEY],  # â­ FIXED
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES if k not in [BATCH_KEY, TISSUE_KEY]]  # â­ FIXED
    )
    
    scanvi_existing_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category="Unknown",
        adata=adata_scanvi_existing, labels_key='scanvi_label_existing_cleaned'
    )
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'batch_size': SCANVI_BATCH_SIZE,
        'early_stopping': True,
        'early_stopping_patience': SCANVI_EARLY_STOPPING_PATIENCE,
        'n_samples_per_label': SCANVI_N_SAMPLES_PER_LABEL,
        'plan_kwargs': {
            'lr': SCANVI_LEARNING_RATE,
            'weight_decay': 0.0
        },
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scanvi_existing_model.train(**train_kwargs)
    print("âœ“ Training complete")
    
    scanvi_existing_model.save(scanvi_existing_model_dir, overwrite=True)
    print("âœ“ Model saved")

# Generate results
print("\nGenerating scANVI-Existing results...")

latent_scanvi_existing_df = pd.DataFrame(
    scanvi_existing_model.get_latent_representation(),
    index=adata_scanvi_existing.obs_names
)
latent_existing_aligned = latent_scanvi_existing_df.reindex(adata.obs_names)
if latent_existing_aligned.isna().any().any():
    raise RuntimeError("scANVI-Existing latent reindex NaNs")
adata.obsm['X_scanvi_existing'] = latent_existing_aligned.to_numpy().astype(np.float32)

predictions_existing_series = pd.Series(
    scanvi_existing_model.predict(),
    index=adata_scanvi_existing.obs_names
)
predictions_existing_aligned = predictions_existing_series.reindex(adata.obs_names)
if predictions_existing_aligned.isna().any():
    raise RuntimeError("scANVI-Existing predictions reindex NaNs")
adata.obs['scanvi_predictions_existing'] = predictions_existing_aligned.astype(str).values

probs_existing = np.asarray(scanvi_existing_model.predict(soft=True))
probs_existing_df = pd.DataFrame(probs_existing, index=adata_scanvi_existing.obs_names)
probs_existing_aligned = probs_existing_df.reindex(adata.obs_names)
if probs_existing_aligned.isna().any().any():
    raise RuntimeError("scANVI-Existing probs reindex NaNs")
adata.obsm['scanvi_probabilities_existing'] = probs_existing_aligned.to_numpy().astype(np.float32)

del adata_scanvi_existing
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scanvi_existing)...")
sc.pp.neighbors(adata, use_rep='X_scanvi_existing', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_existing'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_existing',
             neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print(f"\nâœ“ scANVI-Existing complete")
print(f"  Cell types: {adata.obs['scanvi_predictions_existing'].nunique()}")
print(f"  Top 5:")
for ct, count in adata.obs['scanvi_predictions_existing'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")


# ==============================================================================
# Step 4B: scANVI-CellTypist (PATCH 2B/2C: covariates fix)
# ==============================================================================

print("\n" + "="*70)
print("Step 4B: scANVI-CellTypist (Using CellTypist predictions)")
print("="*70)

print(f"\nâ­ Training scANVI using CellTypist predictions")
print(f"   With confidence filtering (threshold: {LOW_CONFIDENCE_THRESHOLD})")
print(f"   With rare type filtering (threshold: {RARE_TYPE_THRESHOLD} cells)")
print(f"   With n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")
print(f"   With covariates (FIXED)")

# Create semi-supervised labels
print("\nCreating semi-supervised labels...")
labels_celltypist = adata.obs['celltypist_majority_voting'].astype(str).copy()
conf = adata.obs['celltypist_conf_score'].astype(float)

n_original_celltypist = len(labels_celltypist)  # â­ P0-4 FIX: unique variable name
n_low_conf = (conf < LOW_CONFIDENCE_THRESHOLD).sum()
labels_celltypist[conf < LOW_CONFIDENCE_THRESHOLD] = "Unknown"

value_counts = labels_celltypist.value_counts()
rare_types = value_counts[value_counts < RARE_TYPE_THRESHOLD].index
rare_types = rare_types[rare_types != "Unknown"]
n_rare = labels_celltypist.isin(rare_types).sum()
labels_celltypist[labels_celltypist.isin(rare_types)] = "Unknown"

n_unknown_celltypist = (labels_celltypist == "Unknown").sum()  # â­ P0-4 FIX: unique variable name
n_labeled_celltypist = n_original_celltypist - n_unknown_celltypist

print(f"  Total: {n_original_celltypist:,}")
print(f"  Low confidence â†’ Unknown: {n_low_conf:,} ({n_low_conf/n_original_celltypist*100:.1f}%)")
print(f"  Rare types â†’ Unknown: {n_rare:,} ({n_rare/n_original_celltypist*100:.1f}%)")
print(f"  Final labeled: {n_labeled_celltypist:,} ({n_labeled_celltypist/n_original_celltypist*100:.1f}%)")
print(f"  Final Unknown: {n_unknown_celltypist:,} ({n_unknown_celltypist/n_original_celltypist*100:.1f}%)")

adata.obs['scanvi_label_celltypist'] = pd.Categorical(labels_celltypist)

scanvi_celltypist_model_dir = output_dir / "models" / "scanvi_celltypist_model"
scanvi_celltypist_exists = scanvi_celltypist_model_dir.exists() and (scanvi_celltypist_model_dir / "model.pt").exists()

# Prepare data (PATCH 2B: Include covariates)
print("\nPreparing scANVI-CellTypist data...")

if USE_HVG_FOR_SCVI:
    # Check if hvg_indices is available (same as Step 4A)
    if hvg_indices is None:
        if 'highly_variable' in adata.var.columns:
            hvg_mask = adata.var["highly_variable"].values
            hvg_indices = np.flatnonzero(hvg_mask)
        elif hvg_file.exists():
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_saved))
        else:
            raise RuntimeError("USE_HVG_FOR_SCVI=True but hvg_indices is None and cannot be recovered")
    
    adata_scanvi_celltypist = sc.AnnData(
        X=adata.layers["counts"][:, hvg_indices].copy(),
        obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_celltypist']].copy(),  # â­ FIXED
        var=adata.var.iloc[hvg_indices].copy()
    )
    adata_scanvi_celltypist.var_names = adata.var_names[hvg_indices]
else:
    adata_scanvi_celltypist = sc.AnnData(
        X=adata.layers["counts"].copy(),
        obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_celltypist']].copy(),  # â­ FIXED
        var=adata.var.copy()
    )

print(f"âœ“ scANVI-CellTypist AnnData: {adata_scanvi_celltypist.n_obs:,} Ã— {adata_scanvi_celltypist.n_vars:,}")

# Ensure 'Unknown' exists
labels_cat = adata_scanvi_celltypist.obs['scanvi_label_celltypist'].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi_celltypist.obs['scanvi_label_celltypist'] = labels_cat

# Load or train
print("\n" + "-"*70)
print("scANVI-CellTypist Model: Load or Train")
print("-"*70)

if scanvi_celltypist_exists and LOAD_SCANVI_CELLTYPIST_IF_EXISTS:
    print("Loading existing scANVI-CellTypist model...")
    try:
        scanvi_celltypist_model = scvi.model.SCANVI.load(scanvi_celltypist_model_dir, adata=adata_scanvi_celltypist)
        print("âœ“ Model loaded")
    except Exception as e:
        print(f"âš ï¸  Load failed: {e}")
        scanvi_celltypist_exists = False

if not scanvi_celltypist_exists or not LOAD_SCANVI_CELLTYPIST_IF_EXISTS:
    print("Training new scANVI-CellTypist model...")
    
    # PATCH 2C: Setup with covariates (ä¸Ž scVI ä¸€è‡´)
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_celltypist,
        layer=None,
        batch_key=BATCH_KEY,
        labels_key='scanvi_label_celltypist',
        unlabeled_category="Unknown",
        categorical_covariate_keys=[TISSUE_KEY],  # â­ FIXED
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES if k not in [BATCH_KEY, TISSUE_KEY]]  # â­ FIXED
    )
    
    scanvi_celltypist_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, unlabeled_category="Unknown",
        adata=adata_scanvi_celltypist, labels_key='scanvi_label_celltypist'
    )
    
    train_kwargs = {
        'max_epochs': SCANVI_MAX_EPOCHS,
        'batch_size': SCANVI_BATCH_SIZE,
        'early_stopping': True,
        'early_stopping_patience': SCANVI_EARLY_STOPPING_PATIENCE,
        'n_samples_per_label': SCANVI_N_SAMPLES_PER_LABEL,
        'plan_kwargs': {
            'lr': SCANVI_LEARNING_RATE,
            'weight_decay': 0.0
        },
        'enable_progress_bar': True,
    }
    if gpu_available:
        train_kwargs['accelerator'] = 'gpu'
        train_kwargs['devices'] = 'auto'
    
    scanvi_celltypist_model.train(**train_kwargs)
    print("âœ“ Training complete")
    
    scanvi_celltypist_model.save(scanvi_celltypist_model_dir, overwrite=True)
    print("âœ“ Model saved")

# Generate results
print("\nGenerating scANVI-CellTypist results...")

latent_scanvi_celltypist_df = pd.DataFrame(
    scanvi_celltypist_model.get_latent_representation(),
    index=adata_scanvi_celltypist.obs_names
)
latent_celltypist_aligned = latent_scanvi_celltypist_df.reindex(adata.obs_names)
if latent_celltypist_aligned.isna().any().any():
    raise RuntimeError("scANVI-CellTypist latent reindex NaNs")
adata.obsm['X_scanvi_celltypist'] = latent_celltypist_aligned.to_numpy().astype(np.float32)

predictions_celltypist_series = pd.Series(
    scanvi_celltypist_model.predict(),
    index=adata_scanvi_celltypist.obs_names
)
predictions_celltypist_aligned = predictions_celltypist_series.reindex(adata.obs_names)
if predictions_celltypist_aligned.isna().any():
    raise RuntimeError("scANVI-CellTypist predictions reindex NaNs")
adata.obs['scanvi_predictions_celltypist'] = predictions_celltypist_aligned.astype(str).values

probs_celltypist = np.asarray(scanvi_celltypist_model.predict(soft=True))
probs_celltypist_df = pd.DataFrame(probs_celltypist, index=adata_scanvi_celltypist.obs_names)
probs_celltypist_aligned = probs_celltypist_df.reindex(adata.obs_names)
if probs_celltypist_aligned.isna().any().any():
    raise RuntimeError("scANVI-CellTypist probs reindex NaNs")
adata.obsm['scanvi_probabilities_celltypist'] = probs_celltypist_aligned.to_numpy().astype(np.float32)

del adata_scanvi_celltypist
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors (neighbors_scanvi_celltypist)...")
sc.pp.neighbors(adata, use_rep='X_scanvi_celltypist', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_celltypist'] = adata.obsm['X_umap'].copy()

sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_celltypist',
             neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print(f"\nâœ“ scANVI-CellTypist complete")
print(f"  Cell types: {adata.obs['scanvi_predictions_celltypist'].nunique()}")
print(f"  Top 5:")
for ct, count in adata.obs['scanvi_predictions_celltypist'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")


# ==============================================================================
# Step 5: Comparison & Final Export
# ==============================================================================

print("\n" + "="*70)
print("Step 5: Comparison & Final Export")
print("="*70)

# Comparison visualization
print("\nGenerating comparison plots...")

fig, axes = plt.subplots(2, 4, figsize=(24, 12))

# Row 1: scANVI-Existing
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_existing'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, title='scANVI-Existing: Batch')
sc.pl.umap(adata, color='scanvi_label_existing_cleaned', ax=axes[0, 1], show=False, title='scANVI-Existing: Input (cleaned)')
sc.pl.umap(adata, color='scanvi_predictions_existing', ax=axes[0, 2], show=False, title='scANVI-Existing: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_existing', ax=axes[0, 3], show=False, title='scANVI-Existing: Leiden')

# Row 2: scANVI-CellTypist
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0], show=False, title='scANVI-CellTypist: Batch')
sc.pl.umap(adata, color='scanvi_label_celltypist', ax=axes[1, 1], show=False, title='scANVI-CellTypist: Input')
sc.pl.umap(adata, color='scanvi_predictions_celltypist', ax=axes[1, 2], show=False, title='scANVI-CellTypist: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_celltypist', ax=axes[1, 3], show=False, title='scANVI-CellTypist: Leiden')

plt.tight_layout()
plt.savefig(output_dir / "figures" / f"dual_scanvi_comparison.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.close()

print("âœ“ Comparison plots saved")

# Set default UMAP to CellTypist-based
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()

# Save final dataset
output_file = output_dir / "adata_allcells_dual_scanvi_final_v2.3.1_HOTFIX.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"âœ“ Saved ({file_size:.2f} GB)")

# Summary
print("\n" + "="*70)
print("ANALYSIS SUMMARY v2.3.1 HOTFIX")
print("="*70)

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("All-Cells scVI/scANVI Analysis v2.3.1 HOTFIX")
summary_lines.append("="*70)

summary_lines.append("\n[ Critical Hotfixes ]")
summary_lines.append("  âœ“ P0-2: Corrected neighborhood purity calculation")
summary_lines.append("      - Now measures agreement with cell's OWN label")
summary_lines.append("      - Excludes self-loops in neighbor graph")
summary_lines.append("      - Uses efficient CSR matrix indexing")
summary_lines.append("  âœ“ P0-4: Fixed summary variable scoping")
summary_lines.append("      - Separate variables for existing vs celltypist")
summary_lines.append("      - Accurate Unknown counts in final report")

summary_lines.append("\n[ scANVI-Existing (cell_type metadata) ]")
summary_lines.append(f"  Input: cleaned labels (CORRECTED purity + quality)")
summary_lines.append(f"  Total cells: {n_original_existing:,}")
summary_lines.append(f"  Final labeled: {n_labeled_existing:,} ({n_labeled_existing/n_original_existing*100:.1f}%)")
summary_lines.append(f"  Final Unknown: {n_unknown_existing:,} ({n_unknown_existing/n_original_existing*100:.1f}%)")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_existing'].nunique()}")

summary_lines.append("\n[ scANVI-CellTypist (CellTypist predictions) ]")
summary_lines.append(f"  Input: CellTypist predictions")
summary_lines.append(f"  Total cells: {n_original_celltypist:,}")
summary_lines.append(f"  Final labeled: {n_labeled_celltypist:,} ({n_labeled_celltypist/n_original_celltypist*100:.1f}%)")
summary_lines.append(f"  Final Unknown: {n_unknown_celltypist:,} ({n_unknown_celltypist/n_original_celltypist*100:.1f}%)")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_celltypist'].nunique()}")

summary_lines.append("\n[ Key Outputs ]")
summary_lines.append("  scANVI-Existing:")
summary_lines.append(f"    - scanvi_predictions_existing")
summary_lines.append(f"    - scanvi_label_existing_cleaned (with Unknown)")
summary_lines.append(f"    - label_purity_existing (CORRECTED)")
summary_lines.append("  scANVI-CellTypist:")
summary_lines.append(f"    - scanvi_predictions_celltypist")
summary_lines.append(f"    - scanvi_label_celltypist")

summary_lines.append("\n[ Purity Metric Explanation (HOTFIX) ]")
summary_lines.append("  Previous (WRONG): fraction of neighbors in dominant group")
summary_lines.append("    â†’ High purity even if cell mislabeled (if neighbors pure)")
summary_lines.append("  Current (CORRECT): fraction of neighbors with SAME label as cell")
summary_lines.append("    â†’ Low purity = cell label disagrees with neighbors")
summary_lines.append("    â†’ Correctly identifies mislabeled/ambiguous cells")

summary_lines.append("\n" + "="*70)
summary_lines.append("Analysis complete!")
summary_lines.append("="*70)

summary_text = '\n'.join(summary_lines)
print(summary_text)

summary_file = output_dir / "analysis_summary_v2.3.1_HOTFIX.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

print(f"\nâœ“ Summary saved to: {summary_file}")

print("\n" + "="*70)
print("ðŸŽ‰ HOTFIX COMPLETE (v2.3.1 HOTFIX)")
print("="*70)

print(f"\nðŸ“Š Hotfix Summary:")
print(f"  â€¢ P0-2 FIXED: Neighborhood purity now correctly identifies mislabeled cells")
print(f"  â€¢ P0-4 FIXED: Summary statistics now accurate for both pipelines")
print(f"  â€¢ All v2.3.1 features retained (covariates, MT%, normalized scoring)")

print("\nðŸ“ Expected Impact:")
print("  â€¢ More accurate Unknown labeling for existing metadata")
print("  â€¢ Better removal of mislabeled/ambiguous cells")
print("  â€¢ Improved scANVI-Existing model quality")
print("  â€¢ Clearer separation of cell types in UMAP")

print("\n" + "="*70)