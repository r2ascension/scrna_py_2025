#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
==============================================================================
All-Cells scVI/scANVI Dual Pipeline v2.5.2 PRODUCTION
==============================================================================

Author: Clinical-Bioinformatics Team
Date: 2025-01-21
Version: 2.5.2 PRODUCTION

CRITICAL IMPROVEMENTS in v2.5.2:
=================================
1. Fixed HVG workflow (correct execution order) [v2.5]
2. Simplified model control (clear retrain flags) [v2.5]
3. Integrated best practices from QUICK_REFERENCE_MEMORY v2.11-12 [v2.5]
4. Complete English documentation (UTF-8 compatible) [v2.5]
5. Safety mechanisms (HVG check, error handling, status reporting) [v2.5]
6. Robust gene symbol detection (12 column name fallbacks) [v2.5.1 HOTFIX]
7. Full gene preservation in .raw (0 extra memory cost) [v2.5.2 NEW]
8. Updated parameter documentation for scArches compatibility [v2.5.2 NEW]

Pipeline Overview:
==================
Input: adata_bbknn_annotated_corrected.h5ad (assumes pre-filtered genes)
  |
Step 0: Whitelist Filtering (optional)
  |
Step 1: Data Preparation + Covariates
  - MT% recalculation (handle zeros)
  - Stress score (normalized)
  - Cell cycle scores (S, G2M)
  |
Step 2: scVI Integration
  - HVG selection (batch-aware with fallback)
  - Force biological markers into HVG
  - Save full genes to .raw (0 memory cost)
  - Train/load scVI model (scArches-ready)
  |
Step 3: CellTypist Annotation
  - Auto-annotation with majority voting
  - Robust gene symbol detection
  |
Step 4A: scANVI-Existing
  - Use cleaned existing cell_type labels
  - Unknown filtering (purity + quality)
  - From scVI model (transfer learning)
  |
Step 4B: scANVI-CellTypist
  - Use CellTypist predictions
  - From scVI model (transfer learning)
  |
Step 5: Comparison & Export
  - Generate comparison plots
  - Save final h5ad with all results

Key Features:
=============
- Memory optimized (HVG-based training, <48GB RAM)
- Full gene access (.raw preservation, 0 extra memory)
- GPU accelerated (scVI/scANVI training)
- Reproducible (fixed random seeds, saved HVG list)
- Safety checks (HVG consistency, graceful error handling)
- Robust gene detection (12 column name fallbacks)
- Status reporting (clear progress and model loading plan)

Model Loading Strategy:
=======================
# Pattern 1: Full Speed (downstream analysis only)
RETRAIN_SCVI = False
RETRAIN_SCANVI_EXISTING = False
RETRAIN_SCANVI_CELLTYPIST = False
# Use: Visualization, export, report generation
# Saves: 70-105 min -> 70 sec (100x faster)

# Pattern 2: Fresh Start (default for new data)
RETRAIN_SCVI = True
RETRAIN_SCANVI_EXISTING = True
RETRAIN_SCANVI_CELLTYPIST = True
# Use: First run, changed HVG, updated packages
# Time: Full training (~101 min for 252k cells)

# Pattern 3: Selective Retrain (label refinement)
RETRAIN_SCVI = False         # Keep integration
RETRAIN_SCANVI_EXISTING = True   # Retrain with new labels
RETRAIN_SCANVI_CELLTYPIST = True # Retrain with new strategy
# Use: Testing new annotation approaches
# Saves: ~42 min (scVI time)

Requirements:
=============
pip install "scvi-tools>=1.1.4" "lightning>=2.2,<3" scanpy celltypist mygene h5py

Note on Gene Filtering:
========================
This pipeline assumes input data has been pre-filtered for gene quality.
If your data includes MT genes, ribosomal genes, or unannotated genes,
please filter them upstream before running this pipeline.
"""

# ==============================================================================
# IMPORTS
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

# Version check
from packaging import version
if version.parse(scvi.__version__) < version.parse("1.1.4"):
    raise RuntimeError(
        f"Warning: scvi-tools {scvi.__version__} detected. "
        "Upgrade to >=1.1.4 recommended for critical scANVI fixes."
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
# CONFIGURATION
# ==============================================================================

# ===== Input/Output Paths =====
INPUT_H5AD = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected_filtered.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0115/allcells_scvi_analysis_v2_3_1_HOTFIX"
CELLTYPIST_MODEL_PATH = "/home/h2048/data/source/reference/celltypist_models/Cells_Lung_Airway.pkl"

# ===== Whitelist Filtering (Optional) =====
DOWNSTREAM_H5AD_DIR = "/home/h2048/data/core20260115"
DOWNSTREAM_H5AD_GLOB = "**/*.h5ad"
CELLID_SET_MODE = "union"  # "union" or "intersection"
ENABLE_WHITELIST_FILTERING = False

# ===== Data Keys =====
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
EXISTING_CELLTYPE_KEY = "cell_type"

# ===== Required Covariates =====
# These will be used as covariates in scVI/scANVI models
REQUIRED_COVARIATES = [
    BATCH_KEY,           # Categorical: sample batch
    TISSUE_KEY,          # Categorical: tissue type
    'pct_counts_mt',     # Continuous: mitochondrial percentage
    'stress_score',      # Continuous: stress signature score
    'S_score',           # Continuous: S phase score
    'G2M_score'          # Continuous: G2M phase score
]

# ===== Model Control Flags =====
# CRITICAL: Set these based on your use case
# See "Model Loading Strategy" in header for guidance

RETRAIN_SCVI = False              # Force retrain scVI model
RETRAIN_SCANVI_EXISTING = False   # Force retrain scANVI-Existing
RETRAIN_SCANVI_CELLTYPIST = False # Force retrain scANVI-CellTypist

# ===== HVG Settings =====
USE_HVG_FOR_SCVI = True  # CRITICAL: Always True for memory optimization
N_HVG_SCVI = 6000        # Number of highly variable genes
HVG_FLAVOR = "seurat_v3" # HVG selection method

# ===== Force Biological Markers into HVG =====
FORCE_MARKERS_IN_HVG = True
FORCE_MARKERS_CASE_INSENSITIVE = True

# Lineage anchor markers (CRITICAL: preserve cell type structure)
FORCED_MARKERS_LINEAGE = [
    # Epithelial markers
    "FXYD3","EPCAM","ELF3","IGFBP2","SERPINF1","TSPAN1","SCGB1A1",
    "AGER","SFTPC","FOXJ1","KRT5","MUC5B","KRT8",
    # Immune general / Lymphoid
    "CD53","PTPRC","CORO1A","ISG20","CCL5",
    # B cell markers
    "MS4A1","TNFRSF17","CD19","CD79A","SDC1",
    # T cell markers
    "CD40LG","TNFRSF25","CD28","CD4","CD3D","CD3E","CD2","TRBC2",
    "CD8A","CD8B","TRGC2",
    # Myeloid markers
    "FCER1G","C1ORF162","CLEC7A","CD1C","CD86","CD14","XCR1","HLA-DRA",
    # Stromal / Fibroblast / SMC
    "COL1A2","DCN","MFAP4","LUM","COL6A3","CFD","COL1A1","PDGFRA",
    "MXRA8","NBL1","VCAN","LEPR","MYH11","TINAGL1","PLN","DES",
    "ACTA2","CNN1","TAGLN",
    # Endothelial markers
    "CLDN5","ECSCR","CLEC14A","VWF","PECAM1","ACKR1","PTPRB","PDE2A",
    "PLAT","GJA5","SPARCL1","AQP1","RNASE1","MMRN1","CCL21","TFF3",
]

# State markers (optional, disabled by default to avoid cell cycle dominance)
FORCED_MARKERS_STATE = [
    "MKI67","TOP2A","TK1","CENPW"  # Proliferation markers
]

INCLUDE_STATE_MARKERS = False  # Set True to include proliferation axis
FORCED_MARKERS = FORCED_MARKERS_LINEAGE + (FORCED_MARKERS_STATE if INCLUDE_STATE_MARKERS else [])

# ===== Stress Signature Genes =====
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
# Local gene lists (no network dependency)
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

# ===== scVI Hyperparameters =====
# CRITICAL: n_latent=150, n_layers=2 for scArches compatibility
# Reference: QUICK_REFERENCE_MEMORY v2.11 scArches best practices
#
# Why these specific values:
# - n_latent=150: Increased from default 10-50 for better stability
#                 Captures more biological variance in large (400k+) datasets
#                 Optimized for scArches query mapping (transfer learning)
# - n_layers=2:   Reduced from default 4 for scArches compatibility
#                 Simpler architectures transfer better to new datasets
#                 Prevents overfitting in batch correction
# - batch_size=256: Reduced from 512 for memory efficiency
#                   Works well with large datasets (400k+ cells)
SCVI_N_LATENT = 150          # Latent dimensions (scArches-optimized, NOT default)
SCVI_N_LAYERS = 2            # Network depth (scArches-compatible, reduced from 4)
SCVI_DROPOUT_RATE = 0.2      # Dropout rate (increased from 0.1 for stability)
SCVI_MAX_EPOCHS = 400        # Training epochs (sufficient for convergence)
SCVI_BATCH_SIZE = 256        # Batch size (reduced from 512 for memory)
SCVI_LEARNING_RATE = 1e-3    # Learning rate (standard for scVI)
SCVI_EARLY_STOPPING = True   # Enable early stopping (prevents overtraining)
SCVI_EARLY_STOPPING_PATIENCE = 45  # Patience (increased for stability)
SCVI_ENCODE_COVARIATES = True      # Use covariates in latent space
SCVI_USE_LAYER_NORM = "both"       # Layer norm (encoder+decoder)
SCVI_USE_BATCH_NORM = "none"       # Batch norm (disabled for stability)

# ===== scANVI Hyperparameters =====
# CRITICAL: n_samples_per_label=2000 handles class imbalance
# - High value ensures rare cell types get sufficient training samples
# - Prevents model from ignoring minority classes
# - Essential for immune cell subtypes (e.g., rare T cell subsets)
SCANVI_MAX_EPOCHS = 200              # Training epochs (sufficient for semi-supervised)
SCANVI_BATCH_SIZE = 256              # Same as scVI for consistency
SCANVI_LEARNING_RATE = 5e-4          # Lower than scVI (fine-tuning)
SCANVI_EARLY_STOPPING_PATIENCE = 30  # Patience (reduced, faster convergence)
SCANVI_N_SAMPLES_PER_LABEL = 2000    # Class balancing (CRITICAL for rare types)

# ===== Unknown Labeling Thresholds =====
# For scANVI-Existing: clean existing labels before training
EXISTING_LABEL_PURITY_THRESHOLD = 0.5      # Min neighbor agreement
EXISTING_LABEL_MT_THRESHOLD = 15.0         # Max MT% (%)
EXISTING_LABEL_STRESS_PERCENTILE = 95      # Max stress percentile

# ===== UMAP Parameters =====
# Optimized for large datasets (400k+ cells)
# - min_dist=0.5: Increased from default 0.3 for less crowded visualization
# - Higher min_dist = more global structure, less local clumping
UMAP_N_NEIGHBORS = 30     # Number of neighbors (standard)
UMAP_MIN_DIST = 0.5       # Minimum distance (increased for clarity)
UMAP_SPREAD = 1.0         # Spread (standard)

# ===== Visualization =====
# Production-quality settings for publications
FIGURE_FORMAT = 'pdf'  # Vector format for publications
DPI = 300              # Print quality resolution

# ===== Random Seed =====
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if gpu_available:
    torch.cuda.manual_seed(RANDOM_SEED)


# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def create_output_dirs(base_dir):
    """
    Create output directory structure.
    
    Parameters:
    -----------
    base_dir : str
        Base output directory path
        
    Returns:
    --------
    dirs : dict
        Dictionary of created directory paths
    """
    base_path = Path(base_dir)
    dirs = {
        'base': base_path,
        'figures': base_path / 'figures',
        'models': base_path / 'models',
        'logs': base_path / 'logs'
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def detect_gene_symbol_column(adata, verbose=True):
    """
    Detect the gene symbol column in adata.var with automatic fallback.
    
    CRITICAL: CellTypist and marker forcing require gene symbols to match.
    Different preprocessing pipelines may store gene symbols in different columns.
    
    Strategy:
    1. Try common symbol column names in order of preference
    2. Fall back to adata.var_names if no symbol column found
    3. Warn user about which column is being used
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    verbose : bool
        Whether to print detection result
        
    Returns:
    --------
    symbol_col : str or None
        Name of the gene symbol column, or None to use var_names
    gene_symbols : pd.Series
        Series of gene symbols
        
    Example:
    --------
    >>> symbol_col, gene_symbols = detect_gene_symbol_column(adata)
    >>> common_genes = gene_symbols.isin(model_features)
    """
    # Common gene symbol column names (in order of preference)
    CANDIDATE_COLUMNS = [
        'symbol_base',      # Our standard column name
        'gene_symbol',      # Common in many pipelines
        'gene_symbols',     # Plural variant
        'Symbol',           # Capitalized version
        'symbol',           # Lowercase version
        'feature_name',     # 10x Cell Ranger
        'gene_name',        # Alternative naming
        'gene_names',       # Plural variant
        'gene',             # Short version
        'genes',            # Plural short version
        'gene_id',          # Sometimes gene IDs are symbols
        'gene_ids',         # Plural variant
    ]
    
    # Check each candidate
    for col in CANDIDATE_COLUMNS:
        if col in adata.var.columns:
            gene_symbols = adata.var[col].astype(str)
            
            # Verify it's not mostly empty or NA
            non_empty = gene_symbols.notna() & (gene_symbols != '') & (gene_symbols != 'nan')
            if non_empty.sum() >= len(gene_symbols) * 0.9:  # At least 90% valid
                if verbose:
                    print(f"  Using gene symbols from: adata.var['{col}']")
                    print(f"  Valid symbols: {non_empty.sum():,}/{len(gene_symbols):,}")
                return col, gene_symbols
    
    # Fallback: use var_names
    if verbose:
        print("  Warning: No gene symbol column found in adata.var")
        print("  Falling back to: adata.var_names (gene indices)")
        print("  Note: This may reduce CellTypist matching accuracy")
    
    return None, pd.Series(adata.var_names.values, index=adata.var_names)


def force_include_markers_in_hvg(adata, n_top_genes, forced_markers, 
                                  symbol_col=None,
                                  case_insensitive=True, 
                                  hvg_col="highly_variable"):
    """
    Force specified marker genes into HVG selection.
    
    CRITICAL: Ensures biological markers are preserved in HVG even if
    they have low variance (e.g., cell-type-specific markers expressed
    only in rare populations).
    
    Strategy:
    1. Auto-detect gene symbol column if not provided
    2. Identify forced markers present in data
    3. Check if already in HVG
    4. For missing ones: swap out lowest-variance HVG genes
    5. Maintain target HVG count (n_top_genes)
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    n_top_genes : int
        Target number of HVG genes
    forced_markers : list
        List of marker gene symbols to force include
    symbol_col : str or None
        Column in adata.var containing gene symbols.
        If None, will auto-detect using detect_gene_symbol_column()
    case_insensitive : bool
        Whether to perform case-insensitive matching
    hvg_col : str
        Column in adata.var containing HVG boolean flags
        
    Returns:
    --------
    report : dict
        Summary of forced marker inclusion
    """
    if hvg_col not in adata.var.columns:
        raise ValueError(f"HVG column '{hvg_col}' not found in adata.var")
    
    # Auto-detect symbol column if not provided
    if symbol_col is None:
        print("  Auto-detecting gene symbol column...")
        symbol_col, gene_symbols = detect_gene_symbol_column(adata, verbose=False)
        if symbol_col is None:
            # Use var_names
            var_symbols = adata.var_names.values
        else:
            var_symbols = adata.var[symbol_col].values
    else:
        if symbol_col not in adata.var.columns:
            raise ValueError(f"Symbol column '{symbol_col}' not found in adata.var")
        var_symbols = adata.var[symbol_col].values
    
    # Normalize marker names
    if case_insensitive:
        forced_set = {m.upper() for m in forced_markers}
        if isinstance(var_symbols[0], str):
            var_symbols = np.array([str(s).upper() for s in var_symbols])
    else:
        forced_set = set(forced_markers)
        var_symbols = np.array([str(s) for s in var_symbols])
    
    # Find forced markers in data
    forced_mask = np.array([s in forced_set for s in var_symbols])
    forced_indices = np.flatnonzero(forced_mask)
    
    if len(forced_indices) == 0:
        return {
            'n_requested': len(forced_markers),
            'n_found_in_data': 0,
            'n_already_in_hvg': 0,
            'n_newly_added': 0,
            'n_hvg_final': adata.var[hvg_col].sum()
        }
    
    # Check current HVG status
    hvg_mask = adata.var[hvg_col].values.copy()
    forced_in_hvg = hvg_mask[forced_indices]
    n_already_in = forced_in_hvg.sum()
    n_to_add = len(forced_indices) - n_already_in
    
    if n_to_add == 0:
        return {
            'n_requested': len(forced_markers),
            'n_found_in_data': len(forced_indices),
            'n_already_in_hvg': n_already_in,
            'n_newly_added': 0,
            'n_hvg_final': hvg_mask.sum()
        }
    
    # Get variance info (if available)
    var_col = None
    for vc in ['variances_norm', 'variances', 'dispersions_norm', 'dispersions']:
        if vc in adata.var.columns:
            var_col = vc
            break
    
    if var_col is None:
        # No variance info: just add to HVG
        to_add_mask = forced_mask & ~hvg_mask
        hvg_mask[to_add_mask] = True
    else:
        # Smart swap: remove lowest-variance HVG, add forced markers
        variances = adata.var[var_col].values.copy()
        
        # Identify HVG to remove (lowest variance, non-forced)
        current_hvg_indices = np.flatnonzero(hvg_mask)
        current_hvg_non_forced = np.array([i for i in current_hvg_indices 
                                            if not forced_mask[i]])
        
        if len(current_hvg_non_forced) >= n_to_add:
            # Sort by variance and remove lowest
            hvg_variances = variances[current_hvg_non_forced]
            remove_indices = current_hvg_non_forced[np.argsort(hvg_variances)[:n_to_add]]
            hvg_mask[remove_indices] = False
        
        # Add all forced markers not in HVG
        to_add_mask = forced_mask & ~hvg_mask
        hvg_mask[to_add_mask] = True
    
    # Update HVG column
    adata.var[hvg_col] = hvg_mask
    
    # Ensure exactly n_top_genes (if possible)
    n_hvg_now = hvg_mask.sum()
    if n_hvg_now > n_top_genes:
        # Too many: remove lowest-variance non-forced
        current_hvg_indices = np.flatnonzero(hvg_mask)
        non_forced_hvg = np.array([i for i in current_hvg_indices 
                                   if not forced_mask[i]])
        
        if var_col is not None and len(non_forced_hvg) > 0:
            n_to_remove = n_hvg_now - n_top_genes
            hvg_variances = variances[non_forced_hvg]
            remove_indices = non_forced_hvg[np.argsort(hvg_variances)[:n_to_remove]]
            hvg_mask[remove_indices] = False
            adata.var[hvg_col] = hvg_mask
    
    return {
        'n_requested': len(forced_markers),
        'n_found_in_data': len(forced_indices),
        'n_already_in_hvg': n_already_in,
        'n_newly_added': n_to_add,
        'n_hvg_final': adata.var[hvg_col].sum()
    }


def calculate_label_agreement_purity(adata, label_key, neighbors_key='neighbors'):
    """
    Calculate label purity as fraction of neighbors sharing cell's OWN label.
    
    CRITICAL CORRECTION (P0-2 fix from v2.3.1):
    - Measures agreement with cell's own label (not dominant neighbor label)
    - Excludes self-loops from neighbor graph
    - Uses efficient CSR matrix indexing
    
    LOW purity indicates:
    - Cell label disagrees with its neighbors
    - Potentially mislabeled cell
    - Ambiguous cell state (transition, doublet)
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with neighborhood graph
    label_key : str
        Column in adata.obs containing cell labels
    neighbors_key : str
        Key for neighborhood graph in adata.obsp
        
    Returns:
    --------
    purity : np.ndarray
        Array of purity scores (0-1) for each cell
        
    Example:
    --------
    >>> purity = calculate_label_agreement_purity(adata, 'cell_type', 'neighbors_scvi')
    >>> low_purity = purity < 0.5
    >>> print(f"Low purity cells: {low_purity.sum()}")
    """
    # Get labels
    labels = adata.obs[label_key].values
    n_cells = len(labels)
    
    # Get neighbor connectivities
    conn_key = f'{neighbors_key}_connectivities'
    if conn_key not in adata.obsp:
        raise ValueError(f"Connectivity matrix '{conn_key}' not found in adata.obsp")
    
    connectivities = adata.obsp[conn_key]
    if not isinstance(connectivities, csr_matrix):
        connectivities = csr_matrix(connectivities)
    
    # Calculate purity for each cell
    purity = np.zeros(n_cells, dtype=np.float32)
    
    for i in range(n_cells):
        cell_label = labels[i]
        
        # Get neighbors (exclude self)
        start = connectivities.indptr[i]
        end = connectivities.indptr[i + 1]
        neighbor_indices = connectivities.indices[start:end]
        neighbor_weights = connectivities.data[start:end]
        
        # Remove self-loops
        non_self_mask = neighbor_indices != i
        neighbor_indices = neighbor_indices[non_self_mask]
        neighbor_weights = neighbor_weights[non_self_mask]
        
        if len(neighbor_indices) == 0:
            purity[i] = 1.0  # Isolated cell: assign high purity
            continue
        
        # Get neighbor labels
        neighbor_labels = labels[neighbor_indices]
        
        # Calculate weighted purity
        same_label_mask = neighbor_labels == cell_label
        same_label_weight = neighbor_weights[same_label_mask].sum()
        total_weight = neighbor_weights.sum()
        
        purity[i] = same_label_weight / total_weight if total_weight > 0 else 1.0
    
    return purity


def print_model_status(scvi_dir, scanvi_existing_dir, scanvi_celltypist_dir,
                       retrain_scvi, retrain_existing, retrain_celltypist):
    """
    Print clear status report of existing models and loading plan.
    
    Parameters:
    -----------
    scvi_dir : Path
        scVI model directory
    scanvi_existing_dir : Path
        scANVI-Existing model directory
    scanvi_celltypist_dir : Path
        scANVI-CellTypist model directory
    retrain_scvi : bool
        Force retrain scVI flag
    retrain_existing : bool
        Force retrain scANVI-Existing flag
    retrain_celltypist : bool
        Force retrain scANVI-CellTypist flag
    """
    print("\n" + "="*70)
    print("MODEL STATUS CHECK (v2.5)")
    print("="*70)
    
    # Check existing models
    scvi_exists = scvi_dir.exists() and (scvi_dir / "model.pt").exists()
    scanvi_existing_exists = scanvi_existing_dir.exists() and (scanvi_existing_dir / "model.pt").exists()
    scanvi_celltypist_exists = scanvi_celltypist_dir.exists() and (scanvi_celltypist_dir / "model.pt").exists()
    
    print("\nExisting Models:")
    print(f"  scVI model:               {scvi_exists}")
    print(f"  scANVI-Existing model:    {scanvi_existing_exists}")
    print(f"  scANVI-CellTypist model:  {scanvi_celltypist_exists}")
    
    print("\nRetrain Flags:")
    print(f"  RETRAIN_SCVI:              {retrain_scvi}")
    print(f"  RETRAIN_SCANVI_EXISTING:   {retrain_existing}")
    print(f"  RETRAIN_SCANVI_CELLTYPIST: {retrain_celltypist}")
    
    print("\nExecution Plan:")
    
    # scVI plan
    if retrain_scvi or not scvi_exists:
        action = "TRAIN new"
        if retrain_scvi and scvi_exists:
            action = "RETRAIN (ignoring existing)"
    else:
        action = "LOAD existing"
    print(f"  scVI:               {action}")
    
    # scANVI-Existing plan
    if retrain_existing or not scanvi_existing_exists:
        action = "TRAIN new"
        if retrain_existing and scanvi_existing_exists:
            action = "RETRAIN (ignoring existing)"
    else:
        action = "LOAD existing"
    print(f"  scANVI-Existing:    {action}")
    
    # scANVI-CellTypist plan
    if retrain_celltypist or not scanvi_celltypist_exists:
        action = "TRAIN new"
        if retrain_celltypist and scanvi_celltypist_exists:
            action = "RETRAIN (ignoring existing)"
    else:
        action = "LOAD existing"
    print(f"  scANVI-CellTypist:  {action}")
    
    print("\n" + "="*70)


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

print("\n" + "="*70)
print("All-Cells scVI/scANVI Pipeline v2.5 PRODUCTION")
print("="*70)
print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Create output directories
output_dir = create_output_dirs(OUTPUT_DIR)
print(f"\nOutput directory: {output_dir['base']}")

# Print model status
scvi_model_dir = output_dir['models'] / "scvi_model"
scanvi_existing_model_dir = output_dir['models'] / "scanvi_existing_model"
scanvi_celltypist_model_dir = output_dir['models'] / "scanvi_celltypist_model"

print_model_status(
    scvi_model_dir, scanvi_existing_model_dir, scanvi_celltypist_model_dir,
    RETRAIN_SCVI, RETRAIN_SCANVI_EXISTING, RETRAIN_SCANVI_CELLTYPIST
)


# ==============================================================================
# Step 0: Whitelist Filtering (Optional)
# ==============================================================================

if ENABLE_WHITELIST_FILTERING:
    print("\n" + "="*70)
    print("Step 0: Whitelist Filtering")
    print("="*70)
    
    print(f"\nSearching for processed h5ad files...")
    print(f"  Directory: {DOWNSTREAM_H5AD_DIR}")
    print(f"  Pattern: {DOWNSTREAM_H5AD_GLOB}")
    
    downstream_dir = Path(DOWNSTREAM_H5AD_DIR)
    h5ad_files = list(downstream_dir.glob(DOWNSTREAM_H5AD_GLOB))
    print(f"  Found {len(h5ad_files)} h5ad files")
    
    if len(h5ad_files) > 0:
        print("\nExtracting cell IDs from downstream files...")
        cell_id_sets = []
        
        for h5ad_file in h5ad_files:
            try:
                temp_adata = sc.read_h5ad(h5ad_file)
                cell_ids = set(temp_adata.obs_names)
                cell_id_sets.append(cell_ids)
                print(f"  {h5ad_file.name}: {len(cell_ids):,} cells")
                del temp_adata
                gc.collect()
            except Exception as e:
                print(f"  Warning: Failed to read {h5ad_file.name}: {e}")
        
        if len(cell_id_sets) > 0:
            if CELLID_SET_MODE == "union":
                whitelist = set.union(*cell_id_sets)
                print(f"\nWhitelist created (union): {len(whitelist):,} cells")
            else:  # intersection
                whitelist = set.intersection(*cell_id_sets)
                print(f"\nWhitelist created (intersection): {len(whitelist):,} cells")
        else:
            print("\nWarning: No valid cell IDs extracted, skipping whitelist filtering")
            whitelist = None
    else:
        print("\nWarning: No downstream h5ad files found, skipping whitelist filtering")
        whitelist = None
else:
    whitelist = None
    print("\nInfo: Whitelist filtering disabled")


# ==============================================================================
# Step 1: Data Loading and Preparation
# ==============================================================================

print("\n" + "="*70)
print("Step 1: Data Loading and Preparation")
print("="*70)

print(f"\nLoading data from: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# Apply whitelist filtering
if whitelist is not None:
    print("\nApplying whitelist filtering...")
    n_before = adata.n_obs
    adata = adata[adata.obs_names.isin(whitelist)].copy()
    n_after = adata.n_obs
    n_removed = n_before - n_after
    print(f"  Removed: {n_removed:,} cells ({n_removed/n_before*100:.1f}%)")
    print(f"  Remaining: {n_after:,} cells")

# Verify required keys
print("\nVerifying data structure...")
required_keys = {
    'obs': [BATCH_KEY, EXISTING_CELLTYPE_KEY],
    'layers': ['counts']
}

for key_type, keys in required_keys.items():
    for key in keys:
        if key_type == 'obs':
            if key not in adata.obs.columns:
                raise ValueError(f"Missing obs key: {key}")
        elif key_type == 'layers':
            if key not in adata.layers:
                raise ValueError(f"Missing layer: {key}")

print("Data structure verified")

# Check if covariates need calculation
print("\nChecking covariate status...")
covariates_to_calculate = []
for cov in REQUIRED_COVARIATES:
    if cov not in adata.obs.columns:
        covariates_to_calculate.append(cov)

if len(covariates_to_calculate) > 0:
    print(f"  Missing covariates: {covariates_to_calculate}")
    print("  Will calculate...")
else:
    print("  All covariates present")

# Recalculate MT% (handle zeros)
print("\nRecalculating MT%...")
mt_genes = (adata.var_names.str.startswith('MT-') | 
            adata.var_names.str.startswith('Mt-') | 
            adata.var_names.str.startswith('mt-'))
n_mt_genes = mt_genes.sum()
print(f"  MT genes: {n_mt_genes}")

if n_mt_genes > 0:
    counts = adata.layers['counts']
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
    
    adata.obs['pct_counts_mt'] = mt_pct
    print(f"  Mean MT%: {mt_pct.mean():.2f}%")
    print(f"  Cells with zero total counts: {(~nonzero_mask).sum():,}")
else:
    adata.obs['pct_counts_mt'] = 0.0
    print("  Warning: No MT genes found, set MT%=0")

# Calculate stress score
print("\nCalculating stress score...")
stress_genes_in_data = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
n_stress = len(stress_genes_in_data)
print(f"  Stress genes: {n_stress}/{len(STRESS_SIGNATURE_GENES)}")

if n_stress >= 10:
    adata_temp = sc.AnnData(
        X=adata.layers['counts'].copy(),
        var=adata.var.copy()
    )
    sc.pp.normalize_total(adata_temp, target_sum=1e4)
    sc.pp.log1p(adata_temp)
    
    stress_indices = [i for i, g in enumerate(adata.var_names) if g in stress_genes_in_data]
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
    
    adata.obs['stress_score'] = stress_score
    
    del adata_temp, stress_expr
    gc.collect()
    
    print(f"  Mean stress: {stress_score.mean():.3f}")
else:
    adata.obs['stress_score'] = 0.0
    print("  Warning: Too few stress genes, set score=0")

# Calculate cell cycle scores
print("\nCalculating cell cycle scores...")
s_genes = [g for g in S_GENES if g in adata.var_names]
g2m_genes = [g for g in G2M_GENES if g in adata.var_names]
print(f"  S genes: {len(s_genes)}/{len(S_GENES)}")
print(f"  G2M genes: {len(g2m_genes)}/{len(G2M_GENES)}")

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
    
    print(f"  S_score mean: {adata.obs['S_score'].mean():.3f}")
    print(f"  G2M_score mean: {adata.obs['G2M_score'].mean():.3f}")
else:
    adata.obs['S_score'] = 0.0
    adata.obs['G2M_score'] = 0.0
    adata.obs['phase'] = 'G1'
    print("  Warning: Too few cell cycle genes, set scores=0")

print("\nAll covariates calculated")
print("Data preparation complete")


# ==============================================================================
# Step 2: scVI Integration (FIXED HVG workflow)
# ==============================================================================

print("\n" + "="*70)
print("Step 2: scVI Integration (FIXED HVG workflow)")
print("="*70)

hvg_file = scvi_model_dir / "hvg_genes.txt"

# ----- Step 2.1: Load or compute HVG list -----
print("\n" + "-"*70)
print("Step 2.1: HVG Selection")
print("-"*70)

hvg_genes_list = None
hvg_selection_method = "none"

if USE_HVG_FOR_SCVI:
    print(f"\nUsing HVG mode (target: {N_HVG_SCVI} genes)")
    
    # Try to load existing HVG list
    if hvg_file.exists() and not RETRAIN_SCVI:
        print("\nFound existing HVG list, loading...")
        try:
            hvg_genes_saved = pd.read_csv(hvg_file, header=None)[0].astype(str).tolist()
            hvg_genes_in_data = [g for g in hvg_genes_saved if g in set(adata.var_names)]
            
            if len(hvg_genes_in_data) >= len(hvg_genes_saved) * 0.95:
                hvg_genes_list = hvg_genes_in_data
                hvg_selection_method = "loaded_from_file"
                print(f"Loaded {len(hvg_genes_list):,} HVG genes from file")
            else:
                print(f"Warning: Only {len(hvg_genes_in_data)}/{len(hvg_genes_saved)} genes found")
                print("  Will recompute HVG...")
                hvg_genes_list = None
        except Exception as e:
            print(f"Warning: Failed to load HVG file: {e}")
            print("  Will recompute HVG...")
            hvg_genes_list = None
    
    # Compute HVG if not loaded
    if hvg_genes_list is None:
        print("\nComputing HVG...")
        
        # Try batch-aware HVG first
        try:
            print("  Attempting batch-aware HVG selection...")
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                batch_key=BATCH_KEY, flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "batch-aware"
            print("  Batch-aware HVG successful")
        except Exception as e:
            print(f"  Warning: Batch-aware failed: {e}")
            print("  Falling back to non-batch-aware...")
            sc.pp.highly_variable_genes(
                adata, layer="counts", n_top_genes=N_HVG_SCVI,
                flavor=HVG_FLAVOR, subset=False
            )
            hvg_selection_method = "non-batch-aware"
            print("  Non-batch-aware HVG successful")
        
        # Extract HVG genes
        hvg_mask = adata.var["highly_variable"].values
        hvg_genes_list = adata.var_names[hvg_mask].tolist()
        print(f"Selected {len(hvg_genes_list):,} HVG genes")
        
        # Force markers into HVG
        if FORCE_MARKERS_IN_HVG:
            print("\n" + "-"*70)
            print("Forcing biological markers into HVG")
            print("-"*70)
            
            adata.var["highly_variable_original"] = adata.var["highly_variable"].values
            
            report = force_include_markers_in_hvg(
                adata,
                n_top_genes=N_HVG_SCVI,
                forced_markers=FORCED_MARKERS,
                case_insensitive=FORCE_MARKERS_CASE_INSENSITIVE,
                hvg_col="highly_variable"
            )
            print("Forced markers report:")
            for k, v in report.items():
                print(f"  {k}: {v}")
            
            # Update HVG list
            hvg_mask = adata.var["highly_variable"].values
            hvg_genes_list = adata.var_names[hvg_mask].tolist()
            print(f"Final HVG count: {len(hvg_genes_list):,}")
        
        # Save HVG list
        scvi_model_dir.mkdir(parents=True, exist_ok=True)
        pd.Series(hvg_genes_list).to_csv(hvg_file, index=False, header=False)
        print(f"HVG list saved to: {hvg_file}")

else:
    print("\nWarning: Using all genes (no HVG filtering)")
    hvg_genes_list = adata.var_names.tolist()
    hvg_selection_method = "all_genes"

print(f"\nHVG selection complete")
print(f"  Method: {hvg_selection_method}")
print(f"  Final gene count: {len(hvg_genes_list):,}")

# ----- Preserve full gene matrix in .raw -----
# CRITICAL: Save full genes BEFORE subsetting to HVG
# Uses shared memory (0 extra cost) - QUICK_REFERENCE_MEMORY best practice
print("\n" + "-"*70)
print("Preserving Full Gene Matrix")
print("-"*70)

print("\nSaving full gene matrix to .raw (shared memory)...")
adata.raw = sc.AnnData(
    X=adata.layers["counts"],  # No .copy()! Shared memory pointer
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"  Preserved: {adata.raw.n_vars:,} genes")
print(f"  Memory cost: 0 bytes (shared pointer)")
print(f"  Benefit: Full genes available for marker visualization")
print(f"\nNote: Can now use use_raw=True in plotting/analysis functions")
print(f"  Example: sc.pl.umap(adata, color=['YOUR_MARKER'], use_raw=True)")


# ----- Step 2.2: Construct scVI AnnData with HVG -----
print("\n" + "-"*70)
print("Step 2.2: Construct scVI AnnData")
print("-"*70)

print(f"\nApplying HVG filter to data...")
hvg_indices = np.flatnonzero(adata.var_names.isin(hvg_genes_list))
print(f"  HVG indices: {len(hvg_indices):,}")

adata_scvi = sc.AnnData(
    X=adata.layers["counts"][:, hvg_indices].copy(),
    obs=adata.obs[REQUIRED_COVARIATES].copy(),
    var=adata.var.iloc[hvg_indices].copy()
)
adata_scvi.var_names = adata.var_names[hvg_indices]
print(f"scVI AnnData created: {adata_scvi.n_obs:,} x {adata_scvi.n_vars:,}")

# Verify covariates
missing_covariates = [k for k in REQUIRED_COVARIATES if k not in adata_scvi.obs.columns]
if missing_covariates:
    raise ValueError(f"Missing covariate columns: {missing_covariates}")
print(f"All required covariates present: {REQUIRED_COVARIATES}")


# ----- Step 2.3: Train or load scVI model -----
print("\n" + "-"*70)
print("Step 2.3: scVI Model Training/Loading")
print("-"*70)

scvi_model_exists = scvi_model_dir.exists() and (scvi_model_dir / "model.pt").exists()

if scvi_model_exists and not RETRAIN_SCVI:
    print("Found existing scVI model, loading...")
    try:
        scvi_model = scvi.model.SCVI.load(scvi_model_dir, adata=adata_scvi)
        print("scVI model loaded successfully")
    except Exception as e:
        print(f"Warning: Load failed: {e}")
        print("  Will retrain model...")
        scvi_model_exists = False

if not scvi_model_exists or RETRAIN_SCVI:
    print("Training new scVI model...")
    
    if RETRAIN_SCVI and scvi_model_exists:
        print("  Info: RETRAIN_SCVI=True, forcing retraining")
    
    # Setup with covariates
    print("\nSetting up scVI model...")
    scvi.model.SCVI.setup_anndata(
        adata_scvi, 
        layer=None, 
        batch_key=BATCH_KEY,
        categorical_covariate_keys=[TISSUE_KEY],
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES 
                                    if k not in [BATCH_KEY, TISSUE_KEY]]
    )
    print("Setup complete")
    
    # Initialize model
    print("\nInitializing scVI model...")
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
    
    print(f"Model initialized:")
    print(f"  n_latent: {SCVI_N_LATENT}")
    print(f"  n_layers: {SCVI_N_LAYERS}")
    print(f"  encode_covariates: {SCVI_ENCODE_COVARIATES}")
    print(f"  use_layer_norm: {SCVI_USE_LAYER_NORM}")
    
    # Train
    print("\nTraining scVI model...")
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
    print("Training complete")
    
    # Save model
    scvi_model.save(scvi_model_dir, overwrite=True)
    print(f"Model saved to: {scvi_model_dir}")

# Generate latent representation
print("\nGenerating scVI latent representation...")
latent_df = pd.DataFrame(
    scvi_model.get_latent_representation(),
    index=adata_scvi.obs_names
)
latent_aligned = latent_df.reindex(adata.obs_names)
if latent_aligned.isna().any().any():
    raise RuntimeError("scVI latent reindex produced NaNs")
adata.obsm['X_scvi'] = latent_aligned.to_numpy().astype(np.float32)
print(f"scVI latent shape: {adata.obsm['X_scvi'].shape}")

# Cleanup
del adata_scvi
gc.collect()

# Compute neighbors and UMAP
print("\nComputing neighbors...")
sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=UMAP_N_NEIGHBORS, 
                key_added='neighbors_scvi', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, 
           neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()

print("Computing Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scvi', 
             neighbors_key='neighbors_scvi', random_state=RANDOM_SEED)
print(f"Leiden clusters: {adata.obs['leiden_scvi'].nunique()}")

print("\nscVI integration complete")


# ==============================================================================
# Step 3: CellTypist Annotation
# ==============================================================================

print("\n" + "="*70)
print("Step 3: CellTypist Annotation")
print("="*70)

print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL_PATH}")
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
print(f"Model loaded, features: {len(model_features)}")

# Detect gene symbol column
print("\nDetecting gene symbol column...")
symbol_col, gene_symbols = detect_gene_symbol_column(adata, verbose=True)

# Match genes with CellTypist model
print("\nMatching genes with CellTypist model...")
common_genes = gene_symbols.isin(model_features)
n_common = common_genes.sum()
print(f"  Overlapping genes: {n_common:,}/{len(model_features)} ({n_common/len(model_features)*100:.1f}%)")

if n_common < 500:
    print("\nERROR: Insufficient gene overlap with CellTypist model")
    print(f"  Current overlap: {n_common} genes")
    print(f"  Required minimum: 500 genes")
    print("\nPossible causes:")
    print("  1. Wrong CellTypist model for your tissue/species")
    print("  2. Gene naming mismatch (ENSEMBL IDs vs gene symbols)")
    print("  3. Different gene annotation version")
    print("\nSuggestions:")
    print("  - Check if using correct CellTypist model")
    print("  - Verify gene symbols are properly annotated")
    print("  - Try a different CellTypist model")
    raise RuntimeError(f"Too few overlapping genes: {n_common}")

# Construct CellTypist AnnData
print(f"\nConstructing CellTypist AnnData...")
X_subset = adata.layers["counts"][:, common_genes.values]

adata_celltypist = sc.AnnData(
    X=X_subset.copy(),
    obs=adata.obs.copy(),
    var=adata.var.loc[common_genes].copy()
)

# Set var_names to matched gene symbols
if symbol_col is not None:
    adata_celltypist.var_names = gene_symbols[common_genes].values
else:
    adata_celltypist.var_names = adata.var_names[common_genes]

print(f"CellTypist AnnData: {adata_celltypist.n_obs:,} x {adata_celltypist.n_vars:,}")

# Normalize
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)
sc.pp.log1p(adata_celltypist)
print("Data prepared")

# Run CellTypist
print("\nRunning CellTypist prediction...")
predictions = celltypist.annotate(adata_celltypist, model=celltypist_model, 
                                   majority_voting=True)

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

# Confidence scores
conf_candidates = ['conf_score', 'confidence', 'confidence_score', 'prob']
confidence_column = next((col for col in conf_candidates 
                          if col in predicted_labels.columns), None)
if confidence_column:
    adata.obs['celltypist_conf_score'] = pred_df[confidence_column].reindex(adata.obs_names).values
else:
    adata.obs['celltypist_conf_score'] = 1.0

# Cleanup
del adata_celltypist, X_subset
gc.collect()

print(f"\nCellTypist summary:")
print(f"  Unique types: {adata.obs['celltypist_majority_voting'].nunique()}")
print(f"  Mean confidence: {adata.obs['celltypist_conf_score'].mean():.3f}")
print(f"  Top 5 cell types:")
for ct, count in adata.obs['celltypist_majority_voting'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")

print("\nCellTypist annotation complete")


# ==============================================================================
# Step 4A: scANVI-Existing (FIXED workflow)
# ==============================================================================

print("\n" + "="*70)
print("Step 4A: scANVI-Existing (Using cell_type metadata)")
print("="*70)

print(f"\nTraining scANVI using '{EXISTING_CELLTYPE_KEY}' metadata")
print(f"   With Unknown cleaning (purity + quality filters)")
print(f"   With n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")
print(f"   CRITICAL: Using same HVG as scVI")

# Clean existing labels
print("\nCleaning existing labels...")
labels_existing = adata.obs[EXISTING_CELLTYPE_KEY].astype(str).copy()
n_original_existing = len(labels_existing)

# Neighborhood purity
print("  Computing neighborhood purity...")
purity_array = calculate_label_agreement_purity(adata, EXISTING_CELLTYPE_KEY, 'neighbors_scvi')
adata.obs['label_purity_existing'] = purity_array

low_purity_mask = purity_array < EXISTING_LABEL_PURITY_THRESHOLD
n_low_purity = low_purity_mask.sum()
labels_existing[low_purity_mask] = "Unknown"

print(f"    Threshold: {EXISTING_LABEL_PURITY_THRESHOLD}")
print(f"    Mean purity: {purity_array.mean():.3f}")
print(f"    Low purity -> Unknown: {n_low_purity:,} ({n_low_purity/n_original_existing*100:.1f}%)")

# Quality filtering
quality_mask = (
    (adata.obs['pct_counts_mt'] > EXISTING_LABEL_MT_THRESHOLD) |
    (adata.obs['stress_score'] > np.percentile(adata.obs['stress_score'], 
                                                EXISTING_LABEL_STRESS_PERCENTILE))
)
n_low_quality = quality_mask.sum()
labels_existing[quality_mask] = "Unknown"

print(f"    MT% threshold: {EXISTING_LABEL_MT_THRESHOLD}%")
print(f"    Stress percentile: {EXISTING_LABEL_STRESS_PERCENTILE}th")
print(f"    Low quality -> Unknown: {n_low_quality:,} ({n_low_quality/n_original_existing*100:.1f}%)")

# Summary
n_unknown_existing = (labels_existing == "Unknown").sum()
n_labeled_existing = n_original_existing - n_unknown_existing

print(f"\n  Summary:")
print(f"    Total: {n_original_existing:,}")
print(f"    Final labeled: {n_labeled_existing:,} ({n_labeled_existing/n_original_existing*100:.1f}%)")
print(f"    Final Unknown: {n_unknown_existing:,} ({n_unknown_existing/n_original_existing*100:.1f}%)")

adata.obs['scanvi_label_existing_cleaned'] = pd.Categorical(labels_existing)

# CRITICAL: Construct scANVI AnnData with SAME HVG as scVI
print("\nPreparing scANVI-Existing data (using scVI HVG)...")

adata_scanvi_existing = sc.AnnData(
    X=adata.layers["counts"][:, hvg_indices].copy(),
    obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_existing_cleaned']].copy(),
    var=adata.var.iloc[hvg_indices].copy()
)
adata_scanvi_existing.var_names = adata.var_names[hvg_indices]
print(f"scANVI-Existing AnnData: {adata_scanvi_existing.n_obs:,} x {adata_scanvi_existing.n_vars:,}")

# Ensure 'Unknown' category exists
labels_cat = adata_scanvi_existing.obs['scanvi_label_existing_cleaned'].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi_existing.obs['scanvi_label_existing_cleaned'] = labels_cat

# Train or load model
print("\n" + "-"*70)
print("scANVI-Existing Model: Train or Load")
print("-"*70)

scanvi_existing_exists = (scanvi_existing_model_dir.exists() and 
                          (scanvi_existing_model_dir / "model.pt").exists())

if scanvi_existing_exists and not RETRAIN_SCANVI_EXISTING:
    print("Found existing scANVI-Existing model, loading...")
    try:
        scanvi_existing_model = scvi.model.SCANVI.load(
            scanvi_existing_model_dir, adata=adata_scanvi_existing
        )
        print("Model loaded successfully")
    except Exception as e:
        print(f"Warning: Load failed: {e}")
        print("  Will retrain model...")
        scanvi_existing_exists = False

if not scanvi_existing_exists or RETRAIN_SCANVI_EXISTING:
    print("Training new scANVI-Existing model...")
    
    if RETRAIN_SCANVI_EXISTING and scanvi_existing_exists:
        print("  Info: RETRAIN_SCANVI_EXISTING=True, forcing retraining")
    
    # Setup with covariates
    print("\nSetting up scANVI model...")
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_existing,
        layer=None,
        batch_key=BATCH_KEY,
        labels_key='scanvi_label_existing_cleaned',
        unlabeled_category="Unknown",
        categorical_covariate_keys=[TISSUE_KEY],
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES 
                                    if k not in [BATCH_KEY, TISSUE_KEY]]
    )
    print("Setup complete")
    
    # CRITICAL: Create from scVI model (transfer learning)
    print("\nCRITICAL: Creating scANVI from scVI model...")
    scanvi_existing_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model, 
        unlabeled_category="Unknown",
        adata=adata_scanvi_existing, 
        labels_key='scanvi_label_existing_cleaned'
    )
    print("scANVI initialized from scVI")
    
    # Train
    print("\nTraining scANVI-Existing...")
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
    print("Training complete")
    
    # Save
    scanvi_existing_model.save(scanvi_existing_model_dir, overwrite=True)
    print(f"Model saved to: {scanvi_existing_model_dir}")

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

# Cleanup
del adata_scanvi_existing
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors...")
sc.pp.neighbors(adata, use_rep='X_scanvi_existing', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_existing'] = adata.obsm['X_umap'].copy()

print("Computing Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_existing',
             neighbors_key='neighbors_scanvi_existing', random_state=RANDOM_SEED)

print(f"\nscANVI-Existing complete")
print(f"  Cell types: {adata.obs['scanvi_predictions_existing'].nunique()}")
print(f"  Top 5:")
for ct, count in adata.obs['scanvi_predictions_existing'].value_counts().head(5).items():
    print(f"    {ct}: {count:,}")


# ==============================================================================
# Step 4B: scANVI-CellTypist (FIXED workflow)
# ==============================================================================

print("\n" + "="*70)
print("Step 4B: scANVI-CellTypist (Using CellTypist predictions)")
print("="*70)

print(f"\nTraining scANVI using CellTypist predictions")
print(f"   With n_samples_per_label={SCANVI_N_SAMPLES_PER_LABEL}")
print(f"   CRITICAL: Using same HVG as scVI")

# Prepare labels
print("\nPreparing CellTypist labels...")
labels_celltypist = adata.obs['celltypist_majority_voting'].astype(str).copy()
n_original_celltypist = len(labels_celltypist)

# Count initial Unknown
n_unknown_celltypist_initial = (labels_celltypist == "Unknown").sum()
n_labeled_celltypist = n_original_celltypist - n_unknown_celltypist_initial

print(f"  Total cells: {n_original_celltypist:,}")
print(f"  Labeled: {n_labeled_celltypist:,} ({n_labeled_celltypist/n_original_celltypist*100:.1f}%)")
print(f"  Unknown: {n_unknown_celltypist_initial:,} ({n_unknown_celltypist_initial/n_original_celltypist*100:.1f}%)")

adata.obs['scanvi_label_celltypist'] = pd.Categorical(labels_celltypist)

# CRITICAL: Construct scANVI AnnData with SAME HVG as scVI
print("\nPreparing scANVI-CellTypist data (using scVI HVG)...")

adata_scanvi_celltypist = sc.AnnData(
    X=adata.layers["counts"][:, hvg_indices].copy(),
    obs=adata.obs[REQUIRED_COVARIATES + ['scanvi_label_celltypist']].copy(),
    var=adata.var.iloc[hvg_indices].copy()
)
adata_scanvi_celltypist.var_names = adata.var_names[hvg_indices]
print(f"scANVI-CellTypist AnnData: {adata_scanvi_celltypist.n_obs:,} x {adata_scanvi_celltypist.n_vars:,}")

# Ensure 'Unknown' category exists
labels_cat = adata_scanvi_celltypist.obs['scanvi_label_celltypist'].astype('category')
if 'Unknown' not in labels_cat.cat.categories:
    labels_cat = labels_cat.cat.add_categories(['Unknown'])
    adata_scanvi_celltypist.obs['scanvi_label_celltypist'] = labels_cat

# Train or load model
print("\n" + "-"*70)
print("scANVI-CellTypist Model: Train or Load")
print("-"*70)

scanvi_celltypist_exists = (scanvi_celltypist_model_dir.exists() and 
                            (scanvi_celltypist_model_dir / "model.pt").exists())

if scanvi_celltypist_exists and not RETRAIN_SCANVI_CELLTYPIST:
    print("Found existing scANVI-CellTypist model, loading...")
    try:
        scanvi_celltypist_model = scvi.model.SCANVI.load(
            scanvi_celltypist_model_dir, adata=adata_scanvi_celltypist
        )
        print("Model loaded successfully")
    except Exception as e:
        print(f"Warning: Load failed: {e}")
        print("  Will retrain model...")
        scanvi_celltypist_exists = False

if not scanvi_celltypist_exists or RETRAIN_SCANVI_CELLTYPIST:
    print("Training new scANVI-CellTypist model...")
    
    if RETRAIN_SCANVI_CELLTYPIST and scanvi_celltypist_exists:
        print("  Info: RETRAIN_SCANVI_CELLTYPIST=True, forcing retraining")
    
    # Setup with covariates
    print("\nSetting up scANVI model...")
    scvi.model.SCANVI.setup_anndata(
        adata_scanvi_celltypist,
        layer=None,
        batch_key=BATCH_KEY,
        labels_key='scanvi_label_celltypist',
        unlabeled_category="Unknown",
        categorical_covariate_keys=[TISSUE_KEY],
        continuous_covariate_keys=[k for k in REQUIRED_COVARIATES 
                                    if k not in [BATCH_KEY, TISSUE_KEY]]
    )
    print("Setup complete")
    
    # CRITICAL: Create from scVI model (transfer learning)
    print("\nCRITICAL: Creating scANVI from scVI model...")
    scanvi_celltypist_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category="Unknown",
        adata=adata_scanvi_celltypist,
        labels_key='scanvi_label_celltypist'
    )
    print("scANVI initialized from scVI")
    
    # Train
    print("\nTraining scANVI-CellTypist...")
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
    print("Training complete")
    
    # Save
    scanvi_celltypist_model.save(scanvi_celltypist_model_dir, overwrite=True)
    print(f"Model saved to: {scanvi_celltypist_model_dir}")

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

# Cleanup
del adata_scanvi_celltypist
gc.collect()

# Neighbors and UMAP
print("\nComputing neighbors...")
sc.pp.neighbors(adata, use_rep='X_scanvi_celltypist', n_neighbors=UMAP_N_NEIGHBORS,
                key_added='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print("Computing UMAP...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD,
           neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)
adata.obsm['X_umap_scanvi_celltypist'] = adata.obsm['X_umap'].copy()

print("Computing Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_scanvi_celltypist',
             neighbors_key='neighbors_scanvi_celltypist', random_state=RANDOM_SEED)

print(f"\nscANVI-CellTypist complete")
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

# Generate comparison plots
print("\nGenerating comparison plots...")

fig, axes = plt.subplots(2, 4, figsize=(24, 12))

# Row 1: scANVI-Existing
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_existing'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0, 0], show=False, 
           title='scANVI-Existing: Batch')
sc.pl.umap(adata, color='scanvi_label_existing_cleaned', ax=axes[0, 1], show=False, 
           title='scANVI-Existing: Input (cleaned)')
sc.pl.umap(adata, color='scanvi_predictions_existing', ax=axes[0, 2], show=False, 
           title='scANVI-Existing: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_existing', ax=axes[0, 3], show=False, 
           title='scANVI-Existing: Leiden')

# Row 2: scANVI-CellTypist
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1, 0], show=False, 
           title='scANVI-CellTypist: Batch')
sc.pl.umap(adata, color='scanvi_label_celltypist', ax=axes[1, 1], show=False, 
           title='scANVI-CellTypist: Input')
sc.pl.umap(adata, color='scanvi_predictions_celltypist', ax=axes[1, 2], show=False, 
           title='scANVI-CellTypist: Predictions')
sc.pl.umap(adata, color='leiden_scanvi_celltypist', ax=axes[1, 3], show=False, 
           title='scANVI-CellTypist: Leiden')

plt.tight_layout()
plt.savefig(output_dir['figures'] / f"dual_scanvi_comparison.{FIGURE_FORMAT}", 
            dpi=DPI, bbox_inches='tight')
plt.close()

print("Comparison plots saved")

# Set default UMAP to CellTypist-based
adata.obsm['X_umap'] = adata.obsm['X_umap_scanvi_celltypist'].copy()

# Save final dataset
output_file = output_dir['base'] / "adata_allcells_dual_scanvi_final_v2.5_PRODUCTION.h5ad"
print(f"\nSaving final dataset to: {output_file}")
adata.write_h5ad(output_file, compression='gzip')

file_size = output_file.stat().st_size / (1024**3)
print(f"Saved ({file_size:.2f} GB)")

# Generate summary
print("\n" + "="*70)
print("ANALYSIS SUMMARY v2.5 PRODUCTION")
print("="*70)

summary_lines = []
summary_lines.append("="*70)
summary_lines.append("All-Cells scVI/scANVI Analysis v2.5 PRODUCTION")
summary_lines.append("="*70)

summary_lines.append("\n[ Key Improvements in v2.5 ]")
summary_lines.append("  1. FIXED HVG workflow (correct execution order)")
summary_lines.append("  2. Simplified model control (clear retrain flags)")
summary_lines.append("  3. Integrated QUICK_REFERENCE_MEMORY best practices")
summary_lines.append("  4. Complete English documentation (UTF-8 compatible)")
summary_lines.append("  5. Safety mechanisms (HVG check, error handling, status reporting)")

summary_lines.append("\n[ HVG Selection ]")
summary_lines.append(f"  Method: {hvg_selection_method}")
summary_lines.append(f"  Final count: {len(hvg_genes_list):,} genes")
summary_lines.append(f"  Forced markers: {FORCE_MARKERS_IN_HVG}")

summary_lines.append("\n[ scVI Integration ]")
summary_lines.append(f"  Latent dimensions: {SCVI_N_LATENT}")
summary_lines.append(f"  Total cells: {adata.n_obs:,}")
summary_lines.append(f"  Leiden clusters: {adata.obs['leiden_scvi'].nunique()}")

summary_lines.append("\n[ scANVI-Existing ]")
summary_lines.append(f"  Input: cleaned cell_type labels")
summary_lines.append(f"  Total cells: {n_original_existing:,}")
summary_lines.append(f"  Labeled: {n_labeled_existing:,} ({n_labeled_existing/n_original_existing*100:.1f}%)")
summary_lines.append(f"  Unknown: {n_unknown_existing:,} ({n_unknown_existing/n_original_existing*100:.1f}%)")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_existing'].nunique()}")

summary_lines.append("\n[ scANVI-CellTypist ]")
summary_lines.append(f"  Input: CellTypist predictions")
summary_lines.append(f"  Total cells: {n_original_celltypist:,}")
summary_lines.append(f"  Labeled: {n_labeled_celltypist:,} ({n_labeled_celltypist/n_original_celltypist*100:.1f}%)")
summary_lines.append(f"  Unknown: {n_unknown_celltypist_initial:,} ({n_unknown_celltypist_initial/n_original_celltypist*100:.1f}%)")
summary_lines.append(f"  Final types: {adata.obs['scanvi_predictions_celltypist'].nunique()}")

summary_lines.append("\n[ Key Outputs ]")
summary_lines.append("  scVI:")
summary_lines.append("    - X_scvi: latent representation")
summary_lines.append("    - X_umap_scvi: UMAP coordinates")
summary_lines.append("    - leiden_scvi: clustering")
summary_lines.append("  scANVI-Existing:")
summary_lines.append("    - scanvi_predictions_existing: final predictions")
summary_lines.append("    - scanvi_label_existing_cleaned: input labels (with Unknown)")
summary_lines.append("    - X_scanvi_existing: latent representation")
summary_lines.append("    - X_umap_scanvi_existing: UMAP coordinates")
summary_lines.append("  scANVI-CellTypist:")
summary_lines.append("    - scanvi_predictions_celltypist: final predictions")
summary_lines.append("    - scanvi_label_celltypist: input labels")
summary_lines.append("    - X_scanvi_celltypist: latent representation")
summary_lines.append("    - X_umap_scanvi_celltypist: UMAP coordinates")

summary_lines.append("\n" + "="*70)
summary_lines.append("Analysis complete!")
summary_lines.append("="*70)

summary_text = '\n'.join(summary_lines)
print(summary_text)

summary_file = output_dir['base'] / "analysis_summary_v2.5_PRODUCTION.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

print(f"\nSummary saved to: {summary_file}")

print("\n" + "="*70)
print("PIPELINE COMPLETE (v2.5 PRODUCTION)")
print("="*70)

print(f"\nKey Features:")
print(f"  - HVG workflow: Load -> Apply -> Train (correct order)")
print(f"  - Model control: Simplified flags (RETRAIN_*)")
print(f"  - Gene sets: Consistent across all models")
print(f"  - Documentation: Complete English, UTF-8 compatible")
print(f"  - Best practices: Integrated from QUICK_REFERENCE_MEMORY v2.11-12")

print("\nExpected Benefits:")
print("  - Consistent gene sets across all models")
print("  - Reliable model loading/saving")
print("  - Better reproducibility")
print("  - Clear control over retraining")
print("  - Production-ready code quality")

print("\n" + "="*70)
