#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
B Cell: scVI â†’ CellTypist â†’ scANVI Pipeline v3.5.2
===================================================

Cell Types Covered:
- Naive B cells
- Memory B cells (IgM, IgG, IgA)
- Plasma cells / Plasmablasts
- Germinal center B cells
- Activated B cells
- Regulatory B cells (Bregs)

Adapted from T cell pipeline v3.5.1 - Critical Bug Fixes Applied:
- âœ… P0-FIX: scANVI uses dedicated adata_model (labels alignment)
- âœ… P0-FIX: predict(soft=True).max() changed to np.asarray()
- âœ… P1-FIX: CellTypist majority voting properly integrated
- âœ… P1-FIX: CellTypist features validation (hard check)
- âœ… P1-FIX: Pretrained HVG uses var_names (unique, not symbol_base)
- âœ… P2-FIX: Memory optimization (no full X copy, use layers)
- âœ… P2-FIX: Minimal adata_model construction

B Cell-Specific Features:
1. Immunoglobulin isotype markers (IgM, IgG, IgA, IgE)
2. Plasma cell maturation markers
3. Germinal center reaction markers
4. B cell activation and proliferation
5. Type 2 inflammation and IgE production (relevant for CRSwNP/allergic disease)
6. B-T cell interaction markers

Input: Raw B cell data with counts layer
Output: adata_bcell_FINAL.h5ad + count statistics

Author: r2end
Date: 2025-01-11
Version: 3.5.2 - Added counts recovery from external h5ad
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime
# Core libraries
import numpy as np
import pandas as pd
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt
import time
import gc

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
    print("âš ï¸  mygene not available, will use local gene names only")


# ============================================================================
# COUNTS RECOVERY FUNCTION (NEW in v3.5.2)
# ============================================================================

from scipy import sparse

def recover_counts_from_external(adata, external_file_path, counts_layer_name='counts', verbose=True):
    """
    Recover counts layer by matching cell IDs from external h5ad file
    
    Parameters:
        adata: Target AnnData missing counts layer
        external_file_path: Path to h5ad with counts
        counts_layer_name: Name of counts layer (default: 'counts')
        verbose: Print progress info
    
    Returns:
        adata: Updated AnnData with recovered counts layer
    """
    if verbose:
        print(f"\n{'='*80}")
        print("COUNTS RECOVERY FROM EXTERNAL FILE")
        print("="*80)
        print(f"\nLoading external file: {external_file_path}")
    
    adata_ext = sc.read_h5ad(external_file_path)
    
    if verbose:
        print(f"✓ External data loaded:")
        print(f"  Cells: {adata_ext.n_obs:,}")
        print(f"  Genes: {adata_ext.n_vars:,}")
        print(f"  Layers: {list(adata_ext.layers.keys())}")
    
    # Check counts source
    if counts_layer_name in adata_ext.layers:
        counts_source = adata_ext.layers[counts_layer_name]
        ext_gene_names = adata_ext.var_names
        source_name = f"layers['{counts_layer_name}']"
    elif adata_ext.raw is not None:
        counts_source = adata_ext.raw.X
        ext_gene_names = adata_ext.raw.var_names
        source_name = ".raw.X"
    else:
        raise ValueError(f"No '{counts_layer_name}' layer or .raw found")
    
    if verbose:
        print(f"  Using counts from: {source_name}")
    
    # Match cells
    current_cells = set(adata.obs.index)
    external_cells = set(adata_ext.obs.index)
    matching_cells = current_cells & external_cells
    
    if len(matching_cells) == 0:
        raise ValueError("No matching cells found")
    
    if verbose:
        print(f"\nCell matching:")
        print(f"  Current: {len(current_cells):,} cells")
        print(f"  External: {len(external_cells):,} cells")
        print(f"  Matching: {len(matching_cells):,} ({len(matching_cells)/len(current_cells)*100:.1f}%)")
    
    missing = current_cells - external_cells
    if len(missing) > 0:
        print(f"  ⚠️  {len(missing):,} cells not found → will have zero counts")
    
    # Match genes
    if 'symbol_base' in adata.var.columns:
        curr_dict = dict(zip(adata.var.index, adata.var['symbol_base']))
    else:
        curr_dict = dict(zip(adata.var.index, adata.var.index))
    
    if 'symbol_base' in adata_ext.var.columns:
        ext_dict = dict(zip(adata_ext.var.index, adata_ext.var['symbol_base']))
    else:
        ext_dict = dict(zip(adata_ext.var.index, ext_gene_names))
    
    gene_map = {}
    for i, (gid, gname) in enumerate(curr_dict.items()):
        for j, (eid, ename) in enumerate(ext_dict.items()):
            if gname == ename or gid == eid:
                gene_map[i] = j
                break
    
    if verbose:
        print(f"\nGene matching:")
        print(f"  Current: {adata.n_vars:,} genes")
        print(f"  External: {len(ext_gene_names):,} genes")
        print(f"  Matching: {len(gene_map):,} ({len(gene_map)/adata.n_vars*100:.1f}%)")
        if len(gene_map) < adata.n_vars * 0.5:
            print(f"  ⚠️  <50% matched - check gene names")
    
    # Build recovered matrix
    recovered = sparse.lil_matrix((adata.n_obs, adata.n_vars), dtype=np.float32)
    
    cell_idx = {cid: i for i, cid in enumerate(adata.obs.index)}
    ext_idx = {cid: i for i, cid in enumerate(adata_ext.obs.index)}
    
    if verbose:
        print(f"\nCopying counts...")
    
    for cid in matching_cells:
        ci = cell_idx[cid]
        ei = ext_idx[cid]
        
        if sparse.issparse(counts_source):
            counts_vec = counts_source[ei, :].toarray().flatten()
        else:
            counts_vec = counts_source[ei, :]
        
        for curr_g, ext_g in gene_map.items():
            recovered[ci, curr_g] = counts_vec[ext_g]
    
    recovered = recovered.tocsr()
    adata.layers['counts'] = recovered
    
    adata.uns['counts_recovery_stats'] = {
        'external_file': external_file_path,
        'cells_matched': len(matching_cells),
        'genes_matched': len(gene_map),
        'match_rate_cells': len(matching_cells) / len(current_cells),
        'match_rate_genes': len(gene_map) / adata.n_vars
    }
    
    if verbose:
        print(f"\n✓ Counts recovered:")
        print(f"  Shape: {recovered.shape}")
        print(f"  Non-zero: {recovered.nnz:,}")
        print(f"  Total counts: {recovered.sum():,.0f}")
        print(f"  Mean/cell: {recovered.sum()/adata.n_obs:,.0f}")
    
    return adata

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output paths
INPUT_FILE = "/home/h2048/data/R/0104/B/b_final_20260104.h5ad"  # ⭐ UPDATE THIS
OUTPUT_DIR = f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/celltypist_bcell"

# ⭐ External counts source (NEW in v3.5.2)
# If your current data lacks counts layer, provide path to h5ad with raw counts
# Counts will be recovered by matching cell IDs (adata.obs.index)
EXTERNAL_COUNTS_FILE = "/home/h2048/data/py/1208/b_scvi_celltypist_scanvi_v1_1/adata_bcell_scvi_celltypist_scanvi_final_v2.h5ad"  # Path to h5ad with counts layer, or None
# Example: "/path/to/original_full_dataset_with_counts.h5ad"

# â­ Pre-trained model loading
PRETRAINED_SCVI_MODEL = None  # Set path to reuse model, or None to train
# Example: "/path/to/previous_run/models/scvi_model"

# Model hyperparameters
N_LATENT = 50  # Adjust based on dataset complexity (30-75)
N_LAYERS = 2
N_HVG = 4000

# Training parameters
SCVI_MAX_EPOCHS = 400
SCANVI_MAX_EPOCHS = 200
BATCH_SIZE = 2048
LEARNING_RATE = 1e-3
EARLY_STOPPING = True
EARLY_STOPPING_PATIENCE = 50

# CellTypist
CELLTYPIST_MODEL = '/home/h2048/data/source/reference/celltypist_models/Immune_All_Low.pkl'  # Standard immune model covers B cells
CELLTYPIST_MAJORITY_VOTING = True
CELLTYPIST_MIN_CONFIDENCE = 0.5

# â­ Rare type filtering
MIN_CELLS_PER_TYPE = 50  # Types with <50 cells â†’ Unknown

# ============================================================================
# B CELL CORE MARKERS
# ============================================================================

# â­ B CELL CORE MARKERS - Force include in HVG
# These markers are essential for B cell subtype identification
# Will be automatically added to HVG if not selected by variance
BCELL_CORE_MARKERS = {
    # ========== PAN-B CELL MARKERS ==========
    'Pan_B': ['CD19', 'CD79A', 'CD79B', 'MS4A1', 'PTPRC'],  # CD20 = MS4A1
    
    # ========== NAIVE B CELLS ==========
    'Naive_B': ['IGHD', 'IGHM', 'TCL1A', 'IL4R', 'FCER2', 'CD38'],
    
    # ========== MEMORY B CELLS ==========
    'Memory_B': ['CD27', 'TNFRSF13B', 'TNFRSF13C'],  # CD27, TACI, BAFFR
    
    # Immunoglobulin Isotypes (Class-Switched Memory)
    'IgM_Memory': ['IGHM', 'CD27'],
    'IgG_Memory': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4'],
    'IgA_Memory': ['IGHA1', 'IGHA2'],
    'IgE_Memory': ['IGHE'],  # Critical for allergic diseases
    
    # ========== PLASMA CELLS / PLASMABLASTS ==========
    'Plasma_Cells': ['SDC1', 'PRDM1', 'XBP1', 'IRF4', 'MZB1', 'SSR4', 'DERL3'],  # CD138 = SDC1
    'Plasmablasts': ['CD27', 'CD38', 'PRDM1', 'XBP1'],  # CD27+CD38+
    
    # Plasma Cell Immunoglobulin Production
    'IgM_Plasma': ['IGHM', 'PRDM1'],
    'IgG_Plasma': ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'PRDM1'],
    'IgA_Plasma': ['IGHA1', 'IGHA2', 'PRDM1'],
    'IgE_Plasma': ['IGHE', 'PRDM1'],  # Local IgE production in allergic inflammation
    
    # ========== GERMINAL CENTER B CELLS ==========
    'GC_B_Cells': ['BCL6', 'AICDA', 'MME', 'STMN1', 'SUGCT'],  # CD10 = MME
    'Light_Zone_GC': ['CD83', 'EBI3', 'NFKB1'],  # Positive selection
    'Dark_Zone_GC': ['CXCR4', 'AICDA', 'MKI67'],  # Proliferation and SHM
    
    # ========== ACTIVATED B CELLS ==========
    'Activated_B': ['CD69', 'CD86', 'CD40', 'TNFRSF17', 'TNFRSF13B'],  # BCMA, TACI
    
    # ========== REGULATORY B CELLS (Bregs) ==========
    'Regulatory_B': ['IL10', 'CD24', 'CD38', 'TGFB1'],  # IL-10 producing B cells
    
    # ========== B CELL RECEPTOR (BCR) SIGNALING ==========
    'BCR_Signaling': ['CD22', 'CD72', 'FCRL4', 'FCRL5', 'CR2'],  # CD21 = CR2
    
    # ========== CO-STIMULATORY / INHIBITORY ==========
    'Co_Stimulation': ['CD40', 'CD80', 'CD86', 'ICOS', 'ICOSL'],
    'Inhibitory': ['PDCD1', 'FCRL4', 'FCRL5', 'CD22'],
    
    # ========== CHEMOKINE RECEPTORS ==========
    'Chemokine_Receptors': ['CXCR4', 'CXCR5', 'CCR6', 'CCR7'],
    
    # ========== TYPE 2 INFLAMMATION RESPONSE ==========
    'Type2_Response': ['IL4R', 'IL13RA1', 'FCER2', 'STAT6'],  # CD23 = FCER2
    
    # ========== TRANSCRIPTION FACTORS ==========
    'TFs_B_Identity': ['PAX5', 'EBF1', 'IRF4', 'IRF8', 'BCL6', 'PRDM1'],
    
    # ========== PROLIFERATION ==========
    'Proliferation': ['MKI67', 'TOP2A', 'STMN1', 'PCNA', 'CDK1'],
    
    # ========== APOPTOSIS / SURVIVAL ==========
    'Survival': ['BCL2', 'MCL1', 'TNFRSF13B', 'TNFRSF13C', 'TNFRSF17'],  # TACI, BAFFR, BCMA
    
    # ========== CYTOKINES / CYTOKINE RECEPTORS ==========
    'Cytokines': ['IL10', 'IL6', 'TNFA', 'LTA', 'LTB'],  # LTA/LTB = lymphotoxin
    'Cytokine_Receptors': ['IL21R', 'IL4R', 'IL2RA', 'IL6R'],
}

# Flatten to list
ALL_BCELL_MARKERS = []
for category, genes in BCELL_CORE_MARKERS.items():
    ALL_BCELL_MARKERS.extend(genes)
ALL_BCELL_MARKERS = list(set(ALL_BCELL_MARKERS))  # Remove duplicates

print(f"  B cell markers: {len(ALL_BCELL_MARKERS)} unique genes")

# Column names
BATCH_KEY = 'sample'  # â­ UPDATE if your batch column has different name
CELLTYPIST_LABEL_KEY = 'cell_type_celltypist'
SCANVI_LABEL_KEY = 'cell_type_scanvi'

# Random seed
RANDOM_SEED = 42

print("="*80)
print("B CELL PIPELINE v3.5.1")
print("="*80)
print(f"\nCell Types: Naive B, Memory B, Plasma Cells, GC B, Activated B, Bregs")
print(f"\nConfiguration:")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  HVG: {N_HVG}")
print(f"  Latent dims: {N_LATENT}")
print(f"  Batch key: {BATCH_KEY}")
print(f"  Rare type threshold: {MIN_CELLS_PER_TYPE} cells")
print(f"  B cell markers: {len(ALL_BCELL_MARKERS)} genes")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def merge_rare_types(s: pd.Series, min_cells: int = 200, other: str = "Unknown") -> pd.Series:
    """
    Merge rare cell types (< min_cells) into 'other' category.
    
    This improves model stability and visualization clarity.
    """
    s = s.astype(str).copy()
    vc = s.value_counts()
    rare = vc[vc < min_cells].index
    
    if len(rare) > 0:
        print(f"   Merging {len(rare)} rare types (<{min_cells} cells) â†’ {other}")
        for rt in rare[:5]:  # Show first 5
            print(f"      {rt}: {vc[rt]} cells")
        if len(rare) > 5:
            print(f"      ... and {len(rare)-5} more")
        s.loc[s.isin(rare)] = other
    
    return s


def save_celltype_counts(counts_dict: dict, output_dir: Path) -> None:
    """Save cell type count statistics to CSV files"""
    for name, series in counts_dict.items():
        df = pd.DataFrame({
            'cell_type': series.index,
            'count': series.values,
            'percentage': (series.values / series.sum() * 100).round(2)
        })
        df = df.sort_values('count', ascending=False)
        filepath = output_dir / f"{name}.csv"
        df.to_csv(filepath, index=False)
        print(f"   âœ“ Saved: {filepath.name}")


# Set random seeds
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

scvi.settings.seed = RANDOM_SEED

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

# GPU detection
print(f"\n{'='*80}")
print("GPU CHECK")
print("="*80)

GPU_AVAILABLE = torch.cuda.is_available()
print(f"GPU available: {GPU_AVAILABLE}")

if GPU_AVAILABLE:
    gpu_name = torch.cuda.get_device_name(0)
    gpu_memory = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"GPU: {gpu_name}")
    print(f"Memory: {gpu_memory:.1f} GB")
    accelerator = 'gpu'
    devices = 1
else:
    print("âš ï¸  Running on CPU (training will be slow)")
    accelerator = 'cpu'
    devices = 'auto'

# Create output directories
OUTPUT_DIR = Path(OUTPUT_DIR)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FIG_DIR = OUTPUT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

MODEL_DIR = OUTPUT_DIR / "models"
MODEL_DIR.mkdir(exist_ok=True)

print(f"\nâœ“ Output directories created")

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.n_jobs = 48
sc.settings.figdir = FIG_DIR
sc.set_figure_params(dpi=300, facecolor='white', format='pdf')

# ============================================================================
# STAGE 0: DATA LOADING
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 0: DATA LOADING")
print("="*80)

print(f"\nLoading: {INPUT_FILE}")
adata = sc.read_h5ad(INPUT_FILE)

print(f"âœ“ Data loaded:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# Check batch distribution
if BATCH_KEY in adata.obs.columns:
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    print(f"  Batches: {len(batch_counts)}")
else:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found")

# Check data structure
print(f"\nData structure:")
print(f"  Layers: {list(adata.layers.keys())}")
print(f"  .raw: {adata.raw is not None}")

# Ensure counts layer exists
if 'counts' not in adata.layers:
    if adata.raw is not None:
        adata.layers['counts'] = adata.raw.X.copy()
        print(f"  ✓ Created counts layer from .raw")
    elif EXTERNAL_COUNTS_FILE is not None:
        # ⭐ NEW in v3.5.2: Recover from external file
        adata = recover_counts_from_external(
            adata,
            EXTERNAL_COUNTS_FILE,
            counts_layer_name='counts',
            verbose=True
        )
    else:
        raise ValueError(
            "No counts layer found. Options:\n"
            "1. Provide .raw in current file\n"
            "2. Set EXTERNAL_COUNTS_FILE to h5ad with counts"
        )

# ============================================================================
# GENE NAME NORMALIZATION (BEFORE ANY SUBSETTING)
# ============================================================================

print(f"\n{'='*80}")
print("GENE NAME NORMALIZATION")
print("="*80)

def normalize_gene_names_local_first(adata):
    """
    Robust gene name normalization: local columns > mygene > var_names
    Adds 'symbol_base' column to adata.var
    """
    # Try local columns first
    symbol_cols = ['gene_symbol', 'gene_symbols', 'symbol', 'feature_name']
    for col in symbol_cols:
        if col in adata.var.columns:
            vals = adata.var[col].astype(str)
            non_empty = (vals.str.len() > 0).sum() / len(vals)
            if non_empty > 0.5:
                adata.var['symbol_base'] = vals.values
                print(f"âœ“ Using local column: '{col}'")
                return
    
    # Try mygene if ENSEMBL detected
    first_gene = str(adata.var_names[0])
    if MYGENE_AVAILABLE and first_gene.startswith('ENSG'):
        print("Attempting mygene conversion (ENSEMBL detected)...")
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
            
            converted = [symbol_map.get(g, g) for g in adata.var_names]
            frac_converted = np.mean([not str(x).startswith('ENSG') for x in converted])
            
            if frac_converted > 0.3:
                adata.var['symbol_base'] = np.array(converted, dtype=object)
                print(f"âœ“ mygene: {int(frac_converted*len(converted))} genes converted")
                return
        except Exception as e:
            print(f"âš ï¸  mygene failed: {str(e)[:80]}")
    
    # Fallback: use var_names
    adata.var['symbol_base'] = adata.var_names.astype(str).values
    print("âœ“ Using var_names as symbol_base")

normalize_gene_names_local_first(adata)

# â­ CRITICAL: Preserve full genes to .raw BEFORE HVG selection
print(f"\nPreserving full genes to .raw...")
if adata.raw is None or adata.raw.n_vars < adata.n_vars:
    adata.raw = sc.AnnData(
        X=adata.layers['counts'],  # Shared memory reference
        obs=adata.obs.copy(),
        var=adata.var.copy()  # Includes symbol_base
    )
    print(f"âœ“ Created adata.raw: {adata.raw.n_vars:,} genes")
else:
    # Patch symbol_base into existing raw
    if 'symbol_base' in adata.var.columns and 'symbol_base' not in adata.raw.var.columns:
        common = adata.raw.var_names.intersection(adata.var_names)
        if len(common) > 0:
            adata.raw.var.loc[common, 'symbol_base'] = adata.var.loc[common, 'symbol_base'].values
            print(f"âœ“ Patched symbol_base into adata.raw ({len(common):,} genes)")

# ============================================================================
# STAGE 1: HVG SELECTION (NO TRAINING YET)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 1: HVG SELECTION")
print("="*80)

hvg_file = MODEL_DIR / "hvg_genes.txt"
scvi_model_dir = MODEL_DIR / "scvi_model"

if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
    # â­ P1-FIX: Use var_names (unique) for pretrained model, not symbol_base
    print(f"\nâ­ PRETRAINED MODE: Loading HVG list (var_names)...")
    
    pretrained_path = Path(PRETRAINED_SCVI_MODEL)
    gene_list_candidates = [
        pretrained_path.parent / "hvg_genes.txt",
        pretrained_path.parent / "genes_used.txt",
        pretrained_path / "genes.txt"
    ]
    
    genes_for_training = None
    for gf in gene_list_candidates:
        if gf.exists():
            genes_for_training = pd.read_csv(gf, header=None)[0].astype(str).tolist()
            print(f"  âœ“ Loaded: {gf} ({len(genes_for_training)} genes)")
            break
    
    if genes_for_training is None:
        raise FileNotFoundError(f"Cannot find gene list for pretrained model")
    
    # â­ CRITICAL: Match by var_names (unique), not symbol_base
    genes_in_data = [g for g in genes_for_training if g in adata.var_names]
    missing = len(genes_for_training) - len(genes_in_data)
    
    if missing > 0:
        print(f"  âš ï¸  {missing}/{len(genes_for_training)} genes not in data")
        if missing / len(genes_for_training) > 0.2:
            raise ValueError("Too many missing genes for pretrained model")
    
    hvg_mask = adata.var_names.isin(genes_in_data).to_numpy()
    hvg_method = "pretrained_var_names"
    
else:
    # â­ TRAINING MODE: Select HVG
    print(f"\nâ­ TRAINING MODE: Selecting {N_HVG} HVGs...")
    
    try:
        print(f"   Attempting batch-aware HVG...")
        sc.pp.highly_variable_genes(
            adata,
            layer='counts',
            n_top_genes=N_HVG,
            batch_key=BATCH_KEY,
            subset=False,
            flavor='seurat_v3'
        )
        hvg_method = "batch-aware"
        print(f"   âœ“ Success")
    except Exception as e:
        print(f"   âš ï¸  Failed: {str(e)[:80]}")
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
    n_hvg = int(hvg_mask.sum())
    print(f"âœ“ Selected {n_hvg:,} HVGs ({hvg_method})")
    
    # â­ P1-FIX: Save var_names (not symbol_base) for model reuse
    genes_to_save = adata.var_names[hvg_mask].astype(str)
    pd.Series(genes_to_save.values).to_csv(hvg_file, index=False, header=False)
    print(f"âœ“ Gene list saved (var_names): {hvg_file}")

# ============================================================================
# FORCE INCLUDE B CELL MARKERS IN HVG
# ============================================================================

print(f"\n{'='*80}")
print("FORCE INCLUDE B CELL MARKERS IN HVG")
print("="*80)

# Get gene symbol column
gene_symbols = None
symbol_col = None
for col in ['symbol_base', 'gene_symbol', 'gene_symbols', 'feature_name']:
    if col in adata.var.columns:
        gene_symbols = adata.var[col].astype(str)
        symbol_col = col
        break

if gene_symbols is None:
    gene_symbols = adata.var_names.astype(str)
    symbol_col = 'var_names'

print(f"\nUsing gene symbols from: {symbol_col}")
print(f"Total markers to check: {len(ALL_BCELL_MARKERS)}")

# Ensure 'highly_variable' column exists
if 'highly_variable' not in adata.var.columns:
    print(f"âš ï¸  'highly_variable' column missing, creating from hvg_mask...")
    adata.var['highly_variable'] = hvg_mask

# Check each marker
print(f"\n{'='*60}")
print("Checking B cell markers:")
print("="*60)

missing_in_data = []
missing_in_hvg = []
present_in_hvg = []

for category, genes in BCELL_CORE_MARKERS.items():
    print(f"\n{category}:")
    for gene in genes:
        # Find gene (case-insensitive)
        matches = gene_symbols.str.upper() == gene.upper()
        
        if matches.sum() == 0:
            print(f"  âœ— {gene}: NOT IN DATA")
            missing_in_data.append(gene)
        else:
            # Get first match index
            idx = matches[matches].index[0]
            is_hvg = adata.var.loc[idx, 'highly_variable']
            
            if is_hvg:
                print(f"  âœ“ {gene}: in HVG")
                present_in_hvg.append(gene)
            else:
                print(f"  âš ï¸  {gene}: NOT in HVG â†’ adding")
                missing_in_hvg.append(gene)

# Summary
print(f"\n{'='*60}")
print(f"SUMMARY:")
print(f"  Total markers: {len(ALL_BCELL_MARKERS)}")
print(f"  Present in data: {len(present_in_hvg) + len(missing_in_hvg)}")
print(f"  Missing in data: {len(missing_in_data)}")
print(f"  Already in HVG: {len(present_in_hvg)}")
print(f"  TO ADD to HVG: {len(missing_in_hvg)}")
print("="*60)

if len(missing_in_data) > 0:
    print(f"\nâš ï¸  Markers missing in data (will skip):")
    for gene in missing_in_data[:10]:
        print(f"     {gene}")
    if len(missing_in_data) > 10:
        print(f"     ... and {len(missing_in_data)-10} more")

# Force add missing markers to HVG
if len(missing_in_hvg) > 0:
    print(f"\nðŸ”§ Adding {len(missing_in_hvg)} markers to HVG...")
    
    n_added = 0
    for gene in missing_in_hvg:
        matches = gene_symbols.str.upper() == gene.upper()
        if matches.sum() > 0:
            idx = matches[matches].index[0]
            adata.var.loc[idx, 'highly_variable'] = True
            hvg_mask[adata.var_names.get_loc(idx)] = True  # Update mask
            n_added += 1
            print(f"   âœ“ Added {gene}")
    
    new_hvg_count = adata.var['highly_variable'].sum()
    print(f"\nâœ… HVG updated: {N_HVG} â†’ {new_hvg_count} (+{n_added})")
    
    # Update saved gene list
    genes_to_save = adata.var_names[hvg_mask].astype(str)
    pd.Series(genes_to_save.values).to_csv(hvg_file, index=False, header=False)
    print(f"âœ… Gene list re-saved with markers: {hvg_file}")
else:
    print(f"\nâœ… All required markers already in HVG")

# Verify critical markers
critical_markers = ['CD19', 'CD79A', 'MS4A1', 'IGHM', 'IGHG1', 'PRDM1']
print(f"\n{'='*60}")
print(f"VERIFY CRITICAL MARKERS:")
print("="*60)

all_present = True
for gene in critical_markers:
    matches = gene_symbols.str.upper() == gene.upper()
    if matches.sum() > 0:
        idx = matches[matches].index[0]
        is_hvg = adata.var.loc[idx, 'highly_variable']
        status = "âœ… IN HVG" if is_hvg else "âŒ NOT IN HVG"
        print(f"  {gene}: {status}")
        if not is_hvg:
            all_present = False
    else:
        print(f"  {gene}: âŒ NOT IN DATA")
        all_present = False

if not all_present:
    print(f"\nâš ï¸  WARNING: Some critical markers missing!")
    print(f"   B cell classification may fail in downstream analysis")

print("="*60)
print(f"\nâœ… B cell marker check complete")
print(f"   Final HVG count: {adata.var['highly_variable'].sum()}")
print("="*60)

print(f"\nâœ“ Main adata preserved: {adata.n_vars:,} genes (full)")
print(f"âœ“ HVG mask created: {int(hvg_mask.sum()):,} genes for training")

# ============================================================================
# STAGE 2: CellTypist ANNOTATION (ON FULL GENES)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 2: CellTypist ANNOTATION")
print("="*80)

# Load model
print(f"\nLoading CellTypist model: {CELLTYPIST_MODEL}")
try:
    model = models.Model.load(model=CELLTYPIST_MODEL)
    print(f"âœ“ Model loaded")
except:
    print(f"  Downloading...")
    models.download_models(model=CELLTYPIST_MODEL)
    model = models.Model.load(model=CELLTYPIST_MODEL)

# â­ P1-FIX: Hard validation of CellTypist features
print(f"\nâ­ Validating CellTypist model features...")
model_features = None
for attr in ['features', 'genes', 'var_names']:
    if hasattr(model, attr):
        mf = getattr(model, attr)
        if mf is not None and len(mf) > 1000:
            model_features = pd.Index(mf).astype(str)
            print(f"  âœ“ Found {len(model_features)} features via '{attr}'")
            break

if model_features is None:
    raise ValueError(
        "Cannot retrieve CellTypist model features.\n"
        "Expected attributes: features/genes/var_names with >1000 genes."
    )

# Build CellTypist input with proper feature matching
print(f"\nâ­ Building CellTypist input with base symbol matching...")

# Use base symbols (no -1/-2 suffixes) from .raw
sym = pd.Index(adata.raw.var['symbol_base']).astype(str)

# Find overlap: genes in model AND not duplicated in our data
keep_mask = sym.isin(model_features) & ~sym.duplicated(keep='first')
n_overlap = keep_mask.sum()

print(f"  Our data: {len(sym)} genes")
print(f"  Overlap: {n_overlap} genes")
print(f"  Coverage: {n_overlap/len(model_features)*100:.1f}%")

if n_overlap < 1000:
    print(f"\n  âš ï¸  WARNING: Low overlap ({n_overlap} genes)")
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

print(f"âœ“ CellTypist input: {adata_celltypist.n_obs:,} cells Ã— {adata_celltypist.n_vars:,} genes")

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
print(f"âœ“ Prediction complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# â­ P1-FIX: Properly use majority voting if enabled
print(f"\nâ­ Writing results with index alignment...")
pred_df = predictions.predicted_labels
pred_df = pred_df.reindex(adata.obs_names)  # Force alignment

# Choose column: majority_voting if available, else predicted_labels
if CELLTYPIST_MAJORITY_VOTING and 'majority_voting' in pred_df.columns:
    use_col = 'majority_voting'
    print(f"  âœ“ Using majority_voting for scANVI labels")
else:
    use_col = 'predicted_labels'
    print(f"  âœ“ Using predicted_labels for scANVI labels")

adata.obs['cell_type_celltypist_raw'] = pred_df[use_col].astype(str).values

# Also save both versions for reference
if 'predicted_labels' in pred_df.columns:
    adata.obs['celltypist_predicted_labels'] = pred_df['predicted_labels'].astype(str).values
if 'majority_voting' in pred_df.columns:
    adata.obs['celltypist_majority_voting'] = pred_df['majority_voting'].astype(str).values

# Confidence detection
conf_col = None
for name in ['conf_score', 'confidence', 'confidence_score', 'prob']:
    if name in pred_df.columns:
        conf_col = name
        break

if conf_col:
    adata.obs['celltypist_confidence'] = pred_df[conf_col].values
    print(f"  âœ“ Confidence: {conf_col}")
else:
    adata.obs['celltypist_confidence'] = 1.0
    print(f"  âš ï¸  No confidence, using 1.0")

# Clean up
del adata_celltypist, X_ct, var_ct
gc.collect()

# Summary
celltypist_counts_raw = adata.obs['cell_type_celltypist_raw'].value_counts()
mean_conf = adata.obs['celltypist_confidence'].mean()

print(f"\nâœ“ CellTypist raw results:")
print(f"  Unique types: {len(celltypist_counts_raw)}")
print(f"  Mean confidence: {mean_conf:.3f}")

# Save raw counts
print(f"\nSaving CellTypist raw counts...")
save_celltype_counts(
    {'celltypist_counts_raw': celltypist_counts_raw},
    OUTPUT_DIR
)

print(f"\nâœ“ Stage 2 complete")

# ============================================================================
# STAGE 3: PREPARE LABELS FOR scANVI (WITH RARE TYPE FILTERING)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 3: PREPARE LABELS FOR scANVI")
print("="*80)

print(f"\nâ­ Creating scANVI training labels from CellTypist...")
print(f"   Low confidence threshold: {CELLTYPIST_MIN_CONFIDENCE}")
print(f"   Rare type threshold: {MIN_CELLS_PER_TYPE} cells")

# Start with CellTypist predictions (already using majority_voting if enabled)
ct_raw = adata.obs['cell_type_celltypist_raw'].astype(str)

# Step 1: Low confidence â†’ Unknown
print(f"\nStep 1: Handling low confidence cells...")
ct_lc = ct_raw.where(
    adata.obs['celltypist_confidence'] >= CELLTYPIST_MIN_CONFIDENCE,
    'Unknown'
)
n_low_conf = (ct_lc == 'Unknown').sum()
print(f"   Low confidence â†’ Unknown: {n_low_conf:,} cells")

# Step 2: Rare types â†’ Unknown
print(f"\nStep 2: Merging rare types...")
ct_filt = merge_rare_types(
    ct_lc,
    min_cells=MIN_CELLS_PER_TYPE,
    other='Unknown'
)

# Save filtered version
adata.obs['cell_type_celltypist_filt'] = ct_filt

# Convert to categorical for scANVI
labels = ct_filt.astype('category')
if 'Unknown' not in labels.cat.categories:
    labels = labels.cat.add_categories(['Unknown'])

adata.obs['labels_for_scanvi'] = labels

# Statistics
n_unknown = int((labels == 'Unknown').sum())
n_labeled = int(adata.n_obs - n_unknown)
n_unique_types = int(labels.nunique()) - (1 if 'Unknown' in labels.cat.categories else 0)

print(f"\nâœ“ scANVI labels prepared:")
print(f"   Labeled cells: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"   Unknown cells: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(f"   Unique types (excl. Unknown): {n_unique_types}")

# Save filtered counts
print(f"\nSaving CellTypist filtered counts...")
save_celltype_counts(
    {'celltypist_counts_filt': adata.obs['cell_type_celltypist_filt'].value_counts()},
    OUTPUT_DIR
)

# ============================================================================
# STAGE 4: BUILD adata_model FOR scVI/scANVI
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 4: BUILD adata_model (HVG-ONLY FOR TRAINING)")
print("="*80)

print(f"\nâ­â­â­ P0-FIX: Creating dedicated adata_model for scVI/scANVI...")
print(f"   This ensures labels and batch are in the training object")

# â­ P2-FIX: Minimal construction (only necessary columns)
X_hvg = adata.layers['counts'][:, hvg_mask].copy()
obs_hvg = adata.obs[[BATCH_KEY, 'labels_for_scanvi']].copy()
var_hvg = adata.var.loc[hvg_mask, []].copy()  # Minimal var (only index)

adata_model = sc.AnnData(X=X_hvg, obs=obs_hvg, var=var_hvg)
adata_model.layers['counts'] = adata_model.X.copy()

print(f"\nâœ“ adata_model created:")
print(f"   Shape: {adata_model.n_obs:,} cells Ã— {adata_model.n_vars:,} genes")
print(f"   obs columns: {list(adata_model.obs.columns)}")
print(f"   .X: raw counts (for scVI setup)")
print(f"   layers['counts']: raw counts (explicit)")

# Clean up
del X_hvg, obs_hvg, var_hvg
gc.collect()

# ============================================================================
# STAGE 5: scVI BATCH CORRECTION (ON adata_model)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 5: scVI BATCH CORRECTION")
print("="*80)

# Setup for scVI
print(f"\nSetting up scVI on adata_model...")
scvi.model.SCVI.setup_anndata(
    adata_model,
    layer='counts',
    batch_key=BATCH_KEY
)
print(f"âœ“ scVI data setup complete")

# Train or load model
if PRETRAINED_SCVI_MODEL and Path(PRETRAINED_SCVI_MODEL).exists():
    print(f"\nâ­ Loading pretrained scVI: {PRETRAINED_SCVI_MODEL}")
    vae = scvi.model.SCVI.load(str(PRETRAINED_SCVI_MODEL), adata=adata_model)
    print(f"âœ“ Model loaded")
    
else:
    print(f"\nâ­ Training scVI model...")
    vae = scvi.model.SCVI(
        adata_model,
        n_latent=N_LATENT,
        n_layers=N_LAYERS,
        gene_likelihood="nb"
    )
    
    print(f"  Parameters: {sum(p.numel() for p in vae.module.parameters()):,}")
    
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
    print(f"\nâœ“ Training complete: {int(elapsed//60)}m {int(elapsed%60)}s")
    
    # Save
    print(f"\nSaving scVI model: {scvi_model_dir}")
    vae.save(scvi_model_dir, overwrite=True)

# â­ P0-FIX: Index-aligned writeback using pandas reindex
print(f"\nâ­ Writing X_scvi to main adata (index-aligned)...")
latent_scvi = pd.DataFrame(
    vae.get_latent_representation(),
    index=adata_model.obs_names
)
adata.obsm['X_scvi'] = latent_scvi.reindex(adata.obs_names).to_numpy()
print(f"âœ“ X_scvi: {adata.obsm['X_scvi'].shape}")

print(f"\nâœ“ Stage 5 complete")

# ============================================================================
# STAGE 6: scANVI FINE-TUNING (ON adata_model)
# ============================================================================

print(f"\n{'='*80}")
print("STAGE 6: scANVI FINE-TUNING")
print("="*80)

print(f"\nâ­ P0-FIX: Initializing scANVI from scVI with explicit adata...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata_model,  # â­ CRITICAL: Pass adata_model explicitly
    labels_key='labels_for_scanvi',
    unlabeled_category='Unknown'
)
print(f"âœ“ scANVI created")

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
print(f"\nâœ“ Training complete: {int(elapsed//60)}m {int(elapsed%60)}s")

# Save
scanvi_model_dir = MODEL_DIR / "scanvi_model"
print(f"\nSaving scANVI model: {scanvi_model_dir}")
lvae.save(scanvi_model_dir, overwrite=True)

# â­ P0-FIX: Index-aligned predictions using pandas reindex
print(f"\nâ­ Generating scANVI predictions (index-aligned)...")

pred = pd.Series(lvae.predict(), index=adata_model.obs_names)
adata.obs['cell_type_scanvi_raw'] = pred.reindex(adata.obs_names).astype(str).values

# â­ P0-FIX: predict(soft=True) returns ndarray, use np.asarray()
probs = np.asarray(lvae.predict(soft=True))  # ndarray
conf = pd.Series(probs.max(axis=1), index=adata_model.obs_names)
adata.obs['scanvi_confidence'] = conf.reindex(adata.obs_names).values

z = pd.DataFrame(lvae.get_latent_representation(), index=adata_model.obs_names)
adata.obsm['X_scanvi'] = z.reindex(adata.obs_names).to_numpy()

# Filter rare types in scANVI predictions
print(f"\nâ­ Filtering rare types in scANVI predictions...")
adata.obs['cell_type_scanvi_filt'] = merge_rare_types(
    adata.obs['cell_type_scanvi_raw'],
    min_cells=MIN_CELLS_PER_TYPE,
    other='Unknown'
).astype('category')

print(f"âœ“ scANVI predictions complete")
print(f"   X_scanvi: {adata.obsm['X_scanvi'].shape}")
print(f"   Unique types (raw): {adata.obs['cell_type_scanvi_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_scanvi_filt'].nunique()}")

# Save scANVI counts
print(f"\nSaving scANVI counts...")
save_celltype_counts(
    {
        'scanvi_counts_raw': adata.obs['cell_type_scanvi_raw'].value_counts(),
        'scanvi_counts_filt': adata.obs['cell_type_scanvi_filt'].value_counts()
    },
    OUTPUT_DIR
)

# Clean up adata_model
del adata_model
gc.collect()

print(f"\nâœ“ Stage 6 complete")

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
final_counts = adata.obs['cell_type_scanvi_filt'].value_counts()
print(f"\n2. FINAL DISTRIBUTION (filtered)")
print(f"  Cell types: {len(final_counts)}")

# Agreement (using filtered versions)
ct_filt = adata.obs['cell_type_celltypist_filt'].astype('object')
scanvi_filt = adata.obs['cell_type_scanvi_filt'].astype('object')
agreement_mask = (
    (ct_filt == scanvi_filt)
    | (pd.isna(ct_filt) & pd.isna(scanvi_filt))
)
agreements = int(agreement_mask.sum())
agreement_rate = agreements / adata.n_obs * 100

print(f"\n3. AGREEMENT (filtered)")
print(f"  CellTypist-scANVI: {agreement_rate:.1f}%")

# ============================================================================
# VISUALIZATION: DUAL UMAPs (scVI + scANVI)
# ============================================================================

print(f"\n{'='*80}")
print("VISUALIZATION: DUAL UMAPs")
print("="*80)

# â­ P2-FIX: No full X copy - use layer operations
print(f"\nâ­ P2-FIX: Creating log1p layer (no full X copy)...")
if 'log1p' not in adata.layers:
    # Create log1p in-place on counts layer (view, not copy)
    from scipy.sparse import issparse
    if issparse(adata.layers['counts']):
        import scipy.sparse as sp
        # Normalize and log1p in layer
        counts_norm = adata.layers['counts'].copy()
        sc.pp.normalize_total(sc.AnnData(X=counts_norm), target_sum=1e4, inplace=True)
        sc.pp.log1p(sc.AnnData(X=counts_norm))
        adata.layers['log1p'] = counts_norm
        del counts_norm
    else:
        adata.layers['log1p'] = adata.layers['counts'].copy()
        sc.pp.normalize_total(sc.AnnData(X=adata.layers['log1p']), target_sum=1e4, inplace=True)
        sc.pp.log1p(sc.AnnData(X=adata.layers['log1p']))
    
    print(f"âœ“ layers['log1p'] created")

# Set X for UMAP (scanpy requires it)
adata.X = adata.layers['log1p']

# UMAP 1: scVI space
print(f"\nâ­ Computing scVI UMAP (with neighbors_key)...")
sc.pp.neighbors(
    adata,
    use_rep='X_scvi',
    n_neighbors=15,
    key_added='neighbors_scvi'
)
sc.tl.umap(adata, neighbors_key='neighbors_scvi')
adata.obsm['X_umap_scvi'] = adata.obsm['X_umap'].copy()
print(f"âœ“ UMAP(scVI) stored in X_umap_scvi: {adata.obsm['X_umap_scvi'].shape}")

# UMAP 2: scANVI space
print(f"\nâ­ Computing scANVI UMAP (with neighbors_key)...")
sc.pp.neighbors(
    adata,
    use_rep='X_scanvi',
    n_neighbors=15,
    key_added='neighbors_scanvi'
)
sc.tl.umap(adata, neighbors_key='neighbors_scanvi')
adata.obsm['X_umap_scanvi'] = adata.obsm['X_umap'].copy()
print(f"âœ“ UMAP(scANVI) stored in X_umap_scanvi: {adata.obsm['X_umap_scanvi'].shape}")

# Figure A: scVI space
print(f"\nGenerating scVI comparison plot...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='cell_type_celltypist_filt',
    ax=axes[0],
    show=False,
    title='CellTypist (filtered, scVI)',
    size=2,
    legend_loc='right margin',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color='celltypist_confidence',
    ax=axes[1],
    show=False,
    title='CellTypist Confidence (scVI)',
    size=2,
    cmap='viridis',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scvi',
    color=BATCH_KEY,
    ax=axes[2],
    show=False,
    title='Batch (scVI)',
    size=1,
    frameon=False
)

plt.tight_layout()
plt.savefig(FIG_DIR / 'comparison_scvi.pdf', dpi=300, bbox_inches='tight')
plt.close()
print(f"âœ“ Saved: comparison_scvi.pdf")

# Figure B: scANVI space
print(f"\nGenerating scANVI comparison plot...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='cell_type_scanvi_filt',
    ax=axes[0],
    show=False,
    title='scANVI Final (filtered)',
    size=2,
    legend_loc='right margin',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color='scanvi_confidence',
    ax=axes[1],
    show=False,
    title='scANVI Confidence',
    size=2,
    cmap='viridis',
    frameon=False
)

sc.pl.embedding(
    adata,
    basis='umap_scanvi',
    color=BATCH_KEY,
    ax=axes[2],
    show=False,
    title='Batch (scANVI)',
    size=1,
    frameon=False
)

plt.tight_layout()
plt.savefig(FIG_DIR / 'comparison_scanvi.pdf', dpi=300, bbox_inches='tight')
plt.close()
print(f"âœ“ Saved: comparison_scanvi.pdf")

# ============================================================================
# SAVE FINAL RESULTS
# ============================================================================

print(f"\n{'='*80}")
print("SAVE FINAL RESULTS")
print("="*80)

output_file = OUTPUT_DIR / "adata_bcell_FINAL.h5ad"

# Add metadata
adata.uns['pipeline_info'] = {
    'version': '3.5.1-BCELL',
    'date': '2025-01-11',
    'cell_types': 'naive_B, memory_B, plasma_cells, GC_B, activated_B, Bregs',
    'input_file': str(INPUT_FILE),
    'batch_key': BATCH_KEY,
    'n_hvg': int(hvg_mask.sum()),
    'hvg_method': hvg_method,
    'celltypist_model': str(CELLTYPIST_MODEL),
    'celltypist_confidence_threshold': CELLTYPIST_MIN_CONFIDENCE,
    'min_cells_per_type': MIN_CELLS_PER_TYPE,
    'bcell_markers_count': len(ALL_BCELL_MARKERS),
    'critical_fixes_v3_5_1': [
        'P0: scANVI uses dedicated adata_model (labels alignment)',
        'P0: predict(soft=True).max() changed to np.asarray()',
        'P1: CellTypist majority voting properly integrated',
        'P1: CellTypist features validation (hard check)',
        'P1: Pretrained HVG uses var_names (unique)',
        'P2: Memory optimization (layers, no full X copy)',
        'P2: Minimal adata_model construction'
    ],
    'bcell_marker_categories': list(BCELL_CORE_MARKERS.keys())
}

adata.write_h5ad(output_file, compression='gzip')
size_gb = output_file.stat().st_size / 1e9
print(f"\nâœ“ Saved: {output_file}")
print(f"   Size: {size_gb:.2f} GB")

# Export annotations
annotations = adata.obs[[
    'cell_type_celltypist_raw',
    'cell_type_celltypist_filt',
    'celltypist_confidence',
    'cell_type_scanvi_raw',
    'cell_type_scanvi_filt',
    'scanvi_confidence',
    BATCH_KEY
]].copy()
annotations.to_csv(OUTPUT_DIR / "annotations_complete.csv")
print(f"âœ“ Annotations exported")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print("ðŸŽ‰ B CELL PIPELINE COMPLETE - v3.5.1")
print("="*80)

print(f"\nðŸ“Š Analysis Summary:")
print(f"   Cell types: Naive B, Memory B, Plasma Cells, GC B, Activated B, Bregs")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (full): {adata.n_vars:,}")
print(f"   Genes (raw): {adata.raw.n_vars:,}")
print(f"   Genes (HVG): {int(hvg_mask.sum())}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")
print(f"   B cell markers forced: {len(ALL_BCELL_MARKERS)}")

print(f"\nâ­ CellTypist Results:")
print(f"   Unique types (raw): {adata.obs['cell_type_celltypist_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_celltypist_filt'].nunique()}")
print(f"   Mean confidence: {float(np.mean(adata.obs['celltypist_confidence'])):.3f}")

print(f"\nâ­ scANVI Results:")
print(f"   Unique types (raw): {adata.obs['cell_type_scanvi_raw'].nunique()}")
print(f"   Unique types (filtered): {adata.obs['cell_type_scanvi_filt'].nunique()}")
print(f"   Mean confidence: {float(np.mean(adata.obs['scanvi_confidence'])):.3f}")
print(f"   CellTypist-scANVI agreement (filtered): {agreement_rate:.1f}%")

print(f"\nðŸ”§ Critical Fixes Applied:")
print(f"   âœ… P0: Dedicated adata_model for scVI/scANVI")
print(f"   âœ… P0: Index-aligned writeback (pandas reindex)")
print(f"   âœ… P0: predict(soft=True) ndarray handling")
print(f"   âœ… P1: CellTypist majority voting integrated")
print(f"   âœ… P1: Model features validation")
print(f"   âœ… P1: Pretrained HVG uses var_names")
print(f"   âœ… P2: No full X copy (layers-based)")

print(f"\nðŸ”¬ B Cell-Specific Features:")
print(f"   âœ… {len(BCELL_CORE_MARKERS)} marker categories")
print(f"   âœ… Immunoglobulin isotypes (IgM, IgG, IgA, IgE)")
print(f"   âœ… Plasma cell maturation markers")
print(f"   âœ… Germinal center reaction (light/dark zone)")
print(f"   âœ… B cell activation and BCR signaling")
print(f"   âœ… Type 2 inflammation response (IgE relevant for CRSwNP)")

print(f"\nðŸ“ Output Files:")
print(f"   h5ad: {output_file.name}")
print(f"   Annotations: annotations_complete.csv")
print(f"   Cell type counts: celltypist_counts_*.csv, scanvi_counts_*.csv")
print(f"   Models: {scvi_model_dir}, {scanvi_model_dir}")
print(f"   Gene list: hvg_genes.txt (var_names)")

print(f"\nðŸ“ˆ Figures:")
print(f"   comparison_scvi.pdf - scVI UMAP with CellTypist")
print(f"   comparison_scanvi.pdf - scANVI final results")

print("\n" + "="*80)
print("NEXT STEPS:")
print("="*80)
print("1. Examine annotation quality in scANVI UMAPs")
print("2. Validate B cell markers (CD19, MS4A1, IGHM, IGHG)")
print("3. Check immunoglobulin isotype distribution")
print("4. Assess plasma cell frequency and PRDM1 expression")
print("5. Investigate IgE+ cells if studying allergic disease")
print("6. Analyze germinal center B cells if present")
print("7. Compare B cell activation states across conditions")
print("="*80)