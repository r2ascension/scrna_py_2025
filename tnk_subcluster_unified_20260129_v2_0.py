#!/usr/bin/env python3
"""
T/NK Cell Subcluster Unified Pipeline v2.0.1 (PRODUCTION - HOTFIX)
===================================================================
Author: r2end
Date: 2026-01-29
Version: 2.0.1-HOTFIX

Purpose:
--------
All-in-one pipeline that:
1. Performs subcluster analysis on each T/NK cell type
2. Standardizes annotations (L2/L3 format)
3. Merges all celltypes into unified format
4. Computes global UMAP from HVG/PCA
5. Generates comprehensive visualizations and tables

Key Features:
-------------
- Direct output in unified format (no converter needed)
- Memory-optimized workflow
- Production-grade error handling
- Compatible with B cell unified pipeline v2.0

Architecture:
-------------
For each celltype:
  → Subcluster (BBKNN + Leiden)
  → Standardize annotations (L2/L3)
  → Add celltype prefix to obs_names
→ Merge all celltypes
→ Global HVG + PCA + UMAP
→ Save unified h5ad + figures + tables

HOTFIX v2.0.1 Changes:
----------------------
P0 Fixes:
  - Fixed neighbors_within_batch calculation (P0 bug when min_batch_size=3)
  - Fixed categorical loss in leiden transfer
  - Fixed join='outer' unnecessary risk
  - Removed local X_umap_bbknn from merged object
  - Ensured categorical types after merge

P1 Fixes:
  - Added thread control for BLAS libraries (OMP/MKL)
  - Added neighbors/UMAP key isolation
  - Limited marker export to top 200 genes

P2 Fixes:
  - Used local min_cells_per_cluster_ct variable
  - Cleaned unused config

Input:
------
- Pre-processed h5ad with:
  * layers['counts', 'log1p']
  * .raw: Full gene set
  * obs['cell_type_scanvi_cd4cd8_filt']: Cell type annotations

Output:
-------
- results/
  ├── adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad ⭐
  ├── figures/
  │   ├── global_umap_summary.pdf ⭐
  │   ├── umap_level3_highres.pdf ⭐
  │   └── [celltype]_subclustering.pdf (per celltype)
  └── tables/
      ├── subcluster_summary.csv ⭐
      ├── celltype_L2_counts.csv
      ├── celltype_L3_counts.csv
      ├── final_annotations.csv
      └── [celltype]_markers.csv (with L2/L3 columns)

Memory: ~40-60GB total
Runtime: ~2-4 hours
"""

# ===== P1-7 FIX: Control BLAS threading BEFORE imports =====
# Prevent thread over-subscription from MKL/OpenBLAS/OMP
# Note: This limits BLAS libraries, NOT scanpy's n_jobs
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Fix matplotlib backend before any plotting imports
import matplotlib
matplotlib.use('Agg')

import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gc
from pathlib import Path
import warnings
import sys
from datetime import datetime
warnings.filterwarnings('ignore')

print("=" * 80)
print("T/NK CELL SUBCLUSTER UNIFIED PIPELINE v2.0 (ALL-IN-ONE)")
print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)
print("\nKey Features:")
print("  ✓ Subcluster analysis per celltype")
print("  ✓ Standardized L2/L3 annotations")
print("  ✓ Unified format output (no converter needed)")
print("  ✓ Global UMAP from HVG/PCA")
print("  ✓ Production-grade quality control")
print()

# ===== Dependency Check =====
print("🔍 Checking dependencies...")
try:
    import bbknn
    print("  ✓ bbknn available")
except ImportError:
    print("  ❌ bbknn not found. Install with: pip install bbknn")
    sys.exit(1)

try:
    if hasattr(sc.get, 'rank_genes_groups_df'):
        USE_MODERN_MARKER_EXPORT = True
        print("  ✓ scanpy.get.rank_genes_groups_df available")
    else:
        USE_MODERN_MARKER_EXPORT = False
        print("  ⚠️  Using legacy marker export")
except:
    USE_MODERN_MARKER_EXPORT = False
print()

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths
INPUT_PATH = "/home/h2048/data/py/0109/tcell_cd4cd8_scanvi_v1_3/adata_tcell_cd4cd8_FINAL.h5ad"
OUTPUT_BASE = Path("/home/h2048/data/py/0129/tnk_analysis_unified")
OUTPUT_DIR = OUTPUT_BASE / "results" / "subcluster_unified_v2_20260129"
FIG_DIR = OUTPUT_DIR / "figures"
TABLE_DIR = OUTPUT_DIR / "tables"

for d in [OUTPUT_DIR, FIG_DIR, TABLE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Column names
CELLTYPE_COLUMN = 'cell_type_scanvi_cd4cd8_filt'
BATCH_KEY = 'sample'

# Cell type selection
CELLTYPE_WHITELIST = []  # Empty = process all
CELLTYPE_BLACKLIST = []
MIN_CELLS_FOR_CELLTYPE = 100

# Quality control
MIN_CELLS_PER_BATCH = 3
MIN_CELLS_PER_CLUSTER = 10
MIN_CELLS_FOR_MARKER = 20

# Analysis parameters (T/NK optimized)
DEFAULT_N_HVG = 3000
DEFAULT_N_PCS = 50
DEFAULT_RESOLUTION = 0.3
DEFAULT_NEIGHBORS = 5
DEFAULT_HVG_FLAVOR = 'seurat'

# T/NK-specific parameters
TNK_PARAMS = {
    "CD8+ Trm cytotoxic T cells": {
        "n_hvg": 3000, "resolution": 0.4, "min_cells_per_cluster": 15,
        "known_markers": ["CD3D", "CD8A", "CD8B", "GZMB", "PRF1", "ITGAE"]
    },
    "CD8+ Tem/Trm cytotoxic T cells": {
        "n_hvg": 3000, "resolution": 0.4, "min_cells_per_cluster": 15,
        "known_markers": ["CD3D", "CD8A", "GZMK", "CD69", "ITGAE"]
    },
    "CD8+ Tem/Temra cytotoxic T cells": {
        "n_hvg": 3000, "resolution": 0.3, "min_cells_per_cluster": 15,
        "known_markers": ["CD3D", "CD8A", "GZMB", "CX3CR1"]
    },
    "CD4+ Tem/Effector helper T cells": {
        "n_hvg": 3000, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "CD4", "IL7R", "CCL5"]
    },
    "CD4+ Regulatory T cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4"]
    },
    "Regulatory T cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "FOXP3", "IL2RA", "CTLA4"]
    },
    "CD16+ NK cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["GNLY", "NKG7", "FCGR3A", "NCAM1"]
    },
    "NK cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["GNLY", "NKG7", "NCAM1", "KLRD1"]
    },
    "MAIT cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "SLC4A10", "KLRB1", "NCR3"]
    },
    "Follicular helper T cells": {
        "n_hvg": 2500, "resolution": 0.3, "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "CD4", "CXCR5", "PDCD1", "ICOS"]
    }
}

# Global settings
RANDOM_STATE = 42
STOP_ON_ERROR = False
np.random.seed(RANDOM_STATE)

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.n_jobs = 48
sc.settings.set_figure_params(dpi=100, dpi_save=300, frameon=False)

print(f"📂 Input: {INPUT_PATH}")
print(f"📂 Output: {OUTPUT_DIR}")
print(f"Cell type column: {CELLTYPE_COLUMN}")
print(f"Batch key: {BATCH_KEY}")
print()

# ============================================================================
# STEP 1: Load and Validate Input
# ============================================================================

print("=" * 80)
print("STEP 1: LOADING AND VALIDATING INPUT DATA")
print("=" * 80)
print()

print("📂 Loading dataset...")
adata_full = sc.read_h5ad(INPUT_PATH)
print(f"  ✓ Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
print()

print("🔍 Validating preprocessing...")
validation_passed = True
errors = []

# Check layers
for layer in ['counts', 'log1p']:
    if layer not in adata_full.layers:
        errors.append(f"  ❌ Missing layer: '{layer}'")
        validation_passed = False
    else:
        print(f"  ✓ Found layer: '{layer}'")

# Check .raw
if adata_full.raw is None:
    print("  ⚠️  No .raw attribute (optional)")
else:
    print(f"  ✓ Found .raw with {adata_full.raw.n_vars} genes")

# Check columns
if CELLTYPE_COLUMN not in adata_full.obs.columns:
    errors.append(f"  ❌ Missing column: '{CELLTYPE_COLUMN}'")
    validation_passed = False
else:
    print(f"  ✓ Found column: '{CELLTYPE_COLUMN}'")

if BATCH_KEY not in adata_full.obs.columns:
    errors.append(f"  ❌ Missing column: '{BATCH_KEY}'")
    validation_passed = False
else:
    print(f"  ✓ Found column: '{BATCH_KEY}'")

if not validation_passed:
    print("\n❌ VALIDATION FAILED!")
    for error in errors:
        print(error)
    raise ValueError("Input validation failed")

print("\n  ✅ Validation PASSED")
print()

# Show distributions
celltype_counts = adata_full.obs[CELLTYPE_COLUMN].value_counts()
print("📊 Cell type distribution:")
for ct, count in celltype_counts.head(15).items():
    print(f"  {ct}: {count:,} cells")
if len(celltype_counts) > 15:
    print(f"  ... and {len(celltype_counts) - 15} more")
print()

batch_counts = adata_full.obs[BATCH_KEY].value_counts()
print(f"📊 Batch distribution: {len(batch_counts)} batches total")
print()

# ============================================================================
# STEP 2: Determine Cell Types to Process
# ============================================================================

print("=" * 80)
print("STEP 2: DETERMINING CELL TYPES TO PROCESS")
print("=" * 80)
print()

all_celltypes = list(celltype_counts.index)

# Apply filters
if CELLTYPE_WHITELIST:
    celltypes_to_process = [ct for ct in all_celltypes if ct in CELLTYPE_WHITELIST]
else:
    celltypes_to_process = all_celltypes

if CELLTYPE_BLACKLIST:
    celltypes_to_process = [ct for ct in celltypes_to_process if ct not in CELLTYPE_BLACKLIST]

# Apply minimum cell count
celltypes_filtered = []
celltypes_too_small = []
for ct in celltypes_to_process:
    n_cells = celltype_counts[ct]
    if n_cells >= MIN_CELLS_FOR_CELLTYPE:
        celltypes_filtered.append(ct)
    else:
        celltypes_too_small.append((ct, n_cells))

if celltypes_too_small:
    print(f"⚠️  Excluded {len(celltypes_too_small)} cell types with < {MIN_CELLS_FOR_CELLTYPE} cells:")
    for ct, n_cells in celltypes_too_small[:5]:
        print(f"  {ct}: {n_cells} cells")
    if len(celltypes_too_small) > 5:
        print(f"  ... and {len(celltypes_too_small) - 5} more")
    print()

celltypes_to_process = celltypes_filtered

if not celltypes_to_process:
    raise ValueError("No cell types to process!")

print(f"✅ Will process {len(celltypes_to_process)} cell types:")
for ct in celltypes_to_process:
    n_cells = celltype_counts[ct]
    has_params = "✓" if ct in TNK_PARAMS else "⚙️"
    print(f"  {ct}: {n_cells:,} cells {has_params}")
print()

# ============================================================================
# STEP 3: Process Each Cell Type (Subcluster Analysis)
# ============================================================================

print("=" * 80)
print("STEP 3: SUBCLUSTERING EACH CELL TYPE")
print("=" * 80)
print("\n⭐ KEY: Standardizing L2/L3 annotations IN each object")
print()

adata_list = []
celltype_metadata = []
marker_tables_list = []

start_time = datetime.now()

for idx, celltype in enumerate(celltypes_to_process, 1):
    celltype_start = datetime.now()
    
    print("=" * 80)
    print(f"Processing [{idx}/{len(celltypes_to_process)}]: {celltype}")
    print("=" * 80)
    print()
    
    try:
        # ===== P2-9 FIX: Use local variables to avoid global pollution =====
        # Get parameters
        if celltype in TNK_PARAMS:
            params = TNK_PARAMS[celltype]
            N_HVG = params.get("n_hvg", DEFAULT_N_HVG)
            RESOLUTION = params.get("resolution", DEFAULT_RESOLUTION)
            min_cells_per_cluster_ct = params.get("min_cells_per_cluster", 10)  # Local variable
            KNOWN_MARKERS = params.get("known_markers", [])
            print(f"  ✓ Using T/NK-specific parameters")
        else:
            N_HVG = DEFAULT_N_HVG
            RESOLUTION = DEFAULT_RESOLUTION
            min_cells_per_cluster_ct = 10  # Local variable
            KNOWN_MARKERS = []
            print(f"  ⚙️  Using default parameters")
        
        print(f"  Parameters: HVG={N_HVG}, Resolution={RESOLUTION}, MinCluster={min_cells_per_cluster_ct}")
        print()
        
        # Create safe name
        celltype_safe = celltype.lower().replace(' ', '_').replace('/', '_').replace('+', 'plus')
        
        # Filter cells
        print(f"🔬 Filtering cells...")
        cell_mask = adata_full.obs[CELLTYPE_COLUMN] == celltype
        adata_ct = adata_full[cell_mask].copy()
        n_cells_initial = adata_ct.n_obs
        print(f"  ✓ Selected {n_cells_initial:,} cells")
        print()
        
        # Setup data
        print(f"🔬 Setting up data...")
        adata_ct.X = adata_ct.layers['log1p']
        if adata_ct.X.dtype != np.float32:
            adata_ct.X = adata_ct.X.astype(np.float32)
        print(f"  ✓ Using log1p data (float32)")
        print()
        
        # QC: Remove small batches
        print(f"🔍 Quality control...")
        batch_counts_ct = adata_ct.obs[BATCH_KEY].value_counts()
        small_batches = batch_counts_ct[batch_counts_ct < MIN_CELLS_PER_BATCH].index
        if len(small_batches) > 0:
            print(f"  Removing {len(small_batches)} small batches")
            adata_ct = adata_ct[~adata_ct.obs[BATCH_KEY].isin(small_batches)].copy()
            print(f"  ✓ Remaining: {adata_ct.n_obs:,} cells")
        else:
            print(f"  ✓ All batches have ≥{MIN_CELLS_PER_BATCH} cells")
        print()
        
        # HVG selection
        print(f"🎯 HVG selection (n={N_HVG})...")
        try:
            sc.pp.highly_variable_genes(
                adata_ct, layer='log1p', n_top_genes=N_HVG,
                batch_key=BATCH_KEY, flavor=DEFAULT_HVG_FLAVOR, subset=False
            )
            hvg_method = "batch-aware"
        except:
            sc.pp.highly_variable_genes(
                adata_ct, layer='log1p', n_top_genes=N_HVG,
                flavor=DEFAULT_HVG_FLAVOR, subset=False
            )
            hvg_method = "non-batch-aware"
        
        n_hvg = adata_ct.var['highly_variable'].sum()
        print(f"  ✓ Selected {n_hvg} HVGs ({hvg_method})")
        
        # Force include markers
        if KNOWN_MARKERS:
            known_in_data = [m for m in KNOWN_MARKERS if m in adata_ct.var_names]
            if known_in_data:
                adata_ct.var.loc[known_in_data, 'highly_variable'] = True
                n_hvg = adata_ct.var['highly_variable'].sum()
                print(f"  ✓ Force-included {len(known_in_data)} markers (final: {n_hvg} HVGs)")
        print()
        
        # Create HVG subset
        print(f"✂️  Creating HVG subset...")
        adata_hvg = adata_ct[:, adata_ct.var['highly_variable']].copy()
        print(f"  ✓ HVG subset: {adata_hvg.n_obs:,} × {adata_hvg.n_vars}")
        print()
        
        # PCA
        print(f"📊 PCA...")
        sc.tl.pca(adata_hvg, n_comps=DEFAULT_N_PCS, svd_solver='arpack', random_state=RANDOM_STATE)
        print(f"  ✓ Computed {DEFAULT_N_PCS} PCs")
        print()
        
        # BBKNN
        print(f"🔗 BBKNN integration...")
        
        # ===== P0-1 CRITICAL FIX: neighbors_within_batch calculation =====
        # Bug: When min_batch_size=3, old code would give neighbors_within_batch=3
        # which violates bbknn constraint (neighbors < batch_size)
        # Fix: Upper bound must be (min_batch_size - 1), lower bound is 1
        batch_sizes = adata_hvg.obs[BATCH_KEY].value_counts()
        min_batch_size = int(batch_sizes.min())
        neighbors_within_batch = max(1, min(DEFAULT_NEIGHBORS, min_batch_size - 1))
        
        # Sanity check
        if neighbors_within_batch < 1:
            raise ValueError(
                f"Invalid neighbors_within_batch={neighbors_within_batch} "
                f"(min_batch_size={min_batch_size})"
            )
        
        print(f"  Batch sizes: min={min_batch_size}, max={batch_sizes.max()}")
        print(f"  neighbors_within_batch: {neighbors_within_batch}")
        
        sce.pp.bbknn(
            adata_hvg, batch_key=BATCH_KEY,
            neighbors_within_batch=neighbors_within_batch,
            n_pcs=DEFAULT_N_PCS, trim=None
        )
        print(f"  ✓ BBKNN completed")
        print()
        
        # Leiden clustering
        print(f"🎨 Leiden clustering (res={RESOLUTION})...")
        sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=RANDOM_STATE)
        n_clusters_initial = adata_hvg.obs['leiden'].nunique()
        print(f"  ✓ Initial: {n_clusters_initial} subclusters")
        
        # Check cluster sizes using local variable
        cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
        small_clusters = cluster_counts[cluster_counts < min_cells_per_cluster_ct].index.tolist()
        
        if small_clusters:
            print(f"  ⚠️  Removing {len(small_clusters)} small clusters (< {min_cells_per_cluster_ct} cells)")
            adata_hvg = adata_hvg[~adata_hvg.obs['leiden'].isin(small_clusters)].copy()
            
            # Re-run BBKNN and clustering
            sce.pp.bbknn(
                adata_hvg, batch_key=BATCH_KEY,
                neighbors_within_batch=neighbors_within_batch,
                n_pcs=DEFAULT_N_PCS, trim=None
            )
            sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=RANDOM_STATE)
            n_clusters = adata_hvg.obs['leiden'].nunique()
            print(f"  ✓ Final: {adata_hvg.n_obs:,} cells, {n_clusters} subclusters")
        else:
            n_clusters = n_clusters_initial
        
        final_cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
        print(f"  Subcluster sizes: {dict(final_cluster_counts)}")
        print()
        
        # BBKNN UMAP
        print(f"🗺️  BBKNN UMAP...")
        sc.tl.umap(adata_hvg, random_state=RANDOM_STATE)
        print(f"  ✓ UMAP computed")
        print()
        
        # Transfer to full object
        print(f"🔄 Transferring results...")
        adata_ct = adata_ct[adata_hvg.obs_names].copy()
        
        # ===== P0-2 FIX: Preserve categorical type for leiden =====
        # Bug: Using .values strips categorical type, causing .cat.categories to fail
        # Fix: Copy the column directly to preserve dtype, then ensure categorical
        adata_ct.obs['leiden'] = adata_hvg.obs['leiden'].copy()
        adata_ct.obs['leiden'] = adata_ct.obs['leiden'].astype('category')  # Double insurance
        
        adata_ct.obsm['X_umap_bbknn'] = adata_hvg.obsm['X_umap']
        
        # ===== P0-1 FIX: Standardize annotations =====
        print(f"🔧 Standardizing annotations (L2/L3)...")
        adata_ct.obs['cell_type_L2'] = celltype
        adata_ct.obs['cell_type_L3'] = celltype_safe + "_c" + adata_ct.obs['leiden'].astype(str)
        adata_ct.obs['subcluster_id'] = adata_ct.obs['leiden'].astype(str)
        print(f"  L2: {celltype}")
        print(f"  L3: {n_clusters} subclusters ({celltype_safe}_c0 to _c{n_clusters-1})")
        print()
        
        # ===== P0-3 FIX: Add celltype prefix to obs_names =====
        print(f"🔧 Adding celltype prefix to obs_names...")
        original_names = adata_ct.obs_names.tolist()
        adata_ct.obs_names = [f"{celltype_safe}::{x}" for x in original_names]
        print(f"  Example: {adata_ct.obs_names[0]}")
        print()
        
        del adata_hvg
        gc.collect()
        
        # Marker gene analysis
        print(f"🔬 Marker gene analysis...")
        marker_success = False
        try:
            sc.tl.rank_genes_groups(
                adata_ct, groupby='leiden', method='wilcoxon',
                layer='log1p', use_raw=False,
                corr_method='benjamini-hochberg',
                pts=True, tie_correct=True, key_added='rank_genes_groups'
            )
            print(f"  ✓ Markers computed (Wilcoxon + FDR)")
            marker_success = True
            
            # ===== P1-8 FIX: Limit to top 200 genes to reduce file size =====
            # Export markers with L2/L3 info
            if USE_MODERN_MARKER_EXPORT:
                dfs = []
                for g in adata_ct.obs['leiden'].cat.categories:
                    df = sc.get.rank_genes_groups_df(adata_ct, group=g, key='rank_genes_groups')
                    # Limit to top 200 genes
                    df = df.head(200)
                    df['cell_type_L2'] = celltype
                    df['cell_type_L3'] = f"{celltype_safe}_c{g}"
                    df['subcluster_id'] = str(g)
                    dfs.append(df)
                marker_df = pd.concat(dfs, ignore_index=True)
            else:
                # Legacy method with top 200 limit
                result = adata_ct.uns['rank_genes_groups']
                marker_list = []
                for cluster in adata_ct.obs['leiden'].cat.categories:
                    # Get top 200 genes
                    cluster_df = pd.DataFrame({
                        'cell_type_L2': celltype,
                        'cell_type_L3': f"{celltype_safe}_c{cluster}",
                        'subcluster_id': str(cluster),
                        'gene': result['names'][cluster][:200],
                        'log2fc': result['logfoldchanges'][cluster][:200],
                        'pval': result['pvals'][cluster][:200],
                        'pval_adj': result['pvals_adj'][cluster][:200]
                    })
                    marker_list.append(cluster_df)
                marker_df = pd.concat(marker_list, ignore_index=True)
            
            marker_tables_list.append(marker_df)
            print(f"  ✓ Marker table prepared ({len(marker_df)} rows, top 200 genes per cluster)")
            
        except Exception as e:
            print(f"  ⚠️  Marker analysis failed: {e}")
            marker_success = False
        
        print()
        
        # Save per-celltype figure
        print(f"🎨 Generating per-celltype figure...")
        fig, axes = plt.subplots(1, 3, figsize=(21, 6))
        
        # Subclusters on BBKNN UMAP
        sc.pl.embedding(
            adata_ct, basis='umap_bbknn', color='leiden',
            title=f'{celltype}\nSubclusters (BBKNN UMAP)',
            ax=axes[0], show=False, legend_loc='on data'
        )
        
        # Batch distribution
        sc.pl.embedding(
            adata_ct, basis='umap_bbknn', color=BATCH_KEY,
            title=f'{celltype}\nBatch Distribution',
            ax=axes[1], show=False
        )
        
        # L3 labels
        sc.pl.embedding(
            adata_ct, basis='umap_bbknn', color='cell_type_L3',
            title=f'{celltype}\nL3 Labels',
            ax=axes[2], show=False
        )
        
        plt.tight_layout()
        fig_file = FIG_DIR / f"{celltype_safe}_subclustering.pdf"
        plt.savefig(fig_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {fig_file.name}")
        print()
        
        # Store metadata
        celltype_metadata.append({
            'celltype_safe': celltype_safe,
            'celltype_original': celltype,
            'n_cells_initial': n_cells_initial,
            'n_cells_final': adata_ct.n_obs,
            'n_genes': adata_ct.n_vars,
            'n_subclusters': n_clusters,
            'hvg_method': hvg_method,
            'marker_analysis': 'SUCCESS' if marker_success else 'FAILED'
        })
        
        # Add to merge list
        adata_list.append(adata_ct)
        
        runtime = (datetime.now() - celltype_start).total_seconds()
        print(f"✅ {celltype} completed ({runtime:.1f}s)")
        print()
        
    except Exception as e:
        print(f"❌ ERROR processing {celltype}: {e}")
        
        if STOP_ON_ERROR:
            raise
        else:
            print(f"Continuing with next cell type...")
            print()
            continue

if len(adata_list) == 0:
    raise ValueError("No cell types processed successfully!")

print(f"✅ Successfully processed {len(adata_list)}/{len(celltypes_to_process)} cell types")
print()

# ============================================================================
# STEP 4: Merge All Cell Types
# ============================================================================

print("=" * 80)
print("STEP 4: MERGING ALL CELL TYPES")
print("=" * 80)
print()

print(f"Merging {len(adata_list)} AnnData objects...")
print(f"  Using: join='outer', fill_value=0")

print(f"📦 Merging {len(adata_list)} AnnData objects...")

# ===== P0-4 CRITICAL FIX: Validate var_names consistency before merge =====
# Bug: Using join='outer' is unnecessary since all objects come from same adata_full
# Risk: Can trigger gene union + fill_value=0, causing memory/sparse issues
# Fix: Validate var_names are identical, then use join='inner'
print(f"  Validating var_names consistency...")
ref_vars = adata_list[0].var_names
for i, a in enumerate(adata_list[1:], 1):
    if not ref_vars.equals(a.var_names):
        raise ValueError(
            f"var_names mismatch at index {i}. "
            f"This should not happen since all come from same adata_full. "
            f"DO NOT use join='outer' here."
        )
print(f"  ✓ All var_names consistent ({len(ref_vars)} genes)")
print(f"  Using: join='inner' (safe since var_names are identical)")

adata_merged = sc.concat(
    adata_list,
    join='inner',  # P0-4 fix: use 'inner' instead of 'outer'
    label=None,
    index_unique=None
)

print(f"  ✓ Merged: {adata_merged.n_obs:,} cells × {adata_merged.n_vars} genes")
print()

# Clean up
del adata_list
gc.collect()

# ===== P0-3 FIX: Ensure categorical types after merge =====
print(f"🔧 Ensuring categorical types after merge...")
for col in ['cell_type_L2', 'cell_type_L3', 'subcluster_id', BATCH_KEY]:
    if col in adata_merged.obs.columns:
        adata_merged.obs[col] = adata_merged.obs[col].astype('category')
        print(f"  ✓ {col}: {adata_merged.obs[col].dtype}")
print()

# Verify annotations
print(f"📊 Verifying annotations...")
print(f"  L2 cell types: {adata_merged.obs['cell_type_L2'].nunique()}")
print(f"  L3 subclusters: {adata_merged.obs['cell_type_L3'].nunique()}")
print()

# ============================================================================
# STEP 5: Global HVG + PCA + UMAP
# ============================================================================

print("=" * 80)
print("STEP 5: COMPUTING GLOBAL UMAP FROM HVG/PCA")
print("=" * 80)
print("\n⭐ KEY: Recomputing from HVG for comparability (P0-2 fix)")
print()

# Global HVG selection
print(f"🎯 Global HVG selection...")
N_HVG_GLOBAL = 4000
try:
    sc.pp.highly_variable_genes(
        adata_merged, layer='log1p', n_top_genes=N_HVG_GLOBAL,
        batch_key=BATCH_KEY, flavor='seurat', subset=False
    )
    hvg_method = "batch-aware (seurat)"
except:
    sc.pp.highly_variable_genes(
        adata_merged, layer='log1p', n_top_genes=N_HVG_GLOBAL,
        flavor='seurat', subset=False
    )
    hvg_method = "non-batch-aware (seurat)"

n_hvg_global = adata_merged.var['highly_variable'].sum()
print(f"  ✓ Selected {n_hvg_global} global HVGs ({hvg_method})")
print()

# Create HVG subset for global analysis
print(f"✂️  Creating HVG subset...")
adata_hvg_global = adata_merged[:, adata_merged.var['highly_variable']].copy()
print(f"  ✓ {adata_hvg_global.n_obs:,} × {adata_hvg_global.n_vars}")
print()

# Global PCA
print(f"📊 Global PCA...")
sc.tl.pca(adata_hvg_global, n_comps=50, svd_solver='arpack', random_state=RANDOM_STATE)
print(f"  ✓ Computed 50 PCs")
print()

# Global UMAP
print(f"🗺️  Global UMAP...")

# ===== P1-6 FIX: Use key_added to avoid overwriting default neighbors =====
sc.pp.neighbors(
    adata_hvg_global, 
    n_pcs=30, 
    random_state=RANDOM_STATE,
    key_added='neighbors_global'  # Isolate from default 'neighbors'
)

# Try using neighbors_key if supported
try:
    sc.tl.umap(
        adata_hvg_global, 
        random_state=RANDOM_STATE,
        neighbors_key='neighbors_global'
    )
    print(f"  ✓ Global UMAP computed (using key isolation)")
except TypeError:
    # Older scanpy doesn't support neighbors_key
    sc.tl.umap(adata_hvg_global, random_state=RANDOM_STATE)
    print(f"  ✓ Global UMAP computed")
print()

# Transfer to merged object
print(f"🔄 Transferring global UMAP...")
adata_merged.obsm['X_umap_global'] = adata_hvg_global.obsm['X_umap']
print(f"  ✓ Transferred")
print()

del adata_hvg_global
gc.collect()

# ============================================================================
# STEP 6: Generate Visualizations
# ============================================================================

print("=" * 80)
print("STEP 6: GENERATING VISUALIZATIONS")
print("=" * 80)
print()

# Global UMAP summary (2x2 grid)
print(f"🎨 Global UMAP summary...")
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

sc.pl.embedding(
    adata_merged, basis='umap_global', color='cell_type_L2',
    title='Cell Types (L2)', ax=axes[0, 0], show=False
)

sc.pl.embedding(
    adata_merged, basis='umap_global', color='cell_type_L3',
    title='Subclusters (L3)', ax=axes[0, 1], show=False, legend_loc=None
)

sc.pl.embedding(
    adata_merged, basis='umap_global', color=BATCH_KEY,
    title='Batch Distribution', ax=axes[1, 0], show=False
)

# Cell count per L2
l2_counts = adata_merged.obs['cell_type_L2'].value_counts()
axes[1, 1].bar(range(len(l2_counts)), l2_counts.values)
axes[1, 1].set_xticks(range(len(l2_counts)))
axes[1, 1].set_xticklabels(l2_counts.index, rotation=45, ha='right', fontsize=8)
axes[1, 1].set_ylabel('Cell Count')
axes[1, 1].set_title('Cell Count by L2 Type')

plt.tight_layout()
summary_file = FIG_DIR / 'global_umap_summary.pdf'
plt.savefig(summary_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: {summary_file.name}")
print()

# High-res L3 UMAP
print(f"🎨 High-resolution L3 UMAP...")
fig, ax = plt.subplots(figsize=(14, 12))
sc.pl.embedding(
    adata_merged, basis='umap_global', color='cell_type_L3',
    title='T/NK Cell Subclusters (L3 Resolution)',
    ax=ax, show=False, legend_loc='right margin', legend_fontsize=6
)
plt.tight_layout()
highres_file = FIG_DIR / 'umap_level3_highres.pdf'
plt.savefig(highres_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: {highres_file.name}")
print()

# ============================================================================
# STEP 7: Generate Tables
# ============================================================================

print("=" * 80)
print("STEP 7: GENERATING TABLES")
print("=" * 80)
print()

# L2 counts
print(f"📊 Cell type counts...")
l2_counts = adata_merged.obs['cell_type_L2'].value_counts()
l2_df = pd.DataFrame({
    'cell_type': l2_counts.index,
    'count': l2_counts.values,
    'percentage': (l2_counts.values / len(adata_merged) * 100).round(2)
}).sort_values('count', ascending=False)

l2_file = TABLE_DIR / 'celltype_L2_counts.csv'
l2_df.to_csv(l2_file, index=False)
print(f"  ✓ L2 counts: {l2_file.name}")

# L3 counts
l3_counts = adata_merged.obs['cell_type_L3'].value_counts()
l3_df = pd.DataFrame({
    'cell_type': l3_counts.index,
    'count': l3_counts.values,
    'percentage': (l3_counts.values / len(adata_merged) * 100).round(2)
}).sort_values('count', ascending=False)

l3_file = TABLE_DIR / 'celltype_L3_counts.csv'
l3_df.to_csv(l3_file, index=False)
print(f"  ✓ L3 counts: {l3_file.name}")

# Crosstab
crosstab = pd.crosstab(
    adata_merged.obs['cell_type_L2'],
    adata_merged.obs['cell_type_L3']
)
crosstab_file = TABLE_DIR / 'L2_vs_L3_crosstab.csv'
crosstab.to_csv(crosstab_file)
print(f"  ✓ Crosstab: {crosstab_file.name}")

# Subcluster summary
print(f"📊 Subcluster summary...")
subcluster_summary = []

for celltype in adata_merged.obs['cell_type_L2'].cat.categories:
    mask_ct = adata_merged.obs['cell_type_L2'] == celltype
    n_cells_celltype = mask_ct.sum()
    
    for subcluster in adata_merged.obs.loc[mask_ct, 'cell_type_L3'].unique():
        mask_sub = (adata_merged.obs['cell_type_L2'] == celltype) & \
                   (adata_merged.obs['cell_type_L3'] == subcluster)
        n_cells_sub = mask_sub.sum()
        
        sub_ids = adata_merged.obs.loc[mask_sub, 'subcluster_id'].unique()
        subcluster_id = sub_ids[0] if len(sub_ids) > 0 else 'NA'
        
        subcluster_summary.append({
            'celltype_L2': celltype,
            'subcluster_label': subcluster,
            'subcluster_id': subcluster_id,
            'n_cells': n_cells_sub,
            'percentage_within_celltype': round(n_cells_sub / n_cells_celltype * 100, 2),
            'percentage_total': round(n_cells_sub / len(adata_merged) * 100, 2)
        })

summary_df = pd.DataFrame(subcluster_summary)
summary_df = summary_df.sort_values(['celltype_L2', 'subcluster_id'])
summary_file = TABLE_DIR / 'subcluster_summary.csv'
summary_df.to_csv(summary_file, index=False)
print(f"  ✓ Subcluster summary: {summary_file.name}")

# Final annotations
print(f"📊 Final annotations...")
annotation_cols = ['cell_type_L2', 'cell_type_L3', 'subcluster_id']

if BATCH_KEY in adata_merged.obs.columns:
    annotation_cols.insert(0, BATCH_KEY)

for col in ['disease_status', 'tissue', 'patient_id', 'sample']:
    if col in adata_merged.obs.columns and col not in annotation_cols:
        annotation_cols.append(col)

final_annotations = adata_merged.obs[annotation_cols].copy()
annot_file = TABLE_DIR / 'final_annotations.csv'
final_annotations.to_csv(annot_file)
print(f"  ✓ Final annotations: {annot_file.name}")

# Export markers (combined table with L2/L3 info)
if marker_tables_list:
    print(f"📊 Exporting markers...")
    all_markers = pd.concat(marker_tables_list, ignore_index=True)
    marker_file = TABLE_DIR / 'all_markers_combined.csv'
    all_markers.to_csv(marker_file, index=False)
    print(f"  ✓ Combined markers: {marker_file.name} ({len(all_markers)} rows)")
    
    # Also save per-celltype
    for celltype_meta in celltype_metadata:
        ct_safe = celltype_meta['celltype_safe']
        ct_markers = all_markers[all_markers['cell_type_L2'] == celltype_meta['celltype_original']]
        ct_marker_file = TABLE_DIR / f"{ct_safe}_markers.csv"
        ct_markers.to_csv(ct_marker_file, index=False)
    print(f"  ✓ Per-celltype marker tables saved")

print()

# ============================================================================
# STEP 8: Save Final H5AD
# ============================================================================

print("=" * 80)
print("STEP 8: SAVING FINAL H5AD")
print("=" * 80)
print()

# Add comprehensive metadata
adata_merged.uns['subcluster_analysis'] = {
    'date': datetime.now().strftime('%Y-%m-%d'),
    'version': 'unified_v2.0.1_HOTFIX',
    'workflow': 'all_in_one_subcluster_plus_merge',
    'fixes_applied': {
        'P0_critical': {
            'neighbors_within_batch': 'Fixed calculation when min_batch_size=3',
            'categorical_preservation': 'Preserved categorical type for leiden transfer',
            'join_outer_risk': 'Validated var_names consistency, used join=inner',
            'local_umap_removal': 'Removed per-celltype X_umap_bbknn from merged object',
            'categorical_after_merge': 'Ensured categorical types for L2/L3/batch columns'
        },
        'P1_high_priority': {
            'thread_control': 'Limited BLAS threads (OMP/MKL/OpenBLAS/NUMEXPR=1)',
            'key_isolation': 'Used neighbors_global key to avoid overwriting defaults',
            'marker_size_limit': 'Limited marker export to top 200 genes per cluster'
        },
        'P2_maintainability': {
            'local_variables': 'Used min_cells_per_cluster_ct to avoid global pollution',
            'config_cleanup': 'Removed unused PROCESS_ALL_CELLTYPES config'
        }
    },
    'n_celltypes_L2': adata_merged.obs['cell_type_L2'].nunique(),
    'n_subclusters_L3': adata_merged.obs['cell_type_L3'].nunique(),
    'total_cells': adata_merged.n_obs,
    'total_genes': adata_merged.n_vars,
    'batch_key': BATCH_KEY,
    'hvg_method': hvg_method,
    'source_celltypes': [m['celltype_original'] for m in celltype_metadata]
}

# Store celltype metadata
adata_merged.uns['source_metadata'] = {
    m['celltype_safe']: {
        'celltype_original': m['celltype_original'],
        'n_cells_initial': m['n_cells_initial'],
        'n_cells_final': m['n_cells_final'],
        'n_genes': m['n_genes'],
        'n_subclusters': m['n_subclusters'],
        'hvg_method': m['hvg_method'],
        'marker_analysis': m['marker_analysis']
    } for m in celltype_metadata
}

# ===== P0-5 FIX: Remove local X_umap_bbknn to avoid misuse =====
# Bug: X_umap_bbknn is per-celltype local coordinate, meaningless after merge
# Fix: Delete it from merged object (per-celltype PDFs already saved)
if 'X_umap_bbknn' in adata_merged.obsm:
    print(f"🗑️  Removing local X_umap_bbknn (per-celltype coordinate, not comparable after merge)")
    del adata_merged.obsm['X_umap_bbknn']
    print(f"  ✓ Removed (per-celltype figures already saved)")
    print()

# Save
output_file = OUTPUT_DIR / "adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
print(f"💾 Saving to: {output_file}")
adata_merged.write_h5ad(output_file, compression='gzip', compression_opts=9)

size_gb = output_file.stat().st_size / 1e9
print(f"  ✓ File saved: {size_gb:.2f} GB")
print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

total_runtime = (datetime.now() - start_time).total_seconds()

print("=" * 80)
print("✅ ANALYSIS COMPLETE!")
print("=" * 80)

print(f"""
📁 Output Structure
===================

Main Results:
  - {output_file.name} ⭐
  - {annot_file.name}

Figures ({len(list(FIG_DIR.glob('*.pdf')))} files):
  - global_umap_summary.pdf ⭐
  - umap_level3_highres.pdf ⭐
  - [celltype]_subclustering.pdf (per celltype)

Tables ({len(list(TABLE_DIR.glob('*.csv')))} files):
  - subcluster_summary.csv ⭐
  - celltype_L2_counts.csv
  - celltype_L3_counts.csv
  - L2_vs_L3_crosstab.csv
  - final_annotations.csv
  - all_markers_combined.csv (top 200 genes per cluster)
  - [celltype]_markers.csv (per celltype)

📊 Data Summary
===============
  - Total cells: {adata_merged.n_obs:,}
  - Genes: {adata_merged.n_vars:,}
  - Cell types (L2): {adata_merged.obs['cell_type_L2'].nunique()}
  - Subclusters (L3): {adata_merged.obs['cell_type_L3'].nunique()}
  - Source celltypes: {len(celltype_metadata)}

⏱️  Runtime
===========
  - Total: {total_runtime:.1f} seconds ({total_runtime/60:.1f} minutes)
  - Per celltype: {total_runtime/len(celltype_metadata):.1f} seconds average

🔧 HOTFIX v2.0.1 Changes
=========================
  P0 (Critical - Would Cause Errors):
    ✓ Fixed neighbors_within_batch calculation (min_batch_size=3 bug)
    ✓ Preserved categorical types (leiden transfer)
    ✓ Fixed join='outer' risk (validated var_names, used join='inner')
    ✓ Removed local X_umap_bbknn from merged object
    ✓ Ensured categorical types after merge (L2/L3/batch)
  
  P1 (High Priority - Performance/Robustness):
    ✓ Limited BLAS threading (OMP/MKL/OpenBLAS/NUMEXPR=1)
    ✓ Used neighbors_global key isolation
    ✓ Limited marker tables to top 200 genes per cluster
  
  P2 (Maintainability):
    ✓ Used local min_cells_per_cluster_ct variable
    ✓ Cleaned unused config items

📝 Next Steps
=============
  1. Validate global UMAP separates cell types correctly
  2. Review per-celltype subclustering quality
  3. Check marker gene tables for biological sense
  4. Compare with Stromal/B cell unified formats
  5. Perform cross-lineage comparative analysis
  6. Run differential expression between subclusters

⚠️  Important Notes
===================
  - Global UMAP recomputed from HVG/PCA (not per-celltype latent)
  - obs_names have celltype prefix (e.g., "cd8plus_trm::BARCODE")
  - Marker tables limited to top 200 genes (manageable file size)
  - Local X_umap_bbknn removed (meaningless after merge, per-celltype PDFs saved)
  - All figure labels and annotations in English (publication-ready)
  - Thread-safe: BLAS libraries limited to 1 thread, scanpy n_jobs=48
""")

print(f"Output directory: {OUTPUT_DIR}")
print(f"Figures: {FIG_DIR}")
print(f"Tables: {TABLE_DIR}")
print("\n" + "="*80)
print("🎉 PRODUCTION-READY v2.0.1-HOTFIX OUTPUT GENERATED")
print("="*80)
print("\n✅ All P0/P1/P2 fixes applied - Ready for production use")

gc.collect()

# ============================================================================
# END OF UNIFIED PIPELINE v2.0.1 (HOTFIX)
# ============================================================================
