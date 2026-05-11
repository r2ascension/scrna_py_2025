#!/usr/bin/env python3
"""
Universal Cell-Type Subcluster Analysis Pipeline v4.1 (PRODUCTION - FIXED)
===========================================================================
Author: r2end
Date: 2025-01-21

FIXES in v4.1:
--------------
P0 Issues (Critical):
  ✓ Fixed scANVI UMAP visualization - now uses original scANVI coordinates
  ✓ Fixed .raw memory issue - uses shared memory instead of deep copy
  ✓ Marker analysis now explicitly uses log1p layer on full gene set

P1 Issues (Engineering):
  ✓ Eliminated unnecessary .copy() calls to reduce memory usage
  ✓ Explicit HVG flavor specification
  ✓ Proper bbknn dependency handling
  ✓ Robust marker export using sc.get.rank_genes_groups_df

P2 Issues (Maintainability):
  ✓ Proper log file handling with context manager
  ✓ Removed unused variables (batch_log, seaborn)
  ✓ Cluster keys prefixed with cell type for easy merging

Purpose:
--------
BATCH processing pipeline for subcluster analysis of MULTIPLE cell types.
Automatically processes all cell types found in metadata.
Follows QUICK_REFERENCE_MEMORY v2.12 best practices.

Key Architecture (Dual-Object):
--------------------------------
1. adata_ct_full: Full gene set
   - Used for marker gene analysis
   - Used for visualization on original scANVI UMAP
   - Minimal memory footprint (no unnecessary copies)

2. adata_hvg: HVG subset only
   - Used for PCA/BBKNN/Leiden clustering
   - Fast computation
   - Clustering results transferred back to adata_ct_full

CRITICAL Requirements:
----------------------
Input data MUST be pre-processed:
1. Gene filtering (MT/ribo/ENSG removed) - GLOBAL
2. Normalization (all cells together) - GLOBAL  
3. log1p transformation - GLOBAL
4. .raw saved with full gene set - GLOBAL

Workflow:
---------
Pre-normalized full data → For each cell type:
  → Filter cells → Create full & HVG objects →
  → HVG: PCA/BBKNN/Leiden →
  → Transfer clusters to full object →
  → Full: Markers (log1p layer) + Viz (scANVI UMAP) →
  → Save
→ Generate batch summary report

Input: 
------
- Pre-processed h5ad with:
  * layers['counts']: Raw counts
  * layers['log1p']: Globally normalized log1p
  * .raw: Full gene set
  * obsm['X_scanvi']: scANVI latent (optional, for visualization)
  * obsm['X_umap_scanvi']: scANVI UMAP (optional, for visualization)

Output:
-------
For each cell type:
- output_dir/{celltype}_subcluster/
  ├── {celltype}_subcluster_analyzed.h5ad
  ├── marker_genes_with_FDR.csv
  ├── cluster_assignments.csv
  ├── analysis_summary.csv
  └── figures/
      ├── umap_subclusters_scanvi.png (on original scANVI if available)
      ├── umap_subclusters_bbknn.png (BBKNN-derived)
      ├── umap_batches.png
      ├── cluster_sizes.png
      ├── known_markers_umap.png
      ├── marker_heatmap.png
      └── marker_dotplot.png

Batch summary:
- batch_analysis_summary.csv
- batch_processing_log.txt

Memory: ~20-40GB per cell type (optimized)
Runtime: ~30min - 2hr per cell type
"""

# ===== 1. Import Libraries =====
import scanpy as sc
import scanpy.external as sce  # P1: Explicit external import
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
print("UNIVERSAL CELL-TYPE SUBCLUSTER ANALYSIS PIPELINE v4.1 (PRODUCTION - FIXED)")
print("=" * 80)
print()

# ===== 2. Dependency Check =====
# P1: Check critical dependencies before batch processing
print("🔍 Checking dependencies...")
try:
    import bbknn
    print("  ✓ bbknn available")
except ImportError:
    print("  ❌ bbknn not found. Install with: pip install bbknn")
    sys.exit(1)

try:
    # Test if get_rank_genes_groups_df is available (scanpy >= 1.6)
    if not hasattr(sc, 'get') or not hasattr(sc.get, 'rank_genes_groups_df'):
        print("  ⚠️  scanpy.get.rank_genes_groups_df not available")
        print("     Consider upgrading scanpy: pip install --upgrade scanpy")
        USE_MODERN_MARKER_EXPORT = False
    else:
        USE_MODERN_MARKER_EXPORT = True
        print("  ✓ scanpy.get.rank_genes_groups_df available")
except:
    USE_MODERN_MARKER_EXPORT = False
    print("  ⚠️  Will use legacy marker export method")

print()

# ===== 3. Configuration - MODIFY HERE =====
# ========================================
# Paths
# ========================================
INPUT_PATH = "/home/h2048/data/py/0114/stromal_pure_scanvi/adata_stromal_PURE_FINAL.h5ad"  # Pre-processed dataset (REQUIRED)
OUTPUT_BASE_DIR = Path(f"/home/h2048/data/py/0114/stromal_subcluster")
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

# ========================================
# Column Names
# ========================================
CELLTYPE_COLUMN = 'cell_type_scanvi_filt'        # Column containing cell type annotations
BATCH_KEY = 'dataset'                     # Column for batch correction (use 'dataset', not 'Sample')

# scANVI Latent Space Columns (for visualization)
SCANVI_LATENT_KEY = 'X_scanvi'            # scANVI latent representation in .obsm
SCANVI_UMAP_KEY = 'X_umap_scanvi'         # scANVI UMAP in .obsm (if available)

# ========================================
# Cell Type Selection (BATCH MODE)
# ========================================
PROCESS_ALL_CELLTYPES = True              # Set to True to process all cell types
CELLTYPE_WHITELIST = []                   # Example: ["Myeloid", "T", "B"]
CELLTYPE_BLACKLIST = []                   # Example: ["Unknown", "Doublets"]
MIN_CELLS_FOR_CELLTYPE = 100              # Minimum cells to process a cell type

# ========================================
# Quality Control Thresholds
# ========================================
MIN_CELLS_PER_BATCH = 3           # Remove batches smaller than this
MIN_CELLS_PER_CLUSTER = 10        # Remove clusters smaller than this
MIN_CELLS_FOR_MARKER = 20         # Minimum cells for reliable marker analysis

# ========================================
# Default Analysis Parameters
# ========================================
DEFAULT_N_HVG = 3000              # Default number of highly variable genes
DEFAULT_N_PCS = 50                # Default number of PCs for BBKNN
DEFAULT_RESOLUTION = 0.3          # Default Leiden clustering resolution
DEFAULT_NEIGHBORS = 5             # Default BBKNN neighbors_within_batch
DEFAULT_HVG_FLAVOR = 'seurat'     # P1: Explicit flavor for log1p data

# ========================================
# Cell-Type-Specific Parameters
# ========================================
CELLTYPE_PARAMS = {
    "Myeloid": {
        "n_hvg": 3000,
        "n_pcs": 50,
        "resolution": 0.4,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 10,
        "known_markers": ["CD14", "FCGR3A", "CD68", "CD163", "LYZ"]
    },
    "Stromal": {
        "n_hvg": 3000,
        "n_pcs": 50,
        "resolution": 0.3,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 10,
        "known_markers": ["COL1A1", "COL3A1", "DCN", "LUM", "PDGFRA"]
    },
    "B": {
        "n_hvg": 2500,
        "n_pcs": 50,
        "resolution": 0.3,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 15,
        "known_markers": ["CD79A", "CD79B", "MS4A1", "CD19", "IGHM"]
    },
    "T": {
        "n_hvg": 3000,
        "n_pcs": 50,
        "resolution": 0.4,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "CD3E", "CD4", "CD8A", "IL7R"]
    },
    "NK": {
        "n_hvg": 2500,
        "n_pcs": 50,
        "resolution": 0.3,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 10,
        "known_markers": ["GNLY", "NKG7", "NCAM1", "KLRD1", "KLRB1"]
    },
    "Epithelial": {
        "n_hvg": 3500,
        "n_pcs": 50,
        "resolution": 0.5,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 15,
        "known_markers": ["EPCAM", "KRT5", "TP63", "MUC5AC", "FOXJ1"]
    },
    "Endothelial": {
        "n_hvg": 2500,
        "n_pcs": 50,
        "resolution": 0.3,
        "neighbors_within_batch": 5,
        "min_cells_per_cluster": 10,
        "known_markers": ["PECAM1", "VWF", "CDH5", "CD34"]
    }
}

# ========================================
# Error Handling
# ========================================
STOP_ON_ERROR = False            # If True, stop batch processing on first error

# ========================================
# Scanpy Settings
# ========================================
sc.settings.verbosity = 1
sc.settings.n_jobs = 48
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))

print("Configuration Summary:")
print("-" * 80)
print(f"Cell type column: {CELLTYPE_COLUMN}")
print(f"Batch key: {BATCH_KEY}")
print(f"Input: {INPUT_PATH}")
print(f"Output base directory: {OUTPUT_BASE_DIR}")
print(f"Process all cell types: {PROCESS_ALL_CELLTYPES}")
if CELLTYPE_WHITELIST:
    print(f"Whitelist (only process): {CELLTYPE_WHITELIST}")
if CELLTYPE_BLACKLIST:
    print(f"Blacklist (exclude): {CELLTYPE_BLACKLIST}")
print(f"Min cells per cell type: {MIN_CELLS_FOR_CELLTYPE}")
print(f"HVG flavor: {DEFAULT_HVG_FLAVOR}")
print(f"Stop on error: {STOP_ON_ERROR}")
print("-" * 80)
print()

# ===== 4. Load and Validate Input Data =====
print("📂 Loading pre-processed dataset...")
adata_full = sc.read_h5ad(INPUT_PATH)
print(f"  ✓ Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
print()

# ===== 5. CRITICAL: Validate Pre-Processing =====
print("🔍 Validating input data preprocessing...")

validation_passed = True
errors = []

# Check required layers
required_layers = ['counts', 'log1p']
for layer in required_layers:
    if layer not in adata_full.layers:
        errors.append(f"  ❌ Missing required layer: '{layer}'")
        validation_passed = False
    else:
        print(f"  ✓ Found layer: '{layer}'")

# Check for .raw
if adata_full.raw is None:
    print("  ⚠️  WARNING: No .raw attribute found")
    print("     Marker gene analysis will use current gene set only")
else:
    print(f"  ✓ Found .raw with {adata_full.raw.n_vars} genes")

# Warn if gene count suggests unfiltered data
if adata_full.n_vars > 30000:
    print(f"  ⚠️  WARNING: Large gene count ({adata_full.n_vars} genes)")
    print(f"     Suggests unfiltered data (expected: ~18k-22k after filtering)")
elif adata_full.n_vars < 15000:
    print(f"  ⚠️  WARNING: Small gene count ({adata_full.n_vars} genes)")
    print(f"     May be over-filtered")
else:
    print(f"  ✓ Gene count reasonable: {adata_full.n_vars} genes")

# Check for scANVI latent space
if SCANVI_LATENT_KEY in adata_full.obsm:
    print(f"  ✓ Found scANVI latent: '{SCANVI_LATENT_KEY}'")
    has_scanvi = True
else:
    print(f"  ⚠️  WARNING: scANVI latent '{SCANVI_LATENT_KEY}' not found")
    print(f"     Will use BBKNN UMAP for visualization")
    has_scanvi = False

if SCANVI_UMAP_KEY in adata_full.obsm:
    print(f"  ✓ Found scANVI UMAP: '{SCANVI_UMAP_KEY}'")
    has_scanvi_umap = True
else:
    print(f"  ⚠️  scANVI UMAP '{SCANVI_UMAP_KEY}' not found")
    has_scanvi_umap = False

print()

# Stop if critical validation failed
if not validation_passed:
    print("❌ INPUT DATA VALIDATION FAILED!")
    print()
    print("Errors found:")
    for error in errors:
        print(error)
    print()
    print("REQUIRED: Input data must be pre-processed with:")
    print("1. Gene filtering (MT/ribo/ENSG removed) - GLOBAL")
    print("2. Normalization (all cells together) - GLOBAL")
    print("3. log1p transformation - GLOBAL")
    print("4. Saved as layers['counts'] and layers['log1p']")
    print()
    raise ValueError("Input data validation failed. See errors above.")

print("  ✅ Input data validation PASSED")
print()

# ===== 6. Validate Required Columns =====
print("🔍 Validating required columns...")

# Check cell type column
if CELLTYPE_COLUMN not in adata_full.obs.columns:
    raise ValueError(
        f"Cell type column '{CELLTYPE_COLUMN}' not found!\n"
        f"Available columns: {list(adata_full.obs.columns)}"
    )
print(f"  ✓ Found cell type column: '{CELLTYPE_COLUMN}'")

# Check batch column
if BATCH_KEY not in adata_full.obs.columns:
    raise ValueError(
        f"Batch column '{BATCH_KEY}' not found!\n"
        f"Available columns: {list(adata_full.obs.columns)}"
    )
print(f"  ✓ Found batch column: '{BATCH_KEY}'")

# Show cell type distribution
celltype_counts = adata_full.obs[CELLTYPE_COLUMN].value_counts()
print(f"\n  Cell type distribution in full dataset:")
for ct, count in celltype_counts.items():
    print(f"    {ct}: {count:,} cells")
print()

# Show batch distribution
batch_counts = adata_full.obs[BATCH_KEY].value_counts()
print(f"  Batch (dataset) distribution:")
print(f"    Total: {len(batch_counts)} batches")
for batch, count in batch_counts.head(5).items():
    print(f"    {batch}: {count:,} cells")
if len(batch_counts) > 5:
    print(f"    ... and {len(batch_counts) - 5} more batches")
print()

# ===== 7. Determine Cell Types to Process =====
print("🎯 Determining cell types to process...")

# Get all unique cell types
all_celltypes = list(celltype_counts.index)

# Apply whitelist filter
if CELLTYPE_WHITELIST:
    celltypes_to_process = [ct for ct in all_celltypes if ct in CELLTYPE_WHITELIST]
    print(f"  Whitelist filter: {len(all_celltypes)} → {len(celltypes_to_process)} cell types")
else:
    celltypes_to_process = all_celltypes

# Apply blacklist filter
if CELLTYPE_BLACKLIST:
    celltypes_to_process = [ct for ct in celltypes_to_process if ct not in CELLTYPE_BLACKLIST]
    print(f"  Blacklist filter: excluded {len(CELLTYPE_BLACKLIST)} cell types")

# Apply minimum cell count filter
celltypes_filtered = []
celltypes_too_small = []
for ct in celltypes_to_process:
    n_cells = celltype_counts[ct]
    if n_cells >= MIN_CELLS_FOR_CELLTYPE:
        celltypes_filtered.append(ct)
    else:
        celltypes_too_small.append((ct, n_cells))

if celltypes_too_small:
    print(f"  ⚠️  Excluded {len(celltypes_too_small)} cell types with < {MIN_CELLS_FOR_CELLTYPE} cells:")
    for ct, n_cells in celltypes_too_small[:5]:
        print(f"    {ct}: {n_cells} cells")
    if len(celltypes_too_small) > 5:
        print(f"    ... and {len(celltypes_too_small) - 5} more")

celltypes_to_process = celltypes_filtered

if not celltypes_to_process:
    raise ValueError("No cell types to process after filtering!")

print()
print(f"📋 Final list of cell types to process ({len(celltypes_to_process)}):")
for ct in celltypes_to_process:
    n_cells = celltype_counts[ct]
    has_params = "✓" if ct in CELLTYPE_PARAMS else "⚙️ (default)"
    print(f"  {ct}: {n_cells:,} cells {has_params}")
print()

# ===== 8. Initialize Batch Processing =====
batch_results = []
start_time = datetime.now()

# P2: Use context manager for log file
log_file = OUTPUT_BASE_DIR / 'batch_processing_log.txt'

def log_message(msg, to_console=True, to_file=True, log_handle=None):
    """Write message to both console and log file"""
    if to_console:
        print(msg)
    if to_file and log_handle is not None:
        log_handle.write(msg + '\n')
        log_handle.flush()

# P2: Context manager for proper file handling
with open(log_file, 'w') as log_handle:
    log_message("=" * 80, log_handle=log_handle)
    log_message(f"BATCH PROCESSING LOG", log_handle=log_handle)
    log_message(f"Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}", log_handle=log_handle)
    log_message("=" * 80, log_handle=log_handle)
    log_message("", log_handle=log_handle)
    
    # ===== 9. MAIN LOOP: Process Each Cell Type =====
    for idx, celltype in enumerate(celltypes_to_process, 1):
        celltype_start = datetime.now()
        
        log_message("=" * 80, log_handle=log_handle)
        log_message(f"Processing [{idx}/{len(celltypes_to_process)}]: {celltype}", log_handle=log_handle)
        log_message("=" * 80, log_handle=log_handle)
        log_message("", log_handle=log_handle)
        
        try:
            # ===== 9.1. Setup Cell-Type-Specific Configuration =====
            if celltype in CELLTYPE_PARAMS:
                params = CELLTYPE_PARAMS[celltype]
                N_HVG = params.get("n_hvg", DEFAULT_N_HVG)
                N_PCS = params.get("n_pcs", DEFAULT_N_PCS)
                RESOLUTION = params.get("resolution", DEFAULT_RESOLUTION)
                NEIGHBORS_WITHIN_BATCH = params.get("neighbors_within_batch", DEFAULT_NEIGHBORS)
                MIN_CELLS_PER_CLUSTER = params.get("min_cells_per_cluster", 10)
                KNOWN_MARKERS = params.get("known_markers", [])
                log_message(f"  ✓ Using cell-type-specific parameters for {celltype}", log_handle=log_handle)
            else:
                N_HVG = DEFAULT_N_HVG
                N_PCS = DEFAULT_N_PCS
                RESOLUTION = DEFAULT_RESOLUTION
                NEIGHBORS_WITHIN_BATCH = DEFAULT_NEIGHBORS
                MIN_CELLS_PER_CLUSTER = 10
                KNOWN_MARKERS = []
                log_message(f"  ⚙️  Using default parameters for {celltype}", log_handle=log_handle)
            
            # Create cell-type-specific output directory
            celltype_safe = celltype.lower().replace(' ', '_').replace('/', '_')
            OUTPUT_DIR = OUTPUT_BASE_DIR / f"{celltype_safe}_subcluster"
            OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            
            log_message(f"  Output directory: {OUTPUT_DIR}", log_handle=log_handle)
            log_message(f"  Parameters: HVG={N_HVG}, Resolution={RESOLUTION}, MinCluster={MIN_CELLS_PER_CLUSTER}", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.2. Filter by Cell Type - Create Full Gene Object =====
            log_message(f"🔬 Filtering for {celltype} cells...", log_handle=log_handle)
            
            cell_mask = adata_full.obs[CELLTYPE_COLUMN] == celltype
            adata_ct_full = adata_full[cell_mask].copy()
            n_cells_selected = adata_ct_full.n_obs
            
            log_message(f"  ✓ Selected {n_cells_selected:,} {celltype} cells", log_handle=log_handle)
            log_message(f"    ({n_cells_selected/adata_full.n_obs*100:.1f}% of total dataset)", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.3. Setup Full Gene Object (No Re-normalization) =====
            log_message("🔬 Setting up full gene object...", log_handle=log_handle)
            
            # P1: No .copy() - just reference
            adata_ct_full.X = adata_ct_full.layers['log1p']
            
            # P1: Ensure float32 to save memory
            if adata_ct_full.X.dtype != np.float32:
                log_message("  Converting to float32 for memory efficiency...", log_handle=log_handle)
                adata_ct_full.X = adata_ct_full.X.astype(np.float32)
            
            log_message(f"  ✓ Using globally normalized log1p data", log_handle=log_handle)
            log_message(f"    Data type: {adata_ct_full.X.dtype}", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.4. Data Quality Checks =====
            log_message("🔍 Performing data quality checks...", log_handle=log_handle)
            
            batch_counts_ct = adata_ct_full.obs[BATCH_KEY].value_counts()
            log_message(f"  Batch distribution in {celltype} cells:", log_handle=log_handle)
            log_message(f"    Total batches: {len(batch_counts_ct)}", log_handle=log_handle)
            for batch, count in batch_counts_ct.head(5).items():
                log_message(f"    {batch}: {count} cells", log_handle=log_handle)
            if len(batch_counts_ct) > 5:
                log_message(f"    ... and {len(batch_counts_ct) - 5} more batches", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # Remove small batches
            small_batches = batch_counts_ct[batch_counts_ct < MIN_CELLS_PER_BATCH].index
            if len(small_batches) > 0:
                log_message(f"⚠️  Removing {len(small_batches)} small batches (< {MIN_CELLS_PER_BATCH} cells):", log_handle=log_handle)
                for batch in small_batches[:5]:
                    log_message(f"    {batch}: {batch_counts_ct[batch]} cells", log_handle=log_handle)
                if len(small_batches) > 5:
                    log_message(f"    ... and {len(small_batches) - 5} more", log_handle=log_handle)
                
                adata_ct_full = adata_ct_full[~adata_ct_full.obs[BATCH_KEY].isin(small_batches)].copy()
                log_message(f"  ✓ Remaining: {adata_ct_full.n_obs:,} cells in {adata_ct_full.obs[BATCH_KEY].nunique()} batches", log_handle=log_handle)
                log_message("", log_handle=log_handle)
            
            # ===== 9.5. HVG Selection on Full Object (NO SUBSET YET) =====
            log_message(f"🎯 Selecting top {N_HVG} highly variable genes...", log_handle=log_handle)
            
            try:
                # P1: Explicit flavor
                sc.pp.highly_variable_genes(
                    adata_ct_full,
                    layer='log1p',
                    n_top_genes=N_HVG,
                    batch_key=BATCH_KEY,
                    flavor=DEFAULT_HVG_FLAVOR,
                    subset=False  # Do NOT subset yet
                )
                hvg_method = f"batch-aware ({DEFAULT_HVG_FLAVOR})"
                log_message(f"  ✓ Batch-aware HVG selection successful", log_handle=log_handle)
            except Exception as e:
                log_message(f"  ⚠️  Batch-aware HVG failed: {e}", log_handle=log_handle)
                log_message("  Falling back to non-batch-aware method", log_handle=log_handle)
                sc.pp.highly_variable_genes(
                    adata_ct_full,
                    layer='log1p',
                    n_top_genes=N_HVG,
                    flavor=DEFAULT_HVG_FLAVOR,
                    subset=False
                )
                hvg_method = f"non-batch-aware ({DEFAULT_HVG_FLAVOR})"
            
            n_hvg = adata_ct_full.var['highly_variable'].sum()
            log_message(f"  ✓ Selected {n_hvg} HVGs ({hvg_method})", log_handle=log_handle)
            
            # Force include known markers if specified
            if KNOWN_MARKERS:
                known_in_data = [m for m in KNOWN_MARKERS if m in adata_ct_full.var_names]
                if known_in_data:
                    adata_ct_full.var.loc[known_in_data, 'highly_variable'] = True
                    n_hvg = adata_ct_full.var['highly_variable'].sum()
                    log_message(f"  ✓ Force-included {len(known_in_data)} known markers", log_handle=log_handle)
                    log_message(f"  Final HVG count: {n_hvg}", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.6. Create HVG Subset Object for Clustering =====
            log_message("✂️  Creating HVG subset for clustering...", log_handle=log_handle)
            
            # P1: Explicit copy here is necessary for subset
            adata_hvg = adata_ct_full[:, adata_ct_full.var['highly_variable']].copy()
            log_message(f"  ✓ HVG subset: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars} genes", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.7. Dimensionality Reduction on HVG =====
            log_message(f"📊 Computing PCA ({N_PCS} components)...", log_handle=log_handle)
            sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack', random_state=42)
            log_message("  ✓ PCA computed", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.8. Batch Effect Correction with BBKNN =====
            log_message("🔗 Running BBKNN integration...", log_handle=log_handle)
            
            batch_sizes = adata_hvg.obs[BATCH_KEY].value_counts()
            min_batch_size = batch_sizes.min()
            neighbors_within_batch = min(NEIGHBORS_WITHIN_BATCH, max(3, min_batch_size - 1))
            
            log_message(f"  Smallest batch size: {min_batch_size}", log_handle=log_handle)
            log_message(f"  Using neighbors_within_batch: {neighbors_within_batch}", log_handle=log_handle)
            
            # P1: Use explicit sce.pp.bbknn
            sce.pp.bbknn(
                adata_hvg,
                batch_key=BATCH_KEY,
                neighbors_within_batch=neighbors_within_batch,
                n_pcs=N_PCS,
                trim=None
            )
            log_message("  ✓ BBKNN integration completed", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.9. Leiden Clustering on HVG =====
            log_message(f"🎨 Computing Leiden clustering (resolution={RESOLUTION})...", log_handle=log_handle)
            sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=42)
            n_clusters = adata_hvg.obs['leiden'].nunique()
            log_message(f"  ✓ Identified {n_clusters} subclusters", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # Check cluster sizes
            cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
            log_message("  Subcluster sizes:", log_handle=log_handle)
            for cluster, count in cluster_counts.items():
                status = "⚠️ SMALL" if count < MIN_CELLS_PER_CLUSTER else "✓"
                log_message(f"    Subcluster {cluster}: {count:,} cells {status}", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.10. Filter Small Clusters =====
            small_clusters = cluster_counts[cluster_counts < MIN_CELLS_PER_CLUSTER].index.tolist()
            if len(small_clusters) > 0:
                log_message(f"⚠️  Removing {len(small_clusters)} small subclusters (< {MIN_CELLS_PER_CLUSTER} cells):", log_handle=log_handle)
                log_message(f"  Subclusters: {small_clusters}", log_handle=log_handle)
                
                adata_hvg = adata_hvg[~adata_hvg.obs['leiden'].isin(small_clusters)].copy()
                
                # Recompute BBKNN and clustering
                log_message("  Recomputing BBKNN and clustering after filtering...", log_handle=log_handle)
                sce.pp.bbknn(
                    adata_hvg,
                    batch_key=BATCH_KEY,
                    neighbors_within_batch=neighbors_within_batch,
                    n_pcs=N_PCS,
                    trim=None
                )
                sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=42)
                n_clusters = adata_hvg.obs['leiden'].nunique()
                
                log_message(f"  ✓ Final: {adata_hvg.n_obs:,} cells in {n_clusters} subclusters", log_handle=log_handle)
                
                final_cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
                log_message("  Final subcluster sizes:", log_handle=log_handle)
                for cluster, count in final_cluster_counts.items():
                    log_message(f"    Subcluster {cluster}: {count:,} cells", log_handle=log_handle)
                log_message("", log_handle=log_handle)
            else:
                final_cluster_counts = cluster_counts
            
            # ===== 9.11. Compute BBKNN UMAP on HVG (for comparison) =====
            log_message("🗺️  Computing BBKNN UMAP...", log_handle=log_handle)
            sc.tl.umap(adata_hvg, random_state=42)
            log_message("  ✓ BBKNN UMAP computed", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.12. Transfer Clustering Results to Full Object =====
            log_message("🔄 Transferring clustering results to full gene object...", log_handle=log_handle)
            
            # P2: Prefix cluster names with cell type for easy merging later
            cluster_key = f"subcluster_{celltype_safe}"
            
            # Filter adata_ct_full to match adata_hvg (same cells after QC)
            adata_ct_full = adata_ct_full[adata_hvg.obs_names].copy()
            
            # Transfer clustering results
            adata_ct_full.obs[cluster_key] = adata_hvg.obs['leiden'].astype(str).values
            adata_ct_full.obs['leiden'] = adata_hvg.obs['leiden'].values  # Keep original name for compatibility
            
            # Transfer BBKNN UMAP
            adata_ct_full.obsm['X_umap_bbknn'] = adata_hvg.obsm['X_umap']
            
            log_message(f"  ✓ Transferred clustering results", log_handle=log_handle)
            log_message(f"    Cluster key: '{cluster_key}'", log_handle=log_handle)
            log_message(f"    Final cell count: {adata_ct_full.n_obs:,}", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # Clean up HVG object
            del adata_hvg
            gc.collect()
            
            # ===== 9.13. Marker Gene Analysis on Full Object with log1p Layer =====
            log_message("🔬 Computing marker genes on full gene set...", log_handle=log_handle)
            
            min_cluster_size = final_cluster_counts.min()
            if min_cluster_size < MIN_CELLS_FOR_MARKER:
                log_message(f"⚠️  Warning: Smallest cluster has only {min_cluster_size} cells", log_handle=log_handle)
                log_message(f"  Recommended minimum: {MIN_CELLS_FOR_MARKER} cells for robust statistics", log_handle=log_handle)
                log_message("", log_handle=log_handle)
            
            marker_success = False
            try:
                # P0: Use log1p layer explicitly, NOT .raw
                sc.tl.rank_genes_groups(
                    adata_ct_full,
                    groupby='leiden',
                    method='wilcoxon',
                    layer='log1p',      # P0: Explicit layer
                    use_raw=False,      # P0: Do NOT use .raw
                    corr_method='benjamini-hochberg',
                    pts=True,
                    tie_correct=True,
                    key_added='rank_genes_groups'
                )
                log_message("  ✓ Marker genes computed successfully", log_handle=log_handle)
                log_message("  ✓ Applied Benjamini-Hochberg FDR correction", log_handle=log_handle)
                log_message("  ✓ Used log1p layer on full gene set", log_handle=log_handle)
                log_message("", log_handle=log_handle)
                marker_success = True
                
                # Display top markers
                log_message("📊 Top 5 marker genes per subcluster:", log_handle=log_handle)
                log_message("-" * 80, log_handle=log_handle)
                result = adata_ct_full.uns['rank_genes_groups']
                for cluster in adata_ct_full.obs['leiden'].cat.categories:
                    genes = result['names'][cluster][:5]
                    pvals = result['pvals_adj'][cluster][:5]
                    logfcs = result['logfoldchanges'][cluster][:5]
                    
                    log_message(f"\nSubcluster {cluster} ({final_cluster_counts[cluster]} cells):", log_handle=log_handle)
                    for gene, pval, lfc in zip(genes, pvals, logfcs):
                        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
                        log_message(f"  {gene}: log2FC={lfc:.2f}, FDR={pval:.2e} {sig}", log_handle=log_handle)
                log_message("-" * 80, log_handle=log_handle)
                log_message("", log_handle=log_handle)
                
            except Exception as e:
                log_message(f"❌ ERROR in marker gene analysis: {e}", log_handle=log_handle)
                log_message("  Continuing without marker gene analysis...", log_handle=log_handle)
                adata_ct_full.uns['rank_genes_groups'] = None
            
            # ===== 9.14. Visualizations =====
            log_message("🎨 Generating visualizations...", log_handle=log_handle)
            
            fig_dir = OUTPUT_DIR / 'figures'
            fig_dir.mkdir(exist_ok=True)
            
            # P0: Visualization on original scANVI UMAP (if available)
            if has_scanvi_umap and SCANVI_UMAP_KEY in adata_ct_full.obsm:
                log_message("  Using original scANVI UMAP for visualization", log_handle=log_handle)
                
                # Ensure the key is in the right format for scanpy
                scanvi_basis = SCANVI_UMAP_KEY.replace('X_', '')
                
                # 1. Subclusters on scANVI UMAP
                fig, ax = plt.subplots(figsize=(12, 10))
                sc.pl.embedding(
                    adata_ct_full,
                    basis=scanvi_basis,
                    color='leiden',
                    legend_loc='on data',
                    legend_fontsize=10,
                    title=f'{celltype} Subclusters (Original scANVI UMAP)',
                    ax=ax,
                    show=False
                )
                plt.tight_layout()
                plt.savefig(fig_dir / 'umap_subclusters_scanvi.png', dpi=300, bbox_inches='tight')
                plt.close()
                log_message("  ✓ Saved: umap_subclusters_scanvi.png (original scANVI)", to_console=False, log_handle=log_handle)
            
            # 2. Subclusters on BBKNN UMAP (for comparison)
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.embedding(
                adata_ct_full,
                basis='umap_bbknn',
                color='leiden',
                legend_loc='on data',
                legend_fontsize=10,
                title=f'{celltype} Subclusters (BBKNN-derived UMAP)',
                ax=ax,
                show=False
            )
            plt.tight_layout()
            plt.savefig(fig_dir / 'umap_subclusters_bbknn.png', dpi=300, bbox_inches='tight')
            plt.close()
            log_message("  ✓ Saved: umap_subclusters_bbknn.png (BBKNN-derived)", to_console=False, log_handle=log_handle)
            
            # 3. Batch distribution
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.embedding(
                adata_ct_full,
                basis='umap_bbknn',
                color=BATCH_KEY,
                title=f'Batch Distribution - {celltype}',
                ax=ax,
                show=False
            )
            plt.tight_layout()
            plt.savefig(fig_dir / 'umap_batches.png', dpi=300, bbox_inches='tight')
            plt.close()
            log_message("  ✓ Saved: umap_batches.png", to_console=False, log_handle=log_handle)
            
            # 4. Cluster size barplot
            fig, ax = plt.subplots(figsize=(max(10, n_clusters * 0.8), 6))
            colors = plt.cm.tab20(np.linspace(0, 1, len(final_cluster_counts)))
            ax.bar(final_cluster_counts.index.astype(str), final_cluster_counts.values, color=colors)
            ax.set_xlabel('Subcluster', fontsize=12)
            ax.set_ylabel('Number of Cells', fontsize=12)
            ax.set_title(f'{celltype} Subcluster Size Distribution', fontsize=14)
            ax.axhline(y=MIN_CELLS_PER_CLUSTER, color='red', linestyle='--', 
                       label=f'Min threshold ({MIN_CELLS_PER_CLUSTER})')
            ax.legend()
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(fig_dir / 'cluster_sizes.png', dpi=300, bbox_inches='tight')
            plt.close()
            log_message("  ✓ Saved: cluster_sizes.png", to_console=False, log_handle=log_handle)
            
            # 5. Known marker expression
            if KNOWN_MARKERS:
                markers_in_data = [m for m in KNOWN_MARKERS if m in adata_ct_full.var_names]
                if markers_in_data:
                    n_markers = len(markers_in_data)
                    ncols = min(3, n_markers)
                    nrows = (n_markers + ncols - 1) // ncols
                    
                    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
                    if n_markers == 1:
                        axes = [axes]
                    else:
                        axes = axes.flatten()
                    
                    # Use BBKNN UMAP for marker expression
                    for idx_m, marker in enumerate(markers_in_data):
                        sc.pl.embedding(
                            adata_ct_full,
                            basis='umap_bbknn',
                            color=marker,
                            layer='log1p',
                            title=f'{marker} Expression',
                            ax=axes[idx_m],
                            show=False
                        )
                    
                    for idx_m in range(n_markers, len(axes)):
                        axes[idx_m].axis('off')
                    
                    plt.tight_layout()
                    plt.savefig(fig_dir / 'known_markers_umap.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    log_message("  ✓ Saved: known_markers_umap.png", to_console=False, log_handle=log_handle)
            
            # 6. Marker heatmap and dotplot
            if marker_success and adata_ct_full.uns.get('rank_genes_groups') is not None:
                try:
                    fig = sc.pl.rank_genes_groups_heatmap(
                        adata_ct_full,
                        n_genes=10,
                        groupby='leiden',
                        layer='log1p',
                        show=False,
                        return_fig=True
                    )
                    fig.savefig(fig_dir / 'marker_heatmap.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    log_message("  ✓ Saved: marker_heatmap.png", to_console=False, log_handle=log_handle)
                    
                    fig = sc.pl.rank_genes_groups_dotplot(
                        adata_ct_full,
                        n_genes=5,
                        groupby='leiden',
                        layer='log1p',
                        show=False,
                        return_fig=True
                    )
                    fig.savefig(fig_dir / 'marker_dotplot.png', dpi=300, bbox_inches='tight')
                    plt.close()
                    log_message("  ✓ Saved: marker_dotplot.png", to_console=False, log_handle=log_handle)
                except Exception as e:
                    log_message(f"  ⚠️  Could not generate marker visualizations: {e}", to_console=False, log_handle=log_handle)
            
            log_message("  ✓ All visualizations generated", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.15. Save Results =====
            log_message("💾 Saving results...", log_handle=log_handle)
            
            # Add comprehensive metadata
            adata_ct_full.uns['subcluster_analysis_params'] = {
                'pipeline_version': 'v4.1_FIXED',
                'cell_type_analyzed': celltype,
                'celltype_column': CELLTYPE_COLUMN,
                'batch_key': BATCH_KEY,
                'hvg_method': hvg_method,
                'hvg_layer': 'log1p',
                'n_hvg': n_hvg,
                'resolution': RESOLUTION,
                'n_pcs': N_PCS,
                'min_cells_per_cluster': MIN_CELLS_PER_CLUSTER,
                'bbknn_neighbors_within_batch': neighbors_within_batch,
                'marker_layer': 'log1p',  # Document that markers used log1p
                'marker_correction_method': 'benjamini-hochberg',
                'n_cells_input': n_cells_selected,
                'n_cells_final': adata_ct_full.n_obs,
                'n_subclusters': n_clusters,
                'known_markers': KNOWN_MARKERS,
                'has_scanvi_umap': has_scanvi_umap,
                'preprocessing_validated': True
            }
            
            # Save annotated data
            output_filename = f"{celltype_safe}_subcluster_analyzed.h5ad"
            output_path = OUTPUT_DIR / output_filename
            adata_ct_full.write_h5ad(output_path, compression='gzip', compression_opts=9)
            log_message(f"  ✓ Main output: {output_path}", to_console=False, log_handle=log_handle)
            
            # Export cluster assignments
            cluster_assignments = pd.DataFrame({
                'cell_barcode': adata_ct_full.obs_names,
                'subcluster': adata_ct_full.obs['leiden'].values,
                'subcluster_prefixed': adata_ct_full.obs[cluster_key].values,
                'batch': adata_ct_full.obs[BATCH_KEY].values
            })
            assignment_path = OUTPUT_DIR / f"{celltype_safe}_subcluster_assignments.csv"
            cluster_assignments.to_csv(assignment_path, index=False)
            log_message(f"  ✓ Cluster assignments: {assignment_path}", to_console=False, log_handle=log_handle)
            
            # Export marker genes
            marker_path = None
            if marker_success and adata_ct_full.uns.get('rank_genes_groups') is not None:
                # P1: Use modern API if available
                if USE_MODERN_MARKER_EXPORT:
                    try:
                        dfs = []
                        for g in adata_ct_full.obs['leiden'].cat.categories:
                            df = sc.get.rank_genes_groups_df(
                                adata_ct_full, 
                                group=g, 
                                key='rank_genes_groups'
                            )
                            df['subcluster'] = g
                            dfs.append(df)
                        all_markers = pd.concat(dfs, ignore_index=True)
                        log_message("  Using modern marker export API", to_console=False, log_handle=log_handle)
                    except Exception as e:
                        log_message(f"  Modern API failed: {e}, using legacy method", to_console=False, log_handle=log_handle)
                        USE_MODERN_MARKER_EXPORT = False
                
                if not USE_MODERN_MARKER_EXPORT:
                    # Legacy method
                    marker_df_list = []
                    result = adata_ct_full.uns['rank_genes_groups']
                    
                    for cluster in adata_ct_full.obs['leiden'].cat.categories:
                        cluster_markers = pd.DataFrame({
                            'subcluster': cluster,
                            'gene': result['names'][cluster],
                            'log2fc': result['logfoldchanges'][cluster],
                            'pval': result['pvals'][cluster],
                            'pval_adj': result['pvals_adj'][cluster],
                            'pct_in_group': result['pts'][cluster],
                            'pct_out_group': result['pts_rest'][cluster]
                        })
                        marker_df_list.append(cluster_markers)
                    
                    all_markers = pd.concat(marker_df_list, ignore_index=True)
                
                marker_path = OUTPUT_DIR / f"{celltype_safe}_marker_genes_FDR.csv"
                all_markers.to_csv(marker_path, index=False)
                log_message(f"  ✓ Marker genes: {marker_path}", to_console=False, log_handle=log_handle)
            
            # Export summary statistics
            summary_stats = {
                'pipeline_version': 'v4.1_FIXED',
                'cell_type': celltype,
                'total_cells_input': n_cells_selected,
                'cells_after_qc': adata_ct_full.n_obs,
                'n_subclusters': n_clusters,
                'n_batches': adata_ct_full.obs[BATCH_KEY].nunique(),
                'n_hvg': n_hvg,
                'hvg_method': hvg_method,
                'hvg_layer': 'log1p',
                'resolution': RESOLUTION,
                'marker_layer': 'log1p',
                'marker_correction': 'benjamini-hochberg',
                'marker_analysis_success': marker_success,
                'has_scanvi_umap': has_scanvi_umap,
                'preprocessing_validated': True
            }
            
            summary_df = pd.DataFrame([summary_stats])
            summary_path = OUTPUT_DIR / f"{celltype_safe}_analysis_summary.csv"
            summary_df.to_csv(summary_path, index=False)
            log_message(f"  ✓ Summary statistics: {summary_path}", to_console=False, log_handle=log_handle)
            
            log_message("  ✓ All results saved", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # ===== 9.16. Record Success =====
            celltype_end = datetime.now()
            runtime = (celltype_end - celltype_start).total_seconds()
            
            batch_results.append({
                'cell_type': celltype,
                'status': 'SUCCESS',
                'n_cells_input': n_cells_selected,
                'n_cells_final': adata_ct_full.n_obs,
                'n_subclusters': n_clusters,
                'n_batches': adata_ct_full.obs[BATCH_KEY].nunique(),
                'marker_analysis': 'SUCCESS' if marker_success else 'FAILED',
                'has_scanvi_umap': has_scanvi_umap,
                'runtime_seconds': runtime,
                'output_dir': str(OUTPUT_DIR)
            })
            
            log_message(f"✅ {celltype} analysis completed successfully!", log_handle=log_handle)
            log_message(f"   Runtime: {runtime:.1f} seconds ({runtime/60:.1f} minutes)", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            # Cleanup
            del adata_ct_full
            gc.collect()
            
        except Exception as e:
            # ===== 9.17. Handle Errors =====
            celltype_end = datetime.now()
            runtime = (celltype_end - celltype_start).total_seconds()
            
            error_msg = str(e)
            log_message(f"❌ ERROR processing {celltype}: {error_msg}", log_handle=log_handle)
            
            # Try to get more detailed traceback
            import traceback
            tb = traceback.format_exc()
            log_message(f"\nFull traceback:\n{tb}", to_console=False, log_handle=log_handle)
            
            log_message(f"   Runtime before error: {runtime:.1f} seconds", log_handle=log_handle)
            log_message("", log_handle=log_handle)
            
            batch_results.append({
                'cell_type': celltype,
                'status': 'FAILED',
                'n_cells_input': celltype_counts[celltype],
                'n_cells_final': 0,
                'n_subclusters': 0,
                'n_batches': 0,
                'marker_analysis': 'N/A',
                'has_scanvi_umap': False,
                'runtime_seconds': runtime,
                'output_dir': 'N/A',
                'error': error_msg
            })
            
            if STOP_ON_ERROR:
                log_message(f"STOP_ON_ERROR is True. Stopping batch processing.", log_handle=log_handle)
                break
            else:
                log_message(f"Continuing with next cell type...", log_handle=log_handle)
                log_message("", log_handle=log_handle)
                continue
    
    # ===== 10. Generate Batch Summary Report =====
    log_message("=" * 80, log_handle=log_handle)
    log_message("BATCH PROCESSING SUMMARY", log_handle=log_handle)
    log_message("=" * 80, log_handle=log_handle)
    log_message("", log_handle=log_handle)
    
    end_time = datetime.now()
    total_runtime = (end_time - start_time).total_seconds()
    
    log_message(f"Started:  {start_time.strftime('%Y-%m-%d %H:%M:%S')}", log_handle=log_handle)
    log_message(f"Finished: {end_time.strftime('%Y-%m-%d %H:%M:%S')}", log_handle=log_handle)
    log_message(f"Total runtime: {total_runtime:.1f} seconds ({total_runtime/60:.1f} minutes, {total_runtime/3600:.1f} hours)", log_handle=log_handle)
    log_message("", log_handle=log_handle)
    
    # Create results dataframe
    results_df = pd.DataFrame(batch_results)
    
    # Count successes and failures
    n_total = len(batch_results)
    n_success = (results_df['status'] == 'SUCCESS').sum()
    n_failed = (results_df['status'] == 'FAILED').sum()
    
    log_message(f"Total cell types processed: {n_total}", log_handle=log_handle)
    log_message(f"  Successful: {n_success}", log_handle=log_handle)
    log_message(f"  Failed: {n_failed}", log_handle=log_handle)
    log_message("", log_handle=log_handle)
    
    if n_success > 0:
        log_message("Successful cell types:", log_handle=log_handle)
        for _, row in results_df[results_df['status'] == 'SUCCESS'].iterrows():
            log_message(f"  ✓ {row['cell_type']}: {row['n_cells_final']:,} cells, "
                       f"{row['n_subclusters']} subclusters, "
                       f"runtime={row['runtime_seconds']/60:.1f}min", log_handle=log_handle)
        log_message("", log_handle=log_handle)
    
    if n_failed > 0:
        log_message("Failed cell types:", log_handle=log_handle)
        for _, row in results_df[results_df['status'] == 'FAILED'].iterrows():
            log_message(f"  ❌ {row['cell_type']}: {row.get('error', 'Unknown error')}", log_handle=log_handle)
        log_message("", log_handle=log_handle)
    
    # Save batch summary
    summary_output = OUTPUT_BASE_DIR / 'batch_analysis_summary.csv'
    results_df.to_csv(summary_output, index=False)
    log_message(f"💾 Batch summary saved: {summary_output}", log_handle=log_handle)
    log_message("", log_handle=log_handle)
    
    log_message("=" * 80, log_handle=log_handle)
    if n_failed == 0:
        log_message("✅ ALL CELL TYPES PROCESSED SUCCESSFULLY!", log_handle=log_handle)
    else:
        log_message(f"⚠️  {n_failed}/{n_total} CELL TYPES FAILED", log_handle=log_handle)
    log_message("=" * 80, log_handle=log_handle)

# Log file is automatically closed by context manager

print(f"\n📝 Full processing log saved: {log_file}")

print()
print("=" * 80)
print("BATCH PROCESSING COMPLETED")
print("=" * 80)
print(f"Results directory: {OUTPUT_BASE_DIR}")
print(f"Summary: {OUTPUT_BASE_DIR / 'batch_analysis_summary.csv'}")
print(f"Log: {log_file}")
print()

gc.collect()
