# %% [markdown]
# # Epithelial Cell Subcluster Analysis with BBKNN
# 
# **Purpose:**
# - Perform subcluster analysis on epithelial cells
# - Use BBKNN for batch correction (batch_key='dataset')
# - Visualize subclusters on original UMAP
# - Identify marker genes for each subcluster
# 
# **Author:** r2end  
# **Date:** 2025-01-18  
# **Version:** 1.0
# 
# **Key Features:**
# - Memory-optimized workflow with HVG subset
# - Dynamic resolution adjustment based on cluster size
# - BBKNN batch integration
# - Original UMAP preservation
# - Comprehensive marker gene analysis

# %% [markdown]
# ## 📋 Configuration Parameters

# %%
import os
import sys
import gc
import warnings
from pathlib import Path
from datetime import datetime
import time

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

# Suppress warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# ============================================================================
# INPUT/OUTPUT CONFIGURATION
# ============================================================================

# --- Input Data ---
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"

# --- Output Directory ---
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/epithelial_subcluster_bbknn")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# BATCH CORRECTION CONFIGURATION
# ============================================================================

BATCH_KEY = 'dataset'                    # Column name for batch information
NEIGHBORS_WITHIN_BATCH = 5               # BBKNN parameter (3-5 for small batches)
N_PCS = 50                               # Number of PCs for BBKNN

# ============================================================================
# CLUSTERING CONFIGURATION
# ============================================================================

# --- Cell Type Column ---
# Which column to use for major cell type grouping?
CELLTYPE_COLUMN = 'celltypist_pred'      # Options: 'celltypist_pred', 'cell_type_L1', 'Manual_Annotation'

# --- Resolution Settings ---
# Base resolution for subclustering (will be dynamically adjusted)
RESOLUTION_BASE = 0.5

# Dynamic resolution adjustment based on cell counts
def calculate_resolution(n_cells, base_res=0.5):
    """
    Calculate resolution based on cluster size
    
    Principle:
    - Small clusters: lower resolution (avoid over-splitting noise)
    - Large clusters: higher resolution (discover finer subpopulations)
    """
    if n_cells < 200:
        return base_res * 0.4      # 0.2
    elif n_cells < 500:
        return base_res * 0.6      # 0.3
    elif n_cells < 1000:
        return base_res            # 0.5
    elif n_cells < 5000:
        return base_res * 1.6      # 0.8
    else:
        return base_res * 2.0      # 1.0

# --- Minimum Cells for Subclustering ---
MIN_CELLS_FOR_SUBCLUSTER = 100           # Skip clusters with fewer cells

# ============================================================================
# MARKER GENE CONFIGURATION
# ============================================================================

# Marker gene identification parameters
MARKER_METHOD = 'wilcoxon'               # Options: 'wilcoxon', 't-test', 'logreg'
MARKER_MIN_PCT = 0.25                    # Minimum fraction of cells expressing the gene
MARKER_LOGFC_THRESHOLD = 0.5             # Minimum log fold-change
MARKER_PADJ_THRESHOLD = 0.05             # Adjusted p-value threshold

# Top N genes to report/visualize
TOP_N_MARKERS = 20

# ============================================================================
# VISUALIZATION SETTINGS
# ============================================================================

FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"                    # Options: "pdf", "png"

# UMAP plot settings
UMAP_SIZE = 5                            # Point size
UMAP_ALPHA = 0.6                         # Transparency

# ============================================================================
# PERFORMANCE SETTINGS
# ============================================================================

N_JOBS = 48                              # Parallel processing threads
MEMORY_LIMIT_GB = 48                     # Memory constraint

# ============================================================================
# APPLY SETTINGS
# ============================================================================

sc.settings.n_jobs = N_JOBS
sc.settings.figdir = OUTPUT_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

print("="*80)
print("Epithelial Subcluster Analysis Pipeline (BBKNN)")
print("="*80)
print(f"Input:      {INPUT_H5AD}")
print(f"Output:     {OUTPUT_DIR}")
print(f"Batch key:  {BATCH_KEY}")
print(f"Cell type:  {CELLTYPE_COLUMN}")
print(f"Date:       {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("="*80)

# %% [markdown]
# ## 📂 Step 1: Load Data

# %%
print("\n" + "="*80)
print("STEP 1: Loading Data")
print("="*80)

start_time = time.time()

print(f"\nLoading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

elapsed = time.time() - start_time
print(f"✓ Loaded in {elapsed:.1f}s")

print(f"\n📊 Dataset Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# %%
# Check data structure
print(f"\n🔍 Data Structure Check:")
print(f"  .X shape: {adata.X.shape}")
print(f"  .X type: {type(adata.X)}")
print(f"  .layers: {list(adata.layers.keys())}")
print(f"  .raw: {'✓ Present' if adata.raw is not None else '✗ Missing'}")

if adata.raw is not None:
    print(f"  .raw.X shape: {adata.raw.X.shape}")
    print(f"  .raw genes: {adata.raw.n_vars:,}")

# Check for UMAP
print(f"\n🗺️  Embedding Check:")
print(f"  Available: {list(adata.obsm.keys())}")

if 'X_umap' not in adata.obsm:
    print(f"  ⚠️  WARNING: X_umap not found!")
    if 'X_umap_scvi' in adata.obsm:
        print(f"  → Using X_umap_scvi as X_umap")
        adata.obsm['X_umap'] = adata.obsm['X_umap_scvi']
    else:
        print(f"  → Will compute UMAP after BBKNN")

# %%
# Check batch and cell type information
print(f"\n🏷️  Annotation Check:")

if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found in adata.obs!")

print(f"✓ Batch key: '{BATCH_KEY}'")
n_batches = adata.obs[BATCH_KEY].nunique()
print(f"  Number of batches: {n_batches}")
print(f"  Batch distribution:")
batch_counts = adata.obs[BATCH_KEY].value_counts()
for batch, count in batch_counts.items():
    print(f"    {batch}: {count:,} cells")

if CELLTYPE_COLUMN not in adata.obs.columns:
    raise ValueError(f"Cell type column '{CELLTYPE_COLUMN}' not found!")

print(f"\n✓ Cell type column: '{CELLTYPE_COLUMN}'")
celltype_counts = adata.obs[CELLTYPE_COLUMN].value_counts()
print(f"  Number of cell types: {len(celltype_counts)}")
print(f"\n  Cell type distribution:")
for ct, count in celltype_counts.items():
    print(f"    {ct}: {count:,} cells")

# %% [markdown]
# ## 🔍 Step 2: Preserve Original UMAP

# %%
print("\n" + "="*80)
print("STEP 2: Preserving Original UMAP")
print("="*80)

# Save original UMAP coordinates
if 'X_umap' in adata.obsm:
    adata.obsm['X_umap_original'] = adata.obsm['X_umap'].copy()
    print(f"✓ Original UMAP saved to .obsm['X_umap_original']")
else:
    print(f"⚠️  No UMAP coordinates found - will compute after BBKNN")

# Also preserve any other embeddings
for key in ['X_pca', 'X_scvi', 'X_scanvi']:
    if key in adata.obsm:
        adata.obsm[f'{key}_original'] = adata.obsm[key].copy()
        print(f"✓ Saved {key} to {key}_original")

# %% [markdown]
# ## 🧬 Step 3: Data Preparation for Subclustering

# %%
print("\n" + "="*80)
print("STEP 3: Data Preparation")
print("="*80)

# Check if we need to ensure proper data format
print(f"\n🔍 Checking data format...")

# ⭐ CRITICAL: Ensure we have counts layer (required for HVG and marker analysis)
if 'counts' not in adata.layers:
    print(f"⚠️  'counts' layer not found")
    # If .X is raw counts, save to counts layer
    if adata.X.max() > 100:  # Likely raw counts
        print(f"  → .X appears to be raw counts, saving to layer")
        adata.layers['counts'] = adata.X.copy()
    else:
        print(f"  ⚠️  WARNING: .X appears to be normalized")
        print(f"  → Cannot recover raw counts - marker analysis may be affected")
        print(f"  → Consider using input file with counts layer")

# Ensure we have log1p normalized data in .X
if 'log1p' in adata.layers:
    print(f"✓ Found log1p layer, using as .X")
    adata.X = adata.layers['log1p'].copy()
else:
    print(f"⚠️  log1p layer not found")
    # Check if .X is already log-normalized
    if adata.X.max() < 20:  # Likely log-transformed
        print(f"  → .X appears to be log-normalized already")
    else:
        print(f"  → Performing log-normalization")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

print(f"\n✓ Data format verified")

# %% [markdown]
# ## 🔬 Step 4: Subcluster Analysis Loop

# %%
print("\n" + "="*80)
print("STEP 4: Subclustering Each Cell Type")
print("="*80)

# Get unique cell types
celltypes = adata.obs[CELLTYPE_COLUMN].unique()
n_celltypes = len(celltypes)

print(f"\n📊 Will process {n_celltypes} cell types:")
for ct in celltypes:
    n_cells = (adata.obs[CELLTYPE_COLUMN] == ct).sum()
    print(f"  - {ct}: {n_cells:,} cells")

# Initialize subcluster labels
adata.obs['subcluster'] = 'unprocessed'
adata.obs['subcluster_leiden'] = 'unprocessed'

# Store marker genes for each subcluster
all_markers = {}

# %%
# Main subclustering loop
for idx, celltype in enumerate(celltypes, 1):
    print(f"\n{'='*80}")
    print(f"Processing Cell Type {idx}/{n_celltypes}: {celltype}")
    print(f"{'='*80}")
    
    # 1. Subset to current cell type
    mask = adata.obs[CELLTYPE_COLUMN] == celltype
    adata_sub = adata[mask].copy()
    n_cells = adata_sub.n_obs
    
    print(f"\n📊 Subset Statistics:")
    print(f"  Total cells: {n_cells:,}")
    print(f"  Total genes: {adata_sub.n_vars:,}")
    
    # Check batches in subset
    if BATCH_KEY in adata_sub.obs.columns:
        n_batches_sub = adata_sub.obs[BATCH_KEY].nunique()
        print(f"  Batches: {n_batches_sub}")
        batch_sizes = adata_sub.obs[BATCH_KEY].value_counts()
        print(f"  Batch distribution:")
        for batch, count in batch_sizes.items():
            print(f"    {batch}: {count:,} cells")
    
    # 2. Check if enough cells
    if n_cells < MIN_CELLS_FOR_SUBCLUSTER:
        print(f"\n⚠️  Too few cells ({n_cells} < {MIN_CELLS_FOR_SUBCLUSTER}) - skipping")
        adata.obs.loc[mask, 'subcluster'] = f"{celltype}_0"
        adata.obs.loc[mask, 'subcluster_leiden'] = "0"
        
        # Clean up
        del adata_sub
        gc.collect()
        continue
    
    # 3. ⭐ MEMORY OPTIMIZATION: Work with HVG subset
    print(f"\n🧬 Selecting highly variable genes...")
    
    # Select HVG (use batch-aware if multiple batches)
    try:
        if n_batches_sub > 1 and 'counts' in adata_sub.layers:
            sc.pp.highly_variable_genes(
                adata_sub,
                layer='counts',
                n_top_genes=2000,
                batch_key=BATCH_KEY,
                flavor='seurat_v3',  # ⭐ Recommended by QUICK_REFERENCE_MEMORY
                subset=False
            )
            hvg_method = "batch-aware (seurat_v3)"
        else:
            sc.pp.highly_variable_genes(
                adata_sub,
                layer='counts',
                n_top_genes=2000,
                flavor='seurat_v3',  # ⭐ Recommended by QUICK_REFERENCE_MEMORY
                subset=False
            )
            hvg_method = "non-batch-aware (seurat_v3)"
        print(f"✓ HVG selected ({hvg_method})")
    except Exception as e:
        print(f"⚠️  Batch-aware HVG failed: {e}")
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=2000, subset=False)
        hvg_method = "fallback"
    
    n_hvg = adata_sub.var['highly_variable'].sum()
    print(f"  Highly variable genes: {n_hvg:,}")
    
    # 4. ⭐ CRITICAL: Save full genes to .raw before subsetting
    if adata_sub.raw is None:
        print(f"\n💾 Saving full gene matrix to .raw...")
        if 'counts' in adata_sub.layers:
            adata_sub.raw = sc.AnnData(
                X=adata_sub.layers['counts'],
                obs=adata_sub.obs.copy(),
                var=adata_sub.var.copy()
            )
        else:
            adata_sub.raw = sc.AnnData(
                X=adata_sub.X.copy(),
                obs=adata_sub.obs.copy(),
                var=adata_sub.var.copy()
            )
        print(f"✓ Saved {adata_sub.raw.n_vars:,} genes to .raw")
    
    # 5. Subset to HVG for processing
    print(f"\n🔧 Subsetting to HVG for processing...")
    adata_sub = adata_sub[:, adata_sub.var['highly_variable']].copy()
    print(f"✓ Working with {adata_sub.n_vars:,} HVG")
    
    # 6. Preprocessing for BBKNN
    print(f"\n🔬 Preprocessing for BBKNN...")
    
    # ⚠️ IMPORTANT: Do NOT scale before BBKNN/PCA
    # BBKNN works on log1p normalized data, not scaled data
    # Scaling should only be done for visualization, not for graph construction
    
    # ⭐ PCA on log1p normalized data (NOT scaled)
    print(f"  → Running PCA (n_comps={N_PCS})...")
    sc.tl.pca(adata_sub, n_comps=N_PCS, svd_solver='arpack')
    print(f"✓ PCA complete")
    
    # 7. ⭐ BBKNN batch correction
    print(f"\n🔗 Running BBKNN batch correction...")
    
    if n_batches_sub > 1:
        # Multiple batches - use BBKNN
        print(f"  Batches: {n_batches_sub}")
        print(f"  neighbors_within_batch: {NEIGHBORS_WITHIN_BATCH}")
        
        try:
            sc.external.pp.bbknn(
                adata_sub,
                batch_key=BATCH_KEY,
                neighbors_within_batch=NEIGHBORS_WITHIN_BATCH,
                n_pcs=N_PCS,
                trim=None
            )
            print(f"✓ BBKNN complete")
        except Exception as e:
            print(f"⚠️  BBKNN failed: {e}")
            print(f"  → Falling back to standard neighbors")
            sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=N_PCS)
    else:
        # Single batch - use standard neighbors
        print(f"  Single batch detected - using standard neighbors")
        sc.pp.neighbors(adata_sub, n_neighbors=15, n_pcs=N_PCS)
        print(f"✓ Neighbors graph built")
    
    # 8. Calculate dynamic resolution
    resolution = calculate_resolution(n_cells, RESOLUTION_BASE)
    print(f"\n🎯 Clustering with dynamic resolution:")
    print(f"  Cell count: {n_cells:,}")
    print(f"  Resolution: {resolution:.2f}")
    
    # 9. Leiden clustering
    sc.tl.leiden(
        adata_sub,
        resolution=resolution,
        key_added='leiden_sub'
    )
    
    n_subclusters = adata_sub.obs['leiden_sub'].nunique()
    print(f"\n✓ Found {n_subclusters} subclusters")
    
    # Subcluster sizes
    subcluster_counts = adata_sub.obs['leiden_sub'].value_counts()
    print(f"  Subcluster distribution:")
    for sub_id, count in subcluster_counts.sort_index().items():
        pct = count / n_cells * 100
        print(f"    {sub_id}: {count:,} cells ({pct:.1f}%)")
    
    # 10. Create combined labels
    adata_sub.obs['subcluster'] = celltype + "_" + adata_sub.obs['leiden_sub'].astype(str)
    
    # 11. ⭐ Marker gene identification (use .raw for full genes)
    print(f"\n🧬 Identifying marker genes...")
    
    # Verify .raw exists
    if adata_sub.raw is None:
        print(f"  ⚠️  WARNING: .raw not found for {celltype}")
        print(f"  → Marker analysis will use HVG subset only ({adata_sub.n_vars:,} genes)")
        print(f"  → This may miss important markers outside HVG list")
        use_raw = False
    else:
        print(f"  ✓ Using .raw with {adata_sub.raw.n_vars:,} full genes")
        use_raw = True
    
    try:
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby='leiden_sub',
            method=MARKER_METHOD,
            use_raw=use_raw,  # ⭐ CRITICAL: Use full gene set if available
            pts=True,
            tie_correct=True
        )
        print(f"✓ Marker genes identified")
        
        # Extract top markers
        marker_df_list = []
        
        for cluster_id in adata_sub.obs['leiden_sub'].unique():
            cluster_markers = sc.get.rank_genes_groups_df(
                adata_sub,
                group=cluster_id
            )
            
            # Filter by criteria
            cluster_markers = cluster_markers[
                (cluster_markers['pvals_adj'] < MARKER_PADJ_THRESHOLD) &
                (cluster_markers['logfoldchanges'] > MARKER_LOGFC_THRESHOLD)
            ]
            
            # Add metadata
            cluster_markers['celltype'] = celltype
            cluster_markers['subcluster'] = f"{celltype}_{cluster_id}"
            
            marker_df_list.append(cluster_markers)
        
        # Combine all markers for this cell type
        if marker_df_list:
            celltype_markers = pd.concat(marker_df_list, ignore_index=True)
            all_markers[celltype] = celltype_markers
            print(f"✓ Extracted {len(celltype_markers)} significant markers")
        else:
            print(f"⚠️  No significant markers found for {celltype}")
        
    except Exception as e:
        print(f"⚠️  Marker identification failed: {e}")
        print(f"  Continuing without markers for {celltype}")
    
    # 12. Compute UMAP for this subset (optional, for QC visualization)
    # ⚠️ NOTE: This is ONLY for the subset, not the main adata
    if 'X_umap' not in adata_sub.obsm:
        print(f"\n🗺️  Computing UMAP for subset (QC purpose)...")
        sc.tl.umap(adata_sub)
        print(f"✓ UMAP computed for subset")
    else:
        print(f"\n✓ Subset already has UMAP coordinates")
    
    # 13. ⭐ Write results back to main adata
    print(f"\n💾 Writing results back to main adata...")
    
    # Use index-aligned assignment
    cell_indices = adata_sub.obs_names
    adata.obs.loc[cell_indices, 'subcluster'] = adata_sub.obs['subcluster'].values
    adata.obs.loc[cell_indices, 'subcluster_leiden'] = adata_sub.obs['leiden_sub'].values
    
    print(f"✓ Results written for {n_cells:,} cells")
    
    # 14. Save intermediate result for this cell type
    intermediate_path = OUTPUT_DIR / f"subcluster_{celltype.replace(' ', '_')}.h5ad"
    adata_sub.write_h5ad(intermediate_path, compression='gzip')
    print(f"✓ Saved: {intermediate_path.name}")
    
    # 15. Memory cleanup
    del adata_sub
    gc.collect()
    
    print(f"\n✓ Completed processing for {celltype}")

print(f"\n{'='*80}")
print(f"✓ All cell types processed!")
print(f"{'='*80}")

# %% [markdown]
# ## 📊 Step 5: Summary Statistics

# %%
print("\n" + "="*80)
print("STEP 5: Summary Statistics")
print("="*80)

# Overall subcluster distribution
print(f"\n📊 Overall Subcluster Distribution:")

subcluster_counts = adata.obs['subcluster'].value_counts()
n_subclusters_total = len(subcluster_counts)

print(f"  Total subclusters: {n_subclusters_total}")
print(f"\n  Top 20 subclusters:")
for sub, count in subcluster_counts.head(20).items():
    pct = count / adata.n_obs * 100
    print(f"    {sub}: {count:,} cells ({pct:.1f}%)")

# %%
# Per-celltype summary
print(f"\n📊 Subclusters per Cell Type:")

for celltype in celltypes:
    mask = adata.obs[CELLTYPE_COLUMN] == celltype
    n_cells = mask.sum()
    
    if n_cells > 0:
        subclusters = adata.obs.loc[mask, 'subcluster'].unique()
        n_sub = len(subclusters)
        print(f"\n  {celltype}: {n_cells:,} cells → {n_sub} subclusters")
        
        # Show distribution
        sub_counts = adata.obs.loc[mask, 'subcluster'].value_counts()
        for sub, count in sub_counts.items():
            pct = count / n_cells * 100
            print(f"    {sub}: {count:,} ({pct:.1f}%)")

# %% [markdown]
# ## 🎨 Step 6: Visualization on Original UMAP

# %%
print("\n" + "="*80)
print("STEP 6: Visualization")
print("="*80)

# Restore original UMAP if available
if 'X_umap_original' in adata.obsm:
    print(f"\n✓ Using original UMAP coordinates")
    adata.obsm['X_umap'] = adata.obsm['X_umap_original'].copy()
elif 'X_umap' not in adata.obsm:
    print(f"\n⚠️  No UMAP available, computing new UMAP...")
    # Use PCA from one of the subsets (or recompute)
    if 'X_pca_original' in adata.obsm:
        adata.obsm['X_pca'] = adata.obsm['X_pca_original']
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.umap(adata)
    print(f"✓ UMAP computed")

# %%
# Visualization 1: Subclusters on original UMAP
print(f"\n📊 Plotting subclusters on original UMAP...")

fig, ax = plt.subplots(figsize=(14, 12))

sc.pl.umap(
    adata,
    color='subcluster',
    title='Epithelial Subclusters (BBKNN Integration)',
    legend_loc='right margin',
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    frameon=False,
    ax=ax,
    show=False
)

plt.tight_layout()
plt.savefig(
    sc.settings.figdir / f"umap_subclusters_all.{FIGURE_FORMAT}",
    dpi=FIGURE_DPI,
    bbox_inches='tight'
)
plt.show()
print(f"✓ Saved: umap_subclusters_all.{FIGURE_FORMAT}")

# %%
# Visualization 2: Compare original cell types vs subclusters
print(f"\n📊 Comparing original cell types vs subclusters...")

fig, axes = plt.subplots(1, 2, figsize=(24, 10))

# Original cell types
sc.pl.umap(
    adata,
    color=CELLTYPE_COLUMN,
    title='Original Cell Type Annotations',
    legend_loc='right margin',
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    frameon=False,
    ax=axes[0],
    show=False
)

# Subclusters
sc.pl.umap(
    adata,
    color='subcluster',
    title='Refined Subclusters',
    legend_loc='right margin',
    size=UMAP_SIZE,
    alpha=UMAP_ALPHA,
    frameon=False,
    ax=axes[1],
    show=False
)

plt.tight_layout()
plt.savefig(
    sc.settings.figdir / f"umap_comparison_original_vs_subcluster.{FIGURE_FORMAT}",
    dpi=FIGURE_DPI,
    bbox_inches='tight'
)
plt.show()
print(f"✓ Saved: umap_comparison_original_vs_subcluster.{FIGURE_FORMAT}")

# %%
# Visualization 3: Batch distribution on UMAP
if BATCH_KEY in adata.obs.columns:
    print(f"\n📊 Plotting batch distribution...")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(
        adata,
        color=BATCH_KEY,
        title=f'Batch Distribution ({BATCH_KEY})',
        legend_loc='right margin',
        size=UMAP_SIZE,
        alpha=UMAP_ALPHA,
        frameon=False,
        ax=ax,
        show=False
    )
    
    plt.tight_layout()
    plt.savefig(
        sc.settings.figdir / f"umap_batch_distribution.{FIGURE_FORMAT}",
        dpi=FIGURE_DPI,
        bbox_inches='tight'
    )
    plt.show()
    print(f"✓ Saved: umap_batch_distribution.{FIGURE_FORMAT}")

# %%
# Visualization 4: Per-celltype subcluster UMAPs
print(f"\n📊 Creating per-celltype subcluster plots...")

for celltype in celltypes:
    mask = adata.obs[CELLTYPE_COLUMN] == celltype
    n_cells = mask.sum()
    
    if n_cells < MIN_CELLS_FOR_SUBCLUSTER:
        continue
    
    print(f"  Plotting {celltype}...")
    
    # Create a categorical color map for this celltype's subclusters
    adata_temp = adata[mask].copy()
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sc.pl.umap(
        adata_temp,
        color='subcluster',
        title=f'{celltype} Subclusters',
        legend_loc='right margin',
        size=UMAP_SIZE * 1.5,
        alpha=UMAP_ALPHA,
        frameon=False,
        ax=ax,
        show=False
    )
    
    plt.tight_layout()
    
    celltype_safe = celltype.replace(' ', '_').replace('/', '_')
    plt.savefig(
        sc.settings.figdir / f"umap_subcluster_{celltype_safe}.{FIGURE_FORMAT}",
        dpi=FIGURE_DPI,
        bbox_inches='tight'
    )
    plt.close()
    
    del adata_temp
    gc.collect()

print(f"✓ Per-celltype plots saved")

# %% [markdown]
# ## 📈 Step 7: Marker Gene Visualization

# %%
print("\n" + "="*80)
print("STEP 7: Marker Gene Visualization")
print("="*80)

if not all_markers:
    print(f"⚠️  No marker genes available for visualization")
else:
    print(f"\n📊 Visualizing top marker genes for each subcluster...")
    
    for celltype, markers_df in all_markers.items():
        print(f"\n  Processing {celltype}...")
        
        # Get unique subclusters
        subclusters = markers_df['subcluster'].unique()
        n_subclusters = len(subclusters)
        
        if n_subclusters == 0:
            continue
        
        # Top 5 markers per subcluster for heatmap
        top_markers = []
        for subcluster in subclusters:
            sub_markers = markers_df[markers_df['subcluster'] == subcluster]
            top_genes = sub_markers.nsmallest(5, 'pvals_adj')['names'].tolist()
            top_markers.extend(top_genes)
        
        # Remove duplicates while preserving order
        top_markers = list(dict.fromkeys(top_markers))
        
        # Limit to available genes in adata
        available_markers = [g for g in top_markers if g in adata.var_names]
        
        if len(available_markers) < 3:
            print(f"    ⚠️  Too few markers available for {celltype}")
            continue
        
        # Subset to this celltype
        mask = adata.obs[CELLTYPE_COLUMN] == celltype
        adata_temp = adata[mask].copy()
        
        # Create dotplot
        try:
            fig = plt.figure(figsize=(max(10, len(available_markers) * 0.4), max(6, n_subclusters * 0.5)))
            
            sc.pl.dotplot(
                adata_temp,
                var_names=available_markers[:30],  # Limit to 30 genes
                groupby='subcluster',
                dendrogram=False,
                standard_scale='var',
                title=f'{celltype} - Top Marker Genes',
                show=False
            )
            
            plt.tight_layout()
            
            celltype_safe = celltype.replace(' ', '_').replace('/', '_')
            plt.savefig(
                sc.settings.figdir / f"dotplot_markers_{celltype_safe}.{FIGURE_FORMAT}",
                dpi=FIGURE_DPI,
                bbox_inches='tight'
            )
            plt.close()
            
            print(f"    ✓ Saved dotplot for {celltype}")
            
        except Exception as e:
            print(f"    ⚠️  Could not create dotplot: {e}")
        
        # Create heatmap (top 20 markers)
        try:
            if len(available_markers) > 0:
                fig = plt.figure(figsize=(12, max(8, len(available_markers[:20]) * 0.3)))
                
                sc.pl.heatmap(
                    adata_temp,
                    var_names=available_markers[:20],
                    groupby='subcluster',
                    dendrogram_key='dendrogram_subcluster',
                    standard_scale='var',
                    cmap='RdBu_r',
                    show=False
                )
                
                plt.tight_layout()
                
                plt.savefig(
                    sc.settings.figdir / f"heatmap_markers_{celltype_safe}.{FIGURE_FORMAT}",
                    dpi=FIGURE_DPI,
                    bbox_inches='tight'
                )
                plt.close()
                
                print(f"    ✓ Saved heatmap for {celltype}")
        
        except Exception as e:
            print(f"    ⚠️  Could not create heatmap: {e}")
        
        del adata_temp
        gc.collect()

# %%
# Visualization: Top marker genes on UMAP (selected examples)
print(f"\n📊 Visualizing selected markers on UMAP...")

# Get top 3 markers from first celltype (as example)
if all_markers:
    first_celltype = list(all_markers.keys())[0]
    example_markers = all_markers[first_celltype].nsmallest(6, 'pvals_adj')['names'].tolist()
    
    # Filter to available genes
    example_markers = [g for g in example_markers if g in adata.var_names][:6]
    
    if len(example_markers) > 0:
        print(f"  Example markers: {', '.join(example_markers)}")
        
        fig = plt.figure(figsize=(18, 12))
        
        sc.pl.umap(
            adata,
            color=example_markers,
            ncols=3,
            frameon=False,
            cmap='Reds',
            size=UMAP_SIZE * 0.8,
            alpha=UMAP_ALPHA,
            show=False
        )
        
        plt.tight_layout()
        plt.savefig(
            sc.settings.figdir / f"umap_example_markers.{FIGURE_FORMAT}",
            dpi=FIGURE_DPI,
            bbox_inches='tight'
        )
        plt.show()
        print(f"✓ Saved: umap_example_markers.{FIGURE_FORMAT}")

# %% [markdown]
# ## 💾 Step 8: Save Results

# %%
print("\n" + "="*80)
print("STEP 8: Saving Results")
print("="*80)

# 1. Save annotated adata with subclusters
output_h5ad = OUTPUT_DIR / "epithelial_with_subclusters.h5ad"

print(f"\n💾 Saving annotated data with subclusters...")
adata.write_h5ad(output_h5ad, compression='gzip', compression_opts=9)

file_size = output_h5ad.stat().st_size / 1e9
print(f"✓ Saved: {output_h5ad.name}")
print(f"  Size: {file_size:.2f} GB")

# %%
# 2. Export subcluster annotations to CSV
print(f"\n📄 Exporting subcluster annotations...")

subcluster_df = adata.obs[[
    CELLTYPE_COLUMN,
    'subcluster',
    'subcluster_leiden'
]].copy()

if BATCH_KEY in adata.obs.columns:
    subcluster_df[BATCH_KEY] = adata.obs[BATCH_KEY]

subcluster_csv = OUTPUT_DIR / "subcluster_annotations.csv"
subcluster_df.to_csv(subcluster_csv)

print(f"✓ Saved: {subcluster_csv.name}")

# %%
# 3. Save marker genes
if all_markers:
    print(f"\n📊 Saving marker genes...")
    
    # Combine all markers
    all_markers_df = pd.concat(all_markers.values(), ignore_index=True)
    
    # Save all markers
    markers_all_csv = OUTPUT_DIR / "markers_all_subclusters.csv"
    all_markers_df.to_csv(markers_all_csv, index=False)
    print(f"✓ Saved: {markers_all_csv.name}")
    print(f"  Total markers: {len(all_markers_df):,}")
    
    # Save top N markers per subcluster
    top_markers_list = []
    for subcluster in all_markers_df['subcluster'].unique():
        sub_markers = all_markers_df[all_markers_df['subcluster'] == subcluster]
        top_n = sub_markers.nsmallest(TOP_N_MARKERS, 'pvals_adj')
        top_markers_list.append(top_n)
    
    top_markers_df = pd.concat(top_markers_list, ignore_index=True)
    markers_top_csv = OUTPUT_DIR / f"markers_top{TOP_N_MARKERS}_per_subcluster.csv"
    top_markers_df.to_csv(markers_top_csv, index=False)
    print(f"✓ Saved: {markers_top_csv.name}")
    print(f"  Top markers: {len(top_markers_df):,}")
    
    # Save per-celltype marker files
    for celltype, markers_df in all_markers.items():
        celltype_safe = celltype.replace(' ', '_').replace('/', '_')
        markers_ct_csv = OUTPUT_DIR / f"markers_{celltype_safe}.csv"
        markers_df.to_csv(markers_ct_csv, index=False)
    
    print(f"✓ Saved per-celltype marker files")

# %%
# 4. Save summary statistics
print(f"\n📊 Saving summary statistics...")

summary_file = OUTPUT_DIR / "subcluster_summary.txt"

with open(summary_file, 'w') as f:
    f.write("="*80 + "\n")
    f.write("Epithelial Subcluster Analysis Summary (BBKNN)\n")
    f.write("="*80 + "\n\n")
    
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input: {INPUT_H5AD}\n")
    f.write(f"Batch key: {BATCH_KEY}\n")
    f.write(f"Cell type column: {CELLTYPE_COLUMN}\n\n")
    
    f.write(f"Dataset:\n")
    f.write(f"  Total cells: {adata.n_obs:,}\n")
    f.write(f"  Total genes: {adata.n_vars:,}\n")
    f.write(f"  Batches: {adata.obs[BATCH_KEY].nunique()}\n\n")
    
    f.write(f"Subclustering Results:\n")
    f.write(f"  Total subclusters: {adata.obs['subcluster'].nunique()}\n")
    f.write(f"  Base resolution: {RESOLUTION_BASE}\n")
    f.write(f"  BBKNN neighbors_within_batch: {NEIGHBORS_WITHIN_BATCH}\n\n")
    
    f.write(f"Subcluster Distribution:\n")
    for subcluster, count in adata.obs['subcluster'].value_counts().items():
        pct = count / adata.n_obs * 100
        f.write(f"  {subcluster}: {count:,} cells ({pct:.1f}%)\n")
    
    if all_markers:
        f.write(f"\nMarker Genes:\n")
        f.write(f"  Total significant markers: {len(all_markers_df):,}\n")
        f.write(f"  P-value threshold: {MARKER_PADJ_THRESHOLD}\n")
        f.write(f"  Log FC threshold: {MARKER_LOGFC_THRESHOLD}\n")

print(f"✓ Saved: {summary_file.name}")

# %%
# 5. Create quick reference table
print(f"\n📋 Creating quick reference table...")

ref_table = adata.obs.groupby('subcluster').agg({
    CELLTYPE_COLUMN: 'first',
    'subcluster_leiden': 'first'
}).reset_index()

# Add cell counts
ref_table['n_cells'] = adata.obs['subcluster'].value_counts().loc[ref_table['subcluster']].values
ref_table['pct_total'] = ref_table['n_cells'] / adata.n_obs * 100

# Add top 3 markers if available
if all_markers:
    ref_table['top_markers'] = ''
    for idx, row in ref_table.iterrows():
        subcluster = row['subcluster']
        if subcluster in all_markers_df['subcluster'].values:
            top3 = all_markers_df[
                all_markers_df['subcluster'] == subcluster
            ].nsmallest(3, 'pvals_adj')['names'].tolist()
            ref_table.at[idx, 'top_markers'] = ', '.join(top3)

# Sort by cell type and subcluster
ref_table = ref_table.sort_values([CELLTYPE_COLUMN, 'subcluster_leiden'])

ref_csv = OUTPUT_DIR / "subcluster_reference_table.csv"
ref_table.to_csv(ref_csv, index=False)

print(f"✓ Saved: {ref_csv.name}")
print(f"\n📋 Quick Reference Table Preview:")
print(ref_table.head(15).to_string(index=False))

# %% [markdown]
# ## ✅ Pipeline Complete

# %%
print("\n" + "="*80)
print("PIPELINE COMPLETE")
print("="*80)

print(f"\n📂 Output Directory: {OUTPUT_DIR}")

print(f"\n📄 Generated Files:")
print(f"  1. epithelial_with_subclusters.h5ad         (Main results)")
print(f"  2. subcluster_annotations.csv               (Annotation table)")
print(f"  3. subcluster_reference_table.csv           (Quick reference)")
print(f"  4. subcluster_summary.txt                   (Summary statistics)")

if all_markers:
    print(f"  5. markers_all_subclusters.csv              (All markers)")
    print(f"  6. markers_top{TOP_N_MARKERS}_per_subcluster.csv      (Top markers)")
    print(f"  7. markers_[celltype].csv                   (Per-celltype markers)")

print(f"\n📊 Figures (in figures/ subdirectory):")
print(f"  - umap_subclusters_all.{FIGURE_FORMAT}")
print(f"  - umap_comparison_original_vs_subcluster.{FIGURE_FORMAT}")
print(f"  - umap_batch_distribution.{FIGURE_FORMAT}")
print(f"  - umap_subcluster_[celltype].{FIGURE_FORMAT} (per celltype)")
if all_markers:
    print(f"  - dotplot_markers_[celltype].{FIGURE_FORMAT}")
    print(f"  - heatmap_markers_[celltype].{FIGURE_FORMAT}")
    print(f"  - umap_example_markers.{FIGURE_FORMAT}")

print(f"\n📋 Subclustering Summary:")
print(f"  Total cells: {adata.n_obs:,}")
print(f"  Total subclusters: {adata.obs['subcluster'].nunique()}")
print(f"  Cell types processed: {len(celltypes)}")

print(f"\n✅ All tasks completed successfully!")
print(f"\n💡 Next steps:")
print(f"  1. Review UMAP visualizations to validate subcluster quality")
print(f"  2. Check marker genes for biological relevance")
print(f"  3. Use subcluster_reference_table.csv for annotation refinement")
print(f"  4. Consider manual curation of subcluster labels")

# %%
