#!/usr/bin/env python3
"""
T Cell Subcluster Analysis Pipeline v1.0 (PRODUCTION)
======================================================
Author: r2end
Date: 2025-01-26

Purpose:
--------
Subcluster analysis of T cells following scANVI re-training by Sample.
Follows QUICK_REFERENCE_MEMORY v2.12 best practices and 
universal_celltype_subcluster_pipeline_v4.1 architecture.

Key Differences from Universal Pipeline:
-----------------------------------------
1. Input: T cell subset pre-processed with scANVI by Sample
2. scANVI UMAP: Already available from upstream training
3. Known markers: T cell specific (CD3D, CD4, CD8A, etc.)
4. Resolution: 0.4 (higher than default for T cell heterogeneity)

CRITICAL Requirements:
----------------------
Input data MUST have:
1. Gene filtering (MT/ribo/ENSG removed) - GLOBAL
2. Normalization (all cells together) - GLOBAL  
3. log1p transformation - GLOBAL
4. .raw saved with full gene set
5. obsm['X_scanvi']: scANVI latent from Sample-level training
6. obsm['X_umap_scanvi']: scANVI UMAP (optional)

Workflow:
---------
Pre-normalized T cell data → 
  → Filter cells → Create full & HVG objects →
  → HVG: PCA/BBKNN/Leiden →
  → Transfer clusters to full object →
  → Full: Markers (log1p layer) + Viz (scANVI UMAP) →
  → Save

Input: 
------
- T cell h5ad from scANVI re-training (by Sample):
  * layers['counts']: Raw counts
  * layers['log1p']: Globally normalized log1p
  * .raw: Full gene set
  * obsm['X_scanvi']: scANVI latent
  * obsm['X_umap_scanvi']: scANVI UMAP (if available)
  * obs['cell_type_scanvi_by_sample']: scANVI predictions

Output:
-------
- output_dir/tcell_subcluster/
  ├── adata_tcell_subcluster_analyzed.h5ad
  ├── marker_genes_with_FDR.csv
  ├── cluster_assignments.csv
  ├── analysis_summary.csv
  └── figures/
      ├── umap_subclusters_scanvi.png (on original scANVI)
      ├── umap_subclusters_bbknn.png (BBKNN-derived)
      ├── umap_batches.png
      ├── umap_celltypes.png (scANVI predictions)
      ├── cluster_sizes.png
      ├── known_markers_umap.png
      ├── marker_heatmap.png
      └── marker_dotplot.png

Memory: ~20-30GB (optimized)
Runtime: ~30-60min
"""

# ===== 1. Import Libraries =====
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
print("T CELL SUBCLUSTER ANALYSIS PIPELINE v1.0 (PRODUCTION)")
print("=" * 80)
print()

# ===== 2. Dependency Check =====
print("🔍 Checking dependencies...")
try:
    import bbknn
    print("  ✓ bbknn available")
except ImportError:
    print("  ❌ bbknn not found. Install with: pip install bbknn")
    sys.exit(1)

try:
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

# ===== 3. Configuration =====
# ========================================
# Paths
# ========================================
# ⚠️ UPDATE THIS PATH to your T cell scANVI output
INPUT_PATH = "/home/h2048/data/py/0112/cd4cd8_tcell/adata_tcell_RETRAINED_by_sample.h5ad"
OUTPUT_BASE_DIR = Path("/home/h2048/data/py/0126/tcell_subcluster")
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_DIR = OUTPUT_BASE_DIR / "tcell_subcluster"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ========================================
# Column Names
# ========================================
BATCH_KEY = 'dataset'                              # Batch key (lowercase from scANVI training)
CELLTYPE_COLUMN = 'cell_type_scanvi_by_sample'  # scANVI predictions

# scANVI Latent Space Columns
SCANVI_LATENT_KEY = 'X_scanvi'            # scANVI latent from training
SCANVI_UMAP_KEY = 'X_umap_scanvi'         # scANVI UMAP (may need to compute)

# ========================================
# Quality Control Thresholds
# ========================================
MIN_CELLS_PER_BATCH = 3           # Remove batches smaller than this
MIN_CELLS_PER_CLUSTER = 10        # Remove clusters smaller than this
MIN_CELLS_FOR_MARKER = 20         # Minimum cells for reliable marker analysis

# ========================================
# Analysis Parameters (T cell specific)
# ========================================
N_HVG = 3000                      # HVG for clustering
N_PCS = 50                        # PCs for BBKNN
RESOLUTION = 0.4                  # Leiden resolution (higher for T cell heterogeneity)
NEIGHBORS_WITHIN_BATCH = 5        # BBKNN neighbors
HVG_FLAVOR = 'seurat'             # HVG method

# ========================================
# T Cell Known Markers
# ========================================
KNOWN_MARKERS = [
    # Core T cell
    "CD3D", "CD3E", "CD3G",
    # CD4+ T
    "CD4", "IL7R", "CD40LG",
    # CD8+ T
    "CD8A", "CD8B",
    # Naive/Memory
    "CCR7", "SELL", "TCF7", "LEF1",
    # Activation/Effector
    "GZMB", "GZMK", "PRF1", "IFNG",
    # Regulatory T
    "FOXP3", "IL2RA", "IKZF2",
    # Tissue resident
    "CD69", "ITGAE",
    # NK-like
    "GNLY", "NKG7", "KLRD1"
]

# ========================================
# Scanpy Settings
# ========================================
sc.settings.verbosity = 1
sc.settings.n_jobs = 48
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))

print("Configuration Summary:")
print("-" * 80)
print(f"Input: {INPUT_PATH}")
print(f"Output: {OUTPUT_DIR}")
print(f"Batch key: {BATCH_KEY}")
print(f"Cell type column: {CELLTYPE_COLUMN}")
print(f"HVG: {N_HVG}, Resolution: {RESOLUTION}")
print(f"Known markers: {len(KNOWN_MARKERS)} genes")
print("-" * 80)
print()

# ===== 4. Load Data =====
print("📂 Loading T cell data...")
adata_full = sc.read_h5ad(INPUT_PATH)
print(f"  ✓ Loaded: {adata_full.n_obs:,} cells × {adata_full.n_vars} genes")
print()

# ===== 5. Validate Pre-Processing =====
print("🔍 Validating input data...")

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

# Check gene count
if adata_full.n_vars > 30000:
    print(f"  ⚠️  WARNING: Large gene count ({adata_full.n_vars} genes)")
    print(f"     Suggests unfiltered data")
elif adata_full.n_vars < 15000:
    print(f"  ⚠️  WARNING: Small gene count ({adata_full.n_vars} genes)")
    print(f"     May be over-filtered")
else:
    print(f"  ✓ Gene count reasonable: {adata_full.n_vars} genes")

# Check for scANVI latent space
has_scanvi = SCANVI_LATENT_KEY in adata_full.obsm
has_scanvi_umap = SCANVI_UMAP_KEY in adata_full.obsm

if has_scanvi:
    print(f"  ✓ Found scANVI latent: '{SCANVI_LATENT_KEY}'")
else:
    print(f"  ⚠️  WARNING: scANVI latent '{SCANVI_LATENT_KEY}' not found")
    print(f"     Will only use BBKNN UMAP")

if has_scanvi_umap:
    print(f"  ✓ Found scANVI UMAP: '{SCANVI_UMAP_KEY}'")
else:
    print(f"  ⚠️  scANVI UMAP '{SCANVI_UMAP_KEY}' not found")
    print(f"     Will compute from scANVI latent if available")

print()

# Stop if critical validation failed
if not validation_passed:
    print("❌ INPUT DATA VALIDATION FAILED!")
    print()
    print("Errors found:")
    for error in errors:
        print(error)
    print()
    raise ValueError("Input data validation failed. See errors above.")

print("  ✅ Input data validation PASSED")
print()

# ===== 6. Validate Columns =====
print("🔍 Validating required columns...")

# Check batch column
if BATCH_KEY not in adata_full.obs.columns:
    raise ValueError(
        f"Batch column '{BATCH_KEY}' not found!\n"
        f"Available columns: {list(adata_full.obs.columns)}"
    )
print(f"  ✓ Found batch column: '{BATCH_KEY}'")

# Show batch distribution
batch_counts = adata_full.obs[BATCH_KEY].value_counts()
print(f"\n  Batch distribution:")
print(f"    Total: {len(batch_counts)} batches")
for batch, count in batch_counts.head(5).items():
    print(f"    {batch}: {count:,} cells")
if len(batch_counts) > 5:
    print(f"    ... and {len(batch_counts) - 5} more batches")

# Check cell type column
if CELLTYPE_COLUMN in adata_full.obs.columns:
    print(f"\n  ✓ Found cell type column: '{CELLTYPE_COLUMN}'")
    celltype_counts = adata_full.obs[CELLTYPE_COLUMN].value_counts()
    print(f"    Unique types: {celltype_counts.nunique()}")
    print(f"    Top 10 types:")
    for i, (ct, count) in enumerate(celltype_counts.head(10).items(), 1):
        pct = count / adata_full.n_obs * 100
        print(f"      {i:2d}. {ct:30s} {count:6,} ({pct:5.1f}%)")
else:
    print(f"\n  ℹ️  Cell type column '{CELLTYPE_COLUMN}' not found (optional)")

print()

# ===== 7. Setup Full Gene Object =====
print("🔬 Setting up full gene object...")

# Use log1p layer (globally normalized)
adata_full.X = adata_full.layers['log1p']

# Ensure float32
if adata_full.X.dtype != np.float32:
    print("  Converting to float32 for memory efficiency...")
    adata_full.X = adata_full.X.astype(np.float32)

print(f"  ✓ Using globally normalized log1p data")
print(f"    Data type: {adata_full.X.dtype}")
print()

# ===== 8. Remove Small Batches =====
print("🔍 Checking batch sizes...")

batch_counts = adata_full.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index

if len(small_batches) > 0:
    print(f"⚠️  Removing {len(small_batches)} small batches (< {MIN_CELLS_PER_BATCH} cells):")
    for batch in small_batches[:5]:
        print(f"    {batch}: {batch_counts[batch]} cells")
    if len(small_batches) > 5:
        print(f"    ... and {len(small_batches) - 5} more")
    
    adata_full = adata_full[~adata_full.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  ✓ Remaining: {adata_full.n_obs:,} cells in {adata_full.obs[BATCH_KEY].nunique()} batches")
else:
    print(f"  ✓ All batches meet minimum size requirement")

print()

# ===== 9. HVG Selection =====
print(f"🎯 Selecting top {N_HVG} highly variable genes...")

try:
    sc.pp.highly_variable_genes(
        adata_full,
        layer='log1p',
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        flavor=HVG_FLAVOR,
        subset=False
    )
    hvg_method = f"batch-aware ({HVG_FLAVOR})"
    print(f"  ✓ Batch-aware HVG selection successful")
except Exception as e:
    print(f"  ⚠️  Batch-aware HVG failed: {e}")
    print("  Falling back to non-batch-aware method")
    sc.pp.highly_variable_genes(
        adata_full,
        layer='log1p',
        n_top_genes=N_HVG,
        flavor=HVG_FLAVOR,
        subset=False
    )
    hvg_method = f"non-batch-aware ({HVG_FLAVOR})"

n_hvg = adata_full.var['highly_variable'].sum()
print(f"  ✓ Selected {n_hvg} HVGs ({hvg_method})")

# Force include known markers
known_in_data = [m for m in KNOWN_MARKERS if m in adata_full.var_names]
if known_in_data:
    adata_full.var.loc[known_in_data, 'highly_variable'] = True
    n_hvg = adata_full.var['highly_variable'].sum()
    print(f"  ✓ Force-included {len(known_in_data)} known markers")
    print(f"  Final HVG count: {n_hvg}")

print()

# ===== 10. Create HVG Subset =====
print("✂️  Creating HVG subset for clustering...")

adata_hvg = adata_full[:, adata_full.var['highly_variable']].copy()
print(f"  ✓ HVG subset: {adata_hvg.n_obs:,} cells × {adata_hvg.n_vars} genes")
print()

# ===== 11. PCA =====
print(f"📊 Computing PCA ({N_PCS} components)...")
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack', random_state=42)
print("  ✓ PCA computed")
print()

# ===== 12. BBKNN Integration =====
print("🔗 Running BBKNN integration...")

batch_sizes = adata_hvg.obs[BATCH_KEY].value_counts()
min_batch_size = batch_sizes.min()
neighbors_within_batch = min(NEIGHBORS_WITHIN_BATCH, max(3, min_batch_size - 1))

print(f"  Smallest batch size: {min_batch_size}")
print(f"  Using neighbors_within_batch: {neighbors_within_batch}")

sce.pp.bbknn(
    adata_hvg,
    batch_key=BATCH_KEY,
    neighbors_within_batch=neighbors_within_batch,
    n_pcs=N_PCS,
    trim=None
)
print("  ✓ BBKNN integration completed")
print()

# ===== 13. Leiden Clustering =====
print(f"🎨 Computing Leiden clustering (resolution={RESOLUTION})...")
sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=42)
n_clusters = adata_hvg.obs['leiden'].nunique()
print(f"  ✓ Identified {n_clusters} subclusters")
print()

# Check cluster sizes
cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
print("  Subcluster sizes:")
for cluster, count in cluster_counts.items():
    status = "⚠️ SMALL" if count < MIN_CELLS_PER_CLUSTER else "✓"
    print(f"    Subcluster {cluster}: {count:,} cells {status}")
print()

# ⭐ QUALITY CHECK: Record small clusters BEFORE filtering
# This is for QC reporting (50 cell threshold), separate from filtering (10 cell threshold)
SMALL_CLUSTER_THRESHOLD = 50
small_clusters_qc = []  # For QC reporting
for cluster, count in cluster_counts.items():
    if count < SMALL_CLUSTER_THRESHOLD:
        small_clusters_qc.append((f"T_c{cluster}", count))

if small_clusters_qc:
    print(f"  ℹ️  Note: {len(small_clusters_qc)} clusters have <{SMALL_CLUSTER_THRESHOLD} cells")
    print(f"     (will be flagged in QC report)")
print()

# ===== 14. Filter Small Clusters =====
small_clusters = cluster_counts[cluster_counts < MIN_CELLS_PER_CLUSTER].index.tolist()
if len(small_clusters) > 0:
    print(f"⚠️  Removing {len(small_clusters)} small subclusters (< {MIN_CELLS_PER_CLUSTER} cells):")
    print(f"  Subclusters: {small_clusters}")
    
    adata_hvg = adata_hvg[~adata_hvg.obs['leiden'].isin(small_clusters)].copy()
    
    # Recompute
    print("  Recomputing BBKNN and clustering after filtering...")
    sce.pp.bbknn(
        adata_hvg,
        batch_key=BATCH_KEY,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=N_PCS,
        trim=None
    )
    sc.tl.leiden(adata_hvg, resolution=RESOLUTION, key_added='leiden', random_state=42)
    n_clusters = adata_hvg.obs['leiden'].nunique()
    
    print(f"  ✓ Final: {adata_hvg.n_obs:,} cells in {n_clusters} subclusters")
    
    final_cluster_counts = adata_hvg.obs['leiden'].value_counts().sort_index()
    print("  Final subcluster sizes:")
    for cluster, count in final_cluster_counts.items():
        print(f"    Subcluster {cluster}: {count:,} cells")
    print()
else:
    final_cluster_counts = cluster_counts
    print(f"  ✓ All clusters meet minimum size requirement")
    print()

# ===== 15. Compute BBKNN UMAP =====
print("🗺️  Computing BBKNN UMAP...")
sc.tl.umap(adata_hvg, random_state=42)
print("  ✓ BBKNN UMAP computed")
print()

# ===== 16. Transfer Results to Full Object =====
print("🔄 Transferring clustering results to full gene object...")

# Filter adata_full to match adata_hvg
adata_full = adata_full[adata_hvg.obs_names].copy()

# Transfer clustering with standard naming convention (T_c0, T_c1, etc)
adata_full.obs['leiden'] = adata_hvg.obs['leiden'].values
# Standard format: T_c0, T_c1, T_c2, ...
adata_full.obs['cell_type_L3'] = 'T_c' + adata_hvg.obs['leiden'].astype(str)

# Transfer BBKNN UMAP
adata_full.obsm['X_umap_bbknn'] = adata_hvg.obsm['X_umap']

print(f"  ✓ Transferred clustering results")
print(f"    Naming format: T_c0, T_c1, ...")
print(f"    Final cell count: {adata_full.n_obs:,}")
print()

# Clean up
del adata_hvg
gc.collect()

# ===== 17. Compute scANVI UMAP if Needed =====
if has_scanvi and not has_scanvi_umap:
    print("🗺️  Computing UMAP from scANVI latent...")
    sc.pp.neighbors(adata_full, use_rep='X_scanvi', n_neighbors=15, random_state=42)
    sc.tl.umap(adata_full, random_state=42)
    adata_full.obsm['X_umap_scanvi'] = adata_full.obsm['X_umap'].copy()
    has_scanvi_umap = True
    print("  ✓ scANVI UMAP computed")
    print()

# ===== 18. Marker Gene Analysis =====
print("🔬 Computing marker genes on full gene set...")

min_cluster_size = final_cluster_counts.min()
if min_cluster_size < MIN_CELLS_FOR_MARKER:
    print(f"⚠️  Warning: Smallest cluster has only {min_cluster_size} cells")
    print(f"  Recommended minimum: {MIN_CELLS_FOR_MARKER} cells")
    print()

marker_success = False
potential_merges = []  # Initialize here to avoid NameError

try:
    sc.tl.rank_genes_groups(
        adata_full,
        groupby='leiden',
        method='wilcoxon',
        layer='log1p',
        use_raw=False,
        corr_method='benjamini-hochberg',
        pts=True,
        tie_correct=True,
        key_added='rank_genes_groups'
    )
    print("  ✓ Marker genes computed successfully")
    print("  ✓ Applied Benjamini-Hochberg FDR correction")
    print("  ✓ Used log1p layer on full gene set")
    print()
    marker_success = True
    
    # Display top markers
    print("📊 Top 5 marker genes per subcluster:")
    print("-" * 80)
    result = adata_full.uns['rank_genes_groups']
    for cluster in adata_full.obs['leiden'].cat.categories:
        genes = result['names'][cluster][:5]
        pvals = result['pvals_adj'][cluster][:5]
        logfcs = result['logfoldchanges'][cluster][:5]
        
        print(f"\nSubcluster {cluster} ({final_cluster_counts[cluster]} cells):")
        for gene, pval, lfc in zip(genes, pvals, logfcs):
            sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else ""
            print(f"  {gene}: log2FC={lfc:.2f}, FDR={pval:.2e} {sig}")
    print("-" * 80)
    print()
    
except Exception as e:
    print(f"❌ ERROR in marker gene analysis: {e}")
    print("  Continuing without marker gene analysis...")
    adata_full.uns['rank_genes_groups'] = None

# ⭐ QUALITY CHECK: Marker overlap (independent of display success)
if marker_success and adata_full.uns.get('rank_genes_groups') is not None:
    try:
        print("🔍 Quality Check: Marker Gene Overlap")
        print("-" * 80)
        
        n_top_check = 20
        overlap_threshold = 0.5
        result = adata_full.uns['rank_genes_groups']
        cluster_list = list(adata_full.obs['leiden'].cat.categories)
        
        for i, cluster_i in enumerate(cluster_list):
            markers_i = set(result['names'][cluster_i][:n_top_check])
            
            for cluster_j in cluster_list[i+1:]:
                markers_j = set(result['names'][cluster_j][:n_top_check])
                overlap = markers_i & markers_j
                overlap_pct = len(overlap) / n_top_check
                
                if overlap_pct >= overlap_threshold:
                    potential_merges.append({
                        'cluster_1': f'T_c{cluster_i}',
                        'cluster_2': f'T_c{cluster_j}',
                        'overlap_genes': len(overlap),
                        'overlap_pct': f'{overlap_pct*100:.1f}%'
                    })
        
        if potential_merges:
            print(f"  ⚠️  High marker overlap detected (>{overlap_threshold*100:.0f}%):")
            print(f"     These subclusters may be too similar and could be merged:\n")
            for merge in potential_merges:
                print(f"     {merge['cluster_1']} ↔ {merge['cluster_2']}: "
                      f"{merge['overlap_genes']}/{n_top_check} shared markers ({merge['overlap_pct']})")
            print(f"\n     Consider merging after manual validation")
        else:
            print(f"  ✓ All subcluster pairs have distinct markers (<{overlap_threshold*100:.0f}% overlap)")
        
        print("-" * 80)
        print()
        
    except Exception as e:
        print(f"  ⚠️  Overlap check failed: {e}")
        print()

# ===== 19. Visualizations =====
print("🎨 Generating visualizations...")

fig_dir = OUTPUT_DIR / 'figures'
fig_dir.mkdir(exist_ok=True)

# 1. Subclusters on scANVI UMAP (if available)
if has_scanvi_umap:
    print("  Using original scANVI UMAP for visualization")
    
    scanvi_basis = SCANVI_UMAP_KEY.replace('X_', '')
    
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.embedding(
        adata_full,
        basis=scanvi_basis,
        color='leiden',
        legend_loc='on data',
        legend_fontsize=10,
        title='T Cell Subclusters (Original scANVI UMAP)',
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(fig_dir / 'umap_subclusters_scanvi.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: umap_subclusters_scanvi.png")

# 2. Subclusters on BBKNN UMAP
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.embedding(
    adata_full,
    basis='umap_bbknn',
    color='leiden',
    legend_loc='on data',
    legend_fontsize=10,
    title='T Cell Subclusters (BBKNN-derived UMAP)',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(fig_dir / 'umap_subclusters_bbknn.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: umap_subclusters_bbknn.png")

# 3. Batch distribution
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.embedding(
    adata_full,
    basis='umap_bbknn',
    color=BATCH_KEY,
    title='Batch Distribution (Sample)',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(fig_dir / 'umap_batches.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: umap_batches.png")

# 4. Cell types (scANVI predictions)
if CELLTYPE_COLUMN in adata_full.obs.columns:
    fig, ax = plt.subplots(figsize=(14, 10))
    sc.pl.embedding(
        adata_full,
        basis='umap_bbknn',
        color=CELLTYPE_COLUMN,
        title='scANVI Cell Type Predictions',
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(fig_dir / 'umap_celltypes.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: umap_celltypes.png")

# 5. Cluster size barplot
fig, ax = plt.subplots(figsize=(max(10, n_clusters * 0.8), 6))
colors = plt.cm.tab20(np.linspace(0, 1, len(final_cluster_counts)))
ax.bar(final_cluster_counts.index.astype(str), final_cluster_counts.values, color=colors)
ax.set_xlabel('Subcluster', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title('T Cell Subcluster Size Distribution', fontsize=14)
ax.axhline(y=MIN_CELLS_PER_CLUSTER, color='red', linestyle='--', 
           label=f'Min threshold ({MIN_CELLS_PER_CLUSTER})')
ax.legend()
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(fig_dir / 'cluster_sizes.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: cluster_sizes.png")

# 6. Known marker expression
markers_in_data = [m for m in KNOWN_MARKERS if m in adata_full.var_names]
if markers_in_data:
    n_markers = len(markers_in_data)
    ncols = min(3, n_markers)
    nrows = (n_markers + ncols - 1) // ncols
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
    if n_markers == 1:
        axes = [axes]
    elif nrows == 1:
        axes = axes
    else:
        axes = axes.flatten()
    
    for idx_m, marker in enumerate(markers_in_data):
        sc.pl.embedding(
            adata_full,
            basis='umap_bbknn',
            color=marker,
            layer='log1p',
            title=f'{marker} Expression',
            ax=axes[idx_m],
            show=False
        )
    
    if nrows > 1:
        for idx_m in range(n_markers, len(axes)):
            axes[idx_m].axis('off')
    
    plt.tight_layout()
    plt.savefig(fig_dir / 'known_markers_umap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  ✓ Saved: known_markers_umap.png")

# 7. Marker heatmap and dotplot
if marker_success and adata_full.uns.get('rank_genes_groups') is not None:
    try:
        fig = sc.pl.rank_genes_groups_heatmap(
            adata_full,
            n_genes=10,
            groupby='leiden',
            layer='log1p',
            show=False,
            return_fig=True
        )
        fig.savefig(fig_dir / 'marker_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  ✓ Saved: marker_heatmap.png")
        
        fig = sc.pl.rank_genes_groups_dotplot(
            adata_full,
            n_genes=5,
            groupby='leiden',
            layer='log1p',
            show=False,
            return_fig=True
        )
        fig.savefig(fig_dir / 'marker_dotplot.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("  ✓ Saved: marker_dotplot.png")
    except Exception as e:
        print(f"  ⚠️  Could not generate marker visualizations: {e}")

print("  ✓ All visualizations generated")
print()

# ===== 20. Save Results =====
print("💾 Saving results...")

# Add metadata
adata_full.uns['subcluster_analysis_params'] = {
    'pipeline_version': 'v1.0_tcell_subcluster',
    'date': '2025-01-26',
    'input_file': str(INPUT_PATH),
    'batch_key': BATCH_KEY,
    'celltype_column': CELLTYPE_COLUMN,
    'hvg_method': hvg_method,
    'hvg_layer': 'log1p',
    'n_hvg': n_hvg,
    'resolution': RESOLUTION,
    'n_pcs': N_PCS,
    'min_cells_per_cluster': MIN_CELLS_PER_CLUSTER,
    'bbknn_neighbors_within_batch': neighbors_within_batch,
    'marker_layer': 'log1p',
    'marker_correction_method': 'benjamini-hochberg',
    'n_cells_final': adata_full.n_obs,
    'n_subclusters': n_clusters,
    'known_markers': KNOWN_MARKERS,
    'has_scanvi_umap': has_scanvi_umap,
    'preprocessing_validated': True,
    'naming_convention': 'T_c0, T_c1, ... (hierarchical L3 format)',
    'quality_checks': {
        'small_cluster_threshold': 50,
        'marker_overlap_threshold': 0.5,
        'small_clusters_detected': len(small_clusters_qc),
        'potential_merges_detected': len(potential_merges)
    }
}

# Save h5ad
output_filename = "adata_tcell_subcluster_analyzed.h5ad"
output_path = OUTPUT_DIR / output_filename
adata_full.write_h5ad(output_path, compression='gzip', compression_opts=9)
print(f"  ✓ Main output: {output_path}")

# Export cluster assignments
cluster_assignments = pd.DataFrame({
    'cell_barcode': adata_full.obs_names,
    'leiden_id': adata_full.obs['leiden'].values,
    'cell_type_L3': adata_full.obs['cell_type_L3'].values,
    'batch': adata_full.obs[BATCH_KEY].values
})
if CELLTYPE_COLUMN in adata_full.obs.columns:
    cluster_assignments['scanvi_celltype'] = adata_full.obs[CELLTYPE_COLUMN].values

assignment_path = OUTPUT_DIR / "cluster_assignments.csv"
cluster_assignments.to_csv(assignment_path, index=False)
print(f"  ✓ Cluster assignments: {assignment_path}")

# Export marker genes
if marker_success and adata_full.uns.get('rank_genes_groups') is not None:
    if USE_MODERN_MARKER_EXPORT:
        try:
            dfs = []
            for g in adata_full.obs['leiden'].cat.categories:
                df = sc.get.rank_genes_groups_df(
                    adata_full, 
                    group=g, 
                    key='rank_genes_groups'
                )
                df['subcluster'] = g
                dfs.append(df)
            all_markers = pd.concat(dfs, ignore_index=True)
        except:
            USE_MODERN_MARKER_EXPORT = False
    
    if not USE_MODERN_MARKER_EXPORT:
        marker_df_list = []
        result = adata_full.uns['rank_genes_groups']
        
        for cluster in adata_full.obs['leiden'].cat.categories:
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
    
    marker_path = OUTPUT_DIR / "marker_genes_FDR.csv"
    all_markers.to_csv(marker_path, index=False)
    print(f"  ✓ Marker genes: {marker_path}")

# Export summary
summary_stats = {
    'pipeline_version': 'v1.0_tcell_subcluster',
    'date': '2025-01-26',
    'cells_after_qc': adata_full.n_obs,
    'n_subclusters': n_clusters,
    'n_batches': adata_full.obs[BATCH_KEY].nunique(),
    'n_hvg': n_hvg,
    'hvg_method': hvg_method,
    'hvg_layer': 'log1p',
    'resolution': RESOLUTION,
    'marker_layer': 'log1p',
    'marker_correction': 'benjamini-hochberg',
    'marker_analysis_success': marker_success,
    'has_scanvi_umap': has_scanvi_umap
}

summary_df = pd.DataFrame([summary_stats])
summary_path = OUTPUT_DIR / "analysis_summary.csv"
summary_df.to_csv(summary_path, index=False)
print(f"  ✓ Summary statistics: {summary_path}")

print("  ✓ All results saved")
print()

# ===== 21. Final Summary =====
print("=" * 80)
print("✅ T CELL SUBCLUSTER ANALYSIS COMPLETED!")
print("=" * 80)
print()

print(f"📊 Summary:")
print(f"   Final cells: {adata_full.n_obs:,}")
print(f"   Subclusters: {n_clusters}")
print(f"   Naming format: T_c0, T_c1, ... (hierarchical L3)")
print(f"   Batches: {adata_full.obs[BATCH_KEY].nunique()}")
print(f"   Marker analysis: {'SUCCESS' if marker_success else 'FAILED'}")
print(f"   Has scANVI UMAP: {has_scanvi_umap}")
print()

# Quality check summary
if small_clusters_qc:
    print(f"⚠️  Quality Checks:")
    print(f"   Small clusters (<50 cells): {len(small_clusters_qc)}")
    for cluster_name, count in small_clusters_qc:
        print(f"     - {cluster_name}: {count} cells (consider merging)")

if potential_merges:
    if not small_clusters_qc:
        print(f"⚠️  Quality Checks:")
    print(f"   High marker overlap: {len(potential_merges)} cluster pairs")
    for merge in potential_merges[:3]:  # Show first 3
        print(f"     - {merge['cluster_1']} ↔ {merge['cluster_2']}: {merge['overlap_pct']} overlap")
    if len(potential_merges) > 3:
        print(f"     - ... and {len(potential_merges)-3} more pairs")
    print(f"   See console output for details")

if not small_clusters_qc and not potential_merges:
    print(f"✓ Quality Checks:")
    print(f"   All subclusters passed quality thresholds")

print()

print(f"📁 Output Directory:")
print(f"   {OUTPUT_DIR}")
print()

print(f"📄 Key Files:")
print(f"   Main data: {output_filename}")
print(f"   Markers: marker_genes_FDR.csv")
print(f"   Assignments: cluster_assignments.csv")
print(f"   Summary: analysis_summary.csv")
print()

print(f"📈 Figures:")
print(f"   figures/umap_subclusters_scanvi.png (original scANVI)")
print(f"   figures/umap_subclusters_bbknn.png (BBKNN-derived)")
print(f"   figures/umap_batches.png")
if CELLTYPE_COLUMN in adata_full.obs.columns:
    print(f"   figures/umap_celltypes.png")
print(f"   figures/cluster_sizes.png")
print(f"   figures/known_markers_umap.png")
print(f"   figures/marker_heatmap.png")
print(f"   figures/marker_dotplot.png")
print()

print("=" * 80)

gc.collect()
