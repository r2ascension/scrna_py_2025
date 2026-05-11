#!/usr/bin/env python3
"""
Universal Cell-Type Subcluster Analysis Pipeline v3.0
======================================================
Author: r2end
Date: 2025-01-21

Purpose:
--------
Generic pipeline for subcluster analysis of ANY cell type from annotated data.
Only need to modify CELL_TYPE_TO_ANALYZE to analyze different cell types.

Key Features:
-------------
1. Cell-type-specific subsetting from full annotated dataset
2. Robust error handling to prevent ZeroDivisionError
3. Benjamini-Hochberg FDR correction for marker genes
4. Comprehensive quality checks and validation
5. Fully configurable parameters for different cell types

Workflow:
---------
Full annotated data → Filter by cell type → HVG → BBKNN → Leiden → Markers → Save

Usage:
------
Simply modify these parameters in Configuration section:
- CELL_TYPE_TO_ANALYZE: "Myeloid", "Stromal", "B", "T", "Epithelial", etc.
- INPUT_PATH: Path to full annotated h5ad
- OUTPUT_DIR: Where to save results
- CELLTYPE_COLUMN: Column name containing cell type annotations

Input: 
------
- Full annotated dataset (e.g., allcells_scanvi_annotated.h5ad)
- Must contain cell type annotations

Output:
-------
- {celltype}_subcluster_analyzed.h5ad
- marker_genes_with_FDR.csv
- Visualization plots

Memory: ~20-40GB depending on cell type size
Runtime: ~30min - 2hr depending on cell count
"""

# ===== 1. Import Libraries =====
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("=" * 80)
print("UNIVERSAL CELL-TYPE SUBCLUSTER ANALYSIS PIPELINE v3.0")
print("=" * 80)
print()

# ===== 2. Configuration - MODIFY HERE =====
# ========================================
# Cell Type Selection (MAIN PARAMETER TO CHANGE)
# ========================================
CELL_TYPE_TO_ANALYZE = "Myeloid"  # ← Change this: "Myeloid", "Stromal", "B", "T", "Epithelial", etc.

# ========================================
# Paths
# ========================================
INPUT_PATH = "/path/to/allcells_scanvi_annotated.h5ad"  # Full annotated dataset
OUTPUT_DIR = Path(f"/path/to/output/{CELL_TYPE_TO_ANALYZE.lower()}_subcluster")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ========================================
# Cell Type Column
# ========================================
# Column name that contains cell type annotations
# Common options: 'cell_type', 'celltype_major', 'scanvi_predictions', 'celltypist_predictions'
CELLTYPE_COLUMN = 'celltype_major'

# ========================================
# Cell Type Filtering Strategy
# ========================================
# How to match cell types? Options:
# - "exact": Exact match (e.g., only "Myeloid")
# - "contains": Contains substring (e.g., matches "Myeloid", "Myeloid_Macro", etc.)
# - "list": Provide explicit list of cell types to include
FILTER_STRATEGY = "exact"  # or "contains" or "list"

# If FILTER_STRATEGY = "list", specify cell types to include:
CELL_TYPES_TO_INCLUDE = ["Myeloid", "Macrophage", "Monocyte"]  # Used only if FILTER_STRATEGY="list"

# ========================================
# Batch Effect Parameters
# ========================================
BATCH_KEY = 'Sample'  # Column for batch correction

# ========================================
# Quality Control Thresholds
# ========================================
MIN_CELLS_PER_BATCH = 3           # Remove batches smaller than this
MIN_CELLS_PER_CLUSTER = 10        # Remove clusters smaller than this (prevent ZeroDivisionError)
MIN_CELLS_FOR_MARKER = 20         # Minimum cells for reliable marker analysis

# ========================================
# Analysis Parameters
# ========================================
N_HVG = 3000                      # Number of highly variable genes
N_PCS = 50                        # Number of PCs for BBKNN
RESOLUTION = 1.0                  # Leiden clustering resolution
NEIGHBORS_WITHIN_BATCH = 5        # BBKNN parameter (will auto-adjust if needed)

# ========================================
# Cell-Type-Specific Parameters (Optional Fine-tuning)
# ========================================
# Adjust these based on cell type characteristics
CELLTYPE_PARAMS = {
    "Myeloid": {
        "n_hvg": 3000,
        "resolution": 1.0,
        "min_cells_per_cluster": 10,
        "known_markers": ["CD14", "FCGR3A", "CD68", "CD163", "LYZ"]
    },
    "Stromal": {
        "n_hvg": 3000,
        "resolution": 0.8,
        "min_cells_per_cluster": 10,
        "known_markers": ["COL1A1", "COL3A1", "DCN", "LUM", "PDGFRA"]
    },
    "B": {
        "n_hvg": 2500,
        "resolution": 0.8,
        "min_cells_per_cluster": 15,
        "known_markers": ["CD79A", "CD79B", "MS4A1", "CD19", "IGHM"]
    },
    "T": {
        "n_hvg": 3000,
        "resolution": 1.0,
        "min_cells_per_cluster": 10,
        "known_markers": ["CD3D", "CD3E", "CD4", "CD8A", "IL7R"]
    },
    "Epithelial": {
        "n_hvg": 3500,
        "resolution": 1.2,
        "min_cells_per_cluster": 15,
        "known_markers": ["EPCAM", "KRT5", "TP63", "MUC5AC", "FOXJ1"]
    }
}

# Apply cell-type-specific parameters if available
if CELL_TYPE_TO_ANALYZE in CELLTYPE_PARAMS:
    params = CELLTYPE_PARAMS[CELL_TYPE_TO_ANALYZE]
    N_HVG = params.get("n_hvg", N_HVG)
    RESOLUTION = params.get("resolution", RESOLUTION)
    MIN_CELLS_PER_CLUSTER = params.get("min_cells_per_cluster", MIN_CELLS_PER_CLUSTER)
    KNOWN_MARKERS = params.get("known_markers", [])
    print(f"✓ Loaded cell-type-specific parameters for {CELL_TYPE_TO_ANALYZE}")
else:
    KNOWN_MARKERS = []
    print(f"⚠️  Using default parameters for {CELL_TYPE_TO_ANALYZE}")

# ========================================
# Scanpy Settings
# ========================================
sc.settings.verbosity = 1
sc.settings.n_jobs = 48
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))

print()
print("Configuration Summary:")
print("-" * 80)
print(f"Cell type to analyze: {CELL_TYPE_TO_ANALYZE}")
print(f"Cell type column: {CELLTYPE_COLUMN}")
print(f"Filter strategy: {FILTER_STRATEGY}")
print(f"Input: {INPUT_PATH}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Batch key: {BATCH_KEY}")
print(f"HVG count: {N_HVG}")
print(f"Clustering resolution: {RESOLUTION}")
print(f"Min cells per cluster: {MIN_CELLS_PER_CLUSTER}")
if KNOWN_MARKERS:
    print(f"Known markers: {', '.join(KNOWN_MARKERS)}")
print("-" * 80)
print()

# ===== 3. Load Full Annotated Data =====
print("📂 Loading full annotated dataset...")
adata_full = sc.read_h5ad(INPUT_PATH)
print(f"  ✓ Loaded: {adata_full.n_obs} cells × {adata_full.n_vars} genes")
print(f"  Available layers: {list(adata_full.layers.keys())}")
print(f"  Available obs columns: {list(adata_full.obs.columns[:10])}... (showing first 10)")
print()

# ===== 4. Validate Cell Type Column =====
print("🔍 Validating cell type annotations...")
if CELLTYPE_COLUMN not in adata_full.obs.columns:
    raise ValueError(
        f"Cell type column '{CELLTYPE_COLUMN}' not found!\n"
        f"Available columns: {list(adata_full.obs.columns)}"
    )

celltype_counts = adata_full.obs[CELLTYPE_COLUMN].value_counts()
print(f"  Cell type distribution in full dataset:")
for ct, count in celltype_counts.head(10).items():
    print(f"    {ct}: {count:,} cells")
if len(celltype_counts) > 10:
    print(f"    ... and {len(celltype_counts) - 10} more cell types")
print()

# ===== 5. Filter by Cell Type =====
print(f"🔬 Filtering for {CELL_TYPE_TO_ANALYZE} cells...")

if FILTER_STRATEGY == "exact":
    # Exact match
    cell_mask = adata_full.obs[CELLTYPE_COLUMN] == CELL_TYPE_TO_ANALYZE
    
elif FILTER_STRATEGY == "contains":
    # Contains substring (case-insensitive)
    cell_mask = adata_full.obs[CELLTYPE_COLUMN].str.contains(
        CELL_TYPE_TO_ANALYZE, 
        case=False, 
        na=False
    )
    
elif FILTER_STRATEGY == "list":
    # Match any in list
    cell_mask = adata_full.obs[CELLTYPE_COLUMN].isin(CELL_TYPES_TO_INCLUDE)
    print(f"  Matching cell types: {CELL_TYPES_TO_INCLUDE}")
    
else:
    raise ValueError(f"Unknown FILTER_STRATEGY: {FILTER_STRATEGY}")

adata = adata_full[cell_mask].copy()
n_cells_selected = adata.n_obs

if n_cells_selected == 0:
    raise ValueError(
        f"No cells found for {CELL_TYPE_TO_ANALYZE}!\n"
        f"Available cell types: {list(celltype_counts.index)}"
    )

print(f"  ✓ Selected {n_cells_selected:,} cells ({n_cells_selected/adata_full.n_obs*100:.1f}% of total)")

# Show which specific cell types were included
if FILTER_STRATEGY == "contains" or FILTER_STRATEGY == "list":
    included_types = adata.obs[CELLTYPE_COLUMN].value_counts()
    print(f"  Included cell type breakdown:")
    for ct, count in included_types.items():
        print(f"    {ct}: {count:,} cells")
print()

# Clean up full dataset to free memory
del adata_full
gc.collect()

# ===== 6. Data Quality Checks =====
print("🔍 Performing data quality checks...")

# Check for required data structure
if 'counts' not in adata.layers:
    raise ValueError("Missing 'counts' layer in adata!")

# Check batch distribution
batch_counts = adata.obs[BATCH_KEY].value_counts()
print(f"  Batch distribution ({len(batch_counts)} batches):")
for batch, count in batch_counts.head(10).items():
    print(f"    {batch}: {count} cells")
if len(batch_counts) > 10:
    print(f"    ... and {len(batch_counts) - 10} more batches")
print()

# Remove small batches
small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index
if len(small_batches) > 0:
    print(f"⚠️  Removing {len(small_batches)} small batches (< {MIN_CELLS_PER_BATCH} cells):")
    for batch in small_batches[:5]:  # Show first 5
        print(f"    {batch}: {batch_counts[batch]} cells")
    if len(small_batches) > 5:
        print(f"    ... and {len(small_batches) - 5} more")
    
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  ✓ Remaining: {adata.n_obs:,} cells in {adata.obs[BATCH_KEY].nunique()} batches")
    print()

# ===== 7. Preprocessing =====
print("🔬 Preprocessing...")

# Ensure proper data structure
if 'log1p' in adata.layers:
    print("  Using existing log1p layer")
    adata.X = adata.layers['log1p'].copy()
else:
    print("  Computing normalization and log-transformation")
    # Use counts layer for normalization
    adata.X = adata.layers['counts'].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

print("  ✓ Preprocessing completed")
print()

# ===== 8. HVG Selection =====
print(f"🎯 Selecting top {N_HVG} highly variable genes...")

try:
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        batch_key=BATCH_KEY,
        subset=False
    )
    hvg_method = "batch-aware"
    print(f"  ✓ Batch-aware HVG selection successful")
except Exception as e:
    print(f"  ⚠️  Batch-aware HVG failed: {e}")
    print("  Falling back to non-batch-aware method")
    sc.pp.highly_variable_genes(
        adata,
        layer='counts',
        n_top_genes=N_HVG,
        subset=False
    )
    hvg_method = "non-batch-aware"

n_hvg = adata.var['highly_variable'].sum()
print(f"  ✓ Selected {n_hvg} HVGs ({hvg_method})")

# Force include known markers if specified
if KNOWN_MARKERS:
    known_in_data = [m for m in KNOWN_MARKERS if m in adata.var_names]
    if known_in_data:
        adata.var.loc[known_in_data, 'highly_variable'] = True
        n_hvg = adata.var['highly_variable'].sum()
        print(f"  ✓ Force-included {len(known_in_data)} known markers")
        print(f"  Final HVG count: {n_hvg}")
print()

# ===== 9. Save Full Gene Set =====
print("💾 Preserving full gene set for downstream analysis...")
adata.raw = adata.copy()
print("  ✓ Full gene set saved to .raw (shared memory)")
print()

# ===== 10. Subset to HVGs =====
print("✂️  Subsetting to HVGs for integration...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"  ✓ Working with {adata.n_vars} HVGs")
print()

# ===== 11. Dimensionality Reduction =====
print(f"📊 Computing PCA ({N_PCS} components)...")
sc.tl.pca(adata, n_comps=N_PCS, svd_solver='arpack')
print("  ✓ PCA computed")
print()

# ===== 12. Batch Effect Correction with BBKNN =====
print("🔗 Running BBKNN integration...")

# Auto-adjust neighbors_within_batch based on smallest batch
batch_sizes = adata.obs[BATCH_KEY].value_counts()
min_batch_size = batch_sizes.min()
neighbors_within_batch = min(NEIGHBORS_WITHIN_BATCH, max(3, min_batch_size - 1))

print(f"  Smallest batch size: {min_batch_size}")
print(f"  Using neighbors_within_batch: {neighbors_within_batch}")

sc.external.pp.bbknn(
    adata,
    batch_key=BATCH_KEY,
    neighbors_within_batch=neighbors_within_batch,
    n_pcs=N_PCS,
    trim=None
)
print("  ✓ BBKNN integration completed")
print()

# ===== 13. Leiden Clustering =====
print(f"🎨 Computing Leiden clustering (resolution={RESOLUTION})...")
sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden')
n_clusters = adata.obs['leiden'].nunique()
print(f"  ✓ Identified {n_clusters} subclusters")
print()

# Check cluster sizes
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
print("  Subcluster sizes:")
for cluster, count in cluster_counts.items():
    status = "⚠️ SMALL" if count < MIN_CELLS_PER_CLUSTER else "✓"
    print(f"    Subcluster {cluster}: {count:,} cells {status}")
print()

# ===== 14. Filter Small Clusters =====
small_clusters = cluster_counts[cluster_counts < MIN_CELLS_PER_CLUSTER].index.tolist()
if len(small_clusters) > 0:
    print(f"⚠️  Removing {len(small_clusters)} small subclusters (< {MIN_CELLS_PER_CLUSTER} cells):")
    print(f"  Subclusters: {small_clusters}")
    
    adata = adata[~adata.obs['leiden'].isin(small_clusters)].copy()
    
    # Recompute neighbors and clustering after filtering
    print("  Recomputing BBKNN and clustering after filtering...")
    sc.external.pp.bbknn(
        adata,
        batch_key=BATCH_KEY,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=N_PCS,
        trim=None
    )
    sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden')
    n_clusters = adata.obs['leiden'].nunique()
    
    print(f"  ✓ Final: {adata.n_obs:,} cells in {n_clusters} subclusters")
    
    # Show final cluster distribution
    final_cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    print("  Final subcluster sizes:")
    for cluster, count in final_cluster_counts.items():
        print(f"    Subcluster {cluster}: {count:,} cells")
    print()

# ===== 15. UMAP =====
print("🗺️  Computing UMAP...")
sc.tl.umap(adata)
print("  ✓ UMAP computed")
print()

# ===== 16. Marker Gene Analysis with BH Correction =====
print("🔬 Computing marker genes with Benjamini-Hochberg correction...")

# Final validation
final_cluster_counts = adata.obs['leiden'].value_counts()
min_cluster_size = final_cluster_counts.min()

if min_cluster_size < MIN_CELLS_FOR_MARKER:
    print(f"⚠️  Warning: Smallest cluster has only {min_cluster_size} cells")
    print(f"  Recommended minimum: {MIN_CELLS_FOR_MARKER} cells for robust statistics")
    print()

try:
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        use_raw=True,  # Use full gene set from .raw
        corr_method='benjamini-hochberg',  # FDR correction
        pts=True,
        tie_correct=True,
        key_added='rank_genes_groups'
    )
    print("  ✓ Marker genes computed successfully")
    print("  ✓ Applied Benjamini-Hochberg FDR correction")
    print()
    
    # Extract and display top markers
    print("📊 Top 5 marker genes per subcluster:")
    print("-" * 80)
    result = adata.uns['rank_genes_groups']
    for cluster in adata.obs['leiden'].cat.categories:
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
    print("  Possible causes:")
    print("  1. Clusters still too small despite filtering")
    print("  2. Data quality issues")
    print("  3. Insufficient variance in gene expression")
    print("  4. Memory constraints")
    print()
    print("  Continuing without marker gene analysis...")
    adata.uns['rank_genes_groups'] = None

# ===== 17. Visualizations =====
print("🎨 Generating visualizations...")

# Create figure directory
fig_dir = OUTPUT_DIR / 'figures'
fig_dir.mkdir(exist_ok=True)

# 1. UMAP by subcluster
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(
    adata,
    color='leiden',
    legend_loc='on data',
    legend_fontsize=10,
    title=f'{CELL_TYPE_TO_ANALYZE} Subclusters (BBKNN + Leiden)',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(fig_dir / 'umap_subclusters.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. UMAP by batch
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.umap(
    adata,
    color=BATCH_KEY,
    title=f'Batch Distribution - {CELL_TYPE_TO_ANALYZE}',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(fig_dir / 'umap_batches.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. UMAP by original cell type annotation (if multiple were included)
if FILTER_STRATEGY in ["contains", "list"] and adata.obs[CELLTYPE_COLUMN].nunique() > 1:
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(
        adata,
        color=CELLTYPE_COLUMN,
        title=f'Original Cell Type Annotations',
        ax=ax,
        show=False
    )
    plt.tight_layout()
    plt.savefig(fig_dir / 'umap_original_celltypes.png', dpi=300, bbox_inches='tight')
    plt.close()

# 4. Cluster size barplot
fig, ax = plt.subplots(figsize=(max(10, n_clusters * 0.8), 6))
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
colors = plt.cm.tab20(np.linspace(0, 1, len(cluster_counts)))
ax.bar(cluster_counts.index.astype(str), cluster_counts.values, color=colors)
ax.set_xlabel('Subcluster', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title(f'{CELL_TYPE_TO_ANALYZE} Subcluster Size Distribution', fontsize=14)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(fig_dir / 'cluster_sizes.png', dpi=300, bbox_inches='tight')
plt.close()

# 5. Known marker expression (if specified)
if KNOWN_MARKERS:
    markers_in_data = [m for m in KNOWN_MARKERS if m in adata.raw.var_names]
    if markers_in_data:
        n_markers = len(markers_in_data)
        ncols = min(3, n_markers)
        nrows = (n_markers + ncols - 1) // ncols
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows))
        axes = axes.flatten() if n_markers > 1 else [axes]
        
        for idx, marker in enumerate(markers_in_data):
            sc.pl.umap(
                adata,
                color=marker,
                use_raw=True,
                title=f'{marker} Expression',
                ax=axes[idx],
                show=False
            )
        
        # Hide extra subplots
        for idx in range(n_markers, len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        plt.savefig(fig_dir / 'known_markers_umap.png', dpi=300, bbox_inches='tight')
        plt.close()

# 6. Marker heatmap (if markers were computed)
if adata.uns.get('rank_genes_groups') is not None:
    try:
        fig = sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=10,
            groupby='leiden',
            use_raw=True,
            show=False,
            return_fig=True
        )
        fig.savefig(fig_dir / 'marker_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Dotplot
        fig = sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=5,
            groupby='leiden',
            use_raw=True,
            show=False,
            return_fig=True
        )
        fig.savefig(fig_dir / 'marker_dotplot.png', dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"  ⚠️  Could not generate marker heatmap: {e}")

print(f"  ✓ Visualizations saved to {fig_dir}/")
print()

# ===== 18. Save Results =====
print("💾 Saving results...")

# Add comprehensive metadata
adata.uns['subcluster_analysis_params'] = {
    'cell_type_analyzed': CELL_TYPE_TO_ANALYZE,
    'celltype_column': CELLTYPE_COLUMN,
    'filter_strategy': FILTER_STRATEGY,
    'hvg_method': hvg_method,
    'n_hvg': n_hvg,
    'resolution': RESOLUTION,
    'n_pcs': N_PCS,
    'min_cells_per_cluster': MIN_CELLS_PER_CLUSTER,
    'batch_key': BATCH_KEY,
    'bbknn_neighbors_within_batch': neighbors_within_batch,
    'marker_correction_method': 'benjamini-hochberg',
    'n_cells_final': adata.n_obs,
    'n_subclusters': n_clusters,
    'known_markers': KNOWN_MARKERS
}

# Save annotated data
output_filename = f"{CELL_TYPE_TO_ANALYZE.lower().replace(' ', '_')}_subcluster_analyzed.h5ad"
output_path = OUTPUT_DIR / output_filename
adata.write_h5ad(output_path, compression='gzip', compression_opts=9)
print(f"  ✓ Main output: {output_path}")

# Export cluster assignments
cluster_assignments = pd.DataFrame({
    'cell_barcode': adata.obs_names,
    'original_celltype': adata.obs[CELLTYPE_COLUMN].values,
    'subcluster': adata.obs['leiden'].values,
    'batch': adata.obs[BATCH_KEY].values
})
assignment_path = OUTPUT_DIR / f"{CELL_TYPE_TO_ANALYZE.lower()}_subcluster_assignments.csv"
cluster_assignments.to_csv(assignment_path, index=False)
print(f"  ✓ Cluster assignments: {assignment_path}")

# Export marker genes (if available)
if adata.uns.get('rank_genes_groups') is not None:
    marker_df_list = []
    result = adata.uns['rank_genes_groups']
    
    for cluster in adata.obs['leiden'].cat.categories:
        n_genes_cluster = len(result['names'][cluster])
        cluster_markers = pd.DataFrame({
            'subcluster': cluster,
            'gene': result['names'][cluster],
            'log2fc': result['logfoldchanges'][cluster],
            'pval': result['pvals'][cluster],
            'pval_adj': result['pvals_adj'][cluster],  # BH-corrected FDR
            'pct_in_group': result['pts'][cluster],
            'pct_out_group': result['pts_rest'][cluster]
        })
        marker_df_list.append(cluster_markers)
    
    all_markers = pd.concat(marker_df_list, ignore_index=True)
    marker_path = OUTPUT_DIR / f"{CELL_TYPE_TO_ANALYZE.lower()}_marker_genes_FDR.csv"
    all_markers.to_csv(marker_path, index=False)
    print(f"  ✓ Marker genes: {marker_path}")

# Export summary statistics
summary_stats = {
    'cell_type': CELL_TYPE_TO_ANALYZE,
    'total_cells': n_cells_selected,
    'cells_after_qc': adata.n_obs,
    'n_subclusters': n_clusters,
    'n_batches': adata.obs[BATCH_KEY].nunique(),
    'n_hvg': n_hvg,
    'hvg_method': hvg_method,
    'resolution': RESOLUTION,
    'marker_correction': 'benjamini-hochberg'
}

summary_df = pd.DataFrame([summary_stats])
summary_path = OUTPUT_DIR / f"{CELL_TYPE_TO_ANALYZE.lower()}_analysis_summary.csv"
summary_df.to_csv(summary_path, index=False)
print(f"  ✓ Summary statistics: {summary_path}")

print()

# ===== 19. Final Summary =====
print("=" * 80)
print("SUBCLUSTER ANALYSIS SUMMARY")
print("=" * 80)
print(f"Cell type analyzed: {CELL_TYPE_TO_ANALYZE}")
print(f"Cells in final dataset: {adata.n_obs:,}")
print(f"Number of subclusters: {n_clusters}")
print(f"HVG selection: {hvg_method} ({n_hvg} genes)")
print(f"Statistical correction: Benjamini-Hochberg FDR")
print(f"Integration method: BBKNN")
print()
print("Subcluster distribution:")
for cluster, count in cluster_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  Subcluster {cluster}: {count:,} cells ({pct:.1f}%)")
print()
print("Output files:")
print(f"  Main data: {output_path}")
print(f"  Assignments: {assignment_path}")
if adata.uns.get('rank_genes_groups') is not None:
    print(f"  Markers: {marker_path}")
print(f"  Summary: {summary_path}")
print(f"  Figures: {fig_dir}/")
print()
print("✅ Subcluster analysis completed successfully!")
print("=" * 80)

# Cleanup
gc.collect()
