#!/usr/bin/env python3
"""
Stromal Subcluster Analysis Pipeline v2.1 (FIXED)
=================================================
Author: r2end
Date: 2025-01-21

Key Fixes:
- Added minimum cluster size filtering to prevent ZeroDivisionError
- Implemented Benjamini-Hochberg FDR correction
- Enhanced data quality checks before marker gene analysis
- Added robust error handling

Purpose:
- Subcluster analysis of stromal/vascular cells
- BBKNN integration across batches
- Automated marker gene identification with proper statistical correction

Input: stromal_annotated.h5ad
Output: stromal_subcluster_analyzed.h5ad

Memory: ~20GB peak
Runtime: ~1 hour
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

# ===== 2. Configuration =====
# Paths
INPUT_PATH = "/path/to/stromal_annotated.h5ad"
OUTPUT_DIR = Path("/path/to/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Analysis parameters
BATCH_KEY = 'Sample'
MIN_CELLS_PER_CLUSTER = 10  # Critical: prevent division by zero
MIN_CELLS_FOR_MARKER_ANALYSIS = 20  # For robust statistics
N_HVG = 3000
RESOLUTION = 1.0

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.n_jobs = 48
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))

print("=" * 80)
print("STROMAL SUBCLUSTER ANALYSIS PIPELINE v2.1 (FIXED)")
print("=" * 80)
print(f"Input: {INPUT_PATH}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Batch key: {BATCH_KEY}")
print(f"Min cells per cluster: {MIN_CELLS_PER_CLUSTER}")
print(f"Resolution: {RESOLUTION}")
print()

# ===== 3. Load Data =====
print("📂 Loading data...")
adata = sc.read_h5ad(INPUT_PATH)
print(f"  ✓ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Available layers: {list(adata.layers.keys())}")
print()

# ===== 4. Data Quality Checks =====
print("🔍 Performing data quality checks...")

# Check for required data structure
if 'counts' not in adata.layers:
    raise ValueError("Missing 'counts' layer in adata!")

# Check batch distribution
batch_counts = adata.obs[BATCH_KEY].value_counts()
print(f"  Batch distribution:\n{batch_counts}")
print()

# Remove small batches
small_batches = batch_counts[batch_counts < 3].index
if len(small_batches) > 0:
    print(f"⚠️  Removing {len(small_batches)} small batches (< 3 cells):")
    print(f"  {list(small_batches)}")
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  ✓ Remaining: {adata.n_obs} cells")
    print()

# ===== 5. Preprocessing =====
print("🔬 Preprocessing...")

# Ensure proper data structure
if 'log1p' in adata.layers:
    print("  Using existing log1p layer")
    adata.X = adata.layers['log1p'].copy()
else:
    print("  Computing normalization and log-transformation")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

print("  ✓ Preprocessing completed")
print()

# ===== 6. HVG Selection =====
print("🎯 Selecting highly variable genes...")

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
print()

# ===== 7. Save Full Gene Set =====
print("💾 Preserving full gene set for downstream analysis...")
adata.raw = adata.copy()
print("  ✓ Full gene set saved to .raw")
print()

# ===== 8. Subset to HVGs =====
print("✂️  Subsetting to HVGs for integration...")
adata = adata[:, adata.var['highly_variable']].copy()
print(f"  ✓ Working with {adata.n_vars} HVGs")
print()

# ===== 9. Dimensionality Reduction =====
print("📊 Computing PCA...")
sc.tl.pca(adata, n_comps=50, svd_solver='arpack')
print("  ✓ PCA computed")
print()

# ===== 10. Batch Effect Correction with BBKNN =====
print("🔗 Running BBKNN integration...")

# Check batch sizes for BBKNN
batch_sizes = adata.obs[BATCH_KEY].value_counts()
min_batch_size = batch_sizes.min()
neighbors_within_batch = min(5, max(3, min_batch_size - 1))

print(f"  Min batch size: {min_batch_size}")
print(f"  Using neighbors_within_batch: {neighbors_within_batch}")

sc.external.pp.bbknn(
    adata,
    batch_key=BATCH_KEY,
    neighbors_within_batch=neighbors_within_batch,
    n_pcs=50,
    trim=None
)
print("  ✓ BBKNN integration completed")
print()

# ===== 11. Clustering =====
print("🎨 Computing Leiden clustering...")
sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden')
n_clusters = adata.obs['leiden'].nunique()
print(f"  ✓ Identified {n_clusters} clusters at resolution {RESOLUTION}")
print()

# Check cluster sizes
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
print("  Cluster sizes:")
for cluster, count in cluster_counts.items():
    status = "⚠️ SMALL" if count < MIN_CELLS_PER_CLUSTER else "✓"
    print(f"    Cluster {cluster}: {count} cells {status}")
print()

# ===== 12. Filter Small Clusters =====
small_clusters = cluster_counts[cluster_counts < MIN_CELLS_PER_CLUSTER].index.tolist()
if len(small_clusters) > 0:
    print(f"⚠️  Removing {len(small_clusters)} small clusters (< {MIN_CELLS_PER_CLUSTER} cells):")
    print(f"  Clusters: {small_clusters}")
    adata = adata[~adata.obs['leiden'].isin(small_clusters)].copy()
    
    # Recompute neighbors and clustering after filtering
    print("  Recomputing BBKNN after filtering...")
    sc.external.pp.bbknn(
        adata,
        batch_key=BATCH_KEY,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=50,
        trim=None
    )
    sc.tl.leiden(adata, resolution=RESOLUTION, key_added='leiden')
    n_clusters = adata.obs['leiden'].nunique()
    print(f"  ✓ Final: {adata.n_obs} cells in {n_clusters} clusters")
    print()

# ===== 13. UMAP =====
print("🗺️  Computing UMAP...")
sc.tl.umap(adata)
print("  ✓ UMAP computed")
print()

# ===== 14. Marker Gene Analysis with BH Correction =====
print("🔬 Computing marker genes with Benjamini-Hochberg correction...")

# Final cluster size check
final_cluster_counts = adata.obs['leiden'].value_counts()
min_cluster_size = final_cluster_counts.min()

if min_cluster_size < MIN_CELLS_FOR_MARKER_ANALYSIS:
    print(f"⚠️  Warning: Smallest cluster has only {min_cluster_size} cells")
    print(f"  This may affect statistical power for that cluster")
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
    
    # Extract top markers
    print("📊 Top markers per cluster:")
    result = adata.uns['rank_genes_groups']
    for cluster in adata.obs['leiden'].cat.categories:
        genes = result['names'][cluster][:5]
        pvals = result['pvals_adj'][cluster][:5]
        logfcs = result['logfoldchanges'][cluster][:5]
        
        print(f"\n  Cluster {cluster}:")
        for gene, pval, lfc in zip(genes, pvals, logfcs):
            print(f"    {gene}: log2FC={lfc:.2f}, FDR={pval:.2e}")
    print()
    
except Exception as e:
    print(f"❌ ERROR in marker gene analysis: {e}")
    print("  Possible causes:")
    print("  1. Clusters still too small")
    print("  2. Data quality issues")
    print("  3. Memory constraints")
    print()
    # Continue without marker genes
    adata.uns['rank_genes_groups'] = None

# ===== 15. Visualizations =====
print("🎨 Generating visualizations...")

# UMAP by cluster
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color='leiden',
    legend_loc='on data',
    legend_fontsize=8,
    title='Stromal Subclusters (BBKNN + Leiden)',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'umap_clusters.png', dpi=300, bbox_inches='tight')
plt.close()

# UMAP by batch
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color=BATCH_KEY,
    title='Batch Distribution',
    ax=ax,
    show=False
)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'umap_batches.png', dpi=300, bbox_inches='tight')
plt.close()

# Cluster size barplot
fig, ax = plt.subplots(figsize=(10, 6))
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
ax.bar(cluster_counts.index.astype(str), cluster_counts.values)
ax.set_xlabel('Cluster', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_title('Cluster Size Distribution', fontsize=14)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'cluster_sizes.png', dpi=300, bbox_inches='tight')
plt.close()

# Marker heatmap (if markers were computed)
if adata.uns.get('rank_genes_groups') is not None:
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=10,
        groupby='leiden',
        use_raw=True,
        show=False,
        save='_markers.png'
    )

print("  ✓ Visualizations saved")
print()

# ===== 16. Save Results =====
print("💾 Saving results...")

# Add metadata
adata.uns['analysis_params'] = {
    'hvg_method': hvg_method,
    'n_hvg': n_hvg,
    'resolution': RESOLUTION,
    'min_cells_per_cluster': MIN_CELLS_PER_CLUSTER,
    'batch_key': BATCH_KEY,
    'bbknn_neighbors_within_batch': neighbors_within_batch,
    'marker_correction_method': 'benjamini-hochberg'
}

# Save annotated data
output_path = OUTPUT_DIR / 'stromal_subcluster_analyzed.h5ad'
adata.write_h5ad(output_path, compression='gzip', compression_opts=9)
print(f"  ✓ Saved to: {output_path}")

# Export cluster assignments
cluster_assignments = pd.DataFrame({
    'cell_barcode': adata.obs_names,
    'leiden_cluster': adata.obs['leiden'].values,
    'batch': adata.obs[BATCH_KEY].values
})
cluster_assignments.to_csv(
    OUTPUT_DIR / 'cluster_assignments.csv',
    index=False
)
print(f"  ✓ Cluster assignments saved")

# Export marker genes (if available)
if adata.uns.get('rank_genes_groups') is not None:
    marker_df_list = []
    result = adata.uns['rank_genes_groups']
    
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_markers = pd.DataFrame({
            'cluster': cluster,
            'gene': result['names'][cluster],
            'log2fc': result['logfoldchanges'][cluster],
            'pval': result['pvals'][cluster],
            'pval_adj': result['pvals_adj'][cluster],  # BH-corrected
            'pct_in_group': result['pts'][cluster],
            'pct_out_group': result['pts_rest'][cluster]
        })
        marker_df_list.append(cluster_markers)
    
    all_markers = pd.concat(marker_df_list, ignore_index=True)
    all_markers.to_csv(
        OUTPUT_DIR / 'marker_genes_with_FDR.csv',
        index=False
    )
    print(f"  ✓ Marker genes with FDR correction saved")

print()

# ===== 17. Summary =====
print("=" * 80)
print("ANALYSIS SUMMARY")
print("=" * 80)
print(f"Final dataset: {adata.n_obs} cells × {adata.n_vars} HVGs")
print(f"Number of clusters: {n_clusters}")
print(f"HVG selection method: {hvg_method}")
print(f"Statistical correction: Benjamini-Hochberg FDR")
print(f"Integration method: BBKNN")
print()
print("Output files:")
print(f"  - {output_path}")
print(f"  - {OUTPUT_DIR / 'cluster_assignments.csv'}")
if adata.uns.get('rank_genes_groups') is not None:
    print(f"  - {OUTPUT_DIR / 'marker_genes_with_FDR.csv'}")
print(f"  - {OUTPUT_DIR / 'umap_clusters.png'}")
print(f"  - {OUTPUT_DIR / 'umap_batches.png'}")
print(f"  - {OUTPUT_DIR / 'cluster_sizes.png'}")
print()
print("✅ Analysis completed successfully!")
print("=" * 80)

# Cleanup
gc.collect()
