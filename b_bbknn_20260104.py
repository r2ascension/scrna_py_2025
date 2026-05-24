# ===== BBKNN Integration with Data Preprocessing (PRODUCTION VERSION) =====
# Purpose: Batch correction using BBKNN with automatic data normalization
# Author: r2end
# Date: 2025-01-05

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import sparse

# ===== Configuration =====
INPUT_PATH = "/home/h2048/data/R/1228/B/b_filtered_20260104.h5ad"
OUTPUT_DIR = "/home/h2048/data/R/1228/B"
BATCH_KEY = 'study'
N_PCS = 50
N_HVG = 3000  # Number of highly variable genes for PCA

# UMAP parameters for better separation
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.5

# Clustering resolution
LEIDEN_RESOLUTION = 1.0

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Output paths
OUTPUT_H5AD = os.path.join(OUTPUT_DIR, "ciliated_bbknn_integrated.h5ad")
OUTPUT_PDF1 = os.path.join(OUTPUT_DIR, "bbknn_umap_integration.pdf")
OUTPUT_PDF2 = os.path.join(OUTPUT_DIR, "bbknn_overview.pdf")

# Set random seed
np.random.seed(42)
sc.settings.n_jobs = 48
sc.settings.figdir = OUTPUT_DIR

# ===== 1. Load Data =====
adata = sc.read_h5ad(INPUT_PATH)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Batches: {adata.obs[BATCH_KEY].value_counts()}\n")

# ===== 2. Data Structure Check and Preprocessing =====
print("===== Data Preprocessing =====")
print(f"Available layers: {list(adata.layers.keys())}")
print(f"adata.X dtype: {adata.X.dtype}")
print(f"adata.X range: [{adata.X.min():.2f}, {adata.X.max():.2f}]")
print(f"adata.X is sparse: {sparse.issparse(adata.X)}")

# Detect if .X is counts or log-normalized
is_counts = (adata.X.max() > 20) or (adata.X.dtype in [np.int32, np.int64])

if is_counts:
    print("\n⚠️ Detected: adata.X contains RAW COUNTS")
    
    # Step 1: Calculate QC metrics FIRST (before any modification)
    if 'n_genes' not in adata.obs.columns:
        print("   Calculating QC metrics from raw counts...")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        # Rename to standard names
        if 'n_genes_by_counts' in adata.obs.columns:
            adata.obs['n_genes'] = adata.obs['n_genes_by_counts']
        if 'total_counts' in adata.obs.columns:
            adata.obs['n_counts'] = adata.obs['total_counts']
    
    # Step 2: Save raw counts to adata.raw (BEFORE normalization)
    print("   Saving raw counts to adata.raw.X...")
    adata.raw = adata  # This freezes the current state (counts + all genes)
    print(f"   ✅ adata.raw created with shape: {adata.raw.shape}")
    
    # Step 3: Normalize and log-transform
    print("   Normalizing to 10,000 counts per cell...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    print("   Log1p transformation...")
    sc.pp.log1p(adata)
    
    # Step 4: Ensure sparsity
    if not sparse.issparse(adata.X):
        print("   Converting to sparse matrix...")
        adata.X = sparse.csr_matrix(adata.X)
    
    # Step 5: Save log1p data to layers
    adata.layers['log1p'] = adata.X.copy()
    
    print(f"   ✅ Normalized. New .X range: [{adata.X.min():.2f}, {adata.X.max():.2f}]")
    print(f"   ✅ adata.X is sparse: {sparse.issparse(adata.X)}")
    
else:
    print("\n✅ Detected: adata.X already log-normalized")
    
    # Check if adata.raw exists
    if adata.raw is None:
        print("   ⚠️ Warning: No adata.raw found (raw counts not available)")
    else:
        print(f"   ✅ adata.raw exists with shape: {adata.raw.shape}")
    
    # Ensure log1p layer exists
    if 'log1p' not in adata.layers:
        print("   Saving current .X to layer['log1p']")
        adata.layers['log1p'] = adata.X.copy()
    
    # Calculate QC if missing
    if 'n_genes' not in adata.obs.columns:
        print("   Calculating QC metrics...")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        if 'n_genes_by_counts' in adata.obs.columns:
            adata.obs['n_genes'] = adata.obs['n_genes_by_counts']
        if 'total_counts' in adata.obs.columns:
            adata.obs['n_counts'] = adata.obs['total_counts']

# Final check: ensure .X is log-normalized
if 'log1p' in adata.layers:
    adata.X = adata.layers['log1p'].copy()

print(f"\n✅ Data ready for analysis")
print(f"   .X: log1p normalized, sparse={sparse.issparse(adata.X)}")
print(f"   .raw.X: raw counts" if adata.raw is not None else "   ⚠️ .raw: None")
print(f"   .layers: {list(adata.layers.keys())}")

# ===== 3. Filter Small Batches (BEFORE expensive computation) =====
print("\n===== Batch Filtering =====")
batch_counts = adata.obs[BATCH_KEY].value_counts()
min_batch_size = batch_counts.min()
print(f"Smallest batch: {min_batch_size} cells ({batch_counts.idxmin()})")

# ⭐ CRITICAL: N_NEIGHBORS minimum is 4
N_NEIGHBORS = max(4, min(5, min_batch_size // 2))
print(f"Initial neighbors_within_batch = {N_NEIGHBORS}")

# Remove batches smaller than N_NEIGHBORS + 1
min_required = N_NEIGHBORS + 1
small_batches = batch_counts[batch_counts < min_required].index

if len(small_batches) > 0:
    print(f"⚠️ Removing {len(small_batches)} small batches (<{min_required} cells)")
    print(f"   Removed studies:")
    for study in small_batches:
        print(f"      - {study}: {batch_counts[study]} cells")
    
    # Filter adata (raw will be automatically filtered)
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    
    # Recalculate parameters
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    min_batch_size = batch_counts.min()
    N_NEIGHBORS = max(4, min(5, min_batch_size // 2))
    print(f"\n   After filtering:")
    print(f"   - Remaining batches: {len(batch_counts)}")
    print(f"   - Remaining cells: {adata.n_obs:,}")
    print(f"   - Smallest batch: {min_batch_size} cells")
    print(f"   - Updated neighbors_within_batch: {N_NEIGHBORS}")

# ===== 4. Highly Variable Genes Selection =====
print("\n===== HVG Selection =====")
if 'highly_variable' not in adata.var.columns:
    print(f"Selecting {N_HVG} highly variable genes...")
    
    # Try batch-aware HVG first
    try:
        sc.pp.highly_variable_genes(
            adata, 
            n_top_genes=N_HVG,
            batch_key=BATCH_KEY,
            subset=False  # Don't subset yet
        )
        hvg_method = "batch-aware"
        print(f"   ✅ Batch-aware HVG selection successful")
    except Exception as e:
        print(f"   ⚠️ Batch-aware HVG failed: {e}")
        print(f"   Using standard HVG selection...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=N_HVG,
            subset=False
        )
        hvg_method = "standard"
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"   Selected {n_hvg} HVGs ({hvg_method} method)")
else:
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Using existing HVG selection: {n_hvg} genes")

# ===== 5. PCA on HVG =====
print("\n===== PCA =====")
print(f"Running PCA on {n_hvg} HVGs...")

# Subset to HVG for PCA (more efficient)
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')

# Copy PCA results back to full adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
adata.varm['PCs'] = np.zeros((adata.n_vars, N_PCS))
adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']

# Clean up
del adata_hvg

pca_variance = adata.uns['pca']['variance_ratio']
print(f"✅ PCA completed: {N_PCS} components")
print(f"   PC1-10 explain {pca_variance[:10].sum()*100:.1f}% variance")
print(f"   PC1-30 explain {pca_variance[:30].sum()*100:.1f}% variance")

# ===== 6. BBKNN Integration =====
print("\n===== BBKNN Integration =====")
print(f"Parameters:")
print(f"   batch_key: {BATCH_KEY}")
print(f"   neighbors_within_batch: {N_NEIGHBORS}")
print(f"   n_pcs: {N_PCS}")

sc.external.pp.bbknn(
    adata,
    batch_key=BATCH_KEY,
    neighbors_within_batch=N_NEIGHBORS,
    n_pcs=N_PCS,
    trim=None
)
print("✅ BBKNN integration completed")
print("   Results saved in:")
print("   - adata.obsp['connectivities']: neighbor graph connectivity matrix")
print("   - adata.obsp['distances']: neighbor graph distance matrix")
print("   - adata.uns['neighbors']: neighbor graph parameters")

# ===== 7. UMAP =====
print("\n===== UMAP =====")
print(f"Parameters: min_dist={UMAP_MIN_DIST}, spread={UMAP_SPREAD}")
sc.tl.umap(
    adata, 
    min_dist=UMAP_MIN_DIST,
    spread=UMAP_SPREAD
)
print("✅ UMAP completed")
print("   Results saved in:")
print("   - adata.obsm['X_umap']: UMAP coordinates (n_cells × 2)")

# ===== 8. Leiden Clustering =====
print("\n===== Leiden Clustering =====")
sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added='leiden_res1.0')
n_clusters = adata.obs['leiden_res1.0'].nunique()
print(f"✅ Found {n_clusters} clusters at resolution {LEIDEN_RESOLUTION}")
print("   Results saved in:")
print("   - adata.obs['leiden_res1.0']: cluster assignments")

# ===== 9. Visualization =====
print("\n===== Generating Visualizations =====")

# Plot 1: Study vs Clusters
fig, axes = plt.subplots(1, 2, figsize=(18, 7))

sc.pl.umap(
    adata, 
    color=BATCH_KEY,
    ax=axes[0],
    title='UMAP - Colored by Study (Batch Mixing)',
    show=False,
    frameon=False,
    size=20
)

sc.pl.umap(
    adata,
    color='leiden_res1.0',
    ax=axes[1],
    title='UMAP - Colored by Leiden Clusters',
    legend_loc='on data',
    legend_fontsize=10,
    legend_fontweight='bold',
    show=False,
    frameon=False,
    size=20
)

plt.tight_layout()
plt.savefig(OUTPUT_PDF1, dpi=300, bbox_inches='tight')
print(f"✅ Saved: {os.path.basename(OUTPUT_PDF1)}")
plt.close()

# Plot 2: Multi-panel overview
plot_cols = [BATCH_KEY, 'leiden_res1.0']
if 'n_genes' in adata.obs.columns:
    plot_cols.append('n_genes')
if 'n_counts' in adata.obs.columns:
    plot_cols.append('n_counts')

fig = sc.pl.umap(
    adata,
    color=plot_cols,
    ncols=2,
    frameon=False,
    size=15,
    return_fig=True,
    show=False
)
fig.savefig(OUTPUT_PDF2, dpi=300, bbox_inches='tight')
print(f"✅ Saved: {os.path.basename(OUTPUT_PDF2)}")
plt.close()

# ===== 10. Save Results =====
print("\n===== Saving Results =====")

# Store parameters in uns
adata.uns['bbknn_params'] = {
    'batch_key': BATCH_KEY,
    'neighbors_within_batch': N_NEIGHBORS,
    'n_pcs': N_PCS,
    'n_hvg': N_HVG,
    'umap_min_dist': UMAP_MIN_DIST,
    'umap_spread': UMAP_SPREAD,
    'leiden_resolution': LEIDEN_RESOLUTION
}

adata.write_h5ad(OUTPUT_H5AD, compression='gzip', compression_opts=9)
print(f"✅ Saved: {os.path.basename(OUTPUT_H5AD)}")

# ===== 11. Final Summary =====
print("\n" + "="*70)
print("INTEGRATION SUMMARY")
print("="*70)
print(f"Total cells:              {adata.n_obs:,}")
print(f"Total genes:              {adata.n_vars:,}")
print(f"HVG used for PCA:         {n_hvg}")
print(f"Number of studies:        {adata.obs[BATCH_KEY].nunique()}")
print(f"Leiden clusters:          {n_clusters}")
print(f"neighbors_within_batch:   {N_NEIGHBORS}")
print(f"UMAP min_dist:            {UMAP_MIN_DIST}")
print(f"UMAP spread:              {UMAP_SPREAD}")

print("\n" + "="*70)
print("DATA STRUCTURE IN OUTPUT H5AD")
print("="*70)
print("Main data:")
if adata.raw is not None:
    print(f"  .raw.X              : raw counts ({adata.raw.n_obs:,} cells × {adata.raw.n_vars:,} genes)")
else:
    print(f"  .raw                : None")
print(f"  .X                  : log1p normalized (sparse={sparse.issparse(adata.X)})")
print(f"  .layers['log1p']    : log1p normalized (backup)")

print("\nDimensionality reduction:")
print(f"  .obsm['X_pca']      : PCA coordinates ({adata.obsm['X_pca'].shape})")
print(f"  .obsm['X_umap']     : UMAP coordinates ({adata.obsm['X_umap'].shape})")
print(f"  .varm['PCs']        : Principal components loadings ({adata.varm['PCs'].shape})")

print("\nBBKNN neighbor graph:")
print(f"  .obsp['connectivities'] : Neighbor connectivity matrix ({adata.obsp['connectivities'].shape})")
print(f"  .obsp['distances']      : Neighbor distance matrix ({adata.obsp['distances'].shape})")
print(f"  .uns['neighbors']       : Neighbor graph parameters")

print("\nAnnotations:")
print(f"  .obs['leiden_res1.0']   : Leiden clustering results")
print(f"  .obs['study']           : Batch/study information")
print(f"  .obs['n_genes']         : Number of genes per cell")
print(f"  .obs['n_counts']        : Total UMI counts per cell")

print("\nGene annotations:")
print(f"  .var['highly_variable'] : HVG selection ({n_hvg} genes)")
print("="*70)

print("\nCells per study:")
study_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
for study, count in study_counts.items():
    print(f"  {study:30s} {count:6,} cells")

print("\nCells per cluster:")
cluster_counts = adata.obs['leiden_res1.0'].value_counts().sort_values(ascending=False)
for cluster, count in cluster_counts.items():
    pct = count / adata.n_obs * 100
    print(f"  Cluster {cluster:2s}: {count:6,} cells ({pct:5.1f}%)")

print("="*70)

print(f"\n📁 All outputs saved to: {OUTPUT_DIR}")
print("Files generated:")
print(f"  1. {os.path.basename(OUTPUT_H5AD)}")
print(f"  2. {os.path.basename(OUTPUT_PDF1)}")
print(f"  3. {os.path.basename(OUTPUT_PDF2)}")
print("\n✅ Pipeline completed successfully!")