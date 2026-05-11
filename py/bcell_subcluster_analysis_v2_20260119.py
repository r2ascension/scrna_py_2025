#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B Cell Subcluster Analysis Pipeline v2.0 (PRODUCTION)
======================================================

Purpose:
--------
Fine-grained subclustering of B cell populations based on scANVI annotations.
Identifies functional states within each major B cell type (Naive, Memory, 
Plasma, etc.) while maintaining consistent normalization across all cells.

Workflow:
---------
1. Load scANVI-corrected B cell data
2. Global gene filtering (MT/ribosomal/ENSG/unannotated)
3. Global normalization (critical for cross-celltype comparison)
4. Within each celltype: HVG → PCA → BBKNN → clustering → markers
5. Merge subclusters with hierarchical naming (e.g., Memory_B_c0)
6. Global UMAP visualization using scANVI latent space

Input:
------
- adata_bcell_FINAL_corrected_20260114.h5ad
  Required fields:
    - cell_type_scanvi_corrected (scANVI predictions)
    - X_scanvi_corrected (scANVI latent space)
    - Sample (batch key)

Output:
-------
- adata_bcell_subclustered_FINAL_v2_20260119.h5ad
- Per-celltype marker analysis and figures
- Summary UMAP with hierarchical annotations

Author: r2end
Date: 2026-01-19
Memory: < 30GB RAM
Runtime: ~30-45 min

Key Design Principles:
----------------------
✅ GLOBAL normalize BEFORE subsetting (ensures consistent scale)
✅ Filter low-quality genes (MT/ribo/ENSG) at start
✅ Low resolution clustering (avoid over-fragmentation)
✅ BBKNN only if >=3 batches per celltype
✅ Use scANVI latent for final unified UMAP
✅ Hierarchical naming: cell_type_L2 (scANVI) → cell_type_L3 (with subclusters)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gc
import warnings
from scipy import sparse
import re
from datetime import datetime
warnings.filterwarnings('ignore')

# Scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, dpi_save=300, frameon=False)
sc.settings.n_jobs = 48

# ============================================================================
# CELL 1: Configuration
# ============================================================================

print("="*80)
print("B CELL SUBCLUSTER ANALYSIS PIPELINE v2.0")
print("Date:", pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))
print("="*80)

# ===== Paths =====
BASE_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/bcell_analysis") 
INPUT_FILE = Path('/home/h2048/data/py/0111/celltypist_bcell/adata_bcell_FINAL_corrected_20260114.h5ad')
OUTPUT_DIR = BASE_DIR / "results/subcluster_v2_20260119"
FIG_DIR = OUTPUT_DIR / "figures"
TABLE_DIR = OUTPUT_DIR / "tables"

for d in [OUTPUT_DIR, FIG_DIR, TABLE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

print(f"\n📂 Input: {INPUT_FILE}")
print(f"📂 Output: {OUTPUT_DIR}")

# ===== Parameters =====
BATCH_KEY = 'dataset'  # Batch key for BBKNN
CELLTYPE_KEY = 'cell_type_scanvi_corrected'  # scANVI annotations
SCANVI_LATENT_KEY = 'X_scanvi_corrected'  # For global UMAP

# Subclustering parameters (per celltype)
SUBCLUSTER_RESOLUTIONS = {
    'default': 0.3,  # Most celltypes
    'Plasma': 0.2,   # Plasma cells are more homogeneous
    'Naive': 0.3,
    'Memory': 0.4    # Memory B might have more states
}

MIN_CELLS_FOR_SUBCLUSTER = 100  # Skip if too few cells
MIN_BATCHES_FOR_BBKNN = 3  # Require >=3 batches for BBKNN

# Gene filtering
FILTER_MT = True
FILTER_RIBO = True
FILTER_ENSG = True
FILTER_UNANNOTATED = True

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

print(f"\n⚙️  Analysis parameters:")
print(f"  - Batch key: {BATCH_KEY}")
print(f"  - Celltype key: {CELLTYPE_KEY}")
print(f"  - Min cells for subcluster: {MIN_CELLS_FOR_SUBCLUSTER}")
print(f"  - Min batches for BBKNN: {MIN_BATCHES_FOR_BBKNN}")

# ============================================================================
# CELL 2: Load Data and Inspect
# ============================================================================

print("\n" + "="*80)
print("LOADING DATA")
print("="*80)

print(f"\nReading {INPUT_FILE.name}...")
adata = sc.read_h5ad(INPUT_FILE)

print(f"✓ Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"\nData structure:")
print(f"  - .X dtype: {adata.X.dtype}, sparse: {sparse.issparse(adata.X)}")
print(f"  - .layers: {list(adata.layers.keys())}")
print(f"  - .obsm keys: {list(adata.obsm.keys())}")

# Check required fields
required_fields = [CELLTYPE_KEY, BATCH_KEY, SCANVI_LATENT_KEY]
missing = [f for f in required_fields if f not in adata.obs.columns and f not in adata.obsm.keys()]
if missing:
    raise ValueError(f"Missing required fields: {missing}")

print(f"\n✓ Required fields present")

# Inspect cell types
print(f"\nCell type distribution ({CELLTYPE_KEY}):")
celltype_counts = adata.obs[CELLTYPE_KEY].value_counts()
for ct, count in celltype_counts.items():
    pct = count / len(adata) * 100
    print(f"  {ct:30s}: {count:6,} cells ({pct:5.2f}%)")

total_celltypes = len(celltype_counts)
print(f"\nTotal cell types: {total_celltypes}")

# Check batch distribution
print(f"\nBatch distribution:")
batch_counts = adata.obs[BATCH_KEY].value_counts()
print(f"  Total batches: {len(batch_counts)}")
print(f"  Cells per batch: {batch_counts.min():.0f} - {batch_counts.max():.0f}")

gc.collect()

# ============================================================================
# CELL 3: Global Gene Filtering
# ============================================================================

print("\n" + "="*80)
print("GLOBAL GENE FILTERING")
print("="*80)

print(f"\nBefore filtering: {adata.shape[1]:,} genes")

# Annotate gene types
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
adata.var['ensg'] = adata.var_names.str.startswith('ENSG')
adata.var['malat1'] = adata.var_names.str.upper() == 'MALAT1'

# Identify unannotated genes (heuristic: gene names with numbers/underscores but not starting with known patterns)
# This is conservative - adjust based on your data
adata.var['unannotated'] = (
    ~adata.var_names.str.match('^[A-Z][A-Z0-9-]+$') |  # Not standard gene names
    adata.var_names.str.contains('_') |  # Contains underscore
    (adata.var_names.str.len() > 15)  # Suspiciously long names
) & ~adata.var['ensg']  # But not already flagged as ENSG

# Count genes to filter
n_mt = adata.var['mt'].sum()
n_ribo = adata.var['ribo'].sum()
n_ensg = adata.var['ensg'].sum()
n_malat1 = adata.var['malat1'].sum()
n_unannotated = adata.var['unannotated'].sum()

print(f"\nGenes to filter:")
print(f"  - Mitochondrial: {n_mt}")
print(f"  - Ribosomal (RPS*/RPL*): {n_ribo}")
print(f"  - ENSG genes: {n_ensg}")
print(f"  - MALAT1: {n_malat1}")
print(f"  - Unannotated: {n_unannotated}")

# Create keep mask
keep_genes = ~(adata.var['mt'] | adata.var['malat1'])

if FILTER_RIBO:
    keep_genes = keep_genes & ~adata.var['ribo']
    print(f"  ✓ Filtering ribosomal genes")

if FILTER_ENSG:
    keep_genes = keep_genes & ~adata.var['ensg']
    print(f"  ✓ Filtering ENSG genes")

if FILTER_UNANNOTATED:
    keep_genes = keep_genes & ~adata.var['unannotated']
    print(f"  ✓ Filtering unannotated genes")

# Apply filter
print(f"\nFiltering genes...")
adata = adata[:, keep_genes].copy()
print(f"✓ After filtering: {adata.shape[1]:,} genes ({keep_genes.sum() / len(keep_genes) * 100:.1f}% retained)")

# Additional basic filtering (genes in at least 3 cells)
sc.pp.filter_genes(adata, min_cells=3)
print(f"✓ After min_cells filter: {adata.shape[1]:,} genes")

gc.collect()

# ============================================================================
# CELL 4: Global Normalization (CRITICAL!)
# ============================================================================

print("\n" + "="*80)
print("⭐ GLOBAL NORMALIZATION (CRITICAL)")
print("="*80)

print("""
Why global normalization?
-------------------------
If we normalize each celltype separately, they will have different scales,
making cross-celltype comparison invalid and merged UMAP misleading.

Solution: Normalize ALL cells together BEFORE subsetting.
""")

# Check if we have counts
if 'counts' in adata.layers:
    print(f"✓ Found 'counts' layer, using it for normalization")
    counts_layer = 'counts'
elif sparse.issparse(adata.X) and adata.X.dtype in [np.int32, np.int64, np.float32]:
    print(f"⚠️  No 'counts' layer, assuming .X is raw counts")
    adata.layers['counts'] = adata.X.copy()
    counts_layer = 'counts'
else:
    raise ValueError("Cannot find raw counts! Check data structure.")

# Normalize
print(f"\nNormalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata, target_sum=1e4, layer=counts_layer)

print(f"Log-transforming...")
sc.pp.log1p(adata)

# Store normalized data
adata.layers['log1p'] = adata.X.copy()
print(f"✓ Stored normalized data in .layers['log1p']")

# ⭐ CRITICAL: Preserve full gene set in .raw
print(f"\n⭐ Preserving full gene set in .raw...")
adata.raw = sc.AnnData(
    X=adata.layers['counts'],  # Use counts, shared memory
    obs=adata.obs.copy(),
    var=adata.var.copy()
)
print(f"✓ Saved {adata.raw.shape[1]:,} genes to .raw")

print(f"\n{'='*60}")
print(f"GLOBAL NORMALIZATION COMPLETE")
print(f"All cells now on same scale: 10,000 counts/cell, log-transformed")
print(f"{'='*60}")

gc.collect()

# ============================================================================
# CELL 5: Per-Celltype Subclustering Loop
# ============================================================================

print("\n" + "="*80)
print("PER-CELLTYPE SUBCLUSTERING")
print("="*80)

# Get list of celltypes
celltypes = adata.obs[CELLTYPE_KEY].unique()
celltypes = [ct for ct in celltypes if not pd.isna(ct)]  # Remove NaN if any
celltypes = sorted(celltypes)

print(f"\nCell types to analyze: {len(celltypes)}")
for ct in celltypes:
    print(f"  - {ct}")

# Storage for results
subcluster_results = {}
subcluster_labels = pd.Series(index=adata.obs.index, dtype='object')
subcluster_ids = pd.Series(index=adata.obs.index, dtype='object')

# Loop through each celltype
for celltype in celltypes:
    print(f"\n{'='*80}")
    print(f"ANALYZING: {celltype}")
    print(f"{'='*80}")
    
    # Subset cells
    mask = adata.obs[CELLTYPE_KEY] == celltype
    adata_sub = adata[mask].copy()
    
    n_cells = adata_sub.shape[0]
    n_batches = adata_sub.obs[BATCH_KEY].nunique()
    
    print(f"\n📊 Statistics:")
    print(f"  - Cells: {n_cells:,}")
    print(f"  - Batches: {n_batches}")
    print(f"  - Genes: {adata_sub.shape[1]:,}")
    
    # Check if sufficient cells
    if n_cells < MIN_CELLS_FOR_SUBCLUSTER:
        print(f"\n⚠️  Too few cells (<{MIN_CELLS_FOR_SUBCLUSTER}), skipping subclustering")
        print(f"  Keeping as single cluster: {celltype}_c0")
        
        # Assign single cluster
        subcluster_labels[mask] = f"{celltype}_c0"
        subcluster_ids[mask] = "0"
        continue
    
    # HVG selection (celltype-specific)
    n_hvg = min(2000, adata_sub.shape[1] // 2)  # Adaptive HVG count
    print(f"\n🔍 Selecting {n_hvg} HVGs...")
    
    sc.pp.highly_variable_genes(
        adata_sub,
        n_top_genes=n_hvg,
        layer='log1p',
        subset=False
    )
    
    n_hvg_selected = adata_sub.var['highly_variable'].sum()
    print(f"  ✓ Selected {n_hvg_selected} HVGs")
    
    # Subset to HVGs for analysis
    adata_hvg = adata_sub[:, adata_sub.var['highly_variable']].copy()
    print(f"  Working with {adata_hvg.shape[1]} HVGs")
    
    # PCA
    n_pcs = min(30, adata_hvg.shape[1] - 1, n_cells - 1)
    print(f"\n📐 Computing PCA ({n_pcs} components)...")
    sc.tl.pca(adata_hvg, n_comps=n_pcs, svd_solver='arpack', random_state=RANDOM_STATE)
    print(f"  ✓ PCA done")
    
    # Batch correction strategy
    if n_batches >= MIN_BATCHES_FOR_BBKNN:
        print(f"\n🔗 Using BBKNN (sufficient batches: {n_batches})")
        
        # Check for small batches
        batch_counts = adata_hvg.obs[BATCH_KEY].value_counts()
        small_batches = batch_counts[batch_counts < 3].index
        
        if len(small_batches) > 0:
            print(f"  ⚠️  Removing {len(small_batches)} small batches (<3 cells)")
            adata_hvg = adata_hvg[~adata_hvg.obs[BATCH_KEY].isin(small_batches)].copy()
            print(f"  Remaining cells: {adata_hvg.shape[0]}")
        
        # BBKNN
        neighbors_within = 3 if n_batches >= 5 else 5  # Adaptive
        sc.external.pp.bbknn(
            adata_hvg,
            batch_key=BATCH_KEY,
            neighbors_within_batch=neighbors_within,
            n_pcs=n_pcs,
            trim=None
        )
        print(f"  ✓ BBKNN completed (neighbors_within_batch={neighbors_within})")
    else:
        print(f"\n🔗 Using standard neighbors (only {n_batches} batches)")
        sc.pp.neighbors(adata_hvg, n_pcs=n_pcs, random_state=RANDOM_STATE)
        print(f"  ✓ Neighbors computed")
    
    # Leiden clustering
    resolution = SUBCLUSTER_RESOLUTIONS.get(celltype, SUBCLUSTER_RESOLUTIONS['default'])
    print(f"\n🎯 Leiden clustering (resolution={resolution})...")
    
    sc.tl.leiden(adata_hvg, resolution=resolution, key_added='subcluster', 
                 random_state=RANDOM_STATE)
    
    n_subclusters = adata_hvg.obs['subcluster'].nunique()
    print(f"  ✓ Identified {n_subclusters} subclusters")
    
    # Print subcluster sizes
    subcluster_counts = adata_hvg.obs['subcluster'].value_counts().sort_index()
    print(f"\n  Subcluster distribution:")
    for sc_id, count in subcluster_counts.items():
        pct = count / len(adata_hvg) * 100
        print(f"    c{sc_id}: {count:5,} cells ({pct:5.2f}%)")
    
    # UMAP (within celltype)
    print(f"\n🗺️  Computing UMAP...")
    sc.tl.umap(adata_hvg, min_dist=0.3, random_state=RANDOM_STATE)
    print(f"  ✓ UMAP computed")
    
    # Marker genes
    print(f"\n🔬 Computing marker genes...")
    sc.tl.rank_genes_groups(
        adata_hvg,
        groupby='subcluster',
        method='wilcoxon',
        use_raw=False,  # Use adata_hvg (HVGs only)
        key_added='rank_genes_subcluster'
    )
    print(f"  ✓ Marker genes computed")
    
    # Save results
    subcluster_results[celltype] = adata_hvg
    
    # Transfer subcluster labels back to main adata
    # Handle potential cell filtering during BBKNN
    common_cells = adata_hvg.obs.index.intersection(adata_sub.obs.index)
    
    for cell_id in common_cells:
        sc_id = adata_hvg.obs.loc[cell_id, 'subcluster']
        subcluster_labels[cell_id] = f"{celltype}_c{sc_id}"
        subcluster_ids[cell_id] = sc_id
    
    # For cells filtered out during BBKNN, assign to cluster 0
    filtered_cells = adata_sub.obs.index.difference(common_cells)
    if len(filtered_cells) > 0:
        print(f"  ℹ️  {len(filtered_cells)} cells filtered during BBKNN, assigning to c0")
        for cell_id in filtered_cells:
            subcluster_labels[cell_id] = f"{celltype}_c0"
            subcluster_ids[cell_id] = "0"
    
    # Generate per-celltype figure
    print(f"\n📊 Generating figures...")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # UMAP by subcluster
    sc.pl.umap(adata_hvg, color='subcluster', ax=axes[0], show=False,
               title=f'{celltype} Subclusters', size=30, legend_loc='on data')
    
    # UMAP by batch
    sc.pl.umap(adata_hvg, color=BATCH_KEY, ax=axes[1], show=False,
               title='Batch Distribution', size=30)
    
    # Top marker heatmap
    if n_subclusters > 1:
        try:
            sc.pl.rank_genes_groups_heatmap(
                adata_hvg,
                n_genes=5,
                groupby='subcluster',
                key='rank_genes_subcluster',
                show_gene_labels=True,
                cmap='RdBu_r',
                ax=axes[2],
                show=False
            )
        except:
            axes[2].text(0.5, 0.5, 'Marker heatmap\nunavailable',
                        ha='center', va='center', transform=axes[2].transAxes)
    else:
        axes[2].text(0.5, 0.5, 'Only 1 cluster',
                    ha='center', va='center', transform=axes[2].transAxes)
    
    plt.tight_layout()
    fig_file = FIG_DIR / f'{celltype.replace(" ", "_")}_subclustering.pdf'
    plt.savefig(fig_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved: {fig_file.name}")
    
    # Export marker genes
    if n_subclusters > 1:
        result = adata_hvg.uns['rank_genes_subcluster']
        groups = result['names'].dtype.names
        
        marker_df_list = []
        for group in groups:
            group_df = pd.DataFrame({
                'celltype': celltype,
                'subcluster': group,
                'gene': result['names'][group],
                'logfoldchanges': result['logfoldchanges'][group],
                'pvals': result['pvals'][group],
                'pvals_adj': result['pvals_adj'][group]
            })
            marker_df_list.append(group_df)
        
        marker_df = pd.concat(marker_df_list, ignore_index=True)
        marker_df = marker_df[marker_df['pvals_adj'] < 0.05]
        
        marker_file = TABLE_DIR / f'{celltype.replace(" ", "_")}_markers.csv'
        marker_df.to_csv(marker_file, index=False)
        print(f"  ✓ Exported markers: {marker_file.name}")
    
    # Cleanup
    del adata_sub, adata_hvg
    gc.collect()
    
    print(f"\n✅ {celltype} analysis complete")

print(f"\n{'='*80}")
print(f"ALL CELLTYPES ANALYZED")
print(f"{'='*80}")

# ============================================================================
# CELL 6: Add Subcluster Labels to Main AnnData
# ============================================================================

print("\n" + "="*80)
print("MERGING SUBCLUSTER LABELS")
print("="*80)

# Add hierarchical annotations
adata.obs['cell_type_L2'] = adata.obs[CELLTYPE_KEY].copy()  # scANVI level
adata.obs['cell_type_L3'] = subcluster_labels  # With subclusters
adata.obs['subcluster_id'] = subcluster_ids  # Just the cluster ID

print(f"\n📊 Hierarchical annotation levels:")
print(f"  - Level 1: Major lineage (not computed here)")
print(f"  - Level 2 (cell_type_L2): {adata.obs['cell_type_L2'].nunique()} types")
print(f"  - Level 3 (cell_type_L3): {adata.obs['cell_type_L3'].nunique()} subclusters")

print(f"\nLevel 3 distribution:")
l3_counts = adata.obs['cell_type_L3'].value_counts()
for subcluster, count in l3_counts.items():
    pct = count / len(adata) * 100
    print(f"  {subcluster:35s}: {count:6,} ({pct:5.2f}%)")

# ============================================================================
# CELL 7: Global UMAP Visualization
# ============================================================================

print("\n" + "="*80)
print("GLOBAL UMAP VISUALIZATION")
print("="*80)

print(f"\n🗺️  Computing global UMAP using scANVI latent space...")

# Use scANVI latent for unified visualization
if SCANVI_LATENT_KEY in adata.obsm.keys():
    print(f"  Using {SCANVI_LATENT_KEY} for neighbors")
    sc.pp.neighbors(
        adata,
        use_rep=SCANVI_LATENT_KEY,
        n_neighbors=15,
        random_state=RANDOM_STATE
    )
else:
    print(f"  ⚠️  No scANVI latent found, using PCA on HVGs")
    # Fallback: compute HVGs on full dataset
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, layer='log1p', subset=False)
    adata_hvg_global = adata[:, adata.var['highly_variable']].copy()
    sc.tl.pca(adata_hvg_global, n_comps=50, random_state=RANDOM_STATE)
    adata.obsm['X_pca'] = adata_hvg_global.obsm['X_pca']
    sc.pp.neighbors(adata, n_pcs=50, random_state=RANDOM_STATE)
    del adata_hvg_global
    gc.collect()

print(f"  Computing UMAP...")
sc.tl.umap(adata, min_dist=0.3, random_state=RANDOM_STATE)
print(f"  ✓ Global UMAP computed: {adata.obsm['X_umap'].shape}")

# Generate comprehensive UMAP figure
print(f"\n📊 Generating global UMAP figure...")

fig, axes = plt.subplots(2, 3, figsize=(20, 13))
axes = axes.flatten()

# 1. Level 2 (scANVI types)
sc.pl.umap(adata, color='cell_type_L2', ax=axes[0], show=False,
           title='Level 2: scANVI Cell Types', size=10,
           legend_loc='right margin', legend_fontsize=8)

# 2. Level 3 (with subclusters)
sc.pl.umap(adata, color='cell_type_L3', ax=axes[1], show=False,
           title='Level 3: With Subclusters', size=10,
           legend_loc='right margin', legend_fontsize=6)

# 3. Batch
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[2], show=False,
           title='Batch Distribution', size=10)

# 4. scANVI confidence (if available)
if 'scanvi_confidence_corrected' in adata.obs.columns:
    sc.pl.umap(adata, color='scanvi_confidence_corrected', ax=axes[3], show=False,
               title='scANVI Confidence', size=10, cmap='viridis')
else:
    axes[3].text(0.5, 0.5, 'No confidence\ndata available',
                ha='center', va='center', transform=axes[3].transAxes)

# 5. Key B cell marker: CD19
if 'CD19' in adata.raw.var_names:
    sc.pl.umap(adata, color='CD19', ax=axes[4], show=False,
               title='CD19 (Pan B marker)', size=10, use_raw=True, cmap='Reds')
else:
    axes[4].text(0.5, 0.5, 'CD19 not found',
                ha='center', va='center', transform=axes[4].transAxes)

# 6. Key B cell marker: MS4A1 (CD20)
if 'MS4A1' in adata.raw.var_names:
    sc.pl.umap(adata, color='MS4A1', ax=axes[5], show=False,
               title='MS4A1/CD20 (Pan B marker)', size=10, use_raw=True, cmap='Reds')
else:
    axes[5].text(0.5, 0.5, 'MS4A1 not found',
                ha='center', va='center', transform=axes[5].transAxes)

plt.tight_layout()
umap_file = FIG_DIR / 'global_umap_summary.pdf'
plt.savefig(umap_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: {umap_file.name}")

# Generate separate high-res UMAP for Level 3
fig, ax = plt.subplots(figsize=(14, 10))
sc.pl.umap(adata, color='cell_type_L3', ax=ax, show=False,
           title='B Cell Subclusters (Level 3)', size=15,
           legend_loc='right margin', legend_fontsize=10, frameon=False)
plt.tight_layout()
umap_l3_file = FIG_DIR / 'umap_level3_highres.pdf'
plt.savefig(umap_l3_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✓ Saved: {umap_l3_file.name}")

# ============================================================================
# CELL 8: Summary Statistics and Tables
# ============================================================================

print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

# Cell type proportions
print(f"\n📊 Cell type proportions:")

# Level 2
l2_counts = adata.obs['cell_type_L2'].value_counts()
l2_df = pd.DataFrame({
    'cell_type': l2_counts.index,
    'count': l2_counts.values,
    'percentage': (l2_counts.values / len(adata) * 100).round(2)
}).sort_values('count', ascending=False)

l2_file = TABLE_DIR / 'celltype_L2_counts.csv'
l2_df.to_csv(l2_file, index=False)
print(f"  Level 2: {l2_file.name}")

# Level 3
l3_counts = adata.obs['cell_type_L3'].value_counts()
l3_df = pd.DataFrame({
    'cell_type': l3_counts.index,
    'count': l3_counts.values,
    'percentage': (l3_counts.values / len(adata) * 100).round(2)
}).sort_values('count', ascending=False)

l3_file = TABLE_DIR / 'celltype_L3_counts.csv'
l3_df.to_csv(l3_file, index=False)
print(f"  Level 3: {l3_file.name}")

# Cross-tabulation: Level 2 vs Level 3
crosstab = pd.crosstab(adata.obs['cell_type_L2'], adata.obs['cell_type_L3'])
crosstab_file = TABLE_DIR / 'L2_vs_L3_crosstab.csv'
crosstab.to_csv(crosstab_file)
print(f"  Crosstab: {crosstab_file.name}")

# Batch × celltype distribution
if 'disease_status' in adata.obs.columns or 'group' in adata.obs.columns:
    group_col = 'disease_status' if 'disease_status' in adata.obs.columns else 'group'
    
    print(f"\n📊 Cell type proportions by {group_col}:")
    
    # Level 2 by group
    prop_l2 = pd.crosstab(
        adata.obs[group_col],
        adata.obs['cell_type_L2'],
        normalize='index'
    ) * 100
    
    prop_l2_file = TABLE_DIR / f'celltype_L2_by_{group_col}.csv'
    prop_l2.to_csv(prop_l2_file)
    print(f"  Level 2: {prop_l2_file.name}")
    
    # Visualize
    fig, ax = plt.subplots(figsize=(12, 6))
    prop_l2.T.plot(kind='bar', ax=ax, width=0.8)
    ax.set_ylabel('Percentage (%)')
    ax.set_xlabel('Cell Type (Level 2)')
    ax.set_title(f'Cell Type Proportions by {group_col}')
    ax.legend(title=group_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    prop_fig = FIG_DIR / f'celltype_proportions_by_{group_col}.pdf'
    plt.savefig(prop_fig, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Figure: {prop_fig.name}")

# ============================================================================
# CELL 9: Export Per-Celltype Subcluster Results
# ============================================================================

print("\n" + "="*80)
print("EXPORTING PER-CELLTYPE RESULTS")
print("="*80)

# Create a summary table of subclusters
subcluster_summary = []

for celltype in celltypes:
    if celltype in subcluster_results:
        adata_ct = subcluster_results[celltype]
        
        for sc_id in adata_ct.obs['subcluster'].unique():
            mask = adata_ct.obs['subcluster'] == sc_id
            n_cells = mask.sum()
            
            subcluster_summary.append({
                'celltype_L2': celltype,
                'subcluster_id': sc_id,
                'subcluster_label': f"{celltype}_c{sc_id}",
                'n_cells': n_cells,
                'percentage_within_celltype': n_cells / len(adata_ct) * 100,
                'percentage_total': n_cells / len(adata) * 100
            })

summary_df = pd.DataFrame(subcluster_summary)
summary_df = summary_df.sort_values(['celltype_L2', 'subcluster_id'])

summary_file = TABLE_DIR / 'subcluster_summary.csv'
summary_df.to_csv(summary_file, index=False)
print(f"✓ Subcluster summary: {summary_file.name}")

print(f"\nSubcluster summary:")
print(summary_df.to_string(index=False))

# ============================================================================
# CELL 10: Final Save
# ============================================================================

print("\n" + "="*80)
print("SAVING FINAL RESULTS")
print("="*80)

# Final output file
output_file = OUTPUT_DIR / "adata_bcell_subclustered_FINAL_v2_20260119.h5ad"

# Add analysis metadata
adata.uns['subcluster_analysis'] = {
    'date': '2026-01-19',
    'version': 'v2.0',
    'workflow': 'global_normalize_then_subcluster',
    'celltype_key': CELLTYPE_KEY,
    'batch_key': BATCH_KEY,
    'n_celltypes': len(celltypes),
    'n_subclusters_total': adata.obs['cell_type_L3'].nunique(),
    'gene_filtering': {
        'mt': FILTER_MT,
        'ribo': FILTER_RIBO,
        'ensg': FILTER_ENSG,
        'unannotated': FILTER_UNANNOTATED
    },
    'subcluster_resolutions': SUBCLUSTER_RESOLUTIONS,
    'min_cells_for_subcluster': MIN_CELLS_FOR_SUBCLUSTER,
    'min_batches_for_bbknn': MIN_BATCHES_FOR_BBKNN
}

# Save
print(f"\nSaving to: {output_file}")
adata.write_h5ad(output_file, compression='gzip', compression_opts=9)

size_gb = output_file.stat().st_size / 1e9
print(f"✓ File saved: {size_gb:.2f} GB")

# Export final annotations
print(f"\nExporting final annotations...")

annotation_cols = [BATCH_KEY, 'cell_type_L2', 'cell_type_L3', 'subcluster_id']
if 'scanvi_confidence_corrected' in adata.obs.columns:
    annotation_cols.append('scanvi_confidence_corrected')
if 'disease_status' in adata.obs.columns:
    annotation_cols.append('disease_status')

final_annotations = adata.obs[annotation_cols].copy()
annot_file = TABLE_DIR / 'final_annotations.csv'
final_annotations.to_csv(annot_file)
print(f"✓ Annotations exported: {annot_file.name}")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("✅ ANALYSIS COMPLETE!")
print("="*80)

print(f"""
📁 Output Files
===============

Main Results:
  - {output_file.name}
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
  - [celltype]_markers.csv (per celltype)

📊 Analysis Summary
===================
  - Total cells: {adata.shape[0]:,}
  - Genes (after filtering): {adata.shape[1]:,}
  - Cell types (Level 2): {adata.obs['cell_type_L2'].nunique()}
  - Subclusters (Level 3): {adata.obs['cell_type_L3'].nunique()}
  - Batches: {adata.obs[BATCH_KEY].nunique()}

🔑 Key Points
=============
  1. ✅ Global normalization ensures consistent scale
  2. ✅ Gene filtering removed {n_mt + n_ribo + n_ensg} low-quality genes
  3. ✅ Per-celltype subclustering preserves biological context
  4. ✅ BBKNN used when >=3 batches available
  5. ✅ Global UMAP uses scANVI latent space for unified view

📝 Next Steps
=============
  1. Review per-celltype subclustering figures
  2. Examine marker genes for each subcluster
  3. Annotate subclusters with functional states:
     - e.g., Memory_B_c0 → "Resting memory B"
     - e.g., Memory_B_c1 → "Activated memory B"
  4. Perform disease vs healthy comparisons
  5. Investigate subcluster-specific gene programs

""")

print(f"Output directory: {OUTPUT_DIR}")
print(f"Figures: {FIG_DIR}")
print(f"Tables: {TABLE_DIR}")
print("\n" + "="*80)

# ============================================================================
# END OF PIPELINE
# ============================================================================
