# %% [markdown]
"""
# B Cell Population Analysis - Quick Overview
Author: r2end
Date: 2026-01-14
Purpose: Comprehensive analysis of B cell populations from scRNA-seq data

Analysis Pipeline:
1. Data loading and structure validation
2. Quality metrics overview
3. UMAP visualization (batch, cell type, key markers)
4. Cell type composition analysis
5. Marker gene expression patterns
6. (Optional) Differential expression analysis

Input: adata_bcell_FINAL.h5ad
Expected structure:
- .X: log1p normalized
- .layers['counts']: raw counts
- .layers['log1p']: log1p normalized
- .raw: full gene matrix (if available)
- .obs: metadata with cell type annotations
- .obsm: UMAP, PCA coordinates
"""

# %% [markdown]
# ## 0. Setup and Configuration

# %%
# Import libraries
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))
sc.settings.n_jobs = 48  # Multi-core processing

# %%
# Configuration parameters
INPUT_FILE = '/home/h2048/data/py/0111/celltypist_bcell/adata_bcell_FINAL.h5ad'
OUTPUT_DIR = '/home/h2048/data/py/0111/celltypist_bcell/analysis_outputs'
FIG_DIR = f'{OUTPUT_DIR}/figures'

# Create output directories
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Key annotation columns (adjust based on your data)
BATCH_KEY = 'Sample'  # or 'batch', 'donor', 'patient_id'
CELLTYPE_KEY = 'cell_type'  # or 'predicted_labels', 'leiden', 'cluster'
DISEASE_KEY = 'disease_status'  # or 'group', 'condition'

# Key B cell markers
BCELL_MARKERS = {
    'Pan_B': ['CD19', 'MS4A1', 'CD79A', 'CD79B'],  # MS4A1 = CD20
    'Naive_B': ['TCL1A', 'FCER2', 'IL4R'],
    'Memory_B': ['CD27', 'TNFRSF13B'],  # TNFRSF13B = TACI
    'Plasma': ['MZB1', 'JCHAIN', 'IGHA1', 'IGHG1', 'XBP1', 'SDC1'],  # SDC1 = CD138
    'Activated_B': ['CD69', 'CD83', 'BCL2A1'],
    'GC_B': ['AICDA', 'BCL6', 'CXCR4']  # Germinal center B cells
}

# Flatten marker list for quick plotting
ALL_MARKERS = [gene for genes in BCELL_MARKERS.values() for gene in genes]

print(f"✓ Configuration complete")
print(f"  Input: {INPUT_FILE}")
print(f"  Output: {OUTPUT_DIR}")


# %% [markdown]
# ## 1. Data Loading and Structure Validation

# %%
# Load data
print("Loading data...")
adata = sc.read_h5ad(INPUT_FILE)

print(f"\n{'='*80}")
print("DATA STRUCTURE OVERVIEW")
print(f"{'='*80}")
print(f"Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"\nLayers available: {list(adata.layers.keys())}")
print(f"Has .raw: {adata.raw is not None}")
if adata.raw is not None:
    print(f"  .raw shape: {adata.raw.shape[0]:,} cells × {adata.raw.shape[1]:,} genes")

print(f"\n.obsm keys (embeddings): {list(adata.obsm.keys())}")
print(f".obsp keys (graphs): {list(adata.obsp.keys())}")

# %%
# Check metadata columns
print(f"\n{'='*80}")
print("METADATA COLUMNS")
print(f"{'='*80}")
print(f"Total columns: {len(adata.obs.columns)}")
print("\nKey columns:")
for col in adata.obs.columns[:20]:  # Show first 20
    unique_vals = adata.obs[col].nunique()
    print(f"  - {col}: {unique_vals} unique values")

# %%
# Validate expected columns
print(f"\n{'='*80}")
print("COLUMN VALIDATION")
print(f"{'='*80}")

expected_cols = {
    'Batch/Sample': BATCH_KEY,
    'Cell Type': CELLTYPE_KEY,
    'Disease Status': DISEASE_KEY
}

for name, col in expected_cols.items():
    if col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        print(f"✓ {name} ({col}): {n_unique} categories")
        if n_unique <= 20:
            print(f"  Values: {adata.obs[col].unique().tolist()}")
    else:
        print(f"⚠ {name} ({col}): NOT FOUND")
        print(f"  Available columns with similar names:")
        similar = [c for c in adata.obs.columns if any(keyword in c.lower() 
                  for keyword in col.lower().split('_'))]
        print(f"    {similar}")


# %% [markdown]
# ## 2. Quality Metrics Overview

# %%
# Basic QC metrics
print(f"\n{'='*80}")
print("QUALITY CONTROL METRICS")
print(f"{'='*80}")

qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
available_qc = [col for col in qc_metrics if col in adata.obs.columns]

if available_qc:
    for metric in available_qc:
        values = adata.obs[metric]
        print(f"\n{metric}:")
        print(f"  Mean: {values.mean():.2f}")
        print(f"  Median: {values.median():.2f}")
        print(f"  Range: [{values.min():.2f}, {values.max():.2f}]")
else:
    print("⚠ Standard QC metrics not found. Data might be already filtered.")
    print(f"  Available obs columns: {adata.obs.columns.tolist()}")

# %%
# QC visualization (if metrics available)
if available_qc:
    fig, axes = plt.subplots(1, len(available_qc), figsize=(5*len(available_qc), 4))
    if len(available_qc) == 1:
        axes = [axes]
    
    for ax, metric in zip(axes, available_qc):
        sc.pl.violin(adata, metric, ax=ax, show=False)
        ax.set_title(f'{metric} Distribution', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/01_qc_metrics.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/01_qc_metrics.png")
else:
    print("⚠ Skipping QC visualization (metrics not available)")


# %% [markdown]
# ## 3. UMAP Visualization - Batch and Cell Types

# %%
# Check if UMAP exists
if 'X_umap' not in adata.obsm.keys():
    print("⚠ UMAP not found. Computing UMAP...")
    if 'X_pca' not in adata.obsm.keys():
        print("  Computing PCA first...")
        sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    print("✓ UMAP computed")
else:
    print("✓ UMAP coordinates found")

# %%
# UMAP by batch
if BATCH_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, 
               title='B Cells - Batch Distribution', legend_loc='on data')
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/02_umap_batch.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/02_umap_batch.png")
else:
    print(f"⚠ Batch column '{BATCH_KEY}' not found. Skipping batch UMAP.")

# %%
# UMAP by cell type
if CELLTYPE_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.umap(adata, color=CELLTYPE_KEY, ax=ax, show=False,
               title='B Cells - Cell Type Annotation', legend_fontsize=10)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/03_umap_celltype.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/03_umap_celltype.png")
    
    # Print cell type distribution
    print(f"\n{'='*60}")
    print("CELL TYPE DISTRIBUTION")
    print(f"{'='*60}")
    celltype_counts = adata.obs[CELLTYPE_KEY].value_counts()
    for ct, count in celltype_counts.items():
        pct = 100 * count / len(adata)
        print(f"{ct:30s}: {count:6,} cells ({pct:5.2f}%)")
else:
    print(f"⚠ Cell type column '{CELLTYPE_KEY}' not found. Skipping cell type UMAP.")

# %%
# UMAP by disease status (if available)
if DISEASE_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=DISEASE_KEY, ax=ax, show=False,
               title='B Cells - Disease Status', palette='Set2')
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/04_umap_disease.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/04_umap_disease.png")
else:
    print(f"⚠ Disease column '{DISEASE_KEY}' not found. Skipping disease UMAP.")


# %% [markdown]
# ## 4. B Cell Marker Gene Expression

# %%
# Check marker availability
use_raw = adata.raw is not None
gene_names = adata.raw.var_names if use_raw else adata.var_names

available_markers = {}
for category, markers in BCELL_MARKERS.items():
    available = [m for m in markers if m in gene_names]
    if available:
        available_markers[category] = available

print(f"\n{'='*80}")
print("MARKER GENE AVAILABILITY")
print(f"{'='*80}")
print(f"Searching in: {'.raw' if use_raw else '.X'} ({len(gene_names):,} genes)")
for category, markers in BCELL_MARKERS.items():
    available = available_markers.get(category, [])
    missing = set(markers) - set(available)
    print(f"\n{category}:")
    print(f"  Available ({len(available)}/{len(markers)}): {', '.join(available)}")
    if missing:
        print(f"  Missing: {', '.join(missing)}")

# %%
# UMAP colored by key markers (plot available markers)
plot_markers = []
for markers in available_markers.values():
    plot_markers.extend(markers[:2])  # Take first 2 from each category

if plot_markers:
    n_markers = len(plot_markers)
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    axes = axes.flatten() if n_rows * n_cols > 1 else [axes]
    
    for i, marker in enumerate(plot_markers):
        sc.pl.umap(adata, color=marker, ax=axes[i], show=False, 
                   use_raw=use_raw, cmap='viridis', title=marker, 
                   frameon=False, vmin=0)
    
    # Hide unused axes
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/05_umap_markers.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/05_umap_markers.png")
else:
    print("⚠ No B cell markers found in dataset")

# %%
# Dotplot for all available markers (by cell type)
if available_markers and CELLTYPE_KEY in adata.obs.columns:
    all_available = [m for markers in available_markers.values() for m in markers]
    
    # Check if we have enough cells per cell type
    celltype_counts = adata.obs[CELLTYPE_KEY].value_counts()
    valid_celltypes = celltype_counts[celltype_counts >= 10].index.tolist()
    
    if len(valid_celltypes) > 0:
        adata_subset = adata[adata.obs[CELLTYPE_KEY].isin(valid_celltypes)].copy()
        
        fig, ax = plt.subplots(figsize=(12, max(6, len(valid_celltypes)*0.4)))
        sc.pl.dotplot(adata_subset, all_available, groupby=CELLTYPE_KEY, 
                      ax=ax, show=False, use_raw=use_raw, standard_scale='var',
                      title='B Cell Marker Expression Across Cell Types')
        plt.tight_layout()
        plt.savefig(f'{FIG_DIR}/06_dotplot_markers.png', dpi=300, bbox_inches='tight')
        plt.show()
        print(f"✓ Saved: {FIG_DIR}/06_dotplot_markers.png")
    else:
        print("⚠ Not enough cells per cell type for dotplot (need ≥10 cells)")
else:
    print("⚠ Skipping dotplot (no markers or cell types)")


# %% [markdown]
# ## 5. Cell Type Composition Analysis

# %%
# Cell type proportion by batch/sample
if CELLTYPE_KEY in adata.obs.columns and BATCH_KEY in adata.obs.columns:
    # Calculate proportions
    prop_df = adata.obs.groupby([BATCH_KEY, CELLTYPE_KEY]).size().unstack(fill_value=0)
    prop_df = prop_df.div(prop_df.sum(axis=1), axis=0) * 100
    
    # Heatmap
    fig, ax = plt.subplots(figsize=(12, max(6, len(prop_df)*0.3)))
    sns.heatmap(prop_df, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Percentage (%)'}, ax=ax)
    ax.set_title('B Cell Composition Across Samples (%)', fontsize=14)
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel(BATCH_KEY, fontsize=12)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/07_composition_heatmap.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/07_composition_heatmap.png")
    
    # Stacked bar plot
    fig, ax = plt.subplots(figsize=(12, 6))
    prop_df.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_title('B Cell Composition Across Samples', fontsize=14)
    ax.set_xlabel(BATCH_KEY, fontsize=12)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/08_composition_barplot.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/08_composition_barplot.png")

# %%
# Cell type proportion by disease status (if available)
if CELLTYPE_KEY in adata.obs.columns and DISEASE_KEY in adata.obs.columns:
    # Calculate proportions
    prop_disease = adata.obs.groupby([DISEASE_KEY, CELLTYPE_KEY]).size().unstack(fill_value=0)
    prop_disease = prop_disease.div(prop_disease.sum(axis=1), axis=0) * 100
    
    # Bar plot comparison
    fig, ax = plt.subplots(figsize=(12, 6))
    prop_disease.T.plot(kind='bar', ax=ax, width=0.8)
    ax.set_title('B Cell Composition by Disease Status', fontsize=14)
    ax.set_xlabel('Cell Type', fontsize=12)
    ax.set_ylabel('Percentage (%)', fontsize=12)
    ax.legend(title='Disease Status', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/09_composition_disease.png', dpi=300, bbox_inches='tight')
    plt.show()
    print(f"✓ Saved: {FIG_DIR}/09_composition_disease.png")
    
    # Print summary
    print(f"\n{'='*80}")
    print("CELL TYPE PROPORTION BY DISEASE STATUS")
    print(f"{'='*80}")
    print(prop_disease.to_string())


# %% [markdown]
# ## 6. Summary Statistics

# %%
print(f"\n{'='*80}")
print("ANALYSIS SUMMARY")
print(f"{'='*80}")
print(f"Total cells analyzed: {adata.shape[0]:,}")
print(f"Total genes: {adata.shape[1]:,}")
if adata.raw is not None:
    print(f"Full gene matrix: {adata.raw.shape[1]:,} genes")

if BATCH_KEY in adata.obs.columns:
    print(f"\nNumber of batches/samples: {adata.obs[BATCH_KEY].nunique()}")
    
if CELLTYPE_KEY in adata.obs.columns:
    print(f"\nNumber of cell types: {adata.obs[CELLTYPE_KEY].nunique()}")
    print("\nTop 5 cell types:")
    for i, (ct, count) in enumerate(adata.obs[CELLTYPE_KEY].value_counts().head(5).items(), 1):
        pct = 100 * count / len(adata)
        print(f"  {i}. {ct}: {count:,} cells ({pct:.2f}%)")

if DISEASE_KEY in adata.obs.columns:
    print(f"\nDisease status distribution:")
    for status, count in adata.obs[DISEASE_KEY].value_counts().items():
        pct = 100 * count / len(adata)
        print(f"  - {status}: {count:,} cells ({pct:.2f}%)")

print(f"\n{'='*80}")
print(f"All figures saved to: {FIG_DIR}")
print(f"{'='*80}")


# %% [markdown]
# ## 7. Optional: Marker Gene Analysis

# %%
# This section can be expanded for specific analyses:
# - Differential expression between conditions
# - Pseudotime trajectory analysis (if Plasma cells present)
# - BCR clonotype analysis (if BCR data available)
# - Gene set enrichment analysis

# Example: Find top marker genes for each cell type
if CELLTYPE_KEY in adata.obs.columns:
    print("\n[OPTIONAL] Computing marker genes for each cell type...")
    print("This may take a few minutes for large datasets...")
    
    # Check if markers already computed
    if 'rank_genes_groups' not in adata.uns.keys():
        try:
            sc.tl.rank_genes_groups(
                adata, 
                groupby=CELLTYPE_KEY, 
                method='wilcoxon',
                use_raw=use_raw,
                n_genes=100
            )
            print("✓ Marker gene computation complete")
            
            # Save top markers to CSV
            marker_df = sc.get.rank_genes_groups_df(adata, group=None)
            marker_df.to_csv(f'{OUTPUT_DIR}/marker_genes_all_celltypes.csv', index=False)
            print(f"✓ Saved: {OUTPUT_DIR}/marker_genes_all_celltypes.csv")
            
            # Visualize top markers
            fig = sc.pl.rank_genes_groups_dotplot(
                adata, 
                n_genes=5, 
                use_raw=use_raw,
                show=False,
                return_fig=True
            )
            fig.savefig(f'{FIG_DIR}/10_top_markers_dotplot.png', dpi=300, bbox_inches='tight')
            plt.show()
            print(f"✓ Saved: {FIG_DIR}/10_top_markers_dotplot.png")
            
        except Exception as e:
            print(f"⚠ Marker gene computation failed: {e}")
    else:
        print("✓ Marker genes already computed (found in adata.uns)")


# %% [markdown]
# ## 8. Save Processed Data (Optional)

# %%
# If you made any modifications, save the updated AnnData object
SAVE_UPDATED = False  # Set to True if you want to save

if SAVE_UPDATED:
    output_file = f'{OUTPUT_DIR}/adata_bcell_analyzed.h5ad'
    print(f"\nSaving updated AnnData object to: {output_file}")
    adata.write_h5ad(output_file, compression='gzip', compression_opts=9)
    print("✓ Saved successfully")
else:
    print("\n⚠ Data not saved (set SAVE_UPDATED=True to save)")


# %%
print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"Output directory: {OUTPUT_DIR}")
print(f"Figures: {FIG_DIR}")
print("\nGenerated outputs:")
print("  - QC metrics visualization")
print("  - UMAP plots (batch, cell type, disease, markers)")
print("  - Marker gene dotplot")
print("  - Cell composition heatmap and barplots")
print("  - (Optional) Top marker genes per cell type")
print("\nNext steps:")
print("  1. Review UMAP plots for batch effects and cell type separation")
print("  2. Validate marker gene expression patterns")
print("  3. Perform differential expression if comparing conditions")
print("  4. Consider trajectory analysis for B cell maturation")
print("="*80)
