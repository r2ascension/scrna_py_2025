# %% [markdown]
# # T Cell Analysis: BBKNN Integration + Wilcoxon Marker Finding
# 
# **Purpose**: Re-integrate T cells by dataset using BBKNN, then identify marker genes
# 
# **Workflow**:
# 1. Load T cell data (post-cNMF)
# 2. Inspect dataset distribution and filter small datasets
# 3. Re-run preprocessing with batch-aware HVG
# 4. BBKNN batch correction by dataset
# 5. Clustering and visualization
# 6. Visualize T cell canonical markers
# 7. Wilcoxon differential expression
# 8. Export results
# 
# **Author**: r2end  
# **Date**: 2025-01-05

# %% [markdown]
# ## Configuration

# %%
# ============================================================================
# CONFIGURATION PARAMETERS
# ============================================================================

# File paths
INPUT_H5AD = "/home/h2048/data/py/1217/cnmf_batch_production_v1_1_1/T_cells/batch_aware/cnmf_analysis_k40_1/T_cells_with_cnmf_k40.h5ad"
OUTPUT_DIR = "/home/h2048/data/R/0107/t_bbknn_filtered/"

# Epithelial filtering
SCANVI_COLUMN = 'scanvi_predictions'  # Column name for scanvi predictions
# Keywords to identify epithelial cells (case-insensitive matching)
EPITHELIAL_KEYWORDS = ['epithelial', 'basal', 'goblet', 'ciliated', 'secretory', 'club','Epithelial cells']

# Dataset filtering
BATCH_KEY = 'dataset'           # Column name for dataset/batch
MIN_CELLS_PER_DATASET = 20      # Filter datasets with < N cells

# Preprocessing
N_TOP_GENES = 4000              # Number of highly variable genes
N_PCS = 50                      # Number of PCs

# BBKNN parameters
BBKNN_NEIGHBORS_WITHIN_BATCH = 5
BBKNN_N_PCS = 50
BBKNN_TRIM = 35

# UMAP parameters
UMAP_MIN_DIST = 0.4
UMAP_SPREAD = 1.0

# Visualization

# Clustering
LEIDEN_RESOLUTION = 4

# Differential expression
MIN_LOGFC = 0.25
MIN_PCT = 0.1
TOP_N_MARKERS = 50

# Visualization
FIGURE_DPI = 300
FIGURE_FORMAT = 'pdf'
UMAP_SIZE = 3

# Performance
N_JOBS = 8

print("✓ Configuration loaded")

# %% [markdown]
# ## T Cell Canonical Markers
# 
# Define key marker genes for T cell subtypes for validation

# %%
# ============================================================================
# T CELL MARKER GENES
# ============================================================================

# Core T cell markers
TCELL_CORE_MARKERS = {
    'Pan_T': ['CD3D', 'CD3E', 'CD3G'],
    'CD4_T': ['CD4', 'CD40LG'],
    'CD8_T': ['CD8A', 'CD8B'],
    'NK': ['GNLY', 'NKG7', 'KLRD1', 'FCGR3A'],
    'Naive': ['CCR7', 'TCF7', 'LEF1', 'SELL', 'IL7R'],
    'Memory': ['GZMK', 'CD69'],
    'Effector': ['GZMB', 'GZMH', 'PRF1', 'IFNG'],
    'Treg': ['FOXP3', 'IL2RA', 'IKZF2', 'TNFRSF4'],
    'Exhausted': ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT'],
    'Proliferating': ['MKI67', 'TOP2A', 'STMN1'],
    'Th2': ['GATA3', 'IL4', 'IL5', 'IL13'],
    'Th17': ['RORC', 'IL17A', 'IL23R'],
}

# Flat list for quick visualization
TCELL_KEY_MARKERS = [
    'CD3D', 'CD4', 'CD8A',           # Major lineages
    'GNLY', 'NKG7',                  # NK
    'CCR7', 'SELL', 'IL7R',          # Naive
    'GZMK', 'CD69',                  # Memory
    'GZMB', 'PRF1',                  # Effector
    'FOXP3', 'IL2RA',                # Treg
    'PDCD1', 'HAVCR2',               # Exhausted
    'MKI67'                          # Proliferating
]

print("T cell marker genes defined:")
for category, markers in TCELL_CORE_MARKERS.items():
    print(f"  {category}: {', '.join(markers)}")

# %% [markdown]
# ## Imports & Setup

# %%
# ============================================================================
# IMPORTS
# ============================================================================

import scanpy as sc
import scanpy.external as sce
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# Scanpy settings
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=FIGURE_DPI, facecolor='white', frameon=False)
sc.settings.n_jobs = N_JOBS

# Create output directories
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

print("="*80)
print("T Cell BBKNN + Marker Analysis")
print("="*80)
print(f"Output directory: {output_dir}")

# %% [markdown]
# ## Step 1: Load Data & Inspect

# %%
# ============================================================================
# STEP 1: DATA LOADING
# ============================================================================

print("\n" + "="*80)
print("STEP 1: DATA LOADING & INSPECTION")
print("="*80)

print(f"\nLoading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)

print(f"\nData dimensions:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")

# Check metadata
print(f"\nAvailable .obs columns:")
for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    print(f"  - {col}: {n_unique} unique values")

# Check if cNMF results exist
cnmf_columns = [col for col in adata.obs.columns if 'cnmf' in col.lower() or 'usage' in col.lower()]
if len(cnmf_columns) > 0:
    print(f"\n✓ cNMF results detected:")
    for col in cnmf_columns[:5]:  # Show first 5
        print(f"    {col}")
    print(f"  (and {len(cnmf_columns)-5} more...)" if len(cnmf_columns) > 5 else "")
else:
    print(f"\n⚠️  No cNMF results found")

# Check batch key
if BATCH_KEY not in adata.obs.columns:
    raise ValueError(f"Batch key '{BATCH_KEY}' not found!")

# Dataset distribution
print(f"\nDataset distribution ({BATCH_KEY}):")
dataset_counts = adata.obs[BATCH_KEY].value_counts().sort_values(ascending=False)
print(f"  Total datasets: {len(dataset_counts)}")
print(f"\n  Dataset sizes:")
for dataset, count in dataset_counts.items():
    print(f"    {dataset}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

# Identify small datasets
small_datasets = dataset_counts[dataset_counts < MIN_CELLS_PER_DATASET]
if len(small_datasets) > 0:
    print(f"\n  ⚠️  Small datasets (< {MIN_CELLS_PER_DATASET} cells): {len(small_datasets)}")
    for dataset, count in small_datasets.items():
        print(f"      {dataset}: {count} cells (will be removed)")

# %%
# ============================================================================
# STEP 2: FILTER EPITHELIAL CELLS
# ============================================================================

print("\n" + "="*80)
print("STEP 2: EPITHELIAL CONTAMINATION FILTERING")
print("="*80)

# Check if scanvi_predict column exists
if SCANVI_COLUMN not in adata.obs.columns:
    print(f"\n⚠️  Warning: '{SCANVI_COLUMN}' column not found!")
    print(f"Available columns: {', '.join(adata.obs.columns)}")
    print(f"\nSkipping epithelial filtering...")
else:
    # Show cell type distribution before filtering
    print(f"\nCell type distribution (before filtering):")
    celltypes_before = adata.obs[SCANVI_COLUMN].value_counts()
    for ct, count in celltypes_before.head(20).items():
        print(f"  {ct}: {count:,} cells")
    if len(celltypes_before) > 20:
        print(f"  ... and {len(celltypes_before)-20} more cell types")
    
    # Identify epithelial cells (case-insensitive)
    epithelial_mask = adata.obs[SCANVI_COLUMN].str.lower().apply(
        lambda x: any(keyword in str(x).lower() for keyword in EPITHELIAL_KEYWORDS)
    )
    
    n_epithelial = epithelial_mask.sum()
    pct_epithelial = 100 * n_epithelial / len(adata)
    
    print(f"\nEpithelial contamination detected:")
    print(f"  Keywords used: {', '.join(EPITHELIAL_KEYWORDS)}")
    print(f"  Epithelial cells: {n_epithelial:,} ({pct_epithelial:.2f}%)")
    
    if n_epithelial > 0:
        # Show which epithelial types were found
        epithelial_types = adata.obs.loc[epithelial_mask, SCANVI_COLUMN].value_counts()
        print(f"\nEpithelial subtypes to be removed:")
        for etype, count in epithelial_types.items():
            print(f"  - {etype}: {count:,} cells")
        
        # Filter out epithelial cells
        print(f"\nRemoving epithelial cells...")
        adata = adata[~epithelial_mask].copy()
        
        print(f"\n✓ Filtering complete")
        print(f"  Cells after filtering: {adata.n_obs:,}")
        print(f"  Cells removed: {n_epithelial:,} ({pct_epithelial:.2f}%)")
        
        # Show remaining cell type distribution
        print(f"\nRemaining cell type distribution:")
        celltypes_after = adata.obs[SCANVI_COLUMN].value_counts()
        for ct, count in celltypes_after.head(15).items():
            print(f"  {ct}: {count:,} cells")
        if len(celltypes_after) > 15:
            print(f"  ... and {len(celltypes_after)-15} more cell types")
    else:
        print(f"\n✓ No epithelial contamination detected")

# %% [markdown]
# ## Step 2: Filter Small Datasets

# %%
# ============================================================================
# STEP 2: FILTER SMALL DATASETS
# ============================================================================

print("\n" + "="*80)
print("STEP 2: FILTER SMALL DATASETS")
print("="*80)

n_cells_before = adata.n_obs

dataset_counts = adata.obs[BATCH_KEY].value_counts()
valid_datasets = dataset_counts[dataset_counts >= MIN_CELLS_PER_DATASET].index

adata = adata[adata.obs[BATCH_KEY].isin(valid_datasets)].copy()

n_cells_after = adata.n_obs
n_cells_removed = n_cells_before - n_cells_after

print(f"\nFiltering results:")
print(f"  Cells before: {n_cells_before:,}")
print(f"  Cells after: {n_cells_after:,}")
print(f"  Cells removed: {n_cells_removed:,} ({n_cells_removed/n_cells_before*100:.1f}%)")
print(f"  Remaining datasets: {adata.obs[BATCH_KEY].nunique()}")

# %% [markdown]
# ## Step 3: Data Preparation

# %%
# ============================================================================
# STEP 3: DATA STATE CHECK
# ============================================================================

print("\n" + "="*80)
print("STEP 3: DATA STATE VERIFICATION")
print("="*80)

print(f"\nChecking data structure...")

# Check .X
X_min, X_max = adata.X.min(), adata.X.max()
X_mean = adata.X.mean()
print(f"\n.X matrix:")
print(f"  Type: {type(adata.X)}")
print(f"  Range: [{X_min:.4f}, {X_max:.4f}]")
print(f"  Mean: {X_mean:.4f}")

# Check layers
if 'counts' in adata.layers:
    print(f"\n✓ .layers['counts'] exists")
    HAS_COUNTS = True
else:
    print(f"\n⚠️  .layers['counts'] not found")
    HAS_COUNTS = False

# Prepare counts if needed
if not HAS_COUNTS:
    if X_max > 100:
        print(f"  Saving .X to .layers['counts']")
        adata.layers['counts'] = adata.X.copy()
        HAS_COUNTS = True
    else:
        raise ValueError(".X appears normalized but no counts layer found!")

# %%
# ============================================================================
# STEP 7.5: FILTER UNWANTED GENE CATEGORIES
# ============================================================================
# Insert this cell AFTER Step 7 (Visualization) and BEFORE Step 8 (Wilcoxon DE)

print("\n" + "="*80)
print("STEP 7.5: FILTER UNWANTED GENE CATEGORIES")
print("="*80)

print("\nFiltering genes to improve marker quality...")
print("Will remove: MT, ribosomal, histone, pseudogenes, ENSG, unannotated transcripts")

# ===== Configuration =====
REMOVE_MT = True
REMOVE_RIBO = True
REMOVE_HISTONE = True
REMOVE_PSEUDOGENES = True
REMOVE_ENSG = True
REMOVE_UNANNOTATED = True

# ===== Get all gene names from adata.raw =====
if adata.raw is None:
    print("\n⚠️  Warning: adata.raw is None, will filter adata.var instead")
    all_genes = adata.var_names.tolist()
    use_raw = False
else:
    all_genes = adata.raw.var_names.tolist()
    use_raw = True

print(f"\nTotal genes before filtering: {len(all_genes):,}")

# ===== Initialize list of genes to remove =====
genes_to_remove = []

# ===== 1. Mitochondrial genes (MT-) =====
if REMOVE_MT:
    mt_genes = [g for g in all_genes if g.startswith('MT-')]
    genes_to_remove.extend(mt_genes)
    print(f"\n  1. Mitochondrial genes (MT-): {len(mt_genes)} genes")
    if len(mt_genes) > 0:
        print(f"     Examples: {', '.join(mt_genes[:5])}")

# ===== 2. Ribosomal genes (RPS, RPL, MRPS, MRPL) =====
if REMOVE_RIBO:
    import re
    ribo_pattern = re.compile(r'^(RPS|RPL|MRPS|MRPL)')
    ribo_genes = [g for g in all_genes if ribo_pattern.match(g)]
    genes_to_remove.extend(ribo_genes)
    print(f"  2. Ribosomal genes (RPS/RPL/MRPS/MRPL): {len(ribo_genes)} genes")
    if len(ribo_genes) > 0:
        print(f"     Examples: {', '.join(ribo_genes[:5])}")

# ===== 3. Histone genes (H1, H2A, H2B, H3, H4, HIST) =====
if REMOVE_HISTONE:
    histone_pattern = re.compile(r'^(H1|H2A|H2B|H3|H4|HIST)')
    histone_genes = [g for g in all_genes if histone_pattern.match(g)]
    genes_to_remove.extend(histone_genes)
    print(f"  3. Histone genes (H1/H2A/H2B/H3/H4/HIST): {len(histone_genes)} genes")
    if len(histone_genes) > 0:
        print(f"     Examples: {', '.join(histone_genes[:5])}")

# ===== 4. Pseudogenes (e.g., RPS29P1, RPL10P9) =====
if REMOVE_PSEUDOGENES:
    # Pattern: RPS/RPL/MRPS/MRPL + digits + P + digits
    pseudo_pattern = re.compile(r'^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$')
    pseudo_genes = [g for g in all_genes if pseudo_pattern.match(g)]
    genes_to_remove.extend(pseudo_genes)
    print(f"  4. Pseudogenes (*P*): {len(pseudo_genes)} genes")
    if len(pseudo_genes) > 0:
        print(f"     Examples: {', '.join(pseudo_genes[:5])}")

# ===== 5. ENSG unannotated genes =====
if REMOVE_ENSG:
    ensg_pattern = re.compile(r'^ENSG[0-9]+')
    ensg_genes = [g for g in all_genes if ensg_pattern.match(g)]
    genes_to_remove.extend(ensg_genes)
    print(f"  5. ENSG unannotated genes: {len(ensg_genes)} genes")
    if len(ensg_genes) > 0:
        print(f"     Examples: {', '.join(ensg_genes[:5])}")

# ===== 6. Unannotated transcripts =====
if REMOVE_UNANNOTATED:
    # Pattern includes:
    # - AC/AL/AP/BX/Z followed by digits and dot
    # - RP followed by digits and dash
    # - CTD-/CTB-/CTC-
    # - LINC followed by digits
    # - Ending with -AS + digits (antisense)
    # - Ending with -OT + digits (overlapping transcript)
    # - Starting with LOC + digits
    unannotated_pattern = re.compile(
        r'^(AC|AL|AP|BX|Z)[0-9]+\.|'
        r'^RP[0-9]+-|'
        r'^CTD-|^CTB-|^CTC-|'
        r'^LINC[0-9]+|'
        r'-AS[0-9]+$|'
        r'-OT[0-9]+$|'
        r'^LOC[0-9]+'
    )
    unannotated_genes = [g for g in all_genes if unannotated_pattern.search(g)]
    genes_to_remove.extend(unannotated_genes)
    print(f"  6. Unannotated transcripts (AC/AL/RP/CTD/LINC/LOC/etc): {len(unannotated_genes)} genes")
    if len(unannotated_genes) > 0:
        print(f"     Examples: {', '.join(unannotated_genes[:5])}")

# ===== Remove duplicates =====
genes_to_remove = list(set(genes_to_remove))
print(f"\n{'─'*80}")
print(f"Total unique genes to remove: {len(genes_to_remove):,}")

# ===== Determine genes to keep =====
genes_to_keep = [g for g in all_genes if g not in genes_to_remove]
print(f"Genes to keep: {len(genes_to_keep):,}")
print(f"Percentage retained: {len(genes_to_keep)/len(all_genes)*100:.1f}%")

# ===== Filter adata.raw =====
if use_raw:
    print(f"\nFiltering adata.raw...")
    # IMPORTANT: Must convert adata.raw to AnnData, filter, then reassign
    # Direct slicing of adata.raw does not work
    
    # Step 1: Convert raw to full AnnData object
    raw_adata = adata.raw.to_adata()
    
    # Step 2: Filter genes
    raw_adata_filtered = raw_adata[:, genes_to_keep].copy()
    
    # Step 3: Reassign to adata.raw
    adata.raw = raw_adata_filtered
    
    print(f"  ✓ adata.raw filtered: {adata.raw.n_vars:,} genes remain")
    
    # Clean up temporary objects
    del raw_adata, raw_adata_filtered
    
else:
    print(f"\nFiltering adata.var...")
    adata = adata[:, genes_to_keep].copy()
    print(f"  ✓ adata filtered: {adata.n_vars:,} genes remain")

# ===== Summary of removed gene categories =====
print(f"\n{'─'*80}")
print("Summary of removed gene categories:")
if REMOVE_MT:
    print(f"  ✓ Mitochondrial genes (MT-)")
if REMOVE_RIBO:
    print(f"  ✓ Ribosomal genes (RPS/RPL/MRPS/MRPL)")
if REMOVE_HISTONE:
    print(f"  ✓ Histone genes (H1/H2A/H2B/H3/H4/HIST)")
if REMOVE_PSEUDOGENES:
    print(f"  ✓ Pseudogenes (*P*)")
if REMOVE_ENSG:
    print(f"  ✓ ENSG unannotated genes")
if REMOVE_UNANNOTATED:
    print(f"  ✓ Unannotated transcripts (AC/AL/RP/CTD/LINC/LOC/etc)")

print("\n✓ Gene filtering complete")
print("  → Subsequent differential expression will use filtered gene set")

# %% [markdown]
# ## Step 4: Preprocessing

# %%
# ============================================================================
# STEP 4: PREPROCESSING
# ============================================================================

print("\n" + "="*80)
print("STEP 4: PREPROCESSING (HVG + PCA)")
print("="*80)

# Normalize and log
print(f"\nNormalization...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print(f"  ✓ Log-normalized")

# HVG selection
print(f"\nHighly variable genes...")
print(f"  Target: {N_TOP_GENES} genes")
print(f"  Batch key: {BATCH_KEY}")

try:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        batch_key=BATCH_KEY,
        flavor='seurat_v3',
        subset=False
    )
    print(f"  ✓ Batch-aware HVG completed")
except Exception as e:
    print(f"  ⚠️  Batch-aware failed, using non-batch method")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor='seurat_v3',
        subset=False
    )
    print(f"  ✓ Non-batch-aware HVG completed")

n_hvg = adata.var['highly_variable'].sum()
print(f"  Selected HVGs: {n_hvg}")

# Save full data to .raw
print(f"\nPreserving full gene data...")
adata.raw = adata.copy()
print(f"  ✓ adata.raw saved (all {adata.n_vars} genes)")

# Subset to HVGs
adata = adata[:, adata.var['highly_variable']].copy()
print(f"  ✓ Subset to {adata.n_vars} HVGs")

# PCA
print(f"\nPCA...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=N_PCS, svd_solver='arpack')
print(f"  ✓ PCA completed ({N_PCS} components)")

# %% [markdown]
# ## Step 5: BBKNN Integration

# %%
# ============================================================================
# STEP 5: BBKNN BATCH CORRECTION
# ============================================================================

print("\n" + "="*80)
print("STEP 5: BBKNN BATCH CORRECTION")
print("="*80)

print(f"\nBBKNN parameters:")
print(f"  Batch key: {BATCH_KEY}")
print(f"  Neighbors within batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
print(f"  Number of PCs: {BBKNN_N_PCS}")
print(f"  Trim: {BBKNN_TRIM}")  # ← 添加这一行显示trim参数
print(f"  Number of datasets: {adata.obs[BATCH_KEY].nunique()}")

print(f"\nRunning BBKNN...")
sce.pp.bbknn(
    adata,
    batch_key=BATCH_KEY,
    neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
    n_pcs=BBKNN_N_PCS,
    trim=BBKNN_TRIM  # ← 修改这一行，从trim=None改为trim=BBKNN_TRIM
)
print(f"  ✓ BBKNN completed")

# UMAP
print(f"\nComputing UMAP...")
print(f"  Parameters:")
print(f"    - min_dist: {UMAP_MIN_DIST}")
print(f"    - spread: {UMAP_SPREAD}")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD)
print(f"  ✓ UMAP completed")

# %% [markdown]
# ## Step 6: Clustering

# %%
# ============================================================================
# STEP 6: LEIDEN CLUSTERING
# ============================================================================

print("\n" + "="*80)
print("STEP 6: LEIDEN CLUSTERING")
print("="*80)

print(f"\nClustering parameters:")
print(f"  Method: Leiden")
print(f"  Resolution: {LEIDEN_RESOLUTION}")

sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added='leiden')

n_clusters = adata.obs['leiden'].nunique()
print(f"\n✓ Clustering completed: {n_clusters} clusters")

print(f"\nCluster sizes:")
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
for cluster, count in cluster_counts.items():
    print(f"  Cluster {cluster}: {count:,} cells ({count/adata.n_obs*100:.1f}%)")

# %% [markdown]
# ## Step 7: Basic Visualization

# %%
# ============================================================================
# STEP 7: BASIC VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print("STEP 7: VISUALIZATION")
print("="*80)

# UMAP by cluster
print(f"\nGenerating UMAP plots...")

fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata,
    color='leiden',
    ax=ax,
    show=False,
    legend_loc='on data',
    legend_fontsize=8,
    size=UMAP_SIZE,
    title='T Cell Clusters (BBKNN-integrated)'
)
plt.tight_layout()
plt.savefig(fig_dir / f'01_umap_clusters.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
plt.show()
print(f"  ✓ Saved: 01_umap_clusters.{FIGURE_FORMAT}")

# UMAP by dataset
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata,
    color=BATCH_KEY,
    ax=ax,
    show=False,
    size=UMAP_SIZE,
    title=f'Batch Mixing ({BATCH_KEY})'
)
plt.tight_layout()
plt.savefig(fig_dir / f'02_umap_batch.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
plt.show()
print(f"  ✓ Saved: 02_umap_batch.{FIGURE_FORMAT}")


# UMAP by dataset
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(
    adata,
    color=SCANVI_COLUMN,
    ax=ax,
    show=False,
    size=UMAP_SIZE,
    title=f'ScanVI Predictions ({SCANVI_COLUMN})'
)
plt.tight_layout()
plt.savefig(fig_dir / f'03_umap_scanvi.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
plt.show()
print(f"  ✓ Saved: 03_umap_scanvi.{FIGURE_FORMAT}")

# Combined view
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, legend_loc='on data', size=UMAP_SIZE)
sc.pl.umap(adata, color=BATCH_KEY, ax=axes[1], show=False, size=UMAP_SIZE)
axes[0].set_title('Clusters', fontsize=14, weight='bold')
axes[1].set_title('Dataset Batch', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(fig_dir / f'03_umap_combined.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
plt.show()
print(f"  ✓ Saved: 03_umap_combined.{FIGURE_FORMAT}")

# %% [markdown]
# ## Step 8: T Cell Marker Visualization

# %%
# ============================================================================
# STEP 8: T CELL MARKER VISUALIZATION
# ============================================================================

print("\n" + "="*80)
print("STEP 8: T CELL CANONICAL MARKER VISUALIZATION")
print("="*80)

# Check which markers are available
available_markers = [m for m in TCELL_KEY_MARKERS if m in adata.raw.var_names]
missing_markers = [m for m in TCELL_KEY_MARKERS if m not in adata.raw.var_names]

print(f"\nMarker availability:")
print(f"  Available: {len(available_markers)}/{len(TCELL_KEY_MARKERS)}")
if len(missing_markers) > 0:
    print(f"  Missing: {', '.join(missing_markers)}")

if len(available_markers) > 0:
    # UMAP grid
    print(f"\nGenerating marker UMAP grid...")
    n_cols = 4
    n_rows = int(np.ceil(len(available_markers) / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes]
    
    for idx, marker in enumerate(available_markers):
        sc.pl.umap(
            adata,
            color=marker,
            ax=axes[idx],
            show=False,
            use_raw=True,
            vmax='p99',
            frameon=False,
            size=UMAP_SIZE,
            title=marker
        )
    
    # Hide extra subplots
    for idx in range(len(available_markers), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(fig_dir / f'04_tcell_markers_umap.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
    plt.show()
    print(f"  ✓ Saved: 04_tcell_markers_umap.{FIGURE_FORMAT}")
    
    # Dotplot
    print(f"\nGenerating marker dotplot...")
    try:
        fig = sc.pl.dotplot(
            adata,
            var_names=available_markers,
            groupby='leiden',
            use_raw=True,
            show=False,
            figsize=(14, 6),
            standard_scale='var'
        )
        plt.savefig(fig_dir / f'05_tcell_markers_dotplot.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved: 05_tcell_markers_dotplot.{FIGURE_FORMAT}")
    except Exception as e:
        print(f"  ⚠️  Dotplot failed: {e}")

else:
    print(f"\n⚠️  No canonical T cell markers found in dataset")

# %% [markdown]
# ## Step 9: Wilcoxon Differential Expression

# %%
# ============================================================================
# STEP 9: WILCOXON DIFFERENTIAL EXPRESSION
# ============================================================================

print("\n" + "="*80)
print("STEP 9: DIFFERENTIAL EXPRESSION ANALYSIS")
print("="*80)

print(f"\nRunning Wilcoxon rank-sum test...")
print(f"  Method: One-vs-rest")
print(f"  Using: adata.raw (all genes)")
print(f"  Top N genes: {TOP_N_MARKERS}")

sc.tl.rank_genes_groups(
    adata,
    groupby='leiden',
    method='wilcoxon',
    use_raw=True,
    n_genes=TOP_N_MARKERS,
    key_added='rank_genes_wilcox'
)
print(f"  ✓ Differential expression completed")

# %% [markdown]
# ## Step 10: Extract & Filter Markers

# %%
# ============================================================================
# STEP 10: EXTRACT & FILTER MARKERS
# ============================================================================

print("\n" + "="*80)
print("STEP 10: EXTRACTING & FILTERING MARKERS")
print("="*80)

print(f"\nFiltering criteria:")
print(f"  - log2FC > {MIN_LOGFC}")
print(f"  - Adjusted p-value < 0.05")
print(f"  - Expression % > {MIN_PCT*100}%")

all_markers = []

for cluster in adata.obs['leiden'].cat.categories:
    cluster_result = sc.get.rank_genes_groups_df(
        adata,
        group=cluster,
        key='rank_genes_wilcox'
    )
    
    cluster_result_filtered = cluster_result[
        (cluster_result['logfoldchanges'] > MIN_LOGFC) &
        (cluster_result['pvals_adj'] < 0.05)
    ].copy()
    
    if len(cluster_result_filtered) == 0:
        print(f"  Cluster {cluster}: No significant markers")
        continue
    
    # Calculate expression percentage
    cluster_cells = adata.raw.X[adata.obs['leiden'] == cluster]
    other_cells = adata.raw.X[adata.obs['leiden'] != cluster]
    
    pct_in = []
    pct_out = []
    
    for gene in cluster_result_filtered['names']:
        gene_idx = adata.raw.var_names.get_loc(gene)
        
        if sparse.issparse(cluster_cells):
            expr_in = cluster_cells[:, gene_idx].toarray().flatten()
            expr_out = other_cells[:, gene_idx].toarray().flatten()
        else:
            expr_in = cluster_cells[:, gene_idx].flatten()
            expr_out = other_cells[:, gene_idx].flatten()
        
        pct_in.append(np.sum(expr_in > 0) / len(expr_in))
        pct_out.append(np.sum(expr_out > 0) / len(expr_out))
    
    cluster_result_filtered['cluster'] = cluster
    cluster_result_filtered['pct_in_cluster'] = pct_in
    cluster_result_filtered['pct_out_cluster'] = pct_out
    
    cluster_result_filtered = cluster_result_filtered[
        cluster_result_filtered['pct_in_cluster'] > MIN_PCT
    ]
    
    n_markers = len(cluster_result_filtered)
    print(f"  Cluster {cluster}: {n_markers} markers")
    
    all_markers.append(cluster_result_filtered)

# Combine results
if len(all_markers) > 0:
    markers_df = pd.concat(all_markers, ignore_index=True)
    markers_df = markers_df.sort_values(['cluster', 'pvals_adj'])
    
    markers_df.to_csv(output_dir / "tcell_markers_wilcox_all.csv", index=False)
    print(f"\n✓ Total markers: {len(markers_df):,}")
    print(f"✓ Saved: tcell_markers_wilcox_all.csv")
    
    top_markers = markers_df.groupby('cluster').head(10)
    top_markers.to_csv(output_dir / "tcell_markers_wilcox_top10.csv", index=False)
    print(f"✓ Saved: tcell_markers_wilcox_top10.csv")
else:
    print(f"\n⚠️  No markers passed filtering!")
    markers_df = pd.DataFrame()

# %% [markdown]
# ## Step 11: Marker Visualization

# %%
# ============================================================================
# STEP 11: MARKER VISUALIZATION
# ============================================================================

if len(markers_df) > 0:
    print("\n" + "="*80)
    print("STEP 11: MARKER VISUALIZATION")
    print("="*80)
    
    # Heatmap
    print(f"\nGenerating heatmap (top 5 per cluster)...")
    try:
        fig = sc.pl.rank_genes_groups_heatmap(
            adata,
            n_genes=5,
            key='rank_genes_wilcox',
            use_raw=True,
            show=False,
            cmap='RdBu_r',
            figsize=(12, 10),
            vmin=-3,
            vmax=3,
            dendrogram=False
        )
        plt.savefig(fig_dir / f'06_marker_heatmap.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved: 06_marker_heatmap.{FIGURE_FORMAT}")
    except Exception as e:
        print(f"  ⚠️  Heatmap failed: {e}")
    
    # Dotplot
    print(f"\nGenerating dotplot (top 3 per cluster)...")
    try:
        fig = sc.pl.rank_genes_groups_dotplot(
            adata,
            n_genes=3,
            key='rank_genes_wilcox',
            use_raw=True,
            show=False,
            figsize=(14, 6)
        )
        plt.savefig(fig_dir / f'07_marker_dotplot.{FIGURE_FORMAT}', dpi=FIGURE_DPI, bbox_inches='tight')
        plt.show()
        print(f"  ✓ Saved: 07_marker_dotplot.{FIGURE_FORMAT}")
    except Exception as e:
        print(f"  ⚠️  Dotplot failed: {e}")

print("\n✓ Visualization complete")

# %% [markdown]
# ## Step 12: Save Results

# %%
# ============================================================================
# STEP 12: SAVE PROCESSED DATA
# ============================================================================

print("\n" + "="*80)
print("STEP 12: SAVING PROCESSED DATA")
print("="*80)

output_h5ad = output_dir / "tcell_bbknn_processed.h5ad"

print(f"\nSaving AnnData object...")
print(f"  Output: {output_h5ad}")
print(f"  Contents:")
print(f"    - Cells: {adata.n_obs:,}")
print(f"    - HVGs: {adata.n_vars:,}")
print(f"    - Full genes in .raw: {adata.raw.n_vars:,}")
print(f"    - BBKNN neighbors: ✓")
print(f"    - UMAP: ✓")
print(f"    - Leiden clusters: ✓")
print(f"    - Wilcoxon DE: ✓")

adata.write_h5ad(output_h5ad, compression='gzip')
print(f"\n✓ Data saved successfully")

# %% [markdown]
# ## Summary

# %%
# ============================================================================
# ANALYSIS SUMMARY
# ============================================================================

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)

print(f"\nOutput directory: {output_dir}")

print(f"\nGenerated files:")
print(f"  ├── tcell_bbknn_processed.h5ad")
print(f"  ├── tcell_markers_wilcox_all.csv")
print(f"  ├── tcell_markers_wilcox_top10.csv")
print(f"  └── figures/")
print(f"      ├── 01_umap_clusters.{FIGURE_FORMAT}")
print(f"      ├── 02_umap_batch.{FIGURE_FORMAT}")
print(f"      ├── 03_umap_combined.{FIGURE_FORMAT}")
print(f"      ├── 04_tcell_markers_umap.{FIGURE_FORMAT}")
print(f"      ├── 05_tcell_markers_dotplot.{FIGURE_FORMAT}")
print(f"      ├── 06_marker_heatmap.{FIGURE_FORMAT}")
print(f"      └── 07_marker_dotplot.{FIGURE_FORMAT}")

if len(markers_df) > 0:
    print(f"\nMarker summary by cluster:")
    for cluster in sorted(markers_df['cluster'].unique()):
        cluster_markers = markers_df[markers_df['cluster'] == cluster]
        top3 = cluster_markers.head(3)['names'].tolist()
        print(f"  Cluster {cluster}: {len(cluster_markers)} markers")
        print(f"    Top 3: {', '.join(top3)}")

print("\n" + "="*80)
print("Next steps:")
print("  1. Review canonical T cell markers (CD3D, CD4, CD8A, etc.)")
print("  2. Identify potential T cell subtypes:")
print("     - Naive T (CCR7+, SELL+, IL7R+)")
print("     - Memory T (GZMK+, CD69+)")
print("     - Effector T (GZMB+, PRF1+)")
print("     - Treg (FOXP3+, IL2RA+)")
print("     - NK (GNLY+, NKG7+)")
print("  3. Consider trajectory analysis for differentiation dynamics")
print("  4. Compare CD4+ vs CD8+ subsets if both present")
print("="*80)


