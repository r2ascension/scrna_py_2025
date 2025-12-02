# %% [markdown]
# # T Cell Analysis Pipeline with starCAT Integration
# 
# **Author:** Clinical-Bioinformatics Team  
# **Date:** 2025-12-01  
# **Version:** v1.0
# 
# ## Pipeline Overview
# 
# This notebook provides a comprehensive T cell analysis workflow:
# 
# 1. **Data Loading & QC** - Load and quality control
# 2. **Preprocessing** - Normalization and feature selection
# 3. **BBKNN Integration** - Batch effect correction
# 4. **Dimensionality Reduction** - PCA and UMAP
# 5. **Basic Clustering** - Multi-resolution Leiden clustering
# 6. **starCAT Analysis** - Reference-based high-resolution annotation
# 7. **Comprehensive Visualization** - Integrated results visualization
# 
# ## starCAT Method Explanation
# 
# **starCAT (Star Cell Annotation Tool)** is a reference-based cell annotation method:
# 
# ### Core Concept
# - Uses pre-trained **Gene Expression Programs (GEPs)** as references
# - Each GEP represents a specific cellular function/state (e.g., cell cycle, cytotoxicity, exhaustion)
# - Based on **consensus Non-negative Matrix Factorization (cNMF)**
# 
# ### Key Components
# 
# 1. **Reference Spectra**: Gene × Program matrix
#    - Contains weight of each gene in each GEP
#    - For T cells: TCAT.V1 reference with 52 GEPs
# 
# 2. **Usage Matrix**: Cell × Program matrix
#    - Normalized usage of each GEP in each cell (sums to 1 per cell)
#    - Captures cellular heterogeneity at high resolution
# 
# 3. **Scores**: Derived features
#    - **Continuous scores**: e.g., ASA (activation), Proliferation
#    - **Discrete scores**: e.g., ASA_binary, Multinomial_Label (CD4/CD8 subtypes)
# 
# ### Workflow Steps
# 
# ```
# Input: Raw counts (genes × cells)
#   ↓
# Extract local reference (if first run)
#   ↓
# Load Reference (TCAT.V1)
#   ↓
# Match genes with reference
#   ↓
# Compute Usage via NMF projection
#   ↓
# Calculate Scores from Usage
#   ↓
# Output: Usage matrix + Scores dataframe
# ```
# 
# **Note:** This pipeline uses a pre-downloaded local reference file at `/home/h2048/data/source/reference/TCAT.V1.tar.gz` to avoid repeated downloads. On first run, it will be extracted to the cache directory.
# 
# ### Advantages
# - **High resolution**: 52 GEPs capture fine-grained states
# - **Reference-based**: Consistent annotation across datasets
# - **Interpretable**: Each GEP has biological meaning
# - **Complementary**: Works with standard clustering methods

# %% [markdown]
# ---
# ## Configuration Section

# %%
# ==================== Import Libraries ====================
import sys
import os
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import datetime

# Single-cell analysis
import scanpy as sc
import bbknn

# starCAT
from starcat import starCAT

warnings.filterwarnings('ignore')

# Display versions
print(f"scanpy: {sc.__version__}")
print(f"Python: {sys.version}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")

# Set plot parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)
plt.rcParams['figure.figsize'] = (8, 6)

# %%
# ==================== Configuration Parameters ====================

# ========== Input/Output Configuration ==========
INPUT_H5AD_PATH = "/home/h2048/data/py/1128/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1202/bbknn_celltype_analysis/T/tcell_analysis_starcat"
CACHE_DIR = "./starcat_cache"  # Cache directory for starCAT reference

# ========== BBKNN Integration Configuration ==========
BATCH_KEY = "dataset"  # Column name for batch variable
BBKNN_NEIGHBORS_WITHIN_BATCH = 4  # Lower for T cells to enhance separation
BBKNN_N_PCS = 50
BBKNN_METRIC = "correlation"  # Suitable for T cells
BBKNN_TRIM = 35

# ========== Highly Variable Gene Selection ==========
USE_HVG = True
N_TOP_GENES = 4000
HVG_FLAVOR = "seurat_v3"

# ========== Gene Filtering Configuration ==========
MIN_CELLS_PER_GENE = 3

# ========== Normalization Configuration ==========
NORMALIZE_TOTAL = True
TARGET_SUM = 1e4
LOG_TRANSFORM = True
SCALE_DATA = True
MAX_VALUE = 10

# ========== PCA Configuration ==========
N_PCS = 50

# ========== UMAP Configuration ==========
RUN_UMAP = True
UMAP_MIN_DIST = 0.15  # Smaller value for T cells
UMAP_N_NEIGHBORS = 30

# ========== Clustering Configuration ==========
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4, 2.8]
DEFAULT_RESOLUTION = 2.0

# ========== starCAT Configuration ==========
STARCAT_REFERENCE = "TCAT.V1"  # T cell reference
RUN_STARCAT = True
LOCAL_REFERENCE_TAR = "/home/h2048/data/source/reference/TCAT.V1.tar.gz"  # Local reference file path

# ========== Key Marker Genes for T Cells ==========
MARKER_GENES = {
    "Pan_T": ['CD3D', 'CD3E', 'CD3G', 'CD2'],
    "CD4": ['CD4', 'IL7R'],
    "CD8": ['CD8A', 'CD8B'],
    "Naive": ['CCR7', 'SELL', 'LEF1', 'TCF7'],
    "Effector_Memory": ['GZMA', 'GZMB', 'GZMK', 'PRF1', 'NKG7'],
    "Treg": ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
    "Activation": ['IFNG', 'TNF', 'IL2', 'CD69'],
    "Exhaustion": ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'TOX'],
    "Resident_Memory": ['CD69', 'ITGAE', 'CXCR6', 'ZNF683'],
    "Proliferation": ['MKI67', 'TOP2A', 'PCNA']
}

# Flatten marker list
ALL_MARKERS = [gene for genes in MARKER_GENES.values() for gene in genes]
ALL_MARKERS = list(set(ALL_MARKERS))  # Remove duplicates

# ========== Visualization Configuration ==========
DPI = 300
FIGURE_FORMAT = "pdf"

VERBOSE = True

print("Configuration loaded successfully")
print(f"Input file: {INPUT_H5AD_PATH}")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Batch key: {BATCH_KEY}")

# %%
# ==================== Create Output Directories ====================
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

cache_dir = Path(CACHE_DIR)
cache_dir.mkdir(parents=True, exist_ok=True)

print(f"Output directory created: {output_dir}")
print(f"Figure directory: {fig_dir}")
print(f"Cache directory: {cache_dir}")

# %% [markdown]
# ---
# ## Step 1: Load Data

# %%
print("\n" + "="*70)
print("Step 1: Loading Data")
print("="*70)

adata = sc.read_h5ad(INPUT_H5AD_PATH)
print(f"\nLoaded data: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Display basic info
print(f"\nData structure:")
print(f"  obs columns: {len(adata.obs.columns)}")
print(f"  var columns: {len(adata.var.columns)}")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  uns keys: {list(adata.uns.keys())}")

# Check if raw counts are available
if adata.raw is not None:
    print(f"\nRaw data available: {adata.raw.shape[0]:,} cells × {adata.raw.shape[1]:,} genes")
else:
    print("\nNo raw data found, will use .X as raw counts")

# Display batch distribution
if BATCH_KEY in adata.obs.columns:
    print(f"\nBatch distribution (key: {BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        print(f"  {batch}: {count:,} cells ({pct:.1f}%)")
else:
    print(f"\nWarning: Batch key '{BATCH_KEY}' not found in adata.obs")

adata

# %% [markdown]
# ---
# ## Step 2: Check Marker Gene Availability

# %%
print("\n" + "="*70)
print("Step 2: Checking Marker Gene Availability")
print("="*70)

# Collect available gene names
available_genes = set(adata.var_names)
if adata.raw is not None:
    available_genes = available_genes.union(set(adata.raw.var_names))

# Check markers
available_markers = {}
missing_markers = {}

for category, genes in MARKER_GENES.items():
    available = [g for g in genes if g in available_genes]
    missing = [g for g in genes if g not in available_genes]
    available_markers[category] = available
    missing_markers[category] = missing
    
    print(f"\n{category}:")
    print(f"  Available: {len(available)}/{len(genes)} ({len(available)/len(genes)*100:.1f}%)")
    if missing:
        print(f"  Missing: {', '.join(missing)}")

# Summary
total_markers = len(ALL_MARKERS)
total_available = len([g for g in ALL_MARKERS if g in available_genes])
print(f"\nTotal marker availability: {total_available}/{total_markers} ({total_available/total_markers*100:.1f}%)")

# %% [markdown]
# ---
# ## Step 3: Preprocessing for BBKNN Integration

# %%
print("\n" + "="*70)
print("Step 3: Preprocessing for BBKNN Integration")
print("="*70)

# Store raw counts if not already stored
if adata.raw is None:
    print("\nStoring raw counts...")
    adata.raw = adata.copy()
    print("Raw counts stored")

# Create full dataset copy
adata_full = adata.copy()
print(f"\nFull dataset backup created: {adata_full.shape}")

# %%
# Gene filtering
print(f"\nGene filtering (min_cells={MIN_CELLS_PER_GENE})...")
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
print(f"After filtering: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# %%
# Normalization
if NORMALIZE_TOTAL:
    print(f"\nNormalizing to {TARGET_SUM:,.0f} total counts per cell...")
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    print("Normalization complete")

if LOG_TRANSFORM:
    print("\nApplying log1p transformation...")
    sc.pp.log1p(adata)
    print("Log transformation complete")

# %%
# Highly variable gene selection
if USE_HVG:
    print(f"\nSelecting top {N_TOP_GENES} highly variable genes (flavor={HVG_FLAVOR})...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor=HVG_FLAVOR,
        batch_key=BATCH_KEY if BATCH_KEY in adata.obs.columns else None
    )
    
    n_hvg = adata.var['highly_variable'].sum()
    print(f"Selected {n_hvg} highly variable genes")
    
    # Force include marker genes
    print("\nForce including marker genes in HVG set...")
    markers_to_add = [g for g in ALL_MARKERS if g in adata.var_names and not adata.var.loc[g, 'highly_variable']]
    if markers_to_add:
        adata.var.loc[markers_to_add, 'highly_variable'] = True
        print(f"Added {len(markers_to_add)} marker genes to HVG set")
    
    # Filter to HVG
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    print(f"\nFiltered to HVG: {adata_hvg.shape[0]:,} cells × {adata_hvg.shape[1]:,} genes")
else:
    print("\nSkipping HVG selection, using all genes")
    adata_hvg = adata.copy()

# %%
# Scale data
if SCALE_DATA:
    print(f"\nScaling data (max_value={MAX_VALUE})...")
    sc.pp.scale(adata_hvg, max_value=MAX_VALUE)
    print("Scaling complete")

print(f"\nPreprocessing completed")
print(f"Working dataset: {adata_hvg.shape[0]:,} cells × {adata_hvg.shape[1]:,} genes")

# %% [markdown]
# ---
# ## Step 4: PCA Dimensionality Reduction

# %%
print("\n" + "="*70)
print("Step 4: PCA Dimensionality Reduction")
print("="*70)

print(f"\nRunning PCA (n_comps={N_PCS})...")
sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')
print("PCA complete")

# Variance ratio plot
print("\nGenerating variance ratio plot...")
plt.figure(figsize=(10, 4))
sc.pl.pca_variance_ratio(adata_hvg, log=True, n_pcs=50, show=False)
plt.tight_layout()
plt.savefig(fig_dir / f"pca_variance_ratio.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"Variance ratio plot saved")

# %% [markdown]
# ---
# ## Step 5: BBKNN Batch Integration

# %%
print("\n" + "="*70)
print("Step 5: BBKNN Batch Integration")
print("="*70)

print(f"\nBBKNN parameters:")
print(f"  batch_key: {BATCH_KEY}")
print(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
print(f"  n_pcs: {BBKNN_N_PCS}")
print(f"  metric: {BBKNN_METRIC}")
print(f"  trim: {BBKNN_TRIM}")

print("\nRunning BBKNN integration...")
start_time = time.time()

bbknn.bbknn(
    adata_hvg,
    batch_key=BATCH_KEY,
    neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
    n_pcs=BBKNN_N_PCS,
    metric=BBKNN_METRIC,
    trim=BBKNN_TRIM,
    copy=False
)

elapsed = time.time() - start_time
print(f"BBKNN integration completed in {elapsed:.1f}s")
print(f"Neighbor graph constructed: {adata_hvg.obsp['connectivities'].shape}")

# %% [markdown]
# ---
# ## Step 6: UMAP and Clustering

# %%
print("\n" + "="*70)
print("Step 6: UMAP and Clustering")
print("="*70)

# UMAP
if RUN_UMAP:
    print(f"\nRunning UMAP (min_dist={UMAP_MIN_DIST})...")
    sc.tl.umap(adata_hvg, min_dist=UMAP_MIN_DIST)
    print("UMAP completed")

# Multi-resolution clustering
if RUN_CLUSTERING:
    print(f"\nRunning Leiden clustering at multiple resolutions...")
    for res in LEIDEN_RESOLUTIONS:
        cluster_key = f"leiden_bbknn_res{res}"
        print(f"  Resolution {res}...", end=" ")
        sc.tl.leiden(adata_hvg, resolution=res, key_added=cluster_key)
        n_clusters = adata_hvg.obs[cluster_key].nunique()
        print(f"{n_clusters} clusters")
    
    # Add default clustering
    default_key = f"leiden_bbknn_res{DEFAULT_RESOLUTION}"
    if default_key in adata_hvg.obs.columns:
        adata_hvg.obs['leiden_bbknn'] = adata_hvg.obs[default_key]
        print(f"\nDefault clustering: leiden_bbknn (resolution={DEFAULT_RESOLUTION})")

print("\nDimensionality reduction and clustering completed")

# %% [markdown]
# ---
# ## Step 7: Basic Visualization

# %%
print("\n" + "="*70)
print("Step 7: Basic Visualization")
print("="*70)

# UMAP by batch
if BATCH_KEY in adata_hvg.obs.columns:
    print("\nGenerating UMAP colored by batch...")
    fig = sc.pl.umap(adata_hvg, color=BATCH_KEY, return_fig=True, show=False)
    fig.savefig(fig_dir / f"umap_by_batch.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

# UMAP by clustering
if 'leiden_bbknn' in adata_hvg.obs.columns:
    print("\nGenerating UMAP colored by clusters...")
    fig = sc.pl.umap(adata_hvg, color='leiden_bbknn', legend_loc='on data', 
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"umap_by_cluster.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

# Multi-resolution clustering comparison
if RUN_CLUSTERING and len(LEIDEN_RESOLUTIONS) > 1:
    print("\nGenerating multi-resolution clustering comparison...")
    cluster_keys = [f"leiden_bbknn_res{res}" for res in LEIDEN_RESOLUTIONS]
    cluster_keys = [k for k in cluster_keys if k in adata_hvg.obs.columns]
    
    if len(cluster_keys) > 0:
        fig = sc.pl.umap(adata_hvg, color=cluster_keys, ncols=3, 
                         return_fig=True, show=False)
        fig.savefig(fig_dir / f"umap_multiresolution.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
        plt.show()
        plt.close()

print("\nBasic visualization completed")

# %%
# Marker gene expression
print("\nGenerating marker gene expression plots...")
for category, genes in available_markers.items():
    if len(genes) == 0:
        continue
    
    # Filter to genes present in current dataset
    genes_in_data = [g for g in genes if g in adata_hvg.var_names or 
                     (adata_hvg.raw is not None and g in adata_hvg.raw.var_names)]
    
    if len(genes_in_data) == 0:
        continue
    
    print(f"  {category}: {len(genes_in_data)} genes")
    
    # Plot with use_raw=True to access original expression
    fig = sc.pl.umap(adata_hvg, color=genes_in_data, use_raw=True, 
                     ncols=3, vmax='p99', return_fig=True, show=False)
    fig.savefig(fig_dir / f"umap_markers_{category}.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

print("\nMarker gene expression plots saved")

# %% [markdown]
# ---
# ## Step 8: starCAT Analysis
# 
# Now we apply starCAT for high-resolution T cell annotation based on reference GEPs.

# %%
print("\n" + "="*70)
print("Step 8: starCAT Analysis")
print("="*70)

if not RUN_STARCAT:
    print("\nstarCAT analysis is disabled (RUN_STARCAT=False)")
else:
    print(f"\nPreparing starCAT reference: {STARCAT_REFERENCE}")
    print(f"Cache directory: {cache_dir}")
    print(f"Local reference: {LOCAL_REFERENCE_TAR}")
    
    # Check if reference is already extracted
    ref_dir = cache_dir / STARCAT_REFERENCE
    ref_file = ref_dir / f"{STARCAT_REFERENCE}.reference.tsv"
    score_file = ref_dir / f"{STARCAT_REFERENCE}.scores.yaml"
    
    if not ref_file.exists():
        print(f"\nExtracting reference from: {LOCAL_REFERENCE_TAR}")
        
        # Check if tar.gz file exists
        if not Path(LOCAL_REFERENCE_TAR).exists():
            raise FileNotFoundError(f"Reference file not found: {LOCAL_REFERENCE_TAR}")
        
        # Extract tar.gz to cache directory
        import tarfile
        with tarfile.open(LOCAL_REFERENCE_TAR, 'r:gz') as tar:
            tar.extractall(path=cache_dir)
        
        print(f"Reference extracted to: {ref_dir}")
    else:
        print(f"\nUsing cached reference from: {ref_dir}")
    
    # Verify files exist
    if not ref_file.exists():
        raise FileNotFoundError(f"Reference file not found after extraction: {ref_file}")
    if not score_file.exists():
        print(f"Warning: Score file not found: {score_file}")
        print("Will initialize without scores configuration")
    
    # Initialize starCAT with local reference
    print(f"\nInitializing starCAT with local reference...")
    if score_file.exists():
        tcat = starCAT(reference=str(ref_file), score_path=str(score_file))
    else:
        tcat = starCAT(reference=str(ref_file))
    
    print(f"\nReference loaded: {tcat.ref_name}")
    print(f"Reference dimensions: {tcat.ref.shape[0]} genes × {tcat.ref.shape[1]} programs")
    
    # Display first few GEPs
    print("\nFirst 5 GEPs and 5 genes:")
    display(tcat.ref.iloc[:5, :5])
    
    # Display score information
    if hasattr(tcat, 'score_data') and tcat.score_data:
        print("\nstarCAT scores configuration:")
        if 'scores' in tcat.score_data:
            if 'continuous' in tcat.score_data['scores']:
                print(f"  Continuous scores: {[s['name'] for s in tcat.score_data['scores']['continuous']]}")
            if 'discrete' in tcat.score_data['scores']:
                print(f"  Discrete scores: {[s['name'] for s in tcat.score_data['scores']['discrete']]}")
    else:
        print("\nNote: No score configuration loaded")

# %%
# Prepare data for starCAT
# starCAT needs raw counts, so we work with adata_full
print("\nPreparing data for starCAT...")
print(f"Using full dataset: {adata_full.shape}")

# Check if we need to extract raw counts
if adata_full.raw is not None:
    print("Using raw counts from .raw")
    adata_for_starcat = sc.AnnData(
        X=adata_full.raw.X,
        obs=adata_full.obs.copy(),
        var=adata_full.raw.var.copy()
    )
else:
    print("Using .X as raw counts (assuming not log-transformed)")
    adata_for_starcat = adata_full.copy()

print(f"Data for starCAT: {adata_for_starcat.shape}")

# %%
# Run starCAT fit_transform
print("\nRunning starCAT fit_transform...")
print("This may take several minutes depending on dataset size")

start_time = time.time()

usage, scores = tcat.fit_transform(adata_for_starcat)

elapsed = time.time() - start_time
print(f"\nstarCAT completed in {elapsed:.1f}s")

print(f"\nUsage matrix: {usage.shape}")
print(f"Scores dataframe: {scores.shape}")

# Display usage head
print("\nUsage matrix (first 5 cells, first 10 programs):")
display(usage.iloc[:5, :10])

# Display scores head
print("\nScores dataframe (first 5 cells):")
display(scores.head())

# %% [markdown]
# ---
# ## Step 9: Integrate starCAT Results

# %%
print("\n" + "="*70)
print("Step 9: Integrating starCAT Results")
print("="*70)

# Add usage to adata_hvg.obs
print("\nAdding GEP usage to adata...")
for col in usage.columns:
    adata_hvg.obs[f"GEP_{col}"] = usage[col].values
print(f"Added {len(usage.columns)} GEP usage columns")

# Add scores to adata_hvg.obs
print("\nAdding scores to adata...")
for col in scores.columns:
    # Convert binary columns to categorical string for better plotting
    if col.endswith('_binary'):
        adata_hvg.obs[col] = scores[col].astype(str).values
    else:
        adata_hvg.obs[col] = scores[col].values
print(f"Added {len(scores.columns)} score columns")

# Also add to full dataset
print("\nAdding results to full dataset...")
for col in usage.columns:
    adata_full.obs[f"GEP_{col}"] = usage[col].values
for col in scores.columns:
    if col.endswith('_binary'):
        adata_full.obs[col] = scores[col].astype(str).values
    else:
        adata_full.obs[col] = scores[col].values

print("\nstarCAT results integrated successfully")

# %% [markdown]
# ---
# ## Step 10: starCAT Visualization

# %%
print("\n" + "="*70)
print("Step 10: starCAT Visualization")
print("="*70)

# 1. Discrete features from starCAT
print("\nGenerating discrete feature plots...")

discrete_features = [col for col in scores.columns if col.endswith('_binary') or col == 'Multinomial_Label']
if len(discrete_features) > 0:
    print(f"  Discrete features: {discrete_features}")
    fig = sc.pl.umap(adata_hvg, color=discrete_features, ncols=2, 
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"starcat_discrete_features.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

# %%
# 2. Continuous scores from starCAT
print("\nGenerating continuous score plots...")

continuous_features = [col for col in scores.columns if not col.endswith('_binary') and col != 'Multinomial_Label']
if len(continuous_features) > 0:
    print(f"  Continuous features: {continuous_features}")
    fig = sc.pl.umap(adata_hvg, color=continuous_features, ncols=2, vmax='p99',
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"starcat_continuous_scores.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()

# %%
# 3. Selected GEP usage patterns
print("\nGenerating GEP usage plots...")

# Select important GEPs based on their names
important_geps = []
gep_keywords = ['Cytotoxic', 'Exhaustion', 'Activation', 'Tfh', 'Treg', 'Naive', 
                'Memory', 'Proliferation', 'ISG', 'IFN']

for gep in usage.columns:
    for keyword in gep_keywords:
        if keyword.lower() in gep.lower():
            important_geps.append(f"GEP_{gep}")
            break

if len(important_geps) > 0:
    print(f"  Selected {len(important_geps)} important GEPs")
    # Plot in batches
    batch_size = 12
    for i in range(0, len(important_geps), batch_size):
        batch = important_geps[i:i+batch_size]
        batch_num = i // batch_size + 1
        print(f"    Plotting batch {batch_num}: {len(batch)} GEPs")
        
        fig = sc.pl.umap(adata_hvg, color=batch, ncols=3, vmin=0, vmax='p99',
                         return_fig=True, show=False)
        fig.savefig(fig_dir / f"starcat_gep_usage_batch{batch_num}.{FIGURE_FORMAT}", 
                    dpi=DPI, bbox_inches='tight')
        plt.show()
        plt.close()

# %%
# 4. All GEP usage (comprehensive view)
print("\nGenerating comprehensive GEP usage plot...")
print("This will generate a large multi-panel figure with all GEPs")

all_geps = [f"GEP_{col}" for col in usage.columns]
fig = sc.pl.umap(adata_hvg, color=all_geps, ncols=4, vmin=0, vmax='p99',
                 return_fig=True, show=False, size=20)
fig.savefig(fig_dir / f"starcat_all_gep_usage.{FIGURE_FORMAT}", 
            dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print("\nstarCAT visualization completed")

# %% [markdown]
# ---
# ## Step 11: Integrated Analysis - Cluster × starCAT

# %%
print("\n" + "="*70)
print("Step 11: Integrated Analysis - Cluster × starCAT")
print("="*70)

# Compare Leiden clustering with starCAT Multinomial_Label
if 'leiden_bbknn' in adata_hvg.obs.columns and 'Multinomial_Label' in adata_hvg.obs.columns:
    print("\nGenerating cluster vs starCAT label comparison...")
    
    # Create contingency table
    contingency = pd.crosstab(
        adata_hvg.obs['leiden_bbknn'],
        adata_hvg.obs['Multinomial_Label']
    )
    
    print("\nContingency table (Leiden cluster vs starCAT label):")
    display(contingency)
    
    # Save contingency table
    contingency.to_csv(output_dir / "cluster_vs_starcat_label.csv")
    print("\nContingency table saved")
    
    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Heatmap of cell counts
    sns.heatmap(contingency, annot=True, fmt='d', cmap='YlOrRd', ax=axes[0])
    axes[0].set_title('Cell Counts: Leiden Cluster vs starCAT Label')
    axes[0].set_xlabel('starCAT Multinomial Label')
    axes[0].set_ylabel('Leiden Cluster')
    
    # Heatmap of proportions (normalized by cluster)
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0)
    sns.heatmap(contingency_norm, annot=True, fmt='.2f', cmap='YlOrRd', ax=axes[1])
    axes[1].set_title('Proportions: Leiden Cluster vs starCAT Label')
    axes[1].set_xlabel('starCAT Multinomial Label')
    axes[1].set_ylabel('Leiden Cluster')
    
    plt.tight_layout()
    plt.savefig(fig_dir / f"cluster_starcat_comparison.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print("\nCluster comparison plot saved")

# %%
# GEP usage by cluster
if 'leiden_bbknn' in adata_hvg.obs.columns:
    print("\nGenerating GEP usage by cluster heatmap...")
    
    # Calculate mean GEP usage per cluster
    gep_cols = [f"GEP_{col}" for col in usage.columns]
    cluster_gep_mean = adata_hvg.obs.groupby('leiden_bbknn')[gep_cols].mean()
    
    # Transpose for better visualization (clusters as columns)
    cluster_gep_mean_T = cluster_gep_mean.T
    cluster_gep_mean_T.index = [idx.replace('GEP_', '') for idx in cluster_gep_mean_T.index]
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(max(10, len(cluster_gep_mean_T.columns)*0.5), 
                                    max(12, len(cluster_gep_mean_T.index)*0.3)))
    
    sns.heatmap(cluster_gep_mean_T, cmap='RdYlBu_r', center=0, 
                cbar_kws={'label': 'Mean GEP Usage'}, ax=ax)
    ax.set_title('Mean GEP Usage by Leiden Cluster', fontsize=14, fontweight='bold')
    ax.set_xlabel('Leiden Cluster', fontsize=12)
    ax.set_ylabel('Gene Expression Program (GEP)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f"gep_usage_by_cluster_heatmap.{FIGURE_FORMAT}", 
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    
    # Save mean GEP usage table
    cluster_gep_mean_T.to_csv(output_dir / "gep_usage_by_cluster.csv")
    print("\nGEP usage by cluster saved")

# %% [markdown]
# ---
# ## Step 12: Save Results

# %%
print("\n" + "="*70)
print("Step 12: Saving Results")
print("="*70)

# Transfer UMAP and cluster labels to full dataset
print("\nTransferring results to full dataset...")
if 'X_umap' in adata_hvg.obsm:
    adata_full.obsm['X_umap'] = adata_hvg.obsm['X_umap']
if 'X_pca' in adata_hvg.obsm:
    adata_full.obsm['X_pca'] = adata_hvg.obsm['X_pca']

# Transfer clustering results
cluster_cols = [col for col in adata_hvg.obs.columns if col.startswith('leiden_bbknn')]
for col in cluster_cols:
    adata_full.obs[col] = adata_hvg.obs[col]

print(f"Transferred {len(cluster_cols)} clustering columns")

# Save HVG dataset with starCAT results
output_hvg = output_dir / "adata_tcell_bbknn_starcat.h5ad"
print(f"\nSaving HVG dataset to: {output_hvg}")
adata_hvg.write_h5ad(output_hvg)
print(f"Saved: {output_hvg}")

# Save full dataset with all annotations
output_full = output_dir / "adata_tcell_full_annotated.h5ad"
print(f"\nSaving full dataset to: {output_full}")
adata_full.write_h5ad(output_full)
print(f"Saved: {output_full}")

# Save starCAT results separately
print("\nSaving starCAT results...")
usage.to_csv(output_dir / "starcat_usage_matrix.csv")
scores.to_csv(output_dir / "starcat_scores.csv")
print("starCAT results saved")

# %% [markdown]
# ---
# ## Step 13: Analysis Summary

# %%
print("\n" + "="*70)
print("Step 13: Analysis Summary")
print("="*70)

summary = []
summary.append("="*70)
summary.append("T CELL ANALYSIS WITH starCAT - SUMMARY REPORT")
summary.append("="*70)
summary.append(f"\nAnalysis completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary.append("")

summary.append("[ Dataset Information ]")
summary.append(f"  Input: {INPUT_H5AD_PATH}")
summary.append(f"  Cells: {adata_full.shape[0]:,}")
summary.append(f"  Genes (total): {adata_full.shape[1]:,}")
summary.append(f"  Genes (HVG): {adata_hvg.shape[1]:,}")
summary.append("")

summary.append("[ Batch Information ]")
if BATCH_KEY in adata_full.obs.columns:
    batch_counts = adata_full.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata_full.n_obs * 100
        summary.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
summary.append("")

summary.append("[ Integration Method ]")
summary.append(f"  Method: BBKNN")
summary.append(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
summary.append(f"  n_pcs: {BBKNN_N_PCS}")
summary.append(f"  metric: {BBKNN_METRIC}")
summary.append("")

summary.append("[ Clustering Results ]")
if 'leiden_bbknn' in adata_hvg.obs.columns:
    n_clusters = adata_hvg.obs['leiden_bbknn'].nunique()
    summary.append(f"  Default clustering (res={DEFAULT_RESOLUTION}): {n_clusters} clusters")
    summary.append(f"  Resolutions tested: {LEIDEN_RESOLUTIONS}")
summary.append("")

summary.append("[ starCAT Annotation ]")
summary.append(f"  Reference: {STARCAT_REFERENCE}")
summary.append(f"  GEPs: {usage.shape[1]}")
summary.append(f"  Scores: {scores.shape[1]}")
if 'Multinomial_Label' in scores.columns:
    label_counts = adata_hvg.obs['Multinomial_Label'].value_counts()
    summary.append("\n  Cell type distribution (starCAT):")
    for label, count in label_counts.items():
        pct = count / adata_hvg.n_obs * 100
        summary.append(f"    {label}: {count:,} cells ({pct:.1f}%)")
summary.append("")

summary.append("[ Output Files ]")
summary.append(f"  HVG dataset: {output_dir / 'adata_tcell_bbknn_starcat.h5ad'}")
summary.append(f"  Full dataset: {output_dir / 'adata_tcell_full_annotated.h5ad'}")
summary.append(f"  starCAT usage: {output_dir / 'starcat_usage_matrix.csv'}")
summary.append(f"  starCAT scores: {output_dir / 'starcat_scores.csv'}")
summary.append(f"  Figures: {fig_dir / '*.{FIGURE_FORMAT}'}")
summary.append("")

summary.append("="*70)
summary.append("ANALYSIS COMPLETED SUCCESSFULLY")
summary.append("="*70)

summary_text = '\n'.join(summary)

# Save summary
summary_file = output_dir / "analysis_summary.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)

print(summary_text)
print(f"\nSummary saved to: {summary_file}")

# %% [markdown]
# ---
# ## Additional Analysis (Optional)
# 
# Below are some optional analyses you can perform with the annotated dataset.

# %%
# Example: Explore specific GEP patterns
# Select a GEP of interest and examine which clusters express it

gep_of_interest = "Exhaustion"  # Change to any GEP from usage.columns

if gep_of_interest in usage.columns:
    print(f"\nAnalyzing GEP: {gep_of_interest}")
    
    # Distribution by cluster
    if 'leiden_bbknn' in adata_hvg.obs.columns:
        print(f"\nViolin plot: GEP usage by cluster")
        sc.pl.violin(adata_hvg, keys=f"GEP_{gep_of_interest}", 
                     groupby='leiden_bbknn', show=False)
        plt.title(f'GEP Usage: {gep_of_interest}')
        plt.tight_layout()
        plt.show()
        
        print(f"\nUMAP: GEP expression pattern")
        sc.pl.umap(adata_hvg, color=f"GEP_{gep_of_interest}", 
                   vmax='p99', show=False)
        plt.tight_layout()
        plt.show()

# %%
# Example: Compare ASA scores across clusters
if 'ASA' in adata_hvg.obs.columns and 'leiden_bbknn' in adata_hvg.obs.columns:
    print("\nComparing ASA (Activation) scores across clusters")
    
    print("\nViolin plot: ASA score distribution")
    sc.pl.violin(adata_hvg, keys='ASA', groupby='leiden_bbknn', show=False)
    plt.title('ASA Score Distribution by Cluster')
    plt.tight_layout()
    plt.show()
    
    print("\nUMAP: ASA score pattern")
    sc.pl.umap(adata_hvg, color='ASA', vmax='p99', show=False)
    plt.title('ASA Score on UMAP')
    plt.tight_layout()
    plt.show()

# %%
# Example: Correlation between GEPs
# Compute correlation matrix of GEP usage

print("\nComputing GEP correlation matrix...")
gep_corr = usage.corr()

# Plot correlation heatmap
fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(gep_corr, cmap='coolwarm', center=0, 
            cbar_kws={'label': 'Correlation'}, ax=ax)
ax.set_title('GEP Usage Correlation Matrix', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.show()

# Find highly correlated GEPs
print("\nHighly correlated GEP pairs (|r| > 0.7):")
high_corr = []
for i in range(len(gep_corr.columns)):
    for j in range(i+1, len(gep_corr.columns)):
        corr_val = gep_corr.iloc[i, j]
        if abs(corr_val) > 0.7:
            high_corr.append({
                'GEP1': gep_corr.columns[i],
                'GEP2': gep_corr.columns[j],
                'Correlation': corr_val
            })

if high_corr:
    high_corr_df = pd.DataFrame(high_corr).sort_values('Correlation', ascending=False)
    display(high_corr_df)
else:
    print("  No GEP pairs with |r| > 0.7 found")

# %% [markdown]
# ---
# ## End of Analysis
# 
# **Key takeaways:**
# 
# 1. **BBKNN integration** successfully corrected batch effects while preserving biological variation
# 2. **starCAT annotation** provided high-resolution characterization of T cell states using 52 GEPs
# 3. **Complementary approaches**: Leiden clustering captures major cell populations, while starCAT GEPs capture functional states that may span multiple clusters
# 
# **Next steps:**
# - Validate key findings with marker genes
# - Perform differential expression analysis between clusters/cell types
# - Investigate specific GEPs of interest (e.g., Exhaustion, Cytotoxicity)
# - Compare results across different tissue sites or disease conditions
# - Integrate with TCR sequencing data if available


