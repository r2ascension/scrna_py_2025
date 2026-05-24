# -*- coding: utf-8 -*-
# %% [markdown]
# # Epithelial Cell Subcluster Analysis v4.5-HOTFIX (Production Ready)
#
# **Purpose:**
# - Perform subcluster analysis on already-integrated epithelial cells
# - Intelligent small cluster handling (merge or remove)
# - Robust marker gene analysis with proper QC
# - Production-quality error handling and validation
#
# **Author:** r2end
# **Date:** 2025-01-22
# **Version:** 4.5-HOTFIX - All P0/P1 issues fixed
#
# **CRITICAL FIXES in HOTFIX:**
# - P0: All emoji/special chars replaced with ASCII (encoding safety)
# - P0: Write-back uses obs_names alignment (handles remove strategy)
# - P0: Forced categorical conversion for cluster columns
# - P0: Safe neighbors_within_batch calculation
# - P1: Removed restrict_to (simplified logic)
# - P1: Connectivities-based merging (with PCA fallback)
# - P1: Marker analysis skips clusters < MIN_CELLS_FOR_MARKER
# - P1: Outputs both full and top-N marker files
# - P1: Regex uses .contains() for comprehensive filtering
# - P2: Optimized sparse matrix operations
# - P2: Per-celltype immediate disk output

# %% [markdown]
# ## Configuration Parameters

# %%
import os
import sys
import gc
import warnings
from pathlib import Path
from datetime import datetime
import time
import random
import re

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import csr_matrix
from sklearn.metrics.pairwise import cosine_similarity

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# =============================================================================
# PIPELINE METADATA
# =============================================================================
PIPELINE_VERSION = "v4.5-HOTFIX-PRODUCTION"

# =============================================================================
# SMALL CLUSTER HANDLING STRATEGY
# =============================================================================
# Options: "merge", "remove", "keep"
# - "merge": Merge small clusters to nearest neighbor (RECOMMENDED for epithelial)
# - "remove": Remove small clusters and recluster remaining cells
# - "keep": Keep all clusters as-is (no filtering)
SMALL_CLUSTER_STRATEGY = "merge"

# Thresholds
MIN_CELLS_PER_CLUSTER = 30          # Clusters below this are considered "small"
MIN_CELLS_FOR_MARKER = 50           # Minimum cells for reliable marker analysis (STRICTER)
MERGE_MAX_ITERATIONS = 3            # Maximum iterations for merging small clusters

# =============================================================================
# GLOBAL NEIGHBORS GRAPH SELECTION (BBKNN)
# =============================================================================
NEIGHBORS_KEY = "neighbors"

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
try:
    from anndata import settings as anndata_settings
    anndata_settings.seed = RANDOM_SEED
except Exception:
    pass

print(f"[INFO] Random seed set to {RANDOM_SEED}")

# =============================================================================
# INPUT/OUTPUT
# =============================================================================
INPUT_H5AD = "/home/h2048/data/py/0110/celltypist_epithelial/epithelial_celltypist_filtered_final.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%m%d')}/epithelial_subcluster_v4_5_hotfix")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# ANNOTATION KEYS
# =============================================================================
BATCH_KEY = "dataset"
CELLTYPE_COLUMN = "celltypist_pred"

# =============================================================================
# SUBCLUSTERING
# =============================================================================
MIN_CELLS_FOR_SUBCLUSTER = 100
MIN_CELLS_PER_BATCH = 3

RESOLUTION_FIXED = 0.2

# =============================================================================
# MARKER GENE SETTINGS
# =============================================================================
MARKER_METHOD = "wilcoxon"
MARKER_LOGFC_THRESHOLD = 0.5
MARKER_PADJ_THRESHOLD = 0.05
MARKER_MIN_PCT = 0.25
MARKER_MIN_DELTA_PCT = 0.10
TOP_N_MARKERS = 20

# =============================================================================
# MARKER GENE FILTERING
# =============================================================================
FILTER_MT_GENES = True
FILTER_RIBO_GENES = True
FILTER_UNANNOTATED = True
FILTER_HISTONE_GENES = True
FILTER_STRESS_GENES = True
STRESS_MATCH_THRESHOLD = 0.5

STRESS_SIGNATURE_GENES = [
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "ATF3",
    "DUSP1", "DUSP2", "DUSP4", "DUSP5",
    "HSPA1A", "HSPA1B",
    "IER2", "IER3", "IER5",
    "NR4A1", "NR4A2", "NR4A3",
    "ZFP36", "ZFP36L1", "ZFP36L2"
]

# =============================================================================
# GENE UNIVERSE / DE INPUT SELECTION
# =============================================================================
USE_RAW_FOR_GENE_FILTERING = True
ALSO_BUILD_VAR_FILTERS = True

DE_MODE = "auto"
DE_LAYER = None

# =============================================================================
# VISUALIZATION
# =============================================================================
FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"
UMAP_SIZE = 5
UMAP_ALPHA = 0.6

# =============================================================================
# PERFORMANCE
# =============================================================================
N_JOBS = 48
MEMORY_LIMIT_GB = 256

TOP_DE_GENES = 4000

# =============================================================================
# REGEX / SETS
# =============================================================================
STRESS_SET = set(STRESS_SIGNATURE_GENES)
# P1: Use more comprehensive unannotated pattern
UNANNOTATED_REGEX = r"(LINC|RP11-|^RP\d+-|^AC\d+|^AL\d+|CTD-|CTB-|pseudogene)"

sc.settings.n_jobs = N_JOBS
sc.settings.figdir = OUTPUT_DIR / "figures"
sc.settings.figdir.mkdir(exist_ok=True)

print("=" * 80)
print(f"Epithelial Subcluster Analysis Pipeline {PIPELINE_VERSION}")
print("=" * 80)
print(f"Input:              {INPUT_H5AD}")
print(f"Output:             {OUTPUT_DIR}")
print(f"Cell type:          {CELLTYPE_COLUMN}")
print(f"Batch key:          {BATCH_KEY if BATCH_KEY else 'None'}")
print(f"Resolution:         FIXED {RESOLUTION_FIXED}")
print(f"Cluster strategy:   {SMALL_CLUSTER_STRATEGY}")
print(f"Min cluster size:   {MIN_CELLS_PER_CLUSTER}")
print(f"Min for markers:    {MIN_CELLS_FOR_MARKER}")
print(f"DE mode:            {DE_MODE} (layer={DE_LAYER})")
print(f"Date:               {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# %% [markdown]
# ## Helper Functions for Small Cluster Handling

# %%
def merge_small_clusters_connectivities_based(adata, cluster_key='leiden', 
                                             min_cells=30, max_iterations=3,
                                             conn_key='connectivities',
                                             fallback_to_pca=True):
    """
    Merge small clusters to nearest neighbor based on BBKNN connectivities graph.
    Falls back to PCA-based if connectivities not available.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with clustering results
    cluster_key : str
        Column name in adata.obs containing cluster assignments
    min_cells : int
        Minimum number of cells required for a cluster
    max_iterations : int
        Maximum number of merge iterations
    conn_key : str
        Key for connectivities matrix in adata.obsp
    fallback_to_pca : bool
        If True, use PCA-based merging when connectivities unavailable
    
    Returns
    -------
    adata : AnnData
        Modified object with merged clusters
    merge_log : list
        List of (iteration, small_cluster, target_cluster, n_cells, score) tuples
    """
    
    print(f"\n{'='*80}")
    print("SMALL CLUSTER MERGING")
    print(f"{'='*80}")
    print(f"Strategy: Merge clusters with < {min_cells} cells to nearest neighbor")
    print(f"Maximum iterations: {max_iterations}")
    
    # P0: Force categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        print("[INFO] Converting cluster column to categorical")
        adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
    
    # Determine merge method
    use_connectivities = conn_key in adata.obsp
    if use_connectivities:
        print(f"Method: Connectivities-based (using {conn_key})")
        C = adata.obsp[conn_key].tocsr()
    elif fallback_to_pca and 'X_pca' in adata.obsm:
        print("[WARN] Connectivities not found, falling back to PCA-based merging")
        use_connectivities = False
    else:
        print("[WARN] No connectivities or PCA available")
        if 'X_pca' not in adata.obsm:
            print("[INFO] Computing PCA for merging...")
            sc.tl.pca(adata, n_comps=50, random_state=RANDOM_SEED)
        use_connectivities = False
    
    print()
    
    merge_log = []
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        print(f"--- Iteration {iteration} ---")
        
        # Get cluster sizes
        cluster_counts = adata.obs[cluster_key].value_counts()
        small_clusters = cluster_counts[cluster_counts < min_cells].index.tolist()
        
        if len(small_clusters) == 0:
            print(f"[OK] No small clusters remaining. Merging complete.")
            break
            
        print(f"Found {len(small_clusters)} small clusters: {small_clusters}")
        
        clusters = adata.obs[cluster_key].cat.categories
        
        # For each small cluster, find best merge target
        for small_cluster in small_clusters:
            n_cells_small = cluster_counts[small_cluster]
            
            # Get large clusters
            large_clusters = [c for c in clusters 
                            if c != small_cluster and c not in small_clusters]
            
            if len(large_clusters) == 0:
                print(f"  [WARN] No large clusters available to merge {small_cluster}")
                continue
            
            # Compute merge scores
            if use_connectivities:
                # Connectivities-based: sum of edge weights
                idx_small = np.where(adata.obs[cluster_key] == small_cluster)[0]
                
                scores = {}
                for large_cluster in large_clusters:
                    idx_large = np.where(adata.obs[cluster_key] == large_cluster)[0]
                    # Sum of connectivities between small and large cluster
                    score = C[idx_small][:, idx_large].sum()
                    scores[large_cluster] = float(score)
                
                target_cluster = max(scores, key=scores.get)
                max_score = scores[target_cluster]
                score_type = "conn_sum"
                
            else:
                # PCA-based: cosine similarity of centroids
                latent = adata.obsm['X_pca']
                
                small_centroid = latent[adata.obs[cluster_key] == small_cluster].mean(axis=0)
                
                similarities = {}
                for large_cluster in large_clusters:
                    large_centroid = latent[adata.obs[cluster_key] == large_cluster].mean(axis=0)
                    sim = cosine_similarity(
                        small_centroid.reshape(1, -1),
                        large_centroid.reshape(1, -1)
                    )[0, 0]
                    similarities[large_cluster] = sim
                
                target_cluster = max(similarities, key=similarities.get)
                max_score = similarities[target_cluster]
                score_type = "cosine_sim"
            
            # Merge
            mask = adata.obs[cluster_key] == small_cluster
            adata.obs.loc[mask, cluster_key] = target_cluster
            
            print(f"  -> Merged cluster {small_cluster} ({n_cells_small} cells) "
                  f"into {target_cluster} ({score_type}={max_score:.3f})")
            
            merge_log.append((iteration, small_cluster, target_cluster, 
                            n_cells_small, max_score))
        
        # Update categories to remove empty clusters
        adata.obs[cluster_key] = adata.obs[cluster_key].cat.remove_unused_categories()
        
        print()
    
    # Final summary
    print(f"\n{'='*80}")
    print("MERGING SUMMARY")
    print(f"{'='*80}")
    print(f"Total iterations: {iteration}")
    print(f"Total merges: {len(merge_log)}")
    
    if len(merge_log) > 0:
        print("\nMerge history:")
        for iter_num, small, target, n_cells, score in merge_log:
            print(f"  Iter {iter_num}: {small} ({n_cells} cells) -> {target} (score={score:.3f})")
    
    final_counts = adata.obs[cluster_key].value_counts().sort_index()
    print(f"\nFinal cluster distribution:")
    for cluster, count in final_counts.items():
        status = "[OK]" if count >= min_cells else "[WARN] STILL SMALL"
        print(f"  Cluster {cluster}: {count:,} cells {status}")
    
    remaining_small = final_counts[final_counts < min_cells]
    if len(remaining_small) > 0:
        print(f"\n[WARN] {len(remaining_small)} clusters still below threshold")
        print("   Consider increasing max_iterations or using 'remove' strategy")
    
    print(f"{'='*80}\n")
    
    return adata, merge_log


def remove_small_clusters_and_recluster(adata, cluster_key='leiden', min_cells=30,
                                       batch_key='batch', resolution=0.3, n_pcs=50,
                                       conn_key='connectivities'):
    """
    Remove small clusters and recluster remaining cells.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object with clustering results
    cluster_key : str
        Column name in adata.obs containing cluster assignments
    min_cells : int
        Minimum number of cells required for a cluster
    batch_key : str
        Batch key for BBKNN
    resolution : float
        Leiden clustering resolution
    n_pcs : int
        Number of PCs for BBKNN
    conn_key : str
        Key for connectivities in obsp
    
    Returns
    -------
    adata : AnnData
        Filtered and reclustered object
    removed_clusters : list
        List of removed cluster IDs
    """
    
    print(f"\n{'='*80}")
    print("SMALL CLUSTER REMOVAL & RECLUSTERING")
    print(f"{'='*80}")
    print(f"Strategy: Remove clusters with < {min_cells} cells and recluster")
    print()
    
    # P0: Force categorical
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        print("[INFO] Converting cluster column to categorical")
        adata.obs[cluster_key] = adata.obs[cluster_key].astype('category')
    
    # Identify small clusters
    cluster_counts = adata.obs[cluster_key].value_counts()
    small_clusters = cluster_counts[cluster_counts < min_cells].index.tolist()
    
    if len(small_clusters) == 0:
        print("[OK] No small clusters to remove.")
        return adata, []
    
    print(f"Removing {len(small_clusters)} small clusters:")
    for cluster in small_clusters:
        n_cells = cluster_counts[cluster]
        print(f"  - Cluster {cluster}: {n_cells} cells")
    
    n_cells_before = adata.n_obs
    
    # Filter out small clusters
    adata = adata[~adata.obs[cluster_key].isin(small_clusters)].copy()
    
    n_cells_after = adata.n_obs
    n_cells_removed = n_cells_before - n_cells_after
    
    print(f"\nRemoved {n_cells_removed:,} cells ({n_cells_removed/n_cells_before*100:.1f}%)")
    print(f"Remaining: {n_cells_after:,} cells")
    
    # P0: Filter small batches AGAIN in this subset
    print(f"\n[INFO] Re-checking batch sizes in remaining cells...")
    batch_counts = adata.obs[batch_key].value_counts()
    small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index
    
    if len(small_batches) > 0:
        print(f"[WARN] Removing {len(small_batches)} small batches:")
        for batch in small_batches:
            print(f"    {batch}: {batch_counts[batch]} cells")
        
        adata = adata[~adata.obs[batch_key].isin(small_batches)].copy()
        print(f"[OK] Remaining after batch filter: {adata.n_obs:,} cells")
    
    # Recompute BBKNN and clustering
    print("\nRecomputing BBKNN and Leiden clustering...")
    
    # P0: Safe neighbors_within_batch calculation
    batch_sizes = adata.obs[batch_key].value_counts()
    min_batch_size = int(batch_sizes.min())
    
    if min_batch_size < 2:
        raise ValueError(
            f"Smallest batch has {min_batch_size} cells (<2). "
            "BBKNN requires at least 2 cells per batch. "
            "Increase MIN_CELLS_PER_BATCH or filter this cell type."
        )
    
    # P0: Strict formula
    neighbors_within_batch = min(5, min_batch_size - 1)
    neighbors_within_batch = max(1, neighbors_within_batch)
    
    print(f"  Smallest batch size: {min_batch_size}")
    print(f"  Using neighbors_within_batch: {neighbors_within_batch}")
    
    # BBKNN
    sce.pp.bbknn(
        adata,
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=None
    )
    
    # Leiden
    sc.tl.leiden(adata, resolution=resolution, key_added=cluster_key, random_state=42)
    
    n_clusters_new = adata.obs[cluster_key].nunique()
    print(f"[OK] New clustering: {n_clusters_new} clusters")
    
    # Show new cluster sizes
    new_counts = adata.obs[cluster_key].value_counts().sort_index()
    print("\nNew cluster distribution:")
    for cluster, count in new_counts.items():
        status = "[OK]" if count >= min_cells else "[WARN] STILL SMALL"
        print(f"  Cluster {cluster}: {count:,} cells {status}")
    
    print(f"{'='*80}\n")
    
    return adata, small_clusters


# %% [markdown]
# ## Step 1: Load Data

# %%
print("\n" + "=" * 80)
print("STEP 1: Loading Data")
print("=" * 80)

t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time() - t0:.1f}s")

print("\nDataset Summary:")
print(f"  Cells: {adata.n_obs:,}")
print(f"  Genes: {adata.n_vars:,}")
print(f"  Layers: {list(adata.layers.keys()) if adata.layers else 'None'}")
print(f"  Has .raw: {adata.raw is not None}")

if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    print(f"\nBatch distribution ({BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    print(f"  Total batches: {len(batch_counts)}")
    print(f"  Median batch size: {int(batch_counts.median())}")
    print(f"  Range: {batch_counts.min()} - {batch_counts.max()}")
    
    small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH]
    if len(small_batches) > 0:
        print(f"  [WARN] {len(small_batches)} batches have < {MIN_CELLS_PER_BATCH} cells")

if CELLTYPE_COLUMN in adata.obs.columns:
    print(f"\nCell type distribution ({CELLTYPE_COLUMN}):")
    ct_counts = adata.obs[CELLTYPE_COLUMN].value_counts()
    for ct, count in ct_counts.head(10).items():
        print(f"  {ct}: {count:,}")
    if len(ct_counts) > 10:
        print(f"  ... and {len(ct_counts) - 10} more")

# %% [markdown]
# ## Step 2: Data Validation & Preparation

# %%
print("\n" + "=" * 80)
print("STEP 2: Data Validation & Preparation")
print("=" * 80)

# Check connectivity structures
print("\nChecking neighbors/connectivity...")
conn_found = False
CONN_KEY = None

if NEIGHBORS_KEY in adata.uns:
    if 'connectivities_key' in adata.uns[NEIGHBORS_KEY]:
        CONN_KEY = adata.uns[NEIGHBORS_KEY]['connectivities_key']
        if CONN_KEY in adata.obsp:
            print(f"  [OK] Found connectivities via uns['{NEIGHBORS_KEY}']['connectivities_key']: '{CONN_KEY}'")
            conn_found = True
        else:
            print(f"  [WARN] Key '{CONN_KEY}' referenced but not in .obsp")
    else:
        possible_conn = f"{NEIGHBORS_KEY}_connectivities"
        if possible_conn in adata.obsp:
            CONN_KEY = possible_conn
            print(f"  [OK] Found connectivities: '{CONN_KEY}'")
            conn_found = True

if not conn_found:
    for key in ['connectivities', 'bbknn_connectivities', 'neighbors_connectivities']:
        if key in adata.obsp:
            CONN_KEY = key
            print(f"  [OK] Found connectivities: '{CONN_KEY}'")
            conn_found = True
            break

if not conn_found:
    raise ValueError(
        "No connectivities found! Expected BBKNN-generated graph.\n"
        "Available keys in .obsp: " + str(list(adata.obsp.keys())) + "\n"
        "This pipeline requires pre-existing BBKNN integration."
    )

print(f"\nConnectivity matrix: {CONN_KEY}")
print(f"  Shape: {adata.obsp[CONN_KEY].shape}")
print(f"  Stored edges: {adata.obsp[CONN_KEY].nnz:,}")

# Remove small batches
if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    print(f"\nFiltering small batches (< {MIN_CELLS_PER_BATCH} cells)...")
    batch_counts_before = adata.obs[BATCH_KEY].value_counts()
    small_batches = batch_counts_before[batch_counts_before < MIN_CELLS_PER_BATCH].index
    
    if len(small_batches) > 0:
        print(f"  Removing {len(small_batches)} small batches:")
        for batch in small_batches:
            print(f"    {batch}: {batch_counts_before[batch]} cells")
        
        n_before = adata.n_obs
        adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
        n_after = adata.n_obs
        
        print(f"  [OK] Removed {n_before - n_after:,} cells")
        print(f"  Remaining: {n_after:,} cells in {adata.obs[BATCH_KEY].nunique()} batches")
    else:
        print("  [OK] No small batches to remove")

# Ensure sparse matrix
if not isinstance(adata.X, sparse.spmatrix):
    print("\nConverting .X to sparse matrix...")
    adata.X = csr_matrix(adata.X)
    print("  [OK] Converted to sparse CSR")

# %% [markdown]
# ## Step 3: Cell Type Processing Loop

# %%
print("\n" + "=" * 80)
print("STEP 3: Cell Type Processing")
print("=" * 80)

# Get cell types to process
celltypes = adata.obs[CELLTYPE_COLUMN].unique()
celltypes = sorted([ct for ct in celltypes if pd.notna(ct)])

print(f"\nFound {len(celltypes)} cell types:")
for ct in celltypes:
    n_cells = (adata.obs[CELLTYPE_COLUMN] == ct).sum()
    status = "[OK]" if n_cells >= MIN_CELLS_FOR_SUBCLUSTER else "[WARN] TOO SMALL"
    print(f"  {ct}: {n_cells:,} cells {status}")

# Filter cell types by minimum size
celltypes_to_process = [
    ct for ct in celltypes 
    if (adata.obs[CELLTYPE_COLUMN] == ct).sum() >= MIN_CELLS_FOR_SUBCLUSTER
]

print(f"\nProcessing {len(celltypes_to_process)} cell types")

# Initialize results storage
processed_celltypes = []
all_markers = {}
all_markers_filtered = {}
marker_filter_stats = {}

# Gene filter dictionaries
gene_filter_dict_raw = None
gene_filter_dict_var = None

# DE selection
DE_USE_RAW = False
DE_USE_LAYER = None
de_note = ""

# %% [markdown]
# ## Step 4-7: Per-CellType Processing Loop

# %%
for celltype_idx, celltype in enumerate(celltypes_to_process, 1):
    print("\n" + "=" * 80)
    print(f"PROCESSING CELL TYPE [{celltype_idx}/{len(celltypes_to_process)}]: {celltype}")
    print("=" * 80)
    
    ct_start_time = time.time()
    
    # === 4.1. Filter cells for this cell type ===
    print(f"\n[INFO] Filtering cells for {celltype}...")
    ct_mask = adata.obs[CELLTYPE_COLUMN] == celltype
    adata_ct = adata[ct_mask].copy()
    n_cells_ct = adata_ct.n_obs
    
    print(f"  [OK] Selected {n_cells_ct:,} cells ({n_cells_ct/adata.n_obs*100:.1f}% of dataset)")
    
    # === 4.2. Compute subclusters (P1: removed restrict_to) ===
    print(f"\n[INFO] Computing subclusters...")
    print(f"  Resolution: {RESOLUTION_FIXED}")
    print(f"  Using global graph: {CONN_KEY}")
    
    # P1: Simplified - no restrict_to
    sc.tl.leiden(
        adata_ct,
        resolution=RESOLUTION_FIXED,
        key_added='subcluster_leiden',
        adjacency=adata_ct.obsp[CONN_KEY],
        random_state=RANDOM_SEED
    )
    
    print("  [OK] Leiden clustering completed")
    
    # P0: Force categorical
    if not pd.api.types.is_categorical_dtype(adata_ct.obs['subcluster_leiden']):
        adata_ct.obs['subcluster_leiden'] = adata_ct.obs['subcluster_leiden'].astype('category')
    
    # Initial cluster statistics
    cluster_counts_initial = adata_ct.obs['subcluster_leiden'].value_counts().sort_index()
    n_clusters_initial = len(cluster_counts_initial)
    
    print(f"\n  Initial clusters: {n_clusters_initial}")
    print("  Cluster sizes:")
    for cluster, count in cluster_counts_initial.items():
        status = "[OK]" if count >= MIN_CELLS_PER_CLUSTER else "[WARN] SMALL"
        print(f"    Cluster {cluster}: {count:,} cells {status}")
    
    # === 4.3. Handle small clusters ===
    small_cluster_log = []
    removed_clusters = []
    removed_cell_ids = []
    
    if SMALL_CLUSTER_STRATEGY == "merge":
        print(f"\n[INFO] Merging small clusters (strategy: {SMALL_CLUSTER_STRATEGY})...")
        
        adata_ct, merge_log = merge_small_clusters_connectivities_based(
            adata_ct,
            cluster_key='subcluster_leiden',
            min_cells=MIN_CELLS_PER_CLUSTER,
            max_iterations=MERGE_MAX_ITERATIONS,
            conn_key=CONN_KEY,
            fallback_to_pca=True
        )
        
        small_cluster_log = merge_log
        
    elif SMALL_CLUSTER_STRATEGY == "remove":
        print(f"\n[INFO] Removing small clusters (strategy: {SMALL_CLUSTER_STRATEGY})...")
        
        # Store IDs before removal
        cells_before_removal = adata_ct.obs_names.copy()
        
        adata_ct, removed_clusters = remove_small_clusters_and_recluster(
            adata_ct,
            cluster_key='subcluster_leiden',
            min_cells=MIN_CELLS_PER_CLUSTER,
            batch_key=BATCH_KEY,
            resolution=RESOLUTION_FIXED,
            n_pcs=50,
            conn_key=CONN_KEY
        )
        
        # P0: Track removed cell IDs
        removed_cell_ids = cells_before_removal.difference(adata_ct.obs_names).tolist()
        print(f"  [INFO] Tracked {len(removed_cell_ids)} removed cell IDs for write-back")
        
    elif SMALL_CLUSTER_STRATEGY == "keep":
        print(f"\n[INFO] Keeping all clusters (strategy: {SMALL_CLUSTER_STRATEGY})")
        print("  [WARN] Small clusters will be retained as-is")
    
    else:
        raise ValueError(f"Unknown strategy: {SMALL_CLUSTER_STRATEGY}")
    
    # Final cluster statistics
    cluster_counts_final = adata_ct.obs['subcluster_leiden'].value_counts().sort_index()
    n_clusters_final = len(cluster_counts_final)
    
    print(f"\n[INFO] Final cluster distribution:")
    print(f"  Total clusters: {n_clusters_final}")
    for cluster, count in cluster_counts_final.items():
        status = "[OK]" if count >= MIN_CELLS_PER_CLUSTER else "[WARN] SMALL"
        print(f"    Cluster {cluster}: {count:,} cells {status}")
    
    # Create final subcluster label
    adata_ct.obs['subcluster'] = (
        adata_ct.obs[CELLTYPE_COLUMN].astype(str) + "_" + 
        adata_ct.obs['subcluster_leiden'].astype(str)
    )
    
    # === 4.4. Build gene filter dictionaries (only on first cell type) ===
    if celltype_idx == 1:
        print("\n[INFO] Building gene filter dictionaries (first cell type)...")
        
        # Determine gene universe
        if USE_RAW_FOR_GENE_FILTERING and adata_ct.raw is not None:
            gene_names_for_filter = adata_ct.raw.var_names
            filter_source = "raw"
        else:
            gene_names_for_filter = adata_ct.var_names
            filter_source = "var"
        
        print(f"  Gene universe: {len(gene_names_for_filter):,} genes (from {filter_source})")
        
        # Build filter flags
        gene_filter_flags = pd.DataFrame(index=gene_names_for_filter)
        
        # MT genes
        if FILTER_MT_GENES:
            mt_pattern = r"^MT-|^mt-"
            # P1: Use contains for safety
            gene_filter_flags['is_mt'] = gene_filter_flags.index.str.contains(
                mt_pattern, case=False, regex=True, na=False
            )
            n_mt = gene_filter_flags['is_mt'].sum()
            print(f"  MT genes: {n_mt}")
        
        # Ribosomal genes
        if FILTER_RIBO_GENES:
            ribo_pattern = r"^RPL|^RPS|^MRPL|^MRPS"
            gene_filter_flags['is_ribo'] = gene_filter_flags.index.str.contains(
                ribo_pattern, case=False, regex=True, na=False
            )
            n_ribo = gene_filter_flags['is_ribo'].sum()
            print(f"  Ribosomal genes: {n_ribo}")
        
        # Histone genes
        if FILTER_HISTONE_GENES:
            histone_pattern = r"^HIST\d"
            gene_filter_flags['is_histone'] = gene_filter_flags.index.str.contains(
                histone_pattern, case=False, regex=True, na=False
            )
            n_histone = gene_filter_flags['is_histone'].sum()
            print(f"  Histone genes: {n_histone}")
        
        # P1: Unannotated genes - use contains
        if FILTER_UNANNOTATED:
            gene_filter_flags['is_unannotated'] = gene_filter_flags.index.str.contains(
                UNANNOTATED_REGEX, case=False, regex=True, na=False
            )
            n_unannotated = gene_filter_flags['is_unannotated'].sum()
            print(f"  Unannotated genes: {n_unannotated}")
        
        # Stress/IEG genes
        if FILTER_STRESS_GENES and len(STRESS_SET) > 0:
            gene_filter_flags['is_stress'] = gene_filter_flags.index.isin(STRESS_SET)
            n_stress = gene_filter_flags['is_stress'].sum()
            print(f"  Stress/IEG genes: {n_stress}")
        
        # Create filter dictionary
        filter_dict = gene_filter_flags.to_dict('index')
        
        if filter_source == "raw":
            gene_filter_dict_raw = filter_dict
            print(f"  [OK] Created gene_filter_dict_raw with {len(gene_filter_dict_raw):,} entries")
        else:
            gene_filter_dict_var = filter_dict
            print(f"  [OK] Created gene_filter_dict_var with {len(gene_filter_dict_var):,} entries")
        
        # Also build var filter if requested
        if ALSO_BUILD_VAR_FILTERS and filter_source == "raw":
            var_filter_flags = pd.DataFrame(index=adata_ct.var_names)
            if FILTER_MT_GENES:
                var_filter_flags['is_mt'] = var_filter_flags.index.str.contains(
                    r"^MT-|^mt-", case=False, regex=True, na=False
                )
            if FILTER_RIBO_GENES:
                var_filter_flags['is_ribo'] = var_filter_flags.index.str.contains(
                    r"^RPL|^RPS|^MRPL|^MRPS", case=False, regex=True, na=False
                )
            if FILTER_HISTONE_GENES:
                var_filter_flags['is_histone'] = var_filter_flags.index.str.contains(
                    r"^HIST\d", case=False, regex=True, na=False
                )
            if FILTER_UNANNOTATED:
                var_filter_flags['is_unannotated'] = var_filter_flags.index.str.contains(
                    UNANNOTATED_REGEX, case=False, regex=True, na=False
                )
            if FILTER_STRESS_GENES:
                var_filter_flags['is_stress'] = var_filter_flags.index.isin(STRESS_SET)
            
            gene_filter_dict_var = var_filter_flags.to_dict('index')
            print(f"  [OK] Also created gene_filter_dict_var with {len(gene_filter_dict_var):,} entries")
    
    # === 4.5. Determine DE input source ===
    if celltype_idx == 1:
        print("\n[INFO] Determining DE input source...")
        
        if DE_MODE == "auto":
            # Check data characteristics
            has_raw = adata_ct.raw is not None
            has_log1p_layer = 'log1p' in adata_ct.layers
            
            # P2: Use .data for sparse matrices
            if sparse.issparse(adata_ct.X):
                x_min = float(adata_ct.X.data.min()) if adata_ct.X.data.size > 0 else 0
                x_max = float(adata_ct.X.data.max()) if adata_ct.X.data.size > 0 else 0
            else:
                x_min = float(adata_ct.X.min())
                x_max = float(adata_ct.X.max())
            
            x_is_log_like = (x_min >= 0 and x_max < 20)
            
            print(f"  Data characteristics:")
            print(f"    Has .raw: {has_raw}")
            print(f"    Has log1p layer: {has_log1p_layer}")
            print(f"    .X range: [{x_min:.2f}, {x_max:.2f}]")
            print(f"    .X appears log-transformed: {x_is_log_like}")
            
            # Decision logic
            if has_log1p_layer:
                DE_USE_RAW = False
                DE_USE_LAYER = 'log1p'
                de_note = "Using log1p layer (explicit)"
            elif has_raw and x_is_log_like:
                if gene_filter_dict_raw is not None:
                    DE_USE_RAW = True
                    DE_USE_LAYER = None
                    de_note = "Using .raw (assumes log1p, filter_dict matches)"
                else:
                    DE_USE_RAW = False
                    DE_USE_LAYER = None
                    de_note = "Using .X (raw filter not built, .X looks log-like)"
            elif x_is_log_like:
                DE_USE_RAW = False
                DE_USE_LAYER = None
                de_note = "Using .X (appears log-transformed)"
            else:
                print("  [WARN] .X does not appear log-transformed!")
                DE_USE_RAW = False
                DE_USE_LAYER = None
                de_note = "Using .X (WARNING: may not be log1p)"
        
        elif DE_MODE == "raw":
            DE_USE_RAW = True
            DE_USE_LAYER = None
            de_note = "Using .raw (user-specified)"
        
        elif DE_MODE == "X":
            DE_USE_RAW = False
            DE_USE_LAYER = None
            de_note = "Using .X (user-specified)"
        
        elif DE_MODE == "layer":
            DE_USE_RAW = False
            DE_USE_LAYER = DE_LAYER
            de_note = f"Using layer '{DE_LAYER}' (user-specified)"
        
        print(f"\n  [OK] DE decision: {de_note}")
    
    # === 4.6. Marker gene analysis ===
    print(f"\n[INFO] Computing marker genes for {celltype} subclusters...")
    
    # P1: Filter clusters by MIN_CELLS_FOR_MARKER
    valid_clusters = cluster_counts_final[
        cluster_counts_final >= MIN_CELLS_FOR_MARKER
    ].index.tolist()
    
    if len(valid_clusters) == 0:
        print(f"  [WARN] No clusters have >= {MIN_CELLS_FOR_MARKER} cells")
        print("  Skipping marker analysis for this cell type")
    else:
        print(f"  [INFO] Analyzing {len(valid_clusters)}/{n_clusters_final} clusters")
        print(f"  (Skipping {n_clusters_final - len(valid_clusters)} clusters with < {MIN_CELLS_FOR_MARKER} cells)")
        
        skipped_clusters = cluster_counts_final[
            cluster_counts_final < MIN_CELLS_FOR_MARKER
        ].index.tolist()
        if len(skipped_clusters) > 0:
            print(f"  Skipped clusters: {skipped_clusters}")
        
        try:
            # P1: Only analyze valid clusters
            sc.tl.rank_genes_groups(
                adata_ct,
                groupby='subcluster_leiden',
                groups=valid_clusters,  # P1: Only valid clusters
                method=MARKER_METHOD,
                use_raw=DE_USE_RAW,
                layer=DE_USE_LAYER,
                corr_method='benjamini-hochberg',
                pts=True,
                tie_correct=True,
                n_genes=min(TOP_DE_GENES, adata_ct.n_vars),
                key_added='rank_genes_groups'
            )
            
            print(f"  [OK] Marker analysis complete")
            print(f"    Method: {MARKER_METHOD}")
            print(f"    Input: use_raw={DE_USE_RAW}, layer={DE_USE_LAYER}")
            print(f"    Genes analyzed: {min(TOP_DE_GENES, adata_ct.n_vars)}")
            
            # Extract and filter markers
            print("\n  [INFO] Extracting marker genes...")
            
            # Get which filter dict to use
            if DE_USE_RAW and gene_filter_dict_raw is not None:
                active_filter_dict = gene_filter_dict_raw
                filter_source_used = "raw"
            else:
                active_filter_dict = gene_filter_dict_var if gene_filter_dict_var is not None else {}
                filter_source_used = "var"
            
            print(f"    Using filter dictionary: {filter_source_used} ({len(active_filter_dict):,} genes)")
            
            # Extract markers for each cluster
            result = adata_ct.uns['rank_genes_groups']
            
            markers_unfiltered = []
            markers_filtered = []
            
            for cluster in valid_clusters:  # P1: Only valid clusters
                # Get all results for this cluster
                genes = result['names'][cluster]
                pvals = result['pvals'][cluster]
                pvals_adj = result['pvals_adj'][cluster]
                logfcs = result['logfoldchanges'][cluster]
                
                # Get pct values if available
                if 'pts' in result:
                    pts = result['pts'][cluster]
                    pts_rest = result['pts_rest'][cluster]
                else:
                    pts = np.ones_like(pvals)
                    pts_rest = np.zeros_like(pvals)
                
                # Create DataFrame
                cluster_df = pd.DataFrame({
                    'cluster': cluster,
                    'names': genes,
                    'pvals': pvals,
                    'pvals_adj': pvals_adj,
                    'logfoldchanges': logfcs,
                    'pct_in_group': pts,
                    'pct_out_group': pts_rest,
                    'delta_pct': pts - pts_rest
                })
                
                # Basic statistical filter
                cluster_df_stat = cluster_df[
                    (cluster_df['pvals_adj'] <= MARKER_PADJ_THRESHOLD) &
                    (cluster_df['logfoldchanges'] >= MARKER_LOGFC_THRESHOLD) &
                    (cluster_df['pct_in_group'] >= MARKER_MIN_PCT) &
                    (cluster_df['delta_pct'] >= MARKER_MIN_DELTA_PCT)
                ].copy()
                
                markers_unfiltered.append(cluster_df_stat)
                
                # Gene-level filtering
                if active_filter_dict:
                    cluster_df_filtered = cluster_df_stat[
                        ~cluster_df_stat['names'].apply(
                            lambda g: any(active_filter_dict.get(g, {}).values()) 
                            if g in active_filter_dict else False
                        )
                    ].copy()
                    
                    markers_filtered.append(cluster_df_filtered)
                else:
                    markers_filtered.append(cluster_df_stat)
            
            # Combine
            if len(markers_unfiltered) > 0:
                all_markers[celltype] = pd.concat(markers_unfiltered, ignore_index=True)
                print(f"    Unfiltered markers: {len(all_markers[celltype]):,} genes")
            
            if len(markers_filtered) > 0:
                all_markers_filtered[celltype] = pd.concat(markers_filtered, ignore_index=True)
                print(f"    Filtered markers: {len(all_markers_filtered[celltype]):,} genes")
                
                # Calculate filtering stats
                if celltype in all_markers:
                    n_before = len(all_markers[celltype])
                    n_after = len(all_markers_filtered[celltype])
                    n_removed = n_before - n_after
                    pct_removed = (n_removed / n_before * 100) if n_before > 0 else 0
                    
                    marker_filter_stats[celltype] = {
                        'before': n_before,
                        'after': n_after,
                        'removed': n_removed,
                        'pct_removed': pct_removed
                    }
                    
                    print(f"    Filtering removed: {n_removed} genes ({pct_removed:.1f}%)")
            
            # P2: Save per-celltype immediately
            if celltype in all_markers:
                ct_safe = celltype.replace(' ', '_').replace('/', '_')
                
                # Full markers
                markers_all_path = OUTPUT_DIR / f"markers_{ct_safe}_all.csv"
                all_markers[celltype].to_csv(markers_all_path, index=False)
                print(f"  [OK] Saved: {markers_all_path.name}")
                
                # P1: Top N markers
                if celltype in all_markers_filtered:
                    top_markers_list = []
                    for cluster in valid_clusters:
                        cluster_markers = all_markers_filtered[celltype][
                            all_markers_filtered[celltype]['cluster'] == cluster
                        ]
                        top_n = cluster_markers.nsmallest(TOP_N_MARKERS, 'pvals_adj')
                        top_markers_list.append(top_n)
                    
                    if len(top_markers_list) > 0:
                        top_markers = pd.concat(top_markers_list, ignore_index=True)
                        markers_top_path = OUTPUT_DIR / f"markers_{ct_safe}_top{TOP_N_MARKERS}.csv"
                        top_markers.to_csv(markers_top_path, index=False)
                        print(f"  [OK] Saved: {markers_top_path.name} ({len(top_markers)} genes)")
            
        except Exception as e:
            print(f"  [ERROR] Marker gene analysis failed: {e}")
            import traceback
            traceback.print_exc()
    
    # === 4.7. Quick visualization ===
    print(f"\n[INFO] Generating visualizations for {celltype}...")
    
    fig_dir_ct = sc.settings.figdir / celltype.replace(' ', '_').replace('/', '_')
    fig_dir_ct.mkdir(exist_ok=True)
    
    # UMAP of subclusters
    if 'X_umap' in adata_ct.obsm:
        fig, ax = plt.subplots(figsize=(10, 8))
        sc.pl.umap(
            adata_ct,
            color='subcluster_leiden',
            title=f'{celltype} Subclusters',
            legend_loc='right margin',
            ax=ax,
            show=False
        )
        plt.tight_layout()
        plt.savefig(fig_dir_ct / f'umap_subclusters.{FIGURE_FORMAT}', 
                   dpi=FIGURE_DPI, bbox_inches='tight')
        plt.close()
        print(f"  [OK] Saved: {celltype}/umap_subclusters.{FIGURE_FORMAT}")
    
    # Cluster size bar plot
    fig, ax = plt.subplots(figsize=(10, 6))
    cluster_counts_final.plot(kind='bar', ax=ax, color='steelblue')
    ax.axhline(MIN_CELLS_PER_CLUSTER, color='red', linestyle='--', 
               label=f'Min threshold ({MIN_CELLS_PER_CLUSTER})')
    ax.axhline(MIN_CELLS_FOR_MARKER, color='orange', linestyle='--',
               label=f'Marker threshold ({MIN_CELLS_FOR_MARKER})')
    ax.set_xlabel('Subcluster')
    ax.set_ylabel('Number of Cells')
    ax.set_title(f'{celltype} - Cluster Size Distribution')
    ax.legend()
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(fig_dir_ct / f'cluster_sizes.{FIGURE_FORMAT}', 
               dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close()
    print(f"  [OK] Saved: {celltype}/cluster_sizes.{FIGURE_FORMAT}")
    
    # Store results
    processed_celltypes.append(celltype)
    
    # P0: Write back using obs_names alignment (handles remove strategy)
    print(f"\n[INFO] Writing results back to main adata...")
    
    # Write processed cells
    adata.obs.loc[adata_ct.obs_names, 'subcluster'] = adata_ct.obs['subcluster'].astype(str).values
    adata.obs.loc[adata_ct.obs_names, 'subcluster_leiden'] = adata_ct.obs['subcluster_leiden'].astype(str).values
    
    # P0: Mark removed cells (if remove strategy was used)
    if SMALL_CLUSTER_STRATEGY == "remove" and len(removed_cell_ids) > 0:
        adata.obs.loc[removed_cell_ids, 'subcluster'] = f"{celltype}_REMOVED"
        adata.obs.loc[removed_cell_ids, 'subcluster_leiden'] = "REMOVED"
        print(f"  [INFO] Marked {len(removed_cell_ids)} removed cells as 'REMOVED'")
    
    print(f"  [OK] Write-back complete")
    
    ct_elapsed = time.time() - ct_start_time
    print(f"\n[OK] {celltype} processing complete ({ct_elapsed/60:.1f} min)")
    
    # Clean up
    del adata_ct
    gc.collect()

print("\n" + "=" * 80)
print("ALL CELL TYPES PROCESSED")
print("=" * 80)
print(f"Processed: {len(processed_celltypes)} cell types")
print(f"Total subclusters: {adata.obs['subcluster'].nunique()}")

# %% [markdown]
# ## Step 8: Save Results

# %%
print("\n" + "=" * 80)
print("STEP 8: Saving Results")
print("=" * 80)

# Main h5ad
output_h5ad = OUTPUT_DIR / "epithelial_with_subclusters_v4_5_hotfix.h5ad"
print(f"[INFO] Writing h5ad: {output_h5ad}")
adata.write_h5ad(output_h5ad, compression="gzip", compression_opts=9)
file_size_gb = output_h5ad.stat().st_size / 1e9
print(f"[OK] Saved: {output_h5ad.name} ({file_size_gb:.2f} GB)")

# Subcluster annotations
subcluster_df = adata.obs[[CELLTYPE_COLUMN, "subcluster", "subcluster_leiden"]].copy()
if BATCH_KEY and BATCH_KEY in adata.obs.columns:
    subcluster_df[BATCH_KEY] = adata.obs[BATCH_KEY]
subcluster_csv = OUTPUT_DIR / "subcluster_annotations.csv"
subcluster_df.to_csv(subcluster_csv)
print(f"[OK] Saved: {subcluster_csv.name}")

# P1: Combined markers (both full and top-N)
if len(all_markers) > 0:
    all_markers_df = pd.concat(all_markers.values(), ignore_index=True)
    markers_all_csv = OUTPUT_DIR / "markers_all_celltypes_combined.csv"
    all_markers_df.to_csv(markers_all_csv, index=False)
    print(f"[OK] Saved: {markers_all_csv.name} ({len(all_markers_df):,} rows)")

if len(all_markers_filtered) > 0:
    all_markers_filtered_df = pd.concat(all_markers_filtered.values(), ignore_index=True)
    markers_filtered_csv = OUTPUT_DIR / "markers_all_celltypes_filtered.csv"
    all_markers_filtered_df.to_csv(markers_filtered_csv, index=False)
    print(f"[OK] Saved: {markers_filtered_csv.name} ({len(all_markers_filtered_df):,} rows)")
    
    # P1: Top N per cluster (across all celltypes)
    top_list = []
    for celltype in processed_celltypes:
        if celltype in all_markers_filtered:
            ct_df = all_markers_filtered[celltype]
            for cluster in ct_df['cluster'].unique():
                cluster_df = ct_df[ct_df['cluster'] == cluster]
                top_list.append(cluster_df.nsmallest(TOP_N_MARKERS, 'pvals_adj'))
    
    if len(top_list) > 0:
        top_df = pd.concat(top_list, ignore_index=True)
        top_csv = OUTPUT_DIR / f"markers_all_celltypes_top{TOP_N_MARKERS}.csv"
        top_df.to_csv(top_csv, index=False)
        print(f"[OK] Saved: {top_csv.name} ({len(top_df):,} rows)")

# Summary report
summary_file = OUTPUT_DIR / "subcluster_summary.txt"
with open(summary_file, "w") as f:
    f.write("=" * 80 + "\n")
    f.write(f"Epithelial Subcluster Analysis Summary ({PIPELINE_VERSION})\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Input: {INPUT_H5AD}\n")
    f.write(f"Strategy: {SMALL_CLUSTER_STRATEGY}\n")
    f.write(f"Min cluster size: {MIN_CELLS_PER_CLUSTER}\n")
    f.write(f"Min for markers: {MIN_CELLS_FOR_MARKER}\n")
    f.write(f"Random seed: {RANDOM_SEED}\n\n")
    
    f.write("Subclustering:\n")
    f.write(f"  Resolution: FIXED {RESOLUTION_FIXED}\n")
    f.write(f"  Total subclusters: {adata.obs['subcluster'].nunique()}\n")
    f.write(f"  Processed cell types: {len(processed_celltypes)}\n\n")
    
    f.write("Cluster size distribution:\n")
    for sub, cnt in adata.obs["subcluster"].value_counts().items():
        f.write(f"  {sub}: {cnt}\n")
    
    if SMALL_CLUSTER_STRATEGY == "remove":
        n_removed = (adata.obs['subcluster_leiden'] == 'REMOVED').sum()
        if n_removed > 0:
            f.write(f"\nRemoved cells: {n_removed}\n")

print(f"[OK] Saved: {summary_file.name}")

print("\n" + "=" * 80)
print(f"PIPELINE COMPLETE ({PIPELINE_VERSION})")
print("=" * 80)
print(f"Output directory: {OUTPUT_DIR}")
print(f"Strategy used: {SMALL_CLUSTER_STRATEGY}")
print(f"Total subclusters: {adata.obs['subcluster'].nunique()}")
print(f"\nKey outputs:")
print(f"  - epithelial_with_subclusters_v4_5_hotfix.h5ad")
print(f"  - subcluster_annotations.csv")
print(f"  - markers_all_celltypes_combined.csv")
print(f"  - markers_all_celltypes_filtered.csv")
print(f"  - markers_all_celltypes_top{TOP_N_MARKERS}.csv")
print(f"  - markers_[celltype]_all.csv (per celltype)")
print(f"  - markers_[celltype]_top{TOP_N_MARKERS}.csv (per celltype)")

# %%
