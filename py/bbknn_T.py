#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tcell_bbknn_analysis_optimized.py

T cell BBKNN batch integration and analysis pipeline - Unified optimized version
Features:
1. Complete BBKNN batch correction workflow with post-hoc parameter validation
2. Force retention of specified markers (avoid HVG/filtering loss)
3. Differential analysis/plotting reads from .raw full gene matrix (use_raw=True)
4. Multi-resolution clustering; UMAP uses BBKNN neighbor graph
5. Detailed visualization output
6. Consistent data structure and workflow with epithelial cell analysis

Author: Clinical-Bioinformatics Team
Version: v1.1 (unified, fixed)
"""

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
from tqdm import tqdm
import pickle
warnings.filterwarnings('ignore')

# ==================== Configuration Section ====================

# ========== Input/Output Configuration ==========
INPUT_H5AD_PATH = "/home/h2048/data/py/1029/bbknn_celltype_analysis/T/adata_T_bbknn.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1029/bbknn_celltype_analysis/T/output_optimized"
OVERWRITE_EXISTING = True

# ========== BBKNN Integration Configuration ==========
BATCH_KEY = "dataset"                 # Batch variable
BBKNN_NEIGHBORS_WITHIN_BATCH = 4      # Neighbors per batch (T cells tuned lower to increase cluster separation)
BBKNN_N_PCS = 50                      # Number of PCs for BBKNN
BBKNN_METRIC = "correlation"          # Distance metric (similar to cosine, suitable for T cells)
BBKNN_TRIM = 35                       # trim parameter: int or "auto"

# ========== Highly Variable Gene Selection ==========
USE_HVG = True
N_TOP_GENES = 4000
HVG_FLAVOR = "seurat_v3"

# ========== Gene Filtering Configuration ==========
MIN_CELLS_PER_GENE = 3

# ========== Raw Counts Source Configuration ==========
# "auto"|"raw"|"layers"|"X"
RAW_COUNTS_SOURCE = "auto"

# ========== Normalization Configuration ==========
NORMALIZE_TOTAL = True
TARGET_SUM = 1e4
LOG_TRANSFORM = True
SCALE_DATA = True
MAX_VALUE = 10

# ========== PCA Configuration ==========
N_PCS = 50

# ========== Dimensionality Reduction and Visualization Configuration ==========
RUN_UMAP = True
UMAP_MIN_DIST = 0.15  # T cells use smaller value to enhance local structure

# ========== Multi-resolution Clustering Configuration ==========
RUN_CLUSTERING = True
LEIDEN_RESOLUTIONS = [1.2, 1.6, 2.0, 2.4]  # T cell subtypes are numerous, add high resolution
DEFAULT_RESOLUTION = 2.0

# ========== T Cell Marker Genes ==========
MARKER_TCELL = [
    # NK cells
    'KLRD1', 'FCGR3A', 'GNLY', 'TYROBP', 'FCER1G', 'KLRC1', 'FGFBP2', 'SPON2', 'MYOM2',
    # Gamma-delta T
    'TRDC', 'TRGC1', 'TRGC2',
    # Th2
    'GATA3', 'IL5', 'AREG', 'HPGDS', 'IL4', 'IL13',
    # Th17
    'IL23R', 'RORC', 'IL17A', 'CCL20', 'LST1', 'PCDH9', 'TNFSF11', 'KRT86',
    # CD4 T cells
    'CD4', 'CD28', 'CD40LG', 'TRAT1', 'TNFRSF25',
    # CD8 T cells
    'CD8A', 'CD8B',
    # T cell core markers
    'CD3E', 'CD2', 'TRBC2',
    # Naive/Central memory
    'CCR7', 'TCF7', 'LEF1', 'SELL', 'IL7R',
    # Th1
    'TBX21', 'IFNG',
    # Effector/Memory
    'EOMES', 'GZMK', 'KLRG1', 'ITGA1',
    # Cytotoxic
    'GZMB', 'GZMH', 'ZNF683', 'CCL4L2',
    # Exhaustion
    'PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'CTLA4',
    # Regulatory T
    'FOXP3', 'IL2RA', 'IKZF2', 'TNFRSF4',
    # MAIT
    'KLRB1', 'NCR3', 'CEBPD', 'SLC4A10', 'TRAV1-2',
    # Proliferation
    'MKI67', 'TOP2A', 'TK1', 'CENPW',
    # Activation
    'CD69', 'ITGAE', 'CXCR6', 'CD25', 'CD38', 'HLA-DR'
]

KEY_MARKERS = {
    'CD4_Naive': ['CD4', 'CCR7', 'TCF7', 'SELL'],
    'CD4_Memory': ['CD4', 'IL7R', 'CD28'],
    'CD8_Naive': ['CD8A', 'CCR7', 'TCF7', 'SELL'],
    'CD8_Memory': ['CD8A', 'IL7R'],
    'CD8_Effector': ['CD8A', 'GZMB', 'GZMH', 'IFNG'],
    'Treg': ['FOXP3', 'IL2RA', 'IKZF2'],
    'Th1': ['CD4', 'TBX21', 'IFNG'],
    'Th2': ['CD4', 'GATA3', 'IL4', 'IL13'],
    'Th17': ['CD4', 'RORC', 'IL17A'],
    'NK': ['KLRD1', 'FCGR3A', 'GNLY'],
    'NKT': ['CD3E', 'KLRD1', 'FCGR3A'],
    'gdT': ['TRDC', 'TRGC2'],
    'MAIT': ['KLRB1', 'NCR3', 'TRAV1-2'],
    'Proliferating': ['MKI67', 'TOP2A']
}

# ========== Marker Gene Analysis Configuration ==========
RUN_FIND_MARKERS = True
MARKER_MIN_PCT = 0.25
MARKER_LOGFC_THRESHOLD = 0.25
TOP_N_MARKERS = 10

# ========== Performance Optimization Configuration ==========
USE_CACHE = True
CHUNK_SIZE = 100

# ========== Visualization Configuration ==========
GENERATE_DOTPLOT = True
GENERATE_HEATMAP = True
GENERATE_FACET_PLOTS = True
DPI = 300
FIGURE_FORMAT = "png"

VERBOSE = True

# ========= Force retention (not filtered even if not in HVG/not meeting min_cells) =========
FORCE_INCLUDE_MARKERS = sorted({
    *MARKER_TCELL,
    *[g for genes in KEY_MARKERS.values() for g in genes],
})

# ==================== Helper Functions ====================

def log_msg(msg):
    if VERBOSE:
        print(msg)

def log_step(step_num, step_name):
    log_msg("\n" + "="*70)
    log_msg(f"Step {step_num}: {step_name}")
    log_msg("="*70)

def save_checkpoint(data, filename, output_dir):
    if USE_CACHE:
        checkpoint_path = Path(output_dir) / f".cache_{filename}"
        with open(checkpoint_path, 'wb') as f:
            pickle.dump(data, f)
        log_msg(f"   Checkpoint saved: {checkpoint_path}")

def load_checkpoint(filename, output_dir):
    if USE_CACHE:
        checkpoint_path = Path(output_dir) / f".cache_{filename}"
        if checkpoint_path.exists():
            with open(checkpoint_path, 'rb') as f:
                return pickle.load(f)
    return None

# ==================== Main Functionality Functions ====================

def load_anndata(h5ad_path):
    import scanpy as sc
    log_step(1, "Loading Data")
    log_msg(f"\nReading file: {h5ad_path}")
    if not Path(h5ad_path).exists():
        raise FileNotFoundError(f"File not found: {h5ad_path}")

    adata = sc.read_h5ad(h5ad_path)
    log_msg(f"Data loaded successfully")
    log_msg(f"   Cells: {adata.n_obs:,}")
    log_msg(f"   Genes: {adata.n_vars:,}")

    log_msg(f"\nAvailable metadata columns:")
    for col in adata.obs.columns:
        n_unique = adata.obs[col].nunique()
        log_msg(f"   - {col}: {n_unique} unique values")

    if BATCH_KEY not in adata.obs.columns:
        raise ValueError(f"Batch key '{BATCH_KEY}' not found in adata.obs")

    log_msg(f"\nBatch distribution (key: {BATCH_KEY}):")
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        log_msg(f"   {batch}: {count:,} cells ({pct:.1f}%)")
    return adata

def check_marker_genes(adata):
    """Check marker availability (check both var and raw)"""
    log_step(2, "Checking Marker Genes")
    raw_names = set(adata.raw.var_names) if adata.raw is not None else set()
    var_names = set(adata.var_names)

    def present(g):
        return (g in var_names) or (g in raw_names)

    available_markers = [g for g in MARKER_TCELL if present(g)]
    missing_markers   = [g for g in MARKER_TCELL if not present(g)]

    log_msg(f"\nMarker gene availability:")
    log_msg(f"   Total markers: {len(MARKER_TCELL)}")
    log_msg(f"   Available: {len(available_markers)} ({len(available_markers)/len(MARKER_TCELL)*100:.1f}%)")
    log_msg(f"   Missing: {len(missing_markers)} ({len(missing_markers)/len(MARKER_TCELL)*100:.1f}%)")

    if missing_markers:
        log_msg(f"\n   Missing markers: {', '.join(missing_markers[:10])}")
        if len(missing_markers) > 10:
            log_msg(f"   ... and {len(missing_markers)-10} more")

    log_msg(f"\nKey T cell type markers:")
    for celltype, genes in KEY_MARKERS.items():
        available = [g for g in genes if present(g)]
        log_msg(f"   {celltype}: {len(available)}/{len(genes)} available")
    return available_markers

def preprocess_for_bbknn(adata):
    """
    Preprocess data for BBKNN integration
    Key points:
    1. Auto-detect if data is already log-transformed
    2. Force retention of marker genes (even if not meeting min_cells)
    3. Return both full (all genes) and hvg (highly variable genes) versions
    """
    import scanpy as sc
    log_step(3, "Preprocessing for BBKNN")

    # Auto-identify raw count source
    x_data = adata.X
    if issparse(x_data):
        sample = x_data[:1000, :100].toarray()
    else:
        sample = x_data[:1000, :100]
    is_log = np.any((sample > 0) & (sample < 1))

    if RAW_COUNTS_SOURCE == "auto":
        if adata.raw is not None:
            raw = adata.raw.X
            src = ".raw"
        elif 'counts' in adata.layers:
            raw = adata.layers['counts']
            src = "layers['counts']"
        elif not is_log:
            raw = adata.X
            src = ".X"
        else:
            raise ValueError("No raw counts found; set RAW_COUNTS_SOURCE")
    elif RAW_COUNTS_SOURCE == "raw":
        raw = adata.raw.X
        src = ".raw"
    elif RAW_COUNTS_SOURCE == "layers":
        raw = adata.layers['counts']
        src = "layers['counts']"
    elif RAW_COUNTS_SOURCE == "X":
        raw = adata.X
        src = ".X"
    
    log_msg(f"\n   Using {src} as raw counts")
    log_msg(f"   Data type: {'sparse' if issparse(raw) else 'dense'}")

    # Create working copy
    adata_work = adata.copy()
    adata_work.X = raw.copy()
    
    # Convert to appropriate data type
    if issparse(adata_work.X):
        adata_work.X = adata_work.X.tocsr().astype(np.float32)
    else:
        adata_work.X = adata_work.X.astype(np.float32, copy=False)

    log_msg(f"\n   Initial gene count: {adata_work.n_vars:,}")

    # Gene filtering + force retention of markers
    log_msg(f"\n   Filtering genes (min_cells >= {MIN_CELLS_PER_GENE})...")
    if issparse(adata_work.X):
        n_cells_per_gene = np.array((adata_work.X > 0).sum(axis=0)).ravel()
    else:
        n_cells_per_gene = (adata_work.X > 0).sum(axis=0)
    
    keep = n_cells_per_gene >= MIN_CELLS_PER_GENE
    
    # Force retention of markers
    n_protected = 0
    for g in FORCE_INCLUDE_MARKERS:
        if g in adata_work.var_names:
            idx = adata_work.var_names.get_loc(g)
            if not keep[idx]:
                keep[idx] = True
                n_protected += 1
    
    adata_work = adata_work[:, keep].copy()
    log_msg(f"   Genes after filtering: {adata_work.n_vars:,}")
    log_msg(f"   Protected markers: {n_protected}")

    # Normalization
    log_msg(f"\n   Normalizing...")
    if NORMALIZE_TOTAL:
        sc.pp.normalize_total(adata_work, target_sum=TARGET_SUM)
    
    if LOG_TRANSFORM:
        sc.pp.log1p(adata_work)
    
    log_msg(f"   ✓ Normalization and log1p completed")

    # HVG selection (batch-aware)
    if USE_HVG:
        log_msg(f"\n   Selecting highly variable genes...")
        log_msg(f"   Method: {HVG_FLAVOR}")
        log_msg(f"   n_top_genes: {N_TOP_GENES}")
        log_msg(f"   Batch-aware using: {BATCH_KEY}")
        
        sc.pp.highly_variable_genes(
            adata_work,
            n_top_genes=N_TOP_GENES,
            flavor=HVG_FLAVOR,
            batch_key=BATCH_KEY
        )
        
        # Force markers to be marked as highly variable
        n_forced_hvg = 0
        for g in FORCE_INCLUDE_MARKERS:
            if g in adata_work.var_names:
                if not adata_work.var.loc[g, 'highly_variable']:
                    adata_work.var.loc[g, 'highly_variable'] = True
                    n_forced_hvg += 1
        
        n_hvg = adata_work.var['highly_variable'].sum()
        log_msg(f"   HVGs selected: {n_hvg}")
        log_msg(f"   Markers forced into HVGs: {n_forced_hvg}")
        
        adata_hvg = adata_work[:, adata_work.var['highly_variable']].copy()
    else:
        adata_hvg = adata_work.copy()
        log_msg(f"\n   Skipping HVG selection (using all genes)")

    # Scale data
    if SCALE_DATA:
        log_msg(f"\n   Scaling data...")
        log_msg(f"   max_value: {MAX_VALUE}")
        sc.pp.scale(adata_hvg, max_value=MAX_VALUE)
        log_msg(f"   ✓ Scaling completed")

    # Save .raw copy (full gene set, already normalized+log1p)
    log_msg(f"\n   Storing raw copy for downstream analysis...")
    adata_hvg.raw = adata_work.copy()
    
    # Keep full version for final saving
    adata_full = adata_work.copy()
    del adata_work  # Free memory
    
    log_msg(f"\n   Preprocessing summary:")
    log_msg(f"   - adata_hvg shape: {adata_hvg.shape}")
    log_msg(f"   - adata_hvg.raw shape: {adata_hvg.raw.shape}")
    log_msg(f"   - adata_full shape: {adata_full.shape}")

    return adata_full, adata_hvg

def run_pca(adata_hvg):
    """Perform PCA dimensionality reduction"""
    import scanpy as sc
    log_step(3.5, "PCA Analysis")
    
    log_msg(f"\n   Computing PCA...")
    log_msg(f"   n_comps: {N_PCS}")
    
    # Choose solver based on data sparsity
    svd_solver = 'arpack' if issparse(adata_hvg.X) else 'randomized'
    log_msg(f"   svd_solver: {svd_solver}")
    
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver=svd_solver)
    
    # Output variance explained ratio
    var_ratio = adata_hvg.uns['pca']['variance_ratio']
    cumsum_ratio = np.cumsum(var_ratio)
    log_msg(f"\n   PCA variance explained:")
    log_msg(f"   - PC1-10: {cumsum_ratio[9]:.1%}")
    log_msg(f"   - PC1-20: {cumsum_ratio[19]:.1%}")
    log_msg(f"   - PC1-50: {cumsum_ratio[49]:.1%}")
    
    return adata_hvg

def run_bbknn_integration(adata_hvg):
    """Perform BBKNN batch integration"""
    from bbknn import bbknn as bbknn_func
    log_step(4, "BBKNN Integration")
    
    n_batches = adata_hvg.obs[BATCH_KEY].nunique()
    log_msg(f"\n   Number of batches: {n_batches}")
    
    # Auto-calculate trim parameter
    if BBKNN_TRIM == "auto":
        rec_trim = min(20, max(10, int(0.6 * BBKNN_NEIGHBORS_WITHIN_BATCH * n_batches)))
        log_msg(f"   Auto-calculated trim: {rec_trim}")
    else:
        rec_trim = int(BBKNN_TRIM)
        log_msg(f"   User-specified trim: {rec_trim}")
    
    log_msg(f"\n   BBKNN parameters:")
    log_msg(f"   - neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    log_msg(f"   - n_pcs: {BBKNN_N_PCS}")
    log_msg(f"   - metric: {BBKNN_METRIC}")
    log_msg(f"   - trim: {rec_trim}")
    log_msg(f"   - batch_key: {BATCH_KEY}")
    
    log_msg(f"\n   Running BBKNN...")
    bbknn_func(
        adata_hvg,
        batch_key=BATCH_KEY,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
        n_pcs=BBKNN_N_PCS,
        metric=BBKNN_METRIC,
        trim=rec_trim,  # FIXED: Previously commented out
        copy=False
    )
    
    log_msg(f"   ✓ BBKNN integration completed")
    
    # Validate neighbor graph
    if 'connectivities' in adata_hvg.obsp:
        n_neighbors = (adata_hvg.obsp['connectivities'] > 0).sum(axis=1).mean()
        log_msg(f"   Average neighbors per cell: {n_neighbors:.1f}")
    
    return adata_hvg

def run_umap_and_clustering(adata_hvg):
    """Perform UMAP dimensionality reduction and multi-resolution Leiden clustering"""
    import scanpy as sc
    log_step(5, "UMAP and Clustering")
    
    # UMAP
    if RUN_UMAP:
        log_msg(f"\n   Computing UMAP...")
        log_msg(f"   min_dist: {UMAP_MIN_DIST}")
        sc.tl.umap(adata_hvg, min_dist=UMAP_MIN_DIST)
        log_msg(f"   ✓ UMAP completed")
    
    # Multi-resolution clustering
    if RUN_CLUSTERING:
        log_msg(f"\n   Running multi-resolution Leiden clustering...")
        log_msg(f"   Resolutions: {LEIDEN_RESOLUTIONS}")
        
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            sc.tl.leiden(adata_hvg, resolution=res, key_added=key)
            n_clusters = adata_hvg.obs[key].nunique()
            log_msg(f"   Resolution {res}: {n_clusters} clusters")
        
        # Set default resolution
        default_key = f'leiden_bbknn_res{DEFAULT_RESOLUTION}'
        adata_hvg.obs['leiden_bbknn'] = adata_hvg.obs[default_key].astype(str)
        log_msg(f"\n   Default resolution: {DEFAULT_RESOLUTION}")
        log_msg(f"   Default clustering key: 'leiden_bbknn'")
    
    return adata_hvg

def find_all_markers_optimized(adata, cluster_key, output_dir, min_pct=0.25, logfc_threshold=0.25):
    """
    Find cluster marker genes using Wilcoxon rank-sum test
    Read full gene set from .raw for differential analysis
    """
    import scanpy as sc
    log_step(6, "Finding Cluster Markers")
    
    use_raw_flag = adata.raw is not None
    n_genes_for_de = adata.raw.n_vars if use_raw_flag else adata.n_vars
    
    log_msg(f"\n   Cluster key: {cluster_key}")
    log_msg(f"   Using .raw: {use_raw_flag}")
    log_msg(f"   Gene set size: {n_genes_for_de:,}")
    log_msg(f"   Method: wilcoxon")
    log_msg(f"   logFC threshold: {logfc_threshold}")
    
    log_msg(f"\n   Running differential expression...")
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method='wilcoxon',
        use_raw=use_raw_flag,
        n_genes=6000,
        tie_correct=True
    )
    
    # Extract results
    log_msg(f"\n   Extracting significant markers...")
    result = {}
    n_clusters = adata.obs[cluster_key].nunique()
    
    for c in sorted(adata.obs[cluster_key].unique(), key=str):
        names = adata.uns['rank_genes_groups']['names'][c]
        scores = adata.uns['rank_genes_groups']['scores'][c]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][c]
        lfc = adata.uns['rank_genes_groups']['logfoldchanges'][c]
        
        df = pd.DataFrame({
            'gene': names,
            'scores': scores,
            'pvals_adj': pvals_adj,
            'logfoldchanges': lfc
        })
        
        # Filter conditions
        df = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > logfc_threshold)]
        df['cluster'] = str(c)
        result[str(c)] = df
        
        log_msg(f"   Cluster {c}: {len(df)} significant markers")
    
    # Combine all results
    all_markers = pd.concat(result.values(), ignore_index=True)
    
    # Save
    out_path = Path(output_dir) / "cluster_markers.csv"
    all_markers.to_csv(out_path, index=False)
    
    log_msg(f"\n   ✓ Markers saved: {out_path}")
    log_msg(f"   Total significant markers: {len(all_markers):,}")
    
    return all_markers

def generate_visualizations(adata, available_markers, output_dir):
    """Generate visualization plots"""
    import scanpy as sc
    log_step(7, "Generating Visualizations")
    
    fig_dir = Path(output_dir) / "figures"
    fig_dir.mkdir(exist_ok=True)
    
    sc.settings.figdir = fig_dir
    sc.set_figure_params(dpi=DPI, frameon=False, figsize=(7, 7))
    
    log_msg(f"\n   Output directory: {fig_dir}")
    
    # 1) UMAP by clusters
    log_msg("\n   Generating cluster UMAP...")
    try:
        sc.pl.umap(
            adata,
            color='leiden_bbknn',
            legend_loc='on data',
            legend_fontsize=8,
            title='T Cell Clusters (BBKNN)',
            save=f'_clusters.{FIGURE_FORMAT}',
            show=False
        )
        log_msg(f"   ✓ Cluster UMAP saved")
    except Exception as e:
        log_msg(f"   ⚠️ Cluster UMAP failed: {e}")
    
    # 2) UMAP by batch
    if BATCH_KEY in adata.obs.columns:
        log_msg("\n   Generating batch UMAP...")
        try:
            sc.pl.umap(
                adata,
                color=BATCH_KEY,
                title='Batch Distribution',
                save=f'_batch.{FIGURE_FORMAT}',
                show=False
            )
            log_msg(f"   ✓ Batch UMAP saved")
        except Exception as e:
            log_msg(f"   ⚠️ Batch UMAP failed: {e}")
    
    # 3) Marker gene expression
    log_msg("\n   Generating marker expression UMAPs...")
    try:
        # Use first 20 available markers
        markers_to_plot = [m for m in MARKER_TCELL if m in available_markers][:20]
        if markers_to_plot:
            sc.pl.umap(
                adata,
                color=markers_to_plot,
                ncols=4,
                use_raw=True,
                cmap='RdBu_r',
                vmin=-2,
                vmax=2,
                save=f'_markers.{FIGURE_FORMAT}',
                show=False
            )
            log_msg(f"   ✓ Marker UMAPs saved ({len(markers_to_plot)} genes)")
        else:
            log_msg(f"   ⚠️ No markers available for plotting")
    except Exception as e:
        log_msg(f"   ⚠️ Marker UMAP failed: {e}")
    
    # 4) Dotplot for key markers
    if GENERATE_DOTPLOT:
        log_msg("\n   Generating dotplot...")
        try:
            # Select 2-3 representative markers for each cell type
            key_genes = []
            for ct, genes in KEY_MARKERS.items():
                available_in_ct = [g for g in genes if g in available_markers][:2]
                key_genes.extend(available_in_ct)
            
            if key_genes:
                sc.pl.dotplot(
                    adata,
                    var_names=key_genes,
                    groupby='leiden_bbknn',
                    use_raw=True,
                    standard_scale='var',
                    save=f'.{FIGURE_FORMAT}',
                    show=False
                )
                log_msg(f"   ✓ Dotplot saved")
            else:
                log_msg(f"   ⚠️ No key markers available for dotplot")
        except Exception as e:
            log_msg(f"   ⚠️ Dotplot generation failed: {e}")
    
    # 5) Heatmap
    if GENERATE_HEATMAP:
        log_msg("\n   Generating heatmap...")
        try:
            marker_genes = [m for m in MARKER_TCELL if m in available_markers][:30]
            if marker_genes:
                # Calculate mean expression per cluster
                cluster_mean_expr = []
                for cluster in sorted(adata.obs['leiden_bbknn'].unique(), key=str):
                    cluster_cells = adata.obs['leiden_bbknn'] == cluster
                    if adata.raw is not None:
                        cluster_data = adata.raw.X[cluster_cells, :]
                        gene_indices = [adata.raw.var_names.get_loc(g) for g in marker_genes]
                        cluster_expr = cluster_data[:, gene_indices]
                    else:
                        cluster_expr = adata.X[cluster_cells, :][:, [adata.var_names.get_loc(g) for g in marker_genes]]
                    
                    if issparse(cluster_expr):
                        cluster_expr = cluster_expr.toarray()
                    
                    mean_expr = cluster_expr.mean(axis=0)
                    cluster_mean_expr.append(mean_expr)
                
                cluster_mean_expr = pd.DataFrame(
                    cluster_mean_expr,
                    index=sorted(adata.obs['leiden_bbknn'].unique(), key=str),
                    columns=marker_genes
                )
                
                plt.figure(figsize=(12, 8))
                sns.heatmap(
                    cluster_mean_expr.T,
                    cmap='RdYlBu_r',
                    center=0,
                    robust=True,
                    yticklabels=True,
                    xticklabels=True,
                    cbar_kws={'label': 'Mean Expression'}
                )
                plt.title('T Cell Marker Expression Across Clusters')
                plt.xlabel('Cluster')
                plt.ylabel('Marker Genes')
                plt.tight_layout()
                plt.savefig(fig_dir / f'heatmap_marker_expression.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
                plt.close()
                log_msg(f"   ✓ Heatmap saved")
        except Exception as e:
            log_msg(f"   ⚠️ Heatmap generation failed: {e}")
    
    # 6) Batch-split UMAP
    if GENERATE_FACET_PLOTS and BATCH_KEY in adata.obs.columns:
        log_msg("\n   Generating batch-split UMAPs...")
        try:
            batches = sorted(adata.obs[BATCH_KEY].unique())
            n_batches = len(batches)
            n_cols = min(3, n_batches)
            n_rows = (n_batches + n_cols - 1) // n_cols
            
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 5*n_rows))
            axes = [axes] if n_batches == 1 else axes.flatten()
            
            for i, batch in enumerate(batches):
                batch_data = adata[adata.obs[BATCH_KEY] == batch]
                sc.pl.umap(
                    batch_data,
                    color='leiden_bbknn',
                    ax=axes[i],
                    show=False,
                    title=f'{batch}',
                    legend_loc='on data',
                    legend_fontsize=6
                )
            
            for i in range(n_batches, len(axes)):
                axes[i].axis('off')
            
            plt.tight_layout()
            plt.savefig(fig_dir / f'umap_by_batch.{FIGURE_FORMAT}', dpi=DPI, bbox_inches='tight')
            plt.close()
            log_msg(f"   ✓ Batch-split UMAPs saved")
        except Exception as e:
            log_msg(f"   ⚠️ Batch-split UMAP failed: {e}")
    
    log_msg(f"\n   All visualizations saved to: {fig_dir}")

def generate_summary_report(adata, available_markers, output_dir, processing_time=None):
    """Generate analysis summary report"""
    from datetime import datetime
    log_step(8, "Generating Summary Report")
    
    report = []
    report.append("="*70)
    report.append("T Cell Analysis - BBKNN Integration Summary")
    report.append("="*70)
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if processing_time:
        report.append(f"Processing time: {processing_time:.1f} seconds ({processing_time/60:.1f} minutes)")
    
    report.append("")
    report.append("[ Data Overview ]")
    report.append(f"  Cells: {adata.n_obs:,}")
    report.append(f"  Genes: {adata.n_vars:,}")
    report.append("")
    
    if BATCH_KEY in adata.obs.columns:
        report.append(f"[ Batch Distribution (key: {BATCH_KEY}) ]")
        batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
        for batch, count in batch_counts.items():
            pct = count / adata.n_obs * 100
            report.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
        report.append("")
    
    report.append("[ T Cell Marker Genes ]")
    report.append(f"  Total markers: {len(MARKER_TCELL)}")
    report.append(f"  Available: {len(available_markers)} ({len(available_markers)/len(MARKER_TCELL)*100:.1f}%)")
    report.append("")
    
    if RUN_CLUSTERING:
        report.append("[ Multi-resolution Clustering Results ]")
        for res in LEIDEN_RESOLUTIONS:
            key = f'leiden_bbknn_res{res}'
            if key in adata.obs.columns:
                n_clusters = adata.obs[key].nunique()
                default_marker = " (default)" if res == DEFAULT_RESOLUTION else ""
                report.append(f"  Resolution {res}: {n_clusters} clusters{default_marker}")
        report.append("")
    
    report.append("[ BBKNN Configuration ]")
    report.append(f"  neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    report.append(f"  n_pcs: {BBKNN_N_PCS}")
    report.append(f"  metric: {BBKNN_METRIC}")
    report.append(f"  batch_key: {BATCH_KEY}")
    report.append(f"  High variable genes: {USE_HVG} (n={N_TOP_GENES if USE_HVG else 'N/A'})")
    report.append("")
    
    report.append("[ Output Files ]")
    report.append(f"  - Data: {Path(output_dir) / 'tcell_bbknn_integrated.h5ad'}")
    report.append(f"  - Figures: {Path(output_dir) / 'figures/'}*.{FIGURE_FORMAT}")
    if RUN_FIND_MARKERS:
        report.append(f"  - Markers: {Path(output_dir) / 'cluster_markers.csv'}")
    report.append("")
    
    report.append("="*70)
    report.append("Analysis completed successfully")
    report.append("="*70)
    
    report_text = '\n'.join(report)
    report_path = Path(output_dir) / "analysis_summary.txt"
    with open(report_path, 'w') as f:
        f.write(report_text)
    
    log_msg(f"\nReport saved to: {report_path}")
    log_msg("\n" + report_text)

# ==================== Main Program ====================

def main():
    print("\n" + "="*70)
    print("T Cell BBKNN Integration Pipeline (Unified v1.1 - Fixed)")
    print("="*70)

    start_time = time.time()

    log_msg("\nChecking Python environment...")
    try:
        import scanpy as sc   # noqa
        from bbknn import bbknn as bbknn_func  # noqa
        log_msg(f"   scanpy: {sc.__version__}")
        log_msg(f"   bbknn: installed")
    except ImportError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("\nInstallation: pip install scanpy bbknn", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_msg(f"\nOutput directory: {output_dir}")

    # Step 1
    adata = load_anndata(INPUT_H5AD_PATH)

    # Step 2
    available_markers = check_marker_genes(adata)

    # Step 3
    adata_full, adata_hvg = preprocess_for_bbknn(adata)

    # Step 3.5
    adata_hvg = run_pca(adata_hvg)

    # Step 4
    adata_hvg = run_bbknn_integration(adata_hvg)

    # Step 5
    if RUN_UMAP or RUN_CLUSTERING:
        adata_hvg = run_umap_and_clustering(adata_hvg)

    # Step 6
    all_markers = None
    if RUN_FIND_MARKERS:
        all_markers = find_all_markers_optimized(
            adata_hvg,
            'leiden_bbknn',
            output_dir,
            min_pct=MARKER_MIN_PCT,
            logfc_threshold=MARKER_LOGFC_THRESHOLD
        )

    # Step 7
    if RUN_UMAP:
        generate_visualizations(adata_hvg, available_markers, output_dir)

    # Transfer results back to full dataset
    log_msg("\nTransferring results to full dataset...")
    for key in ['X_pca', 'X_umap']:
        if key in adata_hvg.obsm:
            adata_full.obsm[key] = adata_hvg.obsm[key]
    
    for col in adata_hvg.obs.columns:
        if col.startswith('leiden_bbknn'):
            adata_full.obs[col] = adata_hvg.obs[col]
    
    if 'neighbors' in adata_hvg.uns:
        adata_full.uns['neighbors'] = adata_hvg.uns['neighbors']
    
    if 'connectivities' in adata_hvg.obsp:
        adata_full.obsp['connectivities'] = adata_hvg.obsp['connectivities']
    
    if 'distances' in adata_hvg.obsp:
        adata_full.obsp['distances'] = adata_hvg.obsp['distances']

    # Save
    log_msg("\nSaving integrated data...")
    final_path = output_dir / "tcell_bbknn_integrated.h5ad"
    log_msg("   Using gzip compression...")
    adata_full.write_h5ad(final_path, compression='gzip', compression_opts=9)
    file_size = final_path.stat().st_size / (1024**3)
    log_msg(f"   Saved: {final_path} ({file_size:.2f} GB)")

    # Report
    total_time = time.time() - start_time
    generate_summary_report(adata_full, available_markers, output_dir, total_time)

    # Terminal summary
    print("\n" + "="*70)
    print("✅ All analyses completed successfully")
    print("="*70)
    print(f"\nOutput directory: {output_dir}")
    print(f"Data: {final_path} ({file_size:.2f} GB)")
    print(f"Figures: {output_dir / 'figures/'}*.{FIGURE_FORMAT}")
    if RUN_FIND_MARKERS:
        print(f"Markers: {output_dir / 'cluster_markers.csv'}")
    print(f"\nTotal time: {total_time:.1f} seconds ({total_time/60:.1f} minutes)")

    print(f"\n📊 Key Results:")
    print(f"   Cells: {adata_full.n_obs:,}")
    print(f"   Available markers: {len(available_markers)}/{len(MARKER_TCELL)}")
    print(f"   Default clusters: {adata_full.obs['leiden_bbknn'].nunique()}")

    print(f"\n🔧 BBKNN Parameters:")
    print(f"   neighbors_within_batch: {BBKNN_NEIGHBORS_WITHIN_BATCH}")
    print(f"   n_pcs: {BBKNN_N_PCS}")
    print(f"   metric: {BBKNN_METRIC}")
    print(f"   HVG: {N_TOP_GENES}")
    print()

    return adata_full, output_dir

if __name__ == "__main__":
    try:
        adata, output_dir = main()
    except KeyboardInterrupt:
        print("\n\n⚠️ User interrupted", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)