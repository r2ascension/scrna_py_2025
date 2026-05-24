#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
T Cell + Myeloid Subclustering Pipeline (scANVI-based) v1.0
==============================================================================

Purpose:
--------
Subcluster analysis of T cells and Myeloid cells using scANVI latent space
as the foundation for clustering and visualization.

Key Features:
-------------
1. Dual-mode: Process T cells and Myeloid cells together or separately
2. scANVI-based: Uses pre-computed X_scANVI latent space
3. BBKNN refinement: Optional second-round batch correction
4. Leiden clustering: Multiple resolutions for hierarchy
5. Marker analysis: Differential expression per cluster
6. Cell type validation: Compare clusters to original annotations

Input Requirements:
-------------------
- H5AD file with scANVI results:
  * obsm['X_scANVI']: scANVI latent representation
  * obsm['X_umap_scANVI']: scANVI UMAP (optional, will compute if missing)
  * obs['scanvi_pred']: scANVI cell type predictions
  * layers['counts']: Raw counts
  * .raw: Full gene set for markers

Output:
-------
- Subclustered h5ad with multiple Leiden resolutions
- Marker genes per cluster (CSV)
- UMAP visualizations comparing scANVI vs BBKNN
- Cluster annotation suggestions

Author: Claude Code
Date: 2026-02-25
Version: 1.0
==============================================================================
"""

import os
import sys
import warnings
import gc
import json
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import scanpy.external as sce
import bbknn

warnings.filterwarnings('ignore')

print("=" * 80)
print("T Cell + Myeloid Subclustering Pipeline (scANVI-based) v1.0")
print("=" * 80)
print()

# Thread control for HPC
for var in ["OMP", "OPENBLAS", "MKL", "NUMEXPR"]:
    os.environ[f"{var}_NUM_THREADS"] = "1"

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# --- Input/Output ---
INPUT_H5AD = "/path/to/your/scanvi_results.h5ad"  # EDIT THIS
OUTPUT_DIR = Path("/home/h2048/data/py/20260225_tcell_myeloid_subcluster")
OUTPUT_PREFIX = "tcell_myeloid_subcluster"

# --- Cell Type Selection ---
# Which scANVI predictions to include for subclustering
TCELL_PREDICTIONS = ["CD4", "CD8", "T_cell"]  # EDIT based on your data
MYELOID_PREDICTIONS = ["Monocyte", "Macrophage", "DC", "Neutrophil"]  # EDIT

# Or process ALL cells together (set to None to use PREDICTIONS above)
PROCESS_ALL = False

# --- Column Names ---
SCANVI_PRED_COL = "scanvi_pred"       # scANVI prediction column
SCANVI_CONF_COL = "scanvi_confidence" # scANVI confidence column
BATCH_KEY = "sample"                  # Batch key for BBKNN
TISSUE_KEY = "tissue"                 # Tissue/Location key
DATA_SOURCE = "data_source"           # "reference" or "query"

# --- Latent Space Keys ---
SCANVI_LATENT_KEY = "X_scANVI"        # scANVI latent
SCANVI_UMAP_KEY = "X_umap_scANVI"     # scANVI UMAP

# --- Subclustering Parameters ---
USE_BBKNN = True                      # Run BBKNN on top of scANVI?
N_HVG = 3000                          # HVG for BBKNN (if USE_BBKNN)
N_PCS = 50                            # PCs for BBKNN
NEIGHBORS_WITHIN_BATCH = 5            # BBKNN parameter

# Leiden resolutions for multi-level clustering
RESOLUTIONS = [0.2, 0.4, 0.8, 1.2]    # Coarse to fine

# --- Quality Control ---
MIN_CELLS_PER_CLUSTER = 10
MIN_CELLS_FOR_MARKER = 20
CONFIDENCE_THRESHOLD = 0.5            # Filter low-confidence predictions

# --- Known Markers ---
TCELL_MARKERS = {
    "pan_T": ["CD3D", "CD3E", "CD3G", "PTPRC"],
    "CD4": ["CD4", "IL7R", "CD40LG"],
    "CD8": ["CD8A", "CD8B"],
    "Naive": ["CCR7", "SELL", "TCF7", "LEF1"],
    "CM": ["CCR7", "CD27"],
    "EM": ["GZMK", "CXCR3", "CCR5"],
    "Effector": ["GZMB", "PRF1", "IFNG"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "Tissue_resident": ["CD69", "ITGAE"],
    "Proliferating": ["MKI67", "TOP2A"],
}

MYELOID_MARKERS = {
    "pan_Myeloid": ["LYZ", "CD14", "FCGR3A"],
    "Monocyte": ["CD14", "FCGR3A", "S100A8", "S100A9"],
    "Macrophage": ["CD68", "CD163", "MRC1", "MARCO"],
    "cDC": ["CLEC9A", "XCR1", "CLEC10A", "FCER1A"],
    "pDC": ["LILRA4", "CLEC4C", "IRF7"],
    "Neutrophil": ["CSF3R", "FCGR3B", "CXCR2"],
    "M1": ["CD86", "CD80", "IL1B"],
    "M2": ["CD163", "MRC1", "ARG1"],
}

# --- Scanpy Settings ---
sc.settings.verbosity = 2
sc.settings.n_jobs = 8
sc.settings.set_figure_params(dpi=100, facecolor='white')

print("Configuration:")
print(f"  Input: {INPUT_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  BBKNN: {USE_BBKNN}")
print(f"  Resolutions: {RESOLUTIONS}")
print()

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "figures").mkdir(exist_ok=True)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def load_and_validate_data(path):
    """Load and validate input data."""
    print("\n[Step 1] Loading and validating data...")

    adata = sc.read_h5ad(path)
    print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # Check required fields
    required = {
        'obsm': {SCANVI_LATENT_KEY: 'scANVI latent'},
        'layers': {'counts': 'raw counts'},
    }

    missing = []
    for category, items in required.items():
        container = getattr(adata, category)
        for key, desc in items.items():
            if key not in container:
                missing.append(f"  {category}['{key}'] ({desc})")

    if missing:
        print("  ❌ Missing required fields:")
        for m in missing:
            print(m)
        raise ValueError("Input validation failed")

    # Check for .raw
    if adata.raw is None:
        print("  ⚠️  WARNING: No .raw attribute - marker analysis limited")
    else:
        print(f"  ✓ .raw available: {adata.raw.n_vars:,} genes")

    # Check scANVI UMAP
    if SCANVI_UMAP_KEY not in adata.obsm:
        print(f"  ℹ️  {SCANVI_UMAP_KEY} not found, will compute from latent")

    print("  ✓ Validation passed")
    return adata


def subset_cell_types(adata):
    """Subset to T cell and Myeloid predictions."""
    if PROCESS_ALL or SCANVI_PRED_COL not in adata.obs.columns:
        print("\n[Step 2] Processing ALL cells")
        return adata

    print("\n[Step 2] Subsetting to T cell and Myeloid types...")

    target_types = TCELL_PREDICTIONS + MYELOID_PREDICTIONS
    mask = adata.obs[SCANVI_PRED_COL].isin(target_types)

    # Also filter by confidence
    if SCANVI_CONF_COL in adata.obs.columns:
        conf_mask = adata.obs[SCANVI_CONF_COL] >= CONFIDENCE_THRESHOLD
        mask = mask & conf_mask
        print(f"  High confidence (≥{CONFIDENCE_THRESHOLD}): {conf_mask.sum():,} cells")

    adata_subset = adata[mask].copy()
    print(f"  Subset: {adata_subset.n_obs:,} cells")

    # Print composition
    print(f"\n  Cell type composition:")
    print(adata_subset.obs[SCANVI_PRED_COL].value_counts())

    return adata_subset


def compute_umap_if_missing(adata):
    """Compute UMAP from scANVI latent if not present."""
    if SCANVI_UMAP_KEY in adata.obsm:
        print(f"\n[Step 3] Using existing {SCANVI_UMAP_KEY}")
        adata.obsm['X_umap'] = adata.obsm[SCANVI_UMAP_KEY].copy()
        return adata

    print(f"\n[Step 3] Computing UMAP from {SCANVI_LATENT_KEY}...")
    sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, n_neighbors=30)
    sc.tl.umap(adata)
    adata.obsm[SCANVI_UMAP_KEY] = adata.obsm['X_umap'].copy()
    print("  ✓ UMAP computed")
    return adata


def run_bbknn_refinement(adata):
    """Optional BBKNN refinement on gene expression."""
    if not USE_BBKNN:
        print("\n[Step 4] Skipping BBKNN (USE_BBKNN=False)")
        return adata, None

    print("\n[Step 4] Running BBKNN refinement...")

    # Prepare log1p data
    if 'log1p' not in adata.layers:
        print("  Computing log1p layer...")
        adata.layers['counts'] = csr_matrix(adata.layers['counts']) if issparse(adata.layers['counts']) else adata.layers['counts']
        adata.layers['log1p'] = sc.pp.normalize_total(adata, layer='counts', inplace=False)['X']
        adata.layers['log1p'] = np.log1p(adata.layers['log1p'])

    adata.X = adata.layers['log1p']

    # HVG selection
    print(f"  Selecting {N_HVG} HVGs...")
    try:
        sc.pp.highly_variable_genes(
            adata, layer='log1p', n_top_genes=N_HVG,
            batch_key=BATCH_KEY, flavor='seurat_v3'
        )
    except:
        sc.pp.highly_variable_genes(
            adata, layer='log1p', n_top_genes=N_HVG, flavor='seurat_v3'
        )

    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    print(f"  HVG subset: {adata_hvg.n_vars} genes")

    # PCA
    print("  Running PCA...")
    sc.tl.pca(adata_hvg, n_comps=N_PCS, svd_solver='arpack')

    # BBKNN
    print("  Running BBKNN...")
    batch_sizes = adata_hvg.obs[BATCH_KEY].value_counts()
    min_batch = batch_sizes.min()
    neighbors_wb = min(NEIGHBORS_WITHIN_BATCH, max(3, min_batch - 1))

    sce.pp.bbknn(
        adata_hvg,
        batch_key=BATCH_KEY,
        neighbors_within_batch=neighbors_wb,
        n_pcs=N_PCS
    )

    # UMAP
    print("  Computing BBKNN UMAP...")
    sc.tl.umap(adata_hvg)

    # Transfer BBKNN results back
    adata.obsm['X_umap_bbknn'] = adata_hvg.obsm['X_umap'].copy()
    adata.obsp['connectivities'] = adata_hvg.obsp['connectivities'].copy()
    adata.obsp['distances'] = adata_hvg.obsp['distances'].copy()
    adata.uns['neighbors'] = adata_hvg.uns['neighbors'].copy()

    print("  ✓ BBKNN refinement complete")
    return adata, adata_hvg


def run_multi_resolution_clustering(adata):
    """Run Leiden clustering at multiple resolutions."""
    print(f"\n[Step 5] Multi-resolution Leiden clustering...")

    # Use BBKNN neighbors if available, else compute from scANVI
    if 'neighbors' not in adata.uns:
        print("  Computing neighbors from scANVI latent...")
        sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, n_neighbors=30)

    for res in RESOLUTIONS:
        key = f'leiden_res{res}'
        print(f"  Clustering at resolution {res}...")
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        print(f"    {key}: {n_clusters} clusters")

    # Default to middle resolution
    adata.obs['leiden'] = adata.obs[f'leiden_res{RESOLUTIONS[1]}'].copy()

    return adata


def find_marker_genes(adata, cluster_key='leiden'):
    """Find differentially expressed genes per cluster."""
    print(f"\n[Step 6] Finding marker genes ({cluster_key})...")

    # Use log1p data
    if 'log1p' not in adata.layers:
        adata.X = sc.pp.normalize_total(adata, layer='counts', inplace=False)['X']
        adata.X = np.log1p(adata.X)
    else:
        adata.X = adata.layers['log1p']

    # Run rank_genes_groups
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method='wilcoxon',
        use_raw=True if adata.raw else False,
        n_genes=100
    )

    # Extract results
    markers_list = []
    for cluster in adata.obs[cluster_key].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df['cluster'] = cluster
        df['cluster_size'] = (adata.obs[cluster_key] == cluster).sum()
        markers_list.append(df)

    markers_df = pd.concat(markers_list, ignore_index=True)

    # FDR correction
    from scipy import stats
    markers_df['pvals_adj_fdr'] = stats.false_discovery_rate(
        markers_df['pvals'].values
    )[1]

    # Filter significant markers
    markers_sig = markers_df[
        (markers_df['pvals_adj'] < 0.05) &
        (markers_df['logfoldchanges'] > 0.5)
    ].copy()

    output_file = OUTPUT_DIR / f"{OUTPUT_PREFIX}_markers_{cluster_key}.csv"
    markers_sig.to_csv(output_file, index=False)
    print(f"  ✓ Saved: {output_file}")
    print(f"    Significant markers: {len(markers_sig)}")

    return markers_sig


def annotate_clusters(adata, markers_df, cluster_key='leiden'):
    """Annotate clusters based on marker overlap."""
    print(f"\n[Step 7] Annotating {cluster_key} clusters...")

    annotations = {}
    all_markers = {**TCELL_MARKERS, **MYELOID_MARKERS}

    for cluster in adata.obs[cluster_key].cat.categories:
        cluster_markers = markers_df[
            (markers_df['cluster'] == cluster) &
            (markers_df['pvals_adj'] < 0.01)
        ]['names'].head(50).tolist()

        # Score each cell type
        scores = {}
        for cell_type, type_markers in all_markers.items():
            overlap = set(cluster_markers) & set(type_markers)
            scores[cell_type] = len(overlap) / len(type_markers) if type_markers else 0

        # Best match
        best_match = max(scores, key=scores.get) if scores else "Unknown"
        best_score = scores.get(best_match, 0)

        annotations[cluster] = {
            'annotation': best_match if best_score > 0.3 else "Mixed",
            'confidence': best_score,
            'top_markers': ','.join(cluster_markers[:5])
        }

    # Add to adata
    anno_col = f'{cluster_key}_annotation'
    adata.obs[anno_col] = adata.obs[cluster_key].map(
        {k: v['annotation'] for k, v in annotations.items()}
    )

    # Summary
    print(f"  Cluster annotations:")
    for cluster, info in annotations.items():
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        print(f"    {cluster}: {info['annotation']} ({n_cells} cells, conf={info['confidence']:.2f})")

    return adata, annotations


def create_visualizations(adata):
    """Create comprehensive visualizations."""
    print("\n[Step 8] Creating visualizations...")

    fig_dir = OUTPUT_DIR / "figures"

    # Figure 1: scANVI overview
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    sc.pl.umap(adata, color=SCANVI_PRED_COL, ax=axes[0,0], show=False,
               title="scANVI Predictions", legend_loc='on data', s=10)
    sc.pl.umap(adata, color=DATA_SOURCE, ax=axes[0,1], show=False,
               title="Data Source", s=10)
    sc.pl.umap(adata, color=SCANVI_CONF_COL, ax=axes[0,2], show=False,
               title="scANVI Confidence", cmap='viridis', vmin=0, vmax=1, s=10)

    # Leiden clusters at different resolutions
    sc.pl.umap(adata, color='leiden', ax=axes[1,0], show=False,
               title=f"Leiden (res={RESOLUTIONS[1]})", legend_loc='on data', s=10)
    sc.pl.umap(adata, color=f'leiden_res{RESOLUTIONS[0]}', ax=axes[1,1], show=False,
               title=f"Leiden (res={RESOLUTIONS[0]}, coarse)", legend_loc='on data', s=10)
    sc.pl.umap(adata, color=f'leiden_res{RESOLUTIONS[-1]}', ax=axes[1,2], show=False,
               title=f"Leiden (res={RESOLUTIONS[-1]}, fine)", legend_loc='on data', s=10)

    plt.tight_layout()
    plt.savefig(fig_dir / f"{OUTPUT_PREFIX}_overview_scanvi.pdf", dpi=300)
    plt.close()

    # Figure 2: BBKNN comparison (if available)
    if 'X_umap_bbknn' in adata.obsm:
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))

        # scANVI UMAP
        sc.pl.umap(adata, color='leiden', ax=axes[0,0], show=False,
                   title="scANVI UMAP + Leiden", s=10)
        sc.pl.umap(adata, color=SCANVI_PRED_COL, ax=axes[0,1], show=False,
                   title="scANVI UMAP + Predictions", s=10)

        # BBKNN UMAP
        adata_bk = adata.copy()
        adata_bk.obsm['X_umap'] = adata_bk.obsm['X_umap_bbknn']
        sc.pl.umap(adata_bk, color='leiden', ax=axes[1,0], show=False,
                   title="BBKNN UMAP + Leiden", s=10)
        sc.pl.umap(adata_bk, color=SCANVI_PRED_COL, ax=axes[1,1], show=False,
                   title="BBKNN UMAP + Predictions", s=10)

        # Batch distribution
        sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0,2], show=False,
                   title="scANVI: Batches", s=10)
        sc.pl.umap(adata_bk, color=BATCH_KEY, ax=axes[1,2], show=False,
                   title="BBKNN: Batches", s=10)

        del adata_bk

        plt.tight_layout()
        plt.savefig(fig_dir / f"{OUTPUT_PREFIX}_bbknn_comparison.pdf", dpi=300)
        plt.close()

    # Figure 3: Marker gene dotplot
    marker_genes = []
    for m in list(TCELL_MARKERS.values()) + list(MYELOID_MARKERS.values()):
        marker_genes.extend([g for g in m if g in adata.var_names or (adata.raw and g in adata.raw.var_names)])
    marker_genes = list(dict.fromkeys(marker_genes))[:50]  # Top 50 unique

    if marker_genes:
        fig, ax = plt.subplots(figsize=(14, 10))
        sc.pl.dotplot(adata, marker_genes, groupby='leiden',
                      use_raw=True, show=False, ax=ax)
        plt.tight_layout()
        plt.savefig(fig_dir / f"{OUTPUT_PREFIX}_marker_dotplot.pdf", dpi=300)
        plt.close()

    # Figure 4: Key markers on UMAP
    key_markers = ["CD3E", "CD4", "CD8A", "LYZ", "CD14", "CD68"]
    available = [m for m in key_markers if m in adata.var_names or (adata.raw and m in adata.raw.var_names)]

    if available:
        n_rows = (len(available) + 2) // 3
        fig, axes = plt.subplots(n_rows, 3, figsize=(15, 5*n_rows))
        axes = axes.flatten() if n_rows > 1 else [axes] if len(available) == 1 else axes

        for idx, marker in enumerate(available):
            sc.pl.umap(adata, color=marker, ax=axes[idx], show=False,
                       title=marker, cmap='Reds', use_raw=True, s=10)

        for idx in range(len(available), len(axes)):
            axes[idx].axis('off')

        plt.tight_layout()
        plt.savefig(fig_dir / f"{OUTPUT_PREFIX}_key_markers.pdf", dpi=300)
        plt.close()

    print("  ✓ Figures saved")


def create_summary_report(adata, annotations_dict):
    """Create summary statistics."""
    print("\n[Step 9] Creating summary report...")

    summary = {
        'timestamp': datetime.now().isoformat(),
        'input_file': INPUT_H5AD,
        'n_cells': adata.n_obs,
        'n_genes': adata.n_vars,
        'n_batches': adata.obs[BATCH_KEY].nunique(),
        'scanvi_predictions': adata.obs[SCANVI_PRED_COL].value_counts().to_dict(),
    }

    # Cluster stats per resolution
    for res in RESOLUTIONS:
        key = f'leiden_res{res}'
        summary[f'clusters_res{res}'] = adata.obs[key].nunique()

    # Save
    with open(OUTPUT_DIR / f"{OUTPUT_PREFIX}_summary.json", 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    # CSV summary
    cluster_summary = []
    for res in RESOLUTIONS:
        key = f'leiden_res{res}'
        for cluster in adata.obs[key].cat.categories:
            mask = adata.obs[key] == cluster
            cluster_summary.append({
                'resolution': res,
                'cluster': cluster,
                'n_cells': mask.sum(),
                'scanvi_pred_majority': adata.obs.loc[mask, SCANVI_PRED_COL].mode()[0],
                'data_source_ref': (adata.obs.loc[mask, DATA_SOURCE] == 'reference').sum(),
                'data_source_qry': (adata.obs.loc[mask, DATA_SOURCE] == 'query').sum(),
            })

    summary_df = pd.DataFrame(cluster_summary)
    summary_df.to_csv(OUTPUT_DIR / f"{OUTPUT_PREFIX}_cluster_summary.csv", index=False)

    print("  ✓ Summary saved")
    return summary


# ==============================================================================
# MAIN PIPELINE
# ==============================================================================

def main():
    # Step 1: Load data
    adata = load_and_validate_data(INPUT_H5AD)

    # Step 2: Subset cell types
    adata = subset_cell_types(adata)

    # Step 3: Compute UMAP if missing
    adata = compute_umap_if_missing(adata)

    # Step 4: BBKNN refinement (optional)
    adata, adata_hvg = run_bbknn_refinement(adata)

    # Step 5: Multi-resolution clustering
    adata = run_multi_resolution_clustering(adata)

    # Step 6-7: Find markers and annotate (for default resolution)
    markers = find_marker_genes(adata, 'leiden')
    adata, annotations = annotate_clusters(adata, markers, 'leiden')

    # Also for other resolutions
    for res in RESOLUTIONS:
        if res != RESOLUTIONS[1]:
            _ = find_marker_genes(adata, f'leiden_res{res}')

    # Step 8: Visualizations
    create_visualizations(adata)

    # Step 9: Summary
    summary = create_summary_report(adata, annotations)

    # Step 10: Save final data
    print("\n[Step 10] Saving results...")
    output_file = OUTPUT_DIR / f"{OUTPUT_PREFIX}_results.h5ad"
    adata.write_h5ad(output_file, compression='gzip')
    print(f"  ✓ Saved: {output_file}")

    # Clean up
    if adata_hvg is not None:
        del adata_hvg
    gc.collect()

    # Final summary
    print("\n" + "=" * 80)
    print("PIPELINE COMPLETE")
    print("=" * 80)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"\nKey outputs:")
    print(f"  - {OUTPUT_PREFIX}_results.h5ad")
    print(f"  - {OUTPUT_PREFIX}_markers_leiden.csv")
    print(f"  - {OUTPUT_PREFIX}_cluster_summary.csv")
    print(f"  - figures/*.pdf")
    print("\nClustering resolutions available:")
    for res in RESOLUTIONS:
        n_clust = adata.obs[f'leiden_res{res}'].nunique()
        print(f"  - leiden_res{res}: {n_clust} clusters")
    print("=" * 80)


if __name__ == "__main__":
    main()
