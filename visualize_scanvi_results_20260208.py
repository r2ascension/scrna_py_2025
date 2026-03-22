#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
Visualization Script: scANVI Training Results Analysis
==============================================================================

Version: 1.0
Date: 2026-02-07
Author: r2end

Purpose:
    Comprehensive visualization of scANVI training results including:
    - UMAP plots (data source, cell types, batch, confidence)
    - Prediction distribution analysis
    - Marker gene expression panels
    - Quality control metrics
    - Reference vs Query comparison

Usage:
    python visualize_scanvi_results.py \\
        --input /path/to/merged_scanvi_L2_prod.h5ad \\
        --output /path/to/output_dir \\
        --markers CD19,MS4A1,CD27,IGHD,MZB1,SDC1,JCHAIN,XBP1

"""

import sys
import warnings
import argparse
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

import scanpy as sc

warnings.filterwarnings("ignore")

print("=" * 80)
print("scANVI Results Visualization")
print("=" * 80)

# ==============================================================================
# ARGUMENT PARSER
# ==============================================================================

parser = argparse.ArgumentParser(description="Visualize scANVI training results")
parser.add_argument(
    "--input",
    type=str,
    required=True,
    help="Path to trained h5ad file"
)
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="Output directory for figures"
)
parser.add_argument(
    "--markers",
    type=str,
    default="CD19,MS4A1,CD27,IGHD,IGHM,MZB1,SDC1,JCHAIN,XBP1,PRDM1",
    help="Comma-separated marker genes (default: B cell markers)"
)
parser.add_argument(
    "--dpi",
    type=int,
    default=300,
    help="Figure DPI (default: 300)"
)
parser.add_argument(
    "--format",
    type=str,
    default="pdf",
    choices=["pdf", "png"],
    help="Figure format (default: pdf)"
)

args = parser.parse_args()

INPUT_H5AD = args.input
OUTPUT_DIR = args.output
MARKER_GENES = [g.strip() for g in args.markers.split(",")]
DPI = args.dpi
FIG_FORMAT = args.format

print(f"\nInput: {INPUT_H5AD}")
print(f"Output: {OUTPUT_DIR}")
print(f"Markers: {MARKER_GENES}")

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Column names (adjust if your data uses different names)
DATASOURCE_KEY = "data_source"
REF_LABEL_KEY = "Cell_Type_L2"
TRAIN_LABEL_KEY = "Cell_Type_L2_train"
PRED_LABEL_KEY = "L2_scanvi_pred"
CONF_KEY = "L2_scanvi_confidence"
BATCH_KEY = "sample"
TISSUE_KEY = "tissue"

# UMAP keys
UMAP_KEY = "X_umap_scanvi"
LATENT_KEY = "X_scANVI_L2"

# Color schemes
PALETTE_SOURCE = {"reference": "#1f77b4", "query": "#ff7f0e"}
CMAP_CONFIDENCE = "viridis"
CMAP_EXPRESSION = "Reds"

# Plot settings
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIG_FORMAT, vector_friendly=False)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def save_rasterized_figure(fig, path, dpi=300):
    """Save figure with scatter plots rasterized."""
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✅ Saved: {path.name}")

def get_available_markers(adata, markers):
    """Get markers that exist in the data."""
    return [g for g in markers if g in adata.var_names]

def create_celltype_palette(categories):
    """Create fixed color palette for cell types."""
    n = len(categories)
    if n <= 20:
        colors = list(plt.get_cmap("tab20").colors)[:n]
    elif n <= 102:
        colors = sc.pl.palettes.default_102[:n]
    else:
        colors = [plt.cm.hsv(i / n) for i in range(n)]
    return dict(zip(categories, colors))

# ==============================================================================
# LOAD DATA
# ==============================================================================

print("\n" + "=" * 80)
print("LOADING DATA")
print("=" * 80)

print(f"\nLoading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"Shape: {adata.shape}")

# Verify required keys
required_obs = [DATASOURCE_KEY, PRED_LABEL_KEY]
required_obsm = [UMAP_KEY]

missing = []
for key in required_obs:
    if key not in adata.obs.columns:
        missing.append(f"obs['{key}']")
for key in required_obsm:
    if key not in adata.obsm.keys():
        missing.append(f"obsm['{key}']")

if missing:
    raise ValueError(f"Missing required keys: {missing}")

print("✅ All required keys present")

# Create output directory
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"\nOutput directory: {output_dir}")

# ==============================================================================
# DATA SUMMARY
# ==============================================================================

print("\n" + "=" * 80)
print("DATA SUMMARY")
print("=" * 80)

# Data source
if DATASOURCE_KEY in adata.obs.columns:
    print(f"\nData source distribution:")
    for src, count in adata.obs[DATASOURCE_KEY].value_counts().items():
        pct = count / adata.n_obs * 100
        print(f"  {src}: {count:,} ({pct:.1f}%)")

# Predictions
if PRED_LABEL_KEY in adata.obs.columns:
    n_types = adata.obs[PRED_LABEL_KEY].nunique()
    print(f"\nPredicted cell types: {n_types}")
    print(f"\nTop 10 predictions:")
    for ct, count in adata.obs[PRED_LABEL_KEY].value_counts().head(10).items():
        pct = count / adata.n_obs * 100
        print(f"  {ct}: {count:,} ({pct:.1f}%)")

# Confidence
if CONF_KEY in adata.obs.columns:
    conf = adata.obs[CONF_KEY]
    print(f"\nConfidence statistics:")
    print(f"  Mean: {conf.mean():.3f}")
    print(f"  Median: {conf.median():.3f}")
    print(f"  Q1: {conf.quantile(0.25):.3f}")
    print(f"  Q3: {conf.quantile(0.75):.3f}")
    print(f"  High (>0.9): {(conf > 0.9).sum():,} ({(conf > 0.9).mean()*100:.1f}%)")
    print(f"  Medium (0.5-0.9): {((conf >= 0.5) & (conf <= 0.9)).sum():,} ({((conf >= 0.5) & (conf <= 0.9)).mean()*100:.1f}%)")
    print(f"  Low (<0.5): {(conf < 0.5).sum():,} ({(conf < 0.5).mean()*100:.1f}%)")

# Markers
available_markers = get_available_markers(adata, MARKER_GENES)
print(f"\nMarker genes:")
print(f"  Requested: {len(MARKER_GENES)}")
print(f"  Available: {len(available_markers)}")
if len(available_markers) > 0:
    print(f"  Genes: {available_markers}")

# ==============================================================================
# VISUALIZATION 1: Overview (6-panel)
# ==============================================================================

print("\n" + "=" * 80)
print("VISUALIZATION 1: Overview (6-panel)")
print("=" * 80)

fig, axes = plt.subplots(2, 3, figsize=(24, 16))

# Get cell type categories and palette
ct_categories = sorted(adata.obs[PRED_LABEL_KEY].unique())
ct_palette = create_celltype_palette(ct_categories)

# Panel 1: Data source
if DATASOURCE_KEY in adata.obs.columns:
    sc.pl.umap(
        adata,
        color=DATASOURCE_KEY,
        ax=axes[0, 0],
        show=False,
        title="Data Source",
        palette=PALETTE_SOURCE,
        frameon=False,
        s=20
    )

# Panel 2: Predicted cell types
sc.pl.umap(
    adata,
    color=PRED_LABEL_KEY,
    ax=axes[0, 1],
    show=False,
    title="Predicted Cell Types",
    legend_loc="right margin",
    palette=ct_palette,
    frameon=False,
    s=20
)

# Panel 3: Confidence
if CONF_KEY in adata.obs.columns:
    sc.pl.umap(
        adata,
        color=CONF_KEY,
        ax=axes[0, 2],
        show=False,
        title="Mapping Confidence",
        cmap=CMAP_CONFIDENCE,
        vmin=0,
        vmax=1,
        frameon=False,
        s=20
    )
else:
    axes[0, 2].axis('off')

# Panel 4: Batch
if BATCH_KEY in adata.obs.columns:
    sc.pl.umap(
        adata,
        color=BATCH_KEY,
        ax=axes[1, 0],
        show=False,
        title="Batch",
        frameon=False,
        s=20,
        legend_loc=None
    )
else:
    axes[1, 0].axis('off')

# Panel 5: Tissue
if TISSUE_KEY in adata.obs.columns:
    sc.pl.umap(
        adata,
        color=TISSUE_KEY,
        ax=axes[1, 1],
        show=False,
        title="Tissue",
        frameon=False,
        s=20,
        legend_loc="right margin"
    )
else:
    axes[1, 1].axis('off')

# Panel 6: First marker (if available)
if available_markers:
    marker = available_markers[0]
    sc.pl.umap(
        adata,
        color=marker,
        ax=axes[1, 2],
        show=False,
        title=f"{marker} Expression",
        cmap=CMAP_EXPRESSION,
        frameon=False,
        s=20
    )
else:
    axes[1, 2].text(
        0.5, 0.5,
        "No markers\navailable",
        ha='center', va='center',
        transform=axes[1, 2].transAxes,
        fontsize=14
    )
    axes[1, 2].axis('off')

plt.tight_layout()
output_path = output_dir / f"overview_6panel.{FIG_FORMAT}"
save_rasterized_figure(fig, output_path, dpi=DPI)

# ==============================================================================
# VISUALIZATION 2: Reference vs Query Comparison
# ==============================================================================

if DATASOURCE_KEY in adata.obs.columns and REF_LABEL_KEY in adata.obs.columns:
    print("\n" + "=" * 80)
    print("VISUALIZATION 2: Reference vs Query Comparison")
    print("=" * 80)
    
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    # Left: Reference ground truth
    ref_mask = adata.obs[DATASOURCE_KEY] == "reference"
    adata_ref = adata[ref_mask]
    
    sc.pl.umap(
        adata_ref,
        color=REF_LABEL_KEY,
        ax=axes[0],
        show=False,
        title="Reference: Ground Truth Labels",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=30
    )
    
    # Right: Query predictions
    qry_mask = adata.obs[DATASOURCE_KEY] == "query"
    adata_qry = adata[qry_mask]
    
    sc.pl.umap(
        adata_qry,
        color=PRED_LABEL_KEY,
        ax=axes[1],
        show=False,
        title="Query: Predicted Labels",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=30
    )
    
    plt.tight_layout()
    output_path = output_dir / f"reference_vs_query_comparison.{FIG_FORMAT}"
    save_rasterized_figure(fig, output_path, dpi=DPI)

# ==============================================================================
# VISUALIZATION 3: Marker Gene Expression Panel
# ==============================================================================

if available_markers:
    print("\n" + "=" * 80)
    print("VISUALIZATION 3: Marker Gene Expression Panel")
    print("=" * 80)
    
    n_markers = len(available_markers)
    ncols = 4
    nrows = int(np.ceil(n_markers / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 5*nrows))
    if nrows == 1:
        axes = axes.reshape(1, -1)
    
    for i, marker in enumerate(available_markers):
        row = i // ncols
        col = i % ncols
        
        sc.pl.umap(
            adata,
            color=marker,
            ax=axes[row, col],
            show=False,
            title=f"{marker}",
            cmap=CMAP_EXPRESSION,
            frameon=False,
            s=30
        )
    
    # Hide unused subplots
    for i in range(n_markers, nrows * ncols):
        row = i // ncols
        col = i % ncols
        axes[row, col].axis('off')
    
    plt.tight_layout()
    output_path = output_dir / f"marker_expression_panel.{FIG_FORMAT}"
    save_rasterized_figure(fig, output_path, dpi=DPI)

# ==============================================================================
# VISUALIZATION 4: Cell Type Distribution
# ==============================================================================

print("\n" + "=" * 80)
print("VISUALIZATION 4: Cell Type Distribution")
print("=" * 80)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Bar plot
ct_counts = adata.obs[PRED_LABEL_KEY].value_counts()
ct_sorted = ct_counts.sort_values(ascending=True)

axes[0].barh(range(len(ct_sorted)), ct_sorted.values, color='steelblue')
axes[0].set_yticks(range(len(ct_sorted)))
axes[0].set_yticklabels(ct_sorted.index)
axes[0].set_xlabel('Number of Cells', fontsize=12)
axes[0].set_title('Cell Type Distribution', fontsize=14, fontweight='bold')
axes[0].grid(axis='x', alpha=0.3)
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# Pie chart
axes[1].pie(
    ct_counts.values[:10],
    labels=ct_counts.index[:10],
    autopct='%1.1f%%',
    startangle=90
)
axes[1].set_title('Cell Type Proportions (Top 10)', fontsize=14, fontweight='bold')

plt.tight_layout()
output_path = output_dir / f"celltype_distribution.{FIG_FORMAT}"
plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
plt.close(fig)
print(f"  ✅ Saved: {output_path.name}")

# ==============================================================================
# VISUALIZATION 5: Confidence Analysis
# ==============================================================================

if CONF_KEY in adata.obs.columns:
    print("\n" + "=" * 80)
    print("VISUALIZATION 5: Confidence Analysis")
    print("=" * 80)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    conf = adata.obs[CONF_KEY]
    
    # Panel 1: Histogram
    axes[0, 0].hist(conf, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].axvline(0.9, color='red', linestyle='--', linewidth=2, label='High threshold (0.9)')
    axes[0, 0].axvline(0.5, color='orange', linestyle='--', linewidth=2, label='Medium threshold (0.5)')
    axes[0, 0].set_xlabel('Confidence Score', fontsize=12)
    axes[0, 0].set_ylabel('Number of Cells', fontsize=12)
    axes[0, 0].set_title('Confidence Distribution', fontsize=14, fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].spines['top'].set_visible(False)
    axes[0, 0].spines['right'].set_visible(False)
    
    # Panel 2: Confidence by cell type (violin)
    ct_conf_df = pd.DataFrame({
        'cell_type': adata.obs[PRED_LABEL_KEY],
        'confidence': conf
    })
    
    # Sort by median confidence
    ct_order = ct_conf_df.groupby('cell_type')['confidence'].median().sort_values(ascending=False).index[:15]
    ct_conf_df_top = ct_conf_df[ct_conf_df['cell_type'].isin(ct_order)]
    
    sns.violinplot(
        data=ct_conf_df_top,
        x='cell_type',
        y='confidence',
        order=ct_order,
        ax=axes[0, 1],
        palette='Set2'
    )
    axes[0, 1].axhline(0.9, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 1].axhline(0.5, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    axes[0, 1].set_xlabel('Cell Type', fontsize=12)
    axes[0, 1].set_ylabel('Confidence', fontsize=12)
    axes[0, 1].set_title('Confidence by Cell Type (Top 15)', fontsize=14, fontweight='bold')
    axes[0, 1].tick_params(axis='x', rotation=45)
    axes[0, 1].grid(axis='y', alpha=0.3)
    
    # Panel 3: Confidence categories pie
    conf_cats = pd.cut(
        conf,
        bins=[0, 0.5, 0.9, 1.0],
        labels=['Low (<0.5)', 'Medium (0.5-0.9)', 'High (>0.9)']
    )
    conf_cat_counts = conf_cats.value_counts()
    
    axes[1, 0].pie(
        conf_cat_counts.values,
        labels=conf_cat_counts.index,
        autopct='%1.1f%%',
        colors=['#ff6b6b', '#ffd93d', '#6bcf7f'],
        startangle=90
    )
    axes[1, 0].set_title('Confidence Categories', fontsize=14, fontweight='bold')
    
    # Panel 4: Box plot by cell type
    sns.boxplot(
        data=ct_conf_df_top,
        x='cell_type',
        y='confidence',
        order=ct_order,
        ax=axes[1, 1],
        palette='Set2'
    )
    axes[1, 1].axhline(0.9, color='red', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 1].axhline(0.5, color='orange', linestyle='--', linewidth=1, alpha=0.5)
    axes[1, 1].set_xlabel('Cell Type', fontsize=12)
    axes[1, 1].set_ylabel('Confidence', fontsize=12)
    axes[1, 1].set_title('Confidence Distribution by Cell Type (Top 15)', fontsize=14, fontweight='bold')
    axes[1, 1].tick_params(axis='x', rotation=45)
    axes[1, 1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    output_path = output_dir / f"confidence_analysis.{FIG_FORMAT}"
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✅ Saved: {output_path.name}")

# ==============================================================================
# VISUALIZATION 6: Training vs Prediction Comparison
# ==============================================================================

if TRAIN_LABEL_KEY in adata.obs.columns:
    print("\n" + "=" * 80)
    print("VISUALIZATION 6: Training vs Prediction Comparison")
    print("=" * 80)
    
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    
    # Left: Training labels
    sc.pl.umap(
        adata,
        color=TRAIN_LABEL_KEY,
        ax=axes[0],
        show=False,
        title="Training Labels (Ref + High-Conf Query)",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=25
    )
    
    # Right: Predictions
    sc.pl.umap(
        adata,
        color=PRED_LABEL_KEY,
        ax=axes[1],
        show=False,
        title="scANVI Predictions",
        legend_loc="right margin",
        palette=ct_palette,
        frameon=False,
        s=25
    )
    
    plt.tight_layout()
    output_path = output_dir / f"training_vs_prediction.{FIG_FORMAT}"
    save_rasterized_figure(fig, output_path, dpi=DPI)

# ==============================================================================
# VISUALIZATION 7: High-Resolution Individual Exports
# ==============================================================================

print("\n" + "=" * 80)
print("VISUALIZATION 7: High-Resolution Individual Exports")
print("=" * 80)

export_params = {'frameon': False, 's': 50, 'show': False}

# 1. Cell types (on data legend)
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(
    adata,
    color=PRED_LABEL_KEY,
    ax=ax,
    title="",
    legend_loc="on data",
    legend_fontsize=10,
    legend_fontoutline=2,
    palette=ct_palette,
    **export_params
)
output_path = output_dir / f"umap_celltypes_highres.{FIG_FORMAT}"
save_rasterized_figure(fig, output_path, dpi=DPI)

# 2. Confidence
if CONF_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=CONF_KEY,
        ax=ax,
        title="",
        cmap=CMAP_CONFIDENCE,
        vmin=0,
        vmax=1,
        **export_params
    )
    output_path = output_dir / f"umap_confidence_highres.{FIG_FORMAT}"
    save_rasterized_figure(fig, output_path, dpi=DPI)

# 3. Data source
if DATASOURCE_KEY in adata.obs.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(
        adata,
        color=DATASOURCE_KEY,
        ax=ax,
        title="",
        palette=PALETTE_SOURCE,
        **export_params
    )
    output_path = output_dir / f"umap_datasource_highres.{FIG_FORMAT}"
    save_rasterized_figure(fig, output_path, dpi=DPI)

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)

summary_stats = {
    'timestamp': datetime.now().isoformat(),
    'input_file': str(INPUT_H5AD),
    'n_cells': int(adata.n_obs),
    'n_genes': int(adata.n_vars),
    'cell_types': adata.obs[PRED_LABEL_KEY].value_counts().to_dict()
}

if DATASOURCE_KEY in adata.obs.columns:
    summary_stats['data_sources'] = adata.obs[DATASOURCE_KEY].value_counts().to_dict()

if CONF_KEY in adata.obs.columns:
    conf = adata.obs[CONF_KEY]
    summary_stats['confidence'] = {
        'mean': float(conf.mean()),
        'median': float(conf.median()),
        'std': float(conf.std()),
        'min': float(conf.min()),
        'max': float(conf.max()),
        'high_gt_0.9': int((conf > 0.9).sum()),
        'medium_0.5_0.9': int(((conf >= 0.5) & (conf <= 0.9)).sum()),
        'low_lt_0.5': int((conf < 0.5).sum())
    }

# Save summary
import json
summary_path = output_dir / "visualization_summary.json"
with open(summary_path, 'w') as f:
    json.dump(summary_stats, f, indent=2)
print(f"\n✅ Summary saved: {summary_path}")

# Save cell type counts to CSV
ct_counts_df = pd.DataFrame({
    'cell_type': adata.obs[PRED_LABEL_KEY].value_counts().index,
    'count': adata.obs[PRED_LABEL_KEY].value_counts().values,
    'percentage': (adata.obs[PRED_LABEL_KEY].value_counts().values / adata.n_obs * 100)
})
ct_counts_path = output_dir / "celltype_counts.csv"
ct_counts_df.to_csv(ct_counts_path, index=False)
print(f"✅ Cell type counts saved: {ct_counts_path}")

# ==============================================================================
# COMPLETION
# ==============================================================================

print("\n" + "=" * 80)
print("VISUALIZATION COMPLETE!")
print("=" * 80)

print(f"\nOutput directory: {OUTPUT_DIR}")
print(f"\nGenerated files:")

output_files = sorted(output_dir.glob(f"*.{FIG_FORMAT}"))
for f in output_files:
    file_size = f.stat().st_size / 1024  # KB
    print(f"  - {f.name} ({file_size:.1f} KB)")

print(f"\nAdditional files:")
print(f"  - visualization_summary.json")
print(f"  - celltype_counts.csv")

print("\n" + "=" * 80)
print("All visualizations saved successfully!")
print("=" * 80)
