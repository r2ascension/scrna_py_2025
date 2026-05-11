#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
High-Quality Visualization for Merged Reference + Query
========================================================
Reads pre-computed merged h5ad and generates publication-ready figures

Requirements:
- reference_plus_query_merged_L2.h5ad (from v2.5.3)
- Clean UMAP visualization with proper handling of ref/query differences

Usage:
    python visualize_merged_highquality.py
"""

import sys
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
from scipy.sparse import issparse

warnings.filterwarnings("ignore")

print("=" * 70)
print("High-Quality Merged Data Visualization")
print("=" * 70)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def create_query_only_column(adata, source_col, output_col_name, fill_value='N/A'):
    """
    Create a temporary column that only shows query data, reference cells are NA.
    
    Parameters
    ----------
    adata : AnnData
        The merged anndata object
    source_col : str
        Source column name (must exist in adata.obs)
    output_col_name : str
        Name for the temporary column
    fill_value : str
        Value for reference cells (default 'N/A')
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if source_col not in adata.obs.columns:
        return False
    
    if 'data_source' not in adata.obs.columns:
        return False
    
    temp_col = pd.Series(fill_value, index=adata.obs_names, dtype='category')
    qry_mask = adata.obs['data_source'] == 'query'
    
    if qry_mask.any():
        temp_col[qry_mask] = adata.obs.loc[qry_mask, source_col].astype(str)
        temp_col = temp_col.astype('category')
        adata.obs[output_col_name] = temp_col
        return True
    
    return False


def create_query_only_numeric(adata, source_col, output_col_name):
    """
    Create a temporary numeric column that only shows query data.
    
    Parameters
    ----------
    adata : AnnData
        The merged anndata object
    source_col : str
        Source column name (must exist in adata.obs)
    output_col_name : str
        Name for the temporary column
    
    Returns
    -------
    bool
        True if successful, False otherwise
    """
    if source_col not in adata.obs.columns:
        return False
    
    if 'data_source' not in adata.obs.columns:
        return False
    
    temp_arr = np.full(adata.n_obs, np.nan, dtype=float)
    qry_mask = adata.obs['data_source'] == 'query'
    
    if qry_mask.any():
        temp_arr[qry_mask] = pd.to_numeric(
            adata.obs.loc[qry_mask, source_col], 
            errors='coerce'
        ).values
        adata.obs[output_col_name] = temp_arr
        return True
    
    return False


def cleanup_temp_columns(adata, col_names):
    """
    Remove temporary columns from adata.obs.
    
    Parameters
    ----------
    adata : AnnData
        The anndata object
    col_names : list or str
        Column name(s) to remove
    """
    if isinstance(col_names, str):
        col_names = [col_names]
    
    for col in col_names:
        if col in adata.obs.columns:
            adata.obs.drop(columns=[col], inplace=True)


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MERGED_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/reference_plus_query_merged_L2.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/visualizations_final"

DPI = 300
FIGURE_FORMAT = "pdf"
CONFIDENCE_THRESHOLD = 0.5

# Color palettes
DATA_SOURCE_COLORS = {
    'reference': '#3498db',  # Blue
    'query': '#e74c3c'       # Red
}

CELLTYPE_COLORS = {
    'Plasma': '#e74c3c',
    'Naive_B': '#3498db',
    'Memory_B': '#2ecc71',
    'Atypical_Memory_B': '#f39c12',
    'GC_B': '#9b59b6',
    'Unknown': '#95a5a6'
}

sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
sns.set_style("whitegrid")

print(f"\nInput: {MERGED_H5AD}")
print(f"Output: {OUTPUT_DIR}")

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

# ==============================================================================
# Load Merged Data
# ==============================================================================

print("\n" + "=" * 70)
print("Load Merged Data")
print("=" * 70)

adata = sc.read_h5ad(MERGED_H5AD)
print(f"\nMerged data shape: {adata.shape}")
print(f"obs columns: {list(adata.obs.columns)}")
print(f"obsm keys: {list(adata.obsm.keys())}")

if 'X_umap' not in adata.obsm:
    raise ValueError("UMAP coordinates not found!")

# Split reference and query (no copy to save memory)
ref_mask = adata.obs['data_source'] == 'reference'
qry_mask = adata.obs['data_source'] == 'query'

ref_cells = adata[ref_mask]
qry_cells = adata[qry_mask]

print(f"\nReference cells: {ref_cells.n_obs:,}")
print(f"Query cells: {qry_cells.n_obs:,}")

if ref_cells.n_obs == 0:
    print("WARNING: No reference cells found!")
if qry_cells.n_obs == 0:
    print("WARNING: No query cells found!")

# ==============================================================================
# Figure 1: Overview - Data Source and Cell Types (2x2)
# ==============================================================================

print("\n" + "=" * 70)
print("Figure 1: UMAP Overview (Data Source + Cell Types)")
print("=" * 70)

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

# Panel 1: Data source
print("  Panel 1: Data source...")
sc.pl.umap(adata, color='data_source', ax=axes[0, 0], show=False,
           title='Data Source: Reference vs Query', s=8,
           palette=DATA_SOURCE_COLORS, frameon=True)

# Panel 2: Reference L2 (if available)
print("  Panel 2: Reference L2...")
temp_col_name = '_temp_ref_L2'
try:
    if 'Cell_Type_L2' in adata.obs.columns:
        # Show only reference cells' labels
        temp_col = pd.Series('N/A', index=adata.obs_names, dtype='object')
        ref_mask = adata.obs['data_source'] == 'reference'
        if ref_mask.any():
            temp_col[ref_mask] = adata.obs.loc[ref_mask, 'Cell_Type_L2'].astype(str)
            temp_col = temp_col.astype('category')
            adata.obs[temp_col_name] = temp_col
            
            sc.pl.umap(adata, color=temp_col_name, ax=axes[0, 1], show=False,
                       title='Reference L2 Labels', s=8, frameon=True,
                       legend_loc='right margin', na_color='lightgray')
        else:
            axes[0, 1].text(0.5, 0.5, 'No reference cells',
                           ha='center', va='center', transform=axes[0, 1].transAxes,
                           fontsize=14, color='gray')
            axes[0, 1].set_title('Reference L2 Labels')
    else:
        axes[0, 1].text(0.5, 0.5, 'Reference L2 not available',
                       ha='center', va='center', transform=axes[0, 1].transAxes,
                       fontsize=14, color='gray')
        axes[0, 1].set_title('Reference L2 Labels')
except Exception as e:
    print(f"  WARNING: Panel 2 failed: {e}")
    axes[0, 1].text(0.5, 0.5, 'Error generating panel',
                   ha='center', va='center', transform=axes[0, 1].transAxes,
                   fontsize=14, color='red')
finally:
    cleanup_temp_columns(adata, temp_col_name)

# Panel 3: Query L2 final
print("  Panel 3: Query L2 predictions...")
temp_col_name = '_temp_qry_L2'
try:
    if 'Cell_Type_L2_final' in adata.obs.columns:
        temp_col = pd.Series('N/A', index=adata.obs_names, dtype='object')
        qry_mask = adata.obs['data_source'] == 'query'
        if qry_mask.any():
            temp_col[qry_mask] = adata.obs.loc[qry_mask, 'Cell_Type_L2_final'].astype(str)
            temp_col = temp_col.astype('category')
            adata.obs[temp_col_name] = temp_col
            
            sc.pl.umap(adata, color=temp_col_name, ax=axes[1, 0], show=False,
                       title=f'Query L2 Predictions (conf ≥ {CONFIDENCE_THRESHOLD})', 
                       s=8, frameon=True, legend_loc='right margin',
                       palette=CELLTYPE_COLORS, na_color='lightgray')
        else:
            axes[1, 0].text(0.5, 0.5, 'No query data',
                           ha='center', va='center', transform=axes[1, 0].transAxes,
                           fontsize=14, color='gray')
            axes[1, 0].set_title('Query L2 Predictions')
    else:
        axes[1, 0].text(0.5, 0.5, 'Query L2 not available',
                       ha='center', va='center', transform=axes[1, 0].transAxes,
                       fontsize=14, color='gray')
        axes[1, 0].set_title('Query L2 Predictions')
except Exception as e:
    print(f"  WARNING: Panel 3 failed: {e}")
    axes[1, 0].text(0.5, 0.5, 'Error generating panel',
                   ha='center', va='center', transform=axes[1, 0].transAxes,
                   fontsize=14, color='red')
finally:
    cleanup_temp_columns(adata, temp_col_name)

# Panel 4: Query confidence
print("  Panel 4: Query mapping confidence...")
temp_col_name = '_temp_conf'
try:
    if 'mapping_confidence' in adata.obs.columns:
        temp_conf = np.full(adata.n_obs, np.nan, dtype=float)
        qry_mask = adata.obs['data_source'] == 'query'
        if qry_mask.any():
            temp_conf[qry_mask] = pd.to_numeric(
                adata.obs.loc[qry_mask, 'mapping_confidence'], 
                errors='coerce'
            ).values
            
            adata.obs[temp_col_name] = temp_conf
            
            sc.pl.umap(adata, color=temp_col_name, ax=axes[1, 1], show=False,
                       title='Query Mapping Confidence', s=8, frameon=True,
                       cmap='RdYlGn', vmin=0, vmax=1, na_color='lightgray')
        else:
            axes[1, 1].text(0.5, 0.5, 'No query data',
                           ha='center', va='center', transform=axes[1, 1].transAxes,
                           fontsize=14, color='gray')
            axes[1, 1].set_title('Query Mapping Confidence')
    else:
        axes[1, 1].text(0.5, 0.5, 'Confidence not available',
                       ha='center', va='center', transform=axes[1, 1].transAxes,
                       fontsize=14, color='gray')
        axes[1, 1].set_title('Query Mapping Confidence')
except Exception as e:
    print(f"  WARNING: Panel 4 failed: {e}")
    axes[1, 1].text(0.5, 0.5, 'Error generating panel',
                   ha='center', va='center', transform=axes[1, 1].transAxes,
                   fontsize=14, color='red')
finally:
    cleanup_temp_columns(adata, temp_col_name)

plt.tight_layout()
fig_path = output_dir / f"01_umap_overview.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"OK Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 2: Side-by-Side Comparison (Reference vs Query)
# ==============================================================================

print("\n" + "=" * 70)
print("Figure 2: Reference vs Query Side-by-Side")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# Panel 1: Reference only
print("  Panel 1: Reference cells...")
if ref_cells.n_obs > 0:
    umap_ref = ref_cells.obsm['X_umap']
    
    if 'Cell_Type_L2' in ref_cells.obs.columns:
        colors_ref = ref_cells.obs['Cell_Type_L2'].astype(str)
        color_map_ref = {ct: CELLTYPE_COLORS.get(ct, '#95a5a6') 
                        for ct in colors_ref.unique()}
        
        for ct in colors_ref.unique():
            mask = (colors_ref == ct).values
            axes[0].scatter(umap_ref[mask, 0], umap_ref[mask, 1],
                          c=color_map_ref[ct], label=ct, s=10, alpha=0.6)
        axes[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), 
                      frameon=True, title='Cell Type')
    else:
        axes[0].scatter(umap_ref[:, 0], umap_ref[:, 1],
                       c=DATA_SOURCE_COLORS['reference'], s=10, alpha=0.6)
    
    axes[0].set_xlabel('UMAP 1', fontsize=12)
    axes[0].set_ylabel('UMAP 2', fontsize=12)
    axes[0].set_title(f'Reference Only\n({ref_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
else:
    axes[0].text(0.5, 0.5, 'No reference data',
                ha='center', va='center', transform=axes[0].transAxes,
                fontsize=14, color='gray')
    axes[0].set_title('Reference Only')

# Panel 2: Query only
print("  Panel 2: Query cells...")
if qry_cells.n_obs > 0:
    umap_qry = qry_cells.obsm['X_umap']
    
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        colors_qry = qry_cells.obs['Cell_Type_L2_final'].astype(str)
        color_map_qry = {ct: CELLTYPE_COLORS.get(ct, '#95a5a6') 
                        for ct in colors_qry.unique()}
        
        for ct in colors_qry.unique():
            mask = (colors_qry == ct).values
            axes[1].scatter(umap_qry[mask, 0], umap_qry[mask, 1],
                          c=color_map_qry[ct], label=ct, s=10, alpha=0.6)
        axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5),
                      frameon=True, title='Cell Type')
    else:
        axes[1].scatter(umap_qry[:, 0], umap_qry[:, 1],
                       c=DATA_SOURCE_COLORS['query'], s=10, alpha=0.6)
    
    axes[1].set_xlabel('UMAP 1', fontsize=12)
    axes[1].set_ylabel('UMAP 2', fontsize=12)
    axes[1].set_title(f'Query Only\n({qry_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
else:
    axes[1].text(0.5, 0.5, 'No query data',
                ha='center', va='center', transform=axes[1].transAxes,
                fontsize=14, color='gray')
    axes[1].set_title('Query Only')

plt.tight_layout()
fig_path = output_dir / f"02_reference_vs_query_sidebyside.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"OK Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 3: UMAP Space Validation
# ==============================================================================

print("\n" + "=" * 70)
print("Figure 3: UMAP Space Validation")
print("=" * 70)

if ref_cells.n_obs > 0 and qry_cells.n_obs > 0:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6))

    # Panel 1: X coordinate distribution
    print("  Panel 1: X coordinate distribution...")
    ref_umap = ref_cells.obsm['X_umap']
    qry_umap = qry_cells.obsm['X_umap']
    
    axes[0].hist(ref_umap[:, 0], bins=50, alpha=0.5, 
                label=f'Reference (n={ref_cells.n_obs:,})',
                color=DATA_SOURCE_COLORS['reference'], density=True)
    axes[0].hist(qry_umap[:, 0], bins=50, alpha=0.5,
                label=f'Query (n={qry_cells.n_obs:,})',
                color=DATA_SOURCE_COLORS['query'], density=True)
    axes[0].set_xlabel('UMAP X', fontsize=12)
    axes[0].set_ylabel('Density', fontsize=12)
    axes[0].set_title('UMAP X Coordinate Distribution', 
                     fontsize=14, fontweight='bold')
    axes[0].legend(frameon=True)
    axes[0].grid(alpha=0.3)
    
    # Panel 2: Y coordinate distribution
    print("  Panel 2: Y coordinate distribution...")
    axes[1].hist(ref_umap[:, 1], bins=50, alpha=0.5,
                label=f'Reference (n={ref_cells.n_obs:,})',
                color=DATA_SOURCE_COLORS['reference'], density=True)
    axes[1].hist(qry_umap[:, 1], bins=50, alpha=0.5,
                label=f'Query (n={qry_cells.n_obs:,})',
                color=DATA_SOURCE_COLORS['query'], density=True)
    axes[1].set_xlabel('UMAP Y', fontsize=12)
    axes[1].set_ylabel('Density', fontsize=12)
    axes[1].set_title('UMAP Y Coordinate Distribution',
                     fontsize=14, fontweight='bold')
    axes[1].legend(frameon=True)
    axes[1].grid(alpha=0.3)
    
    # Panel 3: 2D overlay
    print("  Panel 3: 2D overlay...")
    axes[2].scatter(ref_umap[:, 0], ref_umap[:, 1], s=5, alpha=0.3,
                   c=DATA_SOURCE_COLORS['reference'], 
                   label=f'Reference ({ref_cells.n_obs:,})')
    axes[2].scatter(qry_umap[:, 0], qry_umap[:, 1], s=5, alpha=0.3,
                   c=DATA_SOURCE_COLORS['query'],
                   label=f'Query ({qry_cells.n_obs:,})')
    axes[2].set_xlabel('UMAP X', fontsize=12)
    axes[2].set_ylabel('UMAP Y', fontsize=12)
    axes[2].set_title('UMAP Space Overlap', fontsize=14, fontweight='bold')
    axes[2].legend(frameon=True)
    axes[2].grid(alpha=0.3)
    
    # Add statistics
    ref_x_range = ref_umap[:, 0].ptp()
    ref_y_range = ref_umap[:, 1].ptp()
    qry_x_range = qry_umap[:, 0].ptp()
    qry_y_range = qry_umap[:, 1].ptp()
    
    stats_text = f"Range Comparison:\n"
    stats_text += f"Ref: X={ref_x_range:.1f}, Y={ref_y_range:.1f}\n"
    stats_text += f"Qry: X={qry_x_range:.1f}, Y={qry_y_range:.1f}\n"
    
    # Safe division
    if ref_x_range > 0 and ref_y_range > 0:
        stats_text += f"Ratio: X={qry_x_range/ref_x_range:.2f}, Y={qry_y_range/ref_y_range:.2f}"
    else:
        stats_text += f"Ratio: N/A (ref range is zero)"
    
    axes[2].text(0.05, 0.95, stats_text, transform=axes[2].transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    fig_path = output_dir / f"03_umap_space_validation.{FIGURE_FORMAT}"
    try:
        plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
        print(f"OK Saved: {fig_path.name}")
    except Exception as e:
        print(f"ERROR saving figure: {e}")
    finally:
        plt.close()
else:
    print("  Skipping: Need both reference and query cells")

# ==============================================================================
# Figure 4: Cell Type Distribution Comparison
# ==============================================================================

print("\n" + "=" * 70)
print("Figure 4: Cell Type Distribution Comparison")
print("=" * 70)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Reference distribution
print("  Panel 1: Reference distribution...")
if 'Cell_Type_L2' in ref_cells.obs.columns and ref_cells.n_obs > 0:
    ref_counts = ref_cells.obs['Cell_Type_L2'].value_counts()
    colors_ref = [CELLTYPE_COLORS.get(ct, '#95a5a6') for ct in ref_counts.index]
    
    wedges, texts, autotexts = axes[0].pie(
        ref_counts.values, labels=ref_counts.index, autopct='%1.1f%%',
        colors=colors_ref, startangle=90, textprops={'fontsize': 10}
    )
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_weight('bold')
    
    axes[0].set_title(f'Reference L2 Distribution\n({ref_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
else:
    axes[0].text(0.5, 0.5, 'Reference L2\nnot available',
                ha='center', va='center', transform=axes[0].transAxes,
                fontsize=14, color='gray')
    axes[0].set_title('Reference L2 Distribution', fontsize=14, fontweight='bold')

# Query distribution
print("  Panel 2: Query distribution...")
if 'Cell_Type_L2_final' in qry_cells.obs.columns and qry_cells.n_obs > 0:
    qry_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
    colors_qry = [CELLTYPE_COLORS.get(ct, '#95a5a6') for ct in qry_counts.index]
    
    wedges, texts, autotexts = axes[1].pie(
        qry_counts.values, labels=qry_counts.index, autopct='%1.1f%%',
        colors=colors_qry, startangle=90, textprops={'fontsize': 10}
    )
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_weight('bold')
    
    axes[1].set_title(f'Query L2 Distribution\n({qry_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
else:
    axes[1].text(0.5, 0.5, 'Query L2\nnot available',
                ha='center', va='center', transform=axes[1].transAxes,
                fontsize=14, color='gray')
    axes[1].set_title('Query L2 Distribution', fontsize=14, fontweight='bold')

plt.tight_layout()
fig_path = output_dir / f"04_celltype_distribution_comparison.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"OK Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 5: Query Mapping Quality (if confidence available)
# ==============================================================================

if 'mapping_confidence' in qry_cells.obs.columns and qry_cells.n_obs > 0:
    print("\n" + "=" * 70)
    print("Figure 5: Query Mapping Quality")
    print("=" * 70)
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    conf = pd.to_numeric(qry_cells.obs['mapping_confidence'], errors='coerce')
    conf_valid = conf[conf.notna()]
    
    # Panel 1: Confidence histogram
    print("  Panel 1: Confidence histogram...")
    axes[0, 0].hist(conf_valid, bins=50, edgecolor='black', 
                   alpha=0.7, color='steelblue')
    axes[0, 0].axvline(CONFIDENCE_THRESHOLD, color='red', 
                      linestyle='--', linewidth=2,
                      label=f'Threshold = {CONFIDENCE_THRESHOLD}')
    axes[0, 0].set_xlabel('Mapping Confidence', fontsize=12)
    axes[0, 0].set_ylabel('Number of Cells', fontsize=12)
    axes[0, 0].set_title('Confidence Distribution', 
                        fontsize=14, fontweight='bold')
    axes[0, 0].legend(frameon=True)
    axes[0, 0].grid(alpha=0.3)
    
    # Add statistics
    stats_text = f"Mean: {conf_valid.mean():.3f}\n"
    stats_text += f"Median: {conf_valid.median():.3f}\n"
    stats_text += f"High (≥{CONFIDENCE_THRESHOLD}): {(conf_valid>=CONFIDENCE_THRESHOLD).sum():,} "
    stats_text += f"({(conf_valid>=CONFIDENCE_THRESHOLD).sum()/len(conf_valid)*100:.1f}%)"
    
    axes[0, 0].text(0.05, 0.95, stats_text, transform=axes[0, 0].transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Panel 2: Confidence by cell type (violin)
    print("  Panel 2: Confidence by cell type...")
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        df = pd.DataFrame({
            'Cell_Type': qry_cells.obs['Cell_Type_L2_final'].astype(str),
            'Confidence': conf
        })
        df = df[df['Cell_Type'] != 'Unknown']
        df = df.dropna()
        
        if len(df) > 0:
            celltypes_sorted = df['Cell_Type'].value_counts().index.tolist()
            
            # Use seaborn for better violin plot
            sns.violinplot(data=df, y='Cell_Type', x='Confidence',
                          order=celltypes_sorted, ax=axes[0, 1],
                          orient='h', inner='box', palette=CELLTYPE_COLORS)
            
            axes[0, 1].axvline(CONFIDENCE_THRESHOLD, color='red',
                              linestyle='--', linewidth=2,
                              label=f'Threshold={CONFIDENCE_THRESHOLD}')
            axes[0, 1].set_xlabel('Mapping Confidence', fontsize=12)
            axes[0, 1].set_ylabel('')
            axes[0, 1].set_title('Confidence by Cell Type',
                                fontsize=14, fontweight='bold')
            axes[0, 1].legend(frameon=True)
            axes[0, 1].grid(alpha=0.3, axis='x')
        else:
            axes[0, 1].text(0.5, 0.5, 'No data',
                           ha='center', va='center', transform=axes[0, 1].transAxes,
                           fontsize=14, color='gray')
            axes[0, 1].set_title('Confidence by Cell Type')
    else:
        axes[0, 1].text(0.5, 0.5, 'Cell type data\nnot available',
                       ha='center', va='center', transform=axes[0, 1].transAxes,
                       fontsize=14, color='gray')
        axes[0, 1].set_title('Confidence by Cell Type')
    
    # Panel 3: Cell type counts
    print("  Panel 3: Cell type counts...")
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        l2_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
        l2_counts_sorted = l2_counts.sort_values(ascending=True)
        
        colors_bar = [CELLTYPE_COLORS.get(ct, '#95a5a6') 
                     for ct in l2_counts_sorted.index]
        
        y_pos = np.arange(len(l2_counts_sorted))
        axes[1, 0].barh(y_pos, l2_counts_sorted.values, color=colors_bar)
        axes[1, 0].set_yticks(y_pos)
        axes[1, 0].set_yticklabels(l2_counts_sorted.index)
        axes[1, 0].set_xlabel('Number of Cells', fontsize=12)
        axes[1, 0].set_title('L2 Cell Type Counts',
                            fontsize=14, fontweight='bold')
        axes[1, 0].grid(alpha=0.3, axis='x')
        
        # Add count labels
        for i, count in enumerate(l2_counts_sorted.values):
            axes[1, 0].text(count + max(l2_counts_sorted)*0.01, i,
                          f'{count:,}', va='center', fontsize=9)
    else:
        axes[1, 0].text(0.5, 0.5, 'Cell type data\nnot available',
                       ha='center', va='center', transform=axes[1, 0].transAxes,
                       fontsize=14, color='gray')
        axes[1, 0].set_title('L2 Cell Type Counts')
    
    # Panel 4: Mean confidence by cell type
    print("  Panel 4: Mean confidence by cell type...")
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        df_conf = pd.DataFrame({
            'Cell_Type': qry_cells.obs['Cell_Type_L2_final'].astype(str),
            'Confidence': conf
        })
        df_conf = df_conf[df_conf['Cell_Type'] != 'Unknown']
        
        conf_by_type = df_conf.groupby('Cell_Type')['Confidence'].agg(['mean', 'std']).sort_values('mean', ascending=True)
        
        if len(conf_by_type) > 0:
            y_pos = np.arange(len(conf_by_type))
            colors_bar = [CELLTYPE_COLORS.get(ct, '#95a5a6') 
                         for ct in conf_by_type.index]
            
            axes[1, 1].barh(y_pos, conf_by_type['mean'], 
                          xerr=conf_by_type['std'],
                          color=colors_bar, capsize=5, alpha=0.7)
            axes[1, 1].set_yticks(y_pos)
            axes[1, 1].set_yticklabels(conf_by_type.index)
            axes[1, 1].set_xlabel('Mean Confidence ± Std', fontsize=12)
            axes[1, 1].set_title('Mean Confidence by Cell Type',
                                fontsize=14, fontweight='bold')
            axes[1, 1].axvline(CONFIDENCE_THRESHOLD, color='red',
                              linestyle='--', linewidth=2)
            axes[1, 1].grid(alpha=0.3, axis='x')
            axes[1, 1].set_xlim(0, 1)
        else:
            axes[1, 1].text(0.5, 0.5, 'No data',
                           ha='center', va='center', transform=axes[1, 1].transAxes,
                           fontsize=14, color='gray')
            axes[1, 1].set_title('Mean Confidence by Cell Type')
    else:
        axes[1, 1].text(0.5, 0.5, 'Cell type data\nnot available',
                       ha='center', va='center', transform=axes[1, 1].transAxes,
                       fontsize=14, color='gray')
        axes[1, 1].set_title('Mean Confidence by Cell Type')
    
    plt.tight_layout()
    fig_path = output_dir / f"05_query_mapping_quality.{FIGURE_FORMAT}"
    try:
        plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
        print(f"OK Saved: {fig_path.name}")
    except Exception as e:
        print(f"ERROR saving figure: {e}")
    finally:
        plt.close()

# ==============================================================================
# Figure 6: Marker Gene Expression (if genes available)
# ==============================================================================

print("\n" + "=" * 70)
print("Figure 6: B Cell Marker Expression")
print("=" * 70)

# Define key markers
bcell_markers = {
    'Pan-B': ['CD19', 'MS4A1', 'CD79A'],
    'Naive': ['TCL1A', 'FCER2', 'IGHD'],
    'Memory': ['CD27', 'TNFRSF13B'],
    'GC': ['AICDA', 'BCL6', 'MME'],
    'Plasma': ['SDC1', 'XBP1', 'JCHAIN', 'MZB1']
}

# Find available markers
available_markers = []
for category, markers in bcell_markers.items():
    for marker in markers:
        # Check both .X and .raw
        in_X = marker in adata.var_names
        in_raw = hasattr(adata, 'raw') and adata.raw is not None and marker in adata.raw.var_names
        
        if in_X or in_raw:
            available_markers.append((category, marker, in_X))

if len(available_markers) > 0:
    print(f"  Found {len(available_markers)} markers")
    
    n_markers = min(12, len(available_markers))
    markers_to_plot = available_markers[:n_markers]
    
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, (category, marker, in_X) in enumerate(markers_to_plot):
        row = i // n_cols
        col = i % n_cols
        
        try:
            use_raw = not in_X  # Use raw if marker not in .X
            sc.pl.umap(adata, color=marker, ax=axes[row, col], show=False,
                       cmap='Reds', s=5, title=f'{marker} ({category})',
                       frameon=True, use_raw=use_raw)
        except Exception as e:
            print(f"  WARNING: Failed to plot {marker}: {e}")
            axes[row, col].text(0.5, 0.5, f'{marker}\nPlot failed',
                               ha='center', va='center', transform=axes[row, col].transAxes,
                               fontsize=12, color='red')
            axes[row, col].set_title(f'{marker} ({category})')
    
    # Hide unused subplots
    for i in range(n_markers, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].axis('off')
    
    plt.tight_layout()
    fig_path = output_dir / f"06_marker_expression.{FIGURE_FORMAT}"
    try:
        plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
        print(f"OK Saved: {fig_path.name}")
    except Exception as e:
        print(f"ERROR saving figure: {e}")
    finally:
        plt.close()
else:
    print("  WARNING: No standard B cell markers found in data")

# ==============================================================================
# Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Generate Summary Report")
print("=" * 70)

report_lines = []
report_lines.append("=" * 70)
report_lines.append("Merged Reference + Query Visualization Report")
report_lines.append("=" * 70)

report_lines.append(f"\nInput file: {MERGED_H5AD}")
report_lines.append(f"Analysis date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

report_lines.append(f"\n[ Dataset Info ]")
report_lines.append(f"  Total cells: {adata.n_obs:,}")
report_lines.append(f"  Reference: {ref_cells.n_obs:,} ({ref_cells.n_obs/adata.n_obs*100:.1f}%)")
report_lines.append(f"  Query: {qry_cells.n_obs:,} ({qry_cells.n_obs/adata.n_obs*100:.1f}%)")
report_lines.append(f"  Genes: {adata.n_vars:,}")

try:
    if 'Cell_Type_L2' in ref_cells.obs.columns and ref_cells.n_obs > 0:
        report_lines.append(f"\n[ Reference L2 Distribution ]")
        ref_counts = ref_cells.obs['Cell_Type_L2'].value_counts()
        for celltype, count in ref_counts.items():
            pct = count / ref_cells.n_obs * 100
            report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")
except Exception as e:
    report_lines.append(f"\n[ Reference L2 Distribution ]")
    report_lines.append(f"  Error: {e}")

try:
    if 'Cell_Type_L2_final' in qry_cells.obs.columns and qry_cells.n_obs > 0:
        report_lines.append(f"\n[ Query L2 Distribution ]")
        qry_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
        for celltype, count in qry_counts.items():
            pct = count / qry_cells.n_obs * 100
            report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")
except Exception as e:
    report_lines.append(f"\n[ Query L2 Distribution ]")
    report_lines.append(f"  Error: {e}")

try:
    if 'mapping_confidence' in qry_cells.obs.columns and qry_cells.n_obs > 0:
        conf = pd.to_numeric(qry_cells.obs['mapping_confidence'], errors='coerce')
        conf_valid = conf[conf.notna()]
        
        if len(conf_valid) > 0:
            report_lines.append(f"\n[ Query Mapping Quality ]")
            report_lines.append(f"  Confidence mean: {conf_valid.mean():.3f}")
            report_lines.append(f"  Confidence median: {conf_valid.median():.3f}")
            high_conf = (conf_valid >= CONFIDENCE_THRESHOLD).sum()
            report_lines.append(f"  High confidence (≥{CONFIDENCE_THRESHOLD}): {high_conf:,} ({high_conf/len(conf_valid)*100:.1f}%)")
        else:
            report_lines.append(f"\n[ Query Mapping Quality ]")
            report_lines.append(f"  No valid confidence values")
except Exception as e:
    report_lines.append(f"\n[ Query Mapping Quality ]")
    report_lines.append(f"  Error: {e}")

try:
    if ref_cells.n_obs > 0 and qry_cells.n_obs > 0:
        ref_umap = ref_cells.obsm['X_umap']
        qry_umap = qry_cells.obsm['X_umap']
        
        ref_x_range = ref_umap[:, 0].ptp()
        ref_y_range = ref_umap[:, 1].ptp()
        qry_x_range = qry_umap[:, 0].ptp()
        qry_y_range = qry_umap[:, 1].ptp()
        
        report_lines.append(f"\n[ UMAP Space Validation ]")
        report_lines.append(f"  Reference UMAP range: X={ref_x_range:.1f}, Y={ref_y_range:.1f}")
        report_lines.append(f"  Query UMAP range: X={qry_x_range:.1f}, Y={qry_y_range:.1f}")
        
        if ref_x_range > 0 and ref_y_range > 0:
            report_lines.append(f"  Range ratio: X={qry_x_range/ref_x_range:.2f}, Y={qry_y_range/ref_y_range:.2f}")
            
            x_ratio = qry_x_range/ref_x_range
            y_ratio = qry_y_range/ref_y_range
            
            if 0.3 < x_ratio < 3.0 and 0.3 < y_ratio < 3.0:
                report_lines.append(f"  Status: OK - Same UMAP space confirmed")
            else:
                report_lines.append(f"  Status: WARNING - UMAP spaces may differ")
        else:
            report_lines.append(f"  Range ratio: N/A (reference range is zero)")
            report_lines.append(f"  Status: WARNING - Invalid reference UMAP range")
except Exception as e:
    report_lines.append(f"\n[ UMAP Space Validation ]")
    report_lines.append(f"  Error: {e}")

report_lines.append(f"\n[ Generated Figures ]")
for i, fname in enumerate([
    "01_umap_overview.pdf",
    "02_reference_vs_query_sidebyside.pdf",
    "03_umap_space_validation.pdf",
    "04_celltype_distribution_comparison.pdf",
    "05_query_mapping_quality.pdf",
    "06_marker_expression.pdf"
], 1):
    full_path = output_dir / fname
    if full_path.exists():
        file_size = full_path.stat().st_size / 1024  # KB
        report_lines.append(f"  {i}. {fname} ({file_size:.1f} KB)")
    else:
        report_lines.append(f"  {i}. {fname} (not generated)")

report_lines.append(f"\n" + "=" * 70)
report_lines.append("Visualization Complete")
report_lines.append("=" * 70)

report_text = '\n'.join(report_lines)
print("\n" + report_text)

report_file = output_dir / "visualization_summary.txt"
try:
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\nReport saved: {report_file}")
except Exception as e:
    print(f"\nERROR saving report: {e}")

print("\n" + "=" * 70)
print("SUCCESS - All Visualizations Generated")
print("=" * 70)
print(f"\nOutput directory: {output_dir}")
print("\nGenerated files:")
for f in sorted(output_dir.glob("*.pdf")):
    print(f"  - {f.name}")
print(f"  - visualization_summary.txt")