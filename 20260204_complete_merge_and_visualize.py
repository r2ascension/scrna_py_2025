#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Complete Reference + Query Merge and Visualization
===================================================
1. Merges reference and query
2. Generates comprehensive visualizations

Usage:
    python complete_merge_and_visualize.py
"""

import sys
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
from scipy.sparse import csr_matrix, issparse

print("=" * 70)
print("Complete Reference + Query Merge and Visualization")
print("=" * 70)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

REF_H5AD = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/reference_with_L2_umap.h5ad"
QUERY_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/query_mapped_L2.h5ad"
OUTPUT_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/reference_plus_query_merged_L2.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/merged_analysis"

DPI = 300
FIGURE_FORMAT = 'pdf'

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)

print(f"\nConfiguration:")
print(f"  Reference: {REF_H5AD}")
print(f"  Query: {QUERY_H5AD}")
print(f"  Output merged: {OUTPUT_H5AD}")
print(f"  Output figures: {OUTPUT_DIR}")

# ==============================================================================
# PART 1: MERGE
# ==============================================================================

print("\n" + "=" * 70)
print("PART 1: Merge Reference and Query")
print("=" * 70)

# Load Reference
print("\nStep 1: Load Reference")
print("-" * 40)
adata_ref = sc.read_h5ad(REF_H5AD)
print(f"  Shape: {adata_ref.shape}")

if 'X_umap' not in adata_ref.obsm:
    print("ERROR: Reference missing X_umap!")
    sys.exit(1)

# Load Query
print("\nStep 2: Load Query")
print("-" * 40)
adata_query = sc.read_h5ad(QUERY_H5AD)
print(f"  Shape: {adata_query.shape}")

if 'X_umap' not in adata_query.obsm:
    print("ERROR: Query missing X_umap!")
    sys.exit(1)

# Pre-Concat Cleanup
print("\nStep 3: Pre-Concat Cleanup")
print("-" * 40)

# Clean var
var_cols_to_drop = ['highly_variable', 'highly_variable_rank',
                    'means', 'dispersions', 'dispersions_norm',
                    'mt', 'n_cells', 'n_counts']

for col in var_cols_to_drop:
    if col in adata_ref.var.columns:
        adata_ref.var.drop(columns=[col], inplace=True)
    if col in adata_query.var.columns:
        adata_query.var.drop(columns=[col], inplace=True)

print("  OK var cleaned")

# Clear uns
adata_ref.uns = {}
adata_query.uns = {}
print("  OK uns cleared")

# Unify sparse format
if issparse(adata_ref.X):
    adata_ref.X = csr_matrix(adata_ref.X)
if issparse(adata_query.X):
    adata_query.X = csr_matrix(adata_query.X)

for layer_key in list(adata_ref.layers.keys()):
    if issparse(adata_ref.layers[layer_key]):
        adata_ref.layers[layer_key] = csr_matrix(adata_ref.layers[layer_key])

for layer_key in list(adata_query.layers.keys()):
    if issparse(adata_query.layers[layer_key]):
        adata_query.layers[layer_key] = csr_matrix(adata_query.layers[layer_key])

print("  OK sparse format unified")

# Align obsm
critical_obsm = {'X_umap', 'X_scANVI_L2'}

for key in list(adata_ref.obsm.keys()):
    if key not in critical_obsm:
        del adata_ref.obsm[key]

for key in list(adata_query.obsm.keys()):
    if key not in critical_obsm:
        del adata_query.obsm[key]

print("  OK obsm aligned")

# Remove obsp
for key in list(adata_ref.obsp.keys()):
    del adata_ref.obsp[key]
for key in list(adata_query.obsp.keys()):
    del adata_query.obsp[key]

print("  OK obsp removed")

# Check uniqueness
if not adata_ref.obs_names.is_unique:
    adata_ref.obs_names_make_unique()
if not adata_query.obs_names.is_unique:
    adata_query.obs_names_make_unique()

print("  OK obs_names unique")

# Add prefix
adata_ref.obs_names = [f"ref::{x}" for x in adata_ref.obs_names]
adata_query.obs_names = [f"qry::{x}" for x in adata_query.obs_names]

print("  OK prefix added")

# Concatenate
print("\nStep 4: Concatenate")
print("-" * 40)

adata_all = sc.concat(
    {"reference": adata_ref, "query": adata_query},
    axis=0,
    join="outer",
    merge="unique",
    fill_value=0,
    label="data_source"
)

print(f"OK Merged shape: {adata_all.shape}")
print(f"  Reference: {(adata_all.obs['data_source']=='reference').sum():,}")
print(f"  Query: {(adata_all.obs['data_source']=='query').sum():,}")

if 'X_umap' not in adata_all.obsm:
    print("WARNING: X_umap lost!")
else:
    print("OK X_umap preserved")

# Save
print("\nStep 5: Save Merged Data")
print("-" * 40)

output_path = Path(OUTPUT_H5AD)
output_path.parent.mkdir(parents=True, exist_ok=True)

try:
    adata_all.write_h5ad(output_path, compression='gzip')
    file_size = output_path.stat().st_size / 1024**3
    print(f"OK Saved ({file_size:.2f} GB)")
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)

# Cleanup reference/query (keep merged)
del adata_ref, adata_query
gc.collect()

print("\nOK Merge complete")

# ==============================================================================
# PART 2: VISUALIZATION
# ==============================================================================

print("\n" + "=" * 70)
print("PART 2: Generate Visualizations")
print("=" * 70)

# Figure 1: UMAP Overview (2x3 grid)
print("\n1. UMAP Overview (6 panels)...")
print("-" * 40)

fig, axes = plt.subplots(2, 3, figsize=(20, 13))

# Panel 1: Data source
sc.pl.umap(adata_all, color='data_source', ax=axes[0, 0], show=False,
           title='Data Source (Ref vs Query)', s=5,
           palette={'reference': '#1f77b4', 'query': '#ff7f0e'})

# Panel 2: Reference L2 (if exists)
if 'Cell_Type_L2' in adata_all.obs.columns:
    sc.pl.umap(adata_all, color='Cell_Type_L2', ax=axes[0, 1], show=False,
               title='Reference L2 Labels', legend_loc='right margin', s=5)
else:
    axes[0, 1].text(0.5, 0.5, 'Reference L2\nnot available',
                    ha='center', va='center', transform=axes[0, 1].transAxes,
                    fontsize=14)
    axes[0, 1].set_title('Reference L2 Labels')

# Panel 3: Query L2 final
if 'Cell_Type_L2_final' in adata_all.obs.columns:
    sc.pl.umap(adata_all, color='Cell_Type_L2_final', ax=axes[0, 2], show=False,
               title='Query L2 Predictions', legend_loc='right margin', s=5)
else:
    axes[0, 2].text(0.5, 0.5, 'Query L2\nnot available',
                    ha='center', va='center', transform=axes[0, 2].transAxes,
                    fontsize=14)
    axes[0, 2].set_title('Query L2 Predictions')

# Panel 4: Mapping confidence (query only)
if 'mapping_confidence' in adata_all.obs.columns:
    query_cells = adata_all[adata_all.obs['data_source'] == 'query']
    
    # Create a full array with NaN for reference cells
    conf_array = np.full(adata_all.n_obs, np.nan)
    query_idx = np.where(adata_all.obs['data_source'] == 'query')[0]
    conf_array[query_idx] = query_cells.obs['mapping_confidence'].values
    
    sc.pl.umap(adata_all, color=conf_array, ax=axes[1, 0], show=False,
               title='Mapping Confidence (Query Only)', cmap='viridis',
               vmin=0, vmax=1, s=5)
else:
    axes[1, 0].text(0.5, 0.5, 'Confidence\nnot available',
                    ha='center', va='center', transform=axes[1, 0].transAxes,
                    fontsize=14)
    axes[1, 0].set_title('Mapping Confidence')

# Panel 5: Pan-B marker (CD19)
if 'CD19' in adata_all.var_names:
    sc.pl.umap(adata_all, color='CD19', ax=axes[1, 1], show=False,
               title='CD19 (Pan-B Marker)', cmap='Reds', s=5)
else:
    axes[1, 1].text(0.5, 0.5, 'CD19\nnot found',
                    ha='center', va='center', transform=axes[1, 1].transAxes,
                    fontsize=14)
    axes[1, 1].set_title('CD19 (Pan-B Marker)')

# Panel 6: Plasma marker (SDC1/CD138)
plasma_marker = None
for marker in ['SDC1', 'CD138']:
    if marker in adata_all.var_names:
        plasma_marker = marker
        break

if plasma_marker:
    sc.pl.umap(adata_all, color=plasma_marker, ax=axes[1, 2], show=False,
               title=f'{plasma_marker} (Plasma Marker)', cmap='Reds', s=5)
else:
    axes[1, 2].text(0.5, 0.5, 'Plasma marker\nnot found',
                    ha='center', va='center', transform=axes[1, 2].transAxes,
                    fontsize=14)
    axes[1, 2].set_title('Plasma Marker')

plt.tight_layout()
fig_path = output_dir / f"merged_umap_overview.{FIGURE_FORMAT}"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_path.name}")

# Figure 2: Cell Type Distribution Comparison
print("\n2. Cell Type Distribution (Ref vs Query)...")
print("-" * 40)

# Get L2 distributions
ref_cells = adata_all[adata_all.obs['data_source'] == 'reference']
qry_cells = adata_all[adata_all.obs['data_source'] == 'query']

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Reference distribution
if 'Cell_Type_L2' in ref_cells.obs.columns:
    ref_counts = ref_cells.obs['Cell_Type_L2'].value_counts()
    colors_ref = plt.cm.Set3(np.linspace(0, 1, len(ref_counts)))
    
    axes[0].pie(ref_counts.values, labels=ref_counts.index, autopct='%1.1f%%',
                colors=colors_ref, startangle=90)
    axes[0].set_title(f'Reference L2 Distribution\n({ref_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
else:
    axes[0].text(0.5, 0.5, 'Reference L2\nnot available',
                ha='center', va='center', transform=axes[0].transAxes, fontsize=14)
    axes[0].set_title('Reference L2 Distribution', fontsize=14, fontweight='bold')

# Query distribution
if 'Cell_Type_L2_final' in qry_cells.obs.columns:
    qry_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
    colors_qry = plt.cm.Set3(np.linspace(0, 1, len(qry_counts)))
    
    axes[1].pie(qry_counts.values, labels=qry_counts.index, autopct='%1.1f%%',
                colors=colors_qry, startangle=90)
    axes[1].set_title(f'Query L2 Distribution\n({qry_cells.n_obs:,} cells)',
                     fontsize=14, fontweight='bold')
else:
    axes[1].text(0.5, 0.5, 'Query L2\nnot available',
                ha='center', va='center', transform=axes[1].transAxes, fontsize=14)
    axes[1].set_title('Query L2 Distribution', fontsize=14, fontweight='bold')

plt.tight_layout()
fig_path = output_dir / f"celltype_comparison.{FIGURE_FORMAT}"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_path.name}")

# Figure 3: UMAP Coordinate Validation (Same Space Check)
print("\n3. UMAP Space Validation...")
print("-" * 40)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# UMAP coordinates
ref_umap = ref_cells.obsm['X_umap']
qry_umap = qry_cells.obsm['X_umap']

# Panel 1: X coordinate distribution
axes[0].hist(ref_umap[:, 0], bins=50, alpha=0.5, label='Reference', color='blue', density=True)
axes[0].hist(qry_umap[:, 0], bins=50, alpha=0.5, label='Query', color='red', density=True)
axes[0].set_xlabel('UMAP X', fontsize=12)
axes[0].set_ylabel('Density', fontsize=12)
axes[0].set_title('UMAP X Coordinate Distribution', fontsize=14, fontweight='bold')
axes[0].legend()
axes[0].grid(alpha=0.3)

# Panel 2: Y coordinate distribution
axes[1].hist(ref_umap[:, 1], bins=50, alpha=0.5, label='Reference', color='blue', density=True)
axes[1].hist(qry_umap[:, 1], bins=50, alpha=0.5, label='Query', color='red', density=True)
axes[1].set_xlabel('UMAP Y', fontsize=12)
axes[1].set_ylabel('Density', fontsize=12)
axes[1].set_title('UMAP Y Coordinate Distribution', fontsize=14, fontweight='bold')
axes[1].legend()
axes[1].grid(alpha=0.3)

# Panel 3: 2D scatter with alpha
axes[2].scatter(ref_umap[:, 0], ref_umap[:, 1], s=5, alpha=0.3, 
               c='blue', label=f'Reference ({ref_cells.n_obs:,})')
axes[2].scatter(qry_umap[:, 0], qry_umap[:, 1], s=5, alpha=0.3, 
               c='red', label=f'Query ({qry_cells.n_obs:,})')
axes[2].set_xlabel('UMAP X', fontsize=12)
axes[2].set_ylabel('UMAP Y', fontsize=12)
axes[2].set_title('UMAP Coordinate Overlap', fontsize=14, fontweight='bold')
axes[2].legend()
axes[2].grid(alpha=0.3)

# Add range statistics
ref_x_range = ref_umap[:, 0].ptp()
ref_y_range = ref_umap[:, 1].ptp()
qry_x_range = qry_umap[:, 0].ptp()
qry_y_range = qry_umap[:, 1].ptp()

stats_text = f"Range comparison:\n"
stats_text += f"Ref X: {ref_x_range:.1f}, Y: {ref_y_range:.1f}\n"
stats_text += f"Qry X: {qry_x_range:.1f}, Y: {qry_y_range:.1f}\n"
stats_text += f"Ratio X: {qry_x_range/ref_x_range:.2f}, Y: {qry_y_range/ref_y_range:.2f}"

axes[2].text(0.05, 0.95, stats_text, transform=axes[2].transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
fig_path = output_dir / f"umap_space_validation.{FIGURE_FORMAT}"
plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
plt.close()
print(f"  Saved: {fig_path.name}")

# Figure 4: Query Mapping Quality
if 'mapping_confidence' in qry_cells.obs.columns:
    print("\n4. Query Mapping Quality...")
    print("-" * 40)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Panel 1: Confidence distribution
    conf = qry_cells.obs['mapping_confidence'].values
    axes[0, 0].hist(conf, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
    axes[0, 0].axvline(0.5, color='red', linestyle='--', linewidth=2, label='Threshold=0.5')
    axes[0, 0].set_xlabel('Mapping Confidence', fontsize=12)
    axes[0, 0].set_ylabel('Number of Cells', fontsize=12)
    axes[0, 0].set_title('Confidence Distribution', fontsize=14, fontweight='bold')
    axes[0, 0].legend()
    axes[0, 0].grid(alpha=0.3)
    
    # Add stats
    stats_text = f"Mean: {conf.mean():.3f}\n"
    stats_text += f"Median: {np.median(conf):.3f}\n"
    stats_text += f"High (>=0.5): {(conf>=0.5).sum():,} ({(conf>=0.5).sum()/len(conf)*100:.1f}%)"
    axes[0, 0].text(0.05, 0.95, stats_text, transform=axes[0, 0].transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Panel 2: Confidence by cell type
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        df = qry_cells.obs[['Cell_Type_L2_final', 'mapping_confidence']].copy()
        df = df[df['Cell_Type_L2_final'] != 'Unknown']
        
        celltypes_sorted = df['Cell_Type_L2_final'].value_counts().index.tolist()
        
        sns.violinplot(data=df, y='Cell_Type_L2_final', x='mapping_confidence',
                      order=celltypes_sorted, ax=axes[0, 1], orient='h')
        axes[0, 1].axvline(0.5, color='red', linestyle='--', linewidth=2)
        axes[0, 1].set_xlabel('Mapping Confidence', fontsize=12)
        axes[0, 1].set_ylabel('L2 Cell Type', fontsize=12)
        axes[0, 1].set_title('Confidence by Cell Type', fontsize=14, fontweight='bold')
        axes[0, 1].grid(axis='x', alpha=0.3)
    
    # Panel 3: Cell type counts
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        l2_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
        l2_counts_sorted = l2_counts.sort_values(ascending=True)
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(l2_counts_sorted)))
        
        axes[1, 0].barh(range(len(l2_counts_sorted)), l2_counts_sorted.values, color=colors)
        axes[1, 0].set_yticks(range(len(l2_counts_sorted)))
        axes[1, 0].set_yticklabels(l2_counts_sorted.index)
        axes[1, 0].set_xlabel('Number of Cells', fontsize=12)
        axes[1, 0].set_title('L2 Cell Type Counts', fontsize=14, fontweight='bold')
        axes[1, 0].grid(axis='x', alpha=0.3)
        
        for i, count in enumerate(l2_counts_sorted.values):
            axes[1, 0].text(count + max(l2_counts_sorted)*0.01, i, f'{count:,}',
                          va='center', fontsize=9)
    
    # Panel 4: Confidence vs cell count scatter
    if 'Cell_Type_L2_final' in qry_cells.obs.columns:
        conf_by_type = qry_cells.obs.groupby('Cell_Type_L2_final')['mapping_confidence'].agg(['mean', 'count'])
        conf_by_type = conf_by_type[conf_by_type.index != 'Unknown']
        
        axes[1, 1].scatter(conf_by_type['count'], conf_by_type['mean'], s=100, alpha=0.6)
        
        for celltype in conf_by_type.index:
            axes[1, 1].text(conf_by_type.loc[celltype, 'count'],
                          conf_by_type.loc[celltype, 'mean'],
                          celltype, fontsize=8, ha='right')
        
        axes[1, 1].set_xlabel('Cell Count', fontsize=12)
        axes[1, 1].set_ylabel('Mean Confidence', fontsize=12)
        axes[1, 1].set_title('Confidence vs Cell Count', fontsize=14, fontweight='bold')
        axes[1, 1].axhline(0.5, color='red', linestyle='--', linewidth=1)
        axes[1, 1].grid(alpha=0.3)
        axes[1, 1].set_xscale('log')
    
    plt.tight_layout()
    fig_path = output_dir / f"query_mapping_quality.{FIGURE_FORMAT}"
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fig_path.name}")

# ==============================================================================
# Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Generate Summary Report")
print("=" * 70)

report_lines = []
report_lines.append("=" * 70)
report_lines.append("Merged Reference + Query Analysis Report")
report_lines.append("=" * 70)

report_lines.append(f"\nAnalysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

report_lines.append(f"\n[ Dataset Info ]")
report_lines.append(f"  Total cells: {adata_all.n_obs:,}")
report_lines.append(f"  Reference: {ref_cells.n_obs:,} ({ref_cells.n_obs/adata_all.n_obs*100:.1f}%)")
report_lines.append(f"  Query: {qry_cells.n_obs:,} ({qry_cells.n_obs/adata_all.n_obs*100:.1f}%)")
report_lines.append(f"  Genes: {adata_all.n_vars:,}")

if 'Cell_Type_L2' in ref_cells.obs.columns:
    report_lines.append(f"\n[ Reference L2 Distribution ]")
    ref_counts = ref_cells.obs['Cell_Type_L2'].value_counts()
    for celltype, count in ref_counts.items():
        pct = count / ref_cells.n_obs * 100
        report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")

if 'Cell_Type_L2_final' in qry_cells.obs.columns:
    report_lines.append(f"\n[ Query L2 Distribution ]")
    qry_counts = qry_cells.obs['Cell_Type_L2_final'].value_counts()
    for celltype, count in qry_counts.items():
        pct = count / qry_cells.n_obs * 100
        report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")

if 'mapping_confidence' in qry_cells.obs.columns:
    conf = qry_cells.obs['mapping_confidence']
    report_lines.append(f"\n[ Query Mapping Quality ]")
    report_lines.append(f"  Confidence mean: {conf.mean():.3f}")
    report_lines.append(f"  Confidence median: {np.median(conf):.3f}")
    high_conf = (conf >= 0.5).sum()
    report_lines.append(f"  High confidence (>=0.5): {high_conf:,} ({high_conf/len(conf)*100:.1f}%)")

report_lines.append(f"\n[ UMAP Space Validation ]")
report_lines.append(f"  Reference UMAP range: X={ref_x_range:.1f}, Y={ref_y_range:.1f}")
report_lines.append(f"  Query UMAP range: X={qry_x_range:.1f}, Y={qry_y_range:.1f}")
report_lines.append(f"  Range ratio: X={qry_x_range/ref_x_range:.2f}, Y={qry_y_range/ref_y_range:.2f}")

if 0.3 < qry_x_range/ref_x_range < 3.0 and 0.3 < qry_y_range/ref_y_range < 3.0:
    report_lines.append(f"  OK Same UMAP space confirmed")
else:
    report_lines.append(f"  WARNING: UMAP spaces may differ")

report_lines.append(f"\n[ Generated Files ]")
report_lines.append(f"  - {OUTPUT_H5AD}")
report_lines.append(f"  - merged_umap_overview.{FIGURE_FORMAT}")
report_lines.append(f"  - celltype_comparison.{FIGURE_FORMAT}")
report_lines.append(f"  - umap_space_validation.{FIGURE_FORMAT}")
if 'mapping_confidence' in qry_cells.obs.columns:
    report_lines.append(f"  - query_mapping_quality.{FIGURE_FORMAT}")

report_lines.append(f"\n" + "=" * 70)
report_lines.append("Analysis Complete")
report_lines.append("=" * 70)

report_text = '\n'.join(report_lines)
print(report_text)

report_file = output_dir / "merged_analysis_summary.txt"
with open(report_file, 'w') as f:
    f.write(report_text)

print(f"\nReport saved: {report_file}")

print("\n" + "=" * 70)
print("SUCCESS - Merge and Visualization Complete")
print("=" * 70)
print(f"\nMerged data: {OUTPUT_H5AD}")
print(f"Figures: {OUTPUT_DIR}")
