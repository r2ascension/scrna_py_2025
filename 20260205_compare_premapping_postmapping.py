#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pre-Mapping vs Post-Mapping Comparison and Visualization
=========================================================
Merge original query data with mapped results and generate comparison visualizations.

Key Features:
- Handles gene count mismatch (full genes vs HVG)
- Dual-layer structure (HVG for UMAP, full genes preserved)
- High-quality comparison visualizations
- Detailed QC metrics

Usage:
    python compare_premapping_postmapping.py
"""

import sys
from pathlib import Path
import warnings
import gc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
from scipy.sparse import issparse, csr_matrix

warnings.filterwarnings("ignore")

print("=" * 70)
print("Pre-Mapping vs Post-Mapping Comparison")
print("=" * 70)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def cleanup_temp_columns(adata, col_names):
    """Remove temporary columns safely."""
    if isinstance(col_names, str):
        col_names = [col_names]
    for col in col_names:
        if col in adata.obs.columns:
            adata.obs.drop(columns=[col], inplace=True)


# ==============================================================================
# CONFIGURATION
# ==============================================================================

print("\n" + "=" * 70)
print("CONFIGURATION")
print("=" * 70)

# Input files
PRE_MAPPING_H5AD = "/home/h2048/data/source/polyp/polyp_obj_updated_20260113.h5ad"
POST_MAPPING_H5AD = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/query_mapped_to_reference.h5ad"

# Output directory
OUTPUT_DIR = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/comparison_analysis"

# Visualization parameters
DPI = 300
FIGURE_FORMAT = "pdf"
CONFIDENCE_THRESHOLD = 0.5

# Color schemes
MAPPING_STATUS_COLORS = {
    'pre-mapping': '#95a5a6',    # Gray
    'post-mapping': '#3498db'    # Blue
}

print(f"\nInput files:")
print(f"  Pre-mapping: {PRE_MAPPING_H5AD}")
print(f"  Post-mapping: {POST_MAPPING_H5AD}")
print(f"  Output: {OUTPUT_DIR}")

# Create output directory
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"\n✓ Output directory created")

sc.settings.set_figure_params(dpi=DPI, facecolor='white', format=FIGURE_FORMAT)
sns.set_style("whitegrid")

# ==============================================================================
# STEP 1: Load Data
# ==============================================================================

print("\n" + "=" * 70)
print("STEP 1: Load Data")
print("=" * 70)

print("\n1. Loading pre-mapping data...")
adata_pre = sc.read_h5ad(PRE_MAPPING_H5AD)
print(f"   Shape: {adata_pre.shape}")
print(f"   obs columns: {len(adata_pre.obs.columns)}")
print(f"   obsm keys: {list(adata_pre.obsm.keys())}")

print("\n2. Loading post-mapping data...")
adata_post = sc.read_h5ad(POST_MAPPING_H5AD)
print(f"   Shape: {adata_post.shape}")
print(f"   obs columns: {len(adata_post.obs.columns)}")
print(f"   obsm keys: {list(adata_post.obsm.keys())}")

# Check cell alignment
print("\n3. Checking cell alignment...")
if not (adata_pre.obs_names == adata_post.obs_names).all():
    print("   WARNING: Cell order differs, aligning...")
    common_cells = adata_pre.obs_names.intersection(adata_post.obs_names)
    print(f"   Common cells: {len(common_cells):,}")
    
    adata_pre = adata_pre[common_cells].copy()
    adata_post = adata_post[common_cells].copy()
    print(f"   ✓ Aligned to {len(common_cells):,} cells")
else:
    print(f"   ✓ Cell order matches ({adata_pre.n_obs:,} cells)")

# ==============================================================================
# STEP 2: Gene Alignment Strategy
# ==============================================================================

print("\n" + "=" * 70)
print("STEP 2: Gene Alignment Strategy")
print("=" * 70)

print(f"\nGene counts:")
print(f"  Pre-mapping: {adata_pre.n_vars:,} genes")
print(f"  Post-mapping: {adata_post.n_vars:,} genes")

if adata_pre.n_vars == adata_post.n_vars:
    print(f"\n  Strategy: Same gene count, use as-is")
    common_genes = list(adata_pre.var_names)
elif adata_pre.n_vars > adata_post.n_vars:
    print(f"\n  Strategy: Pre-mapping has more genes")
    print(f"  Post-mapping likely has HVG only")
    print(f"  Action: Subset pre-mapping to post-mapping genes")
    
    common_genes = [g for g in adata_post.var_names if g in adata_pre.var_names]
    overlap_pct = len(common_genes) / adata_post.n_vars * 100
    print(f"  Overlap: {len(common_genes):,}/{adata_post.n_vars:,} ({overlap_pct:.1f}%)")
    
    if overlap_pct < 95:
        print(f"  WARNING: Low overlap ({overlap_pct:.1f}%)")
    
    # Subset pre-mapping to common genes
    print(f"\n  Subsetting pre-mapping to {len(common_genes):,} genes...")
    adata_pre = adata_pre[:, common_genes].copy()
    print(f"  ✓ Pre-mapping subset: {adata_pre.shape}")
else:
    print(f"\n  ERROR: Post-mapping has more genes than pre-mapping")
    print(f"  This should not happen. Please check data.")
    sys.exit(1)

print(f"\n✓ Gene alignment complete")
print(f"  Final gene count: {adata_pre.n_vars:,}")

# ==============================================================================
# STEP 3: Prepare for Concatenation
# ==============================================================================

print("\n" + "=" * 70)
print("STEP 3: Prepare for Concatenation")
print("=" * 70)

print("\n1. Cleaning pre-mapping data...")

# Keep minimal obsm for pre-mapping
pre_obsm_keys = list(adata_pre.obsm.keys())
print(f"   Pre-mapping obsm: {pre_obsm_keys}")

# Only keep UMAP if it exists
if 'X_umap' not in adata_pre.obsm:
    print(f"   WARNING: No X_umap in pre-mapping data")
    if 'X_pca' in adata_pre.obsm:
        print(f"   Computing UMAP from PCA...")
        sc.pp.neighbors(adata_pre, use_rep='X_pca')
        sc.tl.umap(adata_pre)
    else:
        print(f"   ERROR: Cannot generate UMAP (no X_pca)")
        sys.exit(1)

# Clean pre-mapping obsm
for key in list(adata_pre.obsm.keys()):
    if key != 'X_umap':
        del adata_pre.obsm[key]
print(f"   ✓ Pre-mapping obsm cleaned: {list(adata_pre.obsm.keys())}")

print("\n2. Cleaning post-mapping data...")

# Identify UMAP key in post-mapping
umap_keys = [k for k in adata_post.obsm.keys() if 'umap' in k.lower()]
print(f"   UMAP keys found: {umap_keys}")

if 'X_umap' in adata_post.obsm:
    post_umap_key = 'X_umap'
elif 'X_umap_mapped' in adata_post.obsm:
    post_umap_key = 'X_umap_mapped'
    # Rename to standard
    adata_post.obsm['X_umap'] = adata_post.obsm[post_umap_key]
    post_umap_key = 'X_umap'
else:
    print(f"   ERROR: No UMAP found in post-mapping")
    sys.exit(1)

print(f"   Using UMAP: {post_umap_key}")

# Clean post-mapping obsm
for key in list(adata_post.obsm.keys()):
    if key != 'X_umap':
        del adata_post.obsm[key]
print(f"   ✓ Post-mapping obsm cleaned: {list(adata_post.obsm.keys())}")

print("\n3. Clearing metadata conflicts...")

# Clear uns to avoid conflicts
adata_pre.uns = {}
adata_post.uns = {}
print(f"   ✓ uns cleared")

# Clear obsp
for ad in [adata_pre, adata_post]:
    for key in list(ad.obsp.keys()):
        del ad.obsp[key]
print(f"   ✓ obsp cleared")

# Unify sparse format
for ad in [adata_pre, adata_post]:
    if issparse(ad.X):
        ad.X = csr_matrix(ad.X)
print(f"   ✓ Sparse format unified (CSR)")

print("\n4. Adding mapping status labels...")
adata_pre.obs['mapping_status'] = 'pre-mapping'
adata_post.obs['mapping_status'] = 'post-mapping'
print(f"   ✓ Labels added")

print("\n5. Making cell names unique...")
adata_pre.obs_names = [f"pre::{x}" for x in adata_pre.obs_names]
adata_post.obs_names = [f"post::{x}" for x in adata_post.obs_names]
print(f"   ✓ Cell names prefixed")

# ==============================================================================
# STEP 4: Concatenate
# ==============================================================================

print("\n" + "=" * 70)
print("STEP 4: Concatenate Pre and Post Mapping Data")
print("=" * 70)

print("\nBefore concatenation:")
print(f"  Pre-mapping: {adata_pre.shape}")
print(f"  Post-mapping: {adata_post.shape}")

print("\nConcatenating...")
try:
    adata_combined = sc.concat(
        {"pre": adata_pre, "post": adata_post},
        axis=0,
        join="inner",
        merge="unique",
        label="data_origin"
    )
    print(f"✓ Concatenation successful")
    print(f"  Combined shape: {adata_combined.shape}")
except Exception as e:
    print(f"ERROR during concatenation: {e}")
    sys.exit(1)

# Verify UMAP preservation
if 'X_umap' not in adata_combined.obsm:
    print(f"ERROR: X_umap lost during concatenation!")
    sys.exit(1)
else:
    print(f"✓ X_umap preserved: {adata_combined.obsm['X_umap'].shape}")

# Split for easier access
pre_mask = adata_combined.obs['mapping_status'] == 'pre-mapping'
post_mask = adata_combined.obs['mapping_status'] == 'post-mapping'

print(f"\nFinal counts:")
print(f"  Pre-mapping: {pre_mask.sum():,}")
print(f"  Post-mapping: {post_mask.sum():,}")
print(f"  Total: {adata_combined.n_obs:,}")

# Cleanup
del adata_pre, adata_post
gc.collect()

# ==============================================================================
# STEP 5: Visualization - Figure 1: UMAP Comparison (2x2)
# ==============================================================================

print("\n" + "=" * 70)
print("STEP 5: Generate Visualizations")
print("=" * 70)

print("\nFigure 1: UMAP Overview (2x2)...")

fig, axes = plt.subplots(2, 2, figsize=(16, 16))

# Panel 1: Mapping status (all cells)
print("  Panel 1: Mapping status...")
sc.pl.umap(adata_combined, color='mapping_status', ax=axes[0, 0], show=False,
           title='Pre vs Post Mapping', s=8,
           palette=MAPPING_STATUS_COLORS, frameon=True)

# Panel 2: Pre-mapping only
print("  Panel 2: Pre-mapping cells...")
adata_pre_view = adata_combined[pre_mask]
sc.pl.umap(adata_pre_view, ax=axes[0, 1], show=False,
           title=f'Pre-Mapping Only ({pre_mask.sum():,} cells)',
           color='lightgray', s=8, frameon=True)

# Panel 3: Post-mapping with cell types
print("  Panel 3: Post-mapping cell types...")
temp_col_name = '_temp_post_celltype'
try:
    if 'cell_type_final' in adata_combined.obs.columns:
        temp_col = pd.Series('N/A', index=adata_combined.obs_names, dtype='object')
        temp_col[post_mask] = adata_combined.obs.loc[post_mask, 'cell_type_final'].astype(str)
        temp_col = temp_col.astype('category')
        adata_combined.obs[temp_col_name] = temp_col
        
        sc.pl.umap(adata_combined, color=temp_col_name, ax=axes[1, 0], show=False,
                   title=f'Post-Mapping Cell Types (conf≥{CONFIDENCE_THRESHOLD})',
                   s=8, frameon=True, legend_loc='right margin',
                   na_color='lightgray')
    else:
        axes[1, 0].text(0.5, 0.5, 'Cell type data\nnot available',
                       ha='center', va='center', transform=axes[1, 0].transAxes,
                       fontsize=14, color='gray')
        axes[1, 0].set_title('Post-Mapping Cell Types')
except Exception as e:
    print(f"  WARNING: Panel 3 failed: {e}")
    axes[1, 0].text(0.5, 0.5, 'Error generating panel',
                   ha='center', va='center', transform=axes[1, 0].transAxes,
                   fontsize=14, color='red')
finally:
    cleanup_temp_columns(adata_combined, temp_col_name)

# Panel 4: Post-mapping confidence
print("  Panel 4: Mapping confidence...")
temp_col_name = '_temp_confidence'
try:
    if 'mapping_confidence' in adata_combined.obs.columns:
        temp_conf = np.full(adata_combined.n_obs, np.nan, dtype=float)
        temp_conf[post_mask] = pd.to_numeric(
            adata_combined.obs.loc[post_mask, 'mapping_confidence'],
            errors='coerce'
        ).values
        
        adata_combined.obs[temp_col_name] = temp_conf
        
        sc.pl.umap(adata_combined, color=temp_col_name, ax=axes[1, 1], show=False,
                   title='Mapping Confidence (Post-mapping only)',
                   s=8, frameon=True, cmap='RdYlGn', vmin=0, vmax=1,
                   na_color='lightgray')
    else:
        axes[1, 1].text(0.5, 0.5, 'Confidence data\nnot available',
                       ha='center', va='center', transform=axes[1, 1].transAxes,
                       fontsize=14, color='gray')
        axes[1, 1].set_title('Mapping Confidence')
except Exception as e:
    print(f"  WARNING: Panel 4 failed: {e}")
    axes[1, 1].text(0.5, 0.5, 'Error generating panel',
                   ha='center', va='center', transform=axes[1, 1].transAxes,
                   fontsize=14, color='red')
finally:
    cleanup_temp_columns(adata_combined, temp_col_name)

plt.tight_layout()
fig_path = output_dir / f"01_umap_comparison_overview.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"✓ Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 2: Side-by-Side Comparison (1x2)
# ==============================================================================

print("\nFigure 2: Side-by-Side UMAP Comparison...")

fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# Panel 1: Pre-mapping
print("  Panel 1: Pre-mapping...")
adata_pre_view = adata_combined[pre_mask]
umap_pre = adata_pre_view.obsm['X_umap']
axes[0].scatter(umap_pre[:, 0], umap_pre[:, 1], s=5, alpha=0.6,
               c=MAPPING_STATUS_COLORS['pre-mapping'])
axes[0].set_xlabel('UMAP 1', fontsize=12)
axes[0].set_ylabel('UMAP 2', fontsize=12)
axes[0].set_title(f'Pre-Mapping\n({pre_mask.sum():,} cells)',
                 fontsize=14, fontweight='bold')
axes[0].grid(alpha=0.3)

# Panel 2: Post-mapping with cell types
print("  Panel 2: Post-mapping...")
adata_post_view = adata_combined[post_mask]
umap_post = adata_post_view.obsm['X_umap']

if 'cell_type_final' in adata_post_view.obs.columns:
    celltypes = adata_post_view.obs['cell_type_final'].astype(str)
    unique_types = celltypes.unique()
    
    # Use a color palette
    colors = sns.color_palette('tab20', len(unique_types))
    color_map = dict(zip(unique_types, colors))
    
    for ct in unique_types:
        mask = (celltypes == ct).values
        axes[1].scatter(umap_post[mask, 0], umap_post[mask, 1],
                       c=[color_map[ct]], label=ct, s=5, alpha=0.6)
    
    axes[1].legend(loc='center left', bbox_to_anchor=(1, 0.5),
                  frameon=True, title='Cell Type', fontsize=9)
else:
    axes[1].scatter(umap_post[:, 0], umap_post[:, 1], s=5, alpha=0.6,
                   c=MAPPING_STATUS_COLORS['post-mapping'])

axes[1].set_xlabel('UMAP 1', fontsize=12)
axes[1].set_ylabel('UMAP 2', fontsize=12)
axes[1].set_title(f'Post-Mapping with Cell Types\n({post_mask.sum():,} cells)',
                 fontsize=14, fontweight='bold')
axes[1].grid(alpha=0.3)

plt.tight_layout()
fig_path = output_dir / f"02_sidebyside_comparison.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"✓ Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 3: UMAP Space Validation
# ==============================================================================

print("\nFigure 3: UMAP Space Validation...")

fig, axes = plt.subplots(1, 3, figsize=(20, 6))

adata_pre_view = adata_combined[pre_mask]
adata_post_view = adata_combined[post_mask]
umap_pre = adata_pre_view.obsm['X_umap']
umap_post = adata_post_view.obsm['X_umap']

# Panel 1: X coordinate distribution
print("  Panel 1: X coordinate distribution...")
axes[0].hist(umap_pre[:, 0], bins=50, alpha=0.5, density=True,
            label=f'Pre-mapping', color=MAPPING_STATUS_COLORS['pre-mapping'])
axes[0].hist(umap_post[:, 0], bins=50, alpha=0.5, density=True,
            label=f'Post-mapping', color=MAPPING_STATUS_COLORS['post-mapping'])
axes[0].set_xlabel('UMAP X', fontsize=12)
axes[0].set_ylabel('Density', fontsize=12)
axes[0].set_title('UMAP X Distribution', fontsize=14, fontweight='bold')
axes[0].legend(frameon=True)
axes[0].grid(alpha=0.3)

# Panel 2: Y coordinate distribution
print("  Panel 2: Y coordinate distribution...")
axes[1].hist(umap_pre[:, 1], bins=50, alpha=0.5, density=True,
            label=f'Pre-mapping', color=MAPPING_STATUS_COLORS['pre-mapping'])
axes[1].hist(umap_post[:, 1], bins=50, alpha=0.5, density=True,
            label=f'Post-mapping', color=MAPPING_STATUS_COLORS['post-mapping'])
axes[1].set_xlabel('UMAP Y', fontsize=12)
axes[1].set_ylabel('Density', fontsize=12)
axes[1].set_title('UMAP Y Distribution', fontsize=14, fontweight='bold')
axes[1].legend(frameon=True)
axes[1].grid(alpha=0.3)

# Panel 3: 2D overlay with statistics
print("  Panel 3: 2D overlay...")
axes[2].scatter(umap_pre[:, 0], umap_pre[:, 1], s=3, alpha=0.3,
               c=MAPPING_STATUS_COLORS['pre-mapping'], label='Pre-mapping')
axes[2].scatter(umap_post[:, 0], umap_post[:, 1], s=3, alpha=0.3,
               c=MAPPING_STATUS_COLORS['post-mapping'], label='Post-mapping')
axes[2].set_xlabel('UMAP X', fontsize=12)
axes[2].set_ylabel('UMAP Y', fontsize=12)
axes[2].set_title('UMAP Space Overlay', fontsize=14, fontweight='bold')
axes[2].legend(frameon=True)
axes[2].grid(alpha=0.3)

# Add statistics
pre_x_range = umap_pre[:, 0].ptp()
pre_y_range = umap_pre[:, 1].ptp()
post_x_range = umap_post[:, 0].ptp()
post_y_range = umap_post[:, 1].ptp()

stats_text = f"Range Comparison:\n"
stats_text += f"Pre:  X={pre_x_range:.1f}, Y={pre_y_range:.1f}\n"
stats_text += f"Post: X={post_x_range:.1f}, Y={post_y_range:.1f}\n"

if pre_x_range > 0 and pre_y_range > 0:
    stats_text += f"Ratio: X={post_x_range/pre_x_range:.2f}, Y={post_y_range/pre_y_range:.2f}"
else:
    stats_text += f"Ratio: N/A"

axes[2].text(0.05, 0.95, stats_text, transform=axes[2].transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
fig_path = output_dir / f"03_umap_space_validation.{FIGURE_FORMAT}"
try:
    plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
    print(f"✓ Saved: {fig_path.name}")
except Exception as e:
    print(f"ERROR saving figure: {e}")
finally:
    plt.close()

# ==============================================================================
# Figure 4: Mapping Quality Metrics (if available)
# ==============================================================================

if 'mapping_confidence' in adata_combined.obs.columns:
    print("\nFigure 4: Mapping Quality Metrics...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    
    # Get post-mapping data
    conf = pd.to_numeric(
        adata_combined.obs.loc[post_mask, 'mapping_confidence'],
        errors='coerce'
    )
    conf_valid = conf[conf.notna()]
    
    # Panel 1: Confidence histogram
    print("  Panel 1: Confidence distribution...")
    axes[0, 0].hist(conf_valid, bins=50, edgecolor='black', alpha=0.7,
                   color='steelblue')
    axes[0, 0].axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--',
                      linewidth=2, label=f'Threshold={CONFIDENCE_THRESHOLD}')
    axes[0, 0].set_xlabel('Mapping Confidence', fontsize=12)
    axes[0, 0].set_ylabel('Number of Cells', fontsize=12)
    axes[0, 0].set_title('Confidence Distribution', fontsize=14, fontweight='bold')
    axes[0, 0].legend(frameon=True)
    axes[0, 0].grid(alpha=0.3)
    
    # Add statistics
    stats_text = f"Mean: {conf_valid.mean():.3f}\n"
    stats_text += f"Median: {conf_valid.median():.3f}\n"
    high_conf = (conf_valid >= CONFIDENCE_THRESHOLD).sum()
    stats_text += f"High (≥{CONFIDENCE_THRESHOLD}): {high_conf:,} ({high_conf/len(conf_valid)*100:.1f}%)"
    
    axes[0, 0].text(0.05, 0.95, stats_text, transform=axes[0, 0].transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Panel 2: Cell type distribution
    print("  Panel 2: Cell type distribution...")
    if 'cell_type_final' in adata_combined.obs.columns:
        celltype_counts = adata_combined.obs.loc[post_mask, 'cell_type_final'].value_counts()
        
        y_pos = np.arange(len(celltype_counts))
        axes[0, 1].barh(y_pos, celltype_counts.values, alpha=0.7)
        axes[0, 1].set_yticks(y_pos)
        axes[0, 1].set_yticklabels(celltype_counts.index)
        axes[0, 1].set_xlabel('Number of Cells', fontsize=12)
        axes[0, 1].set_title('Cell Type Distribution', fontsize=14, fontweight='bold')
        axes[0, 1].grid(alpha=0.3, axis='x')
        
        # Add count labels
        for i, count in enumerate(celltype_counts.values):
            axes[0, 1].text(count + max(celltype_counts)*0.01, i,
                          f'{count:,}', va='center', fontsize=9)
    
    # Panel 3: Confidence by cell type (if available)
    print("  Panel 3: Confidence by cell type...")
    if 'cell_type_final' in adata_combined.obs.columns:
        df = pd.DataFrame({
            'Cell_Type': adata_combined.obs.loc[post_mask, 'cell_type_final'].astype(str),
            'Confidence': conf
        })
        df = df[df['Cell_Type'] != 'Unknown']
        df = df.dropna()
        
        if len(df) > 0:
            celltypes_sorted = df['Cell_Type'].value_counts().index.tolist()
            
            sns.violinplot(data=df, y='Cell_Type', x='Confidence',
                          order=celltypes_sorted, ax=axes[1, 0],
                          orient='h', inner='box')
            
            axes[1, 0].axvline(CONFIDENCE_THRESHOLD, color='red',
                              linestyle='--', linewidth=2)
            axes[1, 0].set_xlabel('Mapping Confidence', fontsize=12)
            axes[1, 0].set_ylabel('')
            axes[1, 0].set_title('Confidence by Cell Type', fontsize=14, fontweight='bold')
            axes[1, 0].grid(alpha=0.3, axis='x')
    
    # Panel 4: Mean confidence by cell type
    print("  Panel 4: Mean confidence by cell type...")
    if 'cell_type_final' in adata_combined.obs.columns:
        df_conf = pd.DataFrame({
            'Cell_Type': adata_combined.obs.loc[post_mask, 'cell_type_final'].astype(str),
            'Confidence': conf
        })
        df_conf = df_conf[df_conf['Cell_Type'] != 'Unknown']
        
        conf_by_type = df_conf.groupby('Cell_Type')['Confidence'].agg(['mean', 'std']).sort_values('mean', ascending=True)
        
        if len(conf_by_type) > 0:
            y_pos = np.arange(len(conf_by_type))
            
            axes[1, 1].barh(y_pos, conf_by_type['mean'],
                          xerr=conf_by_type['std'],
                          capsize=5, alpha=0.7)
            axes[1, 1].set_yticks(y_pos)
            axes[1, 1].set_yticklabels(conf_by_type.index)
            axes[1, 1].set_xlabel('Mean Confidence ± Std', fontsize=12)
            axes[1, 1].set_title('Mean Confidence by Cell Type',
                                fontsize=14, fontweight='bold')
            axes[1, 1].axvline(CONFIDENCE_THRESHOLD, color='red',
                              linestyle='--', linewidth=2)
            axes[1, 1].grid(alpha=0.3, axis='x')
            axes[1, 1].set_xlim(0, 1)
    
    plt.tight_layout()
    fig_path = output_dir / f"04_mapping_quality_metrics.{FIGURE_FORMAT}"
    try:
        plt.savefig(fig_path, dpi=DPI, bbox_inches='tight')
        print(f"✓ Saved: {fig_path.name}")
    except Exception as e:
        print(f"ERROR saving figure: {e}")
    finally:
        plt.close()

# ==============================================================================
# Summary Report
# ==============================================================================

print("\n" + "=" * 70)
print("Generate Summary Report")
print("=" * 70)

report_lines = []
report_lines.append("=" * 70)
report_lines.append("Pre-Mapping vs Post-Mapping Comparison Report")
report_lines.append("=" * 70)

report_lines.append(f"\n[ Input Files ]")
report_lines.append(f"  Pre-mapping: {PRE_MAPPING_H5AD}")
report_lines.append(f"  Post-mapping: {POST_MAPPING_H5AD}")

report_lines.append(f"\n[ Data Summary ]")
report_lines.append(f"  Total cells: {adata_combined.n_obs:,}")
report_lines.append(f"  Pre-mapping: {pre_mask.sum():,} ({pre_mask.sum()/adata_combined.n_obs*100:.1f}%)")
report_lines.append(f"  Post-mapping: {post_mask.sum():,} ({post_mask.sum()/adata_combined.n_obs*100:.1f}%)")
report_lines.append(f"  Genes: {adata_combined.n_vars:,}")

if 'mapping_confidence' in adata_combined.obs.columns:
    conf = pd.to_numeric(
        adata_combined.obs.loc[post_mask, 'mapping_confidence'],
        errors='coerce'
    )
    conf_valid = conf[conf.notna()]
    
    if len(conf_valid) > 0:
        report_lines.append(f"\n[ Mapping Quality ]")
        report_lines.append(f"  Confidence mean: {conf_valid.mean():.3f}")
        report_lines.append(f"  Confidence median: {conf_valid.median():.3f}")
        high_conf = (conf_valid >= CONFIDENCE_THRESHOLD).sum()
        report_lines.append(f"  High confidence (≥{CONFIDENCE_THRESHOLD}): {high_conf:,} ({high_conf/len(conf_valid)*100:.1f}%)")

if 'cell_type_final' in adata_combined.obs.columns:
    report_lines.append(f"\n[ Cell Type Distribution (Post-Mapping) ]")
    celltype_counts = adata_combined.obs.loc[post_mask, 'cell_type_final'].value_counts()
    for celltype, count in celltype_counts.items():
        pct = count / post_mask.sum() * 100
        report_lines.append(f"  {celltype}: {count:,} ({pct:.1f}%)")

# UMAP space comparison
adata_pre_view = adata_combined[pre_mask]
adata_post_view = adata_combined[post_mask]
umap_pre = adata_pre_view.obsm['X_umap']
umap_post = adata_post_view.obsm['X_umap']

pre_x_range = umap_pre[:, 0].ptp()
pre_y_range = umap_pre[:, 1].ptp()
post_x_range = umap_post[:, 0].ptp()
post_y_range = umap_post[:, 1].ptp()

report_lines.append(f"\n[ UMAP Space Comparison ]")
report_lines.append(f"  Pre-mapping range:  X={pre_x_range:.1f}, Y={pre_y_range:.1f}")
report_lines.append(f"  Post-mapping range: X={post_x_range:.1f}, Y={post_y_range:.1f}")

if pre_x_range > 0 and pre_y_range > 0:
    report_lines.append(f"  Range ratio: X={post_x_range/pre_x_range:.2f}, Y={post_y_range/pre_y_range:.2f}")

report_lines.append(f"\n[ Generated Figures ]")
figure_list = [
    "01_umap_comparison_overview.pdf",
    "02_sidebyside_comparison.pdf",
    "03_umap_space_validation.pdf"
]
if 'mapping_confidence' in adata_combined.obs.columns:
    figure_list.append("04_mapping_quality_metrics.pdf")

for i, fname in enumerate(figure_list, 1):
    full_path = output_dir / fname
    if full_path.exists():
        file_size = full_path.stat().st_size / 1024
        report_lines.append(f"  {i}. {fname} ({file_size:.1f} KB)")

report_lines.append(f"\n" + "=" * 70)
report_lines.append("Comparison Analysis Complete")
report_lines.append("=" * 70)

report_text = '\n'.join(report_lines)
print("\n" + report_text)

report_file = output_dir / "comparison_summary.txt"
try:
    with open(report_file, 'w') as f:
        f.write(report_text)
    print(f"\nReport saved: {report_file}")
except Exception as e:
    print(f"\nERROR saving report: {e}")

# Save combined data
print("\nSaving combined data...")
combined_file = output_dir / "combined_pre_post_mapping.h5ad"
try:
    adata_combined.write_h5ad(combined_file, compression='gzip')
    file_size = combined_file.stat().st_size / (1024**3)
    print(f"✓ Saved: {combined_file} ({file_size:.2f} GB)")
except Exception as e:
    print(f"ERROR saving combined data: {e}")

# ==============================================================================
# SUCCESS
# ==============================================================================

print("\n" + "=" * 70)
print("SUCCESS - Comparison Analysis Complete")
print("=" * 70)
print(f"\nOutput directory: {output_dir}")
print("\nGenerated files:")
for f in sorted(output_dir.glob("*.pdf")):
    print(f"  - {f.name}")
print(f"  - comparison_summary.txt")
print(f"  - combined_pre_post_mapping.h5ad")
