#!/usr/bin/env python3
"""
T Cell scVI Results - starCAT Analysis Pipeline

This script performs starCAT analysis on scVI-integrated T cell data.
starCAT provides high-resolution annotation using 52 Gene Expression Programs (GEPs).

Author: Clinical-Bioinformatics Team
Date: 2024-12-03
Version: v1.0

Requirements:
- scvi-tools integrated h5ad file with raw counts preserved
- TCAT.V1 reference files in /home/h2048/data/source/reference/TCAT.V1/
"""

import sys
import os
from pathlib import Path
import warnings
import time
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Single-cell analysis
import scanpy as sc

# starCAT
from starcat import starCAT

warnings.filterwarnings('ignore')

# ==================== Configuration ====================
print("="*70)
print("T CELL scVI + starCAT ANALYSIS PIPELINE")
print("="*70)
print(f"\nAnalysis started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Input/Output paths
INPUT_H5AD = "/home/h2048/data/py/1203/T_scvi_integration/adata_tcell_scvi_full.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/1203/T_scvi_integration/starcat_analysis"

# starCAT configuration
STARCAT_REFERENCE = "TCAT.V1"
LOCAL_REFERENCE_DIR = "/home/h2048/data/source/reference/TCAT.V1"

# Clustering key (from scVI)
# Will auto-detect if not found, priority: leiden_3 > leiden_res3.0 > leiden
CLUSTER_KEY = "leiden_3"  # Will try alternatives if not found
BATCH_KEY = "dataset"

# Marker genes for validation
MARKER_GENES = {
    "Pan_T": ['CD3D', 'CD3E', 'CD3G'],
    "NK": ['NCAM1', 'NKG7', 'GNLY', 'KLRD1'],
    "CD4": ['CD4', 'IL7R'],
    "CD8": ['CD8A', 'CD8B'],
    "Naive": ['CCR7', 'SELL', 'LEF1', 'TCF7'],
    "Effector": ['GZMA', 'GZMB', 'GZMK', 'PRF1'],
    "Treg": ['FOXP3', 'IL2RA', 'CTLA4'],
    "Activation": ['IFNG', 'TNF', 'CD69'],
    "Exhaustion": ['PDCD1', 'LAG3', 'TIGIT', 'HAVCR2', 'TOX'],
    "Proliferation": ['MKI67', 'TOP2A', 'PCNA']
}

# Visualization
DPI = 300
FIGURE_FORMAT = "pdf"

# Display configuration
print(f"\nConfiguration:")
print(f"  Input: {INPUT_H5AD}")
print(f"  Output: {OUTPUT_DIR}")
print(f"  starCAT reference: {STARCAT_REFERENCE}")
print(f"  Reference dir: {LOCAL_REFERENCE_DIR}")

# ==================== Create output directories ====================
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)

fig_dir = output_dir / "figures"
fig_dir.mkdir(exist_ok=True)

print(f"\n✓ Output directories created")

# ==================== Step 1: Load scVI results ====================
print("\n" + "="*70)
print("Step 1: Loading scVI-integrated data")
print("="*70)

if not Path(INPUT_H5AD).exists():
    raise FileNotFoundError(f"Input file not found: {INPUT_H5AD}")

adata = sc.read_h5ad(INPUT_H5AD)
print(f"\n✓ Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Display data structure
print(f"\nData structure:")
print(f"  obs columns: {len(adata.obs.columns)}")
print(f"  var columns: {len(adata.var.columns)}")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  layers: {list(adata.layers.keys()) if adata.layers else 'None'}")

# Check for UMAP
if 'X_umap' in adata.obsm:
    print(f"  ✓ UMAP coordinates found")
else:
    print(f"  ⚠️  No UMAP coordinates found")

# Check for clustering
if CLUSTER_KEY in adata.obs.columns:
    n_clusters = adata.obs[CLUSTER_KEY].nunique()
    print(f"  ✓ Clustering found: {CLUSTER_KEY} ({n_clusters} clusters)")
else:
    print(f"  ⚠️  Cluster key '{CLUSTER_KEY}' not found")
    
    # Smart fallback: try common patterns
    leiden_cols = [k for k in adata.obs.columns if 'leiden' in k.lower()]
    
    if leiden_cols:
        print(f"  Available clustering keys: {leiden_cols}")
        
        # Priority order
        priority_patterns = [
            'leiden_3', 'leiden_res3', 'leiden_res3.0',
            'leiden_2', 'leiden_res2', 'leiden_res2.0',
            'leiden_4', 'leiden_res4', 'leiden_res4.0',
            'leiden_5', 'leiden_res5', 'leiden_res5.0',
            'leiden'  # fallback to default
        ]
        
        # Find first match
        for pattern in priority_patterns:
            matches = [k for k in leiden_cols if pattern in k.lower()]
            if matches:
                CLUSTER_KEY = matches[0]
                n_clusters = adata.obs[CLUSTER_KEY].nunique()
                print(f"  → Auto-selected: {CLUSTER_KEY} ({n_clusters} clusters)")
                break
        
        if CLUSTER_KEY not in adata.obs.columns:
            # Just use first available
            CLUSTER_KEY = leiden_cols[0]
            n_clusters = adata.obs[CLUSTER_KEY].nunique()
            print(f"  → Using: {CLUSTER_KEY} ({n_clusters} clusters)")
    else:
        print(f"  ✗ No leiden clustering found in data")

# Check for batch info
if BATCH_KEY in adata.obs.columns:
    n_batches = adata.obs[BATCH_KEY].nunique()
    print(f"  ✓ Batch info found: {BATCH_KEY} ({n_batches} batches)")

# ==================== Step 2: Verify raw counts availability ====================
print("\n" + "="*70)
print("Step 2: Verifying raw counts for starCAT")
print("="*70)

print("\nstarCAT requires RAW COUNTS (not log-transformed)")

# Check data transformation status
raw_available = False
raw_source = None
adata_for_starcat = None

# Option 1: Check .raw
if adata.raw is not None:
    print(f"\n✓ .raw exists: {adata.raw.shape[0]:,} cells × {adata.raw.shape[1]:,} genes")
    
    # Verify it contains raw counts
    x_sample = adata.raw.X[:10000, :10000]
    if issparse(x_sample):
        x_sample = x_sample.toarray()
    x_max = np.max(x_sample)
    x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
    
    print(f"  .raw statistics (sample):")
    print(f"    Max: {x_max:.2f}, Mean (non-zero): {x_mean:.2f}")
    
    if x_max > 100 or (x_max > 20 and x_mean > 1):
        print(f"  ✓ .raw appears to contain raw counts")
        adata_for_starcat = sc.AnnData(
            X=adata.raw.X,
            obs=adata.obs.copy(),
            var=adata.raw.var.copy()
        )
        raw_source = ".raw"
        raw_available = True
    else:
        print(f"  ⚠️  .raw may be log-transformed (low max/mean)")

# Option 2: Check .layers['raw_counts'] (priority)
if not raw_available and 'raw_counts' in adata.layers:
    print(f"\n✓ .layers['raw_counts'] found")
    
    x_sample = adata.layers['raw_counts'][:10000, :10000]
    if issparse(x_sample):
        x_sample = x_sample.toarray()
    x_max = np.max(x_sample)
    x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
    
    print(f"  .layers['raw_counts'] statistics (sample):")
    print(f"    Max: {x_max:.2f}, Mean (non-zero): {x_mean:.2f}")
    
    if x_max > 100 or (x_max > 20 and x_mean > 1):
        print(f"  ✓ .layers['raw_counts'] appears to contain raw counts")
        adata_for_starcat = sc.AnnData(
            X=adata.layers['raw_counts'],
            obs=adata.obs.copy(),
            var=adata.var.copy()
        )
        raw_source = ".layers['raw_counts']"
        raw_available = True
    else:
        print(f"  ⚠️  .layers['raw_counts'] may be log-transformed (low max/mean)")

# Option 3: Check .layers['counts'] (fallback)
if not raw_available and 'counts' in adata.layers:
    print(f"\n✓ .layers['counts'] found")
    
    x_sample = adata.layers['counts'][:10000, :10000]
    if issparse(x_sample):
        x_sample = x_sample.toarray()
    x_max = np.max(x_sample)
    x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
    
    print(f"  .layers['counts'] statistics (sample):")
    print(f"    Max: {x_max:.2f}, Mean (non-zero): {x_mean:.2f}")
    
    if x_max > 100 or (x_max > 20 and x_mean > 1):
        print(f"  ✓ .layers['counts'] appears to contain raw counts")
        adata_for_starcat = sc.AnnData(
            X=adata.layers['counts'],
            obs=adata.obs.copy(),
            var=adata.var.copy()
        )
        raw_source = ".layers['counts']"
        raw_available = True
    else:
        print(f"  ⚠️  .layers['counts'] may be log-transformed (low max/mean)")

# Option 3: Check .X
if not raw_available:
    print(f"\n⚠️  No .raw or .layers['counts'] found, checking .X...")
    
    x_sample = adata.X[:10000, :10000]
    if issparse(x_sample):
        x_sample = x_sample.toarray()
    x_max = np.max(x_sample)
    x_mean = np.mean(x_sample[x_sample > 0]) if np.any(x_sample > 0) else 0
    
    print(f"  .X statistics (sample):")
    print(f"    Max: {x_max:.2f}, Mean (non-zero): {x_mean:.2f}")
    
    if x_max > 100 or (x_max > 20 and x_mean > 1):
        print(f"  ✓ .X appears to contain raw counts")
        adata_for_starcat = adata.copy()
        raw_source = ".X"
        raw_available = True
    else:
        print(f"  ✗ .X appears to be log-transformed")

# Final check
if not raw_available:
    print("\n" + "="*70)
    print("✗ ERROR: No raw counts found!")
    print("="*70)
    print("\nstarCAT requires raw counts (not log-transformed data)")
    print("\nPossible solutions:")
    print("  1. Use the original h5ad file before scVI processing")
    print("  2. Re-run scVI pipeline with .raw preserved")
    print("  3. Load the BBKNN h5ad which has raw counts")
    print("\nExiting...")
    sys.exit(1)

print(f"\n✓ Will use raw counts from: {raw_source}")
print(f"  Shape: {adata_for_starcat.shape}")

# ==================== Step 3: Load starCAT reference ====================
print("\n" + "="*70)
print("Step 3: Loading starCAT reference")
print("="*70)

# Verify reference directory
ref_dir = Path(LOCAL_REFERENCE_DIR)
if not ref_dir.exists():
    raise FileNotFoundError(f"Reference directory not found: {LOCAL_REFERENCE_DIR}")

print(f"\nReference directory: {ref_dir}")

# Check required files
ref_file = ref_dir / f"{STARCAT_REFERENCE}.reference.tsv"
score_file = ref_dir / f"{STARCAT_REFERENCE}.scores.yaml"

print(f"\nVerifying reference files:")
if ref_file.exists():
    size_mb = ref_file.stat().st_size / (1024**2)
    print(f"  ✓ Reference TSV: {ref_file.name} ({size_mb:.2f} MB)")
else:
    raise FileNotFoundError(f"Reference file not found: {ref_file}")

if score_file.exists():
    size_kb = score_file.stat().st_size / 1024
    print(f"  ✓ Scores YAML: {score_file.name} ({size_kb:.2f} KB)")
else:
    raise FileNotFoundError(f"Scores file not found: {score_file}")

# Initialize starCAT
print(f"\nInitializing starCAT...")
tcat = starCAT(reference=str(ref_file), score_path=str(score_file))

print(f"\n✓ starCAT reference loaded: {tcat.ref_name}")
print(f"  Reference dimensions: {tcat.ref.shape[0]} genes × {tcat.ref.shape[1]} GEPs")

# Display first few GEPs and genes
print(f"\nFirst 5 genes and 10 GEPs:")
print(tcat.ref.iloc[:5, :10])

# Display score configuration
if hasattr(tcat, 'score_data') and tcat.score_data:
    print(f"\nstarCAT score configuration:")
    if 'scores' in tcat.score_data:
        if 'continuous' in tcat.score_data['scores']:
            cont_scores = [s['name'] for s in tcat.score_data['scores']['continuous']]
            print(f"  Continuous scores: {cont_scores}")
        if 'discrete' in tcat.score_data['scores']:
            disc_scores = [s['name'] for s in tcat.score_data['scores']['discrete']]
            print(f"  Discrete scores: {disc_scores}")

# ==================== Step 4: Run starCAT analysis ====================
print("\n" + "="*70)
print("Step 4: Running starCAT fit_transform")
print("="*70)

print(f"\nDataset for starCAT:")
print(f"  Cells: {adata_for_starcat.n_obs:,}")
print(f"  Genes: {adata_for_starcat.n_vars:,}")
print(f"  Source: {raw_source}")

print(f"\nRunning starCAT (this may take several minutes)...")
start_time = time.time()

try:
    usage, scores = tcat.fit_transform(adata_for_starcat)
    
    elapsed = time.time() - start_time
    print(f"\n✓ starCAT completed in {elapsed:.1f}s ({elapsed/60:.1f} min)")
    
    print(f"\nResults:")
    print(f"  Usage matrix: {usage.shape} (cells × GEPs)")
    print(f"  Scores dataframe: {scores.shape} (cells × scores)")
    
    # Display samples
    print(f"\nUsage matrix (first 5 cells, first 10 GEPs):")
    print(usage.iloc[:5, :10])
    
    print(f"\nScores dataframe (first 5 cells):")
    print(scores.head())
    
except Exception as e:
    print(f"\n✗ Error during starCAT analysis: {e}")
    print("\nTroubleshooting:")
    print("  1. Verify data contains raw counts")
    print("  2. Check gene names match reference")
    print("  3. Ensure sufficient memory")
    raise

# ==================== Step 5: Integrate results into adata ====================
print("\n" + "="*70)
print("Step 5: Integrating starCAT results")
print("="*70)

# Add GEP usage to adata
print(f"\nAdding GEP usage to adata.obs...")
for col in usage.columns:
    adata.obs[f"GEP_{col}"] = usage[col].values
print(f"✓ Added {len(usage.columns)} GEP usage columns")

# Add scores to adata
print(f"\nAdding scores to adata.obs...")
for col in scores.columns:
    # Convert binary columns to categorical for better plotting
    if col.endswith('_binary'):
        adata.obs[col] = scores[col].astype(str).values
    else:
        adata.obs[col] = scores[col].values
print(f"✓ Added {len(scores.columns)} score columns")

print(f"\n✓ Integration complete")
print(f"  Total new columns: {len(usage.columns) + len(scores.columns)}")

# ==================== Step 6: Basic visualization ====================
print("\n" + "="*70)
print("Step 6: Basic UMAP visualization")
print("="*70)

# Set scanpy parameters
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

# Plot batch and clustering (for reference)
print(f"\nGenerating reference UMAP plots...")

if BATCH_KEY in adata.obs.columns:
    fig = sc.pl.umap(adata, color=BATCH_KEY, return_fig=True, show=False)
    fig.savefig(fig_dir / f"umap_by_batch.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"  ✓ Saved: umap_by_batch.{FIGURE_FORMAT}")

if CLUSTER_KEY in adata.obs.columns:
    fig = sc.pl.umap(adata, color=CLUSTER_KEY, legend_loc='on data',
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"umap_by_cluster.{FIGURE_FORMAT}", dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"  ✓ Saved: umap_by_cluster.{FIGURE_FORMAT}")

# ==================== Step 7: starCAT discrete features ====================
print("\n" + "="*70)
print("Step 7: Visualizing starCAT discrete features")
print("="*70)

discrete_features = [col for col in scores.columns if col.endswith('_binary') or col == 'Multinomial_Label']

if len(discrete_features) > 0:
    print(f"\nDiscrete features: {discrete_features}")
    
    fig = sc.pl.umap(adata, color=discrete_features, ncols=2,
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"starcat_discrete_features.{FIGURE_FORMAT}",
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"✓ Saved: starcat_discrete_features.{FIGURE_FORMAT}")
else:
    print(f"⚠️  No discrete features found")

# ==================== Step 8: starCAT continuous scores ====================
print("\n" + "="*70)
print("Step 8: Visualizing starCAT continuous scores")
print("="*70)

continuous_features = [col for col in scores.columns 
                      if not col.endswith('_binary') and col != 'Multinomial_Label']

if len(continuous_features) > 0:
    print(f"\nContinuous features: {continuous_features}")
    
    fig = sc.pl.umap(adata, color=continuous_features, ncols=2, vmax='p99',
                     return_fig=True, show=False)
    fig.savefig(fig_dir / f"starcat_continuous_scores.{FIGURE_FORMAT}",
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"✓ Saved: starcat_continuous_scores.{FIGURE_FORMAT}")
else:
    print(f"⚠️  No continuous features found")

# ==================== Step 9: Important GEPs ====================
print("\n" + "="*70)
print("Step 9: Visualizing important GEPs")
print("="*70)

# Select important GEPs based on keywords
important_geps = []
gep_keywords = ['Cytotoxic', 'Exhaustion', 'Activation', 'Tfh', 'Treg', 'Naive',
                'Memory', 'Proliferation', 'ISG', 'IFN', 'Effector']

print(f"\nSearching for GEPs with keywords: {gep_keywords}")

for gep in usage.columns:
    for keyword in gep_keywords:
        if keyword.lower() in gep.lower():
            important_geps.append(f"GEP_{gep}")
            break

if len(important_geps) > 0:
    print(f"✓ Found {len(important_geps)} important GEPs")
    print(f"  {[g.replace('GEP_', '') for g in important_geps[:10]]}...")
    
    # Plot in batches
    batch_size = 12
    for i in range(0, len(important_geps), batch_size):
        batch = important_geps[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        print(f"\n  Batch {batch_num}: {len(batch)} GEPs")
        
        fig = sc.pl.umap(adata, color=batch, ncols=3, vmin=0, vmax='p99',
                         return_fig=True, show=False)
        fig.savefig(fig_dir / f"starcat_important_geps_batch{batch_num}.{FIGURE_FORMAT}",
                    dpi=DPI, bbox_inches='tight')
        plt.show()
        plt.close()
    
    print(f"\n✓ Saved {(len(important_geps) - 1) // batch_size + 1} GEP batch plots")
else:
    print(f"⚠️  No important GEPs found with keywords")

# ==================== Step 10: Cluster × starCAT comparison ====================
print("\n" + "="*70)
print("Step 10: Cluster × starCAT label comparison")
print("="*70)

if CLUSTER_KEY in adata.obs.columns and 'Multinomial_Label' in adata.obs.columns:
    print(f"\nComparing {CLUSTER_KEY} vs starCAT Multinomial_Label...")
    
    # Create contingency table
    contingency = pd.crosstab(
        adata.obs[CLUSTER_KEY],
        adata.obs['Multinomial_Label']
    )
    
    print(f"\nContingency table:")
    print(contingency)
    
    # Save table
    contingency.to_csv(output_dir / "cluster_vs_starcat_label.csv")
    
    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Heatmap of counts
    sns.heatmap(contingency, annot=True, fmt='d', cmap='YlOrRd', ax=axes[0])
    axes[0].set_title(f'Cell Counts: {CLUSTER_KEY} vs starCAT Label')
    axes[0].set_xlabel('starCAT Multinomial Label')
    axes[0].set_ylabel(CLUSTER_KEY)
    
    # Heatmap of proportions
    contingency_norm = contingency.div(contingency.sum(axis=1), axis=0)
    sns.heatmap(contingency_norm, annot=True, fmt='.2f', cmap='YlOrRd', ax=axes[1])
    axes[1].set_title(f'Proportions: {CLUSTER_KEY} vs starCAT Label')
    axes[1].set_xlabel('starCAT Multinomial Label')
    axes[1].set_ylabel(CLUSTER_KEY)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f"cluster_starcat_comparison.{FIGURE_FORMAT}",
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    
    print(f"✓ Saved: cluster_starcat_comparison.{FIGURE_FORMAT}")
else:
    print(f"⚠️  Cannot create comparison (missing {CLUSTER_KEY} or Multinomial_Label)")

# ==================== Step 11: GEP usage by cluster ====================
print("\n" + "="*70)
print("Step 11: GEP usage by cluster heatmap")
print("="*70)

if CLUSTER_KEY in adata.obs.columns:
    print(f"\nCalculating mean GEP usage per cluster...")
    
    gep_cols = [f"GEP_{col}" for col in usage.columns]
    cluster_gep_mean = adata.obs.groupby(CLUSTER_KEY)[gep_cols].mean()
    
    # Transpose for better visualization
    cluster_gep_mean_T = cluster_gep_mean.T
    cluster_gep_mean_T.index = [idx.replace('GEP_', '') for idx in cluster_gep_mean_T.index]
    
    # Plot heatmap
    fig, ax = plt.subplots(figsize=(max(12, len(cluster_gep_mean_T.columns)*0.6),
                                    max(15, len(cluster_gep_mean_T.index)*0.25)))
    
    sns.heatmap(cluster_gep_mean_T, cmap='RdYlBu_r', center=0,
                cbar_kws={'label': 'Mean GEP Usage'}, ax=ax)
    ax.set_title(f'Mean GEP Usage by {CLUSTER_KEY}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Cluster', fontsize=12)
    ax.set_ylabel('Gene Expression Program (GEP)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(fig_dir / f"gep_usage_by_cluster_heatmap.{FIGURE_FORMAT}",
                dpi=DPI, bbox_inches='tight')
    plt.show()
    plt.close()
    
    # Save table
    cluster_gep_mean_T.to_csv(output_dir / "gep_usage_by_cluster.csv")
    
    print(f"✓ Saved: gep_usage_by_cluster_heatmap.{FIGURE_FORMAT}")
    print(f"✓ Saved: gep_usage_by_cluster.csv")

# ==================== Step 12: GEP correlation matrix ====================
print("\n" + "="*70)
print("Step 12: GEP correlation matrix")
print("="*70)

print(f"\nComputing GEP correlation matrix...")
gep_corr = usage.corr()

# Plot correlation heatmap
fig, ax = plt.subplots(figsize=(14, 12))
sns.heatmap(gep_corr, cmap='coolwarm', center=0,
            cbar_kws={'label': 'Correlation'}, ax=ax)
ax.set_title('GEP Usage Correlation Matrix', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(fig_dir / f"gep_correlation_matrix.{FIGURE_FORMAT}",
            dpi=DPI, bbox_inches='tight')
plt.show()
plt.close()

print(f"✓ Saved: gep_correlation_matrix.{FIGURE_FORMAT}")

# Find highly correlated GEPs
print(f"\nHighly correlated GEP pairs (|r| > 0.7):")
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
    print(high_corr_df)
    high_corr_df.to_csv(output_dir / "gep_high_correlations.csv", index=False)
    print(f"✓ Saved: gep_high_correlations.csv")
else:
    print(f"  No GEP pairs with |r| > 0.7 found")

# ==================== Step 13: Save annotated data ====================
print("\n" + "="*70)
print("Step 13: Saving annotated data")
print("="*70)

# Save full annotated dataset
output_file = output_dir / "adata_tcell_scvi_starcat_annotated.h5ad"
print(f"\nSaving annotated dataset to: {output_file}")
adata.write_h5ad(output_file)
print(f"✓ Saved: {output_file}")

# Save starCAT results separately
usage.to_csv(output_dir / "starcat_usage_matrix.csv")
scores.to_csv(output_dir / "starcat_scores.csv")
print(f"✓ Saved: starcat_usage_matrix.csv")
print(f"✓ Saved: starcat_scores.csv")

# ==================== Step 14: Analysis summary ====================
print("\n" + "="*70)
print("Step 14: Analysis Summary")
print("="*70)

summary = []
summary.append("="*70)
summary.append("T CELL scVI + starCAT ANALYSIS - SUMMARY")
summary.append("="*70)
summary.append(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
summary.append("")

summary.append("[ Dataset Information ]")
summary.append(f"  Input: {INPUT_H5AD}")
summary.append(f"  Cells: {adata.n_obs:,}")
summary.append(f"  Genes: {adata.n_vars:,}")
summary.append(f"  Raw counts source: {raw_source}")
summary.append("")

summary.append("[ Batch Information ]")
if BATCH_KEY in adata.obs.columns:
    batch_counts = adata.obs[BATCH_KEY].value_counts().sort_index()
    for batch, count in batch_counts.items():
        pct = count / adata.n_obs * 100
        summary.append(f"  {batch}: {count:,} cells ({pct:.1f}%)")
summary.append("")

summary.append("[ Clustering ]")
if CLUSTER_KEY in adata.obs.columns:
    n_clusters = adata.obs[CLUSTER_KEY].nunique()
    summary.append(f"  Key: {CLUSTER_KEY}")
    summary.append(f"  Clusters: {n_clusters}")
summary.append("")

summary.append("[ starCAT Results ]")
summary.append(f"  Reference: {STARCAT_REFERENCE}")
summary.append(f"  GEPs: {usage.shape[1]}")
summary.append(f"  Scores: {scores.shape[1]}")

if 'Multinomial_Label' in adata.obs.columns:
    summary.append("\n  Cell type distribution (starCAT):")
    label_counts = adata.obs['Multinomial_Label'].value_counts()
    for label, count in label_counts.items():
        pct = count / adata.n_obs * 100
        summary.append(f"    {label}: {count:,} cells ({pct:.1f}%)")
summary.append("")

summary.append("[ Output Files ]")
summary.append(f"  Annotated h5ad: {output_dir / 'adata_tcell_scvi_starcat_annotated.h5ad'}")
summary.append(f"  Usage matrix: {output_dir / 'starcat_usage_matrix.csv'}")
summary.append(f"  Scores: {output_dir / 'starcat_scores.csv'}")
summary.append(f"  GEP by cluster: {output_dir / 'gep_usage_by_cluster.csv'}")
summary.append(f"  Figures: {fig_dir / f'*.{FIGURE_FORMAT}'}")
summary.append("")

summary.append("="*70)
summary.append("ANALYSIS COMPLETED SUCCESSFULLY")
summary.append("="*70)

summary_text = '\n'.join(summary)
print(summary_text)

# Save summary
summary_file = output_dir / "analysis_summary.txt"
with open(summary_file, 'w') as f:
    f.write(summary_text)
print(f"\n✓ Summary saved to: {summary_file}")

print("\n" + "="*70)
print("✅ ALL DONE!")
print("="*70)
print(f"\nNext steps:")
print(f"  1. Review figures in: {fig_dir}")
print(f"  2. Load annotated data: {output_file}")
print(f"  3. Explore GEP patterns in specific clusters")
print(f"  4. Compare starCAT labels with marker genes")