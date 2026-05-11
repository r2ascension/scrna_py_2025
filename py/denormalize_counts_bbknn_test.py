#!/usr/bin/env python3
"""
Denormalize Counts and Re-run BBKNN Integration Test (PER-SAMPLE)
==================================================================

Purpose:
--------
1. Detect if layers['counts'] contains normalized (log-transformed) data **PER SAMPLE**
2. Apply denormalization **PER SAMPLE**: raw_counts = exp(log_normalized) - 1
3. Optionally apply size factor correction per sample
4. Re-run BBKNN integration with corrected raw counts
5. Compare results with original integration

Critical Context:
-----------------
- scVI/scANVI REQUIRE raw integer UMI counts in layers['counts']
- BBKNN preprocessing may have stored log-normalized counts in SOME or ALL samples
- **Different samples/batches may have different preprocessing states**
- Detection heuristics (per sample):
  * Raw counts: integers, high values (>100), wide range
  * Normalized counts: floats, low values (<10), narrow range

Key Improvement:
----------------
**PER-SAMPLE PROCESSING**: Each sample is independently checked and denormalized
if needed, rather than treating all samples uniformly.

Author: Claude Code
Date: 2025-12-13
"""

import sys
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

# Configuration
INPUT_FILE = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{datetime.now().strftime('%Y%m%d')}/denormalize_counts_test")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# BBKNN parameters
BBKNN_NEIGHBORS_WITHIN_BATCH = 5
N_PCS = 50

# Logging setup
LOG_FILE = OUTPUT_DIR / f"denormalization_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

def log_print(msg):
    """Print and log simultaneously"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    formatted_msg = f"[{timestamp}] {msg}"
    print(formatted_msg)
    with open(LOG_FILE, 'a') as f:
        f.write(formatted_msg + '\n')

def detect_data_type(data, sample_size=10000):
    """
    Detect if data is raw counts or log-normalized

    Heuristics:
    -----------
    Raw counts:
        - Mostly integers (>95% match floor)
        - High max values (typically >100)
        - Wide range
        - Mean typically 1-50

    Log-normalized (log(counts+1)):
        - Float values
        - Low max values (typically <10)
        - Narrow range
        - Mean typically 0.5-3
    """
    if sparse.issparse(data):
        # Sample random non-zero values
        data_dense = data.data
        if len(data_dense) > sample_size:
            indices = np.random.choice(len(data_dense), sample_size, replace=False)
            sample = data_dense[indices]
        else:
            sample = data_dense
    else:
        flat = data.flatten()
        if len(flat) > sample_size:
            sample = np.random.choice(flat, sample_size, replace=False)
        else:
            sample = flat

    # Remove zeros for analysis
    sample_nonzero = sample[sample > 0]

    if len(sample_nonzero) == 0:
        zero_stats = {
            'mean': 0.0,
            'median': 0.0,
            'max': 0.0,
            'min': 0.0,
            'std': 0.0,
            'pct_integer': 100.0,
            'range': 0.0
        }
        return "ALL_ZERO", zero_stats

    # Calculate statistics
    stats = {
        'mean': float(np.mean(sample_nonzero)),
        'median': float(np.median(sample_nonzero)),
        'max': float(np.max(sample_nonzero)),
        'min': float(np.min(sample_nonzero)),
        'std': float(np.std(sample_nonzero)),
        'pct_integer': float(np.mean(np.abs(sample_nonzero - np.round(sample_nonzero)) < 1e-6) * 100),
        'range': float(np.max(sample_nonzero) - np.min(sample_nonzero))
    }

    # Decision logic
    is_integer = stats['pct_integer'] > 95
    is_high_max = stats['max'] > 50
    is_low_max = stats['max'] < 40
    is_low_mean = stats['mean'] < 20

    if is_integer and is_high_max:
        data_type = "RAW_COUNTS"
    elif not is_integer and is_low_max and is_low_mean:
        data_type = "LOG_NORMALIZED"
    else:
        data_type = "LOG_NORMALIZED"

    return data_type, stats

def denormalize_log_counts(log_data):
    """
    Reverse log-normalization: raw_counts = exp(log_data) - 1
    Assumes original transformation was: log(counts + 1)
    """
    if sparse.issparse(log_data):
        # For sparse matrix, apply element-wise
        raw_counts = log_data.copy()
        raw_counts.data = np.expm1(raw_counts.data)  # expm1(x) = exp(x) - 1
    else:
        raw_counts = np.expm1(log_data)

    # Round to integers (raw counts should be integers)
    if sparse.issparse(raw_counts):
        raw_counts.data = np.round(raw_counts.data)
        raw_counts.eliminate_zeros()
    else:
        raw_counts = np.round(raw_counts)

    return raw_counts

def apply_size_factor_correction(counts, size_factors):
    """
    Apply size factor correction if needed
    Formula: raw_counts_original = raw_counts_normalized * size_factors
    """
    size_factors = np.asarray(size_factors).astype(float)

    if len(size_factors) != counts.shape[0]:
        raise ValueError("Size factor length must match number of rows in counts matrix")

    if sparse.issparse(counts):
        # Scale rows in-place without constructing large dense intermediates
        counts_csr = counts.tocsr(copy=True)
        repeats = np.diff(counts_csr.indptr)
        if repeats.size != counts_csr.shape[0]:
            # indptr structure guarantees len(repeats)==n_rows, but guard anyway
            repeats = np.resize(repeats, counts_csr.shape[0])
        row_scale = np.repeat(size_factors, repeats)
        counts_csr.data = counts_csr.data * row_scale
        corrected = counts_csr
    else:
        corrected = counts * size_factors[:, np.newaxis]

    # Round to integers
    if sparse.issparse(corrected):
        corrected.data = np.round(corrected.data)
        corrected.eliminate_zeros()
    else:
        corrected = np.round(corrected)

    return corrected

def process_sample_denormalization(adata, batch_key='dataset'):
    """
    **PER-SAMPLE DENORMALIZATION**

    For each unique sample/batch:
    1. Extract that sample's data from layers['counts']
    2. Detect if it's normalized or raw
    3. Denormalize if needed
    4. Update the counts layer for that sample

    Returns:
    --------
    per_sample_stats : DataFrame with detection results per sample
    any_denormalized : bool indicating if any sample was denormalized
    """
    log_print("\n" + "="*80)
    log_print("PER-SAMPLE DENORMALIZATION ANALYSIS")
    log_print("="*80)

    if 'counts' not in adata.layers:
        log_print("ERROR: 'counts' layer not found!")
        return None, False

    # Get unique samples
    samples = adata.obs[batch_key].unique()
    n_samples = len(samples)
    log_print(f"\nFound {n_samples} unique samples in '{batch_key}':")
    for i, sample in enumerate(samples[:10]):  # Show first 10
        n_cells = (adata.obs[batch_key] == sample).sum()
        log_print(f"  {i+1}. {sample}: {n_cells} cells")
    if n_samples > 10:
        log_print(f"  ... and {n_samples - 10} more samples")

    # Create backup of original counts
    log_print("\nCreating backup of original counts...")
    adata.layers['counts_original'] = adata.layers['counts'].copy()

    # Process each sample
    per_sample_results = []
    samples_denormalized = []
    samples_raw = []
    samples_uncertain = []

    log_print("\n" + "-"*80)
    log_print("Analyzing each sample...")
    log_print("-"*80)

    for idx, sample in enumerate(samples):
        log_print(f"\n[{idx+1}/{n_samples}] Sample: {sample}")

        # Get cell indices for this sample
        sample_mask = adata.obs[batch_key] == sample
        sample_indices = np.where(sample_mask)[0]
        n_cells = len(sample_indices)

        log_print(f"  Cells: {n_cells}")

        # Extract sample's count data (use integer indices to avoid pandas Series issues)
        sample_counts = adata.layers['counts'][sample_indices, :]

        # Detect data type
        data_type, stats = detect_data_type(sample_counts)

        log_print(f"  Detection: {data_type}")
        log_print(f"    Max: {stats['max']:.2f}, Mean: {stats['mean']:.2f}, % Integer: {stats['pct_integer']:.1f}%")

        # Store results
        result = {
            'sample': sample,
            'n_cells': n_cells,
            'data_type': data_type,
            'denormalized': False,
            **{f'original_{k}': v for k, v in stats.items()}
        }

        # Denormalize if needed
        if data_type == "LOG_NORMALIZED":
            log_print(f"  → DENORMALIZING this sample...")
            samples_denormalized.append(sample)

            # Denormalize
            denorm_counts = denormalize_log_counts(sample_counts)

            # Check for size factors (sample-specific)
            size_factor_col = None
            for col in ['size_factors', 'sizeFactor', 'normalization_factor']:
                if col in adata.obs.columns:
                    size_factor_col = col
                    break

            if size_factor_col is not None:
                sample_size_factors = adata.obs.loc[sample_mask, size_factor_col].values
                sf_mean = np.mean(sample_size_factors)
                sf_std = np.std(sample_size_factors)

                log_print(f"    Applying size factors (mean={sf_mean:.4f}, std={sf_std:.4f})")

                if sf_mean > 0.1 and sf_mean < 10:
                    denorm_counts = apply_size_factor_correction(denorm_counts, sample_size_factors)

            # Update counts for this sample
            adata.layers['counts'][sample_indices, :] = denorm_counts

            # Verify denormalization
            denorm_type, denorm_stats = detect_data_type(denorm_counts)
            log_print(f"  → After denormalization: {denorm_type}")
            log_print(f"    Max: {denorm_stats['max']:.2f}, Mean: {denorm_stats['mean']:.2f}, % Integer: {denorm_stats['pct_integer']:.1f}%")

            result['denormalized'] = True
            result.update({f'denorm_{k}': v for k, v in denorm_stats.items()})

        elif data_type == "RAW_COUNTS":
            log_print(f"  → Already raw counts, no action needed")
            samples_raw.append(sample)
        else:
            log_print(f"  → UNCERTAIN, skipping denormalization")
            samples_uncertain.append(sample)

        per_sample_results.append(result)

    # Create summary DataFrame
    per_sample_df = pd.DataFrame(per_sample_results)

    # Summary
    log_print("\n" + "="*80)
    log_print("PER-SAMPLE ANALYSIS SUMMARY")
    log_print("="*80)
    log_print(f"\nTotal samples: {n_samples}")
    log_print(f"  - Already raw counts: {len(samples_raw)} samples")
    log_print(f"  - Denormalized: {len(samples_denormalized)} samples")
    log_print(f"  - Uncertain: {len(samples_uncertain)} samples")

    if len(samples_denormalized) > 0:
        log_print(f"\nDenormalized samples:")
        for s in samples_denormalized:
            log_print(f"  - {s}")

    if len(samples_uncertain) > 0:
        log_print(f"\nUncertain samples (review manually):")
        for s in samples_uncertain:
            log_print(f"  - {s}")

    # Save per-sample stats
    stats_file = OUTPUT_DIR / 'per_sample_denormalization_stats.csv'
    per_sample_df.to_csv(stats_file, index=False)
    log_print(f"\nPer-sample statistics saved to: {stats_file.name}")

    any_denormalized = len(samples_denormalized) > 0

    return per_sample_df, any_denormalized

def visualize_per_sample_stats(per_sample_df, output_dir):
    """Create per-sample comparison plots"""
    log_print("\nCreating per-sample diagnostic plots...")

    # Separate samples by type
    normalized_samples = per_sample_df[per_sample_df['denormalized'] == True]
    raw_samples = per_sample_df[per_sample_df['data_type'] == 'RAW_COUNTS']

    if len(normalized_samples) == 0:
        log_print("  No samples were denormalized, skipping comparison plots")
        return

    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Max values
    ax = axes[0, 0]
    x_pos = np.arange(len(normalized_samples))
    ax.bar(x_pos - 0.2, normalized_samples['original_max'], width=0.4, label='Original (Normalized)', alpha=0.7)
    ax.bar(x_pos + 0.2, normalized_samples['denorm_max'], width=0.4, label='After Denormalization', alpha=0.7)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Max Value')
    ax.set_title('Max Values: Before vs After Denormalization')
    ax.legend()
    ax.set_xticks(x_pos)
    ax.set_xticklabels(normalized_samples['sample'].values, rotation=45, ha='right', fontsize=8)

    # Plot 2: Mean values
    ax = axes[0, 1]
    ax.bar(x_pos - 0.2, normalized_samples['original_mean'], width=0.4, label='Original (Normalized)', alpha=0.7)
    ax.bar(x_pos + 0.2, normalized_samples['denorm_mean'], width=0.4, label='After Denormalization', alpha=0.7)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Mean Value')
    ax.set_title('Mean Values: Before vs After Denormalization')
    ax.legend()
    ax.set_xticks(x_pos)
    ax.set_xticklabels(normalized_samples['sample'].values, rotation=45, ha='right', fontsize=8)

    # Plot 3: % Integer
    ax = axes[1, 0]
    ax.bar(x_pos - 0.2, normalized_samples['original_pct_integer'], width=0.4, label='Original (Normalized)', alpha=0.7)
    ax.bar(x_pos + 0.2, normalized_samples['denorm_pct_integer'], width=0.4, label='After Denormalization', alpha=0.7)
    ax.set_xlabel('Sample')
    ax.set_ylabel('% Integer Values')
    ax.set_title('% Integer Values: Before vs After Denormalization')
    ax.axhline(y=95, color='r', linestyle='--', label='95% threshold')
    ax.legend()
    ax.set_xticks(x_pos)
    ax.set_xticklabels(normalized_samples['sample'].values, rotation=45, ha='right', fontsize=8)

    # Plot 4: Sample type distribution
    ax = axes[1, 1]
    type_counts = per_sample_df['data_type'].value_counts()
    colors = {'RAW_COUNTS': 'green', 'LOG_NORMALIZED': 'orange', 'UNCERTAIN': 'red'}
    ax.bar(range(len(type_counts)), type_counts.values,
           color=[colors.get(t, 'gray') for t in type_counts.index])
    ax.set_xticks(range(len(type_counts)))
    ax.set_xticklabels(type_counts.index, rotation=45, ha='right')
    ax.set_ylabel('Number of Samples')
    ax.set_title('Sample Data Type Distribution')
    for i, v in enumerate(type_counts.values):
        ax.text(i, v, str(v), ha='center', va='bottom')

    plt.suptitle(f'Per-Sample Denormalization Analysis ({len(per_sample_df)} samples)',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / 'per_sample_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    log_print(f"  Saved: per_sample_comparison.png")

def run_bbknn_integration(adata, batch_key='dataset', n_pcs=50, neighbors_within_batch=5):
    """
    Run BBKNN integration pipeline

    Steps:
    1. Normalize and log-transform for PCA
    2. Run PCA
    3. BBKNN batch correction
    4. UMAP visualization
    """
    import bbknn

    log_print("\n" + "="*60)
    log_print("Running BBKNN Integration")
    log_print("="*60)

    # Check if counts layer exists
    if 'counts' not in adata.layers:
        log_print("ERROR: No 'counts' layer found!")
        return None

    log_print(f"Input shape: {adata.shape}")
    log_print(f"Batch key: {batch_key}")
    log_print(f"Batches: {adata.obs[batch_key].nunique()}")

    # Create working copy
    adata_bbknn = adata.copy()

    # Store raw counts
    adata_bbknn.raw = adata_bbknn

    # Normalize and log-transform for PCA (from counts layer)
    log_print("\nNormalizing and log-transforming for PCA...")
    sc.pp.normalize_total(adata_bbknn, target_sum=1e4)
    sc.pp.log1p(adata_bbknn)

    # Find highly variable genes
    log_print("Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata_bbknn, n_top_genes=2000, batch_key=batch_key)
    log_print(f"  HVG count: {adata_bbknn.var['highly_variable'].sum()}")

    # PCA
    log_print(f"\nRunning PCA (n_comps={n_pcs})...")
    sc.tl.pca(adata_bbknn, n_comps=n_pcs, use_highly_variable=True)

    # BBKNN
    log_print(f"\nRunning BBKNN (neighbors_within_batch={neighbors_within_batch})...")
    bbknn.bbknn(
        adata_bbknn,
        batch_key=batch_key,
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
        trim=None
    )

    # UMAP
    log_print("Computing UMAP...")
    sc.tl.umap(adata_bbknn)

    # Leiden clustering
    log_print("Running Leiden clustering...")
    sc.tl.leiden(adata_bbknn, resolution=1.0)

    log_print("\nBBKNN integration complete!")
    log_print(f"  Final shape: {adata_bbknn.shape}")
    log_print(f"  Clusters: {adata_bbknn.obs['leiden'].nunique()}")

    return adata_bbknn

def main():
    log_print("="*80)
    log_print("PER-SAMPLE DENORMALIZE COUNTS AND BBKNN RE-INTEGRATION TEST")
    log_print("="*80)
    log_print(f"\nInput file: {INPUT_FILE}")
    log_print(f"Output directory: {OUTPUT_DIR}")

    # Load data
    log_print("\n" + "="*60)
    log_print("Loading data...")
    log_print("="*60)
    adata = sc.read_h5ad(INPUT_FILE)
    log_print(f"Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
    log_print(f"Layers: {list(adata.layers.keys())}")
    log_print(f"Obs columns: {list(adata.obs.columns)[:10]}...")

    # Check for batch key
    batch_key = None
    for possible_key in ['dataset', 'sample', 'batch', 'donor']:
        if possible_key in adata.obs.columns:
            batch_key = possible_key
            log_print(f"Found batch key: '{batch_key}' ({adata.obs[batch_key].nunique()} batches)")
            break

    if batch_key is None:
        log_print("WARNING: No batch key found! Using first categorical column...")
        categorical_cols = adata.obs.select_dtypes(include=['category', 'object']).columns
        if len(categorical_cols) > 0:
            batch_key = categorical_cols[0]
            log_print(f"Using: '{batch_key}'")
        else:
            log_print("ERROR: No suitable batch key found!")
            return

    # Check counts layer
    if 'counts' not in adata.layers:
        log_print("\nWARNING: 'counts' layer not found!")
        log_print("Available layers:", list(adata.layers.keys()))

        # Try to find counts elsewhere
        if adata.raw is not None:
            log_print("\nFound .raw.X - copying to layers['counts']")
            adata.layers['counts'] = adata.raw.X.copy()
        elif 'X' in dir(adata):
            log_print("\nWARNING: Using .X as counts (may not be raw!)")
            adata.layers['counts'] = adata.X.copy()
        else:
            log_print("FATAL: Cannot find count data!")
            return

    # **PER-SAMPLE DENORMALIZATION**
    per_sample_df, any_denormalized = process_sample_denormalization(adata, batch_key=batch_key)

    if per_sample_df is None:
        log_print("ERROR: Per-sample processing failed!")
        return

    # Create per-sample visualizations
    visualize_per_sample_stats(per_sample_df, OUTPUT_DIR)

    # Save denormalized data (if any changes were made)
    if any_denormalized:
        log_print("\n" + "="*60)
        log_print("Saving denormalized data...")
        log_print("="*60)
        output_h5ad = OUTPUT_DIR / "adata_denormalized.h5ad"
        adata.write_h5ad(output_h5ad)
        log_print(f"  Saved: {output_h5ad}")
    else:
        log_print("\n✓ All samples already had raw counts - no denormalization needed")

    # Run BBKNN integration with corrected/original counts
    log_print("\n" + "="*60)
    log_print("BBKNN RE-INTEGRATION")
    log_print("="*60)

    adata_integrated = run_bbknn_integration(
        adata,
        batch_key=batch_key,
        n_pcs=N_PCS,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH
    )

    if adata_integrated is None:
        log_print("ERROR: BBKNN integration failed!")
        return

    # Save integrated result
    log_print("\nSaving BBKNN integrated data...")
    output_bbknn = OUTPUT_DIR / "adata_bbknn_reintegrated.h5ad"
    adata_integrated.write_h5ad(output_bbknn)
    log_print(f"  Saved: {output_bbknn}")

    # Create visualization
    log_print("\n" + "="*60)
    log_print("Creating visualization...")
    log_print("="*60)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # UMAP by batch
    sc.pl.umap(adata_integrated, color=batch_key, ax=axes[0], show=False, title='UMAP - Batch')

    # UMAP by cluster
    sc.pl.umap(adata_integrated, color='leiden', ax=axes[1], show=False, title='UMAP - Leiden Clusters')

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'bbknn_umap.png', dpi=150, bbox_inches='tight')
    plt.close()
    log_print("  Saved: bbknn_umap.png")

    # Final report
    log_print("\n" + "="*80)
    log_print("PROCESSING COMPLETE!")
    log_print("="*80)
    log_print(f"\nOutput directory: {OUTPUT_DIR}")
    log_print("\nGenerated files:")
    log_print(f"  1. per_sample_denormalization_stats.csv - Per-sample statistics")
    log_print(f"  2. per_sample_comparison.png - Per-sample comparison plots")
    if any_denormalized:
        log_print(f"  3. adata_denormalized.h5ad - Denormalized counts")
    log_print(f"  4. adata_bbknn_reintegrated.h5ad - Re-integrated BBKNN data")
    log_print(f"  5. bbknn_umap.png - UMAP visualization")
    log_print(f"  6. {LOG_FILE.name} - Full log")

    log_print("\n" + "="*80)
    log_print("RECOMMENDATIONS")
    log_print("="*80)

    if any_denormalized:
        n_denorm = per_sample_df['denormalized'].sum()
        log_print(f"\n✓ SUCCESS: Denormalized {n_denorm} samples")
        log_print("  → Use adata_bbknn_reintegrated.h5ad for downstream scVI/scANVI pipelines")
        log_print("  → Review per_sample_denormalization_stats.csv for details")
    else:
        log_print("\n✓ All samples already had raw counts format")
        log_print("  → Original data is suitable for scVI/scANVI pipelines")

    log_print("\n" + "="*80)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log_print(f"\n{'='*80}")
        log_print(f"FATAL ERROR: {str(e)}")
        log_print(f"{'='*80}")
        import traceback
        log_print("\nFull traceback:")
        log_print(traceback.format_exc())
        sys.exit(1)
