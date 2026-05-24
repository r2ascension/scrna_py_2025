#!/usr/bin/env python3
"""
Quick Counts Layer Inspection (PER-SAMPLE)
===========================================

Quick diagnostic script to check if layers['counts'] contains
raw counts or normalized data on a PER-SAMPLE basis WITHOUT running
full denormalization.

Usage:
------
python quick_counts_inspection.py
"""

import scanpy as sc
import numpy as np
import pandas as pd
from scipy import sparse

INPUT_FILE = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"

def inspect_counts(data, name="counts", sample_size=10000):
    """Quick inspection of count data"""
    # Get sample
    if sparse.issparse(data):
        data_values = data.data
        if len(data_values) > sample_size:
            indices = np.random.choice(len(data_values), sample_size, replace=False)
            sample = data_values[indices]
        else:
            sample = data_values
    else:
        flat = data.flatten()
        if len(flat) > sample_size:
            sample = np.random.choice(flat, sample_size, replace=False)
        else:
            sample = flat

    # Remove zeros
    sample_nonzero = sample[sample > 0]

    if len(sample_nonzero) == 0:
        return None, "UNKNOWN"

    # Statistics
    stats = {
        'shape': data.shape,
        'sparse': sparse.issparse(data),
        'dtype': str(data.dtype),
        'n_nonzero': len(sample_nonzero),
        'min': np.min(sample_nonzero),
        'max': np.max(sample_nonzero),
        'mean': np.mean(sample_nonzero),
        'median': np.median(sample_nonzero),
        'std': np.std(sample_nonzero),
        'pct_integer': np.mean(np.abs(sample_nonzero - np.round(sample_nonzero)) < 1e-6) * 100,
        'sample_values': sample_nonzero[:20]
    }

    # Diagnosis
    is_integer = stats['pct_integer'] > 95
    is_high_max = stats['max'] > 50
    is_low_max = stats['max'] < 15
    is_low_mean = stats['mean'] < 5

    if is_integer and is_high_max:
        diagnosis = "RAW"
    elif not is_integer and is_low_max and is_low_mean:
        diagnosis = "NORMALIZED"
    else:
        diagnosis = "UNCERTAIN"

    return stats, diagnosis

def print_sample_stats(sample_name, stats, diagnosis, n_cells):
    """Print formatted sample statistics"""
    print(f"\n  Sample: {sample_name} ({n_cells} cells)")
    print(f"    Detection: {diagnosis}")
    print(f"    Max: {stats['max']:.2f}, Mean: {stats['mean']:.2f}, % Integer: {stats['pct_integer']:.1f}%")

    if diagnosis == "RAW":
        print(f"    ✓ Raw counts detected")
    elif diagnosis == "NORMALIZED":
        print(f"    ⚠️  Log-normalized detected - NEEDS DENORMALIZATION!")
    else:
        print(f"    ⚠️  Uncertain - review manually")

def main():
    print("="*80)
    print("QUICK COUNTS LAYER INSPECTION (PER-SAMPLE)")
    print("="*80)
    print(f"\nLoading: {INPUT_FILE}")

    adata = sc.read_h5ad(INPUT_FILE)
    print(f"\nDataset: {adata.shape[0]} cells × {adata.shape[1]} genes")
    print(f"Layers: {list(adata.layers.keys())}")

    # Find batch key
    batch_key = None
    for possible_key in ['dataset', 'sample', 'batch', 'donor']:
        if possible_key in adata.obs.columns:
            batch_key = possible_key
            break

    if batch_key is None:
        print("\n⚠️  WARNING: No batch key found, checking global data only")
        # Check global counts layer
        if 'counts' in adata.layers:
            print("\n" + "="*80)
            print("Inspecting layers['counts'] (GLOBAL)")
            print("="*80)
            stats, diagnosis = inspect_counts(adata.layers['counts'], name="layers['counts']")

            if stats:
                print(f"\nShape: {stats['shape']}")
                print(f"Sparse: {stats['sparse']}")
                print(f"Data type: {stats['dtype']}")
                print(f"\nNon-zero value statistics (n={stats['n_nonzero']}):")
                print(f"  Min:      {stats['min']:.6f}")
                print(f"  Max:      {stats['max']:.6f}")
                print(f"  Mean:     {stats['mean']:.6f}")
                print(f"  Median:   {stats['median']:.6f}")
                print(f"  Std:      {stats['std']:.6f}")
                print(f"  % Integer: {stats['pct_integer']:.2f}%")
                print(f"\nSample values (first 20 non-zero):")
                print(stats['sample_values'])

                print(f"\n{'='*80}")
                print("DIAGNOSIS:")
                print(f"{'='*80}")
                if diagnosis == "RAW":
                    print("✓ Likely RAW COUNTS")
                    print("  Evidence: High max values, integer-like")
                elif diagnosis == "NORMALIZED":
                    print("⚠️  Likely LOG-NORMALIZED")
                    print("  Evidence: Low max values, non-integer, low mean")
                    print("  → DENORMALIZATION NEEDED for scVI/scANVI!")
                else:
                    print("⚠️  UNCERTAIN - review statistics")
        else:
            print("\n⚠️  WARNING: 'counts' layer not found!")
            print(f"Available layers: {list(adata.layers.keys())}")
        return

    # Per-sample inspection
    print(f"\n{'='*80}")
    print(f"PER-SAMPLE INSPECTION (batch key: '{batch_key}')")
    print(f"{'='*80}")

    samples = adata.obs[batch_key].unique()
    n_samples = len(samples)
    print(f"\nFound {n_samples} unique samples")

    if 'counts' not in adata.layers:
        print("\n⚠️  ERROR: 'counts' layer not found!")
        print(f"Available layers: {list(adata.layers.keys())}")

        if adata.raw is not None:
            print("\nFound .raw.X, checking that instead...")
            adata.layers['counts'] = adata.raw.X.copy()
        else:
            return

    # Inspect each sample
    results = []
    samples_raw = []
    samples_normalized = []
    samples_uncertain = []

    print(f"\n{'-'*80}")
    print("Analyzing samples...")
    print(f"{'-'*80}")

    # Show progress for large datasets
    show_all = n_samples <= 20
    samples_to_show = samples if show_all else samples[:10]

    for idx, sample in enumerate(samples):
        sample_mask = adata.obs[batch_key] == sample
        n_cells = sample_mask.sum()
        sample_counts = adata.layers['counts'][sample_mask, :]

        stats, diagnosis = inspect_counts(sample_counts, name=f"Sample {sample}")

        if stats is None:
            continue

        results.append({
            'sample': sample,
            'n_cells': n_cells,
            'diagnosis': diagnosis,
            'max': stats['max'],
            'mean': stats['mean'],
            'pct_integer': stats['pct_integer']
        })

        if diagnosis == "RAW":
            samples_raw.append(sample)
        elif diagnosis == "NORMALIZED":
            samples_normalized.append(sample)
        else:
            samples_uncertain.append(sample)

        # Print details for first 10 samples or all if <=20
        if sample in samples_to_show:
            print_sample_stats(sample, stats, diagnosis, n_cells)

    if not show_all:
        print(f"\n  ... and {n_samples - 10} more samples")
        print(f"  (showing first 10, see summary below for all)")

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"\nTotal samples analyzed: {n_samples}")
    print(f"  - Raw counts: {len(samples_raw)} samples ({len(samples_raw)/n_samples*100:.1f}%)")
    print(f"  - Normalized: {len(samples_normalized)} samples ({len(samples_normalized)/n_samples*100:.1f}%)")
    print(f"  - Uncertain: {len(samples_uncertain)} samples ({len(samples_uncertain)/n_samples*100:.1f}%)")

    if len(samples_normalized) > 0:
        print(f"\n⚠️  NORMALIZED SAMPLES ({len(samples_normalized)}):")
        for s in samples_normalized[:20]:  # Show first 20
            print(f"    - {s}")
        if len(samples_normalized) > 20:
            print(f"    ... and {len(samples_normalized) - 20} more")

    if len(samples_uncertain) > 0:
        print(f"\n⚠️  UNCERTAIN SAMPLES ({len(samples_uncertain)}):")
        for s in samples_uncertain[:20]:
            print(f"    - {s}")
        if len(samples_uncertain) > 20:
            print(f"    ... and {len(samples_uncertain) - 20} more")

    # Save summary table
    if len(results) > 0:
        df = pd.DataFrame(results)
        output_file = 'per_sample_inspection_summary.csv'
        df.to_csv(output_file, index=False)
        print(f"\n✓ Detailed results saved to: {output_file}")

    # Recommendations
    print(f"\n{'='*80}")
    print("RECOMMENDATIONS")
    print(f"{'='*80}")

    if len(samples_normalized) > 0:
        print(f"\n⚠️  {len(samples_normalized)} samples contain NORMALIZED data!")
        print("\nNext steps:")
        print("  1. Run: python denormalize_counts_bbknn_test.py")
        print("  2. This will:")
        print("     - Denormalize each sample independently")
        print("     - Formula: raw = exp(log_normalized) - 1")
        print("     - Re-run BBKNN integration")
        print("     - Save corrected data for scVI/scANVI pipelines")
    elif len(samples_uncertain) > 0 and len(samples_raw) == 0:
        print(f"\n⚠️  All {len(samples_uncertain)} samples are UNCERTAIN")
        print("  → Review per_sample_inspection_summary.csv")
        print("  → Consider running denormalize_counts_bbknn_test.py as a test")
    elif len(samples_raw) == n_samples:
        print("\n✓ All samples contain RAW counts")
        print("  → Data is suitable for scVI/scANVI pipelines")
    else:
        print(f"\n⚠️  Mixed data types detected:")
        print(f"    - {len(samples_raw)} raw")
        print(f"    - {len(samples_normalized)} normalized")
        print(f"    - {len(samples_uncertain)} uncertain")
        print("\n  → Run denormalize_counts_bbknn_test.py to correct normalized samples")

    print(f"\n{'='*80}")

if __name__ == "__main__":
    main()
