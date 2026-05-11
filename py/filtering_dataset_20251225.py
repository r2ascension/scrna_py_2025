# ===== Remove Specific Datasets from h5ad =====
# Purpose: Filter out Ji_Hoon_Ahn_2021 and Seon_Pyo_Hong_2023 datasets
# Author: r2end
# Date: 2024-12-25

import scanpy as sc
import numpy as np
import gc

# ===== Configuration =====
INPUT_PATH = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad"
OUTPUT_PATH = "/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected_filtered.h5ad"
DATASET_COL = 'dataset'  # Column name containing dataset information
REMOVE_DATASETS = ['Ji_Hoon_Ahn_2021', 'Seon_Pyo_Hong_2023']

# ===== Load Data =====
print(f"Loading data from {INPUT_PATH}")
adata = sc.read_h5ad(INPUT_PATH)
print(f"Original data: {adata.n_obs} cells, {adata.n_vars} genes")

# Check dataset distribution
print(f"\nDataset distribution before filtering:")
print(adata.obs[DATASET_COL].value_counts())

# ===== Filter Datasets =====
# Create boolean mask for cells to keep
keep_mask = ~adata.obs[DATASET_COL].isin(REMOVE_DATASETS)
removed_count = (~keep_mask).sum()

print(f"\nRemoving {removed_count} cells from datasets: {REMOVE_DATASETS}")

# Subset data
adata_filtered = adata[keep_mask].copy()

# Clean up
del adata
gc.collect()

# ===== Verify Results =====
print(f"\nFiltered data: {adata_filtered.n_obs} cells, {adata_filtered.n_vars} genes")
print(f"\nDataset distribution after filtering:")
print(adata_filtered.obs[DATASET_COL].value_counts())

# Check if removed datasets are completely gone
remaining_removed = adata_filtered.obs[DATASET_COL].isin(REMOVE_DATASETS).sum()
if remaining_removed == 0:
    print(f"\n✓ Successfully removed all cells from {REMOVE_DATASETS}")
else:
    print(f"\n⚠️ Warning: {remaining_removed} cells from removed datasets still remain")

# ===== Save Results =====
print(f"\nSaving filtered data to {OUTPUT_PATH}")
adata_filtered.write_h5ad(OUTPUT_PATH, compression='gzip', compression_opts=9)
print("Done!")