# De-lognormalization Integration Summary

## Overview

This document summarizes the integration of de-lognormalization functions into the QC pipeline script (`rds_folder_qc_smartmerge_20251213_v2_with_denorm.R`).

## What Was Added

### 1. Core De-lognormalization Function

**`restore_raw_counts_log(normalized_data, size_factors)`**

- Reverses log(counts + 1) normalization to recover raw counts
- Formula: `counts = expm1(normalized_data) * size_factors`
- Uses `expm1()` for numerical stability instead of `exp() - 1`

```r
# Example usage:
raw_counts <- restore_raw_counts_log(normalized_data, size_factors)
```

### 2. Seurat-Specific Wrapper Function

**`restore_counts_from_seurat(seurat_obj, assay = "RNA", scale_factor = 10000)`**

- Extracts log-normalized data from Seurat object
- Uses `nCount_RNA` (total UMI per cell) to calculate size factors
- Reverses Seurat's LogNormalize: `log1p(counts / sum(counts) * scale_factor)`
- Returns rounded integer counts in dgCMatrix format

```r
# Example usage:
raw_counts <- restore_counts_from_seurat(seurat_obj, assay = "RNA", scale_factor = 10000)
```

### 3. Modified DecontX Function

**`run_decontx()` - Enhanced with De-lognormalization Support**

Key changes:
- Detects if counts are log-normalized (non-integer)
- If `TRY_DENORMALIZE = TRUE`, attempts to restore raw counts
- Validates restored counts are integer-like before proceeding
- Adds metadata flag `data_type` to track data processing:
  - `"raw_counts"`: Original data was already raw
  - `"restored_counts"`: Successfully restored from log-normalized
  - `"denorm_failed"`: De-lognormalization produced non-integer counts
  - `"denorm_error"`: De-lognormalization threw an error
  - `"non_raw_counts"`: Skipped (TRY_DENORMALIZE = FALSE)

## New Parameters

Added to the PARAMETERS section:

```r
# De-lognormalization (NEW!)
TRY_DENORMALIZE <- TRUE          # Attempt to restore raw counts if data is log-normalized
DENORM_SCALE_FACTOR <- 10000     # Scale factor used in LogNormalize
```

### TRY_DENORMALIZE
- **Default**: `TRUE`
- **Purpose**: Enable automatic de-lognormalization when counts are not integer-like
- **When to set FALSE**: When you want to skip DecontX for log-normalized data

### DENORM_SCALE_FACTOR
- **Default**: `10000`
- **Purpose**: Scale factor used in Seurat's `NormalizeData(scale.factor = 10000)`
- **When to change**: If your data used a different scale factor

## How It Works

### Workflow for DecontX with De-lognormalization:

1. **Extract counts** from Seurat object
2. **Check if integer-like** using `is_integer_like_counts()`
3. **If NOT integer-like AND TRY_DENORMALIZE = TRUE:**
   - Extract log-normalized data from "data" slot
   - Extract nCount_RNA from metadata
   - Calculate size factors: `nCount_RNA / scale_factor`
   - Restore raw counts: `expm1(normalized_data) * size_factors`
   - Round to integers
   - Validate restored counts are integer-like
4. **If validation passes:**
   - Run DecontX on restored counts
   - Store contamination scores
5. **If validation fails:**
   - Skip DecontX
   - Store NA for contamination
   - Flag in metadata

## Mathematical Details

### Seurat's LogNormalize (default)

```r
normalized = log1p(counts / sum(counts) * 10000)
```

Where:
- `counts`: Raw UMI counts per gene
- `sum(counts)`: Total UMI per cell (= nCount_RNA)
- `10000`: Scale factor

### De-lognormalization

To reverse:

```r
counts = expm1(normalized) * sum(counts) / 10000
```

Implementation:
```r
size_factors <- nCount_RNA / 10000
raw_counts <- expm1(normalized_data) * size_factors
raw_counts <- round(raw_counts)  # Convert to integers
```

## Usage Examples

### Example 1: Standard Usage (with de-lognormalization enabled)

```bash
Rscript rds_folder_qc_smartmerge_20251213_v2_with_denorm.R \
  /data/input_rds_folder \
  /data/output_folder
```

- If input data has raw counts → DecontX runs normally
- If input data is log-normalized → Automatically restores raw counts → DecontX runs

### Example 2: Disable De-lognormalization

Edit the script to set:
```r
TRY_DENORMALIZE <- FALSE
```

Then run:
```bash
Rscript rds_folder_qc_smartmerge_20251213_v2_with_denorm.R \
  /data/input_rds_folder \
  /data/output_folder
```

- If input data has raw counts → DecontX runs normally
- If input data is log-normalized → DecontX skipped (contamination = NA)

### Example 3: Custom Scale Factor

If your data was normalized with a different scale factor (e.g., 1e6):

Edit the script:
```r
DENORM_SCALE_FACTOR <- 1e6
```

## Output Metadata

The pipeline adds a new metadata column `data_type` to track processing status:

| Value | Meaning |
|-------|---------|
| `raw_counts` | Original counts were already raw integer counts |
| `restored_counts` | Raw counts successfully restored from log-normalized data |
| `denorm_failed` | De-lognormalization produced non-integer counts (validation failed) |
| `denorm_error` | De-lognormalization threw an error |
| `non_raw_counts` | Data is log-normalized and TRY_DENORMALIZE = FALSE |

### Check Data Types After Processing

```r
# Load merged object
merged <- readRDS("output_folder/merged/merged_seurat_final.rds")

# Check data types
table(merged$data_type)

# Check by sample
table(merged$sample, merged$data_type)
```

## Log Messages

The pipeline now outputs detailed messages during DecontX processing:

**When counts are already raw:**
```
[DecontX] data_type: raw_counts
```

**When de-lognormalization is attempted:**
```
[DecontX] ⚠️  Counts are NOT integer-like (normalized/scaled data detected)
[DecontX] Attempting to restore raw counts via de-lognormalization...
[DecontX] ✓ Successfully restored integer-like counts
```

**When de-lognormalization fails:**
```
[DecontX] ⚠️  Counts are NOT integer-like (normalized/scaled data detected)
[DecontX] Attempting to restore raw counts via de-lognormalization...
[DecontX] ✗ Restored counts are still not integer-like → SKIPPING
```

## Version History

### v2_DENORM (2024-12-14)
- ⭐ Added `restore_raw_counts_log()` function
- ⭐ Added `restore_counts_from_seurat()` function
- ⭐ Modified `run_decontx()` to support de-lognormalization
- ⭐ Added `TRY_DENORMALIZE` and `DENORM_SCALE_FACTOR` parameters
- ⭐ Added `data_type` metadata tracking

### v5_HOTFIX_v6 (2024-12-14)
- DoubletFinder command dependency fix
- ensure_normalize_command function
- Assays() container, dimnames, integer detection fixes

## Validation

To verify the de-lognormalization is working correctly:

```r
# Load a sample
obj <- readRDS("cleaned_samples/sample1_cleaned.rds")

# Check if de-lognormalization was used
table(obj$data_type)

# If data_type = "restored_counts", validation passed:
# 1. Restored counts are integer-like
# 2. DecontX ran successfully
# 3. Contamination scores are available

# Check contamination scores
summary(obj$decontX_contamination)
```

## Notes

1. **Numerical Stability**: Uses `expm1()` instead of `exp() - 1` for better numerical precision
2. **Memory Efficiency**: Maintains sparse matrix format (dgCMatrix) throughout
3. **Validation**: Always validates restored counts before using them
4. **Error Handling**: Gracefully handles failures and continues processing
5. **Metadata Tracking**: `data_type` column tracks processing status for QC

## Comparison with Original Script

| Feature | Original (v5_HOTFIX_v6) | New (v2_DENORM) |
|---------|------------------------|-----------------|
| Handle raw counts | ✓ | ✓ |
| Handle log-normalized data | ✗ (skip DecontX) | ✓ (restore counts) |
| De-lognormalization | ✗ | ✓ |
| Data type tracking | ✗ | ✓ |
| Configurable scale factor | ✗ | ✓ |

## Recommendations

1. **Use TRY_DENORMALIZE = TRUE** (default) for maximum compatibility
2. **Verify DENORM_SCALE_FACTOR** matches your normalization settings
3. **Check data_type** in output to understand processing status
4. **Review QC stats** to see how many samples used de-lognormalization

## Questions?

If you encounter issues:
1. Check log messages for de-lognormalization status
2. Verify `data_type` metadata
3. Check if scale factor matches your normalization
4. Review QC stats CSV for per-sample processing details
