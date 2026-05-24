# Counts Denormalization and BBKNN Re-integration Scripts (PER-SAMPLE)

## Problem Statement

scVI and scANVI models **REQUIRE raw integer UMI counts** in `layers['counts']`, but BBKNN preprocessing may have accidentally stored **log-normalized counts** instead. This causes:

- Poor scVI/scANVI training performance
- Incorrect batch correction
- Invalid cell type annotations

**CRITICAL:** Different samples/batches may have been preprocessed differently! Some samples might have raw counts while others have normalized counts. These scripts check and denormalize **EACH SAMPLE INDEPENDENTLY**.

## Detection Heuristics

### Raw Counts (Correct)
- **Integer values**: >95% match floor
- **High max values**: typically >100
- **Wide range**: 0-10,000+
- **Example**: `[0, 1, 2, 5, 150, 0, 3, ...]`

### Log-Normalized Counts (Incorrect for scVI)
- **Float values**: continuous
- **Low max values**: typically <10
- **Narrow range**: 0-8
- **Example**: `[0, 0.693, 1.099, 1.792, 5.011, 0, ...]`

## Scripts Overview

### 1. Quick Inspection (Recommended First Step)

```bash
python quick_counts_inspection.py
```

**What it does:**
- Loads the h5ad file
- **Iterates through each unique sample/batch**
- For each sample:
  - Samples 10,000 non-zero values from `layers['counts']`
  - Calculates statistics (min, max, mean, % integer)
  - Diagnoses data type (RAW vs LOG-NORMALIZED)
- Generates per-sample summary
- Saves `per_sample_inspection_summary.csv`
- Takes ~1-2 minutes

**Output:**
```
PER-SAMPLE INSPECTION (batch key: 'dataset')
================================================================================
Found 42 unique samples

  Sample: GSM123456 (8234 cells)
    Detection: RAW
    Max: 1234.00, Mean: 12.34, % Integer: 98.5%
    ✓ Raw counts detected

  Sample: GSM123457 (5678 cells)
    Detection: NORMALIZED
    Max: 6.21, Mean: 1.85, % Integer: 12.3%
    ⚠️  Log-normalized detected - NEEDS DENORMALIZATION!

SUMMARY
================================================================================
Total samples analyzed: 42
  - Raw counts: 20 samples (47.6%)
  - Normalized: 18 samples (42.9%)
  - Uncertain: 4 samples (9.5%)
```

### 2. Full Denormalization and BBKNN Re-integration

```bash
python denormalize_counts_bbknn_test.py
```

**What it does:**

1. **Detect PER SAMPLE** if `layers['counts']` contains normalized data
   - Iterates through each unique sample/batch
   - Checks normalization status independently
   - Handles mixed datasets (some samples raw, some normalized)

2. **Denormalize PER SAMPLE** if needed:
   - Formula: `raw_counts = exp(log_normalized) - 1`
   - Optionally apply size factor correction per sample
   - Round to integers
   - **Only denormalizes samples that need it!**

3. **Re-run BBKNN integration** with corrected counts:
   - Normalize/log-transform for PCA (temporary)
   - Run PCA → BBKNN → UMAP → Leiden clustering

4. **Save results**:
   - `per_sample_denormalization_stats.csv` - Detailed per-sample stats
   - `per_sample_comparison.png` - Before/after comparison plots
   - `adata_denormalized.h5ad` - Denormalized data (if any changes)
   - `adata_bbknn_reintegrated.h5ad` - Re-integrated with BBKNN
   - `bbknn_umap.png` - UMAP visualization
   - `denormalization_test_*.log` - Full log

**Expected runtime:** 30-60 minutes for 400k cells

## Denormalization Formula (from R example)

```python
# Original transformation (BBKNN):
log_normalized = log(raw_counts + 1)

# Reverse transformation (this script):
raw_counts = exp(log_normalized) - 1

# With size factor correction (if available):
raw_counts = exp(log_normalized) - 1
raw_counts = raw_counts * size_factors
```

## Usage Workflow

### Step 1: Quick Check
```bash
cd /home/h2048/script/py
python quick_counts_inspection.py
```

### Step 2: If Normalized Data Detected
```bash
python denormalize_counts_bbknn_test.py
```

### Step 3: Use Corrected Data for scVI/scANVI
```python
# In your pipeline scripts (e.g., epithelial_scvi_celltypist_scanvi_*.py)

# OLD (potentially using normalized counts):
# adata = sc.read_h5ad("/home/h2048/data/py/1128/bbknn_annotation_analysis/adata_bbknn_annotated_corrected.h5ad")

# NEW (using denormalized counts):
adata = sc.read_h5ad("/home/h2048/data/py/20251213/denormalize_counts_test/adata_bbknn_reintegrated.h5ad")

# Verify counts are raw integers
assert adata.layers['counts'].data.max() > 50, "Counts seem normalized!"
```

## Output Directory Structure

```
/home/h2048/data/py/20251213/denormalize_counts_test/
├── per_sample_denormalization_stats.csv  # Per-sample detection & denorm results
├── per_sample_comparison.png             # Before/after comparison plots
├── adata_denormalized.h5ad               # Denormalized counts (if changes made)
├── adata_bbknn_reintegrated.h5ad         # FINAL - Use this for scVI/scANVI!
├── bbknn_umap.png                        # UMAP colored by batch & cluster
└── denormalization_test_20251213_*.log   # Detailed log
```

## Key Statistics to Review

After running, check `per_sample_denormalization_stats.csv`:

| Sample | n_cells | data_type | denormalized | original_max | denorm_max | original_pct_integer | denorm_pct_integer |
|--------|---------|-----------|--------------|--------------|------------|---------------------|-------------------|
| GSM123456 | 8234 | RAW_COUNTS | False | 1234.0 | - | 98.5% | - |
| GSM123457 | 5678 | LOG_NORMALIZED | True | 6.2 | 523.0 | 12.3% | 97.8% |

## Integration with Existing Pipelines

### Before (Potentially Wrong)
```python
# epithelial_scvi_celltypist_scanvi_20251212_v3_4.py
INPUT_PATH = "/home/h2048/data/py/1128/bbknn_celltype_analysis/Epithelial/adata_Epithelial_bbknn.h5ad"

# If layers['counts'] is normalized → scVI will fail!
scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='dataset')
```

### After (Corrected)
```python
# Option 1: Use re-integrated data
INPUT_PATH = "/home/h2048/data/py/20251213/denormalize_counts_test/adata_bbknn_reintegrated.h5ad"

# Option 2: Add validation check
adata = sc.read_h5ad(INPUT_PATH)
max_count = adata.layers['counts'].data.max() if sparse.issparse(adata.layers['counts']) else adata.layers['counts'].max()
if max_count < 50:
    raise ValueError(f"layers['counts'] appears normalized (max={max_count})! Run denormalization first.")
```

## Validation Checklist

After denormalization, verify:

- [ ] `per_sample_denormalization_stats.csv` shows successful denormalization for normalized samples
- [ ] Denormalized samples have max value >100 (raw counts range)
- [ ] Denormalized samples have >95% integer-like values
- [ ] UMAP shows good batch mixing (in `bbknn_umap.png`)
- [ ] Leiden clustering is reasonable (not over/under-clustered)
- [ ] `per_sample_comparison.png` shows clear before/after differences

## Common Issues

### Issue 1: "No batch key found"
**Solution:** Script auto-detects 'dataset', 'sample', 'batch', or 'donor'. If none exist, edit line 286:
```python
batch_key = 'your_batch_column_name'
```

### Issue 2: "Size factors seem unusual"
**Solution:** Script skips size factor correction if values are unreasonable. This is safe - denormalization alone usually works.

### Issue 3: "Denormalized data still looks normalized"
**Solution:** Original data may have been doubly-normalized. Try:
```python
# Manual double denormalization
raw = np.expm1(np.expm1(normalized))
```

## Performance Notes

- **Memory usage**: ~3-4x dataset size (keeps backups)
- **Runtime (400k cells)**:
  - Detection: ~1 min
  - Denormalization: ~2-3 min
  - BBKNN: ~20-40 min
  - Total: ~30-60 min
- **GPU**: Not required (BBKNN is CPU-based)

## References

- scVI documentation: https://docs.scvi-tools.org/
- BBKNN paper: Polański et al., Bioinformatics 2020
- R example: User-provided `restore_raw_counts_log()` function

## Contact

For issues or questions about this script, review:
- Full log file: `denormalization_test_*.log`
- CLAUDE.md for pipeline documentation
- Original BBKNN preprocessing scripts

---

**Created:** 2025-12-13
**Purpose:** Fix normalized counts in layers['counts'] for scVI/scANVI compatibility
