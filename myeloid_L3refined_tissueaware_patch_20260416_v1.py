#!/usr/bin/env python3
"""Myeloid L3 refined tissue-aware patch v1.

This script preserves the validated 2026-03-29 patch logic and adds one new rule:
all alveolar macrophage labels outside `lung parenchyma` and `respiratory airway`
are reassigned to `Interstitial macrophages`.
"""

import os
from pathlib import Path
import warnings

import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

INPUT_PATH = Path(
    os.getenv(
        "MYELOID_PATCH_INPUT_H5AD",
        "/home/h2048/data/py/0416/myeloid_validation_optimized_cpu_rerun/adata_myeloid_refined_FINAL.h5ad",
    )
)
OUTPUT_PATH = Path(
    os.getenv(
        "MYELOID_PATCH_OUTPUT_H5AD",
        "/home/h2048/data/py/0416/adata_myeloid_L3refined_tissueaware_patched_v1.h5ad",
    )
)
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

COL_L3 = "cell_type_L3_refined"
COL_L2 = "cell_type_L2"
COL_TISSUE = "tissue"
AUDIT_COL = "alveolar_to_interstitial_tissue_rule_applied"

RESPIRATORY_ALVEOLAR_TISSUES = {"lung parenchyma", "respiratory airway"}
ALVEOLAR_LABELS = {
    "Resident Alveolar macrophages",
    "Resting Alveolar macrophages",
}
KEEP_INTERSTITIAL_LABELS = {
    "CD163L1+ Interstitial macrophages",
    "Inflammatory Interstitial macrophages",
}

print(f"Input : {INPUT_PATH}")
print(f"Output: {OUTPUT_PATH}")
if not INPUT_PATH.exists():
    raise FileNotFoundError(f"Input h5ad not found: {INPUT_PATH}")

adata = sc.read_h5ad(INPUT_PATH)
print(f"Loaded: {adata.shape}")

required_cols = {COL_L3, COL_TISSUE}
missing = sorted(required_cols.difference(adata.obs.columns))
if missing:
    raise KeyError(f"Missing required obs columns: {missing}")

print(f"\nBefore patch — {COL_L3} value counts:")
print(adata.obs[COL_L3].value_counts(dropna=False).to_string())

# ---------------------------------------------------------------------------
# 1. Safety checks inherited from the 2026-03-29 patch
# ---------------------------------------------------------------------------
nan_mask = adata.obs[COL_L3].isna()
if nan_mask.any():
    if COL_L2 not in adata.obs.columns:
        raise KeyError(f"{COL_L2} is required when {COL_L3} contains NaN values")
    l2_of_nan = adata.obs.loc[nan_mask, COL_L2].astype(str).unique().tolist()
    assert l2_of_nan == ["Non-classical monocytes"], (
        f"Unexpected {COL_L2} labels under NaN {COL_L3}: {l2_of_nan}"
    )
    print(f"[OK] All {nan_mask.sum()} NaN cells are Non-classical monocytes in {COL_L2}")
else:
    print("[OK] No NaN labels found in cell_type_L3_refined")

lq_mask = adata.obs[COL_L3].astype(str) == "Low-quality Interstitial macrophages"
print(f"[OK] Low-quality Interstitial macrophages to drop: {int(lq_mask.sum())} cells")

# ---------------------------------------------------------------------------
# 2. Apply patch on mutable object-dtype vectors
# ---------------------------------------------------------------------------
col_l3 = adata.obs[COL_L3].astype(object).copy()
if COL_L2 in adata.obs.columns:
    col_l2 = adata.obs[COL_L2].astype(object).copy()
else:
    col_l2 = None

tissue_series = adata.obs[COL_TISSUE].astype(str)
rule_mask = col_l3.isin(ALVEOLAR_LABELS) & (~tissue_series.isin(RESPIRATORY_ALVEOLAR_TISSUES))

# 2a. Fill NaN exactly as in the validated 2026-03-29 patch
if nan_mask.any():
    col_l3.loc[nan_mask] = "Non-classical monocytes"
    print(f"[2a] Filled {int(nan_mask.sum())} NaN -> 'Non-classical monocytes'")

# 2b. New tissue-aware relabel
col_l3.loc[rule_mask] = "Interstitial macrophages"
if col_l2 is not None:
    col_l2.loc[rule_mask] = "Interstitial macrophages"
print(
    "[2b] Tissue-aware relabel: "
    f"{int(rule_mask.sum())} non-respiratory alveolar macrophages -> Interstitial macrophages"
)

# 2c. Collapse interstitial macrophage subtypes (validated 2026-03-29 behavior)
collapse_mask = (
    pd.Series(col_l3, index=adata.obs_names).astype(str).str.contains("Interstitial macrophages", na=False)
    & (~pd.Series(col_l3, index=adata.obs_names).isin(KEEP_INTERSTITIAL_LABELS))
)
n_collapsed = int(collapse_mask.sum())
col_l3.loc[collapse_mask] = "Interstitial macrophages"
print(f"[2c] Collapsed {n_collapsed} cells -> 'Interstitial macrophages'")

# 2d. Drop low-quality cells using the ORIGINAL low-quality mask
keep_mask = ~lq_mask
print(f"[2d] Dropping {int(lq_mask.sum())} Low-quality Interstitial macrophages cells")

# Write back before subsetting
adata.obs[COL_L3] = pd.Categorical(col_l3)
adata.obs[AUDIT_COL] = pd.Categorical(rule_mask.map({True: "yes", False: "no"}), categories=["no", "yes"])
if col_l2 is not None:
    adata.obs[COL_L2] = pd.Categorical(col_l2)

adata = adata[keep_mask].copy()
print(f"Shape after drop: {adata.shape}")

# ---------------------------------------------------------------------------
# 3. Verify
# ---------------------------------------------------------------------------
print(f"\nAfter patch — {COL_L3} value counts:")
print(adata.obs[COL_L3].value_counts(dropna=False).to_string())

assert adata.obs[COL_L3].isna().sum() == 0, "NaN still present in cell_type_L3_refined"
assert "Low-quality Interstitial macrophages" not in adata.obs[COL_L3].astype(str).values, (
    "Low-quality Interstitial macrophages still present"
)

post_rule_mask = adata.obs[COL_L3].astype(str).isin(ALVEOLAR_LABELS) & (~adata.obs[COL_TISSUE].astype(str).isin(RESPIRATORY_ALVEOLAR_TISSUES))
assert int(post_rule_mask.sum()) == 0, (
    "Found non-respiratory Alveolar macrophages after tissue-aware patch"
)

remaining_interstitial_labels = sorted(
    x for x in adata.obs[COL_L3].astype(str).unique().tolist() if "Interstitial macrophages" in x
)
print("\nRemaining interstitial-related labels:")
for label in remaining_interstitial_labels:
    print(f"  {label}")

print("\nAudit summary (alveolar -> interstitial rule):")
print(adata.obs[AUDIT_COL].value_counts(dropna=False).to_string())

if COL_L2 in adata.obs.columns:
    print(f"\nUpdated {COL_L2} value counts:")
    print(adata.obs[COL_L2].value_counts(dropna=False).to_string())

adata.write_h5ad(OUTPUT_PATH, compression="gzip", compression_opts=9)
print(f"\nSaved: {OUTPUT_PATH}")
print(f"Final shape: {adata.shape}")
