#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Complete Merge with Full Gene Preservation (FIXED v1.2.2)
========================================================
Fixes:
1) Filter QUERY genes with min_cells>=3 (using query.X as raw)
2) Select HVG for merge by overlap with reference HVG (after filter)
3) Concatenate with inner join
4) Standardize query columns:
   - Cell_Type_L2_final
   - mapping_confidence
5) Merged visualization: query panels colored ONLY on query (no NA confusion)

Design:
- .X: HVG only (shared HVGs, ref order)
- .layers['counts']: HVG counts (matches .X)
- .raw: FULL genes aligned to ORIGINAL reference full gene list (ref universe)
"""

import sys
import gc
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix, issparse, vstack as sp_vstack

print("=" * 70)
print("Complete Merge with Full Gene Preservation (FIXED v1.2.2)")
print("=" * 70)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ORIGINAL_REF_H5AD = "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
REF_WITH_UMAP_H5AD = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/reference_with_L2_umap.h5ad"
QUERY_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/query_mapped_L2.h5ad"

OUTPUT_H5AD = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/reference_plus_query_merged_fullgenes_L2_FIXED_v1_2_2.h5ad"
OUTPUT_DIR = "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3/merged_fullgenes_analysis_FIXED_v1_2_2"

# QUERY gene filtering
QUERY_MIN_CELLS = 3

DPI = 300
FIGURE_FORMAT = "pdf"

output_dir = Path(OUTPUT_DIR)
(output_dir / "figures").mkdir(parents=True, exist_ok=True)

sc.settings.set_figure_params(dpi=DPI, facecolor="white", format=FIGURE_FORMAT)

print("\nConfiguration:")
print(f"  Original reference (FULL): {ORIGINAL_REF_H5AD}")
print(f"  Reference with UMAP (HVG): {REF_WITH_UMAP_H5AD}")
print(f"  Query mapped (FULL in X):  {QUERY_H5AD}")
print(f"  Query filter: min_cells>={QUERY_MIN_CELLS}")
print(f"  Output merged:             {OUTPUT_H5AD}")
print(f"  Output dir:                {OUTPUT_DIR}")

# ==============================================================================
# HELPERS
# ==============================================================================

def _to_csr(X):
    if issparse(X):
        return csr_matrix(X)
    return csr_matrix(np.asarray(X))

def filter_genes_min_cells(X_csr: csr_matrix, genes: pd.Index, min_cells: int):
    nnz = np.asarray(X_csr.getnnz(axis=0)).ravel()
    keep = nnz >= int(min_cells)
    kept_genes = genes[keep]
    return keep, kept_genes, nnz

def align_full_to_reference_genes(X_qry_csr: csr_matrix, qry_genes: pd.Index, ref_genes: pd.Index):
    ref_index = pd.Index(ref_genes)
    ref_pos = ref_index.get_indexer(qry_genes)  # -1 if not in ref
    keep_mask = ref_pos >= 0
    n_shared = int(keep_mask.sum())
    if n_shared == 0:
        raise ValueError("No shared genes between query and reference full gene list!")

    X_shared = X_qry_csr[:, np.where(keep_mask)[0]]
    ref_pos_shared = ref_pos[keep_mask]

    order = np.argsort(ref_pos_shared)
    X_shared = X_shared[:, order]
    ref_pos_sorted = ref_pos_shared[order]

    X_coo = X_shared.tocoo()
    new_col = ref_pos_sorted[X_coo.col]
    X_aligned = csr_matrix(
        (X_coo.data, (X_coo.row, new_col)),
        shape=(X_qry_csr.shape[0], len(ref_genes))
    )

    stats = {
        "n_qry_genes": int(len(qry_genes)),
        "n_ref_genes": int(len(ref_genes)),
        "n_shared": n_shared,
        "shared_frac_vs_ref": n_shared / max(1, len(ref_genes)),
        "shared_frac_vs_qry": n_shared / max(1, len(qry_genes)),
    }
    return X_aligned, stats

def pick_first_existing(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None

def standardize_query_columns(adata_qry):
    """
    Ensure:
      - adata_qry.obs['Cell_Type_L2_final'] exists
      - adata_qry.obs['mapping_confidence'] exists
    If missing, copy from best-available candidate.
    """
    # L2 final
    l2_final_candidates = [
        "Cell_Type_L2_final",
        "Cell_Type_L2_pred",
        "Cell_Type_L2",
        "cell_type_L2",
        "predicted_labels",
        "predicted_label",
        "scanvi_pred",
        "scanvi_label",
        "label",
    ]
    src_l2 = pick_first_existing(adata_qry.obs.columns, l2_final_candidates)
    if src_l2 is None:
        l2_like = [c for c in adata_qry.obs.columns if ("L2" in c) or ("label" in c.lower()) or ("pred" in c.lower())]
        raise KeyError(f"Cannot find any query L2 label column. Candidates={l2_final_candidates}. "
                       f"L2/label/pred-like columns={l2_like[:50]}")
    if "Cell_Type_L2_final" not in adata_qry.obs.columns:
        adata_qry.obs["Cell_Type_L2_final"] = adata_qry.obs[src_l2]
    # confidence
    conf_candidates = [
        "mapping_confidence",
        "confidence",
        "pred_confidence",
        "prediction_confidence",
        "scanvi_confidence",
        "probability",
        "max_probability",
        "max_prob",
    ]
    src_conf = pick_first_existing(adata_qry.obs.columns, conf_candidates)
    if src_conf is None:
        # 不强制报错：有些 mapping 文件确实没 confidence
        src_conf = None
    if ("mapping_confidence" not in adata_qry.obs.columns) and (src_conf is not None):
        adata_qry.obs["mapping_confidence"] = adata_qry.obs[src_conf]

    return src_l2, src_conf

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================

print("\n" + "=" * 70)
print("PART 1: Load Data")
print("=" * 70)

print("\nStep 1: Load ORIGINAL reference (FULL genes)")
adata_ref_full = sc.read_h5ad(ORIGINAL_REF_H5AD)
print(f"  ref_full shape: {adata_ref_full.shape}")

print("\nStep 2: Load reference with UMAP (HVG)")
adata_ref_hvg = sc.read_h5ad(REF_WITH_UMAP_H5AD)
print(f"  ref_hvg shape:  {adata_ref_hvg.shape}")
print(f"  ref_hvg obsm:   {list(adata_ref_hvg.obsm.keys())}")

if "X_umap" not in adata_ref_hvg.obsm:
    print("ERROR: ref_hvg missing X_umap")
    sys.exit(1)

print("\nStep 3: Load query mapped (FULL in X per user)")
adata_qry = sc.read_h5ad(QUERY_H5AD)
print(f"  query shape:    {adata_qry.shape}")
print(f"  query obsm:     {list(adata_qry.obsm.keys())}")

if "X_umap" not in adata_qry.obsm:
    print("ERROR: query missing X_umap")
    sys.exit(1)

# ----------------------------
# FIX A1: Use mapped UMAP if available
# ----------------------------
if "X_umap_mapped" in adata_qry.obsm:
    adata_qry.obsm["X_umap"] = adata_qry.obsm["X_umap_mapped"].copy()
    print("  NOTE: query X_umap overwritten by X_umap_mapped (mapped coordinates)")

# ----------------------------
# FIX A2: Standardize query label/confidence columns (no extra plotting)
# ----------------------------
# L2 final: if missing, fallback to a reasonable column
if "Cell_Type_L2_final" not in adata_qry.obs.columns:
    fallback_l2 = None
    for c in ["Cell_Type_L2_pred", "predicted_labels", "Cell_Type_L2", "cell_type_L2", "label"]:
        if c in adata_qry.obs.columns:
            fallback_l2 = c
            break
    if fallback_l2 is None:
        raise KeyError("Query missing Cell_Type_L2_final and no fallback L2 column found.")
    adata_qry.obs["Cell_Type_L2_final"] = adata_qry.obs[fallback_l2]
    print(f"  NOTE: created Cell_Type_L2_final from fallback column: {fallback_l2}")

# mapping confidence: if missing, fallback (optional)
if "mapping_confidence" not in adata_qry.obs.columns:
    fallback_conf = None
    for c in ["confidence", "pred_confidence", "prediction_confidence", "max_probability", "probability"]:
        if c in adata_qry.obs.columns:
            fallback_conf = c
            break
    if fallback_conf is not None:
        adata_qry.obs["mapping_confidence"] = adata_qry.obs[fallback_conf]
        print(f"  NOTE: created mapping_confidence from fallback column: {fallback_conf}")
    else:
        print("  NOTE: mapping_confidence not found (and no fallback). Confidence panel will be N/A.")


# Uniqueness
for ad in (adata_ref_full, adata_ref_hvg, adata_qry):
    if not ad.var_names.is_unique:
        ad.var_names_make_unique()
    if not ad.obs_names.is_unique:
        ad.obs_names_make_unique()

# Standardize query obs columns BEFORE any subsetting
print("\nStep 3b: Standardize query label/confidence columns")
src_l2, src_conf = standardize_query_columns(adata_qry)
print(f"  Query L2 source -> Cell_Type_L2_final: {src_l2}")
print(f"  Query confidence source -> mapping_confidence: {src_conf if src_conf else 'None'}")

# quick debug stats (this is the fastest way to know whether upstream output is all-NA)
na_l2 = adata_qry.obs["Cell_Type_L2_final"].isna().mean()
print(f"  Query Cell_Type_L2_final NA fraction: {na_l2*100:.2f}%")
if "mapping_confidence" in adata_qry.obs.columns:
    na_conf = pd.to_numeric(adata_qry.obs["mapping_confidence"], errors="coerce").isna().mean()
    print(f"  Query mapping_confidence NA fraction: {na_conf*100:.2f}%")
else:
    print("  Query mapping_confidence: NOT AVAILABLE")

# ==============================================================================
# PART 2: Validate + Gene Lists
# ==============================================================================

print("\n" + "=" * 70)
print("PART 2: Validate + Gene Lists")
print("=" * 70)

ref_full_genes = pd.Index(adata_ref_full.var_names)
ref_full_var = adata_ref_full.var.copy()
ref_hvg_genes = pd.Index(adata_ref_hvg.var_names)

print(f"  Master FULL genes: {len(ref_full_genes):,}")
print(f"  Master HVG genes:  {len(ref_hvg_genes):,}")

missing_hvg_in_full = ref_hvg_genes.difference(ref_full_genes)
if len(missing_hvg_in_full) > 0:
    print(f"WARNING: {len(missing_hvg_in_full)} HVG genes not in ref_full. Dropping them for merge.")
    ref_hvg_genes = ref_hvg_genes.intersection(ref_full_genes)
    print(f"  Adjusted HVG genes: {len(ref_hvg_genes):,}")

if not adata_ref_full.obs_names.equals(adata_ref_hvg.obs_names):
    print("ERROR: Reference cell order mismatch between ORIGINAL_REF_H5AD and REF_WITH_UMAP_H5AD")
    sys.exit(1)
print(f"  OK Reference obs_names match ({adata_ref_hvg.n_obs:,} cells)")

for col in adata_ref_full.obs.columns:
    if col not in adata_ref_hvg.obs.columns:
        adata_ref_hvg.obs[col] = adata_ref_full.obs[col].values

# ==============================================================================
# PART 3: Build FULL matrices for .raw (aligned to ref_full_genes)
# ==============================================================================

print("\n" + "=" * 70)
print("PART 3: Build FULL matrices for .raw (aligned to ref_full_genes)")
print("=" * 70)

if "counts" in adata_ref_full.layers:
    ref_full_X = _to_csr(adata_ref_full.layers["counts"])
    ref_full_gene_axis = pd.Index(adata_ref_full.var_names)
    print("\nReference full matrix source:")
    print(f"  picked: layers['counts'], shape={ref_full_X.shape}, gene_axis={len(ref_full_gene_axis):,}")
else:
    ref_full_X = _to_csr(adata_ref_full.X)
    ref_full_gene_axis = pd.Index(adata_ref_full.var_names)
    print("\nReference full matrix source:")
    print(f"  picked: .X, shape={ref_full_X.shape}, gene_axis={len(ref_full_gene_axis):,}")

if not ref_full_gene_axis.equals(ref_full_genes):
    print("  Aligning reference full matrix to master ref_full_genes ...")
    ref_full_X_aligned, st = align_full_to_reference_genes(ref_full_X, ref_full_gene_axis, ref_full_genes)
    print(f"  ref align shared={st['n_shared']:,} ({st['shared_frac_vs_ref']*100:.1f}% of ref)")
else:
    ref_full_X_aligned = ref_full_X

# Query full matrix: FULL raw is in .X
qry_full_X = _to_csr(adata_qry.X)
qry_gene_axis = pd.Index(adata_qry.var_names)
print("\nQuery full matrix source:")
print(f"  picked: .X (user-defined raw), shape={qry_full_X.shape}, gene_axis={len(qry_gene_axis):,}")

print(f"\nFiltering query genes by min_cells>={QUERY_MIN_CELLS} ...")
keep_mask, qry_genes_kept, nnz = filter_genes_min_cells(qry_full_X, qry_gene_axis, QUERY_MIN_CELLS)
print(f"  kept genes: {qry_genes_kept.size:,} / {qry_gene_axis.size:,} "
      f"({qry_genes_kept.size/qry_gene_axis.size*100:.1f}%)")

qry_full_X_filt = qry_full_X[:, keep_mask]
qry_gene_axis_filt = qry_genes_kept

print("  Aligning filtered query full matrix to master ref_full_genes ...")
qry_full_X_aligned, stq = align_full_to_reference_genes(qry_full_X_filt, qry_gene_axis_filt, ref_full_genes)
print(f"  qry align shared={stq['n_shared']:,} "
      f"({stq['shared_frac_vs_ref']*100:.1f}% of ref; {stq['shared_frac_vs_qry']*100:.1f}% of qry_filtered)")

# ==============================================================================
# PART 4: Build HVG objects (.X) via overlap + inner join
# ==============================================================================

print("\n" + "=" * 70)
print("PART 4: Build HVG objects (.X) via overlap + inner join")
print("=" * 70)

qry_gene_set = set(qry_gene_axis_filt.tolist())
hvgs_for_merge = [g for g in ref_hvg_genes.tolist() if g in qry_gene_set]

if len(hvgs_for_merge) == 0:
    print("ERROR: HVG overlap after query filtering is zero.")
    sys.exit(1)

print(f"  HVGs for merge (ref-ordered overlap): {len(hvgs_for_merge):,} / ref_HVG {len(ref_hvg_genes):,}")

adata_ref = adata_ref_hvg[:, hvgs_for_merge].copy()
adata_qry_hvg = adata_qry[:, hvgs_for_merge].copy()

print(f"  ref HVG .X: {adata_ref.X.shape}")
print(f"  qry HVG .X: {adata_qry_hvg.X.shape}")

# keep only critical embeddings
critical_obsm = ["X_umap"]
if ("X_scANVI_L2" in adata_ref.obsm) and ("X_scANVI_L2" in adata_qry_hvg.obsm):
    critical_obsm.append("X_scANVI_L2")

for ad in (adata_ref, adata_qry_hvg):
    for k in list(ad.obsm.keys()):
        if k not in critical_obsm:
            del ad.obsm[k]
    for k in list(ad.obsp.keys()):
        del ad.obsp[k]
    ad.uns = {}

print(f"  critical obsm kept: {critical_obsm}")

# Ensure unique + add prefixes
if not adata_ref.obs_names.is_unique:
    adata_ref.obs_names_make_unique()
if not adata_qry_hvg.obs_names.is_unique:
    adata_qry_hvg.obs_names_make_unique()

adata_ref.obs_names = pd.Index([f"ref::{x}" for x in adata_ref.obs_names])
adata_qry_hvg.obs_names = pd.Index([f"qry::{x}" for x in adata_qry_hvg.obs_names])

# ==============================================================================
# PART 5: CONCATENATE (inner join)
# ==============================================================================

print("\n" + "=" * 70)
print("PART 5: Concatenate (inner join)")
print("=" * 70)

adata_all = sc.concat(
    {"reference": adata_ref, "query": adata_qry_hvg},
    axis=0,
    join="inner",
    merge="unique",
    label="data_source",
)

print(f"  merged .X shape: {adata_all.X.shape}")
print(f"  reference cells: {(adata_all.obs['data_source']=='reference').sum():,}")
print(f"  query cells:     {(adata_all.obs['data_source']=='query').sum():,}")

if "X_umap" not in adata_all.obsm:
    print("ERROR: X_umap lost during concat (unexpected).")
    sys.exit(1)

# ==============================================================================
# PART 6: ATTACH .raw (FULL genes, ref universe) + HVG counts layer
# ==============================================================================

print("\n" + "=" * 70)
print("PART 6: Attach .raw (FULL genes) + HVG counts layer")
print("=" * 70)

X_full_merged = sp_vstack([ref_full_X_aligned, qry_full_X_aligned], format="csr")
print(f"  full merged matrix for raw: {X_full_merged.shape}")

adata_all.raw = AnnData(
    X=X_full_merged,
    obs=adata_all.obs.copy(),
    var=ref_full_var.copy(),
)

assert adata_all.raw.n_obs == adata_all.n_obs
assert adata_all.raw.n_vars == len(ref_full_genes)
assert pd.Index(adata_all.raw.var_names).equals(ref_full_genes)
print("  OK .raw attached (FULL genes aligned to reference gene universe)")

hvg_pos = ref_full_genes.get_indexer(pd.Index(hvgs_for_merge))
if (hvg_pos < 0).any():
    raise RuntimeError("Internal error: some hvgs_for_merge not found in ref_full_genes")

counts_hvg = X_full_merged[:, hvg_pos]
adata_all.layers["counts"] = csr_matrix(counts_hvg)
assert adata_all.layers["counts"].shape == adata_all.X.shape
print("  OK layers['counts'] attached (HVG counts)")

# ==============================================================================
# PART 6b: Add unified L2 column (ref: Cell_Type_L2, qry: Cell_Type_L2_final)
# ==============================================================================

ref_m = (adata_all.obs["data_source"] == "reference")
qry_m = (adata_all.obs["data_source"] == "query")

adata_all.obs["Cell_Type_L2_unified"] = pd.NA  # object dtype first (avoid categorical-setitem error)

if "Cell_Type_L2" in adata_all.obs.columns:
    adata_all.obs.loc[ref_m, "Cell_Type_L2_unified"] = adata_all.obs.loc[ref_m, "Cell_Type_L2"].astype(str).values

if "Cell_Type_L2_final" in adata_all.obs.columns:
    adata_all.obs.loc[qry_m, "Cell_Type_L2_unified"] = adata_all.obs.loc[qry_m, "Cell_Type_L2_final"].astype(str).values

adata_all.obs["Cell_Type_L2_unified"] = adata_all.obs["Cell_Type_L2_unified"].astype("category")

# ==============================================================================
# PART 7: SAVE
# ==============================================================================

print("\n" + "=" * 70)
print("PART 7: Save")
print("=" * 70)

out_path = Path(OUTPUT_H5AD)
out_path.parent.mkdir(parents=True, exist_ok=True)

print(f"Saving: {out_path}")
adata_all.write_h5ad(out_path, compression="gzip")
file_size = out_path.stat().st_size / 1024**3
print(f"OK Saved ({file_size:.2f} GB)")

# ==============================================================================
# PART 8: QUICK VISUALIZATION (FIXED: query panels colored ONLY on query)
# ==============================================================================

print("\n" + "=" * 70)
print("PART 8: Visualization (overview)")
print("=" * 70)

ref_mask = (adata_all.obs["data_source"] == "reference")
qry_mask = (adata_all.obs["data_source"] == "query")

fig, axes = plt.subplots(2, 3, figsize=(20, 13))

# Panel 1
sc.pl.umap(
    adata_all, color="data_source", ax=axes[0, 0], show=False,
    title="Data Source", s=5,
    palette={"reference": "#1f77b4", "query": "#ff7f0e"},
)

# Panel 2: Reference L2 (ref-only)  -- dtype safe (no FutureWarning)
if "Cell_Type_L2" in adata_all.obs.columns:
    adata_all.obs["_ref_L2_only"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
    adata_all.obs.loc[ref_mask, "_ref_L2_only"] = adata_all.obs.loc[ref_mask, "Cell_Type_L2"].astype(str).values
    adata_all.obs["_ref_L2_only"] = adata_all.obs["_ref_L2_only"].astype("category")

    sc.pl.umap(
        adata_all, color="_ref_L2_only", ax=axes[0, 1], show=False,
        title="Reference L2 (Ref Only)", legend_loc="right margin", s=5
    )
    adata_all.obs.drop(columns=["_ref_L2_only"], inplace=True)
else:
    axes[0, 1].text(0.5, 0.5, "Cell_Type_L2\nN/A", ha="center", va="center",
                    transform=axes[0, 1].transAxes)
    axes[0, 1].set_title("Reference L2 (Ref Only)")

# Panel 3: Query L2 FINAL (query-only)  -- dtype safe
if "Cell_Type_L2_final" in adata_all.obs.columns:
    adata_all.obs["_qry_L2_final_only"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
    adata_all.obs.loc[qry_mask, "_qry_L2_final_only"] = adata_all.obs.loc[qry_mask, "Cell_Type_L2_final"].astype(str).values
    adata_all.obs["_qry_L2_final_only"] = adata_all.obs["_qry_L2_final_only"].astype("category")

    sc.pl.umap(
        adata_all, color="_qry_L2_final_only", ax=axes[0, 2], show=False,
        title="Query L2 (Final, Query Only)", legend_loc="right margin", s=5
    )
    adata_all.obs.drop(columns=["_qry_L2_final_only"], inplace=True)
else:
    axes[0, 2].text(0.5, 0.5, "Cell_Type_L2_final\nN/A", ha="center", va="center",
                    transform=axes[0, 2].transAxes)
    axes[0, 2].set_title("Query L2 (Final, Query Only)")

# Panel 4: Query confidence (query-only)  -- dtype safe
if "mapping_confidence" in adata_all.obs.columns:
    adata_all.obs["_qry_conf"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
    adata_all.obs.loc[qry_mask, "_qry_conf"] = pd.to_numeric(
        adata_all.obs.loc[qry_mask, "mapping_confidence"], errors="coerce"
    ).values

    sc.pl.umap(
        adata_all, color="_qry_conf", ax=axes[1, 0], show=False,
        title="Query Confidence (Query Only)", cmap="viridis", vmin=0, vmax=1, s=5
    )
    adata_all.obs.drop(columns=["_qry_conf"], inplace=True)
else:
    axes[1, 0].text(0.5, 0.5, "mapping_confidence\nN/A", ha="center", va="center",
                    transform=axes[1, 0].transAxes)
    axes[1, 0].set_title("Query Confidence (Query Only)")

# Panel 5: CD19 (prefer HVG else raw)
gene = "CD19"
use_raw = (gene not in adata_all.var_names) and (adata_all.raw is not None) and (gene in adata_all.raw.var_names)
if (gene in adata_all.var_names) or use_raw:
    sc.pl.umap(
        adata_all, color=gene, ax=axes[1, 1], show=False,
        title="CD19 (Pan-B)", cmap="Reds", s=5, use_raw=use_raw
    )
else:
    axes[1, 1].text(0.5, 0.5, "CD19\nN/A", ha="center", va="center",
                    transform=axes[1, 1].transAxes)
    axes[1, 1].set_title("CD19 (Pan-B)")

# Panel 6: Plasma marker
plasma_marker = None
for m in ["SDC1", "CD138"]:
    if (m in adata_all.var_names) or (adata_all.raw is not None and m in adata_all.raw.var_names):
        plasma_marker = m
        break

if plasma_marker is not None:
    use_raw_flag = (plasma_marker not in adata_all.var_names)
    sc.pl.umap(
        adata_all, color=plasma_marker, ax=axes[1, 2], show=False,
        title=f"{plasma_marker} (Plasma)", cmap="Reds", s=5, use_raw=use_raw_flag
    )
else:
    axes[1, 2].text(0.5, 0.5, "Plasma\nN/A", ha="center", va="center",
                    transform=axes[1, 2].transAxes)
    axes[1, 2].set_title("Plasma Marker")

plt.tight_layout()
fig_path = output_dir / "figures" / f"merged_umap_overview.{FIGURE_FORMAT}"
plt.savefig(fig_path, dpi=DPI, bbox_inches="tight")
plt.close()
print(f"  Saved: {fig_path}")

print("\n" + "=" * 70)
print("SUCCESS")
print("=" * 70)
print("\nOutput:")
print(f"  Merged h5ad: {out_path}")
print(f"  Figures:     {output_dir / 'figures'}")
print("\nData structure:")
print(f"  .X:                 {adata_all.X.shape} (HVG overlap, ref order)")
print(f"  .layers['counts']:  {adata_all.layers['counts'].shape} (HVG counts)")
print(f"  .raw.X:             {adata_all.raw.X.shape} (FULL genes, ref universe)")
print("  Added obs:          Cell_Type_L2_unified")
