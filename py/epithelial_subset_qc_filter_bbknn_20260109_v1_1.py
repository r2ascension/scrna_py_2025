#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# =============================================================================
# ===== epithelial_subset_qc_filter_bbknn_20260109_v1_1.py =====
#
# Purpose:
#   1) Collect all cells from ALL .h5ad files under SUBSET_DIR (obs-only),
#      extract: cell_id + (cell_type_level_2/3/4, Manual_Annotation if present),
#      add qc_pass_from_subset=1 for collected cells, and save CSV.
#      - Missing columns are allowed (filled with NA).
#      - Duplicated cell_id across subset files are aggregated by first non-null.
#      - Conflicts (same cell_id has >1 distinct non-null value) are detected & logged.
#
#   2) Filter INPUT_H5AD to KEEP ONLY cells whose cell_id appears in (1).
#      Save kept/dropped cell_id lists for reproducibility.
#
#   3) Count cells by Manual_Annotation after filtering and save CSV.
#
#   4) Re-run BBKNN on filtered AnnData using HVG-only preprocessing to avoid
#      scaling/centering all genes; transfer neighbors/umap/leiden/pca back to
#      full-gene object; save final H5AD.
#
# Dependencies:
#   scanpy, anndata, pandas, numpy, bbknn, scipy
#
# Date: 2026-01-09
# =============================================================================

import os
import sys
import glob
import time
import logging
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
from scipy import sparse


# =========================
# User Config (EDIT HERE)
# =========================
SUBSET_DIR = "/home/h2048/data/py/0109/subset"
INPUT_H5AD = "/home/h2048/data/core2/adata_epithelial_with_manual_annotations_raw.h5ad"

OUTDIR = "/home/h2048/data/py/0109/subset_qc_filter_bbknn_20260109_v1_1"
os.makedirs(OUTDIR, exist_ok=True)

# Outputs
OUT_CSV_ALL_SUBSET_CELLS = os.path.join(OUTDIR, "subset_cells_qc_annotations.csv")
OUT_CSV_SUBSET_CONFLICTS = os.path.join(OUTDIR, "subset_cellid_annotation_conflicts.csv")

OUT_TXT_KEPT_CELL_IDS = os.path.join(OUTDIR, "kept_cell_ids.txt")
OUT_TXT_DROPPED_CELL_IDS = os.path.join(OUTDIR, "dropped_cell_ids.txt")

OUT_H5AD_FILTERED_RAW = os.path.join(OUTDIR, "adata_epithelial_filtered_by_subset_raw.h5ad")
OUT_CSV_MANUAL_STATS = os.path.join(OUTDIR, "manual_annotation_counts_after_filter.csv")
OUT_H5AD_BBKNN = os.path.join(OUTDIR, "adata_epithelial_filtered_bbknn.h5ad")

# Columns to extract from subset h5ad obs (optional)
WANTED_COLS = ["cell_type_level_2", "cell_type_level_3", "cell_type_level_4", "Manual_Annotation"]
QC_FLAG_COL = "qc_pass_from_subset"

# BBKNN / embedding parameters
BATCH_KEY = "dataset"  # must exist in adata.obs; edit if needed
N_HVG = 3000
N_PCS = 30
BBKNN_NEIGHBORS_WITHIN_BATCH = 5
BBKNN_TRIM = None  # e.g., 10 or None
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.5
LEIDEN_RESOLUTION = 1.0
RANDOM_SEED = 0

# Preprocess source layer preference (use the first existing layer)
COUNTS_LAYER_CANDIDATES = ["counts", "raw_counts"]


# =========================
# Logging
# =========================
LOG_PATH = os.path.join(OUTDIR, "run.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[logging.FileHandler(LOG_PATH), logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("subset_qc_filter_bbknn")


# =========================
# Helpers
# =========================
def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def list_h5ad_files(folder: str) -> List[str]:
    return sorted(glob.glob(os.path.join(folder, "*.h5ad")))


def pick_counts_layer(a: ad.AnnData) -> Optional[str]:
    for k in COUNTS_LAYER_CANDIDATES:
        if k in a.layers:
            return k
    return None


def detect_conflicts_for_duplicates(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """
    For duplicated cell_id, detect if any column has >1 distinct non-null values.
    Return a long-form table: cell_id, column, n_distinct_nonnull, examples
    """
    if df.empty or "cell_id" not in df.columns:
        return pd.DataFrame(columns=["cell_id", "column", "n_distinct_nonnull", "examples"])

    dup_ids = df["cell_id"][df["cell_id"].duplicated(keep=False)]
    if dup_ids.empty:
        return pd.DataFrame(columns=["cell_id", "column", "n_distinct_nonnull", "examples"])

    df_dup = df[df["cell_id"].isin(dup_ids)].copy()

    rows = []
    for c in cols:
        if c not in df_dup.columns:
            continue

        g = df_dup.groupby("cell_id", sort=False)[c]
        for cell_id, s in g:
            s2 = s.dropna().astype("object")
            if s2.empty:
                continue
            uniq = pd.unique(s2)
            uniq = [x for x in uniq if pd.notna(x)]
            if len(uniq) > 1:
                ex = "; ".join([str(x) for x in uniq[:5]])
                rows.append(
                    {
                        "cell_id": str(cell_id),
                        "column": c,
                        "n_distinct_nonnull": int(len(uniq)),
                        "examples": ex,
                    }
                )

    if not rows:
        return pd.DataFrame(columns=["cell_id", "column", "n_distinct_nonnull", "examples"])
    return pd.DataFrame(rows)


def aggregate_by_cell_id_first_nonnull(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate duplicated cell_id rows by taking first non-null per column.
    """
    if df.empty:
        return df

    def first_nonnull(s: pd.Series):
        s2 = s.dropna()
        if len(s2) == 0:
            return np.nan
        return s2.iloc[0]

    return df.groupby("cell_id", sort=False).agg(first_nonnull).reset_index()


def collect_subset_cells(subset_dir: str, wanted_cols: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read obs from each subset .h5ad (backed), output:
      - subset_df: aggregated unique cell_id table with qc flag and optional columns
      - conflicts_df: conflicts detected prior to aggregation
    """
    files = list_h5ad_files(subset_dir)
    if not files:
        raise FileNotFoundError(f"No .h5ad found under: {subset_dir}")

    logger.info(f"Found {len(files)} subset h5ad files.")
    all_rows = []

    for fp in files:
        logger.info(f"Reading subset obs (backed): {fp}")
        a = ad.read_h5ad(fp, backed="r")
        obs = a.obs

        out = pd.DataFrame({"cell_id": obs.index.astype(str)})
        out[QC_FLAG_COL] = 1

        for c in wanted_cols:
            if c in obs.columns:
                out[c] = obs[c].astype("object").values
            else:
                out[c] = np.nan

        all_rows.append(out)

        try:
            a.file.close()
        except Exception:
            pass

    df_raw = pd.concat(all_rows, axis=0, ignore_index=True)
    logger.info(f"Subset rows collected (pre-aggregation): n={df_raw.shape[0]}")

    conflicts_df = detect_conflicts_for_duplicates(df_raw, cols=[QC_FLAG_COL] + wanted_cols)
    if conflicts_df.shape[0] > 0:
        logger.warning(
            f"Detected conflicts across subset files: {conflicts_df.shape[0]} (cell_id, column) pairs. "
            f"Will aggregate by first non-null; see: {OUT_CSV_SUBSET_CONFLICTS}"
        )
    else:
        logger.info("No annotation conflicts detected across duplicated cell_id in subset files.")

    df_agg = aggregate_by_cell_id_first_nonnull(df_raw)
    if df_agg.shape[0] != df_raw.shape[0]:
        logger.info(f"Aggregated duplicated cell_id: {df_raw.shape[0]} -> {df_agg.shape[0]} unique cell_id")

    return df_agg, conflicts_df


def ensure_obs_columns_from_subset(a: ad.AnnData, subset_df: pd.DataFrame, cols_to_add: List[str]) -> None:
    """
    Left-join subset_df onto a.obs by cell_id (obs_names).
    Ensures safe assignment even if existing obs columns are categorical.
    """
    subset_df2 = subset_df.set_index("cell_id", drop=True)

    # qc flag: all cells in a passed subset filter
    a.obs[QC_FLAG_COL] = 1

    for c in cols_to_add:
        if c not in subset_df2.columns:
            continue
        incoming = subset_df2.reindex(a.obs_names)[c].astype("object")

        if c not in a.obs.columns:
            a.obs[c] = incoming.values
        else:
            # cast to object to avoid categorical new-category assignment crash
            a.obs[c] = a.obs[c].astype("object")
            mask_na = pd.isna(a.obs[c].values)
            if mask_na.any():
                a.obs.loc[mask_na, c] = incoming.values[mask_na]


def manual_annotation_stats(a: ad.AnnData, col: str = "Manual_Annotation") -> pd.DataFrame:
    if col not in a.obs.columns:
        return pd.DataFrame({"Manual_Annotation": ["<MISSING_COLUMN>"], "n_cells": [a.n_obs]})
    vc = a.obs[col].astype("object").fillna("<NA>").value_counts(dropna=False)
    return vc.rename_axis(col).reset_index(name="n_cells")


def _safe_pp_pca(a: ad.AnnData, n_comps: int, zero_center: bool) -> None:
    """
    Prefer sc.pp.pca (supports zero_center) for sparse safety.
    Fallback to sc.tl.pca if signature differs.
    """
    try:
        sc.pp.pca(
            a,
            n_comps=n_comps,
            svd_solver="arpack",
            zero_center=zero_center,
            random_state=RANDOM_SEED,
        )
    except TypeError:
        # older/newer API variations
        sc.tl.pca(a, n_comps=n_comps, svd_solver="arpack", random_state=RANDOM_SEED)


def run_bbknn_pipeline_hvg_only(a_full: ad.AnnData) -> ad.AnnData:
    """
    HVG-only preprocessing to avoid scaling all genes.
    Transfers neighbors / embeddings / clustering back to full-gene object.
    """
    if BATCH_KEY not in a_full.obs.columns:
        candidates = [x for x in ["study", "batch", "dataset", "donor", "sample", "orig.ident"] if x in a_full.obs.columns]
        raise KeyError(
            f"BATCH_KEY='{BATCH_KEY}' not found in adata.obs. "
            f"Available candidates: {candidates}. Please edit BATCH_KEY."
        )

    try:
        import bbknn  # noqa: F401
    except ImportError as e:
        raise ImportError("bbknn is not installed in this environment. Please install bbknn and retry.") from e

    counts_layer = pick_counts_layer(a_full)
    if counts_layer is not None:
        logger.info(f"Using layer='{counts_layer}' as raw counts source for HVG selection/preprocess.")
    else:
        logger.info("No counts layer found; using adata.X as source. Ensure X is counts if you rely on seurat_v3 HVG.")

    # HVG selection on full (no mutation of X)
    logger.info(f"HVG selection: n_top_genes={N_HVG}")
    try:
        sc.pp.highly_variable_genes(
            a_full,
            n_top_genes=N_HVG,
            flavor="seurat_v3",
            layer=counts_layer,
            subset=False,
        )
    except Exception as e:
        logger.warning(f"HVG(seurat_v3) failed: {repr(e)}. Falling back to flavor='cell_ranger'.")
        sc.pp.highly_variable_genes(
            a_full,
            n_top_genes=N_HVG,
            flavor="cell_ranger",
            layer=counts_layer,
            subset=False,
        )

    if "highly_variable" not in a_full.var.columns or int(a_full.var["highly_variable"].sum()) == 0:
        raise RuntimeError("HVG selection resulted in 0 HVGs. Check counts layer / gene filtering.")

    hv_mask = a_full.var["highly_variable"].values
    logger.info(f"HVGs selected: {int(hv_mask.sum())}")

    # Work object: HVG subset
    a = a_full[:, hv_mask].copy()

    # Use counts for preprocessing on the copy only
    if counts_layer is not None and counts_layer in a.layers:
        a.X = a.layers[counts_layer]

    # normalize/log
    logger.info("normalize_total + log1p (HVG-only)")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # scale: keep sparse if possible (avoid zero_center=True on sparse)
    is_sparse = sparse.issparse(a.X)
    if is_sparse:
        zero_center = False
    else:
        # dense path: zero_center only if not too large
        n_entries = a.n_obs * a.n_vars
        zero_center = True if n_entries <= 150_000_000 else False

    # cast to float32 to reduce peak memory for dense
    try:
        if not sparse.issparse(a.X):
            a.X = a.X.astype(np.float32, copy=False)
    except Exception:
        pass

    logger.info(f"scale (HVG-only): zero_center={zero_center}, max_value=10")
    sc.pp.scale(a, max_value=10, zero_center=zero_center)

    # PCA: use zero_center consistent with sparse safety
    # If zero_center=False, scanpy will typically use TruncatedSVD path (sparse-safe).
    logger.info(f"PCA (HVG-only): n_comps={N_PCS}, zero_center={zero_center}")
    _safe_pp_pca(a, n_comps=N_PCS, zero_center=zero_center)

    # BBKNN
    import bbknn

    logger.info(
        f"BBKNN: batch_key={BATCH_KEY}, neighbors_within_batch={BBKNN_NEIGHBORS_WITHIN_BATCH}, "
        f"n_pcs={N_PCS}, trim={BBKNN_TRIM}"
    )
    bbknn.bbknn(
        a,
        batch_key=BATCH_KEY,
        neighbors_within_batch=BBKNN_NEIGHBORS_WITHIN_BATCH,
        n_pcs=N_PCS,
        trim=BBKNN_TRIM,
    )

    # UMAP + Leiden
    logger.info("UMAP")
    sc.tl.umap(a, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=RANDOM_SEED)

    logger.info(f"Leiden: resolution={LEIDEN_RESOLUTION}")
    sc.tl.leiden(a, resolution=LEIDEN_RESOLUTION, random_state=RANDOM_SEED, key_added="leiden_bbknn")

    # Transfer results back to full object
    # neighbors graph
    for k in ["distances", "connectivities"]:
        if k in a.obsp:
            a_full.obsp[k] = a.obsp[k]
    if "neighbors" in a.uns:
        a_full.uns["neighbors"] = a.uns["neighbors"]

    # embeddings
    if "X_pca" in a.obsm:
        a_full.obsm["X_pca"] = a.obsm["X_pca"]
    if "X_umap" in a.obsm:
        a_full.obsm["X_umap"] = a.obsm["X_umap"]

    # clustering
    a_full.obs["leiden_bbknn"] = a.obs["leiden_bbknn"].astype("object").values

    return a_full


def write_id_list(path: str, ids: np.ndarray) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for x in ids:
            f.write(str(x) + "\n")


# =========================
# Main
# =========================
def main() -> None:
    t0 = time.time()
    logger.info("=== START v1.1 ===")
    logger.info(f"Time: {_now()}")
    logger.info(f"SUBSET_DIR: {SUBSET_DIR}")
    logger.info(f"INPUT_H5AD: {INPUT_H5AD}")
    logger.info(f"OUTDIR: {OUTDIR}")
    logger.info(f"scanpy={getattr(sc, '__version__', 'NA')} | anndata={getattr(ad, '__version__', 'NA')} | python={sys.version.split()[0]}")

    # 1) Collect subset cells
    subset_df, conflicts_df = collect_subset_cells(SUBSET_DIR, wanted_cols=WANTED_COLS)

    subset_df.to_csv(OUT_CSV_ALL_SUBSET_CELLS, index=False, encoding="utf-8")
    logger.info(f"Saved subset cell table: {OUT_CSV_ALL_SUBSET_CELLS} (n={subset_df.shape[0]})")

    if conflicts_df.shape[0] > 0:
        conflicts_df.to_csv(OUT_CSV_SUBSET_CONFLICTS, index=False, encoding="utf-8")
        logger.info(f"Saved conflicts table: {OUT_CSV_SUBSET_CONFLICTS}")

        # brief summary per column
        col_summary = conflicts_df.groupby("column", sort=False)["cell_id"].nunique().reset_index(name="n_cell_ids_with_conflict")
        logger.warning("Conflict summary (unique cell_id per column):\n" + col_summary.to_string(index=False))
    else:
        # still write an empty file for pipeline stability
        conflicts_df.to_csv(OUT_CSV_SUBSET_CONFLICTS, index=False, encoding="utf-8")
        logger.info(f"Saved empty conflicts table: {OUT_CSV_SUBSET_CONFLICTS}")

    subset_ids = subset_df["cell_id"].astype(str).values

    # 2) Filter INPUT_H5AD by subset cell IDs
    logger.info("Reading INPUT_H5AD in backed mode for filtering.")
    a_backed = ad.read_h5ad(INPUT_H5AD, backed="r")
    obs_names = a_backed.obs_names.astype(str)

    # vectorized keep mask
    keep_mask = obs_names.isin(subset_ids)
    keep_mask = np.asarray(keep_mask, dtype=bool)

    n_total = int(len(obs_names))
    n_keep = int(keep_mask.sum())
    n_drop = n_total - n_keep

    logger.info(f"Filter by subset IDs: keep {n_keep}/{n_total} cells; drop {n_drop}")

    if n_keep == 0:
        # close handle
        try:
            a_backed.file.close()
        except Exception:
            pass
        raise RuntimeError("After filtering, 0 cells remain. Check cell_id matching between subset and INPUT_H5AD.")

    kept_ids = obs_names[keep_mask].values
    dropped_ids = obs_names[~keep_mask].values

    write_id_list(OUT_TXT_KEPT_CELL_IDS, kept_ids)
    write_id_list(OUT_TXT_DROPPED_CELL_IDS, dropped_ids)
    logger.info(f"Saved kept cell IDs: {OUT_TXT_KEPT_CELL_IDS}")
    logger.info(f"Saved dropped cell IDs: {OUT_TXT_DROPPED_CELL_IDS}")

    # slice then materialize in memory
    a_filt = a_backed[keep_mask].to_memory()

    # close backed
    try:
        a_backed.file.close()
    except Exception:
        pass

    # 2b) Attach subset annotations + qc flag (safe for categorical)
    ensure_obs_columns_from_subset(a_filt, subset_df=subset_df, cols_to_add=WANTED_COLS)

    # Save filtered raw
    logger.info(f"Saving filtered raw h5ad: {OUT_H5AD_FILTERED_RAW}")
    a_filt.write_h5ad(OUT_H5AD_FILTERED_RAW)

    # 3) Manual_Annotation stats after filtering
    stats = manual_annotation_stats(a_filt, col="Manual_Annotation")
    stats.to_csv(OUT_CSV_MANUAL_STATS, index=False, encoding="utf-8")
    logger.info("Manual_Annotation counts after filtering:\n" + stats.to_string(index=False))
    logger.info(f"Saved stats CSV: {OUT_CSV_MANUAL_STATS}")

    # 4) Re-run BBKNN (HVG-only preprocess)
    logger.info("Running BBKNN pipeline (HVG-only preprocessing).")
    a_bbknn = run_bbknn_pipeline_hvg_only(a_filt)

    logger.info(f"Saving BBKNN h5ad: {OUT_H5AD_BBKNN}")
    a_bbknn.write_h5ad(OUT_H5AD_BBKNN)

    dt = time.time() - t0
    logger.info(f"=== DONE in {dt:.1f}s ===")
    logger.info("Outputs:")
    logger.info(f"- {OUT_CSV_ALL_SUBSET_CELLS}")
    logger.info(f"- {OUT_CSV_SUBSET_CONFLICTS}")
    logger.info(f"- {OUT_TXT_KEPT_CELL_IDS}")
    logger.info(f"- {OUT_TXT_DROPPED_CELL_IDS}")
    logger.info(f"- {OUT_H5AD_FILTERED_RAW}")
    logger.info(f"- {OUT_CSV_MANUAL_STATS}")
    logger.info(f"- {OUT_H5AD_BBKNN}")


if __name__ == "__main__":
    main()
