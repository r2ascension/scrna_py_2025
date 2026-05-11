#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
B-cell reference scVI/scANVI full rerun after removing CHOIR clusters 22, 13, and 25.

Goals
-----
1. Start from the original full-gene B-cell reference h5ad used by tissue comparison.
2. Remove the previously excluded CHOIR cluster 22 cells plus the newly requested
   CHOIR clusters 13 and 25 cells.
3. Guarantee both `layers['counts']` and `raw.X` exist before training.
4. Re-run reference scVI/scANVI fully fresh with L3 supervision (`cell_type_expert`).
5. Save a fresh reference h5ad for downstream tissue comparison and overwrite the
   existing 0415 Python output directory in place.

Notes
-----
- Cluster 22 is loaded from the validated 0412 CHOIR assignment table so the prior
  contamination removal is preserved.
- Clusters 13 and 25 are loaded from the current 0415 CHOIR assignment table.
- This script overwrites `cell_type_scanvi_pred` with the fresh rerun L3 predictions
  so the downstream R tissue-comparison pipeline uses the rerun labels instead of
  any stale input column.
- `raw.X` is created from the full-gene counts matrix before HVG subsetting, so the
  saved output preserves full-gene raw counts while the working matrix keeps the HVG
  training view.
"""

from __future__ import annotations

import gc
import json
import os
import shutil
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import umap
from pandas.api.types import is_object_dtype, is_string_dtype
from scipy import sparse
from scipy.sparse import csr_matrix, issparse

warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None

INPUT_REF_H5AD = Path(
    "/home/h2048/data/py/0203/bcell_scarches_v4_1/results/scarches_package/bcell_reference_20260203.h5ad"
)
CHOIR_CLUSTER_SOURCES = (
    {
        "csv": Path(
            "/home/h2048/data/R/0412/bcell_tissue_comparison_v2_6_4_20260412/reports/choir/choir_clusters.csv"
        ),
        "clusters": (22,),
        "source_label": "0412_validated_c22",
    },
    {
        "csv": Path(
            "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/reports/choir/choir_clusters.csv"
        ),
        "clusters": (13, 25),
        "source_label": "0415_current_c13_c25",
    },
)
OUTPUT_DIR = Path("/home/h2048/data/py/0415/bcell_scvi_scanvi_ref_c22drop_l3_20260415")
OUTPUT_H5AD = OUTPUT_DIR / "bcell_reference_c22_c13_c25drop_scanvi_L3_ref_20260417.h5ad"
OUTPUT_SCVI_DIR = OUTPUT_DIR / "scvi_bcell_ref_c22_c13_c25drop_l3_20260417"
OUTPUT_SCANVI_DIR = OUTPUT_DIR / "scanvi_bcell_L3_ref_c22_c13_c25drop_20260417"
OUTPUT_UMAP = OUTPUT_DIR / "umap_operator_scanvi.joblib"
OUTPUT_REMOVED = OUTPUT_DIR / "removed_choir_clusters_c22_c13_c25_cells.csv"
OUTPUT_TARGETS = OUTPUT_DIR / "target_choir_cluster_union_c22_c13_c25.csv"
OUTPUT_CONFIG = OUTPUT_DIR / "training_config.json"
OUTPUT_HVG = OUTPUT_DIR / "hvg_genes.txt"

BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
L3_KEY = "cell_type_expert"
L2_KEY = "Cell_Type_L2"
SCANVI_LABELS_KEY = "scanvi_label_l3_train"
L3_PRED_KEY = "L3_scanvi_pred"
L3_CONFIDENCE_KEY = "L3_scanvi_confidence"
L3_PROBA_KEY = "proba_L3_scanvi"
LEGACY_L3_OUTPUT_KEY = "cell_type_scanvi_pred"
LEGACY_L3_INPUT_BACKUP_KEY = "cell_type_scanvi_pred_input"
LEGACY_L2_INPUT_BACKUP_KEY = "Cell_Type_L2_input"
L3_INPUT_BACKUP_KEY = "cell_type_expert_input"
UNLABELED_CATEGORY = "Unknown"
COUNTS_LAYER = "counts"
NOSE_TISSUE = "nose"
NON_NOSE_GCB_MEMORY_TARGET = "Memory_B"
GC_B_L3_LABELS = {
    "GC_B_Dark_Zone_Centroblast_Cycling",
    "GC_B_Light_Zone_Centrocyte",
    "GC_B_Transitional",
    "GC_B",
}

L3_TO_L2_MAP = {
    "Atypical_Memory_B": "Atypical_Memory_B",
    "IGHEplus_Atypical_Memory_B": "Atypical_Memory_B",
    "GC_B_Light_Zone_Centrocyte": "GC_B",
    "GC_B_Transitional": "GC_B",
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B": "GC_B",
    "Memory_B": "Memory_B",
    "Naive_B": "Naive_B",
    "Plasma_IgA": "Plasma",
    "Plasma_IgG": "Plasma",
    "Plasma": "Plasma",
    "Unknown": "Unknown",
}

STRESS_SIGNATURE_GENES = [
    "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "IER2",
]

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

FORCED_MARKERS_BCELL = [
    "CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "BLK", "FCRL2", "FCRL3", "FCRL5",
    "TCL1A", "FCER2", "IGHD", "IL4R",
    "CD27", "TNFRSF13B", "AIM2",
    "AICDA", "BCL6", "MME", "LMO2", "STMN1", "CXCR4", "RGS13",
    "SDC1", "CD38", "XBP1", "PRDM1", "JCHAIN", "IGHA1", "IGHA2", "IGHG1", "IGHG2",
    "IGHG3", "IGHG4", "MZB1", "DERL3", "SSR4", "ITGAX", "TBX21", "CXCR3",
]

N_HVG = 4000
SCVI_N_LATENT = 50
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_DROPOUT_RATE = 0.1
SCVI_GENE_LIKELIHOOD = "nb"
MAX_EPOCHS_SCVI = 400
MAX_EPOCHS_SCANVI = 60
SCVI_BATCH_SIZE = 256
SCANVI_BATCH_SIZE = 2048
SCANVI_EARLY_STOPPING_PATIENCE = 10
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 0.0
UMAP_N_NEIGHBORS = 30
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
RANDOM_SEED = 42
TORCH_NUM_THREADS = min(56, os.cpu_count() or 8)
TORCH_NUM_INTEROP_THREADS = 4

np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED
torch.set_num_threads(TORCH_NUM_THREADS)
try:
    torch.set_num_interop_threads(TORCH_NUM_INTEROP_THREADS)
except RuntimeError:
    pass


def resolve_training_device() -> tuple[str, int | None, str | None]:
    if not torch.cuda.is_available():
        return "cpu", None, None
    try:
        device_name = torch.cuda.get_device_name(0)
        _ = torch.cuda.get_device_properties(0)
        return "gpu", 1, device_name
    except Exception as exc:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        try:
            torch.cuda.is_available = lambda: False  # type: ignore[assignment]
            torch.cuda.device_count = lambda: 0  # type: ignore[assignment]
        except Exception:
            pass
        print(f"[WARN] CUDA reported as available but is unusable: {exc}")
        print("[WARN] Falling back to CPU training")
        return "cpu", None, None


def print_banner(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def reset_output_dir(output_dir: Path) -> None:
    if output_dir.exists():
        for child in output_dir.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    output_dir.mkdir(parents=True, exist_ok=True)


def load_cluster_rows(cluster_csv: Path, cluster_ids: tuple[int, ...], source_label: str) -> pd.DataFrame:
    df = pd.read_csv(cluster_csv)
    if "choir_cluster" not in df.columns or "cell" not in df.columns:
        raise KeyError(f"{cluster_csv} must contain columns 'cell' and 'choir_cluster'")

    numeric_cluster = pd.to_numeric(df["choir_cluster"], errors="coerce")
    hits = df.loc[numeric_cluster.isin(cluster_ids), ["cell", "choir_cluster"]].copy()
    if hits.empty:
        raise ValueError(f"No cells found for CHOIR clusters {cluster_ids} in {cluster_csv}")

    hits["cell"] = hits["cell"].astype(str)
    hits["choir_cluster"] = pd.to_numeric(hits["choir_cluster"], errors="raise").astype(int)
    hits["cluster_csv"] = str(cluster_csv)
    hits["cluster_source_label"] = source_label
    return hits


def load_union_cluster_rows() -> pd.DataFrame:
    rows = []
    for spec in CHOIR_CLUSTER_SOURCES:
        rows.append(
            load_cluster_rows(
                cluster_csv=spec["csv"],
                cluster_ids=tuple(int(x) for x in spec["clusters"]),
                source_label=str(spec["source_label"]),
            )
        )
    union_df = pd.concat(rows, ignore_index=True)
    union_df = union_df.drop_duplicates(subset=["cell", "choir_cluster", "cluster_csv"])
    return union_df.sort_values(["choir_cluster", "cluster_source_label", "cell"]).reset_index(drop=True)


def describe_cluster_targets(union_df: pd.DataFrame) -> None:
    print("[INFO] CHOIR target cells by source/cluster:")
    summary = (
        union_df.groupby(["cluster_source_label", "choir_cluster"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["cluster_source_label", "choir_cluster"])
    )
    for _, row in summary.iterrows():
        print(
            f"  - {row['cluster_source_label']} | cluster {int(row['choir_cluster'])}: "
            f"{int(row['n_cells']):,} cells"
        )
    print(f"[INFO] Union target cells: {union_df['cell'].nunique():,}")


def ensure_counts_layer(adata: sc.AnnData, counts_layer: str = COUNTS_LAYER) -> None:
    if counts_layer in adata.layers:
        adata.X = adata.layers[counts_layer]
        if issparse(adata.X) and not isinstance(adata.X, csr_matrix):
            adata.X = adata.X.tocsr()
        if issparse(adata.layers[counts_layer]) and not isinstance(adata.layers[counts_layer], csr_matrix):
            adata.layers[counts_layer] = adata.layers[counts_layer].tocsr()
        return

    sample = adata.X.data[:1000] if issparse(adata.X) else np.asarray(adata.X).ravel()[:1000]
    if np.all(np.mod(sample, 1) == 0):
        print(f"[WARN] layers['{counts_layer}'] missing; creating it from integer-like X")
        adata.layers[counts_layer] = adata.X.copy()
        if issparse(adata.layers[counts_layer]) and not isinstance(adata.layers[counts_layer], csr_matrix):
            adata.layers[counts_layer] = adata.layers[counts_layer].tocsr()
        adata.X = adata.layers[counts_layer]
        return

    raise ValueError(
        "Missing raw counts layer and X does not look integer-like. "
        "Cannot continue with scVI/scANVI training."
    )


def ensure_raw_from_counts(adata: sc.AnnData, counts_layer: str = COUNTS_LAYER) -> None:
    if adata.raw is not None:
        print(f"[OK] raw already present: {adata.raw.n_vars:,} genes")
        return

    print("[INFO] raw missing; creating adata.raw from full-gene counts layer")
    raw_counts = adata.layers[counts_layer]
    if issparse(raw_counts) and not isinstance(raw_counts, csr_matrix):
        raw_counts = raw_counts.tocsr()
    raw_adata = sc.AnnData(
        X=raw_counts.copy() if hasattr(raw_counts, "copy") else raw_counts,
        obs=adata.obs.copy(),
        var=adata.var.copy(),
    )
    raw_adata.var_names = adata.var_names.copy()
    adata.raw = raw_adata
    gc.collect()
    print(f"[OK] raw created: {adata.raw.n_vars:,} genes")


def remove_choir_cluster_cells(
    adata: sc.AnnData,
    union_df: pd.DataFrame,
) -> tuple[sc.AnnData, pd.DataFrame]:
    target_cells = pd.Index(union_df["cell"].astype(str).unique())
    obs_names = pd.Index(adata.obs_names.astype(str))
    barcode_idx = pd.Index(adata.obs["barcode"].astype(str)) if "barcode" in adata.obs.columns else pd.Index([])

    hit_obs = obs_names.isin(target_cells)
    hit_barcode = barcode_idx.isin(target_cells) if len(barcode_idx) == adata.n_obs else np.zeros(adata.n_obs, dtype=bool)
    drop_mask = hit_obs | hit_barcode

    n_drop = int(drop_mask.sum())
    if n_drop == 0:
        raise ValueError("Requested CHOIR-cluster cells were not found in the reference h5ad")

    cluster_lookup = (
        union_df.groupby("cell", observed=False)
        .agg(
            drop_choir_clusters=("choir_cluster", lambda x: ",".join(str(int(v)) for v in sorted(set(x)))),
            drop_cluster_sources=("cluster_source_label", lambda x: ",".join(sorted(set(str(v) for v in x)))),
        )
    )

    removed = adata.obs.loc[drop_mask].copy()
    matched_target_cell = np.where(hit_obs[drop_mask], obs_names[drop_mask], barcode_idx[drop_mask] if len(barcode_idx) == adata.n_obs else obs_names[drop_mask])
    matched_target_cell = pd.Index(matched_target_cell.astype(str))
    removed.insert(0, "obs_name", adata.obs_names[drop_mask].astype(str))
    removed.insert(1, "matched_target_cell", matched_target_cell)
    removed.insert(2, "matched_by_obs_name", hit_obs[drop_mask].astype(bool))
    removed.insert(3, "matched_by_barcode", hit_barcode[drop_mask].astype(bool))
    removed.insert(4, "drop_choir_clusters", matched_target_cell.map(cluster_lookup["drop_choir_clusters"]).to_numpy())
    removed.insert(5, "drop_cluster_sources", matched_target_cell.map(cluster_lookup["drop_cluster_sources"]).to_numpy())

    print(f"[OK] removing {n_drop:,} cells from CHOIR targets")
    print("[INFO] removed counts by cluster label:")
    print(removed["drop_choir_clusters"].astype(str).value_counts().to_string())
    if "tissue" in removed.columns:
        print("[INFO] removed tissue distribution:")
        print(removed["tissue"].astype(str).value_counts().to_string())

    keep_mask = ~drop_mask
    adata = adata[keep_mask].copy()
    adata.obs_names_make_unique()
    return adata, removed


def map_l3_to_l2(adata: sc.AnnData) -> None:
    if L3_KEY not in adata.obs.columns:
        raise KeyError(f"Missing L3 source column: {L3_KEY}")

    l3_series = adata.obs[L3_KEY].astype(str)
    adata.obs[L2_KEY] = l3_series.map(L3_TO_L2_MAP)
    n_missing = int(adata.obs[L2_KEY].isna().sum())
    if n_missing > 0:
        unmapped = sorted(pd.unique(l3_series[adata.obs[L2_KEY].isna()]).tolist())
        raise ValueError(f"Unmapped L3 labels after CHOIR-cluster filtering: {unmapped}")

    adata.obs[L2_KEY] = adata.obs[L2_KEY].astype("category")
    if UNLABELED_CATEGORY not in adata.obs[L2_KEY].cat.categories:
        adata.obs[L2_KEY] = adata.obs[L2_KEY].cat.add_categories([UNLABELED_CATEGORY])


def build_non_nose_gcb_mask(labels: pd.Series, tissues: pd.Series) -> pd.Series:
    label_str = labels.astype(str).str.strip()
    tissue_str = tissues.astype(str).str.strip().str.lower()
    return label_str.isin(GC_B_L3_LABELS) & tissue_str.ne(NOSE_TISSUE)


def collapse_non_nose_gcb_to_memory(
    adata: sc.AnnData,
    source_col: str,
    backup_col: str | None = None,
) -> int:
    if source_col not in adata.obs.columns:
        raise KeyError(f"Missing source column for non-nose GC-B remapping: {source_col}")
    if TISSUE_KEY not in adata.obs.columns:
        raise KeyError(f"Missing tissue column for non-nose GC-B remapping: {TISSUE_KEY}")

    labels = pd.Series(adata.obs[source_col], index=adata.obs.index, dtype="object")
    if backup_col is not None and backup_col not in adata.obs.columns:
        adata.obs[backup_col] = labels.astype(str)

    remap_mask = build_non_nose_gcb_mask(labels, adata.obs[TISSUE_KEY])
    n_remapped = int(remap_mask.sum())
    if n_remapped > 0:
        labels.loc[remap_mask] = NON_NOSE_GCB_MEMORY_TARGET
    adata.obs[source_col] = pd.Categorical(labels.astype(str))
    return n_remapped


def collapse_non_nose_gcb_predictions_to_memory(
    l3_pred: pd.Series,
    tissues: pd.Series,
) -> tuple[pd.Series, int]:
    pred = pd.Series(l3_pred, index=l3_pred.index, dtype="object")
    remap_mask = build_non_nose_gcb_mask(pred, tissues)
    n_remapped = int(remap_mask.sum())
    if n_remapped > 0:
        pred.loc[remap_mask] = NON_NOSE_GCB_MEMORY_TARGET
    return pred.astype(str), n_remapped


def prepare_scanvi_l3_training_labels(adata: sc.AnnData) -> None:
    if L3_KEY not in adata.obs.columns:
        raise KeyError(f"Missing L3 source column for scANVI training: {L3_KEY}")

    labels = pd.Series(adata.obs[L3_KEY], index=adata.obs.index)
    labels = labels.where(~labels.isna(), UNLABELED_CATEGORY)
    labels = labels.astype(str).str.strip()
    labels = labels.replace({"": UNLABELED_CATEGORY, "nan": UNLABELED_CATEGORY, "None": UNLABELED_CATEGORY})
    adata.obs[SCANVI_LABELS_KEY] = pd.Categorical(labels)
    if UNLABELED_CATEGORY not in adata.obs[SCANVI_LABELS_KEY].cat.categories:
        adata.obs[SCANVI_LABELS_KEY] = adata.obs[SCANVI_LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])


def derive_l2_from_l3_predictions(l3_pred: pd.Series) -> pd.Series:
    l2_pred = l3_pred.map(L3_TO_L2_MAP)
    if l2_pred.isna().any():
        missing = sorted(pd.unique(l3_pred[l2_pred.isna()]).tolist())
        raise ValueError(f"Unmapped L3 predictions after scANVI rerun: {missing}")
    return l2_pred


def ensure_batch_tissue(adata: sc.AnnData) -> None:
    if BATCH_KEY not in adata.obs.columns:
        adata.obs[BATCH_KEY] = "ref_batch"
    if TISSUE_KEY not in adata.obs.columns:
        adata.obs[TISSUE_KEY] = "ref_tissue"

    adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].astype("category")
    adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].astype("category")
    if "unknown_tissue" not in adata.obs[TISSUE_KEY].cat.categories:
        adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].cat.add_categories(["unknown_tissue"])
    adata.obs[TISSUE_KEY] = adata.obs[TISSUE_KEY].fillna("unknown_tissue")


def ensure_pct_mt(adata: sc.AnnData) -> None:
    for col in ["pct_counts_mt", "percent.mt", "percent_mt"]:
        if col in adata.obs.columns:
            if col != "pct_counts_mt":
                adata.obs["pct_counts_mt"] = adata.obs[col]
            return
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, layer=COUNTS_LAYER)


def ensure_stress_score(adata: sc.AnnData) -> None:
    if "stress_score" in adata.obs.columns:
        return

    genes = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
    if len(genes) < 5:
        adata.obs["stress_score"] = 0.0
        return

    X = adata.layers[COUNTS_LAYER]
    idx = adata.var_names.get_indexer(genes)
    idx = idx[idx >= 0]
    sub = X[:, idx]

    if issparse(X):
        tot = np.asarray(X.sum(axis=1)).ravel()
    else:
        tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1.0

    if issparse(sub):
        scale = sparse.diags(1e4 / tot)
        sub_norm = scale @ sub
        sub_norm = sub_norm.tocsr(copy=True)
        sub_norm.data = np.log1p(sub_norm.data)
        expr = np.asarray(sub_norm.mean(axis=1)).ravel()
    else:
        sub_norm = (sub / tot[:, None]) * 1e4
        expr = np.log1p(sub_norm).mean(axis=1)

    mn, mx = float(expr.min()), float(expr.max())
    adata.obs["stress_score"] = (expr - mn) / (mx - mn) if mx > mn else 0.0


def ensure_cell_cycle_scores(adata: sc.AnnData) -> None:
    if all(k in adata.obs.columns for k in ["S_score", "G2M_score", "phase"]):
        return

    s_in = [g for g in S_GENES if g in adata.var_names]
    g_in = [g for g in G2M_GENES if g in adata.var_names]
    if len(s_in) < 10 or len(g_in) < 10:
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"
        return

    ad_tmp = sc.AnnData(X=adata.layers[COUNTS_LAYER].copy(), var=adata.var.copy())
    sc.pp.normalize_total(ad_tmp, target_sum=1e4)
    sc.pp.log1p(ad_tmp)
    sc.tl.score_genes_cell_cycle(ad_tmp, s_genes=s_in, g2m_genes=g_in)
    adata.obs["S_score"] = ad_tmp.obs["S_score"].values
    adata.obs["G2M_score"] = ad_tmp.obs["G2M_score"].values
    adata.obs["phase"] = ad_tmp.obs["phase"].values
    del ad_tmp
    gc.collect()


def prepare_covariates(adata: sc.AnnData) -> None:
    ensure_counts_layer(adata)
    ensure_batch_tissue(adata)
    ensure_pct_mt(adata)
    ensure_stress_score(adata)
    ensure_cell_cycle_scores(adata)


def select_hvgs(adata: sc.AnnData) -> tuple[sc.AnnData, str]:
    hvg_method = "unknown"
    n_top = min(N_HVG, adata.n_vars)
    n_batches = adata.obs[BATCH_KEY].nunique()

    if n_batches > 1:
        try:
            sc.pp.highly_variable_genes(
                adata,
                layer=COUNTS_LAYER,
                n_top_genes=n_top,
                batch_key=BATCH_KEY,
                flavor="seurat_v3",
                subset=False,
            )
            hvg_method = "batch-aware seurat_v3"
        except Exception as exc:
            print(f"[WARN] batch-aware seurat_v3 failed: {exc}")

    if hvg_method == "unknown":
        try:
            sc.pp.highly_variable_genes(
                adata,
                layer=COUNTS_LAYER,
                n_top_genes=n_top,
                flavor="seurat_v3",
                subset=False,
            )
            hvg_method = "seurat_v3"
        except Exception as exc:
            print(f"[WARN] seurat_v3 failed: {exc}")

    if hvg_method == "unknown":
        sc.pp.highly_variable_genes(
            adata,
            layer=COUNTS_LAYER,
            n_top_genes=n_top,
            flavor="cell_ranger",
            subset=False,
        )
        hvg_method = "cell_ranger"

    before = int(adata.var["highly_variable"].sum())
    for gene in FORCED_MARKERS_BCELL:
        if gene in adata.var_names:
            adata.var.loc[gene, "highly_variable"] = True
    after = int(adata.var["highly_variable"].sum())
    print(f"[OK] HVG method: {hvg_method} | initial={before:,} | final={after:,}")

    hvg_genes = adata.var_names[adata.var["highly_variable"]].tolist()
    OUTPUT_HVG.write_text("\n".join(hvg_genes), encoding="utf-8")
    adata = adata[:, adata.var["highly_variable"]].copy()
    return adata, hvg_method


def _to_legacy_string_object(values: pd.Series) -> pd.Series:
    return pd.Series(
        [None if pd.isna(x) else str(x) for x in values.tolist()],
        index=values.index,
        dtype="object",
    )


def sanitize_dataframe_for_h5ad(df: pd.DataFrame) -> None:
    for col in df.columns:
        s = df[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            if str(s.cat.categories.dtype) == "string":
                categories = [str(x) for x in s.cat.categories.tolist()]
                values = [np.nan if pd.isna(x) else str(x) for x in s.astype(object).tolist()]
                df[col] = pd.Categorical(values, categories=categories, ordered=s.cat.ordered)
            continue
        if is_string_dtype(s):
            df[col] = _to_legacy_string_object(s)
            continue
        if not is_object_dtype(s):
            continue
        non_null = s.dropna()
        if len(non_null) == 0:
            df[col] = pd.Series([None] * len(s), index=s.index, dtype="object")
            continue
        py_types = set(non_null.map(lambda x: type(x).__name__))
        if py_types <= {"bool"}:
            if s.isna().any():
                df[col] = s.map(lambda x: -1 if pd.isna(x) else int(bool(x))).astype(np.int8)
            else:
                df[col] = s.astype(bool)
        elif py_types <= {"str"}:
            df[col] = _to_legacy_string_object(s)
        else:
            df[col] = _to_legacy_string_object(s)


def sanitize_anndata_for_legacy_h5ad(adata: sc.AnnData) -> None:
    sanitize_dataframe_for_h5ad(adata.obs)
    sanitize_dataframe_for_h5ad(adata.var)
    if adata.raw is not None:
        sanitize_dataframe_for_h5ad(adata.raw.var)


def main() -> None:
    accelerator, devices, device_name = resolve_training_device()

    print_banner("B-cell full fresh scVI/scANVI rerun after removing CHOIR clusters 22, 13, and 25")
    print(f"Input reference : {INPUT_REF_H5AD}")
    print(f"Output dir      : {OUTPUT_DIR}")
    print(f"Training device : {accelerator}")
    print(f"Torch threads   : intra={torch.get_num_threads()} interop={TORCH_NUM_INTEROP_THREADS}")
    if device_name is not None:
        print(f"GPU device      : {device_name}")

    print_banner("Step 1: load original reference h5ad + target CHOIR cells")
    adata = sc.read_h5ad(INPUT_REF_H5AD)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    print(f"Loaded shape    : {adata.shape}")
    union_df = load_union_cluster_rows()
    describe_cluster_targets(union_df)

    reset_output_dir(OUTPUT_DIR)
    union_df.to_csv(OUTPUT_TARGETS, index=False)
    print(f"[OK] target manifest saved: {OUTPUT_TARGETS}")

    adata, removed = remove_choir_cluster_cells(adata, union_df)
    removed.to_csv(OUTPUT_REMOVED, index=False)
    print(f"Filtered shape  : {adata.shape}")
    print(f"[OK] removed-cell table saved: {OUTPUT_REMOVED}")

    print_banner("Step 2: ensure counts/raw + historical labels")
    ensure_counts_layer(adata)
    ensure_raw_from_counts(adata)
    n_l3_training_remapped = collapse_non_nose_gcb_to_memory(
        adata,
        source_col=L3_KEY,
        backup_col=L3_INPUT_BACKUP_KEY,
    )
    map_l3_to_l2(adata)
    prepare_scanvi_l3_training_labels(adata)
    prepare_covariates(adata)
    print(
        f"[OK] non-nose GC-B -> {NON_NOSE_GCB_MEMORY_TARGET} applied to {n_l3_training_remapped:,} "
        f"cells in {L3_KEY} before scANVI training"
    )
    print("L3 training-label distribution:")
    print(adata.obs[SCANVI_LABELS_KEY].value_counts().to_string())
    print("\nL2 distribution (derived from expert L3):")
    print(adata.obs[L2_KEY].value_counts().to_string())

    print_banner("Step 3: HVG selection")
    adata_hvg, hvg_method = select_hvgs(adata)
    print(f"Training shape  : {adata_hvg.shape}")
    assert adata_hvg.raw is not None, "raw lost after HVG subsetting"
    assert COUNTS_LAYER in adata_hvg.layers, "counts layer lost after HVG subsetting"

    print_banner("Step 4: fresh scVI training")
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer=COUNTS_LAYER,
        batch_key=BATCH_KEY,
        continuous_covariate_keys=["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
        categorical_covariate_keys=[TISSUE_KEY],
    )
    scvi_model = scvi.model.SCVI(
        adata_hvg,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        n_hidden=SCVI_N_HIDDEN,
        dropout_rate=SCVI_DROPOUT_RATE,
        gene_likelihood=SCVI_GENE_LIKELIHOOD,
    )
    scvi_train_kwargs = {
        "max_epochs": MAX_EPOCHS_SCVI,
        "batch_size": SCVI_BATCH_SIZE,
        "early_stopping": True,
        "early_stopping_patience": 30,
        "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
        "accelerator": accelerator,
        "devices": 1,
    }
    if accelerator == "gpu" and devices is not None:
        scvi_train_kwargs["devices"] = devices
    scvi_model.train(**scvi_train_kwargs)
    scvi_model.save(OUTPUT_SCVI_DIR, overwrite=True)
    print(f"[OK] Fresh scVI model trained + saved: {OUTPUT_SCVI_DIR}")

    print_banner("Step 5: fresh scANVI training")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata_hvg,
        labels_key=SCANVI_LABELS_KEY,
        unlabeled_category=UNLABELED_CATEGORY,
    )
    scanvi_train_kwargs = {
        "max_epochs": MAX_EPOCHS_SCANVI,
        "batch_size": SCANVI_BATCH_SIZE,
        "early_stopping": True,
        "early_stopping_patience": SCANVI_EARLY_STOPPING_PATIENCE,
        "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
        "accelerator": accelerator,
        "devices": 1,
    }
    if accelerator == "gpu" and devices is not None:
        scanvi_train_kwargs["devices"] = devices
    scanvi_model.train(**scanvi_train_kwargs)
    scanvi_model.save(OUTPUT_SCANVI_DIR, overwrite=True)
    print(f"[OK] Fresh scANVI model trained + saved: {OUTPUT_SCANVI_DIR}")

    print_banner("Step 6: export latent + UMAP")
    adata_hvg.obsm["X_scvi"] = scvi_model.get_latent_representation(adata_hvg)
    adata_hvg.obsm["X_scanvi"] = scanvi_model.get_latent_representation(adata_hvg)
    adata_hvg.obsm["X_scanvi_corrected"] = np.asarray(adata_hvg.obsm["X_scanvi"]).copy()
    adata_hvg.obsm["X_scANVI_L3"] = np.asarray(adata_hvg.obsm["X_scanvi"]).copy()
    adata_hvg.obsm["X_scANVI_L2"] = np.asarray(adata_hvg.obsm["X_scanvi"]).copy()
    l3_pred = pd.Series(scanvi_model.predict(adata_hvg), index=adata_hvg.obs_names, dtype="object")
    if LEGACY_L3_OUTPUT_KEY in adata_hvg.obs.columns and LEGACY_L3_INPUT_BACKUP_KEY not in adata_hvg.obs.columns:
        adata_hvg.obs[LEGACY_L3_INPUT_BACKUP_KEY] = adata_hvg.obs[LEGACY_L3_OUTPUT_KEY].astype(str)
    if L2_KEY in adata_hvg.obs.columns and LEGACY_L2_INPUT_BACKUP_KEY not in adata_hvg.obs.columns:
        adata_hvg.obs[LEGACY_L2_INPUT_BACKUP_KEY] = adata_hvg.obs[L2_KEY].astype(str)
    l3_pred, n_prediction_guardrail = collapse_non_nose_gcb_predictions_to_memory(
        l3_pred,
        adata_hvg.obs[TISSUE_KEY],
    )
    adata_hvg.obs[L3_PRED_KEY] = pd.Categorical(l3_pred.astype(str))
    adata_hvg.obs[LEGACY_L3_OUTPUT_KEY] = pd.Categorical(l3_pred.astype(str))
    l2_pred = derive_l2_from_l3_predictions(l3_pred.astype(str))
    adata_hvg.obs["L2_scanvi_pred"] = pd.Categorical(l2_pred.astype(str))
    adata_hvg.obs[L2_KEY] = pd.Categorical(l2_pred.astype(str))
    proba = scanvi_model.predict(adata_hvg, soft=True).astype(np.float32)
    adata_hvg.obsm[L3_PROBA_KEY] = proba
    adata_hvg.obs[L3_CONFIDENCE_KEY] = np.max(proba, axis=1)
    adata_hvg.obs["L2_scanvi_confidence"] = adata_hvg.obs[L3_CONFIDENCE_KEY].values
    print(
        f"[OK] prediction guardrail remapped {n_prediction_guardrail:,} non-nose GC-B predictions "
        f"to {NON_NOSE_GCB_MEMORY_TARGET}"
    )

    umap_op = umap.UMAP(
        n_neighbors=UMAP_N_NEIGHBORS,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=RANDOM_SEED,
    ).fit(np.asarray(adata_hvg.obsm["X_scanvi"]))
    adata_hvg.obsm["X_umap_scanvi"] = np.asarray(umap_op.embedding_).copy()
    adata_hvg.obsm["X_umap_scanvi_corrected"] = np.asarray(umap_op.embedding_).copy()
    adata_hvg.obsm["X_umap"] = np.asarray(umap_op.embedding_).copy()
    joblib.dump(umap_op, OUTPUT_UMAP)
    print(f"[OK] UMAP operator saved: {OUTPUT_UMAP}")

    print_banner("Step 7: save filtered reference h5ad")
    sanitize_anndata_for_legacy_h5ad(adata_hvg)
    if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
        ad.settings.allow_write_nullable_strings = False
    adata_hvg.write_h5ad(OUTPUT_H5AD, compression="gzip")
    print(f"[OK] h5ad saved: {OUTPUT_H5AD}")

    verify = sc.read_h5ad(OUTPUT_H5AD, backed="r")
    has_raw = verify.raw is not None
    has_counts = COUNTS_LAYER in verify.layers.keys()
    print(f"[OK] verify raw={has_raw} | counts_layer={has_counts} | shape={verify.shape}")
    verify.file.close()

    cluster_source_summary = (
        union_df.groupby(["cluster_source_label", "choir_cluster"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .to_dict(orient="records")
    )
    config = {
        "timestamp": datetime.now().isoformat(),
        "input_ref_h5ad": str(INPUT_REF_H5AD),
        "choir_cluster_sources": [
            {
                "csv": str(spec["csv"]),
                "clusters": [int(x) for x in spec["clusters"]],
                "source_label": str(spec["source_label"]),
            }
            for spec in CHOIR_CLUSTER_SOURCES
        ],
        "cluster_source_summary": cluster_source_summary,
        "n_target_cells_union": int(union_df["cell"].nunique()),
        "n_removed_cells": int(removed.shape[0]),
        "n_cells_after_filter": int(adata_hvg.n_obs),
        "n_vars_after_hvg": int(adata_hvg.n_vars),
        "raw_n_vars": int(adata_hvg.raw.n_vars) if adata_hvg.raw is not None else None,
        "hvg_method": hvg_method,
        "scanvi_labels_key": SCANVI_LABELS_KEY,
        "scanvi_l3_prediction_col": L3_PRED_KEY,
        "scanvi_l3_output_col": LEGACY_L3_OUTPUT_KEY,
        "l3_training_distribution": adata_hvg.obs[SCANVI_LABELS_KEY].astype(str).value_counts().to_dict(),
        "l2_distribution": adata_hvg.obs[L2_KEY].astype(str).value_counts().to_dict(),
        "non_nose_gc_b_target": NON_NOSE_GCB_MEMORY_TARGET,
        "non_nose_gc_b_labels": sorted(GC_B_L3_LABELS),
        "non_nose_gc_b_training_relabeled_n": int(n_l3_training_remapped),
        "non_nose_gc_b_prediction_guardrail_n": int(n_prediction_guardrail),
        "scvi_model_dir": str(OUTPUT_SCVI_DIR),
        "scanvi_model_dir": str(OUTPUT_SCANVI_DIR),
        "output_h5ad": str(OUTPUT_H5AD),
        "scvi_batch_size": SCVI_BATCH_SIZE,
        "scanvi_batch_size": SCANVI_BATCH_SIZE,
        "scanvi_early_stopping_patience": SCANVI_EARLY_STOPPING_PATIENCE,
        "torch_num_threads": torch.get_num_threads(),
        "torch_num_interop_threads": TORCH_NUM_INTEROP_THREADS,
        "fresh_scvi": True,
        "fresh_scanvi": True,
    }
    with open(OUTPUT_CONFIG, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)
    print(f"[OK] config saved: {OUTPUT_CONFIG}")

    print_banner("DONE")
    print(f"Removed cells     : {removed.shape[0]:,}")
    print(f"Final h5ad        : {OUTPUT_H5AD}")
    print(f"scVI model        : {OUTPUT_SCVI_DIR}")
    print(f"scANVI model      : {OUTPUT_SCANVI_DIR}")


if __name__ == "__main__":
    main()
