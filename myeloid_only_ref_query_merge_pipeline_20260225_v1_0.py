#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
Myeloid ONLY Reference-Query Merge Pipeline v1.0
==============================================================================

Purpose:
--------
Specialized pipeline for Myeloid cell analysis with focus on monocyte/macrophage/
DC/neutrophil distinction.

Key Features:
-------------
1. Myeloid cell specific marker prioritization (CD14, CD16, CD68, CD1C, etc.)
2. Monocyte subtype scoring (Classical vs Non-classical) for validation
3. Coarse (Myeloid lineage) + Fine (subtypes) dual annotation
4. Direct merge (NOT scArches) - allows discovery of new myeloid subsets
5. Reference cells keep labels, query cells marked as "Unknown"

Myeloid Classification Strategy:
--------------------------------
- Core lineage markers: LYZ, CD14, CD33, PTPRC (CD45)
- Monocyte subtypes: CD14 (Classical), FCGR3A/CD16 (Non-classical)
- Macrophage markers: CD68, CD163, MRC1 (CD206)
- DC markers: CD1C (cDC2), CLEC9A (cDC1), XCR1
- Neutrophil markers: S100A8/A9, FCGR3B, CSF3R
- Mast cells: KIT, TPSAB1

Input:
------
- Reference h5ad: Myeloid cells with annotations (coarse + fine labels preferred)
- Query h5ad: Myeloid cells raw counts, may contain new subsets

Output:
-------
- Merged h5ad with scVI/scANVI latent spaces
- Monocyte subtype scores for validation
- Coarse (lineage) and Fine (subtypes) annotations
- CellTypist predictions (run on full genes, results backfilled)
- UMAP visualizations (scVI, scANVI)
- Novelty detection with query-only Leiden clustering

Based on: tcell_only_ref_query_merge_pipeline_20260225_v1_2.py
Author: Claude Code
Date: 2026-02-25
Version: 1.0
==============================================================================
"""

import sys
import warnings
import json
import gc
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Dict, Any

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
from scipy import sparse
from scipy.stats import entropy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
import scvi
import celltypist
from celltypist import models

import torch

try:
    from scanvi_umap_bundle_helper_20260419_v1 import fit_bundle
except ModuleNotFoundError:
    import importlib.util

    _SCANVI_UMAP_HELPER_PATH = Path(__file__).resolve().with_name(
        "scanvi_umap_bundle_helper_20260419_v1.py"
    )
    _SCANVI_UMAP_HELPER_SPEC = importlib.util.spec_from_file_location(
        "scanvi_umap_bundle_helper_20260419_v1",
        _SCANVI_UMAP_HELPER_PATH,
    )
    if _SCANVI_UMAP_HELPER_SPEC is None or _SCANVI_UMAP_HELPER_SPEC.loader is None:
        raise
    _scanvi_umap_helper = importlib.util.module_from_spec(_SCANVI_UMAP_HELPER_SPEC)
    _SCANVI_UMAP_HELPER_SPEC.loader.exec_module(_scanvi_umap_helper)
    fit_bundle = _scanvi_umap_helper.fit_bundle

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None

print("=" * 80)
print("Myeloid ONLY Ref-Query Merge Pipeline v1.0")
print("=" * 80)

gpu_available = torch.cuda.is_available()
print(f"GPU available: {gpu_available}")
if gpu_available:
    print(f"GPU Device: {torch.cuda.get_device_name(0)}")

scvi.settings.dl_num_workers = 0
print(f"scVI dl_num_workers: {scvi.settings.dl_num_workers}")

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# --- Input Files (EDIT THESE) ---
REFERENCE_H5AD = "/home/h2048/data/py/0129/tnk_analysis_unified/results/subcluster_unified_v2_20260129/adata_tnk_subclustered_FINAL_v2_0_1_20260129.h5ad"
QUERY_H5AD = "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/t_cells.h5ad"

# --- Column Configuration ---
REF_LABEL_COARSE = "cell_type_L2"   # e.g., "Mono", "Macro", "DC", "Neutro"
REF_LABEL_FINE = "cell_type_L3"     # e.g., "Mono_CD14", "Mono_FCGR3A", "cDC1", etc.

BATCH_KEY = "sample"
TISSUE_KEY = "tissue"

# --- Output Configuration ---
OUTPUT_DIR = "/home/h2048/data/py/20260225/myeloid_only_merged_pipeline"
OUTPUT_PREFIX = "myeloid_only_merged"

# --- Cell Type Selection ---
INCLUDE_COARSE_TYPES = None  # e.g., ["Mono", "Macro"]

# --- Pipeline Parameters ---
N_HVG = 4000
FORCE_MARKERS_IN_HVG = True

# scVI Parameters
SCVI_N_LATENT = 100
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_DROPOUT = 0.1
MAX_EPOCHS_SCVI = 400

# scANVI Parameters
MAX_EPOCHS_SCANVI = 200
UNLABELED_CATEGORY = "Unknown"

# Training Parameters
BATCH_SIZE = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 0.0

# CellTypist Parameters
CELLTYPIST_MODEL = "/home/h2048/data/source/reference/celltypist_models/Immune_All_Low.pkl"
CELLTYPIST_MAJORITY_VOTE = True

# Query-only clustering (for novelty detection)
QUERY_LEIDEN_RESOLUTION = 1.0

# --- Myeloid Gene Signatures ---

# Score Thresholds for Monocyte Classification
CLASSICAL_MONO_THRESHOLD = 0.3
NONCLASSICAL_MONO_THRESHOLD = 0.3

# Core Myeloid lineage markers
MYELOID_CORE_MARKERS = ["LYZ", "CD14", "CD33", "PTPRC", "ITGAM", "ITGAX"]

# Monocyte subtype markers
CLASSICAL_MONO_MARKERS = ["CD14", "FCGR1A", "CCR2", "CD36", "SELL"]  # CD14++ CD16-
NONCLASSICAL_MONO_MARKERS = ["FCGR3A", "CD16", "FCER1G", "PICALM", "RHOC"]  # CD14+ CD16++
INTERMEDIATE_MONO_MARKERS = ["CD14", "FCGR3A", "CD86", "HLA-DRA"]  # CD14++ CD16+

# Macrophage markers
MACROPHAGE_MARKERS = ["CD68", "CD163", "MRC1", "MARCO", "MSR1", "FCGR2A", "APOE", "C1QA"]
M1_MACROPHAGE_MARKERS = ["CD86", "CD80", "TNF", "IL1B", "NOS2", "CXCL9", "CXCL10"]
M2_MACROPHAGE_MARKERS = ["CD163", "MRC1", "ARG1", "IL10", "TGFB1", "CCL22", "PPARG"]

# Dendritic Cell markers
CDC1_MARKERS = ["CLEC9A", "XCR1", "WDFY4", "IRF8", "BATF3"]  # Cross-presenting DCs
CDC2_MARKERS = ["CD1C", "FCER1A", "CLEC10A", "CD2", "ESAM"]  # Conventional DCs
PDC_MARKERS = ["LILRA4", "CLEC4C", "NRP1", "TCF4", "IRF7", "GZMB"]  # Plasmacytoid DCs

# Neutrophil markers
NEUTROPHIL_MARKERS = ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CEACAM8", "CXCR2", "FCN1"]

# Mast cell markers
MAST_CELL_MARKERS = ["KIT", "TPSAB1", "TPSB2", "HPGDS", "MS4A2", "FCER1A", "CPA3"]

# Eosinophil markers
EOSINOPHIL_MARKERS = ["EPX", "PRG2", "CLC", "CCL26", "SIGLEC8"]

# Proliferation markers
PROLIF_MARKERS = ["MKI67", "TOP2A", "PCNA"]

# Stress signature
STRESS_SIGNATURE_GENES = [
    "HSPA1A", "HSPA1B", "HSPA8", "HSP90AA1", "HSP90AB1", "DNAJB1",
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "IER2"
]

# Cell cycle genes
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UHRF1",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "HELLS", "RFC2", "RPA2",
    "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7",
    "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1"
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "SMC4", "CCNB1", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E"
]

# All forced markers for HVG selection
FORCED_MARKERS = list(set(
    MYELOID_CORE_MARKERS + 
    CLASSICAL_MONO_MARKERS + NONCLASSICAL_MONO_MARKERS + INTERMEDIATE_MONO_MARKERS +
    MACROPHAGE_MARKERS + M1_MACROPHAGE_MARKERS + M2_MACROPHAGE_MARKERS +
    CDC1_MARKERS + CDC2_MARKERS + PDC_MARKERS +
    NEUTROPHIL_MARKERS + MAST_CELL_MARKERS + EOSINOPHIL_MARKERS +
    PROLIF_MARKERS
))

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.seed = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

def ensure_counts_layer(adata, counts_layer="counts"):
    """Robust counts validation with float32 tolerance."""
    if counts_layer not in (adata.layers or {}):
        print(f"  WARNING: layers['{counts_layer}'] not found, checking .X...")
        if hasattr(adata, 'X') and adata.X is not None:
            if issparse(adata.X):
                X_sample = adata.X[:1000].toarray().flatten()
            else:
                X_sample = adata.X.flat[:1000]
            sample = np.asarray(X_sample, dtype=np.float64)

            if np.any(sample < 0):
                raise ValueError("adata.X contains negative values - not valid raw counts!")

            if np.allclose(sample, np.round(sample), atol=1e-6):
                print(f"  -> Auto-copying .X to layers['{counts_layer}']")
                if issparse(adata.X) and not isinstance(adata.X, csr_matrix):
                    adata.layers[counts_layer] = csr_matrix(adata.X)
                else:
                    adata.layers[counts_layer] = adata.X.copy()
            else:
                raise ValueError(
                    f"CRITICAL ERROR: layers['{counts_layer}'] not found and .X "
                    f"does not look like raw counts (non-integer values). "
                    f"Sample range: [{sample.min():.4f}, {sample.max():.4f}]"
                )
        else:
            raise ValueError(
                f"CRITICAL ERROR: layers['{counts_layer}'] not found and .X is None. "
                "scVI requires raw counts!"
            )

    X_counts = adata.layers[counts_layer]
    if issparse(X_counts):
        sample_data = X_counts.data[:1000]
    else:
        sample_data = X_counts.flat[:1000]

    sample = np.asarray(sample_data, dtype=np.float64)

    if np.any(sample < 0):
        raise ValueError("counts contains negative values!")

    if not np.allclose(sample, np.round(sample), atol=1e-6):
        raise ValueError(
            "counts looks non-integer (possible normalized/log data). "
            f"Sample range: [{sample.min():.4f}, {sample.max():.4f}]"
        )

    if issparse(X_counts) and not isinstance(X_counts, csr_matrix):
        adata.layers[counts_layer] = csr_matrix(X_counts)

    return counts_layer


def ensure_batch_tissue(adata, batch_key, tissue_key):
    """Ensure batch and tissue columns exist and are categorical."""
    if batch_key not in adata.obs.columns:
        print(f"  WARNING: {batch_key} not found, creating placeholder")
        adata.obs[batch_key] = "unknown_batch"
    adata.obs[batch_key] = adata.obs[batch_key].astype("category")

    if tissue_key not in adata.obs.columns:
        print(f"  WARNING: {tissue_key} not found, creating placeholder")
        adata.obs[tissue_key] = "unknown_tissue"
    adata.obs[tissue_key] = adata.obs[tissue_key].astype("category")
    if "unknown_tissue" not in adata.obs[tissue_key].cat.categories:
        adata.obs[tissue_key] = adata.obs[tissue_key].cat.add_categories(["unknown_tissue"])
    adata.obs[tissue_key] = adata.obs[tissue_key].fillna("unknown_tissue")


def compute_module_score_efficient(adata, gene_list, score_name, counts_layer="counts"):
    """
    Memory-efficient module score computation using direct mean-expression.
    """
    genes = [g for g in gene_list if g in adata.var_names]
    if len(genes) < 3:
        adata.obs[score_name] = 0.0
        adata.obs[f"{score_name}_norm"] = 0.0
        return

    gene_idx = adata.var_names.get_indexer(genes)
    gene_idx = gene_idx[gene_idx >= 0]

    if len(gene_idx) < 3:
        adata.obs[score_name] = 0.0
        adata.obs[f"{score_name}_norm"] = 0.0
        return

    X_subset = adata.layers[counts_layer][:, gene_idx]
    var_subset = adata.var.iloc[gene_idx].copy()

    adata_tmp = sc.AnnData(X=X_subset.copy(), var=var_subset)
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)

    X_norm = adata_tmp.X
    if issparse(X_norm):
        mean_expr = np.asarray(X_norm.mean(axis=1)).flatten()
    else:
        mean_expr = np.asarray(X_norm).mean(axis=1).flatten()

    adata.obs[score_name] = mean_expr

    scores = adata.obs[score_name]
    mn, mx = float(scores.min()), float(scores.max())
    if mx > mn:
        adata.obs[f"{score_name}_norm"] = (scores - mn) / (mx - mn)
    else:
        adata.obs[f"{score_name}_norm"] = 0.0

    del adata_tmp, X_subset
    gc.collect()


def compute_monocyte_subtype_scores(adata):
    """
    Compute monocyte subtype scores for validation.
    Classifies cells into Classical (CD14++ CD16-), Non-classical (CD14+ CD16++), 
    or Intermediate based on marker expression.
    """
    print("  -> Computing monocyte subtype scores...")

    compute_module_score_efficient(adata, CLASSICAL_MONO_MARKERS, "Classical_Mono_score")
    compute_module_score_efficient(adata, NONCLASSICAL_MONO_MARKERS, "NonClassical_Mono_score")
    compute_module_score_efficient(adata, INTERMEDIATE_MONO_MARKERS, "Intermediate_Mono_score")

    classical_high = adata.obs["Classical_Mono_score_norm"] > CLASSICAL_MONO_THRESHOLD
    nonclassical_high = adata.obs["NonClassical_Mono_score_norm"] > NONCLASSICAL_MONO_THRESHOLD
    intermediate_high = adata.obs["Intermediate_Mono_score_norm"] > 0.3

    classifications = []
    for i in range(len(adata)):
        c_high = classical_high.iloc[i]
        nc_high = nonclassical_high.iloc[i]
        int_high = intermediate_high.iloc[i]
        
        if c_high and not nc_high:
            classifications.append("Classical")
        elif nc_high and not c_high:
            classifications.append("NonClassical")
        elif c_high and nc_high:
            if int_high:
                classifications.append("Intermediate")
            else:
                classifications.append("DoublePositive")
        else:
            classifications.append("Negative")

    adata.obs["mono_subtype_by_score"] = classifications
    print(f"    Classical: {sum([c == 'Classical' for c in classifications])}")
    print(f"    NonClassical: {sum([c == 'NonClassical' for c in classifications])}")
    print(f"    Intermediate: {sum([c == 'Intermediate' for c in classifications])}")
    print(f"    DoublePositive: {sum([c == 'DoublePositive' for c in classifications])}")
    print(f"    Negative: {sum([c == 'Negative' for c in classifications])}")


def prepare_covariates(adata):
    """Prepare covariates for scVI training."""
    print("  -> Validating counts...")
    ensure_counts_layer(adata, "counts")

    print("  -> Batch/Tissue...")
    ensure_batch_tissue(adata, BATCH_KEY, TISSUE_KEY)

    print("  -> Computing signature scores...")
    if "pct_counts_mt" not in adata.obs.columns:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, layer="counts")

    compute_monocyte_subtype_scores(adata)
    compute_module_score_efficient(adata, MACROPHAGE_MARKERS, "Macrophage_score")
    compute_module_score_efficient(adata, NEUTROPHIL_MARKERS, "Neutrophil_score")
    compute_module_score_efficient(adata, STRESS_SIGNATURE_GENES, "stress_score")

    if not all(k in adata.obs.columns for k in ["S_score", "G2M_score"]):
        s_in = [g for g in S_GENES if g in adata.var_names]
        g_in = [g for g in G2M_GENES if g in adata.var_names]

        if len(s_in) >= 5 and len(g_in) >= 5:
            cc_union = list(dict.fromkeys(s_in + g_in))
            cc_idx = adata.var_names.get_indexer(cc_union)
            X_cc = adata.layers["counts"][:, cc_idx].copy()
            var_cc = adata.var.iloc[cc_idx].copy()

            ad_tmp = sc.AnnData(X=X_cc, var=var_cc)
            sc.pp.normalize_total(ad_tmp, target_sum=1e4)
            sc.pp.log1p(ad_tmp)
            sc.tl.score_genes_cell_cycle(ad_tmp, s_genes=s_in, g2m_genes=g_in)

            adata.obs["S_score"] = ad_tmp.obs["S_score"].values
            adata.obs["G2M_score"] = ad_tmp.obs["G2M_score"].values
            adata.obs["phase"] = ad_tmp.obs["phase"].values

            del ad_tmp
            gc.collect()
        else:
            adata.obs["S_score"] = 0.0
            adata.obs["G2M_score"] = 0.0
            adata.obs["phase"] = "G1"


def run_celltypist_on_full_genes(adata_merged, hvg_mask, model_name=CELLTYPIST_MODEL, majority_vote=True):
    """
    Run CellTypist on full gene matrix (not HVG subset).
    """
    print("\n[CellTypist] Starting annotation on FULL gene matrix...")

    if os.path.isfile(model_name):
        print(f"  -> Loading model from local path: {model_name}")
        model = models.Model.load(model_name)
    else:
        model_basename = os.path.basename(model_name)
        print(f"  -> Local model not found, downloading: {model_basename}")
        models.download_models(model=model_basename)
        model = models.Model.load(model_basename)

    model_genes = set(model.features)

    available_genes = [g for g in adata_merged.var_names if g in model_genes]
    print(f"  -> {len(available_genes)}/{len(model_genes)} model genes present in data")

    if len(available_genes) < 100:
        raise ValueError(f"Too few overlapping genes ({len(available_genes)}) for CellTypist!")

    adata_ct = adata_merged[:, available_genes].copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    predictions = celltypist.annotate(
        adata_ct,
        model=model,
        majority_voting=majority_vote,
        mode='best match'
    )

    pl = predictions.predicted_labels

    if isinstance(pl, pd.DataFrame):
        if "predicted_labels" in pl.columns:
            pred_labels = pl["predicted_labels"].astype(str).values
        else:
            pred_labels = pl.iloc[:, 0].astype(str).values

        conf_values = predictions.probability_matrix.max(axis=1).values

        majority_labels = None
        if majority_vote and "majority_voting" in pl.columns:
            majority_labels = pl["majority_voting"].astype(str).values
    else:
        pred_labels = pl.astype(str).values
        conf_values = predictions.probability_matrix.max(axis=1).values
        majority_labels = None

    adata_merged.obs["celltypist_pred"] = pred_labels
    adata_merged.obs["celltypist_confidence"] = conf_values

    if majority_labels is not None:
        adata_merged.obs["celltypist_majority"] = majority_labels

    print(f"  -> CellTypist predictions added to full matrix")
    print(pd.Series(pred_labels).value_counts().head(10))

    del adata_ct
    gc.collect()
    return predictions


def compute_novelty_scores(adata_merged, proba_df):
    """
    Compute novelty scores for query cells based on prediction entropy.
    """
    print("  -> Computing novelty scores (entropy-based)...")

    proba_array = proba_df.values
    epsilon = 1e-10
    entropies = entropy(proba_array + epsilon, axis=1)

    adata_merged.obs["scanvi_entropy"] = entropies

    ent_min, ent_max = entropies.min(), entropies.max()
    if ent_max > ent_min:
        adata_merged.obs["novelty_score"] = (entropies - ent_min) / (ent_max - ent_min)
    else:
        adata_merged.obs["novelty_score"] = 0.0

    query_mask = adata_merged.obs["data_source"] == "query"
    high_novelty = (adata_merged.obs["novelty_score"] > 0.7) & query_mask
    adata_merged.obs["is_potentially_novel"] = high_novelty

    print(f"    High novelty query cells: {high_novelty.sum()}")


def run_query_only_leiden(adata_merged, resolution=QUERY_LEIDEN_RESOLUTION):
    """
    Run Leiden clustering on query cells only using scVI latent space.
    """
    print(f"\n[Novelty Detection] Running query-only Leiden (resolution={resolution})...")

    query_mask = adata_merged.obs["data_source"] == "query"
    query_cells = adata_merged.obs_names[query_mask]

    if len(query_cells) < 10:
        print("  -> Too few query cells, skipping query-only clustering")
        adata_merged.obs["leiden_query"] = "N/A"
        return

    X_scVI_query = adata_merged.obsm["X_scVI"][query_mask.values]

    adata_qry_tmp = sc.AnnData(
        X=X_scVI_query,
        obs=adata_merged.obs.loc[query_cells].copy()
    )
    adata_qry_tmp.obsm["X_scVI"] = X_scVI_query

    sc.pp.neighbors(adata_qry_tmp, use_rep="X_scVI", n_neighbors=30, random_state=RANDOM_SEED)
    sc.tl.leiden(adata_qry_tmp, resolution=resolution, random_state=RANDOM_SEED)

    leiden_full = pd.Series("N/A", index=adata_merged.obs_names, dtype="object")
    leiden_full.loc[query_cells] = "qry_" + adata_qry_tmp.obs["leiden"].astype(str)
    adata_merged.obs["leiden_query"] = leiden_full.values

    print(f"  -> Found {adata_qry_tmp.obs['leiden'].nunique()} query-only clusters")
    print(adata_merged.obs["leiden_query"].value_counts().head(10))

    del adata_qry_tmp
    gc.collect()


# ==============================================================================
# 3. MAIN PIPELINE
# ==============================================================================

def main():
    """Main pipeline function."""

    output_dir = Path(OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Load Reference and Query ---
    print("\n[Step 1] Loading Reference and Query...")

    print(f"  Loading reference: {REFERENCE_H5AD}")
    adata_ref = sc.read_h5ad(REFERENCE_H5AD)
    adata_ref.var_names_make_unique()
    print(f"    Reference shape: {adata_ref.shape}")

    print(f"  Loading query: {QUERY_H5AD}")
    adata_qry = sc.read_h5ad(QUERY_H5AD)
    adata_qry.var_names_make_unique()
    print(f"    Query shape: {adata_qry.shape}")

    # --- Step 2: Subset Reference if needed ---
    if INCLUDE_COARSE_TYPES is not None and REF_LABEL_COARSE in adata_ref.obs.columns:
        print(f"\n[Step 2] Subsetting reference to: {INCLUDE_COARSE_TYPES}")
        mask = adata_ref.obs[REF_LABEL_COARSE].isin(INCLUDE_COARSE_TYPES)
        adata_ref = adata_ref[mask].copy()
        print(f"    Reference after subset: {adata_ref.shape}")

    # --- Step 3: Prepare labels ---
    print("\n[Step 3] Preparing labels...")

    if REF_LABEL_COARSE in adata_ref.obs.columns:
        adata_ref.obs["cell_type_coarse"] = adata_ref.obs[REF_LABEL_COARSE].astype(str)
        print(f"  -> Coarse labels from: {REF_LABEL_COARSE}")
        print(adata_ref.obs["cell_type_coarse"].value_counts())
    else:
        print(f"  WARNING: {REF_LABEL_COARSE} not found, using 'Myeloid'")
        adata_ref.obs["cell_type_coarse"] = "Myeloid"

    if REF_LABEL_FINE in adata_ref.obs.columns:
        adata_ref.obs["cell_type_fine"] = adata_ref.obs[REF_LABEL_FINE].astype(str)
        print(f"  -> Fine labels from: {REF_LABEL_FINE}")
    else:
        adata_ref.obs["cell_type_fine"] = adata_ref.obs["cell_type_coarse"]
        print(f"  -> Using coarse labels as fine labels")

    adata_qry.obs["cell_type_coarse"] = UNLABELED_CATEGORY
    adata_qry.obs["cell_type_fine"] = UNLABELED_CATEGORY

    # --- Step 4: Find common genes (deterministic order) ---
    print("\n[Step 4] Finding common genes...")

    qry_set = set(adata_qry.var_names)
    common_genes = [g for g in adata_ref.var_names if g in qry_set]

    print(f"  Reference: {adata_ref.n_vars:,}, Query: {adata_qry.n_vars:,}, Common: {len(common_genes):,}")

    if len(common_genes) < 1000:
        raise ValueError(f"Too few common genes ({len(common_genes)}). Check gene naming!")

    adata_ref = adata_ref[:, common_genes].copy()
    adata_qry = adata_qry[:, common_genes].copy()

    # --- Step 5: Validate counts ---
    print("\n[Step 5] Validating counts...")
    ensure_counts_layer(adata_ref, "counts")
    ensure_counts_layer(adata_qry, "counts")

    # --- Step 6: Concatenate ---
    print("\n[Step 6] Concatenating...")

    adata_ref.obs_names = pd.Index([f"ref_{x}" for x in adata_ref.obs_names])
    adata_qry.obs_names = pd.Index([f"qry_{x}" for x in adata_qry.obs_names])

    adata_merged = sc.concat(
        {"reference": adata_ref, "query": adata_qry},
        axis=0,
        join="inner",
        merge="unique",
        label="data_source"
    )

    print(f"  Merged: {adata_merged.shape}")
    print(f"  Reference: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
    print(f"  Query: {(adata_merged.obs['data_source'] == 'query').sum():,}")

    del adata_ref, adata_qry
    gc.collect()

    # --- Step 7: Prepare covariates ---
    print("\n[Step 7] Preparing covariates...")
    prepare_covariates(adata_merged)

    print("  -> Adding symbol_base column...")
    adata_merged.var["symbol_base"] = adata_merged.var_names.str.replace(r"-\d+$", "", regex=True)

    # --- Step 8: HVG Selection ---
    print("\n[Step 8] Selecting HVGs...")
    hvg_method = "unknown"
    try:
        sc.pp.highly_variable_genes(
            adata_merged, layer="counts", n_top_genes=N_HVG,
            batch_key=BATCH_KEY, flavor="seurat_v3", subset=False
        )
        hvg_method = "batch_seurat_v3"
    except Exception as e1:
        print(f"  -> batch-aware failed ({str(e1)[:50]}), trying standard...")
        try:
            sc.pp.highly_variable_genes(
                adata_merged, layer="counts", n_top_genes=N_HVG,
                flavor="seurat_v3", subset=False
            )
            hvg_method = "standard_seurat_v3"
        except Exception as e2:
            print(f"  -> standard failed ({str(e2)[:50]}), fallback to cell_ranger")
            sc.pp.highly_variable_genes(
                adata_merged, layer="counts", n_top_genes=N_HVG,
                flavor="cell_ranger", subset=False
            )
            hvg_method = "cell_ranger"

    print(f"  -> Method: {hvg_method}")

    if FORCE_MARKERS_IN_HVG:
        n_added = 0
        marker_set = set(FORCED_MARKERS)
        for idx, symbol_base in enumerate(adata_merged.var["symbol_base"]):
            if symbol_base in marker_set:
                real_name = adata_merged.var_names[idx]
                if not adata_merged.var.loc[real_name, "highly_variable"]:
                    adata_merged.var.loc[real_name, "highly_variable"] = True
                    n_added += 1
        print(f"  -> Forced {n_added}/{len(FORCED_MARKERS)} markers into HVG")

    n_hvg_final = adata_merged.var["highly_variable"].sum()
    print(f"  -> Final HVG count: {n_hvg_final}")

    hvg_genes = adata_merged.var_names[adata_merged.var["highly_variable"]].tolist()
    with open(output_dir / f"{OUTPUT_PREFIX}_hvg_genes.txt", "w") as f:
        f.write("\n".join(hvg_genes))

    # --- Step 9: Build FULL matrix for .raw ---
    print("\n[Step 9] Building full matrix for .raw...")
    full_counts = adata_merged.layers["counts"]
    if issparse(full_counts) and not isinstance(full_counts, csr_matrix):
        full_counts = csr_matrix(full_counts)
    raw_var = adata_merged.var.copy()
    print(f"  Full matrix shape: {full_counts.shape}")

    # --- Step 10: Create training subset ---
    print("\n[Step 10] Creating training subset...")

    hvg_mask = adata_merged.var["highly_variable"].values
    X_hvg = adata_merged.layers["counts"][:, hvg_mask]

    if issparse(X_hvg) and not isinstance(X_hvg, csr_matrix):
        X_hvg = csr_matrix(X_hvg)

    adata_train = sc.AnnData(
        X=X_hvg.copy(),
        obs=adata_merged.obs.copy(),
        var=adata_merged.var.iloc[hvg_mask].copy()
    )
    adata_train.var_names = adata_merged.var_names[hvg_mask]
    adata_train.layers["counts"] = adata_train.X

    print(f"  Training data: {adata_train.shape}")

    adata_train.obs["scanvi_labels"] = adata_train.obs["cell_type_fine"].astype(str)
    adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")

    if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
        adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])

    print(f"  -> Label distribution:")
    print(adata_train.obs["scanvi_labels"].value_counts())

    gc.collect()

    # --- Step 11: CellTypist (full genes, backfill to training subset) ---
    print("\n[Step 11] Running CellTypist on full gene matrix...")
    try:
        run_celltypist_on_full_genes(adata_merged, hvg_mask)
        adata_train.obs["celltypist_pred"] = adata_merged.obs.loc[adata_train.obs_names, "celltypist_pred"]
        adata_train.obs["celltypist_confidence"] = adata_merged.obs.loc[adata_train.obs_names, "celltypist_confidence"]
        if "celltypist_majority" in adata_merged.obs.columns:
            adata_train.obs["celltypist_majority"] = adata_merged.obs.loc[adata_train.obs_names, "celltypist_majority"]
    except Exception as e:
        print(f"  WARNING: CellTypist failed: {e}")
        import traceback
        traceback.print_exc()

    # --- Step 12: scVI Training ---
    print("\n[Step 12] Training scVI...")

    setup_kwargs = {
        "layer": "counts",
        "batch_key": BATCH_KEY,
        "continuous_covariate_keys": ["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
        "categorical_covariate_keys": [TISSUE_KEY]
    }

    scvi.model.SCVI.setup_anndata(adata_train, **setup_kwargs)

    scvi_model = scvi.model.SCVI(
        adata_train,
        n_latent=SCVI_N_LATENT,
        n_layers=SCVI_N_LAYERS,
        n_hidden=SCVI_N_HIDDEN,
        dropout_rate=SCVI_DROPOUT
    )

    train_kwargs = {
        "max_epochs": MAX_EPOCHS_SCVI,
        "batch_size": BATCH_SIZE,
        "early_stopping": True,
        "early_stopping_patience": 30,
        "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
    }
    if gpu_available:
        train_kwargs["accelerator"] = "gpu"
        train_kwargs["devices"] = 1

    scvi_model.train(**train_kwargs)
    print("  -> scVI complete")

    # --- Step 13: scANVI Training ---
    print("\n[Step 13] Training scANVI...")

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata_train,
        labels_key="scanvi_labels",
        unlabeled_category=UNLABELED_CATEGORY
    )

    scanvi_train_kwargs = {
        "max_epochs": MAX_EPOCHS_SCANVI,
        "batch_size": BATCH_SIZE,
        "early_stopping": True,
        "early_stopping_patience": 20,
        "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
    }
    if gpu_available:
        scanvi_train_kwargs["accelerator"] = "gpu"
        scanvi_train_kwargs["devices"] = 1

    scanvi_model.train(**scanvi_train_kwargs)
    print("  -> scANVI complete")

    # --- Step 14: Export Results ---
    print("\n[Step 14] Exporting Results...")

    latent_scanvi = scanvi_model.get_latent_representation(adata_train)
    latent_df = pd.DataFrame(
        latent_scanvi,
        index=adata_train.obs_names,
        columns=[f"scANVI_{i}" for i in range(latent_scanvi.shape[1])]
    )
    latent_aligned = latent_df.reindex(adata_merged.obs_names)

    if latent_aligned.isna().any().any():
        raise ValueError("CRITICAL: Missing latent representation!")

    adata_merged.obsm["X_scANVI"] = latent_aligned.values

    pred_labels = scanvi_model.predict(adata_train)
    pred_df = pd.Series(pred_labels, index=adata_train.obs_names)
    pred_aligned = pred_df.reindex(adata_merged.obs_names)
    adata_merged.obs["scanvi_pred"] = pred_aligned.values

    proba_raw = scanvi_model.predict(adata_train, soft=True)
    proba = np.asarray(proba_raw, dtype=np.float32)
    label_order = scanvi_model.labels
    proba_df = pd.DataFrame(proba, index=adata_train.obs_names, columns=label_order)
    proba_aligned = proba_df.reindex(adata_merged.obs_names)
    adata_merged.obsm["scanvi_proba"] = proba_aligned.values
    adata_merged.obs["scanvi_confidence"] = proba_aligned.values.max(axis=1)

    adata_merged.uns["scanvi_label_order"] = list(label_order)

    print(f"  -> scANVI predictions:")
    print(adata_merged.obs["scanvi_pred"].value_counts().head(15))

    compute_novelty_scores(adata_merged, proba_aligned)

    # --- Step 15: Attach .raw ---
    print("\n[Step 15] Attaching .raw...")
    from anndata import AnnData

    adata_merged.raw = AnnData(X=full_counts, obs=adata_merged.obs.copy(), var=raw_var)
    print(f"  OK .raw: {adata_merged.raw.n_vars} genes")

    # --- Step 16: Multiple UMAPs (scVI + scANVI) ---
    print("\n[Step 16] Computing Multiple UMAPs...")

    print("  -> Getting scVI latent representation...")
    latent_scvi = scvi_model.get_latent_representation(adata_train)
    latent_scvi_df = pd.DataFrame(
        latent_scvi,
        index=adata_train.obs_names,
        columns=[f"scVI_{i}" for i in range(latent_scvi.shape[1])]
    )
    latent_scvi_aligned = latent_scvi_df.reindex(adata_merged.obs_names)
    if latent_scvi_aligned.isna().any().any():
        raise ValueError("CRITICAL: Missing scVI latent representation after reindex!")
    adata_merged.obsm["X_scVI"] = latent_scvi_aligned.values

    run_query_only_leiden(adata_merged, resolution=QUERY_LEIDEN_RESOLUTION)

    # 16a: UMAP on scVI latent
    print("  -> Computing UMAP on scVI latent...")
    sc.pp.neighbors(adata_merged, use_rep="X_scVI", n_neighbors=30, random_state=RANDOM_SEED,
                    key_added="neighbors_scVI")
    sc.tl.umap(adata_merged, random_state=RANDOM_SEED, neighbors_key="neighbors_scVI")
    adata_merged.obsm["X_umap_scVI"] = adata_merged.obsm["X_umap"].copy()
    print(f"     Saved to X_umap_scVI")

    # 16b: UMAP on scANVI latent (DEFAULT)
    print("  -> Computing UMAP on scANVI latent (DEFAULT)...")
    scanvi_umap_bundle = fit_bundle(
        adata_merged,
        "X_scANVI",
        output_dir,
        primary_umap_key="X_umap_scANVI",
        operator_filename=f"{OUTPUT_PREFIX}_umap_scanvi_operator.joblib",
        manifest_filename=f"{OUTPUT_PREFIX}_umap_scanvi_bundle.json",
        umap_params={
            "n_neighbors": 30,
            "n_components": 2,
            "min_dist": 0.5,
            "spread": 1.0,
            "metric": "euclidean",
            "random_state": RANDOM_SEED,
        },
        set_default_x_umap=True,
        extra_manifest={"pipeline": "Myeloid ONLY Reference-Query Merge Pipeline v1.0"},
    )
    print(f"     Saved to X_umap_scANVI and X_umap (default)")
    print(f"  -> UMAP operator (scANVI) saved: {scanvi_umap_bundle['operator_path']}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(adata_merged.obs["scanvi_confidence"], bins=50, edgecolor="black")
    ax.set_xlabel("Prediction Confidence")
    ax.set_ylabel("Cell Count")
    ax.set_title("scANVI Confidence Distribution")
    plt.tight_layout()
    plt.savefig(output_dir / f"{OUTPUT_PREFIX}_confidence_histogram.png", dpi=150)
    plt.close()

    # --- Step 17: Save Results ---
    print("\n[Step 17] Saving Results...")

    scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
    scvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scvi_model", overwrite=True)

    config = {
        "version": "1.0",
        "timestamp": datetime.now().isoformat(),
        "input": {"reference": REFERENCE_H5AD, "query": QUERY_H5AD},
        "n_hvg": int(n_hvg_final),
        "scanvi_labels": list(label_order),
        "myeloid_specific": {
            "classical_mono_markers": CLASSICAL_MONO_MARKERS,
            "nonclassical_mono_markers": NONCLASSICAL_MONO_MARKERS,
            "macrophage_markers": MACROPHAGE_MARKERS,
            "dc_markers": CDC1_MARKERS + CDC2_MARKERS,
            "neutrophil_markers": NEUTROPHIL_MARKERS,
            "forced_markers_count": len(FORCED_MARKERS)
        },
        "umap_spaces": {
            "X_umap": "DEFAULT - scANVI-based UMAP",
            "X_umap_scVI": "scVI latent space UMAP",
            "X_umap_scANVI": "scANVI latent space UMAP"
        }
    }

    with open(output_dir / f"{OUTPUT_PREFIX}_config.json", "w") as f:
        json.dump(config, f, indent=2)

    output_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"
    adata_merged.write_h5ad(output_h5ad, compression="gzip")
    print(f"  -> Saved: {output_h5ad}")

    train_h5ad = output_dir / f"{OUTPUT_PREFIX}_train_HVG.h5ad"
    adata_train.write_h5ad(train_h5ad, compression="gzip")

    # --- Step 18: Myeloid Specific Visualization ---
    print("\n[Step 18] Creating Myeloid visualizations...")

    fig = plt.figure(figsize=(24, 20))
    gs = fig.add_gridspec(5, 4, hspace=0.3, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, 0])
    sc.pl.umap(adata_merged, color="data_source", ax=ax1, show=False, title="Data Source (scANVI)", s=15)

    ax2 = fig.add_subplot(gs[0, 1])
    adata_merged.obs["_ref_coarse"] = pd.Series(pd.NA, index=adata_merged.obs_names, dtype="object")
    mask_ref = adata_merged.obs["data_source"] == "reference"
    adata_merged.obs.loc[mask_ref, "_ref_coarse"] = adata_merged.obs.loc[mask_ref, "cell_type_coarse"].astype(str).values
    adata_merged.obs["_ref_coarse"] = adata_merged.obs["_ref_coarse"].astype("category")
    sc.pl.umap(adata_merged, color="_ref_coarse", ax=ax2, show=False, title="Reference Coarse Labels", legend_loc="on data", s=15)

    ax3 = fig.add_subplot(gs[0, 2])
    sc.pl.umap(adata_merged, color="scanvi_pred", ax=ax3, show=False, title="scANVI Predictions", legend_loc="on data", s=15)

    ax4 = fig.add_subplot(gs[0, 3])
    sc.pl.umap(adata_merged, color="scanvi_confidence", ax=ax4, show=False, title="Confidence", cmap="viridis", vmin=0, vmax=1, s=15)

    # Row 2: Monocyte subtype scores
    ax5 = fig.add_subplot(gs[1, 0])
    sc.pl.umap(adata_merged, color="Classical_Mono_score", ax=ax5, show=False, title="Classical Mono Score (CD14)", cmap="Reds", s=15)

    ax6 = fig.add_subplot(gs[1, 1])
    sc.pl.umap(adata_merged, color="NonClassical_Mono_score", ax=ax6, show=False, title="NonClassical Mono Score (FCGR3A)", cmap="Blues", s=15)

    ax7 = fig.add_subplot(gs[1, 2])
    sc.pl.umap(adata_merged, color="mono_subtype_by_score", ax=ax7, show=False, title="Monocyte Subtype by Score", legend_loc="on data", s=15)

    ax8 = fig.add_subplot(gs[1, 3])
    query_mask_viz = adata_merged.obs["data_source"] == "query"
    ax8.scatter(
        adata_merged.obs.loc[query_mask_viz, "Classical_Mono_score"],
        adata_merged.obs.loc[query_mask_viz, "NonClassical_Mono_score"],
        c=adata_merged.obs.loc[query_mask_viz, "scanvi_confidence"],
        cmap="viridis", s=5, alpha=0.5
    )
    ax8.set_xlabel("Classical Mono Score")
    ax8.set_ylabel("NonClassical Mono Score")
    ax8.set_title(f"Query: Mono Subtype Scores\n(threshold={CLASSICAL_MONO_THRESHOLD})")
    ax8.axhline(y=NONCLASSICAL_MONO_THRESHOLD, color='k', linestyle='--', alpha=0.3)
    ax8.axvline(x=CLASSICAL_MONO_THRESHOLD, color='k', linestyle='--', alpha=0.3)

    # Row 3: Core Myeloid markers
    ax9 = fig.add_subplot(gs[2, 0])
    if "LYZ" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="LYZ", ax=ax9, show=False, title="LYZ (pan-Myeloid)", cmap="Reds", s=15, use_raw=True)

    ax10 = fig.add_subplot(gs[2, 1])
    if "CD14" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="CD14", ax=ax10, show=False, title="CD14 (Classical Mono)", cmap="Reds", s=15, use_raw=True)

    ax11 = fig.add_subplot(gs[2, 2])
    if "FCGR3A" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="FCGR3A", ax=ax11, show=False, title="FCGR3A/CD16 (NonClassical Mono)", cmap="Reds", s=15, use_raw=True)

    ax12 = fig.add_subplot(gs[2, 3])
    if "CD68" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="CD68", ax=ax12, show=False, title="CD68 (Macrophage)", cmap="Reds", s=15, use_raw=True)

    # Row 4: DC and lineage markers
    ax13 = fig.add_subplot(gs[3, 0])
    if "CD1C" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="CD1C", ax=ax13, show=False, title="CD1C (cDC2)", cmap="Reds", s=15, use_raw=True)

    ax14 = fig.add_subplot(gs[3, 1])
    if "CLEC9A" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="CLEC9A", ax=ax14, show=False, title="CLEC9A (cDC1)", cmap="Reds", s=15, use_raw=True)

    ax15 = fig.add_subplot(gs[3, 2])
    if "S100A8" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="S100A8", ax=ax15, show=False, title="S100A8 (Neutrophil)", cmap="Reds", s=15, use_raw=True)

    ax16 = fig.add_subplot(gs[3, 3])
    if "MKI67" in adata_merged.raw.var_names:
        sc.pl.umap(adata_merged, color="MKI67", ax=ax16, show=False, title="MKI67 (Proliferating)", cmap="Reds", s=15, use_raw=True)

    # Row 5: Query-specific analysis
    ax17 = fig.add_subplot(gs[4, 0])
    ref_dist = adata_merged.obs.loc[adata_merged.obs["data_source"] == "reference", "cell_type_coarse"].value_counts()
    qry_dist = adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "mono_subtype_by_score"].value_counts()
    x = np.arange(len(ref_dist.index))
    width = 0.35
    ax17.bar(x - width/2, ref_dist.values, width, label="Reference (annotated)", alpha=0.8)
    ax17.bar(x + width/2, [qry_dist.get(k, 0) for k in ref_dist.index], width, label="Query (by score)", alpha=0.8)
    ax17.set_xticks(x)
    ax17.set_xticklabels(ref_dist.index, rotation=45, ha="right")
    ax17.set_ylabel("Cell Count")
    ax17.set_title("Myeloid Subtype Distribution")
    ax17.legend()

    ax18 = fig.add_subplot(gs[4, 1])
    sc.pl.umap(adata_merged, color="novelty_score", ax=ax18, show=False,
               title="Novelty Score (Query)", cmap="hot", vmin=0, vmax=1, s=15)

    ax19 = fig.add_subplot(gs[4, 2])
    conf_ref = adata_merged.obs.loc[adata_merged.obs["data_source"] == "reference", "scanvi_confidence"]
    conf_qry = adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "scanvi_confidence"]
    ax19.hist([conf_ref, conf_qry], bins=30, label=["Reference", "Query"], alpha=0.7)
    ax19.set_xlabel("Confidence")
    ax19.set_ylabel("Cell Count")
    ax19.set_title("Confidence Distribution by Source")
    ax19.legend()

    ax20 = fig.add_subplot(gs[4, 3])
    sc.pl.umap(adata_merged, color="leiden_query", ax=ax20, show=False,
               title="Query-only Leiden Clusters", legend_loc="on data", s=15)

    plt.savefig(output_dir / f"{OUTPUT_PREFIX}_myeloid_overview.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  -> Saved: {OUTPUT_PREFIX}_myeloid_overview.pdf")

    # Figure 2: scVI vs scANVI UMAP comparison
    print("  -> Creating scVI vs scANVI UMAP comparison...")
    fig2, axes2 = plt.subplots(2, 3, figsize=(18, 12))

    sc.pl.umap(adata_merged, basis="X_umap_scVI", color="data_source", ax=axes2[0, 0], show=False,
               title="Data Source (scVI UMAP)", s=10)
    sc.pl.umap(adata_merged, basis="X_umap_scVI", color="_ref_coarse", ax=axes2[0, 1], show=False,
               title="Reference Labels (scVI)", legend_loc="on data", s=10)
    sc.pl.umap(adata_merged, basis="X_umap_scVI", color="scanvi_pred", ax=axes2[0, 2], show=False,
               title="scANVI Predictions (scVI)", legend_loc="on data", s=10)

    sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="data_source", ax=axes2[1, 0], show=False,
               title="Data Source (scANVI UMAP)", s=10)
    sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="_ref_coarse", ax=axes2[1, 1], show=False,
               title="Reference Labels (scANVI)", legend_loc="on data", s=10)
    sc.pl.umap(adata_merged, basis="X_umap_scANVI", color="scanvi_pred", ax=axes2[1, 2], show=False,
               title="scANVI Predictions (scANVI)", legend_loc="on data", s=10)

    plt.tight_layout()
    plt.savefig(output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  -> Saved: {OUTPUT_PREFIX}_umap_comparison.pdf")

    # Figure 3: Novelty analysis report
    print("  -> Creating novelty analysis report...")
    fig3, axes3 = plt.subplots(2, 2, figsize=(14, 12))

    qry_mask = adata_merged.obs["data_source"] == "query"

    ax = axes3[0, 0]
    ax.hist(adata_merged.obs.loc[qry_mask, "novelty_score"], bins=50, edgecolor="black", alpha=0.7)
    ax.axvline(x=0.7, color='r', linestyle='--', label='High novelty threshold')
    ax.set_xlabel("Novelty Score")
    ax.set_ylabel("Cell Count")
    ax.set_title("Query Cells: Novelty Score Distribution")
    ax.legend()

    ax = axes3[0, 1]
    scatter = ax.scatter(
        adata_merged.obs.loc[qry_mask, "scanvi_confidence"],
        adata_merged.obs.loc[qry_mask, "novelty_score"],
        c=adata_merged.obs.loc[qry_mask, "scanvi_entropy"],
        cmap="viridis", s=5, alpha=0.5
    )
    ax.set_xlabel("scANVI Confidence")
    ax.set_ylabel("Novelty Score")
    ax.set_title("Query: Confidence vs Novelty (colored by entropy)")
    plt.colorbar(scatter, ax=ax)

    ax = axes3[1, 0]
    leiden_counts = adata_merged.obs.loc[qry_mask, "leiden_query"].value_counts().head(15)
    ax.barh(range(len(leiden_counts)), leiden_counts.values)
    ax.set_yticks(range(len(leiden_counts)))
    ax.set_yticklabels(leiden_counts.index)
    ax.set_xlabel("Cell Count")
    ax.set_title("Query-only Leiden Clusters (Top 15)")

    ax = axes3[1, 1]
    # Macrophage vs Neutrophil score comparison
    ax.scatter(
        adata_merged.obs.loc[qry_mask, "Macrophage_score"],
        adata_merged.obs.loc[qry_mask, "Neutrophil_score"],
        c=adata_merged.obs.loc[qry_mask, "scanvi_confidence"],
        cmap="viridis", s=5, alpha=0.5
    )
    ax.set_xlabel("Macrophage Score")
    ax.set_ylabel("Neutrophil Score")
    ax.set_title("Query: Macrophage vs Neutrophil Score")

    plt.tight_layout()
    plt.savefig(output_dir / f"{OUTPUT_PREFIX}_novelty_analysis.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  -> Saved: {OUTPUT_PREFIX}_novelty_analysis.pdf")

    adata_merged.obs.drop(columns=["_ref_coarse"], inplace=True, errors="ignore")

    # --- Summary ---
    print("\n" + "=" * 80)
    print("MYELOID PIPELINE COMPLETE (v1.0)")
    print("=" * 80)
    print(f"\nOutput: {output_dir}")
    print(f"\nKey observations:")
    print(f"  - Total cells: {adata_merged.n_obs:,}")
    print(f"  - Reference: {(adata_merged.obs['data_source'] == 'reference').sum():,}")
    print(f"  - Query: {(adata_merged.obs['data_source'] == 'query').sum():,}")

    print(f"\nMonocyte Subtype Distribution in Query (by score):")
    print(adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "mono_subtype_by_score"].value_counts())

    print(f"\nTop scANVI predictions in Query:")
    print(adata_merged.obs.loc[adata_merged.obs["data_source"] == "query", "scanvi_pred"].value_counts().head(10))

    print(f"\nNovelty Detection Summary:")
    print(f"  - High novelty cells (score > 0.7): {adata_merged.obs['is_potentially_novel'].sum():,}")
    print(f"  - Query-only Leiden clusters: {adata_merged.obs.loc[qry_mask, 'leiden_query'].nunique()}")

    print(f"\nUMAP spaces available:")
    print(f"  - X_umap: DEFAULT (scANVI-based)")
    print(f"  - X_umap_scVI: scVI latent space")
    print(f"  - X_umap_scANVI: scANVI latent space")

    print("\n" + "=" * 80)
    print("Myeloid Validation Strategy:")
    print("=" * 80)
    print("1. Compare 'scanvi_pred' (from reference labels) with 'mono_subtype_by_score'")
    print(f"   - Score thresholds: Classical >= {CLASSICAL_MONO_THRESHOLD}, NonClassical >= {NONCLASSICAL_MONO_THRESHOLD}")
    print("2. High-novelty query cells (novelty_score > 0.7) may represent new subtypes")
    print("3. Check 'leiden_query' for query-only clusters not present in reference")
    print("4. Use Classical vs NonClassical Mono scores to identify monocyte subtypes")
    print("5. Check core Myeloid markers (LYZ, CD14) to filter out non-Myeloid contaminants")
    print("6. Compare Macrophage_score vs Neutrophil_score to identify lineage")
    print("=" * 80)


if __name__ == "__main__":
    main()
