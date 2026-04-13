#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stromal branch-wise reintegration / reference training v1.5
===========================================================

Purpose
-------
Build three separate stromal reference models after the v1.3 table-driven cleanup:
  - endothelial
  - fibroblast
  - smc   (Smooth_Muscle + Pericyte)

Schwann is intentionally excluded from model training.

Compared with `stromal_reintegration_scvi_scanvi_20260312_v1_4.ipynb`
-------------------------------------------------------------------
1. Keep the same table-driven contamination removal + L1/L2/L3 assignment.
2. After cleanup, split cells into three biologically coherent branches.
3. Train branch-specific scVI/scANVI references independently.
4. Save one branch reference h5ad + model package per branch.
5. Save a root manifest JSON so downstream branch-wise scArches/scHPL can load
   the correct reference automatically.

The script also supports `--smoke-test`, which synthesizes a tiny toy dataset and
runs the full branch loop with 1 epoch to validate code paths end-to-end.
"""

from __future__ import annotations

import argparse
import gc
import importlib
import json
import os
import random
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

for _k in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ.setdefault(_k, "8")

import matplotlib
matplotlib.use("Agg")
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
import matplotlib.pyplot as plt

from anndata import AnnData
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

import sympy

try:
    if not hasattr(sympy, "printing"):
        sympy.printing = importlib.import_module("sympy.printing")
except Exception as e:  # pragma: no cover - import guard
    raise RuntimeError(
        f"Sympy import is broken (sympy module path: {getattr(sympy, '__file__', 'unknown')}). "
        "Please ensure real sympy package is available and no local sympy.py shadows it."
    ) from e


PIPELINE_VERSION = "v1.5-branchwise"
PIPELINE_DATE = "2026-04-07"
PIPELINE_TAG = "v1_5_branchwise"
UNLABELED = "Unknown"
RANDOM_SEED = 42

SRC_L2_COL = "cell_type_L2"
SRC_L3_COL_DEFAULT = "cell_type_L3"
BATCH_KEY = "sample"
OUT_L1 = "cell_type_L1"
OUT_L2 = "cell_type_L2"
OUT_L3 = "cell_type_L3"
TISSUE_KEY = "tissue"
LUNG_TRACHEA_VALUES = {
    "lung", "trachea", "Lung", "Trachea",
    "Lung tissue", "Tracheal tissue",
    "lung_tissue", "trachea_tissue",
}

ANNOTATION_TABLE = [
    ("Endothelia_Lymphatic_c0", "Endothelial", "Endothelia_Lymphatic", "KEEP"),
    ("Endothelia_Lymphatic_c1", "Endothelial", "Endothelia_Lymphatic", "REVIEW"),
    ("Endothelia_Lymphatic_c2", "Endothelial", "Endothelia_Lymphatic", "KEEP"),
    ("Endothelia_vascular_Cap_a_c0", "Contaminants", "Contaminants", "DROP"),
    ("Endothelia_vascular_Cap_a_c1", "Endothelial", "Endothelia_vascular_Cap_a", "KEEP"),
    ("Endothelia_vascular_Cap_a_c2", "Endothelial", "Endothelia_vascular_Cap_a", "KEEP"),
    ("Endothelia_vascular_Cap_g_c0", "Endothelial", "Endothelia_vascular_Cap_g", "REVIEW"),
    ("Endothelia_vascular_Cap_g_c1", "Endothelial", "Endothelia_vascular_Cap_g", "KEEP"),
    ("Endothelia_vascular_Cap_g_c2", "Endothelial", "Endothelia_vascular_Cap_g", "REVIEW"),
    ("Endothelia_vascular_Cap_g_c3", "Endothelial", "Endothelia_vascular_Cap_g", "REVIEW"),
    ("Endothelia_vascular_arterial_pulmonary_c0", "Endothelial", "Endothelia_vascular_arterial_pulmonary", "KEEP"),
    ("Endothelia_vascular_arterial_pulmonary_c1", "Endothelial", "Endothelia_vascular_arterial_pulmonary", "KEEP"),
    ("Endothelia_vascular_arterial_pulmonary_c2", "Endothelial", "Endothelia_vascular_arterial_pulmonary", "REVIEW"),
    ("Endothelia_vascular_arterial_systemic_c0", "Endothelial", "Endothelia_vascular_arterial_systemic", "KEEP"),
    ("Endothelia_vascular_arterial_systemic_c1", "Endothelial", "Endothelia_vascular_arterial_systemic", "REVIEW"),
    ("Endothelia_vascular_venous_pulmonary_c0", "Endothelial", "Endothelia_vascular_venous_pulmonary", "KEEP"),
    ("Endothelia_vascular_venous_pulmonary_c1", "Endothelial", "Endothelia_vascular_venous_pulmonary", "REVIEW"),
    ("Endothelia_vascular_venous_pulmonary_c2", "Endothelial", "Endothelia_vascular_venous_pulmonary", "REVIEW"),
    ("Endothelia_vascular_venous_systemic_c0", "Endothelial", "Endothelia_vascular_venous_systemic", "REVIEW"),
    ("Endothelia_vascular_venous_systemic_c1", "Endothelial", "Endothelia_vascular_venous_systemic", "KEEP"),
    ("Endothelia_vascular_venous_systemic_c2", "Endothelial", "Endothelia_vascular_venous_systemic", "KEEP"),
    ("Endothelia_vascular_venous_systemic_c3", "Contaminants", "Contaminants", "DROP"),
    ("Endothelia_vascular_venous_systemic_c4", "Fibroblast", "Fibro_stress_activated", "REASSIGN"),
    ("Endothelia_vascular_venous_systemic_c5", "Contaminants", "Contaminants", "DROP"),
    ("Fibro_adventitial_c0", "Fibroblast", "Fibro_adventitial", "KEEP"),
    ("Fibro_adventitial_c1", "Fibroblast", "Fibro_adventitial", "KEEP"),
    ("Fibro_adventitial_c2", "Fibroblast", "Fibro_alveolar", "REASSIGN"),
    ("Fibro_alveolar_c0", "Fibroblast", "Fibro_alveolar", "REVIEW"),
    ("Fibro_alveolar_c1", "Fibroblast", "Fibro_alveolar", "KEEP"),
    ("Fibro_alveolar_c2", "Fibroblast", "Fibro_alveolar", "KEEP"),
    ("Fibro_myofibroblast_c0", "Contaminants", "Contaminants", "DROP"),
    ("Fibro_myofibroblast_c1", "Fibroblast", "Fibro_myofibroblast", "REVIEW"),
    ("Fibro_peribronchial_c0", "Fibroblast", "Fibro_peribronchial", "KEEP"),
    ("Fibro_peribronchial_c1", "Fibroblast", "Fibro_peribronchial", "KEEP"),
    ("Fibro_peribronchial_c2", "Fibroblast", "Fibro_peribronchial", "REVIEW"),
    ("Muscle_pericyte_pulmonary_c0", "Pericyte", "Muscle_pericyte_pulmonary", "KEEP"),
    ("Muscle_pericyte_pulmonary_c1", "Pericyte", "Muscle_pericyte_pulmonary", "REVIEW"),
    ("Muscle_pericyte_pulmonary_c2", "Pericyte", "Muscle_pericyte_pulmonary", "KEEP"),
    ("Muscle_pericyte_pulmonary_c3", "Pericyte", "Muscle_pericyte_pulmonary", "KEEP"),
    ("Muscle_pericyte_systemic_c0", "Pericyte", "Muscle_pericyte_systemic", "KEEP"),
    ("Muscle_pericyte_systemic_c1", "Pericyte", "Muscle_pericyte_systemic", "KEEP"),
    ("Muscle_perivascular_immune_recruiting_c0", "Smooth_Muscle", "Muscle_perivascular_immune_recruiting", "KEEP"),
    ("Muscle_perivascular_immune_recruiting_c1", "Pericyte", "Muscle_pericyte_systemic", "REASSIGN"),
    ("Muscle_smooth_arterial_systemic_c0", "Smooth_Muscle", "Muscle_smooth_arterial_systemic", "KEEP"),
    ("Muscle_smooth_arterial_systemic_c1", "Smooth_Muscle", "Muscle_smooth_arterial_systemic", "KEEP"),
    ("Muscle_smooth_pulmonary_c0", "Smooth_Muscle", "Muscle_smooth_pulmonary", "REVIEW"),
    ("Muscle_smooth_pulmonary_c1", "Smooth_Muscle", "Muscle_smooth_pulmonary", "KEEP"),
    ("Schwann_nonmyelinating_c0", "Schwann", "Schwann_nonmyelinating", "KEEP"),
]

ANNOT_DF = pd.DataFrame(
    ANNOTATION_TABLE,
    columns=["cluster", "L2", "L3", "action"],
).set_index("cluster")
DROP_CLUSTERS = ANNOT_DF[ANNOT_DF["action"] == "DROP"].index.tolist()
REASSIGN_CLUSTERS = ANNOT_DF[ANNOT_DF["action"] == "REASSIGN"].index.tolist()
PULMONARY_ALVEOLAR_L3 = sorted(
    ANNOT_DF[ANNOT_DF["L3"].str.contains("pulmonary|alveolar", case=False, regex=True)]["L3"]
    .unique()
    .tolist()
)

BRANCH_CONFIGS: Dict[str, Dict[str, object]] = {
    "endothelial": {
        "display": "Endothelial",
        "allowed_l2": ("Endothelial",),
        "query_subset": "endothelial",
    },
    "fibroblast": {
        "display": "Fibroblast",
        "allowed_l2": ("Fibroblast",),
        "query_subset": "fibroblast",
    },
    "smc": {
        "display": "Smooth_Muscle",
        "allowed_l2": ("Smooth_Muscle", "Pericyte"),
        "query_subset": "smc",
    },
}
IGNORED_L2 = {"Schwann"}
CONTINUOUS_COVARIATES = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]

STRESS_SIGNATURE_GENES = [
    "ALDH18A1","ARFGAP1","ASNS","ATF3","ATF4","ATF6","ATP6V0D1","BAG3","BANF1",
    "CALR","CCL2","CEBPB","CEBPG","CHAC1","CKS1B","CNOT2","CNOT4","CNOT6",
    "CXXC1","DCP1A","DCP2","DCTN1","DDIT4","DDX10","DKC1","DNAJA4","DNAJB9",
    "DNAJC3","EDC4","EDEM1","EEF2","EIF2AK3","EIF2S1","EIF4A1","EIF4A2","EIF4A3",
    "EIF4E","EIF4EBP1","EIF4G1","ERN1","ERO1A","EXOC2","EXOSC1","EXOSC10",
    "EXOSC2","EXOSC4","EXOSC5","EXOSC9","FKBP14","FUS","GEMIN4","GOSR2","H2AX",
    "HERPUD1","HSP90B1","HSPA5","HSPA9","HYOU1","IARS1","IFIT1","IGFBP1","IMP3",
    "KDELR3","KHSRP","KIF5B","LSM1","LSM4","MTHFD2","NFYA","NFYB","NHP2","NOLC1",
    "NOP14","NOP56","NPM1","NABP1","PAIP1","PARN","PDIA5","PDIA6","POP4","PREB",
    "PSAT1","RPS14","RRP9","SDAD1","SEC11A","SEC31A","SERP1","SHC1","MTREX",
    "SLC1A4","SLC30A5","SLC7A5","SPCS1","SPCS3","SRPRA","SRPRB","SSR1","STC2",
    "TARS1","TATDN2","TSPYL2","SKIC3","TUBB2A","VEGFA","WFS1","WIPI1","XBP1",
    "XPOT","YIF1A","YWHAZ","ZBTB17",
]
S_GENES = [
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
    "CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP",
    "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
    "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
    "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8",
]
G2M_GENES = [
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
    "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
    "HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
    "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
    "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA",
]


@dataclass
class PipelineConfig:
    input_h5ad: Path
    output_dir: Path
    n_hvg: int = 4000
    n_latent: int = 75
    n_hidden: int = 128
    n_layers: int = 2
    dropout: float = 0.2
    max_epochs_scvi: int = 400
    max_epochs_scanvi: int = 200
    batch_size: int = 256
    min_cells_batch: int = 3
    purity_k: int = 30
    purity_threshold: float = 0.5
    tissue_key: str = TISSUE_KEY
    source_l3_col: str = SRC_L3_COL_DEFAULT
    use_gpu: bool = torch.cuda.is_available()
    skip_figures: bool = False
    smoke_test: bool = False
    continuous_covariates: Tuple[str, ...] = tuple(CONTINUOUS_COVARIATES)
    encode_covariates: bool = True
    use_layer_norm: str = "both"
    use_batch_norm: str = "none"
    scvi_lr: float = 1e-3
    scanvi_lr: float = 5e-4
    scvi_patience: int = 45
    scanvi_patience: int = 30


def set_seed(seed: int = RANDOM_SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    scvi.settings.seed = seed
    scvi.settings.dl_num_workers = 0
    sc.settings.verbosity = 2
    sc.settings.n_jobs = 16


def log_header(title: str) -> None:
    print("\n" + "=" * 80)
    print(title)
    print("=" * 80)


def _annotation_match_count(values: Iterable[str]) -> int:
    return pd.Index(values).isin(ANNOT_DF.index).sum()


def autodetect_source_l3_col(adata: AnnData, configured_col: str) -> str:
    l3_uniq = adata.obs[configured_col].astype(str).unique()
    n_match = _annotation_match_count(l3_uniq)
    best_col = configured_col

    candidate_cols: List[str] = []
    if "subcluster_id" in adata.obs.columns:
        candidate_cols.append("subcluster_id")
    candidate_cols.extend([
        c for c in adata.obs.columns
        if c not in (configured_col, SRC_L2_COL) and ("l3" in c.lower() or "subcluster" in c.lower())
    ])

    for c in dict.fromkeys(candidate_cols):
        vals = adata.obs[c].astype(str).unique()
        m = _annotation_match_count(vals)
        if m > n_match:
            best_col = c
            l3_uniq = vals
            n_match = m

    if n_match == 0 and SRC_L2_COL in adata.obs.columns and "subcluster_id" in adata.obs.columns:
        reconstructed_col = "__reconstructed_cluster_id"
        reconstructed = (
            adata.obs[SRC_L2_COL].astype(str) + "_c" + adata.obs["subcluster_id"].astype(str)
        )
        adata.obs[reconstructed_col] = pd.Categorical(reconstructed)
        vals = adata.obs[reconstructed_col].astype(str).unique()
        m = _annotation_match_count(vals)
        if m > n_match:
            print(
                f"[INFO] Reconstructed cluster IDs from '{SRC_L2_COL}' + 'subcluster_id' -> "
                f"'{reconstructed_col}' (coverage {m}/{len(vals)})"
            )
            best_col = reconstructed_col
            n_match = m

    if n_match == 0:
        raise ValueError(
            f"[ERROR] No values in '{configured_col}' match ANNOTATION_TABLE. "
            "Check cluster ID format or update the source L3 column."
        )

    if best_col != configured_col:
        print(
            f"[INFO] Auto-switch source L3 column: '{configured_col}' -> '{best_col}' "
            f"(coverage {n_match}/{len(adata.obs[best_col].astype(str).unique())})"
        )

    return best_col


def load_and_clean_input(cfg: PipelineConfig) -> Tuple[AnnData, str, bool, int]:
    log_header("STEP 1-5: LOAD / TABLE ASSIGN / DROP CONTAMINATION / LAYER CHECK / SMALL BATCH FILTER")
    adata = sc.read_h5ad(cfg.input_h5ad)
    n_before = adata.n_obs
    print(f"Loaded input: {cfg.input_h5ad}")
    print(f"Shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    for col in [SRC_L2_COL, BATCH_KEY]:
        if col not in adata.obs.columns:
            raise KeyError(f"[ERROR] Missing required obs column: {col}")
    if cfg.source_l3_col not in adata.obs.columns:
        raise KeyError(f"[ERROR] Missing configured L3 column: {cfg.source_l3_col}")

    src_l3_col = autodetect_source_l3_col(adata, cfg.source_l3_col)
    tissue_key_available = cfg.tissue_key in adata.obs.columns
    if not tissue_key_available:
        similar = [c for c in adata.obs.columns if any(k in c.lower() for k in ("tissue", "organ", "site"))]
        print(f"[WARN] TISSUE_KEY='{cfg.tissue_key}' not found. Candidates: {similar}")
        print("       Tissue-aware Unknown assignment will be skipped.")

    source_cluster_ids = adata.obs[src_l3_col].astype(str).copy()
    if "source_cluster_id" not in adata.obs.columns:
        adata.obs["source_cluster_id"] = pd.Categorical(source_cluster_ids)

    l2_mapped = source_cluster_ids.map(ANNOT_DF["L2"])
    l3_mapped = source_cluster_ids.map(ANNOT_DF["L3"])
    action_mapped = source_cluster_ids.map(ANNOT_DF["action"])
    n_mapped = l2_mapped.notna().sum()
    print(f"Annotation-mapped cells: {n_mapped:,}/{adata.n_obs:,}")
    if n_mapped == 0:
        raise ValueError("[ERROR] Annotation mapping produced 0 matched cells.")

    adata.obs[OUT_L1] = np.where(l2_mapped.fillna("Unknown") == "Contaminants", "Contaminants", "Stromal")
    adata.obs[OUT_L2] = l2_mapped.fillna("Unresolved")
    adata.obs[OUT_L3] = l3_mapped.fillna("Unresolved")
    adata.obs["annot_action"] = action_mapped.fillna("REVIEW")

    drop_mask = source_cluster_ids.isin(DROP_CLUSTERS)
    adata = adata[~drop_mask].copy()
    print(f"Contamination removal: {n_before:,} -> {adata.n_obs:,} cells")

    if "counts" not in adata.layers:
        if adata.raw is None:
            raise ValueError("[ERROR] No 'counts' layer and no .raw -- cannot continue.")
        print("[WARN] 'counts' layer missing -- inferring from .raw.X")
        missing = adata.var_names.difference(adata.raw.var_names)
        if len(missing) > 0:
            raise ValueError(f".raw is missing {len(missing)} current genes; cannot recover counts safely.")
        adata.layers["counts"] = sparse.csr_matrix(adata.raw[:, adata.var_names].X).astype(np.float32)

    if "log1p" not in adata.layers:
        adata.layers["log1p"] = adata.X.copy()

    adata.X = adata.layers["log1p"]
    for lk in ["counts", "log1p"]:
        if lk in adata.layers and not sparse.issparse(adata.layers[lk]):
            adata.layers[lk] = sparse.csr_matrix(adata.layers[lk])
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)

    bc = adata.obs[BATCH_KEY].value_counts()
    small = bc[bc < cfg.min_cells_batch].index.tolist()
    if small:
        n_pre = adata.n_obs
        adata = adata[~adata.obs[BATCH_KEY].isin(small)].copy()
        print(f"Removed {len(small)} small batches: {n_pre:,} -> {adata.n_obs:,}")
    else:
        print(f"[OK] All {adata.obs[BATCH_KEY].nunique()} batches pass minimum size")

    gc.collect()
    return adata, src_l3_col, tissue_key_available, n_before


def subset_branch(adata: AnnData, branch_name: str, branch_cfg: Dict[str, object]) -> AnnData:
    allowed_l2 = set(branch_cfg["allowed_l2"])
    mask = adata.obs[OUT_L2].astype(str).isin(allowed_l2)
    ad = adata[mask].copy()
    ad.obs["reference_branch"] = branch_name
    ad.obs["reference_branch_display"] = branch_cfg["display"]
    return ad


def compute_branch_covariates(adata: AnnData, branch_name: str) -> None:
    print(f"Computing covariates for branch '{branch_name}' on full gene space...")
    if "counts" not in adata.layers:
        raise KeyError("[ERROR] 'counts' layer required for covariate calculation")

    counts = adata.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
        adata.layers["counts"] = counts
    else:
        counts = counts.tocsr()
        adata.layers["counts"] = counts

    gene_names = pd.Index(adata.var_names.astype(str))

    mt_mask = gene_names.str.upper().str.startswith("MT-") | gene_names.str.upper().str.startswith("MT.")
    total_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    if mt_mask.any():
        mt_counts = np.asarray(counts[:, np.where(mt_mask)[0]].sum(axis=1)).ravel().astype(np.float64)
        pct_counts_mt = np.zeros_like(total_counts, dtype=np.float32)
        nz = total_counts > 0
        pct_counts_mt[nz] = (mt_counts[nz] / total_counts[nz] * 100.0).astype(np.float32)
        print(f"  pct_counts_mt: n_genes={int(mt_mask.sum())} mean={pct_counts_mt.mean():.2f}%")
    else:
        pct_counts_mt = np.zeros(adata.n_obs, dtype=np.float32)
        print("  [WARN] No MT genes found; pct_counts_mt set to 0")
    adata.obs["pct_counts_mt"] = pct_counts_mt

    stress_genes_present = [g for g in STRESS_SIGNATURE_GENES if g in gene_names]
    print(f"  Stress genes: {len(stress_genes_present)}/{len(STRESS_SIGNATURE_GENES)} found")
    if len(stress_genes_present) >= 10:
        adata_tmp = AnnData(X=counts.copy(), var=adata.var.copy())
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        stress_idx = gene_names.get_indexer(stress_genes_present)
        stress_expr = adata_tmp.X[:, stress_idx]
        if sparse.issparse(stress_expr):
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
        else:
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
        del stress_expr, adata_tmp
        print(f"  stress_score: mean={float(stress_score.mean()):.4f} max={float(stress_score.max()):.4f}")
    else:
        stress_score = np.zeros(adata.n_obs, dtype=np.float32)
        print("  [WARN] Too few stress genes; stress_score set to 0")
    adata.obs["stress_score"] = stress_score

    s_genes_present = [g for g in S_GENES if g in gene_names]
    g2m_genes_present = [g for g in G2M_GENES if g in gene_names]
    print(
        f"  Cell cycle genes: S={len(s_genes_present)}/{len(S_GENES)} "
        f"G2M={len(g2m_genes_present)}/{len(G2M_GENES)}"
    )
    if len(s_genes_present) >= 10 and len(g2m_genes_present) >= 10:
        adata_tmp = AnnData(X=counts.copy(), obs=adata.obs[[BATCH_KEY]].copy(), var=adata.var.copy())
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)
        sc.tl.score_genes_cell_cycle(adata_tmp, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
        adata.obs["S_score"] = adata_tmp.obs["S_score"].to_numpy(dtype=np.float32)
        adata.obs["G2M_score"] = adata_tmp.obs["G2M_score"].to_numpy(dtype=np.float32)
        adata.obs["phase"] = adata_tmp.obs["phase"].astype(str).to_numpy()
        print(f"  S_score: mean={adata.obs['S_score'].mean():.4f}")
        print(f"  G2M_score: mean={adata.obs['G2M_score'].mean():.4f}")
        del adata_tmp
    else:
        adata.obs["S_score"] = np.zeros(adata.n_obs, dtype=np.float32)
        adata.obs["G2M_score"] = np.zeros(adata.n_obs, dtype=np.float32)
        adata.obs["phase"] = np.repeat("G1", adata.n_obs)
        print("  [WARN] Too few cell-cycle genes; S_score/G2M_score set to 0")

    for cov in CONTINUOUS_COVARIATES:
        if cov not in adata.obs.columns:
            raise KeyError(f"[ERROR] Missing computed covariate: {cov}")

    gc.collect()


def select_hvgs_and_preserve_raw(adata: AnnData, cfg: PipelineConfig, branch_dir: Path) -> str:
    print(f"Selecting up to {cfg.n_hvg} HVGs from {adata.n_vars:,} genes...")
    n_top = min(cfg.n_hvg, adata.n_vars)
    try:
        sc.pp.highly_variable_genes(
            adata,
            layer="counts",
            n_top_genes=n_top,
            batch_key=BATCH_KEY,
            flavor="seurat_v3",
            subset=False,
        )
        hvg_method = "batch-aware seurat_v3"
    except Exception as e:
        print(f"  [WARN] {e}")
        try:
            sc.pp.highly_variable_genes(
                adata,
                layer="counts",
                n_top_genes=n_top,
                flavor="seurat_v3",
                subset=False,
            )
            hvg_method = "seurat_v3 (no batch)"
        except Exception as e2:
            print(f"  [WARN] {e2}")
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, subset=False)
            hvg_method = "default"

    adata.raw = sc.AnnData(X=adata.layers["counts"], obs=adata.obs.copy(), var=adata.var.copy())
    adata._inplace_subset_var(adata.var["highly_variable"].to_numpy())
    if "counts" not in adata.layers:
        adata.layers["counts"] = sparse.csr_matrix(adata.raw[:, adata.var_names].X).astype(np.float32)
    if "log1p" not in adata.layers:
        adata.layers["log1p"] = adata.X.copy()
    adata.X = adata.layers["log1p"]

    pd.Series(adata.var_names.tolist()).to_csv(branch_dir / "hvg_genes.csv", index=False, header=False)
    return hvg_method


def train_scvi(adata: AnnData, cfg: PipelineConfig, model_dir: Path) -> scvi.model.SCVI:
    setup_kwargs = dict(layer="counts", batch_key=BATCH_KEY)
    covariates_present = [cov for cov in cfg.continuous_covariates if cov in adata.obs.columns]
    if covariates_present:
        setup_kwargs["continuous_covariate_keys"] = covariates_present

    scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
    model_scvi = scvi.model.SCVI(
        adata,
        n_latent=cfg.n_latent,
        n_hidden=cfg.n_hidden,
        n_layers=cfg.n_layers,
        dropout_rate=cfg.dropout,
        gene_likelihood="nb",
        dispersion="gene-batch",
        encode_covariates=cfg.encode_covariates,
        use_layer_norm=cfg.use_layer_norm,
        use_batch_norm=cfg.use_batch_norm,
    )
    accelerator = "gpu" if cfg.use_gpu and torch.cuda.is_available() else "cpu"
    print(
        f"Training scVI with accelerator={accelerator}, n_latent={cfg.n_latent}, "
        f"covariates={covariates_present}"
    )
    t0 = time.time()
    model_scvi.train(
        max_epochs=cfg.max_epochs_scvi,
        batch_size=cfg.batch_size,
        early_stopping=True,
        early_stopping_patience=cfg.scvi_patience,
        train_size=0.9,
        accelerator=accelerator,
        devices=1,
        plan_kwargs={"lr": cfg.scvi_lr},
    )
    print(f"[OK] scVI done in {(time.time() - t0) / 60:.1f} min")
    model_scvi.save(str(model_dir), overwrite=True)
    adata.obsm["X_scvi"] = model_scvi.get_latent_representation()
    return model_scvi


def build_scanvi_labels(adata: AnnData, cfg: PipelineConfig, tissue_key_available: bool) -> Tuple[pd.Series, int, int]:
    scanvi_labels = adata.obs[OUT_L3].astype(str).copy()

    if tissue_key_available:
        tissue_str = adata.obs[cfg.tissue_key].astype(str)
        tissue_norm = (
            tissue_str.str.strip().str.lower()
            .str.replace(r"[_-]+", " ", regex=True)
            .str.replace(r"\s+", " ", regex=True)
        )
        lung_trachea_norm = {x.strip().lower().replace("_", " ").replace("-", " ") for x in LUNG_TRACHEA_VALUES}
        exact_tissue_match = tissue_norm.isin(lung_trachea_norm)
        keyword_tissue_match = tissue_norm.str.contains(
            r"\blung\b|\btrachea\b|\bairway\b|\bbronch|\bparenchyma\b|\bpulmon",
            regex=True,
            na=False,
        )
        is_lung_or_trachea = exact_tissue_match | keyword_tissue_match
        is_pulm_alv = scanvi_labels.isin(PULMONARY_ALVEOLAR_L3)
        tissue_mismatch = is_pulm_alv & ~is_lung_or_trachea
        scanvi_labels[tissue_mismatch] = UNLABELED
        print(f"Tissue mismatch -> Unknown: {int(tissue_mismatch.sum()):,}")
    else:
        print("[WARN] Tissue key unavailable -- skipping tissue-aware filtering.")

    print(f"Computing kNN purity (k={cfg.purity_k}) in scVI latent space...")
    knn = NearestNeighbors(n_neighbors=min(cfg.purity_k, max(2, adata.n_obs)), algorithm="auto", n_jobs=16)
    knn.fit(adata.obsm["X_scvi"])
    _, knn_idx = knn.kneighbors(adata.obsm["X_scvi"])
    labels_arr = scanvi_labels.values
    neighbor_labels = labels_arr[knn_idx]
    purity = np.mean(neighbor_labels == labels_arr[:, None], axis=1).astype(np.float32)
    is_unknown_mask = labels_arr == UNLABELED
    purity[is_unknown_mask] = 1.0
    adata.obs["label_purity_scanvi"] = purity
    low_purity = (purity < cfg.purity_threshold) & ~is_unknown_mask
    scanvi_labels[low_purity] = UNLABELED
    print(f"Low purity -> Unknown: {int(low_purity.sum()):,} ({low_purity.sum() / adata.n_obs * 100:.1f}%)")
    del knn, knn_idx, neighbor_labels
    gc.collect()

    adata.obs["scanvi_label"] = scanvi_labels.values
    n_labeled = int((scanvi_labels != UNLABELED).sum())
    n_unknown = int((scanvi_labels == UNLABELED).sum())
    return scanvi_labels, n_labeled, n_unknown


def train_scanvi(
    adata: AnnData,
    model_scvi: scvi.model.SCVI,
    cfg: PipelineConfig,
    scvi_model_dir: Path,
    model_dir: Path,
) -> bool:
    labels = adata.obs["scanvi_label"].astype(str)
    labeled_mask = labels != UNLABELED
    n_labeled_classes = labels[labeled_mask].nunique()
    print(f"Labeled classes (excluding '{UNLABELED}'): {n_labeled_classes}")

    if n_labeled_classes < 2:
        print("[WARN] <2 labeled classes; skipping scANVI and using scVI fallback outputs.")
        adata.obsm["X_scanvi"] = adata.obsm["X_scvi"].copy()
        adata.obs["cell_type_scanvi_pred"] = labels.copy()
        adata.obs["scanvi_uncertainty"] = np.where(labels == UNLABELED, 1.0, 0.0).astype(np.float32)
        return False

    def build_scanvi_from_scvi(scvi_model: scvi.model.SCVI) -> scvi.model.SCANVI:
        return scvi.model.SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category=UNLABELED,
            labels_key="scanvi_label",
        )

    def cleanup_cuda_cache() -> None:
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    def is_cuda_training_error(exc: Exception) -> bool:
        msg = str(exc).lower()
        return any(token in msg for token in ["cuda", "cudnn", "nccl", "device-side", "acceleratorerror"])

    accelerator = "gpu" if cfg.use_gpu and torch.cuda.is_available() else "cpu"
    train_kwargs = dict(
        max_epochs=cfg.max_epochs_scanvi,
        batch_size=cfg.batch_size,
        early_stopping=True,
        early_stopping_patience=cfg.scanvi_patience,
        train_size=0.9,
        devices=1,
        plan_kwargs={"lr": cfg.scanvi_lr, "weight_decay": 0.0},
    )

    model_scanvi = build_scanvi_from_scvi(model_scvi)
    import inspect
    if "n_samples_per_label" in inspect.signature(model_scanvi.train).parameters:
        train_kwargs["n_samples_per_label"] = 100

    try:
        t0 = time.time()
        model_scanvi.train(accelerator=accelerator, **train_kwargs)
        print(f"[OK] scANVI done in {(time.time() - t0) / 60:.1f} min")
    except Exception as exc:
        if accelerator != "gpu" or not is_cuda_training_error(exc):
            raise
        print(f"[WARN] scANVI GPU training failed with CUDA error; retrying on CPU. Error: {exc}")
        del model_scanvi
        cleanup_cuda_cache()
        model_scvi_cpu = scvi.model.SCVI.load(str(scvi_model_dir), adata=adata)
        model_scanvi = build_scanvi_from_scvi(model_scvi_cpu)
        t0 = time.time()
        model_scanvi.train(accelerator="cpu", **train_kwargs)
        print(f"[OK] scANVI CPU fallback done in {(time.time() - t0) / 60:.1f} min")

    del model_scvi
    cleanup_cuda_cache()

    model_scanvi.save(str(model_dir), overwrite=True)
    pd.Series(adata.var_names.tolist()).to_csv(model_dir / "SCANVI_var_names.csv", index=False, header=False)

    adata.obsm["X_scanvi"] = model_scanvi.get_latent_representation()
    adata.obs["cell_type_scanvi_pred"] = model_scanvi.predict()
    soft_pred = model_scanvi.predict(soft=True)
    soft_pred_arr = soft_pred.to_numpy() if hasattr(soft_pred, "to_numpy") else np.asarray(soft_pred)
    adata.obs["scanvi_uncertainty"] = (1 - soft_pred_arr.max(axis=1)).astype(np.float32)
    del model_scanvi
    cleanup_cuda_cache()
    return True


def compute_umap(adata: AnnData, rep_key: str, random_seed: int = RANDOM_SEED) -> None:
    sc.pp.neighbors(adata, use_rep=rep_key, n_neighbors=30, random_state=random_seed, key_added=f"neighbors_{rep_key}")
    sc.tl.umap(adata, neighbors_key=f"neighbors_{rep_key}", random_state=random_seed)


def save_branch_figures(adata: AnnData, branch_name: str, branch_dir: Path, skip_figures: bool) -> None:
    if skip_figures:
        return
    fig_dir = branch_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.vector_friendly = True

    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    sc.pl.embedding(adata, basis="umap", color=OUT_L2, ax=axes[0], show=False, frameon=False, size=4,
                    legend_loc="right margin", legend_fontsize=8, title=f"{branch_name}: L2")
    sc.pl.embedding(adata, basis="umap", color="scanvi_label", ax=axes[1], show=False, frameon=False, size=4,
                    legend_loc="right margin", legend_fontsize=6, title=f"{branch_name}: scANVI training label")
    sc.pl.embedding(adata, basis="umap", color="cell_type_scanvi_pred", ax=axes[2], show=False, frameon=False, size=4,
                    legend_loc="right margin", legend_fontsize=6, title=f"{branch_name}: scANVI prediction")
    plt.tight_layout()
    fig.savefig(fig_dir / f"{branch_name}_scanvi_overview.pdf", dpi=300, bbox_inches="tight")
    plt.close("all")


def sanitize_obs_for_write(adata: AnnData) -> None:
    for df in [adata.obs, adata.var] + ([adata.raw.var] if adata.raw is not None else []):
        for col in df.columns:
            s = df[col]
            if s.dtype != object:
                continue
            non_na = s.dropna()
            if len(non_na) == 0:
                df[col] = s.fillna("").astype(str)
                continue
            if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
                df[col] = s.fillna(False).astype(np.int8)
                continue
            if non_na.map(lambda x: not isinstance(x, str)).any():
                df[col] = s.map(lambda x: "" if pd.isna(x) else str(x))


def branch_train_and_save(
    adata_clean: AnnData,
    branch_name: str,
    branch_cfg: Dict[str, object],
    cfg: PipelineConfig,
    tissue_key_available: bool,
) -> Dict[str, object]:
    log_header(f"BRANCH TRAINING: {branch_name}")
    branch_dir = cfg.output_dir / branch_name
    model_root = branch_dir / "models"
    model_root.mkdir(parents=True, exist_ok=True)

    ad = subset_branch(adata_clean, branch_name, branch_cfg)
    if ad.n_obs == 0:
        raise ValueError(f"[ERROR] Branch '{branch_name}' has 0 cells after filtering")

    print(f"Branch cells: {ad.n_obs:,} x {ad.n_vars:,}")
    print(f"L2 distribution:\n{ad.obs[OUT_L2].value_counts().to_string()}")
    print(f"L3 distribution (top 20):\n{ad.obs[OUT_L3].value_counts().head(20).to_string()}")

    compute_branch_covariates(ad, branch_name)

    hvg_method = select_hvgs_and_preserve_raw(ad, cfg, branch_dir)
    scvi_model_dir = model_root / f"scvi_{branch_name}_{PIPELINE_TAG}"
    scanvi_model_dir = model_root / f"scanvi_{branch_name}_{PIPELINE_TAG}"

    model_scvi = train_scvi(ad, cfg, scvi_model_dir)
    compute_umap(ad, "X_scvi")
    scanvi_labels, n_labeled, n_unknown = build_scanvi_labels(ad, cfg, tissue_key_available)
    scanvi_trained = train_scanvi(ad, model_scvi, cfg, scvi_model_dir, scanvi_model_dir)
    if "X_scanvi" not in ad.obsm:
        ad.obsm["X_scanvi"] = ad.obsm["X_scvi"].copy()
    compute_umap(ad, "X_scanvi")

    ad.uns["branchwise_reference"] = {
        "pipeline_version": PIPELINE_VERSION,
        "branch": branch_name,
        "branch_display": branch_cfg["display"],
        "allowed_l2": list(branch_cfg["allowed_l2"]),
        "query_subset": branch_cfg["query_subset"],
        "label_key": "scanvi_label",
        "latent_key": "X_scanvi",
        "hvg_method": hvg_method,
        "n_cells": int(ad.n_obs),
        "n_hvg": int(ad.n_vars),
        "n_raw_genes": int(ad.raw.n_vars if ad.raw is not None else ad.n_vars),
        "n_labeled": int(n_labeled),
        "n_unknown": int(n_unknown),
        "scanvi_trained": bool(scanvi_trained),
        "continuous_covariates": list(cfg.continuous_covariates),
        "encode_covariates": bool(cfg.encode_covariates),
        "use_layer_norm": cfg.use_layer_norm,
        "use_batch_norm": cfg.use_batch_norm,
    }

    label_summary_cols = [
        BATCH_KEY,
        OUT_L1,
        OUT_L2,
        OUT_L3,
        "annot_action",
        "reference_branch",
        "scanvi_label",
        "label_purity_scanvi",
        "cell_type_scanvi_pred",
        "scanvi_uncertainty",
    ]
    if cfg.tissue_key in ad.obs.columns:
        label_summary_cols.insert(1, cfg.tissue_key)
    ad.obs[label_summary_cols].to_csv(branch_dir / f"{branch_name}_label_summary.csv")

    sanitize_obs_for_write(ad)
    out_h5ad = branch_dir / f"adata_{branch_name}_reference_{PIPELINE_TAG}.h5ad"
    ad.write_h5ad(out_h5ad, compression="gzip", compression_opts=9)
    save_branch_figures(ad, branch_name, branch_dir, cfg.skip_figures)

    metrics = {
        "branch": branch_name,
        "display": branch_cfg["display"],
        "n_cells": int(ad.n_obs),
        "n_hvg": int(ad.n_vars),
        "n_raw_genes": int(ad.raw.n_vars if ad.raw is not None else ad.n_vars),
        "n_l3_labels": int(ad.obs[OUT_L3].nunique()),
        "n_scanvi_labeled": int(n_labeled),
        "n_scanvi_unknown": int(n_unknown),
        "scanvi_trained": bool(scanvi_trained),
        "allowed_l2": list(branch_cfg["allowed_l2"]),
        "query_subset": branch_cfg["query_subset"],
        "reference_h5ad": str(out_h5ad),
        "scvi_model_dir": str(scvi_model_dir),
        "scanvi_model_dir": str(scanvi_model_dir) if scanvi_trained else None,
        "hvg_file": str(branch_dir / "hvg_genes.csv"),
        "label_key": "scanvi_label",
        "latent_key": "X_scanvi",
        "continuous_covariate_keys": list(cfg.continuous_covariates),
        "encode_covariates": bool(cfg.encode_covariates),
        "use_layer_norm": cfg.use_layer_norm,
        "use_batch_norm": cfg.use_batch_norm,
    }
    print(json.dumps(metrics, indent=2, ensure_ascii=False))
    return metrics


def save_manifest(cfg: PipelineConfig, manifest: Dict[str, object]) -> Path:
    manifest_path = cfg.output_dir / "branch_reference_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path


def generate_synthetic_input(smoke_dir: Path) -> Path:
    log_header("SMOKE TEST: GENERATE SYNTHETIC INPUT")
    rng = np.random.default_rng(RANDOM_SEED)
    genes = [f"Gene{i:03d}" for i in range(120)]
    marker_map = {
        "Endothelia_Lymphatic_c0": [0, 1, 2, 3],
        "Endothelia_vascular_Cap_g_c1": [4, 5, 6, 7],
        "Fibro_adventitial_c0": [20, 21, 22, 23],
        "Fibro_alveolar_c1": [24, 25, 26, 27],
        "Muscle_pericyte_pulmonary_c0": [40, 41, 42, 43],
        "Muscle_smooth_pulmonary_c1": [44, 45, 46, 47],
        "Muscle_perivascular_immune_recruiting_c0": [48, 49, 50, 51],
        "Schwann_nonmyelinating_c0": [70, 71, 72, 73],
        "Endothelia_vascular_Cap_a_c0": [90, 91, 92, 93],
    }
    cluster_counts = {
        "Endothelia_Lymphatic_c0": 18,
        "Endothelia_vascular_Cap_g_c1": 18,
        "Fibro_adventitial_c0": 18,
        "Fibro_alveolar_c1": 18,
        "Muscle_pericyte_pulmonary_c0": 18,
        "Muscle_smooth_pulmonary_c1": 18,
        "Muscle_perivascular_immune_recruiting_c0": 12,
        "Schwann_nonmyelinating_c0": 10,
        "Endothelia_vascular_Cap_a_c0": 8,
    }

    rows = []
    X_blocks = []
    sample_levels = ["S1", "S2", "S3"]
    tissue_levels = ["lung", "trachea", "nose"]
    for cluster, n_cells in cluster_counts.items():
        for i in range(n_cells):
            lam = np.full(len(genes), 1.0)
            lam[marker_map[cluster]] = 8.0
            counts = rng.poisson(lam).astype(np.float32)
            X_blocks.append(counts)
            rows.append({
                SRC_L2_COL: cluster.rsplit("_c", 1)[0],
                SRC_L3_COL_DEFAULT: cluster,
                BATCH_KEY: sample_levels[i % len(sample_levels)],
                TISSUE_KEY: tissue_levels[(i + len(cluster)) % len(tissue_levels)],
                "subcluster_id": cluster.split("_c")[-1],
            })

    X = np.vstack(X_blocks)
    obs = pd.DataFrame(rows)
    obs.index = [f"synthetic_cell_{i:04d}" for i in range(len(obs))]
    var = pd.DataFrame(index=genes)
    adata = AnnData(X=sparse.csr_matrix(np.log1p(X)), obs=obs, var=var)
    adata.layers["counts"] = sparse.csr_matrix(X)
    adata.layers["log1p"] = sparse.csr_matrix(np.log1p(X))
    adata.X = adata.layers["log1p"]
    adata.raw = AnnData(X=adata.layers["counts"], obs=adata.obs.copy(), var=adata.var.copy())

    smoke_dir.mkdir(parents=True, exist_ok=True)
    smoke_input = smoke_dir / "synthetic_input_smoke.h5ad"
    adata.write_h5ad(smoke_input)
    print(f"Synthetic smoke input saved: {smoke_input}")
    return smoke_input


def run_pipeline(cfg: PipelineConfig) -> Dict[str, object]:
    set_seed(RANDOM_SEED)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    adata_clean, src_l3_col, tissue_key_available, n_input = load_and_clean_input(cfg)
    ignored_counts = adata_clean.obs[OUT_L2].astype(str).value_counts().to_dict()

    branch_results = {}
    for branch_name, branch_cfg in BRANCH_CONFIGS.items():
        branch_results[branch_name] = branch_train_and_save(
            adata_clean=adata_clean,
            branch_name=branch_name,
            branch_cfg=branch_cfg,
            cfg=cfg,
            tissue_key_available=tissue_key_available,
        )

    summary_df = pd.DataFrame(branch_results.values()).sort_values("branch")
    summary_df.to_csv(cfg.output_dir / "branch_training_summary.csv", index=False)

    manifest = {
        "pipeline_version": PIPELINE_VERSION,
        "pipeline_date": PIPELINE_DATE,
        "pipeline_tag": PIPELINE_TAG,
        "input_h5ad": str(cfg.input_h5ad),
        "output_dir": str(cfg.output_dir),
        "source_l3_col_used": src_l3_col,
        "tissue_key": cfg.tissue_key,
        "tissue_key_available": tissue_key_available,
        "ignored_l2": sorted(IGNORED_L2),
        "n_input_cells": int(n_input),
        "n_clean_cells": int(adata_clean.n_obs),
        "continuous_covariate_keys": list(cfg.continuous_covariates),
        "encode_covariates": bool(cfg.encode_covariates),
        "use_layer_norm": cfg.use_layer_norm,
        "use_batch_norm": cfg.use_batch_norm,
        "clean_l2_distribution": {k: int(v) for k, v in adata_clean.obs[OUT_L2].value_counts().items()},
        "branches": branch_results,
    }
    manifest_path = save_manifest(cfg, manifest)
    print(f"[OK] Manifest saved: {manifest_path}")
    return manifest


def parse_args() -> PipelineConfig:
    parser = argparse.ArgumentParser(description="Train branch-wise stromal scVI/scANVI references")
    parser.add_argument("--input-h5ad", type=Path,
                        default=Path("/home/h2048/data/py/0120/stromal_analysis_unified/results/"
                                     "subcluster_unified_v2_20260128/"
                                     "adata_stromal_subclustered_FINAL_v2_20260128.h5ad"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path(f"/home/h2048/data/py/0407/stromal_reintegration_{PIPELINE_TAG}"))
    parser.add_argument("--n-hvg", type=int, default=4000)
    parser.add_argument("--n-latent", type=int, default=75)
    parser.add_argument("--n-hidden", type=int, default=128)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--max-epochs-scvi", type=int, default=400)
    parser.add_argument("--max-epochs-scanvi", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--min-cells-batch", type=int, default=3)
    parser.add_argument("--purity-k", type=int, default=30)
    parser.add_argument("--purity-threshold", type=float, default=0.5)
    parser.add_argument("--source-l3-col", type=str, default=SRC_L3_COL_DEFAULT)
    parser.add_argument("--cpu-only", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()

    return PipelineConfig(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        n_hvg=args.n_hvg,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        n_layers=args.n_layers,
        dropout=args.dropout,
        max_epochs_scvi=args.max_epochs_scvi,
        max_epochs_scanvi=args.max_epochs_scanvi,
        batch_size=args.batch_size,
        min_cells_batch=args.min_cells_batch,
        purity_k=args.purity_k,
        purity_threshold=args.purity_threshold,
        source_l3_col=args.source_l3_col,
        use_gpu=(not args.cpu_only) and torch.cuda.is_available(),
        skip_figures=args.skip_figures,
        smoke_test=args.smoke_test,
    )


def main() -> None:
    cfg = parse_args()
    if cfg.smoke_test:
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        smoke_input = generate_synthetic_input(cfg.output_dir / "smoke_test")
        cfg.input_h5ad = smoke_input
        if cfg.n_hvg > 64:
            cfg.n_hvg = 64
        if cfg.n_latent > 8:
            cfg.n_latent = 8
        if cfg.n_hidden > 32:
            cfg.n_hidden = 32
        cfg.max_epochs_scvi = min(cfg.max_epochs_scvi, 1)
        cfg.max_epochs_scanvi = min(cfg.max_epochs_scanvi, 1)
        cfg.batch_size = min(cfg.batch_size, 64)
        cfg.skip_figures = True
        cfg.use_gpu = False

    start = time.time()
    manifest = run_pipeline(cfg)
    elapsed = (time.time() - start) / 60
    print("\n" + "=" * 80)
    print(f"PIPELINE COMPLETE ({elapsed:.1f} min)")
    print("=" * 80)
    for branch_name, info in manifest["branches"].items():
        print(
            f"  {branch_name:<12} cells={info['n_cells']:,} "
            f"scanvi_trained={info['scanvi_trained']} reference={info['reference_h5ad']}"
        )


if __name__ == "__main__":
    main()
