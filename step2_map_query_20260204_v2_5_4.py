#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
==============================================================================
Step 2: Map Query to L2 Reference + Merge - v2.5.3 PRODUCTION (FULL FIXED)
==============================================================================

What this version fixes (the real NA root-cause):
1) Step 12 metadata restore now FORCE-overwrites result columns even if they
   already exist in the original query h5ad (previous runs left NA columns).
2) Merged UMAP panels for query use query-only masked temporary columns to avoid
   reference NaN/dtype issues.
3) Merge is done on HVG intersection in reference order (compact + stable).
4) Extra sanity checks + dtype normalization (no Categorical category errors).

Version: 2.5.3 PRODUCTION (FULL FIXED)
Date: 2026-02-04
Author: r2end
"""

import sys
import os
from pathlib import Path
import warnings
import json
import gc

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scanpy as sc
import scvi
import joblib

warnings.filterwarnings("ignore")

print("=" * 70)
print("B Cell Query Mapping + Merge - v2.5.3 PRODUCTION (FULL FIXED)")
print("=" * 70)
print(f"\nscanpy: {sc.__version__}")
print(f"scvi-tools: {scvi.__version__}")
print(f"Python: {sys.version}")

import torch
gpu_available = torch.cuda.is_available()
print(f"\nGPU available: {gpu_available}")
if gpu_available:
    print(f"GPU device: {torch.cuda.get_device_name(0)}")


def _env_path(name: str, default: str) -> str:
    value = os.getenv(name)
    if value is None or not value.strip():
        return default
    return value.strip()

# ==============================================================================
# CONFIGURATION
# ==============================================================================

print("\n" + "=" * 70)
print("CONFIGURATION")
print("=" * 70)

SCANVI_MODEL_DIR = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/scanvi_existing_model"
HVG_FILE = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/hvg_genes.txt"
REF_H5AD_FOR_MERGE = "/home/h2048/data/py/0203/bcell_scarches_v4_1/models/scanvi_bcell_L2_v2_5_3/reference_with_L2_umap.h5ad"

QUERY_H5AD = _env_path(
    "BCELL_QUERY_H5AD",
    "/home/h2048/data/py/0127/scarches_mapping_FIXED_v1_2/subsets/b_cells.h5ad",
)
OUTPUT_DIR = _env_path(
    "BCELL_OUTPUT_DIR",
    "/home/h2048/data/py/0204/scarches_mapping_L2_v2_5_3",
)

BATCH_KEY = "sample"
TISSUE_KEY = "tissue"
L2_KEY = "Cell_Type_L2"
UNLABELED_CATEGORY = "Unknown"

# Covariate gene lists (same as Step1)
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
    "XPOT","YIF1A","YWHAZ","ZBTB17","ARFGEF1","ATG9A","CAND1","CCAR1","CD58",
    "CUL4A","DBI","DDX17","DDX18","DDX21","DYNC1I2","EIF2B4","EIF3A","FAM120A",
    "FKBP1A","GANAB","GNB1","HNRNPA0","HNRNPD","HSPD1","ILF3","ISG15","ITGB1",
    "KRR1","MRPL11","MRPS18C","MVP","NCL","NMD3","NOP2","OGFOD1","PABPC1",
    "PHB","PHB2","POLR1B","PPP1CA","RANBP1","RBM3","RPL18A","RPL3","RPL7",
    "RPS3","RPS6","RTN4","SART1","SET","SSRP1","TCEB2","TMED10","TUBB","VDAC1",
    "VIM","VMP1"
]
S_GENES = [
    'MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1','UNG','GINS2','MCM6',
    'CDCA7','DTL','PRIM1','UHRF1','MLF1IP','HELLS','RFC2','RPA2','NASP',
    'RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7','POLD3','MSH2','ATAD2',
    'RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
    'USP1','CLSPN','POLA1','CHAF1B','BRIP1','E2F8'
]
G2M_GENES = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2',
    'NUF2','CKS1B','MKI67','TMPO','CENPF','TACC3','FAM64A','SMC4','CCNB2',
    'CKAP2L','CKAP2','AURKB','BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B',
    'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1','NCAPD2',
    'DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1','ANLN',
    'LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA'
]

MAX_EPOCHS = 200
WEIGHT_DECAY = 0.0
LEARNING_RATE = 5e-4
BATCH_SIZE = 256
CONFIDENCE_THRESHOLD = 0.5

DPI = 300
FIGURE_FORMAT = "pdf"
sc.settings.set_figure_params(dpi=DPI, facecolor="white", format=FIGURE_FORMAT)

print("\nConfiguration:")
print(f"  Model: {SCANVI_MODEL_DIR}")
print(f"  Query: {QUERY_H5AD}")
print(f"  Reference for merge: {REF_H5AD_FOR_MERGE}")
print(f"  Output: {OUTPUT_DIR}")

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
(output_dir / "figures").mkdir(exist_ok=True)
print(f"\nOK Output directory: {output_dir}")

# ==============================================================================
# Step 1: Load Query Data
# ==============================================================================

print("\n" + "=" * 70)
print("Step 1: Load Query Data")
print("=" * 70)

print(f"\nLoading: {QUERY_H5AD}")
adata_query = sc.read_h5ad(QUERY_H5AD)
adata_query.var_names_make_unique()
adata_query.obs_names_make_unique()
print("OK Gene/cell names made unique")
print(f"Shape: {adata_query.n_obs:,} cells × {adata_query.n_vars:,} genes")

# ==============================================================================
# Step 2: Ensure counts layer exists + use counts as X for covariates
# ==============================================================================

print("\n" + "=" * 70)
print("Step 2: Data Format Validation")
print("=" * 70)

if "counts" not in adata_query.layers:
    adata_query.layers["counts"] = adata_query.X
    print("Created layers['counts'] from .X")
else:
    print("OK layers['counts'] exists")

adata_query.X = adata_query.layers["counts"]
print("Set .X = layers['counts']")

# ==============================================================================
# Step 3: Calculate missing covariates BEFORE prepare_query_anndata
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3: Check and Calculate Covariates (IF MISSING)")
print("=" * 70)

required_covariates = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
missing_covariates = [c for c in required_covariates if c not in adata_query.obs.columns]

if missing_covariates:
    print(f"WARNING  Missing covariates: {missing_covariates}")
else:
    print("OK All required covariates exist (will not recompute)")

# 3.1 MT%
if "pct_counts_mt" not in adata_query.obs.columns:
    adata_query.var["mt"] = adata_query.var_names.str.startswith("MT-")
    if int(adata_query.var["mt"].sum()) > 0:
        sc.pp.calculate_qc_metrics(adata_query, qc_vars=["mt"], inplace=True, layer="counts")
    else:
        adata_query.obs["pct_counts_mt"] = 0.0
    print(f"OK pct_counts_mt mean: {adata_query.obs['pct_counts_mt'].mean():.2f}%")

# 3.2 stress
if "stress_score" not in adata_query.obs.columns:
    counts = adata_query.layers["counts"]
    stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in adata_query.var_names]
    if len(stress_genes) >= 10:
        idx = adata_query.var_names.get_indexer(stress_genes)
        idx = idx[idx >= 0]
        sub = counts[:, idx]
        if issparse(counts):
            tot = np.asarray(counts.sum(axis=1)).ravel()
        else:
            tot = counts.sum(axis=1)
        tot[tot == 0] = 1.0
        from scipy import sparse
        scale = sparse.diags(1e4 / tot)
        sub_norm = scale @ sub
        if issparse(sub_norm):
            sub_norm = sub_norm.tocsr(copy=True)
            sub_norm.data = np.log1p(sub_norm.data)
            expr = np.asarray(sub_norm.mean(axis=1)).ravel()
        else:
            expr = np.log1p(sub_norm).mean(axis=1)
        mn, mx = float(expr.min()), float(expr.max())
        adata_query.obs["stress_score"] = (expr - mn) / (mx - mn) if mx > mn else 0.0
        del sub, sub_norm, scale
        gc.collect()
    else:
        adata_query.obs["stress_score"] = 0.0
    print(f"OK stress_score mean: {adata_query.obs['stress_score'].mean():.3f}")

# 3.3 cell cycle
need_cc = ("S_score" not in adata_query.obs.columns) or ("G2M_score" not in adata_query.obs.columns)
if need_cc:
    s_genes = [g for g in S_GENES if g in adata_query.var_names]
    g2m_genes = [g for g in G2M_GENES if g in adata_query.var_names]
    if len(s_genes) >= 10 and len(g2m_genes) >= 10:
        adata_temp = sc.AnnData(X=adata_query.layers["counts"].copy(), var=adata_query.var.copy())
        sc.pp.normalize_total(adata_temp, target_sum=1e4)
        sc.pp.log1p(adata_temp)
        sc.tl.score_genes_cell_cycle(adata_temp, s_genes=s_genes, g2m_genes=g2m_genes)
        adata_query.obs["S_score"] = adata_temp.obs["S_score"].to_numpy()
        adata_query.obs["G2M_score"] = adata_temp.obs["G2M_score"].to_numpy()
        adata_query.obs["phase"] = adata_temp.obs["phase"].to_numpy()
        del adata_temp
        gc.collect()
    else:
        adata_query.obs["S_score"] = 0.0
        adata_query.obs["G2M_score"] = 0.0
        adata_query.obs["phase"] = "G1"
    print(f"OK S_score mean: {adata_query.obs['S_score'].mean():.3f}")
    print(f"OK G2M_score mean: {adata_query.obs['G2M_score'].mean():.3f}")

# ==============================================================================
# Step 3.5: Validate batch/tissue keys exist
# ==============================================================================

print("\n" + "=" * 70)
print("Step 3.5: Validate Batch and Tissue Keys")
print("=" * 70)

if BATCH_KEY not in adata_query.obs.columns:
    adata_query.obs[BATCH_KEY] = "query_batch"
    print(f"OK Created '{BATCH_KEY}'")
if TISSUE_KEY not in adata_query.obs.columns:
    adata_query.obs[TISSUE_KEY] = "query_tissue"
    print(f"OK Created '{TISSUE_KEY}'")

# ==============================================================================
# Step 4: Validate tissue categories vs reference config
# ==============================================================================

print("\n" + "=" * 70)
print("Step 4: Validate Tissue Categories (P1-1 FIX)")
print("=" * 70)

config_file = Path(SCANVI_MODEL_DIR).parent / "training_config.json"
if not config_file.exists():
    raise FileNotFoundError(f"Training config not found: {config_file}")

with open(config_file, "r") as f:
    ref_config = json.load(f)

ref_tissues = list(ref_config.get("tissue_categories", []))
ref_tissue_set = set(ref_tissues)
print(f"Reference tissue categories: {sorted(ref_tissue_set)}")

t = adata_query.obs[TISSUE_KEY].astype(str).fillna("unknown_tissue")
t = t.where(t.isin(ref_tissue_set), "unknown_tissue")
adata_query.obs[TISSUE_KEY] = pd.Categorical(t, categories=sorted(ref_tissue_set))

print("Query tissue distribution after alignment:")
print(adata_query.obs[TISSUE_KEY].value_counts())

# ==============================================================================
# Step 5: HVG overlap check
# ==============================================================================

print("\n" + "=" * 70)
print("Step 5: Gene Alignment with Reference (HVG overlap)")
print("=" * 70)

hvg_genes_ref = pd.read_csv(HVG_FILE, header=None)[0].astype(str).tolist()
query_genes = set(adata_query.var_names.astype(str))
hvg_in_query = [g for g in hvg_genes_ref if g in query_genes]
overlap_pct = len(hvg_in_query) / max(1, len(hvg_genes_ref)) * 100
print(f"Reference HVG: {len(hvg_genes_ref):,}")
print(f"Overlap: {len(hvg_in_query):,}/{len(hvg_genes_ref):,} ({overlap_pct:.1f}%)")

if overlap_pct < 70:
    raise ValueError(f"Insufficient gene overlap ({overlap_pct:.1f}%).")

# ==============================================================================
# Step 6: Prepare query for scArches (gene pad/reorder)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 6: Prepare Query for scArches")
print("=" * 70)

scvi.model.SCANVI.prepare_query_anndata(adata_query, SCANVI_MODEL_DIR)
print(f"OK Query prepared. Shape after prep: {adata_query.shape}")

# ==============================================================================
# Step 7: Clean obs for model loading (strict types)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 7: Metadata Cleaning (strict categorical typing)")
print("=" * 70)

# Keep a flag so Step12 always runs (even if scvi touches .uns)
adata_query.uns["_restore_original_metadata"] = True

obs_clean = pd.DataFrame(index=adata_query.obs_names)

obs_clean[BATCH_KEY] = pd.Categorical(adata_query.obs[BATCH_KEY].astype(str).fillna("unknown_batch"))
obs_clean[TISSUE_KEY] = adata_query.obs[TISSUE_KEY]

# Labels for scANVI (all Unknown)
obs_clean[L2_KEY] = pd.Categorical(
    [UNLABELED_CATEGORY] * adata_query.n_obs,
    categories=[UNLABELED_CATEGORY],
)

for cov in ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]:
    obs_clean[cov] = adata_query.obs[cov].to_numpy()

adata_query.obs = obs_clean
print("OK Metadata cleaned")

# ==============================================================================
# Step 8: Load query into reference model
# ==============================================================================

print("\n" + "=" * 70)
print("Step 8: Load Query into Reference Model (scArches surgery)")
print("=" * 70)

scanvi_query = scvi.model.SCANVI.load_query_data(adata_query, SCANVI_MODEL_DIR)
print("OK scArches surgery successful")

# ==============================================================================
# Step 9: Fine-tune
# ==============================================================================

print("\n" + "=" * 70)
print("Step 9: Fine-tune Model on Query Data")
print("=" * 70)

train_kwargs = dict(
    max_epochs=MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=True,
    early_stopping_patience=30,
    plan_kwargs={"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
)
if gpu_available:
    train_kwargs["accelerator"] = "gpu"
    train_kwargs["devices"] = "auto"

scanvi_query.train(**train_kwargs)
print("OK Fine-tuning complete")

# ==============================================================================
# Step 10: Extract predictions + latent + confidence
# ==============================================================================

print("\n" + "=" * 70)
print("Step 10: Extract Predictions and Latent")
print("=" * 70)

adata_query.obsm["X_scANVI_L2"] = scanvi_query.get_latent_representation()
print(f"OK X_scANVI_L2: {adata_query.obsm['X_scANVI_L2'].shape}")

pred = scanvi_query.predict()
adata_query.obs["Cell_Type_L2_pred"] = pd.Series(pred, index=adata_query.obs_names).astype(str)

probs = scanvi_query.predict(soft=True)
if hasattr(probs, "to_numpy"):
    probs = probs.to_numpy()
confidence = probs.max(axis=1)
adata_query.obs["mapping_confidence"] = pd.to_numeric(confidence, errors="coerce")

final = np.where(
    adata_query.obs["mapping_confidence"].to_numpy() >= CONFIDENCE_THRESHOLD,
    adata_query.obs["Cell_Type_L2_pred"].to_numpy(),
    UNLABELED_CATEGORY,
)
adata_query.obs["Cell_Type_L2_final"] = pd.Series(final, index=adata_query.obs_names).astype(str)

print(f"OK Predicted L2 types: {adata_query.obs['Cell_Type_L2_pred'].nunique()}")
print(f"OK Confidence median: {np.nanmedian(adata_query.obs['mapping_confidence'].to_numpy()):.3f}")
print("Top Cell_Type_L2_final:")
print(adata_query.obs["Cell_Type_L2_final"].value_counts(dropna=False).head(10))

# ==============================================================================
# Step 11: UMAP projection (same space)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 11: UMAP Projection (same space)")
print("=" * 70)

umap_operator_file = Path(SCANVI_MODEL_DIR).parent / "umap_operator.joblib"
if not umap_operator_file.exists():
    raise FileNotFoundError(f"UMAP operator not found: {umap_operator_file}")

umap_operator = joblib.load(umap_operator_file)
adata_query.obsm["X_umap"] = umap_operator.transform(adata_query.obsm["X_scANVI_L2"])
print(f"OK X_umap: {adata_query.obsm['X_umap'].shape}")

# ==============================================================================
# Step 12: Restore original query metadata (FULL FIX: FORCE OVERWRITE always)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 12: Restore Metadata and Save Query (FORCE OVERWRITE ALWAYS)")
print("=" * 70)

# Always reload original query and overwrite (no conditional)
adata_original = sc.read_h5ad(QUERY_H5AD)
adata_original.var_names_make_unique()
adata_original.obs_names_make_unique()

# align rows
adata_original = adata_original[adata_query.obs_names].copy()

RESULT_COLS_FORCE = [
    "Cell_Type_L2_pred", "mapping_confidence", "Cell_Type_L2_final",
    "pct_counts_mt", "stress_score", "S_score", "G2M_score", "phase",
    BATCH_KEY, TISSUE_KEY,
]

for col in RESULT_COLS_FORCE:
    if col in adata_query.obs.columns:
        # if original has categorical -> object first to avoid "new category"
        if col in adata_original.obs.columns:
            try:
                if pd.api.types.is_categorical_dtype(adata_original.obs[col]):
                    adata_original.obs[col] = adata_original.obs[col].astype("object")
            except Exception:
                pass
        adata_original.obs[col] = pd.Series(adata_query.obs[col].to_numpy(), index=adata_original.obs_names)

# enforce clean dtypes
if "mapping_confidence" in adata_original.obs.columns:
    adata_original.obs["mapping_confidence"] = pd.to_numeric(adata_original.obs["mapping_confidence"], errors="coerce")

if "Cell_Type_L2_final" in adata_original.obs.columns:
    adata_original.obs["Cell_Type_L2_final"] = adata_original.obs["Cell_Type_L2_final"].astype(str).astype("category")

# overwrite embeddings
for k in ["X_scANVI_L2", "X_umap"]:
    if k in adata_query.obsm:
        adata_original.obsm[k] = adata_query.obsm[k]

adata_query = adata_original
del adata_original
gc.collect()

print("\nSanity after restore (MUST NOT be NA):")
print(adata_query.obs["Cell_Type_L2_final"].value_counts(dropna=False).head(10))
print(f"mapping_confidence NA rate = {pd.isna(adata_query.obs['mapping_confidence']).mean():.3f}")

# Save query results
query_output = output_dir / "query_mapped_L2.h5ad"
print(f"\nSaving query results: {query_output}")
adata_query.write_h5ad(query_output, compression="gzip")
print(f"OK Saved ({query_output.stat().st_size / 1024**3:.2f} GB)")

# ==============================================================================
# Step 13: Merge with reference (HVG-only, reference order)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 13: Merge Reference + Query (HVG-only, stable)")
print("=" * 70)

adata_all = None
merged_output = None

if REF_H5AD_FOR_MERGE:
    adata_ref = sc.read_h5ad(REF_H5AD_FOR_MERGE)
    adata_ref.var_names_make_unique()
    adata_ref.obs_names_make_unique()

    if "X_umap" not in adata_ref.obsm:
        raise RuntimeError("Reference missing X_umap in REF_H5AD_FOR_MERGE.")

    # Use reference var_names as authoritative HVG list for merge
    ref_merge_genes = list(adata_ref.var_names)
    common = [g for g in ref_merge_genes if g in adata_query.var_names]
    frac = len(common) / max(1, len(ref_merge_genes))
    print(f"Merge gene overlap vs reference HVG: {len(common):,}/{len(ref_merge_genes):,} ({frac*100:.1f}%)")
    if frac < 0.95:
        print("WARNING: low overlap with reference HVGs; merged expression matrix may be incomplete.")

    adata_ref_m = adata_ref[:, common].copy()
    adata_qry_m = adata_query[:, common].copy()

    # keep minimal obsm
    for ad in (adata_ref_m, adata_qry_m):
        for k in list(ad.obsm.keys()):
            if k not in ["X_umap", "X_scANVI_L2"]:
                del ad.obsm[k]
        ad.uns = {}
        if issparse(ad.X):
            ad.X = csr_matrix(ad.X)

    adata_ref_m.obs_names = [f"ref::{x}" for x in adata_ref_m.obs_names]
    adata_qry_m.obs_names = [f"qry::{x}" for x in adata_qry_m.obs_names]

    adata_all = sc.concat(
        {"reference": adata_ref_m, "query": adata_qry_m},
        axis=0,
        join="inner",
        merge="unique",
        label="data_source",
    )

    if "X_umap" not in adata_all.obsm or adata_all.obsm["X_umap"].shape[0] != adata_all.n_obs:
        raise RuntimeError("UMAP lost/misaligned during merge.")

    merged_output = output_dir / "reference_plus_query_merged_L2.h5ad"
    print(f"Saving merged: {merged_output}")
    adata_all.write_h5ad(merged_output, compression="gzip")
    print(f"OK Saved ({merged_output.stat().st_size / 1024**3:.2f} GB)")

    del adata_ref, adata_ref_m, adata_qry_m
    gc.collect()

# ==============================================================================
# Step 14: Visualizations (query panels FIXED to never be NA)
# ==============================================================================

print("\n" + "=" * 70)
print("Step 14: Generate Visualizations")
print("=" * 70)

# Merged plots
if adata_all is not None:
    print("\nGenerating merged UMAP plots...")
    fig, axes = plt.subplots(2, 3, figsize=(20, 13))

    sc.pl.umap(
        adata_all, color="data_source", ax=axes[0, 0], show=False,
        title="Data Source", s=5,
        palette={"reference": "#1f77b4", "query": "#ff7f0e"},
    )

    if L2_KEY in adata_all.obs.columns:
        sc.pl.umap(adata_all, color=L2_KEY, ax=axes[0, 1], show=False,
                   title="Reference L2 Labels", legend_loc="right margin", s=5)
    else:
        axes[0, 1].text(0.5, 0.5, "N/A", ha="center", va="center", transform=axes[0, 1].transAxes)

    # ---- FIX: query-only masked L2 final
    qry_mask = (adata_all.obs["data_source"] == "query").to_numpy()
    adata_all.obs["_qry_L2_final_only"] = pd.Series(pd.NA, index=adata_all.obs_names, dtype="object")
    if "Cell_Type_L2_final" in adata_all.obs.columns:
        adata_all.obs.loc[qry_mask, "_qry_L2_final_only"] = adata_all.obs.loc[qry_mask, "Cell_Type_L2_final"].astype(str).values
    adata_all.obs["_qry_L2_final_only"] = adata_all.obs["_qry_L2_final_only"].astype("category")

    sc.pl.umap(
        adata_all, color="_qry_L2_final_only", ax=axes[0, 2], show=False,
        title=f"Query L2 Final (Query Only, thr>={CONFIDENCE_THRESHOLD})",
        legend_loc="right margin", s=5
    )
    adata_all.obs.drop(columns=["_qry_L2_final_only"], inplace=True)

    sc.pl.umap(adata_all, color=BATCH_KEY, ax=axes[1, 0], show=False, title="Batch", s=5)

    # ---- FIX: query-only masked confidence
    tmp = np.full(adata_all.n_obs, np.nan, dtype=float)
    if "mapping_confidence" in adata_all.obs.columns:
        tmp[qry_mask] = pd.to_numeric(adata_all.obs.loc[qry_mask, "mapping_confidence"], errors="coerce").values
    adata_all.obs["_qry_conf_only"] = tmp
    sc.pl.umap(
        adata_all, color="_qry_conf_only", ax=axes[1, 1], show=False,
        title="Mapping Confidence (Query Only)", cmap="viridis", vmin=0, vmax=1, s=5
    )
    adata_all.obs.drop(columns=["_qry_conf_only"], inplace=True)

    # marker
    gene = "CD19"
    if gene in adata_all.var_names:
        sc.pl.umap(adata_all, color=gene, ax=axes[1, 2], show=False, cmap="Reds", title="CD19", s=5, use_raw=False)
    else:
        axes[1, 2].text(0.5, 0.5, "CD19\nN/A", ha="center", va="center", transform=axes[1, 2].transAxes)

    plt.tight_layout()
    fig_path = output_dir / "figures" / f"merged_umap_reference_query.{FIGURE_FORMAT}"
    plt.savefig(fig_path, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"OK Saved: {fig_path.name}")

# Query-only plots
print("\nGenerating query-specific plots...")
fig, axes = plt.subplots(2, 2, figsize=(16, 16))

sc.pl.umap(
    adata_query, color="Cell_Type_L2_pred", ax=axes[0, 0], show=False,
    title="L2 Predictions (All)", legend_loc="right margin"
)
sc.pl.umap(
    adata_query, color="Cell_Type_L2_final", ax=axes[0, 1], show=False,
    title=f"L2 Final (thr>={CONFIDENCE_THRESHOLD})", legend_loc="right margin"
)
sc.pl.umap(
    adata_query, color="mapping_confidence", ax=axes[1, 0], show=False,
    title="Mapping Confidence", cmap="viridis", vmin=0, vmax=1
)

conf_vals = pd.to_numeric(adata_query.obs["mapping_confidence"], errors="coerce").to_numpy()
axes[1, 1].hist(conf_vals[np.isfinite(conf_vals)], bins=50, edgecolor="black", alpha=0.7)
axes[1, 1].axvline(CONFIDENCE_THRESHOLD, color="red", linestyle="--", linewidth=2, label=f"thr={CONFIDENCE_THRESHOLD}")
axes[1, 1].set_xlabel("Mapping Confidence")
axes[1, 1].set_ylabel("Number of Cells")
axes[1, 1].set_title("Confidence Distribution")
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
fig_path = output_dir / "figures" / f"query_mapping_summary.{FIGURE_FORMAT}"
plt.savefig(fig_path, dpi=DPI, bbox_inches="tight")
plt.close()
print(f"OK Saved: {fig_path.name}")

# ==============================================================================
# SUMMARY
# ==============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Query saved:  {query_output}")
if merged_output is not None:
    print(f"Merged saved: {merged_output}")
print("Top Cell_Type_L2_final:")
print(adata_query.obs["Cell_Type_L2_final"].value_counts(dropna=False).head(10))
print(f"mapping_confidence NA rate = {pd.isna(adata_query.obs['mapping_confidence']).mean():.3f}")

print("\n" + "=" * 70)
print("SUCCESS")
print("=" * 70)
