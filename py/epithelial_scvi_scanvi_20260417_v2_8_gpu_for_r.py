#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Epithelial scVI/scANVI v2.8-GPU (SELF-for-R)

Purpose
-------
- Re-run epithelial scVI/scANVI with GPU in `scvi_env`
- Export a lean SELF branch h5ad for downstream R tissue comparison
- Write final fine predictions to `cell_type_scanvi_pred`
- Keep curated labels as backup in `cell_type_L3_curated`
- Avoid the extra REF branch and heavy Python-side visualization to shorten restart latency

Date: 2026-04-17
"""

from __future__ import annotations

import gc
import json
import os
import time
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy import sparse
from sklearn.neighbors import NearestNeighbors

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

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1

# -----------------------------------------------------------------------------
# Threading / reproducibility
# -----------------------------------------------------------------------------
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
scvi.settings.seed = RANDOM_SEED

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
INPUT_H5AD = os.environ.get(
    "EPI_SCANVI_INPUT_H5AD",
    "/home/h2048/data/py/0122/epithelial_subcluster_v4_5_2_production/epithelial_with_subclusters_v4_5_2.h5ad",
)
OUTPUT_DIR = Path(
    os.environ.get(
        "EPI_SCANVI_OUTPUT_DIR",
        "/home/h2048/data/py/0417/epithelial_scanvi_v2_8_gpu_rerun",
    )
)
FINAL_H5AD = Path(
    os.environ.get(
        "EPI_SCANVI_FINAL_H5AD",
        str(OUTPUT_DIR / "epithelial_scanvi_v2_8_gpu_SELF_for_R.h5ad"),
    )
)
SCVI_MODEL_DIR = OUTPUT_DIR / "scvi_model"
SCANVI_MODEL_DIR = OUTPUT_DIR / "scanvi_model"
CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints"
SUMMARY_JSON = OUTPUT_DIR / "run_summary.json"
SUMMARY_CSV = OUTPUT_DIR / "scanvi_prediction_summary.csv"

for d in [OUTPUT_DIR, SCVI_MODEL_DIR, SCANVI_MODEL_DIR, CHECKPOINT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Basic config
# -----------------------------------------------------------------------------
BATCH_KEY = "sample"
UNLABELED_CATEGORY = "Unknown"
MIN_CELLS_PER_BATCH = 30
N_HVG = 4000
SCVI_MAX_EPOCHS = 300
SCANVI_MAX_EPOCHS = 120
BATCH_SIZE = 1024
EARLY_STOPPING = True
PURITY_THRESHOLD = 0.50
MT_THRESHOLD = 20.0
STRESS_PERCENTILE = 95
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3

TRAIN_ACCELERATOR = "gpu" if torch.cuda.is_available() else "cpu"
TRAIN_DEVICES = 1

# -----------------------------------------------------------------------------
# Label remap / lineage map
# -----------------------------------------------------------------------------
LABEL_REMAP = {
    "AT2_Canonical": "AT2",
    "AT2_Inflammatory_Repair": "AT2",
    "Epithelial_Cycling": "AT2_Cycling",
    "Goblet_Mucin": "SMG_Mucous",
    "Secretory_Club": "Goblet",
    "Secretory_Club_AT2_Transitional": "Club",
}
EXCLUDE_LABELS = {"Mesenchymal_Contaminant"}

L3_TO_L2_REMAP = {
    "AT1_Canonical": "Alveolar",
    "AT1_MatrixRemodeling": "Alveolar",
    "AT2": "Alveolar",
    "AT2_Cycling": "Alveolar",
    "Basal_Progenitor": "Basal_Lineage",
    "Basal_Cycling": "Basal_Lineage",
    "Basal_Inflammatory": "Basal_Lineage",
    "Basal_EMT_ECM": "Basal_Lineage",
    "Suprabasal_Progenitor": "Basal_Lineage",
    "Suprabasal_Cycling": "Basal_Lineage",
    "Ciliated_Mature": "Ciliated_Lineage",
    "Ciliogenesis_Deuterosomal": "Ciliated_Lineage",
    "Ciliated_Cycling_Immature": "Ciliated_Lineage",
    "Goblet": "Secretory_Lineage",
    "Club": "Secretory_Lineage",
    "SMG_Mucous": "Secretory_Lineage",
    "Goblet_Defense_DUOX2": "Secretory_Lineage",
    "SMG_Serous": "SMG",
    "SMG_Duct_Secretory_Defense": "SMG",
    "Squamous_Metaplasia": "Rare_Specialized",
    "Ionocyte_Brush": "Rare_Specialized",
}

CLUSTER_TO_L3 = {
    "AT1 _0": "AT1_Canonical",
    "AT1 _1": "AT1_MatrixRemodeling",
    "AT1 _2": "AT2_Canonical",
    "AT2_0": "AT2_Canonical",
    "AT2_1": "AT2_Canonical",
    "AT2_2": "AT2_Inflammatory_Repair",
    "AT2_3": "Epithelial_Cycling",
    "Basal_0": "Basal_EMT_ECM",
    "Basal_1": "Basal_Cycling",
    "Basal_2": "Basal_Inflammatory",
    "Basal_3": "Goblet_Mucin",
    "Basal_4": "Goblet_Mucin",
    "Basal_5": "Basal_Progenitor",
    "Ciliated_0": "Ciliated_Mature",
    "Ciliated_1": "Ciliated_Mature",
    "Ciliated_2": "Ciliated_Cycling_Immature",
    "Ciliated_3": "Ciliated_Mature",
    "Ciliated_4": "Ciliated_Mature",
    "Ciliated_5": "Goblet_Mucin",
    "Deuterosomal_0": "Ciliogenesis_Deuterosomal",
    "Deuterosomal_1": "Ciliogenesis_Deuterosomal",
    "Dividing_Basal_0": "Basal_Cycling",
    "Dividing_Basal_1": "Basal_Cycling",
    "Ionocyte_n_Brush_0": "Ionocyte_Brush",
    "Ionocyte_n_Brush_1": "Ionocyte_Brush",
    "SMG_Basal_0": "Basal_EMT_ECM",
    "SMG_Basal_1": "Mesenchymal_Contaminant",
    "SMG_Basal_2": "Basal_EMT_ECM",
    "SMG_Duct_0": "SMG_Duct_Secretory_Defense",
    "SMG_Duct_1": "Squamous_Metaplasia",
    "SMG_Duct_2": "Squamous_Metaplasia",
    "SMG_Duct_3": "Squamous_Metaplasia",
    "SMG_Duct_4": "Squamous_Metaplasia",
    "SMG_Mucous_0": "Goblet_Mucin",
    "SMG_Mucous_1": "Ionocyte_Brush",
    "SMG_Serous_0": "SMG_Serous",
    "SMG_Serous_1": "SMG_Serous",
    "SMG_Serous_2": "SMG_Serous",
    "Secretory_Club_0": "Secretory_Club_AT2_Transitional",
    "Secretory_Goblet_0": "Goblet_Defense_DUOX2",
    "Secretory_Goblet_1": "Secretory_Club",
    "Secretory_Goblet_2": "Basal_Progenitor",
    "Secretory_Goblet_3": "Squamous_Metaplasia",
    "Secretory_Goblet_4": "Ciliated_Cycling_Immature",
    "Suprabasal_0": "Suprabasal_Progenitor",
    "Suprabasal_1": "Suprabasal_Cycling",
    "Suprabasal_2": "Squamous_Metaplasia",
    "Suprabasal_3": "Squamous_Metaplasia",
}

STRESS_SIGNATURE_GENES = [
    "ALDH18A1", "ARFGAP1", "ASNS", "ATF3", "ATF4", "ATF6", "ATP6V0D1", "BAG3", "BANF1",
    "CALR", "CCL2", "CEBPB", "CEBPG", "CHAC1", "CKS1B", "CNOT2", "CNOT4", "CNOT6",
    "CXXC1", "DCP1A", "DCP2", "DCTN1", "DDIT4", "DDX10", "DKC1", "DNAJA4", "DNAJB9",
    "DNAJC3", "EDC4", "EDEM1", "EEF2", "EIF2AK3", "EIF2S1", "EIF4A1", "EIF4A2", "EIF4A3",
    "EIF4E", "EIF4EBP1", "EIF4G1", "ERN1", "ERO1A", "EXOC2", "EXOSC1", "EXOSC10",
    "EXOSC2", "EXOSC4", "EXOSC5", "EXOSC9", "FKBP14", "FUS", "GEMIN4", "GOSR2", "H2AX",
    "HERPUD1", "HSP90B1", "HSPA5", "HSPA9", "HYOU1", "IARS1", "IFIT1", "IGFBP1", "IMP3",
    "KDELR3", "KHSRP", "KIF5B", "LSM1", "LSM4", "MTHFD2", "NFYA", "NFYB", "NHP2", "NOLC1",
    "NOP14", "NOP56", "NPM1", "NABP1", "PAIP1", "PARN", "PDIA5", "PDIA6", "POP4", "PREB",
    "PSAT1", "RPS14", "RRP9", "SDAD1", "SEC11A", "SEC31A", "SERP1", "SHC1", "MTREX",
    "SLC1A4", "SLC30A5", "SLC7A5", "SPCS1", "SPCS3", "SRPRA", "SRPRB", "SSR1", "STC2",
    "TARS1", "TATDN2", "TSPYL2", "SKIC3", "TUBB2A", "VEGFA", "WFS1", "WIPI1", "XBP1",
    "XPOT", "YIF1A", "YWHAZ", "ZBTB17",
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

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%F %T')}] {msg}", flush=True)


def layer_flat_data(x):
    return x.data if sparse.issparse(x) else np.asarray(x).ravel()


def sanitize_sparse_counts(x):
    x = sparse.csr_matrix(x, dtype=np.float32, copy=True)
    bad = ~np.isfinite(x.data)
    n_bad = int(bad.sum())
    if n_bad > 0:
        x.data[bad] = 0.0
        x.eliminate_zeros()
    if np.any(x.data < 0):
        raise ValueError("counts contains negative values")
    return x, n_bad


def map_subclusters_to_l3(adata: ad.AnnData) -> ad.AnnData:
    if "subcluster" not in adata.obs.columns:
        raise ValueError("'subcluster' column not found in adata.obs")
    subcluster_clean = adata.obs["subcluster"].astype(str).str.strip()
    mapping_clean = {str(k).strip(): v for k, v in CLUSTER_TO_L3.items()}
    adata.obs["cell_type_L3_curated"] = subcluster_clean.map(mapping_clean)
    adata.obs.loc[adata.obs["cell_type_L3_curated"].isna(), "cell_type_L3_curated"] = subcluster_clean[
        adata.obs["cell_type_L3_curated"].isna()
    ]
    adata.obs["cell_type_L3_curated"] = adata.obs["cell_type_L3_curated"].replace(LABEL_REMAP)
    exclude_mask = adata.obs["cell_type_L3_curated"].isin(EXCLUDE_LABELS)
    if exclude_mask.any():
        log(f"Removing {int(exclude_mask.sum()):,} excluded cells")
        adata = adata[~exclude_mask].copy()
    adata.obs["cell_type_L3_curated"] = adata.obs["cell_type_L3_curated"].astype("category")
    return adata


def remove_small_batches(adata: ad.AnnData) -> ad.AnnData:
    batch_counts = adata.obs[BATCH_KEY].value_counts()
    small = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index.tolist()
    if small:
        log(f"Removing {len(small)} small batches (<{MIN_CELLS_PER_BATCH} cells)")
        adata = adata[~adata.obs[BATCH_KEY].isin(small)].copy()
    return adata


def ensure_log1p_layer(adata: ad.AnnData) -> None:
    if "log1p" in adata.layers:
        adata.X = adata.layers["log1p"].copy()
        return
    tmp = ad.AnnData(X=adata.layers["counts"].copy(), var=adata.var.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    adata.layers["log1p"] = tmp.X.copy()
    adata.X = adata.layers["log1p"].copy()
    del tmp
    gc.collect()


def compute_covariates(adata: ad.AnnData) -> float:
    if "pct_counts_mt" not in adata.obs.columns or (adata.obs["pct_counts_mt"] == 0).mean() > 0.1:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], layer="counts", percent_top=None, inplace=True)

    tmp = ad.AnnData(X=adata.layers["counts"].copy(), var=adata.var.copy())
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    sg = [g for g in STRESS_SIGNATURE_GENES if g in tmp.var_names]
    if len(sg) > 10:
        idx = [tmp.var_names.get_loc(g) for g in sg]
        x = tmp.X[:, idx]
        adata.obs["stress_score"] = np.asarray(x.mean(axis=1)).ravel() if sparse.issparse(x) else x.mean(axis=1)
    else:
        adata.obs["stress_score"] = 0.0

    s_genes = [g for g in S_GENES if g in tmp.var_names]
    g2m_genes = [g for g in G2M_GENES if g in tmp.var_names]
    if len(s_genes) > 10 and len(g2m_genes) > 10:
        sc.tl.score_genes_cell_cycle(tmp, s_genes=s_genes, g2m_genes=g2m_genes)
        adata.obs["S_score"] = tmp.obs["S_score"].values
        adata.obs["G2M_score"] = tmp.obs["G2M_score"].values
        adata.obs["phase"] = tmp.obs["phase"].values
    else:
        adata.obs["S_score"] = 0.0
        adata.obs["G2M_score"] = 0.0
        adata.obs["phase"] = "G1"

    del tmp
    gc.collect()
    return float(adata.obs["stress_score"].quantile(STRESS_PERCENTILE / 100.0))


def select_hvgs(adata: ad.AnnData) -> str:
    hvg_method = None
    for tier, layer_name, batch_key, flavor in [
        (1, "counts", BATCH_KEY, "seurat_v3"),
        (2, "counts", None, "seurat_v3"),
        (3, "log1p", BATCH_KEY, "seurat"),
        (4, "log1p", None, "seurat"),
    ]:
        try:
            kwargs = dict(layer=layer_name, n_top_genes=N_HVG, flavor=flavor, subset=False)
            if batch_key is not None:
                kwargs["batch_key"] = batch_key
            sc.pp.highly_variable_genes(adata, **kwargs)
            hvg_method = f"tier{tier}:{flavor}({layer_name},batch={'yes' if batch_key else 'no'})"
            break
        except Exception as exc:
            log(f"HVG tier {tier} failed: {exc}")
    if hvg_method is None or int(adata.var["highly_variable"].sum()) == 0:
        raise RuntimeError("No HVGs selected")
    return hvg_method


def build_unknown_mask(adata: ad.AnnData, stress_threshold: float) -> np.ndarray:
    labels = adata.obs["scanvi_labels_fine"].astype(str).values
    knn = NearestNeighbors(n_neighbors=30)
    knn.fit(adata.obsm["X_scvi"])
    indices = knn.kneighbors(return_distance=False)
    purity = []
    for i, nb in enumerate(indices):
        own = labels[i]
        if own == UNLABELED_CATEGORY:
            purity.append(1.0)
        else:
            purity.append(float((labels[nb[1:]] == own).mean()))
    adata.obs["label_purity_fine"] = purity
    low_purity = adata.obs["label_purity_fine"] < PURITY_THRESHOLD
    high_mt = adata.obs["pct_counts_mt"] > MT_THRESHOLD
    high_stress = adata.obs["stress_score"] > stress_threshold
    return (low_purity | high_mt | high_stress).to_numpy()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
pipeline_start = time.time()
log("=" * 80)
log("Epithelial scVI/scANVI v2.8-GPU (SELF-for-R)")
log(f"Input: {INPUT_H5AD}")
log(f"Output dir: {OUTPUT_DIR}")
log(f"Final h5ad: {FINAL_H5AD}")
log(f"torch.cuda.is_available() = {torch.cuda.is_available()}")
if torch.cuda.is_available():
    log(f"CUDA device count = {torch.cuda.device_count()}")
    log(f"CUDA device name = {torch.cuda.get_device_name(0)}")
log(f"Training accelerator = {TRAIN_ACCELERATOR}")
log("=" * 80)

adata = sc.read_h5ad(INPUT_H5AD)
log(f"Loaded input: {adata.shape}")
if "counts" not in adata.layers:
    raise ValueError("Input AnnData missing counts layer")

adata = map_subclusters_to_l3(adata)
adata = remove_small_batches(adata)

adata.layers["counts"], n_bad = sanitize_sparse_counts(adata.layers["counts"])
if n_bad > 0:
    log(f"Sanitized {n_bad} non-finite count entries")

ensure_log1p_layer(adata)
stress_threshold = compute_covariates(adata)
log(f"Global stress threshold (p{STRESS_PERCENTILE}) = {stress_threshold:.4f}")

n_cells_pre, n_genes_pre = adata.shape
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
log(f"filter_genes/filter_cells: {n_cells_pre}x{n_genes_pre} -> {adata.n_obs}x{adata.n_vars}")

hvg_method = select_hvgs(adata)
log(f"HVG method = {hvg_method}; selected = {int(adata.var['highly_variable'].sum())}")

adata.raw = ad.AnnData(X=adata.layers["log1p"], obs=adata.obs.copy(), var=adata.var.copy())
adata = adata[:, adata.var["highly_variable"]].copy()
adata.write_h5ad(CHECKPOINT_DIR / "adata_preprocessed_hvg.h5ad", compression="gzip")
log(f"Training subset saved: {adata.shape}")

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=BATCH_KEY,
    continuous_covariate_keys=["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
)

scvi_model = scvi.model.SCVI(
    adata,
    n_latent=64,
    n_layers=2,
    n_hidden=256,
    dropout_rate=0.2,
    gene_likelihood="nb",
    dispersion="gene-batch",
    encode_covariates=True,
    use_layer_norm="both",
    use_batch_norm="none",
)

train_kwargs = dict(
    max_epochs=SCVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=EARLY_STOPPING,
    train_size=0.9,
    accelerator=TRAIN_ACCELERATOR,
    devices=TRAIN_DEVICES,
    plan_kwargs={"lr": 1e-3},
)
log("Training scVI ...")
scvi_model.train(**train_kwargs)
scvi_model.save(str(SCVI_MODEL_DIR), overwrite=True)
adata.obsm["X_scvi"] = scvi_model.get_latent_representation()
log(f"scVI latent shape = {adata.obsm['X_scvi'].shape}")

adata.obs["scanvi_labels_fine"] = adata.obs["cell_type_L3_curated"].astype("category")
if UNLABELED_CATEGORY not in adata.obs["scanvi_labels_fine"].cat.categories:
    adata.obs["scanvi_labels_fine"] = adata.obs["scanvi_labels_fine"].cat.add_categories([UNLABELED_CATEGORY])

unknown_mask = build_unknown_mask(adata, stress_threshold)
adata.obs.loc[unknown_mask, "scanvi_labels_fine"] = UNLABELED_CATEGORY
log(
    "Unknown labels after purity/QC gate = "
    f"{int((adata.obs['scanvi_labels_fine'] == UNLABELED_CATEGORY).sum()):,} / {adata.n_obs:,}"
)

log("Initializing and training scANVI ...")
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    adata=adata,
    labels_key="scanvi_labels_fine",
    unlabeled_category=UNLABELED_CATEGORY,
)
scanvi_model.train(
    max_epochs=SCANVI_MAX_EPOCHS,
    batch_size=BATCH_SIZE,
    early_stopping=EARLY_STOPPING,
    train_size=0.9,
    n_samples_per_label=2000,
    accelerator=TRAIN_ACCELERATOR,
    devices=TRAIN_DEVICES,
    plan_kwargs={"lr": 1e-3},
)
scanvi_model.save(str(SCANVI_MODEL_DIR), overwrite=True)

adata.obsm["X_scanvi"] = scanvi_model.get_latent_representation()
adata.obsm["X_scanvi_fine"] = adata.obsm["X_scanvi"].copy()
adata.obsm["X_scanvi_major"] = adata.obsm["X_scanvi"].copy()

pred = scanvi_model.predict()
soft = scanvi_model.predict(soft=True)
conf = soft.values.max(axis=1) if isinstance(soft, pd.DataFrame) else np.asarray(soft).max(axis=1)

adata.obs["cell_type_scanvi_pred"] = pd.Categorical(pred)
adata.obs["scanvi_fine_pred"] = pd.Categorical(pred)
adata.obs["scanvi_fine_conf"] = conf
adata.obs["cell_type_scanvi_conf"] = conf
adata.obs["cell_type_L3"] = pd.Categorical(pred)
adata.obs["cell_type_L2"] = pd.Categorical(pd.Series(pred, index=adata.obs_names).map(L3_TO_L2_REMAP))
adata.obs["scanvi_major_pred"] = adata.obs["cell_type_L2"].astype(str)

missing_l2 = adata.obs["cell_type_L2"].isna()
if missing_l2.any():
    bad = sorted(pd.unique(adata.obs.loc[missing_l2, "cell_type_scanvi_pred"].astype(str)))
    raise ValueError(f"Unmapped predicted L3 labels for L2 remap: {bad}")

log("Running neighbors/UMAP on scANVI latent ...")
scanvi_umap_bundle = fit_bundle(
    adata,
    "X_scanvi",
    OUTPUT_DIR,
    primary_umap_key="X_umap_scanvi",
    alias_keys=("X_umap_fine", "X_umap_major"),
    operator_filename="epithelial_umap_scanvi_operator.joblib",
    manifest_filename="epithelial_umap_scanvi_bundle.json",
    umap_params={
        "n_neighbors": 15,
        "n_components": 2,
        "min_dist": 0.5,
        "spread": 1.0,
        "metric": "euclidean",
        "random_state": RANDOM_SEED,
    },
    set_default_x_umap=True,
    extra_manifest={"pipeline_version": "v2.8-GPU-SELF-for-R"},
)
log(f"Saved scanVI UMAP bundle: {scanvi_umap_bundle['operator_path']}")

adata.uns["pipeline_version"] = "v2.8-GPU-SELF-for-R"
adata.uns["hvg_method"] = hvg_method
adata.uns["stress_threshold"] = stress_threshold
adata.uns["train_accelerator"] = TRAIN_ACCELERATOR
adata.uns["train_devices"] = TRAIN_DEVICES

adata.write_h5ad(FINAL_H5AD, compression="gzip")
log(f"Saved final h5ad: {FINAL_H5AD}")

summary_df = (
    adata.obs["cell_type_scanvi_pred"]
    .value_counts()
    .rename_axis("cell_type_scanvi_pred")
    .reset_index(name="n_cells")
)
summary_df["pct"] = 100 * summary_df["n_cells"] / summary_df["n_cells"].sum()
summary_df["cell_type_L2"] = summary_df["cell_type_scanvi_pred"].map(L3_TO_L2_REMAP)
summary_df.to_csv(SUMMARY_CSV, index=False)

summary_payload = {
    "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    "input_h5ad": INPUT_H5AD,
    "final_h5ad": str(FINAL_H5AD),
    "umap_operator_path": str(scanvi_umap_bundle["operator_path"]),
    "umap_manifest_path": str(scanvi_umap_bundle["manifest_path"]),
    "shape": [int(adata.n_obs), int(adata.n_vars)],
    "train_accelerator": TRAIN_ACCELERATOR,
    "train_devices": TRAIN_DEVICES,
    "torch_cuda_available": bool(torch.cuda.is_available()),
    "hvg_method": hvg_method,
    "stress_threshold": stress_threshold,
    "unknown_n": int((adata.obs["scanvi_labels_fine"] == UNLABELED_CATEGORY).sum()),
    "predicted_l3_n": int(adata.obs["cell_type_scanvi_pred"].nunique()),
    "predicted_l2_n": int(adata.obs["cell_type_L2"].nunique()),
    "elapsed_minutes": round((time.time() - pipeline_start) / 60.0, 2),
}
with open(SUMMARY_JSON, "w", encoding="utf-8") as fh:
    json.dump(summary_payload, fh, indent=2)

log(f"Prediction summary saved: {SUMMARY_CSV}")
log(f"Run summary saved: {SUMMARY_JSON}")
log(f"Elapsed minutes = {(time.time() - pipeline_start) / 60.0:.2f}")
