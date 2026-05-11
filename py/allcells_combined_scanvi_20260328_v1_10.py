# %% [markdown]
# # All-Lineage Combined scVI + scANVI Pipeline — v1.10 L1 PRODUCTION
# **Author:** r2end | **Date:** 2026-03-19 | **Version:** v1.10 L1 PRODUCTION
# 
# ## Changes over v1.9
# 
# | Tag | Fix |
# |-----|-----|
# | v1.10-1 | **[P0-2]** Cell 6: overlap audit now checks `cell_type_L2_original` / `cell_type_L3_original`; old L1 audit removed (L1 = lineage names, always unique, audit was vacuous) |
# | v1.10-2 | **[P0-3]** Cell 18: agreement labelled as train-set re-substitution; note added that it cannot serve as independent validation |
# | v1.10-3 | **[P1-2]** `cell_type_final` removed; `cell_type_final_l1` is the sole scANVI output column |
# | v1.10-4 | **[P1-3]** `n_samples_per_label` removed from scANVI training (5-class L1 task; resampling adds no benefit and distorts class priors) |
# 
# ## Changes over v1.8 (v1.9)
# 
# | Tag | Fix |
# |-----|-----|
# | v1.9-1 | **[CRITICAL]** Add stress_score, S_score, G2M_score covariates (Cell 9 expanded) — matching reference pipeline v2.5.2 |
# | v1.9-2 | **[CRITICAL]** scVI: dropout_rate=0.2, use_layer_norm='both', use_batch_norm='none', early_stopping_patience=45 |
# | v1.9-3 | **[CRITICAL]** scVI N_LATENT changed from 175 -> 150 to match reference |
# | v1.9-4 | scANVI: lr=5e-4 (was 1e-3) |
# | v1.9-5 | CONTINUOUS_COVARIATES expanded: pct_counts_mt + stress_score + S_score + G2M_score |
# | v1.9-6 | setup_anndata layer='counts' (consistent with .X=log1p after Cell 14) |

# %% [markdown]
# ## Cell 1 — Logging Setup + Imports

# %%
import os, sys
os.environ["OMP_NUM_THREADS"]      = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"]      = "8"

import matplotlib
matplotlib.use("Agg")

import gc, time, json, anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
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

class _Tee:
    def __init__(self, original, log_fh):
        self._orig = original; self._log = log_fh
    def write(self, data):
        self._orig.write(data); self._orig.flush()
        try: self._log.write(data); self._log.flush()
        except Exception: pass
    def flush(self):
        self._orig.flush()
        try: self._log.flush()
        except Exception: pass
    def fileno(self):  return self._orig.fileno()
    def isatty(self):  return False

_LOG_BASE = Path("/home/h2048/data/py/20260319/allcells_combined_scanvi/logs")
_LOG_BASE.mkdir(parents=True, exist_ok=True)
_TS      = datetime.now().strftime("%Y%m%d_%H%M%S")
_NB_STEM = "allcells_combined_scanvi_20260319_v1_10_PRODUCTION"
LOG_FILE = _LOG_BASE / f"{_NB_STEM}_{_TS}.log"
_log_fh  = open(LOG_FILE, "w", encoding="utf-8", buffering=1)
sys.stdout = _Tee(sys.__stdout__, _log_fh)
sys.stderr = _Tee(sys.__stderr__, _log_fh)

anndata.settings.allow_write_nullable_strings = True
GPU_AVAILABLE = torch.cuda.is_available()
_accelerator  = "gpu" if GPU_AVAILABLE else "cpu"
_devices      = 1     if GPU_AVAILABLE else "auto"
scvi.settings.seed = 42
if GPU_AVAILABLE: torch.cuda.manual_seed_all(42)
np.random.seed(42)
sc.settings.n_jobs    = 8
sc.settings.verbosity = 2
PIPELINE_START = time.time()

print("=" * 70)
print(f"Pipeline   : {_NB_STEM}")
print(f"Log file   : {LOG_FILE}")
print(f"Start      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"scvi-tools : {scvi.__version__}")
print(f"scanpy     : {sc.__version__}")
print(f"Accelerator: {_accelerator}  devices={_devices}")
print("=" * 70)
print(f"  tail -f {LOG_FILE}")

# %% [markdown]
# ## Cell 2 — Configuration

# %%
DATE_TAG   = "20260319"
VERSION    = "1_10_L1"
OUTPUT_DIR = Path(f"/home/h2048/data/py/{DATE_TAG}/allcells_combined_scanvi")
MODEL_DIR  = OUTPUT_DIR / "models"
FIG_DIR    = OUTPUT_DIR / "figures"
for d in [OUTPUT_DIR, MODEL_DIR, FIG_DIR, _LOG_BASE]:
    d.mkdir(parents=True, exist_ok=True)

LINEAGE_INPUTS = {
    "epithelial": {
        "h5ad"           : Path("/home/h2048/data/py/0317/epithelial_v2_7_HOTFIX/SELF/"
                                "epithelial_scanvi_v2_7_HOTFIX_SELF_final.h5ad"),
        "label"          : "cell_type_L3",
        "name"           : "Epithelial",
        "batch_key_local": "sample",
    },
    "tcell": {
        "h5ad"           : Path("/home/h2048/data/py/0318/tnk_subcluster_retrain/"
                                "adata_tnk_scanvi_ref_retrain_v1_2.h5ad"),
        "label"          : "cell_type_L3",
        "name"           : "TNK",
        "batch_key_local": "sample",
    },
    "myeloid": {
        "h5ad"           : Path("/home/h2048/data/py/0209/myeloid_validation_optimized/"
                                "adata_myeloid_refined_FINAL.h5ad"),
        "label"          : "cell_type_L3",
        "name"           : "Myeloid",
        "batch_key_local": "sample",
    },
    "bcell": {
        "h5ad"           : Path("/home/h2048/data/py/0203/bcell_scarches_v4_1/results/"
                                "scarches_package/bcell_reference_20260203.h5ad"),
        "label"          : "cell_type_L3",
        "name"           : "Bcell",
        "batch_key_local": "sample",
    },
    "stromal": {
        "h5ad"           : Path("/home/h2048/data/py/0308/stromal_reintegration_v1_3/"
                                "stromal_reintegrated_scvi_scanvi_v1_3.h5ad"),
        "label"          : "cell_type_L3",
        "name"           : "Stromal",
        "batch_key_local": "sample",
    },
}

LABEL_TIER                 = "L1"
BATCH_KEY                  = "sample"
LABELS_KEY                 = "scanvi_label"
UNLABELED_CATEGORY         = "Unknown"
N_HVG                      = 6000
ALLOW_UNKNOWN              = False
DROP_UNKNOWN_IF_DISALLOWED = True
MIN_CELLS_PER_BATCH        = 3
MIN_CELLS_PER_L1_CLASS     = 10

OVERLAP_WHITELIST = set()

# =====================================================================
# [v1.10] Force biological markers into HVG — copied from reference
# pipeline v2.5.2. Ensures lineage anchor genes are always in HVG even
# if their per-batch variance ranks below the N_HVG threshold.
# Swap strategy: remove lowest-variance non-forced HVGs to make room.
# =====================================================================
FORCE_MARKERS_IN_HVG        = True
FORCE_MARKERS_CASE_INSENSITIVE = True

# Lineage anchor markers (all 5 lineages covered)
FORCED_MARKERS_LINEAGE = [
    # Epithelial
    "FXYD3","EPCAM","ELF3","IGFBP2","SERPINF1","TSPAN1","SCGB1A1",
    "AGER","SFTPC","FOXJ1","KRT5","MUC5B","KRT8",
    # Immune general / Lymphoid
    "CD53","PTPRC","CORO1A","ISG20","CCL5",
    # B cell
    "MS4A1","TNFRSF17","CD19","CD79A","SDC1",
    # T/NK
    "CD40LG","TNFRSF25","CD28","CD4","CD3D","CD3E","CD2","TRBC2",
    "CD8A","CD8B","TRGC2",
    # Myeloid
    "FCER1G","C1ORF162","CLEC7A","CD1C","CD86","CD14","XCR1","HLA-DRA",
    # Stromal / Fibroblast / SMC
    "COL1A2","DCN","MFAP4","LUM","COL6A3","CFD","COL1A1","PDGFRA",
    "MXRA8","NBL1","VCAN","LEPR","MYH11","TINAGL1","PLN","DES",
    "ACTA2","CNN1","TAGLN",
    # Endothelial
    "CLDN5","ECSCR","CLEC14A","VWF","PECAM1","ACKR1","PTPRB","PDE2A",
    "PLAT","GJA5","SPARCL1","AQP1","RNASE1","MMRN1","CCL21","TFF3",
]
# State markers (disabled by default — avoid cell-cycle signal dominance)
FORCED_MARKERS_STATE   = ["MKI67","TOP2A","TK1","CENPW"]
INCLUDE_STATE_MARKERS  = False
FORCED_MARKERS = FORCED_MARKERS_LINEAGE + (FORCED_MARKERS_STATE if INCLUDE_STATE_MARKERS else [])

# =====================================================================
# [v1.9] scVI Hyperparameters — fully aligned with reference pipeline
# v2.5.2 (allcells_scvi_celltypist_scanvi_pipeline_v2_5_2_PRODUCTION.py)
# =====================================================================

# [v1.9-3] N_LATENT: 175 -> 150  (reference standard)
# [v1.9-2] dropout_rate: 0.1 -> 0.2
# [v1.9-2] use_layer_norm: 'both'  (was not set, default='none')
# [v1.9-2] use_batch_norm: 'none'  (was not set, default='both' — contradicts layer_norm)
# [v1.9-2] early_stopping_patience: 45 (was default 45 but now explicit)
SCVI_N_LATENT               = 150
SCVI_N_LAYERS               = 2
SCVI_N_HIDDEN               = 128
SCVI_DROPOUT_RATE           = 0.2
SCVI_MAX_EPOCHS             = 400
SCVI_BATCH_SIZE             = 256
SCVI_LEARNING_RATE          = 1e-3
SCVI_EARLY_STOPPING         = True
SCVI_EARLY_STOPPING_PATIENCE = 45
SCVI_ENCODE_COVARIATES      = True
SCVI_USE_LAYER_NORM         = "both"
SCVI_USE_BATCH_NORM         = "none"

# =====================================================================
# [v1.9] scANVI Hyperparameters
# [v1.9-4] lr: 1e-3 -> 5e-4  (fine-tuning rate, not pretraining rate)
# [v1.9-4] n_samples_per_label=2000  (class balancing for rare subtypes)
# =====================================================================
SCANVI_MAX_EPOCHS             = 200
SCANVI_BATCH_SIZE             = 256
SCANVI_LEARNING_RATE          = 5e-4
SCANVI_EARLY_STOPPING_PATIENCE = 30
# [v1.10-4] n_samples_per_label removed: 5-class L1 task, resampling not needed

# =====================================================================
# [v1.9-5] Covariate registration — must match what Cell 9 computes
# continuous: pct_counts_mt, stress_score, S_score, G2M_score
# categorical: none (no tissue_key in this combined pipeline)
# =====================================================================
CONTINUOUS_COVARIATES  = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
CATEGORICAL_COVARIATES = []

OUTPUT_H5AD   = OUTPUT_DIR / f"allcells_combined_{DATE_TAG}_v{VERSION}.h5ad"
PIPELINE_NAME = f"allcells_combined_v{VERSION}"

for key, cfg in LINEAGE_INPUTS.items():
    p = cfg["h5ad"]
    if not p.exists():
        hits = sorted(Path("/home/h2048/data").rglob(p.name))[:5]
        raise FileNotFoundError(f"[X] {key}: {p}\n    Matches: {hits}")

print(f"Output     : {OUTPUT_H5AD.name}")
print(f"Label tier : {LABEL_TIER}  (major lineage supervision)")
print(f"Batch key  : {BATCH_KEY}  (sample-level unified)")
print(f"ALLOW_UNKNOWN={ALLOW_UNKNOWN}  DROP_UNKNOWN={DROP_UNKNOWN_IF_DISALLOWED}")
print(f"MIN_CELLS_PER_BATCH={MIN_CELLS_PER_BATCH}  MIN_CELLS_PER_L1_CLASS={MIN_CELLS_PER_L1_CLASS}")
print(f"OVERLAP_WHITELIST={OVERLAP_WHITELIST}")
print()
print("[v1.9] scVI hyperparameters:")
print(f"  N_LATENT={SCVI_N_LATENT}  N_LAYERS={SCVI_N_LAYERS}  N_HIDDEN={SCVI_N_HIDDEN}")
print(f"  dropout_rate={SCVI_DROPOUT_RATE}  encode_covariates={SCVI_ENCODE_COVARIATES}")
print(f"  use_layer_norm={SCVI_USE_LAYER_NORM}  use_batch_norm={SCVI_USE_BATCH_NORM}")
print(f"  early_stopping_patience={SCVI_EARLY_STOPPING_PATIENCE}")
print()
print("[v1.9] scANVI hyperparameters:")
print(f"  lr={SCANVI_LEARNING_RATE}  (n_samples_per_label removed in v1.10)")
print(f"  early_stopping_patience={SCANVI_EARLY_STOPPING_PATIENCE}")
print()
print(f"[v1.9] Covariates: continuous={CONTINUOUS_COVARIATES}  categorical={CATEGORICAL_COVARIATES}")
for key, cfg in LINEAGE_INPUTS.items():
    print(f"  {key:12s}: {cfg['h5ad'].name}  ->  L1='{cfg['name']}'")

# %% [markdown]
# ## Cell 3 — Per-Lineage Label Maps + Builder Functions

# %%
_TNK_L3_RENAME = {"CD4 Naive/TCM": "CD4 Naive"}

_TNK_L2_MAP = {
    "CD8 Trm"      : "CD8 T cells", "CD8 Tem"      : "CD8 T cells",
    "CD8 Temra"    : "CD8 T cells", "CD8 Teff"     : "CD8 T cells",
    "CD8 Naive"    : "CD8 T cells",
    "NK"           : "NK cells",    "NK Exhausted"  : "NK cells",
    "CD4 Trm"      : "CD4 T cells", "CD4 Tcm"      : "CD4 T cells",
    "CD4 Tfh"      : "CD4 T cells", "CD4 Tfr"      : "CD4 T cells",
    "CD4 Treg"     : "CD4 T cells", "CD4 Th17"     : "CD4 T cells",
    "CD4 Th1"      : "CD4 T cells", "CD4 Naive"    : "CD4 T cells",
    "gdT"          : "gdT cells",   "MAIT"          : "MAIT cells",
    "ILC3"         : "ILC",
}

_MYELOID_L3_MERGE = {
    "Alveolar Mph CCL3+"        : "Alveolar Mph",
    "Alveolar Mph proliferating": "Alveolar Mph",
    "Alveolar Mph MT-positive"  : "Alveolar Mph",
}
_MYELOID_FILTER_SET = {
    "AT2", "Basal resting", "Club (nasal)", "Club (non-nasal)", "Goblet (nasal)",
    "Adventitial fibroblasts", "Alveolar fibroblasts", "Smooth muscle",
    "Lymphatic EC mature", "B cells", "Plasma cells",
    "CD4 T cells", "CD8 T cells", "Hematopoietic stem cells",
    "Low-quality Interstitial macrophages",
}
_MYELOID_FINEST_L2 = {
    "Alveolar macrophages"          : "Alveolar_Macrophage",
    "Alveolar Mph"                  : "Alveolar_Macrophage",
    "Alveolar Mph CCL3+"            : "Alveolar_Macrophage",
    "Alveolar Mph proliferating"    : "Alveolar_Macrophage",
    "Alveolar Mph MT-positive"      : "Alveolar_Macrophage",
    "Monocyte-derived Mph"          : "Interstitial_Macrophage",
    "Interstitial Mph perivascular" : "Interstitial_Macrophage",
    "Classical monocytes"           : "Monocyte",
    "Non-classical monocytes"       : "Monocyte",
    "DC2"                           : "DC",
    "DC1"                           : "DC",
    "Plasmacytoid DCs"              : "DC",
    "Mast cells"                    : "Mast_cell",
}
_MYELOID_REFINED_L2 = {
    "Resident Alveolar macrophages"                : "Alveolar_Macrophage",
    "Resting Alveolar macrophages"                 : "Alveolar_Macrophage",
    "Inflammatory Interstitial macrophages"        : "Interstitial_Macrophage",
    "M2-like Interstitial macrophages"             : "Interstitial_Macrophage",
    "Atypically activated Interstitial macrophages": "Interstitial_Macrophage",
    "CD163L1+ Interstitial macrophages"            : "Interstitial_Macrophage",
    "Immunoregulatory Interstitial macrophages"    : "Interstitial_Macrophage",
    "Neutrophils"                                  : "Neutrophil",
    "Typical Classical monocytes"                  : "Monocyte",
    "Inflammatory Classical monocytes"             : "Monocyte",
    "Conventional cDC2"                            : "DC",
    "Langerhans-like cDC2"                         : "DC",
    "Migratory DCs"                                : "DC",
    "pDC"                                          : "DC",
    "Mast cells"                                   : "Mast_cell",
}
_MYELOID_ALL_L2 = {**_MYELOID_FINEST_L2, **_MYELOID_REFINED_L2}

_BCELL_L3_MERGE = {"IGHEplus_Atypical_Memory_B": "Atypical_Memory_B"}
_BCELL_L2_MAP   = {
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B_Light_Zone_Centrocyte"        : "GC_B",
    "GC_B_Transitional"                 : "GC_B",
    "Plasma_IgA"                        : "Plasma",
    "Plasma_IgG"                        : "Plasma",
    "Atypical_Memory_B"                 : "Atypical_Memory_B",
    "Memory_B"                          : "Memory_B",
    "Naive_B"                           : "Naive_B",
}


def build_epithelial_labels(adata_lin):
    assert "scanvi_fine_pred"  in adata_lin.obs.columns
    assert "scanvi_major_pred" in adata_lin.obs.columns
    l3 = adata_lin.obs["scanvi_fine_pred"].astype(str).copy()
    l2 = adata_lin.obs["scanvi_major_pred"].astype(str).copy()
    l3 = l3.replace({"Ionocyte_Brush": "Ionocyte"})
    n_overrides = 0
    if "ann_finest_level" in adata_lin.obs.columns:
        mask = adata_lin.obs["ann_finest_level"].astype(str) == "Hillock-like"
        l3[mask] = "Hillock-like"; l2[mask] = "Rare_Specialized"
        n_overrides += int(mask.sum())
    if "ann_coarse_for_GWAS_and_modeling" in adata_lin.obs.columns:
        coarse = adata_lin.obs["ann_coarse_for_GWAS_and_modeling"].astype(str)
        for rare in ["Ionocyte", "Neuroendocrine"]:
            mask = coarse == rare
            l3[mask] = rare; l2[mask] = "Rare_Specialized"
            n_overrides += int(mask.sum())
    adata_lin.obs["cell_type_L3"] = l3
    adata_lin.obs["cell_type_L2"] = l2
    print(f"    -> L3: {l3.nunique()} classes  L2: {l2.nunique()} classes  overrides: {n_overrides:,}")


def build_tcell_labels(adata_lin):
    assert "scanvi_label_refined" in adata_lin.obs.columns
    l3 = adata_lin.obs["scanvi_label_refined"].astype(str).copy().replace(_TNK_L3_RENAME)
    unmapped = set(l3.unique()) - set(_TNK_L2_MAP.keys())
    if unmapped:
        raise ValueError(f"[X] tcell: unmapped L3 -> L2: {sorted(unmapped)}")
    adata_lin.obs["cell_type_L3"] = l3
    adata_lin.obs["cell_type_L2"] = l3.map(_TNK_L2_MAP)
    print(f"    -> L3: {l3.nunique()} classes  L2: {adata_lin.obs['cell_type_L2'].nunique()} classes")


def build_myeloid_labels(adata_lin):
    assert "ann_finest_level"     in adata_lin.obs.columns
    assert "cell_type_L3_refined" in adata_lin.obs.columns
    finest  = adata_lin.obs["ann_finest_level"].astype(str)
    refined = adata_lin.obs["cell_type_L3_refined"].astype(str)
    _bad    = {"Unknown", "nan", "NaN", "None", ""}
    use_finest = ~finest.isin(_bad)
    l3 = finest.where(use_finest, refined)
    l3 = l3.replace({"nan": "Unknown", "NaN": "Unknown", "None": "Unknown"})
    filter_mask = l3.isin(_MYELOID_FILTER_SET)
    if filter_mask.any():
        l3[filter_mask] = "Unknown"
    l3 = l3.replace(_MYELOID_L3_MERGE)
    non_unknown = set(l3[l3 != "Unknown"].unique())
    unmapped    = non_unknown - set(_MYELOID_ALL_L2.keys())
    if unmapped:
        raise ValueError(f"[X] myeloid: unmapped L3 -> L2: {sorted(unmapped)}")
    l2 = l3.map(_MYELOID_ALL_L2).fillna("Unknown")
    adata_lin.obs["cell_type_L3"] = l3
    adata_lin.obs["cell_type_L2"] = l2
    print(f"    -> L3: {l3.nunique()} classes  L2: {l2.nunique()} classes  Unknown: {int((l3=='Unknown').sum()):,}")


def build_bcell_labels(adata_lin):
    assert "cell_type_scanvi_pred" in adata_lin.obs.columns
    l3 = adata_lin.obs["cell_type_scanvi_pred"].astype(str).copy().replace(_BCELL_L3_MERGE)
    unmapped = set(l3.unique()) - set(_BCELL_L2_MAP.keys())
    if unmapped:
        raise ValueError(f"[X] bcell: unmapped L3 -> L2: {sorted(unmapped)}")
    adata_lin.obs["cell_type_L3"] = l3
    adata_lin.obs["cell_type_L2"] = l3.map(_BCELL_L2_MAP)
    print(f"    -> L3: {l3.nunique()} classes  L2: {adata_lin.obs['cell_type_L2'].nunique()} classes")


def build_stromal_labels(adata_lin):
    assert "cell_type_scanvi_pred" in adata_lin.obs.columns
    assert "cell_type_L2"          in adata_lin.obs.columns
    adata_lin.obs["cell_type_L3"] = adata_lin.obs["cell_type_scanvi_pred"].astype(str)
    adata_lin.obs["cell_type_L2"] = adata_lin.obs["cell_type_L2"].astype(str)
    print(f"    -> L3: {adata_lin.obs['cell_type_L3'].nunique()} classes  L2: {adata_lin.obs['cell_type_L2'].nunique()} classes")


_LABEL_BUILDERS = {
    "epithelial": build_epithelial_labels,
    "tcell"     : build_tcell_labels,
    "myeloid"   : build_myeloid_labels,
    "bcell"     : build_bcell_labels,
    "stromal"   : build_stromal_labels,
}

print("[OK] Label maps and builder functions defined")

# %% [markdown]
# ## Cell 4 — Structural Normalization + Counts Resolution Helpers

# %%
_SCVI_INTERNAL_OBS_COLS = {
    "_scvi_batch", "_scvi_labels",
    "_scvi_extra_categorical_covs", "_scvi_extra_continuous_covs",
}
_HVG_AUX_VAR_COLS = {
    "highly_variable", "highly_variable_rank", "means", "variances",
    "variances_norm", "dispersions", "dispersions_norm", "highly_variable_nbatches",
}


def normalize_lineage_structure(adata_lin, key, local_bk, global_bk):
    print(f"  [normalize] {key}")
    if "counts" not in adata_lin.layers:
        if "raw_counts" in adata_lin.layers:
            adata_lin.layers["counts"] = adata_lin.layers["raw_counts"]
        elif adata_lin.raw is not None:
            raw_idx = pd.Index(adata_lin.raw.var_names)
            cur_idx = pd.Index(adata_lin.var_names)
            shared  = raw_idx.intersection(cur_idx)
            if len(shared) == len(cur_idx):
                pos = raw_idx.get_indexer(cur_idx)
                adata_lin.layers["counts"] = sparse.csr_matrix(
                    adata_lin.raw.X[:, pos], dtype=np.float32)
            else:
                raise ValueError(f"[X] {key}: counts absent")
        else:
            raise ValueError(f"[X] {key}: no counts source")
    if not sparse.isspmatrix_csr(adata_lin.layers["counts"]):
        adata_lin.layers["counts"] = sparse.csr_matrix(
            adata_lin.layers["counts"], dtype=np.float32)
    elif adata_lin.layers["counts"].dtype != np.float32:
        adata_lin.layers["counts"] = adata_lin.layers["counts"].astype(np.float32)
    _d = adata_lin.layers["counts"].data
    if len(_d) > 0:
        if not np.isfinite(_d).all():
            raise ValueError(f"[X] {key}: non-finite in counts")
        if np.any(_d < 0):
            raise ValueError(f"[X] {key}: negative in counts")
    if global_bk not in adata_lin.obs.columns:
        if local_bk in adata_lin.obs.columns:
            adata_lin.obs[global_bk] = adata_lin.obs[local_bk].astype(object)
        else:
            cands = [c for c in adata_lin.obs.columns
                     if any(k in c.lower() for k in ("sample","batch","dataset"))]
            raise ValueError(f"[X] {key}: batch col not found. Candidates: {cands}")
    else:
        adata_lin.obs[global_bk] = adata_lin.obs[global_bk].astype(object)
    drop_obs = [c for c in adata_lin.obs.columns if c in _SCVI_INTERNAL_OBS_COLS]
    if drop_obs:
        adata_lin.obs.drop(columns=drop_obs, inplace=True)
    for col in adata_lin.obs.select_dtypes(include=["category"]).columns:
        cats = adata_lin.obs[col].cat.categories
        if hasattr(cats.dtype, "name") and cats.dtype.name in ("string","StringDtype"):
            adata_lin.obs[col] = adata_lin.obs[col].cat.rename_categories(cats.astype(object))
    drop_var = [c for c in adata_lin.var.columns if c in _HVG_AUX_VAR_COLS]
    if drop_var:
        adata_lin.var.drop(columns=drop_var, inplace=True)
    if len(adata_lin.obsm):
        adata_lin.obsm.clear()
    print(f"    counts: CSR float32 {adata_lin.layers['counts'].shape}  batch: {adata_lin.obs[global_bk].nunique()} levels")


def resolve_fullgene_counts_matrix(adata_lin, key):
    candidates = []
    if adata_lin.raw is not None:
        candidates.append(dict(source="adata.raw.X", X=adata_lin.raw.X,
                               var=adata_lin.raw.var.copy(),
                               n_vars=adata_lin.raw.n_vars, rank=0))
    if "raw_counts" in adata_lin.layers:
        candidates.append(dict(source="layers['raw_counts']", X=adata_lin.layers["raw_counts"],
                               var=adata_lin.var.copy(), n_vars=adata_lin.n_vars, rank=1))
    if "counts" in adata_lin.layers:
        candidates.append(dict(source="layers['counts']", X=adata_lin.layers["counts"],
                               var=adata_lin.var.copy(), n_vars=adata_lin.n_vars, rank=2))
    if not candidates:
        raise ValueError(f"[X] {key}: no counts source")
    candidates.sort(key=lambda c: (-c["n_vars"], c["rank"]))
    skipped_hard   = []
    skipped_nonint = []
    nonint_fallback = None
    for cand in candidates:
        X, src = cand["X"], cand["source"]
        _d = X.data if sparse.issparse(X) else np.asarray(X).ravel()
        if len(_d) > 0 and not np.isfinite(_d).all():
            skipped_hard.append(f"    [SKIP] {src}: non-finite"); continue
        if len(_d) > 0 and np.any(_d < 0):
            skipped_hard.append(f"    [SKIP] {src}: negative"); continue
        if len(_d) > 0:
            _s   = _d[:min(50000, len(_d))]
            frac = np.abs(_s - np.round(_s))
            if np.nanmax(frac) > 1e-3:
                skipped_nonint.append(f"    [DEFER] {src} ({cand['n_vars']:,} genes): non-integer")
                if nonint_fallback is None:
                    nonint_fallback = cand
                continue
        if sparse.issparse(X):
            if not sparse.isspmatrix_csr(X): X = X.tocsr()
            if X.dtype != np.float32:        X = X.astype(np.float32)
        else:
            X = sparse.csr_matrix(np.asarray(X, dtype=np.float32))
        print(f"  [OK] {key}: {src} | genes={cand['n_vars']:,}")
        return X, src, cand["var"]
    if nonint_fallback is not None:
        fb = nonint_fallback
        print(f"  [WARNING] {key}: no integer counts found, fallback to {fb['source']}")
        X = fb["X"]
        if sparse.issparse(X):
            if not sparse.isspmatrix_csr(X): X = X.tocsr()
            if X.dtype != np.float32:        X = X.astype(np.float32)
        else:
            X = sparse.csr_matrix(np.asarray(X, dtype=np.float32))
        return X, fb["source"] + " [nonint-fallback]", fb["var"]
    raise ValueError(f"[X] {key}: all candidates failed.\n" + "\n".join(skipped_hard + skipped_nonint))


def _fmt_float(val, fallback="not_set"):
    return f"{val:.4f}" if isinstance(val, (float, np.floating)) else str(fallback)


# ------------------------------------------------------------------ force marker helpers
def detect_gene_symbol_column(adata, verbose=True):
    """
    Detect the gene symbol column in adata.var with automatic fallback.
    Tries 12 common column names; falls back to var_names.
    After outer-join concat, var.columns from 5 lineages may differ;
    this function handles all known naming conventions.
    """
    CANDIDATE_COLUMNS = [
        "symbol_base", "gene_symbol", "gene_symbols", "Symbol",
        "symbol", "feature_name", "gene_name", "gene_names",
        "gene", "genes", "gene_id", "gene_ids",
    ]
    for col in CANDIDATE_COLUMNS:
        if col in adata.var.columns:
            gene_symbols = adata.var[col].astype(str)
            non_empty = gene_symbols.notna() & (gene_symbols != "") & (gene_symbols != "nan")
            if non_empty.sum() >= len(gene_symbols) * 0.9:
                if verbose:
                    print(f"  Gene symbols from: adata.var['{col}']  "
                          f"valid={non_empty.sum():,}/{len(gene_symbols):,}")
                return col, gene_symbols
    if verbose:
        print("  [WARNING] No gene symbol column found; falling back to var_names")
    return None, pd.Series(adata.var_names.values, index=adata.var_names)


def force_include_markers_in_hvg(adata, n_top_genes, forced_markers,
                                  symbol_col=None, case_insensitive=True,
                                  hvg_col="highly_variable"):
    """
    Force specified marker genes into HVG selection.
    For markers not yet in HVG: swap out lowest-variance non-forced HVGs.
    Final HVG count is kept at n_top_genes.
    """
    if hvg_col not in adata.var.columns:
        raise ValueError(f"HVG column '{hvg_col}' not found in adata.var")
    if symbol_col is None:
        symbol_col, _syms = detect_gene_symbol_column(adata, verbose=False)
        var_symbols = adata.var_names.values if symbol_col is None else adata.var[symbol_col].values
    else:
        if symbol_col not in adata.var.columns:
            raise ValueError(f"Symbol column '{symbol_col}' not found")
        var_symbols = adata.var[symbol_col].values
    if case_insensitive:
        forced_set  = {m.upper() for m in forced_markers}
        var_symbols = np.array([str(s).upper() for s in var_symbols])
    else:
        forced_set  = set(forced_markers)
        var_symbols = np.array([str(s) for s in var_symbols])
    forced_mask    = np.array([s in forced_set for s in var_symbols])
    forced_indices = np.flatnonzero(forced_mask)
    if len(forced_indices) == 0:
        return {"n_requested": len(forced_markers), "n_found_in_data": 0,
                "n_already_in_hvg": 0, "n_newly_added": 0,
                "n_hvg_final": int(adata.var[hvg_col].sum())}
    hvg_mask     = adata.var[hvg_col].values.copy()
    n_already_in = int(hvg_mask[forced_indices].sum())
    n_to_add     = len(forced_indices) - n_already_in
    if n_to_add == 0:
        return {"n_requested": len(forced_markers), "n_found_in_data": len(forced_indices),
                "n_already_in_hvg": n_already_in, "n_newly_added": 0,
                "n_hvg_final": int(hvg_mask.sum())}
    var_col = next((vc for vc in ["variances_norm","variances","dispersions_norm","dispersions"]
                    if vc in adata.var.columns), None)
    if var_col is None:
        hvg_mask[forced_mask & ~hvg_mask] = True
    else:
        variances = adata.var[var_col].values.copy()
        non_forced_hvg = np.array([i for i in np.flatnonzero(hvg_mask) if not forced_mask[i]])
        if len(non_forced_hvg) >= n_to_add:
            remove_idx = non_forced_hvg[np.argsort(variances[non_forced_hvg])[:n_to_add]]
            hvg_mask[remove_idx] = False
        hvg_mask[forced_mask & ~hvg_mask] = True
    adata.var[hvg_col] = hvg_mask
    # Trim back if over budget
    n_now = int(hvg_mask.sum())
    if n_now > n_top_genes and var_col is not None:
        non_forced_hvg2 = np.array([i for i in np.flatnonzero(hvg_mask) if not forced_mask[i]])
        if len(non_forced_hvg2) > 0:
            trim = non_forced_hvg2[np.argsort(variances[non_forced_hvg2])[:n_now - n_top_genes]]
            hvg_mask[trim] = False
            adata.var[hvg_col] = hvg_mask
    return {"n_requested": len(forced_markers), "n_found_in_data": len(forced_indices),
            "n_already_in_hvg": n_already_in, "n_newly_added": n_to_add,
            "n_hvg_final": int(adata.var[hvg_col].sum())}


print("[OK] Structural helpers defined")

# %% [markdown]
# ## Cell 5 — Load Lineages + Build L3/L2 Labels + Normalize

# %%
print("\n" + "="*70)
print("[STEP 5] Loading, building labels, normalizing")
print("="*70)

lineage_adatas = []
label_sets     = {}
gene_sets      = {}

for key, cfg in LINEAGE_INPUTS.items():
    print(f"\n[{key}] {cfg['h5ad'].name}")
    adata_lin = sc.read_h5ad(cfg["h5ad"])
    print(f"  shape : {adata_lin.n_obs:,} x {adata_lin.n_vars:,}")

    _LABEL_BUILDERS[key](adata_lin)

    normalize_lineage_structure(
        adata_lin, key,
        local_bk  = cfg.get("batch_key_local", BATCH_KEY),
        global_bk = BATCH_KEY,
    )

    counts_X, count_source, var_df = resolve_fullgene_counts_matrix(adata_lin, key)

    adata_full = sc.AnnData(X=counts_X, obs=adata_lin.obs.copy(), var=var_df)
    adata_full.layers["counts"]             = adata_full.X
    adata_full.obs["lineage_source"]        = cfg["name"]
    adata_full.obs["cell_type_original"]    = cfg["name"]
    adata_full.obs["cell_type_L2_original"] = (
        adata_lin.obs["cell_type_L2"].astype(object).fillna("Unknown")
    )
    adata_full.obs["cell_type_L3_original"] = (
        adata_lin.obs["cell_type_L3"].astype(object).fillna(UNLABELED_CATEGORY)
    )

    label_sets[key] = set(adata_full.obs["cell_type_original"].unique())
    gene_sets[key]  = set(adata_full.var_names.tolist())

    print(f"  -> {adata_full.n_obs:,} x {adata_full.n_vars:,} | src='{count_source}'")

    lineage_adatas.append(adata_full)
    del adata_lin; gc.collect()

print("\n[OK] All lineages loaded")

# %% [markdown]
# ## Cell 6 — Pre-concat Audits

# %%
print("\n" + "="*70)
print("[STEP 6] Pre-concat audits")
print("="*70)

# [v1.10-1] L2/L3 cross-lineage overlap audit.
# L1 labels are lineage names (Epithelial / TNK / ...) — always unique by construction.
# The biologically relevant check is whether any L2 or L3 label string appears
# in >1 lineage, which would indicate label naming collisions that could confuse
# downstream interpretation (even though scANVI only sees L1 here).
# Action: log only; no automatic prefixing (L1 training is not affected).
print("\n[6a] Cross-lineage L2/L3 label overlap audit")

l2_sets = {key: set(ad.obs["cell_type_L2_original"].astype(str).unique()) - {"Unknown", "nan"}
           for key, ad in zip(LINEAGE_INPUTS.keys(), lineage_adatas)}
l3_sets = {key: set(ad.obs["cell_type_L3_original"].astype(str).unique()) - {UNLABELED_CATEGORY, "nan"}
           for key, ad in zip(LINEAGE_INPUTS.keys(), lineage_adatas)}

all_keys = list(l2_sets.keys())
l2_overlap_found = False
l3_overlap_found = False
for i in range(len(all_keys)):
    for j in range(i + 1, len(all_keys)):
        ka, kb = all_keys[i], all_keys[j]
        l2_shared = l2_sets[ka] & l2_sets[kb]
        l3_shared = l3_sets[ka] & l3_sets[kb]
        if l2_shared:
            print(f"  [L2 OVERLAP] {ka} x {kb}: {sorted(l2_shared)}")
            l2_overlap_found = True
        if l3_shared:
            print(f"  [L3 OVERLAP] {ka} x {kb}: {sorted(l3_shared)[:10]}" +
                  (f" ... (+{len(l3_shared)-10} more)" if len(l3_shared) > 10 else ""))
            l3_overlap_found = True

if not l2_overlap_found and not l3_overlap_found:
    print("  [OK] No cross-lineage L2/L3 label string collisions")
else:
    print("  [NOTE] Overlapping label strings logged above. scANVI training uses L1 only")
    print("         and is not affected. Review for downstream L2/L3 interpretation.")

print("\n[6b] Gene universe consistency")
gene_sets = {key: set(ad.var_names.tolist())
             for key, ad in zip(LINEAGE_INPUTS.keys(), lineage_adatas)}
for k in gene_sets: print(f"  {k:12s}: {len(gene_sets[k]):,} genes")
intersection = set.intersection(*gene_sets.values())
union        = set.union(*gene_sets.values())
print(f"  Inner: {len(intersection):,}  Outer: {len(union):,}  Coverage: {len(intersection)/max(len(union),1):.3f}")

# %% [markdown]
# ## Cell 7 — Concatenate (outer join + zero-fill)

# %%
print("\n" + "="*70)
print("[STEP 7] Concatenating with outer join")
print("="*70)

adata = sc.concat(
    lineage_adatas,
    join         = "outer",
    merge        = "unique",
    uns_merge    = "unique",
    label        = "lineage_concat_key",
    keys         = list(LINEAGE_INPUTS.keys()),
    index_unique = "-",
    fill_value   = 0,
)
del lineage_adatas; gc.collect()

print(f"Combined: {adata.n_obs:,} x {adata.n_vars:,}")
print(adata.obs["lineage_source"].value_counts().to_string())

if "counts" in adata.layers:
    X_counts = adata.layers["counts"]
    if not sparse.issparse(X_counts):
        X_arr = np.asarray(X_counts, dtype=np.float32)
        X_arr = np.where(np.isfinite(X_arr), X_arr, 0.0)
        X_arr = np.where(X_arr < 0, 0.0, X_arr)
        adata.layers["counts"] = sparse.csr_matrix(X_arr, dtype=np.float32)
        del X_arr; gc.collect()
    else:
        X_csr = sparse.csr_matrix(X_counts, dtype=np.float32, copy=False)
        if not np.isfinite(X_csr.data).all():
            X_csr = X_csr.copy()
            X_csr.data[~np.isfinite(X_csr.data)] = 0.0
            X_csr.eliminate_zeros()
        if np.any(X_csr.data < 0):
            raise ValueError("[X] Negative values in counts after outer join")
        adata.layers["counts"] = X_csr
else:
    raise ValueError("[X] counts layer missing after concat")

_d = adata.layers["counts"].data
assert np.isfinite(_d).all() and np.all(_d >= 0)
del _d; gc.collect()
print(f"[OK] counts: CSR float32 {adata.layers['counts'].shape}")

# %% [markdown]
# ## Cell 8 — Build `scanvi_label` (L1) + L2 + Normalise + Rare Filter

# %%
print("\n" + "="*70)
print("[STEP 8] Unified labels + normalisation + rare class filter")
print("="*70)

def _norm_label(s):
    return str(s).replace(" ", "_").replace("-", "_")

adata.obs["cell_type_original"]    = adata.obs["cell_type_original"].astype(str).map(_norm_label)
adata.obs["cell_type_L2_original"] = adata.obs["cell_type_L2_original"].astype(str).map(_norm_label)
adata.obs["cell_type_L3_original"] = adata.obs["cell_type_L3_original"].astype(str).map(_norm_label)

adata.obs["cell_type_input_l1"] = adata.obs["cell_type_original"].astype(object)
adata.obs["cell_type_input_l2"] = adata.obs["cell_type_L2_original"].astype(object)
adata.obs["cell_type_input_l3"] = adata.obs["cell_type_L3_original"].astype(object)

adata.obs[LABELS_KEY] = pd.Categorical(
    adata.obs["cell_type_original"].astype(object).fillna(UNLABELED_CATEGORY)
)
if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
    adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

adata.obs["cell_type_L2"] = pd.Categorical(
    adata.obs["cell_type_L2_original"].astype(object).fillna("Unknown")
)

n_labeled = int((adata.obs[LABELS_KEY] != UNLABELED_CATEGORY).sum())
n_unknown  = int((adata.obs[LABELS_KEY] == UNLABELED_CATEGORY).sum())
print(f"L1 Labeled: {n_labeled:,} ({n_labeled/adata.n_obs*100:.1f}%)")
print(f"L1 Unknown: {n_unknown:,} ({n_unknown/adata.n_obs*100:.1f}%)")
print(adata.obs[LABELS_KEY].value_counts().to_string())

label_counts = adata.obs[LABELS_KEY].value_counts()
rare_labels  = label_counts[
    (label_counts < MIN_CELLS_PER_L1_CLASS) &
    (label_counts.index != UNLABELED_CATEGORY)
].index.tolist()

if rare_labels:
    adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].astype(object)
    adata.obs.loc[adata.obs[LABELS_KEY].isin(rare_labels), LABELS_KEY] = UNLABELED_CATEGORY
    adata.obs[LABELS_KEY] = pd.Categorical(adata.obs[LABELS_KEY])
    if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
        adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])
    n_labeled = int((adata.obs[LABELS_KEY] != UNLABELED_CATEGORY).sum())
    n_unknown  = int((adata.obs[LABELS_KEY] == UNLABELED_CATEGORY).sum())

dropped_unknown = 0
if not ALLOW_UNKNOWN and n_unknown > 0:
    if DROP_UNKNOWN_IF_DISALLOWED:
        dropped_unknown = n_unknown
        adata = adata[adata.obs[LABELS_KEY] != UNLABELED_CATEGORY].copy()
        adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.remove_unused_categories()
        n_labeled = adata.n_obs; n_unknown = 0
        print(f"[OK] After drop: {adata.n_obs:,} cells")
    else:
        assert False, f"[X] {n_unknown:,} Unknown cells."

adata.uns["dropped_unknown_cells"] = int(dropped_unknown)
adata.uns["allow_unknown"]         = bool(ALLOW_UNKNOWN)
adata.uns["label_tier"]            = LABEL_TIER
print("[OK] Labels built")

# %% [markdown]
# ## Cell 9 — Covariates: MT% + Stress Score + Cell Cycle Scores
# 
# **[v1.9-1] CRITICAL:** All four covariates must be computed on the **full gene space** before
# any HVG subsetting. The stress score uses log1p(CPM) mean across signature genes (no [0,1]
# normalisation). Cell cycle scores via `sc.tl.score_genes_cell_cycle`. These are then registered
# as `continuous_covariate_keys` in `SCVI.setup_anndata` (Cell 15).

# %%
print("\n" + "="*70)
print("[STEP 9] Covariates (MT% + stress_score + S_score + G2M_score)")
print("="*70)

# ------------------------------------------------------------------ MT%
mt_mask = adata.var_names.str.startswith("MT-")
n_mt    = int(mt_mask.sum())
_total  = np.array(adata.layers["counts"].sum(axis=1)).ravel().astype(np.float64)
if n_mt > 0:
    _mt = np.array(adata.layers["counts"][:, mt_mask].sum(axis=1)).ravel().astype(np.float64)
    adata.obs["pct_counts_mt"] = (_mt / np.maximum(_total, 1) * 100).astype(np.float32)
    del _mt
else:
    adata.obs["pct_counts_mt"] = np.float32(0.0)
    print("[WARNING] No MT- genes found; pct_counts_mt=0")
del _total
print(f"pct_counts_mt: n_genes={n_mt}  mean={adata.obs['pct_counts_mt'].mean():.2f}%")

# ------------------------------------------------------------------ Stress score
# [v1.9-1] Compute on full-gene counts-derived log1p CPM (temporary object).
# Score = mean log1p(CPM) across available stress genes. No [0,1] scaling.
# Reference: QRM Section 15, Key Learning — stress_score definition.
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
    "XPOT","YIF1A","YWHAZ","ZBTB17"
]

stress_genes_in_data = [g for g in STRESS_SIGNATURE_GENES if g in adata.var_names]
n_stress = len(stress_genes_in_data)
print(f"\nStress genes: {n_stress}/{len(STRESS_SIGNATURE_GENES)} found")

if n_stress >= 10:
    # Temporary AnnData for CPM log1p computation (no .copy() on X — new obj)
    adata_tmp = sc.AnnData(
        X   = adata.layers["counts"].copy(),
        var = adata.var.copy()
    )
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    stress_idx = [i for i, g in enumerate(adata.var_names) if g in stress_genes_in_data]
    stress_expr = adata_tmp.X[:, stress_idx]
    if sparse.issparse(stress_expr):
        stress_score_raw = np.array(stress_expr.mean(axis=1)).ravel().astype(np.float32)
    else:
        stress_score_raw = stress_expr.mean(axis=1).astype(np.float32)
    adata.obs["stress_score"] = stress_score_raw
    del adata_tmp, stress_expr; gc.collect()
    print(f"stress_score: mean={adata.obs['stress_score'].mean():.4f}  max={adata.obs['stress_score'].max():.4f}")
else:
    adata.obs["stress_score"] = np.float32(0.0)
    print("[WARNING] Too few stress genes; stress_score=0")

# ------------------------------------------------------------------ Cell cycle scores
# [v1.9-1] S_score and G2M_score via sc.tl.score_genes_cell_cycle on CPM log1p.
# Computed on temporary object; scores transferred back to adata.obs.
S_GENES = [
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
    "CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP",
    "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
    "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
    "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"
]
G2M_GENES = [
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
    "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B",
    "HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2",
    "DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN",
    "LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"
]

s_genes_present   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes_present = [g for g in G2M_GENES if g in adata.var_names]
print(f"\nCell cycle genes: S={len(s_genes_present)}/{len(S_GENES)}  G2M={len(g2m_genes_present)}/{len(G2M_GENES)}")

if len(s_genes_present) >= 10 and len(g2m_genes_present) >= 10:
    adata_tmp = sc.AnnData(
        X   = adata.layers["counts"].copy(),
        var = adata.var.copy()
    )
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.tl.score_genes_cell_cycle(adata_tmp, s_genes=s_genes_present, g2m_genes=g2m_genes_present)
    adata.obs["S_score"]   = adata_tmp.obs["S_score"].values.astype(np.float32)
    adata.obs["G2M_score"] = adata_tmp.obs["G2M_score"].values.astype(np.float32)
    adata.obs["phase"]     = adata_tmp.obs["phase"].values
    del adata_tmp; gc.collect()
    print(f"S_score:   mean={adata.obs['S_score'].mean():.4f}")
    print(f"G2M_score: mean={adata.obs['G2M_score'].mean():.4f}")
else:
    adata.obs["S_score"]   = np.float32(0.0)
    adata.obs["G2M_score"] = np.float32(0.0)
    adata.obs["phase"]     = "G1"
    print("[WARNING] Too few cell cycle genes; S_score=G2M_score=0")

# Verify all registered covariates exist
for cov in CONTINUOUS_COVARIATES:
    assert cov in adata.obs.columns, f"[X] Missing covariate: {cov}"
print(f"\n[OK] All covariates verified: {CONTINUOUS_COVARIATES}")

# %% [markdown]
# ## Cell 10 — Filter Genes, Cells + Drop Small Batches

# %%
print("\n" + "="*70)
print("[STEP 10] Filter genes/cells + drop small batches")
print("="*70)
n_c, n_g = adata.n_obs, adata.n_vars
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=200)
print(f"Genes : {n_g:,} -> {adata.n_vars:,}")
print(f"Cells : {n_c:,} -> {adata.n_obs:,}")
batch_counts  = adata.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index.tolist()
if small_batches:
    n_drop = int(adata.obs[BATCH_KEY].isin(small_batches).sum())
    print(f"Dropping {len(small_batches)} small batches ({n_drop:,} cells)")
    adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
    adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].cat.remove_unused_categories()
else:
    print(f"[OK] All {adata.obs[BATCH_KEY].nunique()} batches pass minimum size")

# %% [markdown]
# ## Cell 11 — HVG Selection (4-tier fallback)

# %%
print("\n" + "="*70)
print("[STEP 11] HVG selection")
print("="*70)
_drop = [c for c in _HVG_AUX_VAR_COLS if c in adata.var.columns]
if _drop: adata.var.drop(columns=_drop, inplace=True)
hvg_method = None
try:
    sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=N_HVG,
                                batch_key=BATCH_KEY, flavor="seurat_v3", subset=False)
    hvg_method = "seurat_v3 batch-aware (counts)"
except Exception as e1:
    print(f"  Tier1 failed: {e1}")
    try:
        sc.pp.highly_variable_genes(adata, layer="counts", n_top_genes=N_HVG,
                                    flavor="seurat_v3", subset=False)
        hvg_method = "seurat_v3 non-batch (counts)"
    except Exception as e2:
        print(f"  Tier2 failed: {e2}")
        _cnt = adata.layers["counts"]
        _tot = np.array(_cnt.sum(axis=1)).ravel().astype(np.float64)
        _sc  = sparse.diags(1e4 / np.maximum(_tot, 1))
        _l1p = (_sc @ _cnt).tocsr(); _l1p.data = np.log1p(_l1p.data).astype(np.float32)
        adata.layers["_log1p_tmp"] = _l1p
        del _cnt, _tot, _sc, _l1p; gc.collect()
        try:
            sc.pp.highly_variable_genes(adata, layer="_log1p_tmp", n_top_genes=N_HVG,
                                        batch_key=BATCH_KEY, flavor="seurat", subset=False)
            hvg_method = "seurat batch-aware (log1p)"
        except Exception as e3:
            print(f"  Tier3 failed: {e3}")
            sc.pp.highly_variable_genes(adata, layer="_log1p_tmp", n_top_genes=N_HVG,
                                        flavor="seurat", subset=False)
            hvg_method = "seurat non-batch (log1p)"
        del adata.layers["_log1p_tmp"]; gc.collect()
n_hvg = int(adata.var["highly_variable"].sum())
assert n_hvg > 0, "[X] Zero HVGs"
adata.uns["hvg_method"] = hvg_method
print(f"[OK] HVG: {n_hvg:,} | {hvg_method}")

# [v1.10] Force lineage anchor markers into HVG (swap lowest-variance non-forced)
if FORCE_MARKERS_IN_HVG:
    print("\n[11b] Forcing biological markers into HVG...")
    adata.var["highly_variable_pre_force"] = adata.var["highly_variable"].copy()
    _force_report = force_include_markers_in_hvg(
        adata,
        n_top_genes      = N_HVG,
        forced_markers   = FORCED_MARKERS,
        case_insensitive = FORCE_MARKERS_CASE_INSENSITIVE,
        hvg_col          = "highly_variable",
    )
    print(f"  requested       : {_force_report['n_requested']}")
    print(f"  found in data   : {_force_report['n_found_in_data']}")
    print(f"  already in HVG  : {_force_report['n_already_in_hvg']}")
    print(f"  newly added     : {_force_report['n_newly_added']}")
    print(f"  final HVG count : {_force_report['n_hvg_final']:,}")
    n_hvg = _force_report["n_hvg_final"]
    adata.uns["force_markers_report"] = _force_report

pd.Series(adata.var_names[adata.var["highly_variable"]].tolist()).to_csv(
    OUTPUT_DIR/"hvg_genes_final.csv", index=False, header=False)
pd.Series(adata.var_names.tolist()).to_csv(
    OUTPUT_DIR/"all_genes_union.csv", index=False, header=False)
print(f"[OK] gene lists saved (HVG: {n_hvg:,}  union: {adata.n_vars:,})")

# %% [markdown]
# ## Cell 12 — Save `.raw` (union gene space)

# %%
print("\n" + "="*70)
print("[STEP 12] Save .raw (union gene space, shared memory)")
print("="*70)
# bool columns in var -> int8 to avoid h5py TypeError
raw_var = adata.var.copy()
bool_cols = raw_var.select_dtypes(include=["bool"]).columns.tolist()
if bool_cols:
    raw_var[bool_cols] = raw_var[bool_cols].astype(np.int8)
# [QRM 2] shared pointer — no .copy() on X
adata.raw = sc.AnnData(
    X   = adata.layers["counts"],
    obs = adata.obs.copy(),
    var = raw_var,
)
print(f"[OK] .raw: {adata.raw.n_vars:,} genes (union outer-join space)")
import psutil
print(f"Memory: {psutil.Process(os.getpid()).memory_info().rss/1e9:.1f} GB")

# %% [markdown]
# ## Cell 13 — Subset to HVG

# %%
print("\n" + "="*70)
print("[STEP 13] HVG subset")
print("="*70)
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"[OK] {adata.n_obs:,} x {adata.n_vars:,}")
print(f"Memory: {psutil.Process(os.getpid()).memory_info().rss/1e9:.1f} GB")

# %% [markdown]
# ## Cell 14 — Normalize `.X` log1p (after HVG subset)

# %%
print("\n" + "="*70)
print("[STEP 14] Normalize .X -> log1p (HVG only)")
print("="*70)
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()
print("[OK] .X = log1p (HVG)")

# %% [markdown]
# ## Cell 15 — scVI Training
# 
# **[v1.9] Key parameter changes vs v1.8:**
# - `layer='counts'` — explicit (`.X` is now log1p; counts layer is raw)  
# - `continuous_covariate_keys`: now 4 covariates (`pct_counts_mt`, `stress_score`, `S_score`, `G2M_score`)  
# - `dropout_rate=0.2` (was 0.1)  
# - `use_layer_norm='both'`, `use_batch_norm='none'` (was default — conflicting)  
# - `n_latent=150` (was 175)  
# - `early_stopping_patience=45` (explicit)

# %%
print("\n" + "="*70)
print("[STEP 15] scVI training")
print("="*70)

# [v1.9-6] layer='counts' because .X is now log1p after Cell 14
# [v1.9-5] full covariate set registered
scvi_setup_kwargs = dict(
    layer     = "counts",
    batch_key = BATCH_KEY,
)
if CONTINUOUS_COVARIATES:
    scvi_setup_kwargs["continuous_covariate_keys"] = CONTINUOUS_COVARIATES
if CATEGORICAL_COVARIATES:
    scvi_setup_kwargs["categorical_covariate_keys"] = CATEGORICAL_COVARIATES

scvi.model.SCVI.setup_anndata(adata, **scvi_setup_kwargs)
print(f"scVI setup: batch='{BATCH_KEY}'")
print(f"  continuous_covariate_keys  = {CONTINUOUS_COVARIATES}")
print(f"  categorical_covariate_keys = {CATEGORICAL_COVARIATES}")

# [v1.9-2] Updated model architecture — matching reference pipeline v2.5.2
model_scvi = scvi.model.SCVI(
    adata,
    n_latent        = SCVI_N_LATENT,
    n_layers        = SCVI_N_LAYERS,
    n_hidden        = SCVI_N_HIDDEN,
    dropout_rate    = SCVI_DROPOUT_RATE,
    gene_likelihood = "nb",
    encode_covariates = SCVI_ENCODE_COVARIATES,
    use_layer_norm  = SCVI_USE_LAYER_NORM,
    use_batch_norm  = SCVI_USE_BATCH_NORM,
)
try:
    n_params = sum(p.numel() for p in model_scvi.module.parameters() if p.requires_grad)
except Exception:
    n_params = -1
print(f"scVI model: n_latent={SCVI_N_LATENT}  n_layers={SCVI_N_LAYERS}  "
      f"dropout={SCVI_DROPOUT_RATE}  params={n_params:,}")
print(f"  use_layer_norm={SCVI_USE_LAYER_NORM}  use_batch_norm={SCVI_USE_BATCH_NORM}")
print(f"  encode_covariates={SCVI_ENCODE_COVARIATES}")

model_scvi.train(
    max_epochs              = SCVI_MAX_EPOCHS,
    batch_size              = SCVI_BATCH_SIZE,
    train_size              = 0.9,
    early_stopping          = SCVI_EARLY_STOPPING,
    early_stopping_patience = SCVI_EARLY_STOPPING_PATIENCE,
    plan_kwargs             = {"lr": SCVI_LEARNING_RATE},
    accelerator             = _accelerator,
    devices                 = _devices,
)
adata.obsm["X_scvi"] = model_scvi.get_latent_representation()

# Record covariate info for manifest
adata.uns["scvi_registered_covariates"] = {
    "continuous"        : CONTINUOUS_COVARIATES,
    "categorical"       : CATEGORICAL_COVARIATES,
    "encode_covariates" : SCVI_ENCODE_COVARIATES,
}

SCVI_MODEL_PATH = MODEL_DIR / "scvi_model"
model_scvi.save(str(SCVI_MODEL_PATH), overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(
    SCVI_MODEL_PATH/"var_names.csv", index=False, header=False)

# scArches dry-run (2k cells)
try:
    _idx = np.random.choice(adata.n_obs, min(2000, adata.n_obs), replace=False)
    _ck  = adata[_idx].copy()
    scvi.model.SCVI.prepare_query_anndata(_ck, str(SCVI_MODEL_PATH))
    del _ck; gc.collect()
    print("[OK] scVI scArches dry-run (2k)")
except Exception as e:
    print(f"[WARNING] scVI dry-run: {e}")

print(f"[OK] scVI saved: {SCVI_MODEL_PATH}")

# %% [markdown]
# ## Cell 16 — UMAP (scVI latent)

# %%
print("\n" + "="*70)
print("[STEP 16] Neighbors + UMAP (scVI)")
print("="*70)
sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)
adata.obsm["X_umap_scvi"] = adata.obsm["X_umap"].copy()
print("[OK] X_umap_scvi saved")

# %% [markdown]
# ## Cell 17 — scANVI Training + Assign `cell_type_final_l1`
# 
# **[v1.9/v1.10] Key parameter notes:**
# - `lr=5e-4` (was 1e-3 in v1.8) — fine-tuning rate lower than scVI pretraining
# - `n_samples_per_label` removed (v1.10-4) — 5-class L1 task does not benefit
# - `weight_decay=0.0` — required for scArches compatibility
# - `early_stopping_patience=30` (explicit)

# %%
print("\n" + "="*70)
print("[STEP 17] scANVI training (L1 fully supervised)")
print("="*70)

assert LABELS_KEY in adata.obs.columns
if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
    adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi,
    unlabeled_category = UNLABELED_CATEGORY,
    labels_key         = LABELS_KEY,
)
print(f"scANVI from scVI: n_latent={SCVI_N_LATENT}  labels='{LABELS_KEY}'")

# [v1.9] lr=5e-4 (fine-tuning); [v1.10-4] n_samples_per_label removed
model_scanvi.train(
    max_epochs              = SCANVI_MAX_EPOCHS,
    batch_size              = SCANVI_BATCH_SIZE,
    train_size              = 0.9,
    early_stopping          = True,
    early_stopping_patience = SCANVI_EARLY_STOPPING_PATIENCE,
    plan_kwargs             = {"lr": SCANVI_LEARNING_RATE, "weight_decay": 0.0},
    accelerator             = _accelerator,
    devices                 = _devices,
)
print(f"  lr={SCANVI_LEARNING_RATE}  weight_decay=0.0")

adata.obsm["X_scanvi"] = model_scanvi.get_latent_representation()

soft_df = model_scanvi.predict(soft=True)
soft_df = soft_df.reindex(adata.obs_names)

adata.obs["cell_type_scanvi_pred_pan"] = soft_df.idxmax(axis=1)
adata.obs["scanvi_confidence"]         = soft_df.max(axis=1).astype(np.float32)
adata.obsm["scanvi_probabilities"]     = soft_df.values.astype(np.float32)
adata.uns["scanvi_celltype_order"]     = soft_df.columns.tolist()

# [v1.10-3] cell_type_final removed — cell_type_final_l1 is the sole output.
# This name is intentionally specific: it reflects L1 lineage-level labels only.
# It is NOT a fine-grained cell type annotation.
adata.obs["cell_type_final_l1"] = adata.obs["cell_type_scanvi_pred_pan"].astype(object)

print(f"[OK] scANVI: {soft_df.shape[1]} L1 classes | "
      f"mean confidence={adata.obs['scanvi_confidence'].mean():.4f}")
print(f"  Primary output: adata.obs['cell_type_final_l1']")

SCANVI_MODEL_PATH = MODEL_DIR / "scanvi_model"
model_scanvi.save(str(SCANVI_MODEL_PATH), overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(
    SCANVI_MODEL_PATH/"var_names.csv", index=False, header=False)

try:
    _idx = np.random.choice(adata.n_obs, min(2000, adata.n_obs), replace=False)
    _ck  = adata[_idx].copy()
    scvi.model.SCANVI.prepare_query_anndata(_ck, str(SCANVI_MODEL_PATH))
    del _ck; gc.collect()
    print("[OK] scANVI scArches dry-run (2k)")
except Exception as e:
    print(f"[WARNING] scANVI dry-run: {e}")

del model_scvi, model_scanvi; gc.collect()
if GPU_AVAILABLE: torch.cuda.empty_cache()
print(f"[OK] scANVI saved: {SCANVI_MODEL_PATH}")

# %% [markdown]
# ## Cell 18 — UMAP (scANVI) + Agreement Check

# %%
print("\n" + "="*70)
print("[STEP 18] Neighbors + UMAP (scANVI) + agreement check")
print("="*70)
scanvi_umap_bundle = fit_bundle(
    adata,
    "X_scanvi",
    OUTPUT_DIR,
    primary_umap_key="X_umap_scanvi",
    operator_filename="umap_scanvi_operator.joblib",
    manifest_filename="scanvi_umap_bundle.json",
    umap_params={
        "n_neighbors": 20,
        "n_components": 2,
        "min_dist": 0.3,
        "spread": 1.0,
        "metric": "euclidean",
        "random_state": 42,
    },
    set_default_x_umap=True,
    extra_manifest={"pipeline": _NB_STEM, "label_tier": LABEL_TIER},
)
print(f"[OK] X_umap_scanvi saved | operator={scanvi_umap_bundle['operator_path']}")

# [v1.10-2 / P0-3] IMPORTANT: this agreement metric is train-set re-substitution.
# All cells were used for training (ALLOW_UNKNOWN=False, no held-out split).
# A high agreement rate here only confirms the model memorised training labels;
# it does NOT constitute independent validation of generalisation.
# Interpretation: use as a sanity check (should be >95% for 5 clean lineage classes).
# Do NOT report this value as classification accuracy in manuscripts.
agree = (adata.obs["cell_type_final_l1"].astype(str) ==
         adata.obs["cell_type_input_l1"].astype(str))
agreement_rate = float(agree.mean())
adata.uns["agreement_rate_overall"]      = agreement_rate
adata.uns["agreement_rate_is_resubstitution"] = True  # explicit flag
print(f"L1 agreement (train-set re-substitution, NOT independent validation): {agreement_rate*100:.2f}%")
for lin in adata.obs["lineage_source"].unique():
    ag = agree[adata.obs["lineage_source"]==lin].mean()
    print(f"  {lin:12s}: {ag*100:.2f}%")

# %% [markdown]
# ## Cell 19 — UMAP Figures

# %%
print("\n" + "="*70)
print("[STEP 19] UMAP figures")
print("="*70)
sc.settings.vector_friendly = True
_cols = ["lineage_source", "cell_type_final_l1", "cell_type_L2", BATCH_KEY, "scanvi_confidence"]
for _basis, _lbl in [("X_umap_scvi","scVI"), ("X_umap_scanvi","scANVI")]:
    fig, axes = plt.subplots(1, len(_cols), figsize=(6*len(_cols), 5))
    for ax, col in zip(axes, _cols):
        sc.pl.embedding(adata, basis=_basis, color=col, ax=ax, show=False, title=col,
                        legend_loc="right margin" if adata.obs[col].nunique()<=30 else "none",
                        frameon=False)
    fig.suptitle(f"{_lbl} UMAP — All Lineages (L1 scANVI refined)", y=1.02, fontsize=14)
    fig.tight_layout()
    fig.savefig(FIG_DIR/f"allcells_{_lbl.lower()}_umap_overview.pdf",
                dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"[OK] {_lbl} UMAP saved")

# %% [markdown]
# ## Cell 20 — Pre-write Cleanup

# %%
print("\n" + "="*70)
print("[STEP 20] Pre-write cleanup")
print("="*70)
for _attr in ("obs","var"):
    _df = getattr(adata, _attr)
    if "_index" in _df.columns:
        setattr(adata, _attr, _df.rename(columns={"_index":"orig_index"}))

OBSM_KEEP = {"X_scvi","X_scanvi","X_umap","X_umap_scvi","X_umap_scanvi","scanvi_probabilities"}
stale = [k for k in list(adata.obsm.keys()) if k not in OBSM_KEEP]
for k in stale: del adata.obsm[k]
if stale: print(f"  Removed stale obsm: {stale}")

for col in adata.obs.select_dtypes(include=["category"]).columns:
    cats = adata.obs[col].cat.categories
    if hasattr(cats.dtype,"name") and cats.dtype.name in ("string","StringDtype"):
        adata.obs[col] = adata.obs[col].cat.rename_categories(cats.astype(object))

def _sanitize_object_cols(df, df_name):
    bool_fixed, str_fixed = [], []
    for col in df.columns:
        s = df[col]
        if s.dtype != object: continue
        non_na = s.dropna()
        if len(non_na) == 0:
            df[col] = s.fillna("").astype(str); str_fixed.append(col); continue
        if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
            df[col] = s.fillna(False).astype(np.int8); bool_fixed.append(col); continue
        if non_na.map(lambda x: not isinstance(x, str)).any():
            df[col] = s.map(lambda x: "" if pd.isna(x) else str(x)); str_fixed.append(col)
    if bool_fixed: print(f"  {df_name} bool->int8: {bool_fixed}")
    if str_fixed:  print(f"  {df_name} mixed->str: {str_fixed}")

_sanitize_object_cols(adata.obs, "obs")
_sanitize_object_cols(adata.var, "var")
if adata.raw is not None:
    _sanitize_object_cols(adata.raw.var, "raw.var")
print("[OK] Pre-write cleanup done")

# %% [markdown]
# ## Cell 21 — Reference Manifest JSON

# %%
print("\n" + "="*70)
print("[STEP 21] Reference manifest JSON")
print("="*70)

_cov = adata.uns.get("scvi_registered_covariates", {
    "continuous"        : CONTINUOUS_COVARIATES,
    "categorical"       : CATEGORICAL_COVARIATES,
    "encode_covariates" : SCVI_ENCODE_COVARIATES,
})

manifest = {
    "reference_h5ad"              : str(OUTPUT_H5AD),
    "scvi_model_dir"              : str(SCVI_MODEL_PATH),
    "scanvi_model_dir"            : str(SCANVI_MODEL_PATH),
    "hvg_var_names_csv"           : str(SCVI_MODEL_PATH/"var_names.csv"),
    "all_genes_union_csv"         : str(OUTPUT_DIR/"all_genes_union.csv"),
    "batch_key"                   : BATCH_KEY,
    "batch_granularity"           : "sample",
    "label_key"                   : LABELS_KEY,
    "label_tier"                  : LABEL_TIER,
    "l2_key"                      : "cell_type_L2",
    "final_label_key"             : "cell_type_final_l1",
    "input_label_key"             : "cell_type_input_l1",
    "input_label_l2_key"          : "cell_type_input_l2",
    "input_label_l3_key"          : "cell_type_input_l3",
    "prediction_key"              : "cell_type_scanvi_pred_pan",
    "confidence_key"              : "scanvi_confidence",
    "unlabeled_category"          : UNLABELED_CATEGORY,
    "dropped_unknown_cells"       : int(adata.uns.get("dropped_unknown_cells", 0)),
    "gene_space"                  : "outer_join_union",
    "n_latent"                    : SCVI_N_LATENT,
    "scvi_dropout_rate"           : SCVI_DROPOUT_RATE,
    "scvi_use_layer_norm"         : SCVI_USE_LAYER_NORM,
    "scvi_use_batch_norm"         : SCVI_USE_BATCH_NORM,
    "scvi_early_stopping_patience": SCVI_EARLY_STOPPING_PATIENCE,
    "scanvi_lr"                   : SCANVI_LEARNING_RATE,
    "scanvi_n_samples_per_label"  : None,   # removed in v1.10-4
    "scanvi_early_stopping_patience": SCANVI_EARLY_STOPPING_PATIENCE,
    "agreement_rate_is_resubstitution": True,   # [v1.10-2] not independent validation
    "continuous_covariate_keys"   : _cov["continuous"],
    "categorical_covariate_keys"  : _cov["categorical"],
    "encode_covariates"           : _cov["encode_covariates"],
    "scarches_compatible"         : True,
    "force_markers_in_hvg"        : FORCE_MARKERS_IN_HVG,
    "force_markers_n_requested"   : len(FORCED_MARKERS),
    "force_markers_report"        : adata.uns.get("force_markers_report", {}),
    "overlap_whitelist"           : sorted(OVERLAP_WHITELIST),
    "agreement_rate"              : float(adata.uns.get("agreement_rate_overall", -1)),
    "version"                     : VERSION,
    "date"                        : datetime.now().strftime("%Y-%m-%d"),
    "lineage_inputs"              : {
        k: {"h5ad": str(v["h5ad"]), "label_built": "cell_type_L3",
            "batch_local": v["batch_key_local"], "l1_label": v["name"]}
        for k, v in LINEAGE_INPUTS.items()
    },
}
with open(OUTPUT_DIR/"reference_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)
print("[OK] reference_manifest.json written")
print(f"  continuous_covariate_keys : {manifest['continuous_covariate_keys']}")
print(f"  categorical_covariate_keys: {manifest['categorical_covariate_keys']}")
print(f"  encode_covariates         : {manifest['encode_covariates']}")
print(f"  scvi_use_layer_norm       : {manifest['scvi_use_layer_norm']}")
print(f"  scvi_use_batch_norm       : {manifest['scvi_use_batch_norm']}")
print(f"  scanvi_lr                 : {manifest['scanvi_lr']}")
print(f"  scanvi_n_samples_per_label: {manifest['scanvi_n_samples_per_label']}")

# %% [markdown]
# ## Cell 22 — Save h5ad

# %%
print("\n" + "="*70)
print(f"[STEP 22] Saving {OUTPUT_H5AD.name}")
print("="*70)
adata.write_h5ad(OUTPUT_H5AD, compression="gzip", compression_opts=9)
print(f"[OK] Saved: {OUTPUT_H5AD}")
print(f"     Size : {OUTPUT_H5AD.stat().st_size/1e9:.2f} GB")

# %% [markdown]
# ## Cell 23 — AnnData Structure JSON Output

# %%
print("\n" + "="*70)
print("[STEP 23] AnnData structure JSON")
print("="*70)
n_lab = int((adata.obs[LABELS_KEY].astype(str) != UNLABELED_CATEGORY).sum())
n_unk = int((adata.obs[LABELS_KEY].astype(str) == UNLABELED_CATEGORY).sum())
def _li(X):
    sp = sparse.issparse(X)
    return {"type": type(X).__name__ if sp else "ndarray",
            "shape": list(X.shape), "dtype": str(X.dtype), "sparse": sp}
def _ci(s):
    d = {"dtype": str(s.dtype), "non_null": int(s.notna().sum())}
    try: d["unique"] = int(s.nunique())
    except: pass
    return d
_cov = adata.uns.get("scvi_registered_covariates", {})
structure = {
    "pipeline"               : PIPELINE_NAME,
    "version"                : VERSION,
    "timestamp"              : datetime.now().isoformat(),
    "output_h5ad"            : str(OUTPUT_H5AD),
    "shape"                  : [adata.n_obs, adata.n_vars],
    "X"                      : _li(adata.X),
    "raw_n_vars"             : int(adata.raw.n_vars) if adata.raw else None,
    "raw_gene_space"         : "post_concat_union (outer join)",
    "layers"                 : {k: _li(v) for k, v in adata.layers.items()},
    "obsm"                   : {k: {"shape": list(v.shape), "dtype": str(v.dtype)}
                                 for k, v in adata.obsm.items()},
    "batch_key"              : BATCH_KEY,
    "n_batches"              : int(adata.obs[BATCH_KEY].nunique()),
    "label_tier"             : LABEL_TIER,
    "labels_key"             : LABELS_KEY,
    "n_L1_classes"           : int(adata.obs[LABELS_KEY].astype(str).nunique()),
    "n_L2_classes"           : int(adata.obs["cell_type_L2"].nunique()),
    "n_labeled"              : n_lab,
    "n_unknown"              : n_unk,
    "dropped_unknown"        : int(adata.uns.get("dropped_unknown_cells", 0)),
    "hvg_method"             : adata.uns.get("hvg_method", "unknown"),
    "n_latent"               : SCVI_N_LATENT,
    "scvi_dropout_rate"      : SCVI_DROPOUT_RATE,
    "scvi_use_layer_norm"    : SCVI_USE_LAYER_NORM,
    "scvi_use_batch_norm"    : SCVI_USE_BATCH_NORM,
    "scanvi_lr"              : SCANVI_LEARNING_RATE,
    "scanvi_n_samples_per_label": None,   # removed in v1.10-4
    "continuous_covariates"  : _cov.get("continuous", []),
    "categorical_covariates" : _cov.get("categorical", []),
    "encode_covariates"      : _cov.get("encode_covariates", True),
    "agreement_rate"         : float(adata.uns.get("agreement_rate_overall", -1)),
    "scanvi_celltype_order"  : adata.uns.get("scanvi_celltype_order", []),
    "key_obs_columns"        : {
        "cell_type_final_l1"        : _ci(adata.obs["cell_type_final_l1"]),
        "cell_type_input_l1"        : _ci(adata.obs["cell_type_input_l1"]),
        "cell_type_input_l2"        : _ci(adata.obs["cell_type_input_l2"]),
        "cell_type_input_l3"        : _ci(adata.obs["cell_type_input_l3"]),
        "cell_type_L2"              : _ci(adata.obs["cell_type_L2"]),
        LABELS_KEY                  : _ci(adata.obs[LABELS_KEY]),
        "cell_type_scanvi_pred_pan" : _ci(adata.obs["cell_type_scanvi_pred_pan"]),
        "scanvi_confidence"         : _ci(adata.obs["scanvi_confidence"]),
        "lineage_source"            : _ci(adata.obs["lineage_source"]),
        BATCH_KEY                   : _ci(adata.obs[BATCH_KEY]),
        "pct_counts_mt"             : _ci(adata.obs["pct_counts_mt"]),
        "stress_score"              : _ci(adata.obs["stress_score"]),
        "S_score"                   : _ci(adata.obs["S_score"]),
        "G2M_score"                 : _ci(adata.obs["G2M_score"]),
    },
    "lineage_distribution"   : adata.obs["lineage_source"].value_counts().to_dict(),
    "L1_final_distribution"  : adata.obs["cell_type_final_l1"].value_counts().to_dict(),
    "confidence_stats"       : {k: float(v)
                                 for k,v in adata.obs["scanvi_confidence"].describe().items()},
    "log_file"               : str(LOG_FILE),
}
struct_path = OUTPUT_DIR / f"anndata_structure_{DATE_TAG}_v{VERSION}.json"
with open(struct_path, "w") as f: json.dump(structure, f, indent=2)
print(f"[OK] Structure JSON: {struct_path}")
print(f"  shape          : {adata.n_obs:,} x {adata.n_vars:,}")
print(f"  .raw           : {adata.raw.n_vars:,} genes")
print(f"  layers         : {list(adata.layers.keys())}")
print(f"  obsm           : {list(adata.obsm.keys())}")
print(f"  covariates     : {structure['continuous_covariates']}")
print(f"  agreement      : {float(adata.uns.get('agreement_rate_overall',-1))*100:.2f}%")

# %% [markdown]
# ## Cell 24 — Pipeline Run Log + Close Log File

# %%
elapsed_min = (time.time() - PIPELINE_START) / 60
n_lab = int((adata.obs[LABELS_KEY].astype(str) != UNLABELED_CATEGORY).sum())
n_unk = int((adata.obs[LABELS_KEY].astype(str) == UNLABELED_CATEGORY).sum())
_cov  = adata.uns.get("scvi_registered_covariates", {})

_log = [
    "="*80, f"{PIPELINE_NAME} Run Log  (v{VERSION})", "="*80,
    f"written_at                    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    f"log_file                      : {LOG_FILE}",
    f"output_h5ad                   : {OUTPUT_H5AD}",
    f"n_cells                       : {adata.n_obs:,}",
    f"n_genes_hvg                   : {adata.n_vars:,}",
    f"n_genes_raw (.raw union space): {adata.raw.n_vars if adata.raw else 'None'}",
    f"n_batches                     : {adata.obs[BATCH_KEY].nunique():,}",
    f"hvg_method                    : {adata.uns.get('hvg_method','unknown')}",
    "",
    "[v1.9] scVI hyperparameters:",
    f"  n_latent                    : {SCVI_N_LATENT}  (was 175 in v1.8)",
    f"  n_layers                    : {SCVI_N_LAYERS}",
    f"  n_hidden                    : {SCVI_N_HIDDEN}",
    f"  dropout_rate                : {SCVI_DROPOUT_RATE}  (was 0.1 in v1.8)",
    f"  use_layer_norm              : {SCVI_USE_LAYER_NORM}  (was default=none in v1.8)",
    f"  use_batch_norm              : {SCVI_USE_BATCH_NORM}  (was default=both in v1.8)",
    f"  encode_covariates           : {SCVI_ENCODE_COVARIATES}",
    f"  early_stopping_patience     : {SCVI_EARLY_STOPPING_PATIENCE}",
    f"  lr                          : {SCVI_LEARNING_RATE}",
    "",
    "[v1.9] scANVI hyperparameters:",
    f"  lr                          : {SCANVI_LEARNING_RATE}  (was 1e-3 in v1.8)",
    f"  n_samples_per_label         : removed in v1.10-4 (5-class L1, not needed)",
    f"  early_stopping_patience     : {SCANVI_EARLY_STOPPING_PATIENCE}",
    f"  weight_decay                : 0.0",
    "",
    "[v1.9] Covariates:",
    f"  continuous  : {_cov.get('continuous', [])}",
    f"  categorical : {_cov.get('categorical', [])}",
    "",
    f"label_tier                    : {LABEL_TIER}",
    f"n_labeled                     : {n_lab:,}",
    f"n_unknown                     : {n_unk:,}",
    f"dropped_unknown               : {int(adata.uns.get('dropped_unknown_cells',0)):,}",
    f"agreement_rate                : {_fmt_float(adata.uns.get('agreement_rate_overall'))}  [train-set re-substitution, NOT independent validation]",
    f"accelerator                   : {_accelerator}",
    f"elapsed_minutes               : {elapsed_min:.2f}",
    "",
    "lineage_source distribution:",
    adata.obs["lineage_source"].value_counts().to_string(),
    "",
    "cell_type_final (L1 scANVI) distribution:",
    adata.obs["cell_type_final"].value_counts().to_string(),
    "",
    "scanvi_confidence percentiles:",
    adata.obs["scanvi_confidence"].describe().to_string(),
]

(OUTPUT_DIR/"pipeline_run_log.txt").write_text("\n".join(_log)+"\n", encoding="utf-8")

print("="*70)
print("[DONE] Pipeline complete")
print(f"Elapsed  : {elapsed_min:.2f} min")
print(f"Output   : {OUTPUT_H5AD}")
print(f"Log      : {LOG_FILE}")
print("="*70)

try:
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
except Exception: pass
try:
    if "_log_fh" in globals() and not getattr(_log_fh,"closed",True):
        _log_fh.flush(); _log_fh.close()
        print(f"[OK] Log closed: {LOG_FILE}")
except Exception as e:
    print(f"[WARNING] Log close: {e}")


