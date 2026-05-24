# %% [markdown]
# # Myeloid Merged Pipeline v2.1 — CellTypist-Filtered Retrain
#
# **Purpose:**
# Downstream of v2.0. Loads the v2.0 preprocessed checkpoint (which already has
# CellTypist labels), removes non-myeloid contaminants based on CellTypist lineage
# classification, then retrains scVI + scANVI on the clean myeloid-only subset.
#
# **Input:**  myeloid_merged_v2_adata_preprocessed.h5ad  (from v2.0 Cell 13)
#             — already contains CellTypist labels, HVG flags, covariates, no .raw yet
#
# **Filter logic:**
#   - Keep cells whose CellTypist label maps to a myeloid lineage
#   - Contaminant labels (T, B, Epithelial, Stromal, NK …) are removed
#   - After filtering: drop batches with < MIN_BATCH_CELLS cells
#
# **Architecture (same as v2.0):**
#   - scANVI labels = CellTypist filtered labels for ALL surviving cells
#   - Low-confidence cells (< threshold) → Unknown
#
# **v2.9 compliance:**
#   - rasterized=True on all UMAP scatter plots, dpi=300

# %% [markdown]
# ## Cell 0 — Imports & Global Settings

# %%
import os
import sys, warnings, json, gc, joblib
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.sparse import issparse, csr_matrix
from scipy.stats import entropy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import torch
import scanpy as sc
import scvi
import celltypist
from umap import UMAP

warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None

print("=" * 80)
print("Myeloid Merged Pipeline v2.1  (CellTypist-Filtered Retrain)")
print("=" * 80)
gpu_available = torch.cuda.is_available()
print(f"GPU available: {gpu_available}")
if gpu_available:
    print(f"GPU: {torch.cuda.get_device_name(0)}")
scvi.settings.dl_num_workers = 0
import time


# %% [markdown]
# ## Cell 1 — Configuration

# %%
# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input: preprocessed checkpoint from v2.0 (has CellTypist labels, HVG flags)
INPUT_H5AD = "/home/h2048/data/py/20260308/myeloid_only_merged_pipeline_v2/myeloid_merged_v2_adata_preprocessed.h5ad"

OUTPUT_DIR    = "/home/h2048/data/py/20260308/myeloid_only_merged_pipeline_v2_1_filtered"
OUTPUT_PREFIX = "myeloid_merged_v2_1"

BATCH_KEY  = "sample"
TISSUE_KEY = "tissue"

# CellTypist column to use for lineage filtering (from v2.0 output)
CELLTYPIST_LABEL_COL  = "celltypist_label_direct"   # unfiltered majority-vote
CELLTYPIST_CONF_KEY   = "celltypist_confidence"
CELLTYPIST_CONF_THRESHOLD = 0.5

# Minimum cells per batch after filtering (batches below this are dropped)
MIN_BATCH_CELLS = 10

# ==============================================================================
# LINEAGE FILTER MAP
# Cells are KEPT if their CellTypist label maps to a MYELOID lineage.
# All other lineages (contaminants) are removed.
# Extend this dict if your CellTypist model produces additional labels.
# ==============================================================================
CELLTYPIST_LINEAGE_MAP = {
    # ── Myeloid: KEEP ──────────────────────────────────────────────────────
    "Classical monocytes":              "Monocyte",
    "Non-classical monocytes":          "Monocyte",
    "Intermediate macrophages":         "Monocyte",
    "Monocytes":                        "Monocyte",
    "Macrophages":                      "Macrophage",
    "Alveolar macrophages":             "Macrophage",
    "Intestinal macrophages":           "Macrophage",
    "Inflammatory macrophages":         "Macrophage",
    "Resident macrophages":             "Macrophage",
    "DC1":                              "DC",
    "DC2":                              "DC",
    "pDC":                              "DC",
    "Plasmacytoid DCs":                 "DC",
    "Migratory DCs":                    "DC",
    "Mast cells":                       "Mast",
    "Basophils":                        "Basophil",
    "Eosinophils":                      "Eosinophil",
    "Neutrophils":                      "Neutrophil",
    "Low-density neutrophils":          "Neutrophil",
    "Erythroid cells":                  "Erythroid",
    "Megakaryocytes":                   "Megakaryocyte",
    # ── Contaminants: REMOVE ──────────────────────────────────────────────
    "Regulatory T cells":               "contaminant_T",
    "Tcm/Naive helper T cells":         "contaminant_T",
    "Tem/Effector helper T cells":      "contaminant_T",
    "Tcm/Naive cytotoxic T cells":      "contaminant_T",
    "Tem/Effector cytotoxic T cells":   "contaminant_T",
    "Tem/Temra cytotoxic T cells":      "contaminant_T",
    "Innate lymphoid cells":            "contaminant_ILC",
    "NK cells":                         "contaminant_NK",
    "NKT cells":                        "contaminant_NK",
    "B cells":                          "contaminant_B",
    "Memory B cells":                   "contaminant_B",
    "Naive B cells":                    "contaminant_B",
    "Plasma cells":                     "contaminant_B",
    "Plasmablasts":                     "contaminant_B",
    "Epithelial cells":                 "contaminant_epi",
    "Fibroblasts":                      "contaminant_stromal",
    "Smooth muscle cells":              "contaminant_stromal",
    "Endothelial cells":                "contaminant_stromal",
    "Pericytes":                        "contaminant_stromal",
}

# Lineage categories that should be KEPT (myeloid)
MYELOID_LINEAGES = {"Monocyte", "Macrophage", "DC", "Mast", "Basophil",
                    "Eosinophil", "Neutrophil", "Erythroid", "Megakaryocyte"}

# scVI / scANVI parameters
N_HVG               = 4000
FORCE_MARKERS_IN_HVG = True

SCVI_N_LATENT   = 100
SCVI_N_LAYERS   = 2
SCVI_N_HIDDEN   = 128
SCVI_DROPOUT    = 0.1
MAX_EPOCHS_SCVI = 400

MAX_EPOCHS_SCANVI  = 200
UNLABELED_CATEGORY = "Unknown"

BATCH_SIZE    = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY  = 0.0

# scANVI label source (should match v2.0 intent)
SCANVI_LABEL_SOURCE          = "filtered"  # use low-conf filtered column
CELLTYPIST_DIRECT_LABEL_KEY  = "celltypist_label_direct"
CELLTYPIST_DIRECT_FILT_KEY   = "celltypist_label_direct_filt"

QUERY_LEIDEN_RESOLUTION = 1.0

MYELOID_CORE_MARKERS      = ["LYZ","CD14","CD33","PTPRC","ITGAM","ITGAX"]
CLASSICAL_MONO_MARKERS    = ["CD14","FCGR1A","CCR2","CD36","SELL"]
NONCLASSICAL_MONO_MARKERS = ["FCGR3A","FCER1G","PICALM","RHOC"]
INTERMEDIATE_MONO_MARKERS = ["CD14","FCGR3A","CD86","HLA-DRA"]
MACROPHAGE_MARKERS        = ["CD68","CD163","MRC1","MARCO","MSR1","FCGR2A","APOE","C1QA"]
M1_MACROPHAGE_MARKERS     = ["CD86","CD80","TNF","IL1B","NOS2","CXCL9","CXCL10"]
M2_MACROPHAGE_MARKERS     = ["CD163","MRC1","ARG1","IL10","TGFB1","CCL22","PPARG"]
CDC1_MARKERS              = ["CLEC9A","XCR1","WDFY4","IRF8","BATF3"]
CDC2_MARKERS              = ["CD1C","FCER1A","CLEC10A","CD2","ESAM"]
PDC_MARKERS               = ["LILRA4","CLEC4C","NRP1","TCF4","IRF7","GZMB"]
NEUTROPHIL_MARKERS        = ["S100A8","S100A9","FCGR3B","CSF3R","CEACAM8","CXCR2","FCN1"]
MAST_CELL_MARKERS         = ["KIT","TPSAB1","TPSB2","HPGDS","MS4A2","FCER1A","CPA3"]
EOSINOPHIL_MARKERS        = ["EPX","PRG2","CLC","CCL26","SIGLEC8"]
PROLIF_MARKERS            = ["MKI67","TOP2A","PCNA"]

FORCED_MARKERS = list(set(
    MYELOID_CORE_MARKERS + CLASSICAL_MONO_MARKERS + NONCLASSICAL_MONO_MARKERS +
    INTERMEDIATE_MONO_MARKERS + MACROPHAGE_MARKERS + M1_MACROPHAGE_MARKERS +
    M2_MACROPHAGE_MARKERS + CDC1_MARKERS + CDC2_MARKERS + PDC_MARKERS +
    NEUTROPHIL_MARKERS + MAST_CELL_MARKERS + EOSINOPHIL_MARKERS + PROLIF_MARKERS
))

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
sc.settings.seed   = RANDOM_SEED
scvi.settings.seed = RANDOM_SEED

PIPELINE_START = time.time()
FIGURE_DPI    = 300
FIGURE_FORMAT = "pdf"


# %% [markdown]
# ## Cell 2 — Helper Functions

# %%
def run_query_only_leiden(adata_merged, resolution=QUERY_LEIDEN_RESOLUTION):
    print(f"\n[Leiden] Query-only clustering (resolution={resolution})...")
    qmask  = adata_merged.obs["data_source"] == "query"
    qcells = adata_merged.obs_names[qmask]
    if len(qcells) < 10:
        adata_merged.obs["leiden_query"] = "N/A"; return
    X_q = adata_merged.obsm["X_scVI"][qmask.values]
    ad_ = sc.AnnData(X=X_q, obs=adata_merged.obs.loc[qcells].copy())
    ad_.obsm["X_scVI"] = X_q
    sc.pp.neighbors(ad_, use_rep="X_scVI", n_neighbors=30, random_state=RANDOM_SEED)
    sc.tl.leiden(ad_, resolution=resolution, random_state=RANDOM_SEED)
    lf = pd.Series("N/A", index=adata_merged.obs_names, dtype="object")
    lf.loc[qcells] = "qry_" + ad_.obs["leiden"].astype(str)
    adata_merged.obs["leiden_query"] = lf.values
    print(f"  -> {ad_.obs['leiden'].nunique()} query-only clusters")
    del ad_; gc.collect()


# %% [markdown]
# ## Cell 3 — Initialization

# %%
output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Output directory: {output_dir}")


# %% [markdown]
# ## Cell 4 — Step 1: Load Preprocessed Checkpoint

# %%
print("\n" + "=" * 80)
print("[Step 1] Loading v2.0 preprocessed checkpoint...")
print("=" * 80)
t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time()-t0:.1f}s  shape: {adata.shape}")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"  layers:      {list(adata.layers.keys())}")

# Verify CellTypist labels are present
assert CELLTYPIST_LABEL_COL in adata.obs.columns, \
    f"CellTypist column '{CELLTYPIST_LABEL_COL}' not found. " \
    f"Available: {list(adata.obs.columns)}"
assert CELLTYPIST_CONF_KEY in adata.obs.columns, \
    f"Confidence column '{CELLTYPIST_CONF_KEY}' not found."
print(f"\n[OK] CellTypist labels present: {adata.obs[CELLTYPIST_LABEL_COL].nunique()} unique types")


# %% [markdown]
# ## Cell 5 — Step 2: Classify Cells by Lineage

# %%
print("\n" + "=" * 80)
print("[Step 2] Classifying cells by CellTypist lineage...")
print("=" * 80)

ct_labels = adata.obs[CELLTYPIST_LABEL_COL].astype(str)

# Map each label to a lineage; unmapped labels get "other_unknown"
lineage_series = ct_labels.map(CELLTYPIST_LINEAGE_MAP).fillna("other_unknown")
adata.obs["celltypist_lineage"] = lineage_series.values

print("\n[INFO] Lineage distribution:")
lineage_counts = adata.obs["celltypist_lineage"].value_counts()
total = len(adata)
for lin, n in lineage_counts.items():
    flag = "  [KEEP]" if lin in MYELOID_LINEAGES else "  [REMOVE]"
    print(f"  {lin:<30} {n:>7,}  ({n/total*100:.1f}%){flag}")

# Warn about unmapped labels (appear as "other_unknown")
unmapped_labels = ct_labels[lineage_series == "other_unknown"].unique()
if len(unmapped_labels) > 0:
    print(f"\n  [WARN] {len(unmapped_labels)} CellTypist labels not in CELLTYPIST_LINEAGE_MAP:")
    for lbl in sorted(unmapped_labels):
        n = int((ct_labels == lbl).sum())
        print(f"    '{lbl}'  n={n:,}  → treated as contaminant (other_unknown)")
    print("  → Add these labels to CELLTYPIST_LINEAGE_MAP if they are myeloid!")


# %% [markdown]
# ## Cell 6 — Step 3: Filter Non-Myeloid Cells

# %%
print("\n" + "=" * 80)
print("[Step 3] Filtering to myeloid cells only...")
print("=" * 80)

n_before = adata.n_obs
keep_mask = adata.obs["celltypist_lineage"].isin(MYELOID_LINEAGES)

# Cells whose lineage is "other_unknown" require manual decision.
# By default they are EXCLUDED (conservative). Change to True to include them.
KEEP_OTHER_UNKNOWN = False
if KEEP_OTHER_UNKNOWN:
    keep_mask |= (adata.obs["celltypist_lineage"] == "other_unknown")
    print("  [NOTE] other_unknown cells are INCLUDED (KEEP_OTHER_UNKNOWN=True)")

adata_filt = adata[keep_mask].copy()
n_after    = adata_filt.n_obs
n_removed  = n_before - n_after

print(f"  Before:  {n_before:,} cells")
print(f"  Removed: {n_removed:,} non-myeloid cells ({n_removed/n_before*100:.1f}%)")
print(f"  After:   {n_after:,} cells")

print("\n  Remaining lineage breakdown:")
for lin, n in adata_filt.obs["celltypist_lineage"].value_counts().items():
    print(f"    {lin:<30} {n:>7,}  ({n/n_after*100:.1f}%)")

print("\n  Remaining data_source breakdown:")
for src, n in adata_filt.obs["data_source"].value_counts().items():
    print(f"    {src:<20} {n:>7,}")

del adata; gc.collect()


# %% [markdown]
# ## Cell 7 — Step 4: Drop Small Batches (Post-Filter)

# %%
print("\n" + "=" * 80)
print(f"[Step 4] Dropping batches with < {MIN_BATCH_CELLS} cells post-filter...")
print("=" * 80)

batch_counts = adata_filt.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < MIN_BATCH_CELLS].index.tolist()

if small_batches:
    print(f"  Dropping {len(small_batches)} small batches:")
    for b in small_batches:
        print(f"    {b}  (n={batch_counts[b]})")
    keep_batch = ~adata_filt.obs[BATCH_KEY].isin(small_batches)
    adata_filt = adata_filt[keep_batch].copy()
    print(f"  Cells after batch drop: {adata_filt.n_obs:,}")
else:
    print(f"  All {batch_counts.shape[0]} batches have >= {MIN_BATCH_CELLS} cells — no drop needed")

# Re-encode BATCH_KEY and TISSUE_KEY as category (categories may have changed)
for key in [BATCH_KEY, TISSUE_KEY, "data_source"]:
    if key in adata_filt.obs.columns:
        adata_filt.obs[key] = adata_filt.obs[key].astype(str).astype("category")

print(f"\n  Final shape: {adata_filt.shape}")
print(f"  Batches: {adata_filt.obs[BATCH_KEY].nunique()}")
gc.collect()


# %% [markdown]
# ## Cell 8 — Step 5: Rebuild CellTypist Filtered Label Column

# %%
print("\n" + "=" * 80)
print("[Step 5] Rebuilding CellTypist filtered label column...")
print("=" * 80)

# Ensure both CellTypist direct and filtered columns exist.
# If the preprocessed checkpoint already has them, just verify; otherwise rebuild.

if CELLTYPIST_DIRECT_FILT_KEY not in adata_filt.obs.columns:
    print(f"  '{CELLTYPIST_DIRECT_FILT_KEY}' not found — rebuilding from confidence...")
    raw_labels = adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].astype(str)
    conf       = adata_filt.obs[CELLTYPIST_CONF_KEY].values
    filt       = np.where(conf >= CELLTYPIST_CONF_THRESHOLD, raw_labels.values, UNLABELED_CATEGORY)
    adata_filt.obs[CELLTYPIST_DIRECT_FILT_KEY] = pd.Categorical(filt)
    n_low = int((conf < CELLTYPIST_CONF_THRESHOLD).sum())
    print(f"  -> {n_low:,} low-conf cells set to '{UNLABELED_CATEGORY}'")
else:
    n_low = int((adata_filt.obs[CELLTYPIST_DIRECT_FILT_KEY].astype(str) == UNLABELED_CATEGORY).sum())
    print(f"  '{CELLTYPIST_DIRECT_FILT_KEY}' already present")
    print(f"  -> {n_low:,} Unknown (low-conf) cells")

# Select scANVI label source
label_source_col = (
    CELLTYPIST_DIRECT_FILT_KEY if SCANVI_LABEL_SOURCE == "filtered"
    else CELLTYPIST_DIRECT_LABEL_KEY
)
print(f"\n  scANVI label source column: {label_source_col}")


# %% [markdown]
# ## Cell 9 — Step 6: HVG Re-Selection on Filtered Data

# %%
print("\n" + "=" * 80)
print("[Step 6] Selecting HVGs on filtered myeloid data...")
print("=" * 80)

# Clear old HVG flags (from v2.0 run on all cells)
if "highly_variable" in adata_filt.var.columns:
    adata_filt.var.drop(columns=["highly_variable"], errors="ignore", inplace=True)

hvg_method = "unknown"
try:
    sc.pp.highly_variable_genes(adata_filt, layer="counts", n_top_genes=N_HVG,
                                batch_key=BATCH_KEY, flavor="seurat_v3", subset=False)
    hvg_method = "batch_seurat_v3"
except Exception as e1:
    print(f"  -> batch-aware failed ({str(e1)[:60]}), trying standard...")
    try:
        sc.pp.highly_variable_genes(adata_filt, layer="counts", n_top_genes=N_HVG,
                                    flavor="seurat_v3", subset=False)
        hvg_method = "standard_seurat_v3"
    except Exception as e2:
        print(f"  -> seurat_v3 failed ({str(e2)[:60]}), fallback to cell_ranger")
        sc.pp.highly_variable_genes(adata_filt, layer="counts", n_top_genes=N_HVG,
                                    flavor="cell_ranger", subset=False)
        hvg_method = "cell_ranger"
print(f"  -> HVG method: {hvg_method}")

# Ensure symbol_base exists for marker forcing
if "symbol_base" not in adata_filt.var.columns:
    adata_filt.var["symbol_base"] = adata_filt.var_names.str.replace(r"-\d+$", "", regex=True)

if FORCE_MARKERS_IN_HVG:
    n_added = 0
    mset = set(FORCED_MARKERS)
    for idx, sb in enumerate(adata_filt.var["symbol_base"]):
        if sb in mset:
            rn = adata_filt.var_names[idx]
            if not adata_filt.var.loc[rn, "highly_variable"]:
                adata_filt.var.loc[rn, "highly_variable"] = True
                n_added += 1
    print(f"  -> Forced {n_added} markers into HVG")

n_hvg_final = int(adata_filt.var["highly_variable"].sum())
print(f"  -> Final HVG count: {n_hvg_final}")

hvg_genes = adata_filt.var_names[adata_filt.var["highly_variable"]].tolist()
(output_dir / f"{OUTPUT_PREFIX}_hvg_genes.txt").write_text("\n".join(hvg_genes))


# %% [markdown]
# ## Cell 10 — Step 7: Capture Full Matrix for .raw

# %%
print("\n" + "=" * 80)
print("[Step 7] Capturing full gene matrix for .raw...")
print("=" * 80)
full_counts = adata_filt.layers["counts"]
if issparse(full_counts) and not isinstance(full_counts, csr_matrix):
    full_counts = csr_matrix(full_counts)
raw_var = adata_filt.var.copy()
print(f"  Full matrix: {full_counts.shape}")


# %% [markdown]
# ## Cell 11 — Step 8: Build Training Subset (HVG)

# %%
print("\n" + "=" * 80)
print("[Step 8] Building training subset (HVG only)...")
print("=" * 80)

hvg_mask    = adata_filt.var["highly_variable"].values
X_hvg       = adata_filt.layers["counts"][:, hvg_mask]
if issparse(X_hvg) and not isinstance(X_hvg, csr_matrix):
    X_hvg = csr_matrix(X_hvg)

adata_train = sc.AnnData(
    X=X_hvg.copy(),
    obs=adata_filt.obs.copy(),
    var=adata_filt.var.iloc[hvg_mask].copy()
)
adata_train.var_names        = adata_filt.var_names[hvg_mask]
adata_train.layers["counts"] = adata_train.X
print(f"  Training data: {adata_train.shape}")

# Assign scANVI training labels from CellTypist filtered column
adata_train.obs["scanvi_labels"] = adata_train.obs[label_source_col].astype(str).astype("category")
if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
    adata_train.obs["scanvi_labels"] = (
        adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])
    )

print("\n  scANVI training label distribution:")
print(adata_train.obs["scanvi_labels"].value_counts())
n_unknown = (adata_train.obs["scanvi_labels"] == UNLABELED_CATEGORY).sum()
print(f"\n  Unknown (low-conf): {n_unknown:,} / {adata_train.n_obs:,} "
      f"({100*n_unknown/adata_train.n_obs:.1f}%)")

gc.collect()


# %% [markdown]
# ## Cell 12 — Step 9: scVI Training

# %%
print("\n" + "=" * 80)
print("[Step 9] Training scVI on filtered myeloid data...")
print("=" * 80)

# Ensure covariates exist; use 0-fallback if not computed previously
required_cont = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
for col in required_cont:
    if col not in adata_train.obs.columns:
        print(f"  [WARN] covariate '{col}' not found — filling with 0.0")
        adata_train.obs[col] = 0.0

scvi.model.SCVI.setup_anndata(
    adata_train,
    layer="counts",
    batch_key=BATCH_KEY,
    continuous_covariate_keys=required_cont,
    categorical_covariate_keys=[TISSUE_KEY]
)

scvi_model = scvi.model.SCVI(
    adata_train,
    n_latent=SCVI_N_LATENT, n_layers=SCVI_N_LAYERS,
    n_hidden=SCVI_N_HIDDEN, dropout_rate=SCVI_DROPOUT
)

train_kwargs = {
    "max_epochs": MAX_EPOCHS_SCVI, "batch_size": BATCH_SIZE,
    "early_stopping": True, "early_stopping_patience": 30,
    "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
}
if gpu_available:
    train_kwargs["accelerator"] = "gpu"; train_kwargs["devices"] = 1

t0 = time.time()
scvi_model.train(**train_kwargs)
print(f"[OK] scVI training complete  ({time.time()-t0:.1f}s)")


# %% [markdown]
# ## Cell 13 — Step 10: scANVI Training

# %%
print("\n" + "=" * 80)
print("[Step 10] Training scANVI (CellTypist labels for all filtered cells)...")
print("=" * 80)

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, adata=adata_train,
    labels_key="scanvi_labels",
    unlabeled_category=UNLABELED_CATEGORY
)

scanvi_train_kwargs = {
    "max_epochs": MAX_EPOCHS_SCANVI, "batch_size": BATCH_SIZE,
    "early_stopping": True, "early_stopping_patience": 20,
    "plan_kwargs": {"lr": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
}
if gpu_available:
    scanvi_train_kwargs["accelerator"] = "gpu"; scanvi_train_kwargs["devices"] = 1

t0 = time.time()
scanvi_model.train(**scanvi_train_kwargs)
print(f"[OK] scANVI training complete  ({time.time()-t0:.1f}s)")


# %% [markdown]
# ## Cell 14 — Step 11: Export Latent & Predictions

# %%
print("\n" + "=" * 80)
print("[Step 11] Exporting latent representations and predictions...")
print("=" * 80)

# scANVI latent — reindex to adata_filt (= adata_train here, but be explicit)
lat    = scanvi_model.get_latent_representation(adata_train)
lat_df = pd.DataFrame(lat, index=adata_train.obs_names,
                      columns=[f"scANVI_{i}" for i in range(lat.shape[1])])
lat_al = lat_df.reindex(adata_filt.obs_names)
if lat_al.isna().any().any():
    raise ValueError("Missing scANVI latent after reindex!")
adata_filt.obsm["X_scANVI"] = lat_al.values

# Hard label predictions
pred_al = (
    pd.Series(scanvi_model.predict(adata_train), index=adata_train.obs_names)
    .reindex(adata_filt.obs_names)
)
adata_filt.obs["scanvi_pred"] = pred_al.values

# Soft probability matrix
proba_raw = scanvi_model.predict(adata_train, soft=True)
if isinstance(proba_raw, pd.DataFrame):
    label_order = list(proba_raw.columns)
    proba       = proba_raw.values.astype(np.float32)
else:
    proba = np.asarray(proba_raw, dtype=np.float32)
    try:
        label_order = list(
            scanvi_model.adata_manager.get_state_registry("labels").categorical_mapping
        )
    except Exception:
        label_order = [f"label_{i}" for i in range(proba.shape[1])]

proba_df = pd.DataFrame(proba, index=adata_train.obs_names, columns=label_order)
proba_al = proba_df.reindex(adata_filt.obs_names)
adata_filt.obsm["scanvi_proba"]      = proba_al.values
adata_filt.obs["scanvi_confidence"]  = proba_al.values.max(axis=1)
adata_filt.uns["scanvi_label_order"] = list(label_order)

print("  scANVI predictions (top 15):")
print(adata_filt.obs["scanvi_pred"].value_counts().head(15))

# Novelty score (entropy-based)
if not np.isfinite(proba_al.values).all():
    raise ValueError("Non-finite values in scANVI probability matrix!")
ents = entropy(proba_al.values + 1e-10, axis=1)
adata_filt.obs["scanvi_entropy"] = ents
emin, emax = ents.min(), ents.max()
adata_filt.obs["novelty_score"] = (ents - emin) / (emax - emin) if emax > emin else 0.0
qmask = adata_filt.obs["data_source"] == "query"
adata_filt.obs["is_potentially_novel"] = (adata_filt.obs["novelty_score"] > 0.7) & qmask
print(f"  High novelty query cells: {adata_filt.obs['is_potentially_novel'].sum()}")


# %% [markdown]
# ## Cell 15 — Step 12: Attach .raw

# %%
print("\n" + "=" * 80)
print("[Step 12] Attaching .raw (full gene matrix)...")
print("=" * 80)
from anndata import AnnData as _AnnData
adata_filt.raw = _AnnData(X=full_counts, obs=adata_filt.obs.copy(), var=raw_var)
print(f"  .raw: {adata_filt.raw.n_vars} genes")


# %% [markdown]
# ## Cell 16 — Step 13: Compute UMAPs

# %%
print("\n" + "=" * 80)
print("[Step 13] Computing UMAPs...")
print("=" * 80)

# Export scVI latent to adata_filt
lat_s    = scvi_model.get_latent_representation(adata_train)
lat_s_df = pd.DataFrame(lat_s, index=adata_train.obs_names,
                         columns=[f"scVI_{i}" for i in range(lat_s.shape[1])])
lat_s_al = lat_s_df.reindex(adata_filt.obs_names)
if lat_s_al.isna().any().any():
    raise ValueError("Missing scVI latent after reindex!")
adata_filt.obsm["X_scVI"] = lat_s_al.values

# Query-only Leiden clustering
run_query_only_leiden(adata_filt, resolution=QUERY_LEIDEN_RESOLUTION)

# scVI UMAP
sc.pp.neighbors(adata_filt, use_rep="X_scVI", n_neighbors=30,
                random_state=RANDOM_SEED, key_added="neighbors_scVI")
sc.tl.umap(adata_filt, random_state=RANDOM_SEED, neighbors_key="neighbors_scVI")
adata_filt.obsm["X_umap_scVI"] = adata_filt.obsm["X_umap"].copy()
print("  -> X_umap_scVI saved")

# scANVI UMAP (default)
sc.pp.neighbors(adata_filt, use_rep="X_scANVI", n_neighbors=30,
                random_state=RANDOM_SEED, key_added="neighbors_scANVI")
sc.tl.umap(adata_filt, random_state=RANDOM_SEED, neighbors_key="neighbors_scANVI")
adata_filt.obsm["X_umap_scANVI"] = adata_filt.obsm["X_umap"].copy()
adata_filt.obsm["X_umap"]        = adata_filt.obsm["X_umap_scANVI"].copy()
print("  -> X_umap_scANVI saved (default X_umap)")

umap_op = UMAP(n_neighbors=30, n_components=2, min_dist=0.5, metric="euclidean",
               random_state=RANDOM_SEED)
umap_op.fit(adata_filt.obsm["X_scANVI"])
joblib.dump(umap_op, output_dir / f"{OUTPUT_PREFIX}_umap_scanvi_operator.joblib")
print("  -> UMAP operator saved")


# %% [markdown]
# ## Cell 17 — Step 14: Run Log

# %%
print("\n" + "=" * 80)
print("[Step 14] Writing run log...")
print("=" * 80)
log_ts   = datetime.now().isoformat()
out_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"

pipeline_log = {
    "version":         "2.1", "timestamp": log_ts,
    "input_h5ad":      INPUT_H5AD,
    "output_dir":      str(output_dir), "output_prefix": OUTPUT_PREFIX,
    "planned_output":  str(out_h5ad),
    "elapsed_min":     round((time.time() - PIPELINE_START) / 60, 2),
    "n_obs_before_filter": n_before,
    "n_obs_after_filter":  adata_filt.n_obs,
    "n_removed":           n_removed,
    "pct_removed":         round(n_removed / n_before * 100, 2),
    "n_vars":          int(adata_filt.n_vars),
    "n_hvg":           n_hvg_final,
    "hvg_method":      hvg_method,
    "gpu_available":   bool(gpu_available),
    "architecture":    "CellTypist-filtered myeloid-only retrain (v2.1)",
    "scanvi_label_source": label_source_col,
    "min_batch_cells": MIN_BATCH_CELLS,
    "small_batches_dropped": small_batches if small_batches else [],
}
anndata_structure = {
    "timestamp": log_ts, "shape": [int(adata_filt.n_obs), int(adata_filt.n_vars)],
    "raw": {"present": True, "n_vars": adata_filt.raw.n_vars},
    "obs_columns": list(adata_filt.obs.columns),
    "obsm_keys":   sorted(adata_filt.obsm.keys()),
    "uns_keys":    sorted(str(k) for k in adata_filt.uns.keys()),
}
adata_filt.uns["pipeline_log"]      = pipeline_log
adata_filt.uns["anndata_structure"] = anndata_structure

(output_dir / f"{OUTPUT_PREFIX}_run_log.txt").write_text(
    "\n".join(f"{k}: {v}" for k, v in pipeline_log.items()), encoding="utf-8"
)
with open(output_dir / f"{OUTPUT_PREFIX}_anndata_structure.json", "w") as f:
    json.dump(anndata_structure, f, indent=2)
print("  -> Log written")


# %% [markdown]
# ## Cell 18 — Step 15: Save

# %%
print("\n" + "=" * 80)
print("[Step 15] Saving results...")
print("=" * 80)

# Cast label columns to category before write_h5ad
_cat_cols = [
    "data_source", "cell_type_fine_ref", "celltypist_lineage",
    "scanvi_labels", "scanvi_pred",
    "mono_subtype_by_score", "leiden_query",
    CELLTYPIST_DIRECT_LABEL_KEY, CELLTYPIST_DIRECT_FILT_KEY,
]
for col in _cat_cols:
    for ad in [adata_filt, adata_train]:
        if col in ad.obs.columns:
            ad.obs[col] = ad.obs[col].astype("category")

scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
scvi_model.save(output_dir   / f"{OUTPUT_PREFIX}_scvi_model",   overwrite=True)

config = {
    "version":   "2.1", "timestamp": datetime.now().isoformat(),
    "input":     INPUT_H5AD,
    "n_hvg":     n_hvg_final, "hvg_method": hvg_method,
    "architecture": {
        "description":         "CellTypist-filtered myeloid retrain — non-myeloid cells removed",
        "filter_policy":       "keep MYELOID_LINEAGES, remove contaminants",
        "myeloid_lineages":    sorted(MYELOID_LINEAGES),
        "keep_other_unknown":  KEEP_OTHER_UNKNOWN,
        "n_removed":           n_removed,
        "pct_removed":         round(n_removed / n_before * 100, 2),
        "scanvi_labels":       label_source_col,
    },
    "scanvi_labels": list(label_order),
    "umap_spaces": {
        "X_umap":        "DEFAULT (scANVI-based)",
        "X_umap_scVI":   "scVI latent UMAP",
        "X_umap_scANVI": "scANVI latent UMAP",
    },
}
with open(output_dir / f"{OUTPUT_PREFIX}_config.json", "w") as f:
    json.dump(config, f, indent=2)

# Statistics CSV
_stats = (
    adata_filt.obs["scanvi_pred"].value_counts()
    .rename_axis("Cell_Type").reset_index(name="Count")
)
_stats["Percentage"] = 100 * _stats["Count"] / _stats["Count"].sum()
_conf = adata_filt.obs.groupby("scanvi_pred")["scanvi_confidence"].agg(["mean","std"])
_stats = _stats.merge(_conf, left_on="Cell_Type", right_index=True, how="left")
_stats.to_csv(output_dir / f"{OUTPUT_PREFIX}_scanvi_statistics.csv", index=False)

# Per-cell annotations CSV
adata_filt.obs[[
    "data_source", "celltypist_lineage", "scanvi_pred",
    CELLTYPIST_DIRECT_LABEL_KEY, CELLTYPIST_DIRECT_FILT_KEY,
    "scanvi_confidence", "novelty_score",
]].to_csv(output_dir / f"{OUTPUT_PREFIX}_annotations.csv")

adata_filt.write_h5ad(out_h5ad, compression="gzip")
print(f"  -> {out_h5ad}")
adata_train.write_h5ad(output_dir / f"{OUTPUT_PREFIX}_train_HVG.h5ad", compression="gzip")
print(f"  -> train HVG saved")


# %% [markdown]
# ## Cell 19 — Step 16: Visualization (v2.9 rasterized)

# %%
print("\n" + "=" * 80)
print("[Step 16] Visualization...")
elapsed = (time.time() - PIPELINE_START) / 60
print(f"  Elapsed so far: {elapsed:.1f} min")
print("=" * 80)


def _umap(adata, color, ax, title, cmap=None, vmin=None, vmax=None,
          basis="X_umap_scANVI", legend_loc="right margin"):
    """Wrapper: sc.pl.embedding with rasterized=True (v2.9 requirement)."""
    kw = dict(
        ax=ax, show=False, title=title, s=8,
        rasterized=True,            # REQUIRED — keeps PDF file size manageable
        legend_loc=legend_loc,
        legend_fontsize=6,
        frameon=False,
    )
    if cmap:  kw["cmap"] = cmap
    if vmin is not None: kw["vmin"] = vmin
    if vmax is not None: kw["vmax"] = vmax
    sc.pl.embedding(adata, basis=basis, color=color, **kw)


# ── Panel 1: Overview 4×4 ────────────────────────────────────────────────
fig = plt.figure(figsize=(24, 20))
gs  = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)

row0_specs = [
    ("data_source",              "Data Source",           None,      None),
    ("celltypist_lineage",       "CellTypist Lineage",    None,      None),
    (CELLTYPIST_DIRECT_LABEL_KEY,"CellTypist Direct",     None,      None),
    (CELLTYPIST_DIRECT_FILT_KEY, "CellTypist Filtered",   None,      None),
]
row1_specs = [
    ("scanvi_pred",              "scANVI Predictions",    None,      None),
    ("scanvi_confidence",        "scANVI Confidence",     "viridis", (0, 1)),
    ("novelty_score",            "Novelty Score",         "hot",     (0, 1)),
    ("mono_subtype_by_score",    "Mono Subtype",          None,      None),
]
for col, (key, title, cmap, vlim) in enumerate(row0_specs):
    if key in adata_filt.obs.columns:
        ax = fig.add_subplot(gs[0, col])
        _umap(adata_filt, key, ax, title, cmap=cmap,
              vmin=vlim[0] if vlim else None, vmax=vlim[1] if vlim else None,
              legend_loc="on data")

for col, (key, title, cmap, vlim) in enumerate(row1_specs):
    if key in adata_filt.obs.columns:
        ax = fig.add_subplot(gs[1, col])
        _umap(adata_filt, key, ax, title, cmap=cmap,
              vmin=vlim[0] if vlim else None, vmax=vlim[1] if vlim else None,
              legend_loc="on data" if not cmap else "right margin")

marker_panels = [
    ("LYZ","LYZ (pan-Myeloid)"),  ("CD14","CD14 (Classical Mono)"),
    ("FCGR3A","FCGR3A (NC-Mono)"),("CD68","CD68 (Macro)"),
    ("CD1C","CD1C (cDC2)"),       ("CLEC9A","CLEC9A (cDC1)"),
    ("S100A8","S100A8 (Neutro)"), ("MKI67","MKI67 (Prolif)"),
]
for i, (gene, title) in enumerate(marker_panels):
    row_i, col_i = 2 + i // 4, i % 4
    ax = fig.add_subplot(gs[row_i, col_i])
    if adata_filt.raw is not None and gene in adata_filt.raw.var_names:
        sc.pl.embedding(
            adata_filt, basis="X_umap_scANVI", color=gene,
            ax=ax, show=False, title=title, cmap="Reds", s=8,
            rasterized=True, frameon=False, use_raw=True,
        )
    else:
        ax.set_title(f"{title} (not found)"); ax.axis("off")

fig.savefig(output_dir / f"{OUTPUT_PREFIX}_overview.pdf",
            dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig)
print(f"  -> {OUTPUT_PREFIX}_overview.pdf")

# ── Panel 2: CellTypist vs scANVI ────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(21, 7))
for ax, (key, title) in zip(axes2, [
    (CELLTYPIST_DIRECT_LABEL_KEY, "CellTypist Direct"),
    (CELLTYPIST_DIRECT_FILT_KEY,  "CellTypist Filtered"),
    ("scanvi_pred",               "scANVI (retrained, myeloid-only)"),
]):
    if key in adata_filt.obs.columns:
        sc.pl.embedding(
            adata_filt, basis="X_umap_scANVI", color=key,
            ax=ax, show=False, title=title,
            legend_loc="on data", legend_fontsize=6, s=8,
            rasterized=True, frameon=False,
        )
plt.tight_layout()
fig2.savefig(output_dir / f"{OUTPUT_PREFIX}_celltypist_vs_scanvi.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig2)
print(f"  -> {OUTPUT_PREFIX}_celltypist_vs_scanvi.pdf")

# ── Panel 3: scVI vs scANVI UMAP comparison ──────────────────────────────
fig3, axes3 = plt.subplots(2, 3, figsize=(18, 12))
for row_i, (basis, label) in enumerate([("X_umap_scVI","scVI"),
                                         ("X_umap_scANVI","scANVI")]):
    for col_i, (color, loc) in enumerate([
        ("data_source",              "right margin"),
        (CELLTYPIST_DIRECT_LABEL_KEY, "on data"),
        ("scanvi_pred",               "on data"),
    ]):
        sc.pl.embedding(
            adata_filt, basis=basis, color=color,
            ax=axes3[row_i, col_i], show=False,
            title=f"{color} ({label})",
            legend_loc=loc, legend_fontsize=6, s=8,
            rasterized=True, frameon=False,
        )
plt.tight_layout()
fig3.savefig(output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig3)
print(f"  -> {OUTPUT_PREFIX}_umap_comparison.pdf")

# ── Panel 4: Filter QC — before vs after ─────────────────────────────────
# Pie chart showing how many cells were removed per lineage
fig4, axes4 = plt.subplots(1, 2, figsize=(14, 6))

# Left: full lineage distribution (from adata_filt + removed cells counts)
# We reconstruct from the lineage_counts we computed before filtering
all_lin = adata_filt.obs["celltypist_lineage"].value_counts()

axes4[0].pie(
    all_lin.values,
    labels=all_lin.index,
    autopct="%1.1f%%",
    startangle=140,
)
axes4[0].set_title(f"Retained Myeloid Lineages\n(n={n_after:,} cells)")

# Right: bar chart of removed contaminants (reconstruct from n_removed per label)
axes4[1].text(0.5, 0.5,
    f"Total removed: {n_removed:,} cells\n"
    f"({n_removed/n_before*100:.1f}% of input)\n\n"
    f"Input:    {n_before:,}\n"
    f"Retained: {n_after:,}\n"
    f"Batches dropped: {len(small_batches)}",
    ha="center", va="center", transform=axes4[1].transAxes, fontsize=14
)
axes4[1].set_title("Filtering Summary")
axes4[1].axis("off")

plt.tight_layout()
fig4.savefig(output_dir / f"{OUTPUT_PREFIX}_filter_qc.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig4)
print(f"  -> {OUTPUT_PREFIX}_filter_qc.pdf")

# ── Final summary ─────────────────────────────────────────────────────────
elapsed = (time.time() - PIPELINE_START) / 60
print("\n" + "=" * 80)
print("MYELOID PIPELINE v2.1 COMPLETE")
print("=" * 80)
print(f"  Input cells:   {n_before:,}")
print(f"  Removed:       {n_removed:,}  ({n_removed/n_before*100:.1f}%)")
print(f"  Retained:      {n_after:,}")
print(f"  Batches used:  {adata_filt.obs[BATCH_KEY].nunique()}")
print(f"  Elapsed:       {elapsed:.1f} min")
print(f"\nscANVI label source: {label_source_col}")
print(f"  Unknown (low-conf): "
      f"{int((adata_filt.obs.get('scanvi_labels', adata_train.obs['scanvi_labels']).astype(str) == UNLABELED_CATEGORY).sum()):,}")
print("\nTop scANVI predictions:")
print(adata_filt.obs["scanvi_pred"].value_counts().head(10))
print("\nTop CellTypist labels (post-filter):")
print(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].value_counts().head(10))
print("=" * 80)
