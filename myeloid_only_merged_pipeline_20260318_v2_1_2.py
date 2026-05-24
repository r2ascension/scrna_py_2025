# %% [markdown]
# # Myeloid Merged Pipeline v2.1.1 HOTFIX — CellTypist-Filtered Retrain
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
#
# **HOTFIX v2.1.1 — all fixes vs v2.1:**
#   P0-1: scanvi_labels written back to adata_filt via reindex (was only in adata_train)
#   P0-2: CELLTYPIST_DIRECT_FILT_KEY always rebuilt — never reuse stale checkpoint column
#   P0-3: Missing covariates raise ValueError instead of silent 0-fill
#   P1-4: Unmapped CellTypist labels raise ValueError — no silent contaminant deletion
#   P1-6: Erythroid / Megakaryocyte removed from MYELOID_LINEAGES (non-myeloid lineages)
#   P2-7: filter_qc.pdf Panel 4 left pie uses full pre-filter lineage_counts_all
#   P2-8: assert counts layer present immediately after load
#   P2-9: old HVG auxiliary columns cleared before re-running HVG selection
#   MISC:  final summary uses adata_filt.obs["scanvi_labels"] directly (no .get() workaround)

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
print("Myeloid Merged Pipeline v2.1.1 HOTFIX  (CellTypist-Filtered Retrain)")
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

# Lineage categories that should be KEPT (true myeloid targets).
# NOTE: Erythroid and Megakaryocyte are intentionally excluded — they are
# hematopoietic but not myeloid in the classical immune sense (monocyte/DC/granulocyte).
# If your study design requires them, add them back explicitly.
MYELOID_LINEAGES = {"Monocyte", "Macrophage", "DC", "Mast", "Basophil",
                    "Eosinophil", "Neutrophil"}

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

# P1-4 guard behaviour:
#   False (default) → raise ValueError if any CellTypist label is absent from
#                     CELLTYPIST_LINEAGE_MAP (safe production default)
#   True            → print warning and continue (temporary override only)
KEEP_OTHER_UNKNOWN = False

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

PIPELINE_VERSION = "2.1.2"
PIPELINE_START   = time.time()
FIGURE_DPI       = 300
FIGURE_FORMAT    = "pdf"

# Rare-label collapse: scANVI labels with fewer cells than this threshold
# are merged into UNLABELED_CATEGORY before training to avoid instability.
MIN_CELLS_PER_SCANVI_LABEL = 20


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

# P2-8: fail fast if counts layer is missing (do not let it surface later)
assert "counts" in adata.layers, \
    "CRITICAL: 'counts' layer not found in input checkpoint. " \
    "Ensure v2.0 preprocessed h5ad was saved before HVG subsetting."

# P1 guard (v2.1.2): verify all required obs columns exist up-front
required_obs = [BATCH_KEY, TISSUE_KEY, "data_source",
                CELLTYPIST_LABEL_COL, CELLTYPIST_CONF_KEY]
missing_obs = [c for c in required_obs if c not in adata.obs.columns]
if missing_obs:
    raise ValueError(f"Missing required obs columns in checkpoint: {missing_obs}")

# P1 guard (v2.1.2): obs_names must be unique — reindex alignment is unsafe otherwise
if not adata.obs_names.is_unique:
    dup_n = int(adata.obs_names.duplicated().sum())
    raise ValueError(
        f"obs_names are not unique ({dup_n} duplicates found). "
        "Fix barcode uniqueness upstream before any reindex-based alignment."
    )
if not adata.var_names.is_unique:
    print("  [WARN] var_names not unique — calling var_names_make_unique()")
    adata.var_names_make_unique()

# Verify CellTypist labels are present
assert CELLTYPIST_LABEL_COL in adata.obs.columns, \
    f"CellTypist column '{CELLTYPIST_LABEL_COL}' not found. " \
    f"Available: {list(adata.obs.columns)}"
assert CELLTYPIST_CONF_KEY in adata.obs.columns, \
    f"Confidence column '{CELLTYPIST_CONF_KEY}' not found."
print(f"\n[OK] CellTypist labels present: {adata.obs[CELLTYPIST_LABEL_COL].nunique()} unique types")
print(f"[OK] obs_names unique: {adata.obs_names.is_unique}")


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
# Store full pre-filter counts for filter_qc.pdf (P2-7)
lineage_counts_all = lineage_counts.copy()
total = len(adata)
for lin, n in lineage_counts.items():
    flag = "  [KEEP]" if lin in MYELOID_LINEAGES else "  [REMOVE]"
    print(f"  {lin:<30} {n:>7,}  ({n/total*100:.1f}%){flag}")

# Warn about unmapped labels (appear as "other_unknown")
# Warn about unmapped labels (appear as "other_unknown")
unmapped_labels = ct_labels[lineage_series == "other_unknown"].unique()
if len(unmapped_labels) > 0:
    msg_lines = [
        f"  Unmapped CellTypist labels detected ({len(unmapped_labels)} unique):",
    ]
    for lbl in sorted(unmapped_labels):
        n = int((ct_labels == lbl).sum())
        msg_lines.append(f"    '{lbl}'  n={n:,}")
    msg_lines += [
        "",
        "  ACTION REQUIRED: Add each label to CELLTYPIST_LINEAGE_MAP before re-running.",
        "  If these labels are myeloid → assign a myeloid lineage.",
        "  If these labels are contaminants → assign a 'contaminant_*' lineage.",
        "  Only set KEEP_OTHER_UNKNOWN=True as a temporary override — it bypasses this guard.",
    ]
    if not KEEP_OTHER_UNKNOWN:
        raise ValueError(
            "P1-4 GUARD: Unmapped CellTypist labels found. "
            "Update CELLTYPIST_LINEAGE_MAP to classify them.\n" +
            "\n".join(msg_lines)
        )
    else:
        print("\n  [WARN] KEEP_OTHER_UNKNOWN=True — skipping unmapped label guard.")
        for line in msg_lines:
            print(line)


# %% [markdown]
# ## Cell 6 — Step 3: Filter Non-Myeloid Cells

# %%
print("\n" + "=" * 80)
print("[Step 3] Filtering to myeloid cells only...")
print("=" * 80)

n_before = adata.n_obs
keep_mask = adata.obs["celltypist_lineage"].isin(MYELOID_LINEAGES)

# KEEP_OTHER_UNKNOWN is set in Cell 1 configuration — do not redefine here.
if KEEP_OTHER_UNKNOWN:
    keep_mask |= (adata.obs["celltypist_lineage"] == "other_unknown")
    print("  [NOTE] other_unknown cells are INCLUDED (KEEP_OTHER_UNKNOWN=True)")

adata_filt        = adata[keep_mask].copy()
n_after_lineage   = adata_filt.n_obs
n_removed_lineage = n_before - n_after_lineage

print(f"  Before:          {n_before:,} cells")
print(f"  Removed (lineage filter): {n_removed_lineage:,} ({n_removed_lineage/n_before*100:.1f}%)")
print(f"  After lineage:   {n_after_lineage:,} cells")

print("\n  Remaining lineage breakdown:")
for lin, n in adata_filt.obs["celltypist_lineage"].value_counts().items():
    print(f"    {lin:<30} {n:>7,}  ({n/n_after_lineage*100:.1f}%)")

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

# v2.1.2: update final cell counts after BOTH lineage filter AND batch drop
n_after_final    = adata_filt.n_obs
n_removed_batch  = n_after_lineage - n_after_final
n_removed_total  = n_before - n_after_final
print(f"\n  Removed by lineage filter:  {n_removed_lineage:,}")
print(f"  Removed by small-batch drop: {n_removed_batch:,}")
print(f"  Total removed:               {n_removed_total:,} ({n_removed_total/n_before*100:.1f}%)")
print(f"  Final retained:              {n_after_final:,}")

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
print("[Step 5] Rebuilding CellTypist filtered label column (always, never reuse checkpoint)...")
print("=" * 80)
# P0-2: Always recompute from current CELLTYPIST_CONF_THRESHOLD.
# Never reuse the column from the checkpoint — the threshold may have changed
# and a stale column would silently produce wrong Unknown assignments.
raw_labels = adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].astype(str)
conf       = adata_filt.obs[CELLTYPIST_CONF_KEY].astype(float).values
filt       = np.where(conf >= CELLTYPIST_CONF_THRESHOLD, raw_labels.values, UNLABELED_CATEGORY)
adata_filt.obs[CELLTYPIST_DIRECT_FILT_KEY] = pd.Series(
    filt, index=adata_filt.obs_names, dtype="object"
).astype("category")
n_low = int((conf < CELLTYPIST_CONF_THRESHOLD).sum())
print(f"  Applied threshold: {CELLTYPIST_CONF_THRESHOLD}")
print(f"  -> {n_low:,} low-conf cells set to '{UNLABELED_CATEGORY}' "
      f"({n_low/adata_filt.n_obs*100:.1f}%)")

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

# P2-9: Clear ALL old HVG auxiliary columns so re-selection starts clean.
# Stale columns from v2.0 (computed on all cells) must not bleed into new results.
_hvg_aux_cols = ["highly_variable", "highly_variable_rank",
                 "highly_variable_nbatches", "highly_variable_intersection",
                 "means", "variances", "variances_norm",
                 "dispersions", "dispersions_norm"]
dropped = [c for c in _hvg_aux_cols if c in adata_filt.var.columns]
if dropped:
    adata_filt.var.drop(columns=dropped, inplace=True)
    print(f"  -> Cleared old HVG columns: {dropped}")

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

# v2.1.2: collapse rare labels into Unknown to stabilise scANVI classifier head.
# Labels with < MIN_CELLS_PER_SCANVI_LABEL cells (excluding Unknown itself) are merged.
label_counts_train = adata_train.obs["scanvi_labels"].value_counts()
rare_labels = label_counts_train[
    (label_counts_train < MIN_CELLS_PER_SCANVI_LABEL) &
    (label_counts_train.index != UNLABELED_CATEGORY)
].index.tolist()
if rare_labels:
    print(f"\n  [v2.1.2] Collapsing {len(rare_labels)} rare labels "
          f"(< {MIN_CELLS_PER_SCANVI_LABEL} cells) → '{UNLABELED_CATEGORY}':")
    for lbl in rare_labels:
        print(f"    '{lbl}'  n={int(label_counts_train[lbl])}")
    adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype(str)
    adata_train.obs.loc[
        adata_train.obs["scanvi_labels"].isin(rare_labels), "scanvi_labels"
    ] = UNLABELED_CATEGORY
    adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")
    if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
        adata_train.obs["scanvi_labels"] = (
            adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])
        )
    print(f"  -> After rare-label collapse: "
          f"{adata_train.obs['scanvi_labels'].nunique()} labels "
          f"({(adata_train.obs['scanvi_labels']==UNLABELED_CATEGORY).sum():,} Unknown)")
else:
    print(f"\n  [v2.1.2] No rare labels (all >= {MIN_CELLS_PER_SCANVI_LABEL} cells)")

# P0-1: Write scanvi_labels back to the final object adata_filt.
# adata_train is a gene-subset of adata_filt sharing the same obs_names,
# so reindex is exact. Without this, adata_filt.write_h5ad() would lack
# the training label column.
adata_filt.obs["scanvi_labels"] = (
    adata_train.obs["scanvi_labels"]
    .reindex(adata_filt.obs_names)
    .astype("category")
)
print(f"  -> scanvi_labels written back to adata_filt "
      f"({adata_filt.obs['scanvi_labels'].nunique()} categories incl. Unknown)")

gc.collect()


# %% [markdown]
# ## Cell 12 — Step 9: scVI Training

# %%
print("\n" + "=" * 80)
print("[Step 9] Training scVI on filtered myeloid data...")
print("=" * 80)

# P0-3: Covariates must exist in the v2.0 checkpoint.
# Filling with 0 is NOT acceptable — missing covariates indicate an upstream
# problem that would silently invalidate scVI's batch/covariate modelling.
required_cont = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
missing_cont  = [col for col in required_cont if col not in adata_train.obs.columns]
if missing_cont:
    raise ValueError(
        f"P0-3 GUARD: Required covariates missing from checkpoint: {missing_cont}. "
        "Re-run v2.0 pipeline to ensure prepare_covariates() was called before saving "
        "the preprocessed checkpoint."
    )
print(f"  [OK] All covariates present: {required_cont}")

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
print("[Step 12] Attaching .raw (full gene matrix — log1p for display)...")
print("=" * 80)
from anndata import AnnData as _AnnData

# v2.1.2 fix (P1): .raw should contain normalized log1p expression, not raw counts.
# Storing counts in .raw and then using use_raw=True in plots produces
# library-size-confounded marker visualizations.
# Priority: log1p layer > compute from counts > fall back to counts with warning.
if "log1p" in adata_filt.layers:
    full_expr = adata_filt.layers["log1p"]
    raw_source = "log1p layer"
elif adata_filt.X is not None:
    import numpy as _np
    _x_sample = adata_filt.X[:5].toarray() if issparse(adata_filt.X) else _np.asarray(adata_filt.X[:5])
    if _x_sample.max() < 30:   # already log1p scale
        full_expr  = adata_filt.X
        raw_source = "adata_filt.X (detected as log1p)"
    else:
        print("  [WARN] No log1p layer found and .X appears to be counts. "
              "Computing log1p from full_counts for .raw.")
        import scipy.sparse as _sp
        full_expr_tmp = full_counts.copy()
        _sc_tmp = sc.AnnData(X=full_expr_tmp)
        sc.pp.normalize_total(_sc_tmp, target_sum=1e4)
        sc.pp.log1p(_sc_tmp)
        full_expr  = _sc_tmp.X
        del _sc_tmp
        raw_source = "log1p computed from counts"
else:
    print("  [WARN] No log1p available — falling back to counts for .raw. "
          "Marker plots will use raw counts (not recommended).")
    full_expr  = full_counts
    raw_source = "raw counts (fallback)"

if issparse(full_expr) and not isinstance(full_expr, csr_matrix):
    full_expr = csr_matrix(full_expr)

adata_filt.raw = _AnnData(X=full_expr, obs=adata_filt.obs.copy(), var=raw_var)
print(f"  .raw source: {raw_source}")
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
    "version":         PIPELINE_VERSION, "timestamp": log_ts,
    "input_h5ad":      INPUT_H5AD,
    "output_dir":      str(output_dir), "output_prefix": OUTPUT_PREFIX,
    "planned_output":  str(out_h5ad),
    "elapsed_min":     round((time.time() - PIPELINE_START) / 60, 2),
    "n_obs_input":            n_before,
    "n_obs_after_lineage":    n_after_lineage,
    "n_obs_final":            n_after_final,
    "n_removed_lineage":      n_removed_lineage,
    "n_removed_batch":        n_removed_batch,
    "n_removed_total":        n_removed_total,
    "pct_removed_total":      round(n_removed_total / n_before * 100, 2),
    "n_vars":          int(adata_filt.n_vars),
    "n_hvg":           n_hvg_final,
    "hvg_method":      hvg_method,
    "gpu_available":   bool(gpu_available),
    "architecture":    f"CellTypist-filtered myeloid-only retrain (v{PIPELINE_VERSION})",
    "scanvi_label_source": label_source_col,
    "min_batch_cells": MIN_BATCH_CELLS,
    "min_cells_per_scanvi_label": MIN_CELLS_PER_SCANVI_LABEL,
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
    "version":   PIPELINE_VERSION, "timestamp": datetime.now().isoformat(),
    "input":     INPUT_H5AD,
    "n_hvg":     n_hvg_final, "hvg_method": hvg_method,
    "architecture": {
        "description":         "CellTypist-filtered myeloid retrain — non-myeloid cells removed",
        "filter_policy":       "keep MYELOID_LINEAGES, remove contaminants",
        "myeloid_lineages":    sorted(MYELOID_LINEAGES),
        "keep_other_unknown":  KEEP_OTHER_UNKNOWN,
        "n_removed_lineage":   n_removed_lineage,
        "n_removed_batch":     n_removed_batch,
        "n_removed_total":     n_removed_total,
        "pct_removed_total":   round(n_removed_total / n_before * 100, 2),
        "scanvi_labels":       label_source_col,
        "min_cells_per_scanvi_label": MIN_CELLS_PER_SCANVI_LABEL,
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


# v3.6 FIX: rasterized=True must NOT be passed directly to sc.pl.embedding.
# scanpy internally passes rasterized via functools.partial, causing:
# TypeError: multiple values for keyword argument 'rasterized'
# Correct approach: set sc.settings.vector_friendly = True once before the block.
sc.settings.vector_friendly = True


def _umap(adata, color, ax, title, cmap=None, vmin=None, vmax=None,
          basis="X_umap_scANVI", legend_loc="right margin"):
    """Wrapper for sc.pl.embedding. Rasterization handled by sc.settings.vector_friendly."""
    kw = dict(
        ax=ax, show=False, title=title, s=8,
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
            frameon=False, use_raw=True,
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
            frameon=False,
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
            frameon=False,
        )
plt.tight_layout()
fig3.savefig(output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig3)
print(f"  -> {OUTPUT_PREFIX}_umap_comparison.pdf")

# ── Panel 4: Filter QC — before vs after (P2-7 fix) ──────────────────────
# Left:  ALL lineages from pre-filter data (includes contaminants)
# Right: retained myeloid lineages only + filtering summary text
fig4, axes4 = plt.subplots(1, 2, figsize=(16, 7))

# Left: full pre-filter lineage breakdown (lineage_counts_all captured at Step 2)
_all_vals   = lineage_counts_all.values
_all_labels = lineage_counts_all.index.tolist()
_colors     = ["#d62728" if "contaminant" in l or l == "other_unknown"
               else "#1f77b4" for l in _all_labels]
axes4[0].pie(
    _all_vals,
    labels=_all_labels,
    autopct="%1.1f%%",
    startangle=140,
    colors=_colors,
    textprops={"fontsize": 7},
)
axes4[0].set_title(
    f"All Lineages — Before Filter\n(n={n_before:,} cells)\n"
    "Blue = myeloid kept  |  Red = contaminants removed",
    fontsize=9
)

# Right: retained myeloid only
_ret_lin = adata_filt.obs["celltypist_lineage"].value_counts()
axes4[1].pie(
    _ret_lin.values,
    labels=_ret_lin.index,
    autopct="%1.1f%%",
    startangle=140,
    textprops={"fontsize": 7},
)
axes4[1].set_title(
    f"Retained Myeloid Lineages — Final\n"
    f"(n={n_after_final:,} cells  |  total removed {n_removed_total:,} = {n_removed_total/n_before*100:.1f}%)\n"
    f"Lineage: {n_removed_lineage:,}  |  Small-batch: {n_removed_batch:,}  |  Batches dropped: {len(small_batches)}",
    fontsize=9
)

plt.tight_layout()
fig4.savefig(output_dir / f"{OUTPUT_PREFIX}_filter_qc.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig4)
print(f"  -> {OUTPUT_PREFIX}_filter_qc.pdf")

# ── Final summary ─────────────────────────────────────────────────────────
elapsed = (time.time() - PIPELINE_START) / 60
print("\n" + "=" * 80)
print(f"MYELOID PIPELINE v{PIPELINE_VERSION} COMPLETE")
print("=" * 80)
print(f"  Input cells:              {n_before:,}")
print(f"  Removed (lineage filter): {n_removed_lineage:,} ({n_removed_lineage/n_before*100:.1f}%)")
print(f"  Removed (small batches):  {n_removed_batch:,} ({n_removed_batch/n_before*100:.1f}%)")
print(f"  Total removed:            {n_removed_total:,} ({n_removed_total/n_before*100:.1f}%)")
print(f"  Final retained:           {n_after_final:,}")
print(f"  Batches used:             {adata_filt.obs[BATCH_KEY].nunique()}")
print(f"  Elapsed:                  {elapsed:.1f} min")
print(f"\nscANVI label source: {label_source_col}")
# P0-1: scanvi_labels is now written back to adata_filt — no .get() fallback needed
n_unknown_final = int(
    (adata_filt.obs["scanvi_labels"].astype(str) == UNLABELED_CATEGORY).sum()
)
print(f"  Unknown (low-conf): {n_unknown_final:,}")
print("\nTop scANVI predictions:")
print(adata_filt.obs["scanvi_pred"].value_counts().head(10))
print("\nTop CellTypist labels (post-filter):")
print(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].value_counts().head(10))
print("=" * 80)
