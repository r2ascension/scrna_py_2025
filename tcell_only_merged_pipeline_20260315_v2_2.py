#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
T/NK Cell Merged Pipeline v2.1.2 — CellTypist-Filtered Retrain
================================================================
Purpose:
  Downstream of v2.1. Loads the v2.1 preprocessed checkpoint (already has
  CellTypist labels, HVG flags, covariates, no .raw yet), removes non-T/NK
  contaminants based on CellTypist lineage classification, then retrains
  scVI + scANVI on the clean T/NK-only subset.

Input:
  tcell_merged_v2_adata_preprocessed.h5ad  (from v2.1 Cell 13 / Step 10)
  — already contains CellTypist labels, HVG flags, covariates, no .raw yet

Filter logic:
  - Keep cells whose CellTypist label maps to a T/NK lineage
  - Contaminants (Epithelial, Endothelial, Mast, B, Myeloid ...) are removed
  - After filtering: drop batches with < MIN_BATCH_CELLS cells

Architecture (same as myeloid v2.1.1):
  - scANVI labels = CellTypist filtered labels for ALL surviving cells
  - Low-confidence cells (< threshold) -> Unknown
  - NK cells NEVER enter CD4/CD8 scoring (GNLY/NKG7 overlap fix)
  - Labels with < MIN_CELLS_PER_LABEL cells -> Unknown before training
  - Post-training label audit saved to CSV

Known CellTypist labels (from v2.0/v2.1 run, Immune_All_Low.pkl):
  CD16+ NK cells, CD16- NK cells, CRTAM+ gamma-delta T cells,
  Cycling T cells, Endothelial cells*, Epithelial cells*, Follicular helper T cells,
  ILC3, MAIT cells, Mast cells*, NK cells, Regulatory T cells,
  Tcm/Naive helper T cells, Tem/Effector helper T cells,
  Tem/Temra cytotoxic T cells, Tem/Trm cytotoxic T cells,
  Trm cytotoxic T cells, Type 1 helper T cells, Type 17 helper T cells,
  gamma-delta T cells
  (* = contaminants, removed)

v2.9 compliance:
  - rasterized=True on all UMAP scatter plots, dpi=300

v2.1.1 HOTFIX — all fixes vs v2.1:
  P0-1: scanvi_labels written back to adata_filt via reindex
  P0-2: CELLTYPIST_DIRECT_FILT_KEY always rebuilt from current threshold
  P0-3: Missing covariates raise ValueError (no silent 0-fill)
  P0-4: NK cells excluded from CD4/CD8 score classification
  P0-5: normalize_label_strings + MIN_CELLS_PER_LABEL guard before training
  P1-6: Unmapped CellTypist labels raise ValueError (KEEP_OTHER_UNKNOWN guard)
  P2-7: filter_qc.pdf uses full pre-filter lineage_counts_all for left pie
  P2-8: assert counts layer present immediately after load
  P2-9: old HVG auxiliary columns cleared before re-running HVG selection
  AUDIT: post-training label comparison -> label_audit.csv

v2.1.2 — additional fixes from code review:
  FIX-A: label normalization moved to Step 2 (before CELLTYPIST_LINEAGE_MAP lookup)
         so whitespace issues never trigger P1-6 guard spuriously
  FIX-B: Cycling NK cells mapped to dedicated "cycling_NK" lineage; NK exclusion
         mask now covers both "NK" and "cycling_NK" — cycling NK cells no longer
         enter CD4/CD8 score classification
  FIX-C: .raw attached BEFORE CD4/CD8 scoring so score_genes() always uses
         log-normalized full-gene expression (not implicit adata_filt.X)
  FIX-D: uns_keys recorded AFTER pipeline_log/anndata_structure are written
         so the JSON matches the actual uns contents
  FIX-E: removed unused fig_dir mkdir and unused _ct_labels variable

Author: r2end
Date: 2026-03-15
Version: 2.1.2
"""

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
print("T/NK Cell Merged Pipeline v2.1.2  (CellTypist-Filtered Retrain)")
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

# Input: preprocessed checkpoint from v2.1 (has CellTypist labels, HVG flags)
INPUT_H5AD = "/home/h2048/data/py/20260310/tcell_only_merged_pipeline_v2/tcell_merged_v2_adata_preprocessed.h5ad"

OUTPUT_DIR    = "/home/h2048/data/py/20260310/tcell_only_merged_pipeline_v2_1_1_filtered"
OUTPUT_PREFIX = "tcell_merged_v2_1_1"

BATCH_KEY  = "sample"
TISSUE_KEY = "tissue"

# CellTypist column to use for lineage filtering (from v2.1 output)
CELLTYPIST_LABEL_COL  = "celltypist_label_direct"   # unfiltered majority-vote
CELLTYPIST_CONF_KEY   = "celltypist_confidence"
CELLTYPIST_CONF_THRESHOLD = 0.5

# Minimum cells per batch after filtering (batches below this are dropped)
MIN_BATCH_CELLS = 10

# ==============================================================================
# LINEAGE FILTER MAP
# Cells are KEPT if their CellTypist label maps to a T/NK lineage.
# All other lineages (contaminants) are removed.
# Labels are populated from the actual v2.0/v2.1 run output
# (Immune_All_Low.pkl on this dataset).
# ==============================================================================
CELLTYPIST_LINEAGE_MAP = {
    # ── CD4 T cells: KEEP ─────────────────────────────────────────────────
    "Regulatory T cells":             "CD4",
    "Tcm/Naive helper T cells":       "CD4",
    "Tem/Effector helper T cells":    "CD4",
    "Follicular helper T cells":      "CD4",
    "Type 1 helper T cells":          "CD4",
    "Type 17 helper T cells":         "CD4",
    # ── CD8 T cells: KEEP ─────────────────────────────────────────────────
    "Tem/Temra cytotoxic T cells":    "CD8",
    "Tem/Trm cytotoxic T cells":      "CD8",
    "Trm cytotoxic T cells":          "CD8",
    "Tcm/Naive cytotoxic T cells":    "CD8",
    # ── NK cells: KEEP ────────────────────────────────────────────────────
    "NK cells":                       "NK",
    "CD16+ NK cells":                 "NK",
    "CD16- NK cells":                 "NK",
    # ── gdT / MAIT / ILC: KEEP ────────────────────────────────────────────
    "gamma-delta T cells":            "gdT",
    "CRTAM+ gamma-delta T cells":     "gdT",
    "MAIT cells":                     "MAIT",
    "ILC3":                           "ILC",
    # ── Cycling T/NK: KEEP ────────────────────────────────────────────────
    "Cycling T cells":                "cycling_T",
    # [v2.1.2 FIX-B] Cycling NK cells get their OWN lineage tag "cycling_NK",
    # NOT "cycling_T". This ensures the NK exclusion mask in Step 12 can
    # correctly identify them via celltypist_lineage. In v2.1.1 they were
    # mapped to "cycling_T", so _nk_mask (lineage == "NK") missed them and
    # they entered CD4/CD8 scoring with high GNLY/NKG7 expression.
    "Cycling NK cells":               "cycling_NK",
    # ── Contaminants: REMOVE ──────────────────────────────────────────────
    "Epithelial cells":               "contaminant_epi",
    "Endothelial cells":              "contaminant_stromal",
    "Fibroblasts":                    "contaminant_stromal",
    "Smooth muscle cells":            "contaminant_stromal",
    "Pericytes":                      "contaminant_stromal",
    "Mast cells":                     "contaminant_mast",
    "Memory B cells":                 "contaminant_B",
    "Naive B cells":                  "contaminant_B",
    "Plasma cells":                   "contaminant_B",
    "Plasmablasts":                   "contaminant_B",
    "Classical monocytes":            "contaminant_myeloid",
    "Non-classical monocytes":        "contaminant_myeloid",
    "Macrophages":                    "contaminant_myeloid",
    "DC1":                            "contaminant_myeloid",
    "DC2":                            "contaminant_myeloid",
    "pDC":                            "contaminant_myeloid",
    "Plasmacytoid DCs":               "contaminant_myeloid",
    "Neutrophils":                    "contaminant_myeloid",
}

# Lineage categories that should be KEPT (true T/NK targets).
# [v2.1.2 FIX-B] cycling_NK is a separate entry from cycling_T so the
# NK exclusion mask can match it independently in Step 12.
TNK_LINEAGES = {"CD4", "CD8", "NK", "cycling_NK", "gdT", "MAIT", "ILC", "cycling_T"}

# P1-6 guard behaviour:
#   False (default) -> raise ValueError if any CellTypist label is absent
#                      from CELLTYPIST_LINEAGE_MAP
#   True            -> print warning and continue (temporary override only)
KEEP_OTHER_UNKNOWN = False

# ==============================================================================
# scVI / scANVI parameters
# ==============================================================================
N_HVG               = 4000
FORCE_MARKERS_IN_HVG = True

# [v2.1] Increased from 100 -> 150 for better T/NK subtype resolution
SCVI_N_LATENT   = 150
SCVI_N_LAYERS   = 2
SCVI_N_HIDDEN   = 128
SCVI_DROPOUT    = 0.1
MAX_EPOCHS_SCVI = 400

MAX_EPOCHS_SCANVI  = 200
UNLABELED_CATEGORY = "Unknown"

BATCH_SIZE    = 256
LEARNING_RATE = 1e-3
WEIGHT_DECAY  = 0.0

# scANVI label source
SCANVI_LABEL_SOURCE          = "filtered"   # use confidence-filtered column
CELLTYPIST_DIRECT_LABEL_KEY  = "celltypist_label_direct"
CELLTYPIST_DIRECT_FILT_KEY   = "celltypist_label_direct_filt"

# [v2.1 FIX-P0] Labels with < MIN_CELLS_PER_LABEL cells -> Unknown before training.
# Prevents scANVI from silently absorbing rare labels into dominant classes.
MIN_CELLS_PER_LABEL = 50

QUERY_LEIDEN_RESOLUTION = 1.0

# ==============================================================================
# MARKER GENES
# ==============================================================================

# CD4/CD8 scoring — used for QC visualization only
CD4_SCORE_GENES = ["CD4","IL7R","MAL","TCF7","SELL","CCR7","FOXP3",
                   "IL2RA","IKZF2","CTLA4","CXCR5","BCL6","PDCD1",
                   "IL21","RORC","CCR6","RORA","IL17A","TBX21","IFNG"]
CD8_SCORE_GENES = ["CD8A","CD8B","GZMB","GZMK","GZMH","GNLY","PRF1",
                   "NKG7","KLRG1","CX3CR1","FGFBP2","FCGR3A",
                   "EOMES","TBX21","ZEB2","PDCD1","LAG3","HAVCR2"]
NK_SCORE_GENES  = ["NCAM1","KLRB1","KLRC1","KLRD1","KLRF1","FCGR3A",
                   "GNLY","NKG7","PRF1","GZMB"]

CD4_SCORE_THRESH = 0.3
CD8_SCORE_THRESH = 0.3

# Forced into HVG
_ILC_MARKERS        = ["RORC","IL22","KIT","IL1RL1","GATA3","IL13","IL5","IL4"]
_TFH_MARKERS        = ["CXCR5","BCL6","PDCD1","IL21","ICOS","SH2D1A","ASCL2"]
_TH17_MARKERS       = ["RORC","IL17A","IL17F","CCR6","RORA","IL23R"]
_TH2_MARKERS        = ["GATA3","IL4","IL5","IL13","CRTH2","IL4R","CCR4"]
_TH1_MARKERS        = ["TBX21","IFNG","TNF","CXCR3","IL12RB2","STAT1"]
_CYTOTOXIC_MARKERS  = ["GZMB","GZMK","GZMH","GZMA","GZMM","PRF1","FASLG","TNFSF10"]
_MAIT_MARKERS       = ["SLC4A10","TRAV1-2","NCR3","KLRB1","IL18R1","CXCR6"]
_GDT_MARKERS        = ["TRDC","TRGC1","TRGC2","TRDV1","TRDV2","TRDV3"]
_NK_SUBTYPE_MARKERS = ["NCAM1","FCGR3A","KLRC1","KLRC2","KLRD1","KLRF1",
                       "KLRB1","KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2",
                       "CD57","B3GAT1","FGFBP2","CX3CR1","S1PR5"]
_TREG_MARKERS       = ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"]
_TRM_MARKERS        = ["ITGAE","CXCR6","CD69","ITGA1","ZNF683"]
_EXHAUSTION_MARKERS = ["HAVCR2","TIGIT","LAG3","PDCD1","TOX","NR4A1"]
_NAIVE_MARKERS      = ["CCR7","SELL","TCF7","LEF1","CD27","CD28"]
_EFFECTOR_MARKERS   = ["GZMB","PRF1","IFNG","TNF","CX3CR1","FGFBP2"]
_ACTIVATION_MARKERS = ["HLA-DRA","HLA-DRB1","CD38","CD69","ICOS","ICOSLG",
                       "TNFRSF4","TNFRSF9","TNFRSF18","ENTPD1","LAYN"]
_MEMORY_MARKERS     = ["IL7R","CD27","CD28","TCF7","LEF1","SELL","CCR7",
                       "CXCR3","CX3CR1","KLRG1","S1PR1","KLF2","ID3"]

FORCED_MARKERS = list(set(
    CD4_SCORE_GENES + CD8_SCORE_GENES + NK_SCORE_GENES +
    _ILC_MARKERS + _TFH_MARKERS + _TH17_MARKERS + _TH2_MARKERS + _TH1_MARKERS +
    _CYTOTOXIC_MARKERS + _MAIT_MARKERS + _GDT_MARKERS + _NK_SUBTYPE_MARKERS +
    _TREG_MARKERS + _TRM_MARKERS + _EXHAUSTION_MARKERS + _NAIVE_MARKERS +
    _EFFECTOR_MARKERS + _ACTIVATION_MARKERS + _MEMORY_MARKERS +
    ["CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2",
     "NCAM1","FCGR3A","CD56",
     "CD4","CD8A","CD8B",
     "FOXP3","IL2RA","IKZF2",
     "ITGAE","CD103","CXCR6","ZNF683",
     "TOX","TOX2","NR4A1","NR4A2","PRDM1","ID2",
     "TCF7","LEF1","KLF2","ID3","MYB",
     "EOMES","TBX21","RUNX3","RUNX1","ZEB2",
     "BCL6","ASCL2","SH2D1A","MAF",
     "GATA3","RORA","RORC",
     "MKI67","TOP2A","PCNA","STMN1","BIRC5",
     "HAVCR2","TIGIT","LAG3","PDCD1","CTLA4",
     "IFNG","TNF","IL2","IL4","IL5","IL13","IL10",
     "GNLY","NKG7","CST7","CCL4","CCL5"]
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
def normalize_label_strings(labels):
    """
    [v2.1.1 P0-5] Strip leading/trailing whitespace and collapse internal
    multiple spaces. Prevents scANVI from treating 'Regulatory T cells' and
    'Regulatory T cells ' as distinct classes during training but identical
    during prediction (phantom merging).
    """
    return pd.Series(labels).astype(str).str.strip().str.replace(r"\s+", " ", regex=True)


def apply_min_cells_filter(labels_series, min_cells, unlabeled_cat):
    """
    [v2.1.1 P0-5] Labels with fewer than min_cells cells are replaced with
    unlabeled_cat. Prevents scANVI from being given classes it cannot learn
    a reliable boundary for — which manifests as those cells being predicted
    as the nearest dominant class (false label merging).

    Returns:
        filtered_series: pd.Series with rare labels replaced
        report_df:       pd.DataFrame documenting which labels were affected
    """
    counts = labels_series.value_counts()
    rare   = [lbl for lbl in counts.index
              if counts[lbl] < min_cells and lbl != unlabeled_cat]
    report_rows = [{"label": lbl, "n_cells": int(counts[lbl]),
                    "action": f"-> {unlabeled_cat}"} for lbl in rare]
    filtered = labels_series.copy()
    if rare:
        filtered = filtered.replace(rare, unlabeled_cat)
    report_df = (pd.DataFrame(report_rows) if report_rows
                 else pd.DataFrame(columns=["label","n_cells","action"]))
    return filtered, report_df


def run_query_only_leiden(adata_filt, resolution=QUERY_LEIDEN_RESOLUTION):
    print(f"\n[Leiden] Query-only clustering (resolution={resolution})...")
    qmask  = adata_filt.obs["data_source"] == "query"
    qcells = adata_filt.obs_names[qmask]
    if len(qcells) < 10:
        adata_filt.obs["leiden_query"] = "N/A"; return
    X_q = adata_filt.obsm["X_scVI"][qmask.values]
    ad_ = sc.AnnData(X=X_q, obs=adata_filt.obs.loc[qcells].copy())
    ad_.obsm["X_scVI"] = X_q
    sc.pp.neighbors(ad_, use_rep="X_scVI", n_neighbors=30, random_state=RANDOM_SEED)
    sc.tl.leiden(ad_, resolution=resolution, random_state=RANDOM_SEED)
    lf = pd.Series("N/A", index=adata_filt.obs_names, dtype="object")
    lf.loc[qcells] = "qry_" + ad_.obs["leiden"].astype(str)
    adata_filt.obs["leiden_query"] = lf.values
    print(f"  -> {ad_.obs['leiden'].nunique()} query-only clusters")
    del ad_; gc.collect()


def _umap(adata, color, ax, title, cmap=None, vmin=None, vmax=None,
          basis="X_umap_scANVI", legend_loc="right margin"):
    """Wrapper: sc.pl.embedding with rasterized=True (v2.9 requirement)."""
    kw = dict(ax=ax, show=False, title=title, s=8,
              rasterized=True, legend_loc=legend_loc,
              legend_fontsize=6, frameon=False)
    if cmap:             kw["cmap"] = cmap
    if vmin is not None: kw["vmin"] = vmin
    if vmax is not None: kw["vmax"] = vmax
    sc.pl.embedding(adata, basis=basis, color=color, **kw)


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
print("[Step 1] Loading v2.1 preprocessed checkpoint...")
print("=" * 80)
t0 = time.time()
adata = sc.read_h5ad(INPUT_H5AD)
print(f"[OK] Loaded in {time.time()-t0:.1f}s  shape: {adata.shape}")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"  layers:      {list(adata.layers.keys())}")

# P2-8: fail fast if counts layer is missing
assert "counts" in adata.layers, (
    "CRITICAL: 'counts' layer not found in input checkpoint. "
    "Ensure v2.1 preprocessed h5ad was saved BEFORE HVG subsetting (Step 10 output)."
)

# Verify CellTypist labels are present
assert CELLTYPIST_LABEL_COL in adata.obs.columns, (
    f"CellTypist column '{CELLTYPIST_LABEL_COL}' not found. "
    f"Available: {list(adata.obs.columns)}"
)
assert CELLTYPIST_CONF_KEY in adata.obs.columns, (
    f"Confidence column '{CELLTYPIST_CONF_KEY}' not found."
)
print(f"\n[OK] CellTypist labels present: {adata.obs[CELLTYPIST_LABEL_COL].nunique()} unique types")
print("\n  All CellTypist labels in checkpoint:")
for lbl, n in adata.obs[CELLTYPIST_LABEL_COL].value_counts().items():
    print(f"    {lbl:<45} {n:>7,}  ({n/adata.n_obs*100:.1f}%)")


# %% [markdown]
# ## Cell 5 — Step 2: Classify Cells by Lineage

# %%
print("\n" + "=" * 80)
print("[Step 2] Classifying cells by CellTypist lineage...")
print("=" * 80)

# [v2.1.2 FIX-A] Normalize label strings BEFORE mapping.
# In v2.1.1 normalization only happened in Step 8 (training subset build).
# If any raw CellTypist label contains leading/trailing/double whitespace,
# map() would return NaN ("other_unknown") and trigger the P1-6 guard with a
# misleading "label not in map" error — even though the label IS in the map,
# just with a whitespace variant. Normalizing here ensures map() sees clean
# strings, and also writes the normalized values back to adata.obs so all
# downstream columns are consistent from the start.
ct_labels = normalize_label_strings(adata.obs[CELLTYPIST_LABEL_COL])
adata.obs[CELLTYPIST_LABEL_COL] = ct_labels.values   # write back normalized

# Map each label to a lineage; unmapped -> "other_unknown"
lineage_series = ct_labels.map(CELLTYPIST_LINEAGE_MAP).fillna("other_unknown")
adata.obs["celltypist_lineage"] = lineage_series.values

print("\n[INFO] Lineage distribution:")
lineage_counts     = adata.obs["celltypist_lineage"].value_counts()
lineage_counts_all = lineage_counts.copy()   # preserve for filter_qc.pdf (P2-7)
total = len(adata)
for lin, n in lineage_counts.items():
    flag = "  [KEEP]" if lin in TNK_LINEAGES else "  [REMOVE]"
    print(f"  {lin:<35} {n:>7,}  ({n/total*100:.1f}%){flag}")

# P1-6: raise on unmapped labels
unmapped_labels = ct_labels[lineage_series == "other_unknown"].unique()
if len(unmapped_labels) > 0:
    msg_lines = [
        f"\n  Unmapped CellTypist labels ({len(unmapped_labels)} unique):",
    ]
    for lbl in sorted(unmapped_labels):
        n = int((ct_labels == lbl).sum())
        msg_lines.append(f"    '{lbl}'  n={n:,}")
    msg_lines += [
        "",
        "  ACTION REQUIRED: Add each label to CELLTYPIST_LINEAGE_MAP.",
        "  T/NK labels -> assign a TNK_LINEAGES value.",
        "  Contaminants -> assign a 'contaminant_*' value.",
        "  Set KEEP_OTHER_UNKNOWN=True only as a temporary override.",
    ]
    if not KEEP_OTHER_UNKNOWN:
        raise ValueError(
            "P1-6 GUARD: Unmapped CellTypist labels found. "
            "Update CELLTYPIST_LINEAGE_MAP.\n" + "\n".join(msg_lines)
        )
    else:
        print("\n  [WARN] KEEP_OTHER_UNKNOWN=True — skipping unmapped label guard.")
        for line in msg_lines:
            print(line)


# %% [markdown]
# ## Cell 6 — Step 3: Filter Non-T/NK Cells

# %%
print("\n" + "=" * 80)
print("[Step 3] Filtering to T/NK cells only...")
print("=" * 80)

n_before  = adata.n_obs
keep_mask = adata.obs["celltypist_lineage"].isin(TNK_LINEAGES)

if KEEP_OTHER_UNKNOWN:
    keep_mask |= (adata.obs["celltypist_lineage"] == "other_unknown")
    print("  [NOTE] other_unknown cells INCLUDED (KEEP_OTHER_UNKNOWN=True)")

adata_filt = adata[keep_mask].copy()
n_after    = adata_filt.n_obs
n_removed  = n_before - n_after

print(f"  Before:  {n_before:,} cells")
print(f"  Removed: {n_removed:,} non-T/NK cells ({n_removed/n_before*100:.1f}%)")
print(f"  After:   {n_after:,} cells")

print("\n  Remaining lineage breakdown:")
for lin, n in adata_filt.obs["celltypist_lineage"].value_counts().items():
    print(f"    {lin:<35} {n:>7,}  ({n/n_after*100:.1f}%)")

print("\n  Remaining data_source breakdown:")
for src, n in adata_filt.obs["data_source"].value_counts().items():
    print(f"    {src:<25} {n:>7,}")

del adata; gc.collect()


# %% [markdown]
# ## Cell 7 — Step 4: Drop Small Batches (Post-Filter)

# %%
print("\n" + "=" * 80)
print(f"[Step 4] Dropping batches with < {MIN_BATCH_CELLS} cells post-filter...")
print("=" * 80)

batch_counts  = adata_filt.obs[BATCH_KEY].value_counts()
small_batches = batch_counts[batch_counts < MIN_BATCH_CELLS].index.tolist()

if small_batches:
    print(f"  Dropping {len(small_batches)} small batches:")
    for b in small_batches:
        print(f"    {b}  (n={batch_counts[b]})")
    keep_batch = ~adata_filt.obs[BATCH_KEY].isin(small_batches)
    adata_filt  = adata_filt[keep_batch].copy()
    print(f"  Cells after batch drop: {adata_filt.n_obs:,}")
else:
    print(f"  All {batch_counts.shape[0]} batches have >= {MIN_BATCH_CELLS} cells — no drop needed")

# Re-encode as category (categories may have changed after filtering)
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
print("[Step 5] Rebuilding CellTypist filtered label column (always fresh, P0-2)...")
print("=" * 80)
# P0-2: Always recompute from current CELLTYPIST_CONF_THRESHOLD.
# Never reuse the column from the checkpoint — threshold may have changed.
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
print("[Step 6] Selecting HVGs on filtered T/NK data...")
print("=" * 80)

# P2-9: Clear ALL old HVG auxiliary columns so re-selection starts clean.
_hvg_aux_cols = ["highly_variable","highly_variable_rank",
                 "means","variances","variances_norm",
                 "dispersions","dispersions_norm"]
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
# ## Cell 11 — Step 8: Build Training Subset (HVG) + Label Fixes

# %%
print("\n" + "=" * 80)
print("[Step 8] Building training subset (HVG only) + applying label fixes...")
print("=" * 80)

hvg_mask = adata_filt.var["highly_variable"].values
X_hvg    = adata_filt.layers["counts"][:, hvg_mask]
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

# ── [v2.1.1 P0-5] Step A: Normalize label strings ────────────────────────────
raw_labels     = normalize_label_strings(adata_train.obs[label_source_col])
n_ws_fixed     = int((raw_labels != adata_train.obs[label_source_col].astype(str)).sum())
if n_ws_fixed > 0:
    print(f"  [P0-5] Normalized whitespace in {n_ws_fixed:,} label strings")

# ── [v2.1.1 P0-5] Step B: Filter rare labels -> Unknown ──────────────────────
filtered_labels, rare_report = apply_min_cells_filter(
    raw_labels, MIN_CELLS_PER_LABEL, UNLABELED_CATEGORY
)

if len(rare_report) > 0:
    print(f"\n  [P0-5] {len(rare_report)} labels below "
          f"MIN_CELLS_PER_LABEL={MIN_CELLS_PER_LABEL}:")
    print(f"  {'Label':<45} {'N cells':>8}  {'Action'}")
    print(f"  {'-'*45} {'-'*8}  {'-'*20}")
    for _, row in rare_report.iterrows():
        print(f"  {row['label']:<45} {row['n_cells']:>8,}  {row['action']}")
    rare_report.to_csv(
        output_dir / f"{OUTPUT_PREFIX}_rare_labels_filtered.csv", index=False
    )
    print(f"\n  [OK] Rare label report: {OUTPUT_PREFIX}_rare_labels_filtered.csv")
else:
    print(f"  [OK] No labels below MIN_CELLS_PER_LABEL={MIN_CELLS_PER_LABEL}")

# ── Assign final scanvi_labels ────────────────────────────────────────────────
adata_train.obs["scanvi_labels"] = filtered_labels.values
adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")
if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
    adata_train.obs["scanvi_labels"] = (
        adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])
    )

# P0-1: Write scanvi_labels back to adata_filt via reindex.
# Without this, adata_filt.write_h5ad() would lack the training label column.
adata_filt.obs["scanvi_labels"] = (
    adata_train.obs["scanvi_labels"]
    .reindex(adata_filt.obs_names)
    .astype("category")
)
print(f"  -> scanvi_labels written back to adata_filt "
      f"({adata_filt.obs['scanvi_labels'].nunique()} categories incl. Unknown)")

# ── [v2.1.1 P0-3] Validate required covariates ───────────────────────────────
required_cont = ["pct_counts_mt","stress_score","S_score","G2M_score"]
missing_cont  = [c for c in required_cont if c not in adata_train.obs.columns]
if missing_cont:
    raise ValueError(
        f"P0-3 GUARD: Required covariates missing from checkpoint: {missing_cont}. "
        "Re-run v2.1 pipeline to ensure prepare_covariates() was called before "
        "saving the preprocessed checkpoint."
    )
print(f"  [OK] All covariates present: {required_cont}")

# ── Pre-training label stats (saved for post-training audit) ─────────────────
print("\n  -> [PRE-TRAINING] scANVI label distribution:")
label_counts_pretrain = adata_train.obs["scanvi_labels"].value_counts()
n_unknown_train = int((adata_train.obs["scanvi_labels"] == UNLABELED_CATEGORY).sum())
n_labeled_types = int((label_counts_pretrain.index != UNLABELED_CATEGORY).sum())
print(f"     Distinct labeled types: {n_labeled_types}")
print(f"     Unknown cells:          {n_unknown_train:,} "
      f"({100*n_unknown_train/adata_train.n_obs:.1f}%)")
print()
for lbl, cnt in label_counts_pretrain.items():
    flag = "  [Unknown]" if lbl == UNLABELED_CATEGORY else ""
    print(f"     {lbl:<45} {cnt:>7,}{flag}")

_pretrain_stats = (label_counts_pretrain.rename("n_cells_pretrain")
                   .reset_index().rename(columns={"index":"label"}))
gc.collect()


# %% [markdown]
# ## Cell 12 — Step 9: scVI Training

# %%
print("\n" + "=" * 80)
print(f"[Step 9] Training scVI (n_latent={SCVI_N_LATENT})...")
print("=" * 80)

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
print("[Step 10] Training scANVI (CellTypist labels for all filtered T/NK cells)...")
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
# ## Cell 14 — Step 11: Export Latent & Predictions + Post-Training Audit

# %%
print("\n" + "=" * 80)
print("[Step 11] Exporting latent representations and predictions...")
print("=" * 80)

# scANVI latent
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

# Novelty score
if not np.isfinite(proba_al.values).all():
    raise ValueError("Non-finite values in scANVI probability matrix!")
ents = entropy(proba_al.values + 1e-10, axis=1)
adata_filt.obs["scanvi_entropy"] = ents
emin, emax = ents.min(), ents.max()
adata_filt.obs["novelty_score"] = (ents - emin) / (emax - emin) if emax > emin else 0.0
qmask = adata_filt.obs["data_source"] == "query"
adata_filt.obs["is_potentially_novel"] = (adata_filt.obs["novelty_score"] > 0.7) & qmask
print(f"  High novelty query cells: {adata_filt.obs['is_potentially_novel'].sum()}")

# ── Post-training label audit ─────────────────────────────────────────────────
print("\n" + "-" * 60)
print("[AUDIT] Post-training label comparison:")
print("-" * 60)

pred_counts = adata_filt.obs["scanvi_pred"].value_counts().rename("n_cells_predicted")
_pred_reset = pred_counts.reset_index()
_pred_reset.columns = ["label","n_cells_predicted"]

_audit = _pretrain_stats.merge(_pred_reset, on="label", how="outer").fillna(0)
_audit["n_cells_pretrain"]   = _audit["n_cells_pretrain"].astype(int)
_audit["n_cells_predicted"]  = _audit["n_cells_predicted"].astype(int)

# Labels in training but absent in predictions -> potential merging
_missing = _audit[
    (_audit["n_cells_pretrain"] > 0) &
    (_audit["label"] != UNLABELED_CATEGORY) &
    (_audit["n_cells_predicted"] == 0)
]
if len(_missing) > 0:
    print(f"\n  [WARNING] {len(_missing)} training labels ABSENT from predictions:")
    print("  These labels may have been absorbed by dominant classes despite the")
    print("  MIN_CELLS_PER_LABEL guard — consider whether they are transcriptomically")
    print("  distinct enough to warrant a separate class.")
    for _, row in _missing.iterrows():
        print(f"    '{row['label']}': trained on {int(row['n_cells_pretrain']):,} cells, 0 predicted")
else:
    print("  [OK] All training labels appear in predictions — no merging detected.")

# Labels in predictions but absent in training (should not happen)
_new = _audit[
    (_audit["n_cells_pretrain"] == 0) &
    (_audit["label"] != UNLABELED_CATEGORY) &
    (_audit["n_cells_predicted"] > 0)
]
if len(_new) > 0:
    print(f"\n  [INFO] {len(_new)} prediction labels not seen in training (unexpected):")
    for _, row in _new.iterrows():
        print(f"    '{row['label']}': {int(row['n_cells_predicted']):,} predicted")

print(f"\n  Training labels:   {int((_audit['n_cells_pretrain'] > 0).sum())}")
print(f"  Prediction labels: {int((_audit['n_cells_predicted'] > 0).sum())}")
_audit.to_csv(output_dir / f"{OUTPUT_PREFIX}_label_audit.csv", index=False)
print(f"  [OK] Audit saved: {OUTPUT_PREFIX}_label_audit.csv")
print("-" * 60)


# %% [markdown]
# ## Cell 15 — Step 12: Attach .raw (Full Gene Matrix)
#
# [v2.1.2 FIX-C] .raw is attached BEFORE CD4/CD8 scoring (Step 13).
# In v2.1.1 .raw was attached after scoring, so score_genes() silently
# fell back to adata_filt.X whose normalization state was uncontrolled
# (depends on checkpoint preprocessing). Attaching .raw first ensures
# scoring always uses log-normalized full-gene expression via use_raw=True.

# %%
print("\n" + "=" * 80)
print("[Step 12] Attaching .raw (full gene matrix — must precede scoring)...")
print("=" * 80)
from anndata import AnnData as _AnnData
adata_filt.raw = _AnnData(X=full_counts, obs=adata_filt.obs.copy(), var=raw_var)
print(f"  .raw: {adata_filt.raw.n_vars} genes")


# %% [markdown]
# ## Cell 16 — Step 13: CD4/CD8 Scoring QC
#
# [v2.1.2 FIX-B] NK exclusion mask now uses the CellTypist direct label
# directly, covering BOTH canonical NK and cycling NK cells.
# In v2.1.1 the mask was: celltypist_lineage == "NK", which missed
# "Cycling NK cells" (mapped to "cycling_T" then, now "cycling_NK").
# Using a direct label set is more robust than relying on the coarse
# lineage tag — any future NK-derived label added to the map is only
# excluded here if it is explicitly listed in NK_DIRECT_LABELS.
#
# [v2.1.2 FIX-C] score_genes() now uses use_raw=True because .raw is
# already attached above (log-normalized full-gene expression guaranteed).

# %%
print("\n" + "=" * 80)
print("[Step 13] CD4/CD8 QC scoring (NK + Cycling NK excluded, FIX-B/FIX-C)...")
print("=" * 80)

# All CellTypist labels that should be treated as NK for scoring purposes.
# Add any future NK-derived labels here explicitly.
NK_DIRECT_LABELS = {
    "NK cells",
    "CD16+ NK cells",
    "CD16- NK cells",
    "Cycling NK cells",
}

# [FIX-C] .raw is now attached — use_raw=True is safe
_use_raw = True   # guaranteed: .raw attached in Step 12
_vnames  = adata_filt.raw.var_names

def _score(genes):
    present = [g for g in genes if g in _vnames]
    if len(present) == 0:
        return None
    sc.tl.score_genes(adata_filt, gene_list=present,
                      score_name="_tmp_score_", use_raw=_use_raw)
    s = adata_filt.obs["_tmp_score_"].values.copy()
    adata_filt.obs.drop(columns=["_tmp_score_"], inplace=True)
    return s

cd4_raw = _score(CD4_SCORE_GENES)
cd8_raw = _score(CD8_SCORE_GENES)
nk_raw  = _score(NK_SCORE_GENES)

adata_filt.obs["cd4_score"] = cd4_raw if cd4_raw is not None else 0.0
adata_filt.obs["cd8_score"] = cd8_raw if cd8_raw is not None else 0.0
adata_filt.obs["nk_score"]  = nk_raw  if nk_raw  is not None else 0.0

# [FIX-B] NK mask: match on direct CellTypist label, not coarse lineage.
# normalize_label_strings() applied for consistency with Step 2 normalization.
_nk_mask = (
    normalize_label_strings(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY])
    .isin(NK_DIRECT_LABELS)
    .values
)
_cd4_pos = adata_filt.obs["cd4_score"].values > CD4_SCORE_THRESH
_cd8_pos = adata_filt.obs["cd8_score"].values > CD8_SCORE_THRESH

# cd4_cd8_by_score categories:
#   NK           -> NK or cycling NK (bypasses CD4/CD8 scoring entirely)
#   CD4_single   -> T cell: CD4+ score only
#   CD8_single   -> T cell: CD8+ score only
#   DP           -> T cell: both scores+ (ambiguous, check for contamination)
#   DN           -> T cell: neither score+ (resting / low-expression)
adata_filt.obs["cd4_cd8_by_score"] = np.select(
    [_nk_mask,
     ~_nk_mask & _cd4_pos & ~_cd8_pos,
     ~_nk_mask & ~_cd4_pos & _cd8_pos,
     ~_nk_mask & _cd4_pos & _cd8_pos,
     ~_nk_mask & ~_cd4_pos & ~_cd8_pos],
    ["NK", "CD4_single", "CD8_single", "DP", "DN"],
    default="other"
).astype(str)
adata_filt.obs["cd4_cd8_by_score"] = adata_filt.obs["cd4_cd8_by_score"].astype("category")

print("[INFO] CD4/CD8 score classification (NK + Cycling NK excluded):")
print(adata_filt.obs["cd4_cd8_by_score"].value_counts().to_string())
print(f"  NK mask total (direct label): {_nk_mask.sum():,} cells")


# %% [markdown]
# ## Cell 17 — Step 14: Compute UMAPs

# %%
print("\n" + "=" * 80)
print("[Step 14] Computing UMAPs...")
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
# ## Cell 18 — Step 15: Run Log

# %%
print("\n" + "=" * 80)
print("[Step 15] Writing run log...")
print("=" * 80)
log_ts   = datetime.now().isoformat()
out_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"

pipeline_log = {
    "version":              "2.1.2",  "timestamp": log_ts,
    "input_h5ad":           INPUT_H5AD,
    "output_dir":           str(output_dir), "output_prefix": OUTPUT_PREFIX,
    "planned_output":       str(out_h5ad),
    "elapsed_min":          round((time.time() - PIPELINE_START) / 60, 2),
    "n_obs_before_filter":  n_before,
    "n_obs_after_filter":   adata_filt.n_obs,
    "n_removed":            n_removed,
    "pct_removed":          round(n_removed / n_before * 100, 2),
    "n_vars":               int(adata_filt.n_vars),
    "n_hvg":                n_hvg_final,
    "hvg_method":           hvg_method,
    "gpu_available":        bool(gpu_available),
    "scvi_n_latent":        SCVI_N_LATENT,
    "min_cells_per_label":  MIN_CELLS_PER_LABEL,
    "n_rare_labels_filtered": len(rare_report),
    "architecture":         "CellTypist-filtered T/NK-only retrain (v2.1.2)",
    "scanvi_label_source":  label_source_col,
    "min_batch_cells":      MIN_BATCH_CELLS,
    "small_batches_dropped": small_batches if small_batches else [],
}
anndata_structure = {
    "timestamp": log_ts, "shape": [int(adata_filt.n_obs), int(adata_filt.n_vars)],
    "raw": {"present": True, "n_vars": adata_filt.raw.n_vars},
    "obs_columns": list(adata_filt.obs.columns),
    "obsm_keys":   sorted(adata_filt.obsm.keys()),
    # [v2.1.2 FIX-D] uns_keys is populated AFTER pipeline_log is written,
    # so the JSON accurately reflects the actual uns contents including
    # "pipeline_log" and "anndata_structure" themselves.
    "uns_keys":    "PLACEHOLDER",
}
# Write pipeline_log first, then finalize uns_keys
adata_filt.uns["pipeline_log"]      = pipeline_log
adata_filt.uns["anndata_structure"] = {}                    # reserve the key
anndata_structure["uns_keys"] = sorted(str(k) for k in adata_filt.uns.keys())
adata_filt.uns["anndata_structure"] = anndata_structure     # now write the real dict
(output_dir / f"{OUTPUT_PREFIX}_run_log.txt").write_text(
    "\n".join(f"{k}: {v}" for k, v in pipeline_log.items()), encoding="utf-8"
)
with open(output_dir / f"{OUTPUT_PREFIX}_anndata_structure.json", "w") as f:
    json.dump(anndata_structure, f, indent=2)
print("  -> Log written")


# %% [markdown]
# ## Cell 19 — Step 16: Save

# %%
print("\n" + "=" * 80)
print("[Step 16] Saving results...")
print("=" * 80)

_cat_cols = [
    "data_source", "cell_type_fine_ref", "celltypist_lineage",
    "scanvi_labels", "scanvi_pred",
    "cd4_cd8_by_score", "leiden_query",
    CELLTYPIST_DIRECT_LABEL_KEY, CELLTYPIST_DIRECT_FILT_KEY,
]
for col in _cat_cols:
    for ad in [adata_filt, adata_train]:
        if col in ad.obs.columns:
            ad.obs[col] = ad.obs[col].astype("category")

scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
scvi_model.save(output_dir   / f"{OUTPUT_PREFIX}_scvi_model",   overwrite=True)

config = {
    "version":   "2.1.2", "timestamp": datetime.now().isoformat(),
    "input":     INPUT_H5AD,
    "n_hvg":     n_hvg_final, "hvg_method": hvg_method,
    "scvi_n_latent": SCVI_N_LATENT,
    "min_cells_per_label": MIN_CELLS_PER_LABEL,
    "architecture": {
        "description":         "CellTypist-filtered T/NK-only retrain — contaminants removed",
        "filter_policy":       "keep TNK_LINEAGES, remove contaminants",
        "tnk_lineages":        sorted(TNK_LINEAGES),
        "keep_other_unknown":  KEEP_OTHER_UNKNOWN,
        "n_removed":           n_removed,
        "pct_removed":         round(n_removed / n_before * 100, 2),
        "scanvi_labels":       label_source_col,
        "fixes": [
            "NK cells excluded from CD4/CD8 score classification (P0-4)",
            "Label string normalization before training (P0-5)",
            f"Labels < {MIN_CELLS_PER_LABEL} cells -> Unknown before training (P0-5)",
            "Post-training label audit -> label_audit.csv",
            "CELLTYPIST_DIRECT_FILT_KEY rebuilt from scratch (P0-2)",
            "scanvi_labels written back to adata_filt via reindex (P0-1)",
            "[v2.1.2] Label normalization moved to Step 2 before lineage map (FIX-A)",
            "[v2.1.2] Cycling NK mapped to cycling_NK; NK mask uses direct labels (FIX-B)",
            "[v2.1.2] .raw attached before CD4/CD8 scoring; use_raw=True enforced (FIX-C)",
            "[v2.1.2] uns_keys recorded after pipeline_log written (FIX-D)",
        ],
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

_stats = (
    adata_filt.obs["scanvi_pred"].value_counts()
    .rename_axis("Cell_Type").reset_index(name="Count")
)
_stats["Percentage"] = 100 * _stats["Count"] / _stats["Count"].sum()
_conf = adata_filt.obs.groupby("scanvi_pred")["scanvi_confidence"].agg(["mean","std"])
_stats = _stats.merge(_conf, left_on="Cell_Type", right_index=True, how="left")
_stats.to_csv(output_dir / f"{OUTPUT_PREFIX}_scanvi_statistics.csv", index=False)

adata_filt.obs[[
    "data_source", "celltypist_lineage", "scanvi_pred",
    CELLTYPIST_DIRECT_LABEL_KEY, CELLTYPIST_DIRECT_FILT_KEY,
    "scanvi_labels", "scanvi_confidence", "novelty_score",
]].to_csv(output_dir / f"{OUTPUT_PREFIX}_annotations.csv")

adata_filt.write_h5ad(out_h5ad, compression="gzip")
print(f"  -> {out_h5ad}")
adata_train.write_h5ad(output_dir / f"{OUTPUT_PREFIX}_train_HVG.h5ad", compression="gzip")
print("  -> train HVG saved")


# %% [markdown]
# ## Cell 20 — Step 17: Visualization (v2.9 rasterized)

# %%
print("\n" + "=" * 80)
print("[Step 17] Visualization...")
elapsed = (time.time() - PIPELINE_START) / 60
print(f"  Elapsed so far: {elapsed:.1f} min")
print("=" * 80)

# ── Panel 1: Overview 4×4 ────────────────────────────────────────────────────
fig = plt.figure(figsize=(24, 20))
gs  = fig.add_gridspec(4, 4, hspace=0.3, wspace=0.3)

row0_specs = [
    ("data_source",              "Data Source"),
    ("celltypist_lineage",       "CellTypist Lineage"),
    (CELLTYPIST_DIRECT_LABEL_KEY,"CellTypist Direct"),
    (CELLTYPIST_DIRECT_FILT_KEY, "CellTypist Filtered"),
]
row1_specs = [
    ("scanvi_pred",       "scANVI Predictions",    None,      None),
    ("scanvi_confidence", "scANVI Confidence",     "viridis", (0, 1)),
    ("novelty_score",     "Novelty Score",         "hot",     (0, 1)),
    ("cd4_cd8_by_score",  "CD4/CD8 by Score (QC)", None,      None),
]
for col_i, (key, title) in enumerate(row0_specs):
    if key in adata_filt.obs.columns:
        ax = fig.add_subplot(gs[0, col_i])
        _umap(adata_filt, key, ax, title, legend_loc="on data")

for col_i, (key, title, cmap, vlim) in enumerate(row1_specs):
    if key in adata_filt.obs.columns:
        ax = fig.add_subplot(gs[1, col_i])
        _umap(adata_filt, key, ax, title, cmap=cmap,
              vmin=vlim[0] if vlim else None, vmax=vlim[1] if vlim else None,
              legend_loc="on data" if not cmap else "right margin")

marker_panels = [
    ("CD3D",  "CD3D (pan-T)"),        ("CD4",   "CD4"),
    ("CD8A",  "CD8A"),                 ("GNLY",  "GNLY (NK)"),
    ("NCAM1", "NCAM1 (NK definitive)"),("FOXP3", "FOXP3 (Treg)"),
    ("HAVCR2","HAVCR2 (Exhaustion)"), ("MKI67", "MKI67 (Prolif)"),
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

# ── Panel 2: CellTypist vs scANVI ─────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(21, 7))
for ax, (key, title) in zip(axes2, [
    (CELLTYPIST_DIRECT_LABEL_KEY, "CellTypist Direct"),
    (CELLTYPIST_DIRECT_FILT_KEY,  "CellTypist Filtered"),
    ("scanvi_pred",               "scANVI (retrained, T/NK-only)"),
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

# ── Panel 3: scVI vs scANVI UMAP comparison ───────────────────────────────────
fig3, axes3 = plt.subplots(2, 3, figsize=(18, 12))
for row_i, (basis, label) in enumerate([("X_umap_scVI","scVI"),
                                         ("X_umap_scANVI","scANVI")]):
    for col_i, (color, loc) in enumerate([
        ("data_source",               "right margin"),
        (CELLTYPIST_DIRECT_LABEL_KEY,  "on data"),
        ("scanvi_pred",                "on data"),
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

# ── Panel 4: NK vs T cell marker separation ───────────────────────────────────
nk_t_genes = ["NCAM1","FCGR3A","KIR2DL1","KIR3DL1",   # NK definitive
               "CD3D","CD3E","TRAC",                    # T definitive
               "CD4","CD8A",                            # lineage split
               "GNLY","NKG7"]                           # shared cytotoxic
_nkt_present = [g for g in nk_t_genes if g in (
    adata_filt.raw.var_names if adata_filt.raw is not None else adata_filt.var_names)]
if _nkt_present:
    ncols = (len(_nkt_present) + 1) // 2
    fig4, axes4 = plt.subplots(2, ncols, figsize=(4 * ncols, 8))
    axes4_flat = np.array(axes4).flatten()
    for ax, gene in zip(axes4_flat, _nkt_present):
        sc.pl.embedding(
            adata_filt, basis="X_umap_scANVI", color=gene,
            ax=ax, show=False, title=gene, cmap="Reds", s=8,
            rasterized=True, frameon=False, use_raw=(adata_filt.raw is not None),
        )
    for ax in axes4_flat[len(_nkt_present):]:
        ax.set_visible(False)
    plt.suptitle("NK vs T cell marker separation", fontsize=14, y=1.01)
    plt.tight_layout()
    fig4.savefig(output_dir / f"{OUTPUT_PREFIX}_nk_t_separation.pdf",
                 dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
    plt.close(fig4)
    print(f"  -> {OUTPUT_PREFIX}_nk_t_separation.pdf")

# ── Panel 5: Filter QC — before vs after (P2-7) ──────────────────────────────
# Left:  ALL lineages from pre-filter data (includes contaminants)
# Right: retained T/NK lineages only + filtering summary
fig5, axes5 = plt.subplots(1, 2, figsize=(16, 7))

_all_vals   = lineage_counts_all.values
_all_labels = lineage_counts_all.index.tolist()
_colors     = ["#d62728" if "contaminant" in l or l == "other_unknown"
               else "#1f77b4" for l in _all_labels]
axes5[0].pie(_all_vals, labels=_all_labels, autopct="%1.1f%%",
             startangle=140, colors=_colors, textprops={"fontsize": 7})
axes5[0].set_title(
    f"All Lineages — Before Filter\n(n={n_before:,} cells)\n"
    "Blue = T/NK kept  |  Red = contaminants removed",
    fontsize=9
)

_ret_lin = adata_filt.obs["celltypist_lineage"].value_counts()
axes5[1].pie(_ret_lin.values, labels=_ret_lin.index, autopct="%1.1f%%",
             startangle=140, textprops={"fontsize": 7})
axes5[1].set_title(
    f"Retained T/NK Lineages — After Filter\n"
    f"(n={n_after:,} cells  |  removed {n_removed:,} = {n_removed/n_before*100:.1f}%)\n"
    f"Batches dropped (< {MIN_BATCH_CELLS} cells): {len(small_batches)}",
    fontsize=9
)
plt.tight_layout()
fig5.savefig(output_dir / f"{OUTPUT_PREFIX}_filter_qc.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig5)
print(f"  -> {OUTPUT_PREFIX}_filter_qc.pdf")

# ── Dotplot: T/NK marker overview ─────────────────────────────────────────────
dotplot_genes = [
    "CD3D","CD3E","TRAC",
    "CD4","IL7R","MAL",
    "CD8A","CD8B",
    "NCAM1","FCGR3A","KLRB1","KLRC1","KLRD1","KLRF1",
    "KIR2DL1","KIR2DL3","KIR3DL1",
    "GNLY","NKG7","GZMB","GZMK","GZMA","GZMH","PRF1",
    "FOXP3","IL2RA","IKZF2","CTLA4","TNFRSF18",
    "CXCR5","BCL6","ICOS","SH2D1A",
    "GATA3","IL4","IL13","CRTH2","CCR4",
    "TBX21","IFNG","CXCR3",
    "RORC","IL17A","CCR6",
    "CCR7","SELL","TCF7","LEF1","KLF2",
    "ITGAE","CD69","CXCR6","ZNF683",
    "PDCD1","HAVCR2","LAG3","TIGIT","TOX","NR4A1",
    "CD38","HLA-DRA","TNFRSF4","TNFRSF9",
    "SLC4A10","TRDC","MKI67","TOP2A","PCNA",
]
_use_raw_dp  = adata_filt.raw is not None
_dp_genes    = [g for g in dotplot_genes if g in (
    adata_filt.raw.var_names if _use_raw_dp else adata_filt.var_names)]
_groupby_col = "scanvi_pred"
if _dp_genes:
    try:
        dp = sc.pl.dotplot(
            adata_filt, var_names=_dp_genes, groupby=_groupby_col,
            use_raw=_use_raw_dp, standard_scale="var",
            cmap="Reds", show=False, return_fig=True,
            figsize=(max(14, len(_dp_genes)*0.55),
                     max(6, adata_filt.obs[_groupby_col].nunique()*0.4)),
            dendrogram=False,
        )
        dp.savefig(str(output_dir / f"{OUTPUT_PREFIX}_dotplot_markers.pdf"),
                   dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close("all")
        print(f"  -> {OUTPUT_PREFIX}_dotplot_markers.pdf")
    except Exception as e:
        print(f"  [WARN] Dotplot failed: {e}")

# ── Final summary ─────────────────────────────────────────────────────────────
elapsed = (time.time() - PIPELINE_START) / 60
print("\n" + "=" * 80)
print("T/NK PIPELINE v2.1.2 COMPLETE")
print("=" * 80)
print(f"  Input cells:   {n_before:,}")
print(f"  Removed:       {n_removed:,}  ({n_removed/n_before*100:.1f}%)")
print(f"  Retained:      {n_after:,}")
print(f"  Batches used:  {adata_filt.obs[BATCH_KEY].nunique()}")
print(f"  Elapsed:       {elapsed:.1f} min")
print(f"\n  scVI n_latent:      {SCVI_N_LATENT}")
print(f"  scANVI labels from: {label_source_col}")
print(f"  MIN_CELLS_PER_LABEL: {MIN_CELLS_PER_LABEL}  "
      f"({len(rare_report)} rare labels filtered -> Unknown)")
# P0-1: scanvi_labels now written back to adata_filt — no fallback needed
n_unknown_final = int(
    (adata_filt.obs["scanvi_labels"].astype(str) == UNLABELED_CATEGORY).sum()
)
print(f"  Unknown (low-conf + rare): {n_unknown_final:,}")
print("\nTop scANVI predictions:")
print(adata_filt.obs["scanvi_pred"].value_counts().head(12))
print("\nTop CellTypist labels (post-filter):")
print(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].value_counts().head(12))
print("\nOutput files:")
for f in sorted(output_dir.glob(f"{OUTPUT_PREFIX}*")):
    print(f"  {f.name}")
print("=" * 80)
