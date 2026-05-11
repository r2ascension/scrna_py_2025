#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
T/NK Cell Merged Pipeline v2.4 — CellTypist-Filtered Retrain
=============================================================
Purpose:
  Downstream of v2.1. Loads the v2.1 preprocessed checkpoint (already has
  CellTypist labels, HVG flags, covariates, no .raw yet), removes non-T/NK
  contaminants based on CellTypist lineage classification, then retrains
  scVI + scANVI on the clean T/NK-only subset.

Input:
  tcell_merged_v2_adata_preprocessed.h5ad  (from v2.1 Cell 13 / Step 10)
  — already contains CellTypist labels, HVG flags, covariates, no .raw yet

Known CellTypist labels (from v2.0/v2.1 run, Immune_All_Low.pkl):
  Trm cytotoxic T cells (23.7%), Tem/Trm cytotoxic T cells (19.3%),
  Regulatory T cells (13.4%), Tem/Effector helper T cells (9.7%),
  CD16+ NK cells (9.7%), Tem/Temra cytotoxic T cells (6.8%),
  CD16- NK cells (4.3%), Tcm/Naive helper T cells (3.7%),
  NK cells (2.3%), Type 17 helper T cells (1.9%), ...
  Endothelial cells*, Epithelial cells*, Mast cells*  (* = contaminants)

v2.9 compliance:
  - rasterized=True on all UMAP scatter plots, dpi=300

v2.1.1 HOTFIX:
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

v2.1.2:
  FIX-A: label normalization moved to Step 2 (before CELLTYPIST_LINEAGE_MAP)
  FIX-B: Cycling NK cells mapped to cycling_NK; NK mask uses direct labels
  FIX-C: .raw attached BEFORE CD4/CD8 scoring; use_raw=True hardcoded
  FIX-D: uns_keys recorded AFTER pipeline_log/anndata_structure are written
  FIX-E: removed unused variables

v2.4 — additional fixes from production review:
  P0-A: Step 13 `or 0.0` -> explicit `is None` check (numpy array truth value
        ambiguity, raises ValueError at runtime)
  P0-B: .raw stores log1p expression, NOT raw counts; score_genes / feature
        plots use log1p (counts-based scoring is depth-confounded and wrong)
  P1-C: LABEL_MERGE_MAP moved BEFORE min-cells filter so merging can rescue
        biologically real but numerically rare labels
  P1-D: apply_max_cells_cap_stratified replaces naive global cap — stratifies
        by data_source to preserve ref/query representation per label
  P1-E: UMAP operator save annotated as "approximate projector"; comment
        clarifies it does NOT reproduce sc.tl.umap exact coordinates
  P2-F: Added obs_names/var_names/BATCH_KEY/TISSUE_KEY uniqueness and notna
        guards in Step 1 (fail-fast before any expensive ops)
  P2-G: RARE_LABEL_WHITELIST skips min-cells filter for biologically important
        rare cell types (MAIT, gdT, ILC variants)
  FIX-F: [P0 crash] _pretrain_stats pandas 3.0 compat — use .columns=[...]
          instead of rename(columns={"index":"label"}) which is a no-op in
          pandas 3.0 (index col is now named after the original Series,
          not "index"). KeyError: 'label' in audit merge.
  FIX-G: [label merging] Add MAX_CELLS_PER_LABEL=15000 to cap dominant
          classes before scANVI training. In the observed run, the top 3
          CD8 TRM-like subtypes totaled ~50k labeled cells (39% of labeled
          data), making scANVI collapse all CD8 subtypes into Trm cytotoxic.
          Capping converts excess cells of dominant classes to Unknown
          (they remain as unlabeled semi-supervised data, not discarded).
  FIX-H: [label merging] Add optional LABEL_MERGE_MAP for biologically
          motivated coarse merging of transcriptomically similar subtypes.
          Defaults to {} (no merging). Populated in Cell 1 Configuration.
  FIX-I: [label merging] Increase SCANVI_PATIENCE from 20 -> 50 to allow
          scANVI more epochs to resolve fine-grained subtype boundaries.
  FIX-J: Add imbalance diagnostics (Top1/Total ratio, Max/Min ratio) before
          training so the user can see whether the dataset is well-balanced.

Author: r2end
Date: 2026-03-18
Version: 2.4
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
print("T/NK Cell Merged Pipeline v2.4  (CellTypist-Filtered Retrain)")
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

INPUT_H5AD = "/home/h2048/data/py/20260310/tcell_only_merged_pipeline_v2/tcell_merged_v2_adata_preprocessed.h5ad"

OUTPUT_DIR    = "/home/h2048/data/py/20260318/tcell_only_merged_pipeline_v2_3_filtered"
OUTPUT_PREFIX = "tcell_merged_v2_3"

BATCH_KEY  = "sample"
TISSUE_KEY = "tissue"

CELLTYPIST_LABEL_COL      = "celltypist_label_direct"
CELLTYPIST_CONF_KEY       = "celltypist_confidence"
CELLTYPIST_CONF_THRESHOLD = 0.5

MIN_BATCH_CELLS = 10

# ==============================================================================
# LINEAGE FILTER MAP
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
    "Cycling NK cells":               "cycling_NK",   # FIX-B: own category
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

TNK_LINEAGES = {"CD4", "CD8", "NK", "cycling_NK", "gdT", "MAIT", "ILC", "cycling_T"}

KEEP_OTHER_UNKNOWN = False

# ==============================================================================
# scVI / scANVI parameters
# ==============================================================================
N_HVG                = 4000
FORCE_MARKERS_IN_HVG = True
SCVI_N_LATENT        = 150
SCVI_N_LAYERS        = 2
SCVI_N_HIDDEN        = 128
SCVI_DROPOUT         = 0.1
MAX_EPOCHS_SCVI      = 400
MAX_EPOCHS_SCANVI    = 200
UNLABELED_CATEGORY   = "Unknown"
BATCH_SIZE           = 256
LEARNING_RATE        = 1e-3
WEIGHT_DECAY         = 0.0

SCANVI_LABEL_SOURCE         = "filtered"
CELLTYPIST_DIRECT_LABEL_KEY = "celltypist_label_direct"
CELLTYPIST_DIRECT_FILT_KEY  = "celltypist_label_direct_filt"

MIN_CELLS_PER_LABEL  = 50    # labels below this -> Unknown before training

# [v2.4 P2-G] Whitelist: these labels are NEVER filtered by min-cells, even if
# their count falls below MIN_CELLS_PER_LABEL. Use for biologically important
# rare cell types where even 10-40 cells carry real signal.
RARE_LABEL_WHITELIST = {
    "MAIT cells",
    "gamma-delta T cells",
    "CRTAM+ gamma-delta T cells",
    "ILC3",
}

# [v2.3 FIX-G] Cap dominant classes to prevent scANVI label collapse.
# In v2.4 this uses stratified capping by data_source (see apply_max_cells_cap_stratified).
# Set to None to disable capping.
MAX_CELLS_PER_LABEL  = 15000

# [v2.3 FIX-H] Optional coarse label merging for transcriptomically similar
# subtypes that scANVI cannot reliably distinguish.
# Populated based on dotplot / biological knowledge — leave {} to skip.
# Applied AFTER MIN_CELLS_PER_LABEL filter and BEFORE MAX_CELLS_PER_LABEL cap.
#
# Example for this dataset:
# LABEL_MERGE_MAP = {
#     "Tem/Trm cytotoxic T cells":  "CD8_Trm",
#     "Trm cytotoxic T cells":      "CD8_Trm",
#     "Tem/Temra cytotoxic T cells":"CD8_Temra",
# }
LABEL_MERGE_MAP = {}

# [v2.3 FIX-I] Increase patience so scANVI has enough epochs to resolve
# fine-grained boundaries between similar subtypes (was 20).
SCANVI_PATIENCE = 50

QUERY_LEIDEN_RESOLUTION = 1.0

# ==============================================================================
# MARKER GENES
# ==============================================================================
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
     "NCAM1","FCGR3A","CD56", "CD4","CD8A","CD8B",
     "FOXP3","IL2RA","IKZF2", "ITGAE","CD103","CXCR6","ZNF683",
     "TOX","TOX2","NR4A1","NR4A2","PRDM1","ID2",
     "TCF7","LEF1","KLF2","ID3","MYB",
     "EOMES","TBX21","RUNX3","RUNX1","ZEB2",
     "BCL6","ASCL2","SH2D1A","MAF", "GATA3","RORA","RORC",
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
    """Strip whitespace and collapse internal spaces (P0-5 / FIX-A)."""
    return pd.Series(labels).astype(str).str.strip().str.replace(r"\s+", " ", regex=True)


def apply_min_cells_filter(labels_series, min_cells, unlabeled_cat):
    """Labels < min_cells cells -> unlabeled_cat (P0-5 / FIX-G)."""
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


def apply_max_cells_cap_stratified(labels_series, strata_series,
                                    max_cells, unlabeled_cat, random_seed=42):
    """
    [v2.4 P1-D] Stratified cap: for each label exceeding max_cells, randomly
    convert surplus cells to unlabeled_cat while preserving the stratum
    (data_source) proportions within that label.

    Naive global capping risks removing nearly all query cells for a dominant
    label, breaking ref/query balance. Stratified capping ensures each
    stratum (ref / query) contributes proportionally to the kept cells.

    Returns:
        capped_series: pd.Series with surplus cells set to unlabeled_cat
        report_df:     pd.DataFrame documenting which labels were capped
    """
    rng    = np.random.default_rng(random_seed)
    capped = labels_series.copy()
    report_rows = []

    df = pd.DataFrame({
        "label":  labels_series.astype(str),
        "strata": strata_series.astype(str),
    }, index=labels_series.index)

    for lbl in df["label"].unique():
        if lbl == unlabeled_cat:
            continue
        idx_lbl = df.index[df["label"] == lbl]
        n_total = len(idx_lbl)
        if n_total <= max_cells:
            continue

        # Proportional keep per stratum
        sub      = df.loc[idx_lbl]
        keep_idx = []
        for _s, idx_s in sub.groupby("strata").groups.items():
            n_s      = len(idx_s)
            n_keep_s = max(1, round(max_cells * n_s / n_total))
            chosen   = rng.choice(list(idx_s), size=min(n_keep_s, n_s), replace=False)
            keep_idx.extend(chosen)

        keep_idx = list(pd.Index(keep_idx).unique())
        # Trim to exactly max_cells if rounding overshot
        if len(keep_idx) > max_cells:
            keep_idx = list(rng.choice(keep_idx, size=max_cells, replace=False))

        cap_idx = idx_lbl.difference(keep_idx)
        capped.loc[cap_idx] = unlabeled_cat

        report_rows.append({
            "label":    lbl,
            "n_before": int(n_total),
            "n_kept":   int(len(keep_idx)),
            "n_capped": int(len(cap_idx)),
        })

    report_df = (pd.DataFrame(report_rows) if report_rows
                 else pd.DataFrame(columns=["label","n_before","n_kept","n_capped"]))
    return capped, report_df


def apply_label_merge(labels_series, merge_map, unlabeled_cat):
    """
    [v2.3 FIX-H] Apply optional coarse label merging before training.
    Cells whose label is in merge_map are relabeled; others are unchanged.
    """
    if not merge_map:
        return labels_series.copy(), pd.DataFrame(columns=["from","to","n_cells"])
    merged = labels_series.copy()
    rows   = []
    for old_lbl, new_lbl in merge_map.items():
        mask = merged == old_lbl
        n    = mask.sum()
        if n > 0:
            merged[mask] = new_lbl
            rows.append({"from": old_lbl, "to": new_lbl, "n_cells": int(n)})
    report_df = (pd.DataFrame(rows) if rows
                 else pd.DataFrame(columns=["from","to","n_cells"]))
    return merged, report_df


def print_imbalance_report(labels_series, unlabeled_cat):
    """
    [v2.3 FIX-J] Print imbalance diagnostics for the labeled cells.
    Rule of thumb: Top1/Total > 30% or Max/Min > 20x suggests problematic
    imbalance that will cause scANVI label merging.
    """
    counts = labels_series.value_counts()
    counts = counts[counts.index != unlabeled_cat]
    if len(counts) == 0:
        print("  [WARN] No labeled cells found!")
        return
    top1_pct  = counts.iloc[0] / counts.sum() * 100
    ratio     = counts.iloc[0] / counts.iloc[-1]
    print(f"\n  [IMBALANCE DIAGNOSTICS]")
    print(f"    Labeled types:   {len(counts)}")
    print(f"    Total labeled:   {counts.sum():,}")
    print(f"    Top 1 / Total:   {top1_pct:.1f}%  "
          f"{'[WARNING >30%]' if top1_pct > 30 else '[OK]'}")
    print(f"    Max / Min ratio: {ratio:.1f}x  "
          f"{'[WARNING >20x]' if ratio > 20 else '[OK]'}")
    print()
    for lbl, cnt in counts.items():
        bar = "#" * min(40, int(cnt / counts.max() * 40))
        print(f"    {lbl:<45} {cnt:>7,}  {bar}")


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
    """sc.pl.embedding wrapper: rasterized=True enforced (v2.9 requirement)."""
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

# P2-8: fail fast if counts layer missing
assert "counts" in adata.layers, (
    "CRITICAL: 'counts' layer not found. "
    "Ensure v2.1 preprocessed h5ad was saved before HVG subsetting."
)
assert CELLTYPIST_LABEL_COL in adata.obs.columns, (
    f"CellTypist column '{CELLTYPIST_LABEL_COL}' not found. "
    f"Available: {list(adata.obs.columns)}"
)
assert CELLTYPIST_CONF_KEY in adata.obs.columns, (
    f"Confidence column '{CELLTYPIST_CONF_KEY}' not found."
)
print(f"\n[OK] CellTypist labels: {adata.obs[CELLTYPIST_LABEL_COL].nunique()} unique types")

# [v2.4 P2-F] Fail-fast guards — all are cheap and catch silent corruption
# that would otherwise surface much later as cryptic errors.
assert adata.obs_names.is_unique, \
    "GUARD: obs_names not unique — reindex operations will silently mis-align data"
assert adata.var_names.is_unique, \
    "GUARD: var_names not unique — HVG forcing and gene lookup will return wrong genes"
assert BATCH_KEY in adata.obs.columns, \
    f"GUARD: BATCH_KEY '{BATCH_KEY}' not in obs.columns"
assert TISSUE_KEY in adata.obs.columns, \
    f"GUARD: TISSUE_KEY '{TISSUE_KEY}' not in obs.columns"
assert CELLTYPIST_DIRECT_LABEL_KEY in adata.obs.columns, \
    f"GUARD: '{CELLTYPIST_DIRECT_LABEL_KEY}' not in obs.columns"
assert adata.obs[CELLTYPIST_CONF_KEY].notna().all(), \
    f"GUARD: NaN values found in '{CELLTYPIST_CONF_KEY}' — confidence filter will be wrong"
assert adata.obs[CELLTYPIST_DIRECT_LABEL_KEY].notna().all(), \
    f"GUARD: NaN values in '{CELLTYPIST_DIRECT_LABEL_KEY}' — label pipeline will fail"
print("[OK] All fail-fast guards passed")
print("\n  Label distribution in checkpoint:")
for lbl, n in adata.obs[CELLTYPIST_LABEL_COL].value_counts().items():
    print(f"    {lbl:<45} {n:>7,}  ({n/adata.n_obs*100:.1f}%)")


# %% [markdown]
# ## Cell 5 — Step 2: Classify Cells by Lineage

# %%
print("\n" + "=" * 80)
print("[Step 2] Classifying cells by CellTypist lineage...")
print("=" * 80)

# [FIX-A] Normalize BEFORE map() to prevent whitespace triggering P1-6 guard
ct_labels = normalize_label_strings(adata.obs[CELLTYPIST_LABEL_COL])
adata.obs[CELLTYPIST_LABEL_COL] = ct_labels.values   # write back normalized

lineage_series = ct_labels.map(CELLTYPIST_LINEAGE_MAP).fillna("other_unknown")
adata.obs["celltypist_lineage"] = lineage_series.values

print("\n[INFO] Lineage distribution:")
lineage_counts     = adata.obs["celltypist_lineage"].value_counts()
lineage_counts_all = lineage_counts.copy()   # preserved for filter_qc.pdf (P2-7)
total = len(adata)
for lin, n in lineage_counts.items():
    flag = "  [KEEP]" if lin in TNK_LINEAGES else "  [REMOVE]"
    print(f"  {lin:<35} {n:>7,}  ({n/total*100:.1f}%){flag}")

# P1-6: raise on unmapped labels
unmapped_labels = ct_labels[lineage_series == "other_unknown"].unique()
if len(unmapped_labels) > 0:
    msg_lines = [f"\n  Unmapped CellTypist labels ({len(unmapped_labels)} unique):"]
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
        print("\n  [WARN] KEEP_OTHER_UNKNOWN=True — skipping guard.")
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

print(f"  Before:  {n_before:,}")
print(f"  Removed: {n_removed:,} non-T/NK cells ({n_removed/n_before*100:.1f}%)")
print(f"  After:   {n_after:,}")

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
    adata_filt = adata_filt[~adata_filt.obs[BATCH_KEY].isin(small_batches)].copy()
    print(f"  Cells after batch drop: {adata_filt.n_obs:,}")
else:
    print(f"  All {batch_counts.shape[0]} batches have >= {MIN_BATCH_CELLS} cells")

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
print("[Step 5] Rebuilding CellTypist filtered label column (P0-2, always fresh)...")
print("=" * 80)
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

# P2-9: Clear old HVG columns before re-selection
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
        print(f"  -> seurat_v3 failed, fallback to cell_ranger")
        sc.pp.highly_variable_genes(adata_filt, layer="counts", n_top_genes=N_HVG,
                                    flavor="cell_ranger", subset=False)
        hvg_method = "cell_ranger"
print(f"  -> HVG method: {hvg_method}")

if "symbol_base" not in adata_filt.var.columns:
    adata_filt.var["symbol_base"] = adata_filt.var_names.str.replace(r"-\d+$", "", regex=True)

if FORCE_MARKERS_IN_HVG:
    n_added, mset = 0, set(FORCED_MARKERS)
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
print("[Step 7] Capturing full-gene log1p matrix for .raw...")
print("=" * 80)
# [v2.4 P0-B] .raw must store log1p-normalized expression, NOT raw counts.
# score_genes() and sc.pl.embedding(use_raw=True) downstream both expect
# log-normalized values: counts-based scoring is confounded by sequencing
# depth, and counts-based feature plots have a different dynamic range.
# Priority: use "log1p" layer if available, otherwise fall back to .X
# (which is log1p normalized in the standard pipeline).
if "log1p" in adata_filt.layers:
    full_log1p = adata_filt.layers["log1p"]
    _raw_src   = "layers['log1p']"
else:
    # .X is log1p normalized in the standard pipeline (from sc.pp.log1p)
    full_log1p = adata_filt.X
    _raw_src   = ".X (assumed log1p)"
if issparse(full_log1p) and not isinstance(full_log1p, csr_matrix):
    full_log1p = csr_matrix(full_log1p)
raw_var = adata_filt.var.copy()
print(f"  Full log1p matrix source: {_raw_src}")
print(f"  Shape: {full_log1p.shape}")


# %% [markdown]
# ## Cell 11 — Step 8: Build Training Subset (HVG) + Label Pipeline
#
# [v2.4] Label pipeline order (changed from v2.3):
#   1. Normalize strings  (FIX-A / P0-5)
#   2. Label merge        (P1-C: moved BEFORE min-cells filter)
#      -> merging can rescue biologically real rare labels
#   3. Min-cells filter   (P0-5 + P2-G whitelist)
#      -> rare labels below threshold -> Unknown, EXCEPT whitelist entries
#   4. Max-cells cap      (P1-D: stratified by data_source)
#      -> dominant labels capped proportionally
#   5. Imbalance report   (FIX-J)

# %%
print("\n" + "=" * 80)
print("[Step 8] Building training subset + label pipeline...")
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

# ── Step A: Normalize label strings (P0-5 / FIX-A) ───────────────────────────
working_labels = normalize_label_strings(adata_train.obs[label_source_col])
n_ws_fixed     = int((working_labels != adata_train.obs[label_source_col].astype(str)).sum())
if n_ws_fixed > 0:
    print(f"  [P0-5] Normalized whitespace in {n_ws_fixed:,} label strings")

# ── Step B: Optional coarse label merge (P1-C: NOW BEFORE min-cells filter) ──
# [v2.4 P1-C] LABEL_MERGE_MAP is applied FIRST so that merging can rescue
# rare labels that are biologically real but numerically small individually.
# Example: two labels each with 30 cells would both be filtered out by
# MIN_CELLS_PER_LABEL=50 if merge is applied after — but if merged first
# into a single label with 60 cells, they survive the filter.
if LABEL_MERGE_MAP:
    working_labels, merge_report = apply_label_merge(
        working_labels, LABEL_MERGE_MAP, UNLABELED_CATEGORY
    )
    print(f"\n  [P1-C] Label merge applied ({len(merge_report)} entries):")
    for _, row in merge_report.iterrows():
        print(f"    '{row['from']}' -> '{row['to']}'  n={row['n_cells']:,}")
    merge_report.to_csv(output_dir / f"{OUTPUT_PREFIX}_label_merges.csv", index=False)
else:
    print("  [P1-C] LABEL_MERGE_MAP={} — no coarse merging applied")

# ── Step C: Min-cells filter with whitelist (P0-5 + P2-G) ────────────────────
# [v2.4 P2-G] Labels in RARE_LABEL_WHITELIST are exempt from the min-cells
# filter even if their count < MIN_CELLS_PER_LABEL.
_counts_now = working_labels.value_counts()
rare_to_filter = [
    lbl for lbl in _counts_now.index
    if _counts_now[lbl] < MIN_CELLS_PER_LABEL
    and lbl != UNLABELED_CATEGORY
    and lbl not in RARE_LABEL_WHITELIST
]
whitelisted_rescued = [
    lbl for lbl in _counts_now.index
    if _counts_now[lbl] < MIN_CELLS_PER_LABEL
    and lbl != UNLABELED_CATEGORY
    and lbl in RARE_LABEL_WHITELIST
]

rare_report_rows = [{"label": lbl, "n_cells": int(_counts_now[lbl]),
                     "action": f"-> {UNLABELED_CATEGORY}"} for lbl in rare_to_filter]
rare_report = (pd.DataFrame(rare_report_rows) if rare_report_rows
               else pd.DataFrame(columns=["label","n_cells","action"]))

if rare_to_filter:
    working_labels = working_labels.replace(rare_to_filter, UNLABELED_CATEGORY)
    print(f"\n  [P0-5] {len(rare_to_filter)} labels below "
          f"MIN_CELLS_PER_LABEL={MIN_CELLS_PER_LABEL} -> Unknown:")
    for _, row in rare_report.iterrows():
        print(f"    {row['label']:<45} n={row['n_cells']:,}  {row['action']}")
    rare_report.to_csv(output_dir / f"{OUTPUT_PREFIX}_rare_labels.csv", index=False)
else:
    print(f"  [OK] No labels below MIN_CELLS_PER_LABEL={MIN_CELLS_PER_LABEL}")

if whitelisted_rescued:
    print(f"\n  [P2-G] {len(whitelisted_rescued)} whitelisted rare labels KEPT despite low count:")
    for lbl in whitelisted_rescued:
        print(f"    {lbl:<45} n={int(_counts_now[lbl]):,}  [WHITELIST — kept]")

# ── Step D: Max-cells cap — stratified by data_source (P1-D) ─────────────────
cap_report = pd.DataFrame(columns=["label","n_before","n_kept","n_capped"])
if MAX_CELLS_PER_LABEL is not None:
    working_labels, cap_report = apply_max_cells_cap_stratified(
        working_labels,
        strata_series=adata_train.obs["data_source"],
        max_cells=MAX_CELLS_PER_LABEL,
        unlabeled_cat=UNLABELED_CATEGORY,
        random_seed=RANDOM_SEED,
    )
    if len(cap_report) > 0:
        print(f"\n  [P1-D] MAX_CELLS_PER_LABEL={MAX_CELLS_PER_LABEL} "
              f"(stratified by data_source): {len(cap_report)} labels capped:")
        for _, row in cap_report.iterrows():
            print(f"    {row['label']:<45} "
                  f"{row['n_before']:>7,} -> {row['n_kept']:>7,} "
                  f"(capped {row['n_capped']:,})")
        cap_report.to_csv(output_dir / f"{OUTPUT_PREFIX}_capped_labels.csv", index=False)
    else:
        print(f"  [P1-D] MAX_CELLS_PER_LABEL={MAX_CELLS_PER_LABEL}: no label exceeded cap")
else:
    print("  [P1-D] MAX_CELLS_PER_LABEL=None — capping disabled")

# ── Assign final scanvi_labels ────────────────────────────────────────────────
adata_train.obs["scanvi_labels"] = working_labels.values
adata_train.obs["scanvi_labels"] = adata_train.obs["scanvi_labels"].astype("category")
if UNLABELED_CATEGORY not in adata_train.obs["scanvi_labels"].cat.categories:
    adata_train.obs["scanvi_labels"] = (
        adata_train.obs["scanvi_labels"].cat.add_categories([UNLABELED_CATEGORY])
    )

# P0-1: Write scanvi_labels back to adata_filt via reindex
adata_filt.obs["scanvi_labels"] = (
    adata_train.obs["scanvi_labels"]
    .reindex(adata_filt.obs_names)
    .astype("category")
)
print(f"\n  -> scanvi_labels written back to adata_filt "
      f"({adata_filt.obs['scanvi_labels'].nunique()} categories incl. Unknown)")

# ── P0-3: Validate required covariates ───────────────────────────────────────
required_cont = ["pct_counts_mt","stress_score","S_score","G2M_score"]
missing_cont  = [c for c in required_cont if c not in adata_train.obs.columns]
if missing_cont:
    raise ValueError(
        f"P0-3 GUARD: Required covariates missing: {missing_cont}. "
        "Re-run v2.1 pipeline."
    )
print(f"  [OK] All covariates present: {required_cont}")

# ── Step E: Pre-training label stats + imbalance report (FIX-J) ──────────────
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

print_imbalance_report(adata_train.obs["scanvi_labels"], UNLABELED_CATEGORY)

# [v2.3 FIX-F] pandas 3.0 compat: use .columns=[...] not rename(columns={"index":...})
# In pandas 3.0, reset_index() names the index column after the original Series
# name (here "scanvi_labels"), NOT "index". rename({"index":"label"}) is a no-op.
_pretrain_stats = label_counts_pretrain.rename("n_cells_pretrain").reset_index()
_pretrain_stats.columns = ["label", "n_cells_pretrain"]   # FIX-F: explicit assignment

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
print(f"[Step 10] Training scANVI (patience={SCANVI_PATIENCE})...")
print("=" * 80)

scanvi_model = scvi.model.SCANVI.from_scvi_model(
    scvi_model, adata=adata_train,
    labels_key="scanvi_labels",
    unlabeled_category=UNLABELED_CATEGORY
)
scanvi_train_kwargs = {
    "max_epochs": MAX_EPOCHS_SCANVI, "batch_size": BATCH_SIZE,
    "early_stopping": True, "early_stopping_patience": SCANVI_PATIENCE,  # FIX-I
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
# [FIX-F] pandas 3.0 compat: explicit column assignment
_pred_reset = pred_counts.reset_index()
_pred_reset.columns = ["label", "n_cells_predicted"]

_audit = _pretrain_stats.merge(_pred_reset, on="label", how="outer").fillna(0)
_audit["n_cells_pretrain"]  = _audit["n_cells_pretrain"].astype(int)
_audit["n_cells_predicted"] = _audit["n_cells_predicted"].astype(int)

_missing = _audit[
    (_audit["n_cells_pretrain"] > 0) &
    (_audit["label"] != UNLABELED_CATEGORY) &
    (_audit["n_cells_predicted"] == 0)
]
if len(_missing) > 0:
    print(f"\n  [WARNING] {len(_missing)} training labels ABSENT from predictions:")
    print("  Labels may have been absorbed by dominant classes.")
    print("  Consider: increasing MAX_CELLS_PER_LABEL cap, merging similar")
    print("  subtypes in LABEL_MERGE_MAP, or reviewing dotplot for distinctiveness.")
    for _, row in _missing.iterrows():
        print(f"    '{row['label']}': trained {int(row['n_cells_pretrain']):,} cells, 0 predicted")
else:
    print("  [OK] All training labels appear in predictions — no merging detected.")

_new = _audit[
    (_audit["n_cells_pretrain"] == 0) &
    (_audit["label"] != UNLABELED_CATEGORY) &
    (_audit["n_cells_predicted"] > 0)
]
if len(_new) > 0:
    print(f"\n  [INFO] {len(_new)} prediction labels not in training (unexpected):")
    for _, row in _new.iterrows():
        print(f"    '{row['label']}': {int(row['n_cells_predicted']):,} predicted")

print(f"\n  Training labels:   {int((_audit['n_cells_pretrain'] > 0).sum())}")
print(f"  Prediction labels: {int((_audit['n_cells_predicted'] > 0).sum())}")
_audit.to_csv(output_dir / f"{OUTPUT_PREFIX}_label_audit.csv", index=False)
print(f"  [OK] Audit: {OUTPUT_PREFIX}_label_audit.csv")
print("-" * 60)


# %% [markdown]
# ## Cell 15 — Step 12: Attach .raw (Must Precede Scoring)
#
# [FIX-C] .raw attached BEFORE CD4/CD8 scoring so use_raw=True is safe.
# [v2.4 P0-B] .raw stores log1p expression (captured in Step 7 as full_log1p).

# %%
print("\n" + "=" * 80)
print("[Step 12] Attaching .raw (log1p expression — must precede scoring)...")
print("=" * 80)
from anndata import AnnData as _AnnData
adata_filt.raw = _AnnData(X=full_log1p, obs=adata_filt.obs.copy(), var=raw_var)
print(f"  .raw source: {_raw_src}")
print(f"  .raw shape:  {adata_filt.raw.n_obs} cells x {adata_filt.raw.n_vars} genes")


# %% [markdown]
# ## Cell 16 — Step 13: CD4/CD8 Scoring QC
#
# [FIX-B] NK mask uses direct label set (covers cycling NK).
# [FIX-C] use_raw=True hardcoded (safe after Step 12).

# %%
print("\n" + "=" * 80)
print("[Step 13] CD4/CD8 QC scoring (NK + Cycling NK excluded)...")
print("=" * 80)

NK_DIRECT_LABELS = {
    "NK cells", "CD16+ NK cells", "CD16- NK cells", "Cycling NK cells",
}

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

# [v2.4 P0-A] Explicit `is None` check — `or 0.0` would trigger
# ValueError: truth value of array with more than one element is ambiguous
cd4_score = _score(CD4_SCORE_GENES)
cd8_score = _score(CD8_SCORE_GENES)
nk_score  = _score(NK_SCORE_GENES)

adata_filt.obs["cd4_score"] = 0.0 if cd4_score is None else cd4_score
adata_filt.obs["cd8_score"] = 0.0 if cd8_score is None else cd8_score
adata_filt.obs["nk_score"]  = 0.0 if nk_score  is None else nk_score

# [FIX-B] NK mask on direct label (covers Cycling NK cells)
_nk_mask = (
    normalize_label_strings(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY])
    .isin(NK_DIRECT_LABELS)
    .values
)
_cd4_pos = adata_filt.obs["cd4_score"].values > CD4_SCORE_THRESH
_cd8_pos = adata_filt.obs["cd8_score"].values > CD8_SCORE_THRESH

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

lat_s    = scvi_model.get_latent_representation(adata_train)
lat_s_df = pd.DataFrame(lat_s, index=adata_train.obs_names,
                         columns=[f"scVI_{i}" for i in range(lat_s.shape[1])])
lat_s_al = lat_s_df.reindex(adata_filt.obs_names)
if lat_s_al.isna().any().any():
    raise ValueError("Missing scVI latent after reindex!")
adata_filt.obsm["X_scVI"] = lat_s_al.values

run_query_only_leiden(adata_filt, resolution=QUERY_LEIDEN_RESOLUTION)

sc.pp.neighbors(adata_filt, use_rep="X_scVI", n_neighbors=30,
                random_state=RANDOM_SEED, key_added="neighbors_scVI")
sc.tl.umap(adata_filt, random_state=RANDOM_SEED, neighbors_key="neighbors_scVI")
adata_filt.obsm["X_umap_scVI"] = adata_filt.obsm["X_umap"].copy()
print("  -> X_umap_scVI saved")

sc.pp.neighbors(adata_filt, use_rep="X_scANVI", n_neighbors=30,
                random_state=RANDOM_SEED, key_added="neighbors_scANVI")
sc.tl.umap(adata_filt, random_state=RANDOM_SEED, neighbors_key="neighbors_scANVI")
adata_filt.obsm["X_umap_scANVI"] = adata_filt.obsm["X_umap"].copy()
adata_filt.obsm["X_umap"]        = adata_filt.obsm["X_umap_scANVI"].copy()
print("  -> X_umap_scANVI saved (default X_umap)")

# [v2.4 P1-E] This UMAP operator is an APPROXIMATE projector for mapping new
# cells into a similar 2D space in the future. It is NOT the same object used
# to produce X_umap_scANVI above (sc.tl.umap runs its own internal UMAP
# instance with specific connectivity parameters). umap_op.transform() will
# give a visually similar but not pixel-identical embedding compared to the
# figures saved in this pipeline.
umap_op = UMAP(n_neighbors=30, n_components=2, min_dist=0.5, metric="euclidean",
               random_state=RANDOM_SEED)
umap_op.fit(adata_filt.obsm["X_scANVI"])
joblib.dump(umap_op, output_dir / f"{OUTPUT_PREFIX}_umap_scanvi_operator.joblib")
print("  -> UMAP approximate projector saved (for future query projection)")


# %% [markdown]
# ## Cell 18 — Step 15: Run Log

# %%
print("\n" + "=" * 80)
print("[Step 15] Writing run log...")
print("=" * 80)
log_ts   = datetime.now().isoformat()
out_h5ad = output_dir / f"{OUTPUT_PREFIX}_results.h5ad"

pipeline_log = {
    "version":               "2.4",  "timestamp": log_ts,
    "input_h5ad":            INPUT_H5AD,
    "output_dir":            str(output_dir), "output_prefix": OUTPUT_PREFIX,
    "planned_output":        str(out_h5ad),
    "elapsed_min":           round((time.time() - PIPELINE_START) / 60, 2),
    "n_obs_before_filter":   n_before,
    "n_obs_after_filter":    adata_filt.n_obs,
    "n_removed":             n_removed,
    "pct_removed":           round(n_removed / n_before * 100, 2),
    "n_vars":                int(adata_filt.n_vars),
    "n_hvg":                 n_hvg_final,
    "hvg_method":            hvg_method,
    "gpu_available":         bool(gpu_available),
    "scvi_n_latent":         SCVI_N_LATENT,
    "min_cells_per_label":   MIN_CELLS_PER_LABEL,
    "max_cells_per_label":   MAX_CELLS_PER_LABEL,
    "n_rare_labels":         len(rare_report),
    "n_capped_labels":       len(cap_report),
    "label_merge_map":       LABEL_MERGE_MAP,
    "scanvi_patience":       SCANVI_PATIENCE,
    "architecture":          "CellTypist-filtered T/NK-only retrain (v2.4)",
    "scanvi_label_source":   label_source_col,
    "min_batch_cells":       MIN_BATCH_CELLS,
    "small_batches_dropped": small_batches if small_batches else [],
}
anndata_structure = {
    "timestamp": log_ts, "shape": [int(adata_filt.n_obs), int(adata_filt.n_vars)],
    "raw": {"present": True, "n_vars": adata_filt.raw.n_vars},
    "obs_columns": list(adata_filt.obs.columns),
    "obsm_keys":   sorted(adata_filt.obsm.keys()),
    "uns_keys":    "PLACEHOLDER",   # FIX-D: filled after uns write
}
# [FIX-D] Write uns first, then record keys
adata_filt.uns["pipeline_log"]      = pipeline_log
adata_filt.uns["anndata_structure"] = {}
anndata_structure["uns_keys"] = sorted(str(k) for k in adata_filt.uns.keys())
adata_filt.uns["anndata_structure"] = anndata_structure

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
    "scanvi_labels", "scanvi_pred", "cd4_cd8_by_score", "leiden_query",
    CELLTYPIST_DIRECT_LABEL_KEY, CELLTYPIST_DIRECT_FILT_KEY,
]
for col in _cat_cols:
    for ad in [adata_filt, adata_train]:
        if col in ad.obs.columns:
            ad.obs[col] = ad.obs[col].astype("category")

scanvi_model.save(output_dir / f"{OUTPUT_PREFIX}_scanvi_model", overwrite=True)
scvi_model.save(output_dir   / f"{OUTPUT_PREFIX}_scvi_model",   overwrite=True)

config = {
    "version":   "2.4", "timestamp": datetime.now().isoformat(),
    "input":     INPUT_H5AD,
    "n_hvg":     n_hvg_final, "hvg_method": hvg_method,
    "scvi_n_latent": SCVI_N_LATENT,
    "label_pipeline": {
        "min_cells_per_label":  MIN_CELLS_PER_LABEL,
        "max_cells_per_label":  MAX_CELLS_PER_LABEL,
        "label_merge_map":      LABEL_MERGE_MAP,
        "scanvi_patience":      SCANVI_PATIENCE,
    },
    "architecture": {
        "description":       "CellTypist-filtered T/NK-only retrain v2.4",
        "filter_policy":     "keep TNK_LINEAGES, remove contaminants",
        "tnk_lineages":      sorted(TNK_LINEAGES),
        "n_removed":         n_removed,
        "pct_removed":       round(n_removed / n_before * 100, 2),
        "scanvi_labels":     label_source_col,
        "fixes": [
            "FIX-F: pandas 3.0 compat — .columns=[...] (KeyError fix)",
            f"FIX-G / P1-D: stratified MAX_CELLS_PER_LABEL={MAX_CELLS_PER_LABEL}",
            "FIX-H / P1-C: LABEL_MERGE_MAP now applied BEFORE min-cells filter",
            f"FIX-I: SCANVI_PATIENCE={SCANVI_PATIENCE}",
            "FIX-J: imbalance diagnostics printed before training",
            "P0-A: numpy `is None` check replaces `or 0.0` crash fix",
            "P0-B: .raw stores log1p (not counts); use_raw=True is semantically correct",
            f"P2-G: RARE_LABEL_WHITELIST={sorted(RARE_LABEL_WHITELIST)} exempt from min-cells",
            "P2-F: obs_names/var_names/BATCH_KEY/TISSUE_KEY/notna guards added",
            "P1-E: UMAP operator annotated as approximate projector",
        ],
    },
    "scanvi_labels": list(label_order),
    "umap_spaces": {
        "X_umap": "DEFAULT (scANVI-based)",
        "X_umap_scVI": "scVI latent UMAP",
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
    ("data_source",               "Data Source"),
    ("celltypist_lineage",        "CellTypist Lineage"),
    (CELLTYPIST_DIRECT_LABEL_KEY, "CellTypist Direct"),
    (CELLTYPIST_DIRECT_FILT_KEY,  "CellTypist Filtered"),
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
    ("CD3D",  "CD3D (pan-T)"),         ("CD4",   "CD4"),
    ("CD8A",  "CD8A"),                  ("GNLY",  "GNLY (NK)"),
    ("NCAM1", "NCAM1 (NK)"),           ("FOXP3", "FOXP3 (Treg)"),
    ("HAVCR2","HAVCR2 (Exhaustion)"), ("MKI67", "MKI67 (Prolif)"),
]
for i, (gene, title) in enumerate(marker_panels):
    row_i, col_i = 2 + i // 4, i % 4
    ax = fig.add_subplot(gs[row_i, col_i])
    if adata_filt.raw is not None and gene in adata_filt.raw.var_names:
        sc.pl.embedding(adata_filt, basis="X_umap_scANVI", color=gene,
                        ax=ax, show=False, title=title, cmap="Reds", s=8,
                        rasterized=True, frameon=False, use_raw=True)
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
    ("scanvi_pred",               "scANVI (T/NK-only retrain v2.4)"),
]):
    if key in adata_filt.obs.columns:
        sc.pl.embedding(adata_filt, basis="X_umap_scANVI", color=key,
                        ax=ax, show=False, title=title,
                        legend_loc="on data", legend_fontsize=6, s=8,
                        rasterized=True, frameon=False)
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
        sc.pl.embedding(adata_filt, basis=basis, color=color,
                        ax=axes3[row_i, col_i], show=False,
                        title=f"{color} ({label})",
                        legend_loc=loc, legend_fontsize=6, s=8,
                        rasterized=True, frameon=False)
plt.tight_layout()
fig3.savefig(output_dir / f"{OUTPUT_PREFIX}_umap_comparison.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig3)
print(f"  -> {OUTPUT_PREFIX}_umap_comparison.pdf")

# ── Panel 4: NK vs T cell marker separation ───────────────────────────────────
nk_t_genes = ["NCAM1","FCGR3A","KIR2DL1","KIR3DL1",
               "CD3D","CD3E","TRAC", "CD4","CD8A", "GNLY","NKG7"]
_nkt_present = [g for g in nk_t_genes if g in (
    adata_filt.raw.var_names if adata_filt.raw is not None else adata_filt.var_names)]
if _nkt_present:
    ncols = (len(_nkt_present) + 1) // 2
    fig4, axes4 = plt.subplots(2, ncols, figsize=(4 * ncols, 8))
    axes4_flat = np.array(axes4).flatten()
    for ax, gene in zip(axes4_flat, _nkt_present):
        sc.pl.embedding(adata_filt, basis="X_umap_scANVI", color=gene,
                        ax=ax, show=False, title=gene, cmap="Reds", s=8,
                        rasterized=True, frameon=False, use_raw=(adata_filt.raw is not None))
    for ax in axes4_flat[len(_nkt_present):]:
        ax.set_visible(False)
    plt.suptitle("NK vs T cell marker separation", fontsize=14, y=1.01)
    plt.tight_layout()
    fig4.savefig(output_dir / f"{OUTPUT_PREFIX}_nk_t_separation.pdf",
                 dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
    plt.close(fig4)
    print(f"  -> {OUTPUT_PREFIX}_nk_t_separation.pdf")

# ── Panel 5: Filter QC — before vs after (P2-7) ──────────────────────────────
fig5, axes5 = plt.subplots(1, 2, figsize=(16, 7))
_all_vals   = lineage_counts_all.values
_all_labels = lineage_counts_all.index.tolist()
_colors     = ["#d62728" if "contaminant" in l or l == "other_unknown"
               else "#1f77b4" for l in _all_labels]
axes5[0].pie(_all_vals, labels=_all_labels, autopct="%1.1f%%",
             startangle=140, colors=_colors, textprops={"fontsize": 7})
axes5[0].set_title(
    f"All Lineages — Before Filter\n(n={n_before:,})\n"
    "Blue = T/NK kept  |  Red = removed",
    fontsize=9
)
_ret_lin = adata_filt.obs["celltypist_lineage"].value_counts()
axes5[1].pie(_ret_lin.values, labels=_ret_lin.index, autopct="%1.1f%%",
             startangle=140, textprops={"fontsize": 7})
axes5[1].set_title(
    f"Retained T/NK — After Filter\n"
    f"(n={n_after:,}  |  removed {n_removed:,} = {n_removed/n_before*100:.1f}%)\n"
    f"Batches dropped (<{MIN_BATCH_CELLS} cells): {len(small_batches)}",
    fontsize=9
)
plt.tight_layout()
fig5.savefig(output_dir / f"{OUTPUT_PREFIX}_filter_qc.pdf",
             dpi=FIGURE_DPI, bbox_inches="tight", format=FIGURE_FORMAT)
plt.close(fig5)
print(f"  -> {OUTPUT_PREFIX}_filter_qc.pdf")

# ── Dotplot: T/NK marker overview ─────────────────────────────────────────────
dotplot_genes = [
    "CD3D","CD3E","TRAC", "CD4","IL7R","MAL", "CD8A","CD8B",
    "NCAM1","FCGR3A","KLRB1","KLRC1","KLRD1","KLRF1","KIR2DL1","KIR2DL3","KIR3DL1",
    "GNLY","NKG7","GZMB","GZMK","GZMA","GZMH","PRF1",
    "FOXP3","IL2RA","IKZF2","CTLA4","TNFRSF18",
    "CXCR5","BCL6","ICOS","SH2D1A",
    "GATA3","IL4","IL13","CRTH2","CCR4",
    "TBX21","IFNG","CXCR3", "RORC","IL17A","CCR6",
    "CCR7","SELL","TCF7","LEF1","KLF2",
    "ITGAE","CD69","CXCR6","ZNF683",
    "PDCD1","HAVCR2","LAG3","TIGIT","TOX","NR4A1",
    "CD38","HLA-DRA","TNFRSF4","TNFRSF9",
    "SLC4A10","TRDC","MKI67","TOP2A","PCNA",
]
_use_raw_dp = adata_filt.raw is not None
_dp_genes   = [g for g in dotplot_genes if g in (
    adata_filt.raw.var_names if _use_raw_dp else adata_filt.var_names)]
if _dp_genes:
    try:
        dp = sc.pl.dotplot(
            adata_filt, var_names=_dp_genes, groupby="scanvi_pred",
            use_raw=_use_raw_dp, standard_scale="var", cmap="Reds",
            show=False, return_fig=True,
            figsize=(max(14, len(_dp_genes)*0.55),
                     max(6, adata_filt.obs["scanvi_pred"].nunique()*0.4)),
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
print("T/NK PIPELINE v2.4 COMPLETE")
print("=" * 80)
print(f"  Input cells:             {n_before:,}")
print(f"  Removed:                 {n_removed:,}  ({n_removed/n_before*100:.1f}%)")
print(f"  Retained:                {n_after:,}")
print(f"  Batches used:            {adata_filt.obs[BATCH_KEY].nunique()}")
print(f"  Elapsed:                 {elapsed:.1f} min")
print(f"\n  scVI n_latent:           {SCVI_N_LATENT}")
print(f"  scANVI patience:         {SCANVI_PATIENCE}")
print(f"  scANVI labels from:      {label_source_col}")
print(f"  MIN_CELLS_PER_LABEL:     {MIN_CELLS_PER_LABEL}  ({len(rare_report)} filtered)")
print(f"  MAX_CELLS_PER_LABEL:     {MAX_CELLS_PER_LABEL}  ({len(cap_report)} capped)")
print(f"  LABEL_MERGE_MAP entries: {len(LABEL_MERGE_MAP)}")
n_unknown_final = int(
    (adata_filt.obs["scanvi_labels"].astype(str) == UNLABELED_CATEGORY).sum()
)
print(f"  Unknown (low-conf + rare + capped): {n_unknown_final:,}")
print("\nTop scANVI predictions:")
print(adata_filt.obs["scanvi_pred"].value_counts().head(12))
print("\nTop CellTypist labels (post-filter):")
print(adata_filt.obs[CELLTYPIST_DIRECT_LABEL_KEY].value_counts().head(12))
print("\nOutput files:")
for f in sorted(output_dir.glob(f"{OUTPUT_PREFIX}*")):
    print(f"  {f.name}")
print("=" * 80)
