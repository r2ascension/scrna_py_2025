"""
patch_runlog_v1_10.py
Completes the run log for allcells_combined_20260319_v1_10_L1.h5ad
after the pipeline crashed at Cell 24 due to residual 'cell_type_final' reference.
The h5ad is already saved correctly. This script regenerates UMAP figures and writes
the missing log files.
"""

import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

# -----------------------------------------------------------------------
OUTPUT_DIR   = Path("/home/h2048/data/py/20260319/allcells_combined_scanvi")
FIG_DIR      = OUTPUT_DIR / "figures"
H5AD_PATH    = OUTPUT_DIR / "allcells_combined_20260319_v1_10_L1.h5ad"
VERSION      = "1_10_L1"
BATCH_KEY    = "sample"
LABELS_KEY   = "scanvi_label"
UNLABELED_CATEGORY = "Unknown"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SCVI_N_LATENT               = 150
SCVI_N_LAYERS               = 2
SCVI_N_HIDDEN               = 128
SCVI_DROPOUT_RATE           = 0.2
SCVI_EARLY_STOPPING_PATIENCE = 45
SCVI_LEARNING_RATE          = 1e-3
SCVI_USE_LAYER_NORM         = "both"
SCVI_USE_BATCH_NORM         = "none"
SCVI_ENCODE_COVARIATES      = True
SCANVI_LEARNING_RATE        = 5e-4
SCANVI_EARLY_STOPPING_PATIENCE = 30
CONTINUOUS_COVARIATES       = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
CATEGORICAL_COVARIATES      = []
LABEL_TIER                  = "L1"
# -----------------------------------------------------------------------

def _fmt_float(val, fallback="not_set"):
    return f"{val:.4f}" if isinstance(val, (float, np.floating)) else str(fallback)

print(f"Loading: {H5AD_PATH}")
adata = sc.read_h5ad(H5AD_PATH)
print(f"Loaded: {adata.n_obs:,} x {adata.n_vars:,}")
print(f"obsm keys: {list(adata.obsm.keys())}")

# -----------------------------------------------------------------------
# Regenerate UMAP figures
# Uses existing X_scvi / X_scanvi latent embeddings already in obsm.
# Recomputes neighbors (n_neighbors=30) and UMAP (min_dist=0.5) to match
# v1.10 parameters and overwrite the figures from the failed run.
# -----------------------------------------------------------------------
sc.settings.vector_friendly = True
sc.settings.n_jobs = 8

_cols = ["lineage_source", "cell_type_final_l1", "cell_type_L2",
         BATCH_KEY, "scanvi_confidence"]

for rep, basis_out, lbl in [
    ("X_scvi",   "X_umap_scvi",   "scVI"),
    ("X_scanvi", "X_umap_scanvi", "scANVI"),
]:
    if rep not in adata.obsm:
        print(f"[WARNING] {rep} not found in obsm, skipping {lbl} UMAP")
        continue

    print(f"\n[{lbl}] Computing neighbors (n_neighbors=30)...")
    sc.pp.neighbors(adata, use_rep=rep, n_neighbors=30)

    print(f"[{lbl}] Computing UMAP (min_dist=0.5)...")
    sc.tl.umap(adata, min_dist=0.5)
    adata.obsm[basis_out] = adata.obsm["X_umap"].copy()

    print(f"[{lbl}] Saving figure...")
    fig, axes = plt.subplots(1, len(_cols), figsize=(6 * len(_cols), 5))
    for ax, col in zip(axes, _cols):
        sc.pl.embedding(
            adata, basis=basis_out, color=col,
            ax=ax, show=False, title=col,
            legend_loc="right margin" if adata.obs[col].nunique() <= 30 else "none",
            frameon=False,
        )
    fig.suptitle(f"{lbl} UMAP — All Lineages (L1 scANVI refined)", y=1.02, fontsize=14)
    fig.tight_layout()
    out_fig = FIG_DIR / f"allcells_{lbl.lower()}_umap_overview.pdf"
    fig.savefig(out_fig, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"[OK] {lbl} UMAP saved: {out_fig}")

# -----------------------------------------------------------------------
# Write run log
# -----------------------------------------------------------------------
n_lab = int((adata.obs[LABELS_KEY].astype(str) != UNLABELED_CATEGORY).sum())
n_unk = int((adata.obs[LABELS_KEY].astype(str) == UNLABELED_CATEGORY).sum())
_cov  = adata.uns.get("scvi_registered_covariates", {
    "continuous": CONTINUOUS_COVARIATES,
    "categorical": CATEGORICAL_COVARIATES,
    "encode_covariates": SCVI_ENCODE_COVARIATES,
})
_pipeline_name = f"allcells_combined_v{VERSION}"

_log = [
    "="*80, f"{_pipeline_name} Run Log  (v{VERSION})", "="*80,
    f"written_at                    : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    f"note                          : log + figures completed by patch_runlog_v1_10.py",
    f"output_h5ad                   : {H5AD_PATH}",
    f"n_cells                       : {adata.n_obs:,}",
    f"n_genes_hvg                   : {adata.n_vars:,}",
    f"n_genes_raw (.raw union space): {adata.raw.n_vars if adata.raw else 'None'}",
    f"n_batches                     : {adata.obs[BATCH_KEY].nunique():,}",
    f"hvg_method                    : {adata.uns.get('hvg_method','unknown')}",
    "",
    "[v1.9/v1.10] scVI hyperparameters:",
    f"  n_latent                    : {SCVI_N_LATENT}",
    f"  n_layers                    : {SCVI_N_LAYERS}",
    f"  n_hidden                    : {SCVI_N_HIDDEN}",
    f"  dropout_rate                : {SCVI_DROPOUT_RATE}",
    f"  use_layer_norm              : {SCVI_USE_LAYER_NORM}",
    f"  use_batch_norm              : {SCVI_USE_BATCH_NORM}",
    f"  encode_covariates           : {SCVI_ENCODE_COVARIATES}",
    f"  early_stopping_patience     : {SCVI_EARLY_STOPPING_PATIENCE}",
    f"  lr                          : {SCVI_LEARNING_RATE}",
    "",
    "[v1.9/v1.10] scANVI hyperparameters:",
    f"  lr                          : {SCANVI_LEARNING_RATE}",
    f"  n_samples_per_label         : removed in v1.10-4",
    f"  early_stopping_patience     : {SCANVI_EARLY_STOPPING_PATIENCE}",
    f"  weight_decay                : 0.0",
    "",
    "[v1.10] UMAP parameters (patch):",
    f"  n_neighbors                 : 30",
    f"  min_dist                    : 0.5",
    "",
    "[v1.10] Covariates:",
    f"  continuous  : {_cov.get('continuous', [])}",
    f"  categorical : {_cov.get('categorical', [])}",
    "",
    f"label_tier                    : {LABEL_TIER}",
    f"n_labeled                     : {n_lab:,}",
    f"n_unknown                     : {n_unk:,}",
    f"dropped_unknown               : {int(adata.uns.get('dropped_unknown_cells', 0)):,}",
    f"agreement_rate                : {_fmt_float(adata.uns.get('agreement_rate_overall'))}  [train-set re-substitution, NOT independent validation]",
    "",
    "lineage_source distribution:",
    adata.obs["lineage_source"].value_counts().to_string(),
    "",
    "cell_type_final_l1 (L1 scANVI) distribution:",
    adata.obs["cell_type_final_l1"].value_counts().to_string(),
    "",
    "cell_type_L2 distribution:",
    adata.obs["cell_type_L2"].value_counts().to_string(),
    "",
    "cell_type_input_l3 (original labels) distribution (top 40):",
    adata.obs["cell_type_input_l3"].value_counts().head(40).to_string(),
    "",
    "scanvi_confidence percentiles:",
    adata.obs["scanvi_confidence"].describe().to_string(),
]

log_path = OUTPUT_DIR / "pipeline_run_log.txt"
log_path.write_text("\n".join(_log) + "\n", encoding="utf-8")
print(f"\n[OK] pipeline_run_log.txt written: {log_path}")

(OUTPUT_DIR / "anndata_structure.txt").write_text(
    f"See anndata_structure_20260319_v{VERSION}.json\n"
    f"n_obs={adata.n_obs:,} n_vars={adata.n_vars:,}\n",
    encoding="utf-8"
)
print(f"[OK] anndata_structure.txt written")
print("\n[DONE] Patch complete.")
