#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TNK refined-label scanvi-only rerun after label corrections.

Purpose
-------
1. Reuse the 2026-04-13 filtered TNK scVI rerun object (after CHOIR-cluster removal).
2. Correct refined supervision labels:
   - CD4 Naive/TCM -> CD4 Naive
   - CD4 Tfr       -> CD4 Tfh
3. Reload the existing filtered scVI model and retrain scANVI only.
4. Save a fresh H5AD / scANVI model for downstream tissue comparison.

Rationale
---------
The user is satisfied with the filtered scVI stage and only wants the label
judgment corrected from the scANVI stage onward. Therefore this script does
NOT retrain scVI from scratch; it reuses the filtered 0413 scVI model and
reruns scANVI with corrected refined labels.
"""

# ============================================================================
# 0. Configuration
# ============================================================================
INPUT_H5AD = (
    "/home/h2048/data/py/0413/tnk_scvi_scanvi_refined_rerun_rm_choir_20260413/"
    "adata_tnk_scanvi_refined_rerun_rm_choir_20260413.h5ad"
)
SCVI_MODEL_DIR = (
    "/home/h2048/data/py/0413/tnk_scvi_scanvi_refined_rerun_rm_choir_20260413/"
    "tnk_scvi_refined_rerun_model"
)
OUTPUT_DIR = \
    "/home/h2048/data/py/0414/tnk_scanvi_relabel_rerun_20260414"

LABELS_KEY = "scanvi_label_refined"
LABELS_BACKUP_KEY = "scanvi_label_refined_pre_20260414_relabel"
PRED_KEY = "scanvi_pred_refined"
PRED_PROB_KEY = "scanvi_pred_prob_refined"
SCANVI_LATENT_KEY = "X_scanvi_refined"
SCANVI_LATENT_ALIAS = "X_scanvi"
SCVI_LATENT_KEY = "X_scvi"
UMAP_KEY = "X_umap"
UMAP_REFINED_KEY = "X_umap_refined"
UNLABELED_CATEGORY = "Unknown"

LABEL_RENAME_MAP = {
    "CD4 Naive/TCM": "CD4 Naive",
    "CD4 Tfr": "CD4 Tfh",
}

BATCH_KEY = "sample"
CANDIDATE_CAT_COVARIATES = ["tissue", "disease_status", "sex"]
CANDIDATE_CONT_COVARIATES = []

SCANVI_EPOCHS = 200
BATCH_SIZE = 256
N_NEIGHBORS = 30
DPI = 300
SEED = 42

# ============================================================================
# 1. Imports & setup
# ============================================================================
import gc
import json
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy import sparse

sc.settings.verbosity = 2
sc.settings.n_jobs = 16
sc.settings.set_figure_params(dpi=DPI, facecolor="white", frameon=False)

np.random.seed(SEED)
scvi.settings.seed = SEED

OUTPUT_PATH = Path(OUTPUT_DIR)
REPORT_DIR = OUTPUT_PATH / "reports"
FIG_DIR = OUTPUT_PATH / "figures"
for directory in (OUTPUT_PATH, REPORT_DIR, FIG_DIR):
    directory.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(FIG_DIR)

FINAL_H5AD_PATH = OUTPUT_PATH / "adata_tnk_scanvi_refined_relabel_rerun_20260414.h5ad"
SCANVI_MODEL_DIR = OUTPUT_PATH / "tnk_scanvi_refined_relabel_model"

print("=" * 80)
print("TNK scanvi-only rerun after refined label correction")
print("=" * 80)
print(f"input h5ad     : {INPUT_H5AD}")
print(f"scVI model dir : {SCVI_MODEL_DIR}")
print(f"output dir     : {OUTPUT_DIR}")
print(f"rename map     : {LABEL_RENAME_MAP}")
print(f"scanpy         : {sc.__version__}")
print(f"scvi-tools     : {scvi.__version__}")
print(f"torch          : {torch.__version__}")
print(f"cuda available : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    try:
        print(f"cuda device    : {torch.cuda.get_device_name(0)}")
    except Exception as exc:
        print(f"[WARN] Unable to query CUDA device name at startup: {exc}")


# ============================================================================
# 2. Helpers
# ============================================================================
def clean_string_series(series, fill_value=UNLABELED_CATEGORY, strip=True):
    series = series.astype(object)
    series = series.where(pd.notna(series), fill_value)
    series = series.astype(str)
    if strip:
        series = series.str.strip()
    series = series.replace({"": fill_value, "nan": fill_value, "None": fill_value})
    return series


def to_clean_category(series, fill_value=UNLABELED_CATEGORY, strip=True):
    return pd.Categorical(clean_string_series(series, fill_value=fill_value, strip=strip))


def save_table(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, index=False)


def write_json(payload, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def plot_loss(history_obj, title, out_path, train_key="train_loss_epoch", val_key=None):
    try:
        fig, ax = plt.subplots(figsize=(8, 4))
        train_hist = history_obj.get(train_key, None)
        if train_hist is not None:
            ax.plot(train_hist.values.flatten(), label="train")
        if val_key is not None:
            val_hist = history_obj.get(val_key, None)
            if val_hist is not None:
                ax.plot(val_hist.values.flatten(), label="validation")
        ax.set_xlabel("Epoch")
        ax.set_ylabel("ELBO / loss")
        ax.set_title(title)
        if len(ax.lines) > 0:
            ax.legend()
        fig.savefig(out_path, bbox_inches="tight", dpi=DPI)
        plt.close(fig)
        print(f"Saved: {out_path}")
    except Exception as exc:
        print(f"[WARN] Loss plot skipped for {title}: {exc}")


def summarize_label_counts(series, label_col=LABELS_KEY):
    return (
        pd.Series(clean_string_series(series), name=label_col)
        .value_counts(dropna=False)
        .rename_axis("label")
        .reset_index(name="n_cells")
        .sort_values(["n_cells", "label"], ascending=[False, True])
        .reset_index(drop=True)
    )


# ============================================================================
# 3. Load current filtered rerun object and rename supervision labels
# ============================================================================
print("\n[Step 1] Loading filtered TNK rerun H5AD...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  shape: {adata.shape}")

if LABELS_KEY not in adata.obs.columns:
    raise KeyError(f"Required label column '{LABELS_KEY}' not found")
if BATCH_KEY not in adata.obs.columns:
    raise KeyError(f"Required batch column '{BATCH_KEY}' not found")
if "counts" not in adata.layers:
    raise KeyError("Expected 'counts' layer for scVI/scANVI")

adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"], dtype=np.float32)
adata.X = sparse.csr_matrix(adata.X, dtype=np.float32)
if "log1p" in adata.layers:
    adata.layers["log1p"] = sparse.csr_matrix(adata.layers["log1p"], dtype=np.float32)

before_labels = clean_string_series(adata.obs[LABELS_KEY])
adata.obs[LABELS_BACKUP_KEY] = before_labels.copy()
after_labels = before_labels.replace(LABEL_RENAME_MAP)
adata.obs[LABELS_KEY] = to_clean_category(after_labels)
if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
    adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

before_counts = summarize_label_counts(before_labels, label_col="label_before")
after_counts = summarize_label_counts(adata.obs[LABELS_KEY], label_col="label_after")
rename_summary = (
    pd.DataFrame({
        "label_before": before_labels,
        "label_after": clean_string_series(adata.obs[LABELS_KEY]),
    })
    .groupby(["label_before", "label_after"], dropna=False)
    .size()
    .rename("n_cells")
    .reset_index()
    .sort_values(["label_before", "n_cells", "label_after"], ascending=[True, False, True])
)

save_table(before_counts, REPORT_DIR / "label_counts_before.tsv")
save_table(after_counts, REPORT_DIR / "label_counts_after.tsv")
save_table(rename_summary, REPORT_DIR / "label_rename_summary.tsv")

print("\nLabel counts before rename:")
print(before_counts.to_string(index=False))
print("\nLabel counts after rename:")
print(after_counts.to_string(index=False))

stored_params = adata.uns.get("rerun_params", {})
cat_covariates = stored_params.get("categorical_covariates", None)
cont_covariates = stored_params.get("continuous_covariates", None)
if cat_covariates is None:
    cat_covariates = [col for col in CANDIDATE_CAT_COVARIATES if col in adata.obs.columns]
else:
    cat_covariates = [col for col in cat_covariates if col in adata.obs.columns]
if cont_covariates is None:
    cont_covariates = [col for col in CANDIDATE_CONT_COVARIATES if col in adata.obs.columns]
else:
    cont_covariates = [col for col in cont_covariates if col in adata.obs.columns]

for col in cat_covariates + [BATCH_KEY]:
    adata.obs[col] = to_clean_category(adata.obs[col], fill_value="Unknown", strip=True)
for col in [LABELS_BACKUP_KEY, "scanvi_label", "scanvi_label_original", "cell_type_L2"]:
    if col in adata.obs.columns:
        adata.obs[col] = to_clean_category(adata.obs[col], fill_value="Unknown", strip=True)

labeled_only = after_counts[after_counts["label"] != UNLABELED_CATEGORY].copy()
if labeled_only.empty:
    raise RuntimeError(f"All cells became '{UNLABELED_CATEGORY}' after relabeling; cannot train scANVI")
if labeled_only.shape[0] < 2:
    raise RuntimeError("Fewer than 2 labeled classes remain; scANVI requires at least 2 supervised classes")
min_class_size = int(labeled_only["n_cells"].min())
n_samples_per_label = min(100, max(2, int(min_class_size * 0.8)))
n_samples_per_label = min(n_samples_per_label, min_class_size)
print(f"\nSmallest labeled class : {labeled_only.loc[labeled_only['n_cells'].idxmin(), 'label']} ({min_class_size} cells)")
print(f"Adaptive n_samples_per_label: {n_samples_per_label}")


# ============================================================================
# 4. Reload filtered scVI model and retrain scANVI only
# ============================================================================
print("\n[Step 2] Reloading filtered scVI model...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=BATCH_KEY,
    categorical_covariate_keys=cat_covariates or None,
    continuous_covariate_keys=cont_covariates or None,
)

accelerator = "gpu" if torch.cuda.is_available() else "cpu"
devices = 1
print(f"  accelerator: {accelerator}")
print(f"  devices    : {devices}")
print(f"  batch key  : {BATCH_KEY}")
print(f"  cat covars : {cat_covariates}")
print(f"  cont covars: {cont_covariates}")

vae = scvi.model.SCVI.load(SCVI_MODEL_DIR, adata=adata)
print(f"  scVI loaded from: {SCVI_MODEL_DIR}")
print(f"  n_latent        : {vae.module.n_latent}")

if SCVI_LATENT_KEY not in adata.obsm:
    adata.obsm[SCVI_LATENT_KEY] = vae.get_latent_representation()

print("\n[Step 3] Building and training scANVI on corrected labels...")
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    unlabeled_category=UNLABELED_CATEGORY,
    labels_key=LABELS_KEY,
)

lvae.train(
    max_epochs=SCANVI_EPOCHS,
    batch_size=BATCH_SIZE,
    train_size=0.9,
    early_stopping=True,
    early_stopping_patience=20,
    n_samples_per_label=n_samples_per_label,
    accelerator=accelerator,
    devices=devices,
)
print("  scANVI training complete")

adata.obsm[SCANVI_LATENT_KEY] = lvae.get_latent_representation()
adata.obsm[SCANVI_LATENT_ALIAS] = np.asarray(adata.obsm[SCANVI_LATENT_KEY]).copy()
adata.obs[PRED_KEY] = lvae.predict().astype(str)
scanvi_proba = lvae.predict(soft=True)
adata.obs[PRED_PROB_KEY] = scanvi_proba.max(axis=1).values.astype(np.float32)

concordance = (
    clean_string_series(adata.obs[PRED_KEY]) == clean_string_series(adata.obs[LABELS_KEY])
).mean()
print(f"  prediction/label concordance: {concordance:.4f}")
print(f"  confidence summary:\n{pd.Series(adata.obs[PRED_PROB_KEY]).describe()}")

plot_loss(
    lvae.history,
    title="TNK scanvi relabel rerun training loss",
    out_path=FIG_DIR / "scanvi_training_loss.pdf",
    train_key="train_loss_epoch",
)

lvae.save(str(SCANVI_MODEL_DIR), overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(SCANVI_MODEL_DIR / "var_names.csv", index=False, header=False)
print(f"  scANVI model saved: {SCANVI_MODEL_DIR}/")


# ============================================================================
# 5. Recompute neighbors + UMAP from new refined scANVI latent
# ============================================================================
print("\n[Step 4] Computing neighbors and UMAP from corrected scanvi latent...")
sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, n_neighbors=N_NEIGHBORS)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
adata.obsm[UMAP_REFINED_KEY] = np.asarray(adata.obsm[UMAP_KEY]).copy()
print(f"  {UMAP_REFINED_KEY}: {adata.obsm[UMAP_REFINED_KEY].shape}")


# ============================================================================
# 6. QC plots
# ============================================================================
print("\n[Step 5] Writing overview plots...")
sc.settings.vector_friendly = True
fig, axes = plt.subplots(2, 3, figsize=(24, 14))

sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=LABELS_BACKUP_KEY,
    title=f"Pre-relabel labels ({LABELS_BACKUP_KEY})",
    ax=axes[0, 0],
    show=False,
    legend_loc="right margin",
    legend_fontsize=6,
)
sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=LABELS_KEY,
    title=f"Corrected labels ({LABELS_KEY})",
    ax=axes[0, 1],
    show=False,
    legend_loc="right margin",
    legend_fontsize=6,
)
sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=PRED_KEY,
    title=f"scANVI predictions ({PRED_KEY})",
    ax=axes[0, 2],
    show=False,
    legend_loc="right margin",
    legend_fontsize=6,
)
sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=PRED_PROB_KEY,
    title="Prediction confidence",
    ax=axes[1, 0],
    show=False,
    color_map="RdYlGn",
    vmin=0,
    vmax=1,
)
sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=BATCH_KEY,
    title=f"Batch ({BATCH_KEY})",
    ax=axes[1, 1],
    show=False,
    legend_loc="right margin",
    legend_fontsize=5,
)
if "cell_type_L2" in adata.obs.columns:
    sc.pl.embedding(
        adata,
        basis="umap_refined",
        color="cell_type_L2",
        title="L2 lineage",
        ax=axes[1, 2],
        show=False,
        legend_loc="right margin",
        legend_fontsize=6,
    )
else:
    axes[1, 2].set_visible(False)

plt.suptitle("TNK corrected-label scanvi rerun overview", fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(FIG_DIR / "umap_scanvi_relabel_overview.pdf", bbox_inches="tight", dpi=DPI)
plt.close(fig)
print(f"Saved: {FIG_DIR / 'umap_scanvi_relabel_overview.pdf'}")


# ============================================================================
# 7. Save metadata + final H5AD
# ============================================================================
print("\n[Step 6] Saving final rerun H5AD...")
label_categories = sorted([c for c in adata.obs[LABELS_KEY].cat.categories if c != UNLABELED_CATEGORY])

adata.uns["scanvi_relabel_rerun_params"] = {
    "script_name": Path(__file__).name,
    "input_h5ad": INPUT_H5AD,
    "scvi_model_dir": SCVI_MODEL_DIR,
    "output_dir": str(OUTPUT_PATH),
    "final_h5ad": str(FINAL_H5AD_PATH),
    "scanvi_model_dir": str(SCANVI_MODEL_DIR),
    "labels_key": LABELS_KEY,
    "labels_backup_key": LABELS_BACKUP_KEY,
    "label_rename_map": LABEL_RENAME_MAP,
    "prediction_key": PRED_KEY,
    "prediction_prob_key": PRED_PROB_KEY,
    "batch_key": BATCH_KEY,
    "categorical_covariates": cat_covariates,
    "continuous_covariates": cont_covariates,
    "scanvi_epochs": int(SCANVI_EPOCHS),
    "batch_size": int(BATCH_SIZE),
    "n_samples_per_label": int(n_samples_per_label),
    "seed": int(SEED),
    "accelerator": accelerator,
    "devices": int(devices),
    "label_categories": label_categories,
    "prediction_label_concordance": float(concordance),
}

summary_json = {
    "input_n_obs": int(adata.n_obs),
    "n_vars": int(adata.n_vars),
    "label_rename_map": LABEL_RENAME_MAP,
    "label_categories": label_categories,
    "min_class_size": int(min_class_size),
    "n_samples_per_label": int(n_samples_per_label),
    "prediction_label_concordance": float(concordance),
}
write_json(summary_json, REPORT_DIR / "scanvi_relabel_summary.json")

for col in [LABELS_KEY, LABELS_BACKUP_KEY, PRED_KEY, BATCH_KEY, "scanvi_label", "scanvi_label_original", "cell_type_L2"]:
    if col in adata.obs.columns:
        adata.obs[col] = to_clean_category(adata.obs[col], fill_value="Unknown", strip=True)

import anndata as _anndata
_anndata.settings.allow_write_nullable_strings = True

for attr in ("obs", "var"):
    df = getattr(adata, attr)
    if "_index" in df.columns:
        setattr(adata, attr, df.rename(columns={"_index": "orig_index"}))

adata.write_h5ad(FINAL_H5AD_PATH, compression="gzip", compression_opts=9)
print(f"Saved: {FINAL_H5AD_PATH}")


del lvae, vae
gc.collect()
if torch.cuda.is_available():
    torch.cuda.empty_cache()

print("\n" + "=" * 80)
print("TNK scanvi relabel rerun complete")
print("=" * 80)
print(f"Final H5AD   : {FINAL_H5AD_PATH}")
print(f"scANVI model : {SCANVI_MODEL_DIR}/")
print(f"Reports      : {REPORT_DIR}/")
print(f"Figures      : {FIG_DIR}/")
