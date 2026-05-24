#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T/NK refined-label scVI/scANVI rerun after removing selected CHOIR clusters.

Purpose
-------
1. Read the current TNK tissue-comparison final H5AD and identify cells in the
   user-selected CHOIR clusters.
2. Remove those cells from the current refined-label TNK scANVI reference H5AD.
3. Re-train scVI from scratch on the filtered object.
4. Re-train scANVI using the existing curated `scanvi_label_refined` labels.
5. Save a fresh H5AD compatible with the TNK tissue-comparison wrapper.

Rationale
---------
The user explicitly asked to remove CHOIR clusters c10, c20, c22, c29, c47,
c2, c17 and then rerun scVI/scANVI "according to the previous labels".
To preserve the validated refined label system while avoiding another round of
subcluster-specific manual remapping, this rerun uses the existing
`scanvi_label_refined` column directly as the supervision label.
"""

# ============================================================================
# 0. Configuration
# ============================================================================
CLUSTER_SOURCE_H5AD = (
    "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/"
    "tnk_tissue_comparison_final.h5ad"
)
INPUT_H5AD = \
    "/home/h2048/data/py/0318/tnk_subcluster_retrain/adata_tnk_scanvi_ref_retrain_v1_2.h5ad"
OUTPUT_DIR = \
    "/home/h2048/data/py/0413/tnk_scvi_scanvi_refined_rerun_rm_choir_20260413"

CHOIR_CLUSTER_COL = "CHOIR_clusters_0.2"
REMOVE_CLUSTERS = ["c10", "c20", "c22", "c29", "c47", "c2", "c17"]

LABELS_KEY = "scanvi_label_refined"
LABELS_BACKUP_KEY = "scanvi_label_refined_pre_rerun"
PRED_KEY = "scanvi_pred_refined"
PRED_PROB_KEY = "scanvi_pred_prob_refined"
SCANVI_LATENT_KEY = "X_scanvi_refined"
SCANVI_LATENT_ALIAS = "X_scanvi"
SCVI_LATENT_KEY = "X_scvi"
UMAP_KEY = "X_umap"
UMAP_REFINED_KEY = "X_umap_refined"
UNLABELED_CATEGORY = "Unknown"

BATCH_KEY = "sample"
CANDIDATE_CAT_COVARIATES = ["tissue", "disease_status", "sex"]
CANDIDATE_CONT_COVARIATES = []

N_LATENT = 75
N_HIDDEN = 128
N_LAYERS = 2
DROPOUT = 0.1
SCVI_EPOCHS = 400
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

FINAL_H5AD_PATH = OUTPUT_PATH / "adata_tnk_scanvi_refined_rerun_rm_choir_20260413.h5ad"
SCVI_MODEL_DIR = OUTPUT_PATH / "tnk_scvi_refined_rerun_model"
SCANVI_MODEL_DIR = OUTPUT_PATH / "tnk_scanvi_refined_rerun_model"

print("=" * 80)
print("TNK refined-label scVI/scANVI rerun after CHOIR-cluster removal")
print("=" * 80)
print(f"cluster source : {CLUSTER_SOURCE_H5AD}")
print(f"input h5ad     : {INPUT_H5AD}")
print(f"output dir     : {OUTPUT_DIR}")
print(f"remove clusters: {REMOVE_CLUSTERS}")
print(f"scanpy         : {sc.__version__}")
print(f"scvi-tools     : {scvi.__version__}")
print(f"torch          : {torch.__version__}")
print(f"cuda available : {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"cuda device    : {torch.cuda.get_device_name(0)}")


# ============================================================================
# 2. Helpers
# ============================================================================
def normalize_cluster_id(value):
    value = str(value).strip()
    if value.lower().startswith("c"):
        value = value[1:]
    return value


def sort_cluster_ids(values):
    def _key(x):
        x = normalize_cluster_id(x)
        return (int(x) if x.isdigit() else x)
    return sorted([normalize_cluster_id(x) for x in values], key=_key)


REMOVE_CLUSTERS_NORM = sort_cluster_ids(REMOVE_CLUSTERS)


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


def write_json(payload, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


# ============================================================================
# 3. Identify cells to remove from current tissue-comparison output
# ============================================================================
print("\n[Step 1] Reading current tissue-comparison final H5AD for CHOIR removal set...")
cluster_adata = ad.read_h5ad(CLUSTER_SOURCE_H5AD, backed="r")
cluster_cols = [
    CHOIR_CLUSTER_COL,
    LABELS_KEY,
    "scanvi_label",
    "scanvi_label_original",
    "cell_type_L2",
]
cluster_cols = [col for col in cluster_cols if col in cluster_adata.obs.columns]
cluster_obs = cluster_adata.obs[cluster_cols].copy()
cluster_obs[CHOIR_CLUSTER_COL] = cluster_obs[CHOIR_CLUSTER_COL].astype(str).map(normalize_cluster_id)
remove_mask = cluster_obs[CHOIR_CLUSTER_COL].isin(REMOVE_CLUSTERS_NORM)
remove_cells = pd.Index(cluster_obs.index[remove_mask])

if len(remove_cells) == 0:
    raise RuntimeError(
        f"No cells matched REMOVE_CLUSTERS={REMOVE_CLUSTERS}. "
        f"Check {CHOIR_CLUSTER_COL} encoding in {CLUSTER_SOURCE_H5AD}."
    )

print(f"  clusters matched : {REMOVE_CLUSTERS_NORM}")
print(f"  cells to remove  : {len(remove_cells):,}")

cluster_size_df = (
    cluster_obs.loc[remove_mask, [CHOIR_CLUSTER_COL]]
    .value_counts()
    .rename("n_cells")
    .reset_index()
    .sort_values(by=CHOIR_CLUSTER_COL, key=lambda s: s.astype(int))
)
cluster_size_df.insert(0, "cluster_id_prefixed", [f"c{x}" for x in cluster_size_df[CHOIR_CLUSTER_COL]])
save_table(cluster_size_df, REPORT_DIR / "removed_cluster_sizes.tsv")

for column in [col for col in [LABELS_KEY, "scanvi_label", "scanvi_label_original", "cell_type_L2"] if col in cluster_obs.columns]:
    tmp = cluster_obs.loc[remove_mask, [CHOIR_CLUSTER_COL, column]].copy()
    tmp[column] = clean_string_series(tmp[column])
    summary = (
        tmp.groupby([CHOIR_CLUSTER_COL, column], dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values([CHOIR_CLUSTER_COL, "n_cells", column], ascending=[True, False, True])
    )
    save_table(summary, REPORT_DIR / f"removed_cluster_by_{column}.tsv")

removed_cells_df = pd.DataFrame({"cell_barcode": remove_cells, "removed_choir_cluster": cluster_obs.loc[remove_cells, CHOIR_CLUSTER_COL].values})
save_table(removed_cells_df, REPORT_DIR / "removed_cells.tsv")
cluster_adata.file.close()


# ============================================================================
# 4. Load refined input H5AD and filter those cells
# ============================================================================
print("\n[Step 2] Loading refined TNK reference H5AD and filtering cells...")
adata = sc.read_h5ad(INPUT_H5AD)
print(f"  original shape : {adata.shape}")

missing_cells = remove_cells.difference(adata.obs_names)
if len(missing_cells) > 0:
    raise RuntimeError(
        f"{len(missing_cells)} cells selected for removal are missing from INPUT_H5AD. "
        f"First few: {missing_cells[:10].tolist()}"
    )

keep_mask = ~adata.obs_names.isin(remove_cells)
adata = adata[keep_mask].copy()
print(f"  filtered shape : {adata.shape}")
print(f"  cells removed  : {(~keep_mask).sum():,}")

if LABELS_KEY not in adata.obs.columns:
    raise KeyError(f"Required label column '{LABELS_KEY}' not found in filtered AnnData")
if BATCH_KEY not in adata.obs.columns:
    raise KeyError(f"Required batch column '{BATCH_KEY}' not found in filtered AnnData")
if "counts" not in adata.layers:
    raise KeyError("Filtered AnnData does not contain a 'counts' layer required for scVI training")

adata.layers["counts"] = sparse.csr_matrix(adata.layers["counts"], dtype=np.float32)
adata.X = sparse.csr_matrix(adata.X, dtype=np.float32)
if "log1p" in adata.layers:
    adata.layers["log1p"] = sparse.csr_matrix(adata.layers["log1p"], dtype=np.float32)

adata.obs[LABELS_BACKUP_KEY] = clean_string_series(adata.obs[LABELS_KEY])
adata.obs[LABELS_KEY] = to_clean_category(adata.obs[LABELS_KEY])
if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
    adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

cat_covariates = [col for col in CANDIDATE_CAT_COVARIATES if col in adata.obs.columns]
cont_covariates = [col for col in CANDIDATE_CONT_COVARIATES if col in adata.obs.columns]
for col in cat_covariates + [BATCH_KEY]:
    adata.obs[col] = to_clean_category(adata.obs[col], fill_value="Unknown", strip=True)

for col in ["scanvi_label", "scanvi_label_original", "cell_type_L2", LABELS_BACKUP_KEY]:
    if col in adata.obs.columns:
        adata.obs[col] = to_clean_category(adata.obs[col], fill_value="Unknown", strip=True)

label_counts = (
    pd.Series(clean_string_series(adata.obs[LABELS_KEY]), name=LABELS_KEY)
    .value_counts()
    .rename_axis("label")
    .reset_index(name="n_cells")
)
save_table(label_counts, REPORT_DIR / "remaining_refined_label_counts.tsv")
print("\nRemaining refined label distribution:")
print(label_counts.to_string(index=False))

labeled_only = label_counts[label_counts["label"] != UNLABELED_CATEGORY].copy()
if labeled_only.empty:
    raise RuntimeError(f"All cells became '{UNLABELED_CATEGORY}' after filtering; cannot train scANVI")
if labeled_only.shape[0] < 2:
    raise RuntimeError("Fewer than 2 labeled classes remain; scANVI requires at least 2 supervised classes")
min_class_size = int(labeled_only["n_cells"].min())
n_samples_per_label = min(100, max(2, int(min_class_size * 0.8)))
n_samples_per_label = min(n_samples_per_label, min_class_size)
print(f"\nSmallest labeled class : {labeled_only.loc[labeled_only['n_cells'].idxmin(), 'label']} ({min_class_size} cells)")
print(f"Adaptive n_samples_per_label: {n_samples_per_label}")

save_table(
    pd.DataFrame(
        {
            "metric": [
                "input_n_obs",
                "filtered_n_obs",
                "n_removed_cells",
                "n_removed_clusters",
                "n_labeled_classes",
                "min_class_size",
                "n_samples_per_label",
            ],
            "value": [
                int(adata.n_obs + len(remove_cells)),
                int(adata.n_obs),
                int(len(remove_cells)),
                int(len(REMOVE_CLUSTERS_NORM)),
                int(labeled_only.shape[0]),
                int(min_class_size),
                int(n_samples_per_label),
            ],
        }
    ),
    REPORT_DIR / "filter_and_label_summary.tsv",
)


# ============================================================================
# 5. Train scVI from scratch on filtered object
# ============================================================================
print("\n[Step 3] Setting up and training scVI...")
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

vae = scvi.model.SCVI(
    adata,
    n_latent=N_LATENT,
    n_hidden=N_HIDDEN,
    n_layers=N_LAYERS,
    dropout_rate=DROPOUT,
    gene_likelihood="nb",
    encode_covariates=True,
    deeply_inject_covariates=False,
    use_layer_norm="both",
    use_batch_norm="none",
)

vae.train(
    max_epochs=SCVI_EPOCHS,
    batch_size=BATCH_SIZE,
    train_size=0.9,
    early_stopping=True,
    early_stopping_patience=30,
    plan_kwargs={"lr": 1e-3},
    accelerator=accelerator,
    devices=devices,
)
print("  scVI training complete")

adata.obsm[SCVI_LATENT_KEY] = vae.get_latent_representation()
print(f"  {SCVI_LATENT_KEY}: {adata.obsm[SCVI_LATENT_KEY].shape}")
plot_loss(
    vae.history,
    title="TNK scVI rerun training loss",
    out_path=FIG_DIR / "scvi_training_loss.pdf",
    train_key="train_loss_epoch",
    val_key="elbo_validation",
)

vae.save(str(SCVI_MODEL_DIR), overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(SCVI_MODEL_DIR / "var_names.csv", index=False, header=False)
print(f"  scVI model saved : {SCVI_MODEL_DIR}/")


# ============================================================================
# 6. Train scANVI using previous refined labels as supervision
# ============================================================================
print("\n[Step 4] Building and training scANVI from the new scVI model...")
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
print(f"  {SCANVI_LATENT_KEY}: {adata.obsm[SCANVI_LATENT_KEY].shape}")

concordance = (
    clean_string_series(adata.obs[PRED_KEY]) == clean_string_series(adata.obs[LABELS_KEY])
).mean()
print(f"  prediction/label concordance: {concordance:.4f}")
print(f"  confidence summary:\n{pd.Series(adata.obs[PRED_PROB_KEY]).describe()}")

plot_loss(
    lvae.history,
    title="TNK scANVI rerun training loss",
    out_path=FIG_DIR / "scanvi_training_loss.pdf",
    train_key="train_loss_epoch",
)

lvae.save(str(SCANVI_MODEL_DIR), overwrite=True)
pd.Series(adata.var_names.tolist()).to_csv(SCANVI_MODEL_DIR / "var_names.csv", index=False, header=False)
print(f"  scANVI model saved: {SCANVI_MODEL_DIR}/")


# ============================================================================
# 7. Recompute neighbors + UMAP from new refined latent space
# ============================================================================
print("\n[Step 5] Computing neighbors and UMAP from new refined scANVI latent...")
umap_bundle = fit_bundle(
    adata,
    SCANVI_LATENT_KEY,
    OUTPUT_PATH,
    primary_umap_key=UMAP_REFINED_KEY,
    operator_filename="tnk_umap_refined_operator.joblib",
    manifest_filename="tnk_umap_refined_bundle.json",
    umap_params={
        "n_neighbors": N_NEIGHBORS,
        "n_components": 2,
        "min_dist": 0.3,
        "spread": 1.0,
        "metric": "euclidean",
        "random_state": SEED,
    },
    set_default_x_umap=True,
    extra_manifest={"removed_choir_clusters": REMOVE_CLUSTERS},
)
print(f"  {UMAP_REFINED_KEY}: {adata.obsm[UMAP_REFINED_KEY].shape}")
print(f"  scanVI UMAP operator: {umap_bundle['operator_path']}")


# ============================================================================
# 8. QC plots
# ============================================================================
print("\n[Step 6] Writing overview plots...")
sc.settings.vector_friendly = True
fig, axes = plt.subplots(2, 3, figsize=(24, 14))

sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=LABELS_KEY,
    title=f"Supervision labels ({LABELS_KEY})",
    ax=axes[0, 0],
    show=False,
    legend_loc="right margin",
    legend_fontsize=6,
)
sc.pl.embedding(
    adata,
    basis="umap_refined",
    color=PRED_KEY,
    title=f"scANVI predictions ({PRED_KEY})",
    ax=axes[0, 1],
    show=False,
    legend_loc="right margin",
    legend_fontsize=6,
)
if LABELS_BACKUP_KEY in adata.obs.columns:
    sc.pl.embedding(
        adata,
        basis="umap_refined",
        color=LABELS_BACKUP_KEY,
        title=f"Pre-rerun labels ({LABELS_BACKUP_KEY})",
        ax=axes[0, 2],
        show=False,
        legend_loc="right margin",
        legend_fontsize=6,
    )
else:
    axes[0, 2].set_visible(False)

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

plt.suptitle("TNK refined-label rerun overview", fontsize=14, y=1.01)
plt.tight_layout()
fig.savefig(FIG_DIR / "umap_refined_rerun_overview.pdf", bbox_inches="tight", dpi=DPI)
plt.close(fig)
print(f"Saved: {FIG_DIR / 'umap_refined_rerun_overview.pdf'}")


# ============================================================================
# 9. Save metadata + final H5AD
# ============================================================================
print("\n[Step 7] Saving final rerun H5AD...")
label_categories = sorted([c for c in adata.obs[LABELS_KEY].cat.categories if c != UNLABELED_CATEGORY])

adata.uns["rerun_params"] = {
    "script_name": Path(__file__).name,
    "cluster_source_h5ad": CLUSTER_SOURCE_H5AD,
    "input_h5ad": INPUT_H5AD,
    "output_dir": str(OUTPUT_PATH),
    "final_h5ad": str(FINAL_H5AD_PATH),
    "scvi_model_dir": str(SCVI_MODEL_DIR),
    "scanvi_model_dir": str(SCANVI_MODEL_DIR),
    "remove_clusters_user": REMOVE_CLUSTERS,
    "remove_clusters_normalized": REMOVE_CLUSTERS_NORM,
    "removed_cell_count": int(len(remove_cells)),
    "labels_key": LABELS_KEY,
    "labels_backup_key": LABELS_BACKUP_KEY,
    "prediction_key": PRED_KEY,
    "prediction_prob_key": PRED_PROB_KEY,
    "batch_key": BATCH_KEY,
    "categorical_covariates": cat_covariates,
    "continuous_covariates": cont_covariates,
    "n_latent": int(N_LATENT),
    "n_hidden": int(N_HIDDEN),
    "n_layers": int(N_LAYERS),
    "dropout": float(DROPOUT),
    "scvi_epochs": int(SCVI_EPOCHS),
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
    "removed_clusters": REMOVE_CLUSTERS,
    "removed_clusters_normalized": REMOVE_CLUSTERS_NORM,
    "removed_cell_count": int(len(remove_cells)),
    "filtered_n_obs": int(adata.n_obs),
    "n_vars": int(adata.n_vars),
    "label_categories": label_categories,
    "min_class_size": int(min_class_size),
    "n_samples_per_label": int(n_samples_per_label),
    "prediction_label_concordance": float(concordance),
}
write_json(summary_json, REPORT_DIR / "rerun_summary.json")

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
print("TNK refined-label rerun complete")
print("=" * 80)
print(f"Final H5AD    : {FINAL_H5AD_PATH}")
print(f"scVI model    : {SCVI_MODEL_DIR}/")
print(f"scANVI model  : {SCANVI_MODEL_DIR}/")
print(f"Reports       : {REPORT_DIR}/")
print(f"Figures       : {FIG_DIR}/")
