#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stromal branch-specific rerun after CHOIR cleanup / relabel.

Branches
--------
- endothelial: remove CHOIR clusters 5, 33 and retrain scVI + scANVI
- fibroblast : remove CHOIR clusters 23, 59 and retrain scVI + scANVI
- smc        : relabel CHOIR cluster 6 to Peripheral_neuron_like and rerun scANVI

This script is intentionally import-safe so lightweight helper tests can import
it without pulling heavy scanpy / scvi dependencies at import time.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import site
import sys
import types
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

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

UNLABELED_CATEGORY = "Unknown"
PIPELINE_DATE = "20260414"
TARGET_SMC_LABEL = "Peripheral_neuron_like"
OUTPUT_ROOT = Path(f"/home/h2048/data/py/0414/stromal_branch_rerun_{PIPELINE_DATE}")

BATCH_KEY = "sample"
CHOIR_CLUSTER_COL = "CHOIR_clusters_0.2"
LABELS_KEY = "scanvi_label"
L2_KEY = "cell_type_L2"
L3_KEY = "cell_type_L3"
PRED_KEY = "cell_type_scanvi_pred"
PRED_PROB_KEY = "scanvi_pred_prob_rerun"
SCVI_LATENT_KEY = "X_scvi"
SCANVI_LATENT_KEY = "X_scanvi"
UMAP_KEY = "X_umap"
UMAP_SCANVI_KEY = "X_umap_scanvi"
UMAP_SCANVI_CORRECTED_KEY = "X_umap_scanvi_corrected"
DEFAULT_CONTINUOUS_COVARIATES = ("pct_counts_mt", "stress_score", "S_score", "G2M_score")

SCVI_N_LATENT = 75
SCVI_N_HIDDEN = 128
SCVI_N_LAYERS = 2
SCVI_DROPOUT = 0.2
SCVI_EPOCHS = 400
SCANVI_EPOCHS = 200
BATCH_SIZE = 256
SCVI_PATIENCE = 45
SCANVI_PATIENCE = 30
SCVI_LR = 1e-3
SCANVI_LR = 5e-4
N_NEIGHBORS = 30
DPI = 300
SEED = 42


@dataclass(frozen=True)
class BranchSpec:
    name: str
    mode: str  # remove or relabel
    cluster_source_h5ad: str
    input_h5ad: str
    output_dir: str
    target_clusters: tuple[str, ...]
    final_h5ad_name: str
    scvi_model_name: str
    scanvi_model_name: str
    label_backup_suffix: str
    existing_scvi_model_dir: str | None = None
    relabel_map: dict[str, str] | None = None


BRANCH_SPECS: dict[str, BranchSpec] = {
    "endothelial": BranchSpec(
        name="endothelial",
        mode="remove",
        cluster_source_h5ad=(
            "/home/h2048/data/R/0414/stromal_endothelial_tissue_comparison_v1_1_0_20260414/"
            "stromal_endothelial_tissue_comparison_final.h5ad"
        ),
        input_h5ad=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/endothelial/"
            "adata_endothelial_reference_v1_5_branchwise.h5ad"
        ),
        output_dir=str(OUTPUT_ROOT / "endothelial"),
        target_clusters=("5", "33"),
        final_h5ad_name="adata_endothelial_reference_rm_choir_5_33_20260414.h5ad",
        scvi_model_name="scvi_endothelial_rm_choir_20260414",
        scanvi_model_name="scanvi_endothelial_rm_choir_20260414",
        label_backup_suffix="pre_rm_choir_20260414",
        existing_scvi_model_dir=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/endothelial/models/"
            "scvi_endothelial_v1_5_branchwise"
        ),
    ),
    "fibroblast": BranchSpec(
        name="fibroblast",
        mode="remove",
        cluster_source_h5ad=(
            "/home/h2048/data/R/0414/stromal_fibroblast_tissue_comparison_v1_1_0_20260414/"
            "stromal_fibroblast_tissue_comparison_final.h5ad"
        ),
        input_h5ad=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/fibroblast/"
            "adata_fibroblast_reference_v1_5_branchwise.h5ad"
        ),
        output_dir=str(OUTPUT_ROOT / "fibroblast"),
        target_clusters=("23", "59"),
        final_h5ad_name="adata_fibroblast_reference_rm_choir_23_59_20260414.h5ad",
        scvi_model_name="scvi_fibroblast_rm_choir_20260414",
        scanvi_model_name="scanvi_fibroblast_rm_choir_20260414",
        label_backup_suffix="pre_rm_choir_20260414",
        existing_scvi_model_dir=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/fibroblast/models/"
            "scvi_fibroblast_v1_5_branchwise"
        ),
    ),
    "smc": BranchSpec(
        name="smc",
        mode="relabel",
        cluster_source_h5ad=(
            "/home/h2048/data/R/0414/stromal_smc_tissue_comparison_v1_1_0_20260414/"
            "stromal_smc_tissue_comparison_final.h5ad"
        ),
        input_h5ad=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/smc/"
            "adata_smc_reference_v1_5_branchwise.h5ad"
        ),
        output_dir=str(OUTPUT_ROOT / "smc"),
        target_clusters=("6",),
        final_h5ad_name="adata_smc_reference_neuronlike_c6_20260414.h5ad",
        scvi_model_name="scvi_smc_reused_20260414",
        scanvi_model_name="scanvi_smc_neuronlike_c6_20260414",
        label_backup_suffix="pre_neuronlike_20260414",
        existing_scvi_model_dir=(
            "/home/h2048/data/py/0407/stromal_reintegration_v1_5_branchwise/smc/models/"
            "scvi_smc_v1_5_branchwise"
        ),
        relabel_map={"6": TARGET_SMC_LABEL},
    ),
}


def normalize_cluster_id(value: object) -> str:
    value = str(value).strip()
    if value.lower().startswith("c"):
        value = value[1:]
    return value


def normalize_cluster_set(values: Iterable[object]) -> set[str]:
    return {normalize_cluster_id(v) for v in values}


def select_target_cells(
    obs_names: Sequence[object],
    choir_values: Sequence[object],
    target_clusters: Iterable[object],
) -> list[str]:
    targets = normalize_cluster_set(target_clusters)
    return [
        str(cell)
        for cell, choir in zip(obs_names, choir_values)
        if normalize_cluster_id(choir) in targets
    ]


def relabel_by_cluster(
    obs_names: Sequence[object],
    choir_values: Sequence[object],
    labels: Sequence[object],
    relabel_map: Mapping[object, str],
) -> list[str]:
    norm_map = {normalize_cluster_id(k): v for k, v in relabel_map.items()}
    relabeled: list[str] = []
    for _cell, choir, label in zip(obs_names, choir_values, labels):
        relabeled.append(norm_map.get(normalize_cluster_id(choir), str(label)))
    return relabeled


def extract_scvi_init_kwargs(
    init_params: Mapping[str, object],
    module_kwargs: Mapping[str, object] | None = None,
) -> dict[str, object]:
    non_kwargs = init_params.get("non_kwargs", {})
    nested_kwargs = init_params.get("kwargs", {})

    flat_kwargs: dict[str, object] = {}
    if isinstance(module_kwargs, Mapping):
        flat_kwargs.update(module_kwargs)

    if isinstance(nested_kwargs, Mapping):
        for group in nested_kwargs.values():
            if isinstance(group, Mapping):
                for key, value in group.items():
                    flat_kwargs.setdefault(key, value)

    merged: dict[str, object] = {}
    if isinstance(non_kwargs, Mapping):
        for key in ("use_observed_lib_size",):
            if key in non_kwargs:
                merged[key] = non_kwargs[key]
    merged.update(flat_kwargs)
    return merged


def configure_runtime() -> None:
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(key, "8")
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    user_site_candidates = []
    try:
        user_site_candidates.append(site.getusersitepackages())
    except Exception:
        pass
    user_site_candidates.extend(
        [
            os.path.expanduser("~/.local/lib/python3.10/site-packages"),
            os.path.expanduser("~/.local/lib/python3.11/site-packages"),
            os.path.expanduser("~/.local/lib/python3.9/site-packages"),
        ]
    )
    user_site_candidates = [p for p in user_site_candidates if isinstance(p, str) and p]
    sys.path[:] = [
        p for p in sys.path
        if not any(str(p).startswith(candidate) for candidate in user_site_candidates)
    ]


def import_runtime_modules():
    configure_runtime()
    import warnings

    warnings.filterwarnings("ignore")

    import anndata as ad  # type: ignore
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    import numpy as np  # type: ignore
    import pandas as pd  # type: ignore
    import scanpy as sc  # type: ignore
    import scvi  # type: ignore
    import torch  # type: ignore
    from scipy import sparse  # type: ignore

    sc.settings.verbosity = 2
    sc.settings.n_jobs = 16
    sc.settings.set_figure_params(dpi=DPI, facecolor="white", frameon=False)
    np.random.seed(SEED)
    scvi.settings.seed = SEED
    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)

    return ad, plt, np, pd, sc, scvi, torch, sparse


def clean_string_series(series, pd_module, fill_value: str = UNLABELED_CATEGORY, strip: bool = True):
    series = series.astype(object)
    series = series.where(pd_module.notna(series), fill_value)
    series = series.astype(str)
    if strip:
        series = series.str.strip()
    return series.replace({"": fill_value, "nan": fill_value, "None": fill_value})


def to_clean_category(series, pd_module, fill_value: str = UNLABELED_CATEGORY, strip: bool = True):
    return pd_module.Categorical(clean_string_series(series, pd_module, fill_value=fill_value, strip=strip))


def save_table(df, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".tsv":
        df.to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, index=False)


def write_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)


def get_output_paths(spec: BranchSpec) -> dict[str, Path]:
    root = Path(spec.output_dir)
    return {
        "root": root,
        "reports": root / "reports",
        "figures": root / "figures",
        "final_h5ad": root / spec.final_h5ad_name,
        "scvi_model": root / spec.scvi_model_name,
        "scanvi_model": root / spec.scanvi_model_name,
    }


def summarize_target_clusters(spec: BranchSpec, pd_module):
    ad, _, _, _, _, _, _, _ = import_runtime_modules()
    cluster_adata = ad.read_h5ad(spec.cluster_source_h5ad, backed="r")
    obs = cluster_adata.obs[[c for c in [CHOIR_CLUSTER_COL, L2_KEY, L3_KEY, PRED_KEY, LABELS_KEY, "source_cluster_id"] if c in cluster_adata.obs.columns]].copy()
    obs[CHOIR_CLUSTER_COL] = obs[CHOIR_CLUSTER_COL].astype(str).map(normalize_cluster_id)
    target_cells = pd_module.Index(select_target_cells(obs.index.tolist(), obs[CHOIR_CLUSTER_COL].tolist(), spec.target_clusters))
    target_obs = obs.loc[target_cells].copy()
    cluster_adata.file.close()

    if len(target_cells) == 0:
        raise RuntimeError(
            f"No cells matched target clusters {spec.target_clusters} in {spec.cluster_source_h5ad}"
        )

    return target_cells, target_obs


def ensure_training_ready(adata, pd_module, sparse_module) -> tuple[list[str], int, int]:
    if "counts" not in adata.layers:
        raise KeyError("Expected `counts` layer in branch reference AnnData")
    adata.layers["counts"] = sparse_module.csr_matrix(adata.layers["counts"], dtype="float32")
    if "log1p" in adata.layers:
        adata.layers["log1p"] = sparse_module.csr_matrix(adata.layers["log1p"], dtype="float32")
        adata.X = adata.layers["log1p"]
    else:
        adata.X = sparse_module.csr_matrix(adata.X, dtype="float32")
    if not sparse_module.issparse(adata.X):
        adata.X = sparse_module.csr_matrix(adata.X, dtype="float32")

    for col in [BATCH_KEY, LABELS_KEY, L2_KEY, L3_KEY, PRED_KEY]:
        if col in adata.obs.columns:
            fill = "Unknown" if col == BATCH_KEY else UNLABELED_CATEGORY
            adata.obs[col] = to_clean_category(adata.obs[col], pd_module, fill_value=fill)
    if LABELS_KEY not in adata.obs.columns:
        raise KeyError(f"Missing required label column: {LABELS_KEY}")
    if BATCH_KEY not in adata.obs.columns:
        raise KeyError(f"Missing required batch column: {BATCH_KEY}")
    if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
        adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

    labeled_counts = (
        pd_module.Series(clean_string_series(adata.obs[LABELS_KEY], pd_module), name=LABELS_KEY)
        .value_counts()
        .rename_axis("label")
        .reset_index(name="n_cells")
    )
    labeled_only = labeled_counts[labeled_counts["label"] != UNLABELED_CATEGORY].copy()
    if labeled_only.empty:
        raise RuntimeError("All supervision labels are Unknown; cannot train scANVI")
    if labeled_only.shape[0] < 2:
        raise RuntimeError("Fewer than 2 labeled classes remain; scANVI requires at least 2 classes")
    min_class_size = int(labeled_only["n_cells"].min())
    n_samples_per_label = min(100, max(2, int(min_class_size * 0.8)))
    n_samples_per_label = min(n_samples_per_label, min_class_size)
    return labeled_counts["label"].tolist(), min_class_size, n_samples_per_label


def get_continuous_covariates(adata) -> list[str]:
    uns_covs = []
    if isinstance(adata.uns.get("branchwise_reference"), dict):
        uns_covs = list(adata.uns["branchwise_reference"].get("continuous_covariates", []))
    if not uns_covs:
        uns_covs = list(DEFAULT_CONTINUOUS_COVARIATES)
    return [cov for cov in uns_covs if cov in adata.obs.columns]


def disable_broken_cuda(torch_module) -> None:
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    try:
        torch_module.cuda.is_available = lambda: False  # type: ignore[assignment]
        torch_module.cuda.device_count = lambda: 0  # type: ignore[assignment]
        torch_module.cuda.get_rng_state_all = lambda: []  # type: ignore[assignment]
        torch_module.cuda.manual_seed_all = lambda *_args, **_kwargs: None  # type: ignore[assignment]
        torch_module.cuda.empty_cache = lambda: None  # type: ignore[assignment]
    except Exception:
        pass


def inject_pandas_pickle_compat(pd_module) -> None:
    if "pandas.core.indexes.numeric" in sys.modules:
        return
    compat_mod = types.ModuleType("pandas.core.indexes.numeric")
    compat_mod.Int64Index = pd_module.Index
    compat_mod.UInt64Index = pd_module.Index
    compat_mod.Float64Index = pd_module.Index
    sys.modules["pandas.core.indexes.numeric"] = compat_mod


def resolve_accelerator(torch_module) -> tuple[str, int | str]:
    if not torch_module.cuda.is_available():
        return "cpu", 1
    try:
        _ = torch_module.cuda.get_device_properties(0)
        return "gpu", 0
    except Exception as exc:
        print(f"[WARN] CUDA reported as available but is unusable; falling back to CPU. Error: {exc}")
        disable_broken_cuda(torch_module)
        return "cpu", 1


def train_scvi_from_scratch(adata, scvi_module, torch_module, continuous_covariates: list[str], model_dir: Path):
    scvi_module.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=BATCH_KEY,
        continuous_covariate_keys=continuous_covariates or None,
    )
    model = scvi_module.model.SCVI(
        adata,
        n_latent=SCVI_N_LATENT,
        n_hidden=SCVI_N_HIDDEN,
        n_layers=SCVI_N_LAYERS,
        dropout_rate=SCVI_DROPOUT,
        gene_likelihood="nb",
        dispersion="gene-batch",
        encode_covariates=True,
        use_layer_norm="both",
        use_batch_norm="none",
    )
    accelerator, _device = resolve_accelerator(torch_module)
    model.train(
        max_epochs=SCVI_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=True,
        early_stopping_patience=SCVI_PATIENCE,
        train_size=0.9,
        accelerator=accelerator,
        devices=1,
        plan_kwargs={"lr": SCVI_LR},
    )
    model.save(str(model_dir), overwrite=True)
    return model, accelerator


def load_existing_scvi(adata, scvi_module, torch_module, continuous_covariates: list[str], model_dir: str | Path):
    import pandas as pd  # local import to keep module import-safe for tests

    def register_current_anndata() -> None:
        scvi_module.model.SCVI.setup_anndata(
            adata,
            layer="counts",
            batch_key=BATCH_KEY,
            continuous_covariate_keys=continuous_covariates or None,
        )

    inject_pandas_pickle_compat(pd)
    register_current_anndata()
    model = scvi_module.model.SCVI.load(
        str(model_dir),
        adata=adata,
        accelerator="cpu",
        device=1,
    )
    align_batch_categories_to_scvi_registry(adata, model, pd)
    # Legacy model loading may re-register the AnnData with an older-style registry.
    # Re-run the current setup so the rehydrated model inherits a modern field registry.
    register_current_anndata()
    init_kwargs = extract_scvi_init_kwargs(
        getattr(model, "init_params_", {}),
        getattr(model, "_module_kwargs", None),
    )
    rehydrated_model = scvi_module.model.SCVI(adata, **init_kwargs)
    state_load = rehydrated_model.module.load_state_dict(model.module.state_dict(), strict=False)
    missing_keys = list(getattr(state_load, "missing_keys", []))
    unexpected_keys = list(getattr(state_load, "unexpected_keys", []))
    if missing_keys or unexpected_keys:
        raise RuntimeError(
            "Failed to rehydrate legacy scVI weights into current registry format. "
            f"missing_keys={missing_keys}, unexpected_keys={unexpected_keys}"
        )

    rehydrated_model.is_trained_ = True
    rehydrated_model.train_indices_ = getattr(model, "train_indices_", None)
    rehydrated_model.test_indices_ = getattr(model, "test_indices_", None)
    rehydrated_model.validation_indices_ = getattr(model, "validation_indices_", None)
    rehydrated_model.history_ = getattr(model, "history_", None)

    del model
    gc.collect()
    return rehydrated_model, "cpu"


def align_batch_categories_to_scvi_registry(adata, scvi_model, pd_module) -> None:
    try:
        # scvi-tools expects the registry key here, not the original obs column name.
        state_registry = scvi_model.adata_manager.get_state_registry("batch")
        categorical_mapping = list(state_registry.categorical_mapping)
    except Exception:
        return

    current = clean_string_series(adata.obs[BATCH_KEY], pd_module, fill_value="Unknown")
    adata.obs[BATCH_KEY] = pd_module.Categorical(current, categories=categorical_mapping)


def train_scanvi(adata, scvi_model, scvi_module, torch_module, n_samples_per_label: int, model_dir: Path, scvi_model_dir: str | Path):
    import pandas as pd  # local import to keep module import-safe for tests

    align_batch_categories_to_scvi_registry(adata, scvi_model, pd)

    def build_from_scvi(base_model):
        return scvi_module.model.SCANVI.from_scvi_model(
            base_model,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key=LABELS_KEY,
        )

    def is_cuda_training_error(exc: Exception) -> bool:
        msg = str(exc).lower()
        return any(token in msg for token in ["cuda", "cudnn", "nccl", "device-side", "acceleratorerror"])

    accelerator, _device = resolve_accelerator(torch_module)
    train_kwargs = dict(
        max_epochs=SCANVI_EPOCHS,
        batch_size=BATCH_SIZE,
        early_stopping=True,
        early_stopping_patience=SCANVI_PATIENCE,
        train_size=0.9,
        devices=1,
        plan_kwargs={"lr": SCANVI_LR, "weight_decay": 0.0},
        n_samples_per_label=n_samples_per_label,
    )

    model = build_from_scvi(scvi_model)
    try:
        model.train(accelerator=accelerator, **train_kwargs)
    except Exception as exc:
        if accelerator != "gpu" or not is_cuda_training_error(exc):
            raise
        print(f"[WARN] SCANVI GPU training failed; retrying on CPU. Error: {exc}")
        del model
        gc.collect()
        if torch_module.cuda.is_available():
            torch_module.cuda.empty_cache()
        if hasattr(scvi_model, "to_device"):
            scvi_model.to_device("cpu")
        elif hasattr(scvi_model, "module"):
            scvi_model.module.to("cpu")
        model = build_from_scvi(scvi_model)
        model.train(accelerator="cpu", **train_kwargs)

    model.save(str(model_dir), overwrite=True)
    return model


def compute_umap(adata, rep_key: str, output_dir: Path):
    bundle_info = fit_bundle(
        adata,
        rep_key,
        output_dir,
        primary_umap_key=UMAP_SCANVI_KEY,
        alias_keys=(UMAP_SCANVI_CORRECTED_KEY,),
        operator_filename="umap_scanvi_operator.joblib",
        manifest_filename="scanvi_umap_bundle.json",
        umap_params={
            "n_neighbors": N_NEIGHBORS,
            "n_components": 2,
            "min_dist": 0.3,
            "spread": 1.0,
            "metric": "euclidean",
            "random_state": SEED,
        },
        set_default_x_umap=True,
        extra_manifest={"branch_pipeline": Path(__file__).name},
    )
    return bundle_info


def plot_overview(adata, spec: BranchSpec, sc_module, plt_module, output_path: Path) -> None:
    sc_module.settings.vector_friendly = True
    fig, axes = plt_module.subplots(2, 2, figsize=(18, 14))
    sc_module.pl.embedding(
        adata, basis="umap", color=LABELS_KEY, ax=axes[0, 0], show=False,
        legend_loc="right margin", legend_fontsize=6, title=f"{spec.name}: supervision"
    )
    sc_module.pl.embedding(
        adata, basis="umap", color=PRED_KEY, ax=axes[0, 1], show=False,
        legend_loc="right margin", legend_fontsize=6, title=f"{spec.name}: scanvi pred"
    )
    sc_module.pl.embedding(
        adata, basis="umap", color=L2_KEY, ax=axes[1, 0], show=False,
        legend_loc="right margin", legend_fontsize=6, title=f"{spec.name}: L2"
    )
    sc_module.pl.embedding(
        adata, basis="umap", color=PRED_PROB_KEY, ax=axes[1, 1], show=False,
        color_map="RdYlGn", vmin=0, vmax=1, title=f"{spec.name}: prediction confidence"
    )
    plt_module.tight_layout()
    fig.savefig(output_path, bbox_inches="tight", dpi=DPI)
    plt_module.close(fig)


def sanitize_for_write(adata, np_module, pd_module) -> None:
    for df in [adata.obs, adata.var] + ([adata.raw.var] if adata.raw is not None else []):
        for col in df.columns:
            s = df[col]
            if s.dtype != object:
                continue
            non_na = s.dropna()
            if len(non_na) == 0:
                df[col] = s.fillna("").astype(str)
                continue
            if non_na.map(lambda x: isinstance(x, (bool, np_module.bool_))).all():
                df[col] = s.fillna(False).astype(np_module.int8)
                continue
            if non_na.map(lambda x: not isinstance(x, str)).any():
                df[col] = s.map(lambda x: "" if pd_module.isna(x) else str(x))


def apply_branch_operation(adata, spec: BranchSpec, target_cells, target_obs, pd_module):
    if spec.mode == "remove":
        missing_cells = pd_module.Index(target_cells).difference(adata.obs_names)
        if len(missing_cells) > 0:
            raise RuntimeError(
                f"{len(missing_cells)} target cells are missing from {spec.input_h5ad}; first few: {missing_cells[:10].tolist()}"
            )
        keep_mask = ~adata.obs_names.isin(target_cells)
        return adata[keep_mask].copy(), target_obs.assign(action="removed")

    if spec.mode == "relabel":
        for col in [L2_KEY, L3_KEY, LABELS_KEY, PRED_KEY]:
            if col in adata.obs.columns:
                adata.obs[f"{col}_{spec.label_backup_suffix}"] = clean_string_series(adata.obs[col], pd_module)

        target_index = pd_module.Index(target_cells).intersection(adata.obs_names)
        if len(target_index) == 0:
            raise RuntimeError(f"No relabel target cells overlapped with {spec.input_h5ad}")

        choir_lookup = target_obs[CHOIR_CLUSTER_COL].to_dict()
        relabel_series = pd_module.Series(index=target_index, dtype=object)
        for cell in target_index:
            cluster_id = normalize_cluster_id(choir_lookup[cell])
            relabel_series.loc[cell] = spec.relabel_map[cluster_id]

        for col in [L2_KEY, L3_KEY, LABELS_KEY, PRED_KEY]:
            if col in adata.obs.columns:
                current = clean_string_series(adata.obs[col], pd_module)
                current.loc[target_index] = relabel_series.astype(str)
                adata.obs[col] = to_clean_category(current, pd_module)
        action_table = target_obs.loc[target_index].copy()
        action_table["new_label"] = relabel_series.astype(str)
        action_table["action"] = "relabeled"
        return adata, action_table

    raise ValueError(f"Unsupported mode: {spec.mode}")


def run_branch(spec: BranchSpec) -> None:
    ad, plt, np, pd, sc, scvi, torch, sparse = import_runtime_modules()
    output_paths = get_output_paths(spec)
    for key in ["root", "reports", "figures"]:
        output_paths[key].mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(output_paths["figures"])

    print("=" * 80)
    print(f"Running branch: {spec.name} ({spec.mode})")
    print("=" * 80)
    print(f"cluster source : {spec.cluster_source_h5ad}")
    print(f"input h5ad     : {spec.input_h5ad}")
    print(f"output dir     : {output_paths['root']}")
    print(f"target clusters: {spec.target_clusters}")
    if spec.mode == "relabel":
        print(f"relabel map    : {spec.relabel_map}")

    target_cells, target_obs = summarize_target_clusters(spec, pd)
    save_table(
        target_obs.reset_index(names="cell_barcode"),
        output_paths["reports"] / "target_cluster_cells.tsv",
    )
    cluster_summary = (
        target_obs[CHOIR_CLUSTER_COL]
        .value_counts()
        .rename_axis(CHOIR_CLUSTER_COL)
        .reset_index(name="n_cells")
        .sort_values(CHOIR_CLUSTER_COL, key=lambda s: s.astype(int))
    )
    save_table(cluster_summary, output_paths["reports"] / "target_cluster_summary.tsv")

    adata = sc.read_h5ad(spec.input_h5ad)
    adata, action_table = apply_branch_operation(adata, spec, target_cells, target_obs, pd)
    save_table(action_table.reset_index(names="cell_barcode"), output_paths["reports"] / "cluster_actions.tsv")

    label_values, min_class_size, n_samples_per_label = ensure_training_ready(adata, pd, sparse)
    save_table(
        pd.DataFrame({LABELS_KEY: label_values}),
        output_paths["reports"] / "label_levels.tsv",
    )
    continuous_covariates = get_continuous_covariates(adata)

    scvi_reused = False
    if spec.existing_scvi_model_dir is not None:
        try:
            scvi_model, accelerator = load_existing_scvi(
                adata,
                scvi,
                torch,
                continuous_covariates,
                spec.existing_scvi_model_dir,
            )
            output_paths["scvi_model"].mkdir(parents=True, exist_ok=True)
            scvi_model.save(str(output_paths["scvi_model"]), overwrite=True)
            write_json(
                {
                    "reused_existing_scvi_model_dir": spec.existing_scvi_model_dir,
                    "reuse_mode": "legacy_weight_rehydration",
                    "saved_rehydrated_scvi_model_dir": str(output_paths["scvi_model"]),
                },
                output_paths["scvi_model"] / "reuse_manifest.json",
            )
            scvi_reused = True
        except Exception as exc:
            print(
                f"[WARN] Failed to reuse existing scVI model for {spec.name}; "
                f"falling back to full scVI retraining. Error: {exc}"
            )
            scvi_model, accelerator = train_scvi_from_scratch(
                adata,
                scvi,
                torch,
                continuous_covariates,
                output_paths["scvi_model"],
            )
    else:
        scvi_model, accelerator = train_scvi_from_scratch(
            adata,
            scvi,
            torch,
            continuous_covariates,
            output_paths["scvi_model"],
        )

    adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
    try:
        scanvi_model = train_scanvi(
            adata,
            scvi_model,
            scvi,
            torch,
            n_samples_per_label,
            output_paths["scanvi_model"],
            output_paths["scvi_model"],
        )
    except Exception as exc:
        if not scvi_reused:
            raise
        print(
            f"[WARN] Reused scVI -> scANVI handoff failed for {spec.name}; "
            f"falling back to full scVI retraining before retrying scANVI. Error: {exc}"
        )
        del scvi_model
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        scvi_model, accelerator = train_scvi_from_scratch(
            adata,
            scvi,
            torch,
            continuous_covariates,
            output_paths["scvi_model"],
        )
        adata.obsm[SCVI_LATENT_KEY] = scvi_model.get_latent_representation()
        scanvi_model = train_scanvi(
            adata,
            scvi_model,
            scvi,
            torch,
            n_samples_per_label,
            output_paths["scanvi_model"],
            output_paths["scvi_model"],
        )
        scvi_reused = False
    adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()
    pred = scanvi_model.predict()
    pred_soft = scanvi_model.predict(soft=True)
    pred_soft_arr = pred_soft.to_numpy() if hasattr(pred_soft, "to_numpy") else np.asarray(pred_soft)
    adata.obs[PRED_KEY] = pd.Categorical(pred.astype(str))
    adata.obs[PRED_PROB_KEY] = pred_soft_arr.max(axis=1).astype(np.float32)

    umap_info = compute_umap(adata, SCANVI_LATENT_KEY, output_paths["root"])
    plot_overview(adata, spec, sc, plt, output_paths["figures"] / f"{spec.name}_rerun_overview.pdf")

    adata.uns["stromal_branch_rerun"] = {
        "script_name": Path(__file__).name,
        "pipeline_date": PIPELINE_DATE,
        "branch": spec.name,
        "mode": spec.mode,
        "cluster_source_h5ad": spec.cluster_source_h5ad,
        "input_h5ad": spec.input_h5ad,
        "output_dir": str(output_paths["root"]),
        "final_h5ad": str(output_paths["final_h5ad"]),
        "scvi_model_dir": str(output_paths["scvi_model"]),
        "scanvi_model_dir": str(output_paths["scanvi_model"]),
        "target_clusters": list(spec.target_clusters),
        "relabel_map": spec.relabel_map,
        "existing_scvi_model_dir": spec.existing_scvi_model_dir,
        "scvi_reused": bool(scvi_reused),
        "continuous_covariates": continuous_covariates,
        "min_class_size": int(min_class_size),
        "n_samples_per_label": int(n_samples_per_label),
        "accelerator": accelerator,
        "umap_operator_path": str(umap_info["operator_path"]),
        "umap_manifest_path": str(umap_info["manifest_path"]),
        "target_cell_count": int(len(target_cells)),
        "n_obs_final": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
    }

    save_table(
        pd.DataFrame({
            "metric": [
                "target_cell_count",
                "n_obs_final",
                "n_vars",
                "min_class_size",
                "n_samples_per_label",
            ],
            "value": [
                int(len(target_cells)),
                int(adata.n_obs),
                int(adata.n_vars),
                int(min_class_size),
                int(n_samples_per_label),
            ],
        }),
        output_paths["reports"] / "rerun_summary.tsv",
    )

    import anndata as _anndata

    _anndata.settings.allow_write_nullable_strings = True
    sanitize_for_write(adata, np, pd)
    adata.write_h5ad(output_paths["final_h5ad"], compression="gzip", compression_opts=9)

    del scanvi_model, scvi_model
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    print(f"[OK] Saved final h5ad: {output_paths['final_h5ad']}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stromal branch rerun after CHOIR cleanup / relabel")
    parser.add_argument(
        "--branch",
        choices=["all", *BRANCH_SPECS.keys()],
        default="all",
        help="Which branch to run",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.branch == "all":
        branches = [BRANCH_SPECS[name] for name in ("endothelial", "fibroblast", "smc")]
    else:
        branches = [BRANCH_SPECS[args.branch]]
    for spec in branches:
        run_branch(spec)


if __name__ == "__main__":
    main()
