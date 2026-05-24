#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Joint reference-query visualization for allcell and B-cell modes.

Outputs per mode:
    - joint_umap_overview.{pdf,png}
    - joint_dotplot.{pdf,png}
    - joint_featureplot.{pdf,png}
    - joint_proportion_violin.{pdf,png}
    - marker_availability.tsv
    - source_label_counts.tsv
    - summary.json
    - lightweight marker-only AnnData for plotting reuse

The allcell path prioritizes the exact 0115 reference-model lineage:
    1128 allcell training input -> 0115 scanvi_existing_model -> 0127 mapped query.

For practicality, continuous covariates are regenerated from counts using the
original training gene sets; cell-cycle scores are approximated from normalized
signature means so we can recover the exact model-matched reference latent space
without rerunning the historical training pipeline.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import warnings
from collections import OrderedDict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex
from matplotlib.lines import Line2D


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "data" / "py" / "20260419" / "joint_ref_query_visualization"

REFERENCE_COLOR = "#4C78A8"
QUERY_COLOR = "#E45756"
UNKNOWN_COLOR = "#B0B0B0"
SOURCE_PALETTE = {"reference": REFERENCE_COLOR, "query": QUERY_COLOR}

EXPR_CMAP = LinearSegmentedColormap.from_list(
    "joint_expr", ["#F2F2F2", "#FDBB84", "#E34A33", "#7F0000"]
)

ALLCELL_LABEL_PALETTE = {
    "Epithelial": "#4E79A7",
    "T": "#F28E2B",
    "B": "#E15759",
    "Myeloid": "#76B7B2",
    "Fibroblast": "#59A14F",
    "Endothelial": "#EDC948",
    "SMC": "#B07AA1",
    "Unknown": UNKNOWN_COLOR,
}

B_CELL_LABEL_PALETTE = {
    "Naive_B": "#4E79A7",
    "Memory_B": "#F28E2B",
    "GC_B": "#E15759",
    "Atypical_Memory_B": "#76B7B2",
    "Plasma": "#59A14F",
    "Unknown": UNKNOWN_COLOR,
    "nan": UNKNOWN_COLOR,
}

ALLCELL_MARKERS = OrderedDict(
    {
        "Epithelial": ["EPCAM", "KRT8", "KRT19", "SCGB1A1", "AGER"],
        "T": ["CD3D", "CD3E", "TRBC1", "IL7R", "LTB"],
        "B": ["CD79A", "MS4A1", "CD74", "HLA-DRA", "BANK1"],
        "Myeloid": ["LYZ", "FCER1G", "LST1", "CTSS", "C1QC"],
        "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5", "EMCN", "KDR"],
        "SMC": ["ACTA2", "TAGLN", "MYH11", "CNN1", "RGS5"],
    }
)

B_CELL_MARKERS = OrderedDict(
    {
        "Naive_B": ["TCL1A", "IGHD", "IGHM", "FCER2", "IL4R"],
        "Memory_B": ["CD27", "BANK1", "TNFRSF13B", "CD44", "AIM2"],
        "GC_B": ["AICDA", "RGS13", "BCL6", "CD83", "MKI67"],
        "Atypical_Memory_B": ["FCRL5", "ITGAX", "TBX21", "CXCR3", "ZEB2"],
        "Plasma": ["MZB1", "JCHAIN", "XBP1", "SDC1", "TNFRSF17"],
    }
)

ALLCELL_FEATURE_GENES = [
    "EPCAM",
    "AGER",
    "CD3D",
    "IL7R",
    "CD79A",
    "MS4A1",
    "LYZ",
    "C1QC",
    "COL1A1",
    "DCN",
    "PECAM1",
    "CLDN5",
    "ACTA2",
    "MYH11",
]

B_CELL_FEATURE_GENES = [
    "TCL1A",
    "IGHD",
    "CD27",
    "BANK1",
    "AICDA",
    "RGS13",
    "FCRL5",
    "ITGAX",
    "MZB1",
    "JCHAIN",
    "XBP1",
    "SDC1",
]

ALLCELL_STRESS_SIGNATURE_GENES = [
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

ALLCELL_S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

ALLCELL_G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


def configure_runtime() -> None:
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

    user_site_fragments = ("/.local/lib/python", "site-packages")
    sys.path[:] = [
        p
        for p in sys.path
        if not (all(fragment in p for fragment in user_site_fragments))
    ]

    try:
        import anndata

        anndata.settings.allow_write_nullable_strings = True
    except Exception:
        pass


def get_mode_config(mode: str) -> dict:
    mode = str(mode).strip().lower()
    configs = {
        "allcell": {
            "mode": "allcell",
            "query_h5ad": PROJECT_ROOT / "data" / "py" / "0127" / "scarches_mapping_FIXED_v1_2" / "query_mapped_to_reference.h5ad",
            "reference_input_h5ad": PROJECT_ROOT / "data" / "py" / "1128" / "bbknn_annotation_analysis" / "adata_bbknn_annotated_corrected_filtered.h5ad",
            "whitelist_dir": PROJECT_ROOT / "data" / "core20260115",
            "scanvi_model_dir": PROJECT_ROOT / "data" / "py" / "0115" / "allcells_scvi_analysis_v2_3_1_HOTFIX" / "models" / "scanvi_existing_model",
            "hvg_file": PROJECT_ROOT / "data" / "py" / "0115" / "allcells_scvi_analysis_v2_3_1_HOTFIX" / "models" / "scvi_model" / "hvg_genes.txt",
            "reference_label_key": "cell_type",
            "query_final_key": "cell_type_final",
            "query_pred_key": "cell_type_mapped",
            "query_conf_key": "mapping_confidence",
            "sample_key": "sample",
            "tissue_key": "tissue",
            "latent_key_query": "X_scANVI_mapped",
            "marker_dict": ALLCELL_MARKERS,
            "feature_genes": ALLCELL_FEATURE_GENES,
            "label_palette": ALLCELL_LABEL_PALETTE,
            "output_dir": OUTPUT_ROOT / "allcell",
        },
        "bcell": {
            "mode": "bcell",
            "candidate_merged_h5ads": [
                PROJECT_ROOT / "data" / "py" / "0419" / "bcell_merge_schpl_v1_3_scanvi_rerun_20260419" / "bcell_reference_plus_query_schpl_v1_3_scanvi_rerun_20260419.h5ad",
                PROJECT_ROOT / "data" / "py" / "0413" / "bcell_merge_schpl_v1_2_scanvi_umap_20260413" / "bcell_reference_plus_query_schpl_v1_2_scanvi_umap_20260413.h5ad",
                PROJECT_ROOT / "data" / "py" / "0412" / "bcell_merge_schpl_unified_v1_0" / "bcell_reference_plus_query_schpl_v1_0_unified_20260412.h5ad",
            ],
            "reference_label_key": "Cell_Type_L2",
            "query_final_key": "Cell_Type_L2_final",
            "query_pred_key": "Cell_Type_L2_pred",
            "query_conf_key": "mapping_confidence",
            "source_key": "data_source",
            "sample_key": "sample",
            "tissue_key": "tissue",
            "umap_key": "X_umap",
            "marker_dict": B_CELL_MARKERS,
            "feature_genes": B_CELL_FEATURE_GENES,
            "label_palette": B_CELL_LABEL_PALETTE,
            "output_dir": OUTPUT_ROOT / "bcell",
        },
    }
    if mode not in configs:
        raise KeyError(f"Unsupported mode: {mode}")
    return dict(configs[mode])


def resolve_umap_key(keys_or_adata, preferred: str | None = None) -> str:
    if hasattr(keys_or_adata, "obsm"):
        available = list(keys_or_adata.obsm.keys())
    else:
        available = list(keys_or_adata)
    candidates = [preferred, "X_umap", "X_umap_mapped", "X_umap_scanvi", "X_umap_scvi", "umap"]
    for candidate in candidates:
        if candidate and candidate in available:
            return candidate
    raise KeyError(f"No UMAP-like key found. Available keys: {available}")


def filter_available_markers(marker_dict: Mapping[str, Sequence[str]], var_names: Iterable[str]):
    var_set = {str(g) for g in var_names}
    filtered: OrderedDict[str, list[str]] = OrderedDict()
    missing: OrderedDict[str, list[str]] = OrderedDict()
    for group, genes in marker_dict.items():
        keep = [str(g) for g in genes if str(g) in var_set]
        drop = [str(g) for g in genes if str(g) not in var_set]
        if keep:
            filtered[group] = keep
        missing[group] = drop
    return filtered, missing


def compute_group_proportions(
    obs: pd.DataFrame,
    sample_key: str,
    label_key: str,
    source_key: str,
) -> pd.DataFrame:
    frame = obs[[sample_key, label_key, source_key]].copy()
    for col in [sample_key, label_key, source_key]:
        frame[col] = _to_clean_string(frame[col])
    counts = (
        frame.groupby([sample_key, source_key, label_key], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = counts.groupby([sample_key, source_key], observed=False)["n_cells"].transform("sum")
    counts["total_cells"] = totals.astype(int)
    counts["proportion"] = counts["n_cells"] / counts["total_cells"].replace(0, np.nan)
    counts["proportion"] = counts["proportion"].fillna(0.0)
    return counts


def build_plot_label(
    obs: pd.DataFrame,
    source_key: str,
    ref_label_key: str,
    query_label_key: str,
    query_fallback_key: str | None = None,
    unknown_label: str = "Unknown",
) -> pd.Series:
    source = _to_clean_string(obs[source_key], unknown=unknown_label).str.lower()
    ref = _to_clean_string(obs.get(ref_label_key), unknown=unknown_label)
    qry = _to_clean_string(obs.get(query_label_key), unknown=unknown_label)
    if query_fallback_key is not None:
        qry = qry.mask(qry.eq(unknown_label), _to_clean_string(obs.get(query_fallback_key), unknown=unknown_label))
    labels = pd.Series(np.where(source.eq("reference"), ref, qry), index=obs.index, dtype=object)
    labels = labels.fillna(unknown_label).replace({"": unknown_label, "nan": unknown_label, "None": unknown_label})
    return labels.astype(str)


def select_best_bcell_candidate(candidates: Sequence[Mapping[str, object]]) -> Mapping[str, object]:
    if not candidates:
        raise ValueError("No B-cell candidates provided")

    def score(candidate: Mapping[str, object]) -> tuple:
        return (
            int(bool(candidate.get("has_data_source"))),
            int(bool(candidate.get("has_umap"))),
            int(bool(candidate.get("has_query_final"))),
            int(candidate.get("marker_hits", 0)),
        )

    return max(candidates, key=score)


def flatten_marker_dict(marker_dict: Mapping[str, Sequence[str]]) -> tuple[list[str], list[tuple[str, int, int]]]:
    genes: list[str] = []
    bounds: list[tuple[str, int, int]] = []
    start = 0
    for group, marker_genes in marker_dict.items():
        unique_genes = [g for g in marker_genes if g not in genes]
        genes.extend(unique_genes)
        end = start + len(unique_genes)
        bounds.append((group, start, end))
        start = end
    return genes, bounds


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


def inject_pandas_pickle_compat(pd_module) -> None:
    import types

    if "pandas.core.indexes.numeric" in sys.modules:
        return
    compat_mod = types.ModuleType("pandas.core.indexes.numeric")
    compat_mod.Int64Index = pd_module.Index
    compat_mod.UInt64Index = pd_module.Index
    compat_mod.Float64Index = pd_module.Index
    sys.modules["pandas.core.indexes.numeric"] = compat_mod


def align_batch_categories_to_scvi_registry(adata: ad.AnnData, scvi_model, pd_module) -> None:
    try:
        state_registry = scvi_model.adata_manager.get_state_registry("batch")
        categorical_mapping = list(state_registry.categorical_mapping)
    except Exception:
        return

    current = _to_clean_string(adata.obs["sample"], unknown="Unknown")
    adata.obs["sample"] = pd_module.Categorical(current, categories=categorical_mapping)


def load_existing_scvi_model(adata_scvi: ad.AnnData, scvi_module, model_dir: Path):
    import gc

    pd_module = pd
    inject_pandas_pickle_compat(pd_module)
    scvi_module.model.SCVI.setup_anndata(
        adata_scvi,
        layer=None,
        batch_key="sample",
        categorical_covariate_keys=["tissue"],
        continuous_covariate_keys=["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
    )
    loaded_model = scvi_module.model.SCVI.load(
        str(model_dir),
        adata=adata_scvi,
        accelerator="cpu",
        device=1,
    )
    align_batch_categories_to_scvi_registry(adata_scvi, loaded_model, pd_module)
    scvi_module.model.SCVI.setup_anndata(
        adata_scvi,
        layer=None,
        batch_key="sample",
        categorical_covariate_keys=["tissue"],
        continuous_covariate_keys=["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
    )
    init_kwargs = extract_scvi_init_kwargs(
        getattr(loaded_model, "init_params_", {}),
        getattr(loaded_model, "_module_kwargs", None),
    )
    rehydrated_model = scvi_module.model.SCVI(adata_scvi, **init_kwargs)
    state_load = rehydrated_model.module.load_state_dict(loaded_model.module.state_dict(), strict=False)
    missing_keys = list(getattr(state_load, "missing_keys", []))
    unexpected_keys = list(getattr(state_load, "unexpected_keys", []))
    if missing_keys or unexpected_keys:
        raise RuntimeError(
            "Failed to rehydrate legacy scVI weights into the current registry format. "
            f"missing_keys={missing_keys}, unexpected_keys={unexpected_keys}"
        )

    rehydrated_model.is_trained_ = True
    rehydrated_model.train_indices_ = getattr(loaded_model, "train_indices_", None)
    rehydrated_model.test_indices_ = getattr(loaded_model, "test_indices_", None)
    rehydrated_model.validation_indices_ = getattr(loaded_model, "validation_indices_", None)
    rehydrated_model.history_ = getattr(loaded_model, "history_", None)

    del loaded_model
    gc.collect()
    return rehydrated_model


def load_scanvi_from_checkpoint(adata_scanvi: ad.AnnData, scvi_model, scvi_module, model_dir: Path):
    import torch

    checkpoint = torch.load(model_dir / "model.pt", map_location="cpu")
    scanvi_model = scvi_module.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata_scanvi,
        labels_key="scanvi_label_existing_cleaned",
        unlabeled_category="Unknown",
    )
    state_load = scanvi_model.module.load_state_dict(checkpoint["model_state_dict"], strict=False)
    missing_keys = list(getattr(state_load, "missing_keys", []))
    unexpected_keys = list(getattr(state_load, "unexpected_keys", []))
    if missing_keys or unexpected_keys:
        warnings.warn(
            "Legacy SCANVI state dict did not match perfectly during rehydration. "
            f"missing_keys={missing_keys}, unexpected_keys={unexpected_keys}"
        )

    attr_dict = checkpoint.get("attr_dict", {})
    scanvi_model.is_trained_ = bool(attr_dict.get("is_trained_", True))
    for key in [
        "train_indices_",
        "test_indices_",
        "validation_indices_",
        "history_",
        "semisupervised_history_",
        "unsupervised_history_",
        "labels_",
        "unlabeled_category_",
    ]:
        if key in attr_dict:
            setattr(scanvi_model, key, attr_dict[key])
    return scanvi_model


def _to_clean_string(values, unknown: str = "Unknown") -> pd.Series:
    if values is None:
        return pd.Series(dtype=object)
    series = pd.Series(values).copy()
    series = series.astype("string")
    series = series.fillna(unknown)
    series = series.replace({"": unknown, "nan": unknown, "None": unknown, "<NA>": unknown})
    return series.astype(str)


def _make_json_safe(value):
    if isinstance(value, pd.DataFrame):
        return {
            "type": "DataFrame",
            "n_rows": int(value.shape[0]),
            "n_cols": int(value.shape[1]),
            "columns": [str(col) for col in value.columns.tolist()],
        }
    if isinstance(value, pd.Series):
        return {
            "type": "Series",
            "n_rows": int(value.shape[0]),
            "name": None if value.name is None else str(value.name),
        }
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {str(k): _make_json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_make_json_safe(v) for v in value]
    return value


def _lazy_import_scanpy():
    import scanpy as sc

    sc.settings.verbosity = 0
    sc.settings.set_figure_params(dpi=120, dpi_save=300, facecolor="white", frameon=False)
    return sc


def _lazy_import_seaborn():
    import seaborn as sns

    return sns


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_counts_layer(adata: ad.AnnData) -> ad.AnnData:
    if "counts" in adata.layers:
        return adata
    if adata.raw is not None:
        adata.layers["counts"] = adata.raw.X.copy()
        return adata
    adata.layers["counts"] = adata.X.copy()
    return adata


def _read_obs_index_from_h5ad(h5ad_path: str) -> list[str]:
    import h5py

    try:
        with h5py.File(h5ad_path, "r") as handle:
            if "obs" in handle and "_index" in handle["obs"]:
                raw = handle["obs"]["_index"][()]
                if getattr(raw, "dtype", None) is not None and raw.dtype.kind in ("S", "O"):
                    return [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in raw]
                return raw.astype(str).tolist()
    except Exception:
        pass

    temp = ad.read_h5ad(h5ad_path, backed="r")
    try:
        return temp.obs_names.astype(str).tolist()
    finally:
        if getattr(temp, "file", None) is not None:
            temp.file.close()


def collect_cell_ids_from_folder(folder: str | Path, pattern: str = "**/*.h5ad", mode: str = "union") -> set[str]:
    folder = Path(folder)
    files = sorted(folder.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No h5ad files found under {folder} with pattern {pattern}")

    cell_set: set[str] | None = None
    for file_path in files:
        ids = set(_read_obs_index_from_h5ad(str(file_path)))
        if cell_set is None:
            cell_set = ids
        elif mode == "intersection":
            cell_set &= ids
        else:
            cell_set |= ids
    return cell_set or set()


def _sum_axis1(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)


def _select_matrix_source(adata: ad.AnnData):
    if "counts" in adata.layers:
        return adata.layers["counts"], adata.var_names.astype(str), "counts"
    if adata.raw is not None:
        return adata.raw.X, adata.raw.var_names.astype(str), "raw"
    return adata.X, adata.var_names.astype(str), "X"


def extract_expression_matrix(
    adata: ad.AnnData,
    genes: Sequence[str],
    prefer_counts: bool = True,
) -> tuple[np.ndarray, dict]:
    genes = [str(g) for g in genes]
    if prefer_counts and "counts" in adata.layers:
        matrix = adata.layers["counts"]
        var_names = adata.var_names.astype(str)
        source = "counts_log1p_cpm"
        normalize = True
    else:
        matrix, var_names, base_source = _select_matrix_source(adata)
        source = base_source
        normalize = False

    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    out = np.zeros((adata.n_obs, len(genes)), dtype=np.float32)
    present: list[str] = []
    missing: list[str] = []
    idxs: list[int] = []
    out_cols: list[int] = []
    for out_idx, gene in enumerate(genes):
        match = gene_to_idx.get(gene)
        if match is None:
            missing.append(gene)
            continue
        present.append(gene)
        idxs.append(match)
        out_cols.append(out_idx)

    if idxs:
        sub = matrix[:, idxs]
        if normalize:
            totals = _sum_axis1(matrix)
            if sparse.issparse(sub):
                sub = sub.tocsr(copy=True).astype(np.float32)
                scale = np.zeros_like(totals, dtype=np.float32)
                nz = totals > 0
                scale[nz] = 1e4 / totals[nz]
                sub = sparse.diags(scale) @ sub
                sub.data = np.log1p(sub.data)
                dense = sub.toarray().astype(np.float32)
            else:
                dense = np.asarray(sub, dtype=np.float32)
                nz = totals > 0
                dense[nz] = np.log1p((dense[nz] / totals[nz, None]) * 1e4)
                dense[~nz] = 0.0
        else:
            dense = sub.toarray().astype(np.float32) if sparse.issparse(sub) else np.asarray(sub, dtype=np.float32)
        out[:, out_cols] = dense

    return out, {"present": present, "missing": missing, "source": source}


def add_allcell_covariates_fast(adata: ad.AnnData) -> ad.AnnData:
    ensure_counts_layer(adata)
    counts = adata.layers["counts"]
    var_names = adata.var_names.astype(str)

    mt_mask = np.array([g.startswith("MT-") or g.startswith("Mt-") or g.startswith("mt-") for g in var_names])
    total = _sum_axis1(counts)
    if mt_mask.any():
        mt_counts = _sum_axis1(counts[:, mt_mask])
        pct_counts_mt = np.zeros_like(total, dtype=np.float32)
        nz = total > 0
        pct_counts_mt[nz] = (mt_counts[nz] / total[nz]) * 100.0
        adata.obs["pct_counts_mt"] = pct_counts_mt
    else:
        adata.obs["pct_counts_mt"] = 0.0

    stress_matrix, _ = extract_expression_matrix(adata, ALLCELL_STRESS_SIGNATURE_GENES, prefer_counts=True)
    adata.obs["stress_score"] = stress_matrix.mean(axis=1).astype(np.float32) if stress_matrix.shape[1] else 0.0

    s_matrix, _ = extract_expression_matrix(adata, ALLCELL_S_GENES, prefer_counts=True)
    g2m_matrix, _ = extract_expression_matrix(adata, ALLCELL_G2M_GENES, prefer_counts=True)
    s_score = s_matrix.mean(axis=1).astype(np.float32) if s_matrix.shape[1] else np.zeros(adata.n_obs, dtype=np.float32)
    g2m_score = g2m_matrix.mean(axis=1).astype(np.float32) if g2m_matrix.shape[1] else np.zeros(adata.n_obs, dtype=np.float32)
    adata.obs["S_score"] = s_score
    adata.obs["G2M_score"] = g2m_score

    phase = np.full(adata.n_obs, "G1", dtype=object)
    phase[s_score > np.maximum(g2m_score, 0)] = "S"
    phase[g2m_score > np.maximum(s_score, 0)] = "G2M"
    adata.obs["phase"] = phase
    return adata


def _estimate_point_size(n_obs: int) -> float:
    return float(np.clip(120000 / max(n_obs, 1), 0.35, 3.5))


def _pick_palette(labels: Sequence[str], preferred: Mapping[str, str] | None = None) -> dict[str, str]:
    unique = list(dict.fromkeys([str(x) for x in labels]))
    palette = dict(preferred or {})
    missing = [label for label in unique if label not in palette]
    if missing:
        cmap = plt.get_cmap("tab20", max(len(missing), 1))
        for idx, label in enumerate(missing):
            palette[label] = to_hex(cmap(idx))
    return {label: palette[label] for label in unique}


def _style_axis(ax, title: str) -> None:
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)


def scatter_categorical(
    ax,
    coords: np.ndarray,
    labels: Sequence[str],
    title: str,
    palette: Mapping[str, str] | None = None,
    point_size: float | None = None,
    order: Sequence[str] | None = None,
):
    labels_series = _to_clean_string(labels)
    categories = list(order) if order is not None else list(dict.fromkeys(labels_series.tolist()))
    categories = [cat for cat in categories if cat in set(labels_series.tolist())]
    palette = _pick_palette(categories, preferred=palette)
    if point_size is None:
        point_size = _estimate_point_size(len(labels_series))

    for category in categories:
        mask = labels_series.eq(category).to_numpy()
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            c=palette[category],
            linewidths=0,
            alpha=0.9,
            rasterized=True,
            label=category,
        )

    _style_axis(ax, title)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", markerfacecolor=palette[cat], markeredgecolor="none", markersize=6, label=cat)
        for cat in categories
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9)


def scatter_query_confidence(
    ax,
    coords: np.ndarray,
    source: Sequence[str],
    confidence: Sequence[float],
    title: str,
    point_size: float | None = None,
):
    source_series = _to_clean_string(source).str.lower()
    confidence_series = pd.to_numeric(pd.Series(confidence), errors="coerce")
    if point_size is None:
        point_size = _estimate_point_size(len(source_series))

    ref_mask = source_series.eq("reference").to_numpy()
    qry_mask = source_series.eq("query").to_numpy()
    ax.scatter(
        coords[ref_mask, 0],
        coords[ref_mask, 1],
        s=point_size,
        c="#DDDDDD",
        linewidths=0,
        alpha=0.7,
        rasterized=True,
    )
    sc = ax.scatter(
        coords[qry_mask, 0],
        coords[qry_mask, 1],
        s=point_size,
        c=confidence_series[qry_mask],
        linewidths=0,
        alpha=0.95,
        rasterized=True,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
    )
    _style_axis(ax, title)
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="Query confidence")


def summarize_dotplot(adata: ad.AnnData, group_key: str, gene_order: Sequence[str]) -> pd.DataFrame:
    frame = adata.obs[[group_key]].copy()
    frame[group_key] = _to_clean_string(frame[group_key])
    X = adata.X
    if sparse.issparse(X):
        X = X.tocsr()
    summaries: list[pd.DataFrame] = []
    for group in frame[group_key].drop_duplicates().tolist():
        idx = np.flatnonzero(frame[group_key].to_numpy() == group)
        sub = X[idx]
        if sparse.issparse(sub):
            mean_expr = np.asarray(sub.mean(axis=0)).ravel()
            pct_expr = np.asarray((sub > 0).mean(axis=0)).ravel()
        else:
            mean_expr = np.asarray(sub.mean(axis=0)).ravel()
            pct_expr = np.asarray((sub > 0).mean(axis=0)).ravel()
        summaries.append(
            pd.DataFrame(
                {
                    "group": group,
                    "gene": list(gene_order),
                    "mean_expr": mean_expr,
                    "pct_expr": pct_expr,
                }
            )
        )
    summary = pd.concat(summaries, ignore_index=True)
    summary["mean_expr_scaled"] = summary.groupby("gene", observed=False)["mean_expr"].transform(_min_max_scale)
    return summary


def _min_max_scale(values: pd.Series) -> pd.Series:
    values = pd.Series(values, copy=False)
    vmin = float(values.min())
    vmax = float(values.max())
    if math.isclose(vmin, vmax):
        return pd.Series(np.full(len(values), 0.5), index=values.index)
    return (values - vmin) / (vmax - vmin)


def plot_dotplot(
    adata: ad.AnnData,
    marker_dict: Mapping[str, Sequence[str]],
    output_dir: Path,
    title: str,
    group_order: Sequence[str] | None = None,
) -> pd.DataFrame:
    output_dir = ensure_dir(output_dir)
    gene_order, bounds = flatten_marker_dict(marker_dict)
    summary = summarize_dotplot(adata, "source_label", gene_order)
    groups = list(group_order) if group_order is not None else summary["group"].drop_duplicates().tolist()
    groups = [group for group in groups if group in set(summary["group"])]

    x_map = {gene: idx for idx, gene in enumerate(gene_order)}
    y_map = {group: idx for idx, group in enumerate(groups)}

    fig_width = max(11, len(gene_order) * 0.55)
    fig_height = max(5.5, len(groups) * 0.45 + 1.8)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    plot_df = summary[summary["group"].isin(groups)].copy()
    x = plot_df["gene"].map(x_map).to_numpy()
    y = plot_df["group"].map(y_map).to_numpy()
    sizes = 30 + plot_df["pct_expr"].to_numpy() * 330
    colors = plot_df["mean_expr_scaled"].to_numpy()
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap="magma", vmin=0.0, vmax=1.0, edgecolors="none")

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=60, ha="right", fontsize=9)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups, fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", color="#EFEFEF", linewidth=0.7)
    ax.set_title(title, fontsize=14, fontweight="bold")
    for spine in ["top", "right", "left", "bottom"]:
        ax.spines[spine].set_visible(False)

    for _, start, end in bounds:
        if start > 0:
            ax.axvline(start - 0.5, color="#DDDDDD", linewidth=1)
    for group_name, start, end in bounds:
        if end <= start:
            continue
        center = (start + end - 1) / 2
        ax.text(center, -1.15, group_name, ha="center", va="bottom", fontsize=10, fontweight="bold")

    colorbar = plt.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Scaled mean expression", fontsize=10)

    size_levels = [0.1, 0.4, 0.7, 1.0]
    size_handles = [
        plt.scatter([], [], s=30 + level * 330, color="#666666", edgecolors="none", label=f"{int(level * 100)}%")
        for level in size_levels
    ]
    ax.legend(
        handles=size_handles,
        title="Pct. expressing",
        loc="center left",
        bbox_to_anchor=(1.14, 0.5),
        frameon=False,
        fontsize=9,
        title_fontsize=10,
    )

    fig.tight_layout()
    for suffix in ["pdf", "png"]:
        fig.savefig(output_dir / f"joint_dotplot.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    summary.to_csv(output_dir / "joint_dotplot_summary.tsv", sep="\t", index=False)
    return summary


def plot_feature_grid(
    adata: ad.AnnData,
    genes: Sequence[str],
    output_dir: Path,
    title: str,
) -> list[str]:
    output_dir = ensure_dir(output_dir)
    available = [gene for gene in genes if gene in set(adata.var_names.astype(str))]
    if not available:
        raise ValueError("No requested feature genes are available for plotting")

    coords = np.asarray(adata.obsm["X_umap"], dtype=np.float32)
    X = adata.X.toarray().astype(np.float32) if sparse.issparse(adata.X) else np.asarray(adata.X, dtype=np.float32)
    point_size = _estimate_point_size(adata.n_obs)
    ncols = 4
    nrows = int(math.ceil(len(available) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.0, nrows * 3.7))
    axes = np.array(axes).reshape(-1)

    for idx, gene in enumerate(available):
        ax = axes[idx]
        expr = X[:, adata.var_names.get_loc(gene)]
        order = np.argsort(expr)
        vmax = float(np.quantile(expr[expr > 0], 0.99)) if np.any(expr > 0) else 1.0
        vmax = max(vmax, 1e-6)
        sc = ax.scatter(
            coords[order, 0],
            coords[order, 1],
            c=expr[order],
            s=point_size,
            linewidths=0,
            cmap=EXPR_CMAP,
            norm=Normalize(vmin=0.0, vmax=vmax),
            rasterized=True,
        )
        _style_axis(ax, gene)
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    for idx in range(len(available), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(title, fontsize=15, fontweight="bold", y=1.01)
    fig.tight_layout()
    for suffix in ["pdf", "png"]:
        fig.savefig(output_dir / f"joint_featureplot.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return available


def plot_umap_overview(
    adata: ad.AnnData,
    output_dir: Path,
    title: str,
    label_palette: Mapping[str, str],
) -> None:
    output_dir = ensure_dir(output_dir)
    coords = np.asarray(adata.obsm["X_umap"], dtype=np.float32)
    point_size = _estimate_point_size(adata.n_obs)
    fig, axes = plt.subplots(1, 3, figsize=(19, 6))

    scatter_categorical(
        axes[0],
        coords,
        adata.obs["data_source"],
        title="Reference vs Query",
        palette=SOURCE_PALETTE,
        point_size=point_size,
        order=["reference", "query"],
    )
    scatter_categorical(
        axes[1],
        coords,
        adata.obs["plot_label"],
        title="Joint label",
        palette=label_palette,
        point_size=point_size,
    )
    scatter_query_confidence(
        axes[2],
        coords,
        adata.obs["data_source"],
        adata.obs.get("mapping_confidence", pd.Series(np.nan, index=adata.obs.index)),
        title="Query confidence",
        point_size=point_size,
    )

    fig.suptitle(title, fontsize=15, fontweight="bold")
    fig.tight_layout()
    for suffix in ["pdf", "png"]:
        fig.savefig(output_dir / f"joint_umap_overview.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_proportion_violin(
    obs: pd.DataFrame,
    output_dir: Path,
    title: str,
    label_order: Sequence[str] | None = None,
) -> pd.DataFrame:
    output_dir = ensure_dir(output_dir)
    props = compute_group_proportions(obs, sample_key="sample", label_key="plot_label", source_key="data_source")
    props = props.copy()
    props["proportion_pct"] = props["proportion"] * 100.0
    if label_order is not None:
        props = props[props["plot_label"].isin(label_order)].copy()
        props["plot_label"] = pd.Categorical(props["plot_label"], categories=list(label_order), ordered=True)

    sns = _lazy_import_seaborn()
    fig_width = max(10, props["plot_label"].nunique() * 1.35)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    violin_kwargs = {
        "data": props,
        "x": "plot_label",
        "y": "proportion_pct",
        "hue": "data_source",
        "palette": SOURCE_PALETTE,
        "cut": 0,
        "inner": None,
        "linewidth": 0.8,
        "ax": ax,
    }
    try:
        sns.violinplot(**violin_kwargs, density_norm="width")
    except TypeError:
        sns.violinplot(**violin_kwargs, scale="width")
    sns.stripplot(
        data=props,
        x="plot_label",
        y="proportion_pct",
        hue="data_source",
        dodge=True,
        palette=SOURCE_PALETTE,
        alpha=0.35,
        size=1.7,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    unique = OrderedDict()
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    ax.legend(unique.values(), unique.keys(), frameon=False, title="data_source")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Cell-type proportion per sample (%)")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", color="#EFEFEF")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    for suffix in ["pdf", "png"]:
        fig.savefig(output_dir / f"joint_proportion_violin.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    props.to_csv(output_dir / "joint_proportion_violin_summary.tsv", sep="\t", index=False)
    return props


def compute_joint_umap_from_latent(reference_latent: np.ndarray, query_latent: np.ndarray) -> np.ndarray:
    sc = _lazy_import_scanpy()
    latent = np.vstack([reference_latent, query_latent]).astype(np.float32)
    temp = ad.AnnData(obs=pd.DataFrame(index=np.arange(latent.shape[0])))
    temp.obsm["X_latent"] = latent
    sc.pp.neighbors(temp, use_rep="X_latent", n_neighbors=50)
    sc.tl.umap(temp, min_dist=0.5, spread=1.0, random_state=42)
    return np.asarray(temp.obsm["X_umap"], dtype=np.float32)


def build_marker_plot_adata(
    ref_adata: ad.AnnData,
    qry_adata: ad.AnnData,
    ref_meta: pd.DataFrame,
    qry_meta: pd.DataFrame,
    ref_coords: np.ndarray,
    qry_coords: np.ndarray,
    marker_dict: Mapping[str, Sequence[str]],
) -> tuple[ad.AnnData, pd.DataFrame]:
    marker_genes, _ = flatten_marker_dict(marker_dict)
    ref_expr, ref_info = extract_expression_matrix(ref_adata, marker_genes, prefer_counts=True)
    qry_expr, qry_info = extract_expression_matrix(qry_adata, marker_genes, prefer_counts=True)

    obs = pd.concat([ref_meta, qry_meta], axis=0)
    obs = obs.copy()
    obs["source_label"] = obs["data_source"].astype(str) + " · " + obs["plot_label"].astype(str)
    X = np.vstack([ref_expr, qry_expr]).astype(np.float32)
    plot_adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(marker_genes, name="gene")))
    plot_adata.obsm["X_umap"] = np.vstack([ref_coords, qry_coords]).astype(np.float32)

    availability_rows = []
    for gene in marker_genes:
        availability_rows.append(
            {
                "gene": gene,
                "in_reference": gene in set(ref_info["present"]),
                "in_query": gene in set(qry_info["present"]),
            }
        )
    availability = pd.DataFrame(availability_rows)
    return plot_adata, availability


def _prepare_reference_metadata(adata: ad.AnnData, label_key: str, sample_key: str, tissue_key: str) -> pd.DataFrame:
    meta = pd.DataFrame(index=adata.obs_names)
    meta["sample"] = _to_clean_string(adata.obs.get(sample_key))
    meta["tissue"] = _to_clean_string(adata.obs.get(tissue_key))
    meta["data_source"] = "reference"
    meta["plot_label"] = _to_clean_string(adata.obs.get(label_key))
    meta["mapping_confidence"] = np.nan
    return meta


def _prepare_query_metadata(
    adata: ad.AnnData,
    sample_key: str,
    tissue_key: str,
    final_key: str,
    pred_key: str,
    conf_key: str,
) -> pd.DataFrame:
    meta = pd.DataFrame(index=adata.obs_names)
    meta["sample"] = _to_clean_string(adata.obs.get(sample_key))
    meta["tissue"] = _to_clean_string(adata.obs.get(tissue_key))
    meta["data_source"] = "query"
    meta["plot_label"] = _to_clean_string(adata.obs.get(final_key)).replace({"Unknown": np.nan})
    meta["plot_label"] = meta["plot_label"].fillna(_to_clean_string(adata.obs.get(pred_key)))
    meta["plot_label"] = _to_clean_string(meta["plot_label"])
    meta["mapping_confidence"] = pd.to_numeric(adata.obs.get(conf_key), errors="coerce")
    return meta


def load_allcell_reference_strict(cfg: Mapping[str, object]) -> tuple[ad.AnnData, pd.DataFrame, np.ndarray, np.ndarray | None, dict]:
    import scvi

    reference = ad.read_h5ad(cfg["reference_input_h5ad"])
    ensure_counts_layer(reference)

    whitelist = collect_cell_ids_from_folder(cfg["whitelist_dir"])
    mask = reference.obs_names.astype(str).isin(whitelist)
    n_whitelist_overlap = int(mask.sum())
    if n_whitelist_overlap > 0:
        reference = reference[mask].copy()
        reference_subset_strategy = "exact_obsname_whitelist_overlap"
    else:
        warnings.warn(
            "No obs_names overlap between the current whitelist AnnData and the historical 1128 reference input. "
            "Falling back to the full 1128 reference input with the 0115 legacy model so the visualization remains tied "
            "to the correct historical model lineage, even though the exact whitelist subset is no longer recoverable "
            "from the current workspace contents."
        )
        reference = reference.copy()
        reference_subset_strategy = "fallback_full_1128_input_zero_whitelist_overlap"

    add_allcell_covariates_fast(reference)

    hvg_genes = pd.read_csv(cfg["hvg_file"], header=None)[0].astype(str).tolist()
    hvg_genes = [gene for gene in hvg_genes if gene in set(reference.var_names.astype(str))]
    if not hvg_genes:
        raise ValueError("No saved HVG genes overlap the reconstructed allcell reference")

    ref_indices = reference.var_names.get_indexer(hvg_genes)
    subset_counts = reference.layers["counts"][:, ref_indices].copy()
    if not sparse.isspmatrix_csr(subset_counts):
        subset_counts = sparse.csr_matrix(subset_counts)

    model_obs = reference.obs[["sample", "tissue", "pct_counts_mt", "stress_score", "S_score", "G2M_score"]].copy()
    model_obs["sample"] = pd.Categorical(_to_clean_string(model_obs["sample"], unknown="Unknown"))
    model_obs["tissue"] = pd.Categorical(_to_clean_string(model_obs["tissue"], unknown="Unknown"))

    obs = model_obs.copy()
    labels = pd.Categorical(reference.obs[cfg["reference_label_key"]].astype(str))
    if "Unknown" not in labels.categories:
        labels = labels.add_categories(["Unknown"])
    obs["scanvi_label_existing_cleaned"] = labels

    scanvi_adata = ad.AnnData(
        X=subset_counts,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(hvg_genes)),
    )

    scvi_adata = ad.AnnData(
        X=subset_counts,
        obs=model_obs,
        var=pd.DataFrame(index=pd.Index(hvg_genes)),
    )

    scvi.model.SCANVI.setup_anndata(
        scanvi_adata,
        layer=None,
        batch_key="sample",
        labels_key="scanvi_label_existing_cleaned",
        unlabeled_category="Unknown",
        categorical_covariate_keys=["tissue"],
        continuous_covariate_keys=["pct_counts_mt", "stress_score", "S_score", "G2M_score"],
    )

    scvi_model = load_existing_scvi_model(scvi_adata, scvi, Path(cfg["hvg_file"]).parent)
    scanvi_model = load_scanvi_from_checkpoint(scanvi_adata, scvi_model, scvi, Path(cfg["scanvi_model_dir"]))
    ref_latent = np.asarray(scanvi_model.get_latent_representation(), dtype=np.float32)
    model_predictions = pd.Series(scanvi_model.predict(), index=reference.obs_names, dtype=object)

    ref_umap = None
    if hasattr(scanvi_model, "umap_op_"):
        try:
            ref_umap = np.asarray(scanvi_model.umap_op_.transform(ref_latent), dtype=np.float32)
        except Exception as exc:
            warnings.warn(f"Reference UMAP operator transform failed, will fall back to recomputing joint UMAP: {exc}")

    meta = _prepare_reference_metadata(reference, cfg["reference_label_key"], cfg["sample_key"], cfg["tissue_key"])
    meta["model_prediction"] = _to_clean_string(model_predictions)
    details = {
        "n_reference_cells": int(reference.n_obs),
        "n_whitelist_cells": int(len(whitelist)),
        "n_whitelist_overlap": n_whitelist_overlap,
        "reference_subset_strategy": reference_subset_strategy,
        "n_hvg": int(len(hvg_genes)),
        "covariate_strategy": "regenerated_from_counts_fast_signature_mode",
        "has_umap_operator": bool(ref_umap is not None),
    }
    return reference, meta, ref_latent, ref_umap, details


def load_allcell_query(cfg: Mapping[str, object]) -> tuple[ad.AnnData, pd.DataFrame, np.ndarray, np.ndarray | None]:
    query = ad.read_h5ad(cfg["query_h5ad"])
    ensure_counts_layer(query)
    latent_key = cfg.get("latent_key_query")
    if latent_key not in query.obsm:
        latent_key = resolve_umap_key(query.obsm.keys(), preferred="X_scANVI_mapped")
    query_latent = np.asarray(query.obsm[latent_key], dtype=np.float32)
    query_umap = None
    if "X_umap_mapped" in query.obsm:
        query_umap = np.asarray(query.obsm["X_umap_mapped"], dtype=np.float32)
    meta = _prepare_query_metadata(
        query,
        sample_key=cfg["sample_key"],
        tissue_key=cfg["tissue_key"],
        final_key=cfg["query_final_key"],
        pred_key=cfg["query_pred_key"],
        conf_key=cfg["query_conf_key"],
    )
    return query, meta, query_latent, query_umap


def inspect_bcell_candidate(path: Path, marker_genes: Sequence[str]) -> dict:
    candidate = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        candidate.update({"has_data_source": False, "has_umap": False, "has_query_final": False, "marker_hits": 0})
        return candidate

    adata = ad.read_h5ad(path, backed="r")
    try:
        obs_cols = set(map(str, adata.obs.columns))
        obsm_keys = set(map(str, adata.obsm.keys()))
        var_names = set(map(str, adata.var_names))
        candidate.update(
            {
                "has_data_source": "data_source" in obs_cols,
                "has_umap": any(key in obsm_keys for key in ["X_umap", "X_umap_scanvi"]),
                "has_query_final": "Cell_Type_L2_final" in obs_cols,
                "marker_hits": int(sum(gene in var_names for gene in marker_genes)),
            }
        )
        return candidate
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def choose_bcell_input(cfg: Mapping[str, object]) -> Path:
    raw_candidates = cfg.get("candidate_merged_h5ads")
    if raw_candidates is None:
        raw_candidates = [cfg["merged_h5ad"], cfg["fallback_merged_h5ad"]]
    candidates = [
        inspect_bcell_candidate(Path(path), flatten_marker_dict(cfg["marker_dict"])[0])
        for path in raw_candidates
    ]
    best = select_best_bcell_candidate(candidates)
    chosen = Path(best["path"])
    if not chosen.exists():
        raise FileNotFoundError(f"No viable B-cell merged h5ad found among candidates: {candidates}")
    return chosen


def load_bcell_mode(cfg: Mapping[str, object]) -> tuple[ad.AnnData, pd.DataFrame, dict]:
    selected = choose_bcell_input(cfg)
    adata = ad.read_h5ad(selected)
    coords_key = resolve_umap_key(adata, preferred=cfg["umap_key"])
    meta = pd.DataFrame(index=adata.obs_names)
    meta["sample"] = _to_clean_string(adata.obs.get(cfg["sample_key"]))
    meta["tissue"] = _to_clean_string(adata.obs.get(cfg["tissue_key"]))
    meta["data_source"] = _to_clean_string(adata.obs.get(cfg["source_key"])).str.lower()
    meta["plot_label"] = build_plot_label(
        adata.obs,
        source_key=cfg["source_key"],
        ref_label_key=cfg["reference_label_key"],
        query_label_key=cfg["query_final_key"],
        query_fallback_key=cfg["query_pred_key"],
    )
    meta["mapping_confidence"] = pd.to_numeric(adata.obs.get(cfg["query_conf_key"]), errors="coerce")
    meta["source_label"] = meta["data_source"].astype(str) + " · " + meta["plot_label"].astype(str)

    marker_genes, _ = flatten_marker_dict(cfg["marker_dict"])
    expr, _ = extract_expression_matrix(adata, marker_genes, prefer_counts=True)
    plot_adata = ad.AnnData(X=expr.astype(np.float32), obs=meta, var=pd.DataFrame(index=pd.Index(marker_genes, name="gene")))
    plot_adata.obsm["X_umap"] = np.asarray(adata.obsm[coords_key], dtype=np.float32)

    availability = pd.DataFrame({"gene": marker_genes, "in_reference": [gene in set(adata.var_names.astype(str)) for gene in marker_genes], "in_query": [gene in set(adata.var_names.astype(str)) for gene in marker_genes]})
    details = {
        "selected_input": str(selected),
        "n_cells": int(adata.n_obs),
        "coords_key": coords_key,
        "marker_availability": availability,
    }
    return plot_adata, availability, details


def prepare_allcell_mode(cfg: Mapping[str, object]) -> tuple[ad.AnnData, pd.DataFrame, dict]:
    ref_adata, ref_meta, ref_latent, ref_umap, ref_details = load_allcell_reference_strict(cfg)
    qry_adata, qry_meta, qry_latent, qry_umap = load_allcell_query(cfg)

    if ref_umap is not None and qry_umap is not None:
        ref_coords = ref_umap
        qry_coords = qry_umap
        joint_strategy = "model_umap_operator_for_reference_plus_existing_query_mapped_umap"
    else:
        joint_coords = compute_joint_umap_from_latent(ref_latent, qry_latent)
        ref_coords = joint_coords[: ref_adata.n_obs]
        qry_coords = joint_coords[ref_adata.n_obs :]
        joint_strategy = "recomputed_from_exact_reference_and_query_scanvi_latent"

    plot_adata, availability = build_marker_plot_adata(
        ref_adata=ref_adata,
        qry_adata=qry_adata,
        ref_meta=ref_meta,
        qry_meta=qry_meta,
        ref_coords=ref_coords,
        qry_coords=qry_coords,
        marker_dict=cfg["marker_dict"],
    )

    details = {
        **ref_details,
        "n_query_cells": int(qry_adata.n_obs),
        "latent_key_query": cfg["latent_key_query"],
        "joint_umap_strategy": joint_strategy,
    }
    return plot_adata, availability, details


def _clean_for_write(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    for col in adata.obs.columns:
        if pd.api.types.is_string_dtype(adata.obs[col].dtype):
            adata.obs[col] = adata.obs[col].astype(object)
    for col in adata.var.columns:
        if pd.api.types.is_string_dtype(adata.var[col].dtype):
            adata.var[col] = adata.var[col].astype(object)
    return adata


def save_lightweight_plot_adata(adata: ad.AnnData, output_dir: Path, mode: str) -> Path:
    output_dir = ensure_dir(output_dir)
    target = output_dir / f"{mode}_joint_marker_plot_data.h5ad"
    writable = _clean_for_write(adata)
    writable.write_h5ad(target, compression="gzip")
    return target


def write_mode_summary(
    mode: str,
    output_dir: Path,
    plot_adata: ad.AnnData,
    availability: pd.DataFrame,
    details: Mapping[str, object],
) -> None:
    output_dir = ensure_dir(output_dir)
    source_counts = (
        plot_adata.obs.groupby(["data_source", "plot_label"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    source_counts.to_csv(output_dir / "source_label_counts.tsv", sep="\t", index=False)
    availability.to_csv(output_dir / "marker_availability.tsv", sep="\t", index=False)

    summary = {
        "mode": mode,
        "n_cells": int(plot_adata.n_obs),
        "n_markers": int(plot_adata.n_vars),
        "n_reference": int((plot_adata.obs["data_source"].astype(str) == "reference").sum()),
        "n_query": int((plot_adata.obs["data_source"].astype(str) == "query").sum()),
        "available_markers": int((availability[["in_reference", "in_query"]].any(axis=1)).sum()),
    }
    summary.update(_make_json_safe(details))
    with open(output_dir / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)


def make_group_order(marker_dict: Mapping[str, Sequence[str]], data_sources: Sequence[str] = ("reference", "query")) -> list[str]:
    order = []
    for label in marker_dict.keys():
        for source in data_sources:
            order.append(f"{source} · {label}")
    return order


def run_mode(mode: str) -> dict:
    configure_runtime()
    cfg = get_mode_config(mode)
    output_dir = ensure_dir(Path(cfg["output_dir"]))

    if mode == "allcell":
        plot_adata, availability, details = prepare_allcell_mode(cfg)
    elif mode == "bcell":
        plot_adata, availability, details = load_bcell_mode(cfg)
    else:
        raise KeyError(f"Unsupported mode: {mode}")

    group_order = [group for group in make_group_order(cfg["marker_dict"]) if group in set(plot_adata.obs["source_label"].astype(str))]
    label_order = [label for label in cfg["marker_dict"].keys() if label in set(plot_adata.obs["plot_label"].astype(str))]

    plot_umap_overview(
        plot_adata,
        output_dir=output_dir,
        title=f"{mode} joint reference-query UMAP",
        label_palette=cfg["label_palette"],
    )
    dotplot_summary = plot_dotplot(
        plot_adata,
        marker_dict=cfg["marker_dict"],
        output_dir=output_dir,
        title=f"{mode} joint dotplot",
        group_order=group_order,
    )
    feature_genes = plot_feature_grid(
        plot_adata,
        genes=cfg["feature_genes"],
        output_dir=output_dir,
        title=f"{mode} joint feature plots",
    )
    props = plot_proportion_violin(
        plot_adata.obs,
        output_dir=output_dir,
        title=f"{mode} sample-level label proportions",
        label_order=label_order,
    )
    plot_data_path = save_lightweight_plot_adata(plot_adata, output_dir=output_dir, mode=mode)
    write_mode_summary(mode, output_dir=output_dir, plot_adata=plot_adata, availability=availability, details=details)

    return {
        "mode": mode,
        "output_dir": str(output_dir),
        "plot_data_path": str(plot_data_path),
        "dotplot_rows": int(len(dotplot_summary)),
        "n_feature_genes": int(len(feature_genes)),
        "proportion_rows": int(len(props)),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Joint reference-query visualization")
    parser.add_argument("--mode", choices=["allcell", "bcell", "both"], default="both")
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    modes = ["allcell", "bcell"] if args.mode == "both" else [args.mode]
    results = []
    for mode in modes:
        print(f"[run] {mode}")
        result = run_mode(mode)
        results.append(result)
        print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()