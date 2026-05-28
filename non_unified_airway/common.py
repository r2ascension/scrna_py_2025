from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

DEFAULT_REPO_ROOT = Path("/home/h2048")
DEFAULT_METHODS_DOC = DEFAULT_REPO_ROOT / "docs/methods/nature_methods_single_cell_respiratory_tract_20260506.md"
DEFAULT_ALLCELLS_REFERENCE_CANDIDATES = [
    DEFAULT_REPO_ROOT / "data/py/20260319/allcells_combined_scanvi/allcells_combined_20260319_v1_10_L1.h5ad",
    DEFAULT_REPO_ROOT / "data/core20260115/adata_cleaned_no_doublets_no_GSE299751.h5ad",
    DEFAULT_REPO_ROOT / "data/R/1215/merge/cleaned_samples_COMPLETE_v4.1/merged_seurat_standardized_filterd.h5ad",
]
DEFAULT_MAPPED_QUERY_CANDIDATES = [
    DEFAULT_REPO_ROOT / "data/py/0127/scarches_mapping_FIXED_v1_2/query_mapped_to_reference.h5ad",
]
CANDIDATE_AUDIT_OBS_COLUMNS = [
    "sample",
    "dataset",
    "study",
    "tissue",
    "tissue_level_2",
    "tissue_sampling_method",
    "condition",
    "disease_level_2",
]
AUDIT_OBS_COLUMNS = [
    "sample",
    "batch",
    "dataset",
    "study",
    "tissue",
    "tissue_level_2",
    "tissue_sampling_method",
    "condition",
    "disease_level_2",
    "lineage_family",
    "lineage_branch",
    "cell_type_L1",
    "cell_type_L2",
    "cell_type_L3",
    "cell_type_scanvi_pred",
    "cell_type",
    "ann_finest_level",
    "schpl_rejected",
    "schpl_novel_candidate",
]
LOW_CONFIDENCE_LABELS = {
    "",
    "na",
    "nan",
    "none",
    "null",
    "unknown",
    "unassigned",
    "unlabeled",
    "rejected",
    "low_confidence",
}
HEALTHY_LABELS = {
    "healthy",
    "control",
    "ctrl",
    "normal",
    "non_disease",
    "non-disease",
    "donor",
}
CASE_LABELS = {
    "case",
    "patient",
    "disease",
    "affected",
}


@dataclass(frozen=True)
class SiteAnnotation:
    site_group: str
    site_axis: str
    analysis_compartment: str


@dataclass(frozen=True)
class PathRecord:
    role: str
    path: str
    exists: bool
    source: str
    label: str


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def slugify_label(value: object) -> str:
    cleaned = _clean_text(value).lower()
    cleaned = re.sub(r"[^0-9a-z]+", "_", cleaned)
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "unknown"


def normalize_disease_group(disease_level_2: object = "", condition: object = "") -> str:
    primary = slugify_label(disease_level_2)
    fallback = slugify_label(condition)
    value = primary if primary != "unknown" else fallback
    if value in HEALTHY_LABELS:
        return "healthy"
    if value in CASE_LABELS:
        return "case"
    return value


def normalize_site(
    tissue: object,
    tissue_level_2: object = "",
    tissue_sampling_method: object = "",
) -> SiteAnnotation:
    tissue_s = slugify_label(tissue)
    joined = " ".join(
        part for part in [_clean_text(tissue_level_2).lower(), _clean_text(tissue_sampling_method).lower()] if part
    )

    if tissue_s in {"nose", "nasal", "nasal_tissue"}:
        return SiteAnnotation("nasal", "upper_airway", "upper_airway")
    if tissue_s == "sinus" or "sinus" in tissue_s:
        return SiteAnnotation("sinus", "upper_airway", "upper_airway")
    if tissue_s == "respiratory_airway":
        if "trache" in joined:
            return SiteAnnotation("trachea", "lower_airway", "lower_airway")
        if "bronch" in joined:
            return SiteAnnotation("bronchial", "lower_airway", "lower_airway")
        return SiteAnnotation("lower_conducting_airway", "lower_airway", "lower_airway")
    if tissue_s in {"lung_parenchyma", "lung"} or "parench" in tissue_s:
        return SiteAnnotation("distal_lung", "lower_airway", "lung_parenchyma")
    return SiteAnnotation("unknown", "unknown", "unknown")


def _resolve_confidence_group(row: Mapping[str, object]) -> str:
    for column in ["cell_type_L3", "ann_finest_level", "cell_type_scanvi_pred", "cell_type_L2", "cell_type"]:
        if column not in row:
            continue
        value = slugify_label(row[column])
        if value == "unknown":
            continue
        if value in LOW_CONFIDENCE_LABELS:
            return "low_confidence"
        return "assigned"
    return "unknown"


def add_canonical_fields(metadata: pd.DataFrame) -> pd.DataFrame:
    out = metadata.copy()
    for column in [
        "sample",
        "batch",
        "dataset",
        "study",
        "tissue",
        "tissue_level_2",
        "tissue_sampling_method",
        "condition",
        "disease_level_2",
    ]:
        if column not in out.columns:
            out[column] = ""

    annotations = [
        normalize_site(tissue, tissue_level_2, tissue_sampling_method)
        for tissue, tissue_level_2, tissue_sampling_method in zip(
            out["tissue"],
            out["tissue_level_2"],
            out["tissue_sampling_method"],
        )
    ]
    out["site_group"] = [annotation.site_group for annotation in annotations]
    out["site_axis"] = [annotation.site_axis for annotation in annotations]
    out["analysis_compartment"] = [annotation.analysis_compartment for annotation in annotations]
    out["disease_group"] = [
        normalize_disease_group(disease_level_2, condition)
        for disease_level_2, condition in zip(out["disease_level_2"], out["condition"])
    ]
    out["is_healthy"] = out["disease_group"].eq("healthy")
    out["is_upper_airway"] = out["analysis_compartment"].eq("upper_airway")
    out["is_lower_airway"] = out["analysis_compartment"].isin(["lower_airway", "lung_parenchyma"])
    out["is_lung_parenchyma"] = out["analysis_compartment"].eq("lung_parenchyma")
    out["confidence_group"] = [
        _resolve_confidence_group(record)
        for record in out.to_dict(orient="records")
    ]
    out["tissue_site_disease"] = [
        "__".join(
            [
                _clean_text(tissue).lower() or "unknown",
                site_group,
                disease_group,
            ]
        )
        for tissue, site_group, disease_group in zip(out["tissue"], out["site_group"], out["disease_group"])
    ]
    out["site_disease"] = out["site_group"].astype(str) + "__" + out["disease_group"].astype(str)
    return out


def build_sample_manifest(metadata: pd.DataFrame) -> pd.DataFrame:
    required = {"sample", "site_group", "disease_group", "tissue_site_disease"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"metadata missing required columns for sample manifest: {missing}")

    out = metadata.copy()
    for column in ["dataset", "study", "batch", "tissue", "tissue_level_2", "tissue_sampling_method"]:
        if column not in out.columns:
            out[column] = ""

    group_cols = [
        "sample",
        "dataset",
        "study",
        "batch",
        "tissue",
        "tissue_level_2",
        "tissue_sampling_method",
        "site_group",
        "site_axis",
        "analysis_compartment",
        "disease_group",
        "tissue_site_disease",
    ]
    group_cols = [column for column in group_cols if column in out.columns]
    agg_spec: dict[str, tuple[str, str]] = {"n_cells": ("sample", "size")}
    if "cell_type_L2" in out.columns:
        agg_spec["n_celltype_l2"] = ("cell_type_L2", "nunique")
    if "cell_type_L3" in out.columns:
        agg_spec["n_celltype_l3"] = ("cell_type_L3", "nunique")
    manifest = (
        out.groupby(group_cols, dropna=False)
        .agg(**agg_spec)
        .reset_index()
        .sort_values([column for column in ["dataset", "study", "sample", "site_group", "disease_group"] if column in group_cols])
        .reset_index(drop=True)
    )
    manifest["is_healthy"] = manifest["disease_group"].eq("healthy")
    manifest["is_upper_airway"] = manifest["analysis_compartment"].eq("upper_airway")
    manifest["is_lower_airway"] = manifest["analysis_compartment"].isin(["lower_airway", "lung_parenchyma"])
    manifest["is_lung_parenchyma"] = manifest["analysis_compartment"].eq("lung_parenchyma")
    return manifest


def _has_meaningful_dataset_column(sample_manifest: pd.DataFrame) -> bool:
    if "dataset" not in sample_manifest.columns:
        return False
    dataset_values = sample_manifest["dataset"].map(slugify_label)
    unique_values = {value for value in dataset_values if value != "unknown"}
    return len(unique_values) >= 2


def _has_within_dataset_variation(sample_manifest: pd.DataFrame, feature: str) -> bool:
    if feature not in sample_manifest.columns:
        return False
    dataset_col = "dataset" if "dataset" in sample_manifest.columns else None
    if dataset_col is None:
        return sample_manifest[feature].nunique(dropna=True) >= 2
    sample_manifest = sample_manifest.copy()
    sample_manifest[dataset_col] = sample_manifest[dataset_col].map(slugify_label)
    for _, sub_df in sample_manifest.groupby(dataset_col, dropna=False):
        if sub_df[feature].nunique(dropna=True) >= 2:
            return True
    return False


def recommend_model_formulas(sample_manifest: pd.DataFrame) -> dict[str, dict[str, str]]:
    required = {"site_group", "disease_group", "tissue_site_disease"}
    missing = sorted(required - set(sample_manifest.columns))
    if missing:
        raise ValueError(f"sample_manifest missing required columns: {missing}")

    result: dict[str, dict[str, str]] = {}
    has_dataset = _has_meaningful_dataset_column(sample_manifest)

    healthy = sample_manifest.loc[sample_manifest["disease_group"] == "healthy"].copy()
    if healthy["site_group"].nunique(dropna=True) < 2:
        result["healthy_mainline"] = {
            "status": "blocked",
            "formula": "",
            "reason": "healthy samples cover fewer than two site groups",
        }
    elif has_dataset and _has_within_dataset_variation(healthy, "site_group"):
        result["healthy_mainline"] = {
            "status": "estimable",
            "formula": "~ dataset + site_group",
            "reason": "at least one dataset contains multiple healthy site groups",
        }
    elif has_dataset:
        result["healthy_mainline"] = {
            "status": "downgraded",
            "formula": "~ tissue_site_disease",
            "reason": "healthy site_group is fully confounded with dataset",
        }
    else:
        result["healthy_mainline"] = {
            "status": "estimable",
            "formula": "~ site_group",
            "reason": "single-dataset or dataset-missing healthy comparison",
        }

    extension = sample_manifest.copy()
    if extension["disease_group"].nunique(dropna=True) < 2:
        result["disease_extension"] = {
            "status": "blocked",
            "formula": "",
            "reason": "fewer than two disease groups available",
        }
    elif has_dataset and _has_within_dataset_variation(extension, "site_group") and _has_within_dataset_variation(extension, "disease_group"):
        result["disease_extension"] = {
            "status": "estimable",
            "formula": "~ dataset + site_group + disease_group",
            "reason": "dataset-adjusted disease and site effects are separately supported",
        }
    elif has_dataset:
        result["disease_extension"] = {
            "status": "downgraded",
            "formula": "~ dataset + tissue_site_disease",
            "reason": "site and disease are not cleanly separable after dataset adjustment",
        }
    else:
        result["disease_extension"] = {
            "status": "estimable",
            "formula": "~ site_group + disease_group",
            "reason": "single-dataset or dataset-missing disease extension",
        }
    return result


def build_recommended_contrasts(sample_manifest: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if sample_manifest.empty:
        return pd.DataFrame(
            columns=[
                "contrast_id",
                "contrast_group",
                "left_sites",
                "right_sites",
                "disease_group",
                "status",
                "n_samples_left",
                "n_samples_right",
            ]
        )

    healthy = sample_manifest.loc[sample_manifest["disease_group"] == "healthy"].copy()
    if not healthy.empty:
        upper_sites = [site for site in ["nasal", "sinus"] if site in set(healthy["site_group"])]
        if upper_sites and "bronchial" in set(healthy["site_group"]):
            rows.append(
                {
                    "contrast_id": "healthy_upper_vs_bronchial",
                    "contrast_group": "healthy_mainline",
                    "left_sites": "+".join(upper_sites),
                    "right_sites": "bronchial",
                    "disease_group": "healthy",
                    "status": "ready",
                    "n_samples_left": int(healthy.loc[healthy["site_group"].isin(upper_sites), "sample"].nunique()),
                    "n_samples_right": int(healthy.loc[healthy["site_group"] == "bronchial", "sample"].nunique()),
                }
            )
        if upper_sites and "distal_lung" in set(healthy["site_group"]):
            rows.append(
                {
                    "contrast_id": "healthy_upper_vs_distal_lung",
                    "contrast_group": "healthy_mainline",
                    "left_sites": "+".join(upper_sites),
                    "right_sites": "distal_lung",
                    "disease_group": "healthy",
                    "status": "ready",
                    "n_samples_left": int(healthy.loc[healthy["site_group"].isin(upper_sites), "sample"].nunique()),
                    "n_samples_right": int(healthy.loc[healthy["site_group"] == "distal_lung", "sample"].nunique()),
                }
            )
        if "bronchial" in set(healthy["site_group"]) and "distal_lung" in set(healthy["site_group"]):
            rows.append(
                {
                    "contrast_id": "healthy_bronchial_vs_distal_lung",
                    "contrast_group": "healthy_mainline",
                    "left_sites": "bronchial",
                    "right_sites": "distal_lung",
                    "disease_group": "healthy",
                    "status": "ready",
                    "n_samples_left": int(healthy.loc[healthy["site_group"] == "bronchial", "sample"].nunique()),
                    "n_samples_right": int(healthy.loc[healthy["site_group"] == "distal_lung", "sample"].nunique()),
                }
            )

    diseased = sample_manifest.loc[sample_manifest["disease_group"] != "healthy"].copy()
    if not diseased.empty:
        upper_airway = diseased.loc[diseased["analysis_compartment"] == "upper_airway"]
        if not upper_airway.empty:
            rows.append(
                {
                    "contrast_id": "upper_airway_disease_specialization",
                    "contrast_group": "disease_extension",
                    "left_sites": "+".join(sorted(set(upper_airway["site_group"]))),
                    "right_sites": "healthy_reference",
                    "disease_group": "+".join(sorted(set(upper_airway["disease_group"]))),
                    "status": "ready",
                    "n_samples_left": int(upper_airway["sample"].nunique()),
                    "n_samples_right": int(sample_manifest.loc[sample_manifest["disease_group"] == "healthy", "sample"].nunique()),
                }
            )
    return pd.DataFrame(rows)


def discover_lineage_output_dirs(methods_text: str, repo_root: Path = DEFAULT_REPO_ROOT) -> list[Path]:
    matches = re.findall(r"`((?:/home/h2048/)?data/R/[^`]+)`", methods_text)
    paths: list[Path] = []
    seen: set[str] = set()
    for match in matches:
        path = Path(match)
        if not path.is_absolute():
            path = repo_root / path
        path_key = str(path)
        if path_key in seen:
            continue
        seen.add(path_key)
        paths.append(path)
    return paths


def discover_input_inventory(
    repo_root: Path = DEFAULT_REPO_ROOT,
    methods_doc_path: Path = DEFAULT_METHODS_DOC,
    allcells_reference_candidates: Sequence[Path] | None = None,
    mapped_query_candidates: Sequence[Path] | None = None,
) -> dict[str, Any]:
    allcells_reference_candidates = list(allcells_reference_candidates or DEFAULT_ALLCELLS_REFERENCE_CANDIDATES)
    mapped_query_candidates = list(mapped_query_candidates or DEFAULT_MAPPED_QUERY_CANDIDATES)
    methods_text = methods_doc_path.read_text(encoding="utf-8")
    lineage_dirs = discover_lineage_output_dirs(methods_text, repo_root=repo_root)
    candidate_audit_df = summarize_allcells_candidates(allcells_reference_candidates)
    selected_allcells_reference = None
    if not candidate_audit_df.empty:
        selected_rows = candidate_audit_df.loc[candidate_audit_df["selected"]]
        if not selected_rows.empty:
            selected_allcells_reference = str(Path(selected_rows.iloc[0]["path"]))

    records: list[PathRecord] = []
    for idx, path in enumerate(allcells_reference_candidates, start=1):
        records.append(
            PathRecord(
                role="allcells_reference",
                path=str(path),
                exists=path.exists(),
                source="hardcoded_priority",
                label=f"reference_candidate_{idx}",
            )
        )
    for idx, path in enumerate(mapped_query_candidates, start=1):
        records.append(
            PathRecord(
                role="mapped_query",
                path=str(path),
                exists=path.exists(),
                source="hardcoded_priority",
                label=f"mapped_query_candidate_{idx}",
            )
        )
    for lineage_dir in lineage_dirs:
        records.append(
            PathRecord(
                role="lineage_output_dir",
                path=str(lineage_dir),
                exists=lineage_dir.exists(),
                source="methods_doc",
                label=lineage_dir.name,
            )
        )

    return {
        "repo_root": str(repo_root),
        "methods_doc_path": str(methods_doc_path),
        "allcells_reference_candidates": [str(path) for path in allcells_reference_candidates],
        "selected_allcells_reference": selected_allcells_reference,
        "allcells_reference_candidate_audit": candidate_audit_df.to_dict(orient="records"),
        "mapped_query_candidates": [str(path) for path in mapped_query_candidates],
        "lineage_output_dirs": [str(path) for path in lineage_dirs],
        "records": [asdict(record) for record in records],
    }


def _decode_array(dataset: h5py.Dataset) -> np.ndarray:
    arr = dataset[()]
    if getattr(arr, "dtype", None) is not None and arr.dtype.kind in {"S", "O"}:
        return np.array(
            [
                value.decode("utf-8", "replace") if isinstance(value, (bytes, bytearray)) else str(value)
                for value in arr
            ],
            dtype=object,
        )
    return arr


def _read_h5ad_column(group: h5py.Group, name: str) -> list[str]:
    obj = group[name]
    if isinstance(obj, h5py.Group):
        if "codes" in obj and "categories" in obj:
            codes = np.asarray(_decode_array(obj["codes"]), dtype=int)
            categories = _decode_array(obj["categories"])
            values: list[str] = []
            for code in codes:
                if code < 0 or code >= len(categories):
                    values.append("")
                else:
                    values.append(str(categories[code]))
            return values
        if "values" in obj:
            return [str(value) for value in _decode_array(obj["values"])]
        raise ValueError(f"Unsupported h5ad encoding for column '{name}': {list(obj.keys())}")
    return [str(value) for value in _decode_array(obj)]


def read_obs_columns_h5py(h5ad_path: Path | str, columns: Sequence[str]) -> pd.DataFrame:
    h5ad_path = Path(h5ad_path)
    with h5py.File(h5ad_path, "r") as handle:
        if "obs" not in handle or "_index" not in handle["obs"]:
            raise ValueError(f"No h5ad obs/_index found in {h5ad_path}")
        obs_group = handle["obs"]
        index = [str(value) for value in _decode_array(obs_group["_index"])]
        data: dict[str, list[str]] = {}
        for column in columns:
            if column in obs_group:
                data[column] = _read_h5ad_column(obs_group, column)
    return pd.DataFrame(data, index=pd.Index(index, name="cell_id"))


def inspect_h5ad_light(h5ad_path: Path | str) -> dict[str, Any]:
    h5ad_path = Path(h5ad_path)
    with h5py.File(h5ad_path, "r") as handle:
        n_obs = len(handle["obs"]["_index"])
        n_vars = len(handle["var"]["_index"])
        obs_columns = sorted(handle["obs"].keys())
        var_columns = sorted(handle["var"].keys())
        layers = sorted(handle["layers"].keys()) if "layers" in handle else []
        obsm = sorted(handle["obsm"].keys()) if "obsm" in handle else []
        uns = sorted(handle["uns"].keys()) if "uns" in handle else []
        raw_var_n = None
        if "raw" in handle and "var" in handle["raw"] and "_index" in handle["raw"]["var"]:
            raw_var_n = len(handle["raw"]["var"]["_index"])
    return {
        "path": str(h5ad_path),
        "n_obs": int(n_obs),
        "n_vars": int(n_vars),
        "obs_columns": obs_columns,
        "var_columns": var_columns,
        "layers": layers,
        "obsm": obsm,
        "uns": uns,
        "raw_var_n": raw_var_n,
        "has_counts_layer": "counts" in layers,
    }


def choose_first_existing_path(candidates: Iterable[str | Path]) -> Path | None:
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return path
    return None


def summarize_allcells_candidate(
    candidate: str | Path,
    *,
    obs_columns: Sequence[str] = CANDIDATE_AUDIT_OBS_COLUMNS,
) -> dict[str, Any]:
    path = Path(candidate)
    summary: dict[str, Any] = {
        "path": str(path),
        "exists": path.exists(),
        "candidate_audit_error": "",
        "healthy_site_groups": "",
    }
    if not path.exists():
        return summary

    try:
        h5ad_info = inspect_h5ad_light(path)
        summary.update(h5ad_info)
        metadata = read_obs_columns_h5py(path, obs_columns)
        annotated = add_canonical_fields(metadata)
        sample_manifest = build_sample_manifest(annotated)
        healthy = sample_manifest.loc[sample_manifest["disease_group"] == "healthy"].copy()
        healthy_sites = sorted(
            {
                site
                for site in healthy["site_group"].dropna().astype(str)
                if site not in {"", "nan", "unknown"}
            }
        )
        summary.update(
            {
                "n_samples": int(sample_manifest["sample"].nunique()) if "sample" in sample_manifest.columns else 0,
                "n_datasets": int(sample_manifest["dataset"].nunique()) if "dataset" in sample_manifest.columns else 0,
                "n_site_groups": int(sample_manifest["site_group"].nunique()) if "site_group" in sample_manifest.columns else 0,
                "n_healthy_samples": int(healthy["sample"].nunique()) if "sample" in healthy.columns else 0,
                "n_healthy_datasets": int(healthy["dataset"].nunique()) if "dataset" in healthy.columns else 0,
                "n_healthy_site_groups": len(healthy_sites),
                "healthy_site_groups": ";".join(healthy_sites),
                "healthy_multisite_ready": len(healthy_sites) >= 2,
            }
        )
    except Exception as exc:  # pragma: no cover - defensive path for large H5AD quirks
        summary["candidate_audit_error"] = f"{type(exc).__name__}: {exc}"
    return summary


def rank_allcells_candidate_summaries(candidate_summaries: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    columns = [
        "candidate_rank",
        "selected",
        "path",
        "exists",
        "has_counts_layer",
        "n_healthy_site_groups",
        "n_healthy_samples",
        "healthy_site_groups",
        "healthy_multisite_ready",
        "n_site_groups",
        "n_samples",
        "n_datasets",
        "n_obs",
        "n_vars",
        "raw_var_n",
        "candidate_audit_error",
        "selection_reason",
    ]
    if not candidate_summaries:
        return pd.DataFrame(columns=columns)

    df = pd.DataFrame(candidate_summaries).copy()
    for bool_col in ["exists", "has_counts_layer", "healthy_multisite_ready"]:
        if bool_col not in df.columns:
            df[bool_col] = False
        df[bool_col] = df[bool_col].fillna(False).astype(bool)
    for text_col in ["healthy_site_groups", "candidate_audit_error"]:
        if text_col not in df.columns:
            df[text_col] = ""
        df[text_col] = df[text_col].fillna("")
    for numeric_col in [
        "n_healthy_site_groups",
        "n_healthy_samples",
        "n_site_groups",
        "n_samples",
        "n_datasets",
        "n_obs",
        "n_vars",
        "raw_var_n",
    ]:
        if numeric_col not in df.columns:
            df[numeric_col] = 0
        df[numeric_col] = pd.to_numeric(df[numeric_col], errors="coerce").fillna(0).astype(int)

    df["candidate_audit_ok"] = df["candidate_audit_error"].eq("")
    df = (
        df.sort_values(
            by=[
                "exists",
                "candidate_audit_ok",
                "n_healthy_site_groups",
                "has_counts_layer",
                "n_healthy_samples",
                "n_site_groups",
                "raw_var_n",
                "n_obs",
            ],
            ascending=[False, False, False, False, False, False, False, False],
        )
        .reset_index(drop=True)
    )
    df["candidate_rank"] = np.arange(1, len(df) + 1)
    df["selected"] = False
    df["selection_reason"] = ""
    selectable_mask = df["exists"] & df["candidate_audit_ok"]
    if selectable_mask.any():
        selected_idx = df.index[selectable_mask][0]
        df.loc[selected_idx, "selected"] = True
        selected_row = df.loc[selected_idx]
        df.loc[selected_idx, "selection_reason"] = (
            "ranked best by healthy-site coverage "
            f"({int(selected_row['n_healthy_site_groups'])} site groups / {int(selected_row['n_healthy_samples'])} samples) "
            f"and counts_layer={bool(selected_row['has_counts_layer'])}"
        )
    elif not df.empty and bool(df.loc[0, "exists"]):
        df.loc[0, "selected"] = True
        df.loc[0, "selection_reason"] = "fallback to first existing candidate because all audits failed"

    return df[[col for col in columns if col in df.columns]]


def summarize_allcells_candidates(
    candidates: Sequence[str | Path],
    *,
    obs_columns: Sequence[str] = CANDIDATE_AUDIT_OBS_COLUMNS,
) -> pd.DataFrame:
    candidate_summaries = [summarize_allcells_candidate(candidate, obs_columns=obs_columns) for candidate in candidates]
    return rank_allcells_candidate_summaries(candidate_summaries)


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def write_tsv(path: Path, table: pd.DataFrame, *, compression: str | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, compression=compression)
