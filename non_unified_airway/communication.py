from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from non_unified_airway.common import (
    DEFAULT_REPO_ROOT,
    LOW_CONFIDENCE_LABELS,
    add_canonical_fields,
    build_sample_manifest,
    normalize_disease_group,
    read_obs_columns_h5py,
    slugify_label,
    write_json,
    write_tsv,
)

DEFAULT_CONFIG_MANIFEST = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/config/config_manifest.json"
DEFAULT_METADATA_TSV = DEFAULT_REPO_ROOT / "data/py/0524/non_unified_airway/metadata_audit/metadata_canonical.tsv.gz"
DEFAULT_COMM_OUT_DIR = DEFAULT_REPO_ROOT / "data/py/0525/non_unified_airway/communication"

COMM_LEVEL_SPECS: dict[str, dict[str, Any]] = {
    "L2": {
        "label": "L2",
        "preferred_columns": ["cell_type_L2", "cell_type_scanvi_pred", "cell_type"],
        "min_cells_per_sample": 20,
        "min_samples_per_group": 3,
    },
    "L3": {
        "label": "L3",
        "preferred_columns": ["cell_type_L3", "ann_finest_level", "cell_type_scanvi_pred", "cell_type"],
        "min_cells_per_sample": 10,
        "min_samples_per_group": 3,
    },
}

METHOD_SCORE_CANDIDATES: dict[str, list[str]] = {
    "liana": ["aggregate_rank", "magnitude_rank", "specificity_rank", "scaled_weight", "lr_means", "score"],
    "cellchat": ["prob", "weight", "count"],
    "cellphonedb": ["mean", "significant_mean", "rank"],
}

LINEAGE_METADATA_COLUMNS: list[str] = [
    "sample",
    "batch",
    "dataset",
    "study",
    "tissue",
    "tissue_level_2",
    "tissue_sampling_method",
    "condition",
    "disease_level_2",
    "COVID_status",
    "cell_type",
    "cell_type_L2",
    "cell_type_L3",
    "cell_type_scanvi_pred",
    "ann_finest_level",
    "barcode",
]

SITE_CONTRAST_DEFINITIONS: list[dict[str, Any]] = [
    {
        "suffix": "upper_vs_bronchial",
        "left_sites": ["nasal", "sinus"],
        "right_sites": ["bronchial"],
        "left_label_base": "upper_airway",
        "right_label_base": "bronchial",
    },
    {
        "suffix": "upper_vs_distal_lung",
        "left_sites": ["nasal", "sinus"],
        "right_sites": ["distal_lung"],
        "left_label_base": "upper_airway",
        "right_label_base": "distal_lung",
    },
    {
        "suffix": "bronchial_vs_distal_lung",
        "left_sites": ["bronchial"],
        "right_sites": ["distal_lung"],
        "left_label_base": "bronchial",
        "right_label_base": "distal_lung",
    },
]


def module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def pick_lineage_expression_h5ad(lineage_output_dir: Path) -> Path:
    fullgene_candidates = sorted(lineage_output_dir.glob("*final_fullgene.h5ad"))
    if fullgene_candidates:
        return fullgene_candidates[0]
    final_candidates = sorted(lineage_output_dir.glob("*final.h5ad"))
    if final_candidates:
        return final_candidates[0]
    raise FileNotFoundError(f"No lineage final h5ad found in {lineage_output_dir}")


def pick_lineage_metadata_h5ad(lineage_output_dir: Path) -> Path:
    final_candidates = sorted(lineage_output_dir.glob("*final.h5ad"))
    if final_candidates:
        return final_candidates[0]
    return pick_lineage_expression_h5ad(lineage_output_dir)


def load_lineage_metadata(config_manifest: Path = DEFAULT_CONFIG_MANIFEST) -> pd.DataFrame:
    config = load_json(config_manifest)
    lineage_dirs = [Path(path) for path in config.get("lineage_output_dirs", [])]
    if not lineage_dirs:
        raise FileNotFoundError(f"No lineage_output_dirs declared in {config_manifest}")

    frames: list[pd.DataFrame] = []
    for lineage_dir in lineage_dirs:
        metadata_h5ad = pick_lineage_metadata_h5ad(lineage_dir)
        expression_h5ad = pick_lineage_expression_h5ad(lineage_dir)
        metadata = read_obs_columns_h5py(metadata_h5ad, LINEAGE_METADATA_COLUMNS)
        metadata.index.name = "cell_id"
        frame = metadata.reset_index()
        frame["source_h5ad"] = str(metadata_h5ad)
        frame["source_expr_h5ad"] = str(expression_h5ad)
        frame["source_lineage"] = lineage_dir.name
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    combined = add_canonical_fields(combined)

    disease_values: list[str] = []
    for row in combined[["disease_level_2", "condition", "COVID_status"]].fillna("").to_dict(orient="records"):
        disease_group = normalize_disease_group(row.get("disease_level_2", ""), row.get("condition", ""))
        if disease_group == "unknown":
            disease_group = normalize_disease_group(row.get("COVID_status", ""), "")
        disease_values.append(disease_group)
    combined["disease_group"] = disease_values
    combined["is_healthy"] = combined["disease_group"].eq("healthy")
    combined["site_disease"] = combined["site_group"].astype(str) + "__" + combined["disease_group"].astype(str)
    combined["tissue_site_disease"] = (
        combined["tissue"].fillna("unknown").astype(str).str.lower()
        + "__"
        + combined["site_group"].astype(str)
        + "__"
        + combined["disease_group"].astype(str)
    )
    return combined


def pick_primary_h5ad(config_manifest: Path = DEFAULT_CONFIG_MANIFEST, override_h5ad: Path | None = None) -> Path:
    if override_h5ad is not None:
        return override_h5ad
    config = load_json(config_manifest)
    for candidate in config.get("allcells_reference_candidates", []):
        path = Path(candidate)
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No existing all-cells reference candidate found in {config_manifest}."
    )


def parse_multi_values(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [item.strip() for item in text.split("+") if item.strip()]


def resolve_annotation_column(metadata: pd.DataFrame, level: str) -> str:
    if level not in COMM_LEVEL_SPECS:
        raise KeyError(f"Unknown communication level: {level}")
    for column in COMM_LEVEL_SPECS[level]["preferred_columns"]:
        if column not in metadata.columns:
            continue
        cleaned = sanitize_comm_celltypes(metadata[column])
        if cleaned.astype(bool).any():
            return column
    raise ValueError(
        f"No usable annotation column found for {level}; tried {COMM_LEVEL_SPECS[level]['preferred_columns']}"
    )


def sanitize_comm_celltypes(values: Sequence[object]) -> pd.Series:
    series = pd.Series(values, copy=False).astype(object)
    cleaned = series.map(lambda x: "" if pd.isna(x) else str(x).strip())
    keep = []
    for value in cleaned:
        slug = slugify_label(value)
        keep.append(bool(value) and slug not in LOW_CONFIDENCE_LABELS and slug != "unknown")
    return pd.Series(np.where(keep, cleaned, ""), index=cleaned.index, dtype=object)


def build_communication_contrasts(
    sample_manifest: pd.DataFrame,
    *,
    min_samples_per_group: int = 3,
    include_disease_extension: bool = True,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def _n_samples(df: pd.DataFrame, mask: pd.Series) -> int:
        if df.empty:
            return 0
        return int(df.loc[mask, "sample"].nunique())

    def _append_site_contrast(
        df: pd.DataFrame,
        contrast_id: str,
        contrast_group: str,
        left_sites: list[str],
        right_sites: list[str],
        left_label: str,
        right_label: str,
        left_disease_groups: str = "",
        right_disease_groups: str = "",
    ) -> None:
        left_mask = df["site_group"].isin(left_sites)
        right_mask = df["site_group"].isin(right_sites)
        n_left = _n_samples(df, left_mask)
        n_right = _n_samples(df, right_mask)
        if min(n_left, n_right) == 0:
            return
        rows.append(
            {
                "contrast_id": contrast_id,
                "contrast_group": contrast_group,
                "comparison_axis": "site_group",
                "left_label": left_label,
                "right_label": right_label,
                "left_sites": "+".join(left_sites),
                "right_sites": "+".join(right_sites),
                "left_disease_groups": left_disease_groups,
                "right_disease_groups": right_disease_groups,
                "analysis_compartments": "upper_airway+lower_airway+lung_parenchyma",
                "status": "ready" if min(n_left, n_right) >= min_samples_per_group else "blocked",
                "n_samples_left": n_left,
                "n_samples_right": n_right,
            }
        )

    healthy = sample_manifest.loc[sample_manifest["disease_group"] == "healthy"].copy()
    healthy_sites = set(healthy["site_group"].dropna().astype(str))
    all_sites = set(sample_manifest["site_group"].dropna().astype(str))

    for spec in SITE_CONTRAST_DEFINITIONS:
        healthy_left_sites = [site for site in spec["left_sites"] if site in healthy_sites]
        healthy_right_sites = [site for site in spec["right_sites"] if site in healthy_sites]
        if healthy_left_sites and healthy_right_sites:
            _append_site_contrast(
                healthy,
                f"healthy_{spec['suffix']}",
                "healthy_mainline",
                healthy_left_sites,
                healthy_right_sites,
                f"{spec['left_label_base']}_healthy",
                f"{spec['right_label_base']}_healthy",
                left_disease_groups="healthy",
                right_disease_groups="healthy",
            )
            continue

        all_left_sites = [site for site in spec["left_sites"] if site in all_sites]
        all_right_sites = [site for site in spec["right_sites"] if site in all_sites]
        if all_left_sites and all_right_sites:
            _append_site_contrast(
                sample_manifest,
                f"all_samples_{spec['suffix']}",
                "site_group_all_samples",
                all_left_sites,
                all_right_sites,
                f"{spec['left_label_base']}_all_samples",
                f"{spec['right_label_base']}_all_samples",
            )

    if include_disease_extension:
        upper_airway = sample_manifest.loc[sample_manifest["analysis_compartment"] == "upper_airway"].copy()
        available_upper_sites = sorted(set(upper_airway["site_group"].dropna().astype(str)) - {"unknown"})
        disease_groups = sorted(
            set(upper_airway["disease_group"].dropna().astype(str)) - {"healthy", "unknown", "", "nan"}
        )
        for disease in disease_groups:
            left_mask = upper_airway["disease_group"] == disease
            right_mask = upper_airway["disease_group"] == "healthy"
            n_left = _n_samples(upper_airway, left_mask)
            n_right = _n_samples(upper_airway, right_mask)
            rows.append(
                {
                    "contrast_id": f"upper_airway_{slugify_label(disease)}_vs_healthy",
                    "contrast_group": "disease_extension",
                    "comparison_axis": "disease_group",
                    "left_label": f"upper_airway_{disease}",
                    "right_label": "upper_airway_healthy",
                    "left_sites": "+".join(available_upper_sites),
                    "right_sites": "+".join(available_upper_sites),
                    "left_disease_groups": disease,
                    "right_disease_groups": "healthy",
                    "analysis_compartments": "upper_airway",
                    "status": "ready" if min(n_left, n_right) >= min_samples_per_group else "blocked",
                    "n_samples_left": n_left,
                    "n_samples_right": n_right,
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(
            columns=[
                "contrast_id",
                "contrast_group",
                "comparison_axis",
                "left_label",
                "right_label",
                "left_sites",
                "right_sites",
                "left_disease_groups",
                "right_disease_groups",
                "analysis_compartments",
                "status",
                "n_samples_left",
                "n_samples_right",
            ]
        )
    return out.sort_values(["contrast_group", "contrast_id"]).reset_index(drop=True)


def apply_communication_contrast(metadata: pd.DataFrame, contrast_row: Mapping[str, Any]) -> pd.DataFrame:
    left_sites = parse_multi_values(contrast_row.get("left_sites", ""))
    right_sites = parse_multi_values(contrast_row.get("right_sites", ""))
    left_diseases = parse_multi_values(contrast_row.get("left_disease_groups", ""))
    right_diseases = parse_multi_values(contrast_row.get("right_disease_groups", ""))
    compartments = parse_multi_values(contrast_row.get("analysis_compartments", ""))

    def _mask(site_values: list[str], disease_values: list[str]) -> pd.Series:
        mask = pd.Series(True, index=metadata.index)
        if site_values:
            mask &= metadata["site_group"].astype(str).isin(site_values)
        if disease_values:
            mask &= metadata["disease_group"].astype(str).isin(disease_values)
        if compartments:
            mask &= metadata["analysis_compartment"].astype(str).isin(compartments)
        return mask

    left_mask = _mask(left_sites, left_diseases)
    right_mask = _mask(right_sites, right_diseases)
    subset = metadata.loc[left_mask | right_mask].copy()
    if subset.empty:
        return subset
    subset["contrast_id"] = str(contrast_row.get("contrast_id", "unknown_contrast"))
    subset["contrast_group"] = str(contrast_row.get("contrast_group", "unknown_group"))
    subset["comparison_axis"] = str(contrast_row.get("comparison_axis", "unknown_axis"))
    subset["contrast_side"] = np.where(left_mask.loc[subset.index], "left", "right")
    subset["contrast_label"] = np.where(
        subset["contrast_side"].eq("left"),
        str(contrast_row.get("left_label", "left")),
        str(contrast_row.get("right_label", "right")),
    )
    return subset


def build_sample_celltype_counts(subset: pd.DataFrame, *, celltype_col: str) -> pd.DataFrame:
    if subset.empty:
        return pd.DataFrame(
            columns=[
                "contrast_id",
                "contrast_side",
                "sample",
                "site_group",
                "disease_group",
                "comm_celltype",
                "n_cells",
            ]
        )
    out = (
        subset.groupby(
            ["contrast_id", "contrast_side", "sample", "site_group", "disease_group", celltype_col],
            dropna=False,
        )
        .size()
        .reset_index(name="n_cells")
        .rename(columns={celltype_col: "comm_celltype"})
        .sort_values(["contrast_side", "comm_celltype", "sample"])
        .reset_index(drop=True)
    )
    return out


def summarize_celltype_support(
    subset: pd.DataFrame,
    *,
    celltype_col: str,
    min_cells_per_sample: int,
    min_samples_per_group: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_counts = build_sample_celltype_counts(subset, celltype_col=celltype_col)
    if sample_counts.empty:
        empty = pd.DataFrame(
            columns=[
                "comm_celltype",
                "n_samples_left_support",
                "n_samples_right_support",
                "n_samples_left_total",
                "n_samples_right_total",
                "n_cells_left_total",
                "n_cells_right_total",
                "eligible",
            ]
        )
        return empty, sample_counts

    sample_counts["sample_meets_threshold"] = sample_counts["n_cells"] >= int(min_cells_per_sample)
    support = (
        sample_counts.groupby(["contrast_side", "comm_celltype"], dropna=False)
        .agg(
            n_samples_total=("sample", "nunique"),
            n_cells_total=("n_cells", "sum"),
        )
        .reset_index()
    )
    support_meeting = (
        sample_counts.loc[sample_counts["sample_meets_threshold"]]
        .groupby(["contrast_side", "comm_celltype"], dropna=False)
        .agg(n_samples_support=("sample", "nunique"))
        .reset_index()
    )
    merged = support.merge(support_meeting, on=["contrast_side", "comm_celltype"], how="left")
    merged["n_samples_support"] = merged["n_samples_support"].fillna(0).astype(int)

    left = merged.loc[merged["contrast_side"] == "left"].drop(columns="contrast_side")
    right = merged.loc[merged["contrast_side"] == "right"].drop(columns="contrast_side")
    summary = left.merge(right, on="comm_celltype", how="outer", suffixes=("_left", "_right")).fillna(0)
    int_cols = [col for col in summary.columns if col != "comm_celltype"]
    summary[int_cols] = summary[int_cols].astype(int)
    summary["eligible"] = (
        summary["n_samples_support_left"].ge(int(min_samples_per_group))
        & summary["n_samples_support_right"].ge(int(min_samples_per_group))
    )
    summary = summary.rename(
        columns={
            "n_samples_support_left": "n_samples_left_support",
            "n_samples_support_right": "n_samples_right_support",
            "n_samples_total_left": "n_samples_left_total",
            "n_samples_total_right": "n_samples_right_total",
            "n_cells_total_left": "n_cells_left_total",
            "n_cells_total_right": "n_cells_right_total",
        }
    )
    summary = summary.sort_values(["eligible", "comm_celltype"], ascending=[False, True]).reset_index(drop=True)
    return summary, sample_counts


def build_sender_receiver_pairs(eligible_celltypes: Sequence[str], *, include_self: bool = True) -> pd.DataFrame:
    values = [str(x) for x in eligible_celltypes if str(x).strip()]
    rows = []
    for sender in values:
        for receiver in values:
            if not include_self and sender == receiver:
                continue
            rows.append(
                {
                    "sender": sender,
                    "receiver": receiver,
                    "pair_id": f"{sender}__to__{receiver}",
                    "self_interaction": sender == receiver,
                }
            )
    return pd.DataFrame(rows)


def prepare_level_subset(
    metadata: pd.DataFrame,
    contrast_row: Mapping[str, Any],
    *,
    level: str,
) -> tuple[pd.DataFrame, str]:
    subset = apply_communication_contrast(metadata, contrast_row)
    if subset.empty:
        return subset, ""
    celltype_col = resolve_annotation_column(subset, level)
    subset = subset.copy()
    subset["comm_celltype"] = sanitize_comm_celltypes(subset[celltype_col])
    subset = subset.loc[subset["comm_celltype"].astype(str).str.len() > 0].copy()
    return subset, celltype_col


def downsample_for_communication_export(
    subset_metadata: pd.DataFrame,
    *,
    max_cells_per_side_celltype: int | None,
    sample_col: str = "sample",
    random_state: int = 0,
) -> pd.DataFrame:
    if subset_metadata.empty or max_cells_per_side_celltype is None:
        return subset_metadata
    max_cells = int(max_cells_per_side_celltype)
    if max_cells <= 0:
        return subset_metadata

    selected_parts: list[pd.DataFrame] = []
    group_cols = ["contrast_side", "comm_celltype"]
    for _, group_df in subset_metadata.groupby(group_cols, sort=False, dropna=False):
        if group_df.shape[0] <= max_cells:
            selected_parts.append(group_df)
            continue

        sample_groups = [sample_df for _, sample_df in group_df.groupby(sample_col, sort=False, dropna=False)]
        if not sample_groups:
            selected_parts.append(group_df.sample(n=max_cells, random_state=random_state))
            continue

        per_sample_quota = max(1, max_cells // len(sample_groups))
        picked_parts: list[pd.DataFrame] = []
        picked_index: set[Any] = set()
        for sample_df in sample_groups:
            take_n = min(sample_df.shape[0], per_sample_quota)
            taken = sample_df if sample_df.shape[0] <= take_n else sample_df.sample(n=take_n, random_state=random_state)
            picked_parts.append(taken)
            picked_index.update(taken.index.tolist())

        picked = pd.concat(picked_parts, axis=0) if picked_parts else group_df.iloc[0:0].copy()
        remaining_budget = max_cells - picked.shape[0]
        if remaining_budget > 0:
            leftovers = group_df.loc[~group_df.index.isin(picked_index)]
            if not leftovers.empty:
                extra_n = min(remaining_budget, leftovers.shape[0])
                extra = leftovers if leftovers.shape[0] <= extra_n else leftovers.sample(n=extra_n, random_state=random_state)
                picked = pd.concat([picked, extra], axis=0)

        selected_parts.append(picked)

    if not selected_parts:
        return subset_metadata.iloc[0:0].copy()
    out = pd.concat(selected_parts, axis=0)
    out = out.loc[~out.index.duplicated(keep="first")].copy()
    return out.sort_index()


def _coerce_dataframe_strings_for_h5ad(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for column in out.columns:
        series = out[column]
        dtype = series.dtype
        if not (
            pd.api.types.is_object_dtype(dtype)
            or pd.api.types.is_string_dtype(dtype)
            or pd.api.types.is_categorical_dtype(dtype)
        ):
            continue
        values = [pd.NA if pd.isna(value) else str(value) for value in series.astype(object)]
        out[column] = pd.Series(values, index=series.index, dtype="string")
    return out


def _prepare_anndata_for_h5ad_write(adata: "ad.AnnData") -> "ad.AnnData":
    adata = adata.copy()
    adata.raw = None
    adata.obs = _coerce_dataframe_strings_for_h5ad(adata.obs)
    adata.var = _coerce_dataframe_strings_for_h5ad(adata.var)
    return adata


def _write_sparse_matrix_market(matrix: Any, path: Path) -> None:
    from scipy import sparse

    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    else:
        matrix = matrix.tocoo()

    try:
        from scipy.io import mmwrite as scipy_mmwrite

        scipy_mmwrite(path, matrix)
        return
    except ImportError:
        pass

    value_kind = "integer" if np.issubdtype(matrix.dtype, np.integer) else "real"
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"%%MatrixMarket matrix coordinate {value_kind} general\n")
        handle.write("%\n")
        handle.write(f"{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n")
        for row, col, value in zip(matrix.row, matrix.col, matrix.data, strict=False):
            if value_kind == "integer":
                value_text = str(int(value))
            else:
                value_text = format(float(value), ".15g")
            handle.write(f"{row + 1} {col + 1} {value_text}\n")


def export_method_ready_subset(
    h5ad_path: Path | None,
    subset_metadata: pd.DataFrame,
    out_dir: Path,
    *,
    counts_layer: str = "counts",
    write_h5ad_subset: bool = True,
    write_mtx: bool = True,
    max_cells_per_side_celltype: int | None = None,
    random_state: int = 0,
) -> dict[str, str]:
    import anndata as ad
    from scipy import sparse

    try:
        ad.settings.allow_write_nullable_strings = True
    except Exception:
        pass

    out_dir.mkdir(parents=True, exist_ok=True)
    subset_metadata = subset_metadata.copy()
    if "cell_id" not in subset_metadata.columns:
        raise ValueError("subset_metadata must contain 'cell_id'")
    subset_metadata = downsample_for_communication_export(
        subset_metadata,
        max_cells_per_side_celltype=max_cells_per_side_celltype,
        random_state=random_state,
    )

    def _write_outputs(subset: ad.AnnData) -> dict[str, str]:
        h5ad_out = out_dir / "adata_subset.h5ad"
        if write_h5ad_subset:
            subset_h5ad = _prepare_anndata_for_h5ad_write(subset)
            subset_h5ad.write_h5ad(h5ad_out, compression="gzip")

        counts_out = out_dir / "counts.mtx"
        genes_out = out_dir / "genes.tsv.gz"
        barcodes_out = out_dir / "barcodes.tsv.gz"
        meta_out = out_dir / "cell_metadata.tsv.gz"
        if write_mtx:
            matrix = subset.layers[counts_layer] if counts_layer in subset.layers else subset.X
            if not sparse.issparse(matrix):
                matrix = sparse.csr_matrix(matrix)
            _write_sparse_matrix_market(matrix.T.tocsr(), counts_out)

            gene_symbol_col = next(
                (col for col in ["gene_symbol", "symbol", "feature_name", "gene_name"] if col in subset.var.columns),
                None,
            )
            genes_df = pd.DataFrame(
                {
                    "gene_id": subset.var_names.astype(str),
                    "gene_symbol": subset.var[gene_symbol_col].astype(str).values if gene_symbol_col else subset.var_names.astype(str),
                }
            )
            genes_df.to_csv(genes_out, sep="\t", index=False, compression="gzip")
            pd.DataFrame({"cell_id": subset.obs_names.astype(str)}).to_csv(
                barcodes_out, sep="\t", index=False, compression="gzip"
            )
            subset.obs.reset_index().rename(columns={"index": "cell_id"}).to_csv(
                meta_out, sep="\t", index=False, compression="gzip"
            )

        manifest = {
            "h5ad_subset": str(h5ad_out) if write_h5ad_subset else "",
            "counts_mtx": str(counts_out) if write_mtx else "",
            "genes_tsv": str(genes_out) if write_mtx else "",
            "barcodes_tsv": str(barcodes_out) if write_mtx else "",
            "cell_metadata_tsv": str(meta_out) if write_mtx else "",
            "n_export_cells": str(subset.n_obs),
        }
        write_json(out_dir / "export_manifest.json", manifest)
        return manifest

    if "source_h5ad" in subset_metadata.columns and subset_metadata["source_h5ad"].astype(str).str.len().gt(0).any():
        adatas: list[ad.AnnData] = []
        obs_order = subset_metadata["cell_id"].astype(str).tolist()
        for _, source_df in subset_metadata.groupby("source_h5ad", sort=False):
            expr_path = Path(source_df.get("source_expr_h5ad", source_df["source_h5ad"]).iloc[0])
            meta_path = Path(source_df["source_h5ad"].iloc[0])
            selected_ids = source_df["cell_id"].astype(str).tolist()
            subset = None
            for candidate in [expr_path, meta_path]:
                adata_backed = ad.read_h5ad(candidate, backed="r")
                candidate_ids = [cell_id for cell_id in selected_ids if cell_id in adata_backed.obs_names]
                if not candidate_ids:
                    continue
                subset = adata_backed[candidate_ids].to_memory()
                break
            if subset is None:
                raise ValueError(f"No overlapping cells found for lineage subset in {expr_path} or {meta_path}")

            lookup = source_df.set_index("cell_id")
            subset.var_names_make_unique()
            for column in lookup.columns:
                subset.obs[column] = lookup.loc[subset.obs_names, column].values
            adatas.append(subset)

        if not adatas:
            raise ValueError("No lineage-backed AnnData subsets could be assembled")
        combined = ad.concat(adatas, join="outer", merge="first", index_unique=None)
        combined = combined[obs_order].copy()
        return _write_outputs(combined)

    lookup = subset_metadata.set_index("cell_id")
    if h5ad_path is None:
        raise ValueError("h5ad_path is required when subset_metadata lacks source_h5ad")
    adata_backed = ad.read_h5ad(h5ad_path, backed="r")
    selected_ids = [cell_id for cell_id in subset_metadata["cell_id"].astype(str).tolist() if cell_id in adata_backed.obs_names]
    if not selected_ids:
        raise ValueError(f"No overlapping cell IDs between subset metadata and {h5ad_path}")
    subset = adata_backed[selected_ids].to_memory()
    for column in lookup.columns:
        subset.obs[column] = lookup.loc[subset.obs_names, column].astype(str).values
    return _write_outputs(subset)


def pick_score_column(df: pd.DataFrame, method: str) -> str | None:
    candidates = METHOD_SCORE_CANDIDATES.get(method.lower(), [])
    for column in candidates:
        if column in df.columns:
            return column
    for column in ["score", "rank", "prob", "mean"]:
        if column in df.columns:
            return column
    return None


def standardize_interaction_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    mapping_candidates = {
        "sender": ["sender", "source", "cell_type_a", "source_celltype"],
        "receiver": ["receiver", "target", "cell_type_b", "target_celltype"],
        "ligand": ["ligand", "ligand_complex", "ligand.complex", "gene_a", "partner_a"],
        "receptor": ["receptor", "receptor_complex", "receptor.complex", "gene_b", "partner_b"],
    }
    for target_col, candidates in mapping_candidates.items():
        if target_col in out.columns:
            continue
        for candidate in candidates:
            if candidate in out.columns:
                out[target_col] = out[candidate].astype(str)
                break
        else:
            out[target_col] = ""
    out["interaction_key"] = (
        out["sender"].astype(str)
        + "|"
        + out["receiver"].astype(str)
        + "|"
        + out["ligand"].astype(str)
        + "|"
        + out["receptor"].astype(str)
    )
    return out


def merge_method_side_results(
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    *,
    method: str,
    left_label: str,
    right_label: str,
    contrast_id: str,
    level: str,
) -> pd.DataFrame:
    left_std = standardize_interaction_columns(left_df)
    right_std = standardize_interaction_columns(right_df)
    score_col_left = pick_score_column(left_std, method)
    score_col_right = pick_score_column(right_std, method)

    left_keep = ["interaction_key", "sender", "receiver", "ligand", "receptor"]
    right_keep = ["interaction_key", "sender", "receiver", "ligand", "receptor"]
    if score_col_left:
        left_keep.append(score_col_left)
    if score_col_right:
        right_keep.append(score_col_right)
    left_keep = [col for col in left_keep if col in left_std.columns]
    right_keep = [col for col in right_keep if col in right_std.columns]

    merged = left_std[left_keep].merge(
        right_std[right_keep],
        on=["interaction_key", "sender", "receiver", "ligand", "receptor"],
        how="outer",
        suffixes=("_left", "_right"),
    )
    left_metric = score_col_left if score_col_left else "score"
    right_metric = score_col_right if score_col_right else left_metric
    if score_col_left:
        merged["left_score"] = pd.to_numeric(merged.get(f"{left_metric}_left"), errors="coerce")
    else:
        merged["left_score"] = np.nan
    if score_col_right:
        merged["right_score"] = pd.to_numeric(merged.get(f"{right_metric}_right"), errors="coerce")
    else:
        merged["right_score"] = np.nan
    merged["delta_score_right_minus_left"] = merged["right_score"].fillna(0) - merged["left_score"].fillna(0)
    merged["method"] = method
    merged["contrast_id"] = contrast_id
    merged["level"] = level
    merged["left_label"] = left_label
    merged["right_label"] = right_label
    merged["score_metric_left"] = score_col_left or ""
    merged["score_metric_right"] = score_col_right or ""
    return merged.sort_values("delta_score_right_minus_left", ascending=False).reset_index(drop=True)


def flatten_cellphonedb_table(df: pd.DataFrame, *, value_name: str) -> pd.DataFrame:
    static_cols = [
        col
        for col in [
            "id_cp_interaction",
            "interacting_pair",
            "partner_a",
            "partner_b",
            "gene_a",
            "gene_b",
            "secreted",
            "receptor_a",
            "receptor_b",
            "annotation_strategy",
            "classification",
        ]
        if col in df.columns
    ]
    dynamic_cols = [col for col in df.columns if col not in static_cols]
    melted = df.melt(id_vars=static_cols, value_vars=dynamic_cols, var_name="cellpair", value_name=value_name)
    pair_split = melted["cellpair"].astype(str).str.split(r"\|", n=1, expand=True)
    if pair_split.shape[1] == 2:
        melted["sender"] = pair_split[0].astype(str)
        melted["receiver"] = pair_split[1].astype(str)
    else:
        melted["sender"] = ""
        melted["receiver"] = ""
    melted["ligand"] = melted["partner_a"].astype(str) if "partner_a" in melted.columns else ""
    melted["receptor"] = melted["partner_b"].astype(str) if "partner_b" in melted.columns else ""
    return melted


def write_status(path: Path, payload: Mapping[str, Any]) -> None:
    write_json(path, dict(payload))


def summarize_prepared_manifest(prepared_rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(prepared_rows)
    if df.empty:
        return {"n_prepared": 0, "n_ready": 0, "n_blocked": 0}
    return {
        "n_prepared": int(df.shape[0]),
        "n_ready": int(df["status"].eq("ready").sum()),
        "n_blocked": int(df["status"].eq("blocked").sum()),
        "n_contrasts": int(df["contrast_id"].nunique()),
    }


def load_metadata_and_manifest(
    metadata_tsv: Path = DEFAULT_METADATA_TSV,
    config_manifest: Path = DEFAULT_CONFIG_MANIFEST,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    try:
        metadata = load_lineage_metadata(config_manifest)
    except Exception:
        metadata = pd.read_csv(metadata_tsv, sep="\t")
    sample_manifest = build_sample_manifest(metadata)
    return metadata, sample_manifest
