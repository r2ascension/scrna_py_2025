#!/usr/bin/env python3
"""Upper/lower airway reproduction helper for the respiratory atlas workspace.

This script intentionally reuses existing epithelial objects and DESeq2 outputs
instead of rerunning the full heavy tissue-comparison pipeline by default.
It creates a compact, auditable bridge from the existing repository outputs to
the upper-vs-lower-airway reproduction plan.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd

DEFAULT_H5AD = Path(
    "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508/"
    "epithelial_tissue_comparison_final_fullgene.h5ad"
)
DEFAULT_EXISTING_R_DIR = Path(
    "/home/h2048/data/R/0508/epithelial_tissue_comparison_v1_3_3_rm_leiden14_17_20260508"
)
DEFAULT_CNMF_DIR = Path("/home/h2048/output/program_full_parallel_methods_20260507/epithelial")
DEFAULT_OUT_DIR = Path("/home/h2048/data/py/0523/upper_lower_airway_reproduction_20260523")

OBS_COLUMNS = [
    "tissue",
    "tissue_sampling_method",
    "tissue_level_2",
    "sample",
    "dataset",
    "condition",
    "cell_type_L1",
    "cell_type_L2",
    "cell_type_L3",
    "cell_type_scanvi_pred",
    "cell_type",
    "ann_finest_level",
]

KEY_GENES = [
    "ACE2",
    "TMPRSS2",
    "TMPRSS4",
    "FURIN",
    "CTSL",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "OAS1",
    "OAS2",
    "OAS3",
    "IFI6",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "ISG15",
    "MX1",
    "MUC5AC",
    "MUC5B",
    "MUC1",
    "FOXJ1",
    "KRT5",
    "KRT15",
    "SCGB1A1",
    "SFTPC",
    "AGER",
    "AQP5",
]

PAIRWISE_SITE_ORDER = {
    "nose": 10,
    "sinus": 20,
    "nasal": 30,
    "trachea": 40,
    "bronchial": 50,
    "lower_conducting_airway": 60,
    "distal_lung": 70,
    "lung": 80,
    "unknown": 999,
}


@dataclass(frozen=True)
class SiteAnnotation:
    site_detail: str
    site_group: str
    site_axis: str
    comparison_group: str
    include_in_main_airway: bool


def _clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def normalize_site(tissue: object, tissue_level_2: object = "") -> SiteAnnotation:
    """Map repository tissue labels to upper/lower airway reproduction labels."""
    tissue_s = _clean_text(tissue).lower()
    level_s = _clean_text(tissue_level_2).lower()

    if tissue_s == "nose":
        return SiteAnnotation("nose", "nasal", "upper_airway", "nasal", True)
    if tissue_s == "sinus":
        return SiteAnnotation("sinus", "sinus", "upper_airway", "nasal", True)
    if tissue_s == "respiratory airway":
        if "trachea" in level_s:
            site_group = "trachea"
        elif "bronch" in level_s or "lobular" in level_s:
            site_group = "bronchial"
        else:
            site_group = "lower_conducting_airway"
        return SiteAnnotation(site_group, site_group, "lower_conducting_airway", "bronchial", True)
    if tissue_s == "lung parenchyma":
        return SiteAnnotation("distal_lung", "distal_lung", "distal_lung", "lung", False)
    return SiteAnnotation("unknown", "unknown", "unknown", "unknown", False)


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
            values = []
            for code in codes:
                if code < 0 or code >= len(categories):
                    values.append(pd.NA)
                else:
                    values.append(str(categories[code]))
            return values
        if "values" in obj:
            return [str(value) for value in _decode_array(obj["values"])]
        raise ValueError(f"Unsupported h5ad obs group encoding for column '{name}': {list(obj.keys())}")
    return [str(value) for value in _decode_array(obj)]


def read_obs_columns_h5py(h5ad_path: Path | str, columns: Sequence[str]) -> pd.DataFrame:
    """Read selected obs columns without loading X/layers from a h5ad file."""
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


def read_var_gene_map_h5py(h5ad_path: Path | str) -> dict[str, str]:
    """Read var index -> gene_symbol mapping without loading expression."""
    h5ad_path = Path(h5ad_path)
    with h5py.File(h5ad_path, "r") as handle:
        if "var" not in handle or "_index" not in handle["var"]:
            return {}
        var_group = handle["var"]
        genes = [str(value) for value in _decode_array(var_group["_index"])]
        if "gene_symbol" in var_group:
            symbols = _read_h5ad_column(var_group, "gene_symbol")
        elif "symbol_base" in var_group:
            symbols = _read_h5ad_column(var_group, "symbol_base")
        else:
            symbols = genes
    mapping: dict[str, str] = {}
    for gene, symbol in zip(genes, symbols):
        symbol_s = _clean_text(symbol) or gene
        mapping[gene] = symbol_s
        mapping.setdefault(symbol_s, symbol_s)
    return mapping


def inspect_h5ad_light(h5ad_path: Path | str) -> dict[str, object]:
    """Return lightweight shape/layer/raw metadata for a h5ad file."""
    h5ad_path = Path(h5ad_path)
    with h5py.File(h5ad_path, "r") as handle:
        n_obs = len(handle["obs"]["_index"])
        n_vars = len(handle["var"]["_index"])
        layers = list(handle["layers"].keys()) if "layers" in handle else []
        raw_var_n = None
        if "raw" in handle and "var" in handle["raw"] and "_index" in handle["raw"]["var"]:
            raw_var_n = len(handle["raw"]["var"]["_index"])
        obsm = list(handle["obsm"].keys()) if "obsm" in handle else []
    return {
        "path": str(h5ad_path),
        "n_obs": int(n_obs),
        "n_vars": int(n_vars),
        "layers": layers,
        "has_counts_layer": "counts" in layers,
        "raw_var_n": raw_var_n,
        "obsm": obsm,
    }


def add_site_annotations(metadata: pd.DataFrame) -> pd.DataFrame:
    """Append site_group/site_axis/comparison_group columns while preserving originals."""
    if "tissue" not in metadata.columns:
        raise ValueError("metadata must contain a 'tissue' column")
    tissue_level = metadata["tissue_level_2"] if "tissue_level_2" in metadata.columns else [""] * len(metadata)
    annotations = [normalize_site(t, lvl) for t, lvl in zip(metadata["tissue"], tissue_level)]
    out = metadata.copy()
    out["site_detail"] = [a.site_detail for a in annotations]
    out["site_group"] = [a.site_group for a in annotations]
    out["site_axis"] = [a.site_axis for a in annotations]
    out["comparison_group"] = [a.comparison_group for a in annotations]
    out["include_in_main_airway"] = [a.include_in_main_airway for a in annotations]
    return out


def summarize_metadata(metadata: pd.DataFrame, celltype_col: str = "cell_type_L3") -> dict[str, pd.DataFrame]:
    """Build composition/QC tables from annotated metadata."""
    required = {"tissue", "site_group", "site_axis", "comparison_group", "sample"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(f"metadata missing required columns: {missing}")
    if celltype_col not in metadata.columns:
        raise ValueError(f"metadata missing cell type column: {celltype_col}")

    meta = metadata.copy()
    for col in ["tissue", "site_group", "site_axis", "comparison_group", "sample", "dataset", celltype_col]:
        if col in meta.columns:
            meta[col] = meta[col].astype("string").fillna("Unknown")

    cell_counts_by_site = (
        meta.groupby(["site_axis", "site_group", "comparison_group", "tissue"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["site_axis", "site_group", "tissue"])
    )

    sample_group_cols = ["sample", "dataset", "site_axis", "site_group", "comparison_group", "tissue"]
    sample_group_cols = [col for col in sample_group_cols if col in meta.columns]
    sample_site_summary = (
        meta.groupby(sample_group_cols, dropna=False)
        .agg(n_cells=(celltype_col, "size"), n_cell_types=(celltype_col, "nunique"))
        .reset_index()
        .sort_values(["comparison_group", "site_group", "sample"])
    )

    celltype_site_summary = (
        meta.groupby([celltype_col, "site_axis", "site_group", "comparison_group", "tissue"], dropna=False)
        .agg(n_cells=(celltype_col, "size"), n_samples=("sample", "nunique"))
        .reset_index()
        .rename(columns={celltype_col: "cell_type"})
        .sort_values(["cell_type", "comparison_group", "site_group"])
    )

    celltype_comparison_summary = (
        meta.groupby([celltype_col, "site_axis", "comparison_group"], dropna=False)
        .agg(n_cells=(celltype_col, "size"), n_samples=("sample", "nunique"))
        .reset_index()
        .rename(columns={celltype_col: "cell_type"})
        .sort_values(["cell_type", "comparison_group"])
    )

    celltype_sample_site_counts = (
        meta.groupby([celltype_col, "sample", "site_axis", "site_group", "comparison_group", "tissue"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .rename(columns={celltype_col: "cell_type"})
        .sort_values(["cell_type", "comparison_group", "sample"])
    )

    return {
        "cell_counts_by_site": cell_counts_by_site,
        "sample_site_summary": sample_site_summary,
        "celltype_comparison_summary": celltype_comparison_summary,
        "celltype_site_summary": celltype_site_summary,
        "celltype_sample_site_counts": celltype_sample_site_counts,
    }


def _split_contrast(contrast: str) -> tuple[str, str]:
    if "_vs_" not in contrast:
        return contrast, ""
    left, right = contrast.split("_vs_", 1)
    return left, right


def _top_genes(df: pd.DataFrame, direction: str, n: int = 10) -> str:
    sub = df[df["direction_resolved"] == direction].copy()
    if sub.empty:
        return ""
    sub["padj_sort"] = pd.to_numeric(sub.get("padj"), errors="coerce").fillna(1.0)
    sub["abs_lfc_sort"] = pd.to_numeric(sub.get("log2FoldChange"), errors="coerce").abs().fillna(0.0)
    sub = sub.sort_values(["padj_sort", "abs_lfc_sort"], ascending=[True, False]).head(n)
    return ";".join(sub["gene_symbol"].astype(str).tolist())


def _contrast_token_to_tissue_label(token: object) -> str:
    return _clean_text(token).replace("_", " ")


def _contrast_token_to_site_detail(token: object) -> str:
    token_s = _clean_text(token).lower().replace(" ", "_")
    mapping = {
        "nose": "nose",
        "sinus": "sinus",
        "respiratory_airway": "bronchial",
        "lung_parenchyma": "distal_lung",
    }
    return mapping.get(token_s, token_s)


def _ordered_pair_label(left: object, right: object) -> str:
    left_s = _clean_text(left)
    right_s = _clean_text(right)
    if not left_s and not right_s:
        return ""
    if not left_s:
        return right_s
    if not right_s:
        return left_s
    left_key = (PAIRWISE_SITE_ORDER.get(left_s, 9999), left_s)
    right_key = (PAIRWISE_SITE_ORDER.get(right_s, 9999), right_s)
    if left_key <= right_key:
        return f"{left_s}_vs_{right_s}"
    return f"{right_s}_vs_{left_s}"


def _pair_sort_columns(pair_label: str) -> tuple[int, int, str]:
    left, right = _split_contrast(pair_label)
    return (
        PAIRWISE_SITE_ORDER.get(left, 9999),
        PAIRWISE_SITE_ORDER.get(right, 9999),
        pair_label,
    )


def summarize_deseq2_tree(
    de_root: Path | str,
    level: str,
    gene_map: Mapping[str, str] | None = None,
    key_genes: Iterable[str] = KEY_GENES,
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize existing DESeq2 result CSV files under pseudobulk_de(_L3)."""
    de_root = Path(de_root)
    gene_map = gene_map or {}
    key_gene_set = {gene.upper() for gene in key_genes}
    summary_records: list[dict[str, object]] = []
    key_records: list[pd.DataFrame] = []

    for result_path in sorted(de_root.glob("*/*/DESeq2_results.csv")):
        cell_type = result_path.parent.parent.name
        contrast = result_path.parent.name
        contrast_left, contrast_right = _split_contrast(contrast)
        try:
            df = pd.read_csv(result_path)
        except Exception as exc:
            summary_records.append(
                {
                    "level": level,
                    "cell_type": cell_type,
                    "contrast": contrast,
                    "status": "read_error",
                    "error": str(exc),
                }
            )
            continue
        if "gene" not in df.columns:
            continue
        df["gene"] = df["gene"].astype(str)
        df["gene_symbol"] = df["gene"].map(lambda gene: gene_map.get(gene, gene))
        padj = pd.to_numeric(df.get("padj"), errors="coerce")
        lfc = pd.to_numeric(df.get("log2FoldChange"), errors="coerce")
        sig_col = df.get("sig")
        if sig_col is not None:
            is_sig = sig_col.astype(str).str.lower().eq("sig") | (padj <= alpha)
        else:
            is_sig = padj <= alpha
        direction_from_lfc = np.where(lfc > 0, "up", np.where(lfc < 0, "down", "flat"))
        if "direction" in df.columns:
            direction = df["direction"].astype(str).str.lower()
            direction = direction.where(direction.isin(["up", "down"]), direction_from_lfc)
        else:
            direction = pd.Series(direction_from_lfc, index=df.index)
        df["is_sig_resolved"] = is_sig.fillna(False)
        df["direction_resolved"] = np.where(df["is_sig_resolved"], direction, "ns")
        sig_df = df[df["is_sig_resolved"]]
        summary_records.append(
            {
                "level": level,
                "cell_type": cell_type,
                "contrast": contrast,
                "contrast_left": contrast_left,
                "contrast_right": contrast_right,
                "status": "ok",
                "result_path": str(result_path),
                "n_genes": int(df.shape[0]),
                "n_sig": int(sig_df.shape[0]),
                "n_up": int((sig_df["direction_resolved"] == "up").sum()),
                "n_down": int((sig_df["direction_resolved"] == "down").sum()),
                "median_abs_log2fc_sig": float(lfc[df["is_sig_resolved"]].abs().median()) if sig_df.shape[0] else 0.0,
                "top_up_genes": _top_genes(df, "up"),
                "top_down_genes": _top_genes(df, "down"),
            }
        )
        key_mask = df["gene_symbol"].astype(str).str.upper().isin(key_gene_set) | df["gene"].astype(str).str.upper().isin(key_gene_set)
        if key_mask.any():
            keep_cols = [
                col
                for col in ["gene", "gene_symbol", "baseMean", "log2FoldChange", "pvalue", "padj", "sig", "direction"]
                if col in df.columns
            ]
            key_df = df.loc[key_mask, keep_cols].copy()
            key_df.insert(0, "level", level)
            key_df.insert(1, "cell_type", cell_type)
            key_df.insert(2, "contrast", contrast)
            key_df.insert(3, "contrast_left", contrast_left)
            key_df.insert(4, "contrast_right", contrast_right)
            key_df["is_key_gene"] = True
            key_df["result_path"] = str(result_path)
            key_records.append(key_df)

    summary = pd.DataFrame(summary_records)
    key_gene_df = pd.concat(key_records, ignore_index=True) if key_records else pd.DataFrame()
    return summary, key_gene_df


def annotate_pairwise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add harmonized tissue/site pair columns to DE summary or key-gene tables."""
    extra_object_cols = [
        "contrast_left_label",
        "contrast_right_label",
        "left_site_detail",
        "right_site_detail",
        "left_site_group",
        "right_site_group",
        "left_site_axis",
        "right_site_axis",
        "left_comparison_group",
        "right_comparison_group",
        "contrast_site_pair",
        "canonical_site_pair",
        "contrast_comparison_pair",
        "canonical_comparison_pair",
    ]
    extra_bool_cols = [
        "left_include_in_main_airway",
        "right_include_in_main_airway",
    ]
    if df.empty:
        out = df.copy()
        for column in extra_object_cols:
            out[column] = pd.Series(dtype="object")
        for column in extra_bool_cols:
            out[column] = pd.Series(dtype="bool")
        return out

    out = df.copy()
    left_tokens = out.get("contrast_left", pd.Series(index=out.index, dtype="object")).fillna("")
    right_tokens = out.get("contrast_right", pd.Series(index=out.index, dtype="object")).fillna("")
    left_labels = left_tokens.map(_contrast_token_to_tissue_label)
    right_labels = right_tokens.map(_contrast_token_to_tissue_label)
    left_site_details = left_tokens.map(_contrast_token_to_site_detail)
    right_site_details = right_tokens.map(_contrast_token_to_site_detail)
    left_annotations = [normalize_site(label) for label in left_labels]
    right_annotations = [normalize_site(label) for label in right_labels]

    out["contrast_left_label"] = left_labels
    out["contrast_right_label"] = right_labels
    out["left_site_detail"] = left_site_details
    out["right_site_detail"] = right_site_details
    out["left_site_group"] = [ann.site_group for ann in left_annotations]
    out["right_site_group"] = [ann.site_group for ann in right_annotations]
    out["left_site_axis"] = [ann.site_axis for ann in left_annotations]
    out["right_site_axis"] = [ann.site_axis for ann in right_annotations]
    out["left_comparison_group"] = [ann.comparison_group for ann in left_annotations]
    out["right_comparison_group"] = [ann.comparison_group for ann in right_annotations]
    out["left_include_in_main_airway"] = [ann.include_in_main_airway for ann in left_annotations]
    out["right_include_in_main_airway"] = [ann.include_in_main_airway for ann in right_annotations]
    out["contrast_site_pair"] = [
        f"{left}_vs_{right}" if left and right else ""
        for left, right in zip(left_site_details, right_site_details)
    ]
    out["canonical_site_pair"] = [
        _ordered_pair_label(left, right)
        for left, right in zip(left_site_details, right_site_details)
    ]
    out["contrast_comparison_pair"] = [
        f"{left.comparison_group}_vs_{right.comparison_group}" if left.comparison_group and right.comparison_group else ""
        for left, right in zip(left_annotations, right_annotations)
    ]
    out["canonical_comparison_pair"] = [
        _ordered_pair_label(left.comparison_group, right.comparison_group)
        for left, right in zip(left_annotations, right_annotations)
    ]
    return out


def summarize_pairwise_catalog(catalog: pd.DataFrame, pair_col: str) -> pd.DataFrame:
    """Collapse annotated DE tables into pairwise availability/count summaries."""
    if catalog.empty or pair_col not in catalog.columns:
        return pd.DataFrame(
            columns=[
                "level",
                "pair_label",
                "n_contrast_files",
                "n_cell_types",
                "cell_types",
                "contrasts",
                "n_sig_total",
                "median_n_sig",
            ]
        )

    data = catalog.copy()
    if "status" in data.columns:
        data = data[data["status"].astype(str) == "ok"].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    if data.empty:
        return pd.DataFrame(
            columns=[
                "level",
                "pair_label",
                "n_contrast_files",
                "n_cell_types",
                "cell_types",
                "contrasts",
                "n_sig_total",
                "median_n_sig",
            ]
        )

    records: list[dict[str, object]] = []
    group_cols = ["level", pair_col] if "level" in data.columns else [pair_col]
    for keys, sub in data.groupby(group_cols, dropna=False):
        if isinstance(keys, tuple):
            level, pair_label = keys
        else:
            level, pair_label = "", keys
        n_sig = pd.to_numeric(sub.get("n_sig"), errors="coerce") if "n_sig" in sub.columns else pd.Series(dtype=float)
        records.append(
            {
                "level": level,
                "pair_label": pair_label,
                "n_contrast_files": int(sub.shape[0]),
                "n_cell_types": int(sub["cell_type"].astype(str).nunique()) if "cell_type" in sub.columns else 0,
                "cell_types": ";".join(sorted(sub["cell_type"].astype(str).unique())) if "cell_type" in sub.columns else "",
                "contrasts": ";".join(sorted(sub["contrast"].astype(str).unique())) if "contrast" in sub.columns else "",
                "n_sig_total": int(n_sig.fillna(0).sum()) if not n_sig.empty else 0,
                "median_n_sig": float(n_sig.median()) if not n_sig.dropna().empty else 0.0,
            }
        )

    summary = pd.DataFrame(records)
    if summary.empty:
        return summary
    sort_parts = summary["pair_label"].astype(str).map(_pair_sort_columns)
    summary["_pair_sort_left"] = sort_parts.map(lambda x: x[0])
    summary["_pair_sort_right"] = sort_parts.map(lambda x: x[1])
    summary = summary.sort_values(["level", "_pair_sort_left", "_pair_sort_right", "pair_label"]).drop(
        columns=["_pair_sort_left", "_pair_sort_right"]
    )
    return summary.reset_index(drop=True)


def _sorted_pair_labels(labels: Iterable[object]) -> list[str]:
    return sorted({_clean_text(label) for label in labels if _clean_text(label)}, key=_pair_sort_columns)


def build_deg_metric_matrix(
    de_summary: pd.DataFrame,
    value_col: str,
    *,
    level: str = "L3",
    pair_col: str = "canonical_site_pair",
    index_col: str = "cell_type",
    aggfunc: str = "sum",
) -> pd.DataFrame:
    """Pivot a DEG summary metric into a cell_type × pair matrix."""
    if de_summary.empty or value_col not in de_summary.columns or pair_col not in de_summary.columns:
        return pd.DataFrame()

    data = de_summary.copy()
    if "status" in data.columns:
        data = data[data["status"].astype(str) == "ok"].copy()
    if "level" in data.columns:
        data = data[data["level"].astype(str) == level].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    if data.empty:
        return pd.DataFrame()

    data[value_col] = pd.to_numeric(data[value_col], errors="coerce").fillna(0.0)
    matrix = data.pivot_table(
        index=index_col,
        columns=pair_col,
        values=value_col,
        aggfunc=aggfunc,
        fill_value=0.0,
    )
    if matrix.empty:
        return matrix

    pair_order = [label for label in _sorted_pair_labels(matrix.columns) if label in matrix.columns]
    matrix = matrix.loc[:, pair_order]
    row_order = matrix.sum(axis=1).sort_values(ascending=False).index.tolist()
    return matrix.loc[row_order]


def build_direction_balance_matrix(
    de_summary: pd.DataFrame,
    *,
    level: str = "L3",
    pair_col: str = "canonical_site_pair",
    index_col: str = "cell_type",
) -> pd.DataFrame:
    """Return a matrix of `(n_up - n_down) / n_sig` for each cell_type × pair."""
    required = {"n_up", "n_down", "n_sig", pair_col}
    if de_summary.empty or not required.issubset(de_summary.columns):
        return pd.DataFrame()

    data = de_summary.copy()
    if "status" in data.columns:
        data = data[data["status"].astype(str) == "ok"].copy()
    if "level" in data.columns:
        data = data[data["level"].astype(str) == level].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    if data.empty:
        return pd.DataFrame()

    n_sig = pd.to_numeric(data["n_sig"], errors="coerce").replace(0, np.nan)
    n_up = pd.to_numeric(data["n_up"], errors="coerce").fillna(0.0)
    n_down = pd.to_numeric(data["n_down"], errors="coerce").fillna(0.0)
    data["direction_balance"] = ((n_up - n_down) / n_sig).replace([np.inf, -np.inf], np.nan).fillna(0.0)

    matrix = data.pivot_table(
        index=index_col,
        columns=pair_col,
        values="direction_balance",
        aggfunc="mean",
        fill_value=0.0,
    )
    if matrix.empty:
        return matrix

    pair_order = [label for label in _sorted_pair_labels(matrix.columns) if label in matrix.columns]
    matrix = matrix.loc[:, pair_order]
    magnitude_order = matrix.abs().sum(axis=1).sort_values(ascending=False).index.tolist()
    return matrix.loc[magnitude_order]


def prepare_pairwise_coverage_table(*pairwise_tables: pd.DataFrame) -> pd.DataFrame:
    """Combine L2/L3 pairwise-count tables into a single long-form table."""
    frames: list[pd.DataFrame] = []
    for table in pairwise_tables:
        if table is None or table.empty:
            continue
        keep_cols = [
            col
            for col in ["level", "pair_label", "n_contrast_files", "n_cell_types", "n_sig_total", "median_n_sig"]
            if col in table.columns
        ]
        if keep_cols:
            frames.append(table.loc[:, keep_cols].copy())
    if not frames:
        return pd.DataFrame(columns=["level", "pair_label", "n_contrast_files", "n_cell_types", "n_sig_total", "median_n_sig"])

    combined = pd.concat(frames, ignore_index=True)
    pair_order = _sorted_pair_labels(combined["pair_label"])
    combined["_pair_sort"] = combined["pair_label"].map({label: idx for idx, label in enumerate(pair_order)})
    combined = combined.sort_values(["_pair_sort", "level", "pair_label"]).drop(columns=["_pair_sort"])
    return combined.reset_index(drop=True)


def summarize_top_gene_recurrence(
    de_summary: pd.DataFrame,
    *,
    direction: str,
    level: str = "L3",
    pair_col: str = "canonical_site_pair",
) -> pd.DataFrame:
    """Count how often top genes recur across cell types for each pairwise contrast."""
    column = "top_up_genes" if direction == "up" else "top_down_genes"
    if de_summary.empty or column not in de_summary.columns or pair_col not in de_summary.columns:
        return pd.DataFrame(
            columns=[
                "level",
                "pair_label",
                "direction",
                "gene_symbol",
                "n_cell_types",
                "n_occurrences",
                "weighted_rank_score",
                "best_rank",
                "cell_types",
            ]
        )

    data = de_summary.copy()
    if "status" in data.columns:
        data = data[data["status"].astype(str) == "ok"].copy()
    if "level" in data.columns:
        data = data[data["level"].astype(str) == level].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    if data.empty:
        return pd.DataFrame()

    records: list[dict[str, object]] = []
    for _, row in data.iterrows():
        genes = [gene.strip() for gene in _clean_text(row.get(column)).split(";") if gene.strip()]
        for rank, gene in enumerate(genes, start=1):
            records.append(
                {
                    "level": level,
                    "pair_label": _clean_text(row.get(pair_col)),
                    "direction": direction,
                    "cell_type": _clean_text(row.get("cell_type")),
                    "gene_symbol": gene,
                    "rank": rank,
                    "rank_score": 1.0 / rank,
                }
            )
    if not records:
        return pd.DataFrame()

    exploded = pd.DataFrame(records)
    summary = (
        exploded.groupby(["level", "pair_label", "direction", "gene_symbol"], dropna=False)
        .agg(
            n_cell_types=("cell_type", "nunique"),
            n_occurrences=("gene_symbol", "size"),
            weighted_rank_score=("rank_score", "sum"),
            best_rank=("rank", "min"),
            cell_types=("cell_type", lambda vals: ";".join(sorted({str(v) for v in vals if _clean_text(v)}))),
        )
        .reset_index()
    )
    pair_order = {label: idx for idx, label in enumerate(_sorted_pair_labels(summary["pair_label"]))}
    summary["_pair_sort"] = summary["pair_label"].map(pair_order)
    summary = summary.sort_values(
        ["_pair_sort", "pair_label", "n_cell_types", "weighted_rank_score", "best_rank", "gene_symbol"],
        ascending=[True, True, False, False, True, True],
    ).drop(columns=["_pair_sort"])
    return summary.reset_index(drop=True)


def prepare_key_gene_panel_matrices(
    key_gene_de: pd.DataFrame,
    *,
    level: str = "L3",
    pair_col: str = "canonical_site_pair",
    celltype_order: Sequence[str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Build per-pair key-gene log2FC heatmap matrices."""
    required = {"gene_symbol", "cell_type", "log2FoldChange", pair_col}
    if key_gene_de.empty or not required.issubset(key_gene_de.columns):
        return {}

    data = key_gene_de.copy()
    if "level" in data.columns:
        data = data[data["level"].astype(str) == level].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    if data.empty:
        return {}

    data["gene_symbol"] = data["gene_symbol"].astype(str)
    data["cell_type"] = data["cell_type"].astype(str)
    data["log2FoldChange"] = pd.to_numeric(data["log2FoldChange"], errors="coerce")
    data["padj_sort"] = pd.to_numeric(data.get("padj"), errors="coerce").fillna(1.0)
    data["abs_lfc_sort"] = data["log2FoldChange"].abs().fillna(0.0)
    data = data.sort_values(["padj_sort", "abs_lfc_sort"], ascending=[True, False])
    data = data.drop_duplicates(subset=[pair_col, "gene_symbol", "cell_type"], keep="first")

    ordered_genes = [gene for gene in KEY_GENES if gene in set(data["gene_symbol"])]
    matrices: dict[str, pd.DataFrame] = {}
    for pair_label in _sorted_pair_labels(data[pair_col]):
        sub = data[data[pair_col].astype(str) == pair_label].copy()
        if sub.empty:
            continue
        matrix = sub.pivot_table(index="gene_symbol", columns="cell_type", values="log2FoldChange", aggfunc="first")
        if ordered_genes:
            matrix = matrix.reindex(index=[gene for gene in ordered_genes if gene in matrix.index])
        else:
            matrix = matrix.sort_index()
        if celltype_order:
            ordered_cols = [cell_type for cell_type in celltype_order if cell_type in matrix.columns]
            extra_cols = sorted([cell_type for cell_type in matrix.columns if cell_type not in ordered_cols])
            matrix = matrix.reindex(columns=ordered_cols + extra_cols)
        else:
            matrix = matrix.reindex(columns=sorted(matrix.columns))
        if not matrix.empty:
            matrices[pair_label] = matrix
    return matrices


def _load_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save_figure(fig, stem: Path) -> tuple[str, str]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    return str(png_path), str(pdf_path)


def plot_pairwise_coverage_by_level(coverage_table: pd.DataFrame, stem: Path) -> tuple[str, str] | None:
    """Plot L2/L3 pairwise coverage as side-by-side barplots."""
    if coverage_table.empty:
        return None

    plt = _load_matplotlib()
    plot_df = coverage_table.copy()
    levels = [level for level in ["L2", "L3"] if level in set(plot_df["level"].astype(str))]
    if not levels:
        levels = sorted(plot_df["level"].astype(str).unique())
    pairs = _sorted_pair_labels(plot_df["pair_label"])
    if not pairs:
        return None

    pivot = plot_df.pivot_table(index="pair_label", columns="level", values="n_cell_types", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(index=pairs, columns=levels, fill_value=0)
    x = np.arange(len(pairs), dtype=float)
    width = 0.38 if len(levels) > 1 else 0.56
    colors = {"L2": "#5B8FF9", "L3": "#D14A61"}

    fig, ax = plt.subplots(figsize=(max(8, len(pairs) * 1.25), 5.8))
    for idx, level_name in enumerate(levels):
        offset = (idx - (len(levels) - 1) / 2.0) * width
        values = pivot[level_name].to_numpy(dtype=float)
        bars = ax.bar(x + offset, values, width=width, color=colors.get(level_name, "#7F7F7F"), label=level_name)
        for bar, value in zip(bars, values):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                value + 0.15,
                f"{int(value)}",
                ha="center",
                va="bottom",
                fontsize=9,
            )

    ax.set_xticks(x)
    ax.set_xticklabels([label.replace("_vs_", "\nvs\n") for label in pairs], fontsize=9)
    ax.set_ylabel("Cell types with reused DE contrasts")
    ax.set_title("DEG pairwise coverage across epithelial L2/L3 summaries")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    ax.set_axisbelow(True)
    fig.tight_layout()
    out = _save_figure(fig, stem)
    plt.close(fig)
    return out


def plot_deg_metric_heatmap(
    matrix: pd.DataFrame,
    *,
    stem: Path,
    title: str,
    cmap: str = "magma",
    annotate: bool = True,
) -> tuple[str, str] | None:
    """Plot a small DEG matrix heatmap."""
    if matrix.empty:
        return None

    plt = _load_matplotlib()
    values = matrix.to_numpy(dtype=float)
    vmax = float(np.nanmax(values)) if np.isfinite(values).any() else 1.0
    fig_w = max(8.0, matrix.shape[1] * 1.2 + 1.8)
    fig_h = max(5.2, matrix.shape[0] * 0.45 + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=0.0, vmax=max(vmax, 1.0))
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels([label.replace("_vs_", "\nvs\n") for label in matrix.columns], fontsize=9)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    ax.set_title(title)
    ax.set_xlabel("Pairwise site comparison")
    ax.set_ylabel("Cell type")

    if annotate and matrix.shape[0] * matrix.shape[1] <= 140:
        threshold = vmax * 0.55 if vmax > 0 else 0.0
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = values[i, j]
                text_color = "white" if value >= threshold else "black"
                ax.text(j, i, f"{int(round(value))}", ha="center", va="center", fontsize=8, color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Significant genes (n)")
    fig.tight_layout()
    out = _save_figure(fig, stem)
    plt.close(fig)
    return out


def plot_direction_balance_heatmap(matrix: pd.DataFrame, *, stem: Path, title: str) -> tuple[str, str] | None:
    """Plot a diverging heatmap of DEG direction balance."""
    if matrix.empty:
        return None

    plt = _load_matplotlib()
    values = matrix.to_numpy(dtype=float)
    vmax = float(np.nanmax(np.abs(values))) if np.isfinite(values).any() else 1.0
    vmax = max(vmax, 1.0)
    fig_w = max(8.0, matrix.shape[1] * 1.2 + 1.8)
    fig_h = max(5.2, matrix.shape[0] * 0.45 + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels([label.replace("_vs_", "\nvs\n") for label in matrix.columns], fontsize=9)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)
    ax.set_title(title)
    ax.set_xlabel("Pairwise site comparison")
    ax.set_ylabel("Cell type")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("(n_up - n_down) / n_sig")
    fig.tight_layout()
    out = _save_figure(fig, stem)
    plt.close(fig)
    return out


def plot_key_gene_panel_heatmaps(
    matrices: Mapping[str, pd.DataFrame],
    *,
    stem: Path,
    title: str,
) -> tuple[str, str] | None:
    """Plot per-pair key-gene log2FC heatmaps in a faceted layout."""
    if not matrices:
        return None

    plt = _load_matplotlib()
    pair_labels = [label for label in _sorted_pair_labels(matrices) if label in matrices]
    if not pair_labels:
        return None

    finite_values = []
    for matrix in matrices.values():
        values = matrix.to_numpy(dtype=float)
        finite_values.extend(np.abs(values[np.isfinite(values)]).tolist())
    vmax = max(finite_values) if finite_values else 1.0
    vmax = max(vmax, 1.0)

    n_panels = len(pair_labels)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols
    max_rows = max(matrix.shape[0] for matrix in matrices.values())
    max_cols = max(matrix.shape[1] for matrix in matrices.values())
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(max(9.5, n_cols * max(max_cols * 0.55 + 2.4, 4.6)), max(5.8, n_rows * max(max_rows * 0.32 + 1.9, 3.6))),
        squeeze=False,
    )

    image = None
    for ax, pair_label in zip(axes.flat, pair_labels):
        matrix = matrices[pair_label]
        values = matrix.to_numpy(dtype=float)
        image = ax.imshow(values, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
        ax.set_title(pair_label.replace("_vs_", " vs "), fontsize=11)
        ax.set_xticks(np.arange(matrix.shape[1]))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(np.arange(matrix.shape[0]))
        ax.set_yticklabels(matrix.index, fontsize=8)
        ax.set_xlabel("Cell type", fontsize=9)
        ax.set_ylabel("Key gene", fontsize=9)

    for ax in axes.flat[n_panels:]:
        ax.axis("off")

    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), fraction=0.018, pad=0.015, label="log2 fold-change")
    fig.suptitle(title, fontsize=13)
    fig.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.9, wspace=0.55, hspace=0.75)
    out = _save_figure(fig, stem)
    plt.close(fig)
    return out


def plot_top_gene_recurrence(
    recurrence: pd.DataFrame,
    *,
    stem: Path,
    title: str,
    top_n: int = 8,
    color: str = "#5B8FF9",
) -> tuple[str, str] | None:
    """Plot recurrent top genes across cell types for each pairwise comparison."""
    if recurrence.empty:
        return None

    plt = _load_matplotlib()
    pair_labels = _sorted_pair_labels(recurrence["pair_label"])
    if not pair_labels:
        return None
    n_panels = len(pair_labels)
    n_cols = min(3, n_panels)
    n_rows = (n_panels + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(max(9.0, n_cols * 5.1), max(5.2, n_rows * 3.8)), squeeze=False)

    for ax, pair_label in zip(axes.flat, pair_labels):
        sub = recurrence[recurrence["pair_label"].astype(str) == pair_label].copy().head(top_n)
        sub = sub.sort_values(["n_cell_types", "weighted_rank_score", "gene_symbol"], ascending=[True, True, False])
        if sub.empty:
            ax.axis("off")
            continue
        bars = ax.barh(sub["gene_symbol"], sub["n_cell_types"], color=color, alpha=0.9)
        for bar, (_, row) in zip(bars, sub.iterrows()):
            ax.text(
                bar.get_width() + 0.08,
                bar.get_y() + bar.get_height() / 2.0,
                f"{int(row['n_cell_types'])} ct | score={row['weighted_rank_score']:.2f}",
                va="center",
                ha="left",
                fontsize=8,
            )
        ax.set_title(pair_label.replace("_vs_", " vs "), fontsize=11)
        ax.set_xlabel("Cell types containing gene in top list")
        ax.grid(axis="x", alpha=0.2)
        ax.set_axisbelow(True)

    for ax in axes.flat[n_panels:]:
        ax.axis("off")

    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = _save_figure(fig, stem)
    plt.close(fig)
    return out


def generate_deg_summary_visuals(
    out_dir: Path,
    *,
    de_summary_l2: pd.DataFrame,
    de_summary_l3: pd.DataFrame,
    de_pairwise_site_l2: pd.DataFrame,
    de_pairwise_site_l3: pd.DataFrame,
    key_gene_de: pd.DataFrame,
) -> tuple[dict[str, str], pd.DataFrame]:
    """Create lightweight DEG summary figures and companion tables."""
    visual_outputs: dict[str, str] = {}
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    try:
        _load_matplotlib()
    except ImportError as exc:
        manifest = pd.DataFrame(
            [
                {
                    "figure_id": "deg_summary_visuals",
                    "status": "skipped",
                    "description": f"matplotlib unavailable: {exc}",
                    "png_path": "",
                    "pdf_path": "",
                }
            ]
        )
        manifest_path = out_dir / "deg_visual_manifest.tsv"
        manifest.to_csv(manifest_path, sep="\t", index=False)
        visual_outputs["deg_visual_manifest"] = str(manifest_path)
        return visual_outputs, manifest

    manifest_records: list[dict[str, object]] = []

    coverage_table = prepare_pairwise_coverage_table(de_pairwise_site_l2, de_pairwise_site_l3)
    coverage_path = out_dir / "deg_pairwise_site_coverage_by_level.tsv"
    coverage_table.to_csv(coverage_path, sep="\t", index=False)
    visual_outputs["deg_pairwise_site_coverage_by_level"] = str(coverage_path)
    coverage_plot = plot_pairwise_coverage_by_level(coverage_table, figures_dir / "deg_pairwise_site_coverage_by_level")
    if coverage_plot:
        manifest_records.append(
            {
                "figure_id": "deg_pairwise_site_coverage_by_level",
                "status": "ok",
                "description": "L2/L3 pairwise DEG coverage barplot across site pairs.",
                "png_path": coverage_plot[0],
                "pdf_path": coverage_plot[1],
            }
        )
        visual_outputs["deg_pairwise_site_coverage_by_level_png"] = coverage_plot[0]
        visual_outputs["deg_pairwise_site_coverage_by_level_pdf"] = coverage_plot[1]

    l3_nsig_matrix = build_deg_metric_matrix(de_summary_l3, "n_sig", level="L3")
    l3_nsig_matrix_path = out_dir / "deg_l3_nsig_matrix.tsv"
    l3_nsig_matrix.to_csv(l3_nsig_matrix_path, sep="\t")
    visual_outputs["deg_l3_nsig_matrix"] = str(l3_nsig_matrix_path)
    nsig_plot = plot_deg_metric_heatmap(
        l3_nsig_matrix,
        stem=figures_dir / "deg_l3_nsig_heatmap",
        title="L3 DEG burden heatmap (significant genes per pair)",
    )
    if nsig_plot:
        manifest_records.append(
            {
                "figure_id": "deg_l3_nsig_heatmap",
                "status": "ok",
                "description": "Heatmap of significant-gene counts for each L3 cell type across pairwise site contrasts.",
                "png_path": nsig_plot[0],
                "pdf_path": nsig_plot[1],
            }
        )
        visual_outputs["deg_l3_nsig_heatmap_png"] = nsig_plot[0]
        visual_outputs["deg_l3_nsig_heatmap_pdf"] = nsig_plot[1]

    balance_matrix = build_direction_balance_matrix(de_summary_l3, level="L3")
    balance_matrix = balance_matrix.reindex(index=l3_nsig_matrix.index, columns=l3_nsig_matrix.columns, fill_value=0.0)
    balance_matrix_path = out_dir / "deg_l3_direction_balance_matrix.tsv"
    balance_matrix.to_csv(balance_matrix_path, sep="\t")
    visual_outputs["deg_l3_direction_balance_matrix"] = str(balance_matrix_path)
    balance_plot = plot_direction_balance_heatmap(
        balance_matrix,
        stem=figures_dir / "deg_l3_direction_balance_heatmap",
        title="L3 DEG direction balance heatmap (up- vs down-weighted)",
    )
    if balance_plot:
        manifest_records.append(
            {
                "figure_id": "deg_l3_direction_balance_heatmap",
                "status": "ok",
                "description": "Diverging heatmap showing whether each L3 contrast is up- or down-dominant.",
                "png_path": balance_plot[0],
                "pdf_path": balance_plot[1],
            }
        )
        visual_outputs["deg_l3_direction_balance_heatmap_png"] = balance_plot[0]
        visual_outputs["deg_l3_direction_balance_heatmap_pdf"] = balance_plot[1]

    key_gene_matrices = prepare_key_gene_panel_matrices(key_gene_de, level="L3", celltype_order=list(l3_nsig_matrix.index))
    key_gene_plot = plot_key_gene_panel_heatmaps(
        key_gene_matrices,
        stem=figures_dir / "deg_l3_key_gene_log2fc_panels",
        title="L3 key-gene log2FC panels across pairwise site contrasts",
    )
    if key_gene_plot:
        manifest_records.append(
            {
                "figure_id": "deg_l3_key_gene_log2fc_panels",
                "status": "ok",
                "description": "Faceted heatmaps of curated key genes across L3 cell types for each pairwise site contrast.",
                "png_path": key_gene_plot[0],
                "pdf_path": key_gene_plot[1],
            }
        )
        visual_outputs["deg_l3_key_gene_log2fc_panels_png"] = key_gene_plot[0]
        visual_outputs["deg_l3_key_gene_log2fc_panels_pdf"] = key_gene_plot[1]

    recurrence_up = summarize_top_gene_recurrence(de_summary_l3, direction="up", level="L3")
    recurrence_down = summarize_top_gene_recurrence(de_summary_l3, direction="down", level="L3")
    recurrence_up_path = out_dir / "deg_l3_top_up_gene_recurrence.tsv"
    recurrence_down_path = out_dir / "deg_l3_top_down_gene_recurrence.tsv"
    recurrence_up.to_csv(recurrence_up_path, sep="\t", index=False)
    recurrence_down.to_csv(recurrence_down_path, sep="\t", index=False)
    visual_outputs["deg_l3_top_up_gene_recurrence"] = str(recurrence_up_path)
    visual_outputs["deg_l3_top_down_gene_recurrence"] = str(recurrence_down_path)

    up_plot = plot_top_gene_recurrence(
        recurrence_up,
        stem=figures_dir / "deg_l3_top_up_gene_recurrence",
        title="L3 recurrent top-up genes across pairwise site contrasts",
        color="#5B8FF9",
    )
    if up_plot:
        manifest_records.append(
            {
                "figure_id": "deg_l3_top_up_gene_recurrence",
                "status": "ok",
                "description": "Per-pair recurrent top up-regulated genes across L3 cell types.",
                "png_path": up_plot[0],
                "pdf_path": up_plot[1],
            }
        )
        visual_outputs["deg_l3_top_up_gene_recurrence_png"] = up_plot[0]
        visual_outputs["deg_l3_top_up_gene_recurrence_pdf"] = up_plot[1]

    down_plot = plot_top_gene_recurrence(
        recurrence_down,
        stem=figures_dir / "deg_l3_top_down_gene_recurrence",
        title="L3 recurrent top-down genes across pairwise site contrasts",
        color="#D14A61",
    )
    if down_plot:
        manifest_records.append(
            {
                "figure_id": "deg_l3_top_down_gene_recurrence",
                "status": "ok",
                "description": "Per-pair recurrent top down-regulated genes across L3 cell types.",
                "png_path": down_plot[0],
                "pdf_path": down_plot[1],
            }
        )
        visual_outputs["deg_l3_top_down_gene_recurrence_png"] = down_plot[0]
        visual_outputs["deg_l3_top_down_gene_recurrence_pdf"] = down_plot[1]

    manifest = pd.DataFrame(manifest_records)
    manifest_path = out_dir / "deg_visual_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    visual_outputs["deg_visual_manifest"] = str(manifest_path)
    return visual_outputs, manifest


def summarize_sig_deg_table(
    de_root: Path | str,
    level: str,
    gene_map: Mapping[str, str] | None = None,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Extract significant DEG rows into a long table suitable for downstream PPI seeding."""
    de_root = Path(de_root)
    gene_map = gene_map or {}
    records: list[pd.DataFrame] = []

    for result_path in sorted(de_root.glob("*/*/DESeq2_results.csv")):
        cell_type = result_path.parent.parent.name
        contrast = result_path.parent.name
        contrast_left, contrast_right = _split_contrast(contrast)
        try:
            df = pd.read_csv(result_path)
        except Exception:
            continue
        if "gene" not in df.columns:
            continue

        df = df.copy()
        df["gene"] = df["gene"].astype(str)
        df["gene_symbol"] = df["gene"].map(lambda gene: gene_map.get(gene, gene))
        padj = pd.to_numeric(df.get("padj"), errors="coerce")
        lfc = pd.to_numeric(df.get("log2FoldChange"), errors="coerce")
        sig_col = df.get("sig")
        if sig_col is not None:
            is_sig = sig_col.astype(str).str.lower().eq("sig") | (padj <= alpha)
        else:
            is_sig = padj <= alpha
        direction_from_lfc = np.where(lfc > 0, "up", np.where(lfc < 0, "down", "flat"))
        if "direction" in df.columns:
            direction = df["direction"].astype(str).str.lower()
            direction = direction.where(direction.isin(["up", "down"]), direction_from_lfc)
        else:
            direction = pd.Series(direction_from_lfc, index=df.index)

        sig_df = df.loc[is_sig.fillna(False)].copy()
        if sig_df.empty:
            continue

        sig_df.insert(0, "level", level)
        sig_df.insert(1, "cell_type", cell_type)
        sig_df.insert(2, "contrast", contrast)
        sig_df.insert(3, "contrast_left", contrast_left)
        sig_df.insert(4, "contrast_right", contrast_right)
        sig_df["result_path"] = str(result_path)
        sig_df["direction_resolved"] = direction.loc[sig_df.index].astype(str)
        sig_df["abs_log2FoldChange"] = pd.to_numeric(sig_df.get("log2FoldChange"), errors="coerce").abs().fillna(0.0)
        sig_df["padj_numeric"] = pd.to_numeric(sig_df.get("padj"), errors="coerce").fillna(1.0).clip(lower=1e-300)
        sig_df["neg_log10_padj"] = -np.log10(sig_df["padj_numeric"])
        sig_df["ppi_priority_score"] = sig_df["abs_log2FoldChange"] * sig_df["neg_log10_padj"]

        keep_cols = [
            "level",
            "cell_type",
            "contrast",
            "contrast_left",
            "contrast_right",
            "gene",
            "gene_symbol",
            "baseMean",
            "log2FoldChange",
            "abs_log2FoldChange",
            "pvalue",
            "padj",
            "padj_numeric",
            "neg_log10_padj",
            "direction_resolved",
            "ppi_priority_score",
            "result_path",
        ]
        keep_cols = [col for col in keep_cols if col in sig_df.columns]
        records.append(sig_df.loc[:, keep_cols].copy())

    if not records:
        return pd.DataFrame(
            columns=[
                "level",
                "cell_type",
                "contrast",
                "contrast_left",
                "contrast_right",
                "gene",
                "gene_symbol",
                "baseMean",
                "log2FoldChange",
                "abs_log2FoldChange",
                "pvalue",
                "padj",
                "padj_numeric",
                "neg_log10_padj",
                "direction_resolved",
                "ppi_priority_score",
                "result_path",
            ]
        )
    return pd.concat(records, ignore_index=True)


def summarize_ppi_seed_candidates(
    ppi_deg: pd.DataFrame,
    *,
    level: str = "L3",
    pair_col: str = "canonical_site_pair",
) -> pd.DataFrame:
    """Aggregate significant DEGs into pairwise PPI seed candidates."""
    template_columns = [
        "level",
        "pair_label",
        "gene_symbol",
        "string_input_gene",
        "n_cell_types",
        "n_occurrences",
        "cell_types",
        "contrasts",
        "mean_log2FoldChange",
        "mean_abs_log2FoldChange",
        "max_abs_log2FoldChange",
        "min_padj",
        "neg_log10_min_padj",
        "up_hits",
        "down_hits",
        "dominant_direction",
        "direction_consistency",
        "ppi_priority_score",
        "seed_rank",
    ]
    required = {pair_col, "gene_symbol", "cell_type", "log2FoldChange", "padj_numeric", "abs_log2FoldChange"}
    if ppi_deg.empty or not required.issubset(ppi_deg.columns):
        return pd.DataFrame(columns=template_columns)

    data = ppi_deg.copy()
    if "level" in data.columns:
        data = data[data["level"].astype(str) == level].copy()
    data = data[data[pair_col].astype(str) != ""].copy()
    data = data[data["gene_symbol"].astype(str) != ""].copy()
    if data.empty:
        return pd.DataFrame(columns=template_columns)

    if "direction_resolved" not in data.columns:
        lfc = pd.to_numeric(data["log2FoldChange"], errors="coerce")
        data["direction_resolved"] = np.where(lfc > 0, "up", np.where(lfc < 0, "down", "flat"))

    summary = (
        data.groupby(["level", pair_col, "gene_symbol"], dropna=False)
        .agg(
            n_cell_types=("cell_type", "nunique"),
            n_occurrences=("gene_symbol", "size"),
            cell_types=("cell_type", lambda vals: ";".join(sorted({str(v) for v in vals if _clean_text(v)}))),
            contrasts=("contrast", lambda vals: ";".join(sorted({str(v) for v in vals if _clean_text(v)}))),
            mean_log2FoldChange=("log2FoldChange", lambda vals: float(pd.to_numeric(vals, errors="coerce").mean())),
            mean_abs_log2FoldChange=("abs_log2FoldChange", lambda vals: float(pd.to_numeric(vals, errors="coerce").mean())),
            max_abs_log2FoldChange=("abs_log2FoldChange", lambda vals: float(pd.to_numeric(vals, errors="coerce").max())),
            min_padj=("padj_numeric", lambda vals: float(pd.to_numeric(vals, errors="coerce").min())),
            up_hits=("direction_resolved", lambda vals: int(pd.Series(vals, dtype="string").eq("up").sum())),
            down_hits=("direction_resolved", lambda vals: int(pd.Series(vals, dtype="string").eq("down").sum())),
        )
        .reset_index()
        .rename(columns={pair_col: "pair_label"})
    )
    if summary.empty:
        return pd.DataFrame(columns=template_columns)

    summary["string_input_gene"] = summary["gene_symbol"].astype(str)
    summary["neg_log10_min_padj"] = -np.log10(pd.to_numeric(summary["min_padj"], errors="coerce").fillna(1.0).clip(lower=1e-300))
    dominant: list[str] = []
    consistency: list[float] = []
    for _, row in summary.iterrows():
        up_hits = int(row["up_hits"])
        down_hits = int(row["down_hits"])
        total = max(int(row["n_occurrences"]), 1)
        if up_hits > down_hits:
            dominant.append("up")
            consistency.append(up_hits / total)
        elif down_hits > up_hits:
            dominant.append("down")
            consistency.append(down_hits / total)
        else:
            dominant.append("mixed")
            consistency.append(max(up_hits, down_hits) / total if total else 0.0)
    summary["dominant_direction"] = dominant
    summary["direction_consistency"] = consistency
    summary["ppi_priority_score"] = (
        pd.to_numeric(summary["n_cell_types"], errors="coerce").fillna(0.0)
        * pd.to_numeric(summary["mean_abs_log2FoldChange"], errors="coerce").fillna(0.0)
        * pd.to_numeric(summary["neg_log10_min_padj"], errors="coerce").fillna(0.0)
        * pd.to_numeric(summary["direction_consistency"], errors="coerce").fillna(0.0)
    )

    pair_order = {label: idx for idx, label in enumerate(_sorted_pair_labels(summary["pair_label"]))}
    summary["_pair_sort"] = summary["pair_label"].map(pair_order)
    summary = summary.sort_values(
        ["_pair_sort", "pair_label", "n_cell_types", "ppi_priority_score", "min_padj", "gene_symbol"],
        ascending=[True, True, False, False, True, True],
    )
    summary["seed_rank"] = summary.groupby("pair_label").cumcount() + 1
    summary = summary.drop(columns=["_pair_sort"]).reset_index(drop=True)
    summary = summary.reindex(columns=template_columns)
    return summary


def write_ppi_gene_list_files(
    summary: pd.DataFrame,
    out_dir: Path,
    *,
    scope_label: str,
    top_n: int,
) -> pd.DataFrame:
    """Write one-gene-per-line files for STRING/Cytoscape import."""
    if summary.empty:
        return pd.DataFrame(columns=["scope", "pair_label", "gene_set_kind", "n_genes", "txt_path"])

    target_dir = out_dir / "ppi_input" / scope_label
    target_dir.mkdir(parents=True, exist_ok=True)
    manifest_records: list[dict[str, object]] = []
    for pair_label in _sorted_pair_labels(summary["pair_label"]):
        sub = summary[summary["pair_label"].astype(str) == pair_label].sort_values("seed_rank").copy()
        if sub.empty:
            continue
        gene_sets = {
            "all_sig": sub["string_input_gene"].astype(str).tolist(),
            f"top{int(top_n)}": sub.loc[sub["seed_rank"] <= int(top_n), "string_input_gene"].astype(str).tolist(),
        }
        for gene_set_kind, genes in gene_sets.items():
            genes = [gene for gene in genes if _clean_text(gene)]
            if not genes:
                continue
            txt_path = target_dir / f"{pair_label}_{gene_set_kind}_genes.txt"
            txt_path.write_text("\n".join(genes) + "\n", encoding="utf-8")
            manifest_records.append(
                {
                    "scope": scope_label,
                    "pair_label": pair_label,
                    "gene_set_kind": gene_set_kind,
                    "n_genes": len(genes),
                    "txt_path": str(txt_path),
                }
            )
    return pd.DataFrame(manifest_records)


def generate_ppi_bridge_outputs(
    out_dir: Path,
    *,
    ppi_deg_l3: pd.DataFrame,
    ppi_top_n: int = 150,
) -> tuple[dict[str, str], pd.DataFrame]:
    """Create PPI-ready DEG tables and importable gene-list files."""
    out_dir.mkdir(parents=True, exist_ok=True)
    output_paths: dict[str, str] = {}

    ppi_deg_path = out_dir / "ppi_ready_deg_L3.tsv"
    ppi_deg_l3.to_csv(ppi_deg_path, sep="\t", index=False)
    output_paths["ppi_ready_deg_L3"] = str(ppi_deg_path)

    site_summary = summarize_ppi_seed_candidates(ppi_deg_l3, level="L3", pair_col="canonical_site_pair")
    group_summary = summarize_ppi_seed_candidates(ppi_deg_l3, level="L3", pair_col="canonical_comparison_pair")
    site_top = site_summary.loc[site_summary["seed_rank"] <= int(ppi_top_n)].copy() if not site_summary.empty else site_summary.copy()
    group_top = group_summary.loc[group_summary["seed_rank"] <= int(ppi_top_n)].copy() if not group_summary.empty else group_summary.copy()

    site_summary_path = out_dir / "ppi_seed_summary_L3_by_site.tsv"
    site_top_path = out_dir / f"ppi_seed_summary_L3_by_site_top{int(ppi_top_n)}.tsv"
    group_summary_path = out_dir / "ppi_seed_summary_L3_by_group.tsv"
    group_top_path = out_dir / f"ppi_seed_summary_L3_by_group_top{int(ppi_top_n)}.tsv"
    site_summary.to_csv(site_summary_path, sep="\t", index=False)
    site_top.to_csv(site_top_path, sep="\t", index=False)
    group_summary.to_csv(group_summary_path, sep="\t", index=False)
    group_top.to_csv(group_top_path, sep="\t", index=False)
    output_paths["ppi_seed_summary_L3_by_site"] = str(site_summary_path)
    output_paths[f"ppi_seed_summary_L3_by_site_top{int(ppi_top_n)}"] = str(site_top_path)
    output_paths["ppi_seed_summary_L3_by_group"] = str(group_summary_path)
    output_paths[f"ppi_seed_summary_L3_by_group_top{int(ppi_top_n)}"] = str(group_top_path)

    site_manifest = write_ppi_gene_list_files(site_summary, out_dir, scope_label="l3_pairwise", top_n=ppi_top_n)
    group_manifest = write_ppi_gene_list_files(group_summary, out_dir, scope_label="l3_grouped", top_n=ppi_top_n)
    ppi_manifest = pd.concat([site_manifest, group_manifest], ignore_index=True) if not site_manifest.empty or not group_manifest.empty else pd.DataFrame(columns=["scope", "pair_label", "gene_set_kind", "n_genes", "txt_path"])
    ppi_manifest_path = out_dir / "ppi_input_manifest_L3.tsv"
    ppi_manifest.to_csv(ppi_manifest_path, sep="\t", index=False)
    output_paths["ppi_input_manifest_L3"] = str(ppi_manifest_path)
    return output_paths, ppi_manifest


def _render_pairwise_lines(pairwise_counts: pd.DataFrame, level: str) -> list[str]:
    if pairwise_counts.empty:
        return [f"- {level} 没有发现可复用的 pairwise DE 对比。"]
    lines = []
    for _, row in pairwise_counts.iterrows():
        lines.append(
            f"- {level} `{row['pair_label']}`: `{int(row['n_cell_types'])}` 个 cell types / "
            f"`{int(row['n_contrast_files'])}` 个 contrast 文件"
        )
    return lines


def summarize_cnmf_status(cnmf_dir: Path | str) -> pd.DataFrame:
    """Collect existing cNMF status JSONs without reading usage matrices."""
    cnmf_dir = Path(cnmf_dir)
    records: list[dict[str, object]] = []
    for status_path in sorted(cnmf_dir.glob("cnmf_by_celltype/*/cnmf_full/cnmf_l3_status.json")):
        try:
            data = json.loads(status_path.read_text())
        except Exception as exc:
            records.append({"status_path": str(status_path), "status": "read_error", "error": str(exc)})
            continue
        records.append(
            {
                "cell_type": data.get("celltype_l3") or status_path.parents[1].name.replace("cnmf_", ""),
                "status": data.get("status"),
                "result_success": data.get("result_success"),
                "n_cells": data.get("n_cells"),
                "n_samples": data.get("n_samples"),
                "n_tissues": data.get("n_tissues"),
                "recommended_k": data.get("recommended_k"),
                "n_score_files": data.get("n_score_files"),
                "n_k_gep_units": data.get("n_k_gep_units"),
                "elapsed_min": data.get("elapsed_min"),
                "output_dir": data.get("output_dir"),
                "status_path": str(status_path),
            }
        )
    return pd.DataFrame(records)


def write_report(
    out_dir: Path,
    h5ad_info: Mapping[str, object],
    summaries: Mapping[str, pd.DataFrame],
    de_summary_l3: pd.DataFrame,
    de_pairwise_site_l3: pd.DataFrame,
    de_pairwise_group_l3: pd.DataFrame,
    key_gene_de: pd.DataFrame,
    cnmf_summary: pd.DataFrame,
    visual_manifest: pd.DataFrame | None = None,
    ppi_deg_l3: pd.DataFrame | None = None,
    ppi_input_manifest: pd.DataFrame | None = None,
    ppi_top_n: int = 150,
) -> Path:
    """Write a compact markdown report for the reproduction handoff."""
    site_counts = summaries["cell_counts_by_site"]
    site_lines = []
    for _, row in site_counts.iterrows():
        site_lines.append(
            f"- `{row['tissue']}` → `{row['site_axis']}/{row['site_group']}` "
            f"(`comparison_group={row['comparison_group']}`): {int(row['n_cells']):,} cells"
        )

    key_gene_n = 0 if key_gene_de.empty else key_gene_de.shape[0]
    cnmf_ok = 0 if cnmf_summary.empty else int((cnmf_summary.get("status") == "ok").sum())
    pairwise_site_lines = _render_pairwise_lines(de_pairwise_site_l3, "L3")
    pairwise_group_lines = _render_pairwise_lines(de_pairwise_group_l3, "L3 grouped")
    ppi_deg_rows = 0 if ppi_deg_l3 is None else int(ppi_deg_l3.shape[0])
    visual_lines: list[str] = []
    if visual_manifest is not None and not visual_manifest.empty:
        for _, row in visual_manifest.iterrows():
            png_name = Path(_clean_text(row.get("png_path"))).name if _clean_text(row.get("png_path")) else ""
            pdf_name = Path(_clean_text(row.get("pdf_path"))).name if _clean_text(row.get("pdf_path")) else ""
            suffix = f" (`{png_name}` / `{pdf_name}`)" if png_name or pdf_name else ""
            visual_lines.append(f"- `{row['figure_id']}`: {row['description']}{suffix}")
    ppi_lines: list[str] = []
    if ppi_input_manifest is not None and not ppi_input_manifest.empty:
        for _, row in ppi_input_manifest.iterrows():
            txt_name = Path(_clean_text(row.get("txt_path"))).name if _clean_text(row.get("txt_path")) else ""
            ppi_lines.append(
                f"- `{row['scope']}` / `{row['pair_label']}` / `{row['gene_set_kind']}`: "
                f"`{int(row['n_genes'])}` genes (`{txt_name}`)"
            )

    text = [
        "# Upper/Lower Airway Reproduction Handoff (2026-05-23)",
        "",
        "## Input objects",
        "",
        f"- Primary epithelial h5ad: `{h5ad_info['path']}`",
        f"- Shape: `{h5ad_info['n_obs']:,} cells × {h5ad_info['n_vars']:,} genes`",
        f"- Layers: `{', '.join(h5ad_info.get('layers', []))}`; counts layer present: `{h5ad_info.get('has_counts_layer')}`",
        "",
        "## Harmonized tissue/site mapping",
        "",
        *site_lines,
        "",
        "## Reused existing evidence",
        "",
        f"- Existing L3 DESeq2 result rows summarized: `{0 if de_summary_l3.empty else de_summary_l3.shape[0]}` cell-type/contrast files.",
        "- L3 detailed pairwise coverage:",
        *pairwise_site_lines,
        "- L3 grouped pairwise coverage:",
        *pairwise_group_lines,
        f"- Key receptor/antiviral/mucosal genes extracted from DE tables: `{key_gene_n}` rows.",
        f"- Existing L3 cNMF status files with `ok`: `{cnmf_ok}`.",
        "",
        "## Interpretation guardrails",
        "",
        "- `nasal` combines `nose` and `sinus` for the high-level upper-airway comparison, while `site_group` preserves `sinus` separately.",
        "- `bronchial` maps the existing `respiratory airway` tissue to lower conducting airway.",
        "- `lung parenchyma` is retained as `distal_lung/lung` and should support divergence claims, not the main nasal-vs-bronchial proxy claim.",
        "- Existing DESeq2 results are reused as prior pipeline output; the new script does not rerun full pseudobulk DE by default.",
        "",
        "## DEG → PPI bridge",
        "",
        f"- PPI-ready L3 DEG rows exported: `{ppi_deg_rows}`.",
        "- `ppi_ready_deg_L3.tsv` keeps one significant DEG row per cell type / contrast / gene with harmonized pair labels and an effect-size-aware `ppi_priority_score`.",
        "- `ppi_seed_summary_L3_by_site.tsv` aggregates genes within each detailed site pair (e.g. `nose_vs_bronchial`) across L3 cell types.",
        "- `ppi_seed_summary_L3_by_group.tsv` repeats the aggregation after collapsing `nose + sinus -> nasal`, so downstream PPI can align with the upper-airway claim.",
        f"- `ppi_seed_summary_*_top{int(ppi_top_n)}.tsv` provides the default top-ranked seed panels for STRING/Cytoscape upload.",
        "- `ppi_input_manifest_L3.tsv` indexes one-gene-per-line `.txt` files that can be pasted directly into STRINGdb or imported into Cytoscape.",
        "- Caveat: RNA differential expression is **not** a direct PPI inference; these outputs only map DEG-encoded proteins onto known interaction databases such as STRING/BioGRID/HINT.",
        "- Recommended handoff: choose a pair -> take the corresponding `topN` gene list -> build the PPI graph in STRING/Cytoscape -> run hub-gene ranking (e.g. degree / MCC) -> return the hubs to the DEG/key-gene tables for expression-side validation.",
        *(ppi_lines if ppi_lines else ["- 本次运行未生成 PPI gene-list files；请检查 `ppi_input_manifest_L3.tsv`。"]),
        "",
        "## DEG summary visuals",
        "",
        *(visual_lines if visual_lines else ["- 本次运行未生成 DEG summary 图；请检查 `deg_visual_manifest.tsv`。"]),
        "",
        "## Generated tables",
        "",
        "- `metadata_upper_lower_airway.tsv.gz`",
        "- `cell_counts_by_site.tsv`",
        "- `sample_site_summary.tsv`",
        "- `celltype_site_summary.tsv`",
        "- `celltype_sample_site_counts.tsv`",
        "- `existing_deseq2_L2_summary.tsv` / `existing_deseq2_L3_summary.tsv`",
        "- `existing_deseq2_L2_pairwise_catalog.tsv` / `existing_deseq2_L3_pairwise_catalog.tsv`",
        "- `existing_deseq2_L2_pairwise_site_counts.tsv` / `existing_deseq2_L3_pairwise_site_counts.tsv`",
        "- `existing_deseq2_L2_pairwise_comparison_counts.tsv` / `existing_deseq2_L3_pairwise_comparison_counts.tsv`",
        "- `deg_pairwise_site_coverage_by_level.tsv`",
        "- `deg_l3_nsig_matrix.tsv`",
        "- `deg_l3_direction_balance_matrix.tsv`",
        "- `deg_l3_top_up_gene_recurrence.tsv` / `deg_l3_top_down_gene_recurrence.tsv`",
        "- `key_gene_existing_deseq2.tsv`",
        "- `ppi_ready_deg_L3.tsv`",
        "- `ppi_seed_summary_L3_by_site.tsv` / `ppi_seed_summary_L3_by_group.tsv`",
        f"- `ppi_seed_summary_L3_by_site_top{int(ppi_top_n)}.tsv` / `ppi_seed_summary_L3_by_group_top{int(ppi_top_n)}.tsv`",
        "- `ppi_input_manifest_L3.tsv`",
        "- `deg_visual_manifest.tsv`",
        "- `cnmf_l3_status_summary.tsv`",
        "- `preflight_summary.json`",
        "",
    ]
    report_path = out_dir / "UPPER_LOWER_AIRWAY_REPRODUCTION_HANDOFF.md"
    report_path.write_text("\n".join(text), encoding="utf-8")
    return report_path


def run_workflow(
    h5ad_path: Path,
    existing_r_dir: Path,
    cnmf_dir: Path,
    out_dir: Path,
    celltype_col: str = "cell_type_L3",
    generate_figures: bool = True,
    ppi_top_n: int = 150,
) -> dict[str, str]:
    """Run the lightweight reproduction bridge and write output files."""
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ad_info = inspect_h5ad_light(h5ad_path)
    metadata = read_obs_columns_h5py(h5ad_path, OBS_COLUMNS)
    metadata = add_site_annotations(metadata)
    summaries = summarize_metadata(metadata, celltype_col=celltype_col)
    gene_map = read_var_gene_map_h5py(h5ad_path)

    output_paths: dict[str, str] = {}
    metadata_path = out_dir / "metadata_upper_lower_airway.tsv.gz"
    metadata.to_csv(metadata_path, sep="\t", compression="gzip")
    output_paths["metadata"] = str(metadata_path)

    for name, table in summaries.items():
        path = out_dir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        output_paths[name] = str(path)

    de_l2_raw, key_l2 = summarize_deseq2_tree(existing_r_dir / "pseudobulk_de", "L2", gene_map=gene_map)
    de_l3_raw, key_l3 = summarize_deseq2_tree(existing_r_dir / "pseudobulk_de_L3", "L3", gene_map=gene_map)
    de_l2 = annotate_pairwise_columns(de_l2_raw)
    de_l3 = annotate_pairwise_columns(de_l3_raw)
    de_l2_path = out_dir / "existing_deseq2_L2_summary.tsv"
    de_l3_path = out_dir / "existing_deseq2_L3_summary.tsv"
    de_l2.to_csv(de_l2_path, sep="\t", index=False)
    de_l3.to_csv(de_l3_path, sep="\t", index=False)
    output_paths["existing_deseq2_L2_summary"] = str(de_l2_path)
    output_paths["existing_deseq2_L3_summary"] = str(de_l3_path)

    de_l2_pairwise_catalog = de_l2.copy()
    de_l3_pairwise_catalog = de_l3.copy()
    de_l2_pairwise_site_counts = summarize_pairwise_catalog(de_l2_pairwise_catalog, "canonical_site_pair")
    de_l3_pairwise_site_counts = summarize_pairwise_catalog(de_l3_pairwise_catalog, "canonical_site_pair")
    de_l2_pairwise_group_counts = summarize_pairwise_catalog(de_l2_pairwise_catalog, "canonical_comparison_pair")
    de_l3_pairwise_group_counts = summarize_pairwise_catalog(de_l3_pairwise_catalog, "canonical_comparison_pair")

    pairwise_tables = {
        "existing_deseq2_L2_pairwise_catalog": de_l2_pairwise_catalog,
        "existing_deseq2_L3_pairwise_catalog": de_l3_pairwise_catalog,
        "existing_deseq2_L2_pairwise_site_counts": de_l2_pairwise_site_counts,
        "existing_deseq2_L3_pairwise_site_counts": de_l3_pairwise_site_counts,
        "existing_deseq2_L2_pairwise_comparison_counts": de_l2_pairwise_group_counts,
        "existing_deseq2_L3_pairwise_comparison_counts": de_l3_pairwise_group_counts,
    }
    for name, table in pairwise_tables.items():
        path = out_dir / f"{name}.tsv"
        table.to_csv(path, sep="\t", index=False)
        output_paths[name] = str(path)

    key_gene_de = pd.concat([key_l2, key_l3], ignore_index=True) if not key_l2.empty or not key_l3.empty else pd.DataFrame()
    key_gene_de = annotate_pairwise_columns(key_gene_de)
    key_gene_path = out_dir / "key_gene_existing_deseq2.tsv"
    key_gene_de.to_csv(key_gene_path, sep="\t", index=False)
    output_paths["key_gene_existing_deseq2"] = str(key_gene_path)

    ppi_deg_l3 = summarize_sig_deg_table(existing_r_dir / "pseudobulk_de_L3", "L3", gene_map=gene_map)
    ppi_deg_l3 = annotate_pairwise_columns(ppi_deg_l3)
    ppi_outputs, ppi_input_manifest = generate_ppi_bridge_outputs(
        out_dir,
        ppi_deg_l3=ppi_deg_l3,
        ppi_top_n=ppi_top_n,
    )
    output_paths.update(ppi_outputs)

    visual_manifest = pd.DataFrame()
    if generate_figures:
        visual_outputs, visual_manifest = generate_deg_summary_visuals(
            out_dir,
            de_summary_l2=de_l2,
            de_summary_l3=de_l3,
            de_pairwise_site_l2=de_l2_pairwise_site_counts,
            de_pairwise_site_l3=de_l3_pairwise_site_counts,
            key_gene_de=key_gene_de,
        )
        output_paths.update(visual_outputs)

    cnmf_summary = summarize_cnmf_status(cnmf_dir)
    cnmf_path = out_dir / "cnmf_l3_status_summary.tsv"
    cnmf_summary.to_csv(cnmf_path, sep="\t", index=False)
    output_paths["cnmf_l3_status_summary"] = str(cnmf_path)

    preflight = {
        "h5ad": h5ad_info,
        "obs_columns_read": [col for col in OBS_COLUMNS if col in metadata.columns],
        "site_counts": summaries["cell_counts_by_site"].to_dict(orient="records"),
        "n_cells_metadata": int(metadata.shape[0]),
        "n_samples": int(metadata["sample"].nunique()) if "sample" in metadata.columns else None,
        "n_cell_types": int(metadata[celltype_col].nunique()) if celltype_col in metadata.columns else None,
        "existing_de_l2_files": int(de_l2.shape[0]),
        "existing_de_l3_files": int(de_l3.shape[0]),
        "existing_de_l3_pairwise_site_counts": de_l3_pairwise_site_counts.to_dict(orient="records"),
        "existing_de_l3_pairwise_comparison_counts": de_l3_pairwise_group_counts.to_dict(orient="records"),
        "key_gene_de_rows": int(key_gene_de.shape[0]),
        "ppi_ready_deg_l3_rows": int(ppi_deg_l3.shape[0]),
        "ppi_input_files": ppi_input_manifest.to_dict(orient="records") if not ppi_input_manifest.empty else [],
        "deg_visual_figures": visual_manifest.to_dict(orient="records") if not visual_manifest.empty else [],
        "cnmf_status_files": int(cnmf_summary.shape[0]),
    }
    preflight_path = out_dir / "preflight_summary.json"
    preflight_path.write_text(json.dumps(preflight, indent=2), encoding="utf-8")
    output_paths["preflight_summary"] = str(preflight_path)

    report_path = write_report(
        out_dir,
        h5ad_info,
        summaries,
        de_l3,
        de_l3_pairwise_site_counts,
        de_l3_pairwise_group_counts,
        key_gene_de,
        cnmf_summary,
        visual_manifest=visual_manifest,
        ppi_deg_l3=ppi_deg_l3,
        ppi_input_manifest=ppi_input_manifest,
        ppi_top_n=ppi_top_n,
    )
    output_paths["report"] = str(report_path)
    return output_paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD, help="Primary epithelial h5ad with full-gene counts layer.")
    parser.add_argument("--existing-r-dir", type=Path, default=DEFAULT_EXISTING_R_DIR, help="Existing epithelial tissue-comparison output directory.")
    parser.add_argument("--cnmf-dir", type=Path, default=DEFAULT_CNMF_DIR, help="Existing epithelial program/cNMF output directory.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for reproduction bridge tables.")
    parser.add_argument("--celltype-col", default="cell_type_L3", help="Cell type column to summarize.")
    parser.add_argument("--ppi-top-n", type=int, default=150, help="Default top-N genes exported per pair for downstream PPI upload.")
    parser.add_argument("--skip-figures", action="store_true", help="Skip DEG summary figure generation.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    outputs = run_workflow(
        h5ad_path=args.h5ad,
        existing_r_dir=args.existing_r_dir,
        cnmf_dir=args.cnmf_dir,
        out_dir=args.out_dir,
        celltype_col=args.celltype_col,
        generate_figures=not args.skip_figures,
        ppi_top_n=args.ppi_top_n,
    )
    print(json.dumps(outputs, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
