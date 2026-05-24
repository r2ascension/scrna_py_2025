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


@dataclass(frozen=True)
class SiteAnnotation:
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
        return SiteAnnotation("nasal", "upper_airway", "nasal", True)
    if tissue_s == "sinus":
        return SiteAnnotation("sinus", "upper_airway", "nasal", True)
    if tissue_s == "respiratory airway":
        if "trachea" in level_s:
            site_group = "trachea"
        elif "bronch" in level_s or "lobular" in level_s:
            site_group = "bronchial"
        else:
            site_group = "lower_conducting_airway"
        return SiteAnnotation(site_group, "lower_conducting_airway", "bronchial", True)
    if tissue_s == "lung parenchyma":
        return SiteAnnotation("distal_lung", "distal_lung", "lung", False)
    return SiteAnnotation("unknown", "unknown", "unknown", False)


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
    key_gene_de: pd.DataFrame,
    cnmf_summary: pd.DataFrame,
) -> Path:
    """Write a compact markdown report for the reproduction handoff."""
    site_counts = summaries["cell_counts_by_site"]
    site_lines = []
    for _, row in site_counts.iterrows():
        site_lines.append(
            f"- `{row['tissue']}` → `{row['site_axis']}/{row['site_group']}` "
            f"(`comparison_group={row['comparison_group']}`): {int(row['n_cells']):,} cells"
        )

    airway_de = de_summary_l3[de_summary_l3.get("contrast", pd.Series(dtype=str)).eq("respiratory_airway_vs_nose")]
    distal_de = de_summary_l3[
        de_summary_l3.get("contrast", pd.Series(dtype=str)).isin(
            ["nose_vs_lung_parenchyma", "respiratory_airway_vs_lung_parenchyma"]
        )
    ]
    key_gene_n = 0 if key_gene_de.empty else key_gene_de.shape[0]
    cnmf_ok = 0 if cnmf_summary.empty else int((cnmf_summary.get("status") == "ok").sum())

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
        f"- L3 `respiratory_airway_vs_nose` contrasts available: `{airway_de.shape[0]}`.",
        f"- L3 distal-lung contrasts available: `{distal_de.shape[0]}`.",
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
        "## Generated tables",
        "",
        "- `metadata_upper_lower_airway.tsv.gz`",
        "- `cell_counts_by_site.tsv`",
        "- `sample_site_summary.tsv`",
        "- `celltype_site_summary.tsv`",
        "- `celltype_sample_site_counts.tsv`",
        "- `existing_deseq2_L2_summary.tsv` / `existing_deseq2_L3_summary.tsv`",
        "- `key_gene_existing_deseq2.tsv`",
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

    de_l2, key_l2 = summarize_deseq2_tree(existing_r_dir / "pseudobulk_de", "L2", gene_map=gene_map)
    de_l3, key_l3 = summarize_deseq2_tree(existing_r_dir / "pseudobulk_de_L3", "L3", gene_map=gene_map)
    de_l2_path = out_dir / "existing_deseq2_L2_summary.tsv"
    de_l3_path = out_dir / "existing_deseq2_L3_summary.tsv"
    de_l2.to_csv(de_l2_path, sep="\t", index=False)
    de_l3.to_csv(de_l3_path, sep="\t", index=False)
    output_paths["existing_deseq2_L2_summary"] = str(de_l2_path)
    output_paths["existing_deseq2_L3_summary"] = str(de_l3_path)

    key_gene_de = pd.concat([key_l2, key_l3], ignore_index=True) if not key_l2.empty or not key_l3.empty else pd.DataFrame()
    key_gene_path = out_dir / "key_gene_existing_deseq2.tsv"
    key_gene_de.to_csv(key_gene_path, sep="\t", index=False)
    output_paths["key_gene_existing_deseq2"] = str(key_gene_path)

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
        "key_gene_de_rows": int(key_gene_de.shape[0]),
        "cnmf_status_files": int(cnmf_summary.shape[0]),
    }
    preflight_path = out_dir / "preflight_summary.json"
    preflight_path.write_text(json.dumps(preflight, indent=2), encoding="utf-8")
    output_paths["preflight_summary"] = str(preflight_path)

    report_path = write_report(out_dir, h5ad_info, summaries, de_l3, key_gene_de, cnmf_summary)
    output_paths["report"] = str(report_path)
    return output_paths


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD, help="Primary epithelial h5ad with full-gene counts layer.")
    parser.add_argument("--existing-r-dir", type=Path, default=DEFAULT_EXISTING_R_DIR, help="Existing epithelial tissue-comparison output directory.")
    parser.add_argument("--cnmf-dir", type=Path, default=DEFAULT_CNMF_DIR, help="Existing epithelial program/cNMF output directory.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for reproduction bridge tables.")
    parser.add_argument("--celltype-col", default="cell_type_L3", help="Cell type column to summarize.")
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
    )
    print(json.dumps(outputs, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
