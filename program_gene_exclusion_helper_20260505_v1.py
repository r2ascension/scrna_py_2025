#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared gene exclusion helper for program-analysis modules."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROGRAM_GENE_EXCLUSION_HELPER_VERSION_20260505_V1 = "20260505_v1"
PROGRAM_GENE_EXCLUSION_SPEC_PATH_20260505_V1 = Path("/home/h2048/script/config/program_gene_exclusion_spec_20260505_v1.json")


def safe_name(x: str, max_len: int = 180) -> str:
    s = str(x)
    for ch in ['/', '\\', ' ', '|', ':', ';', ',', '\t', "'", '"']:
        s = s.replace(ch, '_')
    return s[:max_len]


def safe_mkdir(directory: Path) -> bool:
    try:
        directory.mkdir(parents=True, exist_ok=True)
        return directory.is_dir()
    except Exception:
        return False


def load_gene_exclusion_spec(
    spec_path: Path = PROGRAM_GENE_EXCLUSION_SPEC_PATH_20260505_V1,
    overrides: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    spec.update(overrides or {})
    return spec


def _norm(values: Iterable[Any]) -> List[str]:
    return [str(v).strip() for v in values]


def _upper(values: Iterable[Any]) -> List[str]:
    return [str(v).strip().upper() for v in values]


def _match_prefix_any(values: Iterable[str], prefixes: Iterable[str]) -> np.ndarray:
    value_list = _upper(values)
    prefix_list = _upper(prefixes)
    out = np.zeros(len(value_list), dtype=bool)
    for prefix in prefix_list:
        if not prefix:
            continue
        out |= np.array([val.startswith(prefix) for val in value_list], dtype=bool)
    return out


def _match_regex_any(values: Iterable[str], patterns: Iterable[str]) -> np.ndarray:
    import re

    value_list = _norm(values)
    out = np.zeros(len(value_list), dtype=bool)
    for pattern in _norm(patterns):
        if not pattern:
            continue
        regex = re.compile(pattern, flags=re.IGNORECASE)
        out |= np.array([bool(regex.search(val)) for val in value_list], dtype=bool)
    return out


def _normalize_annotation_df(annotation_df: Optional[pd.DataFrame], spec: Dict[str, Any]) -> Optional[pd.DataFrame]:
    if annotation_df is None or annotation_df.empty:
        return None
    gene_col = next((c for c in spec.get("annotation_gene_column_candidates", []) if c in annotation_df.columns), None)
    biotype_col = next((c for c in spec.get("annotation_column_candidates", []) if c in annotation_df.columns), None)
    if gene_col is None or biotype_col is None:
        return None
    out = annotation_df[[gene_col, biotype_col]].copy()
    out.columns = ["gene_symbol", "annotation_biotype"]
    out["gene_symbol"] = out["gene_symbol"].astype(str).str.strip()
    out["annotation_biotype"] = out["annotation_biotype"].astype(str).str.strip()
    out = out.loc[out["gene_symbol"] != ""].drop_duplicates(subset=["gene_symbol"])
    if out.empty:
        return None
    out["gene_symbol_upper"] = out["gene_symbol"].str.upper()
    out = out.drop_duplicates(subset=["gene_symbol_upper"]).set_index("gene_symbol_upper", drop=False)
    return out


def read_annotation_df(annotation_path: Optional[str], spec: Dict[str, Any]) -> Optional[pd.DataFrame]:
    if not annotation_path:
        return None
    path = Path(str(annotation_path)).expanduser()
    if not path.exists():
        return None
    if path.suffix.lower() in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path)
    return _normalize_annotation_df(df, spec)


def infer_lineage_context(
    adata,
    celltype_col: Optional[str] = None,
    extra_obs_cols: Optional[List[str]] = None,
) -> str:
    cols = []
    if celltype_col:
        cols.append(celltype_col)
    cols.extend(extra_obs_cols or ["lineage", "lineage_tag", "cell_type_L2", "cell_type_L3", "cell_type_final_l2", "cell_type_final_l3"])
    ctx: List[str] = []
    for col in cols:
        if col in adata.obs.columns:
            vals = sorted({str(v).strip() for v in adata.obs[col].dropna().astype(str).tolist() if str(v).strip()})
            ctx.extend(vals)
    if not ctx:
        return ""
    uniq: List[str] = []
    for val in ctx:
        if val not in uniq:
            uniq.append(val)
    return " | ".join(uniq[:25])


def lineage_context_keeps_category(category: str, lineage_context: str, spec: Dict[str, Any]) -> bool:
    category_key = str(category).strip().lower()
    keep_map = spec.get("category_keep_lineage_patterns", {}) or {}
    if category_key == "immunoglobulin":
        policy = str(spec.get("ig_policy", "auto")).strip().lower()
        keep_patterns = keep_map.get("immunoglobulin") or spec.get("ig_keep_lineage_patterns", [])
    elif category_key == "secretory":
        policy = str(spec.get("secretory_policy", "auto")).strip().lower()
        keep_patterns = keep_map.get("secretory") or spec.get("secretory_keep_lineage_patterns", [])
    else:
        return False
    if policy == "keep":
        return True
    if policy == "exclude":
        return False
    if not str(lineage_context).strip():
        return False
    return bool(_match_regex_any([lineage_context], keep_patterns)[0])


def lineage_context_keeps_ig(lineage_context: str, spec: Dict[str, Any]) -> bool:
    return lineage_context_keeps_category("immunoglobulin", lineage_context, spec)


def build_gene_exclusion_packet(
    gene_names: Iterable[str],
    spec: Optional[Dict[str, Any]] = None,
    lineage_context: str = "",
    annotation_df: Optional[pd.DataFrame] = None,
    annotation_path: Optional[str] = None,
) -> Dict[str, Any]:
    spec = load_gene_exclusion_spec(overrides=spec or {}) if not (spec and spec.get("spec_version")) else dict(spec)
    genes = [g for g in _norm(gene_names) if g]
    if not genes:
        raise ValueError("gene_names must contain at least one non-empty symbol")

    ann = _normalize_annotation_df(annotation_df, spec)
    if ann is None:
        ann = read_annotation_df(annotation_path, spec)

    genes_upper = _upper(genes)
    annotation_biotype = np.array([np.nan] * len(genes), dtype=object)
    if ann is not None:
        hits = ann.reindex(genes_upper)
        annotation_biotype = hits["annotation_biotype"].to_numpy(dtype=object)

    mt_match = _match_prefix_any(genes, spec.get("mitochondrial_prefixes", []))
    ribo_match = _match_prefix_any(genes, spec.get("ribosomal_prefixes", []))
    ig_match = _match_regex_any(genes, spec.get("ig_patterns", []))
    ensembl_match = _match_regex_any(genes, spec.get("ensembl_gene_patterns", []))

    lnc_biotypes = {str(x).strip().lower() for x in spec.get("lncrna_biotypes", [])}
    lnc_annotation_match = np.array([
        isinstance(v, str) and v.strip().lower() in lnc_biotypes for v in annotation_biotype
    ], dtype=bool)
    lnc_heuristic_match = _match_regex_any(genes, spec.get("lncrna_symbol_patterns", []))
    lnc_match = lnc_annotation_match | lnc_heuristic_match
    secretory_match = _match_regex_any(genes, spec.get("secretory_symbol_patterns", []))

    keep_ig = lineage_context_keeps_category("immunoglobulin", lineage_context, spec)
    keep_secretory = lineage_context_keeps_category("secretory", lineage_context, spec)
    exclude_mt = bool(spec.get("exclude_mt", True))
    exclude_ribo = bool(spec.get("exclude_ribo", True))
    exclude_lnc = bool(spec.get("exclude_lncRNA", True))
    exclude_ensembl = bool(spec.get("exclude_ensembl_ids", True))
    exclude_secretory = bool(spec.get("exclude_secretory", True))

    category_priority = [
        str(cat).strip().lower()
        for cat in spec.get(
            "category_priority",
            ["mitochondrial", "ribosomal", "immunoglobulin", "ensembl_id", "lncrna", "secretory", "other"],
        )
        if str(cat).strip()
    ]
    category_priority = list(dict.fromkeys(category_priority))
    non_other_priority = [cat for cat in category_priority if cat != "other"]

    category_matches = {
        "mitochondrial": mt_match,
        "ribosomal": ribo_match,
        "immunoglobulin": ig_match,
        "ensembl_id": ensembl_match,
        "lncrna": lnc_match,
        "secretory": secretory_match,
    }
    category_can_exclude = {
        "mitochondrial": exclude_mt,
        "ribosomal": exclude_ribo,
        "immunoglobulin": not keep_ig,
        "ensembl_id": exclude_ensembl,
        "lncrna": exclude_lnc,
        "secretory": exclude_secretory and (not keep_secretory),
    }

    should_exclude = np.zeros(len(genes), dtype=bool)
    exclude_reason = np.array(["kept"] * len(genes), dtype=object)
    rule_source = np.array(["none"] * len(genes), dtype=object)

    for category in non_other_priority:
        update_idx = category_matches.get(category, np.zeros(len(genes), dtype=bool)) & category_can_exclude.get(category, False) & ~should_exclude
        if not np.any(update_idx):
            continue
        should_exclude[update_idx] = True
        exclude_reason[update_idx] = category
        if category == "lncrna":
            rule_source[update_idx] = np.where(lnc_annotation_match[update_idx], "annotation", "heuristic")
        else:
            rule_source[update_idx] = "heuristic"

    kept_ig_idx = ig_match & keep_ig & ~should_exclude
    exclude_reason[kept_ig_idx] = "kept_immunoglobulin_due_to_lineage"
    rule_source[kept_ig_idx] = "heuristic"

    kept_secretory_idx = secretory_match & keep_secretory & ~should_exclude
    exclude_reason[kept_secretory_idx] = "kept_secretory_due_to_lineage"
    rule_source[kept_secretory_idx] = "heuristic"

    matched_categories: List[str] = []
    primary_category: List[str] = []
    for i in range(len(genes)):
        cats: List[str] = []
        for category in non_other_priority:
            match_vec = category_matches.get(category)
            if match_vec is not None and bool(match_vec[i]):
                cats.append(category)
        if not cats:
            cats = ["other"]
        matched_categories.append("|".join(cats))
        primary_category.append(cats[0])

    audit_df = pd.DataFrame(
        {
            "gene_symbol": genes,
            "gene_symbol_upper": genes_upper,
            "primary_category": primary_category,
            "matched_categories": matched_categories,
            "should_exclude": should_exclude,
            "exclude_reason": exclude_reason,
            "rule_source": rule_source,
            "annotation_biotype": annotation_biotype,
        }
    )

    max_examples = int(spec.get("max_example_genes", 12) or 12)
    summary_rows: List[Dict[str, Any]] = []
    excluded_total = int(audit_df["should_exclude"].sum())
    for category in category_priority:
        cat_df = audit_df.loc[audit_df["primary_category"] == category]
        excl_df = cat_df.loc[cat_df["should_exclude"]]
        summary_rows.append(
            {
                "category": category,
                "detected_n": int(len(cat_df)),
                "excluded_n": int(len(excl_df)),
                "kept_n": int(len(cat_df) - len(excl_df)),
                "pct_of_all_features": round(100 * len(cat_df) / max(1, len(audit_df)), 3),
                "pct_of_all_excluded": round(100 * len(excl_df) / max(1, excluded_total), 3),
                "example_genes": ";".join(cat_df["gene_symbol"].head(max_examples).tolist()),
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    keep_genes = audit_df.loc[~audit_df["should_exclude"], "gene_symbol"].tolist()
    exclude_genes = audit_df.loc[audit_df["should_exclude"], "gene_symbol"].tolist()
    manifest = {
        "spec_version": spec.get("spec_version", PROGRAM_GENE_EXCLUSION_HELPER_VERSION_20260505_V1),
        "lineage_context": lineage_context,
        "ig_policy": spec.get("ig_policy", "auto"),
        "keep_ig": keep_ig,
        "secretory_policy": spec.get("secretory_policy", "auto"),
        "keep_secretory": keep_secretory,
        "n_total_features": int(len(audit_df)),
        "n_excluded_features": int(len(exclude_genes)),
        "n_kept_features": int(len(keep_genes)),
        "annotation_used": ann is not None,
    }

    return {
        "audit_df": audit_df,
        "summary_df": summary_df,
        "keep_genes": keep_genes,
        "exclude_genes": exclude_genes,
        "manifest": manifest,
        "spec": spec,
    }


def plot_gene_exclusion_summary(summary_df: pd.DataFrame, output_prefix: Path, title: str) -> Dict[str, str]:
    output_prefix = Path(output_prefix)
    pdf_path = output_prefix.with_suffix(".pdf")
    png_path = output_prefix.with_suffix(".png")
    plot_df = summary_df.copy()
    mat = plot_df[["detected_n", "excluded_n", "kept_n"]].to_numpy(dtype=float).T
    categories = plot_df["category"].tolist()
    rows = ["Detected", "Excluded", "Kept"]

    def draw(save_path: Path, dpi: int = 180) -> None:
        fig, ax = plt.subplots(figsize=(10, 7))
        x = np.arange(len(categories))
        bottom = np.zeros(len(categories), dtype=float)
        colors = ["#9ECAE1", "#FC9272", "#A1D99B"]
        for i, row_name in enumerate(rows):
            vals = mat[i, :]
            ax.bar(x, vals, bottom=bottom, label=row_name, color=colors[i], edgecolor="none")
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(categories, rotation=30, ha="right")
        ax.set_ylabel("Gene count")
        ax.set_title(title)
        for xi, total in zip(x, bottom):
            ax.text(xi, total, f"{int(total)}", ha="center", va="bottom", fontsize=9)
        ax.legend(frameon=False)
        fig.subplots_adjust(left=0.10, right=0.98, bottom=0.23, top=0.90)
        fig.savefig(save_path, dpi=dpi, bbox_inches=None)
        plt.close(fig)

    draw(pdf_path, dpi=180)
    draw(png_path, dpi=180)
    return {"pdf": str(pdf_path), "png": str(png_path)}


def export_gene_exclusion_packet(packet: Dict[str, Any], output_dir: Path, prefix: str = "gene_exclusion") -> Dict[str, Any]:
    output_dir = Path(output_dir)
    if not safe_mkdir(output_dir):
        raise OSError(f"Cannot create output directory: {output_dir}")
    prefix = safe_name(prefix)
    summary_path = output_dir / f"{prefix}_summary.csv"
    audit_path = output_dir / f"{prefix}_audit.csv"
    manifest_path = output_dir / f"{prefix}_manifest.json"
    plot_prefix = output_dir / f"{prefix}_summary"

    packet["summary_df"].to_csv(summary_path, index=False)
    packet["audit_df"].to_csv(audit_path, index=False)
    plot_paths = plot_gene_exclusion_summary(
        packet["summary_df"],
        plot_prefix,
        title=f"Gene exclusion summary ({packet['manifest'].get('lineage_context') or 'lineage:auto'})",
    )
    manifest = dict(packet["manifest"])
    manifest.update(
        {
            "output_dir": str(output_dir),
            "summary_csv": str(summary_path),
            "audit_csv": str(audit_path),
            "summary_plot": plot_paths,
            "manifest_json": str(manifest_path),
            "helper_version": PROGRAM_GENE_EXCLUSION_HELPER_VERSION_20260505_V1,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest


def apply_gene_exclusion_to_adata(
    adata,
    output_dir: Path,
    gene_exclusion_config: Optional[Dict[str, Any]] = None,
    prefix: str = "gene_exclusion",
    celltype_col: Optional[str] = None,
    lineage_context: Optional[str] = None,
):
    config = dict(gene_exclusion_config or {})
    spec = load_gene_exclusion_spec(overrides=config)
    lineage = lineage_context if lineage_context is not None else config.get("lineage_context") or infer_lineage_context(adata, celltype_col=celltype_col)
    packet = build_gene_exclusion_packet(
        gene_names=adata.var_names.tolist(),
        spec=spec,
        lineage_context=lineage,
        annotation_df=config.get("annotation_df"),
        annotation_path=config.get("annotation_path"),
    )
    manifest = export_gene_exclusion_packet(packet, output_dir=Path(output_dir), prefix=prefix)
    keep_mask = adata.var_names.isin(packet["keep_genes"])
    filtered = adata[:, keep_mask].copy()
    filtered.uns.setdefault("program_gene_exclusion", {})
    filtered.uns["program_gene_exclusion"] = {
        "manifest": manifest,
        "summary_table": packet["summary_df"].to_dict(orient="records"),
    }
    return filtered, packet, manifest
