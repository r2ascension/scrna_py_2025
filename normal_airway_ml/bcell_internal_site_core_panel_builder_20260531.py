#!/usr/bin/env python3
"""Build B-cell site-core panels from existing healthy pseudobulk DE outputs.

This script does **not** rerun DE. Instead it reuses the already completed healthy
B-cell `pseudobulk_de` tables and condenses them into stricter, site-anchored
core panels suitable for downstream internal-only query validation.

Key assumptions verified in this repository before implementation:
1. each comparison directory is named as `t2_vs_t1`,
2. DESeq2 positive log2FC means "higher in t2 vs t1",
3. we want site-core differences rather than prediction-centric features, and
4. the disease query currently only supports `nose` and `sinus` tissues, so we
   annotate which healthy panels are directly query-supportable.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

PY_ROOT = Path(__file__).resolve().parents[1]
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from normal_airway_ml_common_20260527 import deep_get, ensure_dir, load_config, write_json  # noqa: E402

CONFIG_KEY = "bcell_internal_site_core_validation"
DEFAULT_CONFIG_PATH = "/home/h2048/script/config/normal_airway_ml_bcell_internal_site_core_validation_20260531.yaml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build strict site-core B-cell panels from healthy pseudobulk DE outputs.")
    parser.add_argument("--config", type=Path, default=Path(DEFAULT_CONFIG_PATH), help="YAML config path.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional explicit run root override.")
    return parser.parse_args()


def slugify_token(text: str) -> str:
    return re.sub(r"[^0-9A-Za-z._-]+", "_", str(text).strip()).strip("_") or "unknown"


def resolve_run_root(cfg: dict[str, Any], output_dir: Path | None) -> Path:
    if output_dir is not None:
        return ensure_dir(output_dir)
    output_root = Path(str(deep_get(cfg, "run", "output_root", default="/home/h2048/data/py/20260531")))
    prefix = str(deep_get(cfg, "run", "run_name_prefix", default="bcell_internal_site_core_validation_20260531"))
    return ensure_dir(output_root / prefix)


def collapse_unique(values: pd.Series) -> str:
    uniq = sorted({str(v) for v in values if pd.notna(v) and str(v) not in {"", "nan", "None", "<NA>"}})
    return "|".join(uniq)


def parse_contrast_name(name: str) -> tuple[str, str]:
    if "_vs_" not in str(name):
        raise ValueError(f"Unsupported contrast directory name: {name}")
    positive_site, negative_site = str(name).split("_vs_", 1)
    return positive_site, negative_site


def load_de_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"gene", "log2FoldChange", "padj"}
    missing = required - set(df.columns)
    if missing:
        raise KeyError(f"{path} is missing required columns: {sorted(missing)}")
    out = df.copy()
    out["gene"] = out["gene"].astype(str).str.strip()
    out = out[out["gene"] != ""].copy()
    out["log2FoldChange"] = pd.to_numeric(out["log2FoldChange"], errors="coerce")
    out["padj"] = pd.to_numeric(out["padj"], errors="coerce")
    out["pvalue"] = pd.to_numeric(out.get("pvalue", np.nan), errors="coerce")
    out["baseMean"] = pd.to_numeric(out.get("baseMean", np.nan), errors="coerce")
    out["abs_log2fc"] = out["log2FoldChange"].abs()
    out = out.dropna(subset=["log2FoldChange", "padj", "abs_log2fc"])
    return out.reset_index(drop=True)


def build_evidence_table(cfg: dict[str, Any]) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    pb_cfg = deep_get(block, "panel_builder", default={}) or {}
    healthy_root = Path(str(deep_get(block, "healthy_de_root", default="")))
    de_filename = str(deep_get(block, "panel_filename", default="DESeq2_results.csv"))
    subtype_whitelist = [str(x) for x in deep_get(pb_cfg, "subtype_whitelist", default=["Naive_B", "Memory_B", "Plasma"])]
    max_padj = float(deep_get(pb_cfg, "max_padj", default=0.05))
    min_abs_log2fc = float(deep_get(pb_cfg, "min_abs_log2fc", default=0.5))

    evidence_frames: list[pd.DataFrame] = []
    discovery_notes: list[dict[str, str]] = []

    for subtype in subtype_whitelist:
        subtype_dir = healthy_root / subtype
        if not subtype_dir.exists():
            discovery_notes.append({"subtype": subtype, "status": "missing_subtype_dir", "path": str(subtype_dir)})
            continue
        de_paths = sorted(subtype_dir.rglob(de_filename))
        if not de_paths:
            discovery_notes.append({"subtype": subtype, "status": "no_de_files", "path": str(subtype_dir)})
            continue

        for de_path in de_paths:
            contrast_name = de_path.parent.name
            try:
                positive_site_from_path, negative_site_from_path = parse_contrast_name(contrast_name)
            except ValueError:
                discovery_notes.append({"subtype": subtype, "status": "skipped_unparsed_contrast", "path": str(de_path)})
                continue

            de = load_de_table(de_path)
            de = de[(de["padj"] <= max_padj) & (de["abs_log2fc"] >= min_abs_log2fc)].copy()
            if de.empty:
                discovery_notes.append({"subtype": subtype, "status": "no_sig_genes_after_filter", "path": str(de_path)})
                continue

            de["subtype"] = subtype
            de["contrast"] = contrast_name
            de["contrast_path"] = str(de_path)
            de["positive_site_from_path"] = positive_site_from_path
            de["negative_site_from_path"] = negative_site_from_path
            de["site_anchor"] = np.where(de["log2FoldChange"] > 0, positive_site_from_path, negative_site_from_path)
            de["opposing_site"] = np.where(de["log2FoldChange"] > 0, negative_site_from_path, positive_site_from_path)
            de["signed_anchor_log2fc"] = np.where(de["log2FoldChange"] > 0, de["log2FoldChange"], -de["log2FoldChange"])
            de["anchor_direction_label"] = np.where(
                de["log2FoldChange"] > 0,
                "higher_in_" + positive_site_from_path + "_vs_" + negative_site_from_path,
                "higher_in_" + negative_site_from_path + "_vs_" + positive_site_from_path,
            )
            evidence_frames.append(
                de[
                    [
                        "subtype",
                        "gene",
                        "contrast",
                        "contrast_path",
                        "site_anchor",
                        "opposing_site",
                        "positive_site_from_path",
                        "negative_site_from_path",
                        "log2FoldChange",
                        "signed_anchor_log2fc",
                        "abs_log2fc",
                        "padj",
                        "pvalue",
                        "baseMean",
                        "anchor_direction_label",
                    ]
                ].copy()
            )

    if not evidence_frames:
        raise RuntimeError("No eligible DE evidence was collected; check healthy_de_root and thresholds.")

    evidence = pd.concat(evidence_frames, axis=0, ignore_index=True)
    evidence = evidence.sort_values(["subtype", "site_anchor", "gene", "padj", "abs_log2fc"], ascending=[True, True, True, True, False]).reset_index(drop=True)
    return evidence, discovery_notes


def build_panel_registry(evidence: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    pb_cfg = deep_get(block, "panel_builder", default={}) or {}
    min_support = int(deep_get(pb_cfg, "min_supporting_contrasts", default=2))
    strict_single_anchor = bool(deep_get(pb_cfg, "strict_single_anchor_site", default=True))
    query_supported_sites = {str(x) for x in deep_get(pb_cfg, "query_supported_sites", default=["nose", "sinus"])}

    grouped = (
        evidence.groupby(["subtype", "gene", "site_anchor"], observed=True)
        .agg(
            supporting_contrasts=("contrast", "nunique"),
            supporting_contrast_list=("contrast", collapse_unique),
            opposing_sites=("opposing_site", collapse_unique),
            median_anchor_log2fc=("signed_anchor_log2fc", "median"),
            median_abs_log2fc=("abs_log2fc", "median"),
            max_abs_log2fc=("abs_log2fc", "max"),
            best_padj=("padj", "min"),
            median_padj=("padj", "median"),
            contrast_paths=("contrast_path", collapse_unique),
        )
        .reset_index()
    )

    gene_totals = (
        grouped.groupby(["subtype", "gene"], observed=True)
        .agg(
            total_supporting_contrasts=("supporting_contrasts", "sum"),
            n_anchor_sites=("site_anchor", "nunique"),
            anchor_sites=("site_anchor", collapse_unique),
            best_supporting_contrasts=("supporting_contrasts", "max"),
        )
        .reset_index()
    )

    registry = (
        grouped.sort_values(
            ["subtype", "gene", "supporting_contrasts", "median_abs_log2fc", "best_padj", "site_anchor"],
            ascending=[True, True, False, False, True, True],
        )
        .drop_duplicates(["subtype", "gene"], keep="first")
        .merge(gene_totals, on=["subtype", "gene"], how="left")
    )

    registry["ambiguous_anchor"] = registry["n_anchor_sites"].fillna(0).astype(int) > 1
    registry["eligible"] = registry["supporting_contrasts"].fillna(0).astype(int) >= min_support
    if strict_single_anchor:
        registry["eligible"] = registry["eligible"] & (~registry["ambiguous_anchor"])
    registry["query_supportable"] = registry["site_anchor"].astype(str).isin(query_supported_sites)
    registry["panel_key"] = registry["subtype"].astype(str) + "__" + registry["site_anchor"].astype(str)
    registry = registry.sort_values(
        ["eligible", "query_supportable", "subtype", "site_anchor", "supporting_contrasts", "median_abs_log2fc", "best_padj", "gene"],
        ascending=[False, False, True, True, False, False, True, True],
    ).reset_index(drop=True)
    return registry


def build_panel_summary(registry: pd.DataFrame) -> pd.DataFrame:
    eligible = registry[registry["eligible"]].copy()
    if eligible.empty:
        return pd.DataFrame(
            columns=[
                "subtype",
                "site_anchor",
                "n_panel_genes",
                "median_supporting_contrasts",
                "median_abs_log2fc",
                "best_padj",
                "query_supportable",
                "top_genes",
            ]
        )

    ranked = eligible.sort_values(
        ["subtype", "site_anchor", "supporting_contrasts", "median_abs_log2fc", "best_padj", "gene"],
        ascending=[True, True, False, False, True, True],
    )
    summary = (
        ranked.groupby(["subtype", "site_anchor", "query_supportable"], observed=True)
        .agg(
            n_panel_genes=("gene", "size"),
            median_supporting_contrasts=("supporting_contrasts", "median"),
            median_abs_log2fc=("median_abs_log2fc", "median"),
            best_padj=("best_padj", "min"),
            top_genes=("gene", lambda s: "|".join(list(pd.Series(s).astype(str).head(8)))),
        )
        .reset_index()
        .sort_values(["subtype", "site_anchor"], ascending=[True, True])
        .reset_index(drop=True)
    )
    return summary


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    run_root = resolve_run_root(cfg, args.output_dir)
    panel_dir = ensure_dir(run_root / "panels")

    evidence, discovery_notes = build_evidence_table(cfg)
    registry = build_panel_registry(evidence, cfg)
    summary = build_panel_summary(registry)

    evidence_path = panel_dir / "panel_evidence.tsv"
    registry_path = panel_dir / "panel_registry.tsv"
    summary_path = panel_dir / "panel_summary.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False)
    registry.to_csv(registry_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    manifest = {
        "config_path": str(args.config),
        "run_root": str(run_root),
        "panel_dir": str(panel_dir),
        "files": {
            "panel_evidence": str(evidence_path),
            "panel_registry": str(registry_path),
            "panel_summary": str(summary_path),
        },
        "counts": {
            "n_evidence_rows": int(evidence.shape[0]),
            "n_registry_rows": int(registry.shape[0]),
            "n_eligible_genes": int(registry[registry["eligible"]].shape[0]),
            "n_query_supportable_genes": int(registry[(registry["eligible"]) & (registry["query_supportable"])].shape[0]),
            "n_panels": int(summary.shape[0]),
        },
        "discovery_notes": discovery_notes,
    }
    write_json(manifest, panel_dir / "panel_manifest.json")

    print(json.dumps(manifest["counts"], ensure_ascii=False, indent=2))
    print(f"[ok] panel registry -> {registry_path}")
    print(f"[ok] panel summary  -> {summary_path}")


if __name__ == "__main__":
    main()
