#!/usr/bin/env python3
"""Score-based internal validation for B-cell site-core panels.

Why score-based here?
- the current disease-side `b_cells.h5ad` is not a raw-count AnnData ready for DESeq2,
- `layers` / `raw` are not available as count-backed sources for this object,
- the user explicitly asked for a repository-grounded initial closure using current data,
- so this pilot validates *retention of healthy site-core signatures* inside the disease
  query by scoring a small set of panel genes with low memory overhead.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

REPO_ROOT = Path(__file__).resolve().parents[1]
PY_ROOT = REPO_ROOT / "py"
if str(PY_ROOT) not in sys.path:
    sys.path.insert(0, str(PY_ROOT))

from normal_airway_ml_common_20260527 import deep_get, ensure_dir, load_config, write_json  # noqa: E402

CONFIG_KEY = "bcell_internal_site_core_validation"
DEFAULT_CONFIG_PATH = "/home/h2048/script/config/normal_airway_ml_bcell_internal_site_core_validation_20260531.yaml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run score-based internal validation for B-cell site-core panels.")
    parser.add_argument("--config", type=Path, default=Path(DEFAULT_CONFIG_PATH), help="YAML config path.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Optional explicit run root override.")
    return parser.parse_args()


def resolve_run_root(cfg: dict[str, Any], output_dir: Path | None) -> Path:
    if output_dir is not None:
        return ensure_dir(output_dir)
    output_root = Path(str(deep_get(cfg, "run", "output_root", default="/home/h2048/data/py/20260531")))
    prefix = str(deep_get(cfg, "run", "run_name_prefix", default="bcell_internal_site_core_validation_20260531"))
    return ensure_dir(output_root / prefix)


def dense_float32(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)


def zscore_panel(matrix: np.ndarray, clip_value: float) -> np.ndarray:
    mean = matrix.mean(axis=0, keepdims=True)
    std = matrix.std(axis=0, keepdims=True)
    std[std <= 1e-8] = 1.0
    z = (matrix - mean) / std
    if np.isfinite(clip_value) and clip_value > 0:
        z = np.clip(z, -clip_value, clip_value)
    return z.astype(np.float32, copy=False)


def build_panel_gene_sets(registry: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    val_cfg = deep_get(block, "validation", default={}) or {}
    query_supported_sites = {str(x) for x in deep_get(val_cfg, "query_supported_sites", default=["nose", "sinus"])}
    min_panel_genes = int(deep_get(val_cfg, "min_panel_genes", default=5))

    eligible = registry[(registry["eligible"]) & (registry["query_supportable"])].copy()
    eligible = eligible[eligible["site_anchor"].astype(str).isin(query_supported_sites)].copy()
    if eligible.empty:
        raise RuntimeError("No eligible + query-supportable panel genes were found in panel_registry.tsv")

    panel_sets = (
        eligible.sort_values(["subtype", "site_anchor", "supporting_contrasts", "median_abs_log2fc", "best_padj", "gene"], ascending=[True, True, False, False, True, True])
        .groupby(["subtype", "site_anchor", "panel_key"], observed=True)
        .agg(
            requested_genes=("gene", lambda s: list(pd.Series(s).astype(str))),
            n_requested_genes=("gene", "size"),
            median_supporting_contrasts=("supporting_contrasts", "median"),
            median_abs_log2fc=("median_abs_log2fc", "median"),
        )
        .reset_index()
    )
    panel_sets = panel_sets[panel_sets["n_requested_genes"] >= min_panel_genes].reset_index(drop=True)
    if panel_sets.empty:
        raise RuntimeError("Panels exist, but none passed min_panel_genes after grouping.")
    return panel_sets


def build_gene_coverage(panel_sets: pd.DataFrame, var_names: pd.Index) -> pd.DataFrame:
    var_name_set = set(var_names.astype(str))
    rows: list[dict[str, Any]] = []
    for row in panel_sets.itertuples(index=False):
        requested = [str(g) for g in row.requested_genes]
        present = [g for g in requested if g in var_name_set]
        missing = [g for g in requested if g not in var_name_set]
        rows.append(
            {
                "panel_key": row.panel_key,
                "subtype": row.subtype,
                "site_anchor": row.site_anchor,
                "n_requested_genes": int(len(requested)),
                "n_present_genes": int(len(present)),
                "n_missing_genes": int(len(missing)),
                "requested_genes": "|".join(requested),
                "present_genes": "|".join(present),
                "missing_genes": "|".join(missing),
            }
        )
    return pd.DataFrame(rows)


def compute_panel_scores(adata: ad.AnnData, label_df: pd.DataFrame, gene_coverage: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    val_cfg = deep_get(block, "validation", default={}) or {}
    query_supported_sites = {str(x) for x in deep_get(val_cfg, "query_supported_sites", default=["nose", "sinus"])}
    min_panel_genes = int(deep_get(val_cfg, "min_panel_genes", default=5))
    zscore_clip = float(deep_get(val_cfg, "zscore_clip", default=5.0))

    cell_rows: list[pd.DataFrame] = []
    base_keep = label_df["keep_for_validation"].astype(bool) & label_df["tissue"].astype(str).isin(query_supported_sites)
    if not bool(base_keep.any()):
        raise RuntimeError("No query cells are marked keep_for_validation inside query-supported tissues.")

    for row in gene_coverage.itertuples(index=False):
        present_genes = [g for g in str(row.present_genes).split("|") if g]
        if len(present_genes) < min_panel_genes:
            continue
        subtype_keep = base_keep & (label_df["transferred_subtype_label"].astype(str) == str(row.subtype))
        subtype_cells = label_df.loc[subtype_keep, "obs_name"].astype(str).tolist()
        if not subtype_cells:
            continue

        matrix = adata[subtype_cells, present_genes].X
        dense = dense_float32(matrix)
        z = zscore_panel(dense, clip_value=zscore_clip)
        score = z.mean(axis=1).astype(np.float32)

        block_df = label_df.loc[subtype_keep, [
            "obs_name",
            "sample",
            "dataset",
            "tissue",
            "condition",
            "disease_level_1",
            "disease_level_2",
            "transferred_subtype_label",
        ]].copy()
        block_df["panel_key"] = str(row.panel_key)
        block_df["subtype"] = str(row.subtype)
        block_df["site_anchor"] = str(row.site_anchor)
        block_df["panel_n_genes"] = int(len(present_genes))
        block_df["panel_score"] = score
        cell_rows.append(block_df)

    if not cell_rows:
        raise RuntimeError("No panel scores were computed; check gene coverage and validation cell filters.")
    return pd.concat(cell_rows, axis=0, ignore_index=True)


def summarize_sample_scores(cell_scores: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    val_cfg = deep_get(block, "validation", default={}) or {}
    min_cells_per_sample_subtype = int(deep_get(val_cfg, "min_cells_per_sample_subtype", default=20))

    sample_scores = (
        cell_scores.groupby(
            ["panel_key", "subtype", "site_anchor", "sample", "dataset", "tissue", "condition", "disease_level_1", "disease_level_2"],
            observed=True,
        )
        .agg(
            n_cells=("obs_name", "size"),
            score_mean=("panel_score", "mean"),
            score_median=("panel_score", "median"),
            score_sd=("panel_score", lambda s: float(np.std(np.asarray(s, dtype=float), ddof=0))),
        )
        .reset_index()
    )
    sample_scores = sample_scores[sample_scores["n_cells"] >= min_cells_per_sample_subtype].reset_index(drop=True)
    return sample_scores


def summarize_retention(sample_scores: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    val_cfg = deep_get(block, "validation", default={}) or {}
    min_samples_per_group = int(deep_get(val_cfg, "min_samples_per_group", default=2))

    rows: list[dict[str, Any]] = []
    for (panel_key, subtype, site_anchor), df in sample_scores.groupby(["panel_key", "subtype", "site_anchor"], observed=True):
        matched = df[df["tissue"].astype(str) == str(site_anchor)].copy()
        other = df[df["tissue"].astype(str) != str(site_anchor)].copy()
        match_scores = matched["score_mean"].astype(float).to_numpy()
        other_scores = other["score_mean"].astype(float).to_numpy()
        mw_p = np.nan
        if len(match_scores) >= min_samples_per_group and len(other_scores) >= min_samples_per_group:
            try:
                mw_p = float(mannwhitneyu(match_scores, other_scores, alternative="two-sided").pvalue)
            except Exception:
                mw_p = np.nan
        rows.append(
            {
                "panel_key": panel_key,
                "subtype": subtype,
                "site_anchor": site_anchor,
                "n_matched_samples": int(len(match_scores)),
                "n_other_samples": int(len(other_scores)),
                "matched_mean": float(np.mean(match_scores)) if len(match_scores) else np.nan,
                "other_mean": float(np.mean(other_scores)) if len(other_scores) else np.nan,
                "matched_median": float(np.median(match_scores)) if len(match_scores) else np.nan,
                "other_median": float(np.median(other_scores)) if len(other_scores) else np.nan,
                "delta_match_minus_other": (float(np.mean(match_scores) - np.mean(other_scores)) if len(match_scores) and len(other_scores) else np.nan),
                "mannwhitney_p": mw_p,
            }
        )
    return pd.DataFrame(rows).sort_values(["subtype", "site_anchor"]).reset_index(drop=True)


def summarize_disease_shift(sample_scores: pd.DataFrame, cfg: dict[str, Any]) -> pd.DataFrame:
    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    grouping_cfg = deep_get(block, "grouping", default={}) or {}
    healthy_reference_label = str(deep_get(grouping_cfg, "healthy_reference_label", default="healthy"))

    matched = sample_scores[sample_scores["tissue"].astype(str) == sample_scores["site_anchor"].astype(str)].copy()
    if matched.empty:
        return pd.DataFrame(
            columns=[
                "panel_key",
                "subtype",
                "site_anchor",
                "disease_level_1",
                "n_samples",
                "score_mean",
                "score_median",
                "score_sd",
                "healthy_baseline_mean",
                "delta_vs_healthy",
            ]
        )

    summary = (
        matched.groupby(["panel_key", "subtype", "site_anchor", "disease_level_1"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            score_mean=("score_mean", "mean"),
            score_median=("score_mean", "median"),
            score_sd=("score_mean", lambda s: float(np.std(np.asarray(s, dtype=float), ddof=0))),
        )
        .reset_index()
    )

    healthy_baseline = summary[summary["disease_level_1"].astype(str) == healthy_reference_label][["panel_key", "score_mean"]].copy()
    healthy_baseline = healthy_baseline.rename(columns={"score_mean": "healthy_baseline_mean"})
    summary = summary.merge(healthy_baseline, on="panel_key", how="left")
    summary["delta_vs_healthy"] = summary["score_mean"] - summary["healthy_baseline_mean"]
    return summary.sort_values(["subtype", "site_anchor", "disease_level_1"]).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    run_root = resolve_run_root(cfg, args.output_dir)
    val_dir = ensure_dir(run_root / "validation")

    block = deep_get(cfg, CONFIG_KEY, default={}) or {}
    query_path = Path(str(deep_get(block, "disease_query_h5ad", default="")))
    label_path = run_root / "labels" / "query_label_transfer.tsv"
    panel_registry_path = run_root / "panels" / "panel_registry.tsv"

    if not label_path.exists():
        raise FileNotFoundError(f"Missing label transfer file: {label_path}")
    if not panel_registry_path.exists():
        raise FileNotFoundError(f"Missing panel registry file: {panel_registry_path}")

    label_df = pd.read_csv(label_path, sep="\t")
    registry = pd.read_csv(panel_registry_path, sep="\t")

    panel_sets = build_panel_gene_sets(registry, cfg)

    qa = ad.read_h5ad(query_path, backed="r")
    try:
        var_names = qa.var_names.astype(str)
        gene_coverage = build_gene_coverage(panel_sets, var_names=var_names)
        gene_coverage_path = val_dir / "panel_gene_coverage.tsv"
        gene_coverage.to_csv(gene_coverage_path, sep="\t", index=False)

        usable_panels = gene_coverage[gene_coverage["n_present_genes"] > 0].copy()
        if usable_panels.empty:
            raise RuntimeError("No panel genes overlap with the disease query var_names.")

        cell_scores = compute_panel_scores(qa, label_df=label_df, gene_coverage=usable_panels, cfg=cfg)
    finally:
        qa.file.close()

    cell_scores_path = val_dir / "cell_scores.tsv"
    cell_scores.to_csv(cell_scores_path, sep="\t", index=False)

    sample_scores = summarize_sample_scores(cell_scores, cfg)
    sample_scores_path = val_dir / "sample_scores.tsv"
    sample_scores.to_csv(sample_scores_path, sep="\t", index=False)

    retention = summarize_retention(sample_scores, cfg)
    retention_path = val_dir / "retention_summary.tsv"
    retention.to_csv(retention_path, sep="\t", index=False)

    disease_shift = summarize_disease_shift(sample_scores, cfg)
    disease_shift_path = val_dir / "disease_shift_summary.tsv"
    disease_shift.to_csv(disease_shift_path, sep="\t", index=False)

    manifest = {
        "config_path": str(args.config),
        "run_root": str(run_root),
        "query_path": str(query_path),
        "files": {
            "panel_gene_coverage": str(gene_coverage_path),
            "cell_scores": str(cell_scores_path),
            "sample_scores": str(sample_scores_path),
            "retention_summary": str(retention_path),
            "disease_shift_summary": str(disease_shift_path),
        },
        "counts": {
            "n_usable_panels": int((usable_panels["n_present_genes"] >= int(deep_get(cfg, CONFIG_KEY, "validation", "min_panel_genes", default=5))).sum()),
            "n_cell_score_rows": int(cell_scores.shape[0]),
            "n_sample_score_rows": int(sample_scores.shape[0]),
            "n_retention_rows": int(retention.shape[0]),
            "n_disease_shift_rows": int(disease_shift.shape[0]),
        },
    }
    write_json(manifest, val_dir / "validation_manifest.json")

    print(json.dumps(manifest["counts"], ensure_ascii=False, indent=2))
    print(f"[ok] sample scores     -> {sample_scores_path}")
    print(f"[ok] retention summary -> {retention_path}")
    print(f"[ok] disease shift     -> {disease_shift_path}")


if __name__ == "__main__":
    main()
