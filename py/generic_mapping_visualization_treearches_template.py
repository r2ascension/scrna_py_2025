#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generic Mapping Visualization Template (treeArches-ready)

Purpose
-------
Provide a reusable visualization scaffold for scArches / scANVI / treeArches
mapping outputs while keeping:
1. project-specific parameters in one config object
2. reusable logic in helper functions
3. non-AnnData outputs in a dedicated artifacts registry

Usage
-----
- Update the paths in `get_config()`.
- Adjust `keys` and `marker_sets` for the current lineage.
- Run the script to generate overview figures, optional treeArches diagnostics,
  and summary exports.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, Optional
import json
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from pandas.api.types import CategoricalDtype

warnings.filterwarnings("ignore")


@dataclass
class MappingVizConfig:
    query_h5ad: str
    merged_h5ad: str
    output_dir: Path
    version: str = "v1.0-template"
    pipeline_name: str = "generic_mapping_visualization"
    run_mode: str = "predict"  # predict | update
    label_level: str = "L2"
    dpi: int = 300
    figure_format: str = "pdf"
    keys: Dict[str, str] = field(default_factory=dict)
    marker_sets: Dict[str, list[str]] = field(default_factory=dict)
    palette_source: Dict[str, str] = field(
        default_factory=lambda: {"reference": "#1f77b4", "query": "#ff7f0e"}
    )
    cmap_confidence: str = "viridis"
    cmap_expression: str = "Reds"
    cmap_novelty: str = "magma"

    @property
    def split_prefix(self) -> str:
        return self.label_level.lower()


def get_config() -> MappingVizConfig:
    """Update this function per project / lineage."""
    return MappingVizConfig(
        query_h5ad="/absolute/path/to/query.h5ad",
        merged_h5ad="/absolute/path/to/reference_plus_query_merged.h5ad",
        output_dir=Path("/absolute/path/to/output_dir"),
        version="v1.0-template",
        run_mode="predict",
        label_level="L2",
        keys={
            "ref_label": "Cell_Type_L2",
            "pred_label": "Cell_Type_L2_pred",
            "final_label": "Cell_Type_L2_final",
            "confidence": "mapping_confidence",
            "datasource": "data_source",
            "batch": "sample",
            "tissue": "tissue",
            "hier_pred": "tree_pred",
            "hier_updated": "tree_updated_label",
            "rejected": "tree_rejected",
            "rejection_reason": "tree_rejection_reason",
            "novelty_score": "tree_novelty_score",
            "novelty_flag": "tree_novelty_flag",
        },
        marker_sets={
            "default": ["CD19", "MS4A1", "CD27", "IGHD", "IGHM", "MZB1", "SDC1", "JCHAIN"],
        },
    )


def setup_plotting(cfg: MappingVizConfig) -> None:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    sc.settings.set_figure_params(
        dpi=cfg.dpi,
        facecolor="white",
        format=cfg.figure_format,
        vector_friendly=False,
    )


def find_first_available_key(adata: sc.AnnData, candidate_keys: Iterable[str], *, required: bool = True) -> Optional[str]:
    for key in candidate_keys:
        if key and key in adata.obs.columns:
            return key
    if required:
        raise KeyError(f"Missing all candidate keys: {list(candidate_keys)}")
    return None


def validate_obs_keys(adata: sc.AnnData, required_keys: Iterable[str], *, adata_name: str) -> None:
    missing = [k for k in required_keys if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"{adata_name} missing required obs columns: {missing}")


def ensure_umap_basis(
    adata: sc.AnnData,
    *,
    adata_name: str,
    preferred_keys: Iterable[str] = ("X_umap", "X_umap_scanvi", "X_umap_scvi"),
) -> str:
    for key in preferred_keys:
        if key in adata.obsm:
            if key != "X_umap":
                adata.obsm["X_umap"] = np.asarray(adata.obsm[key]).copy()
            return key
    raise KeyError(f"{adata_name} missing usable UMAP basis: {tuple(preferred_keys)}")


def build_global_categories(
    adata: sc.AnnData,
    *,
    ref_key: str,
    qry_key: str,
    datasource_key: str,
) -> list[str]:
    ref_mask = adata.obs[datasource_key].eq("reference")
    qry_mask = adata.obs[datasource_key].eq("query")
    ref_vals = adata.obs.loc[ref_mask, ref_key].astype("string")
    qry_vals = adata.obs.loc[qry_mask, qry_key].astype("string")
    return sorted(pd.Index(ref_vals.dropna().unique()).union(pd.Index(qry_vals.dropna().unique())).tolist())


def make_palette(categories: list[str]) -> dict[str, Any]:
    if len(categories) <= 20:
        colors = list(plt.get_cmap("tab20").colors)[: len(categories)]
    else:
        colors = sc.pl.palettes.default_102[: len(categories)]
    return dict(zip(categories, colors))


def prepare_split_obs_columns(
    adata: sc.AnnData,
    *,
    ref_key: str,
    qry_key: str,
    datasource_key: str,
    prefix: str,
    categories: Optional[list[str]] = None,
    confidence_key: Optional[str] = None,
) -> dict[str, Any]:
    ref_only_col = f"{prefix}_reference_only"
    qry_only_col = f"{prefix}_query_only"
    qry_conf_col = f"{prefix}_confidence_query_only"

    ref_mask = adata.obs[datasource_key].eq("reference")
    qry_mask = adata.obs[datasource_key].eq("query")

    adata.obs[ref_only_col] = pd.Series(pd.NA, index=adata.obs_names, dtype="string")
    adata.obs[qry_only_col] = pd.Series(pd.NA, index=adata.obs_names, dtype="string")

    if ref_key in adata.obs.columns:
        adata.obs.loc[ref_mask, ref_only_col] = adata.obs.loc[ref_mask, ref_key].astype("string")
    if qry_key in adata.obs.columns:
        adata.obs.loc[qry_mask, qry_only_col] = adata.obs.loc[qry_mask, qry_key].astype("string")

    if categories is not None:
        dtype = CategoricalDtype(categories=categories, ordered=False)
        adata.obs[ref_only_col] = adata.obs[ref_only_col].astype("string").astype(dtype)
        adata.obs[qry_only_col] = adata.obs[qry_only_col].astype("string").astype(dtype)

    conf_vals = np.full(adata.n_obs, np.nan, dtype=float)
    if confidence_key and confidence_key in adata.obs.columns:
        conf_vals[qry_mask] = pd.to_numeric(adata.obs.loc[qry_mask, confidence_key], errors="coerce").to_numpy()
    adata.obs[qry_conf_col] = conf_vals

    return {
        "ref_only_col": ref_only_col,
        "qry_only_col": qry_only_col,
        "qry_conf_col": qry_conf_col,
        "ref_mask": ref_mask,
        "qry_mask": qry_mask,
    }


def determine_expression_source(adata: sc.AnnData) -> dict[str, Any]:
    if adata.raw is not None and adata.raw.n_vars > adata.n_vars:
        return {"use_raw": True, "layer": None}
    if adata.layers and "log1p" in adata.layers:
        return {"use_raw": False, "layer": "log1p"}
    if adata.layers and "counts" in adata.layers:
        return {"use_raw": False, "layer": "counts"}
    return {"use_raw": False, "layer": None}


def get_available_markers(adata: sc.AnnData, markers: list[str], source_cfg: dict[str, Any]) -> list[str]:
    if source_cfg.get("use_raw") and adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names
    return [gene for gene in markers if gene in var_names]


def save_rasterized_figure(fig: plt.Figure, path: Path, *, dpi: int) -> None:
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_tree_diagnostics(
    adata_query: sc.AnnData,
    adata_merged: sc.AnnData,
    *,
    cfg: MappingVizConfig,
    palette: dict[str, Any],
) -> list[str]:
    created: list[str] = []
    specs = [
        (cfg.keys["rejected"], "Rejected vs Accepted", {False: "#4daf4a", True: "#e41a1c"}),
        (cfg.keys["novelty_score"], "Novelty Score", cfg.cmap_novelty),
        (cfg.keys["hier_updated"], "Hierarchy-updated Label", palette),
    ]
    available = [spec for spec in specs if spec[0] in adata_query.obs.columns or spec[0] in adata_merged.obs.columns]
    if not available:
        return created

    fig, axes = plt.subplots(1, len(available), figsize=(7 * len(available), 6))
    axes = np.atleast_1d(axes)

    for ax, (obs_key, title, style) in zip(axes, available):
        adata_src = adata_query if obs_key in adata_query.obs.columns else adata_merged
        kwargs = {
            "adata": adata_src,
            "color": obs_key,
            "ax": ax,
            "show": False,
            "title": title,
            "frameon": False,
            "s": 25,
        }
        if obs_key == cfg.keys["novelty_score"]:
            kwargs["cmap"] = style
            kwargs["vmin"] = 0
            kwargs["vmax"] = 1
        else:
            kwargs["palette"] = style
            kwargs["na_color"] = "lightgray"
        sc.pl.umap(**kwargs)

    out = cfg.output_dir / f"treearches_ready_diagnostics.{cfg.figure_format}"
    save_rasterized_figure(fig, out, dpi=cfg.dpi)
    created.append(out.name)
    return created


def build_artifacts_registry(**kwargs: Any) -> dict[str, Any]:
    return {
        "adata_query": kwargs.get("adata_query"),
        "adata_ref": kwargs.get("adata_ref"),
        "adata_merged": kwargs.get("adata_merged"),
        "latent_key": kwargs.get("latent_key"),
        "tree": kwargs.get("tree"),
        "classifier": kwargs.get("classifier"),
        "prediction_table": kwargs.get("prediction_table"),
        "novelty_table": kwargs.get("novelty_table"),
        "metadata": kwargs.get("metadata", {}),
    }


def export_summary(
    *,
    cfg: MappingVizConfig,
    adata_query: sc.AnnData,
    adata_merged: sc.AnnData,
    artifacts: dict[str, Any],
    active_query_label_key: str,
    active_merged_query_label_key: str,
    active_reference_label_key: str,
    generated_files: list[str],
) -> None:
    summary_lines = [
        "=" * 80,
        f"{cfg.pipeline_name.upper()} SUMMARY ({cfg.version})",
        "=" * 80,
        f"Run mode: {cfg.run_mode}",
        f"Label level: {cfg.label_level}",
        f"Query cells: {adata_query.n_obs:,}",
        f"Merged cells: {adata_merged.n_obs:,}",
        f"Active query label key: {active_query_label_key}",
        f"Active merged query label key: {active_merged_query_label_key}",
        f"Active reference label key: {active_reference_label_key}",
        "",
        "Generated files:",
        *[f"  - {name}" for name in generated_files],
    ]

    report_path = cfg.output_dir / "visualization_summary_report_generic.txt"
    report_path.write_text("\n".join(summary_lines), encoding="utf-8")

    manifest = {
        "config": asdict(cfg),
        "active_keys": {
            "query_label": active_query_label_key,
            "merged_query_label": active_merged_query_label_key,
            "reference_label": active_reference_label_key,
        },
        "artifacts_metadata": artifacts.get("metadata", {}),
        "generated_files": generated_files,
    }
    manifest_path = cfg.output_dir / "visualization_manifest_generic.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")


def main() -> None:
    cfg = get_config()
    setup_plotting(cfg)

    run_timestamp = datetime.now().isoformat(timespec="seconds")
    print(f"[{run_timestamp}] Loading data...")
    adata_query = sc.read_h5ad(cfg.query_h5ad)
    adata_merged = sc.read_h5ad(cfg.merged_h5ad)

    query_umap_source = ensure_umap_basis(adata_query, adata_name="adata_query")
    merged_umap_source = ensure_umap_basis(adata_merged, adata_name="adata_merged")

    active_query_label_key = find_first_available_key(
        adata_query,
        [cfg.keys["hier_updated"], cfg.keys["final_label"], cfg.keys["pred_label"]],
    )
    active_merged_query_label_key = find_first_available_key(
        adata_merged,
        [cfg.keys["hier_updated"], cfg.keys["final_label"], cfg.keys["pred_label"]],
    )
    active_reference_label_key = find_first_available_key(adata_merged, [cfg.keys["ref_label"]])

    validate_obs_keys(adata_query, [active_query_label_key], adata_name="adata_query")
    validate_obs_keys(
        adata_merged,
        [cfg.keys["datasource"], active_reference_label_key, active_merged_query_label_key],
        adata_name="adata_merged",
    )

    global_categories = build_global_categories(
        adata_merged,
        ref_key=active_reference_label_key,
        qry_key=active_merged_query_label_key,
        datasource_key=cfg.keys["datasource"],
    )
    global_dtype = CategoricalDtype(categories=global_categories, ordered=False)
    palette = make_palette(global_categories)

    split_info = prepare_split_obs_columns(
        adata_merged,
        ref_key=active_reference_label_key,
        qry_key=active_merged_query_label_key,
        datasource_key=cfg.keys["datasource"],
        confidence_key=cfg.keys["confidence"],
        prefix=cfg.split_prefix,
        categories=global_categories,
    )

    for key_name in ["final_label", "pred_label", "hier_pred", "hier_updated"]:
        obs_key = cfg.keys[key_name]
        if obs_key in adata_query.obs.columns:
            adata_query.obs[obs_key] = adata_query.obs[obs_key].astype("string").astype(global_dtype)

    marker_cfg_query = determine_expression_source(adata_query)
    marker_cfg_merged = determine_expression_source(adata_merged)
    markers = cfg.marker_sets.get("default", [])
    available_markers = get_available_markers(adata_query, markers, marker_cfg_query)

    created_files: list[str] = []

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    sc.pl.umap(
        adata_merged,
        color=split_info["ref_only_col"],
        ax=axes[0],
        show=False,
        title=f"Reference {cfg.label_level}",
        legend_loc="right margin",
        palette=palette,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    sc.pl.umap(
        adata_merged,
        color=split_info["qry_only_col"],
        ax=axes[1],
        show=False,
        title=f"Query {cfg.label_level}",
        legend_loc="right margin",
        palette=palette,
        frameon=False,
        s=20,
        na_color="lightgray",
    )
    out = cfg.output_dir / f"reference_vs_query_generic.{cfg.figure_format}"
    save_rasterized_figure(fig, out, dpi=cfg.dpi)
    created_files.append(out.name)

    if available_markers:
        fig, ax = plt.subplots(figsize=(8, 6))
        marker_kwargs = {
            "adata": adata_query,
            "color": available_markers[0],
            "ax": ax,
            "show": False,
            "title": f"Marker: {available_markers[0]}",
            "cmap": cfg.cmap_expression,
            "frameon": False,
            "s": 30,
        }
        if marker_cfg_query["use_raw"]:
            marker_kwargs["use_raw"] = True
        elif marker_cfg_query["layer"] is not None:
            marker_kwargs["layer"] = marker_cfg_query["layer"]
        else:
            marker_kwargs["use_raw"] = False
        sc.pl.umap(**marker_kwargs)
        out = cfg.output_dir / f"marker_preview_generic.{cfg.figure_format}"
        save_rasterized_figure(fig, out, dpi=cfg.dpi)
        created_files.append(out.name)

    created_files.extend(
        plot_tree_diagnostics(
            adata_query,
            adata_merged,
            cfg=cfg,
            palette=palette,
        )
    )

    artifacts = build_artifacts_registry(
        adata_query=adata_query,
        adata_merged=adata_merged,
        metadata={
            "query_umap_source": query_umap_source,
            "merged_umap_source": merged_umap_source,
            "global_categories": global_categories,
            "available_markers": available_markers,
        },
    )

    export_summary(
        cfg=cfg,
        adata_query=adata_query,
        adata_merged=adata_merged,
        artifacts=artifacts,
        active_query_label_key=active_query_label_key,
        active_merged_query_label_key=active_merged_query_label_key,
        active_reference_label_key=active_reference_label_key,
        generated_files=created_files,
    )

    print("Done. Update get_config() with real project paths before production use.")


if __name__ == "__main__":
    main()
