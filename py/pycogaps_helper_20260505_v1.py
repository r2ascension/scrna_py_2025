#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Formal PyCoGAPS helper with shared gene-exclusion sidecars."""

from __future__ import annotations

import json
import logging
import importlib.util
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

_GENE_EXCLUSION_HELPER_PATH = Path(__file__).with_name("program_gene_exclusion_helper_20260505_v1.py")
_gene_exclusion_spec = importlib.util.spec_from_file_location(
    "program_gene_exclusion_helper_20260505_v1",
    _GENE_EXCLUSION_HELPER_PATH,
)
if _gene_exclusion_spec is None or _gene_exclusion_spec.loader is None:
    raise ImportError(f"Cannot load gene exclusion helper: {_GENE_EXCLUSION_HELPER_PATH}")
_gene_exclusion_module = importlib.util.module_from_spec(_gene_exclusion_spec)
sys.modules[_gene_exclusion_spec.name] = _gene_exclusion_module
_gene_exclusion_spec.loader.exec_module(_gene_exclusion_module)
apply_gene_exclusion_to_adata = _gene_exclusion_module.apply_gene_exclusion_to_adata
safe_mkdir = _gene_exclusion_module.safe_mkdir

try:
    import pycogaps  # noqa: F401
    from PyCoGAPS.parameters import CoParams, setParams
    from PyCoGAPS.pycogaps_main import CoGAPS
    PYCOGAPS_AVAILABLE = True
except ImportError:
    PYCOGAPS_AVAILABLE = False
    CoParams = None
    setParams = None
    CoGAPS = None

LOGGER = logging.getLogger("pycogaps_helper")
if not LOGGER.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    LOGGER.addHandler(handler)
LOGGER.setLevel(logging.INFO)
LOGGER.propagate = False

PYCOGAPS_HELPER_VERSION_20260505_V1 = "20260505_v1"


def plot_pycogaps_outputs(
    gene_patterns: pd.DataFrame,
    cell_patterns: pd.DataFrame,
    output_dir: Path,
    top_n_genes: int = 25,
) -> Dict[str, str]:
    """Create lightweight PyCoGAPS method visualizations from exported pattern matrices."""
    output_dir = Path(output_dir)
    viz: Dict[str, str] = {}
    pattern_cols = [c for c in gene_patterns.columns if str(c).startswith("Pattern")]
    if not pattern_cols:
        return viz
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional plotting dependency
        LOGGER.warning("PyCoGAPS plotting skipped: %s", exc)
        return viz

    top_genes: list[str] = []
    for pat in pattern_cols:
        top_genes.extend(gene_patterns[pat].sort_values(ascending=False).head(top_n_genes).index.astype(str).tolist())
    top_genes = list(dict.fromkeys(top_genes))
    if top_genes:
        mat = gene_patterns.loc[top_genes, pattern_cols].to_numpy(dtype=float)
        fig_h = max(5.0, min(16.0, 0.22 * len(top_genes) + 2.0))
        fig_w = max(6.0, 0.65 * len(pattern_cols) + 3.0)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(mat, aspect="auto", cmap="viridis")
        ax.set_xticks(np.arange(len(pattern_cols)))
        ax.set_xticklabels(pattern_cols, rotation=45, ha="right")
        ax.set_yticks(np.arange(len(top_genes)))
        ax.set_yticklabels(top_genes, fontsize=7)
        ax.set_title(f"PyCoGAPS top gene pattern weights (top {top_n_genes}/pattern)")
        ax.set_xlabel("Pattern")
        ax.set_ylabel("Gene")
        fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Weight")
        fig.tight_layout()
        out = output_dir / "pycogaps_gene_pattern_heatmap.png"
        fig.savefig(out, dpi=180)
        plt.close(fig)
        viz["gene_pattern_heatmap_png"] = str(out)

    if not cell_patterns.empty:
        values = [cell_patterns[pat].dropna().to_numpy(dtype=float) for pat in pattern_cols if pat in cell_patterns]
        labels = [pat for pat in pattern_cols if pat in cell_patterns]
        if values:
            fig, ax = plt.subplots(figsize=(max(6.0, 0.7 * len(labels) + 2.0), 5.0))
            ax.boxplot(values, labels=labels, showfliers=False)
            ax.set_title("PyCoGAPS cell pattern score distributions")
            ax.set_ylabel("Pattern score")
            ax.tick_params(axis="x", rotation=45)
            fig.tight_layout()
            out = output_dir / "pycogaps_cell_pattern_score_boxplot.png"
            fig.savefig(out, dpi=180)
            plt.close(fig)
            viz["cell_pattern_score_boxplot_png"] = str(out)

    return viz


def prepare_pycogaps_input(
    adata: sc.AnnData,
    output_dir: Path,
    layer: str = "counts",
    gene_exclusion_config: Optional[Dict[str, Any]] = None,
    celltype_col: Optional[str] = None,
) -> Dict[str, Any]:
    if layer not in adata.layers and layer != "X":
        raise ValueError(f"Layer '{layer}' not present in adata.layers")

    filtered, packet, manifest = apply_gene_exclusion_to_adata(
        adata,
        output_dir=Path(output_dir),
        gene_exclusion_config=gene_exclusion_config,
        prefix="gene_exclusion",
        celltype_col=celltype_col,
    )

    mat = filtered.layers[layer] if layer in filtered.layers else filtered.X
    dense = mat.T.toarray().astype(np.float32, copy=False) if sp.issparse(mat) else np.asarray(mat.T, dtype=np.float32)
    cogaps_adata = ad.AnnData(X=dense)
    cogaps_adata.obs_names = filtered.var_names.copy()
    cogaps_adata.var_names = filtered.obs_names.copy()

    return {
        "filtered_adata": filtered,
        "cogaps_adata": cogaps_adata,
        "gene_exclusion_packet": packet,
        "gene_exclusion_manifest": manifest,
    }


def export_pycogaps_outputs(
    result: ad.AnnData,
    output_dir: Path,
    n_iterations: int,
    run_name: str,
    manifest: Dict[str, Any],
) -> Dict[str, Any]:
    output_dir = Path(output_dir)
    gene_pattern_path = output_dir / "pycogaps_gene_patterns.tsv"
    cell_pattern_path = output_dir / "pycogaps_cell_patterns.tsv"
    top_genes_path = output_dir / "pycogaps_top_genes_by_pattern.json"
    summary_path = output_dir / "pycogaps_summary.json"

    pattern_cols = [c for c in result.obs.columns if str(c).startswith("Pattern")]
    gene_patterns = result.obs[pattern_cols].copy()
    cell_patterns = result.var[pattern_cols].copy()
    gene_patterns.to_csv(gene_pattern_path, sep="\t")
    cell_patterns.to_csv(cell_pattern_path, sep="\t")

    top_genes = {
        pat: gene_patterns[pat].sort_values(ascending=False).head(30).index.tolist()
        for pat in pattern_cols
    }
    top_genes_path.write_text(json.dumps(top_genes, indent=2, ensure_ascii=False), encoding="utf-8")
    visualizations = plot_pycogaps_outputs(gene_patterns, cell_patterns, output_dir)

    summary = {
        "success": True,
        "run_name": run_name,
        "n_genes": int(result.n_obs),
        "n_cells": int(result.n_vars),
        "n_patterns": len(pattern_cols),
        "n_iterations": int(n_iterations),
        "pattern_columns": pattern_cols,
        "gene_pattern_tsv": str(gene_pattern_path),
        "cell_pattern_tsv": str(cell_pattern_path),
        "top_genes_json": str(top_genes_path),
        "visualizations": visualizations,
        "gene_exclusion": manifest,
    }
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def run_pycogaps_full(
    adata: sc.AnnData,
    output_dir: Path,
    run_name: str,
    n_patterns: int = 6,
    n_iterations: int = 50,
    seed: int = 42,
    n_threads: int = 1,
    layer: str = "counts",
    celltype_col: Optional[str] = None,
    gene_exclusion_config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    if not PYCOGAPS_AVAILABLE:
        return {"success": False, "error": "PyCoGAPS not installed"}

    output_dir = Path(output_dir)
    if not safe_mkdir(output_dir):
        return {"success": False, "error": f"Cannot create output directory: {output_dir}"}

    prepared = prepare_pycogaps_input(
        adata=adata,
        output_dir=output_dir,
        layer=layer,
        gene_exclusion_config=gene_exclusion_config,
        celltype_col=celltype_col,
    )
    cogaps_adata = prepared["cogaps_adata"]
    manifest = prepared["gene_exclusion_manifest"]

    params = CoParams(adata=cogaps_adata)
    setParams(
        params,
        {
            "nIterations": int(n_iterations),
            "seed": int(seed),
            "nPatterns": int(n_patterns),
            "useSparseOptimization": False,
        },
    )
    LOGGER.info(
        "Running PyCoGAPS on %d genes x %d cells after exclusion",
        cogaps_adata.n_obs,
        cogaps_adata.n_vars,
    )
    result = CoGAPS(cogaps_adata, params=params, nThreads=int(n_threads), messages=False, outputFrequency=max(10, int(n_iterations // 2)))
    summary = export_pycogaps_outputs(
        result=result,
        output_dir=output_dir,
        n_iterations=n_iterations,
        run_name=run_name,
        manifest=manifest,
    )
    return summary


if __name__ == "__main__":
    print("[INFO] pycogaps_helper_20260505_v1.py loaded")
    print(f"  PyCoGAPS available: {PYCOGAPS_AVAILABLE}")
