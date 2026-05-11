#!/usr/bin/env python3
"""
Harmony2 helper for lineage-level single-cell integration.

This helper intentionally works from existing PCA embeddings when available so it
can run on large lineage h5ad files without loading the full count matrix. It
writes a lightweight h5ad containing Harmony2 embeddings, UMAP coordinates,
Leiden clusters, and summary tables/figures for downstream review.
"""

from __future__ import annotations

import json
import math
import time
import warnings
from pathlib import Path
from typing import Any, Iterable, Optional

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from scipy import sparse
from sklearn.metrics import silhouette_score

warnings.filterwarnings("ignore")

DEFAULT_BATCH_CANDIDATES = ["dataset", "sample", "batch", "orig.ident", "donorID", "donor_id"]
DEFAULT_CELLTYPE_CANDIDATES = ["cell_type_L3", "cell_type", "CellType", "cellType", "ann_level_3"]
DEFAULT_RESOLUTIONS = [0.4, 0.8, 1.2, 1.6, 2.0]


def _safe_mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if math.isfinite(float(value)):
            return float(value)
        return None
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(_json_safe(payload), ensure_ascii=False, indent=2), encoding="utf-8")


def _write_markdown(lines: Iterable[str], path: Path) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _resolve_obs_column(obs: pd.DataFrame, preferred: Optional[str], candidates: list[str], required: bool) -> Optional[str]:
    if preferred and preferred in obs.columns:
        return preferred
    for col in candidates:
        if col in obs.columns:
            return col
    if required:
        raise KeyError(f"None of the candidate columns were found: {candidates}")
    return None


def _clean_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in out.columns:
        if pd.api.types.is_categorical_dtype(out[col]):
            out[col] = out[col].astype(str)
        elif out[col].dtype == object:
            out[col] = out[col].astype(str)
    return out


def _load_existing_basis(h5ad_path: Path, basis_key: str, n_pcs: int) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    """Load obs and an existing PCA-like embedding without materializing X."""
    backed = sc.read_h5ad(h5ad_path, backed="r")
    obs = _clean_obs(backed.obs)
    obsm_keys = list(backed.obsm.keys())
    if basis_key not in obsm_keys:
        raise KeyError(f"basis_key={basis_key!r} not found in obsm. Available: {obsm_keys}")
    basis = np.asarray(backed.obsm[basis_key])
    if basis.ndim != 2 or basis.shape[0] != obs.shape[0]:
        raise ValueError(f"Invalid basis shape for {basis_key}: {basis.shape}; expected ({obs.shape[0]}, n_pcs)")
    use_pcs = min(int(n_pcs), basis.shape[1])
    if use_pcs < 2:
        raise ValueError(f"Need at least 2 PCs for Harmony2; got {use_pcs}")
    basis = np.asarray(basis[:, :use_pcs], dtype=np.float32)
    meta = {
        "input_shape": [int(backed.n_obs), int(backed.n_vars)],
        "input_layers": list(backed.layers.keys()),
        "input_obsm": obsm_keys,
        "basis_key": basis_key,
        "basis_shape": [int(basis.shape[0]), int(basis.shape[1])],
    }
    backed.file.close()
    return obs, basis, meta


def _compute_basis_from_counts(
    h5ad_path: Path,
    batch_col: str,
    n_pcs: int,
    n_top_genes: int,
    seed: int,
) -> tuple[pd.DataFrame, np.ndarray, dict[str, Any]]:
    """Fallback path when no PCA embedding exists. This may load X/layers into memory."""
    adata = sc.read_h5ad(h5ad_path)
    obs = _clean_obs(adata.obs)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    elif adata.raw is not None and adata.raw.X is not None:
        adata.X = adata.raw.X.copy()
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if adata.n_vars > n_top_genes:
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_col, flavor="seurat", subset=True)
        except Exception:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat", subset=True)
    n_comps = min(int(n_pcs), max(2, adata.n_vars - 1), max(2, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=seed)
    basis = np.asarray(adata.obsm["X_pca"][:, :n_comps], dtype=np.float32)
    meta = {
        "input_shape": [int(adata.n_obs), int(adata.n_vars)],
        "input_layers": list(adata.layers.keys()),
        "input_obsm": list(adata.obsm.keys()),
        "basis_key": "computed:X_pca",
        "basis_shape": [int(basis.shape[0]), int(basis.shape[1])],
        "computed_basis_from_counts": True,
        "n_top_genes": int(n_top_genes),
    }
    return obs, basis, meta


def _infer_harmony_nclust(n_cells: int) -> int:
    """Mirror harmonypy's default nclust while keeping tiny smoke tests valid."""
    return max(1, int(min(round(int(n_cells) / 30.0), 100)))


def _run_harmony2(
    obs: pd.DataFrame,
    basis: np.ndarray,
    batch_col: str,
    theta: float,
    lamb: float,
    sigma: float,
    tau: float,
    max_iter_harmony: int,
) -> ad.AnnData:
    work = ad.AnnData(X=sparse.csr_matrix((basis.shape[0], 0)), obs=obs.copy())
    work.obsm["X_pca_harmony2_input"] = basis
    nclust = _infer_harmony_nclust(basis.shape[0])
    sigma_vec = np.repeat(float(sigma), nclust).astype(float)
    sce.pp.harmony_integrate(
        work,
        key=batch_col,
        basis="X_pca_harmony2_input",
        adjusted_basis="X_harmony2",
        theta=float(theta),
        lamb=float(lamb),
        sigma=sigma_vec,
        nclust=nclust,
        tau=float(tau),
        max_iter_harmony=int(max_iter_harmony),
    )
    work.obsm["X_harmony2"] = np.asarray(work.obsm["X_harmony2"], dtype=np.float32)
    return work


def _run_neighbors_umap_clusters(
    work: ad.AnnData,
    resolutions: list[float],
    default_resolution: float,
    n_neighbors: int,
    min_dist: float,
    seed: int,
) -> tuple[str, list[str]]:
    sc.pp.neighbors(work, n_neighbors=int(n_neighbors), use_rep="X_harmony2", random_state=seed)
    sc.tl.umap(work, min_dist=float(min_dist), random_state=seed)
    work.obsm["X_umap_harmony2"] = np.asarray(work.obsm["X_umap"], dtype=np.float32)

    cluster_keys: list[str] = []
    for res in resolutions:
        key = f"leiden_harmony2_res{res:g}"
        try:
            sc.tl.leiden(work, resolution=float(res), key_added=key, flavor="igraph", n_iterations=2, directed=False)
        except TypeError:
            sc.tl.leiden(work, resolution=float(res), key_added=key, n_iterations=2)
        cluster_keys.append(key)
    default_key = f"leiden_harmony2_res{default_resolution:g}"
    if default_key not in work.obs.columns:
        default_key = cluster_keys[min(len(cluster_keys) - 1, len(cluster_keys) // 2)]
    work.obs["leiden_harmony2"] = work.obs[default_key].astype(str)
    return default_key, cluster_keys


def _safe_silhouette(embedding: np.ndarray, labels: pd.Series, sample_size: int, seed: int) -> Optional[float]:
    labels_str = labels.astype(str)
    if labels_str.nunique(dropna=True) < 2:
        return None
    n = embedding.shape[0]
    if n < 3:
        return None
    try:
        return float(
            silhouette_score(
                embedding,
                labels_str.to_numpy(),
                sample_size=min(int(sample_size), n),
                random_state=seed,
                metric="euclidean",
            )
        )
    except Exception:
        return None


def _write_summary_tables(
    work: ad.AnnData,
    output_dir: Path,
    batch_col: str,
    celltype_col: Optional[str],
    default_cluster_key: str,
    seed: int,
    metric_sample_size: int,
) -> dict[str, Any]:
    obs = work.obs.copy()
    obs[batch_col] = obs[batch_col].astype(str)
    obs["leiden_harmony2"] = obs["leiden_harmony2"].astype(str)

    cluster_counts = obs["leiden_harmony2"].value_counts().rename_axis("cluster").reset_index(name="n_cells")
    cluster_counts.to_csv(output_dir / "harmony2_cluster_counts.tsv", sep="\t", index=False)

    batch_by_cluster = pd.crosstab(obs["leiden_harmony2"], obs[batch_col])
    batch_by_cluster.index.name = "cluster"
    batch_by_cluster.to_csv(output_dir / "harmony2_batch_by_cluster_counts.tsv", sep="\t")
    batch_frac = batch_by_cluster.div(batch_by_cluster.sum(axis=1).replace(0, np.nan), axis=0)
    batch_frac.to_csv(output_dir / "harmony2_batch_by_cluster_fraction.tsv", sep="\t")

    table_payload: dict[str, Any] = {
        "cluster_counts_tsv": str(output_dir / "harmony2_cluster_counts.tsv"),
        "batch_by_cluster_counts_tsv": str(output_dir / "harmony2_batch_by_cluster_counts.tsv"),
        "batch_by_cluster_fraction_tsv": str(output_dir / "harmony2_batch_by_cluster_fraction.tsv"),
        "n_clusters_default": int(obs["leiden_harmony2"].nunique()),
        "default_cluster_key": default_cluster_key,
    }

    if celltype_col and celltype_col in obs.columns:
        obs[celltype_col] = obs[celltype_col].astype(str)
        celltype_by_cluster = pd.crosstab(obs["leiden_harmony2"], obs[celltype_col])
        celltype_by_cluster.index.name = "cluster"
        celltype_by_cluster.to_csv(output_dir / "harmony2_celltype_by_cluster_counts.tsv", sep="\t")
        table_payload["celltype_by_cluster_counts_tsv"] = str(output_dir / "harmony2_celltype_by_cluster_counts.tsv")

    umap = np.asarray(work.obsm["X_umap_harmony2"])
    umap_df = pd.DataFrame({
        "cell_id": obs.index.astype(str),
        "UMAP1": umap[:, 0],
        "UMAP2": umap[:, 1],
        "batch": obs[batch_col].to_numpy(),
        "leiden_harmony2": obs["leiden_harmony2"].to_numpy(),
    })
    if celltype_col and celltype_col in obs.columns:
        umap_df[celltype_col] = obs[celltype_col].to_numpy()
    umap_df.to_csv(output_dir / "harmony2_umap_coordinates.tsv", sep="\t", index=False)
    table_payload["umap_coordinates_tsv"] = str(output_dir / "harmony2_umap_coordinates.tsv")

    metrics = {
        "silhouette_batch_harmony2": _safe_silhouette(work.obsm["X_harmony2"], obs[batch_col], metric_sample_size, seed),
        "silhouette_cluster_harmony2": _safe_silhouette(work.obsm["X_harmony2"], obs["leiden_harmony2"], metric_sample_size, seed),
    }
    if celltype_col and celltype_col in obs.columns:
        metrics["silhouette_celltype_harmony2"] = _safe_silhouette(work.obsm["X_harmony2"], obs[celltype_col], metric_sample_size, seed)
    pd.DataFrame([metrics]).to_csv(output_dir / "harmony2_metrics.tsv", sep="\t", index=False)
    table_payload["metrics_tsv"] = str(output_dir / "harmony2_metrics.tsv")
    table_payload["metrics"] = metrics
    return table_payload


def _plot_umap_panel(
    work: ad.AnnData,
    output_dir: Path,
    batch_col: str,
    celltype_col: Optional[str],
    max_plot_cells: int,
    seed: int,
    dpi: int,
) -> Optional[str]:
    if "X_umap_harmony2" not in work.obsm:
        return None
    obs = work.obs
    n = work.n_obs
    rng = np.random.default_rng(seed)
    if n > max_plot_cells:
        idx = np.sort(rng.choice(n, size=int(max_plot_cells), replace=False))
    else:
        idx = np.arange(n)
    umap = np.asarray(work.obsm["X_umap_harmony2"])[idx]

    color_cols = [batch_col, "leiden_harmony2"]
    if celltype_col and celltype_col in obs.columns and obs[celltype_col].astype(str).nunique() <= 50:
        color_cols.append(celltype_col)
    ncols = len(color_cols)
    fig, axes = plt.subplots(1, ncols, figsize=(5.2 * ncols, 4.8), squeeze=False)
    for ax, col in zip(axes[0], color_cols):
        labels = obs.iloc[idx][col].astype(str)
        codes, uniques = pd.factorize(labels, sort=True)
        ax.scatter(umap[:, 0], umap[:, 1], c=codes, cmap="tab20", s=1.0, linewidths=0, alpha=0.75)
        ax.set_title(col)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_xticks([])
        ax.set_yticks([])
        if len(uniques) <= 16:
            handles = [plt.Line2D([0], [0], marker="o", color="w", label=str(lbl), markerfacecolor=plt.cm.tab20(i % 20), markersize=5) for i, lbl in enumerate(uniques)]
            ax.legend(handles=handles, loc="best", fontsize=6, frameon=False, markerscale=1)
    fig.suptitle(f"Harmony2 UMAP (plotted {len(idx):,}/{n:,} cells)", y=1.02)
    fig.tight_layout()
    out = output_dir / "harmony2_umap_overview.png"
    fig.savefig(out, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)
    return str(out)


def run_harmony2_full(
    adata_h5ad_path: str | Path,
    output_dir: str | Path,
    run_name: str = "harmony2",
    batch_col: Optional[str] = None,
    celltype_col: Optional[str] = None,
    basis_key: str = "X_pca",
    n_pcs: int = 50,
    n_top_genes: int = 4000,
    theta: float = 1.0,
    lamb: float = 6.0,
    sigma: float = 0.1,
    tau: float = 0.0,
    max_iter_harmony: int = 20,
    n_neighbors: int = 15,
    resolutions: Optional[list[float]] = None,
    default_resolution: float = 1.2,
    umap_min_dist: float = 0.5,
    seed: int = 42,
    metric_sample_size: int = 10000,
    max_plot_cells: int = 80000,
    dpi: int = 300,
    force_recompute_pca: bool = False,
) -> dict[str, Any]:
    started = time.time()
    h5ad_path = Path(adata_h5ad_path)
    output_dir = _safe_mkdir(Path(output_dir))
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Input h5ad does not exist: {h5ad_path}")

    resolutions = [float(x) for x in (resolutions or DEFAULT_RESOLUTIONS)]
    if len(resolutions) == 0:
        resolutions = DEFAULT_RESOLUTIONS

    try:
        if force_recompute_pca:
            raise KeyError("force_recompute_pca=True")
        obs, basis, meta = _load_existing_basis(h5ad_path, basis_key=basis_key, n_pcs=int(n_pcs))
    except Exception as basis_error:
        # Re-open lightweight obs first to resolve the batch key before full fallback.
        preview = sc.read_h5ad(h5ad_path, backed="r")
        preview_obs = _clean_obs(preview.obs)
        resolved_batch_col = _resolve_obs_column(preview_obs, batch_col, DEFAULT_BATCH_CANDIDATES, required=True)
        preview.file.close()
        obs, basis, meta = _compute_basis_from_counts(
            h5ad_path,
            batch_col=resolved_batch_col,
            n_pcs=int(n_pcs),
            n_top_genes=int(n_top_genes),
            seed=int(seed),
        )
        meta["basis_fallback_reason"] = str(basis_error)

    resolved_batch_col = _resolve_obs_column(obs, batch_col, DEFAULT_BATCH_CANDIDATES, required=True)
    resolved_celltype_col = _resolve_obs_column(obs, celltype_col, DEFAULT_CELLTYPE_CANDIDATES, required=False)

    n_batches = obs[resolved_batch_col].astype(str).nunique(dropna=True)
    if n_batches < 2:
        raise ValueError(f"Harmony2 requires at least 2 batches in {resolved_batch_col}; got {n_batches}")

    work = _run_harmony2(
        obs=obs,
        basis=basis,
        batch_col=resolved_batch_col,
        theta=float(theta),
        lamb=float(lamb),
        sigma=float(sigma),
        tau=float(tau),
        max_iter_harmony=int(max_iter_harmony),
    )
    default_cluster_key, cluster_keys = _run_neighbors_umap_clusters(
        work,
        resolutions=resolutions,
        default_resolution=float(default_resolution),
        n_neighbors=int(n_neighbors),
        min_dist=float(umap_min_dist),
        seed=int(seed),
    )

    table_payload = _write_summary_tables(
        work=work,
        output_dir=output_dir,
        batch_col=resolved_batch_col,
        celltype_col=resolved_celltype_col,
        default_cluster_key=default_cluster_key,
        seed=int(seed),
        metric_sample_size=int(metric_sample_size),
    )
    figure_path = _plot_umap_panel(
        work=work,
        output_dir=output_dir,
        batch_col=resolved_batch_col,
        celltype_col=resolved_celltype_col,
        max_plot_cells=int(max_plot_cells),
        seed=int(seed),
        dpi=int(dpi),
    )

    work.uns["harmony2_config"] = {
        "run_name": run_name,
        "batch_col": resolved_batch_col,
        "celltype_col": resolved_celltype_col,
        "basis_key": meta.get("basis_key", basis_key),
        "n_pcs": int(basis.shape[1]),
        "theta": float(theta),
        "lambda": float(lamb),
        "sigma": float(sigma),
        "tau": float(tau),
        "max_iter_harmony": int(max_iter_harmony),
        "n_neighbors": int(n_neighbors),
        "resolutions": resolutions,
        "default_resolution": float(default_resolution),
        "default_cluster_key": default_cluster_key,
        "seed": int(seed),
    }
    result_h5ad = output_dir / "harmony2_result.h5ad"
    work.write_h5ad(result_h5ad, compression="gzip")

    elapsed_min = (time.time() - started) / 60.0
    result = {
        "success": True,
        "status": "ok",
        "run_name": run_name,
        "input_h5ad": str(h5ad_path),
        "output_dir": str(output_dir),
        "result_h5ad": str(result_h5ad),
        "batch_col": resolved_batch_col,
        "celltype_col": resolved_celltype_col,
        "n_cells": int(work.n_obs),
        "n_pcs": int(work.obsm["X_harmony2"].shape[1]),
        "n_batches": int(n_batches),
        "cluster_keys": cluster_keys,
        "default_cluster_key": default_cluster_key,
        "elapsed_min": elapsed_min,
        "input_meta": meta,
        "tables": table_payload,
        "figure": figure_path,
    }
    _write_json(result, output_dir / "harmony2_summary.json")
    _write_markdown(
        [
            f"# Harmony2 run: {run_name}",
            "",
            f"- Input: `{h5ad_path}`",
            f"- Cells: {work.n_obs:,}",
            f"- Batch column: `{resolved_batch_col}` ({n_batches} batches)",
            f"- Cell type column: `{resolved_celltype_col}`" if resolved_celltype_col else "- Cell type column: not found",
            f"- Corrected embedding: `X_harmony2` ({work.obsm['X_harmony2'].shape[1]} PCs)",
            f"- UMAP embedding: `X_umap_harmony2`",
            f"- Default cluster key: `{default_cluster_key}` copied to `leiden_harmony2`",
            f"- Output h5ad: `{result_h5ad}`",
            f"- Elapsed: {elapsed_min:.2f} min",
            "",
            "## Metrics",
            *(f"- {k}: {v}" for k, v in table_payload.get("metrics", {}).items()),
            "",
            "## Key tables",
            *(f"- {k}: `{v}`" for k, v in table_payload.items() if k.endswith("_tsv")),
        ],
        output_dir / "harmony2_summary.md",
    )
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run Harmony2 integration from a h5ad file")
    parser.add_argument("config_json", help="JSON config written by the R worker wrapper")
    args = parser.parse_args()
    cfg = json.loads(Path(args.config_json).read_text(encoding="utf-8"))
    res = run_harmony2_full(**cfg)
    print(json.dumps(_json_safe(res), ensure_ascii=False, indent=2))
