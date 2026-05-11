#!/usr/bin/env python3
"""Supplemental first-tier trajectory methods for smc_anno.

Runs expression-sufficient Python trajectory/ordering methods when packages are
available. Each method is isolated so one failure does not stop the batch.
"""

from __future__ import annotations

import argparse
import json
import math
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import scipy.sparse as sp


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_json(obj, path: Path) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def dense_sample(X, max_cells: int = 2500):
    if sp.issparse(X):
        X = X.toarray()
    X = np.asarray(X)
    if X.shape[0] > max_cells:
        rng = np.random.default_rng(42)
        idx = np.sort(rng.choice(X.shape[0], max_cells, replace=False))
        return X[idx, :], idx
    return X, np.arange(X.shape[0])


def normalize01(x):
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    out = np.full(x.shape, np.nan, dtype=float)
    if not finite.any():
        return out
    lo, hi = np.nanmin(x[finite]), np.nanmax(x[finite])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        out[finite] = 0.0
    else:
        out[finite] = (x[finite] - lo) / (hi - lo)
    return out


def orient_root_low(pseudo, root_mask):
    pseudo = normalize01(pseudo)
    finite = np.isfinite(pseudo)
    if finite.any() and np.asarray(root_mask).any():
        root_mean = np.nanmean(pseudo[np.asarray(root_mask) & finite])
        other_mean = np.nanmean(pseudo[(~np.asarray(root_mask)) & finite])
        if np.isfinite(root_mean) and np.isfinite(other_mean) and root_mean > other_mean:
            pseudo = 1.0 - pseudo
    return pseudo


def get_root_index(adata, cluster_col: str, root_cluster: str) -> int:
    if cluster_col in adata.obs:
        mask = adata.obs[cluster_col].astype(str).to_numpy() == str(root_cluster)
    else:
        mask = np.zeros(adata.n_obs, dtype=bool)
    if not mask.any():
        return 0
    idx = np.flatnonzero(mask)
    if "X_umap" in adata.obsm and len(idx) > 1:
        emb = np.asarray(adata.obsm["X_umap"])
        center = np.nanmedian(emb[idx, :2], axis=0)
        dist = np.sum((emb[idx, :2] - center) ** 2, axis=1)
        return int(idx[np.nanargmin(dist)])
    return int(idx[0])


def ensure_scanpy_graph(adata, n_neighbors=20):
    import scanpy as sc

    if "X_pca" not in adata.obsm:
        n_comps = max(2, min(50, adata.n_obs - 1, adata.n_vars - 1))
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack", random_state=42)
    nn = max(2, min(int(n_neighbors), adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=nn, use_rep="X_pca", random_state=42)
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata, random_state=42)


def save_pseudotime_table(adata, method: str, pseudo_key: str, out_dir: Path, cluster_col: str, phenotype_col: str | None):
    df = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        "pseudotime": pd.to_numeric(adata.obs[pseudo_key], errors="coerce"),
    })
    if cluster_col in adata.obs:
        df[cluster_col] = adata.obs[cluster_col].astype(str).to_numpy()
        df["cluster_label"] = df[cluster_col]
    if phenotype_col and phenotype_col in adata.obs:
        df[phenotype_col] = adata.obs[phenotype_col].astype(str).to_numpy()
        df["phenotype"] = df[phenotype_col]
    if "sample" in adata.obs:
        df["sample"] = adata.obs["sample"].astype(str).to_numpy()
    df.to_csv(out_dir / f"{method}_pseudotime.tsv", sep="\t", index=False)
    return df


def plot_umap(adata, color_values, out_path: Path, title: str, cmap="viridis"):
    if "X_umap" not in adata.obsm:
        return None
    emb = np.asarray(adata.obsm["X_umap"])
    vals = np.asarray(color_values, dtype=float)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sca = ax.scatter(emb[:, 0], emb[:, 1], c=vals, s=8, cmap=cmap, linewidths=0, alpha=0.88)
    ax.set_title(title)
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.set_aspect("equal", adjustable="datalim")
    cbar = fig.colorbar(sca, ax=ax, shrink=0.82)
    cbar.set_label("pseudotime")
    fig.tight_layout()
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=180)
    plt.close(fig)
    return out_path


def plot_embedding_df(df: pd.DataFrame, x: str, y: str, color: str, out_path: Path, title: str):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    sca = ax.scatter(df[x], df[y], c=df[color], s=8, cmap="viridis", linewidths=0, alpha=0.88)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_aspect("equal", adjustable="datalim")
    fig.colorbar(sca, ax=ax, shrink=0.82).set_label(color)
    fig.tight_layout()
    fig.savefig(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=180)
    plt.close(fig)


def method_status(method: str, status: str, out_dir: Path, message: str = "", **extra):
    payload = {"method": method, "status": status, "message": message, **extra}
    write_json(payload, out_dir / f"{method}_status.json")
    return payload


def run_dpt(adata, out_root: Path, cluster_col: str, phenotype_col: str | None, root_idx: int, root_mask):
    method = "DPT"
    out_dir = ensure_dir(out_root / method)
    try:
        import scanpy as sc
        ensure_scanpy_graph(adata)
        sc.tl.diffmap(adata)
        adata.uns["iroot"] = int(root_idx)
        n_dcs = max(2, min(10, adata.obsm.get("X_diffmap", np.zeros((adata.n_obs, 3))).shape[1] - 1))
        sc.tl.dpt(adata, n_dcs=n_dcs)
        pseudo = orient_root_low(adata.obs["dpt_pseudotime"].to_numpy(), root_mask)
        adata.obs["DPT_pseudotime"] = pseudo
        df = save_pseudotime_table(adata, method, "DPT_pseudotime", out_dir, cluster_col, phenotype_col)
        plot_umap(adata, pseudo, out_dir / "DPT_umap_pseudotime.pdf", "DPT pseudotime (root-oriented)")
        return method_status(method, "ok", out_dir, n_cells=int(df["pseudotime"].notna().sum()))
    except Exception as e:
        (out_dir / "DPT_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return method_status(method, "error", out_dir, str(e))


def run_phate(adata, out_root: Path, cluster_col: str, phenotype_col: str | None, root_mask):
    method = "PHATE"
    out_dir = ensure_dir(out_root / method)
    try:
        import phate
        X, idx = dense_sample(adata.X, max_cells=5000)
        ph = phate.PHATE(n_components=2, knn=max(5, min(20, X.shape[0] - 2)), random_state=42, n_jobs=1, verbose=0)
        emb_sub = ph.fit_transform(X)
        if len(idx) == adata.n_obs:
            emb = emb_sub
        else:
            emb = np.full((adata.n_obs, 2), np.nan)
            emb[idx] = emb_sub
        adata.obsm["X_phate"] = emb
        pseudo = orient_root_low(emb[:, 0], root_mask)
        adata.obs["PHATE_pseudotime"] = pseudo
        df = save_pseudotime_table(adata, method, "PHATE_pseudotime", out_dir, cluster_col, phenotype_col)
        df["PHATE1"] = emb[:, 0]
        df["PHATE2"] = emb[:, 1]
        df.to_csv(out_dir / "PHATE_pseudotime.tsv", sep="\t", index=False)
        plot_embedding_df(df.dropna(subset=["PHATE1", "PHATE2", "pseudotime"]), "PHATE1", "PHATE2", "pseudotime", out_dir / "PHATE_embedding_pseudotime.pdf", "PHATE embedding ordered from root")
        return method_status(method, "ok", out_dir, "PHATE1 root-oriented ordering; PHATE is primarily an embedding method", n_cells=int(df["pseudotime"].notna().sum()))
    except Exception as e:
        (out_dir / "PHATE_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return method_status(method, "error", out_dir, str(e))


def run_palantir(adata, out_root: Path, cluster_col: str, phenotype_col: str | None, root_idx: int, root_mask):
    method = "Palantir"
    out_dir = ensure_dir(out_root / method)
    try:
        import palantir
        X = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        data_df = pd.DataFrame(X, index=adata.obs_names.astype(str), columns=adata.var_names.astype(str))
        dm_res = palantir.utils.run_diffusion_maps(data_df, n_components=max(5, min(10, adata.n_obs - 2)))
        ms_data = palantir.utils.determine_multiscale_space(dm_res)
        root_cell = str(adata.obs_names[root_idx])
        pr_res = palantir.core.run_palantir(ms_data, root_cell, num_waypoints=min(500, max(10, adata.n_obs - 1)))
        pseudo = pd.Series(pr_res.pseudotime).reindex(adata.obs_names.astype(str)).to_numpy(dtype=float)
        pseudo = orient_root_low(pseudo, root_mask)
        adata.obs["Palantir_pseudotime"] = pseudo
        df = save_pseudotime_table(adata, method, "Palantir_pseudotime", out_dir, cluster_col, phenotype_col)
        plot_umap(adata, pseudo, out_dir / "Palantir_umap_pseudotime.pdf", "Palantir pseudotime")
        branch_probs = getattr(pr_res, "branch_probs", None)
        if branch_probs is not None:
            pd.DataFrame(branch_probs).to_csv(out_dir / "Palantir_branch_probabilities.tsv", sep="\t")
        return method_status(method, "ok", out_dir, n_cells=int(df["pseudotime"].notna().sum()))
    except Exception as e:
        (out_dir / "Palantir_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return method_status(method, "error", out_dir, str(e))


def run_via(adata, out_root: Path, cluster_col: str, phenotype_col: str | None, root_idx: int, root_mask):
    method = "VIA"
    out_dir = ensure_dir(out_root / method)
    try:
        from pyVIA.core import VIA
        import pyVIA.core as via_core
        import pyVIA.utils_via as via_utils
        from scipy.sparse import csr_matrix

        def fixed_get_sparse_from_igraph(graph, weight_attr=None):
            edges = graph.get_edgelist()
            shape = (graph.vcount(), graph.vcount())
            if len(edges) == 0:
                return csr_matrix(shape)
            rows, cols = zip(*edges)
            weights = graph.es[weight_attr] if weight_attr is not None and weight_attr in graph.es.attributes() else np.ones(len(edges))
            return csr_matrix((weights, (rows, cols)), shape=shape)

        # pyVIA 0.2.4 uses csr_matrix((weights, zip(*edges))), which is not
        # robust with current SciPy. Patch locally without modifying site-packages.
        via_utils.get_sparse_from_igraph = fixed_get_sparse_from_igraph
        via_core.get_sparse_from_igraph = fixed_get_sparse_from_igraph
        ensure_scanpy_graph(adata)
        X = np.asarray(adata.obsm["X_pca"][:, : min(30, adata.obsm["X_pca"].shape[1])])
        root_candidates = []
        root_indices = np.flatnonzero(np.asarray(root_mask))
        if root_indices.size > 0:
            for q in [0.30, 0.25, 0.50, 0.75, 0.00, 1.00]:
                root_candidates.append(int(root_indices[min(root_indices.size - 1, max(0, int(round(q * (root_indices.size - 1)))))]))
        root_candidates.append(int(root_idx))
        root_candidates = list(dict.fromkeys(root_candidates))
        run_warning = None
        pseudo = None
        via_root_used = None
        last_error = None
        for via_root in root_candidates:
            via = VIA(
                X,
                labels=None,
                root_user=[int(via_root)],
                knn=max(5, min(30, adata.n_obs - 1)),
                random_seed=42,
                num_threads=1,
                do_compute_embedding=False,
                dataset="smc_anno",
                RW2_mode=True,
                num_mcmc_simulations=100,
                x_lazy=0.95,
                alpha_teleport=0.95,
            )
            try:
                via.run_VIA()
            except Exception as run_error:
                last_error = run_error
                run_warning = f"run_VIA raised for root {via_root}: {type(run_error).__name__}: {run_error}"
            for attr in ["single_cell_pt_markov", "single_cell_pt", "pseudotime", "single_cell_bp"]:
                val = getattr(via, attr, None)
                if val is not None:
                    try:
                        arr = np.asarray(val, dtype=float)
                    except Exception:
                        continue
                    if arr.shape[0] == adata.n_obs and np.isfinite(arr).any():
                        pseudo = arr
                        via_root_used = via_root
                        break
            if pseudo is not None:
                break
        if pseudo is None:
            if last_error is not None:
                raise RuntimeError(f"Could not find a VIA single-cell pseudotime vector after trying roots {root_candidates}; last error: {type(last_error).__name__}: {last_error}")
            raise RuntimeError(f"Could not find a VIA single-cell pseudotime vector after trying roots {root_candidates}")
        pseudo = orient_root_low(pseudo, root_mask)
        adata.obs["VIA_pseudotime"] = pseudo
        df = save_pseudotime_table(adata, method, "VIA_pseudotime", out_dir, cluster_col, phenotype_col)
        plot_umap(adata, pseudo, out_dir / "VIA_umap_pseudotime.pdf", "VIA pseudotime")
        return method_status(method, "ok", out_dir, run_warning or "ok", n_cells=int(df["pseudotime"].notna().sum()), via_root_used=via_root_used)
    except Exception as e:
        (out_dir / "VIA_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return method_status(method, "error", out_dir, str(e))


def run_cellrank(adata, out_root: Path, cluster_col: str, phenotype_col: str | None):
    method = "CellRank"
    out_dir = ensure_dir(out_root / method)
    try:
        import cellrank as cr
        ensure_scanpy_graph(adata)
        if "DPT_pseudotime" not in adata.obs:
            raise RuntimeError("CellRank PseudotimeKernel requires DPT_pseudotime from the DPT step")
        adata.obs["cellrank_time"] = pd.to_numeric(adata.obs["DPT_pseudotime"], errors="coerce").fillna(0).to_numpy()
        pk = cr.kernels.PseudotimeKernel(adata, time_key="cellrank_time").compute_transition_matrix()
        ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()
        rows = []
        for kernel_name, kernel in [("PseudotimeKernel", pk), ("ConnectivityKernel", ck)]:
            mat = kernel.transition_matrix
            rows.append({
                "kernel": kernel_name,
                "n_cells": int(mat.shape[0]),
                "nnz": int(mat.nnz) if sp.issparse(mat) else int(np.isfinite(mat).sum()),
                "density": float((mat.nnz if sp.issparse(mat) else np.isfinite(mat).sum()) / (mat.shape[0] * mat.shape[1])),
                "row_sum_min": float(np.nanmin(np.asarray(mat.sum(axis=1)).ravel())),
                "row_sum_max": float(np.nanmax(np.asarray(mat.sum(axis=1)).ravel())),
            })
        df = pd.DataFrame(rows)
        df.to_csv(out_dir / "CellRank_kernel_summary.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.bar(df["kernel"], df["density"], color=["#1F78B4", "#33A02C"])
        ax.set_ylabel("Transition matrix density")
        ax.set_title("CellRank kernel transition matrix density")
        fig.tight_layout()
        fig.savefig(out_dir / "CellRank_kernel_summary.pdf")
        fig.savefig(out_dir / "CellRank_kernel_summary.png", dpi=180)
        plt.close(fig)
        return method_status(method, "ok", out_dir, "Computed PseudotimeKernel and ConnectivityKernel transition matrices", n_kernels=2)
    except Exception as e:
        (out_dir / "CellRank_error.txt").write_text(traceback.format_exc(), encoding="utf-8")
        return method_status(method, "error", out_dir, str(e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--cluster-col", default="celltype")
    ap.add_argument("--phenotype-col", default="group")
    ap.add_argument("--root-cluster", default="Pericyte")
    args = ap.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.h5ad)
    out_root = ensure_dir(Path(args.output_dir))
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name in {"category", "object"}:
            adata.obs[col] = adata.obs[col].astype(str)
    root_idx = get_root_index(adata, args.cluster_col, args.root_cluster)
    root_mask = adata.obs[args.cluster_col].astype(str).to_numpy() == str(args.root_cluster) if args.cluster_col in adata.obs else np.zeros(adata.n_obs, dtype=bool)

    statuses = []
    statuses.append(run_dpt(adata, out_root, args.cluster_col, args.phenotype_col, root_idx, root_mask))
    statuses.append(run_phate(adata, out_root, args.cluster_col, args.phenotype_col, root_mask))
    statuses.append(run_palantir(adata, out_root, args.cluster_col, args.phenotype_col, root_idx, root_mask))
    statuses.append(run_via(adata, out_root, args.cluster_col, args.phenotype_col, root_idx, root_mask))
    statuses.append(run_cellrank(adata, out_root, args.cluster_col, args.phenotype_col))

    pd.DataFrame(statuses).to_csv(out_root / "python_first_tier_method_status.tsv", sep="\t", index=False)
    adata.write_h5ad(out_root / "python_first_tier_methods_with_outputs.h5ad", compression="gzip")
    write_json({"statuses": statuses, "h5ad": str(args.h5ad), "root_cell": str(adata.obs_names[root_idx])}, out_root / "python_first_tier_method_status.json")


if __name__ == "__main__":
    main()
