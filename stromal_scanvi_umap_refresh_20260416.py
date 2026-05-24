#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute diagnostic scanvi UMAP plots from current stromal branch h5ad files.

This script is intentionally lightweight and does not retrain models. It only:
1. loads the existing branch rerun h5ad
2. recomputes neighbors/UMAP from ``X_scanvi``
3. writes PNG/PDF diagnostic plots for quick visual inspection
4. optionally writes explicit ``X_umap_scanvi`` keys back to a sidecar h5ad
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import site
import sys
from pathlib import Path


PIPELINE_DATE = "20260416"
SEED = 42
N_NEIGHBORS = 30
DPI = 300
SCANVI_LATENT_KEY = "X_scanvi"
UMAP_KEY = "X_umap"
UMAP_SCANVI_KEY = "X_umap_scanvi"
UMAP_SCANVI_CORRECTED_KEY = "X_umap_scanvi_corrected"
TISSUE_KEY = "tissue"
L2_KEY = "cell_type_L2"
L3_KEY = "cell_type_L3"
PRED_KEY = "cell_type_scanvi_pred"
PRED_PROB_KEY = "scanvi_pred_prob_rerun"

BRANCH_H5AD = {
    "endothelial": Path(
        "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/endothelial/"
        "adata_endothelial_reference_rm_choir_5_33_20260414.h5ad"
    ),
    "fibroblast": Path(
        "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/fibroblast/"
        "adata_fibroblast_reference_rm_choir_23_59_20260414.h5ad"
    ),
    "smc": Path(
        "/home/h2048/data/py/0414/stromal_branch_rerun_20260414/smc/"
        "adata_smc_reference_neuronlike_c6_20260414.h5ad"
    ),
}


def configure_runtime() -> None:
    for key in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ.setdefault(key, "8")
    os.environ.setdefault("PYTHONNOUSERSITE", "1")
    user_site_candidates = []
    try:
        user_site_candidates.append(site.getusersitepackages())
    except Exception:
        pass
    user_site_candidates.extend(
        [
            os.path.expanduser("~/.local/lib/python3.10/site-packages"),
            os.path.expanduser("~/.local/lib/python3.11/site-packages"),
            os.path.expanduser("~/.local/lib/python3.9/site-packages"),
        ]
    )
    user_site_candidates = [p for p in user_site_candidates if isinstance(p, str) and p]
    sys.path[:] = [
        p for p in sys.path
        if not any(str(p).startswith(candidate) for candidate in user_site_candidates)
    ]


def import_runtime_modules():
    configure_runtime()
    import warnings

    warnings.filterwarnings("ignore")

    import anndata as ad  # type: ignore
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # type: ignore
    import numpy as np  # type: ignore
    import pandas as pd  # type: ignore
    import scanpy as sc  # type: ignore

    sc.settings.verbosity = 2
    sc.settings.n_jobs = 16
    sc.settings.set_figure_params(dpi=DPI, facecolor="white", frameon=False)
    np.random.seed(SEED)

    return ad, plt, np, pd, sc


def clear_stale_graph(adata) -> None:
    for key in ("neighbors", "umap"):
        if key in adata.uns:
            del adata.uns[key]
    for key in ("connectivities", "distances"):
        if key in adata.obsp:
            del adata.obsp[key]


def plot_umap_grid(adata, branch: str, sc_module, plt_module, png_path: Path, pdf_path: Path) -> None:
    panels = [
        (TISSUE_KEY, f"{branch}: tissue"),
        (L2_KEY, f"{branch}: L2"),
        (L3_KEY, f"{branch}: L3"),
        (PRED_KEY, f"{branch}: scanvi pred"),
    ]
    available = [(k, t) for k, t in panels if k in adata.obs.columns]
    n_panels = max(1, len(available))
    n_rows = 2
    n_cols = 2
    fig, axes = plt_module.subplots(n_rows, n_cols, figsize=(18, 14))
    axes = axes.flatten()

    for ax in axes[n_panels:]:
        ax.axis("off")

    for ax, (key, title) in zip(axes, available):
        kwargs = dict(
            adata=adata,
            basis="umap",
            color=key,
            ax=ax,
            show=False,
            title=title,
        )
        if key == PRED_PROB_KEY:
            kwargs.update(color_map="RdYlGn", vmin=0, vmax=1)
        else:
            kwargs.update(legend_loc="right margin", legend_fontsize=6)
        sc_module.pl.embedding(**kwargs)

    plt_module.tight_layout()
    fig.savefig(png_path, bbox_inches="tight", dpi=DPI)
    fig.savefig(pdf_path, bbox_inches="tight", dpi=DPI)
    plt_module.close(fig)


def run_branch(branch: str, h5ad_path: Path, write_sidecar: bool) -> dict[str, object]:
    ad, plt, np, pd, sc = import_runtime_modules()
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing h5ad for {branch}: {h5ad_path}")

    figures_dir = h5ad_path.parent / "figures"
    reports_dir = h5ad_path.parent / "reports"
    figures_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad_path)
    if SCANVI_LATENT_KEY not in adata.obsm:
        raise KeyError(f"{branch} is missing required obsm key: {SCANVI_LATENT_KEY}")

    clear_stale_graph(adata)
    sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY, n_neighbors=N_NEIGHBORS)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)

    coords = np.asarray(adata.obsm[UMAP_KEY], dtype=np.float32)
    adata.obsm[UMAP_KEY] = coords
    adata.obsm[UMAP_SCANVI_KEY] = coords.copy()
    adata.obsm[UMAP_SCANVI_CORRECTED_KEY] = coords.copy()

    png_path = figures_dir / f"{branch}_scanvi_umap_refresh_{PIPELINE_DATE}.png"
    pdf_path = figures_dir / f"{branch}_scanvi_umap_refresh_{PIPELINE_DATE}.pdf"
    plot_umap_grid(adata, branch, sc, plt, png_path, pdf_path)

    sidecar_path = None
    if write_sidecar:
        sidecar_path = h5ad_path.with_name(h5ad_path.stem + f"_scanvi_umap_refresh_{PIPELINE_DATE}.h5ad")
        adata.uns.setdefault("scanvi_umap_refresh", {})
        adata.uns["scanvi_umap_refresh"] = {
            "source_h5ad": str(h5ad_path),
            "pipeline_date": PIPELINE_DATE,
            "neighbors_use_rep": SCANVI_LATENT_KEY,
            "output_png": str(png_path),
            "output_pdf": str(pdf_path),
        }
        import anndata as _anndata

        _anndata.settings.allow_write_nullable_strings = True
        adata.write_h5ad(sidecar_path, compression="gzip", compression_opts=9)

    summary = {
        "branch": branch,
        "source_h5ad": str(h5ad_path),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obsm_keys": list(adata.obsm.keys()),
        "png": str(png_path),
        "pdf": str(pdf_path),
        "sidecar_h5ad": None if sidecar_path is None else str(sidecar_path),
    }
    with open(reports_dir / f"{branch}_scanvi_umap_refresh_{PIPELINE_DATE}.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    del adata
    gc.collect()
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Refresh scanvi UMAP plots for stromal branch h5ad files")
    parser.add_argument("--branch", choices=["all", *BRANCH_H5AD.keys()], default="all")
    parser.add_argument("--write-sidecar-h5ad", action="store_true", help="Write a sidecar h5ad carrying refreshed UMAP keys")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    branches = BRANCH_H5AD.items() if args.branch == "all" else [(args.branch, BRANCH_H5AD[args.branch])]
    summaries = []
    for branch, path in branches:
        print(f"[RUN] {branch}: {path}")
        summary = run_branch(branch, path, write_sidecar=args.write_sidecar_h5ad)
        summaries.append(summary)
        print(f"[OK] {branch}: {summary['png']}")
    print(json.dumps(summaries, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()