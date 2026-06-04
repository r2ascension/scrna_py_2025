#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run pertpy Milo differential abundance analysis across tissues within L2/L3 subclusters.

Notes
-----
- Uses a lightweight AnnData containing only obs + embeddings to avoid loading the
  full gene matrix when the goal is neighbourhood-level DA.
- Runs pertpy Milo separately for each requested cell-type column (for example
  `cell_type_L2` and `cell_type_L3`), testing tissue contrasts using sample as the
  replicate unit.
- Writes Milo's graph-aware `SpatialFDR` together with within-level BH-adjusted
  p-values.
- Produces per-level summaries plus a combined cross-level index for downstream R
  reporting.

Recommended invocation environment
----------------------------------
Use the pertpy environment directly, keep user-site packages disabled, and
prepend the environment `lib/` ahead of system R libraries:

        ENV_PREFIX=/path/to/pertpy_env
        PYTHONNOUSERSITE=1 \
            LD_LIBRARY_PATH=$ENV_PREFIX/lib:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/default-java/lib/server \
            $ENV_PREFIX/bin/python tissue_subcluster_milopy_20260418_v1.py
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from dataclasses import asdict, dataclass
from itertools import combinations
from math import ceil
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

ACTIVE_ENV_PREFIX = str(Path(sys.executable).resolve().parents[1])
R_RUNTIME_LD_LIBRARY_PATH_HINT = (
    f"{ACTIVE_ENV_PREFIX}/lib:/usr/lib/R/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/default-java/lib/server"
)

try:
    import pertpy as pt
except Exception as exc:  # pragma: no cover - environment/runtime failure path
    raise SystemExit(
        "Failed to import pertpy. Run this script with "
        f"PYTHONNOUSERSITE=1 LD_LIBRARY_PATH={R_RUNTIME_LD_LIBRARY_PATH_HINT} "
        f"{ACTIVE_ENV_PREFIX}/bin/python. Original error: {exc}"
    ) from exc


LOGGER = logging.getLogger("tissue_subcluster_milo")


def patch_pertpy_rpy2_setup() -> None:
    """Patch pertpy Milo for rpy2>=3.6 where activate()/deactivate() are deprecated."""

    import pertpy.tools._milo as milo_mod
    from rpy2.robjects import conversion, default_converter, numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr

    def _setup_rpy2(self):
        conversion.set_conversion(default_converter + numpy2ri.converter + pandas2ri.converter)
        edgeR = self._try_import_bioc_library("edgeR")
        limma = self._try_import_bioc_library("limma")
        stats = importr("stats")
        base = importr("base")
        return edgeR, limma, stats, base

    milo_mod.Milo._setup_rpy2 = _setup_rpy2


@dataclass
class MiloRunConfig:
    input_h5ad: str = ""
    output_dir: str = ""
    sample_col: str = "sample"
    tissue_col: str = "tissue"
    celltype_cols: tuple[str, ...] = ("cell_type_L2", "cell_type_L3")
    latent_key_candidates: tuple[str, ...] = (
        "X_scanvi",
        "X_scANVI",
        "X_scanvi_corrected",
        "X_scvi",
        "X_harmony",
        "X_pca",
    )
    umap_key_candidates: tuple[str, ...] = (
        "X_umap_scanvi",
        "X_umap_scANVI",
        "X_umap_scanvi_corrected",
        "X_umap_scvi",
        "X_umap_harmony",
        "X_umap",
    )
    exclude_celltypes: tuple[str, ...] = ("", "Rejected")
    min_cells_per_celltype: int = 200
    min_cells_per_sample: int = 10
    min_samples_per_tissue: int = 3
    n_neighbors: int = 30
    nhood_prop: float = 0.1
    alpha: float = 0.10
    random_seed: int = 0
    max_celltypes_per_level: int | None = None
    make_plots: bool = True
    write_milo_h5ad: bool = False


def parse_args() -> MiloRunConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-col", default=MiloRunConfig.sample_col)
    parser.add_argument("--tissue-col", default=MiloRunConfig.tissue_col)
    parser.add_argument("--celltype-col", action="append", default=None)
    parser.add_argument("--latent-key-candidate", action="append", default=None)
    parser.add_argument("--umap-key-candidate", action="append", default=None)
    parser.add_argument("--exclude-celltypes", nargs="*", default=list(MiloRunConfig.exclude_celltypes))
    parser.add_argument("--min-cells-per-celltype", type=int, default=MiloRunConfig.min_cells_per_celltype)
    parser.add_argument("--min-cells-per-sample", type=int, default=MiloRunConfig.min_cells_per_sample)
    parser.add_argument("--min-samples-per-tissue", type=int, default=MiloRunConfig.min_samples_per_tissue)
    parser.add_argument("--n-neighbors", type=int, default=MiloRunConfig.n_neighbors)
    parser.add_argument("--nhood-prop", type=float, default=MiloRunConfig.nhood_prop)
    parser.add_argument("--alpha", type=float, default=MiloRunConfig.alpha)
    parser.add_argument("--random-seed", type=int, default=MiloRunConfig.random_seed)
    parser.add_argument("--max-celltypes-per-level", type=int, default=None)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--write-milo-h5ad", action="store_true")
    args = parser.parse_args()

    celltype_cols = tuple(args.celltype_col) if args.celltype_col else MiloRunConfig.celltype_cols
    latent_keys = tuple(args.latent_key_candidate) if args.latent_key_candidate else MiloRunConfig.latent_key_candidates
    umap_keys = tuple(args.umap_key_candidate) if args.umap_key_candidate else MiloRunConfig.umap_key_candidates

    return MiloRunConfig(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        sample_col=args.sample_col,
        tissue_col=args.tissue_col,
        celltype_cols=celltype_cols,
        latent_key_candidates=latent_keys,
        umap_key_candidates=umap_keys,
        exclude_celltypes=tuple(args.exclude_celltypes),
        min_cells_per_celltype=args.min_cells_per_celltype,
        min_cells_per_sample=args.min_cells_per_sample,
        min_samples_per_tissue=args.min_samples_per_tissue,
        n_neighbors=args.n_neighbors,
        nhood_prop=args.nhood_prop,
        alpha=args.alpha,
        random_seed=args.random_seed,
        max_celltypes_per_level=args.max_celltypes_per_level,
        make_plots=not args.no_plots,
        write_milo_h5ad=bool(args.write_milo_h5ad),
    )


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def slugify(value: str) -> str:
    value = re.sub(r"[^0-9a-zA-Z]+", "_", str(value).strip())
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        value = "missing"
    if value[0].isdigit():
        value = f"x_{value}"
    return value.lower()


def infer_level_label(celltype_col: str) -> str:
    col_upper = str(celltype_col).upper()
    if col_upper.endswith("_L2") or col_upper == "L2":
        return "L2"
    if col_upper.endswith("_L3") or col_upper == "L3":
        return "L3"
    return slugify(celltype_col).upper()


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def pick_first_existing(candidates: Iterable[str], available: Iterable[str]) -> str | None:
    available_set = set(map(str, available))
    for key in candidates:
        if key in available_set:
            return key
    return None


def lightweight_adata_from_h5ad(input_h5ad: Path, cfg: MiloRunConfig) -> tuple[ad.AnnData, dict[str, str | None]]:
    LOGGER.info("Loading lightweight AnnData from %s", input_h5ad)
    src = ad.read_h5ad(input_h5ad, backed="r")
    try:
        required_cols = [cfg.sample_col, cfg.tissue_col, *cfg.celltype_cols]
        missing_cols = [col for col in required_cols if col not in src.obs.columns]
        if missing_cols:
            raise KeyError(f"Missing required obs columns: {missing_cols}")

        latent_key = pick_first_existing(cfg.latent_key_candidates, src.obsm.keys())
        if latent_key is None:
            raise KeyError(
                "Missing latent embedding in obsm. Checked: "
                f"{list(cfg.latent_key_candidates)} | available={list(src.obsm.keys())}"
            )
        umap_key = pick_first_existing(cfg.umap_key_candidates, src.obsm.keys())

        obs = src.obs.copy()
        dummy_x = sparse.csr_matrix((src.n_obs, 1), dtype=np.float32)
        adata = ad.AnnData(X=dummy_x, obs=obs, var=pd.DataFrame(index=["placeholder_feature"]))
        adata.obsm[latent_key] = np.asarray(src.obsm[latent_key])
        if umap_key is not None:
            adata.obsm[umap_key] = np.asarray(src.obsm[umap_key])
        LOGGER.info("Lightweight AnnData shape: %s | latent=%s | umap=%s", adata.shape, latent_key, umap_key)
        return adata, {"latent_key": latent_key, "umap_key": umap_key}
    finally:
        if getattr(src, "file", None) is not None:
            src.file.close()


def clean_obs_columns(adata: ad.AnnData, cfg: MiloRunConfig) -> None:
    target_cols = [cfg.sample_col, cfg.tissue_col, *cfg.celltype_cols]
    for col in target_cols:
        adata.obs[col] = adata.obs[col].astype(str).str.strip()
        adata.obs[col] = adata.obs[col].replace({"nan": "", "None": ""})


def summarize_celltypes(adata: ad.AnnData, cfg: MiloRunConfig, celltype_col: str, level_label: str) -> pd.DataFrame:
    obs = adata.obs[[celltype_col, cfg.sample_col, cfg.tissue_col]].copy()
    obs = obs.loc[
        ~obs[celltype_col].isin(cfg.exclude_celltypes)
        & (obs[cfg.sample_col] != "")
        & (obs[cfg.tissue_col] != "")
    ].copy()
    rows: list[dict] = []
    for cell_type, sub in obs.groupby(celltype_col, sort=False):
        tissue_sample_counts = sub.groupby(cfg.tissue_col)[cfg.sample_col].nunique().to_dict()
        tissue_cell_counts = sub[cfg.tissue_col].value_counts().to_dict()
        rows.append(
            {
                "level": level_label,
                "celltype_col": celltype_col,
                "cell_type": cell_type,
                "n_cells": int(sub.shape[0]),
                "n_samples": int(sub[cfg.sample_col].nunique()),
                "n_tissues": int(sub[cfg.tissue_col].nunique()),
                "tissue_sample_counts": json.dumps(tissue_sample_counts, ensure_ascii=False, sort_keys=True),
                "tissue_cell_counts": json.dumps(tissue_cell_counts, ensure_ascii=False, sort_keys=True),
            }
        )
    return pd.DataFrame(rows).sort_values(["n_cells", "n_samples"], ascending=[False, False])


def prepare_celltype_subset(
    adata: ad.AnnData,
    cell_type: str,
    celltype_col: str,
    level_label: str,
    cfg: MiloRunConfig,
) -> tuple[ad.AnnData | None, dict]:
    subset = adata[adata.obs[celltype_col] == cell_type].copy()
    info: dict[str, object] = {
        "level": level_label,
        "celltype_col": celltype_col,
        "cell_type": cell_type,
        "n_cells_initial": int(subset.n_obs),
    }
    if subset.n_obs < cfg.min_cells_per_celltype:
        info["status"] = "skipped"
        info["reason"] = f"fewer than {cfg.min_cells_per_celltype} cells"
        return None, info

    sample_tissue = subset.obs[[cfg.sample_col, cfg.tissue_col]].drop_duplicates()
    ambiguous_samples = (
        sample_tissue.groupby(cfg.sample_col)[cfg.tissue_col].nunique().loc[lambda s: s > 1].index.tolist()
    )
    if ambiguous_samples:
        subset = subset[~subset.obs[cfg.sample_col].isin(ambiguous_samples)].copy()
        info["ambiguous_samples_dropped"] = len(ambiguous_samples)

    sample_sizes = subset.obs[cfg.sample_col].value_counts()
    keep_samples = sample_sizes.loc[lambda s: s >= cfg.min_cells_per_sample].index.tolist()
    subset = subset[subset.obs[cfg.sample_col].isin(keep_samples)].copy()
    info["n_cells_after_sample_filter"] = int(subset.n_obs)
    info["n_samples_after_sample_filter"] = int(subset.obs[cfg.sample_col].nunique())

    if subset.n_obs < cfg.min_cells_per_celltype:
        info["status"] = "skipped"
        info["reason"] = f"fewer than {cfg.min_cells_per_celltype} cells after sample filter"
        return None, info

    tissue_sample_counts = subset.obs.groupby(cfg.tissue_col)[cfg.sample_col].nunique().sort_values(ascending=False)
    keep_tissues = tissue_sample_counts.loc[lambda s: s >= cfg.min_samples_per_tissue].index.tolist()
    subset = subset[subset.obs[cfg.tissue_col].isin(keep_tissues)].copy()
    info["n_cells_after_tissue_filter"] = int(subset.n_obs)
    info["tissue_sample_counts"] = json.dumps(
        subset.obs.groupby(cfg.tissue_col)[cfg.sample_col].nunique().to_dict(), ensure_ascii=False, sort_keys=True
    )
    info["tissue_cell_counts"] = json.dumps(
        subset.obs[cfg.tissue_col].value_counts().to_dict(), ensure_ascii=False, sort_keys=True
    )

    if subset.obs[cfg.tissue_col].nunique() < 2:
        info["status"] = "skipped"
        info["reason"] = f"fewer than 2 tissues with >= {cfg.min_samples_per_tissue} samples"
        return None, info

    tissue_map = {t: slugify(t) for t in sorted(subset.obs[cfg.tissue_col].unique())}
    subset.obs["tissue_safe"] = pd.Categorical(subset.obs[cfg.tissue_col].map(tissue_map))
    info["status"] = "ready"
    info["reason"] = ""
    return subset, info


def build_tissue_pairs(subset: ad.AnnData, cfg: MiloRunConfig) -> list[tuple[str, str]]:
    tissues = sorted(subset.obs[cfg.tissue_col].astype(str).unique().tolist())
    return list(combinations(tissues, 2))


def ensure_neighbors_and_umap(subset: ad.AnnData, latent_key: str, umap_key: str, cfg: MiloRunConfig) -> None:
    sc.pp.neighbors(subset, use_rep=latent_key, n_neighbors=cfg.n_neighbors, random_state=cfg.random_seed)
    if umap_key not in subset.obsm:
        sc.tl.umap(subset, random_state=cfg.random_seed)
        if "X_umap" in subset.obsm and umap_key != "X_umap":
            subset.obsm[umap_key] = subset.obsm["X_umap"].copy()


def plot_milo_graph(milo: pt.tl.Milo, mdata, fig_path: Path, title: str, alpha: float) -> None:
    fig = milo.plot_nhood_graph(mdata, alpha=alpha, title=title, return_fig=True)
    if hasattr(fig, "savefig"):
        fig.savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
    else:  # pragma: no cover - defensive
        plt.gcf().savefig(fig_path, dpi=200, bbox_inches="tight")
        plt.close(plt.gcf())


def compute_sample_fraction_table(
    full_adata: ad.AnnData,
    subset: ad.AnnData,
    cell_type: str,
    level_label: str,
    celltype_col: str,
    cfg: MiloRunConfig,
) -> pd.DataFrame:
    sample_meta = subset.obs[[cfg.sample_col, cfg.tissue_col]].drop_duplicates().copy()
    all_counts = full_adata.obs.groupby(cfg.sample_col).size().rename("total_lineage_cells")
    label_counts = subset.obs.groupby(cfg.sample_col).size().rename("label_cells")
    plot_df = sample_meta.copy()
    plot_df["label_cells"] = plot_df[cfg.sample_col].map(label_counts).fillna(0).astype(int)
    plot_df["celltype_cells"] = plot_df["label_cells"]
    plot_df["total_lineage_cells"] = plot_df[cfg.sample_col].map(all_counts).fillna(0).astype(int)
    plot_df = plot_df.loc[plot_df["total_lineage_cells"] > 0].copy()
    plot_df["sample_fraction_pct"] = 100.0 * plot_df["label_cells"] / plot_df["total_lineage_cells"]
    plot_df["sample"] = plot_df[cfg.sample_col].astype(str)
    plot_df["tissue"] = plot_df[cfg.tissue_col].astype(str)
    plot_df["analysis_level"] = level_label
    plot_df["cell_type"] = cell_type
    plot_df["level"] = level_label
    plot_df["celltype_col"] = celltype_col
    plot_df = plot_df.sort_values([cfg.tissue_col, cfg.sample_col]).reset_index(drop=True)
    return plot_df


def pvalue_to_stars(p: float) -> str:
    if pd.isna(p):
        return "na"
    if p <= 1e-4:
        return "****"
    if p <= 1e-3:
        return "***"
    if p <= 1e-2:
        return "**"
    if p <= 5e-2:
        return "*"
    return "ns"


def compute_violin_pairwise_stats(plot_df: pd.DataFrame, cfg: MiloRunConfig) -> tuple[pd.DataFrame, float | None]:
    tissue_order = plot_df[cfg.tissue_col].astype(str).drop_duplicates().tolist()
    rows: list[dict] = []
    grouped = {t: plot_df.loc[plot_df[cfg.tissue_col] == t, "sample_fraction_pct"].astype(float).values for t in tissue_order}

    overall_p: float | None = None
    if len(tissue_order) >= 2 and all(len(grouped[t]) > 0 for t in tissue_order):
        try:
            overall_p = float(kruskal(*[grouped[t] for t in tissue_order]).pvalue)
        except ValueError:
            overall_p = None

    for tissue_a, tissue_b in combinations(tissue_order, 2):
        values_a = grouped[tissue_a]
        values_b = grouped[tissue_b]
        if len(values_a) == 0 or len(values_b) == 0:
            p_raw = np.nan
        else:
            p_raw = float(mannwhitneyu(values_a, values_b, alternative="two-sided").pvalue)
        rows.append(
            {
                "analysis_level": plot_df["analysis_level"].iloc[0],
                "level": plot_df["level"].iloc[0],
                "celltype_col": plot_df["celltype_col"].iloc[0],
                "cell_type": plot_df["cell_type"].iloc[0],
                "tissue_a": tissue_a,
                "tissue_b": tissue_b,
                "n_samples_a": int(len(values_a)),
                "n_samples_b": int(len(values_b)),
                "median_a": float(np.median(values_a)) if len(values_a) else np.nan,
                "median_b": float(np.median(values_b)) if len(values_b) else np.nan,
                "p_raw": p_raw,
                "overall_kruskal_p": overall_p,
            }
        )

    stats_df = pd.DataFrame(rows)
    if stats_df.empty:
        return stats_df, overall_p

    valid = stats_df["p_raw"].notna()
    stats_df["p_bh"] = np.nan
    if valid.any():
        stats_df.loc[valid, "p_bh"] = multipletests(stats_df.loc[valid, "p_raw"], method="fdr_bh")[1]
    stats_df["significant_bh"] = stats_df["p_bh"].lt(cfg.alpha)
    stats_df["label_bh"] = stats_df["p_bh"].map(pvalue_to_stars)
    return stats_df, overall_p


def annotate_violin_significance(
    ax: plt.Axes,
    stats_df: pd.DataFrame,
    tissue_order: list[str],
    y_max: float,
    y_step: float,
) -> None:
    if stats_df.empty:
        ax.text(
            0.99,
            0.98,
            "no pairwise tests",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            color="dimgray",
        )
        return

    y = y_max + y_step
    order_index = {name: idx for idx, name in enumerate(tissue_order)}
    for _, row in stats_df.sort_values(["p_bh", "tissue_a", "tissue_b"], na_position="last").iterrows():
        x1 = order_index[row["tissue_a"]]
        x2 = order_index[row["tissue_b"]]
        if x1 == x2:
            continue
        color = "black" if bool(row.get("significant_bh", False)) else "dimgray"
        ax.plot([x1, x1, x2, x2], [y, y + y_step * 0.25, y + y_step * 0.25, y], lw=1.2, c=color)
        ax.text((x1 + x2) / 2, y + y_step * 0.3, row["label_bh"], ha="center", va="bottom", fontsize=10, color=color)
        y += y_step


def draw_sample_fraction_violin(
    ax: plt.Axes,
    plot_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    cfg: MiloRunConfig,
    *,
    title: str,
    ylabel: str,
    show_legend: bool = False,
) -> None:
    tissue_order = plot_df[cfg.tissue_col].astype(str).drop_duplicates().tolist()
    counts_by_tissue = plot_df.groupby(cfg.tissue_col)[cfg.sample_col].nunique().reindex(tissue_order)
    palette = sns.color_palette("Set2", n_colors=len(tissue_order))

    sns.violinplot(
        data=plot_df,
        x=cfg.tissue_col,
        y="sample_fraction_pct",
        hue=cfg.tissue_col,
        order=tissue_order,
        palette=palette,
        inner=None,
        cut=0,
        linewidth=1,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x=cfg.tissue_col,
        y="sample_fraction_pct",
        order=tissue_order,
        color="black",
        size=4,
        jitter=0.18,
        alpha=0.6,
        ax=ax,
    )

    overall_p = None if stats_df.empty else stats_df["overall_kruskal_p"].dropna().iloc[0] if stats_df["overall_kruskal_p"].notna().any() else None
    title_suffix = f"\nKruskal p={overall_p:.2e}" if overall_p is not None else ""
    ax.set_title(title + title_suffix, fontsize=11)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(tissue_order)))
    ax.set_xticklabels([f"{t}\n(n={int(counts_by_tissue[t])})" for t in tissue_order], rotation=0)

    y_min = float(plot_df["sample_fraction_pct"].min())
    y_max = float(plot_df["sample_fraction_pct"].max())
    y_range = max(y_max - y_min, 1e-3)
    y_step = y_range * 0.10
    annotate_violin_significance(ax, stats_df, tissue_order, y_max=y_max, y_step=y_step)
    extra_levels = max(int(stats_df.shape[0]), 1)
    ax.set_ylim(bottom=min(0.0, y_min - y_range * 0.05), top=y_max + y_step * (extra_levels + 2))
    ax.grid(axis="y", linestyle="--", alpha=0.25)

    if not show_legend and ax.get_legend() is not None:
        ax.get_legend().remove()


def save_violin_outputs(
    full_adata: ad.AnnData,
    subset: ad.AnnData,
    cell_type: str,
    level_label: str,
    celltype_col: str,
    cfg: MiloRunConfig,
    celltype_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    plot_df = compute_sample_fraction_table(full_adata, subset, cell_type, level_label, celltype_col, cfg)
    stats_df, _ = compute_violin_pairwise_stats(plot_df, cfg)
    plot_df.to_csv(celltype_dir / "sample_fraction_per_sample.csv", index=False)
    stats_df.to_csv(celltype_dir / "sample_fraction_pairwise_stats.csv", index=False)

    if cfg.make_plots:
        fig, ax = plt.subplots(figsize=(max(5, 1.6 * plot_df[cfg.tissue_col].nunique()), 5.2))
        draw_sample_fraction_violin(
            ax,
            plot_df,
            stats_df,
            cfg,
            title=f"{level_label} | {cell_type} sample-level abundance across tissues",
            ylabel="Cell type fraction per sample (%)",
        )
        fig.tight_layout()
        fig.savefig(celltype_dir / "sample_fraction_violin.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    return plot_df, stats_df


def plot_violin_overview(
    plot_tables: list[pd.DataFrame],
    stats_tables: list[pd.DataFrame],
    cfg: MiloRunConfig,
    output_dir: Path,
    level_label: str,
) -> None:
    if not plot_tables:
        return
    n_panels = len(plot_tables)
    n_cols = min(3, n_panels)
    n_rows = ceil(n_panels / n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4.8 * n_rows), squeeze=False)

    for ax, plot_df, stats_df in zip(axes.ravel(), plot_tables, stats_tables, strict=False):
        draw_sample_fraction_violin(
            ax,
            plot_df,
            stats_df,
            cfg,
            title=str(plot_df["cell_type"].iloc[0]),
            ylabel="Fraction (%)",
        )

    for ax in axes.ravel()[len(plot_tables):]:
        ax.axis("off")

    fig.suptitle(f"{level_label} subcluster abundance across tissues", fontsize=14, y=0.995)
    fig.tight_layout()
    fig.savefig(output_dir / "sample_fraction_violin_overview.png", dpi=220, bbox_inches="tight")
    fig.savefig(output_dir / f"{level_label.lower()}_sample_fraction_violin_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_da_nhoods_edger_compat(
    milo: pt.tl.Milo,
    mdata,
    *,
    design: str,
    model_contrasts: str | None,
    subset_samples: list[str] | None,
    feature_key: str = "rna",
) -> None:
    """Compatibility implementation of Milo DA using edgeR via rpy2."""

    from rpy2.robjects import conversion, default_converter, numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    sample_adata = mdata["milo"]
    adata = mdata[feature_key]
    covariates = [x.strip(" ") for x in set(re.split("\\+|\\*", design.lstrip("~ ")))]
    sample_col = sample_adata.uns["sample_col"]

    sample_obs = adata.obs[covariates + [sample_col]].drop_duplicates().copy()
    sample_obs.index = sample_obs[sample_col].astype(str)
    sample_adata.obs = sample_obs.loc[sample_adata.obs_names].copy()
    design_df = sample_adata.obs[covariates].copy()

    count_mat = sample_adata.X.T.toarray()
    lib_size = count_mat.sum(0)
    keep_smp = lib_size > 0

    if subset_samples is not None:
        keep_smp = keep_smp & sample_adata.obs_names.isin(subset_samples)
        design_df = design_df.loc[keep_smp].copy()
        for column in design_df.columns:
            if pd.api.types.is_categorical_dtype(design_df[column]):
                design_df[column] = design_df[column].cat.remove_unused_categories()

    keep_nhoods = count_mat[:, keep_smp].sum(1) > 0
    counts = np.asarray(count_mat[keep_nhoods, :][:, keep_smp], dtype=float)
    libs = np.asarray(lib_size[keep_smp], dtype=float)

    edgeR = importr("edgeR")
    limma = importr("limma")
    stats = importr("stats")
    base = importr("base")

    if model_contrasts is not None:
        design = design + " + 0"

    with localconverter(default_converter + pandas2ri.converter):
        r_design_df = conversion.py2rpy(design_df)
    with localconverter(default_converter + numpy2ri.converter):
        r_counts = conversion.py2rpy(counts)
        r_libs = conversion.py2rpy(libs)

    model = stats.model_matrix(object=stats.formula(design), data=r_design_df)
    dge = edgeR.DGEList(counts=r_counts, lib_size=r_libs)
    dge = edgeR.calcNormFactors(dge, method="TMM")
    dge = edgeR.estimateDisp(dge, model)
    fit = edgeR.glmQLFit(dge, model, robust=True)

    if model_contrasts is not None:
        try:
            mod_contrast = limma.makeContrasts(contrasts=model_contrasts, levels=model)
        except Exception as exc:  # pragma: no cover - contrast syntax guard
            raise ValueError(f"Invalid model contrast: {model_contrasts}") from exc
        res_r = base.as_data_frame(edgeR.topTags(edgeR.glmQLFTest(fit, contrast=mod_contrast), sort_by="none", n=counts.shape[0]))
    else:
        n_coef = int(base.ncol(model)[0])
        res_r = base.as_data_frame(edgeR.topTags(edgeR.glmQLFTest(fit, coef=n_coef), sort_by="none", n=counts.shape[0]))

    with localconverter(default_converter + pandas2ri.converter):
        res = conversion.rpy2py(res_r)
    if not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)

    res.index = sample_adata.var_names[keep_nhoods]  # type: ignore[index]
    existing_cols = [col for col in res.columns if col in sample_adata.var.columns]
    if existing_cols:
        sample_adata.var = sample_adata.var.drop(existing_cols, axis=1)
    sample_adata.var = pd.concat([sample_adata.var, res], axis=1)
    milo._graph_spatial_fdr(sample_adata, neighbors_key=adata.uns["nhood_neighbors_key"])


def run_milo_for_celltype(
    full_adata: ad.AnnData,
    subset: ad.AnnData,
    cell_type: str,
    celltype_col: str,
    level_label: str,
    latent_key: str,
    umap_key: str,
    cfg: MiloRunConfig,
    output_dir: Path,
) -> tuple[list[pd.DataFrame], list[dict], pd.DataFrame, pd.DataFrame]:
    LOGGER.info("Running Milo for %s | %s (%d cells)", level_label, cell_type, subset.n_obs)
    celltype_slug = slugify(cell_type)
    celltype_dir = output_dir / celltype_slug
    celltype_dir.mkdir(parents=True, exist_ok=True)

    ensure_neighbors_and_umap(subset, latent_key, umap_key, cfg)
    milo = pt.tl.Milo()
    mdata = milo.load(subset)
    milo.make_nhoods(mdata["rna"], prop=cfg.nhood_prop, seed=cfg.random_seed)
    mdata = milo.count_nhoods(mdata, sample_col=cfg.sample_col)
    milo.annotate_nhoods(mdata, anno_col=cfg.tissue_col)
    milo.build_nhood_graph(mdata, basis=umap_key)
    plot_df, plot_stats_df = save_violin_outputs(full_adata, subset, cell_type, level_label, celltype_col, cfg, celltype_dir)

    pair_summaries: list[dict] = []
    results: list[pd.DataFrame] = []
    tissue_pairs = build_tissue_pairs(subset, cfg)
    if not tissue_pairs:
        pair_summaries.append(
            {
                "level": level_label,
                "celltype_col": celltype_col,
                "cell_type": cell_type,
                "status": "skipped",
                "reason": "no eligible tissue pairs",
            }
        )
        return results, pair_summaries, plot_df, plot_stats_df

    sample_meta = subset.obs[[cfg.sample_col, cfg.tissue_col, "tissue_safe"]].drop_duplicates().copy()
    tissue_map = dict(zip(sample_meta[cfg.tissue_col], sample_meta["tissue_safe"], strict=False))

    for tissue_a, tissue_b in tissue_pairs:
        samples_a = sorted(sample_meta.loc[sample_meta[cfg.tissue_col] == tissue_a, cfg.sample_col].tolist())
        samples_b = sorted(sample_meta.loc[sample_meta[cfg.tissue_col] == tissue_b, cfg.sample_col].tolist())
        if len(samples_a) < cfg.min_samples_per_tissue or len(samples_b) < cfg.min_samples_per_tissue:
            pair_summaries.append(
                {
                    "level": level_label,
                    "celltype_col": celltype_col,
                    "cell_type": cell_type,
                    "tissue_a": tissue_a,
                    "tissue_b": tissue_b,
                    "n_samples_a": len(samples_a),
                    "n_samples_b": len(samples_b),
                    "status": "skipped",
                    "reason": "insufficient samples after filtering",
                }
            )
            continue

        contrast_label = f"{slugify(tissue_b)}_vs_{slugify(tissue_a)}"
        comparison_samples = samples_a + samples_b
        model_contrasts = f"tissue_safe{tissue_map[tissue_b]} - tissue_safe{tissue_map[tissue_a]}"
        LOGGER.info(
            "  Contrast %s | %s (%d samples) vs %s (%d samples)",
            contrast_label,
            tissue_b,
            len(samples_b),
            tissue_a,
            len(samples_a),
        )

        run_da_nhoods_edger_compat(
            milo,
            mdata,
            design="~ tissue_safe",
            model_contrasts=model_contrasts,
            subset_samples=comparison_samples,
            feature_key="rna",
        )

        res = mdata["milo"].var.copy().reset_index(names="nhood")
        valid = res["PValue"].notna()
        res["PValue_BH"] = np.nan
        if valid.any():
            res.loc[valid, "PValue_BH"] = multipletests(res.loc[valid, "PValue"], method="fdr_bh")[1]
        res["level"] = level_label
        res["celltype_col"] = celltype_col
        res["cell_type"] = cell_type
        res["contrast"] = contrast_label
        res["tissue_case"] = tissue_b
        res["tissue_control"] = tissue_a
        res["n_cells_celltype"] = subset.n_obs
        res["n_samples_case"] = len(samples_b)
        res["n_samples_control"] = len(samples_a)
        res["significant_spatialfdr"] = res["SpatialFDR"].lt(cfg.alpha)
        res["significant_bh"] = res["PValue_BH"].lt(cfg.alpha)
        results.append(res)

        contrast_prefix = celltype_dir / contrast_label
        res.to_csv(contrast_prefix.with_suffix(".csv"), index=False)
        res.loc[res["significant_spatialfdr"]].to_csv(
            contrast_prefix.with_name(f"{contrast_label}__significant_spatialfdr.csv"), index=False
        )

        if cfg.write_milo_h5ad:
            try:
                milo_adata = mdata["milo"].copy()
                if "annotation_labels" in milo_adata.uns and isinstance(milo_adata.uns["annotation_labels"], pd.Index):
                    milo_adata.uns["annotation_labels"] = milo_adata.uns["annotation_labels"].astype(str).tolist()
                milo_adata.write_h5ad(contrast_prefix.with_suffix(".h5ad"))
            except Exception as exc:  # pragma: no cover - optional artifact path
                LOGGER.warning("Failed to write Milo h5ad for %s / %s: %s", cell_type, contrast_label, exc)

        if cfg.make_plots:
            try:
                plot_milo_graph(
                    milo,
                    mdata,
                    fig_path=contrast_prefix.with_suffix(".png"),
                    title=f"{level_label} | {cell_type}: {tissue_b} vs {tissue_a}",
                    alpha=cfg.alpha,
                )
            except Exception as exc:  # pragma: no cover - plotting should not kill run
                LOGGER.warning("Plotting failed for %s / %s: %s", cell_type, contrast_label, exc)

        pair_summaries.append(
            {
                "level": level_label,
                "celltype_col": celltype_col,
                "cell_type": cell_type,
                "tissue_a": tissue_a,
                "tissue_b": tissue_b,
                "n_samples_a": len(samples_a),
                "n_samples_b": len(samples_b),
                "n_nhoods": int(res.shape[0]),
                "n_sig_spatialfdr": int(res["significant_spatialfdr"].sum()),
                "n_sig_bh": int(res["significant_bh"].sum()),
                "status": "completed",
                "reason": "",
            }
        )

    return results, pair_summaries, plot_df, plot_stats_df


def apply_global_bh(results: pd.DataFrame, alpha: float) -> pd.DataFrame:
    results = results.copy()
    valid = results["PValue"].notna()
    results["PValue_BH_global"] = np.nan
    if valid.any():
        results.loc[valid, "PValue_BH_global"] = multipletests(results.loc[valid, "PValue"], method="fdr_bh")[1]
    results["significant_bh_global"] = results["PValue_BH_global"].lt(alpha)
    return results


def run_level(
    adata: ad.AnnData,
    celltype_col: str,
    level_label: str,
    latent_key: str,
    umap_key: str,
    cfg: MiloRunConfig,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        output_dir / "run_config.json",
        {
            **asdict(cfg),
            "level": level_label,
            "celltype_col": celltype_col,
            "latent_key": latent_key,
            "umap_key": umap_key,
            "output_dir": str(output_dir),
        },
    )
    coverage_df = summarize_celltypes(adata, cfg, celltype_col, level_label)
    coverage_df.to_csv(output_dir / "celltype_coverage_summary.csv", index=False)

    prefilter_celltypes = coverage_df.loc[
        coverage_df["n_cells"] >= cfg.min_cells_per_celltype, "cell_type"
    ].tolist()
    LOGGER.info("Eligible %s cell types after coverage prefilter: %s", level_label, prefilter_celltypes)
    if cfg.max_celltypes_per_level is not None:
        LOGGER.info(
            "%s will process up to %d ready cell types after sample/tissue filtering.",
            level_label,
            cfg.max_celltypes_per_level,
        )

    level_results: list[pd.DataFrame] = []
    subset_summaries: list[dict] = []
    pair_summaries: list[dict] = []
    violin_tables: list[pd.DataFrame] = []
    violin_stats_tables: list[pd.DataFrame] = []
    processed_celltypes: list[str] = []

    for cell_type in prefilter_celltypes:
        if cfg.max_celltypes_per_level is not None and len(processed_celltypes) >= cfg.max_celltypes_per_level:
            break
        subset, subset_info = prepare_celltype_subset(adata, cell_type, celltype_col, level_label, cfg)
        subset_summaries.append(subset_info)
        if subset is None:
            LOGGER.info("Skipping %s | %s: %s", level_label, cell_type, subset_info.get("reason", "unknown"))
            continue
        processed_celltypes.append(cell_type)
        try:
            results, pair_info, plot_df, plot_stats_df = run_milo_for_celltype(
                adata,
                subset,
                cell_type,
                celltype_col,
                level_label,
                latent_key,
                umap_key,
                cfg,
                output_dir,
            )
            level_results.extend(results)
            pair_summaries.extend(pair_info)
            violin_tables.append(plot_df)
            violin_stats_tables.append(plot_stats_df)
        except Exception as exc:  # pragma: no cover - runtime failure path
            LOGGER.exception("Milo failed for %s | %s", level_label, cell_type)
            pair_summaries.append(
                {
                    "level": level_label,
                    "celltype_col": celltype_col,
                    "cell_type": cell_type,
                    "status": "failed",
                    "reason": str(exc),
                }
            )

    LOGGER.info("Processed %s ready cell types: %s", level_label, processed_celltypes)

    subset_df = pd.DataFrame(subset_summaries)
    pair_df = pd.DataFrame(pair_summaries)
    subset_df.to_csv(output_dir / "celltype_subset_summary.csv", index=False)
    pair_df.to_csv(output_dir / "pairwise_run_summary.csv", index=False)
    if violin_tables:
        violin_plot_all = pd.concat(violin_tables, ignore_index=True)
    else:
        violin_plot_all = pd.DataFrame(
            columns=[
                cfg.sample_col,
                cfg.tissue_col,
                "sample",
                "tissue",
                "label_cells",
                "celltype_cells",
                "total_lineage_cells",
                "sample_fraction_pct",
                "analysis_level",
                "level",
                "celltype_col",
                "cell_type",
            ]
        )
    violin_plot_all.to_csv(output_dir / "sample_fraction_per_sample_all.csv", index=False)

    if violin_stats_tables:
        violin_stats_combined = pd.concat(violin_stats_tables, ignore_index=True)
    else:
        violin_stats_combined = pd.DataFrame(
            columns=[
                "analysis_level",
                "level",
                "celltype_col",
                "cell_type",
                "tissue_a",
                "tissue_b",
                "n_samples_a",
                "n_samples_b",
                "median_a",
                "median_b",
                "p_raw",
                "overall_kruskal_p",
                "p_bh",
                "significant_bh",
                "label_bh",
            ]
        )
    violin_stats_combined.to_csv(output_dir / "sample_fraction_pairwise_stats_all.csv", index=False)
    if cfg.make_plots and violin_tables:
        plot_violin_overview(violin_tables, violin_stats_tables, cfg, output_dir, level_label)

    if level_results:
        results_df = pd.concat(level_results, ignore_index=True)
        results_df = apply_global_bh(results_df, alpha=cfg.alpha)
        results_df.to_csv(output_dir / "all_milo_results.csv", index=False)
        results_df.loc[results_df["significant_spatialfdr"]].to_csv(
            output_dir / "all_milo_results__significant_spatialfdr.csv", index=False
        )
        results_df.loc[results_df["significant_bh_global"]].to_csv(
            output_dir / "all_milo_results__significant_global_bh.csv", index=False
        )
        contrast_summary_df = (
            results_df.groupby(["level", "celltype_col", "cell_type", "contrast", "tissue_case", "tissue_control"], dropna=False)
            .agg(
                n_nhoods=("nhood", "size"),
                n_sig_spatialfdr=("significant_spatialfdr", "sum"),
                n_sig_bh=("significant_bh", "sum"),
                n_sig_bh_global=("significant_bh_global", "sum"),
                min_spatialfdr=("SpatialFDR", "min"),
                min_bh_global=("PValue_BH_global", "min"),
            )
            .reset_index()
            .sort_values(["n_sig_spatialfdr", "n_sig_bh_global", "min_spatialfdr"], ascending=[False, False, True])
        )
        contrast_summary_df.to_csv(output_dir / "contrast_level_summary.csv", index=False)
        status = "completed"
    else:
        results_df = pd.DataFrame()
        contrast_summary_df = pd.DataFrame()
        status = "no_results"
        LOGGER.warning("No Milo results were generated for %s.", level_label)

    level_summary = {
        "level": level_label,
        "celltype_col": celltype_col,
        "latent_key": latent_key,
        "umap_key": umap_key,
        "n_celltypes_coverage": int(coverage_df.shape[0]),
        "n_celltypes_prefilter": int(len(prefilter_celltypes)),
        "n_celltypes_eligible": int(len(processed_celltypes)),
        "n_celltypes_ready": int((subset_df.get("status", pd.Series(dtype=str)) == "ready").sum()) if not subset_df.empty else 0,
        "n_celltypes_with_results": int(results_df["cell_type"].nunique()) if not results_df.empty else 0,
        "n_pairwise_rows": int(pair_df.shape[0]),
        "n_completed_pairwise": int((pair_df.get("status", pd.Series(dtype=str)) == "completed").sum()) if not pair_df.empty else 0,
        "n_contrasts_with_spatialfdr_hits": int((contrast_summary_df.get("n_sig_spatialfdr", pd.Series(dtype=int)) > 0).sum()) if not contrast_summary_df.empty else 0,
        "n_result_rows": int(results_df.shape[0]),
        "status": status,
        "output_dir": str(output_dir),
    }
    return coverage_df, subset_df, pair_df, contrast_summary_df, level_summary


def main() -> None:
    setup_logging()
    cfg = parse_args()
    patch_pertpy_rpy2_setup()
    input_h5ad = Path(cfg.input_h5ad)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "run_config.json", asdict(cfg))

    sc.settings.verbosity = 0
    sc.settings.n_jobs = 1

    adata, embedding_info = lightweight_adata_from_h5ad(input_h5ad, cfg)
    write_json(output_dir / "selected_embeddings.json", embedding_info)
    clean_obs_columns(adata, cfg)

    all_level_summaries: list[dict] = []
    all_coverage: list[pd.DataFrame] = []
    all_subset: list[pd.DataFrame] = []
    all_pairs: list[pd.DataFrame] = []
    all_contrast_summaries: list[pd.DataFrame] = []

    latent_key = str(embedding_info["latent_key"])
    umap_key = str(embedding_info["umap_key"] or "X_umap")

    for celltype_col in cfg.celltype_cols:
        level_label = infer_level_label(celltype_col)
        level_dir = output_dir / level_label
        coverage_df, subset_df, pair_df, contrast_summary_df, level_summary = run_level(
            adata=adata,
            celltype_col=celltype_col,
            level_label=level_label,
            latent_key=latent_key,
            umap_key=umap_key,
            cfg=cfg,
            output_dir=level_dir,
        )
        all_level_summaries.append(level_summary)
        if not coverage_df.empty:
            all_coverage.append(coverage_df)
        if not subset_df.empty:
            all_subset.append(subset_df)
        if not pair_df.empty:
            all_pairs.append(pair_df)
        if not contrast_summary_df.empty:
            all_contrast_summaries.append(contrast_summary_df)

    level_summary_df = pd.DataFrame(all_level_summaries)
    level_summary_df.to_csv(output_dir / "level_run_summary.csv", index=False)

    if all_coverage:
        pd.concat(all_coverage, ignore_index=True).to_csv(output_dir / "combined_celltype_coverage_summary.csv", index=False)
    if all_subset:
        pd.concat(all_subset, ignore_index=True).to_csv(output_dir / "combined_celltype_subset_summary.csv", index=False)
    if all_pairs:
        pd.concat(all_pairs, ignore_index=True).to_csv(output_dir / "combined_pairwise_run_summary.csv", index=False)
    if all_contrast_summaries:
        pd.concat(all_contrast_summaries, ignore_index=True).to_csv(output_dir / "combined_contrast_level_summary.csv", index=False)

    LOGGER.info("Finished Milo run. Output directory: %s", output_dir)
    LOGGER.info("Level summary rows: %d", level_summary_df.shape[0])


if __name__ == "__main__":
    main()
