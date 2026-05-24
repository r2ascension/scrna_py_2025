#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Run pertpy Milo tissue differential abundance analysis for epithelial L2/L3 subtypes.

Notes
-----
- Uses the epithelial scANVI export as input and only keeps obs + selected embeddings,
  so neighbourhood DA does not need to load the full expression matrix.
- Runs Milo separately for each epithelial subtype within both L2 and L3 annotations.
- Uses sample as the replicate unit and tissue as the tested condition.
- Writes per-level, per-subtype, and combined contrast summaries under one output root.

Recommended invocation environment
----------------------------------
Use the working pertpy stack and clear inherited path pollution:

    env -u LD_LIBRARY_PATH -u PYTHONPATH PYTHONNOUSERSITE=1 \
    /home/h2048/miniconda3/envs/scarches_stable_pertpy/bin/python \
      /home/h2048/script/py/epithelial_milopy_tissue_celltype_20260419_v1.py
"""

from __future__ import annotations

import argparse
import json
import logging
import re
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


LOGGER = logging.getLogger("epithelial_milo")

LEVEL_SPECS: dict[str, dict[str, object]] = {
    "L2": {
        "label_col": "cell_type_L2",
        "latent_candidates": ("X_scanvi_major", "X_scanvi", "X_scANVI", "X_scvi"),
        "umap_candidates": ("X_umap_major", "X_umap_scanvi", "X_umap"),
    },
    "L3": {
        "label_col": "cell_type_L3",
        "latent_candidates": ("X_scanvi_fine", "X_scanvi", "X_scANVI", "X_scvi"),
        "umap_candidates": ("X_umap_fine", "X_umap_scanvi", "X_umap"),
    },
}


def get_pertpy():
    try:
        import pertpy as pt
    except Exception as exc:  # pragma: no cover - environment/runtime failure path
        raise SystemExit(
            "Failed to import pertpy. Run this script inside /home/h2048/miniconda3/envs/scarches_stable_pertpy "
            f"with PYTHONNOUSERSITE=1. Original error: {exc}"
        ) from exc
    return pt


def patch_pertpy_rpy2_setup() -> None:
    """Patch pertpy Milo for rpy2>=3.6 where activate()/deactivate() are deprecated."""

    get_pertpy()
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
    input_h5ad: str = "/home/h2048/data/py/0417/epithelial_scanvi_v2_8_gpu_rerun/epithelial_scanvi_v2_8_gpu_SELF_for_R.h5ad"
    output_dir: str = "/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun/pertpy_milo"
    sample_col: str = "sample"
    tissue_col: str = "tissue"
    levels: tuple[str, ...] = ("L2", "L3")
    exclude_labels: tuple[str, ...] = ("", "Unknown", "Rejected")
    min_cells_per_group: int = 80
    min_cells_per_sample: int = 5
    min_samples_per_tissue: int = 2
    n_neighbors: int = 30
    nhood_prop: float = 0.1
    alpha: float = 0.10
    random_seed: int = 0
    max_labels_per_level: int | None = None
    make_plots: bool = True
    write_milo_h5ad: bool = True


def parse_args() -> MiloRunConfig:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", default=MiloRunConfig.input_h5ad)
    parser.add_argument("--output-dir", default=MiloRunConfig.output_dir)
    parser.add_argument("--sample-col", default=MiloRunConfig.sample_col)
    parser.add_argument("--tissue-col", default=MiloRunConfig.tissue_col)
    parser.add_argument("--levels", nargs="+", default=list(MiloRunConfig.levels))
    parser.add_argument("--exclude-labels", nargs="*", default=list(MiloRunConfig.exclude_labels))
    parser.add_argument("--min-cells-per-group", type=int, default=MiloRunConfig.min_cells_per_group)
    parser.add_argument("--min-cells-per-sample", type=int, default=MiloRunConfig.min_cells_per_sample)
    parser.add_argument("--min-samples-per-tissue", type=int, default=MiloRunConfig.min_samples_per_tissue)
    parser.add_argument("--n-neighbors", type=int, default=MiloRunConfig.n_neighbors)
    parser.add_argument("--nhood-prop", type=float, default=MiloRunConfig.nhood_prop)
    parser.add_argument("--alpha", type=float, default=MiloRunConfig.alpha)
    parser.add_argument("--random-seed", type=int, default=MiloRunConfig.random_seed)
    parser.add_argument("--max-labels-per-level", type=int, default=None)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--no-write-milo-h5ad", action="store_true")
    args = parser.parse_args()
    return MiloRunConfig(
        input_h5ad=args.input_h5ad,
        output_dir=args.output_dir,
        sample_col=args.sample_col,
        tissue_col=args.tissue_col,
        levels=tuple(args.levels),
        exclude_labels=tuple(args.exclude_labels),
        min_cells_per_group=args.min_cells_per_group,
        min_cells_per_sample=args.min_cells_per_sample,
        min_samples_per_tissue=args.min_samples_per_tissue,
        n_neighbors=args.n_neighbors,
        nhood_prop=args.nhood_prop,
        alpha=args.alpha,
        random_seed=args.random_seed,
        max_labels_per_level=args.max_labels_per_level,
        make_plots=not args.no_plots,
        write_milo_h5ad=not args.no_write_milo_h5ad,
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


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def choose_available_key(available_keys: Iterable[str], candidates: Iterable[str], kind: str) -> str:
    available = set(available_keys)
    for candidate in candidates:
        if candidate in available:
            return candidate
    raise KeyError(f"Missing {kind}; tried candidates: {list(candidates)}")


def normalize_levels(levels: Iterable[str]) -> tuple[str, ...]:
    normalized = tuple(dict.fromkeys(level.upper() for level in levels))
    invalid = [level for level in normalized if level not in LEVEL_SPECS]
    if invalid:
        raise ValueError(f"Unsupported levels: {invalid}; supported: {sorted(LEVEL_SPECS)}")
    return normalized


def collect_required_obsm_keys(levels: Iterable[str]) -> tuple[str, ...]:
    keys: list[str] = []
    for level in normalize_levels(levels):
        spec = LEVEL_SPECS[level]
        keys.extend(spec["latent_candidates"])
        keys.extend(spec["umap_candidates"])
    return tuple(dict.fromkeys(keys))


def lightweight_adata_from_h5ad(input_h5ad: Path, cfg: MiloRunConfig) -> ad.AnnData:
    LOGGER.info("Loading lightweight AnnData from %s", input_h5ad)
    src = ad.read_h5ad(input_h5ad, backed="r")
    try:
        levels = normalize_levels(cfg.levels)
        required_obs_cols = [cfg.sample_col, cfg.tissue_col] + [LEVEL_SPECS[level]["label_col"] for level in levels]
        missing_cols = [col for col in required_obs_cols if col not in src.obs.columns]
        if missing_cols:
            raise KeyError(f"Missing required obs columns: {missing_cols}")

        needed_obsm_keys = collect_required_obsm_keys(levels)
        found_embedding = any(key in src.obsm.keys() for key in needed_obsm_keys)
        if not found_embedding:
            raise KeyError(f"None of the expected epithelial embeddings were found in obsm: {needed_obsm_keys}")

        obs = src.obs.copy()
        dummy_x = sparse.csr_matrix((src.n_obs, 1), dtype=np.float32)
        adata = ad.AnnData(X=dummy_x, obs=obs, var=pd.DataFrame(index=["placeholder_feature"]))
        for key in needed_obsm_keys:
            if key in src.obsm.keys():
                adata.obsm[key] = np.asarray(src.obsm[key])
        LOGGER.info("Lightweight AnnData shape: %s", adata.shape)
        return adata
    finally:
        if getattr(src, "file", None) is not None:
            src.file.close()


def clean_obs_columns(adata: ad.AnnData, cfg: MiloRunConfig) -> None:
    levels = normalize_levels(cfg.levels)
    label_cols = [LEVEL_SPECS[level]["label_col"] for level in levels]
    for col in [cfg.sample_col, cfg.tissue_col, *label_cols]:
        adata.obs[col] = adata.obs[col].astype(str).str.strip()
        adata.obs[col] = adata.obs[col].replace({"nan": "", "None": "", "NA": ""})


def summarize_level_groups(adata: ad.AnnData, label_col: str, cfg: MiloRunConfig) -> pd.DataFrame:
    obs = adata.obs[[label_col, cfg.sample_col, cfg.tissue_col]].copy()
    obs = obs.loc[
        ~obs[label_col].isin(cfg.exclude_labels)
        & (obs[cfg.sample_col] != "")
        & (obs[cfg.tissue_col] != "")
    ].copy()
    rows: list[dict[str, object]] = []
    for label, sub in obs.groupby(label_col, sort=False):
        tissue_sample_counts = sub.groupby(cfg.tissue_col)[cfg.sample_col].nunique().to_dict()
        tissue_cell_counts = sub[cfg.tissue_col].value_counts().to_dict()
        rows.append(
            {
                "cell_type": label,
                "n_cells": int(sub.shape[0]),
                "n_samples": int(sub[cfg.sample_col].nunique()),
                "n_tissues": int(sub[cfg.tissue_col].nunique()),
                "tissue_sample_counts": json.dumps(tissue_sample_counts, ensure_ascii=False, sort_keys=True),
                "tissue_cell_counts": json.dumps(tissue_cell_counts, ensure_ascii=False, sort_keys=True),
            }
        )
    if not rows:
        return pd.DataFrame(columns=["cell_type", "n_cells", "n_samples", "n_tissues", "tissue_sample_counts", "tissue_cell_counts"])
    return pd.DataFrame(rows).sort_values(["n_cells", "n_samples"], ascending=[False, False]).reset_index(drop=True)


def prepare_level_subset(
    adata: ad.AnnData,
    label_col: str,
    label_value: str,
    cfg: MiloRunConfig,
) -> tuple[ad.AnnData | None, dict[str, object]]:
    subset = adata[adata.obs[label_col] == label_value].copy()
    info: dict[str, object] = {
        "label_col": label_col,
        "cell_type": label_value,
        "n_cells_initial": int(subset.n_obs),
    }
    if subset.n_obs < cfg.min_cells_per_group:
        info["status"] = "skipped"
        info["reason"] = f"fewer than {cfg.min_cells_per_group} cells"
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

    if subset.n_obs < cfg.min_cells_per_group:
        info["status"] = "skipped"
        info["reason"] = f"fewer than {cfg.min_cells_per_group} cells after sample filter"
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


def plot_milo_graph(milo, mdata, fig_path: Path, title: str, alpha: float) -> None:
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
    level: str,
    label_value: str,
    cfg: MiloRunConfig,
) -> pd.DataFrame:
    sample_meta = subset.obs[[cfg.sample_col, cfg.tissue_col]].drop_duplicates().copy()
    all_counts = full_adata.obs.groupby(cfg.sample_col).size().rename("total_lineage_cells")
    label_counts = subset.obs.groupby(cfg.sample_col).size().rename("label_cells")
    plot_df = sample_meta.copy()
    plot_df["label_cells"] = plot_df[cfg.sample_col].map(label_counts).fillna(0).astype(int)
    plot_df["total_lineage_cells"] = plot_df[cfg.sample_col].map(all_counts).fillna(0).astype(int)
    plot_df = plot_df.loc[plot_df["total_lineage_cells"] > 0].copy()
    plot_df["sample_fraction_pct"] = 100.0 * plot_df["label_cells"] / plot_df["total_lineage_cells"]
    plot_df["analysis_level"] = level
    plot_df["cell_type"] = label_value
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


def compute_violin_pairwise_stats(
    plot_df: pd.DataFrame,
    cfg: MiloRunConfig,
) -> tuple[pd.DataFrame, float | None]:
    tissue_order = plot_df[cfg.tissue_col].astype(str).drop_duplicates().tolist()
    rows: list[dict[str, object]] = []
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
        legend=show_legend,
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
    level: str,
    label_value: str,
    cfg: MiloRunConfig,
    label_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    plot_df = compute_sample_fraction_table(full_adata, subset, level, label_value, cfg)
    stats_df, _ = compute_violin_pairwise_stats(plot_df, cfg)
    plot_df.to_csv(label_dir / "sample_fraction_per_sample.csv", index=False)
    stats_df.to_csv(label_dir / "sample_fraction_pairwise_stats.csv", index=False)

    if cfg.make_plots:
        fig, ax = plt.subplots(figsize=(max(5, 1.6 * plot_df[cfg.tissue_col].nunique()), 5.2))
        draw_sample_fraction_violin(
            ax,
            plot_df,
            stats_df,
            cfg,
            title=f"{level} {label_value} sample-level abundance across tissues",
            ylabel="Subtype fraction per sample (%)",
        )
        fig.tight_layout()
        fig.savefig(label_dir / "sample_fraction_violin.png", dpi=220, bbox_inches="tight")
        plt.close(fig)

    return plot_df, stats_df


def plot_violin_overview(
    level: str,
    plot_tables: list[pd.DataFrame],
    stats_tables: list[pd.DataFrame],
    cfg: MiloRunConfig,
    output_dir: Path,
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

    fig.suptitle(f"Epithelial {level} subtype abundance across tissues", fontsize=14, y=0.995)
    fig.tight_layout()
    fig.savefig(output_dir / f"{level.lower()}_sample_fraction_violin_overview.png", dpi=220, bbox_inches="tight")
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


def run_milo_for_label(
    full_adata: ad.AnnData,
    subset: ad.AnnData,
    level: str,
    label_value: str,
    latent_key: str,
    umap_key: str,
    cfg: MiloRunConfig,
    level_dir: Path,
) -> tuple[list[pd.DataFrame], list[dict[str, object]], pd.DataFrame, pd.DataFrame]:
    LOGGER.info("Running Milo for %s %s (%d cells)", level, label_value, subset.n_obs)
    pt = get_pertpy()
    label_slug = slugify(label_value)
    label_dir = level_dir / label_slug
    label_dir.mkdir(parents=True, exist_ok=True)

    ensure_neighbors_and_umap(subset, latent_key, umap_key, cfg)
    milo = pt.tl.Milo()
    mdata = milo.load(subset)
    milo.make_nhoods(mdata["rna"], prop=cfg.nhood_prop, seed=cfg.random_seed)
    mdata = milo.count_nhoods(mdata, sample_col=cfg.sample_col)
    milo.annotate_nhoods(mdata, anno_col=cfg.tissue_col)
    milo.build_nhood_graph(mdata, basis=umap_key)
    plot_df, plot_stats_df = save_violin_outputs(full_adata, subset, level, label_value, cfg, label_dir)

    pair_summaries: list[dict[str, object]] = []
    results: list[pd.DataFrame] = []
    tissue_pairs = build_tissue_pairs(subset, cfg)
    if not tissue_pairs:
        pair_summaries.append({
            "analysis_level": level,
            "cell_type": label_value,
            "status": "skipped",
            "reason": "no eligible tissue pairs",
        })
        return results, pair_summaries, plot_df, plot_stats_df

    sample_meta = subset.obs[[cfg.sample_col, cfg.tissue_col, "tissue_safe"]].drop_duplicates().copy()
    tissue_map = dict(zip(sample_meta[cfg.tissue_col], sample_meta["tissue_safe"], strict=False))

    for tissue_a, tissue_b in tissue_pairs:
        samples_a = sorted(sample_meta.loc[sample_meta[cfg.tissue_col] == tissue_a, cfg.sample_col].tolist())
        samples_b = sorted(sample_meta.loc[sample_meta[cfg.tissue_col] == tissue_b, cfg.sample_col].tolist())
        if len(samples_a) < cfg.min_samples_per_tissue or len(samples_b) < cfg.min_samples_per_tissue:
            pair_summaries.append(
                {
                    "analysis_level": level,
                    "cell_type": label_value,
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
            "  %s %s contrast %s | %s (%d samples) vs %s (%d samples)",
            level,
            label_value,
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
        res["analysis_level"] = level
        res["cell_type"] = label_value
        res["contrast"] = contrast_label
        res["tissue_case"] = tissue_b
        res["tissue_control"] = tissue_a
        res["n_cells_celltype"] = subset.n_obs
        res["n_samples_case"] = len(samples_b)
        res["n_samples_control"] = len(samples_a)
        res["significant_spatialfdr"] = res["SpatialFDR"].lt(cfg.alpha)
        res["significant_bh"] = res["PValue_BH"].lt(cfg.alpha)
        results.append(res)

        contrast_prefix = label_dir / contrast_label
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
                LOGGER.warning("Failed to write Milo h5ad for %s %s / %s: %s", level, label_value, contrast_label, exc)

        if cfg.make_plots:
            try:
                plot_milo_graph(
                    milo,
                    mdata,
                    fig_path=contrast_prefix.with_suffix(".png"),
                    title=f"{level} {label_value}: {tissue_b} vs {tissue_a}",
                    alpha=cfg.alpha,
                )
            except Exception as exc:  # pragma: no cover - plotting should not kill run
                LOGGER.warning("Plotting failed for %s %s / %s: %s", level, label_value, contrast_label, exc)

        pair_summaries.append(
            {
                "analysis_level": level,
                "cell_type": label_value,
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
    level: str,
    cfg: MiloRunConfig,
    output_dir: Path,
) -> tuple[pd.DataFrame | None, pd.DataFrame, pd.DataFrame]:
    spec = LEVEL_SPECS[level]
    label_col = spec["label_col"]
    latent_key = choose_available_key(adata.obsm.keys(), spec["latent_candidates"], kind=f"{level} latent embedding")
    umap_key = choose_available_key(adata.obsm.keys(), spec["umap_candidates"], kind=f"{level} UMAP embedding")
    level_dir = output_dir / level
    level_dir.mkdir(parents=True, exist_ok=True)

    level_cfg_payload = {
        "analysis_level": level,
        "label_col": label_col,
        "latent_key": latent_key,
        "umap_key": umap_key,
    }
    write_json(level_dir / "run_config.json", level_cfg_payload)

    label_summary = summarize_level_groups(adata, label_col, cfg)
    label_summary.to_csv(level_dir / "celltype_coverage_summary.csv", index=False)

    eligible = label_summary.loc[label_summary["n_cells"] >= cfg.min_cells_per_group, "cell_type"].tolist()
    if cfg.max_labels_per_level is not None:
        eligible = eligible[: cfg.max_labels_per_level]
    LOGGER.info("%s eligible labels (pre-filter): %s", level, eligible)

    all_results: list[pd.DataFrame] = []
    subset_summaries: list[dict[str, object]] = []
    pair_summaries: list[dict[str, object]] = []
    violin_tables: list[pd.DataFrame] = []
    violin_stats_tables: list[pd.DataFrame] = []

    for label_value in eligible:
        subset, subset_info = prepare_level_subset(adata, label_col, label_value, cfg)
        subset_info["analysis_level"] = level
        subset_summaries.append(subset_info)
        if subset is None:
            LOGGER.info("Skipping %s %s: %s", level, label_value, subset_info.get("reason", "unknown"))
            continue
        try:
            results, pair_info, plot_df, plot_stats_df = run_milo_for_label(
                adata,
                subset,
                level,
                label_value,
                latent_key,
                umap_key,
                cfg,
                level_dir,
            )
            all_results.extend(results)
            pair_summaries.extend(pair_info)
            violin_tables.append(plot_df)
            violin_stats_tables.append(plot_stats_df)
        except Exception as exc:  # pragma: no cover - runtime failure path
            LOGGER.exception("Milo failed for %s %s", level, label_value)
            pair_summaries.append(
                {
                    "analysis_level": level,
                    "cell_type": label_value,
                    "status": "failed",
                    "reason": str(exc),
                }
            )

    subset_summary_df = pd.DataFrame(subset_summaries)
    pair_summary_df = pd.DataFrame(pair_summaries)
    subset_summary_df.to_csv(level_dir / "celltype_subset_summary.csv", index=False)
    pair_summary_df.to_csv(level_dir / "pairwise_run_summary.csv", index=False)

    if violin_tables:
        pd.concat(violin_tables, ignore_index=True).to_csv(level_dir / "sample_fraction_per_sample_all.csv", index=False)
    if violin_stats_tables:
        violin_stats_combined = pd.concat(violin_stats_tables, ignore_index=True)
        violin_stats_combined.to_csv(level_dir / "sample_fraction_pairwise_stats_all.csv", index=False)
        if cfg.make_plots:
            plot_violin_overview(level, violin_tables, violin_stats_tables, cfg, level_dir)

    if not all_results:
        LOGGER.warning("No Milo results were generated for %s.", level)
        return None, subset_summary_df, pair_summary_df

    combined = pd.concat(all_results, ignore_index=True)
    combined = apply_global_bh(combined, alpha=cfg.alpha)
    combined.to_csv(level_dir / "all_milo_results.csv", index=False)
    combined.loc[combined["significant_spatialfdr"]].to_csv(
        level_dir / "all_milo_results__significant_spatialfdr.csv", index=False
    )
    combined.loc[combined["significant_bh_global"]].to_csv(
        level_dir / "all_milo_results__significant_global_bh.csv", index=False
    )

    summary = (
        combined.groupby(["analysis_level", "cell_type", "contrast", "tissue_case", "tissue_control"], dropna=False)
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
    summary.to_csv(level_dir / "contrast_level_summary.csv", index=False)
    return combined, subset_summary_df, pair_summary_df


def main() -> None:
    setup_logging()
    cfg = parse_args()
    cfg.levels = normalize_levels(cfg.levels)
    patch_pertpy_rpy2_setup()

    input_h5ad = Path(cfg.input_h5ad)
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "run_config.json", asdict(cfg))

    sc.settings.verbosity = 0
    sc.settings.n_jobs = 1

    adata = lightweight_adata_from_h5ad(input_h5ad, cfg)
    clean_obs_columns(adata, cfg)

    root_level_summaries: list[pd.DataFrame] = []
    root_subset_summaries: list[pd.DataFrame] = []
    root_pair_summaries: list[pd.DataFrame] = []
    root_results: list[pd.DataFrame] = []

    for level in cfg.levels:
        level_results, subset_summary_df, pair_summary_df = run_level(adata, level, cfg, output_dir)
        if level_results is not None:
            root_results.append(level_results)
            level_summary_path = output_dir / level / "contrast_level_summary.csv"
            if level_summary_path.exists():
                root_level_summaries.append(pd.read_csv(level_summary_path))
        root_subset_summaries.append(subset_summary_df)
        root_pair_summaries.append(pair_summary_df)

    if root_subset_summaries:
        pd.concat(root_subset_summaries, ignore_index=True).to_csv(output_dir / "all_levels_subset_summary.csv", index=False)
    if root_pair_summaries:
        pd.concat(root_pair_summaries, ignore_index=True).to_csv(output_dir / "all_levels_pairwise_run_summary.csv", index=False)
    if root_level_summaries:
        pd.concat(root_level_summaries, ignore_index=True).to_csv(output_dir / "all_levels_contrast_summary.csv", index=False)
    if root_results:
        pd.concat(root_results, ignore_index=True).to_csv(output_dir / "all_levels_milo_results.csv", index=False)
        LOGGER.info("Finished epithelial Milo run. Output directory: %s", output_dir)
    else:
        LOGGER.warning("Finished epithelial Milo run, but no eligible contrasts produced Milo results.")


if __name__ == "__main__":
    main()
