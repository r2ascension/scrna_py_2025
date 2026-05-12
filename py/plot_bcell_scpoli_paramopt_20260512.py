#!/usr/bin/env python3
"""Visualize B-cell scPoli parameter optimization results."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

COMBINED = Path("/home/h2048/temp/bcell_scpoli_paramopt_combined_20260512")
FAST = Path("/home/h2048/temp/bcell_scpoli_paramopt_fast_20260512")
TRADEOFF = Path("/home/h2048/temp/bcell_scpoli_paramopt_tradeoff_20260512")
BEST_NAME = "batch_mmd_beta0p2_h128_lat16_eta2"
BEST_ROOT = TRADEOFF
BEST_DIR = BEST_ROOT / BEST_NAME
BEST_H5AD = BEST_DIR / "latent_predictions_only.h5ad"
BEST_CONF = BEST_DIR / "holdout_confusion.tsv"
BEST_METRICS = BEST_DIR / "metrics.json"
FIGDIR = COMBINED / "figures"

LABEL_KEY = "cell_type_expert"
PRED_KEY = f"{BEST_NAME}_pred"
UNCERT_KEY = f"{BEST_NAME}_uncert"
LATENT_KEY = f"X_scPoli_{BEST_NAME}"
CONDITION_KEY = "batch"


def short_label(x: str) -> str:
    repl = {
        "GC_B_Dark_Zone_Centroblast_Cycling": "GC DZ cycling",
        "GC_B_Light_Zone_Centrocyte": "GC LZ centrocyte",
        "GC_B_Transitional": "GC transitional",
        "Atypical_Memory_B": "Atypical memory",
        "IGHEplus_Atypical_Memory_B": "IGHE+ atypical",
        "Plasma_IgA": "Plasma IgA",
        "Plasma_IgG": "Plasma IgG",
        "Memory_B": "Memory B",
        "Naive_B": "Naive B",
    }
    return repl.get(str(x), str(x).replace("_", " "))


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk", font_scale=0.75)

    results = pd.read_csv(COMBINED / "combined_results.tsv", sep="\t")
    metrics = json.loads(BEST_METRICS.read_text(encoding="utf-8"))
    adata = ad.read_h5ad(BEST_H5AD)
    conf = pd.read_csv(BEST_CONF, sep="\t", index_col=0)

    # Compact names for plotting.
    name_map = {
        "batch_mmd_beta0p2_h128_lat16_eta2": "batch+MMD β0.2\nh128 lat16 η2",
        "mmd_study_beta0p2_h128_lat16_eta1p5": "study+MMD β0.2\nh128 lat16 η1.5",
        "batch_condition_h128_lat16_eta2": "batch only\nh128 lat16 η2",
        "longer_baseline_study_h128_lat16_eta1p5": "study baseline\nh128 lat16 η1.5",
        "mmd_study_beta0p5_h256_128_lat24_eta2": "study+MMD β0.5\nh256-128 lat24",
        "capacity_study_h256_128_lat24_eta2": "study capacity\nh256-128 lat24",
        "study_mmd_beta0p2_stratified_h128_lat16_eta2": "study+MMD+strat\nh128 lat16 η2",
    }
    results["pretty"] = results["name"].map(name_map).fillna(results["name"])
    results = results.sort_values("composite_score", ascending=False)

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    ax = axes[0, 0]
    melted = results.melt(
        id_vars=["pretty", "name"],
        value_vars=["holdout_accuracy", "holdout_macro_f1", "holdout_balanced_accuracy"],
        var_name="metric",
        value_name="value",
    )
    sns.barplot(data=melted, x="pretty", y="value", hue="metric", ax=ax)
    ax.set_title("A  Holdout classification metrics", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("score")
    ax.set_ylim(0, 0.75)
    ax.tick_params(axis="x", labelrotation=35)
    ax.legend(frameon=False, fontsize=8, title=None)

    ax = axes[0, 1]
    mix_df = results.melt(
        id_vars=["pretty", "name"],
        value_vars=["condition_silhouette", "within_label_condition_asw_abs"],
        var_name="metric",
        value_name="value",
    )
    sns.barplot(data=mix_df, x="pretty", y="value", hue="metric", ax=ax)
    ax.set_title("B  Batch/condition separation (lower is better)", loc="left", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("ASW / silhouette")
    ax.tick_params(axis="x", labelrotation=35)
    ax.legend(frameon=False, fontsize=8, title=None)

    latent = np.asarray(adata.obsm[LATENT_KEY])
    emb = PCA(n_components=2, random_state=0).fit_transform(latent)
    plot_df = pd.DataFrame(
        {
            "PC1": emb[:, 0],
            "PC2": emb[:, 1],
            "label": adata.obs[LABEL_KEY].astype(str).map(short_label).values,
            "pred": adata.obs[PRED_KEY].astype(str).map(short_label).values,
            "split": adata.obs["scpoli_split"].astype(str).values,
            "uncert": adata.obs[UNCERT_KEY].astype(float).values,
        },
        index=adata.obs_names,
    )
    label_levels = sorted(plot_df["label"].unique())
    palette = dict(zip(label_levels, sns.color_palette("tab10", len(label_levels))))
    ax = axes[1, 0]
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="label", style="split", palette=palette, s=18, linewidth=0, alpha=0.8, ax=ax)
    ax.set_title("C  Best latent PCA by expert label", loc="left", fontweight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, title=None)

    ax = axes[1, 1]
    sns.scatterplot(data=plot_df, x="PC1", y="PC2", hue="pred", style="split", palette=palette, s=18, linewidth=0, alpha=0.8, ax=ax)
    ax.set_title("D  Best latent PCA by scPoli prediction", loc="left", fontweight="bold")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7, title=None)

    fig.suptitle(
        "B-cell scPoli parameter optimization\n"
        f"Best: {BEST_NAME}; holdout acc={metrics['holdout_accuracy']:.3f}, "
        f"macro-F1={metrics['holdout_macro_f1']:.3f}, within-label batch ASW={metrics['within_label_condition_asw_abs']:.3f}",
        fontweight="bold",
        fontsize=15,
    )
    fig.savefig(FIGDIR / "bcell_scpoli_paramopt_overview.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGDIR / "bcell_scpoli_paramopt_overview.pdf", bbox_inches="tight")
    plt.close(fig)

    conf_norm = conf.div(conf.sum(axis=1).replace(0, np.nan), axis=0)
    conf_norm.index = [short_label(x) for x in conf_norm.index]
    conf_norm.columns = [short_label(x) for x in conf_norm.columns]
    fig, ax = plt.subplots(figsize=(9, 7.5), constrained_layout=True)
    sns.heatmap(conf_norm, cmap="Blues", vmin=0, vmax=1, linewidths=0.2, linecolor="white", ax=ax, cbar_kws={"label": "row fraction"})
    ax.set_title("Best scPoli holdout confusion matrix", fontweight="bold")
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Expert label")
    ax.tick_params(axis="x", labelrotation=45)
    ax.tick_params(axis="y", labelrotation=0)
    fig.savefig(FIGDIR / "bcell_scpoli_best_holdout_confusion.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # Markdown summary.
    rec = metrics["config"]
    readme = FIGDIR / "README_paramopt_viz.md"
    readme.write_text(
        "# B-cell scPoli parameter optimization\n\n"
        f"Best config: `{BEST_NAME}`\n\n"
        f"- Holdout accuracy: `{metrics['holdout_accuracy']:.3f}`\n"
        f"- Holdout macro-F1: `{metrics['holdout_macro_f1']:.3f}`\n"
        f"- Holdout balanced accuracy: `{metrics['holdout_balanced_accuracy']:.3f}`\n"
        f"- Label silhouette: `{metrics['label_silhouette']:.3f}`\n"
        f"- Condition silhouette: `{metrics['condition_silhouette']:.3f}`\n"
        f"- Within-label condition ASW abs: `{metrics['within_label_condition_asw_abs']:.3f}`\n\n"
        "Recommended parameters:\n\n"
        f"```json\n{json.dumps(rec, indent=2, ensure_ascii=False)}\n```\n\n"
        f"Combined results: `{COMBINED / 'combined_results.tsv'}`\n"
        f"Best model dir: `{BEST_DIR / 'model'}`\n"
        f"Best latent/pred h5ad: `{BEST_H5AD}`\n",
        encoding="utf-8",
    )

    print("Wrote:")
    for p in [FIGDIR / "bcell_scpoli_paramopt_overview.png", FIGDIR / "bcell_scpoli_paramopt_overview.pdf", FIGDIR / "bcell_scpoli_best_holdout_confusion.png", readme]:
        print(p)


if __name__ == "__main__":
    main()
