#!/usr/bin/env python3
"""Brief visualization for the B-cell scPoli smoke-test outputs."""

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

OUTDIR = Path("/home/h2048/temp/bcell_scpoli_smoke_20260511")
FIGDIR = OUTDIR / "figures_20260512"
H5AD = OUTDIR / "bcell_latest_scpoli_smoke.h5ad"
CONFUSION = OUTDIR / "bcell_latest_scpoli_confusion.tsv"
PER_LABEL = OUTDIR / "bcell_latest_scpoli_per_label.tsv"
SUMMARY = OUTDIR / "bcell_latest_scpoli_summary.json"

LABEL_KEY = "cell_type_expert"
PRED_KEY = "scpoli_pred"
UNCERT_KEY = "scpoli_uncert"
LATENT_KEY = "X_scPoli_smoke"


def short_label(x: str) -> str:
    replacements = {
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
    return replacements.get(str(x), str(x).replace("_", " "))


def scatter_panel(ax, emb: np.ndarray, labels: pd.Series, title: str, palette: dict[str, tuple[float, float, float]]) -> None:
    plot_df = pd.DataFrame({"x": emb[:, 0], "y": emb[:, 1], "label": labels.astype(str).map(short_label)})
    sns.scatterplot(
        data=plot_df,
        x="x",
        y="y",
        hue="label",
        palette=palette,
        s=16,
        linewidth=0,
        alpha=0.78,
        ax=ax,
        rasterized=True,
    )
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("scPoli latent PC1")
    ax.set_ylabel("scPoli latent PC2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8, title=None)
    ax.grid(False)


def main() -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk", font_scale=0.72)

    adata = ad.read_h5ad(H5AD)
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    confusion = pd.read_csv(CONFUSION, sep="\t", index_col=0)
    per_label = pd.read_csv(PER_LABEL, sep="\t", index_col=0)

    latent = np.asarray(adata.obsm[LATENT_KEY])
    pca = PCA(n_components=2, random_state=0)
    latent_pca = pca.fit_transform(latent)

    labels_short = adata.obs[LABEL_KEY].astype(str).map(short_label)
    all_levels = sorted(pd.unique(pd.concat([labels_short, adata.obs[PRED_KEY].astype(str).map(short_label)])))
    colors = sns.color_palette("tab10", n_colors=max(10, len(all_levels)))
    palette = dict(zip(all_levels, colors[: len(all_levels)]))

    # Overview figure
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.05, 1.0, 0.85], width_ratios=[1.05, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, 0])
    ax6 = fig.add_subplot(gs[2, 1])

    scatter_panel(ax1, latent_pca, adata.obs[LABEL_KEY], "A  scPoli latent PCA — expert label", palette)
    scatter_panel(ax2, latent_pca, adata.obs[PRED_KEY], "B  scPoli latent PCA — predicted label", palette)

    conf_norm = confusion.div(confusion.sum(axis=1).replace(0, np.nan), axis=0)
    conf_norm.index = [short_label(x) for x in conf_norm.index]
    conf_norm.columns = [short_label(x) for x in conf_norm.columns]
    sns.heatmap(
        conf_norm,
        cmap="Blues",
        vmin=0,
        vmax=1,
        ax=ax3,
        cbar_kws={"label": "row fraction"},
        linewidths=0.2,
        linecolor="white",
    )
    ax3.set_title("C  Prediction confusion matrix", loc="left", fontweight="bold")
    ax3.set_xlabel("Predicted")
    ax3.set_ylabel("Expert label")
    ax3.tick_params(axis="x", labelrotation=45)
    ax3.tick_params(axis="y", labelrotation=0)

    per_label_plot = per_label.copy()
    per_label_plot.index = [short_label(x) for x in per_label_plot.index]
    per_label_plot = per_label_plot.sort_values("accuracy")
    sns.barplot(
        x=per_label_plot["accuracy"],
        y=per_label_plot.index,
        hue=per_label_plot.index,
        dodge=False,
        palette=[palette.get(i, (0.3, 0.3, 0.3)) for i in per_label_plot.index],
        ax=ax4,
        legend=False,
    )
    for i, (acc, n) in enumerate(zip(per_label_plot["accuracy"], per_label_plot["n"])):
        ax4.text(min(float(acc) + 0.02, 0.98), i, f"{acc:.2f} (n={int(n)})", va="center", fontsize=8)
    ax4.set_xlim(0, 1)
    ax4.set_title("D  Per-label self-classification accuracy", loc="left", fontweight="bold")
    ax4.set_xlabel("Accuracy")
    ax4.set_ylabel("")

    uncert_df = adata.obs[[LABEL_KEY, UNCERT_KEY]].copy()
    uncert_df[LABEL_KEY] = uncert_df[LABEL_KEY].astype(str).map(short_label)
    order = uncert_df.groupby(LABEL_KEY, observed=True)[UNCERT_KEY].median().sort_values().index.tolist()
    sns.boxplot(
        data=uncert_df,
        x=UNCERT_KEY,
        y=LABEL_KEY,
        order=order,
        color="#9ecae1",
        fliersize=1.5,
        linewidth=0.8,
        ax=ax5,
    )
    ax5.set_title("E  Classification uncertainty", loc="left", fontweight="bold")
    ax5.set_xlabel("Scaled uncertainty")
    ax5.set_ylabel("")
    ax5.set_xlim(-0.02, 1.02)

    logs = summary.get("trainer_logs", {})
    epochs = np.arange(1, len(logs.get("epoch_cvae_loss", [])) + 1)
    if len(epochs) > 0:
        ax6.plot(epochs, logs.get("epoch_cvae_loss", []), marker="o", label="train CVAE")
        ax6.plot(epochs, logs.get("val_cvae_loss", []), marker="o", label="val CVAE")
        ax6.legend(frameon=False)
    ax6.set_title("F  Training loss trend", loc="left", fontweight="bold")
    ax6.set_xlabel("Epoch")
    ax6.set_ylabel("CVAE loss")

    fig.suptitle(
        "B-cell scPoli smoke-test visualization\n"
        f"subset={summary['shape'][0]} cells × {summary['shape'][1]} genes; "
        f"accuracy={summary['overall_self_classification_accuracy']:.3f}; "
        f"median uncertainty={summary['uncertainty_median']:.3f}",
        fontsize=15,
        fontweight="bold",
    )
    overview_png = FIGDIR / "bcell_scpoli_smoke_overview.png"
    overview_pdf = FIGDIR / "bcell_scpoli_smoke_overview.pdf"
    fig.savefig(overview_png, dpi=220, bbox_inches="tight")
    fig.savefig(overview_pdf, bbox_inches="tight")
    plt.close(fig)

    # Focused split figures for quick embedding in notes.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    scatter_panel(axes[0], latent_pca, adata.obs[LABEL_KEY], "scPoli latent PCA — expert label", palette)
    scatter_panel(axes[1], latent_pca, adata.obs[PRED_KEY], "scPoli latent PCA — predicted label", palette)
    latent_png = FIGDIR / "bcell_scpoli_latent_pca_label_vs_pred.png"
    fig.savefig(latent_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.5, 7), constrained_layout=True)
    sns.heatmap(conf_norm, cmap="Blues", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "row fraction"})
    ax.set_title("B-cell scPoli row-normalized confusion matrix", fontweight="bold")
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Expert label")
    ax.tick_params(axis="x", labelrotation=45)
    heatmap_png = FIGDIR / "bcell_scpoli_confusion_row_fraction.png"
    fig.savefig(heatmap_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    report = FIGDIR / "README_bcell_scpoli_viz.md"
    report.write_text(
        "# B-cell scPoli smoke-test quick visualization\n\n"
        f"- Input h5ad: `{H5AD}`\n"
        f"- Cells × genes: `{summary['shape'][0]} × {summary['shape'][1]}`\n"
        f"- Overall self-classification accuracy: `{summary['overall_self_classification_accuracy']:.3f}`\n"
        f"- Median uncertainty: `{summary['uncertainty_median']:.3f}`\n"
        f"- Overview PNG: `{overview_png}`\n"
        f"- Overview PDF: `{overview_pdf}`\n"
        f"- Latent label/pred PNG: `{latent_png}`\n"
        f"- Confusion heatmap PNG: `{heatmap_png}`\n\n"
        "Visualization choices followed Bizard-style compact biomedical QC: dimensionality-reduction scatter, heatmap, and ranking bar/box summaries.\n",
        encoding="utf-8",
    )

    print("Wrote figures:")
    for p in [overview_png, overview_pdf, latent_png, heatmap_png, report]:
        print(p)


if __name__ == "__main__":
    main()
