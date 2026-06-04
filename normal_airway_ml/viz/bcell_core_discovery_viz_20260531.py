#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

DPI = 300
CELLTYPE_ORDER = ["Memory_B", "Naive_B", "Plasma", "GC_B"]
SITE_ORDER = ["nasal", "sinus", "bronchus", "lung_parenchyma"]
CELLTYPE_COLORS = {
    "Memory_B": "#3b6fb6",
    "Naive_B": "#2c9a57",
    "Plasma": "#d17b0f",
    "GC_B": "#8a5fbf",
}
STATUS_COLORS = {
    "ok": "#4daf4a",
    "failed": "#d73027",
    "skipped": "#9e9e9e",
}
CONTRAST_MARKERS = {
    "pairwise": "o",
    "one_vs_rest": "s",
}
STANDARD_GENE_FILTER_RULES = [
    ("mitochondrial", r"^(MT-|MTRNR)"),
    ("ribosomal", r"^(RPS|RPL|MRPS|MRPL)"),
    ("ribosomal_pseudogene", r"^(RPS|RPL|MRPS|MRPL)[0-9]+P[0-9]+$"),
    ("ensembl_id", r"^ENSG[0-9]"),
    ("technical_ncrna", r"^(MALAT1|NEAT1)$"),
    ("locus_like", r"^(LINC|AC[0-9]+|AL[0-9]+|AP[0-9]+|RP11-|RP[0-9]+-)"),
]
GENE_FILTER_TITLE_NOTE = "standard MT/RP/ENSG/locus-like filter applied"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate PNG/PDF visualizations for the 2026-05-31 B-cell core discovery run.")
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path("/home/h2048/data/py/20260531/normal_airway_ml_bcell_full_run_20260531"),
        help="Root directory of the completed full run.",
    )
    parser.add_argument(
        "--top-site-genes-per-group",
        type=int,
        default=3,
        help="Top genes to retain per (cell_type, enriched_site) for the site-specific dot plot.",
    )
    return parser.parse_args()


def ensure_dirs(paths: Iterable[Path]) -> None:
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, stem: Path) -> list[str]:
    outputs = []
    for suffix in (".pdf", ".png"):
        out = stem.with_suffix(suffix)
        fig.savefig(out, dpi=DPI, bbox_inches="tight")
        outputs.append(str(out))
    plt.close(fig)
    return outputs


def format_contrast_label(row: pd.Series) -> str:
    if row["contrast_type"] == "pairwise":
        contrast = f"{row['positive_label']} vs {row['negative_label']}"
    else:
        contrast = f"{row['positive_label']} vs rest"
    return f"{row['cell_type']} | {contrast}"


def load_tables(summary_dir: Path) -> dict[str, pd.DataFrame]:
    tables = {
        "contrast": pd.read_csv(summary_dir / "contrast_summary.tsv", sep="\t"),
        "site_specific": pd.read_csv(summary_dir / "site_specific_core_genes.tsv", sep="\t"),
        "shared": pd.read_csv(summary_dir / "shared_core_genes.tsv", sep="\t"),
        "celltype": pd.read_csv(summary_dir / "celltype_pseudobulk_summary.tsv", sep="\t"),
        "method_status": pd.read_csv(summary_dir / "stability_method_status.tsv", sep="\t"),
    }
    return tables


def classify_gene_filter_reason(gene: object) -> str | None:
    gene_text = str(gene or "").strip()
    if not gene_text:
        return None
    for reason, pattern in STANDARD_GENE_FILTER_RULES:
        if re.search(pattern, gene_text, flags=re.IGNORECASE):
            return reason
    return None


def apply_standard_gene_filter(df: pd.DataFrame, source_table: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int | str]]:
    annotated = df.copy()
    annotated["gene"] = annotated["gene"].astype(str)
    annotated["filter_reason"] = annotated["gene"].map(classify_gene_filter_reason)

    removed = annotated.loc[annotated["filter_reason"].notna()].copy()
    kept = annotated.loc[annotated["filter_reason"].isna()].drop(columns=["filter_reason"]).copy()

    if not removed.empty:
        removed.insert(0, "source_table", source_table)

    summary = {
        "source_table": source_table,
        "input_rows": int(df.shape[0]),
        "retained_rows": int(kept.shape[0]),
        "filtered_rows": int(removed.shape[0]),
    }
    return kept, removed, summary


def write_gene_filter_artifacts(
    plot_data_dir: Path,
    filter_summaries: list[dict[str, int | str]],
    audit_df: pd.DataFrame,
) -> dict[str, object]:
    summary_df = pd.DataFrame(filter_summaries)
    summary_path = plot_data_dir / "gene_filter_summary.tsv"
    audit_path = plot_data_dir / "gene_filter_audit.tsv"

    example_removed_genes: dict[str, dict[str, list[str]]] = {}
    if not audit_df.empty:
        grouped = audit_df.groupby(["source_table", "filter_reason"], dropna=False, sort=True)
        for (source_table, filter_reason), sub_df in grouped:
            example_removed_genes.setdefault(str(source_table), {})[str(filter_reason)] = sub_df["gene"].astype(str).head(5).tolist()
        summary_df["example_removed_genes"] = summary_df["source_table"].map(
            lambda source: json.dumps(example_removed_genes.get(str(source), {}), ensure_ascii=False)
        )
    else:
        summary_df["example_removed_genes"] = "{}"

    summary_df.to_csv(summary_path, sep="\t", index=False)

    audit_cols = ["source_table", "gene", "filter_reason"]
    if audit_df.empty:
        pd.DataFrame(columns=audit_cols).to_csv(audit_path, sep="\t", index=False)
    else:
        audit_df.loc[:, audit_cols].to_csv(audit_path, sep="\t", index=False)

    return {
        "enabled": True,
        "profile_name": "bcell_standard_top_gene_filter_20260531",
        "rules": [{"reason": reason, "pattern": pattern} for reason, pattern in STANDARD_GENE_FILTER_RULES],
        "summary_table": str(summary_path),
        "audit_table": str(audit_path),
        "source_table_counts": {
            str(row["source_table"]): {
                "input_rows": int(row["input_rows"]),
                "retained_rows": int(row["retained_rows"]),
                "filtered_rows": int(row["filtered_rows"]),
            }
            for row in summary_df.to_dict(orient="records")
        },
        "example_removed_genes": example_removed_genes,
    }


def prep_contrast_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["cell_type"] = pd.Categorical(out["cell_type"], categories=CELLTYPE_ORDER, ordered=True)
    out["display_label"] = out.apply(format_contrast_label, axis=1)
    out["label_wrapped"] = out["display_label"].str.replace(" | ", "\n", regex=False)
    out = out.sort_values(
        ["n_stable_core_genes", "mean_balanced_accuracy", "cell_type", "contrast_type"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    return out


def prep_site_specific_table(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    out = df.copy()
    out["cell_type"] = pd.Categorical(out["cell_type"], categories=CELLTYPE_ORDER, ordered=True)
    out["enriched_site"] = pd.Categorical(out["enriched_site"], categories=SITE_ORDER, ordered=True)
    out = out.sort_values(
        ["cell_type", "enriched_site", "best_consensus_score", "n_supporting_contrasts", "mean_abs_log2_fc"],
        ascending=[True, True, False, False, False],
    )
    out = out.groupby(
        ["cell_type", "enriched_site"],
        as_index=False,
        group_keys=False,
        observed=False,
    ).head(top_n).reset_index(drop=True)
    return out


def prep_shared_table(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["cell_type"] = pd.Categorical(out["cell_type"], categories=CELLTYPE_ORDER, ordered=True)
    out = out.sort_values(["best_consensus_score", "mean_selection_freq"], ascending=[False, False]).reset_index(drop=True)
    out["display_label"] = out["gene"] + " (" + out["cell_type"].astype(str) + ")"
    return out


def prep_method_status(df: pd.DataFrame) -> pd.DataFrame:
    out = (
        df.groupby(["method", "status"], dropna=False)
        .size()
        .reset_index(name="n_records")
    )
    return out


def write_plot_tables(plot_data_dir: Path, contrast: pd.DataFrame, site_specific: pd.DataFrame, shared: pd.DataFrame, method_status: pd.DataFrame, celltype: pd.DataFrame) -> dict[str, str]:
    outputs = {
        "contrast": str(plot_data_dir / "contrast_overview_plot.tsv"),
        "site_specific": str(plot_data_dir / "site_specific_top_plot.tsv"),
        "shared": str(plot_data_dir / "shared_core_plot.tsv"),
        "method_status": str(plot_data_dir / "method_status_counts.tsv"),
        "celltype": str(plot_data_dir / "celltype_coverage_plot.tsv"),
    }
    contrast.to_csv(outputs["contrast"], sep="\t", index=False)
    site_specific.to_csv(outputs["site_specific"], sep="\t", index=False)
    shared.to_csv(outputs["shared"], sep="\t", index=False)
    method_status.to_csv(outputs["method_status"], sep="\t", index=False)
    celltype.to_csv(outputs["celltype"], sep="\t", index=False)
    return outputs


def plot_contrast_overview(df: pd.DataFrame, figure_dir: Path) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(20, 9), gridspec_kw={"width_ratios": [1.2, 1.0]})

    bar_df = df.iloc[::-1].reset_index(drop=True)
    bar_colors = [CELLTYPE_COLORS.get(str(ct), "#777777") for ct in bar_df["cell_type"]]
    axes[0].barh(bar_df["label_wrapped"], bar_df["n_stable_core_genes"], color=bar_colors, edgecolor="black", alpha=0.88)
    axes[0].set_xlabel("Stable core genes")
    axes[0].set_title("Stable core genes by contrast")
    axes[0].grid(axis="x", alpha=0.25, linestyle="--")
    for idx, value in enumerate(bar_df["n_stable_core_genes"]):
        axes[0].text(value + 0.35, idx, str(int(value)), va="center", fontsize=8)

    for contrast_type, marker in CONTRAST_MARKERS.items():
        subset = df[df["contrast_type"] == contrast_type]
        axes[1].scatter(
            subset["mean_balanced_accuracy"],
            subset["n_stable_core_genes"],
            s=40 + subset["n_samples_total"] * 6,
            c=[CELLTYPE_COLORS.get(str(ct), "#777777") for ct in subset["cell_type"]],
            marker=marker,
            edgecolor="black",
            linewidth=0.6,
            alpha=0.9,
            label=contrast_type,
        )

    annotate_df = df.sort_values(["n_stable_core_genes", "mean_balanced_accuracy"], ascending=[False, False]).head(6)
    for _, row in annotate_df.iterrows():
        short_label = f"{row['cell_type']}\n{row['positive_label']}"
        axes[1].annotate(
            short_label,
            (row["mean_balanced_accuracy"], row["n_stable_core_genes"]),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
        )

    axes[1].set_xlabel("Mean balanced accuracy")
    axes[1].set_ylabel("Stable core genes")
    axes[1].set_title("Contrast quality vs core-gene yield")
    axes[1].grid(alpha=0.25, linestyle="--")
    axes[1].set_xlim(0.64, 1.02)
    axes[1].set_ylim(-1, df["n_stable_core_genes"].max() + 4)

    celltype_handles = [
        Patch(facecolor=color, edgecolor="black", label=cell_type)
        for cell_type, color in CELLTYPE_COLORS.items()
        if cell_type in set(df["cell_type"].astype(str))
    ]
    type_handles = [
        Line2D([0], [0], marker=marker, color="white", markerfacecolor="#666666", markeredgecolor="black", markersize=8, linestyle="", label=contrast_type)
        for contrast_type, marker in CONTRAST_MARKERS.items()
    ]
    leg1 = axes[1].legend(handles=celltype_handles, title="Cell type", loc="lower left")
    axes[1].add_artist(leg1)
    axes[1].legend(handles=type_handles, title="Contrast type", loc="lower right")

    fig.suptitle("B-cell core discovery overview (2026-05-31 full run)", fontsize=15, y=1.02)
    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_core_discovery_contrast_overview")


def plot_site_specific_dotplot(df: pd.DataFrame, figure_dir: Path) -> list[str]:
    plot_df = df.copy()
    celltypes = [ct for ct in CELLTYPE_ORDER if ct in set(plot_df["cell_type"].astype(str))]
    fig, axes = plt.subplots(
        len(celltypes),
        1,
        figsize=(14, max(8, len(celltypes) * 3.4)),
        sharex=True,
        constrained_layout=True,
    )
    if len(celltypes) == 1:
        axes = [axes]

    scatter_mappable = None
    for ax, celltype in zip(axes, celltypes):
        sub = plot_df[plot_df["cell_type"].astype(str) == celltype].copy()
        sub = sub.sort_values(["enriched_site", "best_consensus_score"], ascending=[True, False]).reset_index(drop=True)
        sub["site_pos"] = sub["enriched_site"].astype(str).map({site: idx for idx, site in enumerate(SITE_ORDER)})
        sub["gene_label"] = sub["gene"].astype(str)
        y_positions = np.arange(len(sub))[::-1]
        scatter_mappable = ax.scatter(
            sub["site_pos"],
            y_positions,
            c=sub["best_consensus_score"],
            s=90 + sub["n_supporting_contrasts"] * 80,
            cmap="viridis",
            edgecolor="black",
            linewidth=0.5,
            alpha=0.92,
        )
        ax.set_yticks(y_positions)
        ax.set_yticklabels(sub["gene_label"])
        ax.set_title(f"{celltype}: top site-specific core genes", loc="left", fontsize=12)
        ax.grid(axis="x", alpha=0.2, linestyle="--")
        ax.set_xlim(-0.5, len(SITE_ORDER) - 0.5)
        for x in range(len(SITE_ORDER)):
            ax.axvline(x, color="#dddddd", lw=0.5, zorder=0)

    axes[-1].set_xticks(range(len(SITE_ORDER)))
    axes[-1].set_xticklabels(SITE_ORDER, rotation=0)
    fig.suptitle(f"Top site-specific core genes by cell type/site\n({GENE_FILTER_TITLE_NOTE})", fontsize=15, y=1.01)
    cbar = fig.colorbar(scatter_mappable, ax=axes, fraction=0.025, pad=0.015)
    cbar.set_label("Best consensus score")
    return save_figure(fig, figure_dir / "bcell_core_discovery_site_specific_dotplot")


def plot_shared_core(df: pd.DataFrame, figure_dir: Path) -> list[str]:
    fig, ax = plt.subplots(figsize=(12, 4.8))
    y = np.arange(len(df))[::-1]
    colors = [CELLTYPE_COLORS.get(str(ct), "#777777") for ct in df["cell_type"]]
    ax.hlines(y, xmin=0, xmax=df["best_consensus_score"], color="#bbbbbb", linewidth=2)
    ax.scatter(
        df["best_consensus_score"],
        y,
        s=140 + df["n_target_sites"] * 120,
        c=colors,
        edgecolor="black",
        linewidth=0.6,
        alpha=0.9,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(df["display_label"])
    ax.set_xlabel("Best consensus score")
    ax.set_title(f"Shared core genes across multiple target sites\n({GENE_FILTER_TITLE_NOTE})")
    ax.grid(axis="x", alpha=0.25, linestyle="--")
    ax.set_xlim(0, max(0.95, float(df["best_consensus_score"].max()) + 0.08))
    for yi, (_, row) in zip(y, df.iterrows()):
        ax.text(
            float(row["best_consensus_score"]) + 0.015,
            yi,
            f"sites: {row['supported_sites']}",
            va="center",
            fontsize=8,
        )
    handles = [
        Patch(facecolor=color, edgecolor="black", label=cell_type)
        for cell_type, color in CELLTYPE_COLORS.items()
        if cell_type in set(df["cell_type"].astype(str))
    ]
    ax.legend(handles=handles, title="Cell type", loc="lower right")
    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_core_discovery_shared_core")


def plot_run_diagnostics(method_status: pd.DataFrame, celltype_df: pd.DataFrame, figure_dir: Path) -> list[str]:
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [1.2, 1.0]})

    counts = method_status.pivot(index="method", columns="status", values="n_records").fillna(0)
    for status in ["ok", "failed", "skipped"]:
        if status not in counts.columns:
            counts[status] = 0
    counts = counts[["ok", "failed", "skipped"]]
    counts.plot(
        kind="barh",
        stacked=True,
        ax=axes[0],
        color=[STATUS_COLORS[status] for status in counts.columns],
        edgecolor="black",
        alpha=0.9,
    )
    axes[0].set_xlabel("Number of contrast records")
    axes[0].set_ylabel("Method")
    axes[0].set_title("Model-fit status across contrasts")
    axes[0].grid(axis="x", alpha=0.25, linestyle="--")
    axes[0].legend(title="Status", loc="lower right")

    coverage = celltype_df.copy()
    coverage["cell_type"] = pd.Categorical(coverage["cell_type"], categories=CELLTYPE_ORDER, ordered=True)
    coverage = coverage.sort_values("cell_type")
    x = np.arange(len(coverage))
    width = 0.37
    axes[1].bar(x - width / 2, coverage["n_samples"], width=width, color="#4c78a8", edgecolor="black", label="Samples")
    axes[1].bar(x + width / 2, coverage["n_sites"], width=width, color="#f58518", edgecolor="black", label="Sites")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(coverage["cell_type"].astype(str), rotation=20, ha="right")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Pseudobulk coverage by cell type")
    axes[1].grid(axis="y", alpha=0.25, linestyle="--")
    axes[1].legend(loc="upper right")
    for xi, row in enumerate(coverage.itertuples(index=False)):
        axes[1].text(xi - width / 2, row.n_samples + 0.5, str(int(row.n_samples)), ha="center", va="bottom", fontsize=8)
        axes[1].text(xi + width / 2, row.n_sites + 0.2, str(int(row.n_sites)), ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    return save_figure(fig, figure_dir / "bcell_core_discovery_run_diagnostics")


def main() -> None:
    args = parse_args()
    core_dir = args.run_root / "core_discovery"
    summary_dir = core_dir / "summaries"
    figure_dir = core_dir / "figures"
    plot_data_dir = figure_dir / "plot_data"
    ensure_dirs([figure_dir, plot_data_dir])

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "figure.facecolor": "white",
    })

    tables = load_tables(summary_dir)
    contrast_df = prep_contrast_table(tables["contrast"])
    site_specific_source_df, site_specific_removed_df, site_specific_filter_summary = apply_standard_gene_filter(
        tables["site_specific"],
        source_table="site_specific_core_genes",
    )
    shared_source_df, shared_removed_df, shared_filter_summary = apply_standard_gene_filter(
        tables["shared"],
        source_table="shared_core_genes",
    )
    site_specific_df = prep_site_specific_table(site_specific_source_df, args.top_site_genes_per_group)
    shared_df = prep_shared_table(shared_source_df)
    method_status_counts = prep_method_status(tables["method_status"])
    celltype_df = tables["celltype"].copy()

    gene_filter_audit_df = pd.concat(
        [df for df in [site_specific_removed_df, shared_removed_df] if not df.empty],
        ignore_index=True,
    ) if (not site_specific_removed_df.empty or not shared_removed_df.empty) else pd.DataFrame()

    gene_filter_manifest = write_gene_filter_artifacts(
        plot_data_dir=plot_data_dir,
        filter_summaries=[site_specific_filter_summary, shared_filter_summary],
        audit_df=gene_filter_audit_df,
    )
    gene_filter_manifest["plotted_rows"] = {
        "site_specific_dotplot": int(site_specific_df.shape[0]),
        "shared_core": int(shared_df.shape[0]),
    }

    plot_table_paths = write_plot_tables(
        plot_data_dir=plot_data_dir,
        contrast=contrast_df,
        site_specific=site_specific_df,
        shared=shared_df,
        method_status=method_status_counts,
        celltype=celltype_df,
    )

    figure_manifest: list[dict[str, object]] = []
    figure_manifest.append({
        "name": "contrast_overview",
        "outputs": plot_contrast_overview(contrast_df, figure_dir),
        "sources": [plot_table_paths["contrast"]],
    })
    figure_manifest.append({
        "name": "site_specific_dotplot",
        "outputs": plot_site_specific_dotplot(site_specific_df, figure_dir),
        "sources": [plot_table_paths["site_specific"]],
    })
    figure_manifest.append({
        "name": "shared_core",
        "outputs": plot_shared_core(shared_df, figure_dir),
        "sources": [plot_table_paths["shared"]],
    })
    figure_manifest.append({
        "name": "run_diagnostics",
        "outputs": plot_run_diagnostics(method_status_counts, celltype_df, figure_dir),
        "sources": [plot_table_paths["method_status"], plot_table_paths["celltype"]],
    })

    manifest_path = figure_dir / "figure_manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "run_root": str(args.run_root),
                "core_summary_dir": str(summary_dir),
                "plot_data_dir": str(plot_data_dir),
                "gene_filter": gene_filter_manifest,
                "figures": figure_manifest,
            },
            handle,
            indent=2,
            ensure_ascii=False,
        )

    print(f"[done] figures written to: {figure_dir}")
    print(f"[done] plot-data tables written to: {plot_data_dir}")
    print(f"[done] figure manifest: {manifest_path}")
    for item in figure_manifest:
        for output in item["outputs"]:
            print(f"  - {output}")


if __name__ == "__main__":
    main()
