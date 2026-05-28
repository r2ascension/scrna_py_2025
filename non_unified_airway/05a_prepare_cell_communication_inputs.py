#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from non_unified_airway.communication import (
    COMM_LEVEL_SPECS,
    DEFAULT_COMM_OUT_DIR,
    DEFAULT_CONFIG_MANIFEST,
    DEFAULT_METADATA_TSV,
    apply_communication_contrast,
    build_communication_contrasts,
    build_sender_receiver_pairs,
    export_method_ready_subset,
    load_metadata_and_manifest,
    pick_primary_h5ad,
    prepare_level_subset,
    summarize_celltype_support,
    summarize_prepared_manifest,
    write_status,
)
from non_unified_airway.common import write_json, write_tsv


def write_table(path: Path, table: pd.DataFrame) -> None:
    kwargs = {"compression": "gzip"} if str(path).endswith(".gz") else {}
    write_tsv(path, table, **kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare L2/L3 pairwise communication inputs for non-unified airway analysis.")
    parser.add_argument("--config-manifest", type=Path, default=DEFAULT_CONFIG_MANIFEST)
    parser.add_argument("--metadata-tsv", type=Path, default=DEFAULT_METADATA_TSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_COMM_OUT_DIR)
    parser.add_argument("--h5ad", type=Path, default=None, help="Override primary all-cells h5ad path.")
    parser.add_argument("--min-samples-per-group", type=int, default=3)
    parser.add_argument("--healthy-only", action="store_true", help="Only build healthy site-to-site contrasts.")
    parser.add_argument("--export-h5ad", action="store_true", help="Write per-contrast/level AnnData subsets.")
    parser.add_argument("--export-mtx", action="store_true", help="Write MatrixMarket + metadata exports for R/CellPhoneDB.")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--max-cells-per-side-celltype", type=int, default=300)
    parser.add_argument("--random-state", type=int, default=0)
    return parser.parse_args()


def run_prepare_cell_communication_inputs(args: argparse.Namespace) -> dict[str, object]:
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata, sample_manifest = load_metadata_and_manifest(args.metadata_tsv, args.config_manifest)

    contrasts = build_communication_contrasts(
        sample_manifest,
        min_samples_per_group=args.min_samples_per_group,
        include_disease_extension=not args.healthy_only,
    )
    write_table(out_dir / "sample_manifest.tsv", sample_manifest)
    write_table(out_dir / "communication_contrasts.tsv", contrasts)

    export_requested = bool(args.export_h5ad or args.export_mtx)
    h5ad_path = pick_primary_h5ad(args.config_manifest, args.h5ad) if export_requested else None

    prepared_rows: list[dict[str, object]] = []
    for contrast in contrasts.to_dict(orient="records"):
        contrast_id = str(contrast["contrast_id"])
        contrast_dir = out_dir / "prepared" / contrast_id
        contrast_dir.mkdir(parents=True, exist_ok=True)

        contrast_subset = apply_communication_contrast(metadata, contrast)
        write_table(contrast_dir / "cells_all_levels.tsv.gz", contrast_subset)
        write_json(contrast_dir / "contrast_manifest.json", contrast)

        for level, spec in COMM_LEVEL_SPECS.items():
            level_dir = contrast_dir / level
            level_dir.mkdir(parents=True, exist_ok=True)
            level_subset, celltype_col = prepare_level_subset(metadata, contrast, level=level)
            row = {
                "contrast_id": contrast_id,
                "contrast_group": contrast["contrast_group"],
                "level": level,
                "status": "blocked",
                "prepared_dir": str(level_dir),
                "celltype_col": celltype_col,
                "n_cells": int(level_subset.shape[0]),
                "n_eligible_celltypes": 0,
                "left_label": contrast["left_label"],
                "right_label": contrast["right_label"],
                "left_sites": contrast["left_sites"],
                "right_sites": contrast["right_sites"],
                "left_disease_groups": contrast["left_disease_groups"],
                "right_disease_groups": contrast["right_disease_groups"],
                "analysis_compartments": contrast["analysis_compartments"],
                "h5ad_subset": "",
                "counts_mtx": "",
                "genes_tsv": "",
                "barcodes_tsv": "",
                "cell_metadata_tsv": "",
                "n_export_cells": 0,
            }
            if level_subset.empty:
                write_status(
                    level_dir / "prepare_status.json",
                    {
                        **row,
                        "status": "blocked",
                        "reason": "no_cells_after_contrast_and_annotation_filter",
                    },
                )
                prepared_rows.append(row)
                continue

            write_table(level_dir / "cell_assignments.tsv.gz", level_subset)
            support_df, sample_counts = summarize_celltype_support(
                level_subset,
                celltype_col="comm_celltype",
                min_cells_per_sample=spec["min_cells_per_sample"],
                min_samples_per_group=spec["min_samples_per_group"],
            )
            write_table(level_dir / "sample_celltype_counts.tsv.gz", sample_counts)
            write_table(level_dir / "eligible_celltypes.tsv", support_df)

            eligible = support_df.loc[support_df["eligible"], "comm_celltype"].astype(str).tolist()
            pair_df = build_sender_receiver_pairs(eligible)
            write_table(level_dir / "sender_receiver_pairs.tsv", pair_df)
            row["n_eligible_celltypes"] = len(eligible)
            if eligible:
                row["status"] = "ready"
                export_subset = level_subset.loc[level_subset["comm_celltype"].isin(eligible)].copy()
                if export_requested and h5ad_path is not None:
                    export_manifest = export_method_ready_subset(
                        h5ad_path,
                        export_subset,
                        level_dir,
                        counts_layer=args.counts_layer,
                        write_h5ad_subset=args.export_h5ad,
                        write_mtx=args.export_mtx,
                        max_cells_per_side_celltype=args.max_cells_per_side_celltype,
                        random_state=args.random_state,
                    )
                    row.update(export_manifest)
            else:
                row["status"] = "blocked"

            write_status(
                level_dir / "prepare_status.json",
                {
                    **row,
                    "min_cells_per_sample": spec["min_cells_per_sample"],
                    "min_samples_per_group": spec["min_samples_per_group"],
                },
            )
            prepared_rows.append(row)

    prepared_df = pd.DataFrame(prepared_rows)
    write_table(out_dir / "prepared_inputs_manifest.tsv", prepared_df)

    summary = summarize_prepared_manifest(prepared_rows)
    payload = {
        **summary,
        "export_h5ad": bool(args.export_h5ad),
        "export_mtx": bool(args.export_mtx),
        "counts_layer": args.counts_layer,
        "primary_h5ad": str(h5ad_path) if h5ad_path is not None else "",
    }
    write_json(out_dir / "prepare_summary.json", payload)
    return payload


def main() -> None:
    args = parse_args()
    summary = run_prepare_cell_communication_inputs(args)
    print(
        "[non_unified_airway][05a] prepared="
        f"{summary['n_prepared']} ready={summary['n_ready']} blocked={summary['n_blocked']}"
    )


if __name__ == "__main__":
    main()
