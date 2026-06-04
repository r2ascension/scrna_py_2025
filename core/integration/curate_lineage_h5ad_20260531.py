#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import anndata as ad
import pandas as pd
from pandas.api.types import CategoricalDtype

ANNOTATION_COL_PATTERN = re.compile(
    r"(cell.*type|celltype|annotation|annot|label|scanvi|pred|major|minor|level|l[123]\b|final)",
    re.IGNORECASE,
)

DEFAULT_CONTAMINATION_PATTERN = r"contamination"


@dataclass
class CleanupSummary:
    lineage_name: str
    input_h5ad: str
    output_h5ad: str
    copy_strategy: str
    n_obs_input: int
    n_vars_input: int
    n_obs_output: int
    n_vars_output: int
    removed_total: int
    removed_contamination: int
    removed_tissue_label: int
    relabeled_total: int
    annotation_columns_used: list[str]
    contamination_hits_by_column: dict[str, dict[str, int]]
    tissue_label_drop_counts: dict[str, dict[str, int]]
    relabel_counts_by_column: dict[str, dict[str, int]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Curate a lineage main h5ad by removing contamination and lineage-specific subsets.")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--removed-cells-tsv", required=True)
    parser.add_argument("--lineage-name", required=True)
    parser.add_argument("--annotation-columns", default="", help="Comma-separated annotation columns. If omitted, they are auto-detected.")
    parser.add_argument("--contamination-pattern", default=DEFAULT_CONTAMINATION_PATTERN)
    parser.add_argument("--tissue-column", default="")
    parser.add_argument("--drop-tissues", default="", help="Comma-separated tissue names to drop when paired with label hits.")
    parser.add_argument("--drop-label-column", default="")
    parser.add_argument("--drop-label-pattern", default="")
    parser.add_argument("--relabel-columns", default="", help="Comma-separated columns for exact relabeling.")
    parser.add_argument("--relabel-map", action="append", default=[], help="Exact relabel mapping in OLD=NEW form. Repeatable.")
    parser.add_argument("--compression-level", type=int, default=4)
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    if not value.strip():
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_relabel_map(items: Iterable[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --relabel-map entry: {item!r}; expected OLD=NEW")
        old, new = item.split("=", 1)
        old = old.strip()
        new = new.strip()
        if not old:
            raise ValueError(f"Invalid --relabel-map entry with empty OLD value: {item!r}")
        mapping[old] = new
    return mapping


def detect_annotation_columns(obs_columns: Iterable[str]) -> list[str]:
    detected = [col for col in obs_columns if ANNOTATION_COL_PATTERN.search(str(col))]
    return sorted(dict.fromkeys(str(col) for col in detected))


def normalize_string_series(series: pd.Series) -> pd.Series:
    normalized = series.astype("string")
    return normalized.str.strip()


def build_contamination_mask(
    obs: pd.DataFrame,
    annotation_columns: Iterable[str],
    contamination_pattern: str,
) -> tuple[pd.Series, dict[str, dict[str, int]]]:
    pattern = re.compile(contamination_pattern, re.IGNORECASE)
    mask = pd.Series(False, index=obs.index)
    hits: dict[str, dict[str, int]] = {}

    for column in annotation_columns:
        if column not in obs.columns:
            continue
        values = normalize_string_series(obs[column])
        hit = values.str.contains(pattern, na=False)
        if not hit.any():
            continue
        hits[column] = {str(label): int(count) for label, count in values[hit].value_counts().items()}
        mask |= hit

    return mask, hits


def build_tissue_label_mask(
    obs: pd.DataFrame,
    tissue_column: str,
    tissues: Iterable[str],
    label_column: str,
    label_pattern: str,
) -> tuple[pd.Series, dict[str, dict[str, int]]]:
    mask = pd.Series(False, index=obs.index)
    counts: dict[str, dict[str, int]] = {}

    if not tissue_column or not label_column or not label_pattern:
        return mask, counts
    if tissue_column not in obs.columns or label_column not in obs.columns:
        missing = [col for col in (tissue_column, label_column) if col and col not in obs.columns]
        raise KeyError(f"Required cleanup column(s) missing: {missing}")

    tissues_norm = {item.strip() for item in tissues if item.strip()}
    if not tissues_norm:
        return mask, counts

    tissue_values = normalize_string_series(obs[tissue_column])
    label_values = normalize_string_series(obs[label_column])
    tissue_mask = tissue_values.isin(tissues_norm)
    label_mask = label_values.str.contains(re.compile(label_pattern, re.IGNORECASE), na=False)
    mask = tissue_mask & label_mask

    if mask.any():
        ct = pd.crosstab(tissue_values[mask], label_values[mask])
        counts = {
            str(label): {
                str(tissue): int(value)
                for tissue, value in tissue_map.items()
                if int(value) > 0
            }
            for label, tissue_map in ct.to_dict().items()
            if any(int(value) > 0 for value in tissue_map.values())
        }

    return mask, counts


def apply_exact_relabels(obs: pd.DataFrame, columns: Iterable[str], mapping: dict[str, str]) -> dict[str, dict[str, int]]:
    relabel_counts: dict[str, dict[str, int]] = {}
    if not mapping:
        return relabel_counts

    for column in columns:
        if column not in obs.columns:
            continue
        original = obs[column]
        values = normalize_string_series(original)
        mask = values.isin(mapping.keys())
        if not mask.any():
            continue

        relabel_counts[column] = {old: int((values == old).sum()) for old in mapping if int((values == old).sum()) > 0}
        updated = values.replace(mapping)
        if isinstance(original.dtype, CategoricalDtype):
            obs[column] = pd.Categorical(updated)
        else:
            obs[column] = updated.astype(object)

    return relabel_counts


def count_exact_relabels(obs: pd.DataFrame, columns: Iterable[str], mapping: dict[str, str]) -> dict[str, dict[str, int]]:
    relabel_counts: dict[str, dict[str, int]] = {}
    if not mapping:
        return relabel_counts

    for column in columns:
        if column not in obs.columns:
            continue
        values = normalize_string_series(obs[column])
        column_counts = {old: int((values == old).sum()) for old in mapping if int((values == old).sum()) > 0}
        if column_counts:
            relabel_counts[column] = column_counts

    return relabel_counts


def close_backed_file(adata: ad.AnnData | None) -> None:
    if adata is None:
        return
    file_manager = getattr(adata, "file", None)
    if file_manager is None:
        return
    close_fn = getattr(file_manager, "close", None)
    if close_fn is not None:
        close_fn()


def write_filtered_h5ad(
    input_h5ad: Path,
    keep_mask: pd.Series,
    output_h5ad: Path,
    compression_level: int,
) -> str:
    backed = None
    subset = None
    try:
        backed = ad.read_h5ad(input_h5ad, backed="r")
        subset = backed[keep_mask.to_numpy()]
        subset.copy(filename=str(output_h5ad))
        return "backed_copy"
    except Exception as exc:
        print(f"[curate] backed copy unavailable, falling back to in-memory write: {exc}", flush=True)
        full = ad.read_h5ad(input_h5ad)
        filtered = full[keep_mask.to_numpy()].copy()
        filtered.write_h5ad(
            output_h5ad,
            compression="gzip",
            compression_opts=compression_level,
        )
        return "in_memory_write"
    finally:
        close_backed_file(subset)
        close_backed_file(backed)


def relabel_output_h5ad(
    output_h5ad: Path,
    relabel_columns: Iterable[str],
    relabel_map: dict[str, str],
    compression_level: int,
) -> dict[str, dict[str, int]]:
    tmp_output = output_h5ad.with_name(f"{output_h5ad.stem}.tmp{output_h5ad.suffix}")
    adata = ad.read_h5ad(output_h5ad)
    relabel_counts = apply_exact_relabels(adata.obs, relabel_columns, relabel_map)
    adata.write_h5ad(
        tmp_output,
        compression="gzip",
        compression_opts=compression_level,
    )
    tmp_output.replace(output_h5ad)
    return relabel_counts


def write_removed_cells_tsv(
    removed_cells_tsv: Path,
    remove_mask: pd.Series,
    contamination_mask: pd.Series,
    tissue_label_mask: pd.Series,
) -> None:
    removed_df = pd.DataFrame(
        {
            "obs_name": remove_mask.index.astype(str),
            "remove": remove_mask.to_numpy(),
            "remove_contamination": contamination_mask.reindex(remove_mask.index, fill_value=False).to_numpy(),
            "remove_tissue_label": tissue_label_mask.reindex(remove_mask.index, fill_value=False).to_numpy(),
        }
    )
    removed_df = removed_df.loc[removed_df["remove"]].copy()
    removed_df.to_csv(removed_cells_tsv, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    input_h5ad = Path(args.input_h5ad)
    output_h5ad = Path(args.output_h5ad)
    summary_json = Path(args.summary_json)
    removed_cells_tsv = Path(args.removed_cells_tsv)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    removed_cells_tsv.parent.mkdir(parents=True, exist_ok=True)

    annotation_columns = split_csv(args.annotation_columns)
    relabel_columns = split_csv(args.relabel_columns)
    drop_tissues = split_csv(args.drop_tissues)
    relabel_map = parse_relabel_map(args.relabel_map)

    print(f"[curate] loading metadata/obs from {input_h5ad}", flush=True)
    adata_backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        obs = adata_backed.obs.copy()
        obs.columns = obs.columns.astype(str)
        n_obs_input = int(adata_backed.n_obs)
        n_vars_input = int(adata_backed.n_vars)
    finally:
        close_backed_file(adata_backed)
    print(f"[curate] loaded shape={n_obs_input}x{n_vars_input}", flush=True)

    if not annotation_columns:
        annotation_columns = detect_annotation_columns(obs.columns)

    contamination_mask, contamination_hits = build_contamination_mask(
        obs,
        annotation_columns,
        args.contamination_pattern,
    )
    tissue_label_mask, tissue_label_counts = build_tissue_label_mask(
        obs,
        args.tissue_column,
        drop_tissues,
        args.drop_label_column,
        args.drop_label_pattern,
    )

    remove_mask = contamination_mask | tissue_label_mask
    keep_mask = ~remove_mask
    relabel_preview = count_exact_relabels(obs.loc[keep_mask.to_numpy(), :], relabel_columns, relabel_map)
    relabeled_total = int(sum(sum(col_counts.values()) for col_counts in relabel_preview.values()))

    print(
        f"[curate] removed_total={int(remove_mask.sum())} "
        f"removed_contamination={int(contamination_mask.sum())} "
        f"removed_tissue_label={int(tissue_label_mask.sum())} "
        f"relabeled_total={relabeled_total}",
        flush=True,
    )
    print(f"[curate] writing {output_h5ad}", flush=True)
    copy_strategy = write_filtered_h5ad(
        input_h5ad,
        keep_mask,
        output_h5ad,
        args.compression_level,
    )
    if relabeled_total > 0:
        print(f"[curate] applying exact relabels in {output_h5ad}", flush=True)
        relabel_counts = relabel_output_h5ad(
            output_h5ad,
            relabel_columns,
            relabel_map,
            args.compression_level,
        )
    else:
        relabel_counts = relabel_preview
    write_removed_cells_tsv(removed_cells_tsv, remove_mask, contamination_mask, tissue_label_mask)

    summary = CleanupSummary(
        lineage_name=args.lineage_name,
        input_h5ad=str(input_h5ad),
        output_h5ad=str(output_h5ad),
        copy_strategy=copy_strategy,
        n_obs_input=n_obs_input,
        n_vars_input=n_vars_input,
        n_obs_output=int(keep_mask.sum()),
        n_vars_output=n_vars_input,
        removed_total=int(remove_mask.sum()),
        removed_contamination=int(contamination_mask.sum()),
        removed_tissue_label=int(tissue_label_mask.sum()),
        relabeled_total=int(sum(sum(col_counts.values()) for col_counts in relabel_counts.values())),
        annotation_columns_used=annotation_columns,
        contamination_hits_by_column=contamination_hits,
        tissue_label_drop_counts=tissue_label_counts,
        relabel_counts_by_column=relabel_counts,
    )
    summary_json.write_text(json.dumps(asdict(summary), indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[curate] wrote summary {summary_json}", flush=True)
    print(json.dumps(asdict(summary), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
