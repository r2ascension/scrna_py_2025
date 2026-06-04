#!/usr/bin/env python3
"""
curate_query_cleanup_20260604.py
================================
LLM-guided cluster removal from query h5ad files.

Reads the LLM discovery screen TSV produced by the CHOIR/OFA/LLM R pipeline,
identifies clusters flagged as outliers, and removes the corresponding cells
from the input h5ad using the same backed='r' memory-efficient strategy as
curate_lineage_h5ad_20260531.py.

Typical usage:
    python curate_query_cleanup_20260604.py \\
        --input-h5ad /path/to/input.h5ad \\
        --cluster-csv /path/to/reports/choir/choir_clusters.csv \\
        --llm-screen-tsv /path/to/reports/llm_choir_discovery_screen.tsv \\
        --output-h5ad /path/to/cleaned.h5ad \\
        --summary-json /path/to/summary.json \\
        --removed-cells-tsv /path/to/removed_cells.tsv \\
        --lineage-name epithelial \\
        --cluster-column CHOIR_clusters_0.2

Removal rules (configurable via --remove-rules JSON):
    Default: remove clusters where outlier_flag=="yes" OR
             biological_signal_class=="likely_outlier"

Date: 2026-06-04
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd


# ── reuse proven functions from the main cleanup script ──────────────────────
# We inline them here to keep this script self-contained while preserving the
# identical behaviour.  The two scripts share the same backed-copy + fallback
# logic; any future improvements should be applied to both.
# -----------------------------------------------------------------------------

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
    """Write a filtered copy of *input_h5ad* keeping only rows where *keep_mask* is True.

    Returns the strategy used: ``"backed_copy"`` or ``"in_memory_write"``.
    """
    backed = None
    subset = None
    try:
        backed = ad.read_h5ad(input_h5ad, backed="r")
        subset = backed[keep_mask.to_numpy()]
        subset.copy(filename=str(output_h5ad))
        return "backed_copy"
    except Exception as exc:
        print(
            f"[curate_query] backed copy unavailable, falling back to in-memory write: {exc}",
            flush=True,
        )
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


# ── dataclass for the output summary ─────────────────────────────────────────

@dataclass
class QueryCleanupSummary:
    lineage_name: str
    input_h5ad: str
    output_h5ad: str
    copy_strategy: str
    n_obs_input: int
    n_vars_input: int
    n_obs_output: int
    n_vars_output: int
    removed_total: int
    removed_by_cluster: dict[str, int] = field(default_factory=dict)
    clusters_removed: list[str] = field(default_factory=list)
    removal_rules_used: dict[str, str] = field(default_factory=dict)
    llm_screen_tsv: str = ""


# ── argument parsing ─────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="LLM-guided cluster removal from query h5ad files."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--cluster-csv", required=True,
                        help="choir_clusters.csv (cell → cluster_id mapping)")
    parser.add_argument("--llm-screen-tsv", required=True,
                        help="LLM discovery screen TSV from R pipeline")
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--removed-cells-tsv", required=True)
    parser.add_argument("--lineage-name", required=True)
    parser.add_argument("--cluster-column", default="CHOIR_clusters_0.2",
                        help="obs column name for cluster assignments (default: CHOIR_clusters_0.2)")
    parser.add_argument("--remove-rules", default="",
                        help="JSON dict of {column: value} conditions for removal. "
                             "Default: outlier_flag==yes OR biological_signal_class==likely_outlier")
    parser.add_argument("--compression-level", type=int, default=4)
    return parser.parse_args()


# ── helpers ──────────────────────────────────────────────────────────────────

def _parse_remove_rules(raw: str) -> list[dict[str, str]]:
    """Parse --remove-rules into a list of {column: value} dicts.

    If no rules are supplied, return the safe defaults."""
    if not raw.strip():
        return [
            {"outlier_flag": "yes"},
            {"biological_signal_class": "likely_outlier"},
        ]
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        # Fallback: comma-separated COLUMN=VALUE pairs
        parsed = {}
        for pair in raw.split(","):
            pair = pair.strip()
            if "=" not in pair:
                continue
            k, v = pair.split("=", 1)
            parsed[k.strip()] = v.strip()
    if isinstance(parsed, dict):
        # single rule dict → wrap in list
        return [parsed]
    if isinstance(parsed, list):
        return parsed
    raise ValueError(f"Could not parse --remove-rules: {raw!r}")


def _parse_cluster_id(record_id: str) -> str | None:
    """Extract cluster id from a record_id like 'ofa_0012' or 'choir_0012'.

    Returns the leading integer found, or None.
    """
    m = re.search(r"(\d+)", str(record_id))
    return m.group(1) if m else None


def load_cluster_mapping(cluster_csv: Path) -> pd.DataFrame:
    """Load choir_clusters.csv and return a DataFrame with 'cell' and 'cluster_id' columns."""
    df = pd.read_csv(cluster_csv)
    # Auto-detect column names — the R pipeline writes with column names 'cell' and 'cluster_id'
    cols = [str(c).strip().lower() for c in df.columns]
    cell_col = None
    cluster_col = None
    for c in df.columns:
        lc = str(c).strip().lower()
        if lc in ("cell", "barcode", "obs_name", "cell_id", "cell_barcode"):
            cell_col = c
        elif lc in ("cluster_id", "cluster", "choir_cluster", "leiden_cluster"):
            cluster_col = c
    if cell_col is None or cluster_col is None:
        raise ValueError(
            f"Cannot identify cell/cluster columns in {cluster_csv}. "
            f"Columns found: {list(df.columns)}"
        )
    df = df.rename(columns={cell_col: "cell", cluster_col: "cluster_id"})
    df["cell"] = df["cell"].astype(str)
    df["cluster_id"] = df["cluster_id"].astype(str)
    return df


def extract_clusters_to_remove(
    llm_screen_tsv: Path,
    removal_rules: list[dict[str, str]],
) -> tuple[set[str], dict[str, str]]:
    """Parse the LLM discovery screen and return the set of cluster ids to remove.

    Returns (to_remove, reason_map) where reason_map is {cluster_id: reason_string}.
    """
    df = pd.read_csv(llm_screen_tsv, sep="\t")
    # Normalise column names
    df.columns = [str(c).strip() for c in df.columns]

    to_remove: set[str] = set()
    reason_map: dict[str, str] = {}

    for _, row in df.iterrows():
        hit_reasons: list[str] = []
        for rule in removal_rules:
            match = True
            for col, expected_val in rule.items():
                actual = str(row.get(col, "")).strip().lower()
                if actual != str(expected_val).strip().lower():
                    match = False
                    break
            if match:
                hit_reasons.append(" AND ".join(f"{k}={v}" for k, v in rule.items()))

        if not hit_reasons:
            continue

        record_id = str(row.get("record_id", row.get("record_label", "")))
        cluster_id = _parse_cluster_id(record_id)
        if cluster_id is None:
            print(f"[curate_query] WARNING: could not parse cluster id from record_id={record_id!r}",
                  flush=True)
            continue

        to_remove.add(cluster_id)
        reason_map[cluster_id] = "; ".join(hit_reasons)

    return to_remove, reason_map


def build_remove_mask(
    obs_index: pd.Index,
    cluster_df: pd.DataFrame,
    clusters_to_remove: set[str],
) -> pd.Series:
    """Build a boolean mask: True = cell should be removed."""
    # Join obs index with cluster mapping
    cluster_series = obs_index.to_series().map(
        dict(zip(cluster_df["cell"], cluster_df["cluster_id"]))
    )
    # Cells not found in cluster mapping are kept (mask=False)
    mask = cluster_series.isin(clusters_to_remove)
    mask = mask.fillna(False)
    return mask


def write_removed_cells_tsv(
    path: Path,
    remove_mask: pd.Series,
    cluster_df: pd.DataFrame,
    clusters_to_remove: set[str],
    reason_map: dict[str, str],
) -> None:
    """Write a per-cell audit trail of removed cells."""
    obs_index = remove_mask.index
    cluster_series = obs_index.to_series().map(
        dict(zip(cluster_df["cell"], cluster_df["cluster_id"]))
    )
    removed = remove_mask[remove_mask].index
    df = pd.DataFrame({
        "obs_name": removed.astype(str),
        "cluster_id": cluster_series.loc[removed].astype(str),
        "removal_reason": cluster_series.loc[removed].map(
            lambda x: reason_map.get(str(x), "unknown")
        ).fillna("unknown"),
    })
    df.to_csv(path, sep="\t", index=False)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    input_h5ad = Path(args.input_h5ad)
    cluster_csv = Path(args.cluster_csv)
    llm_screen_tsv = Path(args.llm_screen_tsv)
    output_h5ad = Path(args.output_h5ad)
    summary_json = Path(args.summary_json)
    removed_cells_tsv = Path(args.removed_cells_tsv)

    for p in (output_h5ad, summary_json, removed_cells_tsv):
        p.parent.mkdir(parents=True, exist_ok=True)

    # 1. Parse removal rules
    removal_rules = _parse_remove_rules(args.remove_rules)
    print(f"[curate_query] removal rules: {removal_rules}", flush=True)

    # 2. Load cluster mapping
    print(f"[curate_query] loading cluster mapping from {cluster_csv}", flush=True)
    cluster_df = load_cluster_mapping(cluster_csv)
    print(f"[curate_query]   {len(cluster_df)} cells, "
          f"{cluster_df['cluster_id'].nunique()} clusters", flush=True)

    # 3. Identify clusters to remove from LLM screen
    print(f"[curate_query] loading LLM screen from {llm_screen_tsv}", flush=True)
    clusters_to_remove, reason_map = extract_clusters_to_remove(
        llm_screen_tsv, removal_rules
    )
    print(f"[curate_query]   clusters flagged for removal: "
          f"{sorted(clusters_to_remove, key=int) if clusters_to_remove else 'NONE'}", flush=True)

    # 4. Load obs metadata (backed='r' for memory efficiency)
    print(f"[curate_query] loading obs from {input_h5ad}", flush=True)
    adata_backed = ad.read_h5ad(input_h5ad, backed="r")
    try:
        obs_index = adata_backed.obs.index.copy()
        n_obs_input = int(adata_backed.n_obs)
        n_vars_input = int(adata_backed.n_vars)
    finally:
        close_backed_file(adata_backed)
    print(f"[curate_query]   shape={n_obs_input}x{n_vars_input}", flush=True)

    # 5. Build removal mask
    remove_mask = build_remove_mask(obs_index, cluster_df, clusters_to_remove)
    keep_mask = ~remove_mask
    removed_by_cluster: dict[str, int] = {}
    if remove_mask.any():
        cluster_series = obs_index.to_series().map(
            dict(zip(cluster_df["cell"], cluster_df["cluster_id"]))
        )
        vc = cluster_series[remove_mask].value_counts()
        removed_by_cluster = {str(k): int(v) for k, v in vc.items()}

    print(
        f"[curate_query] removed_total={int(remove_mask.sum())} "
        f"removed_by_cluster={removed_by_cluster}",
        flush=True,
    )

    # 6. Write filtered h5ad
    print(f"[curate_query] writing {output_h5ad}", flush=True)
    copy_strategy = write_filtered_h5ad(
        input_h5ad, keep_mask, output_h5ad, args.compression_level
    )
    print(f"[curate_query]   copy_strategy={copy_strategy}", flush=True)

    # 7. Write audit trail
    write_removed_cells_tsv(
        removed_cells_tsv,
        remove_mask,
        cluster_df,
        clusters_to_remove,
        reason_map,
    )
    print(f"[curate_query] wrote removed cells TSV {removed_cells_tsv}", flush=True)

    # 8. Write summary JSON
    summary = QueryCleanupSummary(
        lineage_name=args.lineage_name,
        input_h5ad=str(input_h5ad),
        output_h5ad=str(output_h5ad),
        copy_strategy=copy_strategy,
        n_obs_input=n_obs_input,
        n_vars_input=n_vars_input,
        n_obs_output=int(keep_mask.sum()),
        n_vars_output=n_vars_input,
        removed_total=int(remove_mask.sum()),
        removed_by_cluster=removed_by_cluster,
        clusters_removed=sorted(clusters_to_remove, key=int) if clusters_to_remove else [],
        removal_rules_used={list(r.keys())[0]: list(r.values())[0] for r in removal_rules},
        llm_screen_tsv=str(llm_screen_tsv),
    )
    summary_json.write_text(
        json.dumps(asdict(summary), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"[curate_query] wrote summary {summary_json}", flush=True)
    print(json.dumps(asdict(summary), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
