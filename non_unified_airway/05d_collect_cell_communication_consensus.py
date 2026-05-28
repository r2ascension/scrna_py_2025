#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from non_unified_airway.communication import DEFAULT_COMM_OUT_DIR
from non_unified_airway.common import write_json, write_tsv


def write_table(path: Path, table: pd.DataFrame) -> None:
    kwargs = {"compression": "gzip"} if str(path).endswith(".gz") else {}
    write_tsv(path, table, **kwargs)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Collect LIANA+/CellPhoneDB/CellChat pairwise communication outputs into a consensus table.")
    parser.add_argument("--liana-dir", type=Path, default=DEFAULT_COMM_OUT_DIR / "liana")
    parser.add_argument("--cellphonedb-dir", type=Path, default=DEFAULT_COMM_OUT_DIR / "cellphonedb")
    parser.add_argument("--cellchat-dir", type=Path, default=Path("/home/h2048/data/R/0525/non_unified_airway/communication/cellchat"))
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_COMM_OUT_DIR / "consensus")
    parser.add_argument("--top-n", type=int, default=500)
    return parser.parse_args()


def _load_delta_tables(base_dir: Path, pattern: str, method_name: str, score_col: str = "delta_score_right_minus_left") -> pd.DataFrame:
    rows = []
    if not base_dir.exists():
        return pd.DataFrame()
    for path in sorted(base_dir.glob(pattern)):
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        required = ["contrast_id", "level", "interaction_key", "sender", "receiver", "ligand", "receptor"]
        missing = [col for col in required if col not in df.columns]
        if missing:
            continue
        df = df.copy()
        df["method"] = method_name
        df["method_score"] = pd.to_numeric(df.get(score_col), errors="coerce")
        df["method_result_path"] = str(path)
        rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _summarize_consensus(all_methods: pd.DataFrame) -> pd.DataFrame:
    if all_methods.empty:
        return pd.DataFrame()

    def _method_list(values: pd.Series) -> str:
        uniq = sorted({str(value) for value in values if str(value).strip()})
        return ";".join(uniq)

    summary = (
        all_methods.groupby(
            ["contrast_id", "level", "interaction_key", "sender", "receiver", "ligand", "receptor"],
            dropna=False,
        )
        .agg(
            n_methods=("method", pd.Series.nunique),
            methods=("method", _method_list),
            mean_delta_score=("method_score", "mean"),
            max_delta_score=("method_score", "max"),
            min_delta_score=("method_score", "min"),
        )
        .reset_index()
    )
    summary["consensus_label"] = pd.cut(
        summary["n_methods"],
        bins=[-1, 0, 1, 2, 99],
        labels=["none", "single_method", "two_method_consensus", "three_method_consensus"],
    ).astype(str)
    return summary.sort_values(["n_methods", "mean_delta_score"], ascending=[False, False]).reset_index(drop=True)


def run_collect_consensus(args: argparse.Namespace) -> dict[str, object]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    liana_df = _load_delta_tables(args.liana_dir, "**/liana_delta.tsv.gz", "liana")
    cpdb_df = _load_delta_tables(args.cellphonedb_dir, "**/cellphonedb_delta.tsv.gz", "cellphonedb")
    cellchat_df = _load_delta_tables(args.cellchat_dir, "**/cellchat_delta.tsv.gz", "cellchat")

    all_methods = pd.concat([df for df in [liana_df, cpdb_df, cellchat_df] if not df.empty], ignore_index=True) if any(not df.empty for df in [liana_df, cpdb_df, cellchat_df]) else pd.DataFrame()
    consensus = _summarize_consensus(all_methods)
    top_consensus = consensus.head(args.top_n) if not consensus.empty else consensus

    write_table(args.out_dir / "all_method_delta_rows.tsv.gz", all_methods)
    write_table(args.out_dir / "consensus_interactions.tsv.gz", consensus)
    write_table(args.out_dir / "consensus_interactions_top.tsv", top_consensus)

    summary = {
        "status": "ok",
        "n_liana_rows": int(liana_df.shape[0]),
        "n_cellphonedb_rows": int(cpdb_df.shape[0]),
        "n_cellchat_rows": int(cellchat_df.shape[0]),
        "n_consensus_rows": int(consensus.shape[0]),
        "n_two_method": int(consensus["n_methods"].ge(2).sum()) if not consensus.empty else 0,
        "n_three_method": int(consensus["n_methods"].ge(3).sum()) if not consensus.empty else 0,
    }
    write_json(args.out_dir / "consensus_summary.json", summary)
    return summary


def main() -> None:
    args = parse_args()
    summary = run_collect_consensus(args)
    print(
        "[non_unified_airway][05d] "
        f"consensus_rows={summary.get('n_consensus_rows', 0)} two_method={summary.get('n_two_method', 0)}"
    )


if __name__ == "__main__":
    main()
