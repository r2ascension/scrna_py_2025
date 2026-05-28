#!/usr/bin/env python3
from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from non_unified_airway.communication import DEFAULT_COMM_OUT_DIR, standardize_interaction_columns
from non_unified_airway.common import write_json, write_tsv


COMM_V2_DIR = DEFAULT_COMM_OUT_DIR.parent / "communication_lineage_native_export_v2"
DEFAULT_CELLCHAT_DIR = Path(
	"/home/h2048/data/R/0525/non_unified_airway/communication/cellchat_lineage_native_export_v2"
)
DEFAULT_LIANA_DIR = COMM_V2_DIR / "liana_20260525_v2"
DEFAULT_CELLPHONEDB_DIR = COMM_V2_DIR / "cellphonedb_20260525_v2"
DEFAULT_OUT_DIR = COMM_V2_DIR / "consensus_20260525_v2"

METHODS = ("cellchat", "liana", "cellphonedb")
INTERACTION_KEYS = ["contrast_id", "level", "interaction_key", "sender", "receiver", "ligand", "receptor"]


def write_table(path: Path, table: pd.DataFrame) -> None:
	kwargs = {"compression": "gzip"} if str(path).endswith(".gz") else {}
	write_tsv(path, table, **kwargs)


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description=(
			"Collect LIANA/CellPhoneDB/CellChat pairwise communication outputs into a richer multi-method "
			"consensus package with interaction-, sender-receiver-, and overlap-level summaries."
		)
	)
	parser.add_argument("--liana-dir", type=Path, default=DEFAULT_LIANA_DIR)
	parser.add_argument("--cellphonedb-dir", type=Path, default=DEFAULT_CELLPHONEDB_DIR)
	parser.add_argument("--cellchat-dir", type=Path, default=DEFAULT_CELLCHAT_DIR)
	parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
	parser.add_argument("--top-n", type=int, default=500)
	return parser.parse_args()


def _empty_frame(columns: Iterable[str]) -> pd.DataFrame:
	return pd.DataFrame(columns=list(columns))


def _first_numeric_column(df: pd.DataFrame, candidates: Iterable[str]) -> tuple[pd.Series, str]:
	for column in candidates:
		if column in df.columns:
			return pd.to_numeric(df[column], errors="coerce"), column
	return pd.Series(np.nan, index=df.index, dtype=float), ""


def _direction_label(score: object) -> str:
	value = pd.to_numeric(pd.Series([score]), errors="coerce").iloc[0]
	if pd.isna(value):
		return "missing"
	if value > 0:
		return "right_higher"
	if value < 0:
		return "left_higher"
	return "neutral"


def _build_interaction_label(df: pd.DataFrame) -> pd.Series:
	return (
		df["contrast_id"].astype(str)
		+ " | "
		+ df["level"].astype(str)
		+ " | "
		+ df["sender"].astype(str)
		+ " -> "
		+ df["receiver"].astype(str)
		+ " ["
		+ df["ligand"].astype(str)
		+ "-"
		+ df["receptor"].astype(str)
		+ "]"
	)


def _load_delta_tables(base_dir: Path, pattern: str, method_name: str) -> pd.DataFrame:
	rows: list[pd.DataFrame] = []
	if not base_dir.exists():
		return pd.DataFrame()

	for path in sorted(base_dir.glob(pattern)):
		try:
			raw = pd.read_csv(path, sep="\t")
		except Exception:
			continue
		if raw.empty:
			continue

		df = standardize_interaction_columns(raw)
		if df.empty:
			continue

		for column, fallback in (("contrast_id", path.parents[1].name), ("level", path.parent.name)):
			if column not in df.columns:
				df[column] = fallback
			df[column] = df[column].astype(str)

		required = ["contrast_id", "level", "interaction_key", "sender", "receiver", "ligand", "receptor"]
		if any(column not in df.columns for column in required):
			continue

		delta_score, delta_col = _first_numeric_column(
			df,
			[
				"delta_score_right_minus_left",
				"delta_score",
				"mean_delta_score",
				"score_delta",
			],
		)
		left_score, left_metric = _first_numeric_column(
			df,
			["left_score", "left_mean", "left_prob", "left_weight", "left_count"],
		)
		right_score, right_metric = _first_numeric_column(
			df,
			["right_score", "right_mean", "right_prob", "right_weight", "right_count"],
		)

		df = df.copy()
		df["method"] = method_name
		df["delta_score_right_minus_left"] = delta_score
		df["abs_delta_score"] = df["delta_score_right_minus_left"].abs()
		df["left_score_unified"] = left_score
		df["right_score_unified"] = right_score
		df["left_metric_unified"] = left_metric
		df["right_metric_unified"] = right_metric
		df["delta_metric_unified"] = delta_col or "delta_score_right_minus_left"
		df["direction"] = df["delta_score_right_minus_left"].map(_direction_label)
		df["method_result_path"] = str(path)
		df["interaction_label"] = _build_interaction_label(df)
		rows.append(df)

	if not rows:
		return pd.DataFrame()
	return pd.concat(rows, ignore_index=True, sort=False)


def _sorted_method_list(values: pd.Series) -> str:
	uniq = sorted({str(value) for value in values if str(value).strip()})
	return ";".join(uniq)


def _direction_consistency(summary: pd.DataFrame) -> pd.DataFrame:
	consistent = np.where(
		(summary["n_right_higher"] > 0) & (summary["n_left_higher"] == 0),
		"consistent_right_higher",
		np.where(
			(summary["n_left_higher"] > 0) & (summary["n_right_higher"] == 0),
			"consistent_left_higher",
			np.where(
				(summary["n_right_higher"] == 0) & (summary["n_left_higher"] == 0) & (summary["n_neutral"] > 0),
				"neutral_only",
				"discordant",
			),
		),
	)
	dominant = np.where(
		summary["n_right_higher"] > summary["n_left_higher"],
		"right_higher",
		np.where(
			summary["n_left_higher"] > summary["n_right_higher"],
			"left_higher",
			np.where(summary["n_neutral"] > 0, "neutral", "mixed"),
		),
	)
	summary["direction_consistency"] = consistent
	summary["dominant_direction"] = dominant
	return summary


def _support_tier(n_methods: pd.Series) -> pd.Series:
	return pd.cut(
		n_methods,
		bins=[-1, 0, 1, 2, 99],
		labels=["none", "single_method", "two_method_consensus", "three_method_consensus"],
	).astype(str)


def _summarize_consensus(all_methods: pd.DataFrame) -> pd.DataFrame:
	if all_methods.empty:
		return pd.DataFrame()

	summary = (
		all_methods.groupby(INTERACTION_KEYS, dropna=False)
		.agg(
			n_methods=("method", pd.Series.nunique),
			n_rows=("method", "size"),
			methods=("method", _sorted_method_list),
			mean_delta_score=("delta_score_right_minus_left", "mean"),
			median_delta_score=("delta_score_right_minus_left", "median"),
			mean_abs_delta=("abs_delta_score", "mean"),
			max_abs_delta=("abs_delta_score", "max"),
			min_abs_delta=("abs_delta_score", "min"),
			n_right_higher=("direction", lambda x: int(pd.Series(x).eq("right_higher").sum())),
			n_left_higher=("direction", lambda x: int(pd.Series(x).eq("left_higher").sum())),
			n_neutral=("direction", lambda x: int(pd.Series(x).eq("neutral").sum())),
			n_missing_direction=("direction", lambda x: int(pd.Series(x).eq("missing").sum())),
		)
		.reset_index()
	)

	score_wide = (
		all_methods.pivot_table(
			index=INTERACTION_KEYS,
			columns="method",
			values="delta_score_right_minus_left",
			aggfunc="mean",
		)
		.rename(columns=lambda method: f"{method}_delta_score")
		.reset_index()
	)

	summary = summary.merge(score_wide, on=INTERACTION_KEYS, how="left")
	summary = _direction_consistency(summary)
	summary["support_tier"] = _support_tier(summary["n_methods"])
	summary["interaction_label"] = _build_interaction_label(summary)
	summary = summary.sort_values(
		["n_methods", "direction_consistency", "max_abs_delta", "mean_abs_delta"],
		ascending=[False, True, False, False],
	).reset_index(drop=True)
	return summary


def _summarize_sender_receiver(consensus: pd.DataFrame) -> pd.DataFrame:
	if consensus.empty:
		return pd.DataFrame()

	summary = (
		consensus.groupby(["contrast_id", "level", "sender", "receiver"], dropna=False)
		.agg(
			n_consensus_interactions=("interaction_key", pd.Series.nunique),
			mean_support=("n_methods", "mean"),
			max_support=("n_methods", "max"),
			n_two_method=("n_methods", lambda x: int(pd.Series(x).ge(2).sum())),
			n_three_method=("n_methods", lambda x: int(pd.Series(x).ge(3).sum())),
			n_consistent_direction=("direction_consistency", lambda x: int(pd.Series(x).str.startswith("consistent").sum())),
			sum_abs_mean_delta=("mean_abs_delta", "sum"),
			max_abs_delta=("max_abs_delta", "max"),
		)
		.reset_index()
	)
	summary["sender_receiver_label"] = (
		summary["contrast_id"].astype(str)
		+ " | "
		+ summary["level"].astype(str)
		+ " | "
		+ summary["sender"].astype(str)
		+ " -> "
		+ summary["receiver"].astype(str)
	)
	return summary.sort_values(
		["max_support", "sum_abs_mean_delta", "n_consensus_interactions"],
		ascending=[False, False, False],
	).reset_index(drop=True)


def _pairwise_direction_map(df: pd.DataFrame) -> dict[tuple[str, str], str]:
	if df.empty:
		return {}
	grouped = (
		df.groupby(["method", "interaction_key"], dropna=False)["direction"]
		.agg(lambda x: _sorted_method_list(pd.Series(x)))
		.to_dict()
	)
	return {(str(method), str(key)): str(direction) for (method, key), direction in grouped.items()}


def _summarize_method_overlap(all_methods: pd.DataFrame) -> pd.DataFrame:
	columns = [
		"contrast_id",
		"level",
		"method_a",
		"method_b",
		"n_method_a",
		"n_method_b",
		"n_intersection",
		"n_union",
		"jaccard",
		"overlap_coef",
		"n_direction_compared",
		"n_direction_agree",
		"direction_agreement_rate",
	]
	if all_methods.empty:
		return _empty_frame(columns)

	rows: list[dict[str, object]] = []
	for (contrast_id, level), group_df in all_methods.groupby(["contrast_id", "level"], dropna=False):
		sets = {
			method: set(group_df.loc[group_df["method"].eq(method), "interaction_key"].astype(str))
			for method in METHODS
		}
		direction_map = _pairwise_direction_map(group_df.loc[group_df["direction"].ne("missing")])
		for method_a, method_b in combinations(METHODS, 2):
			set_a = sets.get(method_a, set())
			set_b = sets.get(method_b, set())
			intersection = set_a & set_b
			union = set_a | set_b
			direction_pairs = [
				(
					direction_map.get((method_a, interaction_key), ""),
					direction_map.get((method_b, interaction_key), ""),
				)
				for interaction_key in intersection
			]
			compared = [(left, right) for left, right in direction_pairs if left and right]
			n_agree = sum(1 for left, right in compared if left == right)
			rows.append(
				{
					"contrast_id": str(contrast_id),
					"level": str(level),
					"method_a": method_a,
					"method_b": method_b,
					"n_method_a": int(len(set_a)),
					"n_method_b": int(len(set_b)),
					"n_intersection": int(len(intersection)),
					"n_union": int(len(union)),
					"jaccard": float(len(intersection) / len(union)) if union else np.nan,
					"overlap_coef": float(len(intersection) / min(len(set_a), len(set_b))) if min(len(set_a), len(set_b)) else np.nan,
					"n_direction_compared": int(len(compared)),
					"n_direction_agree": int(n_agree),
					"direction_agreement_rate": float(n_agree / len(compared)) if compared else np.nan,
				}
			)

	if not rows:
		return _empty_frame(columns)
	return pd.DataFrame(rows).sort_values(
		["contrast_id", "level", "jaccard", "n_intersection"],
		ascending=[True, True, False, False],
	).reset_index(drop=True)


def _summarize_support_distribution(consensus: pd.DataFrame) -> pd.DataFrame:
	columns = ["contrast_id", "level", "n_methods", "support_tier", "direction_consistency", "n_interactions"]
	if consensus.empty:
		return _empty_frame(columns)
	return (
		consensus.groupby(["contrast_id", "level", "n_methods", "support_tier", "direction_consistency"], dropna=False)
		.size()
		.reset_index(name="n_interactions")
		.sort_values(["contrast_id", "level", "n_methods", "direction_consistency"], ascending=[True, True, False, True])
		.reset_index(drop=True)
	)


def _summarize_method_rows(all_methods: pd.DataFrame) -> pd.DataFrame:
	columns = [
		"method",
		"contrast_id",
		"level",
		"n_rows",
		"n_unique_interactions",
		"n_sender_receiver_pairs",
		"n_right_higher",
		"n_left_higher",
		"n_neutral",
		"mean_abs_delta",
		"max_abs_delta",
	]
	if all_methods.empty:
		return _empty_frame(columns)
	all_methods = all_methods.copy()
	all_methods["sender_receiver_key"] = (
		all_methods["sender"].astype(str) + "|" + all_methods["receiver"].astype(str)
	)
	summary = (
		all_methods.groupby(["method", "contrast_id", "level"], dropna=False)
		.agg(
			n_rows=("interaction_key", "size"),
			n_unique_interactions=("interaction_key", pd.Series.nunique),
			n_sender_receiver_pairs=("sender_receiver_key", pd.Series.nunique),
			n_right_higher=("direction", lambda x: int(pd.Series(x).eq("right_higher").sum())),
			n_left_higher=("direction", lambda x: int(pd.Series(x).eq("left_higher").sum())),
			n_neutral=("direction", lambda x: int(pd.Series(x).eq("neutral").sum())),
			mean_abs_delta=("abs_delta_score", "mean"),
			max_abs_delta=("abs_delta_score", "max"),
		)
		.reset_index()
	)
	return summary.sort_values(["method", "contrast_id", "level"]).reset_index(drop=True)


def run_collect_consensus(args: argparse.Namespace) -> dict[str, object]:
	args.out_dir.mkdir(parents=True, exist_ok=True)

	liana_df = _load_delta_tables(args.liana_dir, "**/liana_delta.tsv.gz", "liana")
	cellphonedb_df = _load_delta_tables(args.cellphonedb_dir, "**/cellphonedb_delta.tsv.gz", "cellphonedb")
	cellchat_df = _load_delta_tables(args.cellchat_dir, "**/cellchat_delta.tsv.gz", "cellchat")

	all_frames = [frame for frame in [cellchat_df, liana_df, cellphonedb_df] if not frame.empty]
	all_methods = pd.concat(all_frames, ignore_index=True, sort=False) if all_frames else pd.DataFrame()
	if not all_methods.empty:
		all_methods = all_methods.sort_values(["contrast_id", "level", "method", "abs_delta_score"], ascending=[True, True, True, False]).reset_index(drop=True)

	method_summary = _summarize_method_rows(all_methods)
	consensus = _summarize_consensus(all_methods)
	sender_receiver = _summarize_sender_receiver(consensus)
	overlap = _summarize_method_overlap(all_methods)
	support_distribution = _summarize_support_distribution(consensus)

	top_consensus = consensus.head(args.top_n) if not consensus.empty else consensus
	top_sender_receiver = sender_receiver.head(args.top_n) if not sender_receiver.empty else sender_receiver
	top_overlap = overlap.head(args.top_n) if not overlap.empty else overlap

	write_table(args.out_dir / "all_method_delta_rows.tsv.gz", all_methods)
	write_table(args.out_dir / "method_summary.tsv.gz", method_summary)
	write_table(args.out_dir / "consensus_interactions.tsv.gz", consensus)
	write_table(args.out_dir / "consensus_interactions_top.tsv.gz", top_consensus)
	write_table(args.out_dir / "sender_receiver_consensus.tsv.gz", sender_receiver)
	write_table(args.out_dir / "sender_receiver_consensus_top.tsv.gz", top_sender_receiver)
	write_table(args.out_dir / "method_overlap_long.tsv.gz", overlap)
	write_table(args.out_dir / "method_overlap_top.tsv.gz", top_overlap)
	write_table(args.out_dir / "support_distribution.tsv.gz", support_distribution)

	method_row_counts = {
		method: int(all_methods[all_methods["method"].eq(method)].shape[0]) if not all_methods.empty else 0
		for method in METHODS
	}
	summary = {
		"status": "ok",
		"methods_available": sorted([method for method, n_rows in method_row_counts.items() if n_rows > 0]),
		"n_cellchat_rows": method_row_counts["cellchat"],
		"n_liana_rows": method_row_counts["liana"],
		"n_cellphonedb_rows": method_row_counts["cellphonedb"],
		"n_consensus_rows": int(consensus.shape[0]),
		"n_sender_receiver_rows": int(sender_receiver.shape[0]),
		"n_overlap_rows": int(overlap.shape[0]),
		"n_two_method": int(consensus["n_methods"].ge(2).sum()) if not consensus.empty else 0,
		"n_three_method": int(consensus["n_methods"].ge(3).sum()) if not consensus.empty else 0,
		"n_consistent_direction": int(consensus["direction_consistency"].str.startswith("consistent").sum()) if not consensus.empty else 0,
		"n_discordant_direction": int(consensus["direction_consistency"].eq("discordant").sum()) if not consensus.empty else 0,
		"top_n_written": int(args.top_n),
	}
	write_json(args.out_dir / "consensus_summary.json", summary)
	return summary


def main() -> None:
	args = parse_args()
	summary = run_collect_consensus(args)
	print(
		"[non_unified_airway][05d_v2] "
		f"consensus_rows={summary.get('n_consensus_rows', 0)} "
		f"two_method={summary.get('n_two_method', 0)} "
		f"three_method={summary.get('n_three_method', 0)}"
	)


if __name__ == "__main__":
	main()
