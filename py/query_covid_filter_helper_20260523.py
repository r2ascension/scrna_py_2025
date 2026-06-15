from __future__ import annotations

import re
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

DEFAULT_COVID_PATTERNS: tuple[str, ...] = (
    r"\bcovid(?:-?19)?\b",
    r"\bpost[ -]?covid\b",
    r"\bsars[ -]?cov[ -]?2\b",
)
DEFAULT_PRIMARY_COLS: tuple[str, ...] = ("disease",)
DEFAULT_SECONDARY_COLS: tuple[str, ...] = ("COVID_status",)
DEFAULT_DATASET_COLS: tuple[str, ...] = ("dataset",)

_MISSING_LIKE = {
    "",
    "unknown",
    "na",
    "n/a",
    "none",
    "nan",
    "<na>",
    "not available",
    "not_applicable",
    "missing",
}


def _stringify_series(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def _combine_columns(obs: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    present = [col for col in columns if col in obs.columns]
    if not present:
        return pd.Series(pd.array([pd.NA] * len(obs), dtype="string"), index=obs.index)

    combined = pd.Series(pd.array([pd.NA] * len(obs), dtype="string"), index=obs.index)
    for col in present:
        values = _stringify_series(obs[col])
        missing = values.isna() | values.str.lower().isin(_MISSING_LIKE)
        values = values.mask(missing, pd.NA)
        if combined.isna().all():
            combined = values
            continue
        combined = combined.fillna("")
        values = values.fillna("")
        combined = (combined + " | " + values).str.strip(" |")
        combined = combined.replace("", pd.NA).astype("string")
    return combined.astype("string")


def _match_regex(
    series: pd.Series,
    regex: re.Pattern[str],
    *,
    normalize_separators: bool = False,
) -> pd.Series:
    values = series.fillna("")
    if normalize_separators:
        values = values.str.replace(r"[_./]+", " ", regex=True)
    return values.str.contains(regex, regex=True, na=False)


def _missing_like(series: pd.Series) -> pd.Series:
    normalized = series.fillna("").str.strip().str.lower()
    return normalized.isin(_MISSING_LIKE)


def build_covid_filter_debug_frame(
    obs: pd.DataFrame,
    *,
    primary_cols: Sequence[str] = DEFAULT_PRIMARY_COLS,
    secondary_cols: Sequence[str] = DEFAULT_SECONDARY_COLS,
    dataset_cols: Sequence[str] = DEFAULT_DATASET_COLS,
    patterns: Sequence[str] = DEFAULT_COVID_PATTERNS,
    allow_dataset_fallback: bool = False,
) -> pd.DataFrame:
    """Build a row-wise audit frame for removing COVID-related cells.

    Removal policy:
    1. Remove rows whose *primary* metadata (default: ``disease``) matches.
    2. Also remove rows whose *secondary* metadata (default: ``COVID_status``)
       matches. This catches cases like ``disease='normal'`` but
       ``COVID_status='Post-COVID'``.
    3. Optionally use dataset-name fallback only when both primary and secondary
       metadata are missing-like.
    """
    if obs.index.has_duplicates:
        raise ValueError("obs index must be unique when building COVID filter audit frame")

    regex = re.compile("|".join(f"(?:{pattern})" for pattern in patterns), flags=re.IGNORECASE)
    primary_text = _combine_columns(obs, primary_cols)
    secondary_text = _combine_columns(obs, secondary_cols)
    dataset_text = _combine_columns(obs, dataset_cols)

    matched_primary = _match_regex(primary_text, regex)
    matched_secondary = _match_regex(secondary_text, regex)
    matched_dataset = (
        _match_regex(dataset_text, regex, normalize_separators=True)
        if allow_dataset_fallback
        else pd.Series(False, index=obs.index)
    )

    primary_missing = _missing_like(primary_text)
    secondary_missing = _missing_like(secondary_text)
    dataset_fallback = allow_dataset_fallback & primary_missing & secondary_missing & matched_dataset

    remove_mask = matched_primary | matched_secondary | dataset_fallback
    trigger_source = np.where(
        matched_primary,
        "primary",
        np.where(matched_secondary, "secondary", np.where(dataset_fallback, "dataset_fallback", "kept")),
    )

    return pd.DataFrame(
        {
            "primary_text": primary_text,
            "secondary_text": secondary_text,
            "dataset_text": dataset_text,
            "primary_missing": primary_missing.astype(bool),
            "secondary_missing": secondary_missing.astype(bool),
            "matched_primary": matched_primary.astype(bool),
            "matched_secondary": matched_secondary.astype(bool),
            "matched_dataset": matched_dataset.astype(bool),
            "dataset_fallback": pd.Series(dataset_fallback, index=obs.index, dtype=bool),
            "remove_covid_related": remove_mask.astype(bool),
            "trigger_source": pd.Series(trigger_source, index=obs.index, dtype="string"),
        },
        index=obs.index,
    )


def build_covid_related_mask(
    obs: pd.DataFrame,
    **kwargs: object,
) -> pd.Series:
    return build_covid_filter_debug_frame(obs, **kwargs)["remove_covid_related"].astype(bool)


def summarize_covid_filter_debug_frame(debug_frame: pd.DataFrame) -> pd.DataFrame:
    required = {"trigger_source", "remove_covid_related"}
    missing = required - set(debug_frame.columns)
    if missing:
        raise KeyError(f"debug_frame missing required columns: {sorted(missing)}")

    counts = debug_frame["trigger_source"].fillna("kept").astype("string").value_counts(dropna=False)
    total = int(len(debug_frame))
    rows: list[Mapping[str, object]] = []
    for source, n_cells in counts.items():
        rows.append(
            {
                "trigger_source": str(source),
                "n_cells": int(n_cells),
                "fraction": float(n_cells / max(total, 1)),
                "removed": bool(source != "kept"),
            }
        )
    return pd.DataFrame(rows).sort_values(["removed", "n_cells"], ascending=[False, False]).reset_index(drop=True)


def build_removed_kept_value_counts(
    obs: pd.DataFrame,
    mask: Sequence[bool] | pd.Series,
    *,
    columns: Iterable[str],
) -> dict[str, pd.DataFrame]:
    mask_s = pd.Series(mask, index=obs.index, dtype=bool)
    outputs: dict[str, pd.DataFrame] = {}
    for col in columns:
        if col not in obs.columns:
            continue
        values = _stringify_series(obs[col]).fillna("<NA>")
        counts = pd.DataFrame(
            {
                "all": values.value_counts(dropna=False),
                "removed": values[mask_s].value_counts(dropna=False),
                "kept": values[~mask_s].value_counts(dropna=False),
            }
        ).fillna(0).astype(int)
        counts.index.name = col
        outputs[col] = counts.sort_values(["removed", "all"], ascending=[False, False])
    return outputs
