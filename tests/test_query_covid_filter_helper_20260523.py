from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TEST_DIR = Path(__file__).resolve().parent
PY_DIR = TEST_DIR.parent
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))

from query_covid_filter_helper_20260523 import (  # noqa: E402
    build_covid_filter_debug_frame,
    build_removed_kept_value_counts,
)


def _mock_obs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "disease": ["COVID-19", "normal", "Unknown", "MONDO_0004514", "Unknown"],
            "COVID_status": ["Unknown", "Post-COVID", "COVID+", "Healthy", "Unknown"],
            "dataset": [
                "healthy_reference",
                "healthy_reference",
                "healthy_reference",
                "healthy_reference",
                "kerstin_covid_cohort",
            ],
        },
        index=[f"cell_{i}" for i in range(5)],
    )


def test_primary_and_secondary_matches_are_removed() -> None:
    obs = _mock_obs()
    debug = build_covid_filter_debug_frame(obs)

    assert bool(debug.loc["cell_0", "remove_covid_related"]) is True
    assert debug.loc["cell_0", "trigger_source"] == "primary"

    assert bool(debug.loc["cell_1", "remove_covid_related"]) is True
    assert debug.loc["cell_1", "trigger_source"] == "secondary"

    assert bool(debug.loc["cell_2", "remove_covid_related"]) is True
    assert debug.loc["cell_2", "trigger_source"] == "secondary"


def test_non_covid_rows_are_kept_without_dataset_fallback() -> None:
    obs = _mock_obs()
    debug = build_covid_filter_debug_frame(obs, allow_dataset_fallback=False)

    assert bool(debug.loc["cell_3", "remove_covid_related"]) is False
    assert debug.loc["cell_3", "trigger_source"] == "kept"

    assert bool(debug.loc["cell_4", "remove_covid_related"]) is False
    assert debug.loc["cell_4", "trigger_source"] == "kept"


def test_dataset_fallback_only_applies_when_enabled() -> None:
    obs = _mock_obs()
    debug = build_covid_filter_debug_frame(obs, allow_dataset_fallback=True)

    assert bool(debug.loc["cell_4", "remove_covid_related"]) is True
    assert debug.loc["cell_4", "trigger_source"] == "dataset_fallback"


def test_removed_kept_value_counts_are_split_correctly() -> None:
    obs = _mock_obs()
    debug = build_covid_filter_debug_frame(obs, allow_dataset_fallback=True)
    tables = build_removed_kept_value_counts(obs, debug["remove_covid_related"], columns=["disease", "COVID_status"])

    disease_counts = tables["disease"]
    assert int(disease_counts.loc["COVID-19", "removed"]) == 1
    assert int(disease_counts.loc["MONDO_0004514", "kept"]) == 1

    status_counts = tables["COVID_status"]
    assert int(status_counts.loc["Post-COVID", "removed"]) == 1
    assert int(status_counts.loc["Healthy", "kept"]) == 1
