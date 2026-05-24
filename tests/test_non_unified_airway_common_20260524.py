from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, "/home/h2048/script/py")

from non_unified_airway.common import (  # noqa: E402
    add_canonical_fields,
    build_sample_manifest,
    discover_lineage_output_dirs,
    recommend_model_formulas,
)


def test_discover_lineage_output_dirs_parses_methods_table_paths() -> None:
    methods_text = """
| Branch | Final R output directory |
|---|---|
| B cell | `data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415` |
| T/NK | `data/R/0407/tnk_tissue_comparison_v2_6_0` |
| Epithelial | `data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun` |
"""
    paths = discover_lineage_output_dirs(methods_text, repo_root=Path("/home/h2048"))

    assert paths == [
        Path("/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415"),
        Path("/home/h2048/data/R/0407/tnk_tissue_comparison_v2_6_0"),
        Path("/home/h2048/data/R/0415/epithelial_tissue_comparison_v1_3_2_20260415_full_rerun"),
    ]


def test_add_canonical_fields_builds_site_and_disease_groups() -> None:
    meta = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4"],
            "dataset": ["D1", "D1", "D2", "D2"],
            "tissue": ["nose", "respiratory airway", "sinus", "lung parenchyma"],
            "tissue_level_2": ["inferior turbinate", "left main bronchus", "unknown", "parenchyma lower lobe"],
            "condition": ["Control", "healthy", "Case", "Control"],
            "disease_level_2": ["", "Healthy", "CRS", ""],
            "cell_type_L2": ["Epithelial", "Epithelial", "Epithelial", "Epithelial"],
            "cell_type_L3": ["Basal", "Basal", "Secretory", "AT2"],
        },
        index=["c1", "c2", "c3", "c4"],
    )

    annotated = add_canonical_fields(meta)

    assert annotated["site_group"].tolist() == ["nasal", "bronchial", "sinus", "distal_lung"]
    assert annotated["analysis_compartment"].tolist() == [
        "upper_airway",
        "lower_airway",
        "upper_airway",
        "lung_parenchyma",
    ]
    assert annotated["disease_group"].tolist() == ["healthy", "healthy", "crs", "healthy"]
    assert annotated["is_healthy"].tolist() == [True, True, False, True]
    assert annotated["tissue_site_disease"].tolist() == [
        "nose__nasal__healthy",
        "respiratory airway__bronchial__healthy",
        "sinus__sinus__crs",
        "lung parenchyma__distal_lung__healthy",
    ]


def test_recommend_model_formulas_prefers_dataset_adjusted_models_when_estimable() -> None:
    meta = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4", "S5", "S6"],
            "dataset": ["D1", "D1", "D2", "D2", "D1", "D2"],
            "tissue": ["nose", "respiratory airway", "nose", "respiratory airway", "nose", "nose"],
            "tissue_level_2": ["inferior turbinate", "left main bronchus", "inferior turbinate", "left main bronchus", "inferior turbinate", "inferior turbinate"],
            "condition": ["healthy", "healthy", "healthy", "healthy", "case", "case"],
            "disease_level_2": ["healthy", "healthy", "healthy", "healthy", "CRS", "CRS"],
        }
    )
    annotated = add_canonical_fields(meta)
    sample_manifest = build_sample_manifest(annotated)
    plan = recommend_model_formulas(sample_manifest)

    assert plan["healthy_mainline"]["status"] == "estimable"
    assert plan["healthy_mainline"]["formula"] == "~ dataset + site_group"
    assert plan["disease_extension"]["status"] == "estimable"
    assert plan["disease_extension"]["formula"] == "~ dataset + site_group + disease_group"


def test_recommend_model_formulas_downgrades_when_site_is_fully_confounded() -> None:
    meta = pd.DataFrame(
        {
            "sample": ["S1", "S2", "S3", "S4"],
            "dataset": ["D1", "D1", "D2", "D2"],
            "tissue": ["nose", "nose", "respiratory airway", "respiratory airway"],
            "tissue_level_2": ["inferior turbinate", "inferior turbinate", "left main bronchus", "left main bronchus"],
            "condition": ["healthy", "healthy", "healthy", "healthy"],
            "disease_level_2": ["healthy", "healthy", "healthy", "healthy"],
        }
    )
    annotated = add_canonical_fields(meta)
    sample_manifest = build_sample_manifest(annotated)
    plan = recommend_model_formulas(sample_manifest)

    assert plan["healthy_mainline"]["status"] == "downgraded"
    assert plan["healthy_mainline"]["formula"] == "~ tissue_site_disease"


if __name__ == "__main__":
    test_discover_lineage_output_dirs_parses_methods_table_paths()
    test_add_canonical_fields_builds_site_and_disease_groups()
    test_recommend_model_formulas_prefers_dataset_adjusted_models_when_estimable()
    test_recommend_model_formulas_downgrades_when_site_is_fully_confounded()
    print("non-unified airway common tests passed")
