from __future__ import annotations

import pandas as pd

from non_unified_airway.communication import (
    apply_communication_contrast,
    build_communication_contrasts,
    build_sender_receiver_pairs,
    sanitize_comm_celltypes,
    summarize_celltype_support,
)
from non_unified_airway.common import build_sample_manifest


def _mock_metadata() -> pd.DataFrame:
    rows = []
    for site_group, compartment in [("nasal", "upper_airway"), ("bronchial", "lower_airway")]:
        for sample_idx in range(3):
            sample = f"{site_group}_healthy_{sample_idx}"
            for ct in ["Basal", "Secretory"]:
                for cell_idx in range(20):
                    rows.append(
                        {
                            "cell_id": f"{sample}_{ct}_{cell_idx}",
                            "sample": sample,
                            "dataset": "core",
                            "tissue": site_group,
                            "site": site_group,
                            "site_group": site_group,
                            "disease": "healthy",
                            "disease_group": "healthy",
                            "analysis_compartment": compartment,
                            "cell_type_L2": ct,
                            "cell_type_L3": f"{ct}_L3",
                        }
                    )
    for sample_idx in range(3):
        sample = f"nasal_crs_{sample_idx}"
        for ct in ["Basal", "Secretory"]:
            for cell_idx in range(20):
                rows.append(
                    {
                        "cell_id": f"{sample}_{ct}_{cell_idx}",
                        "sample": sample,
                        "dataset": "core",
                        "tissue": "nasal",
                        "site": "nasal",
                        "site_group": "nasal",
                        "disease": "CRS",
                        "disease_group": "CRS",
                        "analysis_compartment": "upper_airway",
                        "cell_type_L2": ct,
                        "cell_type_L3": f"{ct}_L3",
                    }
                )
    metadata = pd.DataFrame(rows)
    metadata["tissue_site_disease"] = (
        metadata["tissue"].astype(str)
        + "__"
        + metadata["site_group"].astype(str)
        + "__"
        + metadata["disease_group"].astype(str)
    )
    return metadata


def test_build_communication_contrasts_creates_site_and_disease_comparisons() -> None:
    metadata = _mock_metadata()
    sample_manifest = build_sample_manifest(metadata)
    contrasts = build_communication_contrasts(sample_manifest)
    contrast_ids = set(contrasts["contrast_id"])
    assert "healthy_upper_vs_bronchial" in contrast_ids
    assert "upper_airway_crs_vs_healthy" in contrast_ids
    assert contrasts.loc[contrasts["contrast_id"] == "healthy_upper_vs_bronchial", "status"].iloc[0] == "ready"


def test_apply_communication_contrast_marks_left_and_right_cells() -> None:
    metadata = _mock_metadata()
    sample_manifest = build_sample_manifest(metadata)
    contrast = build_communication_contrasts(sample_manifest)
    row = contrast.loc[contrast["contrast_id"] == "healthy_upper_vs_bronchial"].iloc[0].to_dict()
    subset = apply_communication_contrast(metadata, row)
    assert subset["contrast_side"].isin(["left", "right"]).all()
    assert subset.loc[subset["contrast_side"] == "left", "site_group"].eq("nasal").all()
    assert subset.loc[subset["contrast_side"] == "right", "site_group"].eq("bronchial").all()


def test_sanitize_and_summarize_celltype_support_filters_low_confidence_labels() -> None:
    metadata = _mock_metadata()
    metadata.loc[metadata.index[:10], "cell_type_L2"] = "Unknown"
    sample_manifest = build_sample_manifest(metadata)
    contrasts = build_communication_contrasts(sample_manifest)
    row = contrasts.loc[contrasts["contrast_id"] == "healthy_upper_vs_bronchial"].iloc[0].to_dict()
    subset = apply_communication_contrast(metadata, row)
    subset["comm_celltype"] = sanitize_comm_celltypes(subset["cell_type_L2"])
    subset = subset.loc[subset["comm_celltype"].astype(str).str.len() > 0].copy()
    summary, sample_counts = summarize_celltype_support(
        subset,
        celltype_col="comm_celltype",
        min_cells_per_sample=10,
        min_samples_per_group=3,
    )
    assert "Unknown" not in set(summary["comm_celltype"])
    assert sample_counts["n_cells"].ge(10).all()
    assert summary["eligible"].any()


def test_build_sender_receiver_pairs_supports_self_pairs() -> None:
    pair_df = build_sender_receiver_pairs(["Basal", "Secretory"], include_self=True)
    assert pair_df.shape[0] == 4
    assert pair_df["self_interaction"].sum() == 2
