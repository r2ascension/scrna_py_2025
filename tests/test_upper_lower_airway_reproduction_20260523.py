from __future__ import annotations

from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/h2048/script/py")

from upper_lower_airway_reproduction_20260523 import (  # noqa: E402
    add_site_annotations,
    normalize_site,
    read_obs_columns_h5py,
    summarize_deseq2_tree,
    summarize_metadata,
)


def _write_categorical(group: h5py.Group, name: str, values: list[str]) -> None:
    categories = list(dict.fromkeys(values))
    codes = np.array([categories.index(v) for v in values], dtype=np.int8)
    cat_group = group.create_group(name)
    cat_group.create_dataset("codes", data=codes)
    cat_group.create_dataset("categories", data=np.array(categories, dtype="S"))


def test_normalize_site_maps_known_respiratory_tissues() -> None:
    assert normalize_site("nose", "inferior turbinate").comparison_group == "nasal"
    assert normalize_site("sinus", "Unknown").site_group == "sinus"
    assert normalize_site("sinus", "Unknown").comparison_group == "nasal"

    bronchial = normalize_site("respiratory airway", "left or right main bronchus")
    assert bronchial.site_axis == "lower_conducting_airway"
    assert bronchial.comparison_group == "bronchial"

    lung = normalize_site("lung parenchyma", "parenchyma lower lobe")
    assert lung.site_axis == "distal_lung"
    assert lung.comparison_group == "lung"

    unknown = normalize_site("blood", "Unknown")
    assert unknown.site_group == "unknown"
    assert unknown.include_in_main_airway is False


def test_read_obs_columns_h5py_decodes_h5ad_categorical_and_string_columns() -> None:
    tmpdir = tempfile.TemporaryDirectory()
    tmp_path = Path(tmpdir.name)
    h5ad_path = tmp_path / "mini.h5ad"
    with h5py.File(h5ad_path, "w") as handle:
        obs = handle.create_group("obs")
        obs.create_dataset("_index", data=np.array(["cell1", "cell2", "cell3"], dtype="S"))
        _write_categorical(obs, "tissue", ["nose", "respiratory airway", "lung parenchyma"])
        obs.create_dataset("sample", data=np.array(["S1", "S2", "S3"], dtype="S"))
        _write_categorical(obs, "cell_type_L3", ["Basal", "Basal", "AT2"])

    meta = read_obs_columns_h5py(h5ad_path, ["tissue", "sample", "cell_type_L3", "missing"])

    assert list(meta.index) == ["cell1", "cell2", "cell3"]
    assert meta["tissue"].tolist() == ["nose", "respiratory airway", "lung parenchyma"]
    assert meta["sample"].tolist() == ["S1", "S2", "S3"]
    assert meta["cell_type_L3"].tolist() == ["Basal", "Basal", "AT2"]
    assert "missing" not in meta.columns
    tmpdir.cleanup()


def test_metadata_annotation_and_summary_keep_upper_lower_groups() -> None:
    meta = pd.DataFrame(
        {
            "tissue": ["nose", "sinus", "respiratory airway", "lung parenchyma"],
            "tissue_level_2": ["inferior turbinate", "Unknown", "trachea", "parenchyma upper lobe"],
            "sample": ["S1", "S1", "S2", "S3"],
            "dataset": ["D1", "D1", "D2", "D3"],
            "cell_type_L3": ["Basal", "Basal", "Basal", "AT2"],
        },
        index=["c1", "c2", "c3", "c4"],
    )

    annotated = add_site_annotations(meta)
    assert annotated["comparison_group"].tolist() == ["nasal", "nasal", "bronchial", "lung"]
    assert annotated["site_axis"].tolist() == [
        "upper_airway",
        "upper_airway",
        "lower_conducting_airway",
        "distal_lung",
    ]

    summaries = summarize_metadata(annotated, celltype_col="cell_type_L3")
    assert set(summaries) == {
        "cell_counts_by_site",
        "celltype_comparison_summary",
        "sample_site_summary",
        "celltype_site_summary",
        "celltype_sample_site_counts",
    }
    celltype_comparison = summaries["celltype_comparison_summary"]
    basal_nasal = celltype_comparison[
        (celltype_comparison["cell_type"] == "Basal")
        & (celltype_comparison["comparison_group"] == "nasal")
    ]
    assert basal_nasal["n_cells"].iloc[0] == 2


def test_summarize_deseq2_tree_extracts_contrasts_and_key_genes() -> None:
    tmpdir = tempfile.TemporaryDirectory()
    tmp_path = Path(tmpdir.name)
    result_dir = tmp_path / "pseudobulk_de_L3" / "Ciliated_Mature" / "respiratory_airway_vs_nose"
    result_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "gene": ["ACE2", "ENSG00000141510", "MUC5AC"],
            "baseMean": [12.0, 200.0, 100.0],
            "log2FoldChange": [1.5, -0.25, -2.2],
            "pvalue": [0.001, 0.5, 0.0001],
            "padj": [0.02, 0.8, 0.001],
            "sig": ["sig", "ns", "sig"],
            "direction": ["up", "ns", "down"],
        }
    ).to_csv(result_dir / "DESeq2_results.csv", index=False)

    summary, key_genes = summarize_deseq2_tree(
        tmp_path / "pseudobulk_de_L3",
        level="L3",
        gene_map={"ENSG00000141510": "TP53"},
        key_genes=["ACE2", "TP53"],
    )

    assert summary.shape[0] == 1
    row = summary.iloc[0]
    assert row["cell_type"] == "Ciliated_Mature"
    assert row["contrast"] == "respiratory_airway_vs_nose"
    assert row["n_sig"] == 2
    assert row["n_up"] == 1
    assert row["n_down"] == 1

    assert set(key_genes["gene_symbol"]) == {"ACE2", "TP53"}
    tp53 = key_genes[key_genes["gene_symbol"] == "TP53"].iloc[0]
    assert tp53["gene"] == "ENSG00000141510"
    assert tp53["is_key_gene"] is True or bool(tp53["is_key_gene"]) is True
    tmpdir.cleanup()


if __name__ == "__main__":
    test_normalize_site_maps_known_respiratory_tissues()
    test_read_obs_columns_h5py_decodes_h5ad_categorical_and_string_columns()
    test_metadata_annotation_and_summary_keep_upper_lower_groups()
    test_summarize_deseq2_tree_extracts_contrasts_and_key_genes()
    print("upper/lower airway reproduction tests passed")
