from __future__ import annotations

from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import pandas as pd

PY_DIR = Path(__file__).resolve().parents[1] / "py"
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))

from upper_lower_airway_reproduction_20260523 import (  # noqa: E402
    annotate_pairwise_columns,
    add_site_annotations,
    build_deg_metric_matrix,
    build_direction_balance_matrix,
    generate_ppi_bridge_outputs,
    generate_deg_summary_visuals,
    normalize_site,
    prepare_key_gene_panel_matrices,
    read_obs_columns_h5py,
    summarize_ppi_seed_candidates,
    summarize_sig_deg_table,
    summarize_top_gene_recurrence,
    summarize_pairwise_catalog,
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
    assert normalize_site("nose", "inferior turbinate").site_detail == "nose"
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


def test_annotate_pairwise_columns_adds_harmonized_site_pairs() -> None:
    summary = pd.DataFrame(
        {
            "level": ["L3", "L3"],
            "cell_type": ["Basal", "Goblet"],
            "contrast": ["respiratory_airway_vs_nose", "sinus_vs_respiratory_airway"],
            "contrast_left": ["respiratory_airway", "sinus"],
            "contrast_right": ["nose", "respiratory_airway"],
            "status": ["ok", "ok"],
            "n_sig": [12, 34],
        }
    )

    annotated = annotate_pairwise_columns(summary)

    first = annotated.iloc[0]
    assert first["left_site_detail"] == "bronchial"
    assert first["right_site_detail"] == "nose"
    assert first["left_comparison_group"] == "bronchial"
    assert first["right_comparison_group"] == "nasal"
    assert first["contrast_site_pair"] == "bronchial_vs_nose"
    assert first["canonical_site_pair"] == "nose_vs_bronchial"
    assert first["canonical_comparison_pair"] == "nasal_vs_bronchial"

    second = annotated.iloc[1]
    assert second["left_site_detail"] == "sinus"
    assert second["right_site_detail"] == "bronchial"
    assert second["canonical_site_pair"] == "sinus_vs_bronchial"


def test_summarize_pairwise_catalog_collapses_all_pairwise_groups() -> None:
    annotated = pd.DataFrame(
        {
            "level": ["L3", "L3", "L3"],
            "cell_type": ["Basal", "Goblet", "AT2"],
            "contrast": [
                "respiratory_airway_vs_nose",
                "sinus_vs_respiratory_airway",
                "nose_vs_lung_parenchyma",
            ],
            "status": ["ok", "ok", "ok"],
            "n_sig": [12, 34, 56],
            "canonical_site_pair": [
                "nose_vs_bronchial",
                "sinus_vs_bronchial",
                "nose_vs_distal_lung",
            ],
            "canonical_comparison_pair": [
                "nasal_vs_bronchial",
                "nasal_vs_bronchial",
                "nasal_vs_lung",
            ],
        }
    )

    detailed = summarize_pairwise_catalog(annotated, "canonical_site_pair")
    grouped = summarize_pairwise_catalog(annotated, "canonical_comparison_pair")

    assert detailed["pair_label"].tolist() == [
        "nose_vs_bronchial",
        "nose_vs_distal_lung",
        "sinus_vs_bronchial",
    ]

    nasal_bronchial = grouped[grouped["pair_label"] == "nasal_vs_bronchial"].iloc[0]
    assert nasal_bronchial["n_contrast_files"] == 2
    assert nasal_bronchial["n_cell_types"] == 2
    assert nasal_bronchial["n_sig_total"] == 46


def test_build_deg_matrices_and_gene_recurrence() -> None:
    de_summary = pd.DataFrame(
        {
            "level": ["L3", "L3", "L3", "L3"],
            "cell_type": ["Basal", "Basal", "Goblet", "Goblet"],
            "status": ["ok", "ok", "ok", "ok"],
            "canonical_site_pair": [
                "nose_vs_bronchial",
                "sinus_vs_bronchial",
                "nose_vs_bronchial",
                "nose_vs_distal_lung",
            ],
            "n_sig": [20, 8, 10, 4],
            "n_up": [15, 2, 3, 1],
            "n_down": [5, 6, 7, 3],
            "top_up_genes": ["KRT5;KRT17", "CXCL8", "MUC1;KRT5", "SFTPA1"],
            "top_down_genes": ["SCGB1A1", "KRT5", "BPIFA1;SCGB1A1", "NAPSA"],
        }
    )

    metric_matrix = build_deg_metric_matrix(de_summary, "n_sig", level="L3")
    assert metric_matrix.index.tolist() == ["Basal", "Goblet"]
    assert metric_matrix.columns.tolist() == ["nose_vs_bronchial", "nose_vs_distal_lung", "sinus_vs_bronchial"]
    assert metric_matrix.loc["Basal", "nose_vs_bronchial"] == 20
    assert metric_matrix.loc["Goblet", "nose_vs_distal_lung"] == 4

    balance_matrix = build_direction_balance_matrix(de_summary, level="L3")
    assert np.isclose(balance_matrix.loc["Basal", "nose_vs_bronchial"], 0.5)
    assert np.isclose(balance_matrix.loc["Basal", "sinus_vs_bronchial"], -0.5)
    assert np.isclose(balance_matrix.loc["Goblet", "nose_vs_bronchial"], -0.4)

    recurrence_up = summarize_top_gene_recurrence(de_summary, direction="up", level="L3")
    krt5 = recurrence_up[(recurrence_up["pair_label"] == "nose_vs_bronchial") & (recurrence_up["gene_symbol"] == "KRT5")].iloc[0]
    assert krt5["n_cell_types"] == 2
    assert np.isclose(krt5["weighted_rank_score"], 1.5)


def test_prepare_key_gene_panel_matrices_and_visual_smoke() -> None:
    key_gene_de = pd.DataFrame(
        {
            "level": ["L3", "L3", "L3", "L3"],
            "cell_type": ["Basal", "Goblet", "Basal", "Goblet"],
            "canonical_site_pair": [
                "nose_vs_bronchial",
                "nose_vs_bronchial",
                "nose_vs_distal_lung",
                "nose_vs_distal_lung",
            ],
            "gene_symbol": ["KRT5", "MUC1", "KRT5", "SFTPA1"],
            "log2FoldChange": [1.2, -0.8, 0.4, -1.6],
            "padj": [0.01, 0.02, 0.2, 0.001],
        }
    )
    matrices = prepare_key_gene_panel_matrices(key_gene_de, level="L3", celltype_order=["Basal", "Goblet"])
    assert sorted(matrices) == ["nose_vs_bronchial", "nose_vs_distal_lung"]
    assert matrices["nose_vs_bronchial"].columns.tolist() == ["Basal", "Goblet"]
    assert np.isclose(matrices["nose_vs_bronchial"].loc["KRT5", "Basal"], 1.2)

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return

    de_summary_l2 = pd.DataFrame(
        {
            "level": ["L2"],
            "cell_type": ["Epithelial"],
            "status": ["ok"],
            "canonical_site_pair": ["nose_vs_bronchial"],
            "n_sig": [12],
            "n_up": [8],
            "n_down": [4],
            "top_up_genes": ["KRT5;KRT17"],
            "top_down_genes": ["SCGB1A1;BPIFA1"],
        }
    )
    de_summary_l3 = pd.DataFrame(
        {
            "level": ["L3", "L3"],
            "cell_type": ["Basal", "Goblet"],
            "status": ["ok", "ok"],
            "canonical_site_pair": ["nose_vs_bronchial", "nose_vs_distal_lung"],
            "n_sig": [20, 9],
            "n_up": [15, 2],
            "n_down": [5, 7],
            "top_up_genes": ["KRT5;KRT17", "MUC1"],
            "top_down_genes": ["SCGB1A1", "BPIFA1;SFTPA1"],
        }
    )
    pairwise_l2 = pd.DataFrame(
        {
            "level": ["L2"],
            "pair_label": ["nose_vs_bronchial"],
            "n_contrast_files": [1],
            "n_cell_types": [1],
            "n_sig_total": [12],
            "median_n_sig": [12.0],
        }
    )
    pairwise_l3 = pd.DataFrame(
        {
            "level": ["L3", "L3"],
            "pair_label": ["nose_vs_bronchial", "nose_vs_distal_lung"],
            "n_contrast_files": [1, 1],
            "n_cell_types": [1, 1],
            "n_sig_total": [20, 9],
            "median_n_sig": [20.0, 9.0],
        }
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        outputs, manifest = generate_deg_summary_visuals(
            Path(tmpdir),
            de_summary_l2=de_summary_l2,
            de_summary_l3=de_summary_l3,
            de_pairwise_site_l2=pairwise_l2,
            de_pairwise_site_l3=pairwise_l3,
            key_gene_de=key_gene_de,
        )
        assert "deg_visual_manifest" in outputs
        assert Path(outputs["deg_visual_manifest"]).exists()
        assert not manifest.empty
        expected_png = Path(tmpdir) / "figures" / "deg_l3_nsig_heatmap.png"
        assert expected_png.exists()


def test_ppi_bridge_outputs_export_ranked_seed_lists() -> None:
    tmpdir = tempfile.TemporaryDirectory()
    tmp_path = Path(tmpdir.name)
    root = tmp_path / "pseudobulk_de_L3"
    basal_dir = root / "Basal" / "respiratory_airway_vs_nose"
    goblet_dir = root / "Goblet" / "respiratory_airway_vs_nose"
    basal_dir.mkdir(parents=True)
    goblet_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "gene": ["KRT5", "ACE2", "ENSG00000141510", "MUC1"],
            "baseMean": [100.0, 40.0, 80.0, 20.0],
            "log2FoldChange": [1.8, 1.2, -1.4, 0.2],
            "pvalue": [1e-5, 2e-3, 5e-4, 0.4],
            "padj": [1e-4, 0.01, 0.002, 0.5],
            "sig": ["sig", "sig", "sig", "ns"],
            "direction": ["up", "up", "down", "ns"],
        }
    ).to_csv(basal_dir / "DESeq2_results.csv", index=False)
    pd.DataFrame(
        {
            "gene": ["KRT5", "MUC1", "BPIFA1"],
            "baseMean": [70.0, 50.0, 45.0],
            "log2FoldChange": [1.1, -1.0, -0.8],
            "pvalue": [3e-3, 2e-4, 5e-3],
            "padj": [0.02, 0.001, 0.03],
            "sig": ["sig", "sig", "sig"],
            "direction": ["up", "down", "down"],
        }
    ).to_csv(goblet_dir / "DESeq2_results.csv", index=False)

    sig_deg = summarize_sig_deg_table(root, "L3", gene_map={"ENSG00000141510": "TP53"})
    annotated = annotate_pairwise_columns(sig_deg)
    assert set(annotated["gene_symbol"]) == {"KRT5", "ACE2", "TP53", "MUC1", "BPIFA1"}
    assert set(annotated["canonical_site_pair"]) == {"nose_vs_bronchial"}
    assert annotated["ppi_priority_score"].gt(0).all()

    site_summary = summarize_ppi_seed_candidates(annotated, level="L3", pair_col="canonical_site_pair")
    krt5 = site_summary[site_summary["gene_symbol"] == "KRT5"].iloc[0]
    assert krt5["n_cell_types"] == 2
    assert krt5["dominant_direction"] == "up"
    assert krt5["seed_rank"] == 1

    out_dir = tmp_path / "ppi_out"
    outputs, manifest = generate_ppi_bridge_outputs(out_dir, ppi_deg_l3=annotated, ppi_top_n=2)
    assert Path(outputs["ppi_ready_deg_L3"]).exists()
    assert Path(outputs["ppi_seed_summary_L3_by_site"]).exists()
    assert Path(outputs["ppi_input_manifest_L3"]).exists()
    assert not manifest.empty

    pairwise_top = manifest[
        (manifest["scope"] == "l3_pairwise")
        & (manifest["pair_label"] == "nose_vs_bronchial")
        & (manifest["gene_set_kind"] == "top2")
    ].iloc[0]
    top_genes = Path(pairwise_top["txt_path"]).read_text(encoding="utf-8").strip().splitlines()
    assert "KRT5" in top_genes
    assert len(top_genes) == 2
    tmpdir.cleanup()


if __name__ == "__main__":
    test_normalize_site_maps_known_respiratory_tissues()
    test_read_obs_columns_h5py_decodes_h5ad_categorical_and_string_columns()
    test_metadata_annotation_and_summary_keep_upper_lower_groups()
    test_summarize_deseq2_tree_extracts_contrasts_and_key_genes()
    test_annotate_pairwise_columns_adds_harmonized_site_pairs()
    test_summarize_pairwise_catalog_collapses_all_pairwise_groups()
    test_build_deg_matrices_and_gene_recurrence()
    test_prepare_key_gene_panel_matrices_and_visual_smoke()
    test_ppi_bridge_outputs_export_ranked_seed_lists()
    print("upper/lower airway reproduction tests passed")
