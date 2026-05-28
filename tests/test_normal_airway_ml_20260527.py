from __future__ import annotations

import copy
import sys
from pathlib import Path

import pandas as pd

PY_DIR = Path(__file__).resolve().parents[1]
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))

from normal_airway_ml_common_20260527 import default_config, make_toy_airway_adata, prepare_obs_contract
from normal_airway_ml_feature_export_20260527 import run_feature_export
from normal_airway_ml_smoke_20260527 import run_smoke
from normal_airway_ml_train_sample_level_20260527 import run_training


def tiny_cfg() -> dict:
    cfg = copy.deepcopy(default_config())
    cfg["smoke"]["toy_healthy_samples_per_site"] = 2
    cfg["smoke"]["toy_disease_samples_per_site"] = 1
    cfg["smoke"]["toy_cells_per_sample"] = 18
    cfg["smoke"]["toy_n_genes"] = 60
    cfg["smoke"]["toy_latent_dim"] = 8
    cfg["feature_export"]["latent_max_dims"] = 8
    cfg["feature_export"]["pseudobulk_top_genes"] = 40
    cfg["feature_export"]["pseudobulk_n_pcs"] = 3
    cfg["train"]["n_splits"] = 3
    cfg["train"]["permutation_repeats"] = 3
    cfg["train"]["model_names"] = ["logistic_regression", "random_forest"]
    return cfg


def test_prepare_obs_contract_normalizes_sites_and_healthy() -> None:
    cfg = tiny_cfg()
    adata = make_toy_airway_adata(cfg, seed=123)
    prepared, contract = prepare_obs_contract(adata.obs.copy(), list(adata.obsm.keys()), cfg)
    healthy = prepared.loc[contract["analysis_mask"]]
    assert contract["sample_col"] == "sample"
    assert contract["condition_col"] == "condition"
    assert contract["cell_type_col"] == "cell_type_L2"
    assert contract["latent_key"] in {"X_scvi", "X_scanvi"}
    assert set(healthy["site_label"].unique()) == {"nasal", "sinus", "bronchus", "lung_parenchyma"}
    assert healthy["condition_resolved"].eq("Healthy").all()


def test_feature_export_and_training_toy(tmp_path: Path) -> None:
    cfg = tiny_cfg()
    adata = make_toy_airway_adata(cfg, seed=321)
    feature_root = tmp_path / "feature_run"
    feature_manifest = run_feature_export(adata=adata, cfg=cfg, output_dir=feature_root)
    assert feature_manifest["n_samples"] > 0
    feature_table = pd.read_csv(feature_root / "features" / "sample_features_wide.tsv.gz", sep="\t")
    assert "site_label" in feature_table.columns
    assert any(col.startswith("composition_fraction__") for col in feature_table.columns)

    train_root = tmp_path / "train_run"
    model_manifest = run_training(feature_dir=feature_root / "features", cfg=cfg, output_dir=train_root)
    assert set(model_manifest["model_names"]) == {"logistic_regression", "random_forest"}
    perf = pd.read_csv(train_root / "model" / "performance_summary.tsv", sep="\t")
    assert not perf.empty
    assert perf["status"].eq("ok").any()


def test_smoke_toy_end_to_end(tmp_path: Path) -> None:
    cfg = tiny_cfg()
    smoke_root = tmp_path / "smoke_run"
    summary = run_smoke(smoke_mode="toy", cfg=cfg, output_dir=smoke_root)
    assert summary["smoke_mode"] == "toy"
    assert Path(summary["audit_summary_json"]).exists()
    assert Path(summary["feature_manifest_json"]).exists()
    assert Path(summary["model_manifest_json"]).exists()
