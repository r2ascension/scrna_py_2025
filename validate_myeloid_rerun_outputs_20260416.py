#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
from pathlib import Path

import anndata as ad

PATCH_OUTPUT = Path(
    os.getenv(
        "MYELOID_VALIDATE_PATCH_OUTPUT",
        "/home/h2048/data/py/0416/adata_myeloid_L3refined_tissueaware_patched_v1.h5ad",
    )
)
R_OUTPUT_DIR = Path(
    os.getenv(
        "MYELOID_VALIDATE_R_OUTPUT_DIR",
        "/home/h2048/data/R/0416/myeloid_tissue_comparison_v1_2_3_20260416",
    )
)
FINAL_RDS = R_OUTPUT_DIR / "myeloid_tissue_comparison_final.rds"
FINAL_H5AD = R_OUTPUT_DIR / "myeloid_tissue_comparison_final.h5ad"
REPORT_MD = R_OUTPUT_DIR / "REPORT.md"
REQUIRED_R_FILES = [FINAL_RDS, FINAL_H5AD, REPORT_MD]
RESPIRATORY_ALVEOLAR_TISSUES = {"lung parenchyma", "respiratory airway"}
ALVEOLAR_LABELS = {"Resident Alveolar macrophages", "Resting Alveolar macrophages"}
SKIP_R_VALIDATION = os.getenv("MYELOID_VALIDATE_SKIP_R", "0").strip().lower() in {"1", "true", "yes"}


def fail(message: str) -> None:
    print(f"[validate] ERROR: {message}", file=sys.stderr, flush=True)
    raise SystemExit(1)


def check_exists_and_nonempty(path: Path) -> None:
    if not path.exists():
        fail(f"Missing expected file: {path}")
    if path.stat().st_size <= 0:
        fail(f"File exists but is empty: {path}")
    print(f"[validate] OK file {path} ({path.stat().st_size} bytes)", flush=True)


def validate_patch_output() -> None:
    check_exists_and_nonempty(PATCH_OUTPUT)
    adata = ad.read_h5ad(PATCH_OUTPUT, backed="r")
    try:
        obs = adata.obs
        if "alveolar_to_interstitial_tissue_rule_applied" not in obs.columns:
            fail("Patch output missing audit column 'alveolar_to_interstitial_tissue_rule_applied'")
        if "cell_type_L3_refined" not in obs.columns:
            fail("Patch output missing 'cell_type_L3_refined'")
        values = obs["alveolar_to_interstitial_tissue_rule_applied"].astype(str)
        unique_values = sorted(set(values.tolist()))
        allowed = {"yes", "no"}
        if not set(unique_values).issubset(allowed):
            fail(f"Unexpected audit values: {unique_values}")
        yes_count = int((values == "yes").sum())
        print(f"[validate] Patch audit column present; yes_count={yes_count}", flush=True)

        if {"cell_type_L3_refined", "tissue"}.issubset(obs.columns):
            l3_values = obs["cell_type_L3_refined"].astype(str)
            tissue_values = obs["tissue"].astype(str)
            invalid_mask = l3_values.isin(ALVEOLAR_LABELS) & (~tissue_values.isin(RESPIRATORY_ALVEOLAR_TISSUES))
            invalid_count = int(invalid_mask.sum())
            if invalid_count != 0:
                fail(
                    "Patch output still contains non-respiratory alveolar macrophage labels "
                    f"(count={invalid_count})"
                )
            print("[validate] No non-respiratory alveolar labels remain in patch output", flush=True)
    finally:
        adata.file.close()


def validate_r_outputs() -> None:
    for path in REQUIRED_R_FILES:
        check_exists_and_nonempty(path)

    report_text = REPORT_MD.read_text(encoding="utf-8", errors="ignore")
    if "Myeloid Tissue Comparison Report" not in report_text:
        fail("REPORT.md does not contain the expected report title")
    if "v1.2.3-MYELOID" not in report_text:
        fail("REPORT.md does not mention v1.2.3-MYELOID")
    print("[validate] REPORT.md title/version checks passed", flush=True)

    adata = ad.read_h5ad(FINAL_H5AD, backed="r")
    try:
        obs_cols = set(adata.obs.columns)
        required_obs_cols = {"cell_type_L2", "cell_type_L3"}
        missing = sorted(required_obs_cols - obs_cols)
        if missing:
            fail(f"Final R h5ad missing obs columns: {missing}")
        print("[validate] Final R h5ad contains cell_type_L2 and cell_type_L3", flush=True)
    finally:
        adata.file.close()


if __name__ == "__main__":
    validate_patch_output()
    if SKIP_R_VALIDATION:
        print("[validate] Skipping R output validation by request.", flush=True)
    else:
        validate_r_outputs()
    print("[validate] All checks passed.", flush=True)
