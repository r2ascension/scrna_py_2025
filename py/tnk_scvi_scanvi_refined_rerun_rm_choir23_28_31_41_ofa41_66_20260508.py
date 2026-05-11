#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TNK 2026-05-08 full scVI/scANVI rerun.

Removes the user-requested TNK outlier clusters before fresh training:
- CHOIR: c23, c28, c31, c41
- OFA:   c41, c66

OFA outliers in the current TNK reports are produced under the CHOIR cluster
backend, so cell removal uses the same `CHOIR_clusters_0.2` assignment column and
records the provenance split in `reports/requested_removal_provenance.tsv`.
"""

from __future__ import annotations

import re
import os
from pathlib import Path

BASE_SCRIPT = Path("/home/h2048/script/py/tnk_scvi_scanvi_refined_rerun_rm_choir_20260413_v1.py")
if not BASE_SCRIPT.exists():
    raise FileNotFoundError(BASE_SCRIPT)

code = BASE_SCRIPT.read_text(encoding="utf-8")

code = re.sub(
    r'CLUSTER_SOURCE_H5AD = \(\n    "/home/h2048/data/R/0413/tnk_tissue_comparison_v2_6_1_20260413/"\n    "tnk_tissue_comparison_final\.h5ad"\n\)',
    'CLUSTER_SOURCE_H5AD = (\n    "/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper/"\n    "tnk_tissue_comparison_final.h5ad"\n)',
    code,
)
code = re.sub(
    r'INPUT_H5AD = \\\n    "/home/h2048/data/py/0318/tnk_subcluster_retrain/adata_tnk_scanvi_ref_retrain_v1_2\.h5ad"',
    'INPUT_H5AD = \\\n    "/home/h2048/data/py/0414/tnk_scanvi_relabel_rerun_20260414/adata_tnk_scanvi_refined_relabel_rerun_20260414.h5ad"',
    code,
)
code = re.sub(
    r'OUTPUT_DIR = \\\n    "/home/h2048/data/py/0413/tnk_scvi_scanvi_refined_rerun_rm_choir_20260413"',
    'OUTPUT_DIR = \\\n    "/home/h2048/data/py/0508/tnk_scvi_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508"',
    code,
)
code = code.replace(
    'REMOVE_CLUSTERS = ["c10", "c20", "c22", "c29", "c47", "c2", "c17"]',
    'REMOVE_CLUSTERS = ["c23", "c28", "c31", "c41", "c66"]\nREQUESTED_CHOIR_CLUSTERS = ["23", "28", "31", "41"]\nREQUESTED_OFA_CLUSTERS = ["41", "66"]',
)
code = code.replace(
    'FINAL_H5AD_PATH = OUTPUT_PATH / "adata_tnk_scanvi_refined_rerun_rm_choir_20260413.h5ad"',
    'FINAL_H5AD_PATH = OUTPUT_PATH / "adata_tnk_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508.h5ad"',
)
code = code.replace(
    'SCVI_MODEL_DIR = OUTPUT_PATH / "tnk_scvi_refined_rerun_model"',
    'SCVI_MODEL_DIR = OUTPUT_PATH / "tnk_scvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_model"',
)
code = code.replace(
    'SCANVI_MODEL_DIR = OUTPUT_PATH / "tnk_scanvi_refined_rerun_model"',
    'SCANVI_MODEL_DIR = OUTPUT_PATH / "tnk_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_model"',
)
code = code.replace(
    'print("TNK refined-label scVI/scANVI rerun after CHOIR-cluster removal")',
    'print("TNK refined-label scVI/scANVI rerun after CHOIR/OFA outlier-cluster removal (20260508)")',
)
code = code.replace(
    'save_table(cluster_size_df, REPORT_DIR / "removed_cluster_sizes.tsv")',
    'cluster_size_df["requested_as_choir"] = cluster_size_df[CHOIR_CLUSTER_COL].astype(str).isin(REQUESTED_CHOIR_CLUSTERS)\n'
    'cluster_size_df["requested_as_ofa"] = cluster_size_df[CHOIR_CLUSTER_COL].astype(str).isin(REQUESTED_OFA_CLUSTERS)\n'
    'save_table(cluster_size_df, REPORT_DIR / "removed_cluster_sizes.tsv")\n\n'
    'requested_provenance_df = pd.DataFrame({\n'
    '    "source": ["CHOIR", "OFA"],\n'
    '    "requested_clusters": [",".join(REQUESTED_CHOIR_CLUSTERS), ",".join(REQUESTED_OFA_CLUSTERS)],\n'
    '    "cluster_assignment_column": [CHOIR_CLUSTER_COL, CHOIR_CLUSTER_COL],\n'
    '    "note": [\n'
    '        "Direct CHOIR outlier clusters requested by user",\n'
    '        "OFA outlier cluster IDs under the CHOIR cluster backend; removed from the same cluster assignment column"\n'
    '    ]\n'
    '})\n'
    'save_table(requested_provenance_df, REPORT_DIR / "requested_removal_provenance.tsv")',
)

required_tokens = [
    '"/home/h2048/data/R/0414/tnk_tissue_comparison_v2_6_3_20260414_relabel_helper/"',
    '"/home/h2048/data/py/0414/tnk_scanvi_relabel_rerun_20260414/adata_tnk_scanvi_refined_relabel_rerun_20260414.h5ad"',
    'REMOVE_CLUSTERS = ["c23", "c28", "c31", "c41", "c66"]',
    'REQUESTED_OFA_CLUSTERS = ["41", "66"]',
    'requested_removal_provenance.tsv',
]
missing_tokens = [token for token in required_tokens if token not in code]
if missing_tokens:
    raise RuntimeError(f"Failed to patch TNK 20260508 rerun script; missing tokens: {missing_tokens}")

if os.environ.get("TNK_20260508_DRY_RUN") == "1":
    for token in ["CLUSTER_SOURCE_H5AD", "INPUT_H5AD", "OUTPUT_DIR", "REMOVE_CLUSTERS"]:
        start = code.index(token)
        end = code.index("\n", start + len(token))
        print(code[start:end])
    raise SystemExit(0)

exec_globals = {"__name__": "__main__", "__file__": str(Path(__file__).resolve())}
exec(compile(code, str(BASE_SCRIPT), "exec"), exec_globals)
