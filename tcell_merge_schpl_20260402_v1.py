#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
T/NK lineage scHPL entrypoint under the unified schpl methodology.

This file uses the **global lineage-wise** mode:
- one named reference ontology,
- one scHPL tree trained on the shared integrated latent,
- one query prediction pass,
- `Rejected` interpreted as a follow-up flag rather than an automatic novel
    cell-type call.
"""

from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from lineage_merge_schpl_core_20260402 import LineageMergeSchPLConfig, run_pipeline


CONFIG = LineageMergeSchPLConfig(
    lineage_name="T/NK",
    lineage_slug="tnk",
    version="v1_0",
    input_merged_h5ad="/home/h2048/data/py/20260318/tcell_only_merged_pipeline_v2_3_filtered/tcell_merged_v2_3_results.h5ad",
    output_dir="/home/h2048/data/py/0402/tcell_merge_schpl_v1_0",
    # NOTE:
    #   `cell_type_fine_ref` is QC-only in the merged artifact and preserves the
    #   old pipeline subcluster names (e.g. *_c0/_c1/_c2), which do not match
    #   the reconciled refined TNK label space used downstream.
    #   `scanvi_labels_reconciled` is the aligned refined label ontology that
    #   best matches the downstream TNK analysis scripts.
    reference_label_key="scanvi_labels_reconciled",
    query_pred_key="scanvi_pred",
    query_final_key="scanvi_labels_reconciled",
    query_conf_key="scanvi_confidence",
    marker_genes=[
        "CD3D", "TRBC2", "CD4", "IL7R", "CCR7", "LTB", "CD8A", "NKG7",
    ],
    latent_key="X_scANVI",
    umap_key="X_umap",
    marker_layer_preferred="log1p",
)


if __name__ == "__main__":
    run_pipeline(CONFIG)
