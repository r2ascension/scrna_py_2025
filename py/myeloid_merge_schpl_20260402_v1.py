#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Myeloid lineage scHPL entrypoint under the unified schpl methodology.

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
    lineage_name="Myeloid",
    lineage_slug="myeloid",
    version="v1_0",
    input_merged_h5ad="/home/h2048/data/py/20260308/myeloid_only_merged_pipeline_v2_1_filtered/myeloid_merged_v2_1_results.h5ad",
    output_dir="/home/h2048/data/py/0402/myeloid_merge_schpl_v1_0",
    # NOTE:
    #   `cell_type_fine_ref` in the merged artifact is a QC-only frozen copy of
    #   the original reference clusters and is not the active training label
    #   space. Using it here yields numeric / cluster-ID-like scHPL outputs.
    #   `scanvi_labels` is the label source written back by the upstream merged
    #   pipeline and matches the intended named myeloid label ontology.
    reference_label_key="scanvi_labels",
    query_pred_key="scanvi_pred",
    query_final_key="scanvi_pred",
    query_conf_key="scanvi_confidence",
    marker_genes=[
        "LYZ", "FCER1G", "CD14", "FCGR3A", "CLEC10A", "XCR1", "FCGR3B", "TPSB2",
    ],
    latent_key="X_scANVI",
    umap_key="X_umap",
    marker_layer_preferred="log1p",
)


if __name__ == "__main__":
    run_pipeline(CONFIG)
