#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""B-cell 2026-05-08 full scVI/scANVI rerun.

Builds on the validated 2026-04-17 B-cell rerun template and removes the
previously retained user-curated outlier CHOIR c14 in addition to the existing
c22/c13/c25 removals before fresh scVI/scANVI training.
"""

from __future__ import annotations

import re
import os
from pathlib import Path

BASE_SCRIPT = Path("/home/h2048/script/py/bcell_scvi_scanvi_ref_c13c25drop_20260417.py")
if not BASE_SCRIPT.exists():
    raise FileNotFoundError(BASE_SCRIPT)

code = BASE_SCRIPT.read_text(encoding="utf-8")

replacement_sources = '''CHOIR_CLUSTER_SOURCES = (
    {
        "csv": Path(
            "/home/h2048/data/R/0412/bcell_tissue_comparison_v2_6_4_20260412/reports/choir/choir_clusters.csv"
        ),
        "clusters": (22,),
        "source_label": "0412_validated_c22",
    },
    {
        "csv": Path(
            "/home/h2048/data/R/0415/bcell_tissue_comparison_v2_6_6_c22drop_l3_20260415/reports/choir/choir_clusters.csv"
        ),
        "clusters": (13, 25, 14),
        "source_label": "0415_current_c13_c25_c14",
    },
)
'''
code = re.sub(
    r'CHOIR_CLUSTER_SOURCES = \(.*?\)\nOUTPUT_DIR = Path',
    replacement_sources + '\nOUTPUT_DIR = Path',
    code,
    flags=re.S,
)
if '"source_label": "0415_current_c13_c25_c14"' not in code:
    raise RuntimeError("Failed to patch CHOIR_CLUSTER_SOURCES with c14 removal")
code = code.replace(
    'OUTPUT_DIR = Path("/home/h2048/data/py/0415/bcell_scvi_scanvi_ref_c22drop_l3_20260415")',
    'OUTPUT_DIR = Path("/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508")',
)
code = code.replace(
    'OUTPUT_H5AD = OUTPUT_DIR / "bcell_reference_c22_c13_c25drop_scanvi_L3_ref_20260417.h5ad"',
    'OUTPUT_H5AD = OUTPUT_DIR / "bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad"',
)
code = code.replace(
    'OUTPUT_SCVI_DIR = OUTPUT_DIR / "scvi_bcell_ref_c22_c13_c25drop_l3_20260417"',
    'OUTPUT_SCVI_DIR = OUTPUT_DIR / "scvi_bcell_ref_c22_c13_c25_c14drop_l3_20260508"',
)
code = code.replace(
    'OUTPUT_SCANVI_DIR = OUTPUT_DIR / "scanvi_bcell_L3_ref_c22_c13_c25drop_20260417"',
    'OUTPUT_SCANVI_DIR = OUTPUT_DIR / "scanvi_bcell_L3_ref_c22_c13_c25_c14drop_20260508"',
)
code = code.replace(
    'OUTPUT_REMOVED = OUTPUT_DIR / "removed_choir_clusters_c22_c13_c25_cells.csv"',
    'OUTPUT_REMOVED = OUTPUT_DIR / "removed_choir_clusters_c22_c13_c25_c14_cells.csv"',
)
code = code.replace(
    'OUTPUT_TARGETS = OUTPUT_DIR / "target_choir_cluster_union_c22_c13_c25.csv"',
    'OUTPUT_TARGETS = OUTPUT_DIR / "target_choir_cluster_union_c22_c13_c25_c14.csv"',
)
code = code.replace(
    'B-cell reference scVI/scANVI full rerun after removing CHOIR clusters 22, 13, and 25.',
    'B-cell reference scVI/scANVI full rerun after removing CHOIR clusters 22, 13, 25, and 14.',
)

if os.environ.get("BCELL_20260508_DRY_RUN") == "1":
    start = code.index("CHOIR_CLUSTER_SOURCES = (")
    end = code.index("OUTPUT_DIR = Path", start)
    print(code[start:end])
    raise SystemExit(0)

exec_globals = {"__name__": "__main__", "__file__": str(Path(__file__).resolve())}
exec(compile(code, str(BASE_SCRIPT), "exec"), exec_globals)
