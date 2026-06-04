from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import pandas as pd

SCRIPT_PATH = Path("/home/h2048/script/py/core/integration/curate_lineage_h5ad_20260531.py")
SPEC = importlib.util.spec_from_file_location("curate_lineage_h5ad_20260531", SCRIPT_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class TestMainLineageH5adCleanup20260531(unittest.TestCase):
    def test_build_contamination_mask_case_insensitive(self) -> None:
        obs = pd.DataFrame(
            {
                "cell_type_L3": pd.Categorical(["AT2", "Contamination", "Basal"]),
                "cell_type_manual": ["AT2", "clean", "possible contamination"],
                "numeric": [1, 2, 3],
            },
            index=["a", "b", "c"],
        )

        mask, hits = MODULE.build_contamination_mask(
            obs,
            ["cell_type_L3", "cell_type_manual"],
            "contamination",
        )

        self.assertEqual(mask.tolist(), [False, True, True])
        self.assertEqual(hits["cell_type_L3"], {"Contamination": 1})
        self.assertEqual(hits["cell_type_manual"], {"possible contamination": 1})

    def test_build_tissue_label_mask_matches_canonical_epithelial_rule(self) -> None:
        obs = pd.DataFrame(
            {
                "tissue": ["nose", "sinus", "lung", "nose", "sinus"],
                "cell_type_scanvi_pred": [
                    "AT1_Canonical",
                    "AT2",
                    "AT2",
                    "Goblet",
                    "AT1_MatrixRemodeling",
                ],
            },
            index=list("abcde"),
        )

        mask, counts = MODULE.build_tissue_label_mask(
            obs,
            "tissue",
            ["nose", "sinus"],
            "cell_type_scanvi_pred",
            r"^(?:AT1(?:_| |$)|AT2(?:_| |$))",
        )

        self.assertEqual(mask.tolist(), [True, True, False, False, True])
        self.assertEqual(
            counts,
            {
                "AT1_Canonical": {"nose": 1},
                "AT1_MatrixRemodeling": {"sinus": 1},
                "AT2": {"sinus": 1},
            },
        )

    def test_apply_exact_relabels_on_categorical_column(self) -> None:
        obs = pd.DataFrame(
            {
                "cell_type_L3": pd.Categorical(["PNEC", "AT2", "PNEC"]),
                "cell_type_scanvi_pred": pd.Categorical(["PNEC", "AT1_Canonical", "Basal"]),
            }
        )

        counts = MODULE.apply_exact_relabels(
            obs,
            ["cell_type_L3", "cell_type_scanvi_pred"],
            {"PNEC": "L3"},
        )

        self.assertEqual(
            counts,
            {
                "cell_type_L3": {"PNEC": 2},
                "cell_type_scanvi_pred": {"PNEC": 1},
            },
        )
        self.assertEqual(list(obs["cell_type_L3"].astype(str)), ["L3", "AT2", "L3"])
        self.assertEqual(list(obs["cell_type_scanvi_pred"].astype(str)), ["L3", "AT1_Canonical", "Basal"])


if __name__ == "__main__":
    unittest.main()
