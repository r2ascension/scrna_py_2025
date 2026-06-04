#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import pandas as pd

MODULE_PATH = Path('/home/h2048/script/py/epithelial/integration/refine_epithelial_l3_labels_20260531.py')
SPEC = importlib.util.spec_from_file_location('refine_epithelial_l3_labels_20260531', MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RefineEpithelialL3LabelsTests(unittest.TestCase):
    def test_build_updated_columns_moves_neuroendocrine_into_l3_and_renames_ionocyte(self) -> None:
        obs = pd.DataFrame(
            {
                'celltypist_pred_original': pd.Categorical(['Neuroendocrine', 'Other', 'Other']),
                'cell_type_scanvi_pred': pd.Categorical(['Ionocyte_Brush', 'Ionocyte_Brush', 'Goblet']),
                'cell_type_L3_curated': pd.Categorical(['Goblet', 'Ionocyte_Brush', 'Basal_Inflammatory']),
            },
            index=['c1', 'c2', 'c3'],
        )
        updated_columns, neuro_mask, changed_rows, summary_columns = MODULE.build_updated_columns(
            obs=obs,
            target_columns=['cell_type_scanvi_pred', 'cell_type_L3_curated'],
            neuro_source_column='celltypist_pred_original',
            neuro_source_label='Neuroendocrine',
            neuro_target_label='Neuroendocrine',
            rename_old='Ionocyte_Brush',
            rename_new='Ionocyte',
        )

        self.assertEqual(int(neuro_mask.sum()), 1)
        self.assertEqual(updated_columns['cell_type_scanvi_pred'].tolist(), ['Neuroendocrine', 'Ionocyte', 'Goblet'])
        self.assertEqual(updated_columns['cell_type_L3_curated'].tolist(), ['Neuroendocrine', 'Ionocyte', 'Basal_Inflammatory'])
        self.assertEqual(int(summary_columns['cell_type_scanvi_pred']['n_neuro_target']), 1)
        self.assertEqual(int(summary_columns['cell_type_scanvi_pred']['n_rename_target']), 1)
        self.assertEqual(changed_rows.shape[0], 2)

    def test_build_categories_replaces_old_name_and_appends_neuroendocrine(self) -> None:
        original = pd.Series(pd.Categorical(['Ionocyte_Brush', 'Goblet', 'Ionocyte_Brush']))
        updated = pd.Series(pd.Categorical(['Neuroendocrine', 'Goblet', 'Ionocyte']))
        categories = MODULE.build_categories(
            original=original,
            updated=updated,
            rename_old='Ionocyte_Brush',
            rename_new='Ionocyte',
            neuro_target_label='Neuroendocrine',
        )
        self.assertEqual(categories, ['Goblet', 'Ionocyte', 'Neuroendocrine'])


if __name__ == '__main__':
    unittest.main()
