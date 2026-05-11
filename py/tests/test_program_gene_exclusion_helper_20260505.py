import importlib.util
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd


MODULE_PATH = pathlib.Path('/home/h2048/script/py/program_gene_exclusion_helper_20260505_v1.py')


def load_module():
    spec = importlib.util.spec_from_file_location('program_gene_exclusion_helper_20260505_v1', MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class ProgramGeneExclusionHelperTest(unittest.TestCase):
    def test_build_gene_exclusion_packet_excludes_expected_categories_outside_b_lineage(self):
        module = load_module()
        packet = module.build_gene_exclusion_packet(
            gene_names=['IGHG1', 'MT-CO1', 'RPS12', 'LINC00152', 'MS4A1'],
            lineage_context='myeloid',
        )
        audit = packet['audit_df'].set_index('gene_symbol')
        self.assertTrue(bool(audit.loc['IGHG1', 'should_exclude']))
        self.assertTrue(bool(audit.loc['MT-CO1', 'should_exclude']))
        self.assertTrue(bool(audit.loc['RPS12', 'should_exclude']))
        self.assertTrue(bool(audit.loc['LINC00152', 'should_exclude']))
        self.assertFalse(bool(audit.loc['MS4A1', 'should_exclude']))

    def test_apply_gene_exclusion_to_adata_keeps_ig_for_b_lineage_and_writes_sidecars(self):
        module = load_module()
        obs = pd.DataFrame({'cell_type': ['Naive_B', 'Naive_B', 'Naive_B']}, index=[f'cell_{i}' for i in range(3)])
        var = pd.DataFrame(index=pd.Index(['IGHG1', 'MT-CO1', 'RPS12', 'LINC00152', 'MS4A1'], name='gene'))
        adata = ad.AnnData(X=np.ones((3, 5), dtype=np.float32), obs=obs, var=var)
        adata.layers['counts'] = adata.X.copy()

        with tempfile.TemporaryDirectory() as tmpdir:
            filtered, packet, manifest = module.apply_gene_exclusion_to_adata(
                adata,
                output_dir=pathlib.Path(tmpdir),
                celltype_col='cell_type',
            )
            self.assertEqual(filtered.var_names.tolist(), ['IGHG1', 'MS4A1'])
            self.assertTrue(pathlib.Path(manifest['summary_csv']).exists())
            self.assertTrue(pathlib.Path(manifest['audit_csv']).exists())
            self.assertTrue(pathlib.Path(manifest['summary_plot']['png']).exists())
            self.assertTrue(pathlib.Path(manifest['manifest_json']).exists())


if __name__ == '__main__':
    unittest.main()
