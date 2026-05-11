import importlib.util
import pathlib
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


MODULE_PATH = pathlib.Path('/home/h2048/script/py/pycogaps_helper_20260505_v1.py')


def load_module():
    spec = importlib.util.spec_from_file_location('pycogaps_helper_20260505_v1', MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class PyCoGAPSHelperTest(unittest.TestCase):
    def test_prepare_pycogaps_input_filters_before_dense_conversion(self):
        module = load_module()
        obs = pd.DataFrame({'cell_type': ['Mono', 'Mono', 'Mono']}, index=[f'cell_{i}' for i in range(3)])
        var = pd.DataFrame(index=pd.Index(['IGHG1', 'MT-CO1', 'RPS12', 'LINC00152', 'MS4A1', 'CD74'], name='gene'))
        x = sp.csr_matrix(np.arange(18, dtype=np.float32).reshape(3, 6))
        adata = ad.AnnData(X=x.copy(), obs=obs, var=var)
        adata.layers['counts'] = x.copy()

        with tempfile.TemporaryDirectory() as tmpdir:
            prepared = module.prepare_pycogaps_input(
                adata=adata,
                output_dir=pathlib.Path(tmpdir),
                layer='counts',
                gene_exclusion_config={'lineage_context': 'myeloid'},
                celltype_col='cell_type',
            )

            filtered = prepared['filtered_adata']
            cogaps_adata = prepared['cogaps_adata']
            self.assertEqual(filtered.var_names.tolist(), ['MS4A1', 'CD74'])
            self.assertEqual(cogaps_adata.obs_names.tolist(), ['MS4A1', 'CD74'])
            self.assertEqual(cogaps_adata.var_names.tolist(), ['cell_0', 'cell_1', 'cell_2'])
            self.assertEqual(cogaps_adata.X.shape, (2, 3))
            self.assertTrue(isinstance(cogaps_adata.X, np.ndarray))
            self.assertTrue(pathlib.Path(prepared['gene_exclusion_manifest']['summary_csv']).exists())


if __name__ == '__main__':
    unittest.main()
