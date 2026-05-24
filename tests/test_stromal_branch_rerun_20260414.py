import importlib.util
import pathlib
import sys
import types
import unittest

import pandas as pd

MODULE_PATH = pathlib.Path("/home/h2048/script/py/stromal_branch_rerun_20260414.py")


def load_module():
    if not MODULE_PATH.exists():
        raise FileNotFoundError(f"Module under test does not exist yet: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location("stromal_branch_rerun_20260414", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class StromalBranchRerunHelpersTest(unittest.TestCase):
    def test_normalize_cluster_id_strips_prefix(self):
        module = load_module()
        self.assertEqual(module.normalize_cluster_id("c33"), "33")
        self.assertEqual(module.normalize_cluster_id("33"), "33")

    def test_extract_scvi_init_kwargs_flattens_nested_groups(self):
        module = load_module()
        init_params = {
            "non_kwargs": {"n_hidden": 128, "n_latent": 75, "use_observed_lib_size": False},
            "kwargs": {
                "model_kwargs": {"dropout_rate": 0.2, "gene_likelihood": "nb"},
                "misc_kwargs": {"dispersion": "gene-batch"},
            },
        }
        merged = module.extract_scvi_init_kwargs(
            init_params,
            module_kwargs={"n_hidden": 256, "n_latent": 30, "latent_distribution": "ln"},
        )
        self.assertEqual(
            merged,
            {
                "use_observed_lib_size": False,
                "n_hidden": 256,
                "n_latent": 30,
                "dropout_rate": 0.2,
                "gene_likelihood": "nb",
                "dispersion": "gene-batch",
                "latent_distribution": "ln",
            },
        )

    def test_select_target_cells_filters_requested_clusters(self):
        module = load_module()
        obs_names = ["cell1", "cell2", "cell3", "cell4"]
        choir_values = ["5", "33", "8", "33"]
        selected = module.select_target_cells(obs_names, choir_values, {"5", "33"})
        self.assertEqual(selected, ["cell1", "cell2", "cell4"])

    def test_configure_runtime_prunes_user_site_paths(self):
        module = load_module()
        original = list(module.sys.path)
        try:
            module.sys.path = [
                "/tmp/fake",
                "/home/h2048/.local/lib/python3.10/site-packages",
                "/another/path",
            ]
            module.configure_runtime()
            self.assertNotIn(
                "/home/h2048/.local/lib/python3.10/site-packages",
                module.sys.path,
            )
            self.assertIn("/tmp/fake", module.sys.path)
            self.assertEqual(module.os.environ.get("PYTHONNOUSERSITE"), "1")
        finally:
            module.sys.path = original

    def test_relabel_cells_updates_only_target_cluster(self):
        module = load_module()
        obs_names = ["cell1", "cell2", "cell3"]
        choir_values = ["6", "1", "6"]
        labels = ["Muscle_smooth_arterial_systemic", "Muscle_pericyte_systemic", "Muscle_smooth_arterial_systemic"]
        relabeled = module.relabel_by_cluster(
            obs_names=obs_names,
            choir_values=choir_values,
            labels=labels,
            relabel_map={"6": "Peripheral_neuron_like"},
        )
        self.assertEqual(
            relabeled,
            ["Peripheral_neuron_like", "Muscle_pericyte_systemic", "Peripheral_neuron_like"],
        )

    def test_load_existing_scvi_forces_cpu_legacy_load(self):
        module = load_module()
        calls = {"setup_calls": []}

        class FakeLoadResult:
            missing_keys = []
            unexpected_keys = []

        class FakeInnerModule:
            def state_dict(self):
                return {"weight": 1}

            def load_state_dict(self, state_dict, strict=False):
                calls["state_dict"] = state_dict
                calls["strict"] = strict
                return FakeLoadResult()

        class LoadedModel:
            def __init__(self):
                self.init_params_ = {
                    "non_kwargs": {"use_observed_lib_size": False},
                    "kwargs": {"model_kwargs": {"dropout_rate": 0.2}},
                }
                self._module_kwargs = {"n_hidden": 128, "n_latent": 30}
                self.module = FakeInnerModule()
                self.train_indices_ = [0, 1]
                self.test_indices_ = [2]
                self.validation_indices_ = [3]
                self.history_ = {"elbo_train": [1.0]}

        class FakeSCVI:
            def __init__(self, adata, **kwargs):
                calls["rehydrated_adata"] = adata
                calls["rehydrated_kwargs"] = kwargs
                self.module = FakeInnerModule()
                self.is_trained_ = False
                self.train_indices_ = None
                self.test_indices_ = None
                self.validation_indices_ = None
                self.history_ = None

            @staticmethod
            def setup_anndata(adata, layer=None, batch_key=None, continuous_covariate_keys=None):
                calls["setup_calls"].append(
                    {
                        "adata": adata,
                        "layer": layer,
                        "batch_key": batch_key,
                        "continuous_covariate_keys": continuous_covariate_keys,
                    }
                )

            @classmethod
            def load(cls, model_dir, adata=None, accelerator=None, device=None):
                calls["load"] = {
                    "model_dir": model_dir,
                    "adata": adata,
                    "accelerator": accelerator,
                    "device": device,
                }
                return LoadedModel()

        fake_scvi_module = types.SimpleNamespace(
            model=types.SimpleNamespace(SCVI=FakeSCVI)
        )
        fake_adata = object()

        rehydrated_model, accelerator = module.load_existing_scvi(
            fake_adata,
            fake_scvi_module,
            torch_module=None,
            continuous_covariates=["pct_counts_mt"],
            model_dir="/tmp/fake_scvi_model",
        )

        self.assertEqual(calls["load"]["accelerator"], "cpu")
        self.assertEqual(calls["load"]["device"], 1)
        self.assertEqual(accelerator, "cpu")
        self.assertEqual(len(calls["setup_calls"]), 2)
        for setup_call in calls["setup_calls"]:
            self.assertEqual(setup_call["layer"], "counts")
            self.assertEqual(setup_call["batch_key"], module.BATCH_KEY)
            self.assertEqual(setup_call["continuous_covariate_keys"], ["pct_counts_mt"])
        self.assertEqual(
            calls["rehydrated_kwargs"],
            {
                "use_observed_lib_size": False,
                "n_hidden": 128,
                "n_latent": 30,
                "dropout_rate": 0.2,
            },
        )
        self.assertTrue(rehydrated_model.is_trained_)
        self.assertEqual(rehydrated_model.train_indices_, [0, 1])
        self.assertEqual(rehydrated_model.test_indices_, [2])
        self.assertEqual(rehydrated_model.validation_indices_, [3])
        self.assertEqual(rehydrated_model.history_, {"elbo_train": [1.0]})

    def test_train_scanvi_gpu_fallback_reuses_in_memory_scvi(self):
        module = load_module()
        original_resolve_accelerator = module.resolve_accelerator
        calls = {"from_scvi": [], "train_accelerators": [], "saved_paths": []}

        class FakeStateRegistry:
            categorical_mapping = ["sample_a"]

        class FakeAdataManager:
            def get_state_registry(self, key):
                self.last_key = key
                return FakeStateRegistry()

        class FakeScviModel:
            def __init__(self):
                self.adata_manager = FakeAdataManager()
                self.to_device_calls = []

            def to_device(self, device):
                self.to_device_calls.append(device)

        class FakeScanviModel:
            def __init__(self, should_fail):
                self.should_fail = should_fail

            def train(self, accelerator=None, **kwargs):
                calls["train_accelerators"].append(accelerator)
                if self.should_fail and accelerator == "gpu":
                    raise RuntimeError("CUDA out of memory")

            def save(self, path, overwrite=False):
                calls["saved_paths"].append((path, overwrite))

        fake_scvi_model = FakeScviModel()
        fake_adata = types.SimpleNamespace(obs=pd.DataFrame({"sample": ["sample_a", "sample_a"]}))

        def fake_from_scvi_model(base_model, unlabeled_category=None, labels_key=None):
            calls["from_scvi"].append(base_model)
            return FakeScanviModel(should_fail=len(calls["from_scvi"]) == 1)

        fake_scvi_module = types.SimpleNamespace(
            model=types.SimpleNamespace(
                SCANVI=types.SimpleNamespace(from_scvi_model=fake_from_scvi_model)
            )
        )
        fake_torch = types.SimpleNamespace(
            cuda=types.SimpleNamespace(
                is_available=lambda: True,
                empty_cache=lambda: calls.setdefault("emptied_cache", True),
            )
        )

        module.resolve_accelerator = lambda _torch: ("gpu", 0)
        try:
            scanvi_model = module.train_scanvi(
                fake_adata,
                fake_scvi_model,
                fake_scvi_module,
                fake_torch,
                n_samples_per_label=8,
                model_dir=pathlib.Path("/tmp/fake_scanvi_model"),
                scvi_model_dir="/tmp/unused_scvi_dir",
            )
        finally:
            module.resolve_accelerator = original_resolve_accelerator

        self.assertEqual(calls["from_scvi"], [fake_scvi_model, fake_scvi_model])
        self.assertEqual(fake_scvi_model.adata_manager.last_key, "batch")
        self.assertEqual(fake_scvi_model.to_device_calls, ["cpu"])
        self.assertEqual(calls["train_accelerators"], ["gpu", "cpu"])
        self.assertEqual(calls["saved_paths"], [("/tmp/fake_scanvi_model", True)])
        self.assertIsNotNone(scanvi_model)


if __name__ == "__main__":
    unittest.main()
