from __future__ import annotations

import gc
import json
import os
import sys
import time
import traceback
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scipy import sparse

warnings.filterwarnings("ignore", category=FutureWarning)
ad.settings.allow_write_nullable_strings = True

PIPELINE_NAME = "allcells_combined_scanvi_coredata0524"
VERSION = "20260524_v1"
DATE_TAG = datetime.now().strftime("%Y%m%d")
SEED = 42

OUTPUT_DIR = Path("/home/h2048/data/coredata0524")
MODEL_DIR = OUTPUT_DIR / "models"
OUTPUT_STEM = f"coredata0524_allcells_scanvi_l1_{DATE_TAG}"
OUTPUT_H5AD = OUTPUT_DIR / f"{OUTPUT_STEM}.h5ad"
LOG_FILE = OUTPUT_DIR / f"{OUTPUT_STEM}.log"
STRUCTURE_JSON = OUTPUT_DIR / f"anndata_structure_{OUTPUT_STEM}.json"
INPUT_TSV = OUTPUT_DIR / f"selected_inputs_{OUTPUT_STEM}.tsv"
ALL_GENES_CSV = OUTPUT_DIR / "all_genes_union.csv"
MANIFEST_JSON = OUTPUT_DIR / "reference_manifest.json"
RUN_LOG = OUTPUT_DIR / "pipeline_run_log.txt"

BATCH_KEY = "sample"
LABELS_KEY = "scanvi_label"
UNLABELED_CATEGORY = "Unknown"
L1_ORDER = ["Bcell", "Epithelial", "Myeloid", "Stromal", "TNK", UNLABELED_CATEGORY]
CONTINUOUS_COVARIATES = ["pct_counts_mt", "stress_score", "S_score", "G2M_score"]
CATEGORICAL_COVARIATES: list[str] = []

MIN_CELLS_PER_BATCH = 20
HVG_TOP_GENES = 6000
MIN_GENES_PER_CELL = 200
MIN_CELLS_PER_GENE = 3

SCVI_N_LATENT = 150
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 256
SCVI_DROPOUT_RATE = 0.2
SCVI_USE_LAYER_NORM = "both"
SCVI_USE_BATCH_NORM = "none"
SCVI_ENCODE_COVARIATES = True
SCVI_MAX_EPOCHS = 40
SCVI_BATCH_SIZE = 2048
SCVI_LEARNING_RATE = 1e-3
SCVI_EARLY_STOPPING = True
SCVI_EARLY_STOPPING_PATIENCE = 10

SCANVI_MAX_EPOCHS = 30
SCANVI_BATCH_SIZE = 2048
SCANVI_LEARNING_RATE = 5e-4
SCANVI_EARLY_STOPPING_PATIENCE = 8

LINEAGE_INPUTS: dict[str, dict[str, object]] = {
    "bcell": {
        "h5ad": Path("/home/h2048/data/py/0508/bcell_scvi_scanvi_ref_c22_c13_c25_c14drop_20260508/bcell_reference_c22_c13_c25_c14drop_scanvi_L3_ref_20260508.h5ad"),
        "l1_label": "Bcell",
        "batch_key_local": "sample",
    },
    "stromal_fibroblast": {
        "h5ad": Path("/home/h2048/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/fibroblast/adata_fibroblast_reference_v1_5_branchwise.h5ad"),
        "l1_label": "Stromal",
        "batch_key_local": "sample",
    },
    "stromal_smc": {
        "h5ad": Path("/home/h2048/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/smc/adata_smc_reference_v1_5_branchwise.h5ad"),
        "l1_label": "Stromal",
        "batch_key_local": "sample",
    },
    "stromal_endothelial": {
        "h5ad": Path("/home/h2048/data/py/0508/stromal_branch_rerun_rm_endothelial6_52_20260508/endothelial/adata_endothelial_reference_v1_5_branchwise.h5ad"),
        "l1_label": "Stromal",
        "batch_key_local": "sample",
    },
    "tnk": {
        "h5ad": Path("/home/h2048/data/py/0508/tnk_scvi_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508/adata_tnk_scanvi_refined_rerun_rm_choir23_28_31_41_ofa41_66_20260508.h5ad"),
        "l1_label": "TNK",
        "batch_key_local": "sample",
    },
    "myeloid": {
        "h5ad": Path("/home/h2048/data/py/0416/adata_myeloid_L3refined_tissueaware_patched_v1.h5ad"),
        "l1_label": "Myeloid",
        "batch_key_local": "sample",
    },
    "epithelial": {
        "h5ad": Path("/home/h2048/data/py/0508/epithelial_scanvi_rm_leiden14_17_20260508/epithelial_scanvi_rm_leiden14_17_SELF_for_R.h5ad"),
        "l1_label": "Epithelial",
        "batch_key_local": "sample",
    },
}

_BCELL_L3_MERGE = {"IGHEplus_Atypical_Memory_B": "Atypical_Memory_B"}
_BCELL_L2_MAP = {
    "GC_B_Dark_Zone_Centroblast_Cycling": "GC_B",
    "GC_B_Light_Zone_Centrocyte": "GC_B",
    "GC_B_Transitional": "GC_B",
    "Plasma_IgA": "Plasma",
    "Plasma_IgG": "Plasma",
    "Atypical_Memory_B": "Atypical_Memory_B",
    "Memory_B": "Memory_B",
    "Naive_B": "Naive_B",
}

_SCVI_INTERNAL_OBS_COLS = {
    "_scvi_batch",
    "_scvi_labels",
    "_scvi_extra_categorical_covs",
    "_scvi_extra_continuous_covs",
}
_HVG_AUX_VAR_COLS = {
    "highly_variable",
    "highly_variable_rank",
    "means",
    "variances",
    "variances_norm",
    "dispersions",
    "dispersions_norm",
    "highly_variable_nbatches",
}

STRESS_SIGNATURE_GENES = [
    "ALDH18A1", "ARFGAP1", "ASNS", "ATF3", "ATF4", "ATF6", "ATP6V0D1", "BAG3", "BANF1",
    "CALR", "CCL2", "CEBPB", "CEBPG", "CHAC1", "CKS1B", "CNOT2", "CNOT4", "CNOT6",
    "CXXC1", "DCP1A", "DCP2", "DCTN1", "DDIT4", "DDX10", "DKC1", "DNAJA4", "DNAJB9",
    "DNAJC3", "EDC4", "EDEM1", "EEF2", "EIF2AK3", "EIF2S1", "EIF4A1", "EIF4A2", "EIF4A3",
    "EIF4E", "EIF4EBP1", "EIF4G1", "ERN1", "ERO1A", "EXOC2", "EXOSC1", "EXOSC10",
    "EXOSC2", "EXOSC4", "EXOSC5", "EXOSC9", "FKBP14", "FUS", "GEMIN4", "GOSR2", "H2AX",
    "HERPUD1", "HSP90B1", "HSPA5", "HSPA9", "HYOU1", "IARS1", "IFIT1", "IGFBP1", "IMP3",
    "KDELR3", "KHSRP", "KIF5B", "LSM1", "LSM4", "MTHFD2", "NFYA", "NFYB", "NHP2", "NOLC1",
    "NOP14", "NOP56", "NPM1", "NABP1", "PAIP1", "PARN", "PDIA5", "PDIA6", "POP4", "PREB",
    "PSAT1", "RPS14", "RRP9", "SDAD1", "SEC11A", "SEC31A", "SERP1", "SHC1", "MTREX",
    "SLC1A4", "SLC30A5", "SLC7A5", "SPCS1", "SPCS3", "SRPRA", "SRPRB", "SSR1", "STC2",
    "TARS1", "TATDN2", "TSPYL2", "SKIC3", "TUBB2A", "VEGFA", "WFS1", "WIPI1", "XBP1",
    "XPOT", "YIF1A", "YWHAZ", "ZBTB17",
]
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B",
    "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
    "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data: str) -> None:
        for stream in self.streams:
            stream.write(data)
            stream.flush()

    def flush(self) -> None:
        for stream in self.streams:
            stream.flush()


LOG_HANDLE = None
PIPELINE_START = time.time()


def setup_logging() -> None:
    global LOG_HANDLE
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    LOG_HANDLE = open(LOG_FILE, "w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, LOG_HANDLE)
    sys.stderr = Tee(sys.__stderr__, LOG_HANDLE)


def close_logging() -> None:
    global LOG_HANDLE
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    if LOG_HANDLE is not None:
        LOG_HANDLE.close()
        LOG_HANDLE = None


def clean_string_series(values, fill_value: str) -> pd.Series:
    if isinstance(values, pd.Series):
        idx = values.index
        s = values.astype(object).copy()
    else:
        idx = None
        s = pd.Series(values, dtype=object)
    if idx is not None:
        s.index = idx
    s = s.where(~pd.isna(s), fill_value)

    def _clean_one(x: object) -> str:
        if pd.isna(x):
            return fill_value
        txt = str(x).strip()
        if txt in {"", "nan", "NaN", "None", "<NA>"}:
            return fill_value
        return txt

    return s.map(_clean_one)


def to_clean_category(values, fill_value: str) -> pd.Categorical:
    return pd.Categorical(clean_string_series(values, fill_value))


def sanitize_category_columns(df: pd.DataFrame) -> None:
    for col in df.select_dtypes(include=["category"]).columns:
        cats = df[col].cat.categories
        if hasattr(cats.dtype, "name") and cats.dtype.name in ("string", "StringDtype"):
            df[col] = df[col].cat.rename_categories(cats.astype(object))


def disable_broken_cuda(torch_module) -> None:
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    try:
        torch_module.cuda.is_available = lambda: False  # type: ignore[assignment]
        torch_module.cuda.device_count = lambda: 0  # type: ignore[assignment]
        torch_module.cuda.get_rng_state_all = lambda: []  # type: ignore[assignment]
        torch_module.cuda.manual_seed_all = lambda *_args, **_kwargs: None  # type: ignore[assignment]
        torch_module.cuda.empty_cache = lambda: None  # type: ignore[assignment]
    except Exception:
        pass


def resolve_accelerator(torch_module) -> tuple[str, int]:
    if not torch_module.cuda.is_available():
        return "cpu", 1
    try:
        _ = torch_module.cuda.get_device_properties(0)
        return "gpu", 1
    except Exception as exc:
        print(f"[WARN] CUDA reported available but failed to initialize, fallback to CPU: {exc}")
        disable_broken_cuda(torch_module)
        return "cpu", 1


def normalize_lineage_structure(adata_lin: ad.AnnData, key: str, local_bk: str, global_bk: str) -> None:
    print(f"  [normalize] {key}")
    if "counts" not in adata_lin.layers:
        if "raw_counts" in adata_lin.layers:
            adata_lin.layers["counts"] = adata_lin.layers["raw_counts"]
        elif adata_lin.raw is not None:
            raw_idx = pd.Index(adata_lin.raw.var_names)
            cur_idx = pd.Index(adata_lin.var_names)
            shared = raw_idx.intersection(cur_idx)
            if len(shared) == len(cur_idx):
                pos = raw_idx.get_indexer(cur_idx)
                adata_lin.layers["counts"] = sparse.csr_matrix(adata_lin.raw.X[:, pos], dtype=np.float32)
            else:
                raise ValueError(f"{key}: counts layer missing and raw does not fully cover current var_names")
        else:
            raise ValueError(f"{key}: no counts source found")

    counts = adata_lin.layers["counts"]
    if sparse.issparse(counts):
        if not sparse.isspmatrix_csr(counts):
            counts = counts.tocsr()
        if counts.dtype != np.float32:
            counts = counts.astype(np.float32)
    else:
        counts = sparse.csr_matrix(np.asarray(counts, dtype=np.float32))
    if counts.nnz > 0:
        if not np.isfinite(counts.data).all():
            raise ValueError(f"{key}: non-finite values in counts")
        if np.any(counts.data < 0):
            raise ValueError(f"{key}: negative values in counts")
    adata_lin.layers["counts"] = counts

    if global_bk not in adata_lin.obs.columns:
        if local_bk in adata_lin.obs.columns:
            adata_lin.obs[global_bk] = clean_string_series(adata_lin.obs[local_bk], "UnknownSample")
        else:
            candidates = [
                c for c in adata_lin.obs.columns
                if any(k in c.lower() for k in ("sample", "batch", "dataset"))
            ]
            raise ValueError(f"{key}: batch column not found. Candidates={candidates}")
    else:
        adata_lin.obs[global_bk] = clean_string_series(adata_lin.obs[global_bk], "UnknownSample")

    drop_obs = [c for c in adata_lin.obs.columns if c in _SCVI_INTERNAL_OBS_COLS]
    if drop_obs:
        adata_lin.obs.drop(columns=drop_obs, inplace=True)
    drop_var = [c for c in adata_lin.var.columns if c in _HVG_AUX_VAR_COLS]
    if drop_var:
        adata_lin.var.drop(columns=drop_var, inplace=True)
    if len(adata_lin.obsm):
        adata_lin.obsm.clear()
    sanitize_category_columns(adata_lin.obs)
    sanitize_category_columns(adata_lin.var)
    print(f"    counts={adata_lin.layers['counts'].shape} batch_levels={adata_lin.obs[global_bk].nunique()}")


def resolve_fullgene_counts_matrix(adata_lin: ad.AnnData, key: str) -> tuple[sparse.csr_matrix, str, pd.DataFrame]:
    candidates: list[dict[str, object]] = []
    if adata_lin.raw is not None:
        candidates.append(
            {
                "source": "adata.raw.X",
                "X": adata_lin.raw.X,
                "var": adata_lin.raw.var.copy(),
                "n_vars": int(adata_lin.raw.n_vars),
                "rank": 0,
            }
        )
    if "raw_counts" in adata_lin.layers:
        candidates.append(
            {
                "source": "layers['raw_counts']",
                "X": adata_lin.layers["raw_counts"],
                "var": adata_lin.var.copy(),
                "n_vars": int(adata_lin.n_vars),
                "rank": 1,
            }
        )
    if "counts" in adata_lin.layers:
        candidates.append(
            {
                "source": "layers['counts']",
                "X": adata_lin.layers["counts"],
                "var": adata_lin.var.copy(),
                "n_vars": int(adata_lin.n_vars),
                "rank": 2,
            }
        )
    if not candidates:
        raise ValueError(f"{key}: no count candidates available")
    candidates.sort(key=lambda c: (-int(c["n_vars"]), int(c["rank"])))

    skipped: list[str] = []
    deferred_nonint: dict[str, object] | None = None
    for cand in candidates:
        X = cand["X"]
        src = str(cand["source"])
        data = X.data if sparse.issparse(X) else np.asarray(X).ravel()
        if len(data) > 0 and not np.isfinite(data).all():
            skipped.append(f"{src}: non-finite")
            continue
        if len(data) > 0 and np.any(data < 0):
            skipped.append(f"{src}: negative")
            continue
        if len(data) > 0:
            sample = data[: min(50000, len(data))]
            frac = np.abs(sample - np.round(sample))
            if np.nanmax(frac) > 1e-3:
                if deferred_nonint is None:
                    deferred_nonint = cand
                skipped.append(f"{src}: non-integer")
                continue
        if sparse.issparse(X):
            if not sparse.isspmatrix_csr(X):
                X = X.tocsr()
            if X.dtype != np.float32:
                X = X.astype(np.float32)
        else:
            X = sparse.csr_matrix(np.asarray(X, dtype=np.float32))
        print(f"  [counts] {key}: {src} | genes={cand['n_vars']:,}")
        return X, src, cand["var"]  # type: ignore[return-value]

    if deferred_nonint is not None:
        X = deferred_nonint["X"]
        if sparse.issparse(X):
            if not sparse.isspmatrix_csr(X):
                X = X.tocsr()
            if X.dtype != np.float32:
                X = X.astype(np.float32)
        else:
            X = sparse.csr_matrix(np.asarray(X, dtype=np.float32))
        src = str(deferred_nonint["source"]) + " [nonint-fallback]"
        print(f"  [WARN] {key}: fallback to {src}")
        return X, src, deferred_nonint["var"]  # type: ignore[return-value]

    raise ValueError(f"{key}: all count candidates failed -> {skipped}")


def first_existing_column(df: pd.DataFrame, columns: list[str]) -> str | None:
    for col in columns:
        if col in df.columns:
            return col
    return None


def build_l2_l3_labels(adata_lin: ad.AnnData, key: str) -> tuple[pd.Series, pd.Series, dict[str, str]]:
    obs = adata_lin.obs
    if key == "bcell":
        l3_col = first_existing_column(
            obs,
            [
                "cell_type_scanvi_pred",
                "L3_scanvi_pred",
                "cell_type_expert",
                "cell_type_expert_input",
                "cell_type_scanvi_pred_input",
                "cell_type_L3",
            ],
        )
        if l3_col is None:
            raise KeyError("bcell: cannot find an L3 label column")
        l3 = clean_string_series(obs[l3_col], UNLABELED_CATEGORY).replace(_BCELL_L3_MERGE)
        l2_col = first_existing_column(obs, ["Cell_Type_L2", "L2_scanvi_pred", "Cell_Type_L2_input"])
        if l2_col is not None:
            l2 = clean_string_series(obs[l2_col], "Unknown")
            return l2, l3, {"l2_source": l2_col, "l3_source": l3_col}
        unmapped = sorted(set(l3.unique()) - set(_BCELL_L2_MAP.keys()))
        if unmapped:
            raise ValueError(f"bcell: unmapped L3 labels for L2 collapse: {unmapped}")
        l2 = l3.map(_BCELL_L2_MAP)
        return l2, l3, {"l2_source": "_BCELL_L2_MAP", "l3_source": l3_col}

    l3_col = first_existing_column(obs, ["cell_type_L3", "scanvi_fine_pred", "cell_type_scanvi_pred", "scanvi_label_refined"])
    l2_col = first_existing_column(obs, ["cell_type_L2", "scanvi_major_pred"])
    if l3_col is None:
        raise KeyError(f"{key}: cannot find an L3 label column")
    l3 = clean_string_series(obs[l3_col], UNLABELED_CATEGORY)
    if l2_col is None:
        print(f"  [WARN] {key}: L2 column missing; fallback to L3 labels")
        l2 = l3.copy()
        l2_source = l3_col + " [fallback-from-L3]"
    else:
        l2 = clean_string_series(obs[l2_col], "Unknown")
        l2_source = l2_col
    return l2, l3, {"l2_source": l2_source, "l3_source": l3_col}


def compute_covariates_from_counts(counts: sparse.csr_matrix, var_names: pd.Index, obs_index: pd.Index) -> pd.DataFrame:
    total = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    upper_names = var_names.astype(str).str.upper()
    mt_mask = upper_names.str.startswith("MT-")
    if int(mt_mask.sum()) > 0:
        mt_sum = np.asarray(counts[:, mt_mask].sum(axis=1)).ravel().astype(np.float64)
        pct_counts_mt = (mt_sum / np.maximum(total, 1.0) * 100.0).astype(np.float32)
    else:
        pct_counts_mt = np.zeros(counts.shape[0], dtype=np.float32)

    tmp = ad.AnnData(X=counts.copy(), var=pd.DataFrame(index=var_names.copy()), obs=pd.DataFrame(index=obs_index.copy()))
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)

    stress_genes = [g for g in STRESS_SIGNATURE_GENES if g in tmp.var_names]
    if len(stress_genes) >= 10:
        stress_idx = [tmp.var_names.get_loc(g) for g in stress_genes]
        stress_expr = tmp.X[:, stress_idx]
        if sparse.issparse(stress_expr):
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
        else:
            stress_score = np.asarray(stress_expr.mean(axis=1)).ravel().astype(np.float32)
    else:
        stress_score = np.zeros(tmp.n_obs, dtype=np.float32)

    s_genes = [g for g in S_GENES if g in tmp.var_names]
    g2m_genes = [g for g in G2M_GENES if g in tmp.var_names]
    if len(s_genes) >= 10 and len(g2m_genes) >= 10:
        sc.tl.score_genes_cell_cycle(tmp, s_genes=s_genes, g2m_genes=g2m_genes)
        s_score = tmp.obs["S_score"].values.astype(np.float32)
        g2m_score = tmp.obs["G2M_score"].values.astype(np.float32)
        phase = tmp.obs["phase"].astype(str).values
    else:
        s_score = np.zeros(tmp.n_obs, dtype=np.float32)
        g2m_score = np.zeros(tmp.n_obs, dtype=np.float32)
        phase = np.repeat("G1", tmp.n_obs)

    del tmp
    gc.collect()
    return pd.DataFrame(
        {
            "pct_counts_mt": pct_counts_mt,
            "stress_score": stress_score,
            "S_score": s_score,
            "G2M_score": g2m_score,
            "phase": phase,
        },
        index=obs_index,
    )


def ensure_covariates(obs: pd.DataFrame, counts: sparse.csr_matrix, var_names: pd.Index, lineage_key: str) -> pd.DataFrame:
    obs = obs.copy()
    numeric_invalid: list[str] = []
    for col in CONTINUOUS_COVARIATES:
        if col in obs.columns:
            numeric = pd.to_numeric(obs[col], errors="coerce")
            if numeric.notna().all():
                obs[col] = numeric.astype(np.float32)
            else:
                numeric_invalid.append(col)
    phase_ok = "phase" in obs.columns
    if phase_ok:
        obs["phase"] = clean_string_series(obs["phase"], "G1")

    missing = [c for c in CONTINUOUS_COVARIATES if c not in obs.columns] + numeric_invalid
    if (not phase_ok) or missing:
        print(f"  [covariates] {lineage_key}: recomputing missing/invalid covariates -> {sorted(set(missing + ([] if phase_ok else ['phase'])))}")
        cov_df = compute_covariates_from_counts(counts, var_names, obs.index)
        for col in CONTINUOUS_COVARIATES:
            if col not in obs.columns or col in numeric_invalid:
                obs[col] = cov_df[col].values
        if not phase_ok:
            obs["phase"] = cov_df["phase"].values

    for col in CONTINUOUS_COVARIATES:
        obs[col] = pd.to_numeric(obs[col], errors="coerce").fillna(0.0).astype(np.float32)
    obs["phase"] = clean_string_series(obs["phase"], "G1")
    return obs


def sanitize_object_cols(df: pd.DataFrame, df_name: str) -> list[str]:
    notes: list[str] = []
    for col in df.columns:
        s = df[col]
        if s.dtype != object:
            continue
        non_na = s.dropna()
        if len(non_na) == 0:
            df[col] = s.fillna("").astype(str)
            notes.append(f"{df_name}.{col}: empty->str")
            continue
        if non_na.map(lambda x: isinstance(x, (bool, np.bool_))).all():
            df[col] = s.fillna(False).astype(np.int8)
            notes.append(f"{df_name}.{col}: bool->int8")
            continue
        if non_na.map(lambda x: not isinstance(x, str)).any():
            df[col] = s.map(lambda x: "" if pd.isna(x) else str(x))
            notes.append(f"{df_name}.{col}: mixed->str")
    return notes


def li_matrix(x) -> dict[str, object]:
    return {
        "type": type(x).__name__ if sparse.issparse(x) else "ndarray",
        "shape": list(x.shape),
        "dtype": str(x.dtype),
        "sparse": bool(sparse.issparse(x)),
    }


def col_info(s: pd.Series) -> dict[str, object]:
    out: dict[str, object] = {"dtype": str(s.dtype), "non_null": int(s.notna().sum())}
    try:
        out["unique"] = int(s.nunique())
    except Exception:
        pass
    return out


def save_reference_artifacts(
    adata: ad.AnnData,
    input_rows: list[dict[str, object]],
    accelerator: str,
    devices: int,
) -> None:
    pd.DataFrame(input_rows).to_csv(INPUT_TSV, sep="\t", index=False)
    pd.Series(adata.raw.var_names if adata.raw is not None else adata.var_names).to_csv(
        ALL_GENES_CSV, index=False, header=False
    )

    agreement = float(adata.uns.get("agreement_rate_overall", -1.0))
    manifest = {
        "pipeline": PIPELINE_NAME,
        "version": VERSION,
        "date": datetime.now().strftime("%Y-%m-%d"),
        "output_h5ad": str(OUTPUT_H5AD),
        "scvi_model_dir": str(MODEL_DIR / "scvi_model"),
        "scanvi_model_dir": str(MODEL_DIR / "scanvi_model"),
        "selected_inputs_tsv": str(INPUT_TSV),
        "all_genes_union_csv": str(ALL_GENES_CSV),
        "batch_key": BATCH_KEY,
        "label_key": LABELS_KEY,
        "final_label_key": "cell_type_final_l1",
        "prediction_key": "cell_type_scanvi_pred_pan",
        "confidence_key": "scanvi_confidence",
        "continuous_covariates": CONTINUOUS_COVARIATES,
        "categorical_covariates": CATEGORICAL_COVARIATES,
        "unlabeled_category": UNLABELED_CATEGORY,
        "n_obs": int(adata.n_obs),
        "n_vars_hvg": int(adata.n_vars),
        "n_vars_raw": int(adata.raw.n_vars) if adata.raw is not None else None,
        "scvi_params": {
            "n_latent": SCVI_N_LATENT,
            "n_layers": SCVI_N_LAYERS,
            "n_hidden": SCVI_N_HIDDEN,
            "dropout_rate": SCVI_DROPOUT_RATE,
            "use_layer_norm": SCVI_USE_LAYER_NORM,
            "use_batch_norm": SCVI_USE_BATCH_NORM,
            "encode_covariates": SCVI_ENCODE_COVARIATES,
            "max_epochs": SCVI_MAX_EPOCHS,
            "batch_size": SCVI_BATCH_SIZE,
            "learning_rate": SCVI_LEARNING_RATE,
            "early_stopping_patience": SCVI_EARLY_STOPPING_PATIENCE,
        },
        "scanvi_params": {
            "max_epochs": SCANVI_MAX_EPOCHS,
            "batch_size": SCANVI_BATCH_SIZE,
            "learning_rate": SCANVI_LEARNING_RATE,
            "early_stopping_patience": SCANVI_EARLY_STOPPING_PATIENCE,
        },
        "accelerator": accelerator,
        "devices": devices,
        "agreement_rate": agreement,
        "agreement_rate_is_resubstitution": True,
        "inputs": input_rows,
    }
    with open(MANIFEST_JSON, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)

    structure = {
        "pipeline": PIPELINE_NAME,
        "version": VERSION,
        "timestamp": datetime.now().isoformat(),
        "output_h5ad": str(OUTPUT_H5AD),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "X": li_matrix(adata.X),
        "raw_n_vars": int(adata.raw.n_vars) if adata.raw is not None else None,
        "layers": {k: li_matrix(v) for k, v in adata.layers.items()},
        "obsm": {k: {"shape": list(v.shape), "dtype": str(v.dtype)} for k, v in adata.obsm.items()},
        "batch_key": BATCH_KEY,
        "label_key": LABELS_KEY,
        "lineage_distribution": adata.obs["lineage_source"].value_counts().to_dict(),
        "lineage_branch_distribution": adata.obs["lineage_branch"].value_counts().to_dict(),
        "L1_final_distribution": adata.obs["cell_type_final_l1"].value_counts().to_dict(),
        "agreement_rate": agreement,
        "key_obs_columns": {
            "lineage_source": col_info(adata.obs["lineage_source"]),
            "lineage_branch": col_info(adata.obs["lineage_branch"]),
            "cell_type_input_l1": col_info(adata.obs["cell_type_input_l1"]),
            "cell_type_input_l2": col_info(adata.obs["cell_type_input_l2"]),
            "cell_type_input_l3": col_info(adata.obs["cell_type_input_l3"]),
            LABELS_KEY: col_info(adata.obs[LABELS_KEY]),
            "cell_type_final_l1": col_info(adata.obs["cell_type_final_l1"]),
            "cell_type_scanvi_pred_pan": col_info(adata.obs["cell_type_scanvi_pred_pan"]),
            "scanvi_confidence": col_info(adata.obs["scanvi_confidence"]),
            BATCH_KEY: col_info(adata.obs[BATCH_KEY]),
            "pct_counts_mt": col_info(adata.obs["pct_counts_mt"]),
            "stress_score": col_info(adata.obs["stress_score"]),
            "S_score": col_info(adata.obs["S_score"]),
            "G2M_score": col_info(adata.obs["G2M_score"]),
        },
        "log_file": str(LOG_FILE),
    }
    with open(STRUCTURE_JSON, "w", encoding="utf-8") as handle:
        json.dump(structure, handle, indent=2, ensure_ascii=False)

    elapsed_min = (time.time() - PIPELINE_START) / 60.0
    lines = [
        f"pipeline={PIPELINE_NAME}",
        f"version={VERSION}",
        f"output_h5ad={OUTPUT_H5AD}",
        f"n_cells={adata.n_obs:,}",
        f"n_genes_hvg={adata.n_vars:,}",
        f"n_genes_raw={adata.raw.n_vars if adata.raw is not None else 'None'}",
        f"batch_key={BATCH_KEY}",
        f"n_batches={adata.obs[BATCH_KEY].nunique()}",
        f"accelerator={accelerator}",
        f"devices={devices}",
        f"agreement_rate={agreement:.4f}",
        f"elapsed_minutes={elapsed_min:.2f}",
        "",
        "lineage_source distribution:",
        adata.obs["lineage_source"].value_counts().to_string(),
        "",
        "cell_type_final_l1 distribution:",
        adata.obs["cell_type_final_l1"].value_counts().to_string(),
        "",
        "scanvi_confidence summary:",
        adata.obs["scanvi_confidence"].describe().to_string(),
    ]
    RUN_LOG.write_text("\n".join(lines) + "\n", encoding="utf-8")



def main() -> None:
    setup_logging()
    try:
        np.random.seed(SEED)
        scvi.settings.seed = SEED
        torch.manual_seed(SEED)
        try:
            torch.set_num_threads(min(16, os.cpu_count() or 1))
        except Exception:
            pass

        accelerator, devices = resolve_accelerator(torch)
        gpu_available = accelerator == "gpu"
        print("=" * 90)
        print(f"[{PIPELINE_NAME}] start | version={VERSION} | date={DATE_TAG}")
        print(f"Output dir   : {OUTPUT_DIR}")
        print(f"Output h5ad  : {OUTPUT_H5AD}")
        print(f"Accelerator  : {accelerator} | devices={devices}")
        print("=" * 90)

        lineage_adatas: list[ad.AnnData] = []
        input_rows: list[dict[str, object]] = []

        for key, cfg in LINEAGE_INPUTS.items():
            h5ad_path = Path(cfg["h5ad"])
            l1_label = str(cfg["l1_label"])
            local_batch = str(cfg.get("batch_key_local", BATCH_KEY))

            print("\n" + "-" * 90)
            print(f"[LOAD] {key} -> {h5ad_path}")
            adata_lin = sc.read_h5ad(h5ad_path)
            print(f"  shape={adata_lin.n_obs:,} x {adata_lin.n_vars:,}")

            l2, l3, label_meta = build_l2_l3_labels(adata_lin, key)
            adata_lin.obs["cell_type_L2"] = l2.values
            adata_lin.obs["cell_type_L3"] = l3.values
            print(
                f"  labels: L1={l1_label} | L2={adata_lin.obs['cell_type_L2'].nunique()} classes | "
                f"L3={adata_lin.obs['cell_type_L3'].nunique()} classes"
            )

            normalize_lineage_structure(adata_lin, key, local_batch, BATCH_KEY)
            counts_x, count_source, raw_var = resolve_fullgene_counts_matrix(adata_lin, key)
            obs_full = ensure_covariates(adata_lin.obs, counts_x, pd.Index(raw_var.index.astype(str)), key)

            adata_full = ad.AnnData(X=counts_x, obs=obs_full.copy(), var=raw_var.copy())
            adata_full.layers["counts"] = adata_full.X
            adata_full.obs["lineage_source"] = l1_label
            adata_full.obs["lineage_branch"] = key
            adata_full.obs["cell_type_input_l1"] = l1_label
            adata_full.obs["cell_type_input_l2"] = clean_string_series(adata_lin.obs["cell_type_L2"], "Unknown").values
            adata_full.obs["cell_type_input_l3"] = clean_string_series(adata_lin.obs["cell_type_L3"], UNLABELED_CATEGORY).values
            adata_full.obs[LABELS_KEY] = adata_full.obs["cell_type_input_l1"]
            adata_full.var_names = adata_full.var_names.astype(str)

            lineage_adatas.append(adata_full)
            input_rows.append(
                {
                    "lineage_branch": key,
                    "lineage_source": l1_label,
                    "h5ad": str(h5ad_path),
                    "n_obs": int(adata_full.n_obs),
                    "n_vars_raw_source": int(adata_full.n_vars),
                    "count_source": count_source,
                    "l2_source": label_meta["l2_source"],
                    "l3_source": label_meta["l3_source"],
                }
            )
            print(f"  prepared={adata_full.n_obs:,} x {adata_full.n_vars:,} | count_source={count_source}")

            del adata_lin
            gc.collect()

        print("\n" + "=" * 90)
        print("[CONCAT] outer-join full raw-count space")
        print("=" * 90)
        adata = ad.concat(
            lineage_adatas,
            axis=0,
            join="outer",
            merge="first",
            uns_merge="first",
            label="concat_key",
            keys=list(LINEAGE_INPUTS.keys()),
            index_unique="__",
            fill_value=0,
        )
        del lineage_adatas
        gc.collect()
        if sparse.issparse(adata.X):
            adata.X = adata.X.tocsr().astype(np.float32)
        else:
            adata.X = sparse.csr_matrix(np.asarray(adata.X, dtype=np.float32))
        adata.layers["counts"] = adata.X
        adata.var_names = adata.var_names.astype(str)
        print(f"  concat shape={adata.n_obs:,} x {adata.n_vars:,}")

        adata.obs[BATCH_KEY] = to_clean_category(adata.obs[BATCH_KEY], "UnknownSample")
        adata.obs["lineage_source"] = to_clean_category(adata.obs["lineage_source"], "UnknownLineage")
        adata.obs["lineage_branch"] = to_clean_category(adata.obs["lineage_branch"], "UnknownBranch")
        adata.obs["cell_type_input_l1"] = pd.Categorical(clean_string_series(adata.obs["cell_type_input_l1"], UNLABELED_CATEGORY), categories=L1_ORDER)
        adata.obs["cell_type_input_l2"] = to_clean_category(adata.obs["cell_type_input_l2"], "Unknown")
        adata.obs["cell_type_input_l3"] = to_clean_category(adata.obs["cell_type_input_l3"], UNLABELED_CATEGORY)
        adata.obs["cell_type_L2"] = adata.obs["cell_type_input_l2"]
        adata.obs["cell_type_L3"] = adata.obs["cell_type_input_l3"]
        adata.obs[LABELS_KEY] = pd.Categorical(
            clean_string_series(adata.obs["cell_type_input_l1"], UNLABELED_CATEGORY),
            categories=L1_ORDER,
        )

        print("\n" + "=" * 90)
        print("[FILTER] genes/cells/batches")
        print("=" * 90)
        n_obs0, n_vars0 = adata.n_obs, adata.n_vars
        sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
        sc.pp.filter_cells(adata, min_genes=MIN_GENES_PER_CELL)
        print(f"  genes: {n_vars0:,} -> {adata.n_vars:,}")
        print(f"  cells: {n_obs0:,} -> {adata.n_obs:,}")
        batch_counts = adata.obs[BATCH_KEY].value_counts()
        small_batches = batch_counts[batch_counts < MIN_CELLS_PER_BATCH].index.tolist()
        if small_batches:
            n_drop = int(adata.obs[BATCH_KEY].isin(small_batches).sum())
            print(f"  dropping {len(small_batches)} small batches | cells={n_drop:,}")
            adata = adata[~adata.obs[BATCH_KEY].isin(small_batches)].copy()
        adata.obs[BATCH_KEY] = pd.Categorical(adata.obs[BATCH_KEY])
        adata.obs["lineage_source"] = pd.Categorical(adata.obs["lineage_source"])
        adata.obs["lineage_branch"] = pd.Categorical(adata.obs["lineage_branch"])
        adata.obs[LABELS_KEY] = pd.Categorical(adata.obs[LABELS_KEY], categories=L1_ORDER)
        if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
            adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])

        print("\n" + "=" * 90)
        print("[RAW] store full-gene raw counts into .raw")
        print("=" * 90)
        raw_var = adata.var.copy()
        bool_cols = raw_var.select_dtypes(include=["bool"]).columns.tolist()
        if bool_cols:
            raw_var[bool_cols] = raw_var[bool_cols].astype(np.int8)
        adata.raw = ad.AnnData(X=adata.layers["counts"], obs=adata.obs.copy(), var=raw_var)
        print(f"  .raw genes={adata.raw.n_vars:,}")

        print("\n" + "=" * 90)
        print("[HVG] select highly variable genes")
        print("=" * 90)
        adata.X = adata.layers["counts"]
        hvg_method = None
        hvg_errors: list[str] = []
        for flavor, batch_key in [("seurat_v3", BATCH_KEY), ("cell_ranger", None), ("seurat", None)]:
            try:
                kwargs = {"n_top_genes": HVG_TOP_GENES, "flavor": flavor, "subset": False, "inplace": True}
                if batch_key is not None:
                    kwargs["batch_key"] = batch_key
                sc.pp.highly_variable_genes(adata, **kwargs)
                if "highly_variable" in adata.var.columns and int(adata.var["highly_variable"].sum()) > 0:
                    hvg_method = f"{flavor}|batch={batch_key if batch_key is not None else 'None'}"
                    break
            except Exception as exc:
                hvg_errors.append(f"{flavor}|batch={batch_key}: {type(exc).__name__}: {exc}")
                drop_cols = [c for c in _HVG_AUX_VAR_COLS if c in adata.var.columns]
                if drop_cols:
                    adata.var.drop(columns=drop_cols, inplace=True)
        if hvg_method is None:
            raise RuntimeError(f"HVG selection failed: {hvg_errors}")
        print(f"  hvg_method={hvg_method} | n_hvg={int(adata.var['highly_variable'].sum()):,}")
        adata.uns["hvg_method"] = hvg_method
        adata = adata[:, adata.var["highly_variable"]].copy()
        print(f"  after HVG subset={adata.n_obs:,} x {adata.n_vars:,}")

        print("\n" + "=" * 90)
        print("[LOG1P] normalize HVG matrix into .X and layers['log1p']")
        print("=" * 90)
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["log1p"] = adata.X.copy()

        print("\n" + "=" * 90)
        print("[scVI] setup and train")
        print("=" * 90)
        scvi.model.SCVI.setup_anndata(
            adata,
            layer="counts",
            batch_key=BATCH_KEY,
            continuous_covariate_keys=CONTINUOUS_COVARIATES,
        )
        model_scvi = scvi.model.SCVI(
            adata,
            n_latent=SCVI_N_LATENT,
            n_layers=SCVI_N_LAYERS,
            n_hidden=SCVI_N_HIDDEN,
            dropout_rate=SCVI_DROPOUT_RATE,
            gene_likelihood="nb",
            encode_covariates=SCVI_ENCODE_COVARIATES,
            use_layer_norm=SCVI_USE_LAYER_NORM,
            use_batch_norm=SCVI_USE_BATCH_NORM,
        )
        model_scvi.train(
            max_epochs=SCVI_MAX_EPOCHS,
            batch_size=SCVI_BATCH_SIZE,
            train_size=0.9,
            early_stopping=SCVI_EARLY_STOPPING,
            early_stopping_patience=SCVI_EARLY_STOPPING_PATIENCE,
            plan_kwargs={"lr": SCVI_LEARNING_RATE},
            accelerator=accelerator,
            devices=devices,
        )
        adata.obsm["X_scvi"] = model_scvi.get_latent_representation().astype(np.float32)
        scvi_model_dir = MODEL_DIR / "scvi_model"
        model_scvi.save(str(scvi_model_dir), overwrite=True)
        pd.Series(adata.var_names).to_csv(scvi_model_dir / "var_names.csv", index=False, header=False)
        print(f"  scVI saved -> {scvi_model_dir}")

        print("\n" + "=" * 90)
        print("[scANVI] L1 fully supervised train")
        print("=" * 90)
        if UNLABELED_CATEGORY not in adata.obs[LABELS_KEY].cat.categories:
            adata.obs[LABELS_KEY] = adata.obs[LABELS_KEY].cat.add_categories([UNLABELED_CATEGORY])
        model_scanvi = scvi.model.SCANVI.from_scvi_model(
            model_scvi,
            unlabeled_category=UNLABELED_CATEGORY,
            labels_key=LABELS_KEY,
        )
        model_scanvi.train(
            max_epochs=SCANVI_MAX_EPOCHS,
            batch_size=SCANVI_BATCH_SIZE,
            train_size=0.9,
            early_stopping=True,
            early_stopping_patience=SCANVI_EARLY_STOPPING_PATIENCE,
            plan_kwargs={"lr": SCANVI_LEARNING_RATE, "weight_decay": 0.0},
            accelerator=accelerator,
            devices=devices,
        )
        adata.obsm["X_scanvi"] = model_scanvi.get_latent_representation().astype(np.float32)
        soft_pred = model_scanvi.predict(soft=True)
        if isinstance(soft_pred, pd.DataFrame):
            soft_df = soft_pred.reindex(adata.obs_names)
        else:
            columns = list(adata.obs[LABELS_KEY].cat.categories[: np.asarray(soft_pred).shape[1]])
            soft_df = pd.DataFrame(np.asarray(soft_pred), index=adata.obs_names, columns=columns)
        adata.obs["cell_type_scanvi_pred_pan"] = soft_df.idxmax(axis=1).astype(str)
        adata.obs["scanvi_confidence"] = soft_df.max(axis=1).astype(np.float32)
        adata.obs["cell_type_final_l1"] = adata.obs["cell_type_scanvi_pred_pan"].astype(object)
        adata.obsm["scanvi_probabilities"] = soft_df.values.astype(np.float32)
        adata.uns["scanvi_celltype_order"] = soft_df.columns.tolist()
        scanvi_model_dir = MODEL_DIR / "scanvi_model"
        model_scanvi.save(str(scanvi_model_dir), overwrite=True)
        pd.Series(adata.var_names).to_csv(scanvi_model_dir / "var_names.csv", index=False, header=False)
        print(f"  scANVI saved -> {scanvi_model_dir}")

        agree = (
            adata.obs["cell_type_final_l1"].astype(str)
            == adata.obs["cell_type_input_l1"].astype(str)
        )
        agreement_rate = float(agree.mean())
        adata.uns["agreement_rate_overall"] = agreement_rate
        adata.uns["agreement_rate_is_resubstitution"] = True
        print(f"  L1 agreement (resubstitution) = {agreement_rate * 100:.2f}%")

        print("\n" + "=" * 90)
        print("[UMAP] scanvi latent")
        print("=" * 90)
        sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=20)
        sc.tl.umap(adata, min_dist=0.3)
        adata.obsm["X_umap_scanvi"] = adata.obsm["X_umap"].copy().astype(np.float32)
        print("  X_umap_scanvi ready")

        print("\n" + "=" * 90)
        print("[CLEANUP] pre-write sanitization")
        print("=" * 90)
        drop_obs = [c for c in adata.obs.columns if c in _SCVI_INTERNAL_OBS_COLS]
        if drop_obs:
            adata.obs.drop(columns=drop_obs, inplace=True)
            print(f"  removed scvi internal obs cols: {drop_obs}")
        keep_obsm = {"X_scvi", "X_scanvi", "X_umap", "X_umap_scanvi", "scanvi_probabilities"}
        stale_obsm = [k for k in list(adata.obsm.keys()) if k not in keep_obsm]
        for key in stale_obsm:
            del adata.obsm[key]
        if stale_obsm:
            print(f"  removed stale obsm: {stale_obsm}")

        sanitize_category_columns(adata.obs)
        sanitize_category_columns(adata.var)
        notes = []
        notes.extend(sanitize_object_cols(adata.obs, "obs"))
        notes.extend(sanitize_object_cols(adata.var, "var"))
        if adata.raw is not None:
            notes.extend(sanitize_object_cols(adata.raw.var, "raw.var"))
        if notes:
            print("  sanitize notes:")
            for note in notes:
                print(f"    - {note}")

        save_reference_artifacts(adata, input_rows, accelerator, devices)

        print("\n" + "=" * 90)
        print(f"[WRITE] {OUTPUT_H5AD}")
        print("=" * 90)
        adata.write_h5ad(OUTPUT_H5AD, compression="gzip", compression_opts=4)
        print(f"  saved h5ad size={OUTPUT_H5AD.stat().st_size / 1e9:.2f} GB")

        if gpu_available:
            try:
                torch.cuda.empty_cache()
            except Exception:
                pass

        elapsed_min = (time.time() - PIPELINE_START) / 60.0
        print("\n" + "=" * 90)
        print("[DONE] pipeline complete")
        print(f"Elapsed minutes : {elapsed_min:.2f}")
        print(f"Output h5ad     : {OUTPUT_H5AD}")
        print(f"Manifest        : {MANIFEST_JSON}")
        print(f"Structure json  : {STRUCTURE_JSON}")
        print(f"Input tsv       : {INPUT_TSV}")
        print(f"Log file        : {LOG_FILE}")
        print("=" * 90)

    except Exception as exc:
        print("\n" + "!" * 90)
        print(f"[FAIL] {type(exc).__name__}: {exc}")
        traceback.print_exc()
        print("!" * 90)
        raise
    finally:
        close_logging()


if __name__ == "__main__":
    main()
