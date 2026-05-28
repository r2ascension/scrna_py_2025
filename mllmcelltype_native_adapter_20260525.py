#!/usr/bin/env python3
"""Native Python adapter for running mLLMCelltype inside the scExtract workflow."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
import scanpy as sc

from scextract_three_way_annotation_workflow_20260524 import (
    deep_get,
    ensure_cluster_key,
    ensure_counts_layer,
    ensure_dir,
    infer_lineage_name_from_source_path,
    normalize_gene_names,
    resolve_cluster_key,
    stringify,
    write_json,
)


LINEAGE_CONTEXT_HINTS: dict[str, str] = {
    "bcell": (
        "This AnnData appears to be a lineage-restricted B/plasma-cell subset from lung/airway tissue. "
        "Expected labels should usually stay within B-cell, plasma-cell, plasmablast, activated/cycling lymphocyte, "
        "or closely related immune states. Only assign unrelated epithelial, stromal, endothelial, or myeloid labels "
        "if the marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "tnk": (
        "This AnnData appears to be a lineage-restricted T/NK-cell subset from lung/airway tissue. "
        "Expected labels should usually stay within T-cell, NK-cell, cytotoxic, exhausted, activated, regulatory, "
        "innate-like, or cycling lymphocyte states. Only assign unrelated non-lymphoid labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "myeloid": (
        "This AnnData appears to be a lineage-restricted myeloid subset from lung/airway tissue. "
        "Expected labels should usually stay within monocyte, macrophage, dendritic, neutrophil, mast-cell, or cycling myeloid states. "
        "Only assign unrelated epithelial, stromal, endothelial, or lymphoid labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "epithelial": (
        "This AnnData appears to be a lineage-restricted epithelial subset from lung/airway tissue. "
        "Expected labels should usually stay within airway/alveolar epithelial programs and their proliferative or transitional states. "
        "Only assign unrelated immune or stromal labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "stromal_endothelial": (
        "This AnnData appears to be a lineage-restricted endothelial/vascular subset from lung/airway tissue. "
        "Expected labels should usually stay within endothelial, vascular, lymphatic, angiogenic, or cycling endothelial states. "
        "Only assign unrelated epithelial or lymphoid labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "stromal_fibroblast": (
        "This AnnData appears to be a lineage-restricted fibroblast/stromal subset from lung/airway tissue. "
        "Expected labels should usually stay within fibroblast, mesenchymal, matrix-remodeling, myofibroblast, or cycling stromal states. "
        "Only assign unrelated epithelial or immune labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
    "stromal_smc": (
        "This AnnData appears to be a lineage-restricted smooth-muscle/pericyte stromal subset from lung/airway tissue. "
        "Expected labels should usually stay within smooth-muscle, pericyte, contractile stromal, vascular support, or cycling stromal states. "
        "Only assign unrelated epithelial or immune labels if marker evidence is overwhelming; otherwise prefer Unknown."
    ),
}


def _normalize_model_spec(spec: Any) -> str | dict[str, str] | None:
    if spec is None:
        return None
    if isinstance(spec, str):
        text = spec.strip()
        return text or None
    if isinstance(spec, dict):
        provider = stringify(spec.get("provider", "")).strip()
        model = stringify(spec.get("model", "")).strip()
        out: dict[str, str] = {}
        if provider:
            out["provider"] = provider
        if model:
            out["model"] = model
        return out or None
    return None


def _resolve_base_urls(config: dict[str, Any], workflow_config: dict[str, Any]) -> dict[str, str] | None:
    provider = stringify(deep_get(config, "provider", default="deepseek")).strip().lower() or "deepseek"
    base_url = stringify(
        deep_get(config, "base_url", default=deep_get(workflow_config, "deepseek", "base_url", default=""))
    ).strip()
    if not base_url:
        return None
    base_url = base_url.rstrip("/")
    if provider == "deepseek" and not base_url.endswith("/chat/completions"):
        if base_url.endswith("/v1"):
            base_url = f"{base_url}/chat/completions"
        elif base_url == "https://api.deepseek.com":
            base_url = f"{base_url}/v1/chat/completions"
        else:
            base_url = f"{base_url}/chat/completions"
    return {provider: base_url}


def _extract_marker_genes(
    adata: ad.AnnData,
    config: dict[str, Any],
    workflow_config: dict[str, Any],
) -> tuple[dict[str, list[str]], dict[str, Any]]:
    counts_layer = stringify(deep_get(workflow_config, "data_contract", "counts_layer", default="counts")).strip() or "counts"
    configured_cluster_key = stringify(deep_get(config, "cluster_key", default="")).strip()
    cluster_key_fallback_allowed = bool(deep_get(config, "allow_cluster_key_fallback", default=True))
    auto_compute_cluster_key = bool(deep_get(workflow_config, "data_contract", "auto_compute_cluster_key", default=True))
    top_n_markers = int(deep_get(config, "top_n_markers", default=10))
    marker_method = stringify(deep_get(config, "marker_method", default="wilcoxon")).strip() or "wilcoxon"

    normalize_gene_names(adata)
    counts_info = ensure_counts_layer(adata, counts_layer)

    if configured_cluster_key and configured_cluster_key.lower() not in {"auto", "null", "none"}:
        if configured_cluster_key in adata.obs.columns:
            cluster_key = configured_cluster_key
            cluster_key_source = "configured"
        elif cluster_key_fallback_allowed:
            cluster_key, resolved_source = resolve_cluster_key(adata, workflow_config)
            cluster_key_source = f"configured_missing_{resolved_source}"
        else:
            cluster_key = configured_cluster_key
            cluster_key_source = "configured_missing_to_be_computed"
    else:
        cluster_key, cluster_key_source = resolve_cluster_key(adata, workflow_config)

    cluster_info = ensure_cluster_key(adata, cluster_key, auto_compute=auto_compute_cluster_key)

    marker_matrix = adata.layers[counts_layer] if counts_layer in adata.layers else adata.X
    marker_obs = adata.obs.loc[:, [cluster_key]].copy()
    if "symbol_base" in adata.var.columns:
        marker_var = adata.var.loc[:, ["symbol_base"]].copy()
    else:
        marker_var = pd.DataFrame(index=adata.var_names.copy())

    marker_adata = ad.AnnData(
        X=marker_matrix.copy(),
        obs=marker_obs,
        var=marker_var,
    )

    if "symbol_base" in marker_adata.var.columns:
        symbols = marker_adata.var["symbol_base"].astype(str).str.strip()
        keep = (symbols != "") & (~symbols.duplicated(keep="first"))
        if int(keep.sum()) > 0:
            marker_adata = marker_adata[:, keep].copy()
            marker_adata.var_names = pd.Index(symbols.loc[keep].astype(str).values)

    if marker_adata.var_names.duplicated().sum() > 0:
        marker_adata.var_names_make_unique()

    marker_adata.obs[cluster_key] = marker_adata.obs[cluster_key].astype(str).astype("category")
    sc.pp.normalize_total(marker_adata, target_sum=1e4)
    sc.pp.log1p(marker_adata)
    sc.tl.rank_genes_groups(marker_adata, groupby=cluster_key, method=marker_method, use_raw=False)

    marker_genes: dict[str, list[str]] = {}
    for cluster in marker_adata.obs[cluster_key].cat.categories:
        df = sc.get.rank_genes_groups_df(marker_adata, group=cluster)
        genes: list[str] = []
        seen: set[str] = set()
        for gene in df["names"].astype(str):
            gene = gene.strip()
            if not gene or gene in seen:
                continue
            genes.append(gene)
            seen.add(gene)
            if len(genes) >= top_n_markers:
                break
        if genes:
            marker_genes[str(cluster)] = genes

    if not marker_genes:
        raise RuntimeError("mLLMCelltype marker extraction produced no usable marker genes")

    summary = {
        "cluster_key": cluster_key,
        "cluster_key_source": cluster_key_source,
        "cluster_key_created": bool(cluster_info.get("created", False)),
        "n_clusters": int(len(marker_genes)),
        "top_n_markers": top_n_markers,
        "marker_method": marker_method,
        "counts_layer_source": counts_info.get("source"),
    }
    return marker_genes, summary


def _resolve_models_for_consensus(config: dict[str, Any]) -> list[str | dict[str, str]]:
    configured = deep_get(config, "models", default=None) or []
    models = [_normalize_model_spec(item) for item in configured]
    models = [item for item in models if item is not None]
    if models:
        return models

    provider = stringify(deep_get(config, "provider", default="deepseek")).strip().lower() or "deepseek"
    primary_model = stringify(deep_get(config, "model", default="deepseek-v4-pro")).strip() or "deepseek-v4-pro"
    fallback_model = stringify(deep_get(config, "fallback_model", default="deepseek-v4-flash")).strip()
    resolved: list[str | dict[str, str]] = [{"provider": provider, "model": primary_model}]
    if fallback_model and fallback_model != primary_model:
        resolved.append({"provider": provider, "model": fallback_model})
    return resolved


def _resolve_lineage_context(
    input_path: Path,
    config: dict[str, Any],
    workflow_config: dict[str, Any],
) -> tuple[str | None, str | None, str | None]:
    configured_lineage = stringify(
        deep_get(
            config,
            "lineage_name",
            default=deep_get(workflow_config, "run", "lineage_name", default=""),
        )
    ).strip()
    configured_lineage = configured_lineage if configured_lineage.lower() not in {"", "auto", "null", "none"} else ""

    if configured_lineage:
        lineage_name = configured_lineage
        lineage_source = "configured"
    else:
        lineage_name = infer_lineage_name_from_source_path(input_path)
        lineage_source = "inferred_from_input_path" if lineage_name else None

    if not lineage_name:
        return None, None, None

    lineage_hint = LINEAGE_CONTEXT_HINTS.get(lineage_name)
    return lineage_name, lineage_source, lineage_hint


def _provider_names_from_models(config: dict[str, Any], models: list[str | dict[str, str]]) -> list[str]:
    providers: list[str] = []
    default_provider = stringify(deep_get(config, "provider", default="deepseek")).strip().lower() or "deepseek"
    for item in models:
        if isinstance(item, dict):
            provider = stringify(item.get("provider", "")).strip().lower()
            if provider:
                providers.append(provider)
        else:
            providers.append(default_provider)
    unique_providers: list[str] = []
    for provider in providers:
        if provider and provider not in unique_providers:
            unique_providers.append(provider)
    return unique_providers


def _resolve_mllm_result(
    result: Any,
    package_mode: str,
) -> tuple[dict[str, str], dict[str, float] | None, dict[str, float] | None, dict[str, Any]]:
    details: dict[str, Any] = {"package_mode": package_mode}

    if package_mode in {"single_model", "annotate_clusters"}:
        if not isinstance(result, dict):
            raise TypeError(f"mLLMCelltype annotate_clusters() returned unexpected type: {type(result)!r}")
        label_map = {str(k): stringify(v).strip() for k, v in result.items() if stringify(v).strip()}
        confidence_map = {cluster: 1.0 for cluster in label_map}
        return label_map, confidence_map, None, details

    if not isinstance(result, dict):
        raise TypeError(f"mLLMCelltype interactive_consensus_annotation() returned unexpected type: {type(result)!r}")

    label_candidates = [result.get("consensus"), result.get("final_annotations"), result.get("annotations")]
    label_map: dict[str, str] | None = None
    for candidate in label_candidates:
        if isinstance(candidate, dict) and candidate:
            label_map = {str(k): stringify(v).strip() for k, v in candidate.items() if stringify(v).strip()}
            break
    if label_map is None:
        scalarish = all(not isinstance(v, (dict, list, tuple, set)) for v in result.values()) if result else False
        if scalarish:
            label_map = {str(k): stringify(v).strip() for k, v in result.items() if stringify(v).strip()}
    if not label_map:
        raise RuntimeError("Could not resolve cluster annotations from mLLMCelltype consensus result")

    confidence_map = None
    confidence_candidate = result.get("consensus_proportion")
    if isinstance(confidence_candidate, dict):
        confidence_map = {str(k): float(v) for k, v in confidence_candidate.items() if v is not None}

    score_map = None
    entropy_candidate = result.get("entropy")
    if isinstance(entropy_candidate, dict):
        score_map = {str(k): float(v) for k, v in entropy_candidate.items() if v is not None}
        details["score_metric"] = "entropy"

    return label_map, confidence_map, score_map, details


def annotate_with_mllmcelltype(
    adata_path: str | None = None,
    input_h5ad: str | None = None,
    output_h5ad: str | None = None,
    output_path: str | None = None,
    branch_dir: str | None = None,
    config: dict[str, Any] | None = None,
    workflow_config: dict[str, Any] | None = None,
    unknown_label: str = "Unknown",
    write_h5ad: bool = True,
    annotation_table_path: str | None = None,
) -> str | dict[str, Any]:
    config = config or {}
    workflow_config = workflow_config or {}

    input_path = Path(input_h5ad or adata_path or "")
    if not input_path.exists():
        raise FileNotFoundError(f"mLLMCelltype adapter input h5ad not found: {input_path}")

    resolved_output = Path(output_h5ad or output_path or "")
    if not resolved_output:
        raise ValueError("mLLMCelltype adapter requires output_h5ad or output_path")
    ensure_dir(resolved_output.parent)

    branch_path = Path(branch_dir) if branch_dir else resolved_output.parent
    ensure_dir(branch_path)

    from mllmcelltype import annotate_clusters, interactive_consensus_annotation

    env_key = stringify(deep_get(config, "env_key", default=deep_get(workflow_config, "deepseek", "env_key", default="DEEPSEEK_API_KEY"))).strip() or "DEEPSEEK_API_KEY"
    api_key = os.getenv(env_key)
    if not api_key:
        raise RuntimeError(f"Missing required environment variable for mLLMCelltype: {env_key}")

    adata = ad.read_h5ad(input_path)
    marker_genes, marker_summary = _extract_marker_genes(adata, config, workflow_config)
    write_json(marker_summary, branch_path / "mllmcelltype_marker_summary.json")
    write_json(marker_genes, branch_path / "mllmcelltype_marker_genes.json")

    package_mode = stringify(deep_get(config, "package_mode", default="consensus")).strip().lower() or "consensus"
    species = stringify(deep_get(config, "species", default="human")).strip() or "human"
    tissue = stringify(deep_get(config, "tissue", default=deep_get(workflow_config, "sctype", "tissue", default=""))).strip() or None
    additional_context = stringify(deep_get(config, "additional_context", default="")).strip() or None
    auto_lineage_context = bool(deep_get(config, "auto_lineage_context", default=True))
    lineage_name = None
    lineage_source = None
    lineage_hint = None
    if auto_lineage_context:
        lineage_name, lineage_source, lineage_hint = _resolve_lineage_context(input_path, config, workflow_config)
        if lineage_hint:
            additional_context = "\n\n".join([part for part in [additional_context, lineage_hint] if part]) or None
    use_cache = bool(deep_get(config, "use_cache", default=True))
    cache_dir_raw = stringify(deep_get(config, "cache_dir", default="")).strip()
    cache_dir = cache_dir_raw or str(branch_path / "cache")
    if use_cache:
        ensure_dir(Path(cache_dir))
    else:
        cache_dir = None
    base_urls = _resolve_base_urls(config, workflow_config)

    if package_mode in {"single_model", "annotate_clusters"}:
        provider = stringify(deep_get(config, "provider", default="deepseek")).strip().lower() or "deepseek"
        model = stringify(deep_get(config, "model", default="deepseek-v4-pro")).strip() or "deepseek-v4-pro"
        log_dir = str(ensure_dir(branch_path / "logs"))
        result = annotate_clusters(
            marker_genes=marker_genes,
            species=species,
            provider=provider,
            model=model,
            api_key=api_key,
            tissue=tissue,
            additional_context=additional_context,
            use_cache=use_cache,
            cache_dir=cache_dir,
            log_dir=log_dir,
            log_level=stringify(deep_get(config, "log_level", default="INFO")).strip() or "INFO",
            base_urls=base_urls,
        )
        invocation = {"package_mode": package_mode, "provider": provider, "model": model}
    else:
        models = _resolve_models_for_consensus(config)
        providers = _provider_names_from_models(config, models)
        api_keys = {provider: api_key for provider in providers}
        consensus_model = _normalize_model_spec(deep_get(config, "consensus_model", default=None))
        clusters_to_analyze = deep_get(config, "clusters_to_analyze", default=None)
        result = interactive_consensus_annotation(
            marker_genes=marker_genes,
            species=species,
            models=models,
            api_keys=api_keys,
            tissue=tissue,
            additional_context=additional_context,
            consensus_threshold=float(deep_get(config, "consensus_threshold", default=0.7)),
            entropy_threshold=float(deep_get(config, "entropy_threshold", default=1.0)),
            max_discussion_rounds=int(deep_get(config, "max_discussion_rounds", default=3)),
            use_cache=use_cache,
            cache_dir=cache_dir,
            verbose=bool(deep_get(config, "verbose", default=False)),
            consensus_model=consensus_model,
            base_urls=base_urls,
            clusters_to_analyze=clusters_to_analyze,
            force_rerun=bool(deep_get(config, "force_rerun", default=False)),
        )
        invocation = {
            "package_mode": package_mode,
            "models": models,
            "consensus_model": consensus_model,
            "providers": providers,
        }

    label_map, confidence_map, score_map, result_summary = _resolve_mllm_result(result, package_mode)
    cluster_key = marker_summary["cluster_key"]
    cluster_ids = adata.obs[cluster_key].astype(str)

    labels = cluster_ids.map(label_map).fillna(unknown_label).astype(str)
    obs_out = pd.DataFrame(index=adata.obs_names.astype(str))
    obs_out["mllmcelltype_annotation"] = labels.values

    if confidence_map is not None:
        obs_out["mllmcelltype_confidence_scextract"] = pd.to_numeric(cluster_ids.map(confidence_map), errors="coerce").values
    if score_map is not None:
        obs_out["mllmcelltype_score"] = pd.to_numeric(cluster_ids.map(score_map), errors="coerce").values

    annotation_table = Path(annotation_table_path) if annotation_table_path else branch_path / "mllmcelltype_annotations.tsv.gz"
    obs_out.to_csv(annotation_table, sep="\t", compression="gzip", index_label="obs_name")

    serializable_result = {
        "label_map": label_map,
        "confidence_map": confidence_map,
        "score_map": score_map,
        "invocation": invocation,
        "result_summary": result_summary,
        "marker_summary": marker_summary,
        "lineage_name": lineage_name,
        "lineage_source": lineage_source,
        "additional_context": additional_context,
        "annotation_table_path": str(annotation_table),
        "write_h5ad": bool(write_h5ad),
    }
    write_json(serializable_result, branch_path / "mllmcelltype_result_summary.json")

    if write_h5ad:
        adata.obs["mllmcelltype_annotation"] = obs_out["mllmcelltype_annotation"].values
        if "mllmcelltype_confidence_scextract" in obs_out.columns:
            adata.obs["mllmcelltype_confidence_scextract"] = pd.to_numeric(obs_out["mllmcelltype_confidence_scextract"], errors="coerce").values
        if "mllmcelltype_score" in obs_out.columns:
            adata.obs["mllmcelltype_score"] = pd.to_numeric(obs_out["mllmcelltype_score"], errors="coerce").values

        adata.uns["mllmcelltype_package_mode"] = package_mode
        adata.uns["mllmcelltype_cluster_key"] = cluster_key
        adata.uns["mllmcelltype_species"] = species
        adata.uns["mllmcelltype_tissue"] = tissue
        if lineage_name:
            adata.uns["mllmcelltype_lineage_name"] = lineage_name
        if lineage_source:
            adata.uns["mllmcelltype_lineage_source"] = lineage_source
        adata.uns["mllmcelltype_invocation"] = json.dumps(invocation, ensure_ascii=False)
        adata.uns["mllmcelltype_result_summary"] = json.dumps(result_summary, ensure_ascii=False)
        adata.write_h5ad(resolved_output, compression="gzip")
        return str(resolved_output)

    return {
        "annotation_table_path": str(annotation_table),
        "label_column": "mllmcelltype_annotation",
        "confidence_column": "mllmcelltype_confidence_scextract" if confidence_map is not None else None,
        "score_column": "mllmcelltype_score" if score_map is not None else None,
        "output_h5ad": None,
    }
