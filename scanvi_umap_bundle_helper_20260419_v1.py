#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Unified scanVI/scVI UMAP bundle helper.

This helper standardizes three related workflows that had diverged across
multiple lineage-specific scripts:

1. ``fit_bundle``: fit a UMAP operator from a latent representation, store the
   resulting coordinates back to AnnData, and persist a reusable operator.
2. ``refresh_from_existing_latent``: recompute a standardized scanVI UMAP from
   an existing AnnData object without retraining the model.
3. ``project_query``: project query latent vectors with a previously saved
   operator so reference/query visualizations stay in the same coordinate space.

The module is intentionally import-safe: heavy optional dependencies such as
``umap`` and ``joblib`` are loaded lazily inside functions.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

DEFAULT_UMAP_PARAMS: dict[str, Any] = {
    "n_neighbors": 30,
    "n_components": 2,
    "min_dist": 0.5,
    "spread": 1.0,
    "metric": "euclidean",
    "random_state": 42,
}

DEFAULT_LATENT_CANDIDATES: tuple[str, ...] = (
    "X_scANVI",
    "X_scanvi",
    "X_scANVI_L2",
    "X_scanvi_refined",
    "X_scvi",
)

DEFAULT_UMAP_CANDIDATES: tuple[str, ...] = (
    "X_umap_scanvi_corrected",
    "X_umap_scanvi",
    "X_umap_scANVI",
    "X_umap_refined",
    "X_umap_scVI",
    "X_umap_scvi",
    "X_umap",
)

DEFAULT_METADATA_KEY = "scanvi_umap_bundle"


def _load_umap_class():
    from umap import UMAP  # type: ignore

    return UMAP


def _load_joblib():
    import joblib  # type: ignore

    return joblib


def _ensure_path(path_like: str | Path) -> Path:
    path = Path(path_like)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _available_keys(keys_or_adata: Any, attr: str = "obsm") -> list[str]:
    if hasattr(keys_or_adata, attr):
        mapping = getattr(keys_or_adata, attr)
        return list(mapping.keys())
    return [str(k) for k in keys_or_adata]


def resolve_latent_key(
    keys_or_adata: Any,
    preferred: str | None = None,
    candidates: Sequence[str] | None = None,
) -> str:
    available = _available_keys(keys_or_adata, attr="obsm")
    search = [preferred] if preferred else []
    search.extend(list(candidates or DEFAULT_LATENT_CANDIDATES))
    for key in search:
        if key and key in available:
            return key
    raise KeyError(f"No latent key found. Available keys: {available}")


def resolve_umap_key(
    keys_or_adata: Any,
    preferred: str | None = None,
    candidates: Sequence[str] | None = None,
) -> str:
    available = _available_keys(keys_or_adata, attr="obsm")
    search = [preferred] if preferred else []
    search.extend(list(candidates or DEFAULT_UMAP_CANDIDATES))
    for key in search:
        if key and key in available:
            return key
    raise KeyError(f"No UMAP key found. Available keys: {available}")


def clear_neighbors_and_umap(adata: Any) -> None:
    for key in ("neighbors", "umap"):
        if key in getattr(adata, "uns", {}):
            del adata.uns[key]
    for key in ("connectivities", "distances"):
        if key in getattr(adata, "obsp", {}):
            del adata.obsp[key]


def _dedupe_keys(keys: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for key in keys:
        key = str(key)
        if not key or key in seen:
            continue
        seen.add(key)
        ordered.append(key)
    return ordered


def enable_h5ad_string_compat() -> None:
    try:
        import anndata  # type: ignore

        anndata.settings.allow_write_nullable_strings = True
    except Exception:
        pass


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return {
            "type": "ndarray",
            "shape": [int(x) for x in value.shape],
            "dtype": str(value.dtype),
        }
    if isinstance(value, Mapping):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def store_umap_coordinates(
    adata: Any,
    coords: np.ndarray,
    *,
    primary_umap_key: str = "X_umap_scanvi",
    alias_keys: Sequence[str] = (),
    set_default_x_umap: bool = True,
) -> dict[str, Any]:
    coords = np.asarray(coords, dtype=np.float32)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f"Expected UMAP coordinates with shape (n_obs, 2); got {coords.shape}")
    if coords.shape[0] != getattr(adata, "n_obs"):
        raise ValueError(
            f"Coordinate rows ({coords.shape[0]}) do not match adata.n_obs ({adata.n_obs})"
        )

    written_keys = _dedupe_keys([primary_umap_key, *alias_keys])
    for key in written_keys:
        adata.obsm[key] = coords.copy()
    if set_default_x_umap:
        adata.obsm["X_umap"] = coords.copy()
        if "X_umap" not in written_keys:
            written_keys.append("X_umap")

    return {
        "primary_umap_key": primary_umap_key,
        "written_umap_keys": written_keys,
        "set_default_x_umap": bool(set_default_x_umap),
        "n_obs": int(coords.shape[0]),
        "n_components": int(coords.shape[1]),
    }


def build_bundle_manifest(
    adata: Any,
    *,
    mode: str,
    latent_key: str,
    primary_umap_key: str,
    operator_path: str | Path | None,
    umap_params: Mapping[str, Any],
    written_umap_keys: Sequence[str],
    source_h5ad: str | Path | None = None,
    metadata_key: str = DEFAULT_METADATA_KEY,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    manifest: dict[str, Any] = {
        "mode": mode,
        "latent_key": latent_key,
        "primary_umap_key": primary_umap_key,
        "written_umap_keys": list(written_umap_keys),
        "operator_path": None if operator_path is None else str(operator_path),
        "umap_params": dict(umap_params),
        "n_obs": int(getattr(adata, "n_obs", -1)),
        "n_vars": int(getattr(adata, "n_vars", -1)),
        "obsm_keys": _available_keys(adata, attr="obsm"),
        "metadata_key": metadata_key,
    }
    if source_h5ad is not None:
        manifest["source_h5ad"] = str(source_h5ad)
    if hasattr(adata, "obs") and "data_source" in adata.obs.columns:
        manifest["data_source_counts"] = {
            str(k): int(v)
            for k, v in adata.obs["data_source"].astype(str).value_counts(dropna=False).items()
        }
    if extra:
        manifest.update({str(k): _json_safe(v) for k, v in extra.items()})
    return manifest


def _write_manifest(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(_json_safe(dict(payload)), handle, indent=2, ensure_ascii=False)


def _fit_umap(coords: np.ndarray, umap_params: Mapping[str, Any]):
    UMAP = _load_umap_class()
    operator = UMAP(**dict(umap_params))
    embedding = operator.fit_transform(coords)
    return operator, np.asarray(embedding, dtype=np.float32)


def fit_bundle(
    adata: Any,
    latent_key: str,
    output_dir: str | Path,
    *,
    primary_umap_key: str = "X_umap_scanvi",
    alias_keys: Sequence[str] = (),
    operator_filename: str = "umap_operator_scanvi.joblib",
    manifest_filename: str = "scanvi_umap_bundle.json",
    umap_params: Mapping[str, Any] | None = None,
    set_default_x_umap: bool = True,
    metadata_key: str = DEFAULT_METADATA_KEY,
    source_h5ad: str | Path | None = None,
    extra_manifest: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if latent_key not in adata.obsm:
        raise KeyError(f"Latent key not found in adata.obsm: {latent_key}")

    clear_neighbors_and_umap(adata)
    latent = np.asarray(adata.obsm[latent_key], dtype=np.float32)
    operator, embedding = _fit_umap(latent, umap_params or DEFAULT_UMAP_PARAMS)
    store_info = store_umap_coordinates(
        adata,
        embedding,
        primary_umap_key=primary_umap_key,
        alias_keys=alias_keys,
        set_default_x_umap=set_default_x_umap,
    )

    operator_path = output_dir / operator_filename
    _load_joblib().dump(operator, operator_path)

    manifest = build_bundle_manifest(
        adata,
        mode="fit_bundle",
        latent_key=latent_key,
        primary_umap_key=primary_umap_key,
        operator_path=operator_path,
        umap_params=umap_params or DEFAULT_UMAP_PARAMS,
        written_umap_keys=store_info["written_umap_keys"],
        source_h5ad=source_h5ad,
        metadata_key=metadata_key,
        extra=extra_manifest,
    )
    adata.uns[metadata_key] = dict(manifest)
    manifest_path = output_dir / manifest_filename
    _write_manifest(manifest_path, manifest)

    return {
        "operator": operator,
        "operator_path": operator_path,
        "manifest": manifest,
        "manifest_path": manifest_path,
        **store_info,
    }


def refresh_from_existing_latent(
    adata: Any,
    *,
    latent_key: str | None = None,
    output_dir: str | Path,
    primary_umap_key: str = "X_umap_scanvi",
    alias_keys: Sequence[str] = (),
    operator_filename: str = "umap_operator_scanvi.joblib",
    manifest_filename: str = "scanvi_umap_refresh.json",
    umap_params: Mapping[str, Any] | None = None,
    set_default_x_umap: bool = True,
    metadata_key: str = "scanvi_umap_refresh",
    source_h5ad: str | Path | None = None,
    write_sidecar_h5ad: bool = False,
    sidecar_h5ad_path: str | Path | None = None,
    sidecar_suffix: str = "_scanvi_umap_refresh",
    extra_manifest: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    resolved_latent_key = resolve_latent_key(adata, preferred=latent_key)
    result = fit_bundle(
        adata,
        resolved_latent_key,
        output_dir,
        primary_umap_key=primary_umap_key,
        alias_keys=alias_keys,
        operator_filename=operator_filename,
        manifest_filename=manifest_filename,
        umap_params=umap_params,
        set_default_x_umap=set_default_x_umap,
        metadata_key=metadata_key,
        source_h5ad=source_h5ad,
        extra_manifest={"refresh_mode": True, **(dict(extra_manifest or {}))},
    )
    result["manifest"]["mode"] = "refresh_from_existing_latent"
    adata.uns[metadata_key]["mode"] = "refresh_from_existing_latent"
    _write_manifest(result["manifest_path"], result["manifest"])

    sidecar_path = None
    if write_sidecar_h5ad:
        if sidecar_h5ad_path is not None:
            sidecar_path = Path(sidecar_h5ad_path)
        elif source_h5ad is not None:
            src = Path(source_h5ad)
            sidecar_path = src.with_name(src.stem + f"{sidecar_suffix}.h5ad")
        else:
            raise ValueError("write_sidecar_h5ad=True requires sidecar_h5ad_path or source_h5ad")
        enable_h5ad_string_compat()
        adata.write_h5ad(sidecar_path, compression="gzip", compression_opts=9)
    result["sidecar_h5ad_path"] = sidecar_path
    return result


def project_query(
    query_latent: np.ndarray,
    *,
    operator: Any | None = None,
    operator_path: str | Path | None = None,
) -> np.ndarray:
    if operator is None and operator_path is None:
        raise ValueError("Provide either operator or operator_path")
    if operator is None:
        operator = _load_joblib().load(operator_path)
    projected = operator.transform(np.asarray(query_latent, dtype=np.float32))
    return np.asarray(projected, dtype=np.float32)
