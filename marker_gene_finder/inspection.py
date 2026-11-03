"""Input inspection and execution planning."""
from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from .input_loading import estimate_dense_expression_bytes
from .version import APP_VERSION
from ._constants import (
    CSV_INSPECTION_CHUNK_ROWS,
    DEFAULT_GENE_CHUNK_SIZE,
    DEFAULT_MAX_DENSE_MEMORY_MB,
    PLAN_SCHEMA_VERSION,
)
from ._validation import (
    _validate_gene_chunk_size,
    _validate_max_dense_memory_mb,
    _validate_workers,
)


def inspect_input(
    input_path: str | Path,
    cluster_column: str = "cluster",
    cell_column: str | None = "cell",
    layer: str | None = None,
    use_raw: bool = False,
    max_dense_memory_mb: float | None = DEFAULT_MAX_DENSE_MEMORY_MB,
) -> dict[str, object]:
    """Inspect dataset dimensions and memory without running marker analysis."""
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()
    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")
    _validate_max_dense_memory_mb(max_dense_memory_mb)
    normalized_dense_limit = (
        None
        if max_dense_memory_mb is None
        else float(max_dense_memory_mb)
    )

    if suffix == ".h5ad":
        adata = ad.read_h5ad(input_path, backed="r")
        try:
            if cluster_column not in adata.obs.columns:
                raise ValueError(
                    f"Missing cluster column in AnnData.obs: {cluster_column!r}."
                )
            if adata.obs[cluster_column].isna().any():
                raise ValueError("Cluster labels cannot contain missing values.")
            if use_raw:
                if adata.raw is None or adata.raw.X is None:
                    raise ValueError(
                        "The AnnData object does not contain raw expression."
                    )
                matrix = adata.raw.X
                source = "raw"
            elif layer is not None:
                if layer not in adata.layers:
                    raise ValueError(f"Missing AnnData layer: {layer!r}.")
                matrix = adata.layers[layer]
                source = f"layer:{layer}"
            else:
                if adata.X is None:
                    raise ValueError(
                        "The AnnData object does not contain an expression matrix in X."
                    )
                matrix = adata.X
                source = "X"
            cells, genes = matrix.shape
            dense_bytes = (
                int(cells) * int(genes) * np.dtype(matrix.dtype).itemsize
            )
            storage = (
                "sparse"
                if sparse.issparse(matrix)
                or getattr(matrix, "format", None) in {"csr", "csc"}
                else "dense"
            )
            clusters = int(adata.obs[cluster_column].nunique())
        finally:
            adata.file.close()
    elif suffix == ".csv":
        if layer is not None or use_raw:
            raise ValueError("Layer and raw selection require H5AD input.")
        header = pd.read_csv(input_path, nrows=0)
        if cluster_column not in header.columns:
            raise ValueError(f"Missing cluster column: {cluster_column!r}.")
        metadata = {cluster_column}
        if cell_column and cell_column in header.columns:
            metadata.add(cell_column)
        gene_columns = [column for column in header.columns if column not in metadata]
        if not gene_columns:
            raise ValueError("No gene-expression columns were found.")
        cells = 0
        dense_bytes = 0
        observed_clusters: set[object] = set()
        for chunk in pd.read_csv(input_path, chunksize=CSV_INSPECTION_CHUNK_ROWS):
            cells += len(chunk)
            if chunk[cluster_column].isna().any():
                raise ValueError("Cluster labels cannot contain missing values.")
            observed_clusters.update(chunk[cluster_column].dropna().unique())
            values = chunk[gene_columns]
            non_numeric = [
                gene for gene in gene_columns
                if not pd.api.types.is_numeric_dtype(values[gene])
            ]
            if non_numeric:
                raise ValueError(
                    "Gene-expression columns must be numeric: "
                    + ", ".join(non_numeric)
                )
            dense_bytes += int(values.memory_usage(index=False).sum())
        genes = len(gene_columns)
        clusters = len(observed_clusters)
        source = "csv"
        storage = "dense"
    else:
        raise ValueError("Input must be a CSV or H5AD file.")

    if cells == 0:
        raise ValueError("The expression dataset contains no cells.")
    if clusters < 2:
        raise ValueError("At least two clusters are required.")

    dense_mib = dense_bytes / (1024 ** 2)
    limit_enabled = normalized_dense_limit not in {None, 0}
    native_within_limit = (
        suffix == ".csv"
        or not limit_enabled
        or dense_mib <= normalized_dense_limit
    )
    recommendation_reasons = []
    if suffix == ".h5ad" and storage == "sparse":
        recommendation_reasons.append(
            "Scanpy preserves the selected sparse H5AD matrix."
        )
    if suffix == ".h5ad" and not native_within_limit:
        recommendation_reasons.append(
            "The estimated dense matrix exceeds the Native safety limit."
        )
    recommended_engine = (
        "scanpy" if recommendation_reasons else "native"
    )

    return {
        "app_version": APP_VERSION,
        "input": str(input_path),
        "input_format": suffix.lstrip("."),
        "expression_source": source,
        "storage": storage,
        "cells": int(cells),
        "genes": int(genes),
        "clusters": clusters,
        "estimated_dense_bytes": dense_bytes,
        "estimated_dense_mib": dense_mib,
        "native_dense_limit_mib": normalized_dense_limit,
        "native_within_limit": native_within_limit,
        "recommended_engine": recommended_engine,
        "recommendation_reasons": recommendation_reasons,
    }


def build_analysis_plan(
    input_path: str | Path,
    requested_engine: str = "native",
    cluster_column: str = "cluster",
    cell_column: str | None = "cell",
    layer: str | None = None,
    use_raw: bool = False,
    max_dense_memory_mb: float | None = DEFAULT_MAX_DENSE_MEMORY_MB,
    gene_chunk_size: int = DEFAULT_GENE_CHUNK_SIZE,
    workers: int = 1,
    registry=None,
) -> dict[str, object]:
    """Build a non-destructive, machine-readable analysis execution plan."""
    from .engines import DEFAULT_ENGINE_REGISTRY, check_capabilities
    from .config import AnalysisConfig
    selected = (registry or DEFAULT_ENGINE_REGISTRY).get(requested_engine)

    _validate_gene_chunk_size(gene_chunk_size)
    _validate_workers(workers)
    inspection = inspect_input(
        input_path,
        cluster_column=cluster_column,
        cell_column=cell_column,
        layer=layer,
        use_raw=use_raw,
        max_dense_memory_mb=max_dense_memory_mb,
    )
    recommended_engine = str(inspection["recommended_engine"])
    resolved_engine = (
        recommended_engine if requested_engine == "auto" else requested_engine
    )
    warnings: list[str] = []
    can_execute = True
    try:
        check_capabilities(selected, Path(input_path).suffix.lower().lstrip("."), AnalysisConfig(
            engine=requested_engine, layer=layer, use_raw=use_raw, workers=workers))
    except ValueError as error:
        warnings.append(str(error))
        can_execute = False

    if requested_engine != "auto" and requested_engine != recommended_engine:
        warnings.append(
            f"Explicit {requested_engine} selection differs from the "
            f"{recommended_engine} resource recommendation."
        )
    if resolved_engine == "scanpy" and inspection["input_format"] != "h5ad":
        warnings.append("Scanpy requires H5AD input.")
        can_execute = False
    if resolved_engine == "native" and not inspection["native_within_limit"]:
        warnings.append(
            "Native dense conversion exceeds the configured memory limit."
        )
        can_execute = False
    if resolved_engine == "scanpy" and workers != 1:
        warnings.append("Scanpy requires workers=1.")
        can_execute = False

    return {
        "schema_version": PLAN_SCHEMA_VERSION,
        **inspection,
        "requested_engine": requested_engine,
        "resolved_engine": resolved_engine,
        "gene_chunk_size": (
            int(gene_chunk_size) if resolved_engine == "native" else None
        ),
        "workers": int(workers),
        "can_execute": can_execute,
        "warnings": warnings,
    }

__all__ = [
    "CSV_INSPECTION_CHUNK_ROWS",
    "PLAN_SCHEMA_VERSION",
    "build_analysis_plan",
    "estimate_dense_expression_bytes",
    "inspect_input",
]
