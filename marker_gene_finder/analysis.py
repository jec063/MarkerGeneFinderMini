"""Single-file analysis orchestration."""
from __future__ import annotations

from ._atomic import write_csv

from collections.abc import Sequence
from pathlib import Path
from time import perf_counter
import anndata as ad
import pandas as pd
from .version import APP_VERSION
from ._constants import (
    DEFAULT_GENE_CHUNK_SIZE, DEFAULT_MAX_DENSE_MEMORY_MB, RANK_BY_CHOICES, RESULT_COLUMNS,
)
from ._validation import _validate_gene_chunk_size, _validate_workers
from .input_loading import load_expression
from .inspection import build_analysis_plan
from .progress import ProgressCallback
from .statistics import _resolve_cluster_arguments, find_markers, find_markers_scanpy


def run(
    input_path: str | Path,
    output_path: str | Path,
    cluster_column: str = "cluster",
    cell_column: str | None = "cell",
    top_n: int = 5,
    pseudocount: float = 0.1,
    min_pct: float = 0.0,
    min_pct_difference: float = 0.0,
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
    engine: str = "native",
    layer: str | None = None,
    use_raw: bool = False,
    target_clusters: Sequence[str] | None = None,
    excluded_genes: Sequence[str] | None = None,
    excluded_prefixes: Sequence[str] | None = None,
    max_pct_out: float = 1.0,
    rank_by: str = "log2FC",
    min_cluster_size: int = 1,
    min_cells_expressing: int = 0,
    min_cluster_mean: float = 0.0,
    min_mean_difference: float = 0.0,
    max_other_mean: float | None = None,
    max_cells_expressing_out: int | None = None,
    max_dense_memory_mb: float | None = DEFAULT_MAX_DENSE_MEMORY_MB,
    gene_chunk_size: int = DEFAULT_GENE_CHUNK_SIZE,
    progress_callback: ProgressCallback | None = None,
    workers: int = 1,
    diagnostics: dict[str, object] | None = None,
) -> pd.DataFrame:
    """Read an expression dataset, find markers and save the results."""

    run_started = perf_counter()
    _validate_gene_chunk_size(gene_chunk_size)
    _validate_workers(workers)

    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")
    if use_raw and Path(input_path).suffix.lower() != ".h5ad":
        raise ValueError("Raw selection requires H5AD input.")

    if engine not in {"native", "scanpy", "auto"}:
        raise ValueError(
            "engine must be 'native', 'scanpy', or 'auto'."
        )
    requested_engine = engine
    if engine == "auto":
        plan = build_analysis_plan(
            input_path,
            requested_engine="auto",
            cluster_column=cluster_column,
            cell_column=cell_column,
            layer=layer,
            use_raw=use_raw,
            max_dense_memory_mb=max_dense_memory_mb,
            gene_chunk_size=gene_chunk_size,
            workers=workers,
        )
        if not plan["can_execute"]:
            raise ValueError("; ".join(plan["warnings"]))
        engine = str(plan["resolved_engine"])
    if engine == "scanpy" and workers != 1:
        raise ValueError(
            "workers greater than 1 are supported by the Native engine only."
        )

    input_path = Path(input_path)
    load_started = perf_counter()

    if engine == "scanpy":
        if input_path.suffix.lower() != ".h5ad":
            raise ValueError(
                "Scanpy engine requires a .h5ad input file."
            )

        adata = ad.read_h5ad(input_path)
        cell_count = int(adata.n_obs)
        if use_raw:
            if adata.raw is None or adata.raw.X is None:
                raise ValueError("The AnnData object does not contain raw expression.")
            gene_count = int(adata.raw.n_vars)
        else:
            gene_count = int(adata.n_vars)
        cluster_count = int(adata.obs[cluster_column].nunique())
        load_seconds = perf_counter() - load_started
        analysis_started = perf_counter()
        selected_clusters = _resolve_cluster_arguments(
            adata.obs[cluster_column], target_clusters
        )

        markers = find_markers_scanpy(
            adata,
            layer=layer,
            use_raw=use_raw,
            cluster_column=cluster_column,
            top_n=top_n,
            min_pct=min_pct,
            max_pct_out=max_pct_out,
            rank_by=rank_by,
            min_pct_difference=min_pct_difference,
            min_log2fc=min_log2fc,
            max_p_adj=max_p_adj,
            target_clusters=selected_clusters,
            excluded_genes=excluded_genes,
            excluded_prefixes=excluded_prefixes,
            min_cluster_size=min_cluster_size,
            min_cells_expressing=min_cells_expressing,
            min_cluster_mean=min_cluster_mean,
            min_mean_difference=min_mean_difference,
            max_other_mean=max_other_mean,
            max_cells_expressing_out=max_cells_expressing_out,
            progress_callback=progress_callback,
        )
        analysis_seconds = perf_counter() - analysis_started
    else:
        expression = load_expression(
            input_path,
            layer=layer,
            use_raw=use_raw,
            cluster_column=cluster_column,
            cell_column=cell_column,
            max_dense_memory_mb=max_dense_memory_mb,
        )
        cell_count = len(expression)
        gene_count = len(expression.columns) - 1 - int(
            cell_column is not None and cell_column in expression.columns
        )
        cluster_count = int(expression[cluster_column].nunique())
        load_seconds = perf_counter() - load_started
        analysis_started = perf_counter()
        selected_clusters = _resolve_cluster_arguments(
            expression[cluster_column], target_clusters
        )

        markers = find_markers(
            expression,
            cluster_column=cluster_column,
            cell_column=cell_column,
            top_n=top_n,
            pseudocount=pseudocount,
            min_pct=min_pct,
            max_pct_out=max_pct_out,
            rank_by=rank_by,
            min_pct_difference=min_pct_difference,
            min_log2fc=min_log2fc,
            max_p_adj=max_p_adj,
            target_clusters=selected_clusters,
            excluded_genes=excluded_genes,
            excluded_prefixes=excluded_prefixes,
            min_cluster_size=min_cluster_size,
            min_cells_expressing=min_cells_expressing,
            min_cluster_mean=min_cluster_mean,
            min_mean_difference=min_mean_difference,
            max_other_mean=max_other_mean,
            max_cells_expressing_out=max_cells_expressing_out,
            gene_chunk_size=gene_chunk_size,
            progress_callback=progress_callback,
            workers=workers,
        )
        analysis_seconds = perf_counter() - analysis_started

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_started = perf_counter()
    write_csv(output_path, markers)
    output_seconds = perf_counter() - output_started

    if diagnostics is not None:
        diagnostics.update(
            {
                "app_version": APP_VERSION,
                "engine": engine,
                "requested_engine": requested_engine,
                "resolved_engine": engine,
                "input_format": input_path.suffix.lower().lstrip("."),
                "input": str(input_path),
                "output": str(output_path),
                "cells": cell_count,
                "genes": gene_count,
                "clusters": cluster_count,
                "marker_rows": len(markers),
                "gene_chunk_size": (
                    gene_chunk_size if engine == "native" else None
                ),
                "workers": workers,
                "load_seconds": load_seconds,
                "analysis_seconds": analysis_seconds,
                "output_seconds": output_seconds,
                "total_seconds": perf_counter() - run_started,
            }
        )

    return markers


def _registered_run(input_path, output_path, *, registry=None, diagnostics, **options):
    if registry is None and options.get("engine", "native") in {"native", "scanpy", "auto"}:
        return run(input_path, output_path, diagnostics=diagnostics, **options)
    from .api import analyze_file
    from .config import AnalysisConfig
    from .engines import DEFAULT_ENGINE_REGISTRY
    values = dict(options)
    for key in ("target_clusters", "excluded_genes", "excluded_prefixes"):
        if key in values:
            values[key] = tuple(values[key] or ())
    started = perf_counter()
    result = analyze_file(input_path, output_path, AnalysisConfig(**values),
                          registry=registry or DEFAULT_ENGINE_REGISTRY)
    result.to_csv(output_path)
    diagnostics.update(result.metadata)
    diagnostics["total_seconds"] = perf_counter() - started
    return result.markers


__all__ = [
    "DEFAULT_GENE_CHUNK_SIZE",
    "DEFAULT_MAX_DENSE_MEMORY_MB",
    "RANK_BY_CHOICES",
    "RESULT_COLUMNS",
    "find_markers",
    "find_markers_scanpy",
    "load_expression",
    "run",
]
