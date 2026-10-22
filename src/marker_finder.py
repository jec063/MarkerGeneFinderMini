"""Find simple marker genes by comparing one cluster with all other cells."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from time import perf_counter

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

if __package__:
    from .version import APP_VERSION
else:
    from version import APP_VERSION

RESULT_COLUMNS = [
    "cluster",
    "gene",
    "rank",
    "cluster_size",
    "reference_size",
    "cluster_mean",
    "other_mean",
    "mean_difference",
    "pct_in",
    "pct_out",
    "cells_expressing_in",
    "cells_expressing_out",
    "pct_difference",
    "log2FC",
    "p_value",
    "p_adj",
]
RANK_BY_CHOICES = (
    "log2FC",
    "cluster_mean",
    "cells_expressing_in",
    "cells_expressing_out",
    "mean_difference",
    "pct_difference",
    "p_adj",
)
DEFAULT_MAX_DENSE_MEMORY_MB = 2048.0
DEFAULT_GENE_CHUNK_SIZE = 1000
ProgressCallback = Callable[[int, int, object], None]


def _validate_workers(workers: int) -> None:
    """Require a positive, non-Boolean worker count."""
    if (
        isinstance(workers, bool)
        or not isinstance(workers, (int, np.integer))
        or workers < 1
    ):
        raise ValueError("workers must be a positive integer.")


def _validate_gene_chunk_size(gene_chunk_size: int) -> None:
    """Require a positive, non-Boolean Native gene chunk size."""
    if (
        isinstance(gene_chunk_size, bool)
        or not isinstance(gene_chunk_size, (int, np.integer))
        or gene_chunk_size < 1
    ):
        raise ValueError("gene_chunk_size must be a positive integer.")


def _rank_markers(markers: pd.DataFrame, rank_by: str, top_n: int) -> pd.DataFrame:
    """Sort filtered markers and assign stable one-based ranks."""
    if rank_by not in RANK_BY_CHOICES:
        choices = ", ".join(RANK_BY_CHOICES)
        raise ValueError(f"rank_by must be one of: {choices}.")

    if rank_by == "p_adj":
        columns = ["p_adj", "log2FC", "gene"]
        ascending = [True, False, True]
    elif rank_by == "cells_expressing_out":
        columns = ["cells_expressing_out", "log2FC", "gene"]
        ascending = [True, False, True]
    else:
        columns = [rank_by, "log2FC", "gene"]
        ascending = [False, False, True]

    ranked = markers.sort_values(columns, ascending=ascending).head(top_n).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    return ranked



def _raw_expression(adata: ad.AnnData) -> ad.AnnData:
    """Copy raw expression with its gene names and current cell metadata."""
    from copy import deepcopy

    if adata.raw is None or adata.raw.X is None:
        raise ValueError("The AnnData object does not contain raw expression.")
    return ad.AnnData(
        X=adata.raw.X.copy(),
        obs=adata.obs.copy(),
        var=adata.raw.var.copy(),
        uns=deepcopy(adata.uns),
    )


def estimate_dense_expression_bytes(
    input_path: str | Path,
    layer: str | None = None,
    use_raw: bool = False,
) -> int:
    """Estimate bytes needed to materialize the selected H5AD matrix."""
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".h5ad":
        raise ValueError("Dense-memory estimation requires H5AD input.")
    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")

    adata = ad.read_h5ad(input_path, backed="r")
    try:
        if use_raw:
            if adata.raw is None or adata.raw.X is None:
                raise ValueError("The AnnData object does not contain raw expression.")
            matrix = adata.raw.X
        elif layer is not None:
            if layer not in adata.layers:
                raise ValueError(f"Missing AnnData layer: {layer!r}.")
            matrix = adata.layers[layer]
        else:
            if adata.X is None:
                raise ValueError(
                    "The AnnData object does not contain an expression matrix in X."
                )
            matrix = adata.X
        rows, columns = matrix.shape
        return int(rows) * int(columns) * np.dtype(matrix.dtype).itemsize
    finally:
        adata.file.close()


def _enforce_dense_memory_limit(
    input_path: str | Path,
    layer: str | None,
    use_raw: bool,
    max_dense_memory_mb: float | None,
) -> None:
    """Reject Native H5AD conversion when its dense matrix exceeds a limit."""
    if max_dense_memory_mb is None:
        return
    if (
        isinstance(max_dense_memory_mb, bool)
        or not isinstance(
            max_dense_memory_mb,
            (int, float, np.integer, np.floating),
        )
        or not np.isfinite(max_dense_memory_mb)
        or max_dense_memory_mb < 0
    ):
        raise ValueError("max_dense_memory_mb must be non-negative or None.")
    if max_dense_memory_mb == 0:
        return

    estimated_bytes = estimate_dense_expression_bytes(
        input_path,
        layer=layer,
        use_raw=use_raw,
    )
    estimated_mb = estimated_bytes / (1024 ** 2)
    if estimated_mb > max_dense_memory_mb:
        raise ValueError(
            "Native dense conversion requires approximately "
            f"{estimated_mb:.1f} MiB, exceeding the "
            f"{max_dense_memory_mb:g} MiB safety limit. Use the Scanpy "
            "engine, raise max_dense_memory_mb, or set it to 0 to disable "
            "the limit."
        )


def load_expression(
    input_path: str | Path,
    cluster_column: str = "cluster",
    cell_column: str | None = "cell_id",
    layer: str | None = None,
    use_raw: bool = False,
    max_dense_memory_mb: float | None = DEFAULT_MAX_DENSE_MEMORY_MB,
) -> pd.DataFrame:
    """Load CSV or AnnData input as an expression table."""

    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")
    if use_raw and Path(input_path).suffix.lower() != ".h5ad":
        raise ValueError("Raw selection requires H5AD input.")
    input_path = Path(input_path)
    suffix = input_path.suffix.lower()

    if suffix == ".csv":
        if layer is not None:
            raise ValueError("Layer selection requires H5AD input.")
        return pd.read_csv(input_path)

    if suffix != ".h5ad":
        raise ValueError(
            "Unsupported input format. Use a .csv or .h5ad file."
        )

    _enforce_dense_memory_limit(
        input_path,
        layer,
        use_raw,
        max_dense_memory_mb,
    )

    adata = ad.read_h5ad(input_path)
    if use_raw:
        adata = _raw_expression(adata)

    if cluster_column not in adata.obs.columns:
        raise ValueError(
            f"Missing cluster column in AnnData.obs: "
            f"{cluster_column!r}."
        )

    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Missing AnnData layer: {layer!r}.")

    if layer is None and adata.X is None:
        raise ValueError(
            "The AnnData object does not contain an "
            "expression matrix in X."
        )

    if not adata.var_names.is_unique:
        raise ValueError(
            "AnnData gene names in var_names must be unique."
        )

    matrix = adata.X if layer is None else adata.layers[layer]
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()

    expression = pd.DataFrame(
        np.asarray(matrix),
        columns=adata.var_names.astype(str),
    )
    expression.insert(
        0,
        cluster_column,
        adata.obs[cluster_column].to_numpy(),
    )

    if cell_column:
        if cell_column == cluster_column:
            raise ValueError(
                "cell_column and cluster_column must be different."
            )
        expression.insert(
            0,
            cell_column,
            adata.obs_names.astype(str),
        )

    return expression


def _benjamini_hochberg(
    p_values: np.ndarray,
) -> np.ndarray:
    """Adjust p-values using Benjamini-Hochberg FDR correction."""

    p_values = np.asarray(p_values, dtype=float)

    if p_values.size == 0:
        return p_values.copy()

    order = np.argsort(p_values, kind="stable")
    sorted_p_values = p_values[order]
    ranks = np.arange(1, p_values.size + 1)

    adjusted_sorted = (
        sorted_p_values * p_values.size / ranks
    )
    adjusted_sorted = np.minimum.accumulate(
        adjusted_sorted[::-1]
    )[::-1]
    adjusted_sorted = np.clip(adjusted_sorted, 0.0, 1.0)

    adjusted = np.empty_like(adjusted_sorted)
    adjusted[order] = adjusted_sorted

    return adjusted


def _selected_clusters(labels, target_clusters, min_cluster_size=1):
    """Validate and return clusters selected for one-versus-rest analysis."""
    if (
        isinstance(min_cluster_size, bool)
        or not isinstance(min_cluster_size, (int, np.integer))
        or min_cluster_size < 1
    ):
        raise ValueError("min_cluster_size must be an integer of at least 1.")

    available = sorted(pd.unique(labels), key=str)
    cluster_sizes = pd.Series(labels).value_counts(dropna=False)
    eligible = [
        cluster for cluster in available
        if int(cluster_sizes.loc[cluster]) >= min_cluster_size
    ]
    if target_clusters is None:
        if not eligible:
            raise ValueError(
                "No clusters meet the minimum cluster size of "
                f"{min_cluster_size}."
            )
        return eligible
    selected = list(dict.fromkeys(target_clusters))
    if not selected:
        raise ValueError("At least one target cluster must be selected.")
    missing = [cluster for cluster in selected if cluster not in available]
    if missing:
        raise ValueError(f"Unknown target clusters: {missing}.")
    undersized = [cluster for cluster in selected if cluster not in eligible]
    if undersized:
        raise ValueError(
            "Target clusters below the minimum cluster size of "
            f"{min_cluster_size}: {undersized}."
        )
    return selected


def _resolve_cluster_arguments(labels, requested):
    """Resolve string CLI cluster values to their original label values."""
    if requested is None:
        return None
    by_text = {}
    for label in pd.unique(labels):
        by_text.setdefault(str(label), []).append(label)
    missing = [value for value in requested if value not in by_text]
    ambiguous = [value for value in requested if len(by_text.get(value, [])) > 1]
    if missing or ambiguous:
        detail = missing or ambiguous
        raise ValueError(f"Unknown or ambiguous target clusters: {detail}.")
    return [by_text[value][0] for value in requested]


def _filter_genes(genes, excluded_genes=None, excluded_prefixes=None):
    """Remove explicitly excluded genes and gene-name prefixes."""
    excluded = set(excluded_genes or [])
    prefixes = tuple(prefix for prefix in (excluded_prefixes or []) if prefix)
    selected = [
        gene for gene in genes
        if str(gene) not in excluded and not str(gene).startswith(prefixes)
    ]
    if not selected:
        raise ValueError("Gene exclusions removed every gene.")
    return selected


def _gene_columns(
    expression: pd.DataFrame,
    cluster_column: str,
    cell_column: str | None,
) -> list[str]:
    """Validate the input table and return the gene-expression columns."""

    if expression.empty:
        raise ValueError("The expression table is empty.")

    if cluster_column not in expression.columns:
        raise ValueError(f"Missing cluster column: {cluster_column!r}.")

    if expression[cluster_column].isna().any():
        raise ValueError("Cluster labels cannot contain missing values.")

    if expression[cluster_column].nunique() < 2:
        raise ValueError("At least two clusters are required.")

    metadata_columns = {cluster_column}

    if cell_column and cell_column in expression.columns:
        metadata_columns.add(cell_column)

    genes = [
        column
        for column in expression.columns
        if column not in metadata_columns
    ]

    if not genes:
        raise ValueError("No gene-expression columns were found.")

    non_numeric = [
        gene
        for gene in genes
        if not pd.api.types.is_numeric_dtype(expression[gene])
    ]

    if non_numeric:
        raise ValueError(
            "Gene-expression columns must be numeric: "
            + ", ".join(non_numeric)
        )

    for start in range(0, len(genes), DEFAULT_GENE_CHUNK_SIZE):
        values = expression[genes[start:start + DEFAULT_GENE_CHUNK_SIZE]]
        if values.isna().any().any():
            raise ValueError(
                "Gene-expression values cannot contain missing values."
            )
        if (values < 0).any().any():
            raise ValueError(
                "Gene-expression values cannot be negative."
            )

    return genes


def find_markers(
    expression: pd.DataFrame,
    cluster_column: str = "cluster",
    cell_column: str | None = "cell",
    top_n: int = 5,
    pseudocount: float = 0.1,
    min_pct: float = 0.0,
    min_pct_difference: float = 0.0,
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
    target_clusters: Sequence[object] | None = None,
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
    gene_chunk_size: int = DEFAULT_GENE_CHUNK_SIZE,
    progress_callback: ProgressCallback | None = None,
    workers: int = 1,
) -> pd.DataFrame:
    """Return top markers with chunked, optionally parallel Native analysis."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1.")

    _validate_gene_chunk_size(gene_chunk_size)
    _validate_workers(workers)

    if (
        isinstance(min_cells_expressing, bool)
        or not isinstance(min_cells_expressing, (int, np.integer))
        or min_cells_expressing < 0
    ):
        raise ValueError(
            "min_cells_expressing must be a non-negative integer."
        )

    if pseudocount <= 0:
        raise ValueError("pseudocount must be greater than 0.")

    if min_cluster_mean < 0:
        raise ValueError("min_cluster_mean must be at least 0.")

    if min_mean_difference < 0:
        raise ValueError("min_mean_difference must be at least 0.")

    if max_other_mean is not None and max_other_mean < 0:
        raise ValueError("max_other_mean must be at least 0 when provided.")

    if max_cells_expressing_out is not None and (
        isinstance(max_cells_expressing_out, bool)
        or not isinstance(max_cells_expressing_out, (int, np.integer))
        or max_cells_expressing_out < 0
    ):
        raise ValueError(
            "max_cells_expressing_out must be a non-negative integer "
            "when provided."
        )

    if not 0 <= min_pct <= 1:
        raise ValueError("min_pct must be between 0 and 1.")

    if not 0 <= max_pct_out <= 1:
        raise ValueError("max_pct_out must be between 0 and 1.")

    if not 0 <= min_pct_difference <= 1:
        raise ValueError("min_pct_difference must be between 0 and 1.")

    if min_log2fc < 0:
        raise ValueError("min_log2fc must be at least 0.")

    if not 0 <= max_p_adj <= 1:
        raise ValueError("max_p_adj must be between 0 and 1.")

    if rank_by not in RANK_BY_CHOICES:
        choices = ", ".join(RANK_BY_CHOICES)
        raise ValueError(f"rank_by must be one of: {choices}.")

    genes = _gene_columns(
        expression,
        cluster_column,
        cell_column,
    )
    genes = _filter_genes(genes, excluded_genes, excluded_prefixes)

    clusters = _selected_clusters(
        expression[cluster_column], target_clusters, min_cluster_size,
    )

    def analyze_cluster(cluster: object) -> pd.DataFrame:
        in_cluster = expression[cluster_column] == cluster
        cluster_size = int(in_cluster.sum())
        reference_size = int((~in_cluster).sum())

        chunk_results = []
        for start in range(0, len(genes), gene_chunk_size):
            chunk = genes[start:start + gene_chunk_size]
            cluster_values = expression.loc[in_cluster, chunk].to_numpy()
            other_values = expression.loc[~in_cluster, chunk].to_numpy()
            cluster_mean = cluster_values.mean(axis=0)
            other_mean = other_values.mean(axis=0)
            cluster_expressed = cluster_values > 0
            other_expressed = other_values > 0
            cells_expressing_in = cluster_expressed.sum(axis=0)
            cells_expressing_out = other_expressed.sum(axis=0)
            pct_in = cells_expressing_in / cluster_size
            pct_out = cells_expressing_out / reference_size
            p_values = np.asarray(
                mannwhitneyu(
                    cluster_values,
                    other_values,
                    alternative="greater",
                    axis=0,
                    method="asymptotic",
                ).pvalue,
                dtype=float,
            )
            chunk_results.append(
                pd.DataFrame(
                    {
                        "cluster": cluster,
                        "gene": chunk,
                        "cluster_size": cluster_size,
                        "reference_size": reference_size,
                        "cluster_mean": cluster_mean,
                        "other_mean": other_mean,
                        "mean_difference": cluster_mean - other_mean,
                        "pct_in": pct_in,
                        "pct_out": pct_out,
                        "cells_expressing_in": cells_expressing_in,
                        "cells_expressing_out": cells_expressing_out,
                        "pct_difference": pct_in - pct_out,
                        "log2FC": np.log2(
                            (cluster_mean + pseudocount)
                            / (other_mean + pseudocount)
                        ),
                        "p_value": p_values,
                    }
                )
            )

        cluster_result = pd.concat(chunk_results, ignore_index=True)
        # Adjust across every tested gene, not independently within chunks.
        cluster_result["p_adj"] = _benjamini_hochberg(
            cluster_result["p_value"].to_numpy()
        )
        other_mean_passes = (
            True if max_other_mean is None
            else cluster_result["other_mean"] <= max_other_mean
        )
        return _rank_markers(
            cluster_result.loc[
                (cluster_result["pct_in"] >= min_pct)
                & (cluster_result["cells_expressing_in"] >= min_cells_expressing)
                & (cluster_result["cluster_mean"] >= min_cluster_mean)
                & (cluster_result["mean_difference"] >= min_mean_difference)
                & other_mean_passes
                & (
                    True if max_cells_expressing_out is None
                    else cluster_result["cells_expressing_out"]
                    <= max_cells_expressing_out
                )
                & (cluster_result["pct_out"] <= max_pct_out)
                & (cluster_result["pct_difference"] >= min_pct_difference)
                & (cluster_result["log2FC"] >= min_log2fc)
                & (cluster_result["p_adj"] <= max_p_adj)
            ],
            rank_by,
            top_n,
        )

    if workers == 1:
        analyzed = map(analyze_cluster, clusters)
        executor = None
    else:
        executor = ThreadPoolExecutor(max_workers=workers)
        analyzed = executor.map(analyze_cluster, clusters)

    try:
        results = []
        for cluster_index, (cluster, cluster_result) in enumerate(
            zip(clusters, analyzed), start=1
        ):
            results.append(cluster_result)
            if progress_callback is not None:
                progress_callback(cluster_index, len(clusters), cluster)
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    if not results:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    return pd.concat(
        results,
        ignore_index=True,
    )[RESULT_COLUMNS]

def find_markers_scanpy(
    adata: ad.AnnData,
    cluster_column: str = "cluster",
    top_n: int = 5,
    min_pct: float = 0.0,
    min_pct_difference: float = 0.0,
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
    layer: str | None = None,
    use_raw: bool = False,
    target_clusters: Sequence[object] | None = None,
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
    progress_callback: ProgressCallback | None = None,
) -> pd.DataFrame:
    """Find marker genes using Scanpy's Wilcoxon implementation."""

    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")

    if use_raw:
        adata = _raw_expression(adata)
    if top_n < 1:
        raise ValueError("top_n must be at least 1.")

    if (
        isinstance(min_cells_expressing, bool)
        or not isinstance(min_cells_expressing, (int, np.integer))
        or min_cells_expressing < 0
    ):
        raise ValueError(
            "min_cells_expressing must be a non-negative integer."
        )

    if not 0 <= min_pct <= 1:
        raise ValueError("min_pct must be between 0 and 1.")

    if min_cluster_mean < 0:
        raise ValueError("min_cluster_mean must be at least 0.")

    if min_mean_difference < 0:
        raise ValueError("min_mean_difference must be at least 0.")

    if max_other_mean is not None and max_other_mean < 0:
        raise ValueError("max_other_mean must be at least 0 when provided.")

    if max_cells_expressing_out is not None and (
        isinstance(max_cells_expressing_out, bool)
        or not isinstance(max_cells_expressing_out, (int, np.integer))
        or max_cells_expressing_out < 0
    ):
        raise ValueError(
            "max_cells_expressing_out must be a non-negative integer "
            "when provided."
        )

    if not 0 <= max_pct_out <= 1:
        raise ValueError("max_pct_out must be between 0 and 1.")

    if not 0 <= min_pct_difference <= 1:
        raise ValueError("min_pct_difference must be between 0 and 1.")

    if min_log2fc < 0:
        raise ValueError("min_log2fc must be at least 0.")

    if not 0 <= max_p_adj <= 1:
        raise ValueError("max_p_adj must be between 0 and 1.")

    if rank_by not in RANK_BY_CHOICES:
        choices = ", ".join(RANK_BY_CHOICES)
        raise ValueError(f"rank_by must be one of: {choices}.")

    if cluster_column not in adata.obs.columns:
        raise ValueError(
            f"Missing cluster column in AnnData.obs: "
            f"{cluster_column!r}."
        )

    if adata.obs[cluster_column].isna().any():
        raise ValueError("Cluster labels cannot contain missing values.")

    if adata.obs[cluster_column].nunique() < 2:
        raise ValueError("At least two clusters are required.")

    if layer is not None and layer not in adata.layers:
        raise ValueError(f"Missing AnnData layer: {layer!r}.")

    if layer is None and adata.X is None:
        raise ValueError(
            "The AnnData object does not contain an "
            "expression matrix in X."
        )

    if adata.n_vars == 0:
        raise ValueError("No genes were found in the AnnData object.")

    if not adata.var_names.is_unique:
        raise ValueError("AnnData gene names must be unique.")

    working = adata.copy()
    if layer is not None:
        working.X = working.layers[layer].copy()
    selected_genes = _filter_genes(
        working.var_names.astype(str).tolist(), excluded_genes, excluded_prefixes
    )
    working = working[:, selected_genes].copy()
    working.obs[cluster_column] = (
        working.obs[cluster_column].astype("category")
    )

    clusters = _selected_clusters(
        working.obs[cluster_column], target_clusters, min_cluster_size
    )

    sc.tl.rank_genes_groups(
        working,
        groupby=cluster_column,
        groups=clusters,
        reference="rest",
        method="wilcoxon",
        corr_method="benjamini-hochberg",
        use_raw=False,
        n_genes=working.n_vars,
        pts=True,
        tie_correct=True,
    )

    ranked = sc.get.rank_genes_groups_df(
        working,
        group=None,
    ).rename(
        columns={
            "group": "cluster",
            "names": "gene",
            "logfoldchanges": "log2FC",
            "pvals": "p_value",
            "pvals_adj": "p_adj",
            "pct_nz_group": "pct_in",
            "pct_nz_reference": "pct_out",
        }
    )
    if "cluster" not in ranked.columns:
        ranked["cluster"] = clusters[0]

    gene_positions = pd.Series(
        np.arange(working.n_vars),
        index=working.var_names.astype(str),
    )

    results = []

    for cluster_index, cluster in enumerate(clusters, start=1):
        in_cluster = np.asarray(
            working.obs[cluster_column] == cluster,
            dtype=bool,
        )
        cluster_size = int(in_cluster.sum())
        reference_size = int((~in_cluster).sum())

        cluster_mean = np.asarray(
            working.X[in_cluster].mean(axis=0)
        ).ravel()
        other_mean = np.asarray(
            working.X[~in_cluster].mean(axis=0)
        ).ravel()
        cells_expressing_in = np.asarray(
            (working.X[in_cluster] > 0).sum(axis=0)
        ).ravel().astype(int)
        cells_expressing_out = np.asarray(
            (working.X[~in_cluster] > 0).sum(axis=0)
        ).ravel().astype(int)

        cluster_result = ranked.loc[
            ranked["cluster"].astype(str) == str(cluster)
        ].copy()

        positions = gene_positions.loc[
            cluster_result["gene"].astype(str)
        ].to_numpy()

        cluster_result["cluster"] = cluster
        cluster_result["cluster_size"] = cluster_size
        cluster_result["reference_size"] = reference_size
        cluster_result["cluster_mean"] = cluster_mean[positions]
        cluster_result["other_mean"] = other_mean[positions]
        cluster_result["mean_difference"] = (
            cluster_result["cluster_mean"] - cluster_result["other_mean"]
        )
        cluster_result["cells_expressing_in"] = cells_expressing_in[positions]
        cluster_result["cells_expressing_out"] = cells_expressing_out[positions]
        cluster_result["pct_difference"] = (
            cluster_result["pct_in"] - cluster_result["pct_out"]
        )

        other_mean_passes = (
            True if max_other_mean is None
            else cluster_result["other_mean"] <= max_other_mean
        )

        cluster_result = _rank_markers(
            cluster_result.loc[
                (cluster_result["pct_in"] >= min_pct)
                & (cluster_result["cells_expressing_in"] >= min_cells_expressing)
                & (cluster_result["cluster_mean"] >= min_cluster_mean)
                & (cluster_result["mean_difference"] >= min_mean_difference)
                & other_mean_passes
                & (
                    True if max_cells_expressing_out is None
                    else cluster_result["cells_expressing_out"]
                    <= max_cells_expressing_out
                )
                & (cluster_result["pct_out"] <= max_pct_out)
                & (cluster_result["pct_difference"] >= min_pct_difference)
                & (cluster_result["log2FC"] >= min_log2fc)
                & (cluster_result["p_adj"] <= max_p_adj)
            ],
            rank_by,
            top_n,
        )

        results.append(cluster_result)
        if progress_callback is not None:
            progress_callback(cluster_index, len(clusters), cluster)

    if not results:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    return pd.concat(
        results,
        ignore_index=True,
    )[RESULT_COLUMNS]


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

    if engine not in {"native", "scanpy"}:
        raise ValueError(
            "engine must be either 'native' or 'scanpy'."
        )
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
    markers.to_csv(output_path, index=False)
    output_seconds = perf_counter() - output_started

    if diagnostics is not None:
        diagnostics.update(
            {
                "app_version": APP_VERSION,
                "engine": engine,
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


def run_batch(
    input_dir: str | Path,
    output_dir: str | Path,
    **analysis_options,
) -> list[dict[str, object]]:
    """Analyze every supported file and write a machine-readable run manifest."""
    batch_started = perf_counter()
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if not input_dir.is_dir():
        raise ValueError(f"Batch input directory does not exist: {input_dir}")

    inputs = sorted(
        path for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".csv", ".h5ad"}
    )
    if not inputs:
        raise ValueError("Batch input directory contains no CSV or H5AD files.")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for input_path in inputs:
        output_path = output_dir / f"{input_path.name}.markers.csv"
        record: dict[str, object] = {
            "input": str(input_path),
            "output": str(output_path),
        }
        run_diagnostics: dict[str, object] = {}
        attempt_started = perf_counter()
        try:
            markers = run(
                input_path,
                output_path,
                diagnostics=run_diagnostics,
                **analysis_options,
            )
        except Exception as error:
            record.update(status="failed", marker_rows=0, error=str(error))
        else:
            record.update(
                status="complete",
                marker_rows=len(markers),
                error=None,
                diagnostics=run_diagnostics,
            )
        record["duration_seconds"] = perf_counter() - attempt_started
        manifest.append(record)

    manifest_path = output_dir / "batch_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "app_version": APP_VERSION,
                "inputs": len(inputs),
                "completed": sum(row["status"] == "complete" for row in manifest),
                "failed": sum(row["status"] == "failed" for row in manifest),
                "duration_seconds": perf_counter() - batch_started,
                "runs": manifest,
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Rank marker genes for each cluster "
            "using log2 fold change."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {APP_VERSION}",
    )

    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--input",
        help="Input expression CSV or H5AD path.",
    )
    inputs.add_argument(
        "--input-dir",
        help="Directory of CSV or H5AD files to analyze as a batch.",
    )

    outputs = parser.add_mutually_exclusive_group(required=True)
    outputs.add_argument(
        "--output",
        help="Output marker CSV path.",
    )
    outputs.add_argument(
        "--output-dir",
        help="Directory for batch marker CSVs and batch_manifest.json.",
    )

    parser.add_argument(
        "--engine",
        choices=["native", "scanpy"],
        default="native",
        help=(
            "Marker-analysis engine. The Scanpy engine "
            "requires H5AD input."
        ),
    )

    parser.add_argument(
        "--cluster-column",
        default="cluster",
        help="Cluster-label column name.",
    )

    parser.add_argument(
        "--cell-column",
        default="cell",
        help="Cell-ID column excluded from gene columns.",
    )
    parser.add_argument(
        "--max-dense-memory-mb",
        type=float,
        default=DEFAULT_MAX_DENSE_MEMORY_MB,
        help=(
            "Native H5AD dense-conversion safety limit in MiB; "
            "use 0 to disable it."
        ),
    )
    parser.add_argument(
        "--gene-chunk-size",
        type=int,
        default=DEFAULT_GENE_CHUNK_SIZE,
        help=(
            "Number of genes processed together by the Native engine "
            f"(default: {DEFAULT_GENE_CHUNK_SIZE})."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Native target clusters analyzed concurrently (default: 1). "
            "Higher values can increase peak memory use."
        ),
    )
    parser.add_argument(
        "--diagnostics-json",
        help=(
            "Write single-file runtime diagnostics as JSON. Batch diagnostics "
            "are included in batch_manifest.json."
        ),
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-cluster progress messages.",
    )

    parser.add_argument(
        "--target-cluster",
        action="append",
        dest="target_clusters",
        help="Cluster to analyze; repeat to select multiple clusters.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=1,
        help="Minimum cell count required for each target cluster.",
    )
    parser.add_argument(
        "--exclude-gene", action="append", default=[],
        help="Exact gene name to exclude; repeat as needed.",
    )
    parser.add_argument(
        "--exclude-prefix", action="append", default=[],
        help="Gene-name prefix to exclude; repeat as needed.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of markers returned per cluster.",
    )

    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.1,
        help="Positive pseudocount used for log2FC.",
    )
    parser.add_argument(
        "--min-pct",
        type=float,
        default=0.0,
        help=(
            "Minimum fraction of cluster cells "
            "expressing a gene."
        ),
    )
    parser.add_argument(
        "--min-cells-expressing",
        type=int,
        default=0,
        help="Minimum target-cluster cells expressing a marker gene.",
    )
    parser.add_argument(
        "--min-cluster-mean",
        type=float,
        default=0.0,
        help="Minimum mean expression required inside the target cluster.",
    )
    parser.add_argument(
        "--min-mean-difference",
        type=float,
        default=0.0,
        help="Minimum target-minus-reference mean-expression difference.",
    )
    parser.add_argument(
        "--max-other-mean",
        type=float,
        default=None,
        help=(
            "Maximum mean expression allowed outside the target cluster; "
            "omit to disable this filter."
        ),
    )
    parser.add_argument(
        "--max-cells-expressing-out",
        type=int,
        default=None,
        help=(
            "Maximum number of reference cells expressing a marker; "
            "omit to disable this filter."
        ),
    )

    parser.add_argument(
        "--min-log2fc",
        type=float,
        default=0.0,
        help=(
            "Minimum log2 fold change required "
            "for a marker gene."
        ),
    )

    parser.add_argument(
        "--max-pct-out",
        type=float,
        default=1.0,
        help="Maximum fraction of non-cluster cells expressing a marker.",
    )

    parser.add_argument(
        "--rank-by",
        choices=RANK_BY_CHOICES,
        default="log2FC",
        help="Column used to rank filtered markers within each cluster.",
    )

    parser.add_argument(
        "--min-pct-difference",
        type=float,
        default=0.0,
        help="Minimum difference between in-cluster and outside expression fractions.",
    )

    parser.add_argument(
        "--max-p-adj",
        type=float,
        default=1.0,
        help=(
            "Maximum Benjamini-Hochberg adjusted "
            "p-value allowed for a marker gene."
        ),
    )

    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--layer",
        default=None,
        help="H5AD layer name; omit to use AnnData.X.",
    )

    source.add_argument(
        "--use-raw",
        action="store_true",
        help="Use AnnData.raw expression and gene names (H5AD only).",
    )

    return parser


def main() -> None:
    """Run the command-line application."""

    args = build_parser().parse_args()

    def report_progress(completed: int, total: int, cluster: object) -> None:
        print(
            f"[{completed}/{total}] Completed cluster {cluster}",
            file=sys.stderr,
            flush=True,
        )

    analysis_options = dict(
        cluster_column=args.cluster_column,
        cell_column=args.cell_column,
        top_n=args.top_n,
        pseudocount=args.pseudocount,
        min_pct=args.min_pct,
        min_cells_expressing=args.min_cells_expressing,
        min_cluster_mean=args.min_cluster_mean,
        min_mean_difference=args.min_mean_difference,
        max_other_mean=args.max_other_mean,
        max_cells_expressing_out=args.max_cells_expressing_out,
        max_pct_out=args.max_pct_out,
        rank_by=args.rank_by,
        min_pct_difference=args.min_pct_difference,
        min_log2fc=args.min_log2fc,
        max_p_adj=args.max_p_adj,
        engine=args.engine,
        layer=args.layer,
        use_raw=args.use_raw,
        target_clusters=args.target_clusters,
        min_cluster_size=args.min_cluster_size,
        excluded_genes=args.exclude_gene,
        excluded_prefixes=args.exclude_prefix,
        max_dense_memory_mb=args.max_dense_memory_mb,
        gene_chunk_size=args.gene_chunk_size,
        progress_callback=None if args.quiet else report_progress,
        workers=args.workers,
    )

    if bool(args.input_dir) != bool(args.output_dir):
        raise SystemExit("--input-dir and --output-dir must be used together.")
    if bool(args.input) != bool(args.output):
        raise SystemExit("--input and --output must be used together.")
    if args.input_dir and args.diagnostics_json:
        raise SystemExit(
            "--diagnostics-json is for single-file analysis; batch diagnostics "
            "are written to batch_manifest.json."
        )

    if args.input_dir:
        manifest = run_batch(args.input_dir, args.output_dir, **analysis_options)
        completed = sum(row["status"] == "complete" for row in manifest)
        failed = len(manifest) - completed
        print(
            f"Batch complete: {completed} succeeded, {failed} failed. "
            f"Manifest: {Path(args.output_dir) / 'batch_manifest.json'}"
        )
        if failed:
            raise SystemExit(1)
        return

    diagnostics: dict[str, object] = {}
    markers = run(
        input_path=args.input,
        output_path=args.output,
        diagnostics=diagnostics,
        **analysis_options,
    )

    if args.diagnostics_json:
        diagnostics_path = Path(args.diagnostics_json)
        diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostics_path.write_text(
            json.dumps(diagnostics, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    print(
        f"Saved {len(markers)} marker rows "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
