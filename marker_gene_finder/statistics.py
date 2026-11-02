"""Native and Scanpy marker statistics with stable ranking."""
from __future__ import annotations

from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu
from ._constants import DEFAULT_GENE_CHUNK_SIZE, RANK_BY_CHOICES, RESULT_COLUMNS
from ._validation import _validate_gene_chunk_size, _validate_workers
from .input_loading import _raw_expression
from .progress import ProgressCallback, report_progress


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
            report_progress(progress_callback, cluster_index, len(clusters), cluster)
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
        report_progress(progress_callback, cluster_index, len(clusters), cluster)

    if not results:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    return pd.concat(
        results,
        ignore_index=True,
    )[RESULT_COLUMNS]
