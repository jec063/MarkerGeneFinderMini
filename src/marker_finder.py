"""Find simple marker genes by comparing one cluster with all other cells."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

RESULT_COLUMNS = [
    "cluster",
    "gene",
    "cluster_mean",
    "other_mean",
    "pct_in",
    "pct_out",
    "log2FC",
    "p_value",
    "p_adj",
]



def load_expression(
    input_path: str | Path,
    cluster_column: str = "cluster",
    cell_column: str | None = "cell_id",
    layer: str | None = None,
) -> pd.DataFrame:
    """Load CSV or AnnData input as an expression table."""
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

    adata = ad.read_h5ad(input_path)

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

    if expression[genes].isna().any().any():
        raise ValueError(
            "Gene-expression values cannot contain missing values."
        )

    if (expression[genes] < 0).any().any():
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
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
) -> pd.DataFrame:
    """Return the top fold-change-ranked genes for every cluster."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1.")

    if pseudocount <= 0:
        raise ValueError("pseudocount must be greater than 0.")

    if not 0 <= min_pct <= 1:
        raise ValueError("min_pct must be between 0 and 1.")

    if min_log2fc < 0:
        raise ValueError("min_log2fc must be at least 0.")

    if not 0 <= max_p_adj <= 1:
        raise ValueError("max_p_adj must be between 0 and 1.")

    genes = _gene_columns(
        expression,
        cluster_column,
        cell_column,
    )

    results = []

    clusters = sorted(
        expression[cluster_column].unique(),
        key=str,
    )

    for cluster in clusters:
        in_cluster = expression[cluster_column] == cluster

        cluster_values = expression.loc[in_cluster, genes]
        other_values = expression.loc[~in_cluster, genes]

        cluster_mean = cluster_values.mean(axis=0)
        other_mean = other_values.mean(axis=0)

        pct_in = cluster_values.gt(0).mean(axis=0)
        pct_out = other_values.gt(0).mean(axis=0)

        log2fc = np.log2(
            (cluster_mean.to_numpy() + pseudocount)
            / (other_mean.to_numpy() + pseudocount)
        )

        test_result = mannwhitneyu(
            cluster_values.to_numpy(),
            other_values.to_numpy(),
            alternative="greater",
            axis=0,
            method="asymptotic",
        )

        p_values = np.asarray(
            test_result.pvalue,
            dtype=float,
        )
        p_adjusted = _benjamini_hochberg(p_values)

        cluster_result = pd.DataFrame(
            {
                "cluster": cluster,
                "gene": genes,
                "cluster_mean": cluster_mean.to_numpy(),
                "other_mean": other_mean.to_numpy(),
                "pct_in": pct_in.to_numpy(),
                "pct_out": pct_out.to_numpy(),
                "log2FC": log2fc,
                "p_value": p_values,
                "p_adj": p_adjusted,
            }
        )

        cluster_result = (
            cluster_result
                        .loc[
                (cluster_result["pct_in"] >= min_pct)
                & (cluster_result["log2FC"] >= min_log2fc)
                & (cluster_result["p_adj"] <= max_p_adj)
            ]
            .sort_values(
                ["log2FC", "gene"],
                ascending=[False, True],
            )
            .head(top_n)
        )

        results.append(cluster_result)

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
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
    layer: str | None = None,
) -> pd.DataFrame:
    """Find marker genes using Scanpy's Wilcoxon implementation."""

    if top_n < 1:
        raise ValueError("top_n must be at least 1.")

    if not 0 <= min_pct <= 1:
        raise ValueError("min_pct must be between 0 and 1.")

    if min_log2fc < 0:
        raise ValueError("min_log2fc must be at least 0.")

    if not 0 <= max_p_adj <= 1:
        raise ValueError("max_p_adj must be between 0 and 1.")

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
    working.obs[cluster_column] = (
        working.obs[cluster_column].astype("category")
    )

    sc.tl.rank_genes_groups(
        working,
        groupby=cluster_column,
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

    gene_positions = pd.Series(
        np.arange(working.n_vars),
        index=working.var_names.astype(str),
    )

    results = []

    for cluster in working.obs[cluster_column].cat.categories:
        in_cluster = np.asarray(
            working.obs[cluster_column] == cluster,
            dtype=bool,
        )

        cluster_mean = np.asarray(
            working.X[in_cluster].mean(axis=0)
        ).ravel()
        other_mean = np.asarray(
            working.X[~in_cluster].mean(axis=0)
        ).ravel()

        cluster_result = ranked.loc[
            ranked["cluster"].astype(str) == str(cluster)
        ].copy()

        positions = gene_positions.loc[
            cluster_result["gene"].astype(str)
        ].to_numpy()

        cluster_result["cluster"] = cluster
        cluster_result["cluster_mean"] = cluster_mean[positions]
        cluster_result["other_mean"] = other_mean[positions]

        cluster_result = (
            cluster_result.loc[
                (cluster_result["pct_in"] >= min_pct)
                & (cluster_result["log2FC"] >= min_log2fc)
                & (cluster_result["p_adj"] <= max_p_adj)
            ]
            .sort_values(
                ["log2FC", "gene"],
                ascending=[False, True],
            )
            .head(top_n)
        )

        results.append(cluster_result)

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
    min_log2fc: float = 0.0,
    max_p_adj: float = 1.0,
    engine: str = "native",
    layer: str | None = None,
) -> pd.DataFrame:
    """Read an expression dataset, find markers and save the results."""

    if engine not in {"native", "scanpy"}:
        raise ValueError(
            "engine must be either 'native' or 'scanpy'."
        )

    input_path = Path(input_path)

    if engine == "scanpy":
        if input_path.suffix.lower() != ".h5ad":
            raise ValueError(
                "Scanpy engine requires a .h5ad input file."
            )

        adata = ad.read_h5ad(input_path)

        markers = find_markers_scanpy(
            adata,
            layer=layer,
            cluster_column=cluster_column,
            top_n=top_n,
            min_pct=min_pct,
            min_log2fc=min_log2fc,
            max_p_adj=max_p_adj,
        )
    else:
        expression = load_expression(
            input_path,
            layer=layer,
            cluster_column=cluster_column,
            cell_column=cell_column,
        )

        markers = find_markers(
            expression,
            cluster_column=cluster_column,
            cell_column=cell_column,
            top_n=top_n,
            pseudocount=pseudocount,
            min_pct=min_pct,
            min_log2fc=min_log2fc,
            max_p_adj=max_p_adj,
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    markers.to_csv(output_path, index=False)

    return markers

def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Rank marker genes for each cluster "
            "using log2 fold change."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input expression CSV or H5AD path.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output marker CSV path.",
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
        "--min-log2fc",
        type=float,
        default=0.0,
        help=(
            "Minimum log2 fold change required "
            "for a marker gene."
        ),
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

    parser.add_argument(
        "--layer",
        default=None,
        help="H5AD layer name; omit to use AnnData.X.",
    )

    return parser


def main() -> None:
    """Run the command-line application."""

    args = build_parser().parse_args()

    markers = run(
        input_path=args.input,
        output_path=args.output,
        cluster_column=args.cluster_column,
        cell_column=args.cell_column,
        top_n=args.top_n,
        pseudocount=args.pseudocount,
        min_pct=args.min_pct,
        min_log2fc=args.min_log2fc,
        max_p_adj=args.max_p_adj,
        engine=args.engine,
        layer=args.layer,
    )

    print(
        f"Saved {len(markers)} marker rows "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
