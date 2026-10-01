"""Summarize selected genes for expression dot plots."""

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def summarize_clusters(labels):
    """Return cell counts and percentages for observed clusters."""
    labels = pd.Series(labels, copy=False)
    if labels.empty or labels.isna().any():
        raise ValueError("Cluster labels must be nonempty and nonmissing.")
    counts = labels.value_counts(sort=False)
    ordered = sorted(counts.index, key=str)
    return pd.DataFrame({
        "cluster": ordered,
        "cells": [int(counts[label]) for label in ordered],
        "percent_of_cells": [100 * counts[label] / len(labels) for label in ordered],
    })


def summarize_expression(
    data,
    genes,
    cluster_column="cluster",
    layer=None,
    use_raw=False,
):
    """Return one row per observed cluster and selected gene."""
    if use_raw and layer is not None:
        raise ValueError("use_raw and layer are mutually exclusive.")

    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("Select at least one gene.")

    if isinstance(data, ad.AnnData):
        labels = data.obs[cluster_column]
        if use_raw:
            if data.raw is None:
                raise ValueError("The AnnData object has no raw expression.")
            matrix = data.raw.X
            names = data.raw.var_names
        else:
            if layer is not None and layer not in data.layers:
                raise ValueError(f"Missing AnnData layer: {layer!r}.")
            matrix = data.X if layer is None else data.layers[layer]
            names = data.var_names
        if matrix is None:
            raise ValueError("The selected expression matrix is missing.")
        if not names.is_unique:
            raise ValueError("Gene names must be unique.")
        positions = names.get_indexer(genes)
        if (positions < 0).any():
            missing = [g for g, p in zip(genes, positions) if p < 0]
            raise ValueError(f"Missing genes: {missing}")
        if sparse.issparse(matrix):
            matrix = matrix.tocsr()
        values = matrix[:, positions]
    elif isinstance(data, pd.DataFrame):
        if layer is not None or use_raw:
            raise ValueError("Select layers or raw before creating the table.")
        if not data.columns.is_unique:
            raise ValueError("Expression columns must be unique.")
        if cluster_column in genes:
            raise ValueError("The cluster column is not a gene.")
        labels = data[cluster_column]
        missing = [gene for gene in genes if gene not in data.columns]
        if missing:
            raise ValueError(f"Missing genes: {missing}")
        values = data.loc[:, genes].to_numpy(dtype=float)
    else:
        raise TypeError("Expected an AnnData object or expression DataFrame.")

    if len(labels) == 0 or labels.isna().any():
        raise ValueError("Cluster labels must be nonempty and nonmissing.")

    values = values.astype(float)
    stored = values.data if sparse.issparse(values) else values
    if not np.isfinite(stored).all() or (stored < 0).any():
        raise ValueError("Expression values must be finite and nonnegative.")

    summaries = []
    for cluster in sorted(labels.unique(), key=str):
        selected = values[np.asarray(labels == cluster)]
        means = np.asarray(selected.mean(axis=0)).ravel()
        fractions = np.asarray((selected > 0).mean(axis=0)).ravel()
        summaries.append(pd.DataFrame({
            "cluster": cluster,
            "gene": genes,
            "mean_expression": means,
            "fraction_expressing": fractions,
        }))
    return pd.concat(summaries, ignore_index=True)
