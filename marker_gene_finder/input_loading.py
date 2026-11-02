"""Expression loading and dense memory safeguards."""
from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from ._constants import DEFAULT_MAX_DENSE_MEMORY_MB
from ._validation import _validate_max_dense_memory_mb


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
    _validate_max_dense_memory_mb(max_dense_memory_mb)
    if max_dense_memory_mb is None:
        return
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
