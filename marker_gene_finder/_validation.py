"""Shared execution option validation."""
from __future__ import annotations

import numpy as np


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


def _validate_max_dense_memory_mb(max_dense_memory_mb: float | None) -> None:
    """Validate the optional Native dense-conversion safety limit."""
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
