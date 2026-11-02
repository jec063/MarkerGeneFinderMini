"""Shared result columns, defaults, and artifact schema versions."""

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
CSV_INSPECTION_CHUNK_ROWS = 10_000
PLAN_SCHEMA_VERSION = 1
BATCH_MANIFEST_SCHEMA_VERSION = 2
BENCHMARK_SCHEMA_VERSION = 1
BENCHMARK_COMPARISON_SCHEMA_VERSION = 1
BATCH_VERIFICATION_SCHEMA_VERSION = 1
BATCH_PLAN_SCHEMA_VERSION = 1
