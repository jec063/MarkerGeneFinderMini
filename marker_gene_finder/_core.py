"""Compatibility facade for historical imports; implementations live in services."""

from ._validation import (
    _validate_workers,
    _validate_gene_chunk_size,
    _validate_max_dense_memory_mb,
)
from .input_loading import (
    _raw_expression,
    estimate_dense_expression_bytes,
    _enforce_dense_memory_limit,
    load_expression,
)
from .inspection import inspect_input, build_analysis_plan
from .statistics import (
    _rank_markers,
    _benjamini_hochberg,
    _selected_clusters,
    _resolve_cluster_arguments,
    _filter_genes,
    _gene_columns,
    find_markers,
    find_markers_scanpy,
)
from .analysis import run, _registered_run
from .batch import (
    _sha256_file,
    _select_batch_inputs,
    _batch_settings_fingerprint,
    build_batch_plan,
    run_batch,
    verify_batch_manifest,
)
from .benchmark import benchmark_analysis, compare_benchmarks
from ._legacy_cli import build_parser, main
from ._constants import (
    RESULT_COLUMNS,
    RANK_BY_CHOICES,
    DEFAULT_MAX_DENSE_MEMORY_MB,
    DEFAULT_GENE_CHUNK_SIZE,
    CSV_INSPECTION_CHUNK_ROWS,
    PLAN_SCHEMA_VERSION,
    BATCH_MANIFEST_SCHEMA_VERSION,
    BENCHMARK_SCHEMA_VERSION,
    BENCHMARK_COMPARISON_SCHEMA_VERSION,
    BATCH_VERIFICATION_SCHEMA_VERSION,
    BATCH_PLAN_SCHEMA_VERSION,
)
from .progress import ProgressCallback
from .version import APP_VERSION


# Historical src.marker_finder aliases this module. Keep patches made through
# that alias visible to service functions that formerly shared this namespace.
import sys as _sys
from types import ModuleType as _ModuleType


class _CompatibilityModule(_ModuleType):
    def __setattr__(self, name, value):
        if name in self.__dict__ and not name.startswith("__"):
            previous = self.__dict__[name]
            for module_name in (
                "_validation", "input_loading", "inspection", "statistics",
                "analysis", "batch", "benchmark", "_legacy_cli", "_constants", "engines",
            ):
                module = _sys.modules.get(f"marker_gene_finder.{module_name}")
                if module is not None and module.__dict__.get(name) is previous:
                    setattr(module, name, value)
        super().__setattr__(name, value)


_sys.modules[__name__].__class__ = _CompatibilityModule

if __name__ == "__main__":
    main()
