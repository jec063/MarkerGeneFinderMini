"""Stable v2 Python API built on the proven v1 analysis engines."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import numpy as np
import json

from .config import AnalysisConfig
from .exceptions import ConfigurationError, InputValidationError
from .engines import DEFAULT_ENGINE_REGISTRY, EngineOutput, EngineRegistry, check_capabilities
from .models import MarkerResult
from .version import APP_VERSION
from .provenance import build_provenance, hash_input


def _validated_output(engine_name: str, output: object) -> EngineOutput:
    """Validate a third-party engine result before exposing it publicly."""
    if not isinstance(output, EngineOutput):
        raise ConfigurationError(
            f"Engine {engine_name!r} returned {type(output).__name__}; "
            "expected EngineOutput."
        )
    if not isinstance(output.markers, pd.DataFrame):
        raise ConfigurationError(
            f"Engine {engine_name!r} returned non-DataFrame markers."
        )
    if not isinstance(output.diagnostics, dict):
        raise ConfigurationError(
            f"Engine {engine_name!r} returned non-dictionary diagnostics."
        )
    try:
        json.dumps(output.diagnostics, allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ConfigurationError("Engine diagnostics must be finite JSON-compatible values.") from error
    table = output.markers
    if not table.columns.is_unique or not {"cluster", "gene"}.issubset(table.columns):
        raise ConfigurationError(f"Engine {engine_name!r} markers require unique cluster and gene columns.")
    if table[["cluster", "gene"]].isna().any().any():
        raise ConfigurationError("Marker cluster and gene values cannot be null.")
    for column in table.columns.intersection(["rank", "cluster_size", "reference_size", "cluster_mean", "other_mean", "mean_difference", "pct_in", "pct_out", "cells_expressing_in", "cells_expressing_out", "pct_difference", "log2FC", "p_value", "p_adj"]):
        values = table[column]
        if len(values) and (not pd.api.types.is_numeric_dtype(values) or pd.api.types.is_bool_dtype(values) or not np.isfinite(values).all()):
            raise ConfigurationError(f"Marker {column} must contain finite numeric values.")
        if column in {"p_value", "p_adj", "pct_in", "pct_out"} and not values.between(0, 1).all():
            raise ConfigurationError(f"Marker {column} must be between 0 and 1.")
        if column in {"rank", "cluster_size", "reference_size", "cells_expressing_in", "cells_expressing_out"} and ((values < (1 if column == "rank" else 0)).any() or (values % 1 != 0).any()):
            raise ConfigurationError(f"Marker {column} must contain valid integer counts.")
    return output


def analyze(
    expression: pd.DataFrame,
    config: AnalysisConfig | None = None,
    *,
    registry: EngineRegistry = DEFAULT_ENGINE_REGISTRY,
) -> MarkerResult:
    """Analyze an in-memory expression table through a registered engine."""
    config = config or AnalysisConfig()
    engine = registry.get(config.engine)
    check_capabilities(engine, "dataframe", config)
    try:
        output = engine.analyze_dataframe(expression, config)
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    output = _validated_output(config.engine, output)
    provenance = build_provenance(
        config, resolved_engine=("native" if config.engine == "auto" else config.engine),
        input_kind="dataframe", engine=engine,
    )
    metadata = {
        **output.diagnostics,
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
        "provenance": provenance,
    }
    return MarkerResult(
        markers=output.markers,
        config=config,
        metadata=metadata,
    )


def analyze_file(
    input_path: str | Path,
    output_path: str | Path,
    config: AnalysisConfig | None = None,
    *,
    registry: EngineRegistry = DEFAULT_ENGINE_REGISTRY,
    hash_input_file: bool = False,
) -> MarkerResult:
    """Analyze a CSV or H5AD file and write its marker table."""
    if not isinstance(hash_input_file, bool):
        raise ConfigurationError("hash_input_file must be boolean.")
    input_digest = hash_input(input_path) if hash_input_file else None
    config = config or AnalysisConfig()
    engine = registry.get(config.engine)
    check_capabilities(engine, Path(input_path).suffix.lower().lstrip("."), config)
    try:
        output = engine.analyze_file(input_path, output_path, config)
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    output = _validated_output(config.engine, output)
    provenance = build_provenance(
        config, resolved_engine=(str(output.diagnostics["resolved_engine"])
            if config.engine == "auto" else config.engine),
        input_kind=Path(input_path).suffix.lower().lstrip("."),
        engine=engine, input_sha256=input_digest,
    )
    metadata = {
        **output.diagnostics,
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
        "provenance": provenance,
    }
    return MarkerResult(
        markers=output.markers, config=config, metadata=metadata
    )
