"""Stable v2 Python API built on the proven v1 analysis engines."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import AnalysisConfig
from .exceptions import InputValidationError
from .engines import DEFAULT_ENGINE_REGISTRY, EngineRegistry
from .models import MarkerResult
from .version import APP_VERSION


def analyze(
    expression: pd.DataFrame,
    config: AnalysisConfig | None = None,
    *,
    registry: EngineRegistry = DEFAULT_ENGINE_REGISTRY,
) -> MarkerResult:
    """Analyze an in-memory expression table through a registered engine."""
    config = config or AnalysisConfig()
    try:
        output = registry.get(config.engine).analyze_dataframe(expression, config)
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    metadata = {
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
        **output.diagnostics,
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
) -> MarkerResult:
    """Analyze a CSV or H5AD file and write its marker table."""
    config = config or AnalysisConfig()
    try:
        output = registry.get(config.engine).analyze_file(
            input_path, output_path, config
        )
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    metadata = {
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
        **output.diagnostics,
    }
    return MarkerResult(
        markers=output.markers, config=config, metadata=metadata
    )
