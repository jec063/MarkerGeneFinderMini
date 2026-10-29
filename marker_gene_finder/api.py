"""Stable v2 Python API built on the proven v1 analysis engines."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import AnalysisConfig
from .exceptions import ConfigurationError, InputValidationError
from .engines import DEFAULT_ENGINE_REGISTRY, EngineOutput, EngineRegistry
from .models import MarkerResult
from .version import APP_VERSION


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
    try:
        output = engine.analyze_dataframe(expression, config)
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    output = _validated_output(config.engine, output)
    metadata = {
        **output.diagnostics,
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
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
    engine = registry.get(config.engine)
    try:
        output = engine.analyze_file(input_path, output_path, config)
    except ValueError as error:
        if isinstance(error, InputValidationError):
            raise
        raise InputValidationError(str(error)) from error
    output = _validated_output(config.engine, output)
    metadata = {
        **output.diagnostics,
        "app_version": APP_VERSION,
        "marker_rows": len(output.markers),
    }
    return MarkerResult(
        markers=output.markers, config=config, metadata=metadata
    )
