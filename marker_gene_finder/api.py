"""Stable v2 Python API built on the proven v1 analysis engines."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.marker_finder import find_markers, run

from .config import AnalysisConfig
from .exceptions import InputValidationError
from .models import MarkerResult
from .version import APP_VERSION


def analyze(
    expression: pd.DataFrame,
    config: AnalysisConfig | None = None,
) -> MarkerResult:
    """Analyze an in-memory expression table with the Native engine."""
    config = config or AnalysisConfig()
    if config.engine not in {"native", "auto"}:
        raise InputValidationError(
            "In-memory DataFrame analysis supports the Native engine only."
        )
    kwargs = config.to_analysis_kwargs()
    for name in ("engine", "layer", "use_raw", "max_dense_memory_mb"):
        kwargs.pop(name)
    try:
        markers = find_markers(expression, **kwargs)
    except ValueError as error:
        raise InputValidationError(str(error)) from error
    return MarkerResult(
        markers=markers,
        config=config,
        metadata={
            "app_version": APP_VERSION,
            "engine": "native",
            "input_kind": "dataframe",
            "marker_rows": len(markers),
        },
    )


def analyze_file(
    input_path: str | Path,
    output_path: str | Path,
    config: AnalysisConfig | None = None,
) -> MarkerResult:
    """Analyze a CSV or H5AD file and write its marker table."""
    config = config or AnalysisConfig()
    diagnostics: dict[str, object] = {}
    try:
        markers = run(
            input_path=input_path,
            output_path=output_path,
            diagnostics=diagnostics,
            **config.to_analysis_kwargs(),
        )
    except ValueError as error:
        raise InputValidationError(str(error)) from error
    return MarkerResult(markers=markers, config=config, metadata=diagnostics)
