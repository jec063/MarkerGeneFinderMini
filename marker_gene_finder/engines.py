"""Formal marker-analysis engine contracts and built-in implementations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

import pandas as pd

from . import _core
from .config import AnalysisConfig, validate_engine_name
from .exceptions import ConfigurationError, InputValidationError


@dataclass(frozen=True, slots=True)
class EngineOutput:
    """Normalized output returned by every analysis engine."""

    markers: pd.DataFrame
    diagnostics: dict[str, object]


@runtime_checkable
class AnalysisEngine(Protocol):
    """Contract implemented by marker-analysis backends."""

    name: str

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        """Analyze an in-memory expression table."""

    def analyze_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        config: AnalysisConfig,
    ) -> EngineOutput:
        """Analyze one input file and write its marker table."""


def _dataframe_kwargs(config: AnalysisConfig) -> dict[str, object]:
    values = config.to_analysis_kwargs()
    for name in ("engine", "layer", "use_raw", "max_dense_memory_mb"):
        values.pop(name)
    return values


def _run_file(
    engine: str,
    input_path: str | Path,
    output_path: str | Path,
    config: AnalysisConfig,
) -> EngineOutput:
    diagnostics: dict[str, object] = {}
    values = config.to_analysis_kwargs()
    values["engine"] = engine
    markers = _core.run(
        input_path=input_path,
        output_path=output_path,
        diagnostics=diagnostics,
        **values,
    )
    return EngineOutput(markers=markers, diagnostics=diagnostics)


class NativeEngine:
    """Built-in engine for native one-versus-rest marker statistics."""

    name = "native"

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        markers = _core.find_markers(expression, **_dataframe_kwargs(config))
        return EngineOutput(
            markers=markers,
            diagnostics={"engine": self.name, "input_kind": "dataframe"},
        )

    def analyze_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        config: AnalysisConfig,
    ) -> EngineOutput:
        return _run_file(self.name, input_path, output_path, config)


class ScanpyEngine:
    """Built-in engine backed by Scanpy rank_genes_groups."""

    name = "scanpy"

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        raise InputValidationError(
            "The Scanpy engine requires an H5AD file; use analyze_file()."
        )

    def analyze_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        config: AnalysisConfig,
    ) -> EngineOutput:
        return _run_file(self.name, input_path, output_path, config)


class AutoEngine:
    """Built-in engine that follows the resource-aware execution plan."""

    name = "auto"

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        return NativeEngine().analyze_dataframe(expression, config)

    def analyze_file(
        self,
        input_path: str | Path,
        output_path: str | Path,
        config: AnalysisConfig,
    ) -> EngineOutput:
        return _run_file(self.name, input_path, output_path, config)


class EngineRegistry:
    """Validated collection of engines used by the public API."""

    def __init__(self, engines: tuple[AnalysisEngine, ...] = ()) -> None:
        self._engines: dict[str, AnalysisEngine] = {}
        for engine in engines:
            self.register(engine)

    def register(self, engine: AnalysisEngine, *, replace: bool = False) -> None:
        """Register an engine, rejecting invalid or duplicate names."""
        name = validate_engine_name(getattr(engine, "name", None))
        if not isinstance(engine, AnalysisEngine):
            raise ConfigurationError(
                f"Engine {name!r} does not implement the AnalysisEngine contract."
            )
        if name in self._engines and not replace:
            raise ConfigurationError(f"Engine {name!r} is already registered.")
        self._engines[name] = engine

    def get(self, name: str) -> AnalysisEngine:
        """Return a registered engine by name."""
        try:
            return self._engines[name]
        except KeyError as error:
            choices = ", ".join(self.names())
            raise ConfigurationError(
                f"Unknown analysis engine {name!r}; choose: {choices}."
            ) from error

    def names(self) -> tuple[str, ...]:
        """Return registered names in deterministic order."""
        return tuple(sorted(self._engines))


DEFAULT_ENGINE_REGISTRY = EngineRegistry(
    (NativeEngine(), ScanpyEngine(), AutoEngine())
)


__all__ = [
    "DEFAULT_ENGINE_REGISTRY",
    "AnalysisEngine",
    "AutoEngine",
    "EngineOutput",
    "EngineRegistry",
    "NativeEngine",
    "ScanpyEngine",
]
