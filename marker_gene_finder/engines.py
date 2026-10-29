"""Formal marker-analysis engine contracts and built-in implementations."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Protocol, runtime_checkable

import pandas as pd

from . import _core
from .config import AnalysisConfig, validate_engine_name
from .exceptions import ConfigurationError, InputValidationError

ENGINE_ENTRY_POINT_GROUP = "marker_gene_finder.engines"


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

    def copy(self) -> "EngineRegistry":
        """Return an independent registry containing the same engines."""
        return EngineRegistry(tuple(self._engines.values()))


DEFAULT_ENGINE_REGISTRY = EngineRegistry(
    (NativeEngine(), ScanpyEngine(), AutoEngine())
)


def _entry_points(group: str) -> Iterable[metadata.EntryPoint]:
    """Return entry points using the supported importlib metadata APIs."""
    discovered = metadata.entry_points()
    if hasattr(discovered, "select"):
        return discovered.select(group=group)
    return discovered.get(group, ())


def _plugin_engine(entry_point: metadata.EntryPoint) -> AnalysisEngine:
    """Load one engine instance from an entry point."""
    try:
        candidate = entry_point.load()
        if isinstance(candidate, type):
            candidate = candidate()
        elif not isinstance(candidate, AnalysisEngine) and isinstance(
            candidate, Callable
        ):
            candidate = candidate()
    except Exception as error:
        raise ConfigurationError(
            f"Cannot load engine plugin {entry_point.name!r}: {error}"
        ) from error
    if not isinstance(candidate, AnalysisEngine):
        raise ConfigurationError(
            f"Engine plugin {entry_point.name!r} does not implement the "
            "AnalysisEngine contract."
        )
    if candidate.name != entry_point.name:
        raise ConfigurationError(
            f"Engine plugin entry point {entry_point.name!r} returned engine "
            f"{candidate.name!r}; the names must match."
        )
    return candidate


def load_engine_plugins(
    registry: EngineRegistry | None = None,
    *,
    group: str = ENGINE_ENTRY_POINT_GROUP,
    replace: bool = False,
) -> EngineRegistry:
    """Load installed engine entry points into an independent registry."""
    result = (registry or DEFAULT_ENGINE_REGISTRY).copy()
    for entry_point in sorted(_entry_points(group), key=lambda item: item.name):
        result.register(_plugin_engine(entry_point), replace=replace)
    return result


__all__ = [
    "DEFAULT_ENGINE_REGISTRY",
    "ENGINE_ENTRY_POINT_GROUP",
    "AnalysisEngine",
    "AutoEngine",
    "EngineOutput",
    "EngineRegistry",
    "NativeEngine",
    "ScanpyEngine",
    "load_engine_plugins",
]
