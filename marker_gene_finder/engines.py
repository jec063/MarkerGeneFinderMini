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
class EngineCapabilities:
    """Optional v2.6 declarations; omitted declarations retain v2 compatibility."""

    input_formats: tuple[str, ...] = ("csv", "h5ad", "dataframe")
    supports_layers: bool = True
    supports_raw: bool = True
    supports_workers: bool = True

    def __post_init__(self):
        if (not isinstance(self.input_formats, tuple) or not self.input_formats
                or any(x not in {"csv", "h5ad", "dataframe"} for x in self.input_formats)):
            raise ConfigurationError("input_formats must be a nonempty tuple of csv, h5ad, dataframe.")
        for name in ("supports_layers", "supports_raw", "supports_workers"):
            if not isinstance(getattr(self, name), bool):
                raise ConfigurationError(f"{name} must be boolean.")


def check_capabilities(engine, input_kind, config):
    """Reject unsupported requests before invoking the engine."""
    caps = getattr(engine, "capabilities", EngineCapabilities())
    if not isinstance(caps, EngineCapabilities):
        raise ConfigurationError("capabilities must be EngineCapabilities.")
    if input_kind not in caps.input_formats:
        raise InputValidationError(f"Engine {engine.name!r} does not support {input_kind} input.")
    if config.layer is not None and (input_kind != "h5ad" or not caps.supports_layers):
        raise InputValidationError(f"Engine {engine.name!r} does not support layers for this input.")
    if config.use_raw and (input_kind != "h5ad" or not caps.supports_raw):
        raise InputValidationError(f"Engine {engine.name!r} does not support raw for this input.")
    if config.workers != 1 and not caps.supports_workers:
        raise InputValidationError(f"Engine {engine.name!r} requires workers=1.")


@dataclass(frozen=True, slots=True)
class EngineOutput:
    """Normalized output returned by every analysis engine."""

    markers: pd.DataFrame
    diagnostics: dict[str, object]


@dataclass(frozen=True, slots=True)
class EngineInfo:
    """Serializable discovery details for one built-in or plugin engine."""

    name: str
    source: str
    status: str
    entry_point: str | None = None
    distribution: str | None = None
    version: str | None = None
    error: str | None = None

    def to_dict(self) -> dict[str, str | None]:
        """Return a stable JSON-compatible representation."""
        return {
            "name": self.name,
            "source": self.source,
            "status": self.status,
            "entry_point": self.entry_point,
            "distribution": self.distribution,
            "version": self.version,
            "error": self.error,
        }


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
    capabilities = EngineCapabilities()

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        options = _dataframe_kwargs(config)
        if config.target_clusters and config.cluster_column in expression.columns:
            options["target_clusters"] = _core._resolve_cluster_arguments(
                expression[config.cluster_column], list(config.target_clusters)
            )
        markers = _core.find_markers(expression, **options)
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
    capabilities = EngineCapabilities(input_formats=("h5ad",), supports_workers=False)

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
    capabilities = EngineCapabilities()

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


def validate_engine(
    engine: object, *, expected_name: str | None = None
) -> AnalysisEngine:
    """Validate an engine contract and optionally its advertised name."""
    name = validate_engine_name(getattr(engine, "name", None))
    if not isinstance(getattr(engine, "capabilities", EngineCapabilities()), EngineCapabilities):
        raise ConfigurationError("capabilities must be EngineCapabilities.")
    if not isinstance(engine, AnalysisEngine):
        raise ConfigurationError(
            f"Engine {name!r} does not implement the AnalysisEngine contract."
        )
    if expected_name is not None and name != expected_name:
        raise ConfigurationError(
            f"Engine plugin entry point {expected_name!r} returned engine "
            f"{name!r}; the names must match."
        )
    return engine


class EngineRegistry:
    """Validated collection of engines used by the public API."""

    def __init__(self, engines: tuple[AnalysisEngine, ...] = ()) -> None:
        self._engines: dict[str, AnalysisEngine] = {}
        for engine in engines:
            self.register(engine)

    def register(self, engine: AnalysisEngine, *, replace: bool = False) -> None:
        """Register an engine, rejecting invalid or duplicate names."""
        engine = validate_engine(engine)
        name = engine.name
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
    return validate_engine(candidate, expected_name=entry_point.name)


def _entry_point_distribution(
    entry_point: metadata.EntryPoint,
) -> tuple[str | None, str | None]:
    distribution = getattr(entry_point, "dist", None)
    if distribution is None:
        return None, None
    return (
        getattr(distribution, "name", None),
        getattr(distribution, "version", None),
    )


def inspect_engines(
    *, load_plugins: bool = False, group: str = ENGINE_ENTRY_POINT_GROUP
) -> tuple[EngineInfo, ...]:
    """Inspect built-ins and plugins without aborting on broken plugins."""
    records = [
        EngineInfo(name=name, source="built-in", status="available")
        for name in DEFAULT_ENGINE_REGISTRY.names()
    ]
    if not load_plugins:
        return tuple(records)
    registered = set(DEFAULT_ENGINE_REGISTRY.names())
    for entry_point in sorted(_entry_points(group), key=lambda item: item.name):
        distribution, version = _entry_point_distribution(entry_point)
        try:
            engine = _plugin_engine(entry_point)
            if engine.name in registered:
                raise ConfigurationError(
                    f"Engine {engine.name!r} is already registered."
                )
            registered.add(engine.name)
            records.append(
                EngineInfo(
                    name=engine.name,
                    source="plugin",
                    status="available",
                    entry_point=getattr(entry_point, "value", None),
                    distribution=distribution,
                    version=version,
                )
            )
        except ConfigurationError as error:
            records.append(
                EngineInfo(
                    name=entry_point.name,
                    source="plugin",
                    status="error",
                    entry_point=getattr(entry_point, "value", None),
                    distribution=distribution,
                    version=version,
                    error=str(error),
                )
            )
    return tuple(records)


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
    "EngineInfo",
    "EngineOutput",
    "EngineCapabilities",
    "EngineRegistry",
    "NativeEngine",
    "ScanpyEngine",
    "inspect_engines",
    "load_engine_plugins",
    "validate_engine",
]
