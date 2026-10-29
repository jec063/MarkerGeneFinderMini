"""Stable public API for Marker Gene Finder v2."""

from .api import analyze, analyze_file
from .config import AnalysisConfig, CONFIG_SCHEMA_VERSION
from .exceptions import (
    ConfigurationError,
    InputValidationError,
    MarkerFinderError,
    SchemaError,
)
from .models import MarkerResult
from .engines import (
    DEFAULT_ENGINE_REGISTRY,
    ENGINE_ENTRY_POINT_GROUP,
    AnalysisEngine,
    AutoEngine,
    EngineOutput,
    EngineRegistry,
    NativeEngine,
    ScanpyEngine,
    load_engine_plugins,
)
from .settings import (
    SETTINGS_SCHEMA_NAME,
    SETTINGS_SCHEMA_VERSION,
    build_settings_document,
    load_analysis_settings,
    migrate_settings_document,
)
from .version import APP_VERSION

__all__ = [
    "load_engine_plugins",
    "APP_VERSION",
    "CONFIG_SCHEMA_VERSION",
    "AnalysisConfig",
    "ConfigurationError",
    "DEFAULT_ENGINE_REGISTRY",
    "ENGINE_ENTRY_POINT_GROUP",
    "AnalysisEngine",
    "AutoEngine",
    "EngineOutput",
    "EngineRegistry",
    "InputValidationError",
    "MarkerFinderError",
    "MarkerResult",
    "NativeEngine",
    "ScanpyEngine",
    "SchemaError",
    "SETTINGS_SCHEMA_NAME",
    "SETTINGS_SCHEMA_VERSION",
    "analyze",
    "analyze_file",
    "build_settings_document",
    "load_analysis_settings",
    "migrate_settings_document",
]
