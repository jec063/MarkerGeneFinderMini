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
from .settings import (
    SETTINGS_SCHEMA_NAME,
    SETTINGS_SCHEMA_VERSION,
    build_settings_document,
    load_analysis_settings,
    migrate_settings_document,
)
from .version import APP_VERSION

__all__ = [
    "APP_VERSION",
    "CONFIG_SCHEMA_VERSION",
    "AnalysisConfig",
    "ConfigurationError",
    "InputValidationError",
    "MarkerFinderError",
    "MarkerResult",
    "SchemaError",
    "SETTINGS_SCHEMA_NAME",
    "SETTINGS_SCHEMA_VERSION",
    "analyze",
    "analyze_file",
    "build_settings_document",
    "load_analysis_settings",
    "migrate_settings_document",
]
