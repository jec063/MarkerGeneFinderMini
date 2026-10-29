"""Contract tests for the formal v2.3 analysis-engine interface."""

from __future__ import annotations

import subprocess
import sys

import pandas as pd
import pytest

from marker_gene_finder import (
    AnalysisConfig,
    ConfigurationError,
    EngineOutput,
    EngineRegistry,
    analyze,
)
from marker_gene_finder.engines import DEFAULT_ENGINE_REGISTRY


def expression_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell": ["a", "b"],
            "cluster": ["A", "B"],
            "gene": [2.0, 0.0],
        }
    )


class RecordingEngine:
    name = "native"

    def __init__(self) -> None:
        self.config: AnalysisConfig | None = None

    def analyze_dataframe(
        self, expression: pd.DataFrame, config: AnalysisConfig
    ) -> EngineOutput:
        self.config = config
        return EngineOutput(
            markers=pd.DataFrame({"cluster": ["A"], "gene": ["gene"]}),
            diagnostics={"engine": "recording", "input_kind": "dataframe"},
        )

    def analyze_file(self, input_path, output_path, config) -> EngineOutput:
        raise AssertionError("file analysis was not requested")


def test_public_api_dispatches_through_engine_registry() -> None:
    engine = RecordingEngine()
    registry = EngineRegistry((engine,))
    config = AnalysisConfig(top_n=1)

    result = analyze(expression_table(), config, registry=registry)

    assert engine.config is config
    assert result.metadata["engine"] == "recording"
    assert result.metadata["marker_rows"] == 1


def test_default_registry_exposes_all_builtin_engines() -> None:
    assert DEFAULT_ENGINE_REGISTRY.names() == ("auto", "native", "scanpy")


def test_registry_rejects_duplicates_and_incomplete_engines() -> None:
    registry = EngineRegistry((RecordingEngine(),))

    with pytest.raises(ConfigurationError, match="already registered"):
        registry.register(RecordingEngine())
    with pytest.raises(ConfigurationError, match="contract"):
        registry.register(type("Incomplete", (), {"name": "broken"})())


def test_legacy_import_emits_deprecation_notice() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-W",
            "always::DeprecationWarning",
            "-c",
            "import src.filter_presets",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "src.filter_presets is deprecated" in completed.stderr
    assert "removed in Marker Gene Finder 3.0" in completed.stderr
