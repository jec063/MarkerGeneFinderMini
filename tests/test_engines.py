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
    inspect_engines,
    load_engine_plugins,
    validate_engine,
)
from marker_gene_finder.engines import DEFAULT_ENGINE_REGISTRY, metadata


def expression_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell": ["a", "b"],
            "cluster": ["A", "B"],
            "gene": [2.0, 0.0],
        }
    )


class RecordingEngine:
    name = "recording"

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
    config = AnalysisConfig(engine="recording", top_n=1)

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


def test_registry_rejects_invalid_names_and_reports_unknown_engines() -> None:
    registry = EngineRegistry()

    with pytest.raises(ConfigurationError, match="lowercase"):
        registry.register(type("BadName", (), {"name": "Not Valid"})())
    with pytest.raises(ConfigurationError, match="Unknown analysis engine"):
        analyze(
            expression_table(),
            AnalysisConfig(engine="missing"),
            registry=registry,
        )


def test_public_api_rejects_invalid_custom_engine_output() -> None:
    class InvalidOutputEngine(RecordingEngine):
        name = "invalid"

        def analyze_dataframe(self, expression, config):
            return expression

    with pytest.raises(ConfigurationError, match="expected EngineOutput"):
        analyze(
            expression_table(),
            AnalysisConfig(engine="invalid"),
            registry=EngineRegistry((InvalidOutputEngine(),)),
        )


def test_public_api_keeps_authoritative_result_metadata() -> None:
    class MisleadingEngine(RecordingEngine):
        name = "misleading"

        def analyze_dataframe(self, expression, config):
            return EngineOutput(
                markers=pd.DataFrame({"gene": ["real"]}),
                diagnostics={"app_version": "0.0.0", "marker_rows": 99},
            )

    result = analyze(
        expression_table(),
        AnalysisConfig(engine="misleading"),
        registry=EngineRegistry((MisleadingEngine(),)),
    )

    assert result.metadata["app_version"] == "2.5.0"
    assert result.metadata["marker_rows"] == 1


class FakeEntryPoint:
    def __init__(self, name, value=None, error=None):
        self.name = name
        self.value = value
        self.error = error

    def load(self):
        if self.error is not None:
            raise self.error
        return self.value


class FakeEntryPoints(list):
    def select(self, *, group):
        assert group == "marker_gene_finder.engines"
        return self


def test_engine_plugins_are_loaded_without_mutating_defaults(monkeypatch) -> None:
    class PluginEngine(RecordingEngine):
        name = "plugin"

    monkeypatch.setattr(
        metadata,
        "entry_points",
        lambda: FakeEntryPoints([FakeEntryPoint("plugin", PluginEngine)]),
    )

    registry = load_engine_plugins()

    assert registry.names() == ("auto", "native", "plugin", "scanpy")
    assert DEFAULT_ENGINE_REGISTRY.names() == ("auto", "native", "scanpy")
    result = analyze(
        expression_table(),
        AnalysisConfig(engine="plugin"),
        registry=registry,
    )
    assert result.metadata["engine"] == "recording"


def test_engine_plugin_names_must_match(monkeypatch) -> None:
    monkeypatch.setattr(
        metadata,
        "entry_points",
        lambda: FakeEntryPoints(
            [FakeEntryPoint("different", RecordingEngine())]
        ),
    )

    with pytest.raises(ConfigurationError, match="names must match"):
        load_engine_plugins()


def test_engine_plugin_load_failures_are_public_errors(monkeypatch) -> None:
    monkeypatch.setattr(
        metadata,
        "entry_points",
        lambda: FakeEntryPoints(
            [FakeEntryPoint("broken", error=RuntimeError("unavailable"))]
        ),
    )

    with pytest.raises(ConfigurationError, match="broken.*unavailable"):
        load_engine_plugins()


def test_engine_inspection_reports_plugins_and_failures(monkeypatch) -> None:
    class PluginEngine(RecordingEngine):
        name = "plugin"

    monkeypatch.setattr(
        metadata,
        "entry_points",
        lambda: FakeEntryPoints(
            [
                FakeEntryPoint("plugin", PluginEngine),
                FakeEntryPoint("broken", error=RuntimeError("unavailable")),
            ]
        ),
    )

    records = inspect_engines(load_plugins=True)

    assert [record.name for record in records] == [
        "auto",
        "native",
        "scanpy",
        "broken",
        "plugin",
    ]
    assert records[3].status == "error"
    assert records[3].error is not None
    assert "unavailable" in records[3].error
    assert records[4].status == "available"


def test_public_engine_validator_checks_expected_name() -> None:
    assert validate_engine(RecordingEngine()).name == "recording"
    with pytest.raises(ConfigurationError, match="names must match"):
        validate_engine(RecordingEngine(), expected_name="different")


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
