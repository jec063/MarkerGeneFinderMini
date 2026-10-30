"""Reference plugin exercises v2.6 services without changing built-in statistics."""
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from marker_gene_finder import (
    AnalysisConfig, EngineCapabilities, EngineOutput, EngineRegistry,
    NativeEngine, ConfigurationError, InputValidationError, analyze, analyze_file,
)
from marker_gene_finder.batch import run_batch, build_batch_plan
from marker_gene_finder.benchmark import benchmark_analysis
from marker_gene_finder.inspection import build_analysis_plan
from marker_gene_finder.cli import main
from marker_gene_finder import cli

EXAMPLE = Path(__file__).resolve().parents[1] / "data/example_expression.csv"


class ReferenceEngine(NativeEngine):
    name = "reference"
    capabilities = EngineCapabilities(input_formats=("csv", "dataframe"),
                                      supports_layers=False, supports_raw=False,
                                      supports_workers=False)

    def analyze_file(self, input_path, output_path, config):
        return NativeEngine().analyze_file(input_path, output_path,
                                          replace(config, engine="native"))


@pytest.fixture
def registry():
    return EngineRegistry((ReferenceEngine(),))


def test_reference_plugin_all_workflows(tmp_path, registry):
    config = AnalysisConfig(engine="reference")
    expected = analyze(pd.read_csv(EXAMPLE)).markers
    pd.testing.assert_frame_equal(analyze(pd.read_csv(EXAMPLE), config, registry=registry).markers, expected)
    result = analyze_file(EXAMPLE, tmp_path / "single.csv", config, registry=registry)
    pd.testing.assert_frame_equal(result.markers, expected)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "sample.csv").write_bytes(EXAMPLE.read_bytes())
    output = tmp_path / "batch"
    assert build_analysis_plan(EXAMPLE, requested_engine="reference", registry=registry)["can_execute"]
    assert build_batch_plan(inputs, output, engine="reference", registry=registry)["entries"][0]["can_execute"]
    assert run_batch(inputs, output, engine="reference", registry=registry)[0]["status"] == "complete"
    assert run_batch(inputs, output, engine="reference", registry=registry, resume=True)[0]["status"] == "skipped"
    report = benchmark_analysis(EXAMPLE, repeats=2, engine="reference", registry=registry)
    assert len(report["runs"]) == 2
    assert report["median_seconds"] > 0


def test_capabilities_reject_before_execution(tmp_path, registry):
    with pytest.raises(InputValidationError, match="workers=1"):
        analyze_file(EXAMPLE, tmp_path / "out.csv", AnalysisConfig(engine="reference", workers=2), registry=registry)
    assert not (tmp_path / "out.csv").exists()
    plan = build_analysis_plan(EXAMPLE, requested_engine="reference", workers=2, registry=registry)
    assert not plan["can_execute"]
    assert any("workers=1" in warning for warning in plan["warnings"])
    with pytest.raises(InputValidationError, match="h5ad"):
        analyze_file(tmp_path / "missing.h5ad", tmp_path / "out.csv", AnalysisConfig(engine="reference"), registry=registry)
    with pytest.raises(InputValidationError, match="layers"):
        analyze(pd.read_csv(EXAMPLE), AnalysisConfig(engine="reference", layer="counts"), registry=registry)


@pytest.mark.parametrize("column,value", [("p_adj", 1.1), ("p_value", float("nan")),
    ("log2FC", "bad"), ("rank", 0), ("cells_expressing_in", 1.5), ("gene", None)])
def test_invalid_plugin_values(column, value):
    class Invalid(ReferenceEngine):
        def analyze_dataframe(self, expression, config):
            table = pd.DataFrame({"cluster": ["A"], "gene": ["G"]})
            table[column] = value
            return EngineOutput(table, {})
    with pytest.raises(ConfigurationError):
        analyze(pd.DataFrame(), AnalysisConfig(engine="reference"), registry=EngineRegistry((Invalid(),)))


def test_missing_columns_and_invalid_capabilities():
    class Invalid(ReferenceEngine):
        def analyze_dataframe(self, expression, config):
            return EngineOutput(pd.DataFrame({"gene": ["G"]}), {})
    with pytest.raises(ConfigurationError, match="cluster and gene"):
        analyze(pd.DataFrame(), AnalysisConfig(engine="reference"), registry=EngineRegistry((Invalid(),)))
    with pytest.raises(ConfigurationError):
        EngineCapabilities(input_formats=("invalid",))


def test_cli_plugin_workflows(tmp_path, monkeypatch, registry):
    monkeypatch.setattr(cli, "load_engine_plugins", lambda: registry)
    config = tmp_path / "config.json"
    AnalysisConfig(engine="reference").to_json(config)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "sample.csv").write_bytes(EXAMPLE.read_bytes())
    for args in [
        ["analyze", str(EXAMPLE), str(tmp_path / "single.csv")],
        ["plan", str(EXAMPLE), "--output", str(tmp_path / "plan.json")],
        ["batch-plan", str(inputs), str(tmp_path / "batch"), "--output", str(tmp_path / "batch-plan.json")],
        ["batch", str(inputs), str(tmp_path / "batch")],
        ["benchmark", str(EXAMPLE), "--repeats", "1", "--output", str(tmp_path / "benchmark.json")],
    ]:
        assert main(args + ["--config", str(config), "--load-engine-plugins"]) == 0
