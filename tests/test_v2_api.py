"""Contract tests for the stable v2 package and subcommand CLI."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from marker_gene_finder import (
    APP_VERSION,
    AnalysisConfig,
    ConfigurationError,
    InputValidationError,
    MarkerResult,
    SchemaError,
    analyze,
)
from marker_gene_finder.cli import main


def expression_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell": ["a", "b", "c", "d"],
            "cluster": ["A", "A", "B", "B"],
            "gene_a": [8.0, 7.0, 0.0, 0.0],
            "gene_b": [0.0, 0.0, 8.0, 7.0],
        }
    )


def test_v2_config_round_trip_and_public_result(tmp_path) -> None:
    config = AnalysisConfig(top_n=1, target_clusters=("A",))
    path = tmp_path / "analysis-config.json"
    config.to_json(path)

    restored = AnalysisConfig.from_json(path)
    result = analyze(expression_table(), restored)

    assert APP_VERSION == "2.3.1"
    assert restored == config
    assert isinstance(result, MarkerResult)
    assert result.markers["cluster"].tolist() == ["A"]
    assert result.metadata["app_version"] == APP_VERSION
    assert json.loads(path.read_text())["schema_version"] == 2


def test_v2_config_rejects_invalid_and_unknown_values() -> None:
    with pytest.raises(ConfigurationError, match="engine"):
        AnalysisConfig(engine="Not Valid")
    with pytest.raises(ConfigurationError, match="cluster_column"):
        AnalysisConfig(cluster_column=1)  # type: ignore[arg-type]
    with pytest.raises(ConfigurationError, match="only strings"):
        AnalysisConfig(excluded_genes=(1,))  # type: ignore[arg-type]
    with pytest.raises(SchemaError, match="Unknown"):
        AnalysisConfig.from_dict({"schema_version": 2, "surprise": True})
    with pytest.raises(SchemaError, match="schema version"):
        AnalysisConfig.from_dict({"schema_version": 99})


def test_public_api_wraps_input_validation_errors() -> None:
    with pytest.raises(InputValidationError, match="cluster"):
        analyze(expression_table().drop(columns="cluster"))


def test_v2_cli_analyze_inspect_and_plan(tmp_path, capsys) -> None:
    input_path = tmp_path / "expression.csv"
    output_path = tmp_path / "markers.csv"
    inspection_path = tmp_path / "inspection.json"
    plan_path = tmp_path / "plan.json"
    expression_table().to_csv(input_path, index=False)

    assert main(["analyze", str(input_path), str(output_path)]) == 0
    assert output_path.is_file()
    assert "Wrote" in capsys.readouterr().out

    assert main(["inspect", str(input_path), "--output", str(inspection_path)]) == 0
    assert json.loads(inspection_path.read_text())["genes"] == 2

    assert main(["plan", str(input_path), "--output", str(plan_path)]) == 0
    plan = json.loads(plan_path.read_text())
    assert plan["can_execute"] is True
    assert plan["resolved_engine"] == "native"


def test_v2_cli_reports_bad_config_without_traceback(tmp_path, capsys) -> None:
    config_path = tmp_path / "bad.json"
    config_path.write_text('{"schema_version": 99}')

    status = main(
        ["inspect", "missing.csv", "--output", str(tmp_path / "out.json"),
         "--config", str(config_path)]
    )

    assert status == 2
    assert "Unsupported analysis configuration schema version" in capsys.readouterr().err
