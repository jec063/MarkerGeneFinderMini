"""Malformed portable settings fail with public, field-specific errors."""
import json

import pytest

from marker_gene_finder import AnalysisConfig, ConfigurationError, SchemaError
from marker_gene_finder.cli import main


@pytest.mark.parametrize("field", [
    "min_pct", "min_pct_difference", "max_p_adj", "max_pct_out",
    "pseudocount", "min_log2fc", "min_cluster_mean", "min_mean_difference",
])
@pytest.mark.parametrize("value", [None, True, "bad", [], float("nan"), float("inf")])
def test_numeric_fields_reject_non_numbers_and_nonfinite_values(field, value):
    with pytest.raises(ConfigurationError, match=field):
        AnalysisConfig.from_dict({field: value})


@pytest.mark.parametrize("field", ["max_other_mean", "max_dense_memory_mb"])
@pytest.mark.parametrize("value", [True, "bad", {}, float("nan")])
def test_optional_numbers_still_validate_present_values(field, value):
    with pytest.raises(ConfigurationError, match=field):
        AnalysisConfig(**{field: value})


@pytest.mark.parametrize("field", ["min_cells_expressing", "max_cells_expressing_out"])
@pytest.mark.parametrize("value", [True, 1.5, "1", [], -1])
def test_counts_require_nonnegative_integers(field, value):
    with pytest.raises(ConfigurationError, match=field):
        AnalysisConfig.from_dict({field: value})


@pytest.mark.parametrize("field,value", [
    ("rank_by", []), ("rank_by", {}), ("use_raw", "false"), ("use_raw", 1),
    ("use_raw", None), ("min_cells_expressing", None), ("pseudocount", 10 ** 400),
])
def test_bad_types_do_not_escape_as_python_errors(field, value):
    with pytest.raises(ConfigurationError, match=field):
        AnalysisConfig.from_dict({field: value})


def test_valid_values_keep_their_meanings_and_round_trip():
    config = AnalysisConfig(
        min_pct=1, max_pct_out=0, pseudocount=1, min_log2fc=-2,
        min_cluster_mean=-1, min_mean_difference=-1, max_other_mean=None,
        max_dense_memory_mb=None, min_cells_expressing=0, max_cells_expressing_out=0,
        use_raw=False,
    )
    assert AnalysisConfig.from_dict(json.loads(json.dumps(config.to_dict()))) == config


@pytest.mark.parametrize("schema", [True, 2.0, "2", None])
def test_schema_version_requires_an_integer(schema):
    with pytest.raises(SchemaError, match="schema version"):
        AnalysisConfig.from_dict({"schema_version": schema})


def test_nonstring_python_dictionary_keys_are_schema_errors():
    with pytest.raises(SchemaError, match="field names"):
        AnalysisConfig.from_dict({1: "bad", "unknown": True})


def test_invalid_utf8_configuration_is_a_schema_error(tmp_path):
    path = tmp_path / "bad.json"
    path.write_bytes(b"\xff")
    with pytest.raises(SchemaError, match="Cannot load analysis configuration"):
        AnalysisConfig.from_json(path)


@pytest.mark.parametrize("document,field", [
    ({"min_pct": "bad"}, "min_pct"), ({"rank_by": []}, "rank_by"),
    ({"use_raw": "false"}, "use_raw"), ({"min_cells_expressing": 1.5}, "min_cells_expressing"),
])
def test_cli_reports_invalid_types_without_touching_output(tmp_path, capsys, document, field):
    config = tmp_path / "config.json"
    config.write_text(json.dumps(document))
    output = tmp_path / "report.json"
    output.write_text("previous report")
    assert main(["inspect", "missing.csv", "--config", str(config), "--output", str(output)]) == 2
    error = capsys.readouterr().err
    assert field in error and "Traceback" not in error
    assert output.read_text() == "previous report"
