"""Tests for importing exported analysis settings."""

import json

import pytest

from src.settings_import import load_analysis_settings


def encode(settings: dict) -> bytes:
    """Encode settings as an uploaded JSON file."""
    return json.dumps(settings).encode("utf-8")


def test_load_analysis_settings_restores_portable_values() -> None:
    restored = load_analysis_settings(
        encode(
            {
                "schema_version": 1,
                "top_n": 12,
                "filter_preset": "Custom",
                "rank_by": "pct_difference",
                "min_cluster_size": 4,
                "min_pct": 0.25,
                "min_cells_expressing": 3,
                "min_cluster_mean": 1.5,
                "min_mean_difference": 0.5,
                "max_other_mean": 0.4,
                "max_cells_expressing_out": 2,
                "max_pct_out": 0.2,
                "min_pct_difference": 0.1,
                "min_log2fc": 1.0,
                "max_p_adj": 0.05,
                "pseudocount": 0.2,
                "excluded_genes": ["MALAT1"],
                "excluded_prefixes": ["MT-", "RPL"],
            }
        )
    )

    assert restored["top_n"] == 12
    assert restored["rank_by"] == "pct_difference"
    assert restored["filter_preset"] == "Custom"
    assert restored["min_cluster_size"] == 4
    assert restored["limit_other_mean"] is True
    assert restored["max_other_mean"] == 0.4
    assert restored["limit_cells_expressing_out"] is True
    assert restored["max_cells_expressing_out"] == 2
    assert restored["excluded_gene_text"] == "MALAT1"
    assert restored["excluded_prefix_text"] == "MT-, RPL"


def test_null_optional_limits_disable_their_controls() -> None:
    restored = load_analysis_settings(
        encode(
            {
                "schema_version": 1,
                "max_other_mean": None,
                "max_cells_expressing_out": None,
            }
        )
    )

    assert restored["limit_other_mean"] is False
    assert restored["limit_cells_expressing_out"] is False


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"schema_version": 2}, "schema version 1"),
        ({"schema_version": 1, "top_n": 0}, "top_n"),
        ({"schema_version": 1, "min_pct": 1.1}, "min_pct"),
        ({"schema_version": 1, "rank_by": "unknown"}, "rank_by"),
        ({"schema_version": 1, "min_cluster_mean": float("nan")}, "finite"),
        ({"schema_version": 1, "excluded_genes": "MALAT1"}, "excluded_genes"),
    ],
)
def test_invalid_analysis_settings_are_rejected(settings: dict, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        load_analysis_settings(encode(settings))


def test_invalid_json_is_rejected() -> None:
    with pytest.raises(ValueError, match="valid UTF-8 JSON"):
        load_analysis_settings(b"not-json")
