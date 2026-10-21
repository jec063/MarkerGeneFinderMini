"""Tests for named marker-filter presets."""

import pytest

from src.filter_presets import FILTER_DEFAULTS, FILTER_PRESETS, get_filter_preset


def test_permissive_preset_matches_backward_compatible_defaults() -> None:
    assert FILTER_DEFAULTS == {
        "min_pct": 0.0,
        "max_pct_out": 1.0,
        "min_log2fc": 0.0,
        "min_pct_difference": 0.0,
        "max_p_adj": 1.0,
    }


def test_presets_become_progressively_more_selective() -> None:
    permissive = FILTER_PRESETS["Permissive"]
    balanced = FILTER_PRESETS["Balanced"]
    specific = FILTER_PRESETS["High specificity"]

    assert permissive["min_pct"] < balanced["min_pct"] < specific["min_pct"]
    assert permissive["max_pct_out"] > balanced["max_pct_out"] > specific["max_pct_out"]
    assert permissive["max_p_adj"] > balanced["max_p_adj"] > specific["max_p_adj"]


def test_get_filter_preset_returns_an_independent_copy() -> None:
    values = get_filter_preset("Balanced")
    values["min_pct"] = 0.9

    assert FILTER_PRESETS["Balanced"]["min_pct"] == 0.25


def test_unknown_filter_preset_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown filter preset"):
        get_filter_preset("Unknown")
