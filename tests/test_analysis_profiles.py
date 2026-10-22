"""Tests for saved analysis profiles."""

import json

from src.analysis_profiles import (
    PROFILE_DEFAULTS,
    build_analysis_profile,
    compare_profile_state,
    is_newer_app_version,
    reset_profile_state,
)
from src.version import APP_VERSION


def test_build_analysis_profile_serializes_portable_state() -> None:
    state = PROFILE_DEFAULTS.copy()
    state.update(
        {
            "top_n": 9,
            "limit_other_mean": True,
            "max_other_mean": 0.25,
            "excluded_gene_text": "MALAT1, GAPDH",
        }
    )

    profile = json.loads(build_analysis_profile(state, "Focused markers"))

    assert profile["schema_version"] == 1
    assert profile["app_version"] == APP_VERSION
    assert profile["profile_name"] == "Focused markers"
    assert profile["top_n"] == 9
    assert profile["max_other_mean"] == 0.25
    assert profile["excluded_genes"] == ["MALAT1", "GAPDH"]
    assert "input_file" not in profile


def test_profile_comparison_only_lists_changes() -> None:
    changes = compare_profile_state(
        {"top_n": 5, "rank_by": "log2FC"},
        {"top_n": 10, "rank_by": "log2FC"},
    )

    assert changes == [{"setting": "top_n", "current": 5, "imported": 10}]


def test_reset_profile_state_restores_independent_defaults() -> None:
    state = {"top_n": 99, "unrelated": "retained"}

    reset_profile_state(state)

    assert state["top_n"] == 5
    assert state["unrelated"] == "retained"


def test_newer_version_detection() -> None:
    assert is_newer_app_version("99.0.0") is True
    assert is_newer_app_version(APP_VERSION) is False
    assert is_newer_app_version("not-a-version") is False
