"""Portable, named analysis profiles for the Streamlit controls."""

from __future__ import annotations

import json
from collections.abc import Mapping

from .filter_presets import FILTER_DEFAULTS
from .version import APP_VERSION
from .settings import build_settings_document


PROFILE_DEFAULTS: dict[str, object] = {
    "top_n": 5,
    "rank_by": "log2FC",
    "min_cluster_size": 1,
    "min_cells_expressing": 0,
    "min_cluster_mean": 0.0,
    "min_mean_difference": 0.0,
    "limit_other_mean": False,
    "max_other_mean": 0.0,
    "limit_cells_expressing_out": False,
    "max_cells_expressing_out": 0,
    "pseudocount": 0.1,
    "excluded_gene_text": "",
    "excluded_prefix_text": "",
    "filter_preset": "Permissive",
    **FILTER_DEFAULTS,
}

PROFILE_STATE_KEYS = tuple(PROFILE_DEFAULTS)


def reset_profile_state(state: dict) -> None:
    """Reset portable controls to their stable application defaults."""
    state.update(PROFILE_DEFAULTS)


def build_analysis_profile(state: Mapping[str, object], name: str) -> bytes:
    """Serialize the current portable controls without requiring an analysis."""
    profile_name = name.strip() or "Untitled profile"
    settings = {
        "top_n": int(state["top_n"]),
        "filter_preset": state["filter_preset"],
        "rank_by": state["rank_by"],
        "min_cluster_size": int(state["min_cluster_size"]),
        "min_pct": float(state["min_pct"]),
        "min_cells_expressing": int(state["min_cells_expressing"]),
        "min_cluster_mean": float(state["min_cluster_mean"]),
        "min_mean_difference": float(state["min_mean_difference"]),
        "max_other_mean": (
            float(state["max_other_mean"]) if state["limit_other_mean"] else None
        ),
        "max_cells_expressing_out": (
            int(state["max_cells_expressing_out"])
            if state["limit_cells_expressing_out"] else None
        ),
        "max_pct_out": float(state["max_pct_out"]),
        "min_pct_difference": float(state["min_pct_difference"]),
        "min_log2fc": float(state["min_log2fc"]),
        "max_p_adj": float(state["max_p_adj"]),
        "pseudocount": float(state["pseudocount"]),
        "excluded_genes": _split_values(str(state["excluded_gene_text"])),
        "excluded_prefixes": _split_values(str(state["excluded_prefix_text"])),
    }
    profile = build_settings_document(
        settings, profile_name=profile_name, app_version=APP_VERSION
    )
    return (json.dumps(profile, indent=2, ensure_ascii=False) + "\n").encode()


def compare_profile_state(
    state: Mapping[str, object], restored: Mapping[str, object]
) -> list[dict[str, object]]:
    """Return portable settings whose imported values differ from current state."""
    return [
        {"setting": key, "current": state.get(key), "imported": value}
        for key, value in restored.items()
        if state.get(key) != value
    ]


def is_newer_app_version(imported: str | None) -> bool:
    """Return whether an imported semantic version is newer than this app."""
    if not imported:
        return False
    try:
        imported_parts = tuple(int(part) for part in imported.split("."))
        current_parts = tuple(int(part) for part in APP_VERSION.split("."))
    except ValueError:
        return False
    return imported_parts > current_parts


def _split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]
