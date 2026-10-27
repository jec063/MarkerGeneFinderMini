"""Named marker-filter presets shared by the Streamlit interface and tests."""

FILTER_PRESETS = {
    "Permissive": {
        "min_pct": 0.0,
        "max_pct_out": 1.0,
        "min_log2fc": 0.0,
        "min_pct_difference": 0.0,
        "max_p_adj": 1.0,
    },
    "Balanced": {
        "min_pct": 0.25,
        "max_pct_out": 0.5,
        "min_log2fc": 0.5,
        "min_pct_difference": 0.1,
        "max_p_adj": 0.05,
    },
    "High specificity": {
        "min_pct": 0.5,
        "max_pct_out": 0.25,
        "min_log2fc": 1.0,
        "min_pct_difference": 0.25,
        "max_p_adj": 0.01,
    },
}

FILTER_DEFAULTS = FILTER_PRESETS["Permissive"]


def get_filter_preset(name: str) -> dict[str, float]:
    """Return an independent copy of a named preset's filter values."""
    try:
        return FILTER_PRESETS[name].copy()
    except KeyError as error:
        choices = ", ".join(FILTER_PRESETS)
        raise ValueError(f"Unknown filter preset {name!r}; choose: {choices}.") from error
