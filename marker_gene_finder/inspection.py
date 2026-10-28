"""Input inspection and non-destructive execution planning."""

from ._core import (
    PLAN_SCHEMA_VERSION,
    build_analysis_plan,
    estimate_dense_expression_bytes,
    inspect_input,
)

__all__ = [
    "PLAN_SCHEMA_VERSION",
    "build_analysis_plan",
    "estimate_dense_expression_bytes",
    "inspect_input",
]
