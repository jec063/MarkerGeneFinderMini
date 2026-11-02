"""Repeatable analysis benchmarking and comparison."""
from __future__ import annotations

import json
import platform
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
import numpy as np
from .version import APP_VERSION
from ._constants import BENCHMARK_COMPARISON_SCHEMA_VERSION, BENCHMARK_SCHEMA_VERSION
from .analysis import _registered_run


def benchmark_analysis(
    input_path: str | Path,
    repeats: int = 3,
    registry=None,
    **analysis_options,
) -> dict[str, object]:
    """Benchmark repeated analyses and return machine-readable timings."""
    if isinstance(repeats, bool) or not isinstance(repeats, int) or repeats < 1:
        raise ValueError("repeats must be a positive integer.")
    runs: list[dict[str, object]] = []
    with TemporaryDirectory() as temporary_directory:
        for index in range(repeats):
            diagnostics: dict[str, object] = {}
            output = Path(temporary_directory) / f"run-{index + 1}.csv"
            _registered_run(input_path, output, registry=registry, diagnostics=diagnostics, **analysis_options)
            runs.append({"repeat": index + 1, **diagnostics})
    totals = [float(item["total_seconds"]) for item in runs]
    return {
        "schema_version": BENCHMARK_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "input": str(input_path),
        "repeats": repeats,
        "min_seconds": min(totals),
        "median_seconds": median(totals),
        "max_seconds": max(totals),
        "runs": runs,
    }


def compare_benchmarks(
    baseline_path: str | Path,
    candidate_path: str | Path,
    max_regression_percent: float = 10.0,
) -> dict[str, object]:
    """Compare benchmark medians against an allowed regression threshold."""
    if (
        isinstance(max_regression_percent, bool)
        or not np.isfinite(max_regression_percent)
        or max_regression_percent < 0
    ):
        raise ValueError("max_regression_percent must be finite and non-negative.")

    def load(path: str | Path) -> tuple[Path, dict[str, object]]:
        normalized = Path(path)
        report = json.loads(normalized.read_text(encoding="utf-8"))
        if report.get("schema_version") != BENCHMARK_SCHEMA_VERSION:
            raise ValueError(f"Unsupported benchmark schema: {normalized}.")
        value = report.get("median_seconds")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Benchmark median is missing: {normalized}.")
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"Benchmark median must be positive: {normalized}.")
        return normalized, report

    baseline_path, baseline = load(baseline_path)
    candidate_path, candidate = load(candidate_path)
    baseline_median = float(baseline["median_seconds"])
    candidate_median = float(candidate["median_seconds"])
    ratio = candidate_median / baseline_median
    change_percent = (ratio - 1.0) * 100.0
    warnings = [
        f"{field}_mismatch"
        for field in ("platform", "python_version", "input")
        if baseline.get(field) != candidate.get(field)
    ]
    regression = change_percent > float(max_regression_percent)
    return {
        "schema_version": BENCHMARK_COMPARISON_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "baseline": str(baseline_path),
        "candidate": str(candidate_path),
        "baseline_median_seconds": baseline_median,
        "candidate_median_seconds": candidate_median,
        "ratio": ratio,
        "change_percent": change_percent,
        "max_regression_percent": float(max_regression_percent),
        "regression": regression,
        "passed": not regression,
        "compatibility_warnings": warnings,
    }

__all__ = [
    "BENCHMARK_COMPARISON_SCHEMA_VERSION",
    "BENCHMARK_SCHEMA_VERSION",
    "benchmark_analysis",
    "compare_benchmarks",
]
