from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.marker_finder import (
    _benjamini_hochberg,
    build_parser,
    find_markers,
    run,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = (
    PROJECT_ROOT / "data" / "example_expression.csv"
)


def test_expected_markers_are_ranked_first() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=2,
    )

    expected = {
        0: {"CD3D", "CD3E"},
        1: {"MS4A1", "CD79A"},
        2: {"LST1", "S100A8"},
        3: {"PECAM1", "VWF"},
    }

    observed = {
        cluster: set(group["gene"])
        for cluster, group in markers.groupby("cluster")
    }

    assert observed == expected


def test_output_contains_expression_prevalence() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=2,
    )

    assert markers["pct_in"].between(0, 1).all()
    assert markers["pct_out"].between(0, 1).all()
    assert (markers["log2FC"] > 0).all()

def test_output_contains_statistical_significance() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=2,
    )

    assert markers["p_value"].between(0, 1).all()
    assert markers["p_adj"].between(0, 1).all()
    assert (
        markers["p_adj"] >= markers["p_value"]
    ).all()


def test_benjamini_hochberg_adjustment() -> None:
    p_values = np.array([0.01, 0.04, 0.03])

    adjusted = _benjamini_hochberg(p_values)

    np.testing.assert_allclose(
        adjusted,
        np.array([0.03, 0.04, 0.04]),
    )


def test_default_filter_excludes_negative_log2fc() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=100,
    )

    assert not markers.empty
    assert (markers["log2FC"] >= 0).all()

def test_min_pct_filters_rare_genes() -> None:
    expression = pd.DataFrame(
        {
            "cell": ["a", "b", "c", "d"],
            "cluster": [0, 0, 1, 1],
            "rare": [5.0, 0.0, 0.0, 0.0],
            "common": [2.0, 2.0, 0.0, 0.0],
        }
    )

    markers = find_markers(
        expression,
        top_n=2,
        min_pct=0.75,
    )

    cluster_zero_genes = set(
        markers.loc[
            markers["cluster"] == 0,
            "gene",
        ]
    )

    assert "common" in cluster_zero_genes
    assert "rare" not in cluster_zero_genes


def test_run_writes_output_file(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "markers.csv"

    markers = run(
        EXAMPLE_DATA,
        output_path,
        top_n=2,
    )

    saved = pd.read_csv(output_path)

    pd.testing.assert_frame_equal(
        saved,
        markers,
    )

def test_run_applies_min_log2fc(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "filtered_markers.csv"

    markers = run(
        EXAMPLE_DATA,
        output_path,
        top_n=100,
        min_log2fc=5.3,
    )

    assert not markers.empty
    assert (markers["log2FC"] >= 5.3).all()

def test_parser_accepts_min_log2fc() -> None:
    args = build_parser().parse_args(
        [
            "--input",
            "input.csv",
            "--output",
            "output.csv",
            "--min-log2fc",
            "0.5",
        ]
    )

    assert args.min_log2fc == pytest.approx(0.5)

def test_negative_min_log2fc_is_rejected() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    with pytest.raises(
        ValueError,
        match="min_log2fc must be at least 0",
    ):
        find_markers(
            expression,
            min_log2fc=-0.1,
        )

def test_single_cluster_is_rejected() -> None:
    expression = pd.DataFrame(
        {
            "cell": ["a", "b"],
            "cluster": [0, 0],
            "CD3D": [1.0, 2.0],
        }
    )

    with pytest.raises(
        ValueError,
        match="At least two clusters",
    ):
        find_markers(expression)
