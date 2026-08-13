from pathlib import Path

import pandas as pd
import pytest

from src.marker_finder import find_markers, run


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