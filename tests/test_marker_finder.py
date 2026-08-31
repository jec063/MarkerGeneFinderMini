from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.marker_finder import (
    _benjamini_hochberg,
    build_parser,
    find_markers,
    load_expression,
    run,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_DATA = (
    PROJECT_ROOT / "data" / "example_expression.csv"
)


def test_load_expression_reads_h5ad(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "expression.h5ad"
    adata = ad.AnnData(
        X=np.array(
            [
                [5.0, 0.0],
                [4.0, 0.1],
                [0.0, 5.0],
                [0.2, 4.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": ["A", "A", "B", "B"]},
            index=["cell_1", "cell_2", "cell_3", "cell_4"],
        ),
        var=pd.DataFrame(
            index=["GENE_A", "GENE_B"],
        ),
    )
    adata.write_h5ad(input_path)

    expression = load_expression(
        input_path,
        cluster_column="cluster",
        cell_column="cell_id",
    )

    expected = pd.DataFrame(
        {
            "cell_id": [
                "cell_1",
                "cell_2",
                "cell_3",
                "cell_4",
            ],
            "cluster": ["A", "A", "B", "B"],
            "GENE_A": [5.0, 4.0, 0.0, 0.2],
            "GENE_B": [0.0, 0.1, 5.0, 4.0],
        }
    )

    pd.testing.assert_frame_equal(
        expression,
        expected,
        check_dtype=False,
    )


def test_run_reads_sparse_h5ad(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "sparse_expression.h5ad"
    output_path = tmp_path / "markers.csv"

    adata = ad.AnnData(
        X=sparse.csr_matrix(
            [
                [5.0, 0.0],
                [4.0, 0.1],
                [0.0, 5.0],
                [0.2, 4.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": ["A", "A", "B", "B"]},
            index=["cell_1", "cell_2", "cell_3", "cell_4"],
        ),
        var=pd.DataFrame(
            index=["GENE_A", "GENE_B"],
        ),
    )
    adata.write_h5ad(input_path)

    markers = run(
        input_path,
        output_path,
        top_n=1,
    )
    saved = pd.read_csv(output_path)

    observed = dict(
        zip(markers["cluster"], markers["gene"])
    )
    assert observed == {
        "A": "GENE_A",
        "B": "GENE_B",
    }
    pd.testing.assert_frame_equal(saved, markers)


def test_h5ad_missing_cluster_column_is_rejected(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "missing_cluster.h5ad"
    adata = ad.AnnData(
        X=np.array(
            [
                [1.0],
                [2.0],
            ]
        ),
        obs=pd.DataFrame(
            index=["cell_1", "cell_2"],
        ),
        var=pd.DataFrame(
            index=["GENE_A"],
        ),
    )
    adata.write_h5ad(input_path)

    with pytest.raises(
        ValueError,
        match="Missing cluster column in AnnData.obs",
    ):
        load_expression(input_path)


def test_unsupported_input_format_is_rejected(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "expression.txt"

    with pytest.raises(
        ValueError,
        match="Unsupported input format",
    ):
        load_expression(input_path)


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


def test_run_applies_max_p_adj(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "significant_markers.csv"

    markers = run(
        EXAMPLE_DATA,
        output_path,
        top_n=100,
        max_p_adj=0.02,
    )

    saved = pd.read_csv(output_path)

    assert not markers.empty
    assert (markers["p_adj"] <= 0.02).all()
    pd.testing.assert_frame_equal(saved, markers)


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


def test_max_p_adj_filters_nonsignificant_genes() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=100,
        max_p_adj=0.02,
    )

    assert not markers.empty
    assert (markers["p_adj"] <= 0.02).all()


def test_parser_accepts_max_p_adj() -> None:
    args = build_parser().parse_args(
        [
            "--input",
            "input.csv",
            "--output",
            "output.csv",
            "--max-p-adj",
            "0.05",
        ]
    )

    assert args.max_p_adj == pytest.approx(0.05)


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


def test_invalid_max_p_adj_is_rejected() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    for invalid_value in (-0.01, 1.01):
        with pytest.raises(
            ValueError,
            match="max_p_adj must be between 0 and 1",
        ):
            find_markers(
                expression,
                max_p_adj=invalid_value,
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
