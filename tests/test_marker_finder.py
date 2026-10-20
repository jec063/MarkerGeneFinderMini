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
    assert markers.groupby("cluster")["rank"].apply(list).to_dict() == {
        cluster: [1, 2] for cluster in expected
    }


def test_rank_reflects_filtered_marker_order() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=3,
        min_pct=0.5,
    )

    for _, cluster_markers in markers.groupby("cluster", sort=False):
        assert cluster_markers["rank"].tolist() == list(
            range(1, len(cluster_markers) + 1)
        )
        assert cluster_markers["log2FC"].is_monotonic_decreasing


def test_output_contains_expression_prevalence() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression,
        top_n=2,
    )

    assert markers["pct_in"].between(0, 1).all()
    assert markers["pct_out"].between(0, 1).all()
    assert (markers["cells_expressing_in"] >= 0).all()
    assert (markers["cells_expressing_out"] >= 0).all()
    assert pd.api.types.is_integer_dtype(markers["cells_expressing_in"])
    assert pd.api.types.is_integer_dtype(markers["cells_expressing_out"])
    assert (markers["log2FC"] > 0).all()


def test_output_contains_comparison_group_sizes() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "B", "B", "B"],
            "gene": [5.0, 4.0, 0.0, 0.0, 0.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=1,
    )

    assert markers["cluster_size"].tolist() == [2]
    assert markers["reference_size"].tolist() == [3]
    assert pd.api.types.is_integer_dtype(markers["cluster_size"])
    assert pd.api.types.is_integer_dtype(markers["reference_size"])


def test_min_cells_expressing_filters_rare_markers() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A"] * 4 + ["B"] * 4,
            "rare": [10, 0, 0, 0, 0, 0, 0, 0],
            "supported": [2, 2, 2, 0, 0, 0, 0, 0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        min_cells_expressing=2,
    )

    assert markers["gene"].tolist() == ["supported"]
    assert markers["cells_expressing_in"].tolist() == [3]


@pytest.mark.parametrize("value", [-1, 1.5, True])
def test_invalid_min_cells_expressing_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="min_cells_expressing"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            min_cells_expressing=value,
        )


def test_parser_accepts_min_cells_expressing() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cells-expressing", "3",
        ]
    )

    assert args.min_cells_expressing == 3


def test_min_cluster_mean_filters_low_expression_markers() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A"] * 4 + ["B"] * 4,
            "low_mean": [0.2, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0, 0.0],
            "high_mean": [2.0, 2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        min_cluster_mean=1.0,
    )

    assert markers["gene"].tolist() == ["high_mean"]
    assert (markers["cluster_mean"] >= 1.0).all()


def test_negative_min_cluster_mean_is_rejected() -> None:
    with pytest.raises(ValueError, match="min_cluster_mean"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            min_cluster_mean=-0.1,
        )


def test_min_mean_difference_filters_weakly_enriched_markers() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "B", "B"],
            "weak": [1.5, 1.5, 1.0, 1.0],
            "strong": [3.0, 3.0, 1.0, 1.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        min_mean_difference=1.0,
    )

    assert markers["gene"].tolist() == ["strong"]
    assert markers["mean_difference"].tolist() == pytest.approx([2.0])


def test_negative_min_mean_difference_is_rejected() -> None:
    with pytest.raises(ValueError, match="min_mean_difference"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            min_mean_difference=-0.1,
        )


def test_max_other_mean_filters_broad_reference_expression() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "B", "B"],
            "specific": [3.0, 3.0, 0.2, 0.2],
            "broad": [5.0, 5.0, 2.0, 2.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        max_other_mean=0.5,
    )

    assert markers["gene"].tolist() == ["specific"]
    assert markers["other_mean"].tolist() == pytest.approx([0.2])


def test_negative_max_other_mean_is_rejected() -> None:
    with pytest.raises(ValueError, match="max_other_mean"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            max_other_mean=-0.1,
        )


def test_max_reference_expressing_cells_filters_broad_markers() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A"] * 4 + ["B"] * 4,
            "specific": [3.0] * 4 + [1.0, 0.0, 0.0, 0.0],
            "broad": [5.0] * 4 + [1.0, 1.0, 1.0, 0.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        max_cells_expressing_out=1,
    )

    assert markers["gene"].tolist() == ["specific"]
    assert markers["cells_expressing_out"].tolist() == [1]


@pytest.mark.parametrize("value", [-1, 1.5, True])
def test_invalid_max_reference_expressing_cells_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="max_cells_expressing_out"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            max_cells_expressing_out=value,
        )


def test_parser_accepts_min_cluster_mean() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cluster-mean", "1.5",
        ]
    )

    assert args.min_cluster_mean == pytest.approx(1.5)


def test_parser_accepts_min_mean_difference() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-mean-difference", "0.75",
        ]
    )

    assert args.min_mean_difference == pytest.approx(0.75)


def test_parser_accepts_max_other_mean() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--max-other-mean", "0.5",
        ]
    )

    assert args.max_other_mean == pytest.approx(0.5)


def test_parser_accepts_max_reference_expressing_cells() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--max-cells-expressing-out", "3",
        ]
    )

    assert args.max_cells_expressing_out == 3

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


def test_min_pct_difference_filters_nonspecific_genes() -> None:
    expression = pd.DataFrame({
        "cluster": ["A", "A", "B", "B"],
        "specific": [2.0, 2.0, 0.0, 0.0],
        "shared": [2.0, 2.0, 1.0, 1.0],
    })

    markers = find_markers(expression, cell_column=None, top_n=10,
                           min_pct_difference=0.5)

    assert set(markers.loc[markers["cluster"] == "A", "gene"]) == {"specific"}
    assert (markers["pct_difference"] >= 0.5).all()


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_invalid_min_pct_difference_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="min_pct_difference"):
        find_markers(pd.read_csv(EXAMPLE_DATA), min_pct_difference=value)


def test_max_pct_out_filters_broadly_expressed_genes() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "A", "A", "B", "B", "B", "B"],
            "specific": [5, 4, 5, 4, 0, 0, 0, 0],
            "broad": [5, 5, 4, 4, 1, 1, 1, 0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=10,
        max_pct_out=0.5,
    )

    assert markers["gene"].tolist() == ["specific"]
    assert (markers["pct_out"] <= 0.5).all()


@pytest.mark.parametrize("value", [-0.01, 1.01])
def test_invalid_max_pct_out_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="max_pct_out"):
        find_markers(pd.read_csv(EXAMPLE_DATA), max_pct_out=value)


def test_rank_by_pct_difference_changes_marker_order() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A"] * 4 + ["B"] * 4,
            "high_fc": [10, 0, 0, 0, 0, 0, 0, 0],
            "specific": [2, 2, 2, 2, 1, 0, 0, 0],
        }
    )

    fold_change = find_markers(
        expression, cell_column=None, target_clusters=["A"], top_n=2,
    )
    specificity = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=2,
        rank_by="pct_difference",
    )

    assert fold_change["gene"].tolist() == ["high_fc", "specific"]
    assert specificity["gene"].tolist() == ["specific", "high_fc"]
    assert specificity["rank"].tolist() == [1, 2]


def test_rank_by_mean_difference_changes_marker_order() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "B", "B"],
            "high_ratio": [2.0, 2.0, 0.0, 0.0],
            "large_difference": [10.0, 10.0, 7.0, 7.0],
        }
    )

    fold_change = find_markers(
        expression, cell_column=None, target_clusters=["A"], top_n=2,
    )
    mean_difference = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=2,
        rank_by="mean_difference",
    )

    assert fold_change["gene"].tolist() == ["high_ratio", "large_difference"]
    assert mean_difference["gene"].tolist() == [
        "large_difference", "high_ratio",
    ]
    assert mean_difference["mean_difference"].tolist() == [3.0, 2.0]
    assert mean_difference["rank"].tolist() == [1, 2]


def test_rank_by_cluster_mean_changes_marker_order() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["A", "A", "B", "B"],
            "high_ratio": [2.0, 2.0, 0.0, 0.0],
            "high_mean": [10.0, 10.0, 7.0, 7.0],
        }
    )

    fold_change = find_markers(
        expression, cell_column=None, target_clusters=["A"], top_n=2,
    )
    cluster_mean = find_markers(
        expression,
        cell_column=None,
        target_clusters=["A"],
        top_n=2,
        rank_by="cluster_mean",
    )

    assert fold_change["gene"].tolist() == ["high_ratio", "high_mean"]
    assert cluster_mean["gene"].tolist() == ["high_mean", "high_ratio"]
    assert cluster_mean["cluster_mean"].tolist() == [10.0, 2.0]
    assert cluster_mean["rank"].tolist() == [1, 2]


def test_invalid_rank_by_is_rejected() -> None:
    with pytest.raises(ValueError, match="rank_by must be one of"):
        find_markers(pd.read_csv(EXAMPLE_DATA), rank_by="unknown")


def test_rank_by_cells_expressing_in_changes_marker_order() -> None:
    expression = pd.DataFrame(
        {
            "cell": ["a1", "a2", "a3", "b1", "b2", "b3"],
            "cluster": ["A", "A", "A", "B", "B", "B"],
            "broad": [1.0, 1.0, 1.0, 0.8, 0.8, 0.8],
            "sparse": [8.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        }
    )

    markers = find_markers(
        expression,
        target_clusters=["A"],
        top_n=2,
        rank_by="cells_expressing_in",
    )

    assert markers["gene"].tolist() == ["broad", "sparse"]
    assert markers["cells_expressing_in"].tolist() == [3, 1]
    assert markers["rank"].tolist() == [1, 2]


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


def test_target_clusters_limit_output_but_compare_with_rest() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(expression, top_n=2, target_clusters=[1])

    assert markers["cluster"].unique().tolist() == [1]
    assert set(markers["gene"]) == {"MS4A1", "CD79A"}


def test_unknown_target_cluster_is_rejected() -> None:
    with pytest.raises(ValueError, match="Unknown target clusters"):
        find_markers(pd.read_csv(EXAMPLE_DATA), target_clusters=["missing"])


def test_min_cluster_size_skips_small_target_clusters() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["small", "large", "large", "large"],
            "small_gene": [5.0, 0.0, 0.0, 0.0],
            "large_gene": [0.0, 5.0, 4.0, 5.0],
        }
    )

    markers = find_markers(
        expression,
        cell_column=None,
        min_cluster_size=2,
    )

    assert markers["cluster"].unique().tolist() == ["large"]


def test_selected_cluster_below_minimum_size_is_rejected() -> None:
    expression = pd.DataFrame(
        {
            "cluster": ["small", "large", "large"],
            "gene": [5.0, 0.0, 0.0],
        }
    )

    with pytest.raises(ValueError, match="below the minimum cluster size"):
        find_markers(
            expression,
            cell_column=None,
            target_clusters=["small"],
            min_cluster_size=2,
        )


@pytest.mark.parametrize("value", [0, -1, 1.5, True])
def test_invalid_min_cluster_size_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="min_cluster_size"):
        find_markers(
            pd.read_csv(EXAMPLE_DATA),
            min_cluster_size=value,
        )


def test_parser_accepts_min_cluster_size() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cluster-size", "10",
        ]
    )

    assert args.min_cluster_size == 10


def test_gene_and_prefix_exclusions_are_applied() -> None:
    expression = pd.read_csv(EXAMPLE_DATA)

    markers = find_markers(
        expression, top_n=100,
        excluded_genes=["CD3D"], excluded_prefixes=["MS4"],
    )

    assert "CD3D" not in set(markers["gene"])
    assert not markers["gene"].str.startswith("MS4").any()


def test_excluding_every_gene_is_rejected() -> None:
    expression = pd.DataFrame({"cluster": ["A", "B"], "G": [1, 0]})
    with pytest.raises(ValueError, match="removed every gene"):
        find_markers(expression, cell_column=None, excluded_genes=["G"])
