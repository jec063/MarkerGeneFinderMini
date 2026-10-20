from pathlib import Path
import pytest
import anndata as ad
import numpy as np
import pandas as pd

from src.marker_finder import (
    RESULT_COLUMNS,
    build_parser,
    find_markers_scanpy,
    run,
)

def test_find_markers_scanpy_ranks_expected_genes() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [5.0, 0.0],
                [4.0, 0.0],
                [6.0, 0.0],
                [5.0, 0.0],
                [0.0, 5.0],
                [0.0, 4.0],
                [0.0, 6.0],
                [0.0, 5.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "A", "A", "B", "B", "B", "B"]
                )
            },
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["GENE_A", "GENE_B"]),
    )

    markers = find_markers_scanpy(
        adata,
        cluster_column="cluster",
        top_n=1,
    )

    observed = dict(
        zip(markers["cluster"], markers["gene"])
    )

    assert markers.columns.tolist() == RESULT_COLUMNS
    assert markers.groupby("cluster")["rank"].apply(list).to_dict() == {
        "A": [1],
        "B": [1],
    }
    assert observed == {
        "A": "GENE_A",
        "B": "GENE_B",
    }
    assert markers["cluster_size"].tolist() == [4, 4]
    assert markers["reference_size"].tolist() == [4, 4]

def test_parser_accepts_scanpy_engine() -> None:
    args = build_parser().parse_args(
        [
            "--input",
            "expression.h5ad",
            "--output",
            "markers.csv",
            "--engine",
            "scanpy",
        ]
    )

    assert args.engine == "scanpy"


def test_scanpy_max_pct_out_filters_broadly_expressed_genes() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [5, 5], [4, 5], [5, 4], [4, 4],
                [0, 1], [0, 1], [0, 1], [0, 0],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["specific", "broad"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=10,
        max_pct_out=0.5,
    )

    assert markers["gene"].tolist() == ["specific"]


def test_parser_accepts_max_pct_out() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--max-pct-out", "0.25",
        ]
    )

    assert args.max_pct_out == pytest.approx(0.25)


def test_parser_accepts_rank_by() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--rank-by", "pct_difference",
        ]
    )

    assert args.rank_by == "pct_difference"


def test_parser_accepts_mean_difference_ranking() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--rank-by", "mean_difference",
        ]
    )

    assert args.rank_by == "mean_difference"


def test_parser_accepts_cluster_mean_ranking() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--rank-by", "cluster_mean",
        ]
    )

    assert args.rank_by == "cluster_mean"


def test_scanpy_min_cluster_size_skips_small_targets() -> None:
    adata = ad.AnnData(
        X=np.array([[5.0, 0.0], [0.0, 4.0], [0.0, 5.0]]),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["small", "large", "large"])},
            index=["cell_1", "cell_2", "cell_3"],
        ),
        var=pd.DataFrame(index=["SMALL_GENE", "LARGE_GENE"]),
    )

    markers = find_markers_scanpy(
        adata,
        top_n=1,
        min_cluster_size=2,
    )

    assert markers["cluster"].unique().tolist() == ["large"]


def test_scanpy_min_cells_expressing_filters_rare_markers() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [10.0, 2.0], [0.0, 2.0], [0.0, 2.0], [0.0, 0.0],
                [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["rare", "supported"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=10,
        min_cells_expressing=2,
    )

    assert markers["gene"].tolist() == ["supported"]
    assert markers["cells_expressing_in"].tolist() == [3]


def test_scanpy_min_cluster_mean_filters_low_expression_markers() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [0.2, 2.0], [0.2, 2.0], [0.2, 2.0], [0.2, 2.0],
                [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["low_mean", "high_mean"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=10,
        min_cluster_mean=1.0,
    )

    assert markers["gene"].tolist() == ["high_mean"]
    assert (markers["cluster_mean"] >= 1.0).all()


def test_scanpy_min_mean_difference_filters_weak_markers() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [1.5, 3.0], [1.5, 3.0], [1.5, 3.0], [1.5, 3.0],
                [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["weak", "strong"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=10,
        min_mean_difference=1.0,
    )

    assert markers["gene"].tolist() == ["strong"]
    assert markers["mean_difference"].tolist() == pytest.approx([2.0])


def test_scanpy_can_rank_by_mean_difference() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [2.0, 10.0], [2.0, 10.0], [2.0, 10.0], [2.0, 10.0],
                [0.0, 7.0], [0.0, 7.0], [0.0, 7.0], [0.0, 7.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["high_ratio", "large_difference"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=2,
        rank_by="mean_difference",
    )

    assert markers["gene"].tolist() == ["large_difference", "high_ratio"]
    assert markers["mean_difference"].tolist() == pytest.approx([3.0, 2.0])
    assert markers["rank"].tolist() == [1, 2]


def test_scanpy_can_rank_by_cluster_mean() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [2.0, 10.0], [2.0, 10.0], [2.0, 10.0], [2.0, 10.0],
                [0.0, 7.0], [0.0, 7.0], [0.0, 7.0], [0.0, 7.0],
            ]
        ),
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["high_ratio", "high_mean"]),
    )

    markers = find_markers_scanpy(
        adata,
        target_clusters=["A"],
        top_n=2,
        rank_by="cluster_mean",
    )

    assert markers["gene"].tolist() == ["high_mean", "high_ratio"]
    assert markers["cluster_mean"].tolist() == pytest.approx([10.0, 2.0])
    assert markers["rank"].tolist() == [1, 2]


def test_run_uses_scanpy_engine(tmp_path: Path) -> None:
    input_path = tmp_path / "expression.h5ad"
    output_path = tmp_path / "markers.csv"

    adata = ad.AnnData(
        X=np.array(
            [
                [5.0, 0.0],
                [4.0, 0.0],
                [0.0, 5.0],
                [0.0, 4.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "B", "B"]
                )
            },
            index=["cell_1", "cell_2", "cell_3", "cell_4"],
        ),
        var=pd.DataFrame(index=["GENE_A", "GENE_B"]),
    )
    adata.write_h5ad(input_path)

    markers = run(
        input_path,
        output_path,
        cluster_column="cluster",
        top_n=1,
        engine="scanpy",
    )

    observed = dict(
        zip(markers["cluster"], markers["gene"])
    )

    assert observed == {
        "A": "GENE_A",
        "B": "GENE_B",
    }
    assert output_path.exists()


def test_scanpy_engine_rejects_csv(tmp_path: Path) -> None:
    input_path = tmp_path / "expression.csv"
    output_path = tmp_path / "markers.csv"

    pd.DataFrame(
        {
            "cell": ["cell_1", "cell_2"],
            "cluster": ["A", "B"],
            "GENE_A": [5.0, 0.0],
        }
    ).to_csv(input_path, index=False)

    with pytest.raises(
        ValueError,
        match=r"Scanpy engine requires a \.h5ad input file",
    ):
        run(
            input_path,
            output_path,
            engine="scanpy",
        )
