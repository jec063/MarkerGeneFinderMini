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
