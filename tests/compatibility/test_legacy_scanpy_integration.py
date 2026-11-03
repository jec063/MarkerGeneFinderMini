"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
from pathlib import Path
import pytest
from src.marker_finder import build_parser

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


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

def test_parser_accepts_expressing_cell_count_ranking() -> None:
    args = build_parser().parse_args(
        [
            "--input", "input.h5ad",
            "--output", "markers.csv",
            "--rank-by", "cells_expressing_in",
        ]
    )

    assert args.rank_by == "cells_expressing_in"

def test_parser_accepts_reference_expressing_cell_count_ranking() -> None:
    args = build_parser().parse_args(
        [
            "--input", "input.h5ad",
            "--output", "markers.csv",
            "--rank-by", "cells_expressing_out",
        ]
    )

    assert args.rank_by == "cells_expressing_out"
