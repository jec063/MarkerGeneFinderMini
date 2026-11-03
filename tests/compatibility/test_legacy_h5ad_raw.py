"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
import pytest
from src.marker_finder import build_parser

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


def test_parser_accepts_raw():
    args = build_parser().parse_args([
        "--input", "input.h5ad", "--output", "output.csv", "--use-raw",
    ])
    assert args.use_raw is True
    assert args.layer is None

def test_parser_rejects_raw_with_layer():
    with pytest.raises(SystemExit):
        build_parser().parse_args([
            "--input", "input.h5ad", "--output", "output.csv",
            "--use-raw", "--layer", "alternate",
        ])
