"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
from src.marker_finder import build_parser

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


def test_parser_accepts_layer():
    args = build_parser().parse_args([
        "--input", "expression.h5ad",
        "--output", "markers.csv",
        "--layer", "alternate",
    ])
    assert args.layer == "alternate"
