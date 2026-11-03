"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
from __future__ import annotations
import subprocess
import sys

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


def test_legacy_import_emits_deprecation_notice() -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-W",
            "always::DeprecationWarning",
            "-c",
            "import src.filter_presets",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert "src.filter_presets is deprecated" in completed.stderr
    assert "removed in Marker Gene Finder 3.0" in completed.stderr
