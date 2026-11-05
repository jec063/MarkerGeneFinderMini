"""Validate the actual v3 removals rather than simulating blocked imports."""
import importlib.util
from pathlib import Path

import pytest

from marker_gene_finder.cli import main

ROOT = Path(__file__).resolve().parents[1]


def test_legacy_modules_are_physically_absent():
    assert not (ROOT / "src").exists()
    for module in ("src", "marker_gene_finder._core", "marker_gene_finder._legacy_cli"):
        assert importlib.util.find_spec(module) is None, module
    for filename in ("_core.py", "_legacy_cli.py"):
        assert not (ROOT / "marker_gene_finder" / filename).exists()
    assert '"src"' not in (ROOT / "pyproject.toml").read_text()


def test_flat_cli_is_rejected_without_writing_output(tmp_path):
    output = tmp_path / "markers.csv"
    with pytest.raises(SystemExit) as error:
        main(["--input", "expression.csv", "--output", str(output)])
    assert error.value.code == 2
    assert not output.exists()
