"""Tests for the installable project metadata."""

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_pyproject_defines_version_and_cli_entry_point() -> None:
    contents = (PROJECT_ROOT / "pyproject.toml").read_text(encoding="utf-8")

    assert 'version = {attr = "src.version.APP_VERSION"}' in contents
    assert 'marker-gene-finder = "src.marker_finder:main"' in contents
    assert 'requires-python = ">=3.10"' in contents
