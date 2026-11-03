"""Regression tests for the v2.1 internal package boundary."""

from __future__ import annotations

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_public_package_does_not_import_legacy_src_modules() -> None:
    package = PROJECT_ROOT / "marker_gene_finder"
    offenders = {
        path.relative_to(PROJECT_ROOT): line
        for path in package.glob("*.py")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith("from src") or line.startswith("import src")
    }

    assert offenders == {}

    app_source = (PROJECT_ROOT / "app.py").read_text(encoding="utf-8")
    assert "from src" not in app_source
    assert "import src" not in app_source
