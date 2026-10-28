"""Regression tests for the v2.1 internal package boundary."""

from __future__ import annotations

from pathlib import Path

import marker_gene_finder._core as core
import src.marker_finder as legacy


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


def test_legacy_module_aliases_v2_core_for_monkeypatch_compatibility() -> None:
    assert legacy is core
    assert legacy.find_markers is core.find_markers
    assert legacy.run_batch is core.run_batch


def test_service_modules_are_owned_by_public_package() -> None:
    from marker_gene_finder import analysis, batch, benchmark, inspection
    from marker_gene_finder import expression_summary, filter_presets
    from marker_gene_finder import profiles, report_bundle, settings

    assert analysis.run is core.run
    assert batch.run_batch is core.run_batch
    assert benchmark.benchmark_analysis is core.benchmark_analysis
    assert inspection.inspect_input is core.inspect_input
    assert callable(expression_summary.summarize_expression)
    assert filter_presets.FILTER_PRESETS
    assert callable(profiles.build_analysis_profile)
    assert callable(report_bundle.build_analysis_bundle)
    assert settings.SETTINGS_SCHEMA_VERSION == 2
