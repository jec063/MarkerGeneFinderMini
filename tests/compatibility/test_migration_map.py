"""Every shipped legacy entry point has a documented public replacement."""
import importlib
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from marker_gene_finder import AnalysisConfig, analyze_file

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("old,new", [
    ("analysis_profiles", "profiles"), ("expression_summary", "expression_summary"),
    ("filter_presets", "filter_presets"), ("report_bundle", "report_bundle"),
    ("settings_import", "settings"), ("version", "version"),
])
def test_helper_shims_keep_public_names_and_have_documented_replacements(old, new):
    legacy = importlib.import_module("src." + old)
    public = importlib.import_module("marker_gene_finder." + new)
    for name, value in vars(public).items():
        if not name.startswith("_") and (name.isupper() or getattr(value, "__module__", None) == public.__name__):
            assert getattr(legacy, name) is value
    guide = (ROOT / "docs/migration.md").read_text()
    assert "src." + old in guide
    assert "marker_gene_finder." + new in guide


def test_all_flat_options_and_marker_exports_are_documented():
    from src.marker_finder import build_parser
    import src.marker_finder as legacy
    guide = (ROOT / "docs/migration.md").read_text()
    for action in build_parser()._actions:
        for flag in action.option_strings:
            assert f"`{flag}`" in guide, flag
    for name in vars(legacy):
        if not name.startswith("_"):
            assert f"`{name}`" in guide, name


def test_flat_analysis_and_public_api_produce_same_rows(tmp_path):
    source = ROOT / "data/example_expression.csv"
    legacy_output = tmp_path / "legacy.csv"
    subprocess.run([sys.executable, "-m", "src.marker_finder", "--input", str(source),
                    "--output", str(legacy_output), "--top-n", "2", "--quiet"], check=True)
    result = analyze_file(source, tmp_path / "public.csv", AnalysisConfig(top_n=2))
    pd.testing.assert_frame_equal(pd.read_csv(legacy_output), result.markers,
                                  check_dtype=False, rtol=1e-6, atol=1e-10)
