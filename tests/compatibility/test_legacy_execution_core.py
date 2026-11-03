"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
from pathlib import Path
import pandas as pd
from marker_gene_finder import analysis

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


def test_legacy_monkeypatch_reaches_extracted_implementation(tmp_path, monkeypatch):
    from marker_gene_finder import _core
    original = analysis.find_markers
    expected = pd.read_csv(FIXTURES / "native.csv")
    with monkeypatch.context() as patch:
        patch.setattr(_core, "find_markers", lambda *args, **kwargs: expected)
        actual = _core.run(FIXTURES / "expression.csv", tmp_path / "markers.csv")
        assert actual is expected
    assert analysis.find_markers is original
