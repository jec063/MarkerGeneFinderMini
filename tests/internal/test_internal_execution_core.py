import ast
from pathlib import Path
import pandas as pd
import pytest
from marker_gene_finder import analysis, batch, benchmark, inspection
from marker_gene_finder._atomic import atomic_path, write_bytes, write_text
from marker_gene_finder.config import AnalysisConfig
from marker_gene_finder.models import MarkerResult
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures/v2_7"


@pytest.mark.parametrize("existing", [False, True])
@pytest.mark.parametrize("exception", [OSError, KeyboardInterrupt])
def test_interrupted_write_never_publishes_partial_output(tmp_path, existing, exception):
    destination = tmp_path / "markers.csv"
    if existing:
        destination.write_bytes(b"previous complete output")
    with pytest.raises(exception):
        with atomic_path(destination) as staging:
            staging.write_bytes(b"partial")
            raise exception("interrupted")
    if existing:
        assert destination.read_bytes() == b"previous complete output"
    else:
        assert not destination.exists()
    assert list(tmp_path.glob(".*.tmp")) == []


def test_atomic_publish_failure_cleans_up_all_export_types(tmp_path, monkeypatch):
    from marker_gene_finder import _atomic
    def fail(*args):
        raise OSError("replace failed")
    monkeypatch.setattr(_atomic.os, "replace", fail)
    config = AnalysisConfig()
    result = MarkerResult(pd.read_csv(FIXTURES / "native.csv"), config, {})
    for export in [lambda p: write_text(p, "text"), lambda p: write_bytes(p, b"zip"),
                   config.to_json, result.to_csv]:
        destination = tmp_path / "output"
        destination.write_text("old")
        with pytest.raises(OSError, match="replace failed"):
            export(destination)
        assert destination.read_text() == "old"
        assert list(tmp_path.iterdir()) == [destination]


def test_services_own_implementations_and_do_not_depend_on_core():
    for module, function in [(analysis, analysis.run), (batch, batch.run_batch),
                             (benchmark, benchmark.benchmark_analysis),
                             (inspection, inspection.inspect_input)]:
        assert function.__module__ == module.__name__
    package = Path(analysis.__file__).parent
    for path in package.glob("*.py"):
        if path.name == "_core.py":
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom):
                assert node.module not in {"_core", "marker_gene_finder._core"}
                assert all(alias.name != "_core" for alias in node.names)
