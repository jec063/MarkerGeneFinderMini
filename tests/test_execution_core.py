"""v2.7 scientific fixtures and v2.8 execution guarantees."""
import ast
from pathlib import Path
import threading

import anndata as ad
import pandas as pd
import pytest
from scipy import sparse

from marker_gene_finder import analysis, batch, benchmark, inspection
from marker_gene_finder._atomic import atomic_path, write_bytes, write_text
from marker_gene_finder.config import AnalysisConfig
from marker_gene_finder.models import MarkerResult

FIXTURES = Path(__file__).parent / "fixtures" / "v2_7"


@pytest.mark.parametrize("engine", ["native", "scanpy"])
@pytest.mark.parametrize("source", ["X", "layer", "raw"])
def test_engine_matches_its_v27_fixture(tmp_path, engine, source):
    frame = pd.read_csv(FIXTURES / "expression.csv")
    data = ad.AnnData(sparse.csr_matrix(frame.drop(columns="cluster").to_numpy()),
                      obs=pd.DataFrame({"cluster": pd.Categorical(frame.cluster)},
                                       index=frame.index.astype(str)),
                      var=pd.DataFrame(index=frame.columns[1:]))
    options = {}
    if source == "layer":
        data.layers["counts"] = data.X.copy()
        data.X = data.X * 0
        options["layer"] = "counts"
    elif source == "raw":
        data.raw = data.copy()
        data.X = data.X * 0
        options["use_raw"] = True
    input_path = tmp_path / "input.h5ad"
    data.write_h5ad(input_path)
    actual = analysis.run(input_path, tmp_path / "markers.csv", engine=engine,
                          top_n=8, **options)
    expected = pd.read_csv(FIXTURES / f"{engine}.csv")
    pd.testing.assert_frame_equal(actual, expected, check_dtype=False,
                                  rtol=1e-6, atol=1e-10)


@pytest.mark.parametrize("workers", [1, 2])
def test_csv_fixture_and_ordered_progress(tmp_path, workers):
    events = []
    caller = threading.get_ident()
    actual = analysis.run(FIXTURES / "expression.csv", tmp_path / "markers.csv",
                          workers=workers, gene_chunk_size=3, top_n=8,
                          progress_callback=lambda *event: events.append(
                              (*event, threading.get_ident())))
    pd.testing.assert_frame_equal(actual, pd.read_csv(FIXTURES / "native.csv"),
                                  check_dtype=False, rtol=1e-6, atol=1e-10)
    assert events == [(1, 3, "A", caller), (2, 3, "B", caller), (3, 3, "C", caller)]


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


def test_failed_csv_serialization_preserves_existing_output(tmp_path, monkeypatch):
    destination = tmp_path / "markers.csv"
    destination.write_text("old")
    def interrupted(self, path, **kwargs):
        Path(path).write_text("partial")
        raise KeyboardInterrupt
    monkeypatch.setattr(pd.DataFrame, "to_csv", interrupted)
    with pytest.raises(KeyboardInterrupt):
        analysis.run(FIXTURES / "expression.csv", destination)
    assert destination.read_text() == "old"
    assert list(tmp_path.iterdir()) == [destination]


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


def test_invalid_plugin_output_does_not_replace_previous_file(tmp_path):
    from marker_gene_finder import (
        analyze_file, EngineRegistry, NativeEngine, EngineOutput, ConfigurationError,
    )
    class InvalidEngine(NativeEngine):
        name = "invalid"

        def analyze_file(self, input_path, output_path, config):
            Path(output_path).write_text("partial plugin output")
            return EngineOutput(pd.DataFrame({"gene": ["G"]}), {})

    destination = tmp_path / "markers.csv"
    destination.write_text("previous results")
    with pytest.raises(ConfigurationError, match="cluster and gene"):
        analyze_file(FIXTURES / "expression.csv", destination,
                     AnalysisConfig(engine="invalid"),
                     registry=EngineRegistry((InvalidEngine(),)))
    assert destination.read_text() == "previous results"
    assert list(tmp_path.iterdir()) == [destination]


def test_interrupted_batch_preserves_previous_manifest_and_output(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "sample.csv").write_bytes((FIXTURES / "expression.csv").read_bytes())
    outputs = tmp_path / "outputs"
    batch.run_batch(inputs, outputs)
    before = {path.name: path.read_bytes() for path in outputs.iterdir()}

    def interrupted(input_path, output_path, **kwargs):
        Path(output_path).write_text("partial")
        raise KeyboardInterrupt

    monkeypatch.setattr(batch, "_registered_run", interrupted)
    with pytest.raises(KeyboardInterrupt):
        batch.run_batch(inputs, outputs)
    assert {path.name: path.read_bytes() for path in outputs.iterdir()} == before


def test_bundle_replace_failure_preserves_previous_bundle(tmp_path, monkeypatch):
    from marker_gene_finder import analyze, _atomic
    result = analyze(pd.read_csv(FIXTURES / "expression.csv"))
    destination = tmp_path / "bundle.zip"
    destination.write_bytes(b"previous bundle")

    def interrupted(*args):
        raise OSError("replace failed")

    monkeypatch.setattr(_atomic.os, "replace", interrupted)
    with pytest.raises(OSError, match="replace failed"):
        result.to_bundle(destination)
    assert destination.read_bytes() == b"previous bundle"
    assert list(tmp_path.iterdir()) == [destination]


def test_legacy_monkeypatch_reaches_extracted_implementation(tmp_path, monkeypatch):
    from marker_gene_finder import _core
    original = analysis.find_markers
    expected = pd.read_csv(FIXTURES / "native.csv")
    with monkeypatch.context() as patch:
        patch.setattr(_core, "find_markers", lambda *args, **kwargs: expected)
        actual = _core.run(FIXTURES / "expression.csv", tmp_path / "markers.csv")
        assert actual is expected
    assert analysis.find_markers is original
