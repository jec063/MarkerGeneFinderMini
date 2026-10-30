"""Reproducibility exports preserve configuration and actual execution context."""
import hashlib
import io
import json
from pathlib import Path
from zipfile import ZipFile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from marker_gene_finder import AnalysisConfig, analyze, analyze_file
from marker_gene_finder.cli import main
from marker_gene_finder.provenance import build_provenance
from marker_gene_finder.settings import load_analysis_settings_upload

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "data/example_expression.csv"


def test_api_bundle_round_trip_and_hash(tmp_path):
    config = AnalysisConfig(engine="auto", top_n=2)
    result = analyze_file(EXAMPLE, tmp_path / "markers.csv", config, hash_input_file=True)
    provenance = result.metadata["provenance"]
    assert provenance["requested_engine"] == "auto"
    assert provenance["resolved_engine"] == "native"
    assert provenance["engine_version"] == "2.7.0"
    assert provenance["input_sha256"] == hashlib.sha256(EXAMPLE.read_bytes()).hexdigest()
    assert provenance["expression_source"] == "CSV"
    assert all(provenance["dependency_versions"].values())
    contents = result.to_bundle(tmp_path / "result.zip")
    assert contents == (tmp_path / "result.zip").read_bytes()
    assert contents == result.to_bundle()
    with ZipFile(io.BytesIO(contents)) as archive:
        restored = AnalysisConfig.from_dict(json.loads(archive.read("analysis_config.json")))
        assert restored == config
        assert json.loads(archive.read("provenance.json"))["resolved_engine"] == "native"
        pd.testing.assert_frame_equal(pd.read_csv(archive.open("marker_gene_results.csv")),
                                      pd.read_csv(tmp_path / "markers.csv"))
    assert load_analysis_settings_upload(contents, "result.zip")["top_n"] == 2


def test_hashing_is_opt_in(tmp_path, monkeypatch):
    from marker_gene_finder import api
    def forbidden(path):
        raise AssertionError("Unexpected input hashing")
    monkeypatch.setattr(api, "hash_input", forbidden)
    result = analyze_file(EXAMPLE, tmp_path / "markers.csv")
    assert result.metadata["provenance"]["input_sha256"] is None


def test_cli_and_api_provenance_agree(tmp_path):
    api_result = analyze_file(EXAMPLE, tmp_path / "api.csv", hash_input_file=True)
    assert main(["analyze", str(EXAMPLE), str(tmp_path / "cli.csv"),
                 "--bundle", str(tmp_path / "cli.zip"),
                 "--metadata", str(tmp_path / "metadata.json"), "--hash-input"]) == 0
    saved = json.loads((tmp_path / "metadata.json").read_text())
    assert saved["provenance"] == json.loads(json.dumps(api_result.metadata["provenance"]))
    with ZipFile(tmp_path / "cli.zip") as archive:
        assert json.loads(archive.read("provenance.json")) == saved["provenance"]


def test_cli_rejects_colliding_paths_before_analysis(tmp_path):
    source = tmp_path / "input.csv"
    source.write_bytes(EXAMPLE.read_bytes())
    before = source.read_bytes()
    assert main(["analyze", str(source), str(tmp_path / "out.csv"),
                 "--bundle", str(source)]) == 2
    assert source.read_bytes() == before
    assert not (tmp_path / "out.csv").exists()


def test_web_and_api_effective_settings_agree():
    app = AppTest.from_file(ROOT / "app.py").run()
    app.checkbox(key="use_example_data").check().run()
    next(button for button in app.button if button.label == "Find marker genes").click().run(timeout=30)
    assert not app.exception
    saved = app.session_state["analysis_settings"]
    config = AnalysisConfig.from_dict(saved["analysis_config"])
    result = analyze(pd.read_csv(EXAMPLE), config)
    assert saved["provenance"]["effective_config"] == result.metadata["provenance"]["effective_config"]
    assert saved["provenance"]["resolved_engine"] == "native"
    assert saved["provenance"]["expression_source"] == "CSV"


@pytest.mark.parametrize("use_raw,layer,source", [
    (False, None, "AnnData.X"), (True, None, "AnnData.raw.X"),
    (False, "counts", "AnnData.layers"),
])
def test_h5ad_expression_sources(tmp_path, use_raw, layer, source):
    data = ad.AnnData(np.array([[3., 0.], [4., 0.], [0., 4.], [0., 3.]]),
                      obs=pd.DataFrame({"cluster": ["a", "a", "b", "b"]}, index=list("abcd")))
    data.raw = data.copy()
    data.layers["counts"] = data.X.copy()
    path = tmp_path / "input.h5ad"
    data.write_h5ad(path)
    result = analyze_file(path, tmp_path / "markers.csv", AnalysisConfig(use_raw=use_raw, layer=layer))
    assert result.metadata["provenance"]["expression_source"] == source
    assert result.metadata["provenance"]["layer"] == layer
    assert result.metadata["provenance"]["effective_config"]["cell_column"] is None


def test_plugin_versions_and_scanpy_unused_settings():
    class Plugin:
        version = "1.2.3"
    config = AnalysisConfig(engine="plugin")
    assert build_provenance(config, resolved_engine="plugin", input_kind="csv", engine=Plugin())["engine_version"] == "1.2.3"
    assert build_provenance(config, resolved_engine="plugin", input_kind="csv", engine=object())["engine_version"] is None
    p = build_provenance(AnalysisConfig(engine="scanpy"), resolved_engine="scanpy", input_kind="h5ad")
    assert p["engine_version"] == p["dependency_versions"]["scanpy"]
    assert p["effective_config"]["pseudocount"] is None


def test_dataframe_cluster_resolution_preserves_labels_and_input():
    table = pd.read_csv(EXAMPLE)
    original = table.copy(deep=True)
    result = analyze(table, AnalysisConfig(target_clusters=("0",)))
    assert set(result.markers["cluster"]) == {0}
    pd.testing.assert_frame_equal(table, original)


def test_dataframe_ambiguous_cluster_labels_are_rejected():
    from marker_gene_finder import InputValidationError
    table = pd.DataFrame({"cluster": [0, "0"], "gene": [2.0, 1.0]})
    with pytest.raises(InputValidationError, match="ambiguous"):
        analyze(table, AnalysisConfig(target_clusters=("0",)))
