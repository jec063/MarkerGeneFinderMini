"""Read and resume artifacts captured from the unmodified v2.9 release."""
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from marker_gene_finder import AnalysisConfig, SchemaError, analyze_file
from marker_gene_finder.batch import run_batch, verify_batch_manifest
from marker_gene_finder.benchmark import compare_benchmarks
from marker_gene_finder.settings import load_analysis_profile_upload, load_analysis_settings_upload

FIXTURES = Path(__file__).resolve().parent / "fixtures/v2_9"


def test_v29_config_replays_same_markers(tmp_path):
    config = AnalysisConfig.from_json(FIXTURES / "analysis_config.json")
    result = analyze_file(FIXTURES / "sample.csv", tmp_path / "markers.csv", config)
    pd.testing.assert_frame_equal(result.markers,
                                  pd.read_csv(FIXTURES / "sample.csv.markers.csv"),
                                  check_dtype=False, rtol=1e-6, atol=1e-10)


def test_v29_settings_profiles_and_benchmarks_remain_readable(tmp_path):
    restored = load_analysis_settings_upload(
        (FIXTURES / "analysis_settings.json").read_bytes(), "settings.json")
    assert restored["top_n"] == 2
    state, metadata = load_analysis_profile_upload(
        (FIXTURES / "profile.json").read_bytes(), "profile.json")
    assert state["top_n"] == 2 and metadata["app_version"] == "2.9.0"
    assert compare_benchmarks(FIXTURES / "benchmark.json", FIXTURES / "benchmark.json")["passed"]


def test_v29_batch_manifest_verifies_and_resumes(tmp_path, monkeypatch):
    manifest = json.loads((FIXTURES / "batch_manifest.json").read_text())
    assert manifest["app_version"] == "2.9.0"
    (tmp_path / "inputs").mkdir()
    (tmp_path / "outputs").mkdir()
    source = tmp_path / "inputs/sample.csv"
    source.write_bytes((FIXTURES / "sample.csv").read_bytes())
    mtime = manifest["runs"][0]["input_fingerprint"]["mtime_ns"]
    os.utime(source, ns=(mtime, mtime))
    for name in ("batch_manifest.json", "sample.csv.markers.csv"):
        (tmp_path / "outputs" / name).write_bytes((FIXTURES / name).read_bytes())
    monkeypatch.chdir(tmp_path)
    assert verify_batch_manifest("outputs/batch_manifest.json")["valid"]
    assert run_batch("inputs", "outputs", top_n=2, resume=True)[0]["status"] == "skipped"


def test_future_artifact_schemas_still_receive_precise_errors(tmp_path):
    with pytest.raises(SchemaError, match="Unsupported analysis configuration schema version"):
        AnalysisConfig.from_dict({"schema_version": 99})
    with pytest.raises(ValueError, match="supported versions are 1 and 2"):
        load_analysis_settings_upload(b'{"schema_version":99}', "settings.json")
    path = tmp_path / "manifest.json"
    path.write_text('{"schema_version":99}')
    with pytest.raises(ValueError, match="Unsupported batch manifest schema version"):
        verify_batch_manifest(path)
