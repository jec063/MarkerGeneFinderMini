"""Reporting must preserve the input files it describes."""
import json
from pathlib import Path

import pytest

from marker_gene_finder.cli import main

FIXTURES = Path(__file__).parent / "fixtures/v2_9"


@pytest.mark.parametrize("command", ["inspect", "plan", "benchmark", "verify"])
def test_report_cannot_replace_its_input(tmp_path, capsys, command):
    source = tmp_path / "input.csv"
    source.write_bytes(b"original input")
    assert main([command, str(source), "--output", str(source)]) == 2
    assert "distinct from source file" in capsys.readouterr().err
    assert source.read_bytes() == b"original input"


@pytest.mark.parametrize("which", ["baseline", "candidate"])
def test_comparison_cannot_replace_either_input(tmp_path, which):
    paths = {name: tmp_path / f"{name}.json" for name in ("baseline", "candidate")}
    for path in paths.values():
        path.write_text("original benchmark")
    assert main(["compare-benchmarks", str(paths["baseline"]), str(paths["candidate"]),
                 "--output", str(paths[which])]) == 2
    assert all(path.read_text() == "original benchmark" for path in paths.values())


@pytest.mark.parametrize("command", ["inspect", "plan", "benchmark"])
def test_report_cannot_replace_configuration(tmp_path, command):
    config = tmp_path / "config.json"
    config.write_text('{"top_n":2}')
    assert main([command, "missing.csv", "--config", str(config), "--output", str(config)]) == 2
    assert config.read_text() == '{"top_n":2}'


def test_report_rejects_symlink_alias(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("original input")
    alias = tmp_path / "alias.json"
    alias.symlink_to(source)
    assert main(["inspect", str(source), "--output", str(alias)]) == 2
    assert alias.is_symlink() and source.read_text() == "original input"


@pytest.mark.parametrize("destination", ["input", "markers", "manifest"])
def test_batch_plan_preserves_data_and_recovery_files(tmp_path, destination):
    inputs, outputs = tmp_path / "inputs", tmp_path / "outputs"
    inputs.mkdir()
    outputs.mkdir()
    paths = {"input": inputs / "sample.csv", "markers": outputs / "sample.csv.markers.csv",
             "manifest": outputs / "batch_manifest.json"}
    paths["input"].write_bytes((FIXTURES / "sample.csv").read_bytes())
    for name in ("markers", "manifest"):
        paths[name].write_text("previous result")
    before = {name: path.read_bytes() for name, path in paths.items()}
    assert main(["batch-plan", str(inputs), str(outputs), "--output", str(paths[destination])]) == 2
    assert {name: path.read_bytes() for name, path in paths.items()} == before


def test_comparing_a_report_with_itself_remains_supported(tmp_path):
    source = FIXTURES / "benchmark.json"
    output = tmp_path / "comparison.json"
    assert main(["compare-benchmarks", str(source), str(source), "--output", str(output)]) == 0
    assert json.loads(output.read_text())["passed"]
