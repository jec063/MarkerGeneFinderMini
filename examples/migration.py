"""Executable migration examples using only the supported public interfaces.

Run: python examples/migration.py /tmp/marker-gene-finder-migration
"""
from dataclasses import replace
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from marker_gene_finder import (
    AnalysisConfig, EngineCapabilities, EngineRegistry, NativeEngine,
    analyze, analyze_file, migrate_settings_document,
)
from marker_gene_finder.cli import main
from marker_gene_finder.profiles import PROFILE_DEFAULTS, build_analysis_profile
from marker_gene_finder.settings import load_analysis_profile_upload


class ReferenceEngine(NativeEngine):
    name = "reference"
    capabilities = EngineCapabilities()

    def analyze_file(self, input_path, output_path, config):
        return NativeEngine().analyze_file(
            input_path, output_path, replace(config, engine="native")
        )


def run_demo(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    inputs = directory / "inputs"
    inputs.mkdir(exist_ok=True)
    expression = pd.DataFrame({
        "cluster": ["A"] * 4 + ["B"] * 4,
        "GENE_A": [5., 4., 6., 5., 0., 0., 0., 0.],
        "GENE_B": [0., 0., 0., 0., 5., 4., 6., 5.],
    })
    source = inputs / "expression.csv"
    expression.to_csv(source, index=False)
    config = AnalysisConfig(top_n=1)
    config_path = directory / "config.json"
    config.to_json(config_path)
    assert AnalysisConfig.from_json(config_path) == config

    in_memory = analyze(expression, config)
    result = analyze_file(source, directory / "markers.csv", config)
    pd.testing.assert_frame_equal(result.markers, in_memory.markers)
    result.to_bundle(directory / "bundle.zip")

    # Source choice is configuration, independent of the API/CLI entry point.
    matrix = np.log1p(expression.drop(columns="cluster").to_numpy())
    data = ad.AnnData(
        sparse.csr_matrix(matrix),
        obs=pd.DataFrame({"cluster": pd.Categorical(expression.cluster)},
                         index=[str(i) for i in range(len(expression))]),
        var=pd.DataFrame(index=["GENE_A", "GENE_B"]),
    )
    data.layers["lognormalized"] = data.X.copy()
    data.raw = data.copy()
    h5ad = inputs / "expression.h5ad"
    data.write_h5ad(h5ad)
    for engine in ("native", "scanpy"):
        for source_name, options in (
            ("X", {}), ("layer", {"layer": "lognormalized"}), ("raw", {"use_raw": True}),
        ):
            selected = replace(config, engine=engine, **options)
            selected_result = analyze_file(
                h5ad, directory / f"{engine}-{source_name}.csv", selected
            )
            assert selected_result.markers.gene.tolist() == ["GENE_A", "GENE_B"]

    plugin = analyze_file(source, directory / "plugin.csv",
                          replace(config, engine="reference"),
                          registry=EngineRegistry((ReferenceEngine(),)))
    pd.testing.assert_frame_equal(plugin.markers, result.markers)

    # Settings/profile schemas are independent of application versions.
    migrated = migrate_settings_document({"schema_version": 1, "top_n": 1})
    assert migrated["settings"]["top_n"] == 1
    profile = build_analysis_profile({**PROFILE_DEFAULTS, "top_n": 1}, "Migration")
    restored, metadata = load_analysis_profile_upload(profile, "profile.json")
    assert restored["top_n"] == 1 and metadata["profile_name"] == "Migration"

    commands = [
        ["list-engines", "--output", str(directory / "engines.json")],
        ["analyze", str(source), str(directory / "cli.csv"), "--config", str(config_path),
         "--bundle", str(directory / "cli.zip"), "--metadata", str(directory / "metadata.json"), "--hash-input"],
        ["inspect", str(h5ad), "--output", str(directory / "inspection.json")],
        ["plan", str(h5ad), "--output", str(directory / "plan.json")],
        ["batch-plan", str(inputs), str(directory / "batch"), "--output", str(directory / "batch-plan.json")],
        ["batch", str(inputs), str(directory / "batch")],
        ["batch", str(inputs), str(directory / "batch"), "--resume"],
        ["verify", str(directory / "batch/batch_manifest.json"), "--output", str(directory / "verification.json")],
        ["benchmark", str(source), "--output", str(directory / "benchmark.json"), "--repeats", "1"],
        ["compare-benchmarks", str(directory / "benchmark.json"), str(directory / "benchmark.json"), "--output", str(directory / "comparison.json")],
    ]
    for command in commands:
        assert main(command) == 0, command
    manifest = json.loads((directory / "batch/batch_manifest.json").read_text())
    assert manifest["skipped"] == 2


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("directory", type=Path)
    run_demo(parser.parse_args().directory)
    print("Migration examples passed.")
