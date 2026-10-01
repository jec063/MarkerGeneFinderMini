import io
import json
import zipfile

import pandas as pd

from src.report_bundle import build_analysis_bundle


def test_bundle_contains_results_settings_and_optional_summary():
    markers = pd.DataFrame({"cluster": ["A"], "gene": ["G1"], "rank": [1]})
    summary = pd.DataFrame({"cluster": ["A"], "gene": ["G1"],
                            "mean_expression": [2.0]})

    contents = build_analysis_bundle(markers, {"engine": "native"}, summary)

    with zipfile.ZipFile(io.BytesIO(contents)) as archive:
        assert set(archive.namelist()) == {
            "marker_gene_results.csv", "analysis_settings.json",
            "expression_summary.csv", "README.md",
        }
        saved_markers = pd.read_csv(archive.open("marker_gene_results.csv"))
        pd.testing.assert_frame_equal(saved_markers, markers)
        assert json.load(archive.open("analysis_settings.json")) == {
            "engine": "native"
        }


def test_bundle_omits_unavailable_expression_summary():
    contents = build_analysis_bundle(pd.DataFrame(), {})

    with zipfile.ZipFile(io.BytesIO(contents)) as archive:
        assert "expression_summary.csv" not in archive.namelist()
