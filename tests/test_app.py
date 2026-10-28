"""Tests for the Streamlit interface."""

from pathlib import Path

from streamlit.testing.v1 import AppTest

from src.version import APP_VERSION


APP_PATH = Path(__file__).resolve().parents[1] / "app.py"


def test_app_loads_without_errors():
    """The app should load and request an expression-data upload."""
    app = AppTest.from_file(APP_PATH).run()

    assert len(app.exception) == 0
    assert "CellCoPilot Marker Gene Finder" in app.title[0].value
    assert any(f"Version {APP_VERSION}" in item.value for item in app.caption)
    assert app.info[0].value == "Upload a CSV or H5AD file to begin."

def test_example_data_flow():
    """Example data can be analyzed, and switching inputs clears results."""
    app = AppTest.from_file(APP_PATH).run()
    app.checkbox(key="use_example_data").check().run()

    assert len(app.exception) == 0
    assert len(app.error) == 0
    assert any("Using the bundled example CSV" in item.value for item in app.info)

    button = next(
        item for item in app.button if item.label == "Find marker genes"
    )
    button.click().run(timeout=30)

    assert len(app.exception) == 0
    assert len(app.error) == 0
    markers = app.session_state["marker_results"]
    assert not markers.empty
    assert markers["cluster"].nunique() == 4
    assert {"CD3D", "MS4A1", "LST1", "PECAM1"}.issubset(set(markers["gene"]))

    settings = app.session_state["analysis_settings"]
    assert settings["schema_version"] == 2
    assert settings["settings"]["input_file"] == "example_expression.csv"
    assert settings["app_version"] == APP_VERSION
    assert settings["settings"]["example_data"] is True
    assert settings["settings"]["engine"] == "native"
    assert settings["settings"]["expression_source"] == "CSV"
    assert settings["settings"]["layer"] is None
    assert settings["settings"]["use_raw"] is False
    assert settings["settings"]["cluster_column"] == "cluster"
    assert settings["settings"]["pseudocount"] == 0.1
    assert settings["settings"]["clusters"] == 4
    assert settings["settings"]["min_cluster_size"] == 1
    assert settings["settings"]["min_cells_expressing"] == 0
    assert settings["settings"]["min_cluster_mean"] == 0.0
    assert settings["settings"]["min_mean_difference"] == 0.0
    assert settings["settings"]["max_other_mean"] is None
    assert settings["settings"]["max_cells_expressing_out"] is None
    assert settings["settings"]["marker_rows"] == len(markers)

    top_n_widget = next(
        item for item in app.number_input if item.label == "Markers per cluster"
    )
    assert settings["settings"]["top_n"] == int(top_n_widget.value)
    assert settings["settings"]["filter_preset"] == "Permissive"
    top_n_widget.set_value(1).run()
    assert len(app.exception) == 0
    assert "marker_results" not in app.session_state
    assert "analysis_settings" not in app.session_state

    next(
        item for item in app.button if item.label == "Find marker genes"
    ).click().run(timeout=30)
    assert len(app.exception) == 0
    assert len(app.error) == 0
    assert app.session_state["analysis_settings"]["settings"]["top_n"] == 1

    app.checkbox(key="use_example_data").uncheck().run()

    assert len(app.exception) == 0
    assert "marker_results" not in app.session_state
    assert "analysis_settings" not in app.session_state
    assert app.info[0].value == "Upload a CSV or H5AD file to begin."


def test_filter_preset_updates_core_controls() -> None:
    """Selecting a preset updates its filters before analysis."""
    app = AppTest.from_file(APP_PATH).run()
    app.checkbox(key="use_example_data").check().run()

    preset = next(item for item in app.selectbox if item.label == "Filter preset")
    preset.select("Balanced").run()

    assert app.session_state["filter_preset"] == "Balanced"
    assert app.session_state["min_pct"] == 0.25
    assert app.session_state["max_pct_out"] == 0.5
    assert app.session_state["min_log2fc"] == 0.5
    assert app.session_state["min_pct_difference"] == 0.1
    assert app.session_state["max_p_adj"] == 0.05


def test_settings_import_controls_are_available() -> None:
    """The app should expose settings upload and application controls."""
    app = AppTest.from_file(APP_PATH).run()
    app.checkbox(key="use_example_data").check().run()

    assert any(
        item.label == "Import analysis settings" for item in app.file_uploader
    )
    apply_button = next(
        item for item in app.button if item.label == "Apply imported settings"
    )
    assert apply_button.disabled is True
    assert any(
        item.label == "Download current profile" for item in app.get("download_button")
    )
    assert any(
        item.label == "Reset profile to defaults" for item in app.button
    )
