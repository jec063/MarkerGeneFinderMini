"""Tests for the Streamlit interface."""

from pathlib import Path

from streamlit.testing.v1 import AppTest


APP_PATH = Path(__file__).resolve().parents[1] / "app.py"


def test_app_loads_without_errors():
    """The app should load and request an expression-data upload."""
    app = AppTest.from_file(APP_PATH).run()

    assert len(app.exception) == 0
    assert "CellCoPilot Marker Gene Finder" in app.title[0].value
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

    app.checkbox(key="use_example_data").uncheck().run()

    assert len(app.exception) == 0
    assert "marker_results" not in app.session_state
    assert app.info[0].value == "Upload a CSV or H5AD file to begin."
