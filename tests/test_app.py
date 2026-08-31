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
