"""Build a portable archive of marker-analysis outputs."""

import io
import json
import zipfile

import pandas as pd


def build_analysis_bundle(
    markers: pd.DataFrame,
    settings: dict,
    expression_summary: pd.DataFrame | None = None,
) -> bytes:
    """Return a ZIP containing results, settings, and an optional summary."""
    buffer = io.BytesIO()
    files = {
        "marker_gene_results.csv": markers.to_csv(index=False),
        "analysis_settings.json": json.dumps(
            settings, indent=2, ensure_ascii=False, allow_nan=False
        ) + "\n",
        "README.md": (
            "# Marker Gene Finder analysis bundle\n\n"
            "This archive contains marker results and the settings used to "
            "produce them. If present, `expression_summary.csv` contains the "
            "genes selected in the expression plots.\n"
        ),
    }
    if expression_summary is not None:
        files["expression_summary.csv"] = expression_summary.to_csv(index=False)

    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, contents in files.items():
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, contents.encode("utf-8"))
    return buffer.getvalue()
