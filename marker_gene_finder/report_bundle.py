"""Build a portable archive of marker-analysis outputs."""

import io
import json
import zipfile

import pandas as pd

from .settings import build_settings_document, migrate_settings_document


def build_analysis_bundle(
    markers: pd.DataFrame,
    settings: dict,
    expression_summary: pd.DataFrame | None = None,
    cluster_summary: pd.DataFrame | None = None,
    *,
    provenance: dict | None = None,
    config: dict | None = None,
) -> bytes:
    """Return a ZIP containing results, settings, and an optional summary."""
    buffer = io.BytesIO()
    if "schema_version" not in settings:
        settings = build_settings_document(settings)
    else:
        settings = migrate_settings_document(settings)
    provenance = provenance if provenance is not None else settings.get("provenance")
    config = config if config is not None else settings.get("analysis_config")
    app_version = settings.get("app_version", "unknown")
    files = {
        "marker_gene_results.csv": markers.to_csv(index=False),
        "analysis_settings.json": json.dumps(
            settings, indent=2, ensure_ascii=False, allow_nan=False
        ) + "\n",
        "README.md": (
            "# Marker Gene Finder analysis bundle\n\n"
            f"Application version: {app_version}\n\n"
            "This archive contains marker results and the settings used to "
            "produce them. If present, `expression_summary.csv` contains the "
            "genes selected in the expression plots. `cluster_summary.csv` "
            "contains cluster cell counts and dataset percentages.\n"
        ),
    }
    if provenance is not None:
        files["provenance.json"] = json.dumps(provenance, indent=2, allow_nan=False) + "\n"
    if config is not None:
        files["analysis_config.json"] = json.dumps(config, indent=2, allow_nan=False) + "\n"
    if provenance is not None or config is not None:
        files["README.md"] += (
            "\n`provenance.json` records the execution environment and effective settings. "
            "`analysis_config.json` is the requested configuration, reusable with --config. "
            "Input hashes are optional; the expression matrix is not included.\n"
        )
    if expression_summary is not None:
        files["expression_summary.csv"] = expression_summary.to_csv(index=False)
    if cluster_summary is not None:
        files["cluster_summary.csv"] = cluster_summary.to_csv(index=False)

    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, contents in files.items():
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, contents.encode("utf-8"))
    return buffer.getvalue()
