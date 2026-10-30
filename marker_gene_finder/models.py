"""Typed result containers for the public v2 API."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .config import AnalysisConfig


@dataclass(frozen=True, slots=True)
class MarkerResult:
    """Marker table plus stable analysis metadata."""

    markers: pd.DataFrame
    config: AnalysisConfig
    metadata: dict[str, Any]

    def to_csv(self, path: str | Path) -> None:
        """Write marker rows without a DataFrame index."""
        self.markers.to_csv(path, index=False)

    def to_bundle(self, path: str | Path | None = None) -> bytes:
        """Export results, replayable config, settings, and captured provenance."""
        from .report_bundle import build_analysis_bundle
        from .settings import build_settings_document

        settings = dict(self.metadata["provenance"]["effective_config"])
        settings.pop("schema_name", None)
        settings.pop("schema_version", None)
        settings["expression_source"] = self.metadata["provenance"]["expression_source"]
        contents = build_analysis_bundle(
            self.markers, build_settings_document(settings),
            provenance=self.metadata["provenance"], config=self.config.to_dict(),
        )
        if path is not None:
            Path(path).write_bytes(contents)
        return contents
