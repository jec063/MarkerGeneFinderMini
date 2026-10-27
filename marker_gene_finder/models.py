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
