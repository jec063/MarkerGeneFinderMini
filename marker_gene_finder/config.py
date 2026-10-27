"""Typed and serializable v2 analysis configuration."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

from .exceptions import ConfigurationError, SchemaError

CONFIG_SCHEMA_VERSION = 2
RANK_BY_CHOICES = {
    "log2FC",
    "cluster_mean",
    "cells_expressing_in",
    "cells_expressing_out",
    "mean_difference",
    "pct_difference",
    "p_adj",
}


@dataclass(frozen=True, slots=True)
class AnalysisConfig:
    """Portable configuration shared by the v2 API and CLI."""

    cluster_column: str = "cluster"
    cell_column: str | None = "cell"
    engine: str = "native"
    top_n: int = 5
    pseudocount: float = 0.1
    min_pct: float = 0.0
    min_pct_difference: float = 0.0
    min_log2fc: float = 0.0
    max_p_adj: float = 1.0
    max_pct_out: float = 1.0
    rank_by: str = "log2FC"
    min_cluster_size: int = 1
    min_cells_expressing: int = 0
    min_cluster_mean: float = 0.0
    min_mean_difference: float = 0.0
    max_other_mean: float | None = None
    max_cells_expressing_out: int | None = None
    layer: str | None = None
    use_raw: bool = False
    target_clusters: tuple[str, ...] = ()
    excluded_genes: tuple[str, ...] = ()
    excluded_prefixes: tuple[str, ...] = ()
    max_dense_memory_mb: float | None = 2048.0
    gene_chunk_size: int = 1000
    workers: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.cluster_column, str) or not self.cluster_column.strip():
            raise ConfigurationError("cluster_column cannot be empty.")
        if self.cell_column is not None and not isinstance(self.cell_column, str):
            raise ConfigurationError("cell_column must be a string or null.")
        if self.layer is not None and not isinstance(self.layer, str):
            raise ConfigurationError("layer must be a string or null.")
        if self.engine not in {"native", "scanpy", "auto"}:
            raise ConfigurationError("engine must be native, scanpy, or auto.")
        if self.rank_by not in RANK_BY_CHOICES:
            raise ConfigurationError(f"Unsupported rank_by value: {self.rank_by}.")
        for name in ("top_n", "min_cluster_size", "gene_chunk_size", "workers"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ConfigurationError(f"{name} must be a positive integer.")
        if self.min_cells_expressing < 0:
            raise ConfigurationError("min_cells_expressing cannot be negative.")
        if (
            self.max_cells_expressing_out is not None
            and self.max_cells_expressing_out < 0
        ):
            raise ConfigurationError("max_cells_expressing_out cannot be negative.")
        for name in ("min_pct", "min_pct_difference", "max_p_adj", "max_pct_out"):
            value = getattr(self, name)
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ConfigurationError(f"{name} must be between 0 and 1.")
        for name in (
            "pseudocount", "min_log2fc", "min_cluster_mean",
            "min_mean_difference",
        ):
            value = getattr(self, name)
            if not math.isfinite(value):
                raise ConfigurationError(f"{name} must be finite.")
        if self.pseudocount <= 0:
            raise ConfigurationError("pseudocount must be positive.")
        if self.max_other_mean is not None and not math.isfinite(self.max_other_mean):
            raise ConfigurationError("max_other_mean must be finite when provided.")
        if (
            self.max_dense_memory_mb is not None
            and (
                not math.isfinite(self.max_dense_memory_mb)
                or self.max_dense_memory_mb <= 0
            )
        ):
            raise ConfigurationError(
                "max_dense_memory_mb must be positive or null."
            )
        if self.use_raw and self.layer is not None:
            raise ConfigurationError("use_raw and layer are mutually exclusive.")
        for name in ("target_clusters", "excluded_genes", "excluded_prefixes"):
            values = getattr(self, name)
            if not isinstance(values, tuple) or not all(
                isinstance(value, str) for value in values
            ):
                raise ConfigurationError(f"{name} must contain only strings.")

    def to_dict(self) -> dict[str, Any]:
        """Return a versioned JSON-compatible representation."""
        document = asdict(self)
        document["schema_name"] = "marker-gene-finder.analysis-config"
        document["schema_version"] = CONFIG_SCHEMA_VERSION
        return document

    def to_json(self, path: str | Path) -> None:
        """Write this configuration as formatted JSON."""
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    def to_analysis_kwargs(self) -> dict[str, Any]:
        """Return keyword arguments accepted by the internal v1 engine layer."""
        values = asdict(self)
        for name in ("target_clusters", "excluded_genes", "excluded_prefixes"):
            values[name] = list(values[name]) or None
        return values

    @classmethod
    def from_dict(cls, document: dict[str, Any]) -> "AnalysisConfig":
        """Load a v2 configuration document and reject unknown fields."""
        if not isinstance(document, dict):
            raise SchemaError("Analysis configuration must be a JSON object.")
        schema = document.get("schema_version", CONFIG_SCHEMA_VERSION)
        if schema != CONFIG_SCHEMA_VERSION:
            raise SchemaError(
                f"Unsupported analysis configuration schema version: {schema}."
            )
        schema_name = document.get(
            "schema_name", "marker-gene-finder.analysis-config"
        )
        if schema_name != "marker-gene-finder.analysis-config":
            raise SchemaError(f"Unsupported schema name: {schema_name}.")
        allowed = {field.name for field in fields(cls)}
        metadata = {"schema_name", "schema_version", "app_version"}
        unknown = set(document) - allowed - metadata
        if unknown:
            raise SchemaError(
                "Unknown analysis configuration fields: "
                + ", ".join(sorted(unknown))
                + "."
            )
        values = {key: document[key] for key in allowed if key in document}
        for name in ("target_clusters", "excluded_genes", "excluded_prefixes"):
            if name in values:
                if not isinstance(values[name], (list, tuple)):
                    raise SchemaError(f"{name} must be a JSON array.")
                values[name] = tuple(values[name])
        return cls(**values)

    @classmethod
    def from_json(cls, path: str | Path) -> "AnalysisConfig":
        """Load and validate a configuration JSON file."""
        try:
            document = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise SchemaError(f"Cannot load analysis configuration: {error}") from error
        return cls.from_dict(document)
