"""Versioned, JSON-safe provenance shared by API, CLI, and web exports."""
from __future__ import annotations

import hashlib
import platform
from importlib import metadata
from pathlib import Path

from .config import AnalysisConfig
from .version import APP_VERSION

PROVENANCE_SCHEMA_VERSION = 1


def dependency_versions() -> dict[str, str | None]:
    """Report installed scientific dependencies without requiring extra packages."""
    result: dict[str, str | None] = {"python": platform.python_version()}
    for name in ("numpy", "pandas", "scipy", "anndata", "scanpy"):
        try:
            result[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            result[name] = None
    return result


def hash_input(path: str | Path) -> str:
    """Hash a file in bounded chunks, only when explicitly requested."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_provenance(
    config: AnalysisConfig, *, resolved_engine: str, input_kind: str,
    engine=None, input_sha256: str | None = None,
) -> dict:
    """Record effective settings separately from the replayable requested config.

    Unknown plugin versions remain null; they are never inferred from the app
    version. Plugins can expose a string ``version`` attribute.
    """
    versions = dependency_versions()
    engine_version = getattr(engine, "version", None)
    if not isinstance(engine_version, str):
        engine_version = None
    if engine is None or type(engine).__module__ == "marker_gene_finder.engines":
        if resolved_engine == "native":
            engine_version = APP_VERSION
        elif resolved_engine == "scanpy":
            engine_version = versions["scanpy"]
    effective = config.to_dict()
    effective["engine"] = resolved_engine
    if input_kind == "h5ad":
        effective["cell_column"] = None
    if resolved_engine == "scanpy":
        for key in ("pseudocount", "gene_chunk_size", "max_dense_memory_mb"):
            effective[key] = None
    source = (
        "AnnData.raw.X" if config.use_raw else
        "AnnData.layers" if config.layer is not None else "AnnData.X"
    ) if input_kind == "h5ad" else "CSV" if input_kind == "csv" else "dataframe"
    return {
        "schema_name": "marker-gene-finder.provenance",
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "requested_engine": config.engine,
        "resolved_engine": resolved_engine,
        "engine_version": engine_version,
        "effective_config": effective,
        "expression_source": source,
        "layer": config.layer,
        "dependency_versions": versions,
        "input_sha256": input_sha256,
    }
