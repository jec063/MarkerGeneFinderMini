"""Validate portable analysis settings for the Streamlit interface."""

from __future__ import annotations

import json
import math
from io import BytesIO
from pathlib import Path
from zipfile import BadZipFile, ZipFile

from .filter_presets import FILTER_PRESETS
from .analysis import RANK_BY_CHOICES


MAX_SETTINGS_BYTES = 1_000_000
SETTINGS_BUNDLE_PATH = "analysis_settings.json"
SETTINGS_SCHEMA_NAME = "marker-gene-finder-analysis-settings"
SETTINGS_SCHEMA_VERSION = 2


def build_settings_document(
    settings: dict[str, object],
    *,
    profile_name: str | None = None,
    app_version: str | None = None,
) -> dict[str, object]:
    """Wrap settings in the stable v2 interchange envelope."""
    from .version import APP_VERSION

    values = dict(settings)
    embedded_version = values.pop("app_version", None)
    document: dict[str, object] = {
        "schema_name": SETTINGS_SCHEMA_NAME,
        "schema_version": SETTINGS_SCHEMA_VERSION,
        "app_version": app_version or embedded_version or APP_VERSION,
        "settings": values,
    }
    if profile_name is not None:
        document["profile_name"] = profile_name
    return document


def migrate_settings_document(document: dict) -> dict[str, object]:
    """Return a v2 document, migrating the supported flat v1 format."""
    schema_version = document.get("schema_version")
    if schema_version == 1:
        metadata = {"schema_version", "app_version", "profile_name"}
        return build_settings_document(
            {key: value for key, value in document.items() if key not in metadata},
            profile_name=document.get("profile_name"),
            app_version=document.get("app_version"),
        )
    if schema_version != SETTINGS_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported analysis settings schema version {schema_version!r}; "
            f"supported versions are 1 and {SETTINGS_SCHEMA_VERSION}."
        )
    if document.get("schema_name") != SETTINGS_SCHEMA_NAME:
        raise ValueError(f"schema_name must be {SETTINGS_SCHEMA_NAME!r}.")
    if not isinstance(document.get("settings"), dict):
        raise ValueError("settings must be an object.")
    return dict(document)


def extract_analysis_settings_upload(contents: bytes, filename: str) -> bytes:
    """Extract the settings JSON bytes from JSON or an analysis-bundle ZIP."""
    suffix = Path(filename).suffix.lower()
    if suffix == ".json":
        if len(contents) > MAX_SETTINGS_BYTES:
            raise ValueError("Analysis settings JSON must be no larger than 1 MB.")
        return contents
    if suffix != ".zip":
        raise ValueError("Analysis settings must be uploaded as JSON or ZIP.")

    try:
        with ZipFile(BytesIO(contents)) as archive:
            try:
                entry = archive.getinfo(SETTINGS_BUNDLE_PATH)
            except KeyError as error:
                raise ValueError(
                    "Analysis bundle does not contain analysis_settings.json."
                ) from error
            if entry.flag_bits & 0x1:
                raise ValueError("Encrypted analysis settings are not supported.")
            if entry.file_size > MAX_SETTINGS_BYTES:
                raise ValueError("Analysis settings JSON must be no larger than 1 MB.")
            with archive.open(entry) as settings_file:
                settings_contents = settings_file.read(MAX_SETTINGS_BYTES + 1)
    except BadZipFile as error:
        raise ValueError("Analysis settings ZIP is invalid.") from error

    if len(settings_contents) > MAX_SETTINGS_BYTES:
        raise ValueError("Analysis settings JSON must be no larger than 1 MB.")
    return settings_contents


def load_analysis_settings_upload(contents: bytes, filename: str) -> dict[str, object]:
    """Load widget state from an exported JSON file or analysis-bundle ZIP."""
    return load_analysis_settings(extract_analysis_settings_upload(contents, filename))


def load_analysis_profile_upload(
    contents: bytes, filename: str
) -> tuple[dict[str, object], dict[str, str | None]]:
    """Load portable state and display metadata from a settings upload."""
    settings_contents = extract_analysis_settings_upload(contents, filename)
    restored = load_analysis_settings(settings_contents)
    settings = migrate_settings_document(json.loads(settings_contents.decode("utf-8")))
    metadata = {
        "profile_name": settings.get("profile_name"),
        "app_version": settings.get("app_version"),
    }
    for key, value in metadata.items():
        if value is not None and not isinstance(value, str):
            raise ValueError(f"{key} must be a string when provided.")
    return restored, metadata


def _number(settings: dict, name: str, minimum: float = 0.0) -> float | None:
    """Return a validated numeric setting when it is present."""
    if name not in settings:
        return None
    value = settings[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric.")
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite.")
    if value < minimum:
        raise ValueError(f"{name} must be at least {minimum}.")
    return float(value)


def _fraction(settings: dict, name: str) -> float | None:
    """Return a validated zero-to-one setting when it is present."""
    value = _number(settings, name)
    if value is not None and value > 1:
        raise ValueError(f"{name} must be between 0 and 1.")
    return value


def _integer(settings: dict, name: str, minimum: int = 0) -> int | None:
    """Return a validated integer setting when it is present."""
    if name not in settings:
        return None
    value = settings[name]
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer of at least {minimum}.")
    return value


def _string_list(settings: dict, name: str) -> list[str] | None:
    """Return a validated list of strings when it is present."""
    if name not in settings:
        return None
    value = settings[name]
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"{name} must be a list of strings.")
    return value


def load_analysis_settings(contents: bytes) -> dict[str, object]:
    """Return validated widget state from exported analysis-settings JSON."""
    try:
        settings = json.loads(contents.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("Analysis settings must be valid UTF-8 JSON.") from error
    if not isinstance(settings, dict):
        raise ValueError("Analysis settings JSON must contain an object.")
    document = migrate_settings_document(settings)
    settings = document["settings"]

    restored: dict[str, object] = {}
    top_n = _integer(settings, "top_n", 1)
    if top_n is not None:
        if top_n > 100:
            raise ValueError("top_n must be at most 100.")
        restored["top_n"] = top_n

    rank_by = settings.get("rank_by")
    if rank_by is not None:
        if rank_by not in RANK_BY_CHOICES:
            raise ValueError("rank_by is not supported by this version.")
        restored["rank_by"] = rank_by

    preset = settings.get("filter_preset", "Custom")
    if preset not in {*FILTER_PRESETS, "Custom"}:
        raise ValueError("filter_preset is not supported by this version.")
    restored["filter_preset"] = preset

    for name in ("min_pct", "max_pct_out", "min_pct_difference", "max_p_adj"):
        value = _fraction(settings, name)
        if value is not None:
            restored[name] = value
    for name in ("min_cluster_mean", "min_mean_difference", "min_log2fc"):
        value = _number(settings, name)
        if value is not None:
            restored[name] = value
    for name in ("min_cells_expressing",):
        value = _integer(settings, name)
        if value is not None:
            restored[name] = value
    min_cluster_size = _integer(settings, "min_cluster_size", 1)
    if min_cluster_size is not None:
        restored["min_cluster_size"] = min_cluster_size

    max_other_mean = settings.get("max_other_mean")
    restored["limit_other_mean"] = max_other_mean is not None
    if max_other_mean is not None:
        restored["max_other_mean"] = _number(settings, "max_other_mean")

    max_cells_out = settings.get("max_cells_expressing_out")
    restored["limit_cells_expressing_out"] = max_cells_out is not None
    if max_cells_out is not None:
        restored["max_cells_expressing_out"] = _integer(
            settings, "max_cells_expressing_out"
        )

    pseudocount = settings.get("pseudocount")
    if pseudocount is not None:
        restored["pseudocount"] = _number(settings, "pseudocount", 0.0001)

    for name, state_name in (
        ("excluded_genes", "excluded_gene_text"),
        ("excluded_prefixes", "excluded_prefix_text"),
    ):
        values = _string_list(settings, name)
        if values is not None:
            restored[state_name] = ", ".join(values)

    return restored
