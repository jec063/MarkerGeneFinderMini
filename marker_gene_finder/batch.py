"""Batch selection, execution, recovery, and verification."""
from __future__ import annotations

from ._atomic import atomic_path, write_text

import fnmatch
import hashlib
import json
from collections.abc import Sequence
from pathlib import Path
from time import perf_counter
from .version import APP_VERSION
from ._constants import (
    BATCH_MANIFEST_SCHEMA_VERSION,
    BATCH_PLAN_SCHEMA_VERSION,
    BATCH_VERIFICATION_SCHEMA_VERSION,
)
from .analysis import _registered_run
from .inspection import build_analysis_plan


def _sha256_file(path: str | Path) -> str:
    """Return a streaming SHA-256 digest for a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _select_batch_inputs(
    input_dir: Path,
    include_patterns: Sequence[str] = (),
    exclude_patterns: Sequence[str] = (),
) -> list[Path]:
    """Select supported batch inputs using ordered filename filtering."""
    inputs = sorted(
        path for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".csv", ".h5ad"}
    )
    if include_patterns:
        inputs = [
            path for path in inputs
            if any(fnmatch.fnmatchcase(path.name, pattern) for pattern in include_patterns)
        ]
    if exclude_patterns:
        inputs = [
            path for path in inputs
            if not any(
                fnmatch.fnmatchcase(path.name, pattern)
                for pattern in exclude_patterns
            )
        ]
    return inputs


def _batch_settings_fingerprint(analysis_options: dict[str, object]) -> str:
    """Hash analysis settings while excluding non-serializable callbacks."""
    fingerprint_options = {
        key: value for key, value in analysis_options.items()
        if key != "progress_callback"
    }
    settings_json = json.dumps(
        fingerprint_options,
        sort_keys=True,
        default=str,
        allow_nan=False,
    )
    return hashlib.sha256(settings_json.encode()).hexdigest()


def build_batch_plan(
    input_dir: str | Path,
    output_dir: str | Path,
    resume: bool = False,
    retry_failed: bool = False,
    include_patterns: Sequence[str] = (),
    exclude_patterns: Sequence[str] = (),
    registry=None,
    **analysis_options,
) -> dict[str, object]:
    """Describe which batch inputs would run or be skipped without analysis."""
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if not input_dir.is_dir():
        raise ValueError(f"Batch input directory does not exist: {input_dir}")
    inputs = _select_batch_inputs(input_dir, include_patterns, exclude_patterns)
    if not inputs:
        raise ValueError("Batch selection contains no CSV or H5AD files.")

    manifest_path = output_dir / "batch_manifest.json"
    previous_runs: dict[str, dict[str, object]] = {}
    if (resume or retry_failed) and manifest_path.is_file():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        previous_runs = {
            str(row["input"]): row for row in previous.get("runs", [])
        }
    settings_fingerprint = _batch_settings_fingerprint(analysis_options)
    entries: list[dict[str, object]] = []
    for input_path in inputs:
        output_path = output_dir / f"{input_path.name}.markers.csv"
        stat = input_path.stat()
        input_fingerprint = {
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": _sha256_file(input_path),
        }
        previous = previous_runs.get(str(input_path))
        action = "run"
        reason = "resume_disabled"
        if resume or retry_failed:
            if previous is None:
                reason = "new_input"
            elif previous.get("status") not in {"complete", "skipped"}:
                reason = f"previous_{previous.get('status', 'unknown')}"
            elif previous.get("input_fingerprint") != input_fingerprint:
                reason = "input_changed"
            elif previous.get("settings_fingerprint") != settings_fingerprint:
                reason = "settings_changed"
            elif not output_path.is_file():
                reason = "output_missing"
            elif not previous.get("output_sha256"):
                reason = "output_hash_missing"
            elif previous["output_sha256"] != _sha256_file(output_path):
                reason = "output_changed"
            else:
                action = "skip"
                reason = "unchanged_success"
        execution = build_analysis_plan(input_path, requested_engine=analysis_options.get("engine", "native"), registry=registry,
            **{key: value for key, value in analysis_options.items() if key in {"cluster_column", "cell_column", "layer", "use_raw", "workers", "gene_chunk_size", "max_dense_memory_mb"}})
        entries.append({
            "can_execute": execution["can_execute"],
            "warnings": execution["warnings"],
            "input": str(input_path),
            "output": str(output_path),
            "action": action,
            "reason": reason,
            "input_fingerprint": input_fingerprint,
        })
    return {
        "schema_version": BATCH_PLAN_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "settings_fingerprint": settings_fingerprint,
        "include_patterns": list(include_patterns),
        "exclude_patterns": list(exclude_patterns),
        "selected": len(entries),
        "run": sum(entry["action"] == "run" for entry in entries),
        "skip": sum(entry["action"] == "skip" for entry in entries),
        "entries": entries,
    }


def run_batch(
    input_dir: str | Path,
    output_dir: str | Path,
    resume: bool = False,
    retry_failed: bool = False,
    fail_fast: bool = False,
    include_patterns: Sequence[str] = (),
    exclude_patterns: Sequence[str] = (),
    registry=None,
    **analysis_options,
) -> list[dict[str, object]]:
    """Analyze supported files with resumable, atomic manifest tracking."""
    batch_started = perf_counter()
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    if not input_dir.is_dir():
        raise ValueError(f"Batch input directory does not exist: {input_dir}")

    inputs = _select_batch_inputs(input_dir, include_patterns, exclude_patterns)
    if not inputs:
        raise ValueError("Batch selection contains no CSV or H5AD files.")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "batch_manifest.json"
    previous_runs: dict[str, dict[str, object]] = {}
    if (resume or retry_failed) and manifest_path.is_file():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        previous_runs = {
            str(row["input"]): row for row in previous.get("runs", [])
        }

    settings_fingerprint = _batch_settings_fingerprint(analysis_options)
    manifest: list[dict[str, object]] = []
    stopped = False
    for input_path in inputs:
        output_path = output_dir / f"{input_path.name}.markers.csv"
        stat = input_path.stat()
        input_fingerprint = {
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": _sha256_file(input_path),
        }
        record: dict[str, object] = {
            "input": str(input_path),
            "output": str(output_path),
            "input_fingerprint": input_fingerprint,
            "settings_fingerprint": settings_fingerprint,
        }
        if stopped:
            record.update(
                status="pending", marker_rows=0,
                error="Not attempted because fail_fast stopped the batch.",
                duration_seconds=0.0,
            )
            manifest.append(record)
            continue
        previous = previous_runs.get(str(input_path))
        output_matches = (
            output_path.is_file()
            and previous is not None
            and previous.get("output_sha256") is not None
            and previous.get("output_sha256") == _sha256_file(output_path)
        )
        if (
            previous
            and previous.get("status") in {"complete", "skipped"}
            and previous.get("input_fingerprint") == input_fingerprint
            and previous.get("settings_fingerprint") == settings_fingerprint
            and output_matches
        ):
            record.update(
                status="skipped",
                marker_rows=previous.get("marker_rows", 0),
                error=None,
                diagnostics=previous.get("diagnostics", {}),
                output_sha256=previous["output_sha256"],
                duration_seconds=0.0,
            )
            manifest.append(record)
            continue
        run_diagnostics: dict[str, object] = {}
        attempt_started = perf_counter()
        try:
            with atomic_path(output_path) as temporary_output:
                markers = _registered_run(
                    input_path,
                    temporary_output,
                    registry=registry,
                    diagnostics=run_diagnostics,
                    **analysis_options,
                )
        except Exception as error:
            record.update(status="failed", marker_rows=0, error=str(error))
            stopped = fail_fast
        else:
            run_diagnostics["output"] = str(output_path)
            record.update(
                status="complete",
                marker_rows=len(markers),
                error=None,
                diagnostics=run_diagnostics,
                output_sha256=_sha256_file(output_path),
            )
        record["duration_seconds"] = perf_counter() - attempt_started
        manifest.append(record)

    document = {
        "schema_version": BATCH_MANIFEST_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "settings_fingerprint": settings_fingerprint,
        "include_patterns": list(include_patterns),
        "exclude_patterns": list(exclude_patterns),
        "inputs": len(inputs),
        "completed": sum(
            row["status"] in {"complete", "skipped"} for row in manifest
        ),
        "skipped": sum(row["status"] == "skipped" for row in manifest),
        "failed": sum(row["status"] == "failed" for row in manifest),
        "pending": sum(row["status"] == "pending" for row in manifest),
        "duration_seconds": perf_counter() - batch_started,
        "runs": manifest,
    }
    write_text(manifest_path,
        json.dumps(document, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify_batch_manifest(manifest_path: str | Path) -> dict[str, object]:
    """Verify completed batch inputs and outputs against recorded hashes."""
    manifest_path = Path(manifest_path)
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    if document.get("schema_version") != BATCH_MANIFEST_SCHEMA_VERSION:
        raise ValueError("Unsupported batch manifest schema version.")
    if not isinstance(document.get("runs"), list):
        raise ValueError("Batch manifest runs must be a list.")

    entries: list[dict[str, object]] = []
    for record in document["runs"]:
        status = record.get("status")
        entry: dict[str, object] = {
            "input": record.get("input"),
            "output": record.get("output"),
            "batch_status": status,
            "valid": True,
            "issues": [],
        }
        issues = entry["issues"]
        if status in {"complete", "skipped"}:
            input_path = Path(str(record.get("input", "")))
            output_path = Path(str(record.get("output", "")))
            fingerprint = record.get("input_fingerprint", {})
            if not input_path.is_file():
                issues.append("input_missing")
            elif not isinstance(fingerprint, dict):
                issues.append("input_fingerprint_missing")
            else:
                stat = input_path.stat()
                if fingerprint.get("size") != stat.st_size:
                    issues.append("input_size_mismatch")
                if fingerprint.get("mtime_ns") != stat.st_mtime_ns:
                    issues.append("input_mtime_mismatch")
                expected_input_hash = fingerprint.get("sha256")
                if not expected_input_hash:
                    issues.append("input_sha256_missing")
                elif expected_input_hash != _sha256_file(input_path):
                    issues.append("input_sha256_mismatch")
            if not output_path.is_file():
                issues.append("output_missing")
            else:
                expected_output_hash = record.get("output_sha256")
                if not expected_output_hash:
                    issues.append("output_sha256_missing")
                elif expected_output_hash != _sha256_file(output_path):
                    issues.append("output_sha256_mismatch")
        entry["valid"] = not issues
        entries.append(entry)

    invalid = sum(not bool(entry["valid"]) for entry in entries)
    return {
        "schema_version": BATCH_VERIFICATION_SCHEMA_VERSION,
        "app_version": APP_VERSION,
        "manifest": str(manifest_path),
        "valid": invalid == 0,
        "checked": len(entries),
        "invalid": invalid,
        "entries": entries,
    }

__all__ = [
    "BATCH_MANIFEST_SCHEMA_VERSION",
    "BATCH_PLAN_SCHEMA_VERSION",
    "BATCH_VERIFICATION_SCHEMA_VERSION",
    "build_batch_plan",
    "run_batch",
    "verify_batch_manifest",
]
