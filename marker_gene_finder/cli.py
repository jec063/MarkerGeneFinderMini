"""Version 2 command-line interface with explicit subcommands."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .api import analyze_file
from .batch import build_batch_plan, run_batch, verify_batch_manifest
from .benchmark import benchmark_analysis, compare_benchmarks
from .inspection import build_analysis_plan, inspect_input
from .engines import inspect_engines, load_engine_plugins
from .config import AnalysisConfig
from .exceptions import MarkerFinderError
from .version import APP_VERSION


def _add_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        type=Path,
        help="Versioned v2 analysis configuration JSON (defaults if omitted).",
    )


def _add_batch_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    _add_config_argument(parser)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--retry-failed", action="store_true")
    parser.add_argument("--include-pattern", action="append", default=[])
    parser.add_argument("--exclude-pattern", action="append", default=[])


def build_parser() -> argparse.ArgumentParser:
    """Create the stable v2 subcommand parser."""
    parser = argparse.ArgumentParser(
        prog="marker-gene-finder",
        description="Find marker genes with a typed, reproducible configuration.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {APP_VERSION}"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    analyze = commands.add_parser("analyze", help="Analyze one CSV or H5AD file.")
    analyze.add_argument("input", type=Path)
    analyze.add_argument("output", type=Path)
    _add_config_argument(analyze)
    analyze.add_argument("--bundle", type=Path, help="Write a reproducibility ZIP.")
    analyze.add_argument("--metadata", type=Path, help="Write result metadata as JSON.")
    analyze.add_argument("--hash-input", action="store_true", help="Include the input SHA-256 digest.")
    analyze.add_argument(
        "--load-engine-plugins",
        action="store_true",
        help="Load installed marker_gene_finder.engines entry points.",
    )

    list_engines = commands.add_parser(
        "list-engines", help="List built-in engines and inspect installed plugins."
    )
    list_engines.add_argument("--load-plugins", action="store_true")
    list_engines.add_argument("--output", type=Path)

    inspect = commands.add_parser("inspect", help="Inspect an input without analysis.")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--output", type=Path, required=True)
    _add_config_argument(inspect)

    plan = commands.add_parser("plan", help="Create a non-destructive execution plan.")
    plan.add_argument("input", type=Path)
    plan.add_argument("--output", type=Path, required=True)
    _add_config_argument(plan)

    batch = commands.add_parser("batch", help="Analyze a directory of inputs.")
    _add_batch_arguments(batch)
    batch.add_argument("--fail-fast", action="store_true")

    batch_plan = commands.add_parser(
        "batch-plan", help="Preview batch actions without analysis."
    )
    _add_batch_arguments(batch_plan)
    batch_plan.add_argument("--output", type=Path, required=True)

    verify = commands.add_parser("verify", help="Verify a batch manifest.")
    verify.add_argument("manifest", type=Path)
    verify.add_argument("--output", type=Path, required=True)

    benchmark = commands.add_parser("benchmark", help="Benchmark one analysis.")
    benchmark.add_argument("input", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)
    benchmark.add_argument("--repeats", type=int, default=3)
    _add_config_argument(benchmark)

    compare = commands.add_parser(
        "compare-benchmarks", help="Compare two benchmark reports."
    )
    compare.add_argument("baseline", type=Path)
    compare.add_argument("candidate", type=Path)
    compare.add_argument("--output", type=Path, required=True)
    compare.add_argument("--max-regression-percent", type=float, default=10.0)
    for command_parser in (plan, batch, batch_plan, benchmark):
        command_parser.add_argument("--load-engine-plugins", action="store_true")
    return parser


def _config(path: Path | None) -> AnalysisConfig:
    return AnalysisConfig.from_json(path) if path else AnalysisConfig()


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def _inspection_options(config: AnalysisConfig) -> dict[str, Any]:
    return {
        "cluster_column": config.cluster_column,
        "cell_column": config.cell_column,
        "layer": config.layer,
        "use_raw": config.use_raw,
        "max_dense_memory_mb": config.max_dense_memory_mb,
    }


def execute(args: argparse.Namespace) -> int:
    """Execute one parsed v2 command and return its process status."""
    command = args.command
    registry = load_engine_plugins() if getattr(args, "load_engine_plugins", False) else None
    if command == "analyze":
        options = {"registry": registry} if registry is not None else {}
        paths = [args.input, args.output, args.bundle, args.metadata, args.config]
        resolved = [path.resolve() for path in paths if path is not None]
        if len(resolved) != len(set(resolved)):
            raise ValueError("Input, config, output, bundle and metadata paths must be distinct.")
        result = analyze_file(
            args.input, args.output, _config(args.config), hash_input_file=args.hash_input, **options
        )
        if args.bundle is not None:
            result.to_bundle(args.bundle)
        if args.metadata is not None:
            _write_json(args.metadata, result.metadata)
        print(f"Wrote {len(result.markers)} marker rows to {args.output}.")
        return 0
    if command == "list-engines":
        records = [
            record.to_dict()
            for record in inspect_engines(load_plugins=args.load_plugins)
        ]
        document = {"app_version": APP_VERSION, "engines": records}
        if args.output is not None:
            _write_json(args.output, document)
        else:
            print(json.dumps(document, indent=2, allow_nan=False))
        return 1 if any(row["status"] == "error" for row in records) else 0
    if command == "inspect":
        document = inspect_input(args.input, **_inspection_options(_config(args.config)))
        _write_json(args.output, document)
        return 0
    if command == "plan":
        config = _config(args.config)
        document = build_analysis_plan(
            args.input,
            requested_engine=config.engine,
            registry=registry,
            gene_chunk_size=config.gene_chunk_size,
            workers=config.workers,
            **_inspection_options(config),
        )
        _write_json(args.output, document)
        return 0 if document["can_execute"] else 1
    if command in {"batch", "batch-plan"}:
        config = _config(args.config)
        options = {
            "registry": registry,
            "resume": args.resume,
            "retry_failed": args.retry_failed,
            "include_patterns": args.include_pattern,
            "exclude_patterns": args.exclude_pattern,
            **config.to_analysis_kwargs(),
        }
        if command == "batch-plan":
            document = build_batch_plan(args.input_dir, args.output_dir, **options)
            _write_json(args.output, document)
            return 0 if all(row["can_execute"] for row in document["entries"]) else 1
        records = run_batch(
            args.input_dir, args.output_dir, fail_fast=args.fail_fast, **options
        )
        return 1 if any(row["status"] == "failed" for row in records) else 0
    if command == "verify":
        document = verify_batch_manifest(args.manifest)
        _write_json(args.output, document)
        return 0 if document["valid"] else 1
    if command == "benchmark":
        document = benchmark_analysis(
            args.input, repeats=args.repeats, registry=registry, **_config(args.config).to_analysis_kwargs()
        )
        _write_json(args.output, document)
        return 0
    if command == "compare-benchmarks":
        document = compare_benchmarks(
            args.baseline, args.candidate, args.max_regression_percent
        )
        _write_json(args.output, document)
        return 0 if document["passed"] else 1
    raise AssertionError(f"Unhandled command: {command}")


def main(argv: list[str] | None = None) -> int:
    """Run the v2 CLI and present expected failures without a traceback."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return execute(args)
    except (MarkerFinderError, ValueError, OSError, json.JSONDecodeError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
