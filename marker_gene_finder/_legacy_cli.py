"""Deprecated flat command interface retained through v2."""
from __future__ import annotations

from .progress import console_progress
from ._atomic import write_text

import argparse
import json
import sys
from pathlib import Path
from .version import APP_VERSION
from ._constants import (
    DEFAULT_GENE_CHUNK_SIZE,
    DEFAULT_MAX_DENSE_MEMORY_MB,
    RANK_BY_CHOICES,
)
from .analysis import run
from .batch import build_batch_plan, run_batch, verify_batch_manifest
from .benchmark import benchmark_analysis, compare_benchmarks
from .inspection import build_analysis_plan, inspect_input


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Rank marker genes for each cluster "
            "using log2 fold change."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {APP_VERSION}",
    )

    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "--input",
        help="Input expression CSV or H5AD path.",
    )
    inputs.add_argument(
        "--input-dir",
        help="Directory of CSV or H5AD files to analyze as a batch.",
    )
    inputs.add_argument(
        "--verify-batch-manifest",
        help="Verify batch input and output checksums and exit.",
    )
    inputs.add_argument(
        "--compare-benchmarks",
        nargs=2,
        metavar=("BASELINE", "CANDIDATE"),
        help="Compare two versioned benchmark JSON reports and exit.",
    )

    outputs = parser.add_mutually_exclusive_group()
    outputs.add_argument(
        "--output",
        help="Output marker CSV path.",
    )
    outputs.add_argument(
        "--output-dir",
        help="Directory for batch marker CSVs and batch_manifest.json.",
    )
    outputs.add_argument(
        "--verification-json",
        help="Output JSON path for --verify-batch-manifest.",
    )
    outputs.add_argument(
        "--comparison-json",
        help="Output JSON path for --compare-benchmarks.",
    )

    parser.add_argument(
        "--engine",
        choices=["native", "scanpy", "auto"],
        default="native",
        help=(
            "Marker-analysis engine. Auto follows the preflight resource "
            "recommendation; Scanpy requires H5AD input."
        ),
    )

    parser.add_argument(
        "--cluster-column",
        default="cluster",
        help="Cluster-label column name.",
    )

    parser.add_argument(
        "--cell-column",
        default="cell",
        help="Cell-ID column excluded from gene columns.",
    )
    parser.add_argument(
        "--max-dense-memory-mb",
        type=float,
        default=DEFAULT_MAX_DENSE_MEMORY_MB,
        help=(
            "Native H5AD dense-conversion safety limit in MiB; "
            "use 0 to disable it."
        ),
    )
    parser.add_argument(
        "--gene-chunk-size",
        type=int,
        default=DEFAULT_GENE_CHUNK_SIZE,
        help=(
            "Number of genes processed together by the Native engine "
            f"(default: {DEFAULT_GENE_CHUNK_SIZE})."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help=(
            "Native target clusters analyzed concurrently (default: 1). "
            "Higher values can increase peak memory use."
        ),
    )
    parser.add_argument(
        "--diagnostics-json",
        help=(
            "Write single-file runtime diagnostics as JSON. Batch diagnostics "
            "are included in batch_manifest.json."
        ),
    )
    parser.add_argument(
        "--inspect-json",
        help=(
            "Inspect input dimensions and estimated dense memory, write JSON, "
            "and exit without marker analysis."
        ),
    )
    parser.add_argument(
        "--plan-json",
        help=(
            "Write a versioned execution plan and exit without marker "
            "analysis. Requires a single --input."
        ),
    )
    parser.add_argument(
        "--benchmark-json",
        help=(
            "Run repeated single-input analyses, write benchmark JSON, and "
            "exit without retaining marker CSV output."
        ),
    )
    parser.add_argument(
        "--benchmark-repeats",
        type=int,
        default=3,
        help="Number of benchmark repetitions (default: 3).",
    )
    parser.add_argument(
        "--max-regression-percent",
        type=float,
        default=10.0,
        help="Maximum allowed benchmark median regression (default: 10).",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume a batch by skipping unchanged completed inputs.",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Retry failed batch inputs while skipping unchanged successes.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop a batch after its first failed input.",
    )
    parser.add_argument(
        "--include-pattern",
        action="append",
        default=[],
        help="Batch filename glob to include; repeat to combine patterns.",
    )
    parser.add_argument(
        "--exclude-pattern",
        action="append",
        default=[],
        help="Batch filename glob to exclude; repeat to combine patterns.",
    )
    parser.add_argument(
        "--batch-plan-json",
        help="Write a batch run/skip plan without running analyses.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-cluster progress messages.",
    )

    parser.add_argument(
        "--target-cluster",
        action="append",
        dest="target_clusters",
        help="Cluster to analyze; repeat to select multiple clusters.",
    )
    parser.add_argument(
        "--min-cluster-size",
        type=int,
        default=1,
        help="Minimum cell count required for each target cluster.",
    )
    parser.add_argument(
        "--exclude-gene", action="append", default=[],
        help="Exact gene name to exclude; repeat as needed.",
    )
    parser.add_argument(
        "--exclude-prefix", action="append", default=[],
        help="Gene-name prefix to exclude; repeat as needed.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="Number of markers returned per cluster.",
    )

    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.1,
        help="Positive pseudocount used for log2FC.",
    )
    parser.add_argument(
        "--min-pct",
        type=float,
        default=0.0,
        help=(
            "Minimum fraction of cluster cells "
            "expressing a gene."
        ),
    )
    parser.add_argument(
        "--min-cells-expressing",
        type=int,
        default=0,
        help="Minimum target-cluster cells expressing a marker gene.",
    )
    parser.add_argument(
        "--min-cluster-mean",
        type=float,
        default=0.0,
        help="Minimum mean expression required inside the target cluster.",
    )
    parser.add_argument(
        "--min-mean-difference",
        type=float,
        default=0.0,
        help="Minimum target-minus-reference mean-expression difference.",
    )
    parser.add_argument(
        "--max-other-mean",
        type=float,
        default=None,
        help=(
            "Maximum mean expression allowed outside the target cluster; "
            "omit to disable this filter."
        ),
    )
    parser.add_argument(
        "--max-cells-expressing-out",
        type=int,
        default=None,
        help=(
            "Maximum number of reference cells expressing a marker; "
            "omit to disable this filter."
        ),
    )

    parser.add_argument(
        "--min-log2fc",
        type=float,
        default=0.0,
        help=(
            "Minimum log2 fold change required "
            "for a marker gene."
        ),
    )

    parser.add_argument(
        "--max-pct-out",
        type=float,
        default=1.0,
        help="Maximum fraction of non-cluster cells expressing a marker.",
    )

    parser.add_argument(
        "--rank-by",
        choices=RANK_BY_CHOICES,
        default="log2FC",
        help="Column used to rank filtered markers within each cluster.",
    )

    parser.add_argument(
        "--min-pct-difference",
        type=float,
        default=0.0,
        help="Minimum difference between in-cluster and outside expression fractions.",
    )

    parser.add_argument(
        "--max-p-adj",
        type=float,
        default=1.0,
        help=(
            "Maximum Benjamini-Hochberg adjusted "
            "p-value allowed for a marker gene."
        ),
    )

    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--layer",
        default=None,
        help="H5AD layer name; omit to use AnnData.X.",
    )

    source.add_argument(
        "--use-raw",
        action="store_true",
        help="Use AnnData.raw expression and gene names (H5AD only).",
    )

    return parser


def main() -> None:
    """Run the command-line application."""

    args = build_parser().parse_args()

    batch_only_options = (
        args.resume,
        args.retry_failed,
        args.fail_fast,
        bool(args.include_pattern),
        bool(args.exclude_pattern),
        bool(args.batch_plan_json),
    )
    if any(batch_only_options) and not args.input_dir:
        raise SystemExit(
            "Batch recovery, selection, and planning options require --input-dir."
        )

    if args.verify_batch_manifest:
        if not args.verification_json:
            raise SystemExit(
                "--verify-batch-manifest requires --verification-json."
            )
        verification = verify_batch_manifest(args.verify_batch_manifest)
        verification_path = Path(args.verification_json)
        verification_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(verification_path,
            json.dumps(verification, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved batch verification to {verification_path}")
        if not verification["valid"]:
            raise SystemExit(1)
        return

    if args.compare_benchmarks:
        if not args.comparison_json:
            raise SystemExit("--compare-benchmarks requires --comparison-json.")
        comparison = compare_benchmarks(
            *args.compare_benchmarks,
            max_regression_percent=args.max_regression_percent,
        )
        comparison_path = Path(args.comparison_json)
        comparison_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(comparison_path,
            json.dumps(comparison, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved benchmark comparison to {comparison_path}")
        if not comparison["passed"]:
            raise SystemExit(1)
        return

    exclusive_modes = [args.inspect_json, args.plan_json, args.benchmark_json]
    if sum(value is not None for value in exclusive_modes) > 1:
        raise SystemExit(
            "--inspect-json, --plan-json, and --benchmark-json are mutually "
            "exclusive."
        )

    if args.plan_json:
        if not args.input or args.output or args.input_dir or args.output_dir:
            raise SystemExit(
                "--plan-json requires --input and cannot be combined with "
                "analysis outputs or batch input."
            )
        if args.diagnostics_json:
            raise SystemExit(
                "--plan-json cannot be combined with --diagnostics-json."
            )
        plan = build_analysis_plan(
            args.input,
            requested_engine=args.engine,
            cluster_column=args.cluster_column,
            cell_column=args.cell_column,
            layer=args.layer,
            use_raw=args.use_raw,
            max_dense_memory_mb=args.max_dense_memory_mb,
            gene_chunk_size=args.gene_chunk_size,
            workers=args.workers,
        )
        plan_path = Path(args.plan_json)
        plan_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(plan_path,
            json.dumps(plan, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved analysis plan to {plan_path}")
        return

    if args.inspect_json:
        if not args.input or args.output or args.input_dir or args.output_dir:
            raise SystemExit(
                "--inspect-json requires --input and cannot be combined with "
                "analysis outputs or batch input."
            )
        if args.diagnostics_json:
            raise SystemExit(
                "--inspect-json cannot be combined with --diagnostics-json."
            )
        inspection = inspect_input(
            args.input,
            cluster_column=args.cluster_column,
            cell_column=args.cell_column,
            layer=args.layer,
            use_raw=args.use_raw,
            max_dense_memory_mb=args.max_dense_memory_mb,
        )
        inspection_path = Path(args.inspect_json)
        inspection_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(inspection_path,
            json.dumps(inspection, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved input inspection to {inspection_path}")
        return

    analysis_options = dict(
        cluster_column=args.cluster_column,
        cell_column=args.cell_column,
        top_n=args.top_n,
        pseudocount=args.pseudocount,
        min_pct=args.min_pct,
        min_cells_expressing=args.min_cells_expressing,
        min_cluster_mean=args.min_cluster_mean,
        min_mean_difference=args.min_mean_difference,
        max_other_mean=args.max_other_mean,
        max_cells_expressing_out=args.max_cells_expressing_out,
        max_pct_out=args.max_pct_out,
        rank_by=args.rank_by,
        min_pct_difference=args.min_pct_difference,
        min_log2fc=args.min_log2fc,
        max_p_adj=args.max_p_adj,
        engine=args.engine,
        layer=args.layer,
        use_raw=args.use_raw,
        target_clusters=args.target_clusters,
        min_cluster_size=args.min_cluster_size,
        excluded_genes=args.exclude_gene,
        excluded_prefixes=args.exclude_prefix,
        max_dense_memory_mb=args.max_dense_memory_mb,
        gene_chunk_size=args.gene_chunk_size,
        progress_callback=None if args.quiet else console_progress,
        workers=args.workers,
    )

    if args.batch_plan_json:
        if not args.input_dir or not args.output_dir or args.input or args.output:
            raise SystemExit(
                "--batch-plan-json requires --input-dir and --output-dir."
            )
        if args.diagnostics_json:
            raise SystemExit(
                "--batch-plan-json cannot be combined with --diagnostics-json."
            )
        if args.inspect_json or args.plan_json or args.benchmark_json:
            raise SystemExit(
                "--batch-plan-json cannot be combined with another planning, "
                "inspection, or benchmark mode."
            )
        batch_plan = build_batch_plan(
            args.input_dir,
            args.output_dir,
            resume=args.resume,
            retry_failed=args.retry_failed,
            include_patterns=args.include_pattern,
            exclude_patterns=args.exclude_pattern,
            **analysis_options,
        )
        batch_plan_path = Path(args.batch_plan_json)
        batch_plan_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(batch_plan_path,
            json.dumps(batch_plan, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved batch plan to {batch_plan_path}")
        return

    if args.benchmark_json:
        if not args.input or args.output or args.input_dir or args.output_dir:
            raise SystemExit(
                "--benchmark-json requires --input and cannot be combined "
                "with analysis outputs or batch input."
            )
        if args.diagnostics_json:
            raise SystemExit(
                "--benchmark-json cannot be combined with --diagnostics-json."
            )
        benchmark = benchmark_analysis(
            args.input,
            repeats=args.benchmark_repeats,
            **analysis_options,
        )
        benchmark_path = Path(args.benchmark_json)
        benchmark_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(benchmark_path,
            json.dumps(benchmark, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(f"Saved benchmark to {benchmark_path}")
        return

    if bool(args.input_dir) != bool(args.output_dir):
        raise SystemExit("--input-dir and --output-dir must be used together.")
    if bool(args.input) != bool(args.output):
        raise SystemExit("--input and --output must be used together.")
    if not args.input and not args.input_dir:
        raise SystemExit(
            "Provide --input with --output, or --input-dir with --output-dir."
        )
    if args.input_dir and args.diagnostics_json:
        raise SystemExit(
            "--diagnostics-json is for single-file analysis; batch diagnostics "
            "are written to batch_manifest.json."
        )

    if args.input_dir:
        manifest = run_batch(
            args.input_dir,
            args.output_dir,
            resume=args.resume,
            retry_failed=args.retry_failed,
            fail_fast=args.fail_fast,
            include_patterns=args.include_pattern,
            exclude_patterns=args.exclude_pattern,
            **analysis_options,
        )
        completed = sum(
            row["status"] in {"complete", "skipped"} for row in manifest
        )
        failed = sum(row["status"] == "failed" for row in manifest)
        print(
            f"Batch complete: {completed} succeeded, {failed} failed. "
            f"Manifest: {Path(args.output_dir) / 'batch_manifest.json'}"
        )
        if failed:
            raise SystemExit(1)
        return

    diagnostics: dict[str, object] = {}
    markers = run(
        input_path=args.input,
        output_path=args.output,
        diagnostics=diagnostics,
        **analysis_options,
    )

    if args.diagnostics_json:
        diagnostics_path = Path(args.diagnostics_json)
        diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
        write_text(diagnostics_path,
            json.dumps(diagnostics, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    print(
        f"Saved {len(markers)} marker rows "
        f"to {args.output}"
    )
