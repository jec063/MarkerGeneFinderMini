"""Dedicated compatibility coverage; removable with the v3 legacy shims."""
from pathlib import Path
import pytest
from src.marker_finder import build_parser
from src.version import APP_VERSION

PROJECT_ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_DATA = PROJECT_ROOT / "data/example_expression.csv"
FIXTURES = PROJECT_ROOT / "tests/fixtures/v2_7"


def test_parser_reports_application_version(capsys) -> None:
    with pytest.raises(SystemExit, match="0"):
        build_parser().parse_args(["--version"])

    assert capsys.readouterr().out.strip().endswith(APP_VERSION)

def test_parser_accepts_batch_directories() -> None:
    args = build_parser().parse_args(
        ["--input-dir", "inputs", "--output-dir", "outputs"]
    )

    assert args.input_dir == "inputs"
    assert args.output_dir == "outputs"
    assert args.input is None
    assert args.output is None

def test_parser_accepts_diagnostics_path() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--diagnostics-json", "diagnostics.json",
        ]
    )

    assert args.diagnostics_json == "diagnostics.json"

def test_parser_accepts_native_gene_chunk_size() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--gene-chunk-size", "17",
        ]
    )

    assert args.gene_chunk_size == 17

def test_parser_accepts_quiet_progress_mode() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--quiet",
        ]
    )

    assert args.quiet is True

def test_parser_accepts_native_workers() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--workers", "3",
        ]
    )

    assert args.workers == 3

def test_parser_accepts_input_inspection_without_analysis_output() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--inspect-json", "inspection.json",
        ]
    )

    assert args.output is None
    assert args.inspect_json == "inspection.json"

def test_parser_accepts_analysis_plan_and_auto_engine() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.h5ad",
            "--engine", "auto",
            "--plan-json", "plan.json",
        ]
    )

    assert args.output is None
    assert args.engine == "auto"
    assert args.plan_json == "plan.json"

def test_parser_accepts_batch_recovery_and_benchmark_options() -> None:
    batch = build_parser().parse_args(
        [
            "--input-dir", "inputs", "--output-dir", "outputs",
            "--resume", "--retry-failed", "--fail-fast",
        ]
    )
    benchmark = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--benchmark-json", "benchmark.json",
            "--benchmark-repeats", "2",
        ]
    )

    assert batch.resume and batch.retry_failed and batch.fail_fast
    assert benchmark.benchmark_json == "benchmark.json"
    assert benchmark.benchmark_repeats == 2

def test_parser_accepts_batch_selection_and_dry_run_plan() -> None:
    args = build_parser().parse_args(
        [
            "--input-dir", "inputs",
            "--output-dir", "outputs",
            "--include-pattern", "sample-*.csv",
            "--include-pattern", "*.h5ad",
            "--exclude-pattern", "*-draft.csv",
            "--resume",
            "--batch-plan-json", "batch-plan.json",
        ]
    )

    assert args.include_pattern == ["sample-*.csv", "*.h5ad"]
    assert args.exclude_pattern == ["*-draft.csv"]
    assert args.batch_plan_json == "batch-plan.json"

def test_main_rejects_batch_selection_for_single_input(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "marker-gene-finder",
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--include-pattern", "*.csv",
        ],
    )

    from src.marker_finder import main

    with pytest.raises(SystemExit, match="require --input-dir"):
        main()

def test_parser_accepts_verification_and_benchmark_comparison() -> None:
    verification = build_parser().parse_args(
        [
            "--verify-batch-manifest", "batch_manifest.json",
            "--verification-json", "verification.json",
        ]
    )
    comparison = build_parser().parse_args(
        [
            "--compare-benchmarks", "baseline.json", "candidate.json",
            "--comparison-json", "comparison.json",
            "--max-regression-percent", "15",
        ]
    )

    assert verification.verify_batch_manifest == "batch_manifest.json"
    assert verification.verification_json == "verification.json"
    assert comparison.compare_benchmarks == ["baseline.json", "candidate.json"]
    assert comparison.comparison_json == "comparison.json"
    assert comparison.max_regression_percent == 15

def test_parser_accepts_dense_memory_limit() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.h5ad",
            "--output", "markers.csv",
            "--max-dense-memory-mb", "512",
        ]
    )

    assert args.max_dense_memory_mb == 512

def test_parser_accepts_min_cells_expressing() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cells-expressing", "3",
        ]
    )

    assert args.min_cells_expressing == 3

def test_parser_accepts_min_cluster_mean() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cluster-mean", "1.5",
        ]
    )

    assert args.min_cluster_mean == pytest.approx(1.5)

def test_parser_accepts_min_mean_difference() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-mean-difference", "0.75",
        ]
    )

    assert args.min_mean_difference == pytest.approx(0.75)

def test_parser_accepts_max_other_mean() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--max-other-mean", "0.5",
        ]
    )

    assert args.max_other_mean == pytest.approx(0.5)

def test_parser_accepts_max_reference_expressing_cells() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--max-cells-expressing-out", "3",
        ]
    )

    assert args.max_cells_expressing_out == 3

def test_parser_accepts_min_log2fc() -> None:
    args = build_parser().parse_args(
        [
            "--input",
            "input.csv",
            "--output",
            "output.csv",
            "--min-log2fc",
            "0.5",
        ]
    )

    assert args.min_log2fc == pytest.approx(0.5)

def test_parser_accepts_max_p_adj() -> None:
    args = build_parser().parse_args(
        [
            "--input",
            "input.csv",
            "--output",
            "output.csv",
            "--max-p-adj",
            "0.05",
        ]
    )

    assert args.max_p_adj == pytest.approx(0.05)

def test_parser_accepts_min_cluster_size() -> None:
    args = build_parser().parse_args(
        [
            "--input", "expression.csv",
            "--output", "markers.csv",
            "--min-cluster-size", "10",
        ]
    )

    assert args.min_cluster_size == 10
