# v2.8 execution-core validation

The baseline is v2.7.0, commit `e6e1ab76736045336013832a3a31b8f76a1e729e`.
Reference CSVs in `tests/fixtures/v2_7` were generated from that unmodified
implementation before extraction, using `dataset(30, 18)` in
`execution_core.py`, seed 280, and `top_n=8`. Native and Scanpy have independent
expected outputs. Tests preserve exact row/column order, labels, counts, and
ranks; floating-point results use `rtol=1e-6`, `atol=1e-10` to accommodate
supported dependency/platform variation. X, layers, raw, CSV, chunking, and
Native workers are exercised. The existing scientific tests remain in place.

## Runtime and memory

Run from a checkout with its dependencies installed:

```sh
python benchmarks/execution_core.py --output small.json
python benchmarks/execution_core.py --cells 3000 --genes 1500 --output large.json
```

Run the same script with `PYTHONPATH` pointing to a v2.7.0 source checkout for
the baseline. Both versions must use the same Python and dependency environment.
Each case uses a fresh subprocess, one untimed warmup, five timed runs, and
one thread for BLAS/OpenMP/Numba. Timings include loading, analysis, and output;
data generation is excluded. Linux peak RSS includes dependencies, generation,
warmup, and all runs; it is a process high-water mark, not an allocation delta.

Inputs are seeded Poisson matrices with three clusters and enriched marker
genes. Cases cover Native CSV, Native sparse H5AD converted to dense, Scanpy
sparse H5AD, and a two-file Native batch. Small inputs have 600 cells × 300
genes; large inputs have 3,000 cells × 1,500 genes. Raw timing samples and
dependency versions are retained in the adjacent JSON files. These synthetic
workloads check the refactor; they do not establish performance for every atlas.

The release check flags median runtime or peak RSS increases above 10% for
investigation. Small differences include filesystem/cache and shared-host noise;
atomic staging also adds file operations. No statistical algorithm changed.

## Output integrity

Tests inject serialization, replacement, and interruption failures and confirm
that prior files survive, new destinations stay absent, and temporary files are
cleaned up on Python exceptions, including `KeyboardInterrupt`. Tests also
cover rejected plugin results, interrupted batch execution, and bundle exports.
Atomic replacement is per file; a process kill may leave hidden staging files.
There is no batch-wide transaction or power-loss durability promise.

## Architecture and compatibility

`input_loading`, `statistics`, `inspection`, `analysis`, `batch`, and `benchmark`
own their implementations. `_constants`, `_validation`, `_atomic`, and `progress`
provide shared contracts. `_core` re-exports historical names, forwards legacy
monkeypatch assignments, and delegates the old command interface to `_legacy_cli`.
Normal service execution does not import `_core`. Existing public module export
lists, legacy imports, CLI behavior, and artifact schemas remain supported.

## Recorded comparison

| Cells × genes | Case | v2.7 median (s) | v2.8 median (s) | Runtime change | Peak RSS change |
| --- | --- | ---: | ---: | ---: | ---: |
| 600 × 300 | csv | 0.0771 | 0.0839 | +8.76% | +0.01% |
| 600 × 300 | h5ad | 0.0792 | 0.0826 | +4.30% | +0.07% |
| 600 × 300 | scanpy | 0.0824 | 0.0886 | +7.56% | +0.01% |
| 600 × 300 | batch | 0.1561 | 0.1466 | -6.06% | -0.01% |
| 3,000 × 1,500 | csv | 1.3435 | 1.4711 | +9.50% | +2.36% |
| 3,000 × 1,500 | h5ad | 1.0109 | 1.1810 | +16.83% | -0.40% |
| 3,000 × 1,500 | scanpy | 0.5901 | 0.6003 | +1.74% | -0.02% |
| 3,000 × 1,500 | batch | 2.5547 | 2.5839 | +1.14% | -1.09% |

The initial large Native H5AD timing exceeded the threshold (+16.83%).
Three additional pairs of fresh processes alternated candidate/baseline order
with other local compute checks paused. Across 15 measured runs per version,
the pooled medians were 1.11597 s (v2.7) and 1.09566 s (v2.8), a **1.82%
improvement**. The slowdown did not reproduce under this controlled follow-up,
consistent with shared-host timing variation. Original results are retained
above, and all follow-up samples are in `h5ad-investigation.json`. No persistent
runtime regression was observed; all initial peak RSS changes were below 2.4%.
