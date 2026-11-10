# RC2 promotion validation — passed after investigation

Date: 2026-11-10. Candidate: published v3.0.0rc2, commit `0d11acda1a3758450d8dd8688ec5b4f1ce0bbe9c`. Reference: v2.9.0.

## Decision

RC2 is ready for stable release preparation. No reproducible application regression was found. The validation harness errors and network access blocker were resolved. Performance flags were investigated, with original failures retained below. No application-code fix was warranted. The final versioned v3.0.0 build still needs the normal release build and CI checks before publication; this assessment validates RC2, not an unbuilt stable artifact.

## Correctness, installation, migration, and compatibility

- All 342 repository tests passed; dependency consistency passed. The suite includes independent Native and Scanpy scientific fixtures, the frozen public contract, removed-interface checks, plugin workflows, legacy saved artifacts, atomic-write interruption checks, and migration examples.
- Downloaded published wheel and source archive passed archive validation. SHA-256 digests match the original tested release (environment.json).
- The published wheel was installed outside the checkout. Migration workflows, supported CLI, public Python API, reference plugin, H5AD layers/raw, profile/bundle exports, batch resume/verification, and Streamlit example passed (wheel-smoke.log).
- Original release CI verified successful on Python 3.10, 3.11 and 3.12, including builds and installed-wheel workflows (ci.json). This is verified existing CI, not a newly triggered matrix.
- Reviewed docs/migration.md against executable examples and compatibility tests; no missing documented replacement was found.

## Real biological data

Used Scanpy's bundled 10x Genomics PBMC68k reduced dataset: 700 cells, 765 genes, 10 annotated cell types. Dataset source URL and SHA-256 are recorded in biological-report.json. Existing non-negative raw expression was cast to float64; scaled negative X was not used and no new normalization was applied.

Native CSV and Native/Scanpy H5AD X/layer/raw produced 50 markers per analysis and matched their own v2.9.0 reference results at rtol=1e-6, atol=1e-10, preserving ordering. CLI/API agreement, configuration round trips, ZIP integrity, batch resume, and manifest verification passed for both engines. The initial harness incorrectly requested unsupported Scanpy CSV; it now explicitly tests rejection and restricts Scanpy batch inputs to H5AD. These are workflow/regression checks on real measurements, not a claim of biological discovery accuracy.

## Performance investigation

The initial fresh-process test flagged small Native H5AD at +101.58% elapsed time. Within-version timings varied substantially. The executed analysis, input-loading, statistics, batch, and engine modules are unchanged between versions.

A follow-up used persistent separate workers, three warmups, 21 individual measurements per version/case, alternating execution order, one thread, and identical dependencies. Garbage collection was requested before each timed sample, while normal automatic collection remained enabled. Both elapsed and process CPU time were retained. RSS is the process high-water mark, including input generation/dependencies. No other local compute validation ran concurrently.

| Cells × genes | Case | Elapsed change | CPU change | Peak RSS change |
| --- | --- | ---: | ---: | ---: |
| 600 × 300 | csv | +22.43% | +22.43% | +0.00% |
| 600 × 300 | h5ad | +8.54% | +8.54% | +0.04% |
| 600 × 300 | scanpy | -20.61% | -20.59% | -0.02% |
| 600 × 300 | batch | +5.14% | +4.58% | +0.00% |
| 3000 × 1500 | csv | +1.48% | +1.62% | +0.31% |
| 3000 × 1500 | h5ad | -3.65% | -3.48% | +0.73% |
| 3000 × 1500 | scanpy | -7.51% | -6.81% | -0.00% |
| 3000 × 1500 | batch | -2.69% | -2.87% | +0.50% |

All large cases passed the +10% investigation threshold. Small H5AD passed, but small CSV flagged +22.43%, inconsistent with its initial +2.49%. To investigate process variability, both small cases were measured again with PYTHONHASHSEED=0 and both workers pinned to the same available CPU, retaining alternating individual samples.

| Small case | Fixed-seed/CPU elapsed change | CPU change | Peak RSS change |
| --- | ---: | ---: | ---: |
| csv | -5.22% | -5.22% | +0.00% |
| h5ad | -7.51% | -7.51% | -0.00% |

Both passed. An A/A control then ran v2.9.0 in both worker slots with the same protocol; elapsed differences were CSV -6.99% and H5AD -2.70%, demonstrating variation even with identical code. The original 102% slowdown did not reproduce in either H5AD follow-up or the large workload. The evidence supports benchmark/process variability rather than a persistent version regression; it does not identify a unique host-level cause or promise speedups. All initial and follow-up samples are retained.

## Reproduction and artifacts

- biological.py: run with baseline or candidate; select baseline source or installed wheel via PYTHONPATH.
- execution_core.py, performance.py: original benchmark and flagged measurements.
- controlled_worker.py, controlled.py: alternating measurements, all eight cases.
- fixed.py: fixed seed and CPU; MGF_AA=1 selects the identical-code control.
- initial-performance.json, controlled-performance.json, fixed-performance.json, aa-performance.json: raw samples.
- environment.json: dependency versions and published artifact hashes.

Application repository is unchanged. No stable tag, commit, PR, or release was created during validation.

This report records the completed RC2 assessment before the stable version bump.
Executable validation scripts and logs are retained in the release evidence archive.
