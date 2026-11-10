# v3 release validation

v3.0.0 promotes the validated RC2 implementation. The breaking change is removal of `src.*`,
the `_core` compatibility facade, and the flat CLI implementation. Supported
commands and Python services follow the contract frozen in v2.9. Statistical
defaults, result columns, engine capabilities, and data schemas are retained.

Profile version comparisons now understand PEP 440 release candidates, so a
profile from rc2 or stable 3.0.0 is correctly recognized as newer than rc1.

RC2 adds field-specific errors for malformed configuration types and invalid UTF-8.
Numeric strings, booleans used as numbers, fractional counts, and non-boolean
`use_raw` values are rejected. Valid configurations and scientific defaults are
unchanged. CLI report destinations cannot replace their source files, including
symlink aliases; batch plans also protect selected marker outputs and manifests.

## Required checks

1. Run `python -m pytest -q -x` on the source tree with removed modules physically
   absent. Tests assert that old modules cannot be found and flat arguments are
   rejected, and compare the public contract with the unchanged v2.9 snapshot.
2. Run the independent Native and Scanpy scientific fixtures at `rtol=1e-6`,
   `atol=1e-10`, preserving marker order, labels, counts, and ranks. Do not require
   the engines to agree with each other.
3. Replay v2.9 configs/settings/profiles, benchmark reports, and hashed batch
   manifests from `tests/fixtures/v2_9`. Verify/resume unchanged outputs and
   ensure unsupported schemas produce precise errors. Existing schema-1 settings
   and legacy-bundle reader tests remain in the suite.
4. Build distributions and run
   `python scripts/check_distribution.py dist/*.whl dist/*.tar.gz` to reject
   removed modules in wheel and source archives.
5. Install the wheel outside the checkout and run the copied
   `scripts/check_installed_wheel.py` with its app/data/example files. It verifies
   the installed package location and actual absence of removed modules, runs
   the CSV/API/CLI examples, H5AD X/layers/raw in both engines, a plugin, profile
   and bundle exports, planning, batch resume/verification, benchmark comparison,
   and Streamlit's example analysis. No import-blocking shim simulates removal.
6. Run the test/build/installed-wheel matrix on Python 3.10, 3.11, and 3.12 before
   publishing the prerelease. Repeat it on the merged release commit.

## Stable promotion

RC2 completed the [promotion validation](validation/v3/STATUS.md), including
real PBMC data and investigation of noisy performance measurements. Original
and follow-up samples are retained alongside the report. No scientific code or
defaults changed during promotion. Build the actual versioned v3.0.0 artifacts,
run all source and installed-wheel checks, and require passing PR and merged-main
CI before publishing stable v3.0.0 as the latest release.
