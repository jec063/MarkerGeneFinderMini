# v3.0.0rc1 validation

This is a prerelease, not stable v3. The breaking change is removal of `src.*`,
the `_core` compatibility facade, and the flat CLI implementation. Supported
commands and Python services follow the contract frozen in v2.9. Statistical
defaults, result columns, engine capabilities, and data schemas are retained.

Profile version comparisons now understand PEP 440 release candidates, so a
profile from rc2 or stable 3.0.0 is correctly recognized as newer than rc1.

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

## Promotion to stable

Keep v2.9.0 as the latest stable release while rc1 is exercised. Review this
checklist, CI results, and any RC feedback before releasing stable 3.0.0. Changes
after rc1 require the affected checks to be rerun. No new statistical methods,
distributed execution, or platform expansion are included in this transition.
