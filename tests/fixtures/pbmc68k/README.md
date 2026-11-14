# Biological regression references

These independent Native and Scanpy marker tables were generated with unmodified
v2.9.0 during v3 promotion validation. Never compare the engines against each
other. Row ordering is preserved and numeric tolerances are rtol=1e-6, atol=1e-10.

The test consumes Scanpy's bundled `10x_pbmc68k_reduced.h5ad`; it does not download
or redistribute the expression dataset. `provenance.json` records its source,
container and scientific-content SHA-256 fingerprints, dimensions, cell-type counts, reference version, and settings. The
original measurements are from the 10x Genomics PBMC68k dataset. Scanpy's reduced
fixture contains 700 cells and 765 genes. Tests use existing non-negative raw
expression, cast to float64, with no new normalization. Scaled X is not used.

Scanpy 1.11.5 and 1.12.6 have different HDF5 container bytes but identical raw
expression, cell IDs, gene IDs, and labels (verified element by element). The test
pins this ordered scientific content with a documented canonical encoding.

A scientific-content checksum change requires scientific review, not automatic reference updates.
Run `python -m pytest tests/test_biological_regression.py -q` to exercise both
engines, H5AD X/layer/raw, Native CSV, CLI/API agreement, exports and batch resume.
The normal CI matrix runs these tests without external dataset access.
