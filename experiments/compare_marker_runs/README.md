# CellCoPilot prototype: compare existing marker exports

One static, export-only comparison experiment. Stable app, CLI, API, defaults and
export formats are unchanged. No source expression is needed to compare exports;
`generate.py` uses PBMC only to create the sample and independently check it.

## Run

From repository root, using the existing requirements/virtual environment:

```sh
mkdir -p /workspace/compare-marker-input
curl -fL https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad \
  -o /workspace/compare-marker-input/pbmc3k.h5ad
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mgf-mpl PYTHONPATH=. \
  .venv/bin/python experiments/compare_marker_runs/generate.py \
  /workspace/compare-marker-input/pbmc3k.h5ad
```

The public dataset is the same [PBMC source-selection QA file](../pbmc_h5ad_source_check/README.md),
SHA-256 `0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38`.
The generator fails on different bytes. Keep the H5AD outside Git. The recorded
run used Scanpy 1.12.5; both exports record the full scientific dependency versions.

The generator runs the public `analyze_file` API twice with input hashing enabled,
then extracts files from the existing `MarkerResult.to_bundle()` format. Both use
`.raw`, Scanpy, `louvain`, B cells versus rest (342 versus 2,296 cells), top_n=20,
rank_by=log2FC and otherwise default settings. **Only min_pct changes: 0.05 to
0.10** (5% to 10% target-cluster detection). No normalization is performed here.
The matrix is the file's selected raw.X; raw is storage, not a scale guarantee.

Compare any two current-format bundle ZIPs or extracted bundle directories:

```sh
.venv/bin/python experiments/compare_marker_runs/compare.py \
  experiments/compare_marker_runs/sample/A \
  experiments/compare_marker_runs/sample/B /tmp/marker-comparison
```

Open `/tmp/marker-comparison/report.html`, or the committed
[sample report](sample/report.html). The [no-change sample](sample/identical.html)
compares A with itself. Complete settings/provenance are expandable. The report
joins on (cluster, gene), distinguishes absent rows from zero detection, and lists
added/removed/reranked/shared markers with statistics and the relevant threshold.
`sample/A` and `sample/B` contain unmodified current-format marker CSV, settings,
config and provenance files; no new export schema is introduced. SHA-256s of the
input and marker CSVs establish which files are being compared.

## What this example established

| Change | Count |
| --- | ---: |
| Added in B | 7 |
| Removed from A | 7 |
| Shared, reranked | 12 |
| Shared, same rank | 1 |

All seven removed genes have exported A detection between 0.05 and 0.10. For
example MACROD2 (rank 2, pct_in=0.055556) and BACE2 (rank 12, pct_in=0.084795)
pass A's threshold and fail B's. Provenance confirms the same input, matrix,
engine version and dependencies, with min_pct the sole requested/effective
setting change. This supports a threshold-crossing explanation for these removals.
MS4A1 moves from rank 17 to 10 with unchanged statistics; CD79A moves from 5 to 4.
The report distinguishes these list-position changes from statistical changes.

Added markers have no A row. The report explicitly leaves their prior rank/filter
outcome unresolved, rather than pretending absence means a failed filter. A
separate full direct-Scanpy check found all seven passed A's filter but ranked
below its exported top 20 (recorded in sample/verification.json). That information
is intentionally NOT fed into the export-only report. Similarly, if a gene lacks
exported detection or comparable input/source/engine provenance, no causal removal
explanation is asserted. Settings outside the supported monotone top_n/min_pct cases are reported as
differences without assigning a cause.

## Verification and comparison with side-by-side files

The generator independently calls direct Scanpy rank_genes_groups with explicit
raw selection, Wilcoxon, tie correction and BH adjustment. It matches permissive
filters and ordering against both exported top-20 tables, checking detection,
adjusted p-values and fold changes. It verifies every positive threshold crossing,
all added markers' full A ranks, identical-export no-change, and refusal to infer
a cause with missing provenance or missing detection statistics. Assertions fail
on disagreement. Compact results and verification records are committed.

Opening the two CSVs and settings JSON files side by side exposes the same
information but requires manually finding the changed setting, joining genes,
tracking shifted ranks and comparing detection with thresholds. This sample report
puts all seven supported crossings, 12 rank movements and seven evidence gaps in
one view, with source/configuration checks alongside them. It made those checks
explicit and auditable; no reviewer time savings or user preference was measured.
It cannot recover truncated rows, infer preprocessing, or adjudicate biology.

## Extended checks: top N and multiple settings

The original single-filter comparison remains a regression check. New runs keep
B cells, raw, Scanpy, input hash and all other settings fixed:

| Comparison (A → B) | Added | Removed | Reranked | Explainable changed entries |
| --- | ---: | ---: | ---: | ---: |
| Original: min_pct .05 → .10, top20 fixed | 7 | 7 | 12 | 19/26 |
| Top N only: top20 → top10, min_pct .05 fixed | 0 | 10 | 0 | 10/10 |
| Combined: top20/.05 → top10/.10 | 0 | 10 | 9 | 19/19 |

Counts include rank changes, not just membership. Membership-only explanations
are 7/14, 10/10, and 10/10 respectively. The combined case attributes seven removals
to a failing detection threshold and three to truncation after filtering. For
example, MS4A1 passes the new .10 threshold but its post-filter rank is 10; FAM129C
passes too but moves to 11, so top10 truncates it. The nine shared rank movements
are established by counting surviving earlier entries in the exported A prefix.

This inference requires matching context, min_pct staying fixed or increasing,
only top_n/min_pct changing, requested/effective settings agreeing, contiguous
export ranks, matching expected B prefix, and unchanged shared statistics.
The code does not fill in unseen rows. In the original comparison, seven additions
still lack A evidence; their previous ranks remain unknown to the report.

The synthetic `extended/context_mismatch/` example keeps the CSVs but changes
B's declared input hash. It is explicitly a metadata-safety test, not a biological
run. The HTML exposes the mismatch and assigns **zero causal explanations**, even
though the values look like familiar threshold crossings. Additional automated
checks cover matrix, engine, engine-version and grouping mismatches, missing
metadata, and missing detection statistics. Differences in context are observations,
not grounds for attributing a gene's change to a threshold.

Generate all examples after the original sample:

```sh
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mgf-mpl PYTHONPATH=. \
  .venv/bin/python experiments/compare_marker_runs/extend.py \
  /workspace/compare-marker-input/pbmc3k.h5ad
.venv/bin/python -m pytest experiments/compare_marker_runs/test_compare.py -q
```

**One command to try with two existing bundle ZIPs:**

```sh
.venv/bin/python experiments/compare_marker_runs/compare.py \
  experiments/compare_marker_runs/extended/A.zip \
  experiments/compare_marker_runs/extended/combined/B.zip /tmp/marker-comparison
```

Substitute your own bundle paths. It writes `/tmp/marker-comparison/report.html`
and `/tmp/marker-comparison/comparison.json`; no dataset is opened or analysis run.
See [top-N report](extended/top_n/report.html),
[combined report](extended/combined/report.html), and
[context mismatch report](extended/context_mismatch/report.html).
Machine-readable summaries/checks: `extended/verification.json`.

**Recommendation: the evidence supports a small, read-only experimental Streamlit
integration for two uploaded bundles**, limited to inspecting settings/context and
showing supported explanations or explicit unknowns. Stable app/CLI code remains
untouched here. It makes ranking shifts, cutoff failures and context mismatch
checks auditable on this example; it does not establish user time savings or
broad scientific validity. Avoid changing analysis defaults or claiming arbitrary
causal explanations across changed datasets/settings.

The report does not mostly return unknown in these three examples. The remaining
gap is specific: a top-N export omits prior rows. The smallest useful future export
addition would be a per-gene filter outcome and post-filter/pre-top-N rank for
candidate genes outside the shortlist (plus the measured statistic for a claimed
failed threshold). This would resolve additions such as PNOC's missing A rank;
it would not justify comparing different biological contexts. No export change
is made in this experiment.

To regenerate outside the committed sample, add `--output /tmp/compare-marker-rerun`
to generate.py. To inspect reports in browsers that block file URLs, serve the
output with `python -m http.server 8767 --bind 127.0.0.1 --directory /tmp/marker-comparison`.


Recorded checks: all generator assertions passed; both static reports rendered in
Chromium with seven threshold explanations and seven explicit evidence gaps;
an independent rerun reproduced every sample export, report and verification file
byte for byte. No user study or timing benchmark was performed.

For an independent extended rerun, use `extend.py --output /tmp/compare-extended-rerun`.
All generated extended exports/reports can be compared byte for byte with the
committed examples. The standalone comparison command was also exercised on ZIPs
and checked for parity with extracted-directory inputs.

Extended verification completed: five export-only regressions passed via
`python -m pytest experiments/compare_marker_runs/test_compare.py -q`; the existing
repository suite passed 345 tests. Original generator checks were rerun, including
identical exports, missing provenance/statistics and full direct Scanpy validation.
The extended generator checked direct full results for both new shortlists and
all inferred post-filter ranks. All three new HTML reports rendered in Chromium.
A second extended run reproduced every file (including bundle ZIPs) byte for byte.
The documented two-ZIP command produced an identical machine-readable diff to the
sample. A future upload UI should stay experimental and add bounded archive
validation before accepting arbitrary public bundles.
