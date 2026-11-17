# Compare marker exports: standalone Streamlit experiment

A read-only page for two existing MarkerGeneFinderMini export bundles. It shows
changed settings and recorded matrix/engine context, then additions, removals and
rank changes **within exported shortlists**. It runs separately from the stable
app; the stable app, CLI, API, export formats and analysis defaults are unchanged.
The page makes no annotation verdicts. Separate Community Cloud deployment is
prepared but blocked on account sign-in; no hosted demo URL is verified. See
[deployment settings and clean-install instructions](DEPLOYMENT.md).

## Try it without downloading expression data

From repository root, using the existing requirements/virtual environment:

```sh
PYTHONPATH=. MPLCONFIGDIR=/tmp/mgf-mpl .venv/bin/python -m streamlit run \
  experiments/compare_marker_runs/streamlit_app.py --server.port 8507
```

Open http://localhost:8507 and click **Try sample comparison** to load the Muraro
bundles immediately. Use **Comparison inputs** to switch between the sample and
**Upload my bundles**; changing inputs clears the previous results and cluster
selection. For uploads, choose an Earlier run and a Later run:

| Example | Earlier run | Later run | Expected change |
| --- | --- | --- | --- |
| PBMC3k: detection threshold | [extended/A.zip](extended/A.zip) | [sample/B.zip](sample/B.zip) | min_pct .05 → .10; 7 added, 7 removed, 12 reranked |
| PBMC3k: top N | [extended/A.zip](extended/A.zip) | [extended/top_n/B.zip](extended/top_n/B.zip) | top20 → top10; 10 removed |
| PBMC3k: both settings | [extended/A.zip](extended/A.zip) | [extended/combined/B.zip](extended/combined/B.zip) | 10 removed, 9 reranked |
| Muraro pancreas | [pancreas/results/earlier.zip](pancreas/results/earlier.zip) | [pancreas/results/later.zip](pancreas/results/later.zip) | min_pct .10 → .25; choose alpha, beta or duct |
| Identical runs | Any bundle above | The same bundle | No shortlist changes |

In Muraro, beta changes most: **28 entries** (9 added, 9 removed, 10 reranked),
versus duct 25 and alpha 23. Beta **ZNF385D** has detection 19.64%, so it fails the
new 25% threshold. **WSCD2** moves from rank 3 to 2 with unchanged statistics.
**GLP1R** appears at rank 12, but its earlier rank is unavailable in the exports;
the page leaves its cause unknown. Independent full results place it at earlier
rank 22; that information is deliberately not supplied to the page.

The overview distinguishes membership changes (added + removed) from all changed
entries (also reranked); neither measures biological change. The JSON download
contains all clusters; the selector only filters the view.

## Review map and retained artifacts

1. `compare.py`: existing export-only comparison and conservative explanation rules.
2. `upload.py`: bounded ZIP validation and adapter; `streamlit_app.py`: separate UI.
3. `test_compare.py`, `test_page.py`: 17 expression-data-free experiment tests, including sample/upload state changes.
4. [PBMC page verification](PAGE_VERIFICATION.md) and [Muraro QA](pancreas/QA.md):
   browser observations, independent checks, limitations and reproduction commands.

| Retained files | Reason |
| --- | --- |
| ZIPs linked above | Immediately uploadable examples; the Muraro ZIPs are the canonical complete exports |
| `sample/A`, `sample/B`, `extended/*/B` | Small readable PBMC CSV/settings/config/provenance fixtures used by tests, browser variants and side-by-side checks |
| `sample/verification.json`, `extended/verification.json` | Full direct-Scanpy checks, counts and hidden-rank validation for PBMC |
| `pancreas/results/independent_verification.json` | Muraro source hash, counts, preprocessing, direct-Scanpy checks and independent earlier ranks |
| `page-browser-verification.json`, `pancreas/browser_after.json` | One aggregate browser record per dataset; downloaded JSONs matched standalone comparisons |

Generated HTML reports and full comparison JSONs repeat evidence derivable from
the retained exports, so they are no longer tracked. Extracted Muraro exports
repeat the ZIP contents; unzip a bundle to inspect CSV/settings side by side.
Duplicate PBMC browser and pre-fix Muraro records were removed; the QA note retains
the pre-fix observations. ZIP copies of PBMC fixtures are intentionally retained
for immediate uploads and ZIP/directory parity checks. Regenerated outputs are
ignored locally; downloaded/derived H5AD files must stay outside Git.

## Generate reports from retained exports (no dataset required)

```sh
.venv/bin/python experiments/compare_marker_runs/compare.py \
  experiments/compare_marker_runs/extended/A.zip \
  experiments/compare_marker_runs/sample/B.zip /tmp/compare-pbmc
.venv/bin/python experiments/compare_marker_runs/compare.py \
  experiments/compare_marker_runs/pancreas/results/earlier.zip \
  experiments/compare_marker_runs/pancreas/results/later.zip /tmp/compare-muraro
```

Each command writes `report.html` and `comparison.json`. Substitute the top-N or
combined later bundle for those reports. Use the same bundle twice for the
no-change report. For the synthetic context-mismatch report, use `sample/A` and
`extended/context_mismatch/B`; its changed input hash must suppress all causes.
This is a metadata-safety check, not a biological run.

## Regenerate analyses and independent checks

The public PBMC file is the [source-selection QA dataset](../pbmc_h5ad_source_check/README.md),
SHA-256 `0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38`.
Both runs select `.raw`, Scanpy, B cells versus rest (342 versus 2,296 cells),
top20 and log2FC ranking. Only min_pct changes .05 → .10; no preprocessing is
performed. `.raw` is a storage location, not a scale guarantee.

```sh
mkdir -p /workspace/compare-marker-input
curl -fL https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad \
  -o /workspace/compare-marker-input/pbmc3k.h5ad
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mgf-mpl PYTHONPATH=. \
  .venv/bin/python experiments/compare_marker_runs/generate.py \
  /workspace/compare-marker-input/pbmc3k.h5ad --output /tmp/compare-sample
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mgf-mpl PYTHONPATH=. \
  .venv/bin/python experiments/compare_marker_runs/extend.py \
  /workspace/compare-marker-input/pbmc3k.h5ad --sample /tmp/compare-sample \
  --output /tmp/compare-extended
```

These commands regenerate the extracted exports, ZIPs, reports and verification
records, including the original identical-run check. Scientific dependency
versions are recorded in each export. The [Muraro QA](pancreas/QA.md) gives its
versioned dataset URL, checksum, explicit named-layer preparation and exact
`pancreas/run.py` command to regenerate all Muraro outputs. Both generators verify
shortlists and explanations independently against full direct Scanpy results.

## Tests and browser reproduction

```sh
MPLCONFIGDIR=/tmp/mgf-mpl .venv/bin/python -m pytest -v \
  experiments/compare_marker_runs/test_compare.py \
  experiments/compare_marker_runs/test_page.py
```

The original passing run [44063938885](https://github.com/jec063/MarkerGeneFinderMini/actions/runs/44063938885)
ran only 345 stable tests: `pytest.ini` defaults to `tests/`. The existing Actions
workflow explicitly invokes these two experiment files with verbose output
on Python 3.10, 3.11 and 3.12, separately from the stable suite and wheel checks.
It needs neither expression downloads nor browser installation.

For real Chromium upload checks, run the local server above, then use a Python
with Playwright and `/usr/bin/chromium` installed:

```sh
python experiments/compare_marker_runs/browser_demo.py \
  --url http://127.0.0.1:8507 --output /tmp/comparison-demo-browser
python experiments/compare_marker_runs/browser_page.py \
  --url http://127.0.0.1:8507 --output /tmp/compare-page-verification
python experiments/compare_marker_runs/pancreas/browser.py \
  --url http://127.0.0.1:8507 --expect-overview --output /tmp/compare-pancreas-browser
```

The page rejects invalid/incomplete bundles, clears stale results after upload
changes and handles identical/no-shared-cluster cases. It shows context differences
and leaves causes unavailable when metadata is missing or incompatible. Matching
recorded input hashes are not an independent rehash of the original dataset.

The examples made settings, threshold crossings, rank shifts and evidence gaps
easier to verify than manually joining CSV rows. They do not measure human time
savings, annotation accuracy or broad scientific validity. Top-N exports cannot
recover unseen prior rows. A future export extension could provide per-gene filter
outcomes and pre-top-N ranks, but this experiment changes no export schema.
