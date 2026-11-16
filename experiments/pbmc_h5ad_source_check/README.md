# PBMC H5AD source-selection QA

Independent of the annotation-evidence pilot. Base: latest main at
`bc86b18ddc14d059d37fccf83d3bcc53c371dbba`. Branch:
`experiment/pbmc-h5ad-source-check`. No analysis defaults changed.

## Source and reproduction

Public example from CZ CELLxGENE:
https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad

SHA-256: `0db367b991dd95809732b218539ede489bea99113807f62ebd7ccc970025fe38`.
The URL tracks main; verification pins the actual bytes and fails if they change.
Downloaded H5AD, screenshots, exports and large generated files remain outside Git.

Commands run from repository root in the Cloud workspace:

```sh
mkdir -p /workspace/pbmc-source-check
curl -fL https://raw.githubusercontent.com/chanzuckerberg/cellxgene/main/example-dataset/pbmc3k.h5ad \
  -o /workspace/pbmc-source-check/pbmc3k.h5ad
MPLCONFIGDIR=/tmp/mgf-mpl .venv/bin/python -m streamlit run app.py \
  --server.headless true --server.address 127.0.0.1 --server.port 8506 \
  --browser.gatherUsageStats false
```

In another terminal, with Python Playwright and `/usr/bin/chromium` available:

```sh
python experiments/pbmc_h5ad_source_check/browser.py /workspace/pbmc-source-check
MPLCONFIGDIR=/tmp/mgf-mpl OPENBLAS_NUM_THREADS=1 PYTHONPATH=. \
  .venv/bin/python experiments/pbmc_h5ad_source_check/verify.py /workspace/pbmc-source-check
MPLCONFIGDIR=/tmp/mgf-mpl XDG_CACHE_HOME=/tmp/mgf-cache \
  .venv/bin/python -m pytest -q -x
```

The browser script is the post-fix regression: it asserts the new warning, uploads
the actual file, selects louvain, tests X, chooses raw/Scanpy, keeps only B cells,
runs analysis, downloads four exports, and switches back to X to check clearing.
Cloud socket restrictions require approved execution outside the sandbox for the
server/browser/network commands; retain the configured network proxy and trust.
The static verification.json records this run's environment and results.

## Dataset inspection

AnnData opened the file successfully. There are 2,638 cells, eight louvain groups,
1,838 genes in X and 13,714 in raw.X, with no named layers. B cells: 342; rest:
2,296. X ranges from -2.8491054 to 10; raw.X from 0 to 6.0402546. The supplied
example is described as scaled X and log-normalized raw.X. Its uns does NOT
contain log1p preprocessing metadata; this QA verifies storage, dimensions and
values, not the original normalization recipe. Preprocessing state is unspecified
by metadata and must not be inferred merely from the storage name raw.

## First-time browser observations before the fix

Used local Chromium against the actual Streamlit app, not mocked uploads.

- Upload accepted. The initial cluster-column default was n_genes, giving many
  apparent groups; explicitly selecting louvain yielded the expected eight.
  This observation was recorded but not changed in this focused fix.
- Expression choices were AnnData.X (default) and AnnData.raw; no layers.
  The help correctly said source selection does not normalize/log-transform and
  that raw does not indicate preprocessing. No visible warning identified
  negative X before analysis.
- Leaving X and clicking Find marker genes produced only
  **“Gene-expression values cannot be negative.”** No source or scaled/centered
  explanation accompanied that error. This is the verified usability issue.
- With Scanpy and raw selected, the preview reported 13,714 genes. Selecting only
  B cells completed analysis and produced a five-row marker table, fold-change
  bar chart, expression dot plot and heatmap. The bar caption explained log2 fold
  change rather than significance. Dot/heatmap captions explicitly identified
  AnnData.raw, mean expression including zeros on the selected matrix's existing
  scale, and detection >0. They did not falsely claim raw meant unprocessed counts.
- CSV results, analysis settings JSON, bundle ZIP and expression-summary CSV
  downloads worked. The bundle also contains provenance, replay config and cluster
  summaries. Settings recorded expression_source=AnnData.raw.X, use_raw=true,
  layer=null, engine=scanpy, target_clusters=["B cells"]. Separate JSON and bundled
  settings matched. Switching back to X removed prior results and downloads.

## One verified fix

Added an upfront warning beside the expression selector whenever the selected
matrix contains negative values. It explains that the marker filters and
expression summaries assume nonnegative expression and count values >0 as
expressing, that scaling/centering can produce negatives, and
that a suitable matrix and preprocessing documentation are needed. It explicitly
avoids treating raw as a guarantee of preprocessing. It does not assert that all
negative matrices were scaled or all nonnegative matrices are normalized.

A cached helper inspects the chosen X/raw/layer in 256-row blocks using backed
AnnData; it does not densify an entire sparse matrix or transform values. The
normal validation/error remains in place. The warning disappears for this file's
nonnegative raw selection. No default selection, filter, ranking or statistical
method changed. No additional source/provenance fix was needed.

## Independent Scanpy, CLI and plot checks

Called rank_genes_groups directly on a freshly loaded file, explicitly use_raw=True,
B cells versus rest, method=wilcoxon, tie_correct=True,
corr_method=benjamini-hochberg, pts=True, n_genes=13714. Stored rankings in the
example file were not reused. Applied the app's permissive filters (nonnegative
log2FC, mean and detection differences; adjusted p<=1) and sorted log2FC descending,
gene ascending, top five. Statistics matched at rtol=1e-6, atol=1e-12, and gene
order matched exactly. CLI analyze with an explicit versioned raw/Scanpy config
produced the same marker table as the browser export.

| Marker | Scanpy log2FC | Detection in B cells | Detection in rest |
| --- | ---: | ---: | ---: |
| FAM177B | 24.974190 | 0.046784 | 0 |
| RP11-16E12.2 | 24.879610 | 0.043860 | 0 |
| PKHD1L1 | 23.768335 | 0.020468 | 0 |
| HRK | 23.544477 | 0.017544 | 0 |
| RP11-297B17.3 | 23.440454 | 0.014620 | 0 |

These are rare detections with zero outside detection. Scanpy's approximate fold
change and the app's permissive log2FC ranking place them above familiar B-cell
markers; their large fold changes are not a source mismatch or evidence of high
prevalence. This QA preserves the defaults rather than changing the biological
ranking question. P-values, adjusted p-values and detection fractions also
matched direct Scanpy for all five markers.

All 40 exported plot rows (five genes × eight groups) were checked against direct
raw.X means and fractions >0 and the existing summarize_expression API. Thus dot
area is a fraction/percentage of cells, color is mean selected expression including
zeros, and neither is a fold change or p-value. No re-normalization was applied.

## Post-fix verification and conclusion

Repeated the entire browser workflow after the change: the warning appeared
before analysis on X, raw/Scanpy analysis and all downloads succeeded, and source
switching cleared results. Re-ran direct Scanpy/CLI/summary verification on these
new downloads, checked ZIP integrity and identical bundled/standalone settings.
Added a small synthetic H5AD AppTest proving the warning is source-specific and
appears before results. Full suite: **345 passed**.

The result justified one small interpretation fix: earlier, source-specific
negative-input guidance. Analysis results, plot source use and export provenance
were already consistent. This is source-selection QA, not validation of the
file's biological annotation or a guarantee that any matrix named raw is suitable.
