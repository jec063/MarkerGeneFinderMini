# Multi-cluster pancreas workflow check

Starting commit: `66bb1658f08cf16ef8aac9f59c84f96a478bd4af`.
**Outcome:** the page answers the three questions correctly after two small
presentation fixes. This supports a draft PR for the read-only experimental page,
not a claim of measured user productivity or broader scientific validity.

## Dataset and two realistic runs

Public Muraro et al. human pancreas atlas, Cell Systems 2016,
https://doi.org/10.1016/j.cels.2016.09.002; GEO GSE85241. CELLxGENE version URL:
https://datasets.cellxgene.cziscience.com/ac56150b-add4-4336-9059-6d3d3ce17f3b.h5ad
SHA-256: `183673651cfa8c473a26641d42011d43be44eb2fea44e6e6ab8e2b0065d07483`.

2,126 cells and 15,643 features; published `cellular_classification` targets
alpha (812), beta (448), and duct (246). Every source cell, including other cell
types, remains in each target-versus-rest comparison. We normalize the source's
nonnegative integer raw counts to 10,000 per cell and log1p, without scaling or
imputation, and store that explicitly prepared expression in `log1p_cp10k`.
Both bundles use that same named layer, Scanpy, top20 ranked by log2FC, and default
filters except min_pct: **0.10 → 0.25**. This is a plausible exploration of how
requiring broader target-cell detection changes marker candidates. The analysis
uses the existing public engine; no statistics or export format is modified.

## Answers checked in the actual page

| Annotated cluster | Added | Removed | Reranked | Membership changes | Total changed entries |
| --- | ---: | ---: | ---: | ---: | ---: |
| beta | 9 | 9 | 10 | 18 | 28 |
| duct | 5 | 5 | 15 | 10 | 25 |
| alpha | 4 | 4 | 15 | 8 | 23 |

1. **Setting and largest shortlist change:** min_pct changed from 10% to 25%.
   Beta changed most under both membership-only and membership-plus-rank counts.
   These are changes in exported entries, not evidence that beta expression
   changed biologically between runs.
2. **Why a marker disappeared/moved:** beta ZNF385D, earlier rank 2, has detection
   0.1964286 (88/448 cells), so it fails the new 0.25 threshold. Its log2FC is
   6.977443. WSCD2 retains detection 0.5691964 and log2FC 6.9106436 but moves
   from rank 3 to 2 after the earlier entry is filtered out. CSV/settings checks
   confirmed those exact values and that all other settings match.
3. **What remains unknown:** GLP1R appears at beta rank 12 but has no earlier
   exported row. The page correctly says its earlier rank/cause is unavailable,
   not that it was absent or failed the earlier threshold. Separate full Scanpy
   validation finds that it passed earlier filters at rank 22; this extra evidence
   is not passed into the comparator. Across all groups, 18 additions remain
   unresolved from exports. All 18 removals and 40 rank movements are supported.

## Side-by-side check and interpretation findings

Opening the two CSVs and settings JSON files side by side gives the same evidence,
but finding beta's largest change required joining genes and counting rank shifts
for three groups. The overview now gives the 28/25/23 comparison directly, and
ZNF385D's old detection sits next to the new threshold and an explicit missing
later rank. This made the cutoff check and rank-versus-statistic distinction
concrete to verify; no researcher timing study was performed.

The pre-fix page showed only one cluster's counts at a time. It did not make a
false claim, but its initial alpha view (23 changed entries) did not answer which
cluster changed most without manually recording each selection. Added a sorted
shared-cluster overview that distinguishes membership changes from all changed
entries. Also, the matrix display previously said only AnnData.layers; the layer
name was available only in expanded settings. It now displays `log1p_cp10k`
directly. Neither fix alters the comparator or underlying diff.

No incorrect data interpretation or causal explanation was observed in this
workflow. The page explicitly retained unknowns for additions rather than implying
new expression. These are observations from one dataset and an agent-operated
workflow check, not evidence of adoption, annotation accuracy or human speed gains.

## Verification and reproduction

`run.py` checks the source checksum/counts and direct full Scanpy Wilcoxon results
(tie correction, BH adjustment, same layer/targets/rest), applies matching filters
and ranking, and verifies both exported shortlists and every reported cutoff/rank
explanation. It also verifies the true earlier ranks of additions independently.
Checks passed: **14 experimental tests and 345 stable tests**. The original PBMC
min_pct-only, top-N-only, combined, identical-run and mismatched-context browser
checks passed, including invalid bundles, missing identity, no shared clusters and
stale-result clearing. Compact browser evidence is in `browser_before.json`,
`browser_after.json` and `pbmc_browser_regression.json`.

Current CSVs/config/settings/provenance and small bundle ZIPs are in `results/`;
raw and derived H5AD files stay outside Git. Browser JSON downloads equal the
standalone full diff for each of the three selected clusters, before/after the
presentation fixes. A fresh local server was used after changing an imported
helper because the prior Streamlit process retained that module in memory. Browser
automation initially timed out selecting an off-screen custom dropdown; explicitly
scrolling the selector into view before opening it resolved the harness issue.

From repository root, with existing requirements installed:

```sh
mkdir -p /workspace/compare-pancreas
curl -fL https://datasets.cellxgene.cziscience.com/ac56150b-add4-4336-9059-6d3d3ce17f3b.h5ad \
  -o /workspace/compare-pancreas/muraro.h5ad
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mgf-mpl PYTHONPATH=. \
  .venv/bin/python experiments/compare_marker_runs/pancreas/run.py \
  /workspace/compare-pancreas/muraro.h5ad
PYTHONPATH=. MPLCONFIGDIR=/tmp/mgf-mpl .venv/bin/python -m streamlit run \
  experiments/compare_marker_runs/streamlit_app.py --server.port 8509
```

Upload `results/earlier.zip` and `results/later.zip`. For automated browser checks,
use system Python with Playwright and `/usr/bin/chromium`:

```sh
python experiments/compare_marker_runs/pancreas/browser.py \
  --url http://127.0.0.1:8509 --expect-overview --output /tmp/compare-pancreas-after
MPLCONFIGDIR=/tmp/mgf-mpl .venv/bin/python -m pytest \
  experiments/compare_marker_runs/test_compare.py \
  experiments/compare_marker_runs/test_page.py -q
```

The dataset source and commands are sufficient to regenerate the bundles; the
derived input hash can depend on HDF5 serialization, so use the recorded source
hash and compare scientific results rather than promising identical file bytes
across environments. The page remains local and undeployed; stable app/CLI and
main are untouched.
