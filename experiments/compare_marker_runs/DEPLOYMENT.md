# Separate experimental comparison demo

**Working demo:** https://cellcopilot-compare-runs.streamlit.app/

Public sample and upload workflows verified on November 17, 2026. This is a
separate experimental app; the existing marker-analysis app is unchanged.

## Local launch and clean installation

From a checkout, create a fresh Python 3.12 environment (3.10–3.12 are tested):

```sh
python -m venv /tmp/comparison-demo-venv
/tmp/comparison-demo-venv/bin/python -m pip install -r requirements.txt
/tmp/comparison-demo-venv/bin/python -m streamlit run \
  experiments/compare_marker_runs/streamlit_app.py --server.port 8510
```

Open http://localhost:8510 and select **Try sample comparison**. The Muraro bundles
are committed and loaded relative to the entrypoint, without an expression data
download or analysis rerun. The sample uses the same bounded upload validation and
comparison logic as **Upload my bundles**. Changing modes or uploads clears the
previous comparison and cluster selection. The JSON download covers all clusters.

The entrypoint also works from a different working directory using an absolute
script path; no editable package install or `PYTHONPATH` is required. Keep the
repository checkout intact, including `marker_gene_finder/` and the two ZIPs in
`experiments/compare_marker_runs/pancreas/results/`. The root `requirements.txt`
is the deployment dependency source; no extra service or secret is required by
the comparison page. Do not repoint the existing marker-analysis app at this file.

## Streamlit Community Cloud configuration

Deployment target settings (use these for recreating the separate app):

| Field | Value |
| --- | --- |
| Repository | `jec063/MarkerGeneFinderMini` |
| Branch | `main` |
| Main file path | `experiments/compare_marker_runs/streamlit_app.py` |
| Python version | `3.12` |
| App URL | https://cellcopilot-compare-runs.streamlit.app/ |
| Dependencies | Root `requirements.txt` |
| Secrets | None required by this demo |

**Deployment status:** published and publicly verified. The earlier account-access
blocker was resolved by the owner deploying the app. Public verification requires
no Community Cloud account. Private management settings and the server's exact Git
revision are not exposed by this check; results are compared against repository
reference `504cb5174cf496d4536e9233998c2152d6fb9b28` and its committed Muraro bundles.
Keep the existing `cellcopilot-marker-finder.streamlit.app` app's settings unchanged.

## Browser acceptance check

With Python Playwright and `/usr/bin/chromium` available:

```sh
python experiments/compare_marker_runs/browser_demo.py \
  --url https://cellcopilot-compare-runs.streamlit.app/ --output /tmp/comparison-demo-browser
python experiments/compare_marker_runs/browser_page.py \
  --url https://cellcopilot-compare-runs.streamlit.app/ --output /tmp/comparison-upload-browser
```

For local checks, substitute `http://127.0.0.1:8510`. The iframe-aware runners locate
the rendered application inside Community Cloud and wait for expected result state
instead of relying on local response times. The first check loads the sample,
selects beta and duct, verifies settings/overview/counts, compares JSON downloads
with the same uploaded Muraro bundles, switches input modes, rejects an invalid
replacement without showing stale sample results, and restores the sample. The
second exercises the original PBMC comparisons, incompatible contexts, incomplete
bundles, missing dataset identity and stale-result handling. Outputs and screenshots
remain outside Git.

Expected Muraro evidence: min_pct 10% → 25%; beta most changed at 28 entries.
ZNF385D's 19.64% detection fails 25%; GLP1R's earlier rank remains unavailable
because the earlier export lacks its row. These checks demonstrate evidence for
exported shortlist changes, not biological differences or measured researcher time
savings. Those usability benefits remain pending.

## Verified in the cloud workspace

A fresh virtual environment installed only root `requirements.txt`. From `/tmp`,
with `PYTHONPATH` unset, the absolute entrypoint started successfully. In that
installation, **17 experiment tests and 345 stable tests passed**. Chromium verified
sample beta/duct selection, setting/overview counts, all four sample/upload diff
downloads, mode changes clearing visible results, invalid uploads not reusing the
sample, and the sample action resetting the cluster after an invalid upload.
Original PBMC browser checks also passed, including context mismatch, missing
identity, identical runs, incomplete bundles and stale-result clearing.

Compact results: [demo-verification.json](demo-verification.json). The sample
JSON downloads match the export-only comparator exactly. Sample files, statistical
logic and export formats are unchanged. This original local record predates public
deployment; the current hosted evidence is recorded below. Researcher time savings
remain unmeasured.

## Hosted verification

The normal public URL opened the standalone comparison page. Selecting **Try sample
comparison** showed min_pct **0.10 → 0.25**, the same named `log1p_cp10k` layer and
Scanpy, and beta as the most changed cluster: **28 exported entries** (9 added,
9 removed, 10 reranked), versus duct 25 and alpha 23.

- Beta **ZNF385D**: earlier rank 2, detection **0.1964286 = 19.64%**, below the new
  **25%** threshold. Removal is from the exported shortlist, not absence of expression.
- **GLP1R**: later rank 12; earlier rank **Unavailable** because the earlier export
  contains no row. Its earlier rank or filter outcome cannot be inferred from these
  exports. No full-analysis rank is inserted into the page.

Selected beta and duct in the sample, then uploaded the two committed Muraro ZIPs
through **Upload my bundles** and selected beta. Replaced an upload with invalid
bytes, confirmed error feedback and
removal of stale results/download controls, and returned to a fresh sample. Both
complete hosted runs passed; all **eight Muraro JSON downloads** matched the same
committed comparison exactly. Original PBMC checks also passed for min_pct-only,
top-N-only, combined settings, identical runs, mismatched contexts, incomplete
bundles, missing identity and no shared clusters.

[Fix PR #78](https://github.com/jec063/MarkerGeneFinderMini/pull/78) explicitly resets
the selector for new comparisons after an intermittent recovery retained the old
cluster. It also makes the browser checks handle Cloud's iframe and await settled
results. The hosted checks were repeated after the fix merged. Comparison logic,
statistics and exported evidence did not change.

Compact evidence: [hosted-verification.json](hosted-verification.json), including the
repository reference, canonical diff hash and CI links. Downloaded diffs and browser
screenshots remain outside Git. This is an agent-operated functional check, not a
researcher usability study; claims about researcher time savings remain pending.
