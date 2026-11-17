# Separate experimental comparison demo

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

Create a **new app** from https://share.streamlit.io/ with these settings after
merging the PR:

| Field | Value |
| --- | --- |
| Repository | `jec063/MarkerGeneFinderMini` |
| Branch | `main` |
| Main file path | `experiments/compare_marker_runs/streamlit_app.py` |
| Python version | `3.12` |
| App URL slug | Request `cellcopilot-compare-runs`, or another available separate slug |
| Dependencies | Root `requirements.txt` |
| Secrets | None required by this demo |

**Deployment status:** prepared, not published. In the Codex workspace,
Community Cloud opens at **Continue to sign-in**. No authenticated Chrome profile
or Streamlit deployment credential is configured. GitHub repository access does
not provide a Community Cloud login. The exact remaining action is for an account
with deployment access to sign in, create the new app using the settings above,
and choose an available URL slug. No working demo URL has been verified; the
requested slug is not a claim that the URL exists. Keep the existing
`cellcopilot-marker-finder.streamlit.app` app's settings unchanged.

Once deployed, record the actual working URL here and in the experiment README,
then repeat the checks below against that URL. Hosted sample loading, upload,
cluster selection and download checks remain unverified until account access is
available. The local server is not a public deployment.

## Browser acceptance check

With Python Playwright and `/usr/bin/chromium` available:

```sh
python experiments/compare_marker_runs/browser_demo.py \
  --url http://127.0.0.1:8510 --output /tmp/comparison-demo-browser
python experiments/compare_marker_runs/browser_page.py \
  --url http://127.0.0.1:8510 --output /tmp/comparison-upload-browser
```

Substitute the hosted HTTPS URL after deployment. The first check loads the sample,
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
logic and export formats are unchanged. Hosted verification is still blocked as
described above; no working hosted URL or researcher time saving is claimed.
