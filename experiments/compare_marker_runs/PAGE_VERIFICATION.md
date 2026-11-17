# Read-only Streamlit experiment verification

Starting reference: `b6eeee47b0a2c703319197cb74040490b608441c`.
The separate `streamlit_app.py` entrypoint runs locally only; it is not part of
the stable app's navigation and has not been deployed. The comparator, statistical
engines, CLI and export formats are unchanged.

## Browser observations

Chromium uploaded actual current-format ZIPs created from the committed PBMC
examples. No expression dataset was uploaded or analysis rerun by the page.

| Example | Added | Removed | Reranked | Unchanged | JSON parity |
| --- | ---: | ---: | ---: | ---: | --- |
| Original min_pct | 7 | 7 | 12 | 1 | Exact standalone diff |
| Top N only | 0 | 10 | 0 | 10 | Exact standalone diff |
| Top N + min_pct | 0 | 10 | 9 | 1 | Exact standalone diff |
| Identical runs | 0 | 0 | 0 | 20 | Exact standalone diff |
| Synthetic context mismatch | 7 | 7 | 12 | 1 | Exact standalone diff; no causes assigned |

Changed settings appear before cluster counts. The shared-cluster selector shows
B cells. Added/removed/reranked tables use plain-language evidence labels; missing
ranks are explicitly unavailable. The original seven additions remain unexplained.
The combined case maps the existing engine's seven threshold crossings, three
truncations and nine rank shifts without adding any new causal inference.

The identical case shows no-change feedback, including when equivalent exports
are repacked with different ZIP timestamps. The context-mismatch case shows the
recorded input-hash difference and suppresses all causes. Missing identity warns
that dataset identity cannot be verified; matching hashes are described as recorded
claims, not an independent inspection of the original dataset.

Invalid ZIP and incomplete bundle uploads showed errors, with prior results and
download controls removed. No shared exported clusters showed clear feedback and
no cluster selector. Replacing either upload with invalid content cleared stale
results. All five downloaded JSON files matched the standalone comparator exactly;
the cluster selector does not alter the all-cluster download.

## Validation and limits

`test_compare.py` and `test_page.py` verify upload/standalone parity, source and
context suppression, bounded flat-archive validation, source-specific evidence
labels, identical/no-shared cases, and stale state clearing. All 13 experiment tests passed; the existing repository suite passed 345 tests. The original export-only generation/Scanpy evidence remains the
basis for explanations; no statistical algorithm or expression values changed.

The upload adapter allows <=20 MiB compressed bytes, <=50 MiB declared expanded
size, <=32 flat entries and <=100,000 marker rows. It rejects duplicate filenames,
nested paths, encryption and symlinks; only known CSV/JSON members are written to
a temporary directory. It reuses the existing loader and typed config validation,
checks numeric columns, and does not retain upload files after comparison. This
is an experimental adapter, not a general archive-processing service.

Differences from the standalone report are presentation only: Earlier/Later
labels, plain-language evidence, shared-cluster filtering, and upload feedback.
Unsupported/incomplete bundles are rejected before comparison. No user time savings,
biological accuracy, public deployment or production security guarantee is claimed.

Aggregate browser checks are recorded in `page-browser-verification.json`.
Screenshots and downloaded diffs remain in /tmp and are not committed.
