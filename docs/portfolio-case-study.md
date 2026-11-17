# Making expression-source assumptions visible

**MarkerGeneFinderMini · Python, AnnData, Scanpy, Streamlit · v3.0.2**

I reproduced a real-data failure in the [hosted marker-gene app](https://cellcopilot-marker-finder.streamlit.app/)
and shipped a focused warning with reproducible scientific checks.

## The problem and decision

Uploading CZ CELLxGENE’s PBMC3k H5AD and selecting `louvain` exposed 2,638 cells
in eight groups. The file offered scaled `X` with 1,838 genes and `.raw.X` with
13,714 genes. Leaving `X` selected produced only “Gene-expression values cannot
be negative” after analysis started.

The app counts values above zero as expressing. After centering or scaling,
zero can represent a reference level rather than absence, making that detection
interpretation unsuitable. I added an upfront warning for negatives in the
selected matrix, checked in cached row blocks, without changing analysis defaults.
It presents scaling and centering as possible causes and directs users to suitable
expression and preprocessing documentation. `.raw` is a storage location whose
contents depend on the file. Missing preprocessing metadata remained unspecified.

## Verification across interfaces

- Ran `.raw`/Scanpy for **342 B cells versus 2,296 other cells**. Direct
  `rank_genes_groups`, explicitly selecting `.raw`, Wilcoxon, tie correction and
  Benjamini–Hochberg adjustment, matched all five exported markers’ statistics
  within documented tolerances and exact order after matching filters/ranking.
  CLI output reproduced the browser table.
- Checked **all 40 exported plot-summary rows** against direct `.raw.X` means
  including zeros and detection fractions above zero, verifying the source behind
  dot colors and sizes. Switching sources cleared previous results and downloads.
- Confirmed hosted **v3.0.2** and exported settings identifying
  `expression_source=AnnData.raw.X`, `use_raw=true`, and `engine=scanpy`.
  Hosted downloads passed the same Scanpy, CLI and summary checks.

The patch shipped in [PR #73](https://github.com/jec063/MarkerGeneFinderMini/pull/73)
and [v3.0.2](https://github.com/jec063/MarkerGeneFinderMini/releases/tag/v3.0.2),
with **345 passing tests** and Python **3.10–3.12 CI**. The
[QA report and reproduction scripts](../experiments/pbmc_h5ad_source_check/README.md)
record the dataset URL, SHA-256 and numerical checks. The outcome is clearer input
guidance and verified matrix-to-export provenance; user adoption and annotation
accuracy were not measured.
