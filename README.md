# Marker Gene Finder Mini

A small Python project for identifying candidate marker genes from clustered gene-expression data. This is the first working component of the broader CellCoPilot project.

Current application version: **2.6.0**. The version is shown in the web
interface, available from the CLI with `--version`, and recorded in new
analysis-settings exports and reproducibility bundles.

**Live app:** [Launch the CellCoPilot Marker Gene Finder](https://cellcopilot-marker-finder.streamlit.app/)

The program compares gene expression inside each cluster with all other cells and ranks genes using log2 fold change.

## Input

The program accepts `.csv` and `.h5ad` expression files.

A CSV file contains:

- one row per cell
- an optional cell-identifier column
- a cluster-label column
- numeric gene-expression columns

See `data/example_expression.csv` for an example.

An H5AD file uses:

- `AnnData.X` by default, a named layer, or `AnnData.raw.X` for the cell-by-gene expression matrix
- `AnnData.obs_names` for cell identifiers
- a selected `AnnData.obs` column for cluster labels
- `AnnData.var_names` for X and layers, or `AnnData.raw.var_names` when selecting raw

Dense and sparse expression matrices are supported. The native
engine converts the expression matrix to an in-memory dense table.
The Scanpy engine works directly with the AnnData object and preserves
sparse matrices, so it is recommended for larger sparse H5AD datasets.

Before the Native engine materializes an H5AD expression matrix, it estimates
the dense matrix size from its shape and numeric data type. Conversion stops
when the estimate exceeds the default 2,048 MiB safety limit, avoiding a large
unexpected allocation. Use Scanpy to preserve sparse storage, change the CLI
limit with `--max-dense-memory-mb`, or explicitly pass `0` to disable the
guard. CSV inputs and Scanpy analysis are not subject to this dense-conversion
limit. Negative, non-finite, Boolean, and non-numeric limits are rejected.

The Native engine also calculates gene statistics in chunks of 1,000 genes to
limit temporary per-cluster arrays. Use `--gene-chunk-size` to tune this for
the available memory. Benjamini-Hochberg correction is still applied across
all tested genes in a cluster, so changing the chunk size does not change
marker values, ordering, or adjusted p-values.

Command-line analyses report completed clusters to standard error, including
the current and total cluster counts. This keeps CSV output and diagnostics
machine-readable while making long Native or Scanpy runs observable. Pass
`--quiet` to suppress progress messages. Library callers can supply a
`progress_callback` to `find_markers()`, `find_markers_scanpy()`, or `run()`;
the callback receives the completed count, total count, and cluster label.

Native command-line analysis can process independent target clusters
concurrently with `--workers`. The default is one worker, which preserves the
lowest peak-memory behavior. Multiple workers preserve cluster ordering and
marker results, but each active worker needs its own temporary gene-chunk
arrays, so increase the value only when additional memory is available. Scanpy
already performs its ranking in one engine call and therefore requires one
worker.

Before committing to an analysis, use `--inspect-json` to scan an input and
write a machine-readable preflight report without calculating markers. The
report includes the application version, expression source, cell, gene, and
cluster counts, storage type, estimated dense expression size, configured
Native safety limit, and an operational engine recommendation with its
reasons. Sparse H5AD inputs and H5AD inputs that exceed the Native limit
recommend Scanpy; CSV and dense H5AD inputs within the limit recommend Native.
This is resource guidance, not a claim that the engines use identical
statistics. H5AD inspection opens the file in backed mode; CSV inspection reads bounded
10,000-row chunks rather than loading the entire table solely for preflight.
Inspection reports normalize the configured memory limit to a regular JSON
number, including when library callers provide a NumPy numeric scalar.

Use `--plan-json` to turn the same preflight data into a versioned execution
plan without calculating markers. Plans keep `requested_engine`,
`recommended_engine`, and `resolved_engine` separate, record resolved Native
chunk and worker settings, and report conflicts through `can_execute` and
`warnings`. Explicit engine choices are never silently replaced.

The opt-in `--engine auto` mode follows the preflight resource recommendation
for each input: CSV uses Native, sparse H5AD uses Scanpy, and dense H5AD uses
Native unless its estimated dense representation exceeds the configured
safety limit. The default remains `native`, and automatic selection is
resource-based rather than a guarantee of statistically identical results.

## Analysis engines

MarkerGeneFinderMini provides two analysis engines:

- **Native:** the original CellCoPilot implementation for CSV and H5AD input.
- **Scanpy:** an H5AD-only engine using `scanpy.tl.rank_genes_groups()` with the Wilcoxon method and Benjamini-Hochberg correction.

The native engine remains the default. Scanpy works directly with AnnData and reports its own approximate log2 fold changes.

Select one or more **Clusters to analyze** in the web interface to limit the
reported groups. Each selected cluster is still compared with every other cell.
On the CLI, repeat `--target-cluster`, for example
`--target-cluster 0 --target-cluster 2`. Omitting it analyzes every cluster.
Set **Minimum cells per target cluster** or use `--min-cluster-size` to avoid
testing undersized groups. Small clusters remain part of the one-versus-rest
reference population; they are only excluded as analysis targets. The default
of `1` preserves the previous behavior.

Exact genes and gene-name prefixes can be excluded before testing. In the web
interface, enter comma-separated values under **Exclude genes** or **Exclude
gene prefixes**. On the CLI, repeat `--exclude-gene` or `--exclude-prefix`.
No exclusions are applied by default; organism-specific naming conventions
remain under the user's control.

The web interface provides **Permissive**, **Balanced**, and **High
specificity** filter presets. A preset supplies starting values for the core
expression-fraction, fold-change, and adjusted-p-value filters without changing
cluster selection, gene exclusions, or the remaining advanced controls. Editing
one of those core thresholds marks the configuration as **Custom**. The selected
preset name and the resulting individual thresholds are retained in exported
analysis settings.

## Method

For each cluster and gene, the program calculates:

1. Number of cells inside the target cluster and in the reference population
2. Mean expression inside the cluster
3. Mean expression outside the cluster
4. Difference between target-cluster and reference mean expression
5. Fraction of cells expressing the gene inside and outside the cluster
6. Number of cells expressing the gene inside and outside the cluster
7. Difference between the inside and outside expression fractions
8. Log2 fold change
9. One-sided Wilcoxon-Mann-Whitney p-value
10. Benjamini-Hochberg adjusted p-value

```text
log2FC = log2(
    (cluster_mean + pseudocount)
    /
    (other_mean + pseudocount)
)
```

The native engine default pseudocount is `0.1`.

The `cluster_size` and `reference_size` result columns record the sample sizes
used for each one-versus-rest comparison. They are included in the results CSV
and analysis bundle so statistical context remains attached to every marker.

By default, genes with negative log2 fold change are excluded. A higher minimum log2 fold-change threshold can be used to require stronger enrichment.

Statistical testing uses a one-sided Mann-Whitney U test, which is the Wilcoxon rank-sum test for independent samples. The alternative hypothesis is that expression is higher inside the cluster than outside it. For each cluster, p-values for all genes are adjusted together using the Benjamini-Hochberg false-discovery-rate procedure.

Results are ranked after applying the selected thresholds. Choose log2 fold
change (the default), target-cluster mean expression, target-cluster expressing
cell count, outside-cluster expressing-cell count, mean-expression difference,
expression-fraction difference, or adjusted p-value in the web interface or
with `--rank-by`. Ranking by
`cells_expressing_in` prioritizes markers supported by the largest absolute
number of target cells. Ranking by `cells_expressing_out` prioritizes genes
detected in the fewest reference-population cells. The existing fold-change
and gene-name tie-breakers keep both orderings deterministic. The output
includes raw p-values in `p_value` and adjusted p-values in `p_adj`.
The one-based `rank` column records each gene's position within its cluster
after filtering, so ranks restart at 1 for every cluster.
Use the minimum expression-fraction difference to require markers to be
expressed in a larger share of target-cluster cells than other cells.
Use the maximum outside-cluster expression fraction to exclude broadly
expressed genes even when they are enriched in the target cluster. Its default
of `1.0` preserves the previous behavior.
Use the minimum target-cluster expressing-cell count to require evidence from
a specific number of cells in addition to the fractional expression threshold.
Set **Minimum target-cluster cells expressing gene** in the web interface or
use `--min-cells-expressing` on the CLI. Its default of `0` preserves the
previous behavior.
Use the minimum target-cluster mean-expression threshold to remove genes whose
average expression inside the target cluster is too low on the input matrix's
existing scale. Set **Minimum target-cluster mean expression** in the web
interface or use `--min-cluster-mean` on the CLI. Its default of `0.0`
preserves the previous behavior.
Use the optional maximum outside-cluster mean-expression threshold to exclude
genes with excessive average expression in the reference population. Enable
**Limit outside-cluster mean expression** in the web interface or pass
`--max-other-mean` on the CLI. The value uses the input matrix's existing
scale, and omitting it preserves the previous behavior.
Use the optional maximum outside-cluster expressing-cell count to cap the
absolute number of reference cells in which a marker is detected. Enable
**Limit outside-cluster cells expressing gene** in the web interface or pass
`--max-cells-expressing-out` on the CLI. Omitting it preserves the previous
behavior.
Use the minimum mean-expression difference to require a marker's target-cluster
mean to exceed its reference mean by a chosen absolute amount on the input
matrix's existing scale. The `mean_difference` output column records this
target-minus-reference value. Set **Minimum mean-expression difference** in the
web interface or use `--min-mean-difference` on the CLI. Its default of `0.0`
preserves the previous behavior.

The default maximum adjusted p-value is `1.0`, which preserves all otherwise eligible markers. Set a lower value, such as `0.05`, to retain only markers that pass the selected false-discovery-rate threshold.

## Installation

```bash
git clone https://github.com/jec063/MarkerGeneFinderMini.git
cd MarkerGeneFinderMini

python3 -m venv .venv
source .venv/bin/activate

python -m pip install -r requirements.txt
```

To install the analysis library and command-line entry point directly from the
repository, use:

```bash
python -m pip install .
marker-gene-finder --version
```

Install the optional Streamlit interface with `python -m pip install '.[web]'`
or the test dependencies with `python -m pip install '.[test]'`.

## Version 2 API and CLI

Version 2 provides a stable `marker_gene_finder` Python package, a validated
`AnalysisConfig`, a typed `MarkerResult`, documented public exceptions, and an
explicit subcommand CLI. The established analysis engines remain available
behind this interface, so the release changes orchestration without changing
the marker calculations.

Version 2.1 moves the analysis, inspection, batch, and benchmark services into
the `marker_gene_finder` package. The historical `src.marker_finder` module is
now a compatibility entry point, so existing imports and its flat CLI continue
to work during the v2 migration window without remaining a dependency of the
public package.

Version 2.2 moves profiles, presets, expression summaries, settings import,
and report bundles into the public package. New exports use analysis-settings
schema 2, while schema-1 JSON and ZIP bundles are migrated automatically on
import. The Streamlit application now imports only public
`marker_gene_finder` services.

Version 2.3 adds a formal `AnalysisEngine` protocol, normalized `EngineOutput`,
and a validated `EngineRegistry`. The public API dispatches Native, Scanpy, and
automatic file analysis through this interface, and callers can inject a
registry for controlled integrations and tests. Legacy `src.*` imports now
emit deprecation notices and are scheduled for removal in version 3.0.

Version 2.3.1 completes the custom-engine contract: portable configurations
may name registered third-party engines, unknown engines remain configuration
errors, and public API calls validate custom `EngineOutput` objects before
building a result. Engine identifiers use lowercase letters, digits,
underscores, and hyphens.

Version 2.4 adds opt-in engine plugin discovery through Python entry points.
Plugins register under the `marker_gene_finder.engines` group and must expose
an engine whose name matches its entry-point name. Library callers use
`load_engine_plugins()`; the `analyze` command enables discovery with
`--load-engine-plugins`. Discovery returns an independent registry and never
mutates the built-in global registry.

Plugin packages declare their engine class or zero-argument factory in
`pyproject.toml`:

```toml
[project.entry-points."marker_gene_finder.engines"]
my-engine = "my_package.engine:MyEngine"
```

The object must implement `AnalysisEngine` and report `name = "my-engine"`.
Use it with an analysis configuration whose `engine` is `my-engine`, then run
`marker-gene-finder analyze input.csv output.csv --config config.json
--load-engine-plugins`.

Version 2.5 adds plugin diagnostics and contract validation. Run
`marker-gene-finder list-engines` to list built-ins, or add `--load-plugins`
to inspect installed entry points without stopping at the first broken plugin.
The command prints JSON by default and accepts `--output engines.json` for
automation. Plugin authors can use the public `validate_engine()` helper to
check the engine protocol and its advertised name.

```python
import pandas as pd

from marker_gene_finder import AnalysisConfig, analyze

expression = pd.read_csv("expression.csv")
result = analyze(expression, AnalysisConfig(top_n=10, min_log2fc=1.0))
result.to_csv("markers.csv")
print(result.metadata)
```

Run a file analysis with the installed command:

```bash
marker-gene-finder analyze expression.csv markers.csv \
  --config analysis-config.json
```

The versioned configuration schema is shared by the API and every analysis
subcommand. Unknown fields and unsupported schema versions are rejected rather
than silently ignored. A minimal configuration is:

```json
{
  "schema_name": "marker-gene-finder.analysis-config",
  "schema_version": 2,
  "engine": "native",
  "cluster_column": "cluster",
  "cell_column": "cell",
  "top_n": 10
}
```

The v2 commands are `analyze`, `inspect`, `plan`, `batch`, `batch-plan`,
`verify`, `benchmark`, and `compare-benchmarks`. For example:

```bash
marker-gene-finder inspect expression.h5ad --output inspection.json \
  --config analysis-config.json
marker-gene-finder plan expression.h5ad --output plan.json \
  --config analysis-config.json
marker-gene-finder batch inputs outputs --config analysis-config.json --resume
```

### Migrating from v1

The installed command now requires a subcommand and moves analysis settings
into a reusable JSON configuration. For example, replace
`marker-gene-finder --input expression.csv --output markers.csv --top-n 10`
with the `analyze` command above and set `top_n` to `10` in the configuration.
The original flat parser remains temporarily available through
`python src/marker_finder.py` for migration, but new automation should use the
v2 subcommands and public package.

Analysis profiles and bundle settings now use a named envelope:

```json
{
  "schema_name": "marker-gene-finder-analysis-settings",
  "schema_version": 2,
  "app_version": "2.4.0",
  "settings": {"top_n": 10, "min_pct": 0.25}
}
```

Imports still accept the flat schema-1 document and migrate it in memory.
Unsupported future schemas and malformed schema-2 envelopes are rejected.

## Legacy v1 command-line reference

The examples below document the transitional flat parser invoked with
`python src/marker_finder.py`. Prefer the v2 commands above for new workflows.

Inspect a dataset without creating marker results:

```bash
python src/marker_finder.py \
  --input expression.h5ad \
  --cluster-column leiden \
  --inspect-json inspection.json
```

The inspection respects `--layer`, `--use-raw`, and
`--max-dense-memory-mb`. The report records `storage`,
`native_within_limit`, `recommended_engine`, and `recommendation_reasons` so
automation can select an engine before analysis. It cannot be combined with an
analysis output, batch input, or runtime diagnostics.

Create a reviewable execution plan without running marker analysis:

```bash
python src/marker_finder.py \
  --input expression.h5ad \
  --engine auto \
  --cluster-column leiden \
  --plan-json analysis-plan.json
```

After reviewing the plan, run with `--engine auto` to use its resource-based
engine resolution. Runtime diagnostics record both the requested and resolved
engine. For batch input, auto resolution is performed independently per file.

Use `--input-dir` and `--output-dir` to apply one configuration to every CSV
or H5AD file directly inside a directory:

```bash
python src/marker_finder.py \
  --input-dir batch_inputs \
  --output-dir batch_outputs \
  --cluster-column cluster \
  --min-log2fc 1.0 \
  --max-p-adj 0.05
```

Each input receives its own `<input-name>.markers.csv` output. The output
directory also contains `batch_manifest.json`, which records the application
version, output path, marker-row count, and success or error for every input.
Inputs are processed in filename order. A failed input is recorded without
preventing later files from running, and the command exits with a nonzero
status after the batch finishes if any input failed. Batch mode intentionally
uses the same analysis options as single-file mode.

Successful batch records also include cells, genes, clusters, marker rows, and
load, analysis, output, and total durations. The manifest records the overall
batch duration. For a single input, pass `--diagnostics-json diagnostics.json`
to write the same runtime and dataset metadata to a standalone JSON file:

```bash
python src/marker_finder.py \
  --input expression.h5ad \
  --output markers.csv \
  --engine native \
  --cluster-column leiden \
  --gene-chunk-size 500 \
  --workers 2 \
  --diagnostics-json diagnostics.json
```

For automation that should emit only its final summary, add `--quiet`.

Durations use a monotonic high-resolution clock and are intended for runtime
comparison and troubleshooting, not as stable performance assertions. Native
diagnostics also record the resolved gene chunk size.

Batch output is written atomically. Use `--resume` to skip unchanged completed
inputs, `--retry-failed` to retry failures while retaining unchanged successes,
or `--fail-fast` to stop after the first error. The versioned manifest records
input size and modification-time fingerprints, a settings fingerprint, and
`complete`, `skipped`, `failed`, or `pending` status for every input. Changing
an input, analysis setting, or missing output invalidates the saved success and
causes that input to run again.

Benchmark a complete single-input analysis without retaining marker output:

```bash
python src/marker_finder.py \
  --input expression.h5ad \
  --engine auto \
  --benchmark-repeats 3 \
  --benchmark-json benchmark.json
```

The versioned benchmark report records the platform, Python and application
versions, per-repeat runtime diagnostics, and minimum, median, and maximum
total durations. Benchmarks are intended for comparing configurations on the
same machine; shared-runner timings should not be treated as strict performance
guarantees.
They also record the configured worker count.

Compare a candidate benchmark with a baseline and fail when its median runtime
regresses by more than an allowed percentage:

```bash
python src/marker_finder.py \
  --compare-benchmarks baseline.json candidate.json \
  --max-regression-percent 10 \
  --comparison-json comparison.json
```

The comparison records the median ratio and percentage change. It also warns
when the platform, Python version, or input differs, since those differences
can make timings unsuitable for direct comparison. A detected regression is
written to JSON before the command exits with a nonzero status.

Batch manifests now record SHA-256 digests for every input and successful
marker output. Verify that a completed batch still matches those recorded
artifacts without rerunning the analysis:

```bash
python src/marker_finder.py \
  --verify-batch-manifest batch_outputs/batch_manifest.json \
  --verification-json batch-verification.json
```

Verification reports missing or modified inputs and outputs and exits nonzero
after writing its JSON report when any completed record is invalid. Resume mode
also reruns an input when its prior marker output no longer matches its digest.

Select a batch subset with repeatable filename globs. Include patterns are
combined with OR logic, then matching exclude patterns are removed:

```bash
python src/marker_finder.py \
  --input-dir batch_inputs \
  --output-dir batch_outputs \
  --include-pattern 'sample-*.csv' \
  --include-pattern '*.h5ad' \
  --exclude-pattern '*-draft.csv'
```

Preview the exact selected files and whether resume mode would run or skip each
one without executing marker analysis or changing the batch manifest:

```bash
python src/marker_finder.py \
  --input-dir batch_inputs \
  --output-dir batch_outputs \
  --include-pattern 'sample-*.csv' \
  --resume \
  --batch-plan-json batch-plan.json
```

The versioned plan records input fingerprints, the settings fingerprint, run or
skip actions, and explicit reasons such as `new_input`, `input_changed`,
`settings_changed`, `output_missing`, `output_changed`, or `unchanged_success`.

## Web interface

Launch the Streamlit interface:

```bash
python -m streamlit run app.py
```

To capture a local screenshot, install Playwright and its Chromium browser,
start Streamlit, and run the screenshot helper from a second terminal:

```bash
python -m pip install playwright
python -m playwright install chromium
python screenshot_app.py --output marker-gene-finder.png
```

Use `--url` when Streamlit is running somewhere other than
`http://localhost:8501`.

For H5AD uploads, select either **Native** or **Scanpy** from the analysis-engine control. CSV uploads use the native engine.

## Try example data

Select **Use example data** in the Streamlit app to load the bundled
example CSV without uploading a file. Click **Find marker genes** to
explore the results table, marker chart, dot plot, and heatmap.

The example uses the Native engine. Uncheck **Use example data** to
return to uploaded-file analysis. Switching between example data and
uploaded data clears previous marker results.

## Analysis progress

Before running analysis, expand **Cluster sizes** to review each cluster's
cell count and percentage of the dataset. The app warns when a cluster has
fewer than 10 cells because marker statistics may be unstable. Use **Download
cluster summary as CSV** to save this table before or after analysis.

When marker analysis starts, the Streamlit app displays the selected
engine and the number of cells, genes, and clusters being analyzed.
The status panel remains active while marker testing is running.

After completion, it reports the elapsed time and the number of marker
rows retained after filtering. If analysis fails validation, the panel
changes to an error state and displays the error message.

## Marker visualization

After finding markers in Streamlit, use **Cluster to visualize**
to display a horizontal bar chart for one cluster.
The chart shows up to 20 markers with the highest log2 fold changes
from the filtered results and identifies the analysis engine used.

Bar length represents log2 fold change, not statistical significance.
Switching the displayed cluster does not rerun the analysis.
The full results table and CSV download remain available.

## Expression dot plot

After finding markers, select up to 20 genes from the results using
**Genes to compare**. The dot plot compares them across all clusters.

- Dot area represents the fraction of cells with expression greater than zero.
- Color represents mean expression across all cells in the cluster, including zeros.
- No dot is drawn when the fraction expressing is zero.
- Hover over a dot to inspect its expression summary.

The plot supports CSV input and H5AD X, named layers, and raw expression
with both Native and Scanpy analysis. It uses the selected expression
matrix without additional normalization, log transformation, or scaling.
Changing the gene selection does not rerun marker analysis.

## Expression heatmap

The heatmap shares the **Genes to compare** selection with the dot plot
and displays mean expression for each selected gene across all clusters.

- Color represents mean expression across all cells in each cluster, including zeros.
- Hover over a tile to inspect the gene, cluster, mean expression, and fraction expressing.
- Values use the selected matrix's existing scale without additional normalization,
  log transformation, or gene-wise scaling.

The heatmap supports CSV input and H5AD X, named layers, and raw expression
with both Native and Scanpy analysis. Changing the gene selection updates
both expression plots without rerunning marker analysis.

## Download expression summaries

After finding markers, choose genes under **Genes to compare**, then click
**Download expression summary as CSV** below the heatmap.

The CSV contains one row per selected gene and observed cluster:

- `cluster`: cluster label.
- `gene`: selected gene name.
- `mean_expression`: mean across all cells in the cluster, including zeros.
- `fraction_expressing`: fraction of cells with expression greater than zero, from 0 to 1.

The export uses the same expression summary as the dot plot and heatmap,
on the selected matrix's existing scale. Changing the gene selection
updates the export. The download is hidden when no genes are selected.

## Download analysis settings

After a successful analysis, click **Download analysis settings as JSON**
to save the settings associated with the results.

The export records the input filename, example-data status, analysis
engine, expression source, layer or raw selection, cluster and cell-ID
settings, filtering thresholds, dataset dimensions, and marker-row count.
The pseudocount is recorded for Native analysis and is null for Scanpy.

Settings are captured when analysis succeeds, including runs where no
markers pass filtering. Changing an input or analysis setting clears
both the previous results and their settings download.

The JSON records analysis settings; it does not contain the expression
data. Use **Import analysis settings** in the Streamlit sidebar to restore
portable thresholds, ranking, exclusions, and preset selection from either the
standalone JSON file or a previously downloaded analysis-bundle ZIP.
Dataset-specific choices—including the input file, expression source, cluster
column, target clusters, and cell-ID column—are intentionally left unchanged
and must match the current dataset. Review the restored controls before running
the analysis.

## Saved analysis profiles

Use **Download current profile** to save the portable controls at any time after
loading a dataset; running marker analysis is not required. Give the profile a
descriptive name before downloading it. Profiles use the same validated JSON
schema as exported analysis settings and can be imported from a standalone JSON
file or an analysis-bundle ZIP.

Selecting a profile displays a before-and-after table of settings that would
change. The app warns when a profile was created by a newer application version,
and no values change until **Apply imported settings** is clicked. Use **Reset
profile to defaults** to restore every portable control to the stable defaults.
Dataset-specific selections are deliberately excluded from profiles.

## Download an analysis bundle

After a successful analysis, click **Download analysis bundle as ZIP** to save
the marker-results CSV, analysis-settings JSON, and a short README together.
The bundle includes `cluster_summary.csv` with cluster counts and percentages.
When genes are selected for the expression plots, it also contains the current
`expression_summary.csv`. The archive does not contain the uploaded expression
matrix.

## Run

```bash
python src/marker_finder.py \
  --input data/example_expression.csv \
  --output outputs/top_markers.csv \
  --cluster-column cluster \
  --top-n 2 \
  --min-cluster-size 10 \
  --min-cells-expressing 3 \
  --min-cluster-mean 1.0 \
  --min-mean-difference 0.5 \
  --max-other-mean 0.5 \
  --max-cells-expressing-out 5 \
  --min-pct 0.5 \
  --max-pct-out 0.25 \
  --rank-by pct_difference \
  --min-log2fc 0.0 \
  --max-p-adj 0.05
```

For H5AD input, select the `AnnData.obs` column containing
cluster labels:

```bash
python src/marker_finder.py --input path/to/expression.h5ad --output outputs/top_markers.csv --cluster-column leiden --engine scanpy
```
The example should identify:

- Cluster 0: `CD3D`, `CD3E`
- Cluster 1: `MS4A1`, `CD79A`
- Cluster 2: `LST1`, `S100A8`
- Cluster 3: `PECAM1`, `VWF`

## H5AD expression layers

Both Native and Scanpy use `AnnData.X` by default. To analyze a named
matrix in `AnnData.layers`, pass `--layer`:

```bash
python src/marker_finder.py --input expression.h5ad --output outputs/markers.csv --cluster-column leiden --engine scanpy --layer lognormalized
```

Replace `lognormalized` with a layer name present in your file.
The same option works with `--engine native`. CSV input does not support
layer selection, and missing layer names produce an error.

In Streamlit, use the **Expression matrix** dropdown after uploading an
H5AD file. Changing the selected expression matrix clears previous results.

Layer selection does not normalize or log-transform expression. For
Scanpy, select appropriately normalized, log-transformed data. Layer
names alone do not establish how the data were processed.


## AnnData.raw expression

Both engines support expression stored in `AnnData.raw`:

```bash
python src/marker_finder.py --input expression.h5ad --output outputs/raw_markers.csv --cluster-column leiden --engine scanpy --use-raw
```

The same option works with `--engine native`. Raw selection uses
`AnnData.raw.X` and its own gene names, which may differ from X.
Cell identifiers and cluster labels use the current observations.
The input file is not modified.

`--use-raw` and `--layer` are mutually exclusive. Raw selection
requires H5AD input and reports an error if raw expression is missing.
Omitting both options continues to use `AnnData.X`.

In Streamlit, choose **AnnData.raw** from **Expression matrix**.
This option appears only when raw expression is available.
Changing the selection clears previous results. For Scanpy,
the displayed gene count reflects the selected matrix.

Raw selection does not normalize or log-transform expression.
The name `raw` does not establish how the data were processed.
For Scanpy, use appropriately normalized, log-transformed expression.

## Test

```bash
pytest -q
```

## Version 2.6 engine workflows

Plugins can now be used with `analyze`, `plan`, `batch-plan`, `batch`, and
`benchmark`: pass `--load-engine-plugins` and a `--config` file naming the
engine. Python callers can pass `registry=...` to `analyze`, `analyze_file`,
`build_analysis_plan`, `build_batch_plan`, `run_batch`, and `benchmark_analysis`.
Registries remain local to the call; discovery never modifies built-in defaults.

Engines may declare an optional `capabilities = EngineCapabilities(...)`:

| Capability | Default | Meaning |
| --- | --- | --- |
| `input_formats` | `("csv", "h5ad", "dataframe")` | Accepted input kinds |
| `supports_layers` | `True` | Named H5AD layers |
| `supports_raw` | `True` | H5AD raw expression |
| `supports_workers` | `True` | Worker counts greater than one |

Import `EngineCapabilities` from `marker_gene_finder`. Scanpy declares H5AD-only
input and one worker. Unsupported declared combinations are rejected before
engine execution. Plans report capability conflicts through `can_execute` and
`warnings`; a CLI batch plan exits unsuccessfully if any input cannot execute.
Resource recommendations still describe Native versus Scanpy, not plugin
performance. Plugins without declarations retain the permissive v2 contract;
they must validate any additional engine-specific restrictions themselves.

Plugin marker tables must have unique columns including non-null `cluster` and
`gene`. Additional statistical columns remain optional for v2 compatibility,
but known numeric columns must be finite, probabilities and expression fractions
must be in [0, 1], and counts/ranks must be valid integers. Diagnostics must be
JSON-compatible with finite numbers. Batch execution validates results before
promoting temporary output to a successful artifact, and benchmarks measure
elapsed time independently of plugin diagnostics. Built-in statistical defaults
are unchanged.
