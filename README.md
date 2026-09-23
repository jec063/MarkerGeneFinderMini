# Marker Gene Finder Mini

A small Python project for identifying candidate marker genes from clustered gene-expression data. This is the first working component of the broader CellCoPilot project.

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

## Analysis engines

MarkerGeneFinderMini provides two analysis engines:

- **Native:** the original CellCoPilot implementation for CSV and H5AD input.
- **Scanpy:** an H5AD-only engine using `scanpy.tl.rank_genes_groups()` with the Wilcoxon method and Benjamini-Hochberg correction.

The native engine remains the default. Scanpy works directly with AnnData and reports its own approximate log2 fold changes.

Select one or more **Clusters to analyze** in the web interface to limit the
reported groups. Each selected cluster is still compared with every other cell.
On the CLI, repeat `--target-cluster`, for example
`--target-cluster 0 --target-cluster 2`. Omitting it analyzes every cluster.

Exact genes and gene-name prefixes can be excluded before testing. In the web
interface, enter comma-separated values under **Exclude genes** or **Exclude
gene prefixes**. On the CLI, repeat `--exclude-gene` or `--exclude-prefix`.
No exclusions are applied by default; organism-specific naming conventions
remain under the user's control.

## Method

For each cluster and gene, the program calculates:

1. Mean expression inside the cluster
2. Mean expression outside the cluster
3. Fraction of cells expressing the gene inside and outside the cluster
4. Difference between the inside and outside expression fractions
5. Log2 fold change
6. One-sided Wilcoxon-Mann-Whitney p-value
7. Benjamini-Hochberg adjusted p-value

```text
log2FC = log2(
    (cluster_mean + pseudocount)
    /
    (other_mean + pseudocount)
)
```

The native engine default pseudocount is `0.1`.

By default, genes with negative log2 fold change are excluded. A higher minimum log2 fold-change threshold can be used to require stronger enrichment.

Statistical testing uses a one-sided Mann-Whitney U test, which is the Wilcoxon rank-sum test for independent samples. The alternative hypothesis is that expression is higher inside the cluster than outside it. For each cluster, p-values for all genes are adjusted together using the Benjamini-Hochberg false-discovery-rate procedure.

Results are ranked by log2 fold change after applying the selected thresholds. The output includes raw p-values in `p_value` and adjusted p-values in `p_adj`.
The one-based `rank` column records each gene's position within its cluster
after filtering, so ranks restart at 1 for every cluster.
Use the minimum expression-fraction difference to require markers to be
expressed in a larger share of target-cluster cells than other cells.

The default maximum adjusted p-value is `1.0`, which preserves all otherwise eligible markers. Set a lower value, such as `0.05`, to retain only markers that pass the selected false-discovery-rate threshold.

## Installation

```bash
git clone https://github.com/jec063/MarkerGeneFinderMini.git
cd MarkerGeneFinderMini

python3 -m venv .venv
source .venv/bin/activate

python -m pip install -r requirements.txt
```

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
fewer than 10 cells because marker statistics may be unstable.

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
data or automatically restore a session.

## Download an analysis bundle

After a successful analysis, click **Download analysis bundle as ZIP** to save
the marker-results CSV, analysis-settings JSON, and a short README together.
When genes are selected for the expression plots, the bundle also contains the
current `expression_summary.csv`. The archive does not contain the uploaded
expression matrix.

## Run

```bash
python src/marker_finder.py \
  --input data/example_expression.csv \
  --output outputs/top_markers.csv \
  --cluster-column cluster \
  --top-n 2 \
  --min-pct 0.5 \
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
