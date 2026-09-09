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

## Method

For each cluster and gene, the program calculates:

1. Mean expression inside the cluster
2. Mean expression outside the cluster
3. Fraction of cells expressing the gene inside and outside the cluster
4. Log2 fold change
5. One-sided Wilcoxon-Mann-Whitney p-value
6. Benjamini-Hochberg adjusted p-value

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

For H5AD uploads, select either **Native** or **Scanpy** from the analysis-engine control. CSV uploads use the native engine.

## Marker visualization

After finding markers in Streamlit, use **Cluster to visualize**
to display a horizontal bar chart for one cluster.
The chart shows up to 20 markers with the highest log2 fold changes
from the filtered results and identifies the analysis engine used.

Bar length represents log2 fold change, not statistical significance.
Switching the displayed cluster does not rerun the analysis.
The full results table and CSV download remain available.

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

## Next milestone

Planned improvements:

- Add expression dot plots and heatmaps
- Improve progress reporting for larger datasets
