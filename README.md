# Marker Gene Finder Mini

A small Python project for identifying candidate marker genes from clustered gene-expression data. This is the first working component of the broader CellCoPilot project.

**Live app:** [Launch the CellCoPilot Marker Gene Finder](https://cellcopilot-marker-finder.streamlit.app/)

The program compares gene expression inside each cluster with all other cells and ranks genes using log2 fold change.

## Input

The input CSV contains:

- one row per cell
- a `cell` identifier column
- a `cluster` label column
- numeric gene-expression columns

See `data/example_expression.csv` for an example.

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

The default pseudocount is `0.1`.

By default, genes with negative log2 fold change are excluded. A higher minimum log2 fold-change threshold can be used to require stronger enrichment.

Statistical testing uses a one-sided Mann-Whitney U test, which is the Wilcoxon rank-sum test for independent samples. The alternative hypothesis is that expression is higher inside the cluster than outside it. For each cluster, p-values for all genes are adjusted together using the Benjamini-Hochberg false-discovery-rate procedure.

Results are ranked by log2 fold change after applying the expression and fold-change thresholds. The output includes raw p-values in `p_value` and adjusted p-values in `p_adj`.

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

## Run

```bash
python src/marker_finder.py \
  --input data/example_expression.csv \
  --output outputs/top_markers.csv \
  --cluster-column cluster \
  --top-n 2 \
  --min-pct 0.5 \
  --min-log2fc 0.0
```

The example should identify:

- Cluster 0: `CD3D`, `CD3E`
- Cluster 1: `MS4A1`, `CD79A`
- Cluster 2: `LST1`, `S100A8`
- Cluster 3: `PECAM1`, `VWF`

## Test

```bash
pytest -q
```

## Next milestone

Planned improvements:

- adjusted p-value filtering
- support for `.h5ad` files
- Scanpy integration
