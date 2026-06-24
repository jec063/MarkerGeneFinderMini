# Marker Gene Finder Mini

A small Python project for identifying marker genes from clustered gene expression data.

This project compares average gene expression inside one cluster versus outside that cluster, ranks genes by log2 fold change, and outputs the top marker genes for each cluster.

## Goal

The goal is to practice a basic single-cell analysis workflow:

```text
expression matrix + cluster labels
↓
marker gene detection
↓
top markers per cluster
```

## Input Format

The input is a CSV file where:

* each row is one cell
* one column contains the cluster label
* the remaining columns are gene expression values

Example:

```csv
cell,cluster,CD3D,CD3E,MS4A1,CD79A,LST1,S100A8,PECAM1,VWF
cell1,0,5,4.8,0,0,0.2,0.1,0,0
cell2,0,4.7,5,0.1,0,0.1,0.2,0,0
cell3,1,0.1,0,5,4.9,0.2,0.1,0,0
cell4,1,0,0.1,5.2,5,0.1,0.1,0,0
```

## Method

For each cluster and each gene, the script calculates:

1. Mean expression inside the cluster
2. Mean expression outside the cluster
3. Log2 fold change

Formula:

```text
log2FC = log2((cluster_mean + 0.1) / (other_mean + 0.1))
```

The `0.1` is a pseudocount that prevents division by zero.

Genes with the highest log2 fold change are selected as marker genes.

## Project Structure

```text
MarkerGeneFinderMini/
  README.md
  data/
    example_expression.csv
  src/
    marker_finder.py
  outputs/
    top_markers.csv
  requirements.txt
  .gitignore
```

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/MarkerGeneFinderMini.git
cd MarkerGeneFinderMini
```

Install the required packages:

```bash
pip install -r requirements.txt
```

## Run

```bash
python src/marker_finder.py
```

## Output

The script saves the result to:

```text
outputs/top_markers.csv
```

Example output:

| cluster | gene  | cluster_mean | other_mean | log2FC |
| ------- | ----- | -----------: | ---------: | -----: |
| 0       | CD3D  |         4.85 |       0.07 |   5.39 |
| 0       | CD3E  |         4.90 |       0.05 |   5.64 |
| 1       | MS4A1 |         5.10 |       0.03 |   5.98 |

## Expected Result

With the example dataset, the script should identify:

* Cluster 0: T-cell markers such as `CD3D` and `CD3E`
* Cluster 1: B-cell markers such as `MS4A1` and `CD79A`
* Cluster 2: Monocyte/macrophage markers such as `LST1` and `S100A8`
* Cluster 3: Endothelial markers such as `PECAM1` and `VWF`

## Future Improvements

Possible next steps:

* Add filtering by minimum expression
* Add statistical testing
* Add adjusted p-values
* Add a simple web interface
* Add cell type prediction using known marker genes
* Add plain-English explanations for each cluster

## Why This Project Matters

Marker gene detection is one of the first steps in understanding clustered single-cell expression data. This mini project helps practice the basic logic behind comparing clusters and identifying genes that define each group of cells.
