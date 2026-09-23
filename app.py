"""Streamlit interface for the CellCoPilot Marker Gene Finder."""

from __future__ import annotations

import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd
import streamlit as st

from src.expression_summary import summarize_clusters, summarize_expression

from src.marker_finder import (
    find_markers,
    find_markers_scanpy,
    load_expression,
)
from src.report_bundle import build_analysis_bundle
from time import perf_counter


st.set_page_config(
    page_title="CellCoPilot Marker Finder",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 CellCoPilot Marker Gene Finder")
st.write(
    "Upload clustered gene-expression data and identify candidate "
    "marker genes for each cluster."
)


@st.cache_data(show_spinner=False)
def load_csv(file_contents: bytes) -> pd.DataFrame:
    """Load an uploaded CSV file."""
    return pd.read_csv(io.BytesIO(file_contents))


@st.cache_data(show_spinner=False)
def h5ad_obs_columns(
    file_contents: bytes,
    file_name: str,
) -> list[str]:
    """Return available observation columns from an h5ad file."""
    with TemporaryDirectory() as directory:
        input_path = Path(directory) / Path(file_name).name
        input_path.write_bytes(file_contents)

        adata = ad.read_h5ad(input_path, backed="r")
        try:
            return [
                str(column)
                for column in adata.obs.columns
            ]
        finally:
            adata.file.close()


@st.cache_data(show_spinner=False)
def load_uploaded_h5ad(
    file_contents: bytes,
    file_name: str,
    cluster_column: str,
    layer: str | None = None,
    use_raw: bool = False,
) -> pd.DataFrame:
    """Convert an uploaded h5ad file to an expression table."""
    with TemporaryDirectory() as directory:
        input_path = Path(directory) / Path(file_name).name
        input_path.write_bytes(file_contents)

        return load_expression(
            input_path,
            cluster_column=cluster_column,
            cell_column="cell_id",
            layer=layer,
            use_raw=use_raw,
        )


@st.cache_data(show_spinner=False)
def load_uploaded_anndata(
    file_contents: bytes,
    file_name: str,
) -> ad.AnnData:
    """Load an uploaded H5AD file without converting it to a table."""
    with TemporaryDirectory() as directory:
        input_path = Path(directory) / Path(file_name).name
        input_path.write_bytes(file_contents)

        return ad.read_h5ad(input_path)


@st.cache_data(show_spinner=False)
def h5ad_expression_sources(
    file_contents: bytes,
    file_name: str,
) -> list[tuple[str, str | None]]:
    """List available expression matrices without densifying them."""
    with TemporaryDirectory() as directory:
        input_path = Path(directory) / Path(file_name).name
        input_path.write_bytes(file_contents)
        adata = ad.read_h5ad(input_path, backed="r")
        try:
            sources = [("X", None)]
            if adata.raw is not None and adata.raw.X is not None:
                sources.append(("raw", None))
            sources.extend(("layer", name) for name in adata.layers.keys())
            return sources
        finally:
            adata.file.close()


def clear_marker_results() -> None:
    """Clear results when the input or analysis settings change."""
    st.session_state.pop("marker_results", None)
    st.session_state.pop("analysis_settings", None)


use_example_data = st.checkbox(
    "Use example data",
    key="use_example_data",
    on_change=clear_marker_results,
    help="Try the app with the bundled example CSV.",
)

uploaded_file = st.file_uploader(
    "Upload an expression CSV or H5AD file",
    on_change=clear_marker_results,
    type=["csv", "h5ad"],
    help=(
        "CSV files should contain one row per cell, a cluster "
        "column, an optional cell-ID column, and numeric gene "
        "columns. H5AD expression can come from X, layers, or raw, "
        "and cluster labels in AnnData.obs."
    ),
)

if use_example_data:
    example_path = Path(__file__).resolve().parent / "data" / "example_expression.csv"
    try:
        file_contents = example_path.read_bytes()
    except OSError as error:
        st.error(f"Could not load example data: {error}")
        st.stop()
    file_name = example_path.name
    st.info(
        "Using the bundled example CSV. Uncheck Use example data "
        "to analyze an uploaded file."
    )
else:
    if uploaded_file is None:
        st.info("Upload a CSV or H5AD file to begin.")
        st.stop()
    file_contents = uploaded_file.getvalue()
    file_name = uploaded_file.name

file_suffix = Path(file_name).suffix.lower()

if file_suffix == ".csv":
    try:
        expression = load_csv(file_contents)
    except (
        pd.errors.ParserError,
        pd.errors.EmptyDataError,
        UnicodeDecodeError,
    ) as error:
        st.error(f"Could not read the uploaded CSV: {error}")
        st.stop()

    if expression.empty:
        st.error("The uploaded CSV is empty.")
        st.stop()

    columns = expression.columns.tolist()
    cluster_options = columns
    cluster_index = (
        columns.index("cluster")
        if "cluster" in columns
        else 0
    )
else:
    try:
        cluster_options = h5ad_obs_columns(
            file_contents,
            file_name,
        )
    except (KeyError, OSError, ValueError) as error:
        st.error(f"Could not inspect the uploaded H5AD file: {error}")
        st.stop()

    if not cluster_options:
        st.error(
            "The uploaded H5AD file has no columns in AnnData.obs."
        )
        st.stop()

    preferred_cluster_columns = [
        "cluster",
        "cell_type",
        "leiden",
        "customclassif",
    ]
    default_cluster = next(
        (
            column
            for column in preferred_cluster_columns
            if column in cluster_options
        ),
        cluster_options[0],
    )
    cluster_index = cluster_options.index(default_cluster)

with st.sidebar:
    st.header("Analysis settings")

    cluster_column = st.selectbox(
        "Cluster column",
        on_change=clear_marker_results,
        options=cluster_options,
        index=cluster_index,
    )

    if file_suffix == ".h5ad":
        engine_label = st.selectbox(
            "Analysis engine",
            on_change=clear_marker_results,
            options=["Native", "Scanpy"],
            help=(
                "Native preserves the existing CellCoPilot method. "
                "Scanpy uses scanpy.tl.rank_genes_groups with the "
                "Wilcoxon method."
            ),
        )
        analysis_engine = engine_label.lower()

        try:
            expression_sources = h5ad_expression_sources(
                file_contents,
                file_name,
            )
        except (KeyError, OSError, ValueError) as error:
            st.error(f"Could not inspect H5AD expression matrices: {error}")
            st.stop()

        selected_source = st.selectbox(
            "Expression matrix",
            options=expression_sources,
            format_func=lambda source: (
                "AnnData.X (default)" if source[0] == "X"
                else "AnnData.raw" if source[0] == "raw"
                else f"Layer: {source[1]}"
            ),
            on_change=clear_marker_results,
            help=(
                "Choose X, a named layer, or raw when available. "
                "This does not normalize or log-transform the data. "
                "For Scanpy, select normalized, log-transformed expression. "
                "The name raw does not indicate how the data were processed."
            ),
        )
        selected_layer = (
            selected_source[1] if selected_source[0] == "layer" else None
        )
        use_raw = selected_source[0] == "raw"
    else:
        analysis_engine = "native"
        selected_layer = None
        use_raw = False

    if file_suffix == ".csv":
        cell_options = [
            "None",
            *[
                column
                for column in columns
                if column != cluster_column
            ],
        ]
        cell_index = (
            cell_options.index("cell")
            if "cell" in cell_options
            else 0
        )
        selected_cell_column = st.selectbox(
            "Cell-ID column",
            on_change=clear_marker_results,
            options=cell_options,
            index=cell_index,
        )
        cell_column = (
            None
            if selected_cell_column == "None"
            else selected_cell_column
        )
    else:
        st.caption(
            "Cell IDs are read from AnnData.obs_names."
        )
        cell_column = "cell_id"

adata = None

if file_suffix == ".h5ad":
    try:
        if analysis_engine == "scanpy":
            adata = load_uploaded_anndata(
                file_contents,
                file_name,
            )
            expression = None
        else:
            expression = load_uploaded_h5ad(
                file_contents,
                file_name,
                cluster_column,
                layer=selected_layer,
                use_raw=use_raw,
            )
    except (KeyError, OSError, ValueError) as error:
        st.error(f"Could not read the uploaded H5AD file: {error}")
        st.stop()

if analysis_engine == "scanpy":
    if adata is None or adata.n_obs == 0:
        st.error("The uploaded expression dataset is empty.")
        st.stop()
else:
    if expression is None or expression.empty:
        st.error("The uploaded expression dataset is empty.")
        st.stop()

st.subheader("Data preview")

if analysis_engine == "scanpy":
    preview = adata.obs.head(20).copy()
    preview.insert(0, "cell_id", preview.index.astype(str))
    st.caption(
        "Showing AnnData observations without converting the "
        "expression matrix to a dense table."
    )
    st.dataframe(preview, hide_index=True, width="stretch")
else:
    st.dataframe(expression.head(20), width="stretch")

with st.sidebar:

    top_n = st.number_input(
        "Markers per cluster",
        on_change=clear_marker_results,
        min_value=1,
        max_value=100,
        value=5,
        step=1,
    )

    observed_clusters = sorted(
        (
            adata.obs[cluster_column].unique()
            if analysis_engine == "scanpy"
            else expression[cluster_column].unique()
        ),
        key=str,
    )
    target_clusters = st.multiselect(
        "Clusters to analyze",
        options=observed_clusters,
        default=observed_clusters,
        format_func=str,
        on_change=clear_marker_results,
        help="Each selected cluster is compared with all other cells.",
    )
    excluded_gene_text = st.text_input(
        "Exclude genes",
        on_change=clear_marker_results,
        help="Comma-separated exact gene names.",
    )
    excluded_prefix_text = st.text_input(
        "Exclude gene prefixes",
        on_change=clear_marker_results,
        help="Comma-separated prefixes, such as MT-, RPL, RPS.",
    )
    excluded_genes = [
        value.strip() for value in excluded_gene_text.split(",") if value.strip()
    ]
    excluded_prefixes = [
        value.strip() for value in excluded_prefix_text.split(",") if value.strip()
    ]

    min_pct = st.slider(
        "Minimum fraction expressing gene",
        on_change=clear_marker_results,
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.05,
    )

    min_log2fc = st.number_input(
        "Minimum log2 fold change",
        on_change=clear_marker_results,
        min_value=0.0,
        value=0.0,
        step=0.1,
        format="%.2f",
    )

    min_pct_difference = st.slider(
        "Minimum expression-fraction difference",
        on_change=clear_marker_results,
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.05,
        help="Minimum pct_in minus pct_out required for a marker.",
    )

    max_p_adj = st.number_input(
        "Maximum adjusted p-value",
        on_change=clear_marker_results,
        min_value=0.0,
        max_value=1.0,
        value=1.0,
        step=0.01,
        format="%.3f",
        help=(
            "Use 0.05 to keep markers with "
            "Benjamini-Hochberg adjusted p-values "
            "at or below 0.05."
        ),
    )

    if analysis_engine == "native":
        pseudocount = st.number_input(
            "Pseudocount",
            on_change=clear_marker_results,
            min_value=0.0001,
            value=0.1,
            step=0.1,
            format="%.4f",
        )
    else:
        pseudocount = 0.1
        st.caption(
            "Scanpy calculates its own approximate log2 fold changes."
        )

metric_columns = st.columns(3)

if analysis_engine == "scanpy":
    metric_columns[0].metric("Cells", adata.n_obs)
    metric_columns[1].metric(
        "Genes", adata.raw.n_vars if use_raw else adata.n_vars
    )
    metric_columns[2].metric(
        "Clusters",
        adata.obs[cluster_column].nunique(),
    )
else:
    metric_columns[0].metric("Cells", len(expression))
    metric_columns[1].metric("Columns", len(expression.columns))
    metric_columns[2].metric(
        "Clusters",
        expression[cluster_column].nunique(),
    )

cluster_labels = (
    adata.obs[cluster_column]
    if analysis_engine == "scanpy"
    else expression[cluster_column]
)
cluster_summary = summarize_clusters(cluster_labels)
with st.expander("Cluster sizes"):
    st.dataframe(cluster_summary, hide_index=True, width="stretch")
    small_clusters = cluster_summary.loc[cluster_summary["cells"] < 10, "cluster"]
    if not small_clusters.empty:
        st.warning(
            "Clusters with fewer than 10 cells may produce unstable marker "
            "statistics: " + ", ".join(map(str, small_clusters))
        )

if st.button("Find marker genes", type="primary"):
    clear_marker_results()
    engine_display = analysis_engine.title()

    if analysis_engine == "scanpy":
        cell_count = int(adata.n_obs)
        gene_count = int(adata.raw.n_vars if use_raw else adata.n_vars)
        cluster_count = int(adata.obs[cluster_column].nunique())
    else:
        cell_count = int(len(expression))
        gene_count = int(
            len(
                [
                    column
                    for column in expression.columns
                    if column not in {cluster_column, cell_column}
                ]
            )
        )
        cluster_count = int(expression[cluster_column].nunique())

    analysis_status = st.status(
        f"Preparing {engine_display} analysis...",
        expanded=True,
    )
    analysis_status.write(
        f"Input: {cell_count:,} cells, {gene_count:,} genes, "
        f"{cluster_count:,} clusters."
    )
    analysis_status.write(
        "Testing genes and ranking markers. "
        "Large datasets may take several minutes."
    )
    analysis_status.update(
        label=f"Running {engine_display} marker analysis...",
        state="running",
    )

    started_at = perf_counter()

    try:
        if not target_clusters:
            raise ValueError("Select at least one cluster to analyze.")
        if analysis_engine == "scanpy":
            markers = find_markers_scanpy(
                adata,
                layer=selected_layer,
                use_raw=use_raw,
                cluster_column=cluster_column,
                top_n=int(top_n),
                min_pct=float(min_pct),
                min_pct_difference=float(min_pct_difference),
                min_log2fc=float(min_log2fc),
                max_p_adj=float(max_p_adj),
                target_clusters=target_clusters,
                excluded_genes=excluded_genes,
                excluded_prefixes=excluded_prefixes,
            )
        else:
            markers = find_markers(
                expression,
                cluster_column=cluster_column,
                cell_column=cell_column,
                top_n=int(top_n),
                pseudocount=float(pseudocount),
                min_pct=float(min_pct),
                min_pct_difference=float(min_pct_difference),
                min_log2fc=float(min_log2fc),
                max_p_adj=float(max_p_adj),
                target_clusters=target_clusters,
                excluded_genes=excluded_genes,
                excluded_prefixes=excluded_prefixes,
            )
    except ValueError as error:
        elapsed = perf_counter() - started_at
        analysis_status.update(
            label=f"Analysis failed after {elapsed:.1f} seconds.",
            state="error",
            expanded=True,
        )
        st.error(str(error))
    else:
        elapsed = perf_counter() - started_at
        st.session_state["marker_results"] = markers
        st.session_state["analysis_settings"] = {
            "schema_version": 1,
            "input_file": Path(file_name).name,
            "example_data": bool(use_example_data),
            "engine": analysis_engine,
            "expression_source": (
                "CSV" if file_suffix == ".csv"
                else "AnnData.raw.X" if use_raw
                else "AnnData.layers" if selected_layer is not None
                else "AnnData.X"
            ),
            "layer": selected_layer,
            "use_raw": bool(use_raw),
            "cluster_column": cluster_column,
            "cell_column": cell_column if file_suffix == ".csv" else None,
            "cell_id_source": (
                cell_column if file_suffix == ".csv" else "AnnData.obs_names"
            ),
            "top_n": int(top_n),
            "min_pct": float(min_pct),
            "min_pct_difference": float(min_pct_difference),
            "min_log2fc": float(min_log2fc),
            "max_p_adj": float(max_p_adj),
            "pseudocount": (
                float(pseudocount) if analysis_engine == "native" else None
            ),
            "cells": cell_count,
            "genes": gene_count,
            "clusters": cluster_count,
            "target_clusters": [str(cluster) for cluster in target_clusters],
            "excluded_genes": excluded_genes,
            "excluded_prefixes": excluded_prefixes,
            "marker_rows": int(len(markers)),
        }
        analysis_status.write(
            f"Produced {len(markers):,} marker rows after filtering."
        )
        analysis_status.update(
            label=f"Analysis complete in {elapsed:.1f} seconds.",
            state="complete",
            expanded=False,
        )

if "marker_results" in st.session_state:
    markers = st.session_state["marker_results"]
    expression_summary_for_bundle = None

    st.subheader("Marker-gene results")
    if "analysis_settings" in st.session_state:
        settings_json = (
            json.dumps(
                st.session_state["analysis_settings"],
                indent=2,
                ensure_ascii=False,
                allow_nan=False,
            ) + "\n"
        ).encode("utf-8")
        st.download_button(
            "Download analysis settings as JSON",
            data=settings_json,
            file_name="analysis_settings.json",
            mime="application/json",
            key="analysis_settings_download",
        )

    if markers.empty:
        st.warning(
            "No genes passed the selected filtering settings."
        )
    else:
        st.dataframe(
            markers,
            hide_index=True,
            width="stretch",
        )

        st.subheader("Marker visualization")
        chart_cluster = st.selectbox(
            "Cluster to visualize",
            options=sorted(markers["cluster"].unique().tolist(), key=str),
            format_func=str,
        )
        chart_data = (
            markers.loc[
                markers["cluster"] == chart_cluster,
                ["gene", "log2FC"],
            ]
            .sort_values(
                ["log2FC", "gene"],
                ascending=[False, True],
            )
            .head(20)
            .copy()
        )
        engine_name = (
            "Scanpy" if analysis_engine == "scanpy" else "Native"
        )
        st.caption(
            f"{engine_name} engine · Cluster {chart_cluster} · "
            f"{len(chart_data)} markers shown. "
            "Shows up to 20 markers from the filtered results. "
            "Bar length represents log2 fold change, not statistical significance."
        )
        st.bar_chart(
            chart_data,
            x="gene",
            y="log2FC",
            horizontal=True,
            x_label="Gene",
            y_label="log2 fold change",
            color="#2b6ca3",
            height=max(300, 28 * len(chart_data) + 80),
        )

        st.subheader("Expression dot plot")
        available_genes = list(dict.fromkeys(markers["gene"].tolist()))
        dot_genes = st.multiselect(
            "Genes to compare",
            options=available_genes,
            default=available_genes[:10],
            max_selections=20,
            help="Choose up to 20 genes from the marker results.",
        )
        if not dot_genes:
            st.info("Select at least one gene to display the expression plots.")
        else:
            try:
                if analysis_engine == "scanpy":
                    dot_data = summarize_expression(
                        adata,
                        dot_genes,
                        cluster_column=cluster_column,
                        layer=selected_layer,
                        use_raw=use_raw,
                    )
                else:
                    dot_data = summarize_expression(
                        expression,
                        dot_genes,
                        cluster_column=cluster_column,
                    )
            except (ValueError, KeyError) as error:
                st.error(f"Could not create expression plots: {error}")
            else:
                dot_data["cluster"] = dot_data["cluster"].astype(str)
                cluster_order = dot_data["cluster"].drop_duplicates().tolist()
                source_name = (
                    "CSV expression" if file_suffix == ".csv"
                    else "AnnData.raw" if use_raw
                    else f"Layer: {selected_layer}" if selected_layer is not None
                    else "AnnData.X"
                )
                st.caption(
                    f"{engine_name} engine · {source_name}. "
                    "Dot area shows the percentage of cells with expression > 0. "
                    "Color shows mean expression across all cells in each cluster, "
                    "including zeros, on the selected matrix's existing scale. "
                    "A missing dot means 0% expressing."
                )
                st.vega_lite_chart(
                    dot_data,
                    spec={
                        "mark": {"type": "circle", "opacity": 1},
                        "height": max(300, 30 * len(cluster_order)),
                        "encoding": {
                            "x": {
                                "field": "gene",
                                "type": "nominal",
                                "sort": dot_genes,
                                "title": "Gene",
                                "axis": {"labelAngle": -45},
                            },
                            "y": {
                                "field": "cluster",
                                "type": "nominal",
                                "sort": cluster_order,
                                "title": "Cluster",
                            },
                            "size": {
                                "field": "fraction_expressing",
                                "type": "quantitative",
                                "title": "Cells expressing",
                                "scale": {
                                    "type": "linear",
                                    "domain": [0, 1],
                                    "range": [0, 450],
                                },
                                "legend": {
                                    "format": ".0%",
                                    "values": [0.25, 0.5, 0.75, 1],
                                    "orient": "bottom",
                                    "direction": "horizontal",
                                    "columns": 4,
                                },
                            },
                            "color": {
                                "field": "mean_expression",
                                "type": "quantitative",
                                "title": "Mean expression",
                                "scale": {"scheme": "blues", "zero": True},
                            },
                            "tooltip": [
                                {"field": "cluster", "type": "nominal"},
                                {"field": "gene", "type": "nominal"},
                                {
                                    "field": "mean_expression",
                                    "type": "quantitative",
                                    "format": ".3f",
                                },
                                {
                                    "field": "fraction_expressing",
                                    "type": "quantitative",
                                    "format": ".1%",
                                },
                            ],
                        },
                    },
                    use_container_width=True,
                )

                st.subheader("Expression heatmap")
                st.caption(
                    f"{engine_name} engine · {source_name}. "
                    "Color shows mean expression across all cells in each cluster, "
                    "including zeros. Values use the selected matrix's existing "
                    "scale without additional transformation or gene-wise scaling."
                )
                st.vega_lite_chart(
                    dot_data,
                    spec={
                        "mark": {
                            "type": "rect",
                            "stroke": "white",
                            "strokeWidth": 1,
                        },
                        "height": max(180, 30 * len(cluster_order)),
                        "encoding": {
                            "x": {
                                "field": "gene",
                                "type": "nominal",
                                "sort": dot_genes,
                                "title": "Gene",
                                "axis": {"labelAngle": -45},
                            },
                            "y": {
                                "field": "cluster",
                                "type": "nominal",
                                "sort": cluster_order,
                                "title": "Cluster",
                            },
                            "color": {
                                "field": "mean_expression",
                                "type": "quantitative",
                                "title": "Mean expression",
                                "scale": {"scheme": "blues", "zero": True},
                            },
                            "tooltip": [
                                {"field": "cluster", "type": "nominal"},
                                {"field": "gene", "type": "nominal"},
                                {
                                    "field": "mean_expression",
                                    "type": "quantitative",
                                    "format": ".3f",
                                },
                                {
                                    "field": "fraction_expressing",
                                    "type": "quantitative",
                                    "format": ".1%",
                                },
                            ],
                        },
                    },
                    use_container_width=True,
                )

                summary_csv = dot_data.to_csv(index=False).encode("utf-8")
                expression_summary_for_bundle = dot_data
                st.download_button(
                    "Download expression summary as CSV",
                    data=summary_csv,
                    file_name="expression_summary.csv",
                    mime="text/csv",
                    key="expression_summary_download",
                )
                st.caption(
                    "Exports the selected genes across all clusters. "
                    "fraction_expressing ranges from 0 to 1; "
                    "mean_expression includes cells with zero expression."
                )

        result_csv = markers.to_csv(index=False).encode("utf-8")

        st.download_button(
            "Download results as CSV",
            data=result_csv,
            file_name="marker_gene_results.csv",
            mime="text/csv",
        )

        if "analysis_settings" in st.session_state:
            analysis_bundle = build_analysis_bundle(
                markers,
                st.session_state["analysis_settings"],
                expression_summary_for_bundle,
            )
            st.download_button(
                "Download analysis bundle as ZIP",
                data=analysis_bundle,
                file_name="marker_gene_analysis.zip",
                mime="application/zip",
                key="analysis_bundle_download",
            )

st.caption(
    "This tool ranks markers using log2 fold change and reports "
    "engine-specific Wilcoxon p-values with "
    "Benjamini-Hochberg adjustment."
)
