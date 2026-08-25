"""Streamlit interface for the CellCoPilot Marker Gene Finder."""

from __future__ import annotations

import io

import pandas as pd
import streamlit as st

from src.marker_finder import find_markers


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


@st.cache_data
def load_csv(file_contents: bytes) -> pd.DataFrame:
    """Load an uploaded CSV file."""
    return pd.read_csv(io.BytesIO(file_contents))


uploaded_file = st.file_uploader(
    "Upload an expression CSV",
    type=["csv"],
    help=(
        "The file should contain one row per cell, a cluster column, "
        "an optional cell-ID column, and numeric gene-expression columns."
    ),
)

if uploaded_file is None:
    st.info("Upload a CSV file to begin.")
    st.stop()

try:
    expression = load_csv(uploaded_file.getvalue())
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

st.subheader("Data preview")
st.dataframe(expression.head(20), width="stretch")

cluster_index = columns.index("cluster") if "cluster" in columns else 0

with st.sidebar:
    st.header("Analysis settings")

    cluster_column = st.selectbox(
        "Cluster column",
        options=columns,
        index=cluster_index,
    )

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
        options=cell_options,
        index=cell_index,
    )

    cell_column = (
        None
        if selected_cell_column == "None"
        else selected_cell_column
    )

    top_n = st.number_input(
        "Markers per cluster",
        min_value=1,
        max_value=100,
        value=5,
        step=1,
    )

    min_pct = st.slider(
        "Minimum fraction expressing gene",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.05,
    )

    min_log2fc = st.number_input(
        "Minimum log2 fold change",
        min_value=0.0,
        value=0.0,
        step=0.1,
        format="%.2f",
    )

    max_p_adj = st.number_input(
        "Maximum adjusted p-value",
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

    pseudocount = st.number_input(
        "Pseudocount",
        min_value=0.0001,
        value=0.1,
        step=0.1,
        format="%.4f",
    )

metric_columns = st.columns(3)
metric_columns[0].metric("Cells", len(expression))
metric_columns[1].metric("Columns", len(expression.columns))
metric_columns[2].metric(
    "Clusters",
    expression[cluster_column].nunique(),
)

if st.button("Find marker genes", type="primary"):
    try:
        with st.spinner("Finding marker genes..."):
            markers = find_markers(
                expression,
                cluster_column=cluster_column,
                cell_column=cell_column,
                top_n=int(top_n),
                pseudocount=float(pseudocount),
                min_pct=float(min_pct),
                min_log2fc=float(min_log2fc),
                max_p_adj=float(max_p_adj),
            )

        st.session_state["marker_results"] = markers

    except ValueError as error:
        st.error(str(error))

if "marker_results" in st.session_state:
    markers = st.session_state["marker_results"]

    st.subheader("Marker-gene results")

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

        result_csv = markers.to_csv(index=False).encode("utf-8")

        st.download_button(
            "Download results as CSV",
            data=result_csv,
            file_name="marker_gene_results.csv",
            mime="text/csv",
        )

st.caption(
    "This tool ranks markers using log2 fold change and reports "
    "one-sided Wilcoxon-Mann-Whitney p-values with "
    "Benjamini-Hochberg adjustment."
)
