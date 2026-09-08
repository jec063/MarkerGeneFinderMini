"""Streamlit interface for the CellCoPilot Marker Gene Finder."""

from __future__ import annotations

import io
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd
import streamlit as st

from src.marker_finder import (
    find_markers,
    find_markers_scanpy,
    load_expression,
)


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

if uploaded_file is None:
    st.info("Upload a CSV or H5AD file to begin.")
    st.stop()

file_contents = uploaded_file.getvalue()
file_suffix = Path(uploaded_file.name).suffix.lower()

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
            uploaded_file.name,
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
                uploaded_file.name,
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
                uploaded_file.name,
            )
            expression = None
        else:
            expression = load_uploaded_h5ad(
                file_contents,
                uploaded_file.name,
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

if st.button("Find marker genes", type="primary"):
    try:
        with st.spinner("Finding marker genes..."):
            if analysis_engine == "scanpy":
                markers = find_markers_scanpy(
                    adata,
                    layer=selected_layer,
                    use_raw=use_raw,
                    cluster_column=cluster_column,
                    top_n=int(top_n),
                    min_pct=float(min_pct),
                    min_log2fc=float(min_log2fc),
                    max_p_adj=float(max_p_adj),
                )
            else:
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
    "engine-specific Wilcoxon p-values with "
    "Benjamini-Hochberg adjustment."
)
