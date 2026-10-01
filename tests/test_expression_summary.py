import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.expression_summary import summarize_clusters, summarize_expression


def test_cluster_summary_reports_counts_and_percentages():
    result = summarize_clusters(pd.Series(["B", "A", "B", "B"]))

    assert result["cluster"].tolist() == ["A", "B"]
    assert result["cells"].tolist() == [1, 3]
    np.testing.assert_allclose(result["percent_of_cells"], [25, 75])


@pytest.mark.parametrize("source", ["table", "dense", "sparse", "layer", "raw"])
def test_summary_uses_selected_matrix(source):
    values = np.array([[0., 2.], [4., 0.], [8., 6.], [0., 2.]])
    obs = pd.DataFrame(
        {"cluster": ["A", "A", "B", "B"]},
        index=["c0", "c1", "c2", "c3"],
    )
    genes = ["G1", "G2"]
    options = {}

    if source == "table":
        data = pd.DataFrame(values, columns=genes)
        data["cluster"] = obs["cluster"].to_numpy()
    else:
        data = ad.AnnData(
            X=values.copy(), obs=obs,
            var=pd.DataFrame(index=genes),
        )
        if source == "sparse":
            data.X = sparse.csr_matrix(values)
        elif source == "layer":
            data.X = np.zeros_like(values)
            data.layers["chosen"] = sparse.csr_matrix(values)
            options["layer"] = "chosen"
        elif source == "raw":
            genes = ["RAW1", "RAW2"]
            data.raw = ad.AnnData(
                X=sparse.csr_matrix(
                    np.column_stack([values, np.ones(4)])
                ),
                obs=obs.copy(),
                var=pd.DataFrame(index=genes + ["EXTRA"]),
            )
            data.X = np.zeros_like(values)
            options["use_raw"] = True

    result = summarize_expression(data, genes, **options)
    assert result["cluster"].tolist() == ["A", "A", "B", "B"]
    assert result["gene"].tolist() == genes * 2
    np.testing.assert_allclose(result["mean_expression"], [2, 1, 4, 4])
    np.testing.assert_allclose(
        result["fraction_expressing"], [.5, .5, .5, 1.]
    )


def test_conflicting_sources_rejected():
    with pytest.raises(ValueError, match="mutually exclusive"):
        summarize_expression(
            ad.AnnData(np.ones((2, 1))),
            ["G1"], layer="chosen", use_raw=True,
        )
