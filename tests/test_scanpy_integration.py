import anndata as ad
import numpy as np
import pandas as pd

from src.marker_finder import (
    RESULT_COLUMNS,
    find_markers_scanpy,
)


def test_find_markers_scanpy_ranks_expected_genes() -> None:
    adata = ad.AnnData(
        X=np.array(
            [
                [5.0, 0.0],
                [4.0, 0.0],
                [6.0, 0.0],
                [5.0, 0.0],
                [0.0, 5.0],
                [0.0, 4.0],
                [0.0, 6.0],
                [0.0, 5.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "cluster": pd.Categorical(
                    ["A", "A", "A", "A", "B", "B", "B", "B"]
                )
            },
            index=[f"cell_{index}" for index in range(8)],
        ),
        var=pd.DataFrame(index=["GENE_A", "GENE_B"]),
    )

    markers = find_markers_scanpy(
        adata,
        cluster_column="cluster",
        top_n=1,
    )

    observed = dict(
        zip(markers["cluster"], markers["gene"])
    )

    assert markers.columns.tolist() == RESULT_COLUMNS
    assert observed == {
        "A": "GENE_A",
        "B": "GENE_B",
    }
