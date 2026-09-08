import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.marker_finder import build_parser, run


@pytest.fixture
def layered_input(tmp_path):
    matrix = np.log1p([
        [5., 0.],
        [4., 0.],
        [6., 0.],
        [5., 0.],
        [0., 5.],
        [0., 4.],
        [0., 6.],
        [0., 5.],
    ])
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
            index=[f"cell_{i}" for i in range(8)],
        ),
        var=pd.DataFrame(index=["GENE_A", "GENE_B"]),
    )
    # Reverse the markers in a sparse layer.
    adata.layers["alternate"] = sparse.csr_matrix(matrix[:, ::-1])
    path = tmp_path / "layered.h5ad"
    adata.write_h5ad(path)
    return path


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_layer_selection_changes_markers(layered_input, tmp_path, engine):
    default = run(
        layered_input,
        tmp_path / "default.csv",
        engine=engine,
        top_n=1,
    )
    selected = run(
        layered_input,
        tmp_path / "selected.csv",
        engine=engine,
        layer="alternate",
        top_n=1,
    )

    assert dict(zip(default["cluster"], default["gene"])) == {
        "A": "GENE_A", "B": "GENE_B",
    }
    assert dict(zip(selected["cluster"], selected["gene"])) == {
        "A": "GENE_B", "B": "GENE_A",
    }
    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "selected.csv"),
        selected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_missing_layer_is_rejected(layered_input, tmp_path, engine):
    output = tmp_path / "missing.csv"
    with pytest.raises(ValueError, match="Missing AnnData layer"):
        run(
            layered_input,
            output,
            engine=engine,
            layer="missing",
        )
    assert not output.exists()


def test_csv_rejects_layer_selection(tmp_path):
    source = tmp_path / "expression.csv"
    source.write_text("cluster,GENE_A\nA,5\nB,0\n")
    output = tmp_path / "markers.csv"

    with pytest.raises(ValueError, match="Layer selection requires H5AD"):
        run(source, output, layer="alternate")
    assert not output.exists()


def test_parser_accepts_layer():
    args = build_parser().parse_args([
        "--input", "expression.h5ad",
        "--output", "markers.csv",
        "--layer", "alternate",
    ])
    assert args.layer == "alternate"
