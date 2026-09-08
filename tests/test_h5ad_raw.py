import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.marker_finder import build_parser, find_markers_scanpy, run


@pytest.fixture
def raw_input(tmp_path):
    obs = pd.DataFrame(
        {"cluster": pd.Categorical(["A"] * 4 + ["B"] * 4)},
        index=[f"cell_{i}" for i in range(8)],
    )
    values = np.log1p(
        [[5, 0], [4, 0], [6, 0], [5, 0],
         [0, 5], [0, 4], [0, 6], [0, 5]]
    )
    data = ad.AnnData(
        X=values.copy(),
        obs=obs,
        var=pd.DataFrame(index=["X_A", "X_B"]),
    )
    raw_values = np.column_stack([values, np.log1p([1, 2, 1, 2] * 2)])
    data.raw = ad.AnnData(
        X=sparse.csr_matrix(raw_values),
        obs=obs.copy(),
        var=pd.DataFrame(index=["RAW_A", "RAW_B", "RAW_EXTRA"]),
    )
    data.layers["alternate"] = values.copy()
    path = tmp_path / "raw.h5ad"
    data.write_h5ad(path)
    return path


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_raw_uses_its_own_genes(raw_input, tmp_path, engine):
    default = run(
        raw_input, tmp_path / "default.csv", engine=engine, top_n=1
    )
    selected = run(
        raw_input, tmp_path / "raw.csv",
        engine=engine, use_raw=True, top_n=1,
    )
    assert dict(zip(default["cluster"], default["gene"])) == {
        "A": "X_A", "B": "X_B",
    }
    assert dict(zip(selected["cluster"], selected["gene"])) == {
        "A": "RAW_A", "B": "RAW_B",
    }
    np.testing.assert_allclose(
        selected["cluster_mean"],
        np.log1p([5, 4, 6, 5]).mean(),
    )
    np.testing.assert_allclose(selected["other_mean"], 0)
    pd.testing.assert_frame_equal(
        pd.read_csv(tmp_path / "raw.csv"),
        selected,
        check_dtype=False,
        check_exact=False,
        rtol=1e-6,
        atol=1e-8,
    )


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_missing_raw_rejected(raw_input, tmp_path, engine):
    data = ad.read_h5ad(raw_input)
    data.raw = None
    data.write_h5ad(raw_input)
    output = tmp_path / "missing.csv"
    with pytest.raises(ValueError, match="does not contain raw"):
        run(raw_input, output, engine=engine, use_raw=True)
    assert not output.exists()


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_raw_and_layer_rejected(raw_input, tmp_path, engine):
    output = tmp_path / "conflict.csv"
    with pytest.raises(ValueError, match="mutually exclusive"):
        run(
            raw_input, output, engine=engine,
            use_raw=True, layer="alternate",
        )
    assert not output.exists()


@pytest.mark.parametrize("engine", ["native", "scanpy"])
def test_raw_rejects_csv(tmp_path, engine):
    source = tmp_path / "input.csv"
    source.write_text("cluster,GENE\nA,5\nB,0\n")
    output = tmp_path / "output.csv"
    with pytest.raises(ValueError, match="Raw selection requires H5AD"):
        run(source, output, engine=engine, use_raw=True)
    assert not output.exists()


def test_scanpy_raw_preserves_input(raw_input):
    data = ad.read_h5ad(raw_input)
    original_x = data.X.copy()
    original_raw = data.raw.X.copy()
    original_obs = data.obs.copy()
    original_var = data.var.copy()

    find_markers_scanpy(data, use_raw=True, top_n=1)

    np.testing.assert_array_equal(data.X, original_x)
    assert sparse.issparse(data.raw.X)
    assert (data.raw.X != original_raw).nnz == 0
    pd.testing.assert_frame_equal(data.obs, original_obs)
    pd.testing.assert_frame_equal(data.var, original_var)
    assert "rank_genes_groups" not in data.uns


def test_parser_accepts_raw():
    args = build_parser().parse_args([
        "--input", "input.h5ad", "--output", "output.csv", "--use-raw",
    ])
    assert args.use_raw is True
    assert args.layer is None


def test_parser_rejects_raw_with_layer():
    with pytest.raises(SystemExit):
        build_parser().parse_args([
            "--input", "input.h5ad", "--output", "output.csv",
            "--use-raw", "--layer", "alternate",
        ])
