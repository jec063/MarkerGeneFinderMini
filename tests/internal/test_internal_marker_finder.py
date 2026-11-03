import numpy as np
from marker_gene_finder.statistics import _benjamini_hochberg


def test_benjamini_hochberg_adjustment() -> None:
    p_values = np.array([0.01, 0.04, 0.03])

    adjusted = _benjamini_hochberg(p_values)

    np.testing.assert_allclose(
        adjusted,
        np.array([0.03, 0.04, 0.04]),
    )
