import numpy as np
import pytest

from sklearn.neighbors._kd_tree import KDTree

DIMENSION = 3

METRICS = {'euclidean': {},
           'manhattan': {},
           'chebyshev': {},
           'minkowski': dict(p=3)}


def test_array_object_type():
    """Check that we do not accept object dtype array."""
    X = np.array([(1, 2, 3), (2, 5), (5, 5, 1, 2)], dtype=object)
    with pytest.raises(
        ValueError,
        match="setting an array element with a sequence"
    ):
        KDTree(X)
