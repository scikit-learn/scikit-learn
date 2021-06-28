import pytest

import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._feature_names import _get_feature_names


def test_get_feature_names_pandas():
    pd = pytest.importorskip("pandas")
    column_names = [f"col_{i}" for i in range(3)]
    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)
    X = pd.DataFrame(X, columns=column_names)

    names = _get_feature_names(X)
    assert_array_equal(names, column_names)


def test_get_feature_names_numpy():
    X = np.array([[1, 2, 3], [4, 5, 6]], dtype=float)
    names = _get_feature_names(X)
    assert names is None
