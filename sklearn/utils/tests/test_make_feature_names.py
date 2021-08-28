import pytest

import numpy as np
from numpy.testing import assert_array_equal
from sklearn.utils._feature_names import _make_feature_names


@pytest.mark.parametrize(
    "n_features, prefix, input_features, expected_names",
    [
        (3, "x", None, ["x0", "x1", "x2"]),
        (4, "x", ["cat", "dog", "snake"], ["cat", "dog", "snake"]),
        (4, "pca", None, ["pca0", "pca1", "pca2", "pca3"]),
    ],
)
def test_make_feature_names(n_features, prefix, input_features, expected_names):
    feature_names = _make_feature_names(
        n_features=n_features, prefix=prefix, input_features=input_features
    )
    assert isinstance(feature_names, np.ndarray)
    assert_array_equal(expected_names, feature_names)
