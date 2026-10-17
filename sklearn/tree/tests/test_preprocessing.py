import re

import numpy as np
import pytest

from sklearn.tree._preprocessing import _check_categorical_features
from sklearn.utils._testing import _convert_container, assert_allclose


@pytest.mark.parametrize(
    ["categorical_features", "expected_msg"],
    [
        (
            [b"hello", b"world"],
            re.escape(
                "categorical_features must be an array-like of bool, int or str, "
                "got: bytes40."
            ),
        ),
        (
            np.array([b"hello", 1.3], dtype=object),
            re.escape(
                "categorical_features must be an array-like of bool, int or str, "
                "got: bytes, float."
            ),
        ),
        (
            [0, -1],
            re.escape(
                "categorical_features set as integer indices must be in "
                "[0, n_features - 1]"
            ),
        ),
        (
            [True, True, False, False, True],
            re.escape(
                "categorical_features set as a boolean mask must have shape "
                "(n_features,)"
            ),
        ),
    ],
)
def test_check_categorical_features_raises(categorical_features, expected_msg):
    """Test that check_categorical_features raises expected errors."""
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 10
    X = rng.randint(0, 3, size=(n_samples, n_features))

    with pytest.raises(ValueError, match=expected_msg):
        _check_categorical_features(X, categorical_features)


@pytest.mark.parametrize(
    ["categorical_features", "on_array"],
    [
        ([False, True, True, False], True),
        ([1, 2], True),
        (["b", "c"], False),
        ("from_dtype", False),
    ],
)
@pytest.mark.parametrize("constructor_name", ["array", "pandas", "polars"])
def test_check_categorical_features(categorical_features, on_array, constructor_name):
    """Test that check_categorical_features returns as expected on simple data."""
    rng = np.random.RandomState(0)
    n_samples, n_features = 30, 4
    X = rng.randint(0, 3, size=(n_samples, n_features))
    if constructor_name == "array" and not on_array:
        return
    elif constructor_name == "polars":
        X = X.astype(str)
    X = _convert_container(
        X,
        constructor_name,
        column_names=["a", "b", "c", "d"],
        categorical_feature_names=["b", "c"],
    )

    result = _check_categorical_features(X, categorical_features)
    assert_allclose(result, [False, True, True, False])
