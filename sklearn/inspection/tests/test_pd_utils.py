import numpy as np
import pytest

from sklearn.inspection._pd_utils import (
    _check_feature_names,
    _get_feature_index,
    _get_is_categorical,
)
from sklearn.utils._testing import _convert_container


@pytest.mark.parametrize(
    "feature_names, array_type, expected_feature_names",
    [
        (None, "array", ["x0", "x1", "x2"]),
        (None, "pandas", ["a", "b", "c"]),
        (np.array(["a", "b", "c"]), "array", ["a", "b", "c"]),
    ],
)
def test_check_feature_names(feature_names, array_type, expected_feature_names):
    X = np.random.randn(10, 3)
    column_names = ["a", "b", "c"]
    X = _convert_container(X, constructor_name=array_type, column_names=column_names)
    feature_names_validated = _check_feature_names(X, feature_names)
    assert feature_names_validated == expected_feature_names


def test_check_feature_names_error():
    X = np.random.randn(10, 3)
    feature_names = ["a", "b", "c", "a"]
    msg = "feature_names should not contain duplicates."
    with pytest.raises(ValueError, match=msg):
        _check_feature_names(X, feature_names)


@pytest.mark.parametrize("fx, idx", [(0, 0), (1, 1), ("a", 0), ("b", 1), ("c", 2)])
def test_get_feature_index(fx, idx):
    feature_names = ["a", "b", "c"]
    assert _get_feature_index(fx, feature_names) == idx


@pytest.mark.parametrize(
    "fx, feature_names, err_msg",
    [
        ("a", None, "Cannot plot partial dependence for feature 'a'"),
        ("d", ["a", "b", "c"], "Feature 'd' not in feature_names"),
    ],
)
def test_get_feature_names_error(fx, feature_names, err_msg):
    with pytest.raises(ValueError, match=err_msg):
        _get_feature_index(fx, feature_names)


@pytest.mark.parametrize(
    "categorical_features, expected",
    [
        (None, [False, False]),
        ([True, False, True], [False, True]),
        ([2], [False, True]),
        (["c"], [False, True]),
    ],
)
def test_get_is_categorical(categorical_features, expected):
    feature_names = ["a", "b", "c"]
    is_categorical = _get_is_categorical(
        categorical_features, [1, 2], feature_names, n_features=3
    )
    assert is_categorical == expected


@pytest.mark.parametrize(
    "categorical_features, err_msg",
    [
        ([], "Passing an empty list"),
        ([True, False], "the array should be of shape"),
        ([1.5], "Expected `categorical_features` to be an array-like"),
    ],
)
def test_get_is_categorical_error(categorical_features, err_msg):
    with pytest.raises(ValueError, match=err_msg):
        _get_is_categorical(categorical_features, [0], ["a", "b", "c"], n_features=3)
