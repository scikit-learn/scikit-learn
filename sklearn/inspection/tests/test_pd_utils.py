import numpy as np
import pytest
from numpy.testing import assert_allclose

from sklearn.inspection._pd_utils import (
    _check_feature_names,
    _get_feature_index,
    _nanpercentile,
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


def test_nanpercentile():
    X = np.array([[0.0], [1.0], [2.0], [3.0], [4.0]])
    percentiles = _nanpercentile(X, 0, (0.0, 0.5, 1.0))
    assert_allclose(percentiles, [0.0, 2.0, 4.0])


def test_nanpercentile_ignores_nan():
    """Non-regression test for gh-34928: a single `np.nan` used to make every
    computed percentile `np.nan` when using `np.quantile` naively.
    """
    X = np.array([[np.nan], [1.0], [2.0], [3.0], [4.0]])
    percentiles = _nanpercentile(X, 0, (0.0, 0.5, 1.0))
    assert not np.isnan(percentiles).any()
    assert_allclose(percentiles, [1.0, 2.5, 4.0])


def test_nanpercentile_boolean_column():
    """Non-regression test for gh-34928: `np.nanquantile` raises a `TypeError`
    on boolean arrays.
    """
    X = np.array([[True], [False], [True], [False]])
    percentiles = _nanpercentile(X, 0, (0.0, 1.0))
    assert_allclose(percentiles, [0.0, 1.0])
