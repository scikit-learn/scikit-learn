import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse import csr_matrix

from sklearn import config_context
from sklearn.utils._array_out import _make_array_out
from sklearn.utils._array_out import _get_feature_names


@pytest.mark.parametrize("X", [
   np.array([[1, 2, 3], [4, 5, 6]]),
   [[1, 2, 3], [4, 5, 6]],
   csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'list', 'sparse'])
def test_feature_names_no_names(X):
    assert _get_feature_names(X) is None


@pytest.mark.parametrize("column_name", [
    "columns", "my_special_name"
])
def test_feature_names_xarray(column_name):
    # the column names will always be the second axes
    xr = pytest.importorskip("xarray")
    X = np.array([[1, 2, 3], [4, 5, 6]])
    feature_names = [f"feature_{i}" for i in range(3)]
    X = xr.DataArray(X, dims=("index", column_name),
                     coords={column_name: feature_names})

    names = _get_feature_names(X)
    assert_array_equal(names, feature_names)


def test_feature_names_pandas():
    pd = pytest.importorskip("pandas")
    X = np.array([[1, 2, 3], [4, 5, 6]])
    feature_names = [f"feature_{i}" for i in range(3)]
    X = pd.DataFrame(X, columns=feature_names)

    names = _get_feature_names(X)
    assert_array_equal(names, feature_names)


@pytest.mark.parametrize("X_out", [
    np.array([[1, 2, 3], [2, 3, 4]]),
    csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'sparse'])
def test_make_array_out_default(X_out):
    # array_out='default' is an noop
    X_orig = np.ones((2, 10))
    with config_context(array_out='default'):
        out = _make_array_out(X_out, X_orig,
                              lambda names: names)
        assert out is X_out


def test_make_array_out_pandas():
    pd = pytest.importorskip("pandas")


def test_make_array_out_xarray():
    xr = pytest.importorskip("xarray")
