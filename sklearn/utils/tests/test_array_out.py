import pytest
import numpy as np
from numpy.testing import assert_array_equal
from scipy.sparse import csr_matrix

from sklearn import config_context
from sklearn.utils._array_out import _make_array_out
from sklearn.utils._array_out import _get_feature_names
from sklearn.utils._testing import assert_allclose_dense_sparse


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
        out = _make_array_out(X_out, X_orig, lambda names: names)
        assert out is X_out


@pytest.mark.parametrize("X_out", [
    np.array([[1, 2, 3], [2, 3, 4]]),
    csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'sparse'])
def test_make_array_out_error(X_out):
    X_orig = np.ones((2, 10))
    msg = "array_out must be 'default', 'pandas' or 'xarray'"
    with config_context(array_out='bad'):
        with pytest.raises(ValueError, match=msg):
            _make_array_out(X_out, X_orig, lambda names: names)


@pytest.mark.parametrize("is_sparse", [True, False])
@pytest.mark.parametrize("out_features", [
    ['feat_1', 'feat_2'], None
])
def test_make_array_out_pandas(is_sparse, out_features):
    pd = pytest.importorskip("pandas")

    index = [2, 3]
    X_orig = pd.DataFrame(np.array([[1, 2, 3], [3, 4, 5]]),
                          columns=[f'col_{i}' for i in range(3)],
                          index=index)

    X_out = np.array([[0, 1], [1, 0]])
    if is_sparse:
        X_out = csr_matrix(X_out)

    with config_context(array_out="pandas"):
        df_out = _make_array_out(X_out, X_orig, lambda: out_features)

    assert isinstance(df_out, pd.DataFrame)
    assert_array_equal(df_out.index, index)

    if out_features is None:
        # default output feature names
        assert_array_equal(df_out.columns, ["X0", "X1"])
    else:
        assert_array_equal(df_out.columns, out_features)

    if is_sparse:
        unwrapped = df_out.sparse.to_coo()
    else:
        unwrapped = df_out.to_numpy()

    assert_allclose_dense_sparse(X_out, unwrapped)


def test_make_array_out_pandas_zero_features():
    pd = pytest.importorskip("pandas")

    index = [2, 3]
    X_orig = pd.DataFrame(np.array([[1, 2, 3], [3, 4, 5]]),
                          columns=[f'col_{i}' for i in range(3)],
                          index=index)

    X_out = np.zeros((2, 0))
    with config_context(array_out="pandas"):
        df_out = _make_array_out(X_out, X_orig, lambda: [])
    assert isinstance(df_out, pd.DataFrame)
    assert_array_equal(df_out.index, index)


@pytest.mark.parametrize("is_sparse", [True, False])
@pytest.mark.parametrize("out_features", [
    ['feat_1', 'feat_2'], None
])
def test_make_array_out_xarray(is_sparse, out_features):
    xr = pytest.importorskip("xarray")
    if is_sparse:
        pytest.importorskip("sparse")

    index = [2, 3]
    X_orig = xr.DataArray(np.array([[1, 2, 3], [3, 4, 5]]),
                          dims=("index", "columns"),
                          coords={"index": index,
                                  "columns": [f"col_{i}" for i in range(3)]})

    X_out = np.array([[0, 1], [1, 0]])
    if is_sparse:
        X_out = csr_matrix(X_out)

    with config_context(array_out="xarray"):
        df_out = _make_array_out(X_out, X_orig, lambda: out_features)

    assert isinstance(df_out, xr.DataArray)
    assert_array_equal(df_out.coords["index"], index)

    if out_features is None:
        # default output feature names
        assert_array_equal(df_out.coords["columns"], ["X0", "X1"])
    else:
        assert_array_equal(df_out.coords["columns"], out_features)

    unwrapped = df_out.data
    if is_sparse:
        unwrapped = unwrapped.to_scipy_sparse()

    assert_allclose_dense_sparse(X_out, unwrapped)


def test_make_array_out_xarray_zero_features():
    xr = pytest.importorskip("xarray")

    index = [2, 3]
    X_orig = xr.DataArray(np.array([[1, 2, 3], [3, 4, 5]]),
                          dims=("index", "columns"),
                          coords={"index": index,
                                  "columns": [f"col_{i}" for i in range(3)]})

    X_out = np.zeros((2, 0))
    with config_context(array_out="xarray"):
        df_out = _make_array_out(X_out, X_orig, lambda: [])
    assert isinstance(df_out, xr.DataArray)
    assert_array_equal(df_out.coords["index"], index)
