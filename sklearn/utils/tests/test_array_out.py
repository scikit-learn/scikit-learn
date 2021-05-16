import numpy as np
import pytest
from scipy.sparse import csr_matrix
from numpy.testing import assert_array_equal

from sklearn.utils._array_out import _get_feature_names
from sklearn.utils._array_out import _get_index
from sklearn.utils._array_out import _make_array_out
from sklearn.utils._testing import assert_allclose_dense_sparse


@pytest.mark.parametrize("X", [
   np.array([[1, 2, 3], [4, 5, 6]]),
   [[1, 2, 3], [4, 5, 6]],
   csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'list', 'sparse'])
def test_feature_names_no_names(X):
    assert _get_feature_names(X) is None


@pytest.mark.parametrize("X", [
   np.array([[1, 2, 3], [4, 5, 6]]),
   [[1, 2, 3], [4, 5, 6]],
   csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'list', 'sparse'])
def test_get_index(X):
    assert _get_index(X) is None


@pytest.mark.parametrize("feature_names", [
    ["feat_0", "feat_1", "feat_2"],
    [1, 0, 2],
])
def test_feature_names_pandas(feature_names):
    pd = pytest.importorskip("pandas")
    X = np.array([[1, 2, 3], [4, 5, 6]])
    X = pd.DataFrame(X, columns=feature_names)

    names = _get_feature_names(X)
    assert_array_equal(names, feature_names)


@pytest.mark.parametrize("index", [
    [0, 1], ["a", "b"]
])
def test_get_feature_names_pandas(index):
    pd = pytest.importorskip("pandas")
    X = np.array([[1, 2, 3], [4, 5, 6]])
    X = pd.DataFrame(X, index=index)

    index_names = _get_index(X)
    assert_array_equal(index_names, index)


@pytest.mark.parametrize("X_out", [
    np.array([[1, 2, 3], [2, 3, 4]]),
    csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'sparse'])
def test_make_array_out_default(X_out):
    out = _make_array_out(X_out, None, lambda: None, array_out="default")
    assert out is X_out


@pytest.mark.parametrize("X_out", [
    np.array([[1, 2, 3], [2, 3, 4]]),
    csr_matrix([[1, 0, 0], [0, 0, 1]])
], ids=['ndarray', 'sparse'])
def test_make_array_out_error(X_out):
    msg = "array_out must be 'default' or 'pandas'"
    with pytest.raises(ValueError, match=msg):
        _make_array_out(X_out, None, lambda: None, array_out="bad")


@pytest.mark.parametrize("is_sparse", [True, False])
@pytest.mark.parametrize("out_features, expected_columns", [
    (['feat_1', 'feat_2'], ["feat_1", "feat_2"]),
    ([0, 1], [0, 1]),
    (None, ["X0", "X1"]),
])
@pytest.mark.parametrize("index, expected_index", [
    ([2, 3, 1], [2, 3, 1]),
    (["a", "c", "d"], ["a", "c", "d"]),
    (None, [0, 1, 2]),
])
def test_make_array_out_pandas(is_sparse, out_features, expected_columns,
                               index, expected_index):
    pd = pytest.importorskip("pandas")

    X_out = np.array([[0, 1], [1, 0], [2, 3]])
    if is_sparse:
        X_out = csr_matrix(X_out)

    df_out = _make_array_out(X_out, index, lambda: out_features,
                             array_out="pandas")

    assert isinstance(df_out, pd.DataFrame)
    assert_array_equal(df_out.index, expected_index)
    assert_array_equal(df_out.columns, expected_columns)

    if is_sparse:
        unwrapped = df_out.sparse.to_coo()
    else:
        unwrapped = df_out.to_numpy()

    assert_allclose_dense_sparse(X_out, unwrapped)
