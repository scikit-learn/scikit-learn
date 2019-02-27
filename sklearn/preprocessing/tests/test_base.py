import numpy as np
import pytest
from scipy import sparse

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.preprocessing._encoders import _transform_selected
from sklearn.preprocessing.data import Binarizer


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def _check_transform_selected(X, X_expected, dtype, sel):
    for M in (X, sparse.csr_matrix(X)):
        Xtr = _transform_selected(M, Binarizer().transform, dtype, sel)
        assert_array_equal(toarray(Xtr), X_expected)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_transform_selected(output_dtype, input_dtype):
    X = np.asarray([[3, 2, 1], [0, 1, 1]], dtype=input_dtype)

    X_expected = np.asarray([[1, 2, 1], [0, 1, 1]], dtype=output_dtype)
    _check_transform_selected(X, X_expected, output_dtype, [0])
    _check_transform_selected(X, X_expected, output_dtype,
                              [True, False, False])

    X_expected = np.asarray([[1, 1, 1], [0, 1, 1]], dtype=output_dtype)
    _check_transform_selected(X, X_expected, output_dtype, [0, 1, 2])
    _check_transform_selected(X, X_expected, output_dtype, [True, True, True])
    _check_transform_selected(X, X_expected, output_dtype, "all")

    _check_transform_selected(X, X, output_dtype, [])
    _check_transform_selected(X, X, output_dtype, [False, False, False])


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_transform_selected_copy_arg(output_dtype, input_dtype):
    # transformer that alters X
    def _mutating_transformer(X):
        X[0, 0] = X[0, 0] + 1
        return X

    original_X = np.asarray([[1, 2], [3, 4]], dtype=input_dtype)
    expected_Xtr = np.asarray([[2, 2], [3, 4]], dtype=output_dtype)

    X = original_X.copy()
    Xtr = _transform_selected(X, _mutating_transformer, output_dtype,
                              copy=True, selected='all')

    assert_array_equal(toarray(X), toarray(original_X))
    assert_array_equal(toarray(Xtr), expected_Xtr)


def test_transform_selected_retain_order():
    X = [[-1, 1], [2, -2]]

    assert_raise_message(ValueError,
                         "The retain_order option can only be set to True "
                         "for dense matrices.",
                         _transform_selected, sparse.csr_matrix(X),
                         Binarizer().transform, dtype=np.int, selected=[0],
                         retain_order=True)

    def transform(X):
        return np.hstack((X, [[0], [0]]))

    assert_raise_message(ValueError,
                         "The retain_order option can only be set to True "
                         "if the dimensions of the input array match the "
                         "dimensions of the transformed array.",
                         _transform_selected, X, transform, dtype=np.int,
                         selected=[0], retain_order=True)

    X_expected = [[-1, 1], [2, 0]]
    Xtr = _transform_selected(X, Binarizer().transform, dtype=np.int,
                              selected=[1], retain_order=True)
    assert_array_equal(toarray(Xtr), X_expected)

    X_expected = [[0, 1], [1, -2]]
    Xtr = _transform_selected(X, Binarizer().transform, dtype=np.int,
                              selected=[0], retain_order=True)
    assert_array_equal(toarray(Xtr), X_expected)
