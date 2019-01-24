from sklearn.utils.testing import (assert_array_equal, assert_equal,
                                   assert_raises)

from scipy.sparse import bsr_matrix, csc_matrix, csr_matrix
import numpy as np

from sklearn.feature_selection import VarianceThreshold

data = [[0, 1, 2, 3, 4],
        [0, 2, 2, 3, 5],
        [1, 1, 2, 4, 0]]


def test_zero_variance():
    # Test VarianceThreshold with default setting, zero variance.

    for X in [data, csr_matrix(data), csc_matrix(data), bsr_matrix(data)]:
        sel = VarianceThreshold().fit(X)
        assert_array_equal([0, 1, 3, 4], sel.get_support(indices=True))

    assert_raises(ValueError, VarianceThreshold().fit, [[0, 1, 2, 3]])
    assert_raises(ValueError, VarianceThreshold().fit, [[0, 1], [0, 1]])


def test_variance_threshold():
    # Test VarianceThreshold with custom variance.
    for X in [data, csr_matrix(data)]:
        X = VarianceThreshold(threshold=.4).fit_transform(X)
        assert_equal((len(data), 1), X.shape)


def test_variance_nan():
    arr = np.array(data, dtype=np.float64)
    # add single NaN and feature should still be included
    arr[0, 0] = np.NaN
    # make all values in feature NaN and feature should be rejected
    arr[:, 1] = np.NaN

    for X in [arr, csr_matrix(arr), csc_matrix(arr), bsr_matrix(arr)]:
        sel = VarianceThreshold().fit(X)
        assert_array_equal([0, 3, 4], sel.get_support(indices=True))
