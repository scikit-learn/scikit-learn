import numpy as np
import pytest

from sklearn.utils.testing import (assert_array_equal, assert_equal,
                                   assert_raises)

from scipy.sparse import bsr_matrix, csc_matrix, csr_matrix

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
        assert (len(data), 1) == X.shape


def test_zero_variance_floating_point_error():
    # Test that VarianceThreshold(0.0).fit eliminates features that have
    # the same value in every sample, even when floating point errors
    # cause np.var not to be 0 for the feature.
    # See #13691

    data = [[-0.13725701]] * 10
    if np.var(data) == 0:
        pytest.skip('This test is not valid for this platform, as it relies '
                    'on numerical instabilities.')
    for X in [data, csr_matrix(data), csc_matrix(data), bsr_matrix(data)]:
        msg = "No feature in X meets the variance threshold 0.00000"
        with pytest.raises(ValueError, match=msg):
            VarianceThreshold().fit(X)
