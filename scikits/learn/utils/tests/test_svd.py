# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np
from scipy import sparse
from scipy import linalg

from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal

from scikits.learn.utils.extmath import fast_svd
from scikits.learn.datasets.samples_generator import low_rank_fat_tail


def test_fast_svd_low_rank():
    """Check that extmath.fast_svd is consistent with linalg.svd"""
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # generate a matrix X of approximate effective rank `rank` and no noise
    # component (very structured signal):
    X = low_rank_fat_tail(n_samples, n_features, effective_rank=rank,
                          tail_strength=0.0, seed=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    U, s, V = linalg.svd(X, full_matrices=False)

    # compute the singular values of X using the fast approximate method
    Ua, sa, Va = fast_svd(X, k)
    assert_equal(Ua.shape, (n_samples, k))
    assert_equal(sa.shape, (k,))
    assert_equal(Va.shape, (k, n_features))

    # ensure that the singular values of both methods are equal up to the real
    # rank of the matrix
    assert_almost_equal(s[:k], sa)

    # check the singular vectors too (while not checking the sign)
    assert_almost_equal(np.dot(U[:, :k], V[:k, :]), np.dot(Ua, Va))

    # check the sparse matrix representation
    X = sparse.csr_matrix(X)

    # compute the singular values of X using the fast approximate method
    Ua, sa, Va = fast_svd(X, k)
    assert_almost_equal(s[:rank], sa[:rank])


def test_fast_svd_low_rank_with_noise():
    """Check that extmath.fast_svd can handle noisy matrices"""
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # generate a matrix X wity structure approximate rank `rank` and an
    # important noisy component
    X = low_rank_fat_tail(n_samples, n_features, effective_rank=rank,
                          tail_strength=0.5, seed=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    _, s, _ = linalg.svd(X, full_matrices=False)

    # compute the singular values of X using the fast approximate method without
    # the iterated power method
    _, sa, _ = fast_svd(X, k, q=0)

    # the approximation does not tolerate the noise:
    assert np.abs(s[:k] - sa).max() > 0.05

    # compute the singular values of X using the fast approximate method with
    # iterated power method
    _, sap, _ = fast_svd(X, k, q=5)

    # the iterated power method is helping getting rid of the noise:
    assert_almost_equal(s[:k], sap, decimal=3)


def test_fast_svd_infinite_rank():
    """Check that extmath.fast_svd can handle noisy matrices"""
    n_samples = 100
    n_features = 500
    rank = 5
    k = 10

    # let us try again without 'low_rank component': just regularly but slowly
    # decreasing singular values: the rank of the data matrix is infinite
    X = low_rank_fat_tail(n_samples, n_features, effective_rank=rank,
                          tail_strength=1.0, seed=0)
    assert_equal(X.shape, (n_samples, n_features))

    # compute the singular values of X using the slow exact method
    _, s, _ = linalg.svd(X, full_matrices=False)

    # compute the singular values of X using the fast approximate method without
    # the iterated power method
    _, sa, _ = fast_svd(X, k, q=0)

    # the approximation does not tolerate the noise:
    assert np.abs(s[:k] - sa).max() > 0.1

    # compute the singular values of X using the fast approximate method with
    # iterated power method
    _, sap, _ = fast_svd(X, k, q=5)

    # the iterated power method is still managing to get most of the structure
    # at the requested rank
    assert_almost_equal(s[:k], sap, decimal=3)


def test_fast_svd_transpose_consistency():
    """Check that transposing the design matrix has limit impact"""
    n_samples = 100
    n_features = 500
    rank = 4
    k = 10

    X = low_rank_fat_tail(n_samples, n_features, effective_rank=rank,
                          tail_strength=0.5, seed=0)
    assert_equal(X.shape, (n_samples, n_features))

    U1, s1, V1 = fast_svd(X, k, q=3, transpose=False, random_state=0)
    U2, s2, V2 = fast_svd(X, k, q=3, transpose=True, random_state=0)
    U3, s3, V3 = fast_svd(X, k, q=3, transpose='auto', random_state=0)
    U4, s4, V4 = linalg.svd(X, full_matrices=False)

    assert_almost_equal(s1, s4[:k], decimal=3)
    assert_almost_equal(s2, s4[:k], decimal=3)
    assert_almost_equal(s3, s4[:k], decimal=3)

    assert_almost_equal(np.dot(U1, V1), np.dot(U4[:, :k], V4[:k, :]),
                        decimal=2)
    assert_almost_equal(np.dot(U2, V2), np.dot(U4[:, :k], V4[:k, :]),
                        decimal=2)


    # in this case 'auto' is equivalent to transpose
    assert_almost_equal(s2, s3)


if __name__ == '__main__':
    import nose
    nose.runmodule()
