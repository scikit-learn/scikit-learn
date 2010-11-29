# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np
from scipy import sparse
from scipy import linalg

from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal

from scikits.learn.utils.extmath import fast_svd
from scikits.learn.datasets.samples_generator import low_rank_fat_tail


def test_fast_svd():
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
    assert_almost_equal(s[:rank], sa[:rank])

    # check the singular vectors too (while not checking the sign)
    assert_almost_equal(np.dot(U[:, :rank], V[:rank, :]),
                        np.dot(Ua[:, :rank], Va[:rank, :]))

    # check the sparse matrix representation
    X = sparse.csr_matrix(X)

    # compute the singular values of X using the fast approximate method
    Ua, sa, Va = fast_svd(X, k)
    assert_almost_equal(s[:rank], sa[:rank])


def test_fast_svd_with_noise():
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
    assert np.abs(s[:rank] - sa[:rank]).max() > 0.1

    # compute the singular values of X using the fast approximate method with
    # iterated power method
    _, sap, _ = fast_svd(X, k, q=3)

    # the iterated power method is helping getting rid of the noise:
    assert_almost_equal(s[:rank], sap[:rank], decimal=5)


