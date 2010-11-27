# Author: Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD

import numpy as np
from scipy import sparse
from scipy import linalg

from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal

from ..extmath import fast_svd


def test_fast_svd():
    """Check that extmath.fast_svd is consistent with linalg.svd"""
    n_samples = 100
    n_features = 500
    rank = 5
    k = 100

    # generate a matrix X of rank `rank`
    np.random.seed(42)
    X = np.dot(np.random.randn(n_samples, rank),
               np.random.randn(rank, n_features))
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


