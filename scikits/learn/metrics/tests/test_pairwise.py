import numpy as np
from numpy import linalg
from numpy.testing import assert_array_almost_equal, assert_array_equal

from ..pairwise import euclidean_distances, linear_kernel, polynomial_kernel, \
                       rbf_kernel

np.random.seed(0)

def test_euclidean_distances():
    """Check that the pairwise euclidian distances computation"""
    X = [[0]]
    Y = [[1], [2]]
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])


def test_kernel_symmetry():
    """valid kernels should be symmetric"""
    X = np.random.random((5, 4))
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel):
        K = kernel(X, X)
        assert_array_equal(K, K.T)


def test_linear_kernel():
    X = np.random.random((5, 4))
    K = linear_kernel(X, X)
    # the diagonal elements of a linear kernel are their squared norm
    assert_array_almost_equal(K.flat[::6], [linalg.norm(x) ** 2 for x in X])


def test_rbf_kernel():
    X = np.random.random((5, 4))
    K = rbf_kernel(X, X)
    # the diagonal elements of a rbf kernel are 1
    assert_array_almost_equal(K.flat[::6], np.ones(5))
