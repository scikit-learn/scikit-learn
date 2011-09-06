import numpy as np
from numpy import linalg
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises
from nose.tools import assert_true
from scipy.sparse import csr_matrix

from ..pairwise import euclidean_distances, linear_kernel, polynomial_kernel, \
                       rbf_kernel, sigmoid_kernel
from .. import pairwise_distances, pairwise_kernels
from ..pairwise import pairwise_distance_functions, pairwise_kernel_functions

np.random.seed(0)


def test_pairwise_distances():
    """ Test the pairwise_distance helper function. """
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    S = pairwise_distances(X, metric="euclidean")
    S2 = euclidean_distances(X)
    assert_array_equal(S, S2)

    Y = rng.random_sample((2, 4))
    S = pairwise_distances(X, Y, metric="euclidean")
    S2 = euclidean_distances(X, Y)
    assert_array_equal(S, S2)

    S = pairwise_distances(X, metric="cityblock")
    assert_equal(S.shape[0], S.shape[1])
    assert_equal(S.shape[0], X.shape[0])

    S = pairwise_distances(X, Y, metric="cityblock")
    assert_equal(S.shape[0], X.shape[0])
    assert_equal(S.shape[1], Y.shape[0])

    S = np.dot(X, X.T)
    S2 = pairwise_distances(S, metric="precomputed")
    assert_true(S is S2)
    assert_raises(ValueError, pairwise_distances, X, None, "precomputed")


def test_pairwise_kernels():
    """ Test the pairwise_kernels helper function. """
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((2, 4))
    # Test with all metrics that should be in pairwise_kernel_functions.
    test_metrics = ["rbf", "sigmoid", "polynomial", "linear"]
    for metric in test_metrics:
        function = pairwise_kernel_functions[metric]
        # Test with Y=None
        K1 = pairwise_kernels(X, metric=metric)
        K2 = function(X)
        assert_equal(K1, K2)
        # Test with Y=Y
        K1 = pairwise_kernels(X, Y=Y, metric=metric)
        K2 = function(X, Y=Y)
        assert_equal(K1, K2)
    # Test with a callable function, with given keywords.
    metric = callable_rbf_kernel
    kwds = {}
    kwds['gamma'] = 0.5
    K1 = pairwise_kernels(X, Y=Y, metric=metric, **kwds)
    K2 = rbf_kernel(X, Y=Y, **kwds)    
    assert_equal(K1, K2)


def callable_rbf_kernel(x, y, **kwds):
    """ Callable version of pairwise.rbf_kernel. """
    return rbf_kernel(np.atleast_2d(x), np.atleast_2d(y), **kwds)


def test_euclidean_distances():
    """Check the pairwise Euclidean distances computation"""
    X = [[0]]
    Y = [[1], [2]]
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])

    X = csr_matrix(X)
    Y = csr_matrix(Y)
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])


def test_kernel_symmetry():
    """valid kernels should be symmetric"""
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   sigmoid_kernel):
        K = kernel(X, X)
        assert_array_almost_equal(K, K.T, 15)


def test_kernel_sparse():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    X_sparse = csr_matrix(X)
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   sigmoid_kernel):
        K = kernel(X, X)
        K2 = kernel(X_sparse, X_sparse)
        assert_array_almost_equal(K, K2)


def test_linear_kernel():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    K = linear_kernel(X, X)
    # the diagonal elements of a linear kernel are their squared norm
    assert_array_almost_equal(K.flat[::6], [linalg.norm(x) ** 2 for x in X])


def test_rbf_kernel():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    K = rbf_kernel(X, X)
    # the diagonal elements of a rbf kernel are 1
    assert_array_almost_equal(K.flat[::6], np.ones(5))
