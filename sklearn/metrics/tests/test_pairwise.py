import numpy as np
from numpy import linalg
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises
from nose.tools import assert_true
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cosine, cityblock, minkowski

from ..pairwise import (euclidean_distances, linear_kernel, polynomial_kernel,
                        rbf_kernel, sigmoid_kernel)
from .. import pairwise_distances, pairwise_kernels
from ..pairwise import pairwise_kernel_functions
from ..pairwise import check_pairwise_arrays

np.random.seed(0)


def test_pairwise_distances():
    """ Test the pairwise_distance helper function. """
    rng = np.random.RandomState(0)
    # Euclidean distance should be equivalent to calling the function.
    X = rng.random_sample((5, 4))
    S = pairwise_distances(X, metric="euclidean")
    S2 = euclidean_distances(X)
    assert_array_almost_equal(S, S2)
    # Euclidean distance, with Y != X.
    Y = rng.random_sample((2, 4))
    S = pairwise_distances(X, Y, metric="euclidean")
    S2 = euclidean_distances(X, Y)
    assert_array_almost_equal(S, S2)
    # Test with tuples as X and Y
    X_tuples = tuple([tuple([v for v in row]) for row in X])
    Y_tuples = tuple([tuple([v for v in row]) for row in Y])
    S2 = pairwise_distances(X_tuples, Y_tuples, metric="euclidean")
    assert_array_almost_equal(S, S2)
    # "cityblock" uses sklearn metric, cityblock (function) is scipy.spatial.
    S = pairwise_distances(X, metric="cityblock")
    S2 = pairwise_distances(X, metric=cityblock)
    assert_equal(S.shape[0], S.shape[1])
    assert_equal(S.shape[0], X.shape[0])
    assert_array_almost_equal(S, S2)
    # The manhattan metric should be equivalent to cityblock.
    S = pairwise_distances(X, Y, metric="manhattan")
    S2 = pairwise_distances(X, Y, metric=cityblock)
    assert_equal(S.shape[0], X.shape[0])
    assert_equal(S.shape[1], Y.shape[0])
    assert_array_almost_equal(S, S2)
    # Test cosine as a string metric versus cosine callable
    S = pairwise_distances(X, Y, metric="cosine")
    S2 = pairwise_distances(X, Y, metric=cosine)
    assert_equal(S.shape[0], X.shape[0])
    assert_equal(S.shape[1], Y.shape[0])
    assert_array_almost_equal(S, S2)
    # Tests that precomputed metric returns pointer to, and not copy of, X.
    S = np.dot(X, X.T)
    S2 = pairwise_distances(S, metric="precomputed")
    assert_true(S is S2)
    assert_raises(ValueError, pairwise_distances, X, None, "precomputed")
    # Test with sparse X and Y
    X_sparse = csr_matrix(X)
    Y_sparse = csr_matrix(Y)
    S = pairwise_distances(X_sparse, Y_sparse, metric="euclidean")
    S2 = euclidean_distances(X_sparse, Y_sparse)
    assert_array_almost_equal(S, S2)
    # Test with scipy.spatial.distance metric, with a kwd
    kwds = {"p": 2.0}
    S = pairwise_distances(X, Y, metric="minkowski", **kwds)
    S2 = pairwise_distances(X, Y, metric=minkowski, **kwds)
    assert_array_almost_equal(S, S2)
    # Test that scipy distance metrics throw an error if sparse matrix given
    assert_raises(TypeError, pairwise_distances, X_sparse, metric="minkowski")
    assert_raises(TypeError, pairwise_distances, X, Y_sparse,
                  metric="minkowski")


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
        assert_array_almost_equal(K1, K2)
        # Test with Y=Y
        K1 = pairwise_kernels(X, Y=Y, metric=metric)
        K2 = function(X, Y=Y)
        assert_array_almost_equal(K1, K2)
        # Test with tuples as X and Y
        X_tuples = tuple([tuple([v for v in row]) for row in X])
        Y_tuples = tuple([tuple([v for v in row]) for row in Y])
        K2 = pairwise_kernels(X_tuples, Y_tuples, metric=metric)
        assert_array_almost_equal(K1, K2)
        # Test with sparse X and Y
        X_sparse = csr_matrix(X)
        Y_sparse = csr_matrix(Y)
        K1 = pairwise_kernels(X_sparse, Y=Y_sparse, metric=metric)
        assert_array_almost_equal(K1, K2)
    # Test with a callable function, with given keywords.
    metric = callable_rbf_kernel
    kwds = {}
    kwds['gamma'] = 0.
    K1 = pairwise_kernels(X, Y=Y, metric=metric, **kwds)
    K2 = rbf_kernel(X, Y=Y, **kwds)
    assert_array_almost_equal(K1, K2)


def callable_rbf_kernel(x, y, **kwds):
    """ Callable version of pairwise.rbf_kernel. """
    K = rbf_kernel(np.atleast_2d(x), np.atleast_2d(y), **kwds)
    return K


def test_euclidean_distances():
    """ Check the pairwise Euclidean distances computation"""
    X = [[0]]
    Y = [[1], [2]]
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])

    X = csr_matrix(X)
    Y = csr_matrix(Y)
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])


def test_kernel_symmetry():
    """ Valid kernels should be symmetric"""
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


def test_check_dense_matrices():
    """ Ensure that pairwise array check works for dense matrices."""
    # Check that if XB is None, XB is returned as reference to XA
    XA = np.resize(np.arange(40), (5, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, None)
    assert_true(XA_checked is XB_checked)
    assert_equal(XA, XA_checked)


def test_check_XB_returned():
    """ Ensure that if XA and XB are given correctly, they return as equal."""
    # Check that if XB is not None, it is returned equal.
    # Note that the second dimension of XB is the same as XA.
    XA = np.resize(np.arange(40), (5, 8))
    XB = np.resize(np.arange(32), (4, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, XB)
    assert_equal(XA, XA_checked)
    assert_equal(XB, XB_checked)


def test_check_different_dimensions():
    """ Ensure an error is raised if the dimensions are different. """
    XA = np.resize(np.arange(45), (5, 9))
    XB = np.resize(np.arange(32), (4, 8))
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)


def test_check_invalid_dimensions():
    """ Ensure an error is raised on 1D input arrays. """
    XA = np.arange(45)
    XB = np.resize(np.arange(32), (4, 8))
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)
    XA = np.resize(np.arange(45), (5, 9))
    XB = np.arange(32)
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)


def test_check_sparse_arrays():
    """ Ensures that checks return valid sparse matrices. """
    rng = np.random.RandomState(0)
    XA = rng.random_sample((5, 4))
    XA_sparse = csr_matrix(XA)
    XB = rng.random_sample((5, 4))
    XB_sparse = csr_matrix(XB)
    XA_checked, XB_checked = check_pairwise_arrays(XA_sparse, XB_sparse)
    assert_equal(XA_sparse, XA_checked)
    assert_equal(XB_sparse, XB_checked)


def tuplify(X):
    """ Turns a numpy matrix (any n-dimensional array) into tuples."""
    s = X.shape
    if len(s) > 1:
        # Tuplify each sub-array in the input.
        return tuple(tuplify(row) for row in X)
    else:
        # Single dimension input, just return tuple of contents.
        return tuple(r for r in X)


def test_check_tuple_input():
    """ Ensures that checks return valid tuples. """
    rng = np.random.RandomState(0)
    XA = rng.random_sample((5, 4))
    XA_tuples = tuplify(XA)
    XB = rng.random_sample((5, 4))
    XB_tuples = tuplify(XB)
    XA_checked, XB_checked = check_pairwise_arrays(XA_tuples, XB_tuples)
    assert_equal(XA_tuples, XA_checked)
    assert_equal(XB_tuples, XB_checked)
