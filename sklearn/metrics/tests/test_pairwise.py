import numpy as np
from numpy import linalg

from scipy.sparse import csr_matrix
from scipy.spatial.distance import cosine, cityblock, minkowski

from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_true

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import linear_kernel
from sklearn.metrics.pairwise import chi2_kernel, additive_chi2_kernel
from sklearn.metrics.pairwise import polynomial_kernel
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.metrics.pairwise import sigmoid_kernel
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.metrics.pairwise import pairwise_kernel_functions
from sklearn.metrics.pairwise import check_pairwise_arrays
from sklearn.metrics.pairwise import _parallel_pairwise
from sklearn.preprocessing import normalize


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
    # manhattan does not support sparse matrices atm.
    assert_raises(ValueError, pairwise_distances, csr_matrix(X),
                  metric="manhattan")
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
    # same with Y = None
    kwds = {"p": 2.0}
    S = pairwise_distances(X, metric="minkowski", **kwds)
    S2 = pairwise_distances(X, metric=minkowski, **kwds)
    assert_array_almost_equal(S, S2)
    # Test that scipy distance metrics throw an error if sparse matrix given
    assert_raises(TypeError, pairwise_distances, X_sparse, metric="minkowski")
    assert_raises(TypeError, pairwise_distances, X, Y_sparse,
                  metric="minkowski")


def test_pairwise_parallel():
    rng = np.random.RandomState(0)
    for func in (np.array, csr_matrix):
        X = func(rng.random_sample((5, 4)))
        Y = func(rng.random_sample((3, 4)))

        S = euclidean_distances(X)
        S2 = _parallel_pairwise(X, None, euclidean_distances, n_jobs=-1)
        assert_array_almost_equal(S, S2)

        S = euclidean_distances(X, Y)
        S2 = _parallel_pairwise(X, Y, euclidean_distances, n_jobs=-1)
        assert_array_almost_equal(S, S2)


def test_pairwise_kernels():
    """ Test the pairwise_kernels helper function. """

    def callable_rbf_kernel(x, y, **kwds):
        """ Callable version of pairwise.rbf_kernel. """
        K = rbf_kernel(np.atleast_2d(x), np.atleast_2d(y), **kwds)
        return K

    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((2, 4))
    # Test with all metrics that should be in pairwise_kernel_functions.
    test_metrics = ["rbf", "sigmoid", "polynomial", "linear", "chi2",
                    "additive_chi2"]
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
        if metric in ["chi2", "additive_chi2"]:
            # these don't support sparse matrices yet
            assert_raises(ValueError, pairwise_kernels,
                          X_sparse, Y=Y_sparse, metric=metric)
            continue
        K1 = pairwise_kernels(X_sparse, Y=Y_sparse, metric=metric)
        assert_array_almost_equal(K1, K2)
    # Test with a callable function, with given keywords.
    metric = callable_rbf_kernel
    kwds = {}
    kwds['gamma'] = 0.1
    K1 = pairwise_kernels(X, Y=Y, metric=metric, **kwds)
    K2 = rbf_kernel(X, Y=Y, **kwds)
    assert_array_almost_equal(K1, K2)

    # callable function, X=Y
    K1 = pairwise_kernels(X, Y=X, metric=metric, **kwds)
    K2 = rbf_kernel(X, Y=X, **kwds)
    assert_array_almost_equal(K1, K2)


def test_pairwise_kernels_filter_param():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((2, 4))
    K = rbf_kernel(X, Y, gamma=0.1)
    params = {"gamma": 0.1, "blabla": ":)"}
    K2 = pairwise_kernels(X, Y, metric="rbf", filter_params=True, **params)
    assert_array_almost_equal(K, K2)

    assert_raises(TypeError, pairwise_kernels, X, Y, "rbf", **params)


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


def test_chi_square_kernel():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((10, 4))
    K_add = additive_chi2_kernel(X, Y)
    gamma = 0.1
    K = chi2_kernel(X, Y, gamma=gamma)
    assert_equal(K.dtype, np.float)
    for i, x in enumerate(X):
        for j, y in enumerate(Y):
            chi2 = -np.sum((x - y) ** 2 / (x + y))
            chi2_exp = np.exp(gamma * chi2)
            assert_almost_equal(K_add[i, j], chi2)
            assert_almost_equal(K[i, j], chi2_exp)

    # check diagonal is ones for data with itself
    K = chi2_kernel(Y)
    assert_array_equal(np.diag(K), 1)
    # check off-diagonal is < 1 but > 0:
    assert_true(np.all(K > 0))
    assert_true(np.all(K - np.diag(np.diag(K)) < 1))
    # check that float32 is preserved
    X = rng.random_sample((5, 4)).astype(np.float32)
    Y = rng.random_sample((10, 4)).astype(np.float32)
    K = chi2_kernel(X, Y)
    assert_equal(K.dtype, np.float32)

    # check integer type gets converted,
    # check that zeros are handled
    X = rng.random_sample((10, 4)).astype(np.int32)
    K = chi2_kernel(X, X)
    assert_true(np.isfinite(K).all())
    assert_equal(K.dtype, np.float)

    # check that kernel of similar things is greater than dissimilar ones
    X = [[.3, .7], [1., 0]]
    Y = [[0,   1], [.9, .1]]
    K = chi2_kernel(X, Y)
    assert_greater(K[0, 0], K[0, 1])
    assert_greater(K[1, 1], K[1, 0])

    # test negative input
    assert_raises(ValueError, chi2_kernel, [[0, -1]])
    assert_raises(ValueError, chi2_kernel, [[0, -1]], [[-1, -1]])
    assert_raises(ValueError, chi2_kernel, [[0, 1]], [[-1, -1]])

    # different n_features in X and Y
    assert_raises(ValueError, chi2_kernel, [[0, 1]], [[.2, .2, .6]])

    # sparse matrices
    assert_raises(ValueError, chi2_kernel, csr_matrix(X), csr_matrix(Y))
    assert_raises(ValueError, additive_chi2_kernel,
                  csr_matrix(X), csr_matrix(Y))


def test_kernel_symmetry():
    """ Valid kernels should be symmetric"""
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   sigmoid_kernel, cosine_similarity):
        K = kernel(X, X)
        assert_array_almost_equal(K, K.T, 15)


def test_kernel_sparse():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    X_sparse = csr_matrix(X)
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   sigmoid_kernel, cosine_similarity):
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


def test_cosine_similarity():
    """ Test the cosine_similarity. """

    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((3, 4))
    Xcsr = csr_matrix(X)
    Ycsr = csr_matrix(Y)

    for X_, Y_ in ((X, None), (X, Y),
                   (Xcsr, None), (Xcsr, Ycsr)):
        # Test that the cosine is kernel is equal to a linear kernel when data
        # has been previously normalized by L2-norm.
        K1 = pairwise_kernels(X_, Y=Y_, metric="cosine")
        X_ = normalize(X_)
        if Y_ is not None:
            Y_ = normalize(Y_)
        K2 = pairwise_kernels(X_, Y=Y_, metric="linear")
        assert_array_almost_equal(K1, K2)


def test_check_dense_matrices():
    """ Ensure that pairwise array check works for dense matrices."""
    # Check that if XB is None, XB is returned as reference to XA
    XA = np.resize(np.arange(40), (5, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, None)
    assert_true(XA_checked is XB_checked)
    assert_array_equal(XA, XA_checked)


def test_check_XB_returned():
    """ Ensure that if XA and XB are given correctly, they return as equal."""
    # Check that if XB is not None, it is returned equal.
    # Note that the second dimension of XB is the same as XA.
    XA = np.resize(np.arange(40), (5, 8))
    XB = np.resize(np.arange(32), (4, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, XB)
    assert_array_equal(XA, XA_checked)
    assert_array_equal(XB, XB_checked)


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
    assert_array_equal(XA_tuples, XA_checked)
    assert_array_equal(XB_tuples, XB_checked)
