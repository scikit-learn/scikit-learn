import numpy as np
from numpy import linalg

from scipy.sparse import dok_matrix, csr_matrix, issparse
from scipy.spatial.distance import cosine, cityblock, minkowski, wminkowski

from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import ignore_warnings

from sklearn.externals.six import iteritems

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import manhattan_distances
from sklearn.metrics.pairwise import linear_kernel
from sklearn.metrics.pairwise import chi2_kernel, additive_chi2_kernel
from sklearn.metrics.pairwise import polynomial_kernel
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.metrics.pairwise import laplacian_kernel
from sklearn.metrics.pairwise import sigmoid_kernel
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics.pairwise import cosine_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import pairwise_distances_argmin_min
from sklearn.metrics.pairwise import pairwise_distances_argmin
from sklearn.metrics.pairwise import pairwise_kernels
from sklearn.metrics.pairwise import PAIRWISE_KERNEL_FUNCTIONS
from sklearn.metrics.pairwise import PAIRWISE_DISTANCE_FUNCTIONS
from sklearn.metrics.pairwise import PAIRWISE_BOOLEAN_FUNCTIONS
from sklearn.metrics.pairwise import PAIRED_DISTANCES
from sklearn.metrics.pairwise import check_pairwise_arrays
from sklearn.metrics.pairwise import check_paired_arrays
from sklearn.metrics.pairwise import paired_distances
from sklearn.metrics.pairwise import paired_euclidean_distances
from sklearn.metrics.pairwise import paired_manhattan_distances
from sklearn.preprocessing import normalize
from sklearn.exceptions import DataConversionWarning


def test_pairwise_distances():
    # Test the pairwise_distance helper function.
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
    # "cityblock" uses scikit-learn metric, cityblock (function) is
    # scipy.spatial.
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
    # Low-level function for manhattan can divide in blocks to avoid
    # using too much memory during the broadcasting
    S3 = manhattan_distances(X, Y, size_threshold=10)
    assert_array_almost_equal(S, S3)
    # Test cosine as a string metric versus cosine callable
    # The string "cosine" uses sklearn.metric,
    # while the function cosine is scipy.spatial
    S = pairwise_distances(X, Y, metric="cosine")
    S2 = pairwise_distances(X, Y, metric=cosine)
    assert_equal(S.shape[0], X.shape[0])
    assert_equal(S.shape[1], Y.shape[0])
    assert_array_almost_equal(S, S2)
    # Test with sparse X and Y,
    # currently only supported for Euclidean, L1 and cosine.
    X_sparse = csr_matrix(X)
    Y_sparse = csr_matrix(Y)
    S = pairwise_distances(X_sparse, Y_sparse, metric="euclidean")
    S2 = euclidean_distances(X_sparse, Y_sparse)
    assert_array_almost_equal(S, S2)
    S = pairwise_distances(X_sparse, Y_sparse, metric="cosine")
    S2 = cosine_distances(X_sparse, Y_sparse)
    assert_array_almost_equal(S, S2)
    S = pairwise_distances(X_sparse, Y_sparse.tocsc(), metric="manhattan")
    S2 = manhattan_distances(X_sparse.tobsr(), Y_sparse.tocoo())
    assert_array_almost_equal(S, S2)
    S2 = manhattan_distances(X, Y)
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

    # Test that a value error is raised if the metric is unknown
    assert_raises(ValueError, pairwise_distances, X, Y, metric="blah")


# ignore conversion to boolean in pairwise_distances
@ignore_warnings(category=DataConversionWarning)
def test_pairwise_boolean_distance():
    # test that we convert to boolean arrays for boolean distances
    rng = np.random.RandomState(0)
    X = rng.randn(5, 4)
    Y = X.copy()
    Y[0, 0] = 1 - Y[0, 0]

    for metric in PAIRWISE_BOOLEAN_FUNCTIONS:
        for Z in [Y, None]:
            res = pairwise_distances(X, Z, metric=metric)
            res[np.isnan(res)] = 0
            assert_true(np.sum(res != 0) == 0)


def test_pairwise_precomputed():
    for func in [pairwise_distances, pairwise_kernels]:
        # Test correct shape
        assert_raises_regexp(ValueError, '.* shape .*',
                             func, np.zeros((5, 3)), metric='precomputed')
        # with two args
        assert_raises_regexp(ValueError, '.* shape .*',
                             func, np.zeros((5, 3)), np.zeros((4, 4)),
                             metric='precomputed')
        # even if shape[1] agrees (although thus second arg is spurious)
        assert_raises_regexp(ValueError, '.* shape .*',
                             func, np.zeros((5, 3)), np.zeros((4, 3)),
                             metric='precomputed')

        # Test not copied (if appropriate dtype)
        S = np.zeros((5, 5))
        S2 = func(S, metric="precomputed")
        assert_true(S is S2)
        # with two args
        S = np.zeros((5, 3))
        S2 = func(S, np.zeros((3, 3)), metric="precomputed")
        assert_true(S is S2)

        # Test always returns float dtype
        S = func(np.array([[1]], dtype='int'), metric='precomputed')
        assert_equal('f', S.dtype.kind)

        # Test converts list to array-like
        S = func([[1.]], metric='precomputed')
        assert_true(isinstance(S, np.ndarray))


def check_pairwise_parallel(func, metric, kwds):
    rng = np.random.RandomState(0)
    for make_data in (np.array, csr_matrix):
        X = make_data(rng.random_sample((5, 4)))
        Y = make_data(rng.random_sample((3, 4)))

        try:
            S = func(X, metric=metric, n_jobs=1, **kwds)
        except (TypeError, ValueError) as exc:
            # Not all metrics support sparse input
            # ValueError may be triggered by bad callable
            if make_data is csr_matrix:
                assert_raises(type(exc), func, X, metric=metric,
                              n_jobs=2, **kwds)
                continue
            else:
                raise
        S2 = func(X, metric=metric, n_jobs=2, **kwds)
        assert_array_almost_equal(S, S2)

        S = func(X, Y, metric=metric, n_jobs=1, **kwds)
        S2 = func(X, Y, metric=metric, n_jobs=2, **kwds)
        assert_array_almost_equal(S, S2)


def test_pairwise_parallel():
    wminkowski_kwds = {'w': np.arange(1, 5).astype('double'), 'p': 1}
    metrics = [(pairwise_distances, 'euclidean', {}),
               (pairwise_distances, wminkowski, wminkowski_kwds),
               (pairwise_distances, 'wminkowski', wminkowski_kwds),
               (pairwise_kernels, 'polynomial', {'degree': 1}),
               (pairwise_kernels, callable_rbf_kernel, {'gamma': .1}),
               ]
    for func, metric, kwds in metrics:
        yield check_pairwise_parallel, func, metric, kwds


def test_pairwise_callable_nonstrict_metric():
    # paired_distances should allow callable metric where metric(x, x) != 0
    # Knowing that the callable is a strict metric would allow the diagonal to
    # be left uncalculated and set to 0.
    assert_equal(pairwise_distances([[1.]], metric=lambda x, y: 5)[0, 0], 5)


def callable_rbf_kernel(x, y, **kwds):
    # Callable version of pairwise.rbf_kernel.
    K = rbf_kernel(np.atleast_2d(x), np.atleast_2d(y), **kwds)
    return K


def test_pairwise_kernels():    # Test the pairwise_kernels helper function.

    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((2, 4))
    # Test with all metrics that should be in PAIRWISE_KERNEL_FUNCTIONS.
    test_metrics = ["rbf", "laplacian", "sigmoid", "polynomial", "linear",
                    "chi2", "additive_chi2"]
    for metric in test_metrics:
        function = PAIRWISE_KERNEL_FUNCTIONS[metric]
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
    kwds = {'gamma': 0.1}
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


def test_paired_distances():
    # Test the pairwise_distance helper function.
    rng = np.random.RandomState(0)
    # Euclidean distance should be equivalent to calling the function.
    X = rng.random_sample((5, 4))
    # Euclidean distance, with Y != X.
    Y = rng.random_sample((5, 4))
    for metric, func in iteritems(PAIRED_DISTANCES):
        S = paired_distances(X, Y, metric=metric)
        S2 = func(X, Y)
        assert_array_almost_equal(S, S2)
        S3 = func(csr_matrix(X), csr_matrix(Y))
        assert_array_almost_equal(S, S3)
        if metric in PAIRWISE_DISTANCE_FUNCTIONS:
            # Check the pairwise_distances implementation
            # gives the same value
            distances = PAIRWISE_DISTANCE_FUNCTIONS[metric](X, Y)
            distances = np.diag(distances)
            assert_array_almost_equal(distances, S)

    # Check the callable implementation
    S = paired_distances(X, Y, metric='manhattan')
    S2 = paired_distances(X, Y, metric=lambda x, y: np.abs(x - y).sum(axis=0))
    assert_array_almost_equal(S, S2)

    # Test that a value error is raised when the lengths of X and Y should not
    # differ
    Y = rng.random_sample((3, 4))
    assert_raises(ValueError, paired_distances, X, Y)


def test_pairwise_distances_argmin_min():
    # Check pairwise minimum distances computation for any metric
    X = [[0], [1]]
    Y = [[-1], [2]]

    Xsp = dok_matrix(X)
    Ysp = csr_matrix(Y, dtype=np.float32)

    # euclidean metric
    D, E = pairwise_distances_argmin_min(X, Y, metric="euclidean")
    D2 = pairwise_distances_argmin(X, Y, metric="euclidean")
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(D2, [0, 1])
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(E, [1., 1.])

    # sparse matrix case
    Dsp, Esp = pairwise_distances_argmin_min(Xsp, Ysp, metric="euclidean")
    assert_array_equal(Dsp, D)
    assert_array_equal(Esp, E)
    # We don't want np.matrix here
    assert_equal(type(Dsp), np.ndarray)
    assert_equal(type(Esp), np.ndarray)

    # Non-euclidean scikit-learn metric
    D, E = pairwise_distances_argmin_min(X, Y, metric="manhattan")
    D2 = pairwise_distances_argmin(X, Y, metric="manhattan")
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(D2, [0, 1])
    assert_array_almost_equal(E, [1., 1.])
    D, E = pairwise_distances_argmin_min(Xsp, Ysp, metric="manhattan")
    D2 = pairwise_distances_argmin(Xsp, Ysp, metric="manhattan")
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(E, [1., 1.])

    # Non-euclidean Scipy distance (callable)
    D, E = pairwise_distances_argmin_min(X, Y, metric=minkowski,
                                         metric_kwargs={"p": 2})
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(E, [1., 1.])

    # Non-euclidean Scipy distance (string)
    D, E = pairwise_distances_argmin_min(X, Y, metric="minkowski",
                                         metric_kwargs={"p": 2})
    assert_array_almost_equal(D, [0, 1])
    assert_array_almost_equal(E, [1., 1.])

    # Compare with naive implementation
    rng = np.random.RandomState(0)
    X = rng.randn(97, 149)
    Y = rng.randn(111, 149)

    dist = pairwise_distances(X, Y, metric="manhattan")
    dist_orig_ind = dist.argmin(axis=0)
    dist_orig_val = dist[dist_orig_ind, range(len(dist_orig_ind))]

    dist_chunked_ind, dist_chunked_val = pairwise_distances_argmin_min(
        X, Y, axis=0, metric="manhattan", batch_size=50)
    np.testing.assert_almost_equal(dist_orig_ind, dist_chunked_ind, decimal=7)
    np.testing.assert_almost_equal(dist_orig_val, dist_chunked_val, decimal=7)


def test_euclidean_distances():
    # Check the pairwise Euclidean distances computation
    X = [[0]]
    Y = [[1], [2]]
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])

    X = csr_matrix(X)
    Y = csr_matrix(Y)
    D = euclidean_distances(X, Y)
    assert_array_almost_equal(D, [[1., 2.]])

    rng = np.random.RandomState(0)
    X = rng.random_sample((10, 4))
    Y = rng.random_sample((20, 4))
    X_norm_sq = (X ** 2).sum(axis=1).reshape(1, -1)
    Y_norm_sq = (Y ** 2).sum(axis=1).reshape(1, -1)

    # check that we still get the right answers with {X,Y}_norm_squared
    D1 = euclidean_distances(X, Y)
    D2 = euclidean_distances(X, Y, X_norm_squared=X_norm_sq)
    D3 = euclidean_distances(X, Y, Y_norm_squared=Y_norm_sq)
    D4 = euclidean_distances(X, Y, X_norm_squared=X_norm_sq,
                             Y_norm_squared=Y_norm_sq)
    assert_array_almost_equal(D2, D1)
    assert_array_almost_equal(D3, D1)
    assert_array_almost_equal(D4, D1)

    # check we get the wrong answer with wrong {X,Y}_norm_squared
    X_norm_sq *= 0.5
    Y_norm_sq *= 0.5
    wrong_D = euclidean_distances(X, Y,
                                  X_norm_squared=np.zeros_like(X_norm_sq),
                                  Y_norm_squared=np.zeros_like(Y_norm_sq))
    assert_greater(np.max(np.abs(wrong_D - D1)), .01)


def test_cosine_distances():
    # Check the pairwise Cosine distances computation
    rng = np.random.RandomState(1337)
    x = np.abs(rng.rand(910))
    XA = np.vstack([x, x])
    D = cosine_distances(XA)
    assert_array_almost_equal(D, [[0., 0.], [0., 0.]])
    # check that all elements are in [0, 2]
    assert_true(np.all(D >= 0.))
    assert_true(np.all(D <= 2.))
    # check that diagonal elements are equal to 0
    assert_array_almost_equal(D[np.diag_indices_from(D)], [0., 0.])

    XB = np.vstack([x, -x])
    D2 = cosine_distances(XB)
    # check that all elements are in [0, 2]
    assert_true(np.all(D2 >= 0.))
    assert_true(np.all(D2 <= 2.))
    # check that diagonal elements are equal to 0 and non diagonal to 2
    assert_array_almost_equal(D2, [[0., 2.], [2., 0.]])

    # check large random matrix
    X = np.abs(rng.rand(1000, 5000))
    D = cosine_distances(X)
    # check that diagonal elements are equal to 0
    assert_array_almost_equal(D[np.diag_indices_from(D)], [0.] * D.shape[0])
    assert_true(np.all(D >= 0.))
    assert_true(np.all(D <= 2.))


# Paired distances

def test_paired_euclidean_distances():
    # Check the paired Euclidean distances computation
    X = [[0], [0]]
    Y = [[1], [2]]
    D = paired_euclidean_distances(X, Y)
    assert_array_almost_equal(D, [1., 2.])


def test_paired_manhattan_distances():
    # Check the paired manhattan distances computation
    X = [[0], [0]]
    Y = [[1], [2]]
    D = paired_manhattan_distances(X, Y)
    assert_array_almost_equal(D, [1., 2.])


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
    Y = [[0, 1], [.9, .1]]
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
    # Valid kernels should be symmetric
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   laplacian_kernel, sigmoid_kernel, cosine_similarity):
        K = kernel(X, X)
        assert_array_almost_equal(K, K.T, 15)


def test_kernel_sparse():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    X_sparse = csr_matrix(X)
    for kernel in (linear_kernel, polynomial_kernel, rbf_kernel,
                   laplacian_kernel, sigmoid_kernel, cosine_similarity):
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


def test_laplacian_kernel():
    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    K = laplacian_kernel(X, X)
    # the diagonal elements of a laplacian kernel are 1
    assert_array_almost_equal(np.diag(K), np.ones(5))

    # off-diagonal elements are < 1 but > 0:
    assert_true(np.all(K > 0))
    assert_true(np.all(K - np.diag(np.diag(K)) < 1))


def test_cosine_similarity_sparse_output():
    # Test if cosine_similarity correctly produces sparse output.

    rng = np.random.RandomState(0)
    X = rng.random_sample((5, 4))
    Y = rng.random_sample((3, 4))
    Xcsr = csr_matrix(X)
    Ycsr = csr_matrix(Y)

    K1 = cosine_similarity(Xcsr, Ycsr, dense_output=False)
    assert_true(issparse(K1))

    K2 = pairwise_kernels(Xcsr, Y=Ycsr, metric="cosine")
    assert_array_almost_equal(K1.todense(), K2)


def test_cosine_similarity():
    # Test the cosine_similarity.

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
    # Ensure that pairwise array check works for dense matrices.
    # Check that if XB is None, XB is returned as reference to XA
    XA = np.resize(np.arange(40), (5, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, None)
    assert_true(XA_checked is XB_checked)
    assert_array_equal(XA, XA_checked)


def test_check_XB_returned():
    # Ensure that if XA and XB are given correctly, they return as equal.
    # Check that if XB is not None, it is returned equal.
    # Note that the second dimension of XB is the same as XA.
    XA = np.resize(np.arange(40), (5, 8))
    XB = np.resize(np.arange(32), (4, 8))
    XA_checked, XB_checked = check_pairwise_arrays(XA, XB)
    assert_array_equal(XA, XA_checked)
    assert_array_equal(XB, XB_checked)

    XB = np.resize(np.arange(40), (5, 8))
    XA_checked, XB_checked = check_paired_arrays(XA, XB)
    assert_array_equal(XA, XA_checked)
    assert_array_equal(XB, XB_checked)


def test_check_different_dimensions():
    # Ensure an error is raised if the dimensions are different.
    XA = np.resize(np.arange(45), (5, 9))
    XB = np.resize(np.arange(32), (4, 8))
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)

    XB = np.resize(np.arange(4 * 9), (4, 9))
    assert_raises(ValueError, check_paired_arrays, XA, XB)


def test_check_invalid_dimensions():
    # Ensure an error is raised on 1D input arrays.
    # The modified tests are not 1D. In the old test, the array was internally
    # converted to 2D anyways
    XA = np.arange(45).reshape(9, 5)
    XB = np.arange(32).reshape(4, 8)
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)
    XA = np.arange(45).reshape(9, 5)
    XB = np.arange(32).reshape(4, 8)
    assert_raises(ValueError, check_pairwise_arrays, XA, XB)


def test_check_sparse_arrays():
    # Ensures that checks return valid sparse matrices.
    rng = np.random.RandomState(0)
    XA = rng.random_sample((5, 4))
    XA_sparse = csr_matrix(XA)
    XB = rng.random_sample((5, 4))
    XB_sparse = csr_matrix(XB)
    XA_checked, XB_checked = check_pairwise_arrays(XA_sparse, XB_sparse)
    # compare their difference because testing csr matrices for
    # equality with '==' does not work as expected.
    assert_true(issparse(XA_checked))
    assert_equal(abs(XA_sparse - XA_checked).sum(), 0)
    assert_true(issparse(XB_checked))
    assert_equal(abs(XB_sparse - XB_checked).sum(), 0)

    XA_checked, XA_2_checked = check_pairwise_arrays(XA_sparse, XA_sparse)
    assert_true(issparse(XA_checked))
    assert_equal(abs(XA_sparse - XA_checked).sum(), 0)
    assert_true(issparse(XA_2_checked))
    assert_equal(abs(XA_2_checked - XA_checked).sum(), 0)


def tuplify(X):
    # Turns a numpy matrix (any n-dimensional array) into tuples.
    s = X.shape
    if len(s) > 1:
        # Tuplify each sub-array in the input.
        return tuple(tuplify(row) for row in X)
    else:
        # Single dimension input, just return tuple of contents.
        return tuple(r for r in X)


def test_check_tuple_input():
    # Ensures that checks return valid tuples.
    rng = np.random.RandomState(0)
    XA = rng.random_sample((5, 4))
    XA_tuples = tuplify(XA)
    XB = rng.random_sample((5, 4))
    XB_tuples = tuplify(XB)
    XA_checked, XB_checked = check_pairwise_arrays(XA_tuples, XB_tuples)
    assert_array_equal(XA_tuples, XA_checked)
    assert_array_equal(XB_tuples, XB_checked)


def test_check_preserve_type():
    # Ensures that type float32 is preserved.
    XA = np.resize(np.arange(40), (5, 8)).astype(np.float32)
    XB = np.resize(np.arange(40), (5, 8)).astype(np.float32)

    XA_checked, XB_checked = check_pairwise_arrays(XA, None)
    assert_equal(XA_checked.dtype, np.float32)

    # both float32
    XA_checked, XB_checked = check_pairwise_arrays(XA, XB)
    assert_equal(XA_checked.dtype, np.float32)
    assert_equal(XB_checked.dtype, np.float32)

    # mismatched A
    XA_checked, XB_checked = check_pairwise_arrays(XA.astype(np.float),
                                                   XB)
    assert_equal(XA_checked.dtype, np.float)
    assert_equal(XB_checked.dtype, np.float)

    # mismatched B
    XA_checked, XB_checked = check_pairwise_arrays(XA,
                                                   XB.astype(np.float))
    assert_equal(XA_checked.dtype, np.float)
    assert_equal(XB_checked.dtype, np.float)
