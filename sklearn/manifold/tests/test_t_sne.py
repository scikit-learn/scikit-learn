import sys
from io import StringIO
import numpy as np
import scipy.sparse as sp

import pytest

from sklearn.neighbors import BallTree
from sklearn.neighbors import NearestNeighbors
from sklearn.utils.testing import assert_less_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import skip_if_32bit
from sklearn.utils import check_random_state
from sklearn.manifold.t_sne import _joint_probabilities
from sklearn.manifold.t_sne import _joint_probabilities_nn
from sklearn.manifold.t_sne import _kl_divergence
from sklearn.manifold.t_sne import _kl_divergence_bh
from sklearn.manifold.t_sne import _gradient_descent
from sklearn.manifold.t_sne import trustworthiness
from sklearn.manifold.t_sne import TSNE
from sklearn.manifold import _barnes_hut_tsne
from sklearn.manifold._utils import _binary_search_perplexity
from sklearn.datasets import make_blobs
from scipy.optimize import check_grad
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import manhattan_distances
from sklearn.metrics.pairwise import cosine_distances


x = np.linspace(0, 1, 10)
xx, yy = np.meshgrid(x, x)
X_2d_grid = np.hstack([
    xx.ravel().reshape(-1, 1),
    yy.ravel().reshape(-1, 1),
])


def test_gradient_descent_stops():
    # Test stopping conditions of gradient descent.
    class ObjectiveSmallGradient:
        def __init__(self):
            self.it = -1

        def __call__(self, _, compute_error=True):
            self.it += 1
            return (10 - self.it) / 10.0, np.array([1e-5])

    def flat_function(_, compute_error=True):
        return 0.0, np.ones(1)

    # Gradient norm
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            ObjectiveSmallGradient(), np.zeros(1), 0, n_iter=100,
            n_iter_without_progress=100, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=1e-5, verbose=2)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout
    assert_equal(error, 1.0)
    assert_equal(it, 0)
    assert("gradient norm" in out)

    # Maximum number of iterations without improvement
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            flat_function, np.zeros(1), 0, n_iter=100,
            n_iter_without_progress=10, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=0.0, verbose=2)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout
    assert_equal(error, 0.0)
    assert_equal(it, 11)
    assert("did not make any progress" in out)

    # Maximum number of iterations
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            ObjectiveSmallGradient(), np.zeros(1), 0, n_iter=11,
            n_iter_without_progress=100, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=0.0, verbose=2)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout
    assert_equal(error, 0.0)
    assert_equal(it, 10)
    assert("Iteration 10" in out)


def test_binary_search():
    # Test if the binary search finds Gaussians with desired perplexity.
    random_state = check_random_state(0)
    distances = random_state.randn(50, 2).astype(np.float32)
    # Distances shouldn't be negative
    distances = np.abs(distances.dot(distances.T))
    np.fill_diagonal(distances, 0.0)
    desired_perplexity = 25.0
    P = _binary_search_perplexity(distances, None, desired_perplexity,
                                  verbose=0)
    P = np.maximum(P, np.finfo(np.double).eps)
    mean_perplexity = np.mean([np.exp(-np.sum(P[i] * np.log(P[i])))
                               for i in range(P.shape[0])])
    assert_almost_equal(mean_perplexity, desired_perplexity, decimal=3)


def test_binary_search_neighbors():
    # Binary perplexity search approximation.
    # Should be approximately equal to the slow method when we use
    # all points as neighbors.
    n_samples = 500
    desired_perplexity = 25.0
    random_state = check_random_state(0)
    distances = random_state.randn(n_samples, 2).astype(np.float32)
    # Distances shouldn't be negative
    distances = np.abs(distances.dot(distances.T))
    np.fill_diagonal(distances, 0.0)
    P1 = _binary_search_perplexity(distances, None, desired_perplexity,
                                   verbose=0)

    # Test that when we use all the neighbors the results are identical
    k = n_samples
    neighbors_nn = np.argsort(distances, axis=1)[:, 1:k].astype(np.int64,
                                                                copy=False)
    distances_nn = np.array([distances[k, neighbors_nn[k]]
                            for k in range(n_samples)])
    P2 = _binary_search_perplexity(distances_nn, neighbors_nn,
                                   desired_perplexity, verbose=0)
    P_nn = np.array([P1[k, neighbors_nn[k]] for k in range(n_samples)])
    assert_array_almost_equal(P_nn, P2, decimal=4)

    # Test that the highest P_ij are the same when few neighbors are used
    for k in np.linspace(80, n_samples, 5):
        k = int(k)
        topn = k * 10  # check the top 10 *k entries out of k * k entries
        neighbors_nn = np.argsort(distances, axis=1)[:, :k].astype(np.int64,
                                                                   copy=False)
        distances_nn = np.array([distances[k, neighbors_nn[k]]
                                for k in range(n_samples)])
        P2k = _binary_search_perplexity(distances_nn, neighbors_nn,
                                        desired_perplexity, verbose=0)
        idx = np.argsort(P1.ravel())[::-1]
        P1top = P1.ravel()[idx][:topn]
        idx = np.argsort(P2k.ravel())[::-1]
        P2top = P2k.ravel()[idx][:topn]
        assert_array_almost_equal(P1top, P2top, decimal=2)


def test_binary_perplexity_stability():
    # Binary perplexity search should be stable.
    # The binary_search_perplexity had a bug wherein the P array
    # was uninitialized, leading to sporadically failing tests.
    k = 10
    n_samples = 100
    random_state = check_random_state(0)
    distances = random_state.randn(n_samples, 2).astype(np.float32)
    # Distances shouldn't be negative
    distances = np.abs(distances.dot(distances.T))
    np.fill_diagonal(distances, 0.0)
    last_P = None
    neighbors_nn = np.argsort(distances, axis=1)[:, :k].astype(np.int64,
                                                               copy=False)
    for _ in range(100):
        P = _binary_search_perplexity(distances.copy(), neighbors_nn.copy(),
                                      3, verbose=0)
        P1 = _joint_probabilities_nn(distances, neighbors_nn, 3, verbose=0)
        # Convert the sparse matrix to a dense one for testing
        P1 = P1.toarray()
        if last_P is None:
            last_P = P
            last_P1 = P1
        else:
            assert_array_almost_equal(P, last_P, decimal=4)
            assert_array_almost_equal(P1, last_P1, decimal=4)


def test_gradient():
    # Test gradient of Kullback-Leibler divergence.
    random_state = check_random_state(0)

    n_samples = 50
    n_features = 2
    n_components = 2
    alpha = 1.0

    distances = random_state.randn(n_samples, n_features).astype(np.float32)
    distances = np.abs(distances.dot(distances.T))
    np.fill_diagonal(distances, 0.0)
    X_embedded = random_state.randn(n_samples, n_components).astype(np.float32)

    P = _joint_probabilities(distances, desired_perplexity=25.0,
                             verbose=0)

    def fun(params):
        return _kl_divergence(params, P, alpha, n_samples, n_components)[0]

    def grad(params):
        return _kl_divergence(params, P, alpha, n_samples, n_components)[1]

    assert_almost_equal(check_grad(fun, grad, X_embedded.ravel()), 0.0,
                        decimal=5)


def test_trustworthiness():
    # Test trustworthiness score.
    random_state = check_random_state(0)

    # Affine transformation
    X = random_state.randn(100, 2)
    assert_equal(trustworthiness(X, 5.0 + X / 10.0), 1.0)

    # Randomly shuffled
    X = np.arange(100).reshape(-1, 1)
    X_embedded = X.copy()
    random_state.shuffle(X_embedded)
    assert_less(trustworthiness(X, X_embedded), 0.6)

    # Completely different
    X = np.arange(5).reshape(-1, 1)
    X_embedded = np.array([[0], [2], [4], [1], [3]])
    assert_almost_equal(trustworthiness(X, X_embedded, n_neighbors=1), 0.2)


def test_preserve_trustworthiness_approximately():
    # Nearest neighbors should be preserved approximately.
    random_state = check_random_state(0)
    n_components = 2
    methods = ['exact', 'barnes_hut']
    X = random_state.randn(50, n_components).astype(np.float32)
    for init in ('random', 'pca'):
        for method in methods:
            tsne = TSNE(n_components=n_components, init=init, random_state=0,
                        method=method)
            X_embedded = tsne.fit_transform(X)
            t = trustworthiness(X, X_embedded, n_neighbors=1)
            assert_greater(t, 0.85, msg='Trustworthiness={:0.3f} < 0.85 '
                                        'for method={} and '
                                        'init={}'.format(t, method, init))


def test_optimization_minimizes_kl_divergence():
    """t-SNE should give a lower KL divergence with more iterations."""
    random_state = check_random_state(0)
    X, _ = make_blobs(n_features=3, random_state=random_state)
    kl_divergences = []
    for n_iter in [250, 300, 350]:
        tsne = TSNE(n_components=2, perplexity=10, learning_rate=100.0,
                    n_iter=n_iter, random_state=0)
        tsne.fit_transform(X)
        kl_divergences.append(tsne.kl_divergence_)
    assert_less_equal(kl_divergences[1], kl_divergences[0])
    assert_less_equal(kl_divergences[2], kl_divergences[1])


def test_fit_csr_matrix():
    # X can be a sparse matrix.
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    X[(np.random.randint(0, 100, 50), np.random.randint(0, 2, 50))] = 0.0
    X_csr = sp.csr_matrix(X)
    tsne = TSNE(n_components=2, perplexity=10, learning_rate=100.0,
                random_state=0, method='exact')
    X_embedded = tsne.fit_transform(X_csr)
    assert_almost_equal(trustworthiness(X_csr, X_embedded, n_neighbors=1), 1.0,
                        decimal=1)


def test_preserve_trustworthiness_approximately_with_precomputed_distances():
    # Nearest neighbors should be preserved approximately.
    random_state = check_random_state(0)
    for i in range(3):
        X = random_state.randn(100, 2)
        D = squareform(pdist(X), "sqeuclidean")
        tsne = TSNE(n_components=2, perplexity=2, learning_rate=100.0,
                    early_exaggeration=2.0, metric="precomputed",
                    random_state=i, verbose=0)
        X_embedded = tsne.fit_transform(D)
        t = trustworthiness(D, X_embedded, n_neighbors=1, metric="precomputed")
        assert t > .95


def test_trustworthiness_precomputed_deprecation():
    # FIXME: Remove this test in v0.23

    # Use of the flag `precomputed` in trustworthiness parameters has been
    # deprecated, but will still work until v0.23.
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    assert_equal(assert_warns(DeprecationWarning, trustworthiness,
                              pairwise_distances(X), X, precomputed=True), 1.)
    assert_equal(assert_warns(DeprecationWarning, trustworthiness,
                              pairwise_distances(X), X, metric='precomputed',
                              precomputed=True), 1.)
    assert_raises(ValueError, assert_warns, DeprecationWarning,
                  trustworthiness, X, X, metric='euclidean', precomputed=True)
    assert_equal(assert_warns(DeprecationWarning, trustworthiness,
                              pairwise_distances(X), X, metric='euclidean',
                              precomputed=True), 1.)


def test_trustworthiness_not_euclidean_metric():
    # Test trustworthiness with a metric different from 'euclidean' and
    # 'precomputed'
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    assert_equal(trustworthiness(X, X, metric='cosine'),
                 trustworthiness(pairwise_distances(X, metric='cosine'), X,
                                 metric='precomputed'))


def test_early_exaggeration_too_small():
    # Early exaggeration factor must be >= 1.
    tsne = TSNE(early_exaggeration=0.99)
    assert_raises_regexp(ValueError, "early_exaggeration .*",
                         tsne.fit_transform, np.array([[0.0], [0.0]]))


def test_too_few_iterations():
    # Number of gradient descent iterations must be at least 200.
    tsne = TSNE(n_iter=199)
    assert_raises_regexp(ValueError, "n_iter .*", tsne.fit_transform,
                         np.array([[0.0], [0.0]]))


def test_non_square_precomputed_distances():
    # Precomputed distance matrices must be square matrices.
    tsne = TSNE(metric="precomputed")
    assert_raises_regexp(ValueError, ".* square distance matrix",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_non_positive_precomputed_distances():
    # Precomputed distance matrices must be positive.
    bad_dist = np.array([[0., -1.], [1., 0.]])
    for method in ['barnes_hut', 'exact']:
        tsne = TSNE(metric="precomputed", method=method)
        assert_raises_regexp(ValueError, "All distances .*precomputed.*",
                             tsne.fit_transform, bad_dist)


def test_non_positive_computed_distances():
    # Computed distance matrices must be positive.
    def metric(x, y):
        return -1

    tsne = TSNE(metric=metric, method='exact')
    X = np.array([[0.0, 0.0], [1.0, 1.0]])
    assert_raises_regexp(ValueError, "All distances .*metric given.*",
                         tsne.fit_transform, X)


def test_init_not_available():
    # 'init' must be 'pca', 'random', or numpy array.
    tsne = TSNE(init="not available")
    m = "'init' must be 'pca', 'random', or a numpy array"
    assert_raises_regexp(ValueError, m, tsne.fit_transform,
                         np.array([[0.0], [1.0]]))


def test_init_ndarray():
    # Initialize TSNE with ndarray and test fit
    tsne = TSNE(init=np.zeros((100, 2)))
    X_embedded = tsne.fit_transform(np.ones((100, 5)))
    assert_array_equal(np.zeros((100, 2)), X_embedded)


def test_init_ndarray_precomputed():
    # Initialize TSNE with ndarray and metric 'precomputed'
    # Make sure no FutureWarning is thrown from _fit
    tsne = TSNE(init=np.zeros((100, 2)), metric="precomputed")
    tsne.fit(np.zeros((100, 100)))


def test_distance_not_available():
    # 'metric' must be valid.
    tsne = TSNE(metric="not available", method='exact')
    assert_raises_regexp(ValueError, "Unknown metric not available.*",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))

    tsne = TSNE(metric="not available", method='barnes_hut')
    assert_raises_regexp(ValueError, "Metric 'not available' not valid.*",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_method_not_available():
    # 'nethod' must be 'barnes_hut' or 'exact'
    tsne = TSNE(method='not available')
    assert_raises_regexp(ValueError, "'method' must be 'barnes_hut' or ",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_angle_out_of_range_checks():
    # check the angle parameter range
    for angle in [-1, -1e-6, 1 + 1e-6, 2]:
        tsne = TSNE(angle=angle)
        assert_raises_regexp(ValueError, "'angle' must be between 0.0 - 1.0",
                             tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_pca_initialization_not_compatible_with_precomputed_kernel():
    # Precomputed distance matrices must be square matrices.
    tsne = TSNE(metric="precomputed", init="pca")
    assert_raises_regexp(ValueError, "The parameter init=\"pca\" cannot be "
                         "used with metric=\"precomputed\".",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_n_components_range():
    # barnes_hut method should only be used with n_components <= 3
    tsne = TSNE(n_components=4, method="barnes_hut")
    assert_raises_regexp(ValueError, "'n_components' should be .*",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_early_exaggeration_used():
    # check that the ``early_exaggeration`` parameter has an effect
    random_state = check_random_state(0)
    n_components = 2
    methods = ['exact', 'barnes_hut']
    X = random_state.randn(25, n_components).astype(np.float32)
    for method in methods:
        tsne = TSNE(n_components=n_components, perplexity=1,
                    learning_rate=100.0, init="pca", random_state=0,
                    method=method, early_exaggeration=1.0)
        X_embedded1 = tsne.fit_transform(X)
        tsne = TSNE(n_components=n_components, perplexity=1,
                    learning_rate=100.0, init="pca", random_state=0,
                    method=method, early_exaggeration=10.0)
        X_embedded2 = tsne.fit_transform(X)

        assert not np.allclose(X_embedded1, X_embedded2)


def test_n_iter_used():
    # check that the ``n_iter`` parameter has an effect
    random_state = check_random_state(0)
    n_components = 2
    methods = ['exact', 'barnes_hut']
    X = random_state.randn(25, n_components).astype(np.float32)
    for method in methods:
        for n_iter in [251, 500]:
            tsne = TSNE(n_components=n_components, perplexity=1,
                        learning_rate=0.5, init="random", random_state=0,
                        method=method, early_exaggeration=1.0, n_iter=n_iter)
            tsne.fit_transform(X)

            assert tsne.n_iter_ == n_iter - 1


def test_answer_gradient_two_points():
    # Test the tree with only a single set of children.
    #
    # These tests & answers have been checked against the reference
    # implementation by LvdM.
    pos_input = np.array([[1.0, 0.0], [0.0, 1.0]])
    pos_output = np.array([[-4.961291e-05, -1.072243e-04],
                           [9.259460e-05, 2.702024e-04]])
    neighbors = np.array([[1],
                          [0]])
    grad_output = np.array([[-2.37012478e-05, -6.29044398e-05],
                            [2.37012478e-05, 6.29044398e-05]])
    _run_answer_test(pos_input, pos_output, neighbors, grad_output)


def test_answer_gradient_four_points():
    # Four points tests the tree with multiple levels of children.
    #
    # These tests & answers have been checked against the reference
    # implementation by LvdM.
    pos_input = np.array([[1.0, 0.0], [0.0, 1.0],
                          [5.0, 2.0], [7.3, 2.2]])
    pos_output = np.array([[6.080564e-05, -7.120823e-05],
                           [-1.718945e-04, -4.000536e-05],
                           [-2.271720e-04, 8.663310e-05],
                           [-1.032577e-04, -3.582033e-05]])
    neighbors = np.array([[1, 2, 3],
                          [0, 2, 3],
                          [1, 0, 3],
                          [1, 2, 0]])
    grad_output = np.array([[5.81128448e-05, -7.78033454e-06],
                            [-5.81526851e-05, 7.80976444e-06],
                            [4.24275173e-08, -3.69569698e-08],
                            [-2.58720939e-09, 7.52706374e-09]])
    _run_answer_test(pos_input, pos_output, neighbors, grad_output)


def test_skip_num_points_gradient():
    # Test the kwargs option skip_num_points.
    #
    # Skip num points should make it such that the Barnes_hut gradient
    # is not calculated for indices below skip_num_point.
    # Aside from skip_num_points=2 and the first two gradient rows
    # being set to zero, these data points are the same as in
    # test_answer_gradient_four_points()
    pos_input = np.array([[1.0, 0.0], [0.0, 1.0],
                          [5.0, 2.0], [7.3, 2.2]])
    pos_output = np.array([[6.080564e-05, -7.120823e-05],
                           [-1.718945e-04, -4.000536e-05],
                           [-2.271720e-04, 8.663310e-05],
                           [-1.032577e-04, -3.582033e-05]])
    neighbors = np.array([[1, 2, 3],
                          [0, 2, 3],
                          [1, 0, 3],
                          [1, 2, 0]])
    grad_output = np.array([[0.0, 0.0],
                            [0.0, 0.0],
                            [4.24275173e-08, -3.69569698e-08],
                            [-2.58720939e-09, 7.52706374e-09]])
    _run_answer_test(pos_input, pos_output, neighbors, grad_output,
                     False, 0.1, 2)


def _run_answer_test(pos_input, pos_output, neighbors, grad_output,
                     verbose=False, perplexity=0.1, skip_num_points=0):
    distances = pairwise_distances(pos_input).astype(np.float32)
    args = distances, perplexity, verbose
    pos_output = pos_output.astype(np.float32)
    neighbors = neighbors.astype(np.int64, copy=False)
    pij_input = _joint_probabilities(*args)
    pij_input = squareform(pij_input).astype(np.float32)
    grad_bh = np.zeros(pos_output.shape, dtype=np.float32)

    from scipy.sparse import csr_matrix
    P = csr_matrix(pij_input)

    neighbors = P.indices.astype(np.int64)
    indptr = P.indptr.astype(np.int64)

    _barnes_hut_tsne.gradient(P.data, pos_output, neighbors, indptr,
                              grad_bh, 0.5, 2, 1, skip_num_points=0)
    assert_array_almost_equal(grad_bh, grad_output, decimal=4)


def test_verbose():
    # Verbose options write to stdout.
    random_state = check_random_state(0)
    tsne = TSNE(verbose=2)
    X = random_state.randn(5, 2)

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        tsne.fit_transform(X)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    assert("[t-SNE]" in out)
    assert("nearest neighbors..." in out)
    assert("Computed conditional probabilities" in out)
    assert("Mean sigma" in out)
    assert("early exaggeration" in out)


def test_chebyshev_metric():
    # t-SNE should allow metrics that cannot be squared (issue #3526).
    random_state = check_random_state(0)
    tsne = TSNE(metric="chebyshev")
    X = random_state.randn(5, 2)
    tsne.fit_transform(X)


def test_reduction_to_one_component():
    # t-SNE should allow reduction to one component (issue #4154).
    random_state = check_random_state(0)
    tsne = TSNE(n_components=1)
    X = random_state.randn(5, 2)
    X_embedded = tsne.fit(X).embedding_
    assert(np.all(np.isfinite(X_embedded)))


def test_no_sparse_on_barnes_hut():
    # No sparse matrices allowed on Barnes-Hut.
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    X[(np.random.randint(0, 100, 50), np.random.randint(0, 2, 50))] = 0.0
    X_csr = sp.csr_matrix(X)
    tsne = TSNE(n_iter=199, method='barnes_hut')
    assert_raises_regexp(TypeError, "A sparse matrix was.*",
                         tsne.fit_transform, X_csr)


@pytest.mark.parametrize('method', ['barnes_hut', 'exact'])
@pytest.mark.parametrize('dt', [np.float32, np.float64])
def test_64bit(method, dt):
    # Ensure 64bit arrays are handled correctly.
    random_state = check_random_state(0)

    X = random_state.randn(50, 2).astype(dt, copy=False)
    tsne = TSNE(n_components=2, perplexity=2, learning_rate=100.0,
                random_state=0, method=method, verbose=0)
    X_embedded = tsne.fit_transform(X)
    effective_type = X_embedded.dtype

    # tsne cython code is only single precision, so the output will
    # always be single precision, irrespectively of the input dtype
    assert effective_type == np.float32


@pytest.mark.parametrize('method', ['barnes_hut', 'exact'])
def test_kl_divergence_not_nan(method):
    # Ensure kl_divergence_ is computed at last iteration
    # even though n_iter % n_iter_check != 0, i.e. 1003 % 50 != 0
    random_state = check_random_state(0)

    X = random_state.randn(50, 2)
    tsne = TSNE(n_components=2, perplexity=2, learning_rate=100.0,
                random_state=0, method=method, verbose=0, n_iter=1003)
    tsne.fit_transform(X)

    assert not np.isnan(tsne.kl_divergence_)


def test_barnes_hut_angle():
    # When Barnes-Hut's angle=0 this corresponds to the exact method.
    angle = 0.0
    perplexity = 10
    n_samples = 100
    for n_components in [2, 3]:
        n_features = 5
        degrees_of_freedom = float(n_components - 1.0)

        random_state = check_random_state(0)
        distances = random_state.randn(n_samples, n_features)
        distances = distances.astype(np.float32)
        distances = abs(distances.dot(distances.T))
        np.fill_diagonal(distances, 0.0)
        params = random_state.randn(n_samples, n_components)
        P = _joint_probabilities(distances, perplexity, verbose=0)
        kl_exact, grad_exact = _kl_divergence(params, P, degrees_of_freedom,
                                              n_samples, n_components)

        k = n_samples - 1
        bt = BallTree(distances)
        distances_nn, neighbors_nn = bt.query(distances, k=k + 1)
        neighbors_nn = neighbors_nn[:, 1:]
        distances_nn = np.array([distances[i, neighbors_nn[i]]
                                 for i in range(n_samples)])
        assert np.all(distances[0, neighbors_nn[0]] == distances_nn[0]),\
            abs(distances[0, neighbors_nn[0]] - distances_nn[0])
        P_bh = _joint_probabilities_nn(distances_nn, neighbors_nn,
                                       perplexity, verbose=0)
        kl_bh, grad_bh = _kl_divergence_bh(params, P_bh, degrees_of_freedom,
                                           n_samples, n_components,
                                           angle=angle, skip_num_points=0,
                                           verbose=0)

        P = squareform(P)
        P_bh = P_bh.toarray()
        assert_array_almost_equal(P_bh, P, decimal=5)
        assert_almost_equal(kl_exact, kl_bh, decimal=3)


@skip_if_32bit
def test_n_iter_without_progress():
    # Use a dummy negative n_iter_without_progress and check output on stdout
    random_state = check_random_state(0)
    X = random_state.randn(100, 10)
    for method in ["barnes_hut", "exact"]:
        tsne = TSNE(n_iter_without_progress=-1, verbose=2, learning_rate=1e8,
                    random_state=0, method=method, n_iter=351, init="random")
        tsne._N_ITER_CHECK = 1
        tsne._EXPLORATION_N_ITER = 0

        old_stdout = sys.stdout
        sys.stdout = StringIO()
        try:
            tsne.fit_transform(X)
        finally:
            out = sys.stdout.getvalue()
            sys.stdout.close()
            sys.stdout = old_stdout

        # The output needs to contain the value of n_iter_without_progress
        assert_in("did not make any progress during the "
                  "last -1 episodes. Finished.", out)


def test_min_grad_norm():
    # Make sure that the parameter min_grad_norm is used correctly
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    min_grad_norm = 0.002
    tsne = TSNE(min_grad_norm=min_grad_norm, verbose=2,
                random_state=0, method='exact')

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        tsne.fit_transform(X)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    lines_out = out.split('\n')

    # extract the gradient norm from the verbose output
    gradient_norm_values = []
    for line in lines_out:
        # When the computation is Finished just an old gradient norm value
        # is repeated that we do not need to store
        if 'Finished' in line:
            break

        start_grad_norm = line.find('gradient norm')
        if start_grad_norm >= 0:
            line = line[start_grad_norm:]
            line = line.replace('gradient norm = ', '').split(' ')[0]
            gradient_norm_values.append(float(line))

    # Compute how often the gradient norm is smaller than min_grad_norm
    gradient_norm_values = np.array(gradient_norm_values)
    n_smaller_gradient_norms = \
        len(gradient_norm_values[gradient_norm_values <= min_grad_norm])

    # The gradient norm can be smaller than min_grad_norm at most once,
    # because in the moment it becomes smaller the optimization stops
    assert_less_equal(n_smaller_gradient_norms, 1)


def test_accessible_kl_divergence():
    # Ensures that the accessible kl_divergence matches the computed value
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    tsne = TSNE(n_iter_without_progress=2, verbose=2,
                random_state=0, method='exact')

    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        tsne.fit_transform(X)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # The output needs to contain the accessible kl_divergence as the error at
    # the last iteration
    for line in out.split('\n')[::-1]:
        if 'Iteration' in line:
            _, _, error = line.partition('error = ')
            if error:
                error, _, _ = error.partition(',')
                break
    assert_almost_equal(tsne.kl_divergence_, float(error), decimal=5)


def check_uniform_grid(method, seeds=[0, 1, 2], n_iter=1000):
    """Make sure that TSNE can approximately recover a uniform 2D grid

    Due to ties in distances between point in X_2d_grid, this test is platform
    dependent for ``method='barnes_hut'`` due to numerical imprecision.

    Also, t-SNE is not assured to converge to the right solution because bad
    initialization can lead to convergence to bad local minimum (the
    optimization problem is non-convex). To avoid breaking the test too often,
    we re-run t-SNE from the final point when the convergence is not good
    enough.
    """
    for seed in seeds:
        tsne = TSNE(n_components=2, init='random', random_state=seed,
                    perplexity=20, n_iter=n_iter, method=method)
        Y = tsne.fit_transform(X_2d_grid)

        try_name = "{}_{}".format(method, seed)
        try:
            assert_uniform_grid(Y, try_name)
        except AssertionError:
            # If the test fails a first time, re-run with init=Y to see if
            # this was caused by a bad initialization. Note that this will
            # also run an early_exaggeration step.
            try_name += ":rerun"
            tsne.init = Y
            Y = tsne.fit_transform(X_2d_grid)
            assert_uniform_grid(Y, try_name)


def assert_uniform_grid(Y, try_name=None):
    # Ensure that the resulting embedding leads to approximately
    # uniformly spaced points: the distance to the closest neighbors
    # should be non-zero and approximately constant.
    nn = NearestNeighbors(n_neighbors=1).fit(Y)
    dist_to_nn = nn.kneighbors(return_distance=True)[0].ravel()
    assert dist_to_nn.min() > 0.1

    smallest_to_mean = dist_to_nn.min() / np.mean(dist_to_nn)
    largest_to_mean = dist_to_nn.max() / np.mean(dist_to_nn)

    assert_greater(smallest_to_mean, .5, msg=try_name)
    assert_less(largest_to_mean, 2, msg=try_name)


@pytest.mark.parametrize('method', ['barnes_hut', 'exact'])
def test_uniform_grid(method):
    check_uniform_grid(method)


def test_bh_match_exact():
    # check that the ``barnes_hut`` method match the exact one when
    # ``angle = 0`` and ``perplexity > n_samples / 3``
    random_state = check_random_state(0)
    n_features = 10
    X = random_state.randn(30, n_features).astype(np.float32)
    X_embeddeds = {}
    n_iter = {}
    for method in ['exact', 'barnes_hut']:
        tsne = TSNE(n_components=2, method=method, learning_rate=1.0,
                    init="random", random_state=0, n_iter=251,
                    perplexity=30.0, angle=0)
        # Kill the early_exaggeration
        tsne._EXPLORATION_N_ITER = 0
        X_embeddeds[method] = tsne.fit_transform(X)
        n_iter[method] = tsne.n_iter_

    assert n_iter['exact'] == n_iter['barnes_hut']
    assert_array_almost_equal(X_embeddeds['exact'], X_embeddeds['barnes_hut'],
                              decimal=3)


def test_tsne_with_different_distance_metrics():
    """Make sure that TSNE works for different distance metrics"""
    random_state = check_random_state(0)
    n_components_original = 3
    n_components_embedding = 2
    X = random_state.randn(50, n_components_original).astype(np.float32)
    metrics = ['manhattan', 'cosine']
    dist_funcs = [manhattan_distances, cosine_distances]
    for metric, dist_func in zip(metrics, dist_funcs):
        X_transformed_tsne = TSNE(
            metric=metric, n_components=n_components_embedding,
            random_state=0).fit_transform(X)
        X_transformed_tsne_precomputed = TSNE(
            metric='precomputed', n_components=n_components_embedding,
            random_state=0).fit_transform(dist_func(X))
        assert_array_equal(X_transformed_tsne, X_transformed_tsne_precomputed)
