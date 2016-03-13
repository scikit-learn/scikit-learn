import sys
from sklearn.externals.six.moves import cStringIO as StringIO
import numpy as np
import scipy.sparse as sp

from sklearn.neighbors import BallTree
from sklearn.utils.testing import assert_less_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises_regexp
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


def test_gradient_descent_stops():
    # Test stopping conditions of gradient descent.
    class ObjectiveSmallGradient:
        def __init__(self):
            self.it = -1

        def __call__(self, _):
            self.it += 1
            return (10 - self.it) / 10.0, np.array([1e-5])

    def flat_function(_):
        return 0.0, np.ones(1)

    # Gradient norm
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            ObjectiveSmallGradient(), np.zeros(1), 0, n_iter=100,
            n_iter_without_progress=100, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=1e-5, min_error_diff=0.0, verbose=2)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout
    assert_equal(error, 1.0)
    assert_equal(it, 0)
    assert("gradient norm" in out)

    # Error difference
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            ObjectiveSmallGradient(), np.zeros(1), 0, n_iter=100,
            n_iter_without_progress=100, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=0.0, min_error_diff=0.2, verbose=2)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout
    assert_equal(error, 0.9)
    assert_equal(it, 1)
    assert("error difference" in out)

    # Maximum number of iterations without improvement
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    try:
        _, error, it = _gradient_descent(
            flat_function, np.zeros(1), 0, n_iter=100,
            n_iter_without_progress=10, momentum=0.0, learning_rate=0.0,
            min_gain=0.0, min_grad_norm=0.0, min_error_diff=-1.0, verbose=2)
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
            min_gain=0.0, min_grad_norm=0.0, min_error_diff=0.0, verbose=2)
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
    neighbors_nn = np.argsort(distances, axis=1)[:, :k].astype(np.int64)
    P2 = _binary_search_perplexity(distances, neighbors_nn,
                                   desired_perplexity, verbose=0)
    assert_array_almost_equal(P1, P2, decimal=4)

    # Test that the highest P_ij are the same when few neighbors are used
    for k in np.linspace(80, n_samples, 10):
        k = int(k)
        topn = k * 10  # check the top 10 *k entries out of k * k entries
        neighbors_nn = np.argsort(distances, axis=1)[:, :k].astype(np.int64)
        P2k = _binary_search_perplexity(distances, neighbors_nn,
                                        desired_perplexity, verbose=0)
        idx = np.argsort(P1.ravel())[::-1]
        P1top = P1.ravel()[idx][:topn]
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
    neighbors_nn = np.argsort(distances, axis=1)[:, :k].astype(np.int64)
    for _ in range(100):
        P = _binary_search_perplexity(distances.copy(), neighbors_nn.copy(),
                                      3, verbose=0)
        P1 = _joint_probabilities_nn(distances, neighbors_nn, 3, verbose=0)
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
    distances = distances.dot(distances.T)
    np.fill_diagonal(distances, 0.0)
    X_embedded = random_state.randn(n_samples, n_components)

    P = _joint_probabilities(distances, desired_perplexity=25.0,
                             verbose=0)
    fun = lambda params: _kl_divergence(params, P, alpha, n_samples,
                                        n_components)[0]
    grad = lambda params: _kl_divergence(params, P, alpha, n_samples,
                                         n_components)[1]
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
    # The Barnes-Hut approximation uses a different method to estimate
    # P_ij using only a number of nearest neighbors instead of all
    # points (so that k = 3 * perplexity). As a result we set the
    # perplexity=5, so that the number of neighbors is 5%.
    n_components = 2
    methods = ['exact', 'barnes_hut']
    X = random_state.randn(100, n_components).astype(np.float32)
    for init in ('random', 'pca'):
        for method in methods:
            tsne = TSNE(n_components=n_components, perplexity=50,
                        learning_rate=100.0, init=init, random_state=0,
                        method=method)
            X_embedded = tsne.fit_transform(X)
            T = trustworthiness(X, X_embedded, n_neighbors=1)
            assert_almost_equal(T, 1.0, decimal=1)


def test_optimization_minimizes_kl_divergence():
    """t-SNE should give a lower KL divergence with more iterations."""
    random_state = check_random_state(0)
    X, _ = make_blobs(n_features=3, random_state=random_state)
    kl_divergences = []
    for n_iter in [200, 250, 300]:
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
    X = random_state.randn(100, 2)
    D = squareform(pdist(X), "sqeuclidean")
    tsne = TSNE(n_components=2, perplexity=2, learning_rate=100.0,
                metric="precomputed", random_state=0, verbose=0)
    X_embedded = tsne.fit_transform(D)
    assert_almost_equal(trustworthiness(D, X_embedded, n_neighbors=1,
                                        precomputed=True), 1.0, decimal=1)


def test_early_exaggeration_too_small():
    # Early exaggeration factor must be >= 1.
    tsne = TSNE(early_exaggeration=0.99)
    assert_raises_regexp(ValueError, "early_exaggeration .*",
                         tsne.fit_transform, np.array([[0.0]]))


def test_too_few_iterations():
    # Number of gradient descent iterations must be at least 200.
    tsne = TSNE(n_iter=199)
    assert_raises_regexp(ValueError, "n_iter .*", tsne.fit_transform,
                         np.array([[0.0]]))


def test_non_square_precomputed_distances():
    # Precomputed distance matrices must be square matrices.
    tsne = TSNE(metric="precomputed")
    assert_raises_regexp(ValueError, ".* square distance matrix",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_init_not_available():
    # 'init' must be 'pca' or 'random'.
    m = "'init' must be 'pca', 'random' or a NumPy array"
    assert_raises_regexp(ValueError, m, TSNE, init="not available")


def test_distance_not_available():
    # 'metric' must be valid.
    tsne = TSNE(metric="not available")
    assert_raises_regexp(ValueError, "Unknown metric not available.*",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_pca_initialization_not_compatible_with_precomputed_kernel():
    # Precomputed distance matrices must be square matrices.
    tsne = TSNE(metric="precomputed", init="pca")
    assert_raises_regexp(ValueError, "The parameter init=\"pca\" cannot be "
                         "used with metric=\"precomputed\".",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


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
    neighbors = neighbors.astype(np.int64)
    pij_input = _joint_probabilities(*args)
    pij_input = squareform(pij_input).astype(np.float32)
    grad_bh = np.zeros(pos_output.shape, dtype=np.float32)

    _barnes_hut_tsne.gradient(pij_input, pos_output, neighbors,
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
    assert("Computing pairwise distances" in out)
    assert("Computed conditional probabilities" in out)
    assert("Mean sigma" in out)
    assert("Finished" in out)
    assert("early exaggeration" in out)
    assert("Finished" in out)


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


def test_64bit():
    # Ensure 64bit arrays are handled correctly.
    random_state = check_random_state(0)
    methods = ['barnes_hut', 'exact']
    for method in methods:
        for dt in [np.float32, np.float64]:
            X = random_state.randn(100, 2).astype(dt)
            tsne = TSNE(n_components=2, perplexity=2, learning_rate=100.0,
                        random_state=0, method=method)
            tsne.fit_transform(X)


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
        distances = distances.dot(distances.T)
        np.fill_diagonal(distances, 0.0)
        params = random_state.randn(n_samples, n_components)
        P = _joint_probabilities(distances, perplexity, False)
        kl, gradex = _kl_divergence(params, P, degrees_of_freedom, n_samples,
                                    n_components)

        k = n_samples - 1
        bt = BallTree(distances)
        distances_nn, neighbors_nn = bt.query(distances, k=k + 1)
        neighbors_nn = neighbors_nn[:, 1:]
        Pbh = _joint_probabilities_nn(distances, neighbors_nn,
                                      perplexity, False)
        kl, gradbh = _kl_divergence_bh(params, Pbh, neighbors_nn,
                                       degrees_of_freedom, n_samples,
                                       n_components, angle=angle,
                                       skip_num_points=0, verbose=False)
        assert_array_almost_equal(Pbh, P, decimal=5)
        assert_array_almost_equal(gradex, gradbh, decimal=5)


def test_quadtree_similar_point():
    # Introduce a point into a quad tree where a similar point already exists.
    # Test will hang if it doesn't complete.
    Xs = []

    # check the case where points are actually different
    Xs.append(np.array([[1, 2], [3, 4]], dtype=np.float32))
    # check the case where points are the same on X axis
    Xs.append(np.array([[1.0, 2.0], [1.0, 3.0]], dtype=np.float32))
    # check the case where points are arbitrarily close on X axis
    Xs.append(np.array([[1.00001, 2.0], [1.00002, 3.0]], dtype=np.float32))
    # check the case where points are the same on Y axis
    Xs.append(np.array([[1.0, 2.0], [3.0, 2.0]], dtype=np.float32))
    # check the case where points are arbitrarily close on Y axis
    Xs.append(np.array([[1.0, 2.00001], [3.0, 2.00002]], dtype=np.float32))
    # check the case where points are arbitrarily close on both axes
    Xs.append(np.array([[1.00001, 2.00001], [1.00002, 2.00002]],
              dtype=np.float32))

    # check the case where points are arbitrarily close on both axes
    # close to machine epsilon - x axis
    Xs.append(np.array([[1, 0.0003817754041], [2, 0.0003817753750]],
              dtype=np.float32))

    # check the case where points are arbitrarily close on both axes
    # close to machine epsilon - y axis
    Xs.append(np.array([[0.0003817754041, 1.0], [0.0003817753750, 2.0]],
              dtype=np.float32))

    for X in Xs:
        counts = np.zeros(3, dtype='int64')
        _barnes_hut_tsne.check_quadtree(X, counts)
        m = "Tree consistency failed: unexpected number of points at root node"
        assert_equal(counts[0], counts[1], m)
        m = "Tree consistency failed: unexpected number of points on the tree"
        assert_equal(counts[0], counts[2], m)


def test_index_offset():
    # Make sure translating between 1D and N-D indices are preserved
    assert_equal(_barnes_hut_tsne.test_index2offset(), 1)
    assert_equal(_barnes_hut_tsne.test_index_offset(), 1)
