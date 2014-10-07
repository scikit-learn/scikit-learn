import sys
from sklearn.externals.six.moves import cStringIO as StringIO
import numpy as np
import scipy.sparse as sp
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils import check_random_state
from sklearn.manifold.t_sne import _joint_probabilities
from sklearn.manifold.t_sne import _kl_divergence
from sklearn.manifold.t_sne import _gradient_descent
from sklearn.manifold.t_sne import trustworthiness
from sklearn.manifold.t_sne import TSNE
from sklearn.manifold._utils import _binary_search_perplexity
from scipy.optimize import check_grad
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform


def test_gradient_descent_stops():
    """Test stopping conditions of gradient descent."""
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
    """Test if the binary search finds Gaussians with desired perplexity."""
    random_state = check_random_state(0)
    distances = random_state.randn(50, 2)
    distances = distances.dot(distances.T)
    np.fill_diagonal(distances, 0.0)
    desired_perplexity = 25.0
    P = _binary_search_perplexity(distances, desired_perplexity, verbose=0)
    P = np.maximum(P, np.finfo(np.double).eps)
    mean_perplexity = np.mean([np.exp(-np.sum(P[i] * np.log(P[i])))
                               for i in range(P.shape[0])])
    assert_almost_equal(mean_perplexity, desired_perplexity, decimal=3)


def test_gradient():
    """Test gradient of Kullback-Leibler divergence."""
    random_state = check_random_state(0)

    n_samples = 50
    n_features = 2
    n_components = 2
    alpha = 1.0

    distances = random_state.randn(n_samples, n_features)
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
    """Test trustworthiness score."""
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
    """Nearest neighbors should be preserved approximately."""
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    for init in ('random', 'pca'):
        tsne = TSNE(n_components=2, perplexity=10, learning_rate=100.0,
                    init=init, random_state=0)
        X_embedded = tsne.fit_transform(X)
        assert_almost_equal(trustworthiness(X, X_embedded, n_neighbors=1), 1.0,
                            decimal=1)


def test_fit_csr_matrix():
    """X can be a sparse matrix."""
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    X[(np.random.randint(0, 100, 50), np.random.randint(0, 2, 50))] = 0.0
    X_csr = sp.csr_matrix(X)
    tsne = TSNE(n_components=2, perplexity=10, learning_rate=100.0,
                random_state=0)
    X_embedded = tsne.fit_transform(X_csr)
    assert_almost_equal(trustworthiness(X_csr, X_embedded, n_neighbors=1), 1.0,
                        decimal=1)


def test_preserve_trustworthiness_approximately_with_precomputed_distances():
    """Nearest neighbors should be preserved approximately."""
    random_state = check_random_state(0)
    X = random_state.randn(100, 2)
    D = squareform(pdist(X), "sqeuclidean")
    tsne = TSNE(n_components=2, perplexity=10, learning_rate=100.0,
                metric="precomputed", random_state=0)
    X_embedded = tsne.fit_transform(D)
    assert_almost_equal(trustworthiness(D, X_embedded, n_neighbors=1,
                                        precomputed=True), 1.0, decimal=1)


def test_early_exaggeration_too_small():
    """Early exaggeration factor must be >= 1."""
    tsne = TSNE(early_exaggeration=0.99)
    assert_raises_regexp(ValueError, "early_exaggeration .*",
                         tsne.fit_transform, np.array([[0.0]]))


def test_too_few_iterations():
    """Number of gradient descent iterations must be at least 200."""
    tsne = TSNE(n_iter=199)
    assert_raises_regexp(ValueError, "n_iter .*", tsne.fit_transform,
                         np.array([[0.0]]))


def test_non_square_precomputed_distances():
    """Precomputed distance matrices must be square matrices."""
    tsne = TSNE(metric="precomputed")
    assert_raises_regexp(ValueError, ".* square distance matrix",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_init_not_available():
    """'init' must be 'pca' or 'random'."""
    assert_raises_regexp(ValueError, "'init' must be either 'pca' or 'random'",
                         TSNE, init="not available")


def test_distance_not_available():
    """'metric' must be valid."""
    tsne = TSNE(metric="not available")
    assert_raises_regexp(ValueError, "Unknown metric not available.*",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_pca_initialization_not_compatible_with_precomputed_kernel():
    """Precomputed distance matrices must be square matrices."""
    tsne = TSNE(metric="precomputed", init="pca")
    assert_raises_regexp(ValueError, "The parameter init=\"pca\" cannot be "
                         "used with metric=\"precomputed\".",
                         tsne.fit_transform, np.array([[0.0], [1.0]]))


def test_verbose():
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
