import numpy as np
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.manifold.t_sne import _joint_probabilities
from sklearn.manifold.t_sne import _kl_divergence
from sklearn.manifold.t_sne import _gradient_descent
from sklearn.manifold.t_sne import trustworthiness
from sklearn.manifold.t_sne import TSNE
from sklearn.manifold._binary_search import _binary_search_perplexity
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
    _, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=1e-5,
        min_error_diff=0.0)
    assert_equal(error, 1.0)
    assert_equal(it, 0)

    # Error difference
    _, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=0.2)
    assert_equal(error, 0.9)
    assert_equal(it, 1)

    # Maximum number of iterations without improvement
    _, error, it = _gradient_descent(flat_function, np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=10, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=-1.0)
    assert_equal(error, 0.0)
    assert_equal(it, 11)

    # Maximum number of iterations
    _, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=11, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=0.0)
    assert_equal(error, 0.0)
    assert_equal(it, 10)


def test_binary_search():
    """Test if the binary search finds Gaussians with desired perplexity."""
    affinities = np.random.randn(50, 2)
    affinities = affinities.dot(affinities.T)
    np.fill_diagonal(affinities, 0.0)
    desired_perplexity = 25.0
    P = _binary_search_perplexity(affinities, desired_perplexity, verbose=1)
    P = np.maximum(P, np.finfo(np.double).eps)
    mean_perplexity = np.mean([np.exp(-np.sum(P[i] * np.log(P[i])))
                               for i in range(P.shape[0])])
    assert_almost_equal(mean_perplexity, desired_perplexity, decimal=3)


def test_gradient():
    """Test gradient of Kullback-Leibler divergence."""
    n_samples = 50
    n_features = 2
    n_components = 2
    alpha = 1.0

    affinities = np.random.randn(n_samples, n_features)
    affinities = affinities.dot(affinities.T)
    np.fill_diagonal(affinities, 0.0)
    X_embedded = np.random.randn(n_samples, n_components)

    P = _joint_probabilities(affinities, desired_perplexity=25.0,
                             verbose=0)
    fun = lambda params: _kl_divergence(params, P, alpha, n_samples,
                                        n_components)[0]
    grad = lambda params: _kl_divergence(params, P, alpha, n_samples,
                                         n_components)[1]
    assert_almost_equal(check_grad(fun, grad, X_embedded.ravel()), 0.0,
                        decimal=5)

def test_trustworthiness():
    """Test trustworthiness score."""
    # Affine transformation
    X = np.random.randn(100, 2)
    assert_equal(trustworthiness(X, 5.0 + X / 10.0), 1.0)

    # Randomly shuffled
    X = np.arange(100).reshape(-1, 1)
    X_embedded = X.copy()
    np.random.shuffle(X_embedded)
    assert_less(trustworthiness(X, X_embedded), 0.6)

    # Completely different
    X = np.arange(5).reshape(-1, 1)
    X_embedded = np.array([[0], [2], [4], [1], [3]])
    assert_almost_equal(trustworthiness(X, X_embedded, n_neighbors=1), 0.2)


def test_preserve_trustworthiness_approximately():
    """Nearest neighbors should be preserved approximately."""
    X = np.random.randn(100, 2)
    tsne = TSNE(n_components=2, perplexity=10, random_state=0)
    X_embedded = tsne.fit_transform(X)
    assert_almost_equal(tsne.score(X, n_neighbors=1), 1.0, decimal=1)


def test_preserve_trustworthiness_approximately_with_precomputed_affinities():
    """Nearest neighbors should be preserved approximately."""
    X = np.random.randn(100, 2)
    D = squareform(pdist(X), "sqeuclidean")
    tsne = TSNE(n_components=2, perplexity=10, affinity="precomputed",
                random_state=0)
    tsne.fit(D)
    assert_almost_equal(tsne.score(D, n_neighbors=1), 1.0, decimal=1)


def test_transform_inverse_transform():
    """Training data can be reconstructed approximately."""
    X = np.random.randn(5, 2)
    tsne = TSNE(n_components=2, fit_inverse_transform=True, random_state=0)
    tsne.fit(X)
    X_reconstructed = tsne.inverse_transform(tsne.transform(X))
    assert_array_almost_equal(X, X_reconstructed)


def test_precomputed_affinities_inverse_transform():
    """Cannot fit inverse tranform with precomputed affinities."""
    assert_raises_regexp(ValueError, "Cannot fit_inverse_transform",
        TSNE, affinity="precomputed", fit_inverse_transform=True)


def test_early_exaggeration_too_small():
    """Early exaggeration factor must be >= 1."""
    tsne = TSNE(early_exaggeration=0.99)
    assert_raises_regexp(ValueError, "early_exaggeration", tsne.fit,
                         np.array([[0.0]]))


def test_too_few_iterations():
    """Number of gradient descent iterations must be at least 200."""
    tsne = TSNE(n_iter=199)
    assert_raises_regexp(ValueError, "n_iter", tsne.fit, np.array([[0.0]]))


def test_non_square_precomputed_affinities():
    """Precomputed affinity matrices must be sqaure matrices."""
    tsne = TSNE(affinity="precomputed")
    assert_raises_regexp(ValueError, "square affinity matrix", tsne.fit,
                         np.array([[0.0], [1.0]]))


def test_inverse_transform_not_fitted():
    """The inverse transform must be activated."""
    X = np.ones((2, 2))
    tsne = TSNE(fit_inverse_transform=False)
    tsne.fit(X)
    assert_raises_regexp(ValueError, "not fitted", tsne.inverse_transform, X)


def test_cannot_generalize_with_precomputed_affinities():
    """Cannot generalize with precomputed affinity matrix."""
    X_train = np.ones((2, 2))
    X_test = np.zeros((2, 2))
    tsne = TSNE(affinity="precomputed")
    tsne.fit(X_train)
    assert_raises_regexp(ValueError, "Can only transform training data",
                         tsne.transform, X_test)

# TODO tests with verbose = 1 or 2
