import numpy as np
from nose.tools import assert_equals, assert_almost_equal
from sklearn.manifold.t_sne import _gradient_descent
from sklearn.manifold._binary_search import _binary_search_perplexity


def test_gradient_descent_stops():
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
    assert_equals(error, 1.0)
    assert_equals(it, 0)

    # Error difference
    _, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=0.2)
    assert_equals(error, 0.9)
    assert_equals(it, 1)

    # Maximum number of iterations without improvement
    _, error, it = _gradient_descent(flat_function, np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=10, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=-1.0)
    assert_equals(error, 0.0)
    assert_equals(it, 11)

    # Maximum number of iterations
    _, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=11, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=0.0,
        min_error_diff=0.0)
    assert_equals(error, 0.0)
    assert_equals(it, 10)


def test_binary_search():
    dist = np.random.randn(100, 2)
    dist = dist.dot(dist.T)
    np.fill_diagonal(dist, 0.0)
    desired_perplexity = 50.0
    P = _binary_search_perplexity(dist, desired_perplexity, verbose=1)
    P = np.maximum(P, np.finfo(np.double).eps)
    mean_perplexity = np.mean([np.exp(-np.sum(P[i] * np.log(P[i])))
                               for i in range(P.shape[0])])
    assert_almost_equal(mean_perplexity, desired_perplexity, places=4)
