import numpy as np
from nose.tools import assert_equals
from sklearn.manifold.t_sne import _gradient_descent


def test_gradient_descent_stops():
    class ObjectiveSmallGradient:
        def __init__(self):
            self.it = -1

        def __call__(self, params):
            self.it += 1
            return 10 - self.it, np.array([1e-5])

    p, error, it = _gradient_descent(ObjectiveSmallGradient(), np.zeros(1), 0,
        n_iter=100, n_iter_without_progress=100, momentum=0.0,
        learning_rate=0.0, min_gain=0.0, min_grad_norm=1e-5,
        min_error_diff=0.0)
    assert_equals(error, 10)
    assert_equals(it, 0)
