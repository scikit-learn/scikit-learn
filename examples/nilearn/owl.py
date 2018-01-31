import numpy as np
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.isotonic import isotonic_regression

rng_global = check_random_state(0)


def invert_permutation(w, perm):
    inv = np.empty(w.size, w.dtype)
    inv[perm] = w
    return inv


def owl(w, weights, increasing=False):
    absw = np.sort(np.abs(w))
    if increasing:
        absw = absw[:--1]
    return np.dot(absw, weights)


def prox_owl(w, lambd, weights, increasing=False):
    """Proximal operator of the OWL norm using isotonic regression.
    The time complexity is O(n log n), where n = len(w).
    """
    sample_weights = np.copy(weights)
    sample_weights *= lambd
    signs = np.sign(w)
    w = np.abs(w)
    perm = np.argsort(w)
    if not increasing:
        perm = perm[::-1]
    w = w[perm]
    w -= sample_weights
    w = isotonic_regression(w, y_min=0., increasing=increasing)
    w = invert_permutation(w, perm)
    w *= signs
    return w


def _make_oscar_weights(alpha, n, increasing=False):
    assert 0. <= alpha <= 1. / (n - 1.)
    weights = np.ones(n, dtype=np.float)
    weights -= alpha * np.arange(n)
    if increasing:
        weights = weights[::-1]
    return weights


def prox_oscar(w, lambd, alpha, increasing=False):
    weights = _make_oscar_weights(alpha, len(w), increasing=increasing)
    return prox_owl(w, lambd, weights, increasing=increasing)


def test_owl_generalizes_l1():
    from nilearn.decoding.proximal_operators import _prox_l1
    w = np.arange(100000.)
    weights = np.ones(len(w), dtype=float)
    assert_array_almost_equal(prox_owl(w.copy(), 5., weights),
                              _prox_l1(w.copy(), 5.), decimal=10)

    w = rng_global.randn(len(w))
    assert_array_almost_equal(prox_owl(w.copy(), 1., weights),
                              _prox_l1(w.copy(), 1.), decimal=10)
