# Author: Vlad Niculae
# License: BSD style

import numpy as np
from nose.tools import assert_raises
from numpy.testing import assert_equal, assert_array_almost_equal
from nose.plugins.skip import SkipTest

from .. import orthogonal_mp, orthogonal_mp_gram, OrthogonalMatchingPursuit
from ...utils.fixes import count_nonzero
from ...utils import check_random_state
from ...datasets import generate_sparse_coded_signal

n_samples, n_features, n_nonzero_coefs = 10, 15, 6
y, X, gamma = generate_sparse_coded_signal(3, n_features, n_samples,
                                           n_nonzero_coefs, random_state=0)

# this makes X (n_samples, n_features)
#        and y (n_samples, 3)

def generate_data(n_samples, n_features, random_state=42):
    rng = check_random_state(random_state)
    X = rng.randn(n_samples, n_features)
    X /= np.sqrt(np.sum((X ** 2), axis=0))
    gamma = rng.randn(n_features)
    return X, np.dot(X, gamma)

# XXX: change samples_generator to the transpose problem, makes more sense


def test_correct_shapes():
    assert_equal(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=6).shape,
                 (n_features,))
    assert_equal(orthogonal_mp(X, y, n_nonzero_coefs=6).shape,
                 (n_features, 3))


def test_correct_shapes_gram():
    G, Xy = np.dot(X.T, X), np.dot(X.T, y[:, 0])
    assert_equal(orthogonal_mp_gram(G, Xy, n_nonzero_coefs=6).shape, 
                 (n_features,))
    Xy = np.dot(X.T, y)
    assert_equal(orthogonal_mp_gram(G, Xy, n_nonzero_coefs=6).shape,
                 (n_features, 3))


def test_n_nonzero_coefs():
    assert count_nonzero(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=6)) <= 6
    assert count_nonzero(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=6,
                                       compute_gram=True)) <= 6


def test_eps():
    eps = 0.5
    gamma = orthogonal_mp(X, y[:, 0], eps=eps)
    gamma_gram = orthogonal_mp(X, y[:, 0], eps=eps, compute_gram=True)
    assert np.sum((y[:, 0] - np.dot(X, gamma)) ** 2) <= eps
    assert np.sum((y[:, 0] - np.dot(X, gamma_gram)) ** 2) <= eps


def test_with_without_gram():
    assert_array_almost_equal(orthogonal_mp(X, y, n_nonzero_coefs=6),
                              orthogonal_mp(X, y, n_nonzero_coefs=6,
                                            compute_gram=True))


def test_with_without_gram_eps():
    assert_array_almost_equal(orthogonal_mp(X, y, eps=1.),
                              orthogonal_mp(X, y, eps=1., compute_gram=True))


def test_unreachable_accuracy():
    assert_array_almost_equal(orthogonal_mp(X, y, eps=0),
                              orthogonal_mp(X, y, n_nonzero_coefs=n_features))


def test_bad_input():
    G, Xy = np.dot(X.T, X), np.dot(X.T, y)
    assert_raises(ValueError, orthogonal_mp, X, y, eps=-1)
    assert_raises(ValueError, orthogonal_mp, X, y, n_nonzero_coefs=-1)
    assert_raises(ValueError, orthogonal_mp, X, y,
                  n_nonzero_coefs=n_features + 1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy, eps=-1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy, n_nonzero_coefs=-1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy,
                  n_nonzero_coefs=n_features + 1)


def test_perfect_signal_recovery():
    # XXX: use signal generator
    G, Xy = np.dot(X.T, X), np.dot(X.T, y[:, 0])
    idx, = gamma[:, 0].nonzero()
    gamma_rec = orthogonal_mp(X, y[:, 0], 6)
    gamma_gram = orthogonal_mp_gram(G, Xy, 6)
    assert_equal(idx, np.flatnonzero(gamma_rec))
    assert_equal(idx, np.flatnonzero(gamma_gram))
    assert_array_almost_equal(gamma, gamma_rec, decimal=2)
    assert_array_almost_equal(gamma, gamma_gram, decimal=2)


def test_estimator_shapes():
    OMP = OrthogonalMatchingPursuit(n_nonzero_coefs=6)
    OMP.fit(X, y[:, 0])
    assert_equal(OMP.coef_.shape, (15, ))
    assert_equal(OMP.intercept_.shape, ())
    assert count_nonzero(OMP.coef_) <= 6

    OMP.fit(X, y)
    assert_equal(OMP.coef_.shape, (3, 15))
    assert_equal(OMP.intercept_.shape, (3, ))
    assert count_nonzero(OMP.coef_) <= 3 * 6
