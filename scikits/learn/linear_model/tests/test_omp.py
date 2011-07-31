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

def generate_data(n_samples, n_features, random_state=42):
    rng = check_random_state(random_state)
    X = rng.randn(n_samples, n_features)
    X /= np.sqrt(np.sum((X ** 2), axis=0))
    gamma = rng.randn(n_features)
    return X, np.dot(X, gamma)

# XXX: change samples_generator to the transpose problem, makes more sense


def test_correct_shapes():
    y, _, X = generate_sparse_coded_signal(3, n_features, n_samples,
                                           n_nonzero_coefs, random_state=0)
    assert_equal(orthogonal_mp(X.T, y[0], n_nonzero_coefs=6).shape,
                 (n_features,))

    assert_equal(orthogonal_mp(X.T, y.T, n_nonzero_coefs=6).shape,
                 (n_features, 3))


def test_correct_shapes_gram():
    y, _, X = generate_sparse_coded_signal(3, n_features, n_samples,
                                           n_nonzero_coefs, random_state=0)
    G, Xy = np.dot(X, X.T), np.dot(X, y[0])
    assert_equal(orthogonal_mp_gram(G, Xy, n_nonzero_coefs=6).shape, 
                 (n_features,))

    Xy = np.dot(X, y.T)
    assert_equal(orthogonal_mp_gram(G, Xy, n_nonzero_coefs=6).shape,
                 (n_features, 3))


def test_n_nonzero_coefs():
    n_samples, n_features = 10, 15
    X, y = generate_data(n_samples, n_features)
    assert_equal(6, count_nonzero(orthogonal_mp(X, y, n_nonzero_coefs=6)))
    assert_equal(6, count_nonzero(orthogonal_mp(X, y, n_nonzero_coefs=6,
                                                   compute_gram=True)))


def test_eps():
    n_samples, n_features = 10, 15
    eps = 0.5
    X, y = generate_data(n_samples, n_features)
    gamma = orthogonal_mp(X, y, eps=eps)
    gamma_gram = orthogonal_mp(X, y, eps=eps, compute_gram=True)
    assert(np.sum((y - np.dot(X, gamma)) ** 2) <= eps)
    assert(np.sum((y - np.dot(X, gamma_gram)) ** 2) <= eps)


def test_with_without_gram():
    n_samples, n_features = 10, 15
    X, y = generate_data(n_samples, n_features)
    assert_array_almost_equal(orthogonal_mp(X, y, n_nonzero_coefs=6),
                              orthogonal_mp(X, y, n_nonzero_coefs=6,
                                            compute_gram=True))


def test_with_without_gram_eps():
    n_samples, n_features = 10, 15
    X, y = generate_data(n_samples, n_features)
    assert_array_almost_equal(orthogonal_mp(X, y, eps=1.),
                              orthogonal_mp(X, y, eps=1., compute_gram=True))


def test_unreachable_accuracy():
    n_samples, n_features = 10, 15
    X, y = generate_data(n_samples, n_features)
    assert_array_almost_equal(orthogonal_mp(X, y, eps=0),
                              orthogonal_mp(X, y, n_nonzero_coefs=n_features))


def test_bad_input():
    n_samples, n_features = 10, 15
    X, y = generate_data(n_samples, n_features)
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
    rng = np.random.RandomState(0)
    n_samples, n_features = 10, 15
    n_atoms = 4
    X, _ = generate_data(n_samples, n_features, rng)
    gamma = np.zeros(n_features)
    idx = np.arange(n_features)
    rng.shuffle(idx)
    idx = idx[:n_atoms]
    idx = np.sort(idx)
    gamma[idx] = rng.randn(n_atoms)
    y = np.dot(X, gamma)
    G, Xy = np.dot(X.T, X), np.dot(X.T, y)

    gamma_rec = orthogonal_mp(X, y, n_atoms)
    gamma_gram = orthogonal_mp_gram(G, Xy, n_atoms)
    assert_equal(idx, np.flatnonzero(gamma_rec))
    assert_equal(idx, np.flatnonzero(gamma_gram))
    assert_array_almost_equal(gamma, gamma_rec, decimal=2)
    assert_array_almost_equal(gamma, gamma_gram, decimal=2)


def test_estimator_shapes():
    X = np.random.randn(10, 15)
    y1 = np.random.randn(10)
    y2 = np.random.randn(10, 3)

    OMP = OrthogonalMatchingPursuit(n_nonzero_coefs=5)
    OMP.fit(X, y1)
    assert_equal((15,), OMP.coef_.shape)
    assert_equal((), OMP.intercept_.shape)
    assert_equal(5, count_nonzero(OMP.coef_))

    OMP.fit(X, y2)
    assert_equal((3, 15), OMP.coef_.shape)
    assert_equal((3,), OMP.intercept_.shape)
    assert_equal(3 * 5, count_nonzero(OMP.coef_))
