# Author: Vlad Niculae
# License: BSD style

import warnings
from sys import version_info

import numpy as np

from nose import SkipTest
from nose.tools import assert_raises
from numpy.testing import assert_equal, assert_array_almost_equal

from .. import orthogonal_mp, orthogonal_mp_gram, OrthogonalMatchingPursuit
from ...utils.fixes import count_nonzero
from ...datasets import make_sparse_coded_signal

n_samples, n_features, n_nonzero_coefs, n_targets = 20, 30, 5, 3
y, X, gamma = make_sparse_coded_signal(n_targets, n_features, n_samples,
                                       n_nonzero_coefs, random_state=0)
G, Xy = np.dot(X.T, X), np.dot(X.T, y)
# this makes X (n_samples, n_features)
# and y (n_samples, 3)


def check_warnings():
    if version_info < (2, 6):
        raise SkipTest("Testing for warnings is not supported in versions \
        older than Python 2.6")


def test_correct_shapes():
    assert_equal(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=5).shape,
                 (n_features,))
    assert_equal(orthogonal_mp(X, y, n_nonzero_coefs=5).shape,
                 (n_features, 3))


def test_correct_shapes_gram():
    assert_equal(orthogonal_mp_gram(G, Xy[:, 0], n_nonzero_coefs=5).shape,
                 (n_features,))
    assert_equal(orthogonal_mp_gram(G, Xy, n_nonzero_coefs=5).shape,
                 (n_features, 3))


def test_n_nonzero_coefs():
    assert count_nonzero(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=5)) <= 5
    assert count_nonzero(orthogonal_mp(X, y[:, 0], n_nonzero_coefs=5,
                                       precompute_gram=True)) <= 5


def test_tol():
    tol = 0.5
    gamma = orthogonal_mp(X, y[:, 0], tol=tol)
    gamma_gram = orthogonal_mp(X, y[:, 0], tol=tol, precompute_gram=True)
    assert np.sum((y[:, 0] - np.dot(X, gamma)) ** 2) <= tol
    assert np.sum((y[:, 0] - np.dot(X, gamma_gram)) ** 2) <= tol


def test_with_without_gram():
    assert_array_almost_equal(
        orthogonal_mp(X, y, n_nonzero_coefs=5),
        orthogonal_mp(X, y, n_nonzero_coefs=5, precompute_gram=True))


def test_with_without_gram_tol():
    assert_array_almost_equal(
        orthogonal_mp(X, y, tol=1.),
        orthogonal_mp(X, y, tol=1., precompute_gram=True))


def test_unreachable_accuracy():
    check_warnings()  # Skip if unsupported Python version
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        assert_array_almost_equal(
            orthogonal_mp(X, y, tol=0),
            orthogonal_mp(X, y, n_nonzero_coefs=n_features))

        assert_array_almost_equal(
            orthogonal_mp(X, y, tol=0, precompute_gram=True),
            orthogonal_mp(X, y, precompute_gram=True,
                          n_nonzero_coefs=n_features))
        assert len(w) > 0  # warnings should be raised


def test_bad_input():
    assert_raises(ValueError, orthogonal_mp, X, y, tol=-1)
    assert_raises(ValueError, orthogonal_mp, X, y, n_nonzero_coefs=-1)
    assert_raises(ValueError, orthogonal_mp, X, y,
                  n_nonzero_coefs=n_features + 1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy, tol=-1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy, n_nonzero_coefs=-1)
    assert_raises(ValueError, orthogonal_mp_gram, G, Xy,
                  n_nonzero_coefs=n_features + 1)


def test_perfect_signal_recovery():
    # XXX: use signal generator
    idx, = gamma[:, 0].nonzero()
    gamma_rec = orthogonal_mp(X, y[:, 0], 5)
    gamma_gram = orthogonal_mp_gram(G, Xy[:, 0], 5)
    assert_equal(idx, np.flatnonzero(gamma_rec))
    assert_equal(idx, np.flatnonzero(gamma_gram))
    assert_array_almost_equal(gamma[:, 0], gamma_rec, decimal=2)
    assert_array_almost_equal(gamma[:, 0], gamma_gram, decimal=2)


def test_estimator_shapes():
    omp = OrthogonalMatchingPursuit(n_nonzero_coefs=n_nonzero_coefs)
    omp.fit(X, y[:, 0])
    assert_equal(omp.coef_.shape, (n_features,))
    assert_equal(omp.intercept_.shape, ())
    assert count_nonzero(omp.coef_) <= n_nonzero_coefs

    omp.fit(X, y)
    assert_equal(omp.coef_.shape, (n_targets, n_features))
    assert_equal(omp.intercept_.shape, (n_targets,))
    assert count_nonzero(omp.coef_) <= n_targets * n_nonzero_coefs

    omp.fit(X, y[:, 0], Gram=G, Xy=Xy[:, 0])
    assert_equal(omp.coef_.shape, (n_features,))
    assert_equal(omp.intercept_.shape, ())
    assert count_nonzero(omp.coef_) <= n_nonzero_coefs

    omp.fit(X, y, Gram=G, Xy=Xy)
    assert_equal(omp.coef_.shape, (n_targets, n_features))
    assert_equal(omp.intercept_.shape, (n_targets,))
    assert count_nonzero(omp.coef_) <= n_targets * n_nonzero_coefs


def test_identical_regressors():
    check_warnings()  # Skip if unsupported Python version
    newX = X.copy()
    newX[:, 1] = newX[:, 0]
    gamma = np.zeros(n_features)
    gamma[0] = gamma[1] = 1.
    newy = np.dot(newX, gamma)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        orthogonal_mp(newX, newy, 2)
        assert len(w) == 1


def test_swapped_regressors():
    gamma = np.zeros(n_features)
    # X[:, 21] should be selected first, then X[:, 0] selected second,
    # which will take X[:, 21]'s place in case the algorithm does
    # column swapping for optimization (which is the case at the moment)
    gamma[21] = 1.0
    gamma[0] = 0.5
    new_y = np.dot(X, gamma)
    new_Xy = np.dot(X.T, new_y)
    gamma_hat = orthogonal_mp(X, new_y, 2)
    gamma_hat_gram = orthogonal_mp_gram(G, new_Xy, 2)
    assert_equal(np.flatnonzero(gamma_hat), [0, 21])
    assert_equal(np.flatnonzero(gamma_hat_gram), [0, 21])


def test_no_atoms():
    y_empty = np.zeros_like(y)
    Xy_empty = np.dot(X.T, y_empty)
    gamma_empty = orthogonal_mp(X, y_empty, 1)
    gamma_empty_gram = orthogonal_mp_gram(G, Xy_empty, 1)
    assert_equal(np.all(gamma_empty == 0), True)
    assert_equal(np.all(gamma_empty_gram == 0), True)
