# Authors: Jaques Grobler <jaques.grobler@inria.fr>
# License: BSD Style.

import warnings
from sys import version_info

import numpy as np

from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater

#from sklearn import datasets
from sklearn.linear_model.nng import NonNegativeGarrote
from sklearn.linear_model import LinearRegression
from sklearn.linear_model import non_negative_garrote_path
from sklearn import datasets

diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target


def test_nng_zero():
    """Check that the nng can handle zero data without crashing"""
    X = [[0], [0], [0]]
    y = [0, 0, 0]
    clf = NonNegativeGarrote(alpha=0.1).fit(X, y)
    pred = clf.predict([[1], [2], [3]])
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(pred, [0, 0, 0])


def test_nng_toy():
    """
    Test nng on a toy example for various values of alpha.
    """

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # a straight line
    T = [[2], [3], [4]]  # test sample

    clf = NonNegativeGarrote(alpha=1e-8)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])

    clf = NonNegativeGarrote(alpha=0.1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.85])
    assert_array_almost_equal(pred, [1.7, 2.55, 3.4])

    clf = NonNegativeGarrote(alpha=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.25])
    assert_array_almost_equal(pred, [0.5, 0.75, 1.])

    clf = NonNegativeGarrote(alpha=1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.0])
    assert_array_almost_equal(pred, [0, 0, 0])


def test_nng_alpha_warning():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        X = [[-1], [0], [1]]
        Y = [-1, 0, 1]       # just a straight line

        clf = NonNegativeGarrote(alpha=0)
        clf.fit(X, Y)

        assert_greater(len(w), 0)  # warnings should be raised


def test_nng_positive_constraint():
    """
    Check that the nng performs it's positive constraint on the
    shrinkage coefficients
    """
    X = [[-1], [0], [1]]
    y = [1, 0, -1]       # just a straight line with negative slope

    nngarrote = NonNegativeGarrote(alpha=0.1, max_iter=1000)
    nngarrote.fit(X, y)
    assert_true(min(nngarrote.shrink_coef_) >= 0)

    nngarrote = NonNegativeGarrote(alpha=0.1, max_iter=1000,
                                   precompute=True)
    nngarrote.fit(X, y)
    assert_true(min(nngarrote.shrink_coef_) >= 0)


def test_small_alpha():
    """
    Test to see that if alpha goes small, the coefs will be the same
    as OLS
    """
    X = [[-1], [0], [1]]
    y = [-1, 0, 1]       # a straight line

    nng_coef = NonNegativeGarrote(alpha=1e-8).fit(X, y).coef_
    ols_coef = LinearRegression().fit(X, y).coef_
    assert_array_almost_equal(nng_coef, ols_coef)


def build_dataset(n_samples=50, n_features=200, n_informative_features=10,
                  n_targets=1):
    """
    build an ill-posed linear regression problem with many noisy features and
    comparatively few samples
    """
    random_state = np.random.RandomState(0)
    if n_targets > 1:
        w = random_state.randn(n_features, n_targets)
    else:
        w = random_state.randn(n_features)
    w[n_informative_features:] = 0.0
    X = random_state.randn(n_samples, n_features)
    y = np.dot(X, w)
    X_test = random_state.randn(n_samples, n_features)
    y_test = np.dot(X_test, w)
    return X, y, X_test, y_test


def test_positive_well_conditioned():
    pass


def test_less_sample_than_dimentions():
    pass


def test_lasso_gives_lstsq_solution():
    """
    Test that Non-Negative Garrote gives least square solution at the end
    of the path
    """
    coef_path_, _ = non_negative_garrote_path(X, y)
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq, coef_path_[:, -1])


def test_singular_matrix():
    """
    Test when input is a singular matrix
    """
    X1 = np.array([[1, 1.], [1., 1.]])
    y1 = np.array([1, 1])
    coef_path, scp = non_negative_garrote_path(X1, y1)
    assert_array_almost_equal(coef_path.T, [[0, 0], [1, 0]])


def test_lars_add_features():
    """
    assure that at least some features get added if necessary

    test for 6d2b4c
    """
    # Hilbert matrix
    n = 5
    H = 1. / (np.arange(1, n + 1) + np.arange(n)[:, np.newaxis])
    clf = NonNegativeGarrote(fit_intercept=False).fit(
        H, np.arange(n))
    assert_true(np.all(np.isfinite(clf.coef_)))


if __name__ == '__main__':
    import nose
    nose.runmodule()
