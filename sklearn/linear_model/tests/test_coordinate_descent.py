# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

import warnings
from sys import version_info

import numpy as np
from scipy import interpolate

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater

from sklearn.linear_model.coordinate_descent import Lasso, \
    LassoCV, ElasticNet, ElasticNetCV, MultiTaskLasso, MultiTaskElasticNet, \
    lasso_path
from sklearn.linear_model import LassoLarsCV, lars_path


def check_warnings():
    if version_info < (2, 6):
        raise SkipTest("Testing for warnings is not supported in versions \
        older than Python 2.6")


def test_lasso_zero():
    """Check that the lasso can handle zero data without crashing"""
    X = [[0], [0], [0]]
    y = [0, 0, 0]
    clf = Lasso(alpha=0.1).fit(X, y)
    pred = clf.predict([[1], [2], [3]])
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_, 0)


def test_lasso_toy():
    """
    Test Lasso on a toy example for various values of alpha.

    When validating this against glmnet notice that glmnet divides it
    against nobs.
    """

    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line
    T = [[2], [3], [4]]  # test sample

    clf = Lasso(alpha=1e-8)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=0.1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.85])
    assert_array_almost_equal(pred, [1.7, 2.55, 3.4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.25])
    assert_array_almost_equal(pred, [0.5, 0.75, 1.])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = Lasso(alpha=1)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [.0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_, 0)


def test_enet_toy():
    """
    Test ElasticNet for various parameters of alpha and l1_ratio.

    Actually, the parameters alpha = 0 should not be allowed. However,
    we test it as a border case.

    ElasticNet is tested with and without precomputed Gram matrix
    """

    X = np.array([[-1.], [0.], [1.]])
    Y = [-1, 0, 1]       # just a straight line
    T = [[2.], [3.], [4.]]  # test sample

    # this should be the same as lasso
    clf = ElasticNet(alpha=1e-8, l1_ratio=1.0)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(pred, [2, 3, 4])
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, l1_ratio=0.3, max_iter=100,
                     precompute=False)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=True)
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf.set_params(max_iter=100, precompute=np.dot(X.T, X))
    clf.fit(X, Y)  # with Gram
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.50819], decimal=3)
    assert_array_almost_equal(pred, [1.0163, 1.5245, 2.0327], decimal=3)
    assert_almost_equal(clf.dual_gap_, 0)

    clf = ElasticNet(alpha=0.5, l1_ratio=0.5)
    clf.fit(X, Y)
    pred = clf.predict(T)
    assert_array_almost_equal(clf.coef_, [0.45454], 3)
    assert_array_almost_equal(pred, [0.9090, 1.3636, 1.8181], 3)
    assert_almost_equal(clf.dual_gap_, 0)


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


def test_lasso_cv():
    X, y, X_test, y_test = build_dataset()
    max_iter = 150
    clf = LassoCV(n_alphas=10, eps=1e-3, max_iter=max_iter).fit(X, y)
    assert_almost_equal(clf.alpha_, 0.056, 2)

    clf = LassoCV(n_alphas=10, eps=1e-3, max_iter=max_iter, precompute=True)
    clf.fit(X, y)
    assert_almost_equal(clf.alpha_, 0.056, 2)

    # Check that the lars and the coordinate descent implementation
    # select a similar alpha
    lars = LassoLarsCV(normalize=False, max_iter=30).fit(X, y)
    # for this we check that they don't fall in the grid of
    # clf.alphas further than 1
    assert_true(np.abs(
        np.searchsorted(clf.alphas_[::-1], lars.alpha_)
        - np.searchsorted(clf.alphas_[::-1], clf.alpha_)) <= 1)
    # check that they also give a similar MSE
    mse_lars = interpolate.interp1d(lars.cv_alphas_, lars.cv_mse_path_.T)
    np.testing.assert_approx_equal(mse_lars(clf.alphas_[5]).mean(),
                                   clf.mse_path_[5].mean(), significant=2)

    # test set
    assert_greater(clf.score(X_test, y_test), 0.99)


def test_lasso_path_return_models_vs_new_return_gives_same_coefficients():
    # Test that lasso_path with lars_path style output gives the
    # same result

    # Some toy data
    X = np.array([[1, 2, 3.1], [2.3, 5.4, 4.3]]).T
    y = np.array([1, 2, 3.1])
    alphas = [5., 1., .5]
    # Compute the lasso_path
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        coef_path = [e.coef_ for e in lasso_path(X, y, alphas=alphas,
                                                 fit_intercept=False)]

    # Use lars_path and lasso_path(new output) with 1D linear interpolation
    # to compute the the same path
    alphas_lars, _, coef_path_lars = lars_path(X, y, method='lasso')
    coef_path_cont_lars = interpolate.interp1d(alphas_lars[::-1],
                                               coef_path_lars[:, ::-1])
    alphas_lasso2, coef_path_lasso2, _ = lasso_path(X, y, alphas=alphas,
                                                 fit_intercept=False,
                                                 return_models=False)
    coef_path_cont_lasso = interpolate.interp1d(alphas_lasso2[::-1],
                                                coef_path_lasso2[:, ::-1])

    np.testing.assert_array_almost_equal(coef_path_cont_lasso(alphas),
                                         np.asarray(coef_path).T, decimal=1)
    np.testing.assert_array_almost_equal(coef_path_cont_lasso(alphas),
                                         coef_path_cont_lars(alphas),
                                         decimal=1)


def test_enet_path():
    # We use a large number of samples and of informative features so that
    # the l1_ratio selected is more toward ridge than lasso
    X, y, X_test, y_test = build_dataset(n_samples=200, n_features=100,
                                         n_informative_features=100)
    max_iter = 150

    with warnings.catch_warnings():
        # Here we have a small number of iterations, and thus the
        # ElasticNet might not converge. This is to speed up tests
        warnings.simplefilter("ignore", UserWarning)
        clf = ElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7], cv=3,
                           max_iter=max_iter)
        clf.fit(X, y)
        # Well-conditionned settings, we should have selected our
        # smallest penalty
        assert_almost_equal(clf.alpha_, min(clf.alphas_))
        # Non-sparse ground truth: we should have seleted an elastic-net
        # that is closer to ridge than to lasso
        assert_equal(clf.l1_ratio_, min(clf.l1_ratio))

        clf = ElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7], cv=3,
                           max_iter=max_iter, precompute=True)
        clf.fit(X, y)

    # Well-conditionned settings, we should have selected our
    # smallest penalty
    assert_almost_equal(clf.alpha_, min(clf.alphas_))
    # Non-sparse ground truth: we should have seleted an elastic-net
    # that is closer to ridge than to lasso
    assert_equal(clf.l1_ratio_, min(clf.l1_ratio))

    # We are in well-conditionned settings with low noise: we should
    # have a good test-set performance
    assert_greater(clf.score(X_test, y_test), 0.99)


def test_path_parameters():
    X, y, _, _ = build_dataset()
    max_iter = 50

    clf = ElasticNetCV(n_alphas=50, eps=1e-3, max_iter=max_iter,
                       l1_ratio=0.5)
    clf.fit(X, y)  # new params
    assert_almost_equal(0.5, clf.l1_ratio)
    assert_equal(50, clf.n_alphas)
    assert_equal(50, len(clf.alphas_))


def test_warm_start():
    X, y, _, _ = build_dataset()
    clf = ElasticNet(alpha=1.0, max_iter=50, warm_start=True)
    clf.fit(X, y)
    clf.set_params(alpha=0.1)
    clf.fit(X, y)

    clf2 = ElasticNet(alpha=0.1, max_iter=500)
    clf2.fit(X, y)
    assert_array_almost_equal(clf2.coef_, clf.coef_)


def test_lasso_alpha_warning():
    check_warnings()  # Skip if unsupported Python version
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        X = [[-1], [0], [1]]
        Y = [-1, 0, 1]       # just a straight line

        clf = Lasso(alpha=0)
        clf.fit(X, Y)

        assert_greater(len(w), 0)  # warnings should be raised


def test_lasso_positive_constraint():
    X = [[-1], [0], [1]]
    y = [1, 0, -1]       # just a straight line with negative slope

    lasso = Lasso(alpha=0.1, max_iter=1000, positive=True)
    lasso.fit(X, y)
    assert_true(min(lasso.coef_) >= 0)

    lasso = Lasso(alpha=0.1, max_iter=1000, precompute=True, positive=True)
    lasso.fit(X, y)
    assert_true(min(lasso.coef_) >= 0)


def test_enet_positive_constraint():
    X = [[-1], [0], [1]]
    y = [1, 0, -1]       # just a straight line with negative slope

    enet = ElasticNet(alpha=0.1, max_iter=1000, positive=True)
    enet.fit(X, y)
    assert_true(min(enet.coef_) >= 0)


def test_multi_task_lasso_and_enet():
    X, y, X_test, y_test = build_dataset()
    Y = np.c_[y, y]
    #Y_test = np.c_[y_test, y_test]
    clf = MultiTaskLasso(alpha=1, tol=1e-8).fit(X, Y)
    assert_true(0 < clf.dual_gap_ < 1e-5)
    assert_array_almost_equal(clf.coef_[0], clf.coef_[1])

    clf = MultiTaskElasticNet(alpha=1, tol=1e-8).fit(X, Y)
    assert_true(0 < clf.dual_gap_ < 1e-5)
    assert_array_almost_equal(clf.coef_[0], clf.coef_[1])


def test_enet_multitarget():
    n_targets = 3
    X, y, _, _ = build_dataset(n_samples=10, n_features=8,
                               n_informative_features=10, n_targets=n_targets)
    estimator = ElasticNet(alpha=0.01, fit_intercept=True)
    estimator.fit(X, y)
    coef, intercept, dual_gap = (estimator.coef_, estimator.intercept_,
                                 estimator.dual_gap_)

    for k in range(n_targets):
        estimator.fit(X, y[:, k])
        assert_array_almost_equal(coef[k, :], estimator.coef_)
        assert_array_almost_equal(intercept[k], estimator.intercept_)
        assert_array_almost_equal(dual_gap[k], estimator.dual_gap_)


if __name__ == '__main__':
    import nose
    nose.runmodule()
