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

def check_warnings():
    if version_info < (2, 6):
        raise SkipTest("Testing for warnings is not supported in versions \
        older than Python 2.6")

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
    check_warnings()  # Skip if unsupported Python version
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

#def test_simple():
#    """
#    Principle of Lars is to keep covariances tied and decreasing
#    """
#
#    alphas_, active, coef_path_ = linear_model.lars_path(
#        diabetes.data, diabetes.target, method="lar")
#
#    for (i, coef_) in enumerate(coef_path_.T):
#        res = y - np.dot(X, coef_)
#        cov = np.dot(X.T, res)
#        C = np.max(abs(cov))
#        eps = 1e-3
#        ocur = len(cov[C - eps < abs(cov)])
#        if i < X.shape[1]:
#            assert_true(ocur == i + 1)
#        else:
#            # no more than max_pred variables can go into the active set
#            assert_true(ocur == X.shape[1])

#def test_collinearity():
#    """Check that lars_path is robust to collinearity in input"""
#    X = np.array([[3., 3., 1.],
#                  [2., 2., 0.],
#                  [1., 1., 0]])
#    y = np.array([1., 0., 0])

#    _, _, coef_path_ = linear_model.lars_path(X, y, alpha_min=0.01)
#    assert_true(not np.isnan(coef_path_).any())
#    residual = np.dot(X, coef_path_[:, -1]) - y
#    assert_less((residual ** 2).sum(), 1.)  # just make sure it's bounded

#    n_samples = 10
#    X = np.random.rand(n_samples, 5)
#    y = np.zeros(n_samples)
#    _, _, coef_path_ = linear_model.lars_path(X, y, Gram='auto', copy_X=False,
#            copy_Gram=False, alpha_min=0., method='lasso', verbose=0,
#            max_iter=500)
#    assert_array_almost_equal(coef_path_, np.zeros_like(coef_path_))


#def test_lasso_lars_vs_lasso_cd(verbose=False):
#    """
#    Test that LassoLars and Lasso using coordinate descent give the
#    same results
#    """
#    X = 3 * diabetes.data

#    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
#    lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
#    for c, a in zip(lasso_path.T, alphas):
#        if a == 0:
#            continue
#        lasso_cd.alpha = a
#        lasso_cd.fit(X, y)
#        error = np.linalg.norm(c - lasso_cd.coef_)
#        assert_less(error, 0.01)

    # similar test, with the classifiers
#    for alpha in np.linspace(1e-2, 1 - 1e-2):
#        clf1 = linear_model.LassoLars(alpha=alpha, normalize=False).fit(X, y)
#        clf2 = linear_model.Lasso(alpha=alpha, tol=1e-8,
#                                  normalize=False).fit(X, y)
#        err = np.linalg.norm(clf1.coef_ - clf2.coef_)
#        assert_less(err, 1e-3)

    # same test, with normalized data
#    X = diabetes.data
#    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
#    lasso_cd = linear_model.Lasso(fit_intercept=False, normalize=True,
#                                  tol=1e-8)
#    for c, a in zip(lasso_path.T, alphas):
#        if a == 0:
#            continue
#        lasso_cd.alpha = a
#        lasso_cd.fit(X, y)
#        error = np.linalg.norm(c - lasso_cd.coef_)
#        assert_less(error, 0.01)


#def test_lasso_lars_vs_lasso_cd_early_stopping(verbose=False):
#    """
#    Test that LassoLars and Lasso using coordinate descent give the
#    same results when early stopping is used.
#    (test : before, in the middle, and in the last part of the path)
#    """
#    alphas_min = [10, 0.9, 1e-4]
#    for alphas_min in alphas_min:
#        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
#                                                    alpha_min=0.9)
#        lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
#        lasso_cd.alpha = alphas[-1]
#        lasso_cd.fit(X, y)
#        error = np.linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
#        assert_less(error, 0.01)

#    alphas_min = [10, 0.9, 1e-4]
#    # same test, with normalization
#    for alphas_min in alphas_min:
#        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
#                                                    alpha_min=0.9)
#        lasso_cd = linear_model.Lasso(fit_intercept=True, normalize=True,
#                                      tol=1e-8)
#        lasso_cd.alpha = alphas[-1]
#        lasso_cd.fit(X, y)
#        error = np.linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
#        assert_less(error, 0.01)



#def test_multitarget():
#    """
#    Assure that estimators receiving multidimensional y do the right thing
#    """
#    X = diabetes.data
#    Y = np.vstack([diabetes.target, diabetes.target ** 2]).T
#    n_targets = Y.shape[1]

#    for estimator in (linear_model.LassoLars(), linear_model.Lars()):
#        estimator.fit(X, Y)
#        alphas, active, coef, path = (estimator.alphas_, estimator.active_,
#                                      estimator.coef_, estimator.coef_path_)
#        for k in xrange(n_targets):
#            estimator.fit(X, Y[:, k])
#            assert_array_almost_equal(alphas[k, :], estimator.alphas_)
#            assert_array_almost_equal(active[k, :], estimator.active_)
#            assert_array_almost_equal(coef[k, :], estimator.coef_)
#            assert_array_almost_equal(path[k, :, :], estimator.coef_path_)



if __name__ == '__main__':
    import nose
    nose.runmodule()
