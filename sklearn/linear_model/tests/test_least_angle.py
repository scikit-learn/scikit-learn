from nose.tools import assert_equal

import numpy as np
from scipy import linalg

from sklearn.cross_validation import train_test_split
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_no_warnings, assert_warns
from sklearn.utils.testing import TempMemmap
from sklearn.utils import ConvergenceWarning
from sklearn import linear_model, datasets
from sklearn.linear_model.least_angle import _lars_path_residues

diabetes = datasets.load_diabetes()
X, y = diabetes.data, diabetes.target

# TODO: use another dataset that has multiple drops


def test_simple():
    # Principle of Lars is to keep covariances tied and decreasing

    # also test verbose output
    from sklearn.externals.six.moves import cStringIO as StringIO
    import sys
    old_stdout = sys.stdout
    try:
        sys.stdout = StringIO()

        alphas_, active, coef_path_ = linear_model.lars_path(
            diabetes.data, diabetes.target, method="lar", verbose=10)

        sys.stdout = old_stdout

        for (i, coef_) in enumerate(coef_path_.T):
            res = y - np.dot(X, coef_)
            cov = np.dot(X.T, res)
            C = np.max(abs(cov))
            eps = 1e-3
            ocur = len(cov[C - eps < abs(cov)])
            if i < X.shape[1]:
                assert_true(ocur == i + 1)
            else:
                # no more than max_pred variables can go into the active set
                assert_true(ocur == X.shape[1])
    finally:
        sys.stdout = old_stdout


def test_simple_precomputed():
    # The same, with precomputed Gram matrix

    G = np.dot(diabetes.data.T, diabetes.data)
    alphas_, active, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, Gram=G, method="lar")

    for i, coef_ in enumerate(coef_path_.T):
        res = y - np.dot(X, coef_)
        cov = np.dot(X.T, res)
        C = np.max(abs(cov))
        eps = 1e-3
        ocur = len(cov[C - eps < abs(cov)])
        if i < X.shape[1]:
            assert_true(ocur == i + 1)
        else:
            # no more than max_pred variables can go into the active set
            assert_true(ocur == X.shape[1])


def test_all_precomputed():
    # Test that lars_path with precomputed Gram and Xy gives the right answer
    X, y = diabetes.data, diabetes.target
    G = np.dot(X.T, X)
    Xy = np.dot(X.T, y)
    for method in 'lar', 'lasso':
        output = linear_model.lars_path(X, y, method=method)
        output_pre = linear_model.lars_path(X, y, Gram=G, Xy=Xy, method=method)
        for expected, got in zip(output, output_pre):
            assert_array_almost_equal(expected, got)


def test_lars_lstsq():
    # Test that Lars gives least square solution at the end
    # of the path
    X1 = 3 * diabetes.data  # use un-normalized dataset
    clf = linear_model.LassoLars(alpha=0.)
    clf.fit(X1, y)
    coef_lstsq = np.linalg.lstsq(X1, y)[0]
    assert_array_almost_equal(clf.coef_, coef_lstsq)


def test_lasso_gives_lstsq_solution():
    # Test that Lars Lasso gives least square solution at the end
    # of the path
    alphas_, active, coef_path_ = linear_model.lars_path(X, y, method="lasso")
    coef_lstsq = np.linalg.lstsq(X, y)[0]
    assert_array_almost_equal(coef_lstsq, coef_path_[:, -1])


def test_collinearity():
    # Check that lars_path is robust to collinearity in input
    X = np.array([[3., 3., 1.],
                  [2., 2., 0.],
                  [1., 1., 0]])
    y = np.array([1., 0., 0])

    f = ignore_warnings
    _, _, coef_path_ = f(linear_model.lars_path)(X, y, alpha_min=0.01)
    assert_true(not np.isnan(coef_path_).any())
    residual = np.dot(X, coef_path_[:, -1]) - y
    assert_less((residual ** 2).sum(), 1.)  # just make sure it's bounded

    n_samples = 10
    X = np.random.rand(n_samples, 5)
    y = np.zeros(n_samples)
    _, _, coef_path_ = linear_model.lars_path(X, y, Gram='auto', copy_X=False,
                                              copy_Gram=False, alpha_min=0.,
                                              method='lasso', verbose=0,
                                              max_iter=500)
    assert_array_almost_equal(coef_path_, np.zeros_like(coef_path_))


def test_no_path():
    # Test that the ``return_path=False`` option returns the correct output

    alphas_, active_, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar")
    alpha_, active, coef = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar", return_path=False)

    assert_array_almost_equal(coef, coef_path_[:, -1])
    assert_true(alpha_ == alphas_[-1])


def test_no_path_precomputed():
    # Test that the ``return_path=False`` option with Gram remains correct

    G = np.dot(diabetes.data.T, diabetes.data)

    alphas_, active_, coef_path_ = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar", Gram=G)
    alpha_, active, coef = linear_model.lars_path(
        diabetes.data, diabetes.target, method="lar", Gram=G,
        return_path=False)

    assert_array_almost_equal(coef, coef_path_[:, -1])
    assert_true(alpha_ == alphas_[-1])


def test_no_path_all_precomputed():
    # Test that the ``return_path=False`` option with Gram and Xy remains
    # correct
    X, y = 3 * diabetes.data, diabetes.target
    G = np.dot(X.T, X)
    Xy = np.dot(X.T, y)

    alphas_, active_, coef_path_ = linear_model.lars_path(
        X, y, method="lasso", Gram=G, Xy=Xy, alpha_min=0.9)
    print("---")
    alpha_, active, coef = linear_model.lars_path(
        X, y, method="lasso", Gram=G, Xy=Xy, alpha_min=0.9, return_path=False)

    assert_array_almost_equal(coef, coef_path_[:, -1])
    assert_true(alpha_ == alphas_[-1])


def test_singular_matrix():
    # Test when input is a singular matrix
    X1 = np.array([[1, 1.], [1., 1.]])
    y1 = np.array([1, 1])
    alphas, active, coef_path = linear_model.lars_path(X1, y1)
    assert_array_almost_equal(coef_path.T, [[0, 0], [1, 0]])


def test_rank_deficient_design():
    # consistency test that checks that LARS Lasso is handling rank
    # deficient input data (with n_features < rank) in the same way
    # as coordinate descent Lasso
    y = [5, 0, 5]
    for X in ([[5, 0],
               [0, 5],
               [10, 10]],

              [[10, 10, 0],
               [1e-32, 0, 0],
               [0, 0, 1]],
              ):
        # To be able to use the coefs to compute the objective function,
        # we need to turn off normalization
        lars = linear_model.LassoLars(.1, normalize=False)
        coef_lars_ = lars.fit(X, y).coef_
        obj_lars = (1. / (2. * 3.)
                    * linalg.norm(y - np.dot(X, coef_lars_)) ** 2
                    + .1 * linalg.norm(coef_lars_, 1))
        coord_descent = linear_model.Lasso(.1, tol=1e-6, normalize=False)
        coef_cd_ = coord_descent.fit(X, y).coef_
        obj_cd = ((1. / (2. * 3.)) * linalg.norm(y - np.dot(X, coef_cd_)) ** 2
                  + .1 * linalg.norm(coef_cd_, 1))
        assert_less(obj_lars, obj_cd * (1. + 1e-8))


def test_lasso_lars_vs_lasso_cd(verbose=False):
    # Test that LassoLars and Lasso using coordinate descent give the
    # same results.
    X = 3 * diabetes.data

    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
    lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
    for c, a in zip(lasso_path.T, alphas):
        if a == 0:
            continue
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
        error = linalg.norm(c - lasso_cd.coef_)
        assert_less(error, 0.01)

    # similar test, with the classifiers
    for alpha in np.linspace(1e-2, 1 - 1e-2, 20):
        clf1 = linear_model.LassoLars(alpha=alpha, normalize=False).fit(X, y)
        clf2 = linear_model.Lasso(alpha=alpha, tol=1e-8,
                                  normalize=False).fit(X, y)
        err = linalg.norm(clf1.coef_ - clf2.coef_)
        assert_less(err, 1e-3)

    # same test, with normalized data
    X = diabetes.data
    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso')
    lasso_cd = linear_model.Lasso(fit_intercept=False, normalize=True,
                                  tol=1e-8)
    for c, a in zip(lasso_path.T, alphas):
        if a == 0:
            continue
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
        error = linalg.norm(c - lasso_cd.coef_)
        assert_less(error, 0.01)


def test_lasso_lars_vs_lasso_cd_early_stopping(verbose=False):
    # Test that LassoLars and Lasso using coordinate descent give the
    # same results when early stopping is used.
    # (test : before, in the middle, and in the last part of the path)
    alphas_min = [10, 0.9, 1e-4]
    for alphas_min in alphas_min:
        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                       alpha_min=0.9)
        lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8)
        lasso_cd.alpha = alphas[-1]
        lasso_cd.fit(X, y)
        error = linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
        assert_less(error, 0.01)

    alphas_min = [10, 0.9, 1e-4]
    # same test, with normalization
    for alphas_min in alphas_min:
        alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                       alpha_min=0.9)
        lasso_cd = linear_model.Lasso(fit_intercept=True, normalize=True,
                                      tol=1e-8)
        lasso_cd.alpha = alphas[-1]
        lasso_cd.fit(X, y)
        error = linalg.norm(lasso_path[:, -1] - lasso_cd.coef_)
        assert_less(error, 0.01)


def test_lasso_lars_path_length():
    # Test that the path length of the LassoLars is right
    lasso = linear_model.LassoLars()
    lasso.fit(X, y)
    lasso2 = linear_model.LassoLars(alpha=lasso.alphas_[2])
    lasso2.fit(X, y)
    assert_array_almost_equal(lasso.alphas_[:3], lasso2.alphas_)
    # Also check that the sequence of alphas is always decreasing
    assert_true(np.all(np.diff(lasso.alphas_) < 0))


def test_lasso_lars_vs_lasso_cd_ill_conditioned():
    # Test lasso lars on a very ill-conditioned design, and check that
    # it does not blow up, and stays somewhat close to a solution given
    # by the coordinate descent solver
    # Also test that lasso_path (using lars_path output style) gives
    # the same result as lars_path and previous lasso output style
    # under these conditions.
    rng = np.random.RandomState(42)

    # Generate data
    n, m = 70, 100
    k = 5
    X = rng.randn(n, m)
    w = np.zeros((m, 1))
    i = np.arange(0, m)
    rng.shuffle(i)
    supp = i[:k]
    w[supp] = np.sign(rng.randn(k, 1)) * (rng.rand(k, 1) + 1)
    y = np.dot(X, w)
    sigma = 0.2
    y += sigma * rng.rand(*y.shape)
    y = y.squeeze()
    lars_alphas, _, lars_coef = linear_model.lars_path(X, y, method='lasso')

    _, lasso_coef2, _ = linear_model.lasso_path(X, y,
                                                alphas=lars_alphas,
                                                tol=1e-6,
                                                fit_intercept=False)

    assert_array_almost_equal(lars_coef, lasso_coef2, decimal=1)


def test_lasso_lars_vs_lasso_cd_ill_conditioned2():
    # Create an ill-conditioned situation in which the LARS has to go
    # far in the path to converge, and check that LARS and coordinate
    # descent give the same answers
    # Note it used to be the case that Lars had to use the drop for good
    # strategy for this but this is no longer the case with the
    # equality_tolerance checks
    X = [[1e20, 1e20, 0],
         [-1e-32, 0, 0],
         [1, 1, 1]]
    y = [10, 10, 1]
    alpha = .0001

    def objective_function(coef):
        return (1. / (2. * len(X)) * linalg.norm(y - np.dot(X, coef)) ** 2
                + alpha * linalg.norm(coef, 1))

    lars = linear_model.LassoLars(alpha=alpha, normalize=False)
    assert_warns(ConvergenceWarning, lars.fit, X, y)
    lars_coef_ = lars.coef_
    lars_obj = objective_function(lars_coef_)

    coord_descent = linear_model.Lasso(alpha=alpha, tol=1e-4, normalize=False)
    cd_coef_ = coord_descent.fit(X, y).coef_
    cd_obj = objective_function(cd_coef_)

    assert_less(lars_obj, cd_obj * (1. + 1e-8))


def test_lars_add_features():
    # assure that at least some features get added if necessary
    # test for 6d2b4c
    # Hilbert matrix
    n = 5
    H = 1. / (np.arange(1, n + 1) + np.arange(n)[:, np.newaxis])
    clf = linear_model.Lars(fit_intercept=False).fit(
        H, np.arange(n))
    assert_true(np.all(np.isfinite(clf.coef_)))


def test_lars_n_nonzero_coefs(verbose=False):
    lars = linear_model.Lars(n_nonzero_coefs=6, verbose=verbose)
    lars.fit(X, y)
    assert_equal(len(lars.coef_.nonzero()[0]), 6)
    # The path should be of length 6 + 1 in a Lars going down to 6
    # non-zero coefs
    assert_equal(len(lars.alphas_), 7)


@ignore_warnings
def test_multitarget():
    # Assure that estimators receiving multidimensional y do the right thing
    X = diabetes.data
    Y = np.vstack([diabetes.target, diabetes.target ** 2]).T
    n_targets = Y.shape[1]

    for estimator in (linear_model.LassoLars(), linear_model.Lars()):
        estimator.fit(X, Y)
        Y_pred = estimator.predict(X)
        Y_dec = assert_warns(DeprecationWarning, estimator.decision_function, X)
        assert_array_almost_equal(Y_pred, Y_dec)
        alphas, active, coef, path = (estimator.alphas_, estimator.active_,
                                      estimator.coef_, estimator.coef_path_)
        for k in range(n_targets):
            estimator.fit(X, Y[:, k])
            y_pred = estimator.predict(X)
            assert_array_almost_equal(alphas[k], estimator.alphas_)
            assert_array_almost_equal(active[k], estimator.active_)
            assert_array_almost_equal(coef[k], estimator.coef_)
            assert_array_almost_equal(path[k], estimator.coef_path_)
            assert_array_almost_equal(Y_pred[:, k], y_pred)


def test_lars_cv():
    # Test the LassoLarsCV object by checking that the optimal alpha
    # increases as the number of samples increases.
    # This property is not actually garantied in general and is just a
    # property of the given dataset, with the given steps chosen.
    old_alpha = 0
    lars_cv = linear_model.LassoLarsCV()
    for length in (400, 200, 100):
        X = diabetes.data[:length]
        y = diabetes.target[:length]
        lars_cv.fit(X, y)
        np.testing.assert_array_less(old_alpha, lars_cv.alpha_)
        old_alpha = lars_cv.alpha_


def test_lasso_lars_ic():
    # Test the LassoLarsIC object by checking that
    # - some good features are selected.
    # - alpha_bic > alpha_aic
    # - n_nonzero_bic < n_nonzero_aic
    lars_bic = linear_model.LassoLarsIC('bic')
    lars_aic = linear_model.LassoLarsIC('aic')
    rng = np.random.RandomState(42)
    X = diabetes.data
    y = diabetes.target
    X = np.c_[X, rng.randn(X.shape[0], 4)]  # add 4 bad features
    lars_bic.fit(X, y)
    lars_aic.fit(X, y)
    nonzero_bic = np.where(lars_bic.coef_)[0]
    nonzero_aic = np.where(lars_aic.coef_)[0]
    assert_greater(lars_bic.alpha_, lars_aic.alpha_)
    assert_less(len(nonzero_bic), len(nonzero_aic))
    assert_less(np.max(nonzero_bic), diabetes.data.shape[1])

    # test error on unknown IC
    lars_broken = linear_model.LassoLarsIC('<unknown>')
    assert_raises(ValueError, lars_broken.fit, X, y)


def test_no_warning_for_zero_mse():
    # LassoLarsIC should not warn for log of zero MSE.
    y = np.arange(10, dtype=float)
    X = y.reshape(-1, 1)
    lars = linear_model.LassoLarsIC(normalize=False)
    assert_no_warnings(lars.fit, X, y)
    assert_true(np.any(np.isinf(lars.criterion_)))


def test_lars_path_readonly_data():
    # When using automated memory mapping on large input, the
    # fold data is in read-only mode
    # This is a non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/4597
    splitted_data = train_test_split(X, y, random_state=42)
    with TempMemmap(splitted_data) as (X_train, X_test, y_train, y_test):
        # The following should not fail despite copy=False
        _lars_path_residues(X_train, y_train, X_test, y_test, copy=False)


def test_lars_path_positive_constraint():
    # this is the main test for the positive parameter on the lars_path method
    # the estimator classes just make use of this function

    # we do the test on the diabetes dataset

    # ensure that we get negative coefficients when positive=False
    # and all positive when positive=True
    # for method 'lar' (default) and lasso
    for method in ['lar', 'lasso']:
        alpha, active, coefs = \
            linear_model.lars_path(diabetes['data'], diabetes['target'],
                                   return_path=True, method=method,
                                   positive=False)
        assert_true(coefs.min() < 0)

        alpha, active, coefs = \
            linear_model.lars_path(diabetes['data'], diabetes['target'],
                                   return_path=True, method=method,
                                   positive=True)
        assert_true(coefs.min() >= 0)


# now we gonna test the positive option for all estimator classes

default_parameter = {'fit_intercept': False}

estimator_parameter_map = {'Lars': {'n_nonzero_coefs': 5},
                           'LassoLars': {'alpha': 0.1},
                           'LarsCV': {},
                           'LassoLarsCV': {},
                           'LassoLarsIC': {}}


def test_estimatorclasses_positive_constraint():
    # testing the transmissibility for the positive option of all estimator
    # classes in this same function here

    for estname in estimator_parameter_map:
        params = default_parameter.copy()
        params.update(estimator_parameter_map[estname])
        estimator = getattr(linear_model, estname)(positive=False, **params)
        estimator.fit(diabetes['data'], diabetes['target'])
        assert_true(estimator.coef_.min() < 0)
        estimator = getattr(linear_model, estname)(positive=True, **params)
        estimator.fit(diabetes['data'], diabetes['target'])
        assert_true(min(estimator.coef_) >= 0)


def test_lasso_lars_vs_lasso_cd_positive(verbose=False):
    # Test that LassoLars and Lasso using coordinate descent give the
    # same results when using the positive option

    # This test is basically a copy of the above with additional positive
    # option. However for the middle part, the comparison of coefficient values
    # for a range of alphas, we had to make an adaptations. See below.

    # not normalized data
    X = 3 * diabetes.data

    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                   positive=True)
    lasso_cd = linear_model.Lasso(fit_intercept=False, tol=1e-8, positive=True)
    for c, a in zip(lasso_path.T, alphas):
        if a == 0:
            continue
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
        error = linalg.norm(c - lasso_cd.coef_)
        assert_less(error, 0.01)

    # The range of alphas chosen for coefficient comparison here is restricted
    # as compared with the above test without the positive option. This is due
    # to the circumstance that the Lars-Lasso algorithm does not converge to
    # the least-squares-solution for small alphas, see 'Least Angle Regression'
    # by Efron et al 2004. The coefficients are typically in congruence up to
    # the smallest alpha reached by the Lars-Lasso algorithm and start to
    # diverge thereafter.  See
    # https://gist.github.com/michigraber/7e7d7c75eca694c7a6ff

    for alpha in np.linspace(6e-1, 1 - 1e-2, 20):
        clf1 = linear_model.LassoLars(fit_intercept=False, alpha=alpha,
                                      normalize=False, positive=True).fit(X, y)
        clf2 = linear_model.Lasso(fit_intercept=False, alpha=alpha, tol=1e-8,
                                  normalize=False, positive=True).fit(X, y)
        err = linalg.norm(clf1.coef_ - clf2.coef_)
        assert_less(err, 1e-3)

    # normalized data
    X = diabetes.data
    alphas, _, lasso_path = linear_model.lars_path(X, y, method='lasso',
                                                   positive=True)
    lasso_cd = linear_model.Lasso(fit_intercept=False, normalize=True,
                                  tol=1e-8, positive=True)
    for c, a in zip(lasso_path.T[:-1], alphas[:-1]):  # don't include alpha=0
        lasso_cd.alpha = a
        lasso_cd.fit(X, y)
        error = linalg.norm(c - lasso_cd.coef_)
        assert_less(error, 0.01)
