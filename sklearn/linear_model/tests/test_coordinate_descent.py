# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

from sys import version_info

import numpy as np
from scipy import interpolate, sparse
from copy import deepcopy

from sklearn.datasets import load_boston
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import TempMemmap

from sklearn.linear_model.coordinate_descent import Lasso, \
    LassoCV, ElasticNet, ElasticNetCV, MultiTaskLasso, MultiTaskElasticNet, \
    MultiTaskElasticNetCV, MultiTaskLassoCV, lasso_path, enet_path
from sklearn.linear_model import LassoLarsCV, lars_path
from sklearn.utils import check_array


def test_lasso_zero():
    # Check that the lasso can handle zero data without crashing
    X = [[0], [0], [0]]
    y = [0, 0, 0]
    clf = Lasso(alpha=0.1).fit(X, y)
    pred = clf.predict([[1], [2], [3]])
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(pred, [0, 0, 0])
    assert_almost_equal(clf.dual_gap_, 0)


def test_lasso_toy():
    # Test Lasso on a toy example for various values of alpha.
    # When validating this against glmnet notice that glmnet divides it
    # against nobs.

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
    # Test ElasticNet for various parameters of alpha and l1_ratio.
    # Actually, the parameters alpha = 0 should not be allowed. However,
    # we test it as a border case.
    # ElasticNet is tested with and without precomputed Gram matrix

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
        np.searchsorted(clf.alphas_[::-1], lars.alpha_) -
        np.searchsorted(clf.alphas_[::-1], clf.alpha_)) <= 1)
    # check that they also give a similar MSE
    mse_lars = interpolate.interp1d(lars.cv_alphas_, lars.mse_path_.T)
    np.testing.assert_approx_equal(mse_lars(clf.alphas_[5]).mean(),
                                   clf.mse_path_[5].mean(), significant=2)

    # test set
    assert_greater(clf.score(X_test, y_test), 0.99)


def test_lasso_cv_positive_constraint():
    X, y, X_test, y_test = build_dataset()
    max_iter = 500

    # Ensure the unconstrained fit has a negative coefficient
    clf_unconstrained = LassoCV(n_alphas=3, eps=1e-1, max_iter=max_iter, cv=2,
                                n_jobs=1)
    clf_unconstrained.fit(X, y)
    assert_true(min(clf_unconstrained.coef_) < 0)

    # On same data, constrained fit has non-negative coefficients
    clf_constrained = LassoCV(n_alphas=3, eps=1e-1, max_iter=max_iter,
                              positive=True, cv=2, n_jobs=1)
    clf_constrained.fit(X, y)
    assert_true(min(clf_constrained.coef_) >= 0)


def test_lasso_path_return_models_vs_new_return_gives_same_coefficients():
    # Test that lasso_path with lars_path style output gives the
    # same result

    # Some toy data
    X = np.array([[1, 2, 3.1], [2.3, 5.4, 4.3]]).T
    y = np.array([1, 2, 3.1])
    alphas = [5., 1., .5]

    # Use lars_path and lasso_path(new output) with 1D linear interpolation
    # to compute the same path
    alphas_lars, _, coef_path_lars = lars_path(X, y, method='lasso')
    coef_path_cont_lars = interpolate.interp1d(alphas_lars[::-1],
                                               coef_path_lars[:, ::-1])
    alphas_lasso2, coef_path_lasso2, _ = lasso_path(X, y, alphas=alphas,
                                                    return_models=False)
    coef_path_cont_lasso = interpolate.interp1d(alphas_lasso2[::-1],
                                                coef_path_lasso2[:, ::-1])

    assert_array_almost_equal(
        coef_path_cont_lasso(alphas), coef_path_cont_lars(alphas),
        decimal=1)


def test_enet_path():
    # We use a large number of samples and of informative features so that
    # the l1_ratio selected is more toward ridge than lasso
    X, y, X_test, y_test = build_dataset(n_samples=200, n_features=100,
                                         n_informative_features=100)
    max_iter = 150

    # Here we have a small number of iterations, and thus the
    # ElasticNet might not converge. This is to speed up tests
    clf = ElasticNetCV(alphas=[0.01, 0.05, 0.1], eps=2e-3,
                       l1_ratio=[0.5, 0.7], cv=3,
                       max_iter=max_iter)
    ignore_warnings(clf.fit)(X, y)
    # Well-conditioned settings, we should have selected our
    # smallest penalty
    assert_almost_equal(clf.alpha_, min(clf.alphas_))
    # Non-sparse ground truth: we should have selected an elastic-net
    # that is closer to ridge than to lasso
    assert_equal(clf.l1_ratio_, min(clf.l1_ratio))

    clf = ElasticNetCV(alphas=[0.01, 0.05, 0.1], eps=2e-3,
                       l1_ratio=[0.5, 0.7], cv=3,
                       max_iter=max_iter, precompute=True)
    ignore_warnings(clf.fit)(X, y)

    # Well-conditioned settings, we should have selected our
    # smallest penalty
    assert_almost_equal(clf.alpha_, min(clf.alphas_))
    # Non-sparse ground truth: we should have selected an elastic-net
    # that is closer to ridge than to lasso
    assert_equal(clf.l1_ratio_, min(clf.l1_ratio))

    # We are in well-conditioned settings with low noise: we should
    # have a good test-set performance
    assert_greater(clf.score(X_test, y_test), 0.99)

    # Multi-output/target case
    X, y, X_test, y_test = build_dataset(n_features=10, n_targets=3)
    clf = MultiTaskElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7],
                                cv=3, max_iter=max_iter)
    ignore_warnings(clf.fit)(X, y)
    # We are in well-conditioned settings with low noise: we should
    # have a good test-set performance
    assert_greater(clf.score(X_test, y_test), 0.99)
    assert_equal(clf.coef_.shape, (3, 10))

    # Mono-output should have same cross-validated alpha_ and l1_ratio_
    # in both cases.
    X, y, _, _ = build_dataset(n_features=10)
    clf1 = ElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7])
    clf1.fit(X, y)
    clf2 = MultiTaskElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7])
    clf2.fit(X, y[:, np.newaxis])
    assert_almost_equal(clf1.l1_ratio_, clf2.l1_ratio_)
    assert_almost_equal(clf1.alpha_, clf2.alpha_)


def test_path_parameters():
    X, y, _, _ = build_dataset()
    max_iter = 100

    clf = ElasticNetCV(n_alphas=50, eps=1e-3, max_iter=max_iter,
                       l1_ratio=0.5, tol=1e-3)
    clf.fit(X, y)  # new params
    assert_almost_equal(0.5, clf.l1_ratio)
    assert_equal(50, clf.n_alphas)
    assert_equal(50, len(clf.alphas_))


def test_warm_start():
    X, y, _, _ = build_dataset()
    clf = ElasticNet(alpha=0.1, max_iter=5, warm_start=True)
    ignore_warnings(clf.fit)(X, y)
    ignore_warnings(clf.fit)(X, y)  # do a second round with 5 iterations

    clf2 = ElasticNet(alpha=0.1, max_iter=10)
    ignore_warnings(clf2.fit)(X, y)
    assert_array_almost_equal(clf2.coef_, clf.coef_)


def test_lasso_alpha_warning():
    X = [[-1], [0], [1]]
    Y = [-1, 0, 1]       # just a straight line

    clf = Lasso(alpha=0)
    assert_warns(UserWarning, clf.fit, X, Y)


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


def test_enet_cv_positive_constraint():
    X, y, X_test, y_test = build_dataset()
    max_iter = 500

    # Ensure the unconstrained fit has a negative coefficient
    enetcv_unconstrained = ElasticNetCV(n_alphas=3, eps=1e-1,
                                        max_iter=max_iter,
                                        cv=2, n_jobs=1)
    enetcv_unconstrained.fit(X, y)
    assert_true(min(enetcv_unconstrained.coef_) < 0)

    # On same data, constrained fit has non-negative coefficients
    enetcv_constrained = ElasticNetCV(n_alphas=3, eps=1e-1, max_iter=max_iter,
                                      cv=2, positive=True, n_jobs=1)
    enetcv_constrained.fit(X, y)
    assert_true(min(enetcv_constrained.coef_) >= 0)


def test_uniform_targets():
    enet = ElasticNetCV(fit_intercept=True, n_alphas=3)
    m_enet = MultiTaskElasticNetCV(fit_intercept=True, n_alphas=3)
    lasso = LassoCV(fit_intercept=True, n_alphas=3)
    m_lasso = MultiTaskLassoCV(fit_intercept=True, n_alphas=3)

    models_single_task = (enet, lasso)
    models_multi_task = (m_enet, m_lasso)

    rng = np.random.RandomState(0)

    X_train = rng.random_sample(size=(10, 3))
    X_test = rng.random_sample(size=(10, 3))

    y1 = np.empty(10)
    y2 = np.empty((10, 2))

    for model in models_single_task:
        for y_values in (0, 5):
            y1.fill(y_values)
            assert_array_equal(model.fit(X_train, y1).predict(X_test), y1)
            assert_array_equal(model.alphas_, [np.finfo(float).resolution]*3)

    for model in models_multi_task:
        for y_values in (0, 5):
            y2[:, 0].fill(y_values)
            y2[:, 1].fill(2 * y_values)
            assert_array_equal(model.fit(X_train, y2).predict(X_test), y2)
            assert_array_equal(model.alphas_, [np.finfo(float).resolution]*3)


def test_multi_task_lasso_and_enet():
    X, y, X_test, y_test = build_dataset()
    Y = np.c_[y, y]
    # Y_test = np.c_[y_test, y_test]
    clf = MultiTaskLasso(alpha=1, tol=1e-8).fit(X, Y)
    assert_true(0 < clf.dual_gap_ < 1e-5)
    assert_array_almost_equal(clf.coef_[0], clf.coef_[1])

    clf = MultiTaskElasticNet(alpha=1, tol=1e-8).fit(X, Y)
    assert_true(0 < clf.dual_gap_ < 1e-5)
    assert_array_almost_equal(clf.coef_[0], clf.coef_[1])

    clf = MultiTaskElasticNet(alpha=1.0, tol=1e-8, max_iter=1)
    assert_warns_message(ConvergenceWarning, 'did not converge', clf.fit, X, Y)


def test_lasso_readonly_data():
    X = np.array([[-1], [0], [1]])
    Y = np.array([-1, 0, 1])   # just a straight line
    T = np.array([[2], [3], [4]])  # test sample
    with TempMemmap((X, Y)) as (X, Y):
        clf = Lasso(alpha=0.5)
        clf.fit(X, Y)
        pred = clf.predict(T)
        assert_array_almost_equal(clf.coef_, [.25])
        assert_array_almost_equal(pred, [0.5, 0.75, 1.])
        assert_almost_equal(clf.dual_gap_, 0)


def test_multi_task_lasso_readonly_data():
    X, y, X_test, y_test = build_dataset()
    Y = np.c_[y, y]
    with TempMemmap((X, Y)) as (X, Y):
        Y = np.c_[y, y]
        clf = MultiTaskLasso(alpha=1, tol=1e-8).fit(X, Y)
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


def test_multioutput_enetcv_error():
    rng = np.random.RandomState(0)
    X = rng.randn(10, 2)
    y = rng.randn(10, 2)
    clf = ElasticNetCV()
    assert_raises(ValueError, clf.fit, X, y)


def test_multitask_enet_and_lasso_cv():
    X, y, _, _ = build_dataset(n_features=50, n_targets=3)
    clf = MultiTaskElasticNetCV().fit(X, y)
    assert_almost_equal(clf.alpha_, 0.00556, 3)
    clf = MultiTaskLassoCV().fit(X, y)
    assert_almost_equal(clf.alpha_, 0.00278, 3)

    X, y, _, _ = build_dataset(n_targets=3)
    clf = MultiTaskElasticNetCV(n_alphas=10, eps=1e-3, max_iter=100,
                                l1_ratio=[0.3, 0.5], tol=1e-3)
    clf.fit(X, y)
    assert_equal(0.5, clf.l1_ratio_)
    assert_equal((3, X.shape[1]), clf.coef_.shape)
    assert_equal((3, ), clf.intercept_.shape)
    assert_equal((2, 10, 3), clf.mse_path_.shape)
    assert_equal((2, 10), clf.alphas_.shape)

    X, y, _, _ = build_dataset(n_targets=3)
    clf = MultiTaskLassoCV(n_alphas=10, eps=1e-3, max_iter=100, tol=1e-3)
    clf.fit(X, y)
    assert_equal((3, X.shape[1]), clf.coef_.shape)
    assert_equal((3, ), clf.intercept_.shape)
    assert_equal((10, 3), clf.mse_path_.shape)
    assert_equal(10, len(clf.alphas_))


def test_1d_multioutput_enet_and_multitask_enet_cv():
    X, y, _, _ = build_dataset(n_features=10)
    y = y[:, np.newaxis]
    clf = ElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7])
    clf.fit(X, y[:, 0])
    clf1 = MultiTaskElasticNetCV(n_alphas=5, eps=2e-3, l1_ratio=[0.5, 0.7])
    clf1.fit(X, y)
    assert_almost_equal(clf.l1_ratio_, clf1.l1_ratio_)
    assert_almost_equal(clf.alpha_, clf1.alpha_)
    assert_almost_equal(clf.coef_, clf1.coef_[0])
    assert_almost_equal(clf.intercept_, clf1.intercept_[0])


def test_1d_multioutput_lasso_and_multitask_lasso_cv():
    X, y, _, _ = build_dataset(n_features=10)
    y = y[:, np.newaxis]
    clf = LassoCV(n_alphas=5, eps=2e-3)
    clf.fit(X, y[:, 0])
    clf1 = MultiTaskLassoCV(n_alphas=5, eps=2e-3)
    clf1.fit(X, y)
    assert_almost_equal(clf.alpha_, clf1.alpha_)
    assert_almost_equal(clf.coef_, clf1.coef_[0])
    assert_almost_equal(clf.intercept_, clf1.intercept_[0])


def test_sparse_input_dtype_enet_and_lassocv():
    X, y, _, _ = build_dataset(n_features=10)
    clf = ElasticNetCV(n_alphas=5)
    clf.fit(sparse.csr_matrix(X), y)
    clf1 = ElasticNetCV(n_alphas=5)
    clf1.fit(sparse.csr_matrix(X, dtype=np.float32), y)
    assert_almost_equal(clf.alpha_, clf1.alpha_, decimal=6)
    assert_almost_equal(clf.coef_, clf1.coef_, decimal=6)

    clf = LassoCV(n_alphas=5)
    clf.fit(sparse.csr_matrix(X), y)
    clf1 = LassoCV(n_alphas=5)
    clf1.fit(sparse.csr_matrix(X, dtype=np.float32), y)
    assert_almost_equal(clf.alpha_, clf1.alpha_, decimal=6)
    assert_almost_equal(clf.coef_, clf1.coef_, decimal=6)


def test_precompute_invalid_argument():
    X, y, _, _ = build_dataset()
    for clf in [ElasticNetCV(precompute="invalid"),
                LassoCV(precompute="invalid")]:
        assert_raises_regex(ValueError, ".*should be.*True.*False.*auto.*"
                            "array-like.*Got 'invalid'", clf.fit, X, y)

    # Precompute = 'auto' is not supported for ElasticNet
    assert_raises_regex(ValueError, ".*should be.*True.*False.*array-like.*"
                        "Got 'auto'", ElasticNet(precompute='auto').fit, X, y)


def test_warm_start_convergence():
    X, y, _, _ = build_dataset()
    model = ElasticNet(alpha=1e-3, tol=1e-3).fit(X, y)
    n_iter_reference = model.n_iter_

    # This dataset is not trivial enough for the model to converge in one pass.
    assert_greater(n_iter_reference, 2)

    # Check that n_iter_ is invariant to multiple calls to fit
    # when warm_start=False, all else being equal.
    model.fit(X, y)
    n_iter_cold_start = model.n_iter_
    assert_equal(n_iter_cold_start, n_iter_reference)

    # Fit the same model again, using a warm start: the optimizer just performs
    # a single pass before checking that it has already converged
    model.set_params(warm_start=True)
    model.fit(X, y)
    n_iter_warm_start = model.n_iter_
    assert_equal(n_iter_warm_start, 1)


def test_warm_start_convergence_with_regularizer_decrement():
    boston = load_boston()
    X, y = boston.data, boston.target

    # Train a model to converge on a lightly regularized problem
    final_alpha = 1e-5
    low_reg_model = ElasticNet(alpha=final_alpha).fit(X, y)

    # Fitting a new model on a more regularized version of the same problem.
    # Fitting with high regularization is easier it should converge faster
    # in general.
    high_reg_model = ElasticNet(alpha=final_alpha * 10).fit(X, y)
    assert_greater(low_reg_model.n_iter_, high_reg_model.n_iter_)

    # Fit the solution to the original, less regularized version of the
    # problem but from the solution of the highly regularized variant of
    # the problem as a better starting point. This should also converge
    # faster than the original model that starts from zero.
    warm_low_reg_model = deepcopy(high_reg_model)
    warm_low_reg_model.set_params(warm_start=True, alpha=final_alpha)
    warm_low_reg_model.fit(X, y)
    assert_greater(low_reg_model.n_iter_, warm_low_reg_model.n_iter_)


def test_random_descent():
    # Test that both random and cyclic selection give the same results.
    # Ensure that the test models fully converge and check a wide
    # range of conditions.

    # This uses the coordinate descent algo using the gram trick.
    X, y, _, _ = build_dataset(n_samples=50, n_features=20)
    clf_cyclic = ElasticNet(selection='cyclic', tol=1e-8)
    clf_cyclic.fit(X, y)
    clf_random = ElasticNet(selection='random', tol=1e-8, random_state=42)
    clf_random.fit(X, y)
    assert_array_almost_equal(clf_cyclic.coef_, clf_random.coef_)
    assert_almost_equal(clf_cyclic.intercept_, clf_random.intercept_)

    # This uses the descent algo without the gram trick
    clf_cyclic = ElasticNet(selection='cyclic', tol=1e-8)
    clf_cyclic.fit(X.T, y[:20])
    clf_random = ElasticNet(selection='random', tol=1e-8, random_state=42)
    clf_random.fit(X.T, y[:20])
    assert_array_almost_equal(clf_cyclic.coef_, clf_random.coef_)
    assert_almost_equal(clf_cyclic.intercept_, clf_random.intercept_)

    # Sparse Case
    clf_cyclic = ElasticNet(selection='cyclic', tol=1e-8)
    clf_cyclic.fit(sparse.csr_matrix(X), y)
    clf_random = ElasticNet(selection='random', tol=1e-8, random_state=42)
    clf_random.fit(sparse.csr_matrix(X), y)
    assert_array_almost_equal(clf_cyclic.coef_, clf_random.coef_)
    assert_almost_equal(clf_cyclic.intercept_, clf_random.intercept_)

    # Multioutput case.
    new_y = np.hstack((y[:, np.newaxis], y[:, np.newaxis]))
    clf_cyclic = MultiTaskElasticNet(selection='cyclic', tol=1e-8)
    clf_cyclic.fit(X, new_y)
    clf_random = MultiTaskElasticNet(selection='random', tol=1e-8,
                                     random_state=42)
    clf_random.fit(X, new_y)
    assert_array_almost_equal(clf_cyclic.coef_, clf_random.coef_)
    assert_almost_equal(clf_cyclic.intercept_, clf_random.intercept_)

    # Raise error when selection is not in cyclic or random.
    clf_random = ElasticNet(selection='invalid')
    assert_raises(ValueError, clf_random.fit, X, y)


def test_enet_path_positive():
    # Test that the coefs returned by positive=True in enet_path are positive

    X, y, _, _ = build_dataset(n_samples=50, n_features=50)
    for path in [enet_path, lasso_path]:
        pos_path_coef = path(X, y, positive=True)[1]
        assert_true(np.all(pos_path_coef >= 0))


def test_sparse_dense_descent_paths():
    # Test that dense and sparse input give the same input for descent paths.
    X, y, _, _ = build_dataset(n_samples=50, n_features=20)
    csr = sparse.csr_matrix(X)
    for path in [enet_path, lasso_path]:
        _, coefs, _ = path(X, y, fit_intercept=False)
        _, sparse_coefs, _ = path(csr, y, fit_intercept=False)
        assert_array_almost_equal(coefs, sparse_coefs)


def test_check_input_false():
    X, y, _, _ = build_dataset(n_samples=20, n_features=10)
    X = check_array(X, order='F', dtype='float64')
    y = check_array(X, order='F', dtype='float64')
    clf = ElasticNet(selection='cyclic', tol=1e-8)
    # Check that no error is raised if data is provided in the right format
    clf.fit(X, y, check_input=False)
    X = check_array(X, order='F', dtype='float32')
    clf.fit(X, y, check_input=True)
    # Check that an error is raised if data is provided in the wrong dtype,
    # because of check bypassing
    assert_raises(ValueError, clf.fit, X, y, check_input=False)

    # With no input checking, providing X in C order should result in false
    # computation
    X = check_array(X, order='C', dtype='float64')
    assert_raises(ValueError, clf.fit, X, y, check_input=False)


def test_overrided_gram_matrix():
    X, y, _, _ = build_dataset(n_samples=20, n_features=10)
    Gram = X.T.dot(X)
    clf = ElasticNet(selection='cyclic', tol=1e-8, precompute=Gram,
                     fit_intercept=True)
    assert_warns_message(UserWarning,
                         "Gram matrix was provided but X was centered"
                         " to fit intercept, "
                         "or X was normalized : recomputing Gram matrix.",
                         clf.fit, X, y)


def test_lasso_non_float_y():
    X = [[0, 0], [1, 1], [-1, -1]]
    y = [0, 1, 2]
    y_float = [0.0, 1.0, 2.0]

    for model in [ElasticNet, Lasso]:
        clf = model(fit_intercept=False)
        clf.fit(X, y)
        clf_float = model(fit_intercept=False)
        clf_float.fit(X, y_float)
        assert_array_equal(clf.coef_, clf_float.coef_)


def test_enet_float_precision():
    # Generate dataset
    X, y, X_test, y_test = build_dataset(n_samples=20, n_features=10)
    # Here we have a small number of iterations, and thus the
    # ElasticNet might not converge. This is to speed up tests

    for normalize in [True, False]:
        for fit_intercept in [True, False]:
            coef = {}
            intercept = {}
            for dtype in [np.float64, np.float32]:
                clf = ElasticNet(alpha=0.5, max_iter=100, precompute=False,
                                 fit_intercept=fit_intercept,
                                 normalize=normalize)

                X = dtype(X)
                y = dtype(y)
                ignore_warnings(clf.fit)(X, y)

                coef[('simple', dtype)] = clf.coef_
                intercept[('simple', dtype)] = clf.intercept_

                assert_equal(clf.coef_.dtype, dtype)

                # test precompute Gram array
                Gram = X.T.dot(X)
                clf_precompute = ElasticNet(alpha=0.5, max_iter=100,
                                            precompute=Gram,
                                            fit_intercept=fit_intercept,
                                            normalize=normalize)
                ignore_warnings(clf_precompute.fit)(X, y)
                assert_array_almost_equal(clf.coef_, clf_precompute.coef_)
                assert_array_almost_equal(clf.intercept_,
                                          clf_precompute.intercept_)

                # test multi task enet
                multi_y = np.hstack((y[:, np.newaxis], y[:, np.newaxis]))
                clf_multioutput = MultiTaskElasticNet(
                    alpha=0.5, max_iter=100, fit_intercept=fit_intercept,
                    normalize=normalize)
                clf_multioutput.fit(X, multi_y)
                coef[('multi', dtype)] = clf_multioutput.coef_
                intercept[('multi', dtype)] = clf_multioutput.intercept_
                assert_equal(clf.coef_.dtype, dtype)

            for v in ['simple', 'multi']:
                assert_array_almost_equal(coef[(v, np.float32)],
                                          coef[(v, np.float64)],
                                          decimal=4)
                assert_array_almost_equal(intercept[(v, np.float32)],
                                          intercept[(v, np.float64)],
                                          decimal=4)


def test_enet_l1_ratio():
    # Test that an error message is raised if an estimator that
    # uses _alpha_grid is called with l1_ratio=0
    msg = ("Automatic alpha grid generation is not supported for l1_ratio=0. "
           "Please supply a grid by providing your estimator with the "
           "appropriate `alphas=` argument.")
    X = np.array([[1, 2, 4, 5, 8], [3, 5, 7, 7, 8]]).T
    y = np.array([12, 10, 11, 21, 5])

    assert_raise_message(ValueError, msg, ElasticNetCV(
        l1_ratio=0, random_state=42).fit, X, y)
    assert_raise_message(ValueError, msg, MultiTaskElasticNetCV(
        l1_ratio=0, random_state=42).fit, X, y[:, None])

    # Test that l1_ratio=0 is allowed if we supply a grid manually
    alphas = [0.1, 10]
    estkwds = {'alphas': alphas, 'random_state': 42}
    est_desired = ElasticNetCV(l1_ratio=0.00001, **estkwds)
    est = ElasticNetCV(l1_ratio=0, **estkwds)
    with ignore_warnings():
        est_desired.fit(X, y)
        est.fit(X, y)
    assert_array_almost_equal(est.coef_, est_desired.coef_, decimal=5)

    est_desired = MultiTaskElasticNetCV(l1_ratio=0.00001, **estkwds)
    est = MultiTaskElasticNetCV(l1_ratio=0, **estkwds)
    with ignore_warnings():
        est.fit(X, y[:, None])
        est_desired.fit(X, y[:, None])
    assert_array_almost_equal(est.coef_, est_desired.coef_, decimal=5)
