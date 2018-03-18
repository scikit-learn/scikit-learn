# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause

import numpy as np
from sklearn.utils.testing import assert_almost_equal, \
    assert_raises, assert_array_almost_equal
from sklearn.datasets import make_regression
from sklearn.linear_model import HuberRegressor, QuantileRegressor

import warnings


def test_quantile_equals_huber_for_low_epsilon():
    X, y = make_regression(n_samples=100, n_features=20, random_state=0,
                           noise=1.0)
    huber = HuberRegressor(epsilon=1+1e-4, alpha=1e-4).fit(X, y)
    quant = QuantileRegressor(alpha=1e-4).fit(X, y)
    assert_almost_equal(huber.intercept_, quant.intercept_, 1)
    assert_almost_equal(huber.coef_, quant.coef_, 1)


def test_quantile_estimates_fraction():
    # Test that model estimates percentage of points below the prediction
    X, y = make_regression(n_samples=1000, n_features=20, random_state=0,
                           noise=1.0)
    for q in [0.5, 0.9, 0.05]:
        quant = QuantileRegressor(quantile=q, alpha=0).fit(X, y)
        fraction_below = np.mean(y < quant.predict(X))
        assert_almost_equal(fraction_below, q, 2)


def test_quantile_is_approximately_sparse():
    # Now most of coefficients are not exact zero,
    # but with large n_samples they are close enough
    X, y = make_regression(n_samples=3000, n_features=100, n_informative=10,
                           random_state=0, noise=1.0)
    q = QuantileRegressor(l1_ratio=1, alpha=0.1).fit(X, y)
    share_zeros = np.mean(np.abs(q.coef_) > 1e-1)
    assert_almost_equal(share_zeros, 0.1, 2)


def test_quantile_without_intercept():
    X, y = make_regression(n_samples=300, n_features=20, random_state=0,
                           noise=1.0)
    quant = QuantileRegressor(alpha=1e-4, fit_intercept=False).fit(X, y)
    # check that result is similar to Huber
    huber = HuberRegressor(epsilon=1 + 1e-4, alpha=1e-4, fit_intercept=False
                           ).fit(X, y)
    assert_almost_equal(huber.intercept_, quant.intercept_, 1)
    assert_almost_equal(huber.coef_, quant.coef_, 1)
    # check that we still predict fraction
    fraction_below = np.mean(y < quant.predict(X))
    assert_almost_equal(fraction_below, 0.5, 1)


def test_quantile_sample_weight():
    # test that with unequal sample weights we still estimate weighted fraction
    n = 1000
    X, y = make_regression(n_samples=n, n_features=10, random_state=0,
                           noise=10.0)
    weight = np.ones(n)
    # when we increase weight of upper observaions,
    # estimate of quantile should go up
    weight[y > y.mean()] = 100
    quant = QuantileRegressor(quantile=0.5, alpha=1e-4)
    quant.fit(X, y, sample_weight=weight)
    fraction_below = np.mean(y < quant.predict(X))
    assert fraction_below > 0.5
    weighted_fraction_below = np.sum((y < quant.predict(X)) * weight) \
        / np.sum(weight)
    assert_almost_equal(weighted_fraction_below, 0.5, 2)


def test_quantile_incorrect_quantile():
    X, y = make_regression(n_samples=10, n_features=1, random_state=0, noise=1)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=2.0).fit(X, y)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=1.0).fit(X, y)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=0.0).fit(X, y)


def test_quantile_warm_start():
    X, y = make_regression()
    warm = QuantileRegressor(fit_intercept=True, alpha=1.0, max_iter=10000,
                             warm_start=True, gtol=1e-1)
    warm.fit(X, y)
    warm_coef = warm.coef_.copy()
    warm.fit(X, y)

    # SciPy performs the tol check after doing the coef updates, so
    # these would be almost same but not equal.
    assert_array_almost_equal(warm.coef_, warm_coef, 1)


def test_quantile_convergence():
    # Quantile loss may not converge to unique solution
    # if there is no regularization
    # need to check that warning is not thrown if model has converged.
    X, y = make_regression(n_samples=300, n_features=20, random_state=0,
                           noise=1.0)

    # check that for small n_iter, warning is thrown
    with warnings.catch_warnings(record=True) as w:
        QuantileRegressor(max_iter=10).fit(X, y)
        assert len(w) == 1
        assert 'QuantileRegressor convergence failed' in str(w[-1].message)

    # check that for large n_iter, it is not thrown
    with warnings.catch_warnings(record=True) as w:
        QuantileRegressor(max_iter=10000).fit(X, y)
        assert len(w) == 0
