# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause

import numpy as np

from sklearn.utils._testing import assert_allclose, assert_raises
from sklearn.datasets import make_regression
from sklearn.linear_model import HuberRegressor, QuantileRegressor


def test_quantile_toy_example():
    # test how different parameters affect a small intuitive example
    X = [[0], [1], [1]]
    y = [1, 2, 11]
    # for 50% quantile w/o regularization, any slope in [1, 10] is okay
    model = QuantileRegressor(quantile=0.5, alpha=0).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert model.coef_[0] >= 1
    assert model.coef_[0] <= 10

    # if positive error costs more, the slope is maximal
    model = QuantileRegressor(quantile=0.51, alpha=0).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert_allclose(model.coef_[0], 10, atol=1e-2)

    # if negative error costs more, the slope is minimal
    model = QuantileRegressor(quantile=0.49, alpha=0).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert_allclose(model.coef_[0], 1, atol=1e-2)

    # for a small lasso penalty, the slope is also minimal
    model = QuantileRegressor(quantile=0.5, alpha=0.01).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert_allclose(model.coef_[0], 1, atol=1e-2)

    # for a large lasso penalty, the model predicts constant median
    model = QuantileRegressor(quantile=0.5, alpha=100).fit(X, y)
    assert_allclose(model.intercept_, 2, atol=1e-2)
    assert_allclose(model.coef_[0], 0, atol=1e-2)


def test_quantile_equals_huber_for_low_epsilon():
    X, y = make_regression(n_samples=100, n_features=20, random_state=0,
                           noise=1.0)
    huber = HuberRegressor(epsilon=1+1e-4, alpha=1e-4).fit(X, y)
    quant = QuantileRegressor(alpha=1e-4).fit(X, y)
    assert_allclose(huber.intercept_, quant.intercept_, atol=1e-1)
    assert_allclose(huber.coef_, quant.coef_, atol=1e-1)


def test_quantile_estimates_fraction():
    # Test that model estimates percentage of points below the prediction
    X, y = make_regression(n_samples=1000, n_features=20, random_state=0,
                           noise=1.0)
    for q in [0.5, 0.9, 0.05]:
        quant = QuantileRegressor(quantile=q, alpha=0).fit(X, y)
        fraction_below = np.mean(y < quant.predict(X))
        assert_allclose(fraction_below, q, atol=1e-2)


def test_quantile_without_intercept():
    X, y = make_regression(n_samples=300, n_features=20, random_state=0,
                           noise=1.0)
    quant = QuantileRegressor(alpha=1e-4, fit_intercept=False).fit(X, y)
    # check that result is similar to Huber
    huber = HuberRegressor(epsilon=1 + 1e-4, alpha=1e-4, fit_intercept=False
                           ).fit(X, y)
    assert_allclose(huber.intercept_, quant.intercept_, atol=1e-1)
    assert_allclose(huber.coef_, quant.coef_, atol=1e-1)
    # check that we still predict fraction
    fraction_below = np.mean(y < quant.predict(X))
    assert_allclose(fraction_below, 0.5, atol=1e-1)


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
    assert_allclose(weighted_fraction_below, 0.5, atol=1e-2)


def test_quantile_incorrect_quantile():
    X, y = make_regression(n_samples=10, n_features=1, random_state=0, noise=1)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=2.0).fit(X, y)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=1.0).fit(X, y)
    with assert_raises(ValueError):
        QuantileRegressor(quantile=0.0).fit(X, y)
