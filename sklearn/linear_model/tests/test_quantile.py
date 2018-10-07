# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause

import pytest
import numpy as np
from sklearn.utils.testing import assert_allclose, assert_raises
from sklearn.datasets import make_regression
from sklearn.linear_model import HuberRegressor, QuantileRegressor
from sklearn.model_selection import cross_val_score


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

    # for a small ridge penalty, the slope is also minimal
    model = QuantileRegressor(quantile=0.5, alpha=0.01).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert_allclose(model.coef_[0], 1, atol=1e-2)

    # for a small lasso penalty, the slope is also minimal
    model = QuantileRegressor(quantile=0.5, alpha=0.01, l1_ratio=1).fit(X, y)
    assert_allclose(model.intercept_, 1, atol=1e-2)
    assert_allclose(model.coef_[0], 1, atol=1e-2)

    # for a large ridge penalty, the model no longer minimizes MAE
    # (1.75, 0.25) minimizes c^2 + 0.5 (abs(1-b) + abs(2-b-c) + abs(11-b-c))
    model = QuantileRegressor(quantile=0.5, alpha=1).fit(X, y)
    assert_allclose(model.intercept_, 1.75, atol=1e-2)
    assert_allclose(model.coef_[0], 0.25, atol=1e-2)


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


def test_quantile_is_approximately_sparse():
    # Now most of coefficients are not exact zero,
    # but with large n_samples they are close enough
    X, y = make_regression(n_samples=3000, n_features=100, n_informative=10,
                           random_state=0, noise=1.0)
    q = QuantileRegressor(l1_ratio=1, alpha=0.1).fit(X, y)
    share_zeros = np.mean(np.abs(q.coef_) > 1e-1)
    assert_allclose(share_zeros, 0.1, atol=1e-2)


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


def test_normalize():
    # test that normalization works ok if features have different scales
    X, y = make_regression(n_samples=1000, n_features=20, random_state=0,
                           noise=10.0)
    rng = np.random.RandomState(0)
    X += rng.normal(size=X.shape[1], scale=3)
    X *= rng.normal(size=X.shape[1], scale=3)
    y = y * 10 + 100
    model1 = QuantileRegressor(alpha=1e-6, normalize=False, max_iter=10000)
    model2 = QuantileRegressor(alpha=1e-6, normalize=True, max_iter=10000)
    cvs1 = cross_val_score(model1, X, y, cv=3).mean()
    cvs2 = cross_val_score(model2, X, y, cv=3).mean()
    assert cvs1 > 0.99
    assert cvs2 > 0.99


def test_quantile_warm_start():
    # test that warm restart leads to the same point
    X, y = make_regression(random_state=0, n_samples=1000)
    warm = QuantileRegressor(fit_intercept=True, alpha=1.0, max_iter=10000,
                             warm_start=True, gamma=1e-10,
                             xtol=1e-10, gtol=1e-10)
    warm.fit(X, y)
    warm_coef = warm.coef_.copy()
    warm_iter = sum(warm.total_iter_)
    warm.fit(X, y)

    # SciPy performs the tol check after doing the coef updates, so
    # these would be almost same but not necessarily equal.
    assert_allclose(warm.coef_, warm_coef, atol=1e-1)
    # assert a smaller number of iterations than the first fit
    assert sum(warm.total_iter_) < warm_iter


def test_quantile_convergence():
    # Quantile loss may not converge to unique solution
    # if there is no regularization
    # need to check that warning is not thrown if model has converged.
    X, y = make_regression(n_samples=300, n_features=20, random_state=0,
                           noise=1.0)

    # check that for small n_iter, warning is thrown
    with pytest.warns(None) as record:
        QuantileRegressor(max_iter=1).fit(X, y)
        assert len(record) == 1
        assert 'QuantileRegressor did not converge' in str(record[-1].message)

    # check that for large n_iter, it is not thrown
    with pytest.warns(None) as record:
        QuantileRegressor(max_iter=10000).fit(X, y)
        assert len(record) == 0
