# Authors: David Dale dale.david@mail.ru
# License: BSD 3 clause

import numpy as np
import pytest
import scipy.stats

from sklearn.utils import check_random_state
from sklearn.utils._testing import assert_allclose
from sklearn.datasets import make_regression
from sklearn.linear_model import HuberRegressor, QuantileRegressor
from sklearn.utils.fixes import sp_version, parse_version

_SCIPY_TOO_OLD = "requires at least scipy 1.0.0"


@pytest.mark.parametrize(
    'quantile, alpha, intercept, coef',
    [
        # for 50% quantile w/o regularization, any slope in [1, 10] is okay
        [0.5, 0, 1, None],
        # if positive error costs more, the slope is maximal
        [0.51, 0, 1, 10],
        # if negative error costs more, the slope is minimal
        [0.49, 0, 1, 1],
        # for a small lasso penalty, the slope is also minimal
        [0.5, 0.01, 1, 1],
        # for a large lasso penalty, the model predicts the constant median
        [0.5, 100, 2, 0],
    ]
)
@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
def test_quantile_toy_example(quantile, alpha, intercept, coef):
    # test how different parameters affect a small intuitive example
    X = [[0], [1], [1]]
    y = [1, 2, 11]
    model = QuantileRegressor(quantile=quantile, alpha=alpha).fit(X, y)
    assert_allclose(model.intercept_, intercept, atol=1e-2)
    if coef is not None:
        assert_allclose(model.coef_[0], coef, atol=1e-2)
    if alpha < 100:
        assert model.coef_[0] >= 1
    assert model.coef_[0] <= 10


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
def test_quantile_equals_huber_for_low_epsilon():
    X, y = make_regression(n_samples=100, n_features=20, random_state=0,
                           noise=1.0)
    huber = HuberRegressor(epsilon=1+1e-4, alpha=1e-4).fit(X, y)
    quant = QuantileRegressor(alpha=1e-4).fit(X, y)
    assert_allclose(huber.intercept_, quant.intercept_, atol=1e-1)
    assert_allclose(huber.coef_, quant.coef_, atol=1e-1)


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
@pytest.mark.parametrize("q", [0.5, 0.9, 0.05])
def test_quantile_estimates_calibration(q):
    # Test that model estimates percentage of points below the prediction
    X, y = make_regression(n_samples=1000, n_features=20, random_state=0,
                           noise=1.0)
    quant = QuantileRegressor(quantile=q, alpha=0).fit(X, y)
    fraction_below = np.mean(y < quant.predict(X))
    assert_allclose(fraction_below, q, atol=1e-2)


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
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


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
def test_quantile_sample_weight():
    # test that with unequal sample weights we still estimate weighted fraction
    n = 1000
    X, y = make_regression(n_samples=n, n_features=10, random_state=0,
                           noise=10.0)
    weight = np.ones(n)
    # when we increase weight of upper observations,
    # estimate of quantile should go up
    weight[y > y.mean()] = 100
    quant = QuantileRegressor(quantile=0.5, alpha=1e-8)
    quant.fit(X, y, sample_weight=weight)
    fraction_below = np.mean(y < quant.predict(X))
    assert fraction_below > 0.5
    weighted_fraction_below = np.sum((y < quant.predict(X)) * weight) \
        / np.sum(weight)
    assert_allclose(weighted_fraction_below, 0.5, atol=1e-2)


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
@pytest.mark.parametrize('quantile', [2.0, 1.0, 0.0, -1])
def test_quantile_incorrect_quantile(quantile):
    X, y = make_regression(n_samples=10, n_features=1, random_state=0, noise=1)
    with pytest.raises(
        ValueError,
        match="Quantile should be strictly between 0.0 and 1.0"
    ):
        QuantileRegressor(quantile=quantile).fit(X, y)


@pytest.mark.skipif(sp_version < parse_version("1.0.0"), reason=_SCIPY_TOO_OLD)
@pytest.mark.parametrize('quantile', [0.1, 0.5, 0.9])
def test_asymmetric_error(quantile):
    n_samples = 1000
    n_features = 3
    bias = 12.3
    param = 1.0
    generator = check_random_state(42)
    X = generator.randn(n_samples, n_features)
    ground_truth = generator.rand(n_features)
    dist = scipy.stats.expon(param)
    noise = dist.rvs(size=n_samples, random_state=42)
    y = np.dot(X, ground_truth) + bias + noise
    model = QuantileRegressor(quantile=quantile).fit(X, y)

    assert_allclose(model.intercept_, bias + dist.ppf(quantile), atol=0.1)
    assert_allclose(model.coef_, ground_truth, atol=0.3)
    assert_allclose(np.mean(model.predict(X) > y), quantile, atol=0.003)
