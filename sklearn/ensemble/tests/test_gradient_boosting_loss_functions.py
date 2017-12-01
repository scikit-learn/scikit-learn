"""
Testing for the gradient boosting loss functions and initial estimators.
"""

import pytest

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_raises_regex

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises
from sklearn.ensemble.gradient_boosting import BinomialDeviance
from sklearn.ensemble.gradient_boosting import LogOddsEstimator
from sklearn.ensemble.gradient_boosting import LeastSquaresError
from sklearn.ensemble.gradient_boosting import MultinomialDeviance
from sklearn.ensemble.gradient_boosting import RegressionLossFunction
from sklearn.ensemble.gradient_boosting import LOSS_FUNCTIONS
from sklearn.ensemble.gradient_boosting import _weighted_percentile
from sklearn.ensemble.gradient_boosting import QuantileLossFunction


def test_binomial_deviance():
    # Check binomial deviance loss.
    # Check against alternative definitions in ESLII.
    bd = BinomialDeviance(2)

    # pred has the same BD for y in {0, 1}
    assert_equal(bd(np.array([0.0]), np.array([0.0])),
                 bd(np.array([1.0]), np.array([0.0])))

    assert_almost_equal(bd(np.array([1.0, 1.0, 1.0]),
                           np.array([100.0, 100.0, 100.0])),
                        0.0)
    assert_almost_equal(bd(np.array([1.0, 0.0, 0.0]),
                           np.array([100.0, -100.0, -100.0])), 0)

    # check if same results as alternative definition of deviance (from ESLII)
    alt_dev = lambda y, pred: np.mean(np.logaddexp(0.0, -2.0 *
                                                   (2.0 * y - 1) * pred))
    test_data = [(np.array([1.0, 1.0, 1.0]), np.array([100.0, 100.0, 100.0])),
                 (np.array([0.0, 0.0, 0.0]), np.array([100.0, 100.0, 100.0])),
                 (np.array([0.0, 0.0, 0.0]),
                  np.array([-100.0, -100.0, -100.0])),
                 (np.array([1.0, 1.0, 1.0]),
                  np.array([-100.0, -100.0, -100.0]))]

    for datum in test_data:
        assert_almost_equal(bd(*datum), alt_dev(*datum))

    # check the gradient against the
    alt_ng = lambda y, pred: (2 * y - 1) / (1 + np.exp(2 * (2 * y - 1) * pred))
    for datum in test_data:
        assert_almost_equal(bd.negative_gradient(*datum), alt_ng(*datum))


def test_log_odds_estimator():
    # Check log odds estimator.
    est = LogOddsEstimator()
    assert_raises(ValueError, est.fit, None, np.array([1]))

    est.fit(None, np.array([1.0, 0.0]))
    assert_equal(est.prior, 0.0)
    assert_array_equal(est.predict(np.array([[1.0], [1.0]])),
                       np.array([[0.0], [0.0]]))


def test_sample_weight_smoke():
    rng = check_random_state(13)
    y = rng.rand(100)
    pred = rng.rand(100)

    # least squares
    loss = LeastSquaresError(1)
    loss_wo_sw = loss(y, pred)
    loss_w_sw = loss(y, pred, np.ones(pred.shape[0], dtype=np.float32))
    assert_almost_equal(loss_wo_sw, loss_w_sw)


def test_sample_weight_init_estimators():
    # Smoke test for init estimators with sample weights.
    rng = check_random_state(13)
    X = rng.rand(100, 2)
    sample_weight = np.ones(100)
    reg_y = rng.rand(100)

    clf_y = rng.randint(0, 2, size=100)

    for Loss in LOSS_FUNCTIONS.values():
        if Loss is None:
            continue
        if issubclass(Loss, RegressionLossFunction):
            k = 1
            y = reg_y
        else:
            k = 2
            y = clf_y
            if Loss.is_multi_class:
                # skip multiclass
                continue

        loss = Loss(k)
        init_est = loss.init_estimator()
        init_est.fit(X, y)
        out = init_est.predict(X)
        assert_equal(out.shape, (y.shape[0], 1))

        sw_init_est = loss.init_estimator()
        sw_init_est.fit(X, y, sample_weight=sample_weight)
        sw_out = init_est.predict(X)
        assert_equal(sw_out.shape, (y.shape[0], 1))

        # check if predictions match
        assert_array_equal(out, sw_out)


def test_weighted_percentile():
    y = np.empty(102, dtype=np.float64)
    y[:50] = 0
    y[-51:] = 2
    y[-1] = 100000
    y[50] = 1
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert score == 1


def test_weighted_percentile_equal():
    y = np.empty(102, dtype=np.float64)
    y.fill(0.0)
    sw = np.ones(102, dtype=np.float64)
    sw[-1] = 0.0
    score = _weighted_percentile(y, sw, 50)
    assert score == 0


def test_weighted_percentile_zero_weight():
    y = np.empty(102, dtype=np.float64)
    y.fill(1.0)
    sw = np.ones(102, dtype=np.float64)
    sw.fill(0.0)
    score = _weighted_percentile(y, sw, 50)
    assert score == 1.0


def test_quantile_loss_function():
    # Non regression test for the QuantileLossFunction object
    # There was a sign problem when evaluating the function
    # for negative values of 'ytrue - ypred'
    x = np.asarray([-1.0, 0.0, 1.0])
    y_found = QuantileLossFunction(1, 0.9)(x, np.zeros_like(x))
    y_expected = np.asarray([0.1, 0.0, 0.9]).mean()
    np.testing.assert_allclose(y_found, y_expected)


def test_sample_weight_deviance():
    # Test if deviance supports sample weights.
    rng = check_random_state(13)
    X = rng.rand(100, 2)
    sample_weight = np.ones(100)
    reg_y = rng.rand(100)
    clf_y = rng.randint(0, 2, size=100)
    mclf_y = rng.randint(0, 3, size=100)

    for Loss in LOSS_FUNCTIONS.values():
        if Loss is None:
            continue
        if issubclass(Loss, RegressionLossFunction):
            k = 1
            y = reg_y
            p = reg_y
        else:
            k = 2
            y = clf_y
            p = clf_y
            if Loss.is_multi_class:
                k = 3
                y = mclf_y
                # one-hot encoding
                p = np.zeros((y.shape[0], k), dtype=np.float64)
                for i in range(k):
                    p[:, i] = y == i

        loss = Loss(k)
        deviance_w_w = loss(y, p, sample_weight)
        deviance_wo_w = loss(y, p)
        assert deviance_wo_w == deviance_w_w


@pytest.mark.parametrize('n_classes, n_samples',
                         [(3, 100), (5, 57), (7, 13)])
def test_multinomial_deviance(n_classes, n_samples):
    # Check multinomial deviance with and without sample weights.
    rng = check_random_state(13)
    sample_weight = np.ones(n_samples)
    y = rng.randint(0, n_classes, size=sample_weight.shape[0])
    p = np.zeros((y.shape[0], n_classes), dtype=np.float64)
    for i in range(p.shape[1]):
        p[:, i] = y == i

    loss = MultinomialDeviance(n_classes)
    loss_wo_sw = loss(y, p)
    assert loss_wo_sw > 0
    loss_w_sw = loss(y, p, np.ones(p.shape[0], dtype=np.float32))
    assert_almost_equal(loss_wo_sw, loss_w_sw)
    loss_w_sw = loss(y, p, 0.5 * np.ones(p.shape[0], dtype=np.float32))
    assert_almost_equal(loss_wo_sw, loss_w_sw)


@pytest.mark.parametrize('pred, y, expected_loss',
                         [(np.array([[1.0, 0, 0],
                                     [0, 0.5, 0.5]]),
                           np.array([0, 1]),
                           0.75473)])
def test_mdl_computation_unweighted(pred, y, expected_loss):
    # MultinomialDeviance loss computation with uniform/without weights.
    weights = np.array([1, 1])
    loss = MultinomialDeviance(3)
    assert_almost_equal(loss(y, pred, weights), expected_loss, decimal=4)
    assert_almost_equal(loss(y, pred, None), expected_loss, decimal=4)
    assert_almost_equal(loss(y, pred), expected_loss, decimal=4)


@pytest.mark.parametrize('pred, y, weights, expected_loss',
                         [(np.array([[1.0, 0, 0],
                                     [0, 0.5, 0.5]]),
                           np.array([0, 1]),
                           np.array([1, 3]),
                           0.85637)])
def test_mdl_computation_weighted(pred, y, weights, expected_loss):
    # MultinomialDeviance loss computation with weights.
    loss = MultinomialDeviance(3)
    assert_almost_equal(loss(y, pred, weights), expected_loss, decimal=4)


@pytest.mark.parametrize('n', [0, 1, 2])
def test_mdl_exception(n):
    # Check that MultinomialDeviance throws an exception when n_classes <= 2
    err_msg = 'MultinomialDeviance requires more than 2 classes.'
    assert_raises_regex(ValueError, err_msg, MultinomialDeviance, n)
