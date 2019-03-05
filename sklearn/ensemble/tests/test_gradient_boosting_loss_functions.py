"""
Testing for the gradient boosting loss functions and initial estimators.
"""

import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_allclose
from numpy.testing import assert_equal

from sklearn.utils import check_random_state
from sklearn.utils.stats import _weighted_percentile
from sklearn.ensemble._gb_losses import RegressionLossFunction
from sklearn.ensemble._gb_losses import LeastSquaresError
from sklearn.ensemble._gb_losses import LeastAbsoluteError
from sklearn.ensemble._gb_losses import HuberLossFunction
from sklearn.ensemble._gb_losses import QuantileLossFunction
from sklearn.ensemble._gb_losses import BinomialDeviance
from sklearn.ensemble._gb_losses import MultinomialDeviance
from sklearn.ensemble._gb_losses import ExponentialLoss
from sklearn.ensemble._gb_losses import LOSS_FUNCTIONS


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
        out = loss.get_init_raw_predictions(X, init_est)
        assert_equal(out.shape, (y.shape[0], 1))

        sw_init_est = loss.init_estimator()
        sw_init_est.fit(X, y, sample_weight=sample_weight)
        sw_out = loss.get_init_raw_predictions(X, sw_init_est)
        assert_equal(sw_out.shape, (y.shape[0], 1))

        # check if predictions match
        assert_allclose(out, sw_out, rtol=1e-2)


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


def test_init_raw_predictions_shapes():
    # Make sure get_init_raw_predictions returns float64 arrays with shape
    # (n_samples, K) where K is 1 for binary classification and regression, and
    # K = n_classes for multiclass classification
    rng = np.random.RandomState(0)

    n_samples = 100
    X = rng.normal(size=(n_samples, 5))
    y = rng.normal(size=n_samples)
    for loss in (LeastSquaresError(n_classes=1),
                 LeastAbsoluteError(n_classes=1),
                 QuantileLossFunction(n_classes=1),
                 HuberLossFunction(n_classes=1)):
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        assert raw_predictions.shape == (n_samples, 1)
        assert raw_predictions.dtype == np.float64

    y = rng.randint(0, 2, size=n_samples)
    for loss in (BinomialDeviance(n_classes=2),
                 ExponentialLoss(n_classes=2)):
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        assert raw_predictions.shape == (n_samples, 1)
        assert raw_predictions.dtype == np.float64

    for n_classes in range(3, 5):
        y = rng.randint(0, n_classes, size=n_samples)
        loss = MultinomialDeviance(n_classes=n_classes)
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        assert raw_predictions.shape == (n_samples, n_classes)
        assert raw_predictions.dtype == np.float64


def test_init_raw_predictions_values():
    # Make sure the get_init_raw_predictions() returns the expected values for
    # each loss.
    rng = np.random.RandomState(0)

    n_samples = 100
    X = rng.normal(size=(n_samples, 5))
    y = rng.normal(size=n_samples)

    # Least squares loss
    loss = LeastSquaresError(n_classes=1)
    init_estimator = loss.init_estimator().fit(X, y)
    raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
    # Make sure baseline prediction is the mean of all targets
    assert_almost_equal(raw_predictions, y.mean())

    # Least absolute and huber loss
    for Loss in (LeastAbsoluteError, HuberLossFunction):
        loss = Loss(n_classes=1)
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        # Make sure baseline prediction is the median of all targets
        assert_almost_equal(raw_predictions, np.median(y))

    # Quantile loss
    for alpha in (.1, .5, .9):
        loss = QuantileLossFunction(n_classes=1, alpha=alpha)
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        # Make sure baseline prediction is the alpha-quantile of all targets
        assert_almost_equal(raw_predictions, np.percentile(y, alpha * 100))

    y = rng.randint(0, 2, size=n_samples)

    # Binomial deviance
    loss = BinomialDeviance(n_classes=2)
    init_estimator = loss.init_estimator().fit(X, y)
    # Make sure baseline prediction is equal to link_function(p), where p
    # is the proba of the positive class. We want predict_proba() to return p,
    # and by definition
    # p = inverse_link_function(raw_prediction) = sigmoid(raw_prediction)
    # So we want raw_prediction = link_function(p) = log(p / (1 - p))
    raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
    p = y.mean()
    assert_almost_equal(raw_predictions, np.log(p / (1 - p)))

    # Exponential loss
    loss = ExponentialLoss(n_classes=2)
    init_estimator = loss.init_estimator().fit(X, y)
    raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
    p = y.mean()
    assert_almost_equal(raw_predictions, .5 * np.log(p / (1 - p)))

    # Multinomial deviance loss
    for n_classes in range(3, 5):
        y = rng.randint(0, n_classes, size=n_samples)
        loss = MultinomialDeviance(n_classes=n_classes)
        init_estimator = loss.init_estimator().fit(X, y)
        raw_predictions = loss.get_init_raw_predictions(y, init_estimator)
        for k in range(n_classes):
            p = (y == k).mean()
        assert_almost_equal(raw_predictions[:, k], np.log(p))
