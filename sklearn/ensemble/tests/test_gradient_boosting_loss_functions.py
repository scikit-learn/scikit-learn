"""
Testing for the gradient boosting loss functions and initial estimators.
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises
from sklearn.ensemble.gradient_boosting import BinomialDeviance
from sklearn.ensemble.gradient_boosting import LogOddsEstimator
from sklearn.ensemble.gradient_boosting import LeastSquaresError
from sklearn.ensemble.gradient_boosting import RegressionLossFunction
from sklearn.ensemble.gradient_boosting import LOSS_FUNCTIONS
from sklearn.ensemble.gradient_boosting import _weighted_percentile


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
