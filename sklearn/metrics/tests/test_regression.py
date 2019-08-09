
import numpy as np
from numpy.testing import assert_allclose
from itertools import product
import pytest

from sklearn.utils.testing import assert_raises, assert_raises_regex
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal

from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_squared_log_error
from sklearn.metrics import median_absolute_error
from sklearn.metrics import max_error
from sklearn.metrics import r2_score
from sklearn.metrics import mean_tweedie_deviance

from sklearn.metrics.regression import _check_reg_targets

from ...exceptions import UndefinedMetricWarning


def test_regression_metrics(n_samples=50):
    y_true = np.arange(n_samples)
    y_pred = y_true + 1

    assert_almost_equal(mean_squared_error(y_true, y_pred), 1.)
    assert_almost_equal(mean_squared_log_error(y_true, y_pred),
                        mean_squared_error(np.log(1 + y_true),
                                           np.log(1 + y_pred)))
    assert_almost_equal(mean_absolute_error(y_true, y_pred), 1.)
    assert_almost_equal(median_absolute_error(y_true, y_pred), 1.)
    assert_almost_equal(max_error(y_true, y_pred), 1.)
    assert_almost_equal(r2_score(y_true, y_pred),  0.995, 2)
    assert_almost_equal(explained_variance_score(y_true, y_pred), 1.)
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=0),
                        mean_squared_error(y_true, y_pred))

    # Tweedie deviance needs positive y_pred, except for p=0,
    # p>=2 needs positive y_true
    # results evaluated by sympy
    y_true = np.arange(1, 1 + n_samples)
    y_pred = 2 * y_true
    n = n_samples
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=-1),
                        5/12 * n * (n**2 + 2 * n + 1))
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=1),
                        (n + 1) * (1 - np.log(2)))
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=2),
                        2 * np.log(2) - 1)
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=3/2),
                        ((6 * np.sqrt(2) - 8) / n) * np.sqrt(y_true).sum())
    assert_almost_equal(mean_tweedie_deviance(y_true, y_pred, p=3),
                        np.sum(1 / y_true) / (4 * n))


def test_multioutput_regression():
    y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
    y_pred = np.array([[0, 0, 0, 1], [1, 0, 1, 1], [0, 0, 0, 1]])

    error = mean_squared_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    error = mean_squared_error(y_true, y_pred, squared=False)
    assert_almost_equal(error, 0.645, decimal=2)

    error = mean_squared_log_error(y_true, y_pred)
    assert_almost_equal(error, 0.200, decimal=2)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    error = mean_absolute_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    error = r2_score(y_true, y_pred, multioutput='variance_weighted')
    assert_almost_equal(error, 1. - 5. / 2)
    error = r2_score(y_true, y_pred, multioutput='uniform_average')
    assert_almost_equal(error, -.875)


def test_regression_metrics_at_limits():
    assert_almost_equal(mean_squared_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(mean_squared_error([0.], [0.], squared=False), 0.00, 2)
    assert_almost_equal(mean_squared_log_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(mean_absolute_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(median_absolute_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(max_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(explained_variance_score([0.], [0.]), 1.00, 2)
    assert_almost_equal(r2_score([0., 1], [0., 1]), 1.00, 2)
    assert_raises_regex(ValueError, "Mean Squared Logarithmic Error cannot be "
                        "used when targets contain negative values.",
                        mean_squared_log_error, [-1.], [-1.])
    assert_raises_regex(ValueError, "Mean Squared Logarithmic Error cannot be "
                        "used when targets contain negative values.",
                        mean_squared_log_error, [1., 2., 3.], [1., -2., 3.])
    assert_raises_regex(ValueError, "Mean Squared Logarithmic Error cannot be "
                        "used when targets contain negative values.",
                        mean_squared_log_error, [1., -2., 3.], [1., 2., 3.])

    # Tweedie deviance error
    p = -1.2
    assert_allclose(mean_tweedie_deviance([0], [1.], p=p),
                    2./(2.-p), rtol=1e-3)
    with pytest.raises(ValueError,
                       match="can only be used on strictly positive y_pred."):
        mean_tweedie_deviance([0.], [0.], p=p)
    assert_almost_equal(mean_tweedie_deviance([0.], [0.], p=0), 0.00, 2)

    msg = "only be used on non-negative y_true and strictly positive y_pred."
    with pytest.raises(ValueError, match=msg):
        mean_tweedie_deviance([0.], [0.], p=1.0)

    p = 1.5
    assert_allclose(mean_tweedie_deviance([0.], [1.], p=p), 2./(2.-p))
    msg = "only be used on non-negative y_true and strictly positive y_pred."
    with pytest.raises(ValueError, match=msg):
        mean_tweedie_deviance([0.], [0.], p=p)
    p = 2.
    assert_allclose(mean_tweedie_deviance([1.], [1.], p=p), 0.00,
                    atol=1e-8)
    msg = "can only be used on strictly positive y_true and y_pred."
    with pytest.raises(ValueError, match=msg):
        mean_tweedie_deviance([0.], [0.], p=p)
    p = 3.
    assert_allclose(mean_tweedie_deviance([1.], [1.], p=p),
                    0.00, atol=1e-8)

    msg = "can only be used on strictly positive y_true and y_pred."
    with pytest.raises(ValueError, match=msg):
        mean_tweedie_deviance([0.], [0.], p=p)

    with pytest.raises(ValueError,
                       match="deviance is only defined for p<=0 and p>=1."):
        mean_tweedie_deviance([0.], [0.], p=0.5)


def test__check_reg_targets():
    # All of length 3
    EXAMPLES = [
        ("continuous", [1, 2, 3], 1),
        ("continuous", [[1], [2], [3]], 1),
        ("continuous-multioutput", [[1, 1], [2, 2], [3, 1]], 2),
        ("continuous-multioutput", [[5, 1], [4, 2], [3, 1]], 2),
        ("continuous-multioutput", [[1, 3, 4], [2, 2, 2], [3, 1, 1]], 3),
    ]

    for (type1, y1, n_out1), (type2, y2, n_out2) in product(EXAMPLES,
                                                            repeat=2):

        if type1 == type2 and n_out1 == n_out2:
            y_type, y_check1, y_check2, multioutput = _check_reg_targets(
                y1, y2, None)
            assert type1 == y_type
            if type1 == 'continuous':
                assert_array_equal(y_check1, np.reshape(y1, (-1, 1)))
                assert_array_equal(y_check2, np.reshape(y2, (-1, 1)))
            else:
                assert_array_equal(y_check1, y1)
                assert_array_equal(y_check2, y2)
        else:
            assert_raises(ValueError, _check_reg_targets, y1, y2, None)


def test__check_reg_targets_exception():
    invalid_multioutput = 'this_value_is_not_valid'
    expected_message = ("Allowed 'multioutput' string values are.+"
                        "You provided multioutput={!r}".format(
                            invalid_multioutput))
    assert_raises_regex(ValueError, expected_message,
                        _check_reg_targets,
                        [1, 2, 3],
                        [[1], [2], [3]],
                        invalid_multioutput)


def test_regression_multioutput_array():
    y_true = [[1, 2], [2.5, -1], [4.5, 3], [5, 7]]
    y_pred = [[1, 1], [2, -1], [5, 4], [5, 6.5]]

    mse = mean_squared_error(y_true, y_pred, multioutput='raw_values')
    mae = mean_absolute_error(y_true, y_pred, multioutput='raw_values')
    r = r2_score(y_true, y_pred, multioutput='raw_values')
    evs = explained_variance_score(y_true, y_pred, multioutput='raw_values')

    assert_array_almost_equal(mse, [0.125, 0.5625], decimal=2)
    assert_array_almost_equal(mae, [0.25, 0.625], decimal=2)
    assert_array_almost_equal(r, [0.95, 0.93], decimal=2)
    assert_array_almost_equal(evs, [0.95, 0.93], decimal=2)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    y_true = [[0, 0]]*4
    y_pred = [[1, 1]]*4
    mse = mean_squared_error(y_true, y_pred, multioutput='raw_values')
    mae = mean_absolute_error(y_true, y_pred, multioutput='raw_values')
    r = r2_score(y_true, y_pred, multioutput='raw_values')
    assert_array_almost_equal(mse, [1., 1.], decimal=2)
    assert_array_almost_equal(mae, [1., 1.], decimal=2)
    assert_array_almost_equal(r, [0., 0.], decimal=2)

    r = r2_score([[0, -1], [0, 1]], [[2, 2], [1, 1]], multioutput='raw_values')
    assert_array_almost_equal(r, [0, -3.5], decimal=2)
    assert np.mean(r) == r2_score([[0, -1], [0, 1]], [[2, 2], [1, 1]],
                                  multioutput='uniform_average')
    evs = explained_variance_score([[0, -1], [0, 1]], [[2, 2], [1, 1]],
                                   multioutput='raw_values')
    assert_array_almost_equal(evs, [0, -1.25], decimal=2)

    # Checking for the condition in which both numerator and denominator is
    # zero.
    y_true = [[1, 3], [-1, 2]]
    y_pred = [[1, 4], [-1, 1]]
    r2 = r2_score(y_true, y_pred, multioutput='raw_values')
    assert_array_almost_equal(r2, [1., -3.], decimal=2)
    assert np.mean(r2) == r2_score(y_true, y_pred,
                                   multioutput='uniform_average')
    evs = explained_variance_score(y_true, y_pred, multioutput='raw_values')
    assert_array_almost_equal(evs, [1., -3.], decimal=2)
    assert np.mean(evs) == explained_variance_score(y_true, y_pred)

    # Handling msle separately as it does not accept negative inputs.
    y_true = np.array([[0.5, 1], [1, 2], [7, 6]])
    y_pred = np.array([[0.5, 2], [1, 2.5], [8, 8]])
    msle = mean_squared_log_error(y_true, y_pred, multioutput='raw_values')
    msle2 = mean_squared_error(np.log(1 + y_true), np.log(1 + y_pred),
                               multioutput='raw_values')
    assert_array_almost_equal(msle, msle2, decimal=2)


def test_regression_custom_weights():
    y_true = [[1, 2], [2.5, -1], [4.5, 3], [5, 7]]
    y_pred = [[1, 1], [2, -1], [5, 4], [5, 6.5]]

    msew = mean_squared_error(y_true, y_pred, multioutput=[0.4, 0.6])
    rmsew = mean_squared_error(y_true, y_pred, multioutput=[0.4, 0.6],
                               squared=False)
    maew = mean_absolute_error(y_true, y_pred, multioutput=[0.4, 0.6])
    rw = r2_score(y_true, y_pred, multioutput=[0.4, 0.6])
    evsw = explained_variance_score(y_true, y_pred, multioutput=[0.4, 0.6])

    assert_almost_equal(msew, 0.39, decimal=2)
    assert_almost_equal(rmsew, 0.62, decimal=2)
    assert_almost_equal(maew, 0.475, decimal=3)
    assert_almost_equal(rw, 0.94, decimal=2)
    assert_almost_equal(evsw, 0.94, decimal=2)

    # Handling msle separately as it does not accept negative inputs.
    y_true = np.array([[0.5, 1], [1, 2], [7, 6]])
    y_pred = np.array([[0.5, 2], [1, 2.5], [8, 8]])
    msle = mean_squared_log_error(y_true, y_pred, multioutput=[0.3, 0.7])
    msle2 = mean_squared_error(np.log(1 + y_true), np.log(1 + y_pred),
                               multioutput=[0.3, 0.7])
    assert_almost_equal(msle, msle2, decimal=2)


@pytest.mark.parametrize('metric', [r2_score])
def test_regression_single_sample(metric):
    y_true = [0]
    y_pred = [1]
    warning_msg = 'not well-defined with less than two samples.'

    # Trigger the warning
    with pytest.warns(UndefinedMetricWarning, match=warning_msg):
        score = metric(y_true, y_pred)
        assert np.isnan(score)


def test_tweedie_deviance_continuity():
    n_samples = 100

    y_true = np.random.RandomState(0).rand(n_samples) + 0.1
    y_pred = np.random.RandomState(1).rand(n_samples) + 0.1

    assert_allclose(mean_tweedie_deviance(y_true, y_pred, p=0 - 1e-10),
                    mean_tweedie_deviance(y_true, y_pred, p=0))

    # Ws we get closer to the limit, with 1e-12 difference the absolute
    # tolerance to pass the below check increases. There are likely
    # numerical precision issues on the edges of different definition
    # regions.
    assert_allclose(mean_tweedie_deviance(y_true, y_pred, p=1 + 1e-10),
                    mean_tweedie_deviance(y_true, y_pred, p=1),
                    atol=1e-6)

    assert_allclose(mean_tweedie_deviance(y_true, y_pred, p=2 - 1e-10),
                    mean_tweedie_deviance(y_true, y_pred, p=2),
                    atol=1e-6)

    assert_allclose(mean_tweedie_deviance(y_true, y_pred, p=2 + 1e-10),
                    mean_tweedie_deviance(y_true, y_pred, p=2),
                    atol=1e-6)
