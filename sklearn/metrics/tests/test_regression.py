from __future__ import division, print_function

import numpy as np
from itertools import product

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal

from sklearn.metrics import explained_variance_score
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
from sklearn.metrics import median_absolute_error
from sklearn.metrics import r2_score

from sklearn.metrics.regression import _check_reg_targets


def test_regression_metrics(n_samples=50):
    y_true = np.arange(n_samples)
    y_pred = y_true + 1

    assert_almost_equal(mean_squared_error(y_true, y_pred), 1.)
    assert_almost_equal(mean_absolute_error(y_true, y_pred), 1.)
    assert_almost_equal(median_absolute_error(y_true, y_pred), 1.)
    assert_almost_equal(r2_score(y_true, y_pred),  0.995, 2)
    assert_almost_equal(explained_variance_score(y_true, y_pred), 1.)


def test_multioutput_regression():
    y_true = np.array([[1, 0, 0, 1], [0, 1, 1, 1], [1, 1, 0, 1]])
    y_pred = np.array([[0, 0, 0, 1], [1, 0, 1, 1], [0, 0, 0, 1]])

    error = mean_squared_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    # mean_absolute_error and mean_squared_error are equal because
    # it is a binary problem.
    error = mean_absolute_error(y_true, y_pred)
    assert_almost_equal(error, (1. / 3 + 2. / 3 + 2. / 3) / 4.)

    error = r2_score(y_true, y_pred)
    assert_almost_equal(error, 1 - 5. / 2)


def test_regression_metrics_at_limits():
    assert_almost_equal(mean_squared_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(mean_absolute_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(median_absolute_error([0.], [0.]), 0.00, 2)
    assert_almost_equal(explained_variance_score([0.], [0.]), 1.00, 2)
    assert_almost_equal(r2_score([0., 1], [0., 1]), 1.00, 2)


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
            y_type, y_check1, y_check2 = _check_reg_targets(y1, y2)
            assert_equal(type1, y_type)
            if type1 == 'continuous':
                assert_array_equal(y_check1, np.reshape(y1, (-1, 1)))
                assert_array_equal(y_check2, np.reshape(y2, (-1, 1)))
            else:
                assert_array_equal(y_check1, y1)
                assert_array_equal(y_check2, y2)
        else:
            assert_raises(ValueError, _check_reg_targets, y1, y2)
