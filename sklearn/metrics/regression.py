"""Metrics to assess performance on regression task

Functions named as ``*_score`` return a scalar value to maximize: the higher
the better

Function named as ``*_error`` or ``*_loss`` return a scalar value to minimize:
the lower the better
"""

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Arnaud Joly <a.joly@ulg.ac.be>
#          Jochen Wersdorfer <jochen@wersdoerfer.de>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
#          Joel Nothman <joel.nothman@gmail.com>
#          Noel Dawe <noel@dawe.me>
# License: BSD 3 clause

from __future__ import division

import numpy as np

from ..utils.validation import check_array, check_consistent_length
from ..utils.validation import column_or_1d

__ALL__ = [
    "mean_absolute_error",
    "mean_squared_error",
    "r2_score",
    "explained_variance_score"
]


def _check_reg_targets(y_true, y_pred):
    """Check that y_true and y_pred belong to the same regression task

    Parameters
    ----------
    y_true : array-like,

    y_pred : array-like,

    Returns
    -------
    type_true : one of {'continuous', continuous-multioutput'}
        The type of the true target data, as output by
        ``utils.multiclass.type_of_target``

    y_true : array-like of shape = [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples, n_outputs]
        Estimated target values.
    """
    check_consistent_length(y_true, y_pred)
    y_true = check_array(y_true, ensure_2d=False)
    y_pred = check_array(y_pred, ensure_2d=False)

    if y_true.ndim == 1:
        y_true = y_true.reshape((-1, 1))

    if y_pred.ndim == 1:
        y_pred = y_pred.reshape((-1, 1))

    if y_true.shape[1] != y_pred.shape[1]:
        raise ValueError("y_true and y_pred have different number of output "
                         "({0}!={1})".format(y_true.shape[1], y_pred.shape[1]))

    y_type = 'continuous' if y_true.shape[1] == 1 else 'continuous-multioutput'

    return y_type, y_true, y_pred


def _average_and_variance(values, sample_weight=None):
    """
    Compute the (weighted) average and variance.

    Parameters
    ----------
    values : array-like of shape = [n_samples] or [n_samples, n_outputs]

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    average : float
        The weighted average

    variance : float
        The weighted variance
    """
    values = np.asarray(values)
    if values.ndim == 1:
        values = values.reshape((-1, 1))
    if sample_weight is not None:
        sample_weight = np.asarray(sample_weight)
        if sample_weight.ndim == 1:
            sample_weight = sample_weight.reshape((-1, 1))
    average = np.average(values, weights=sample_weight)
    variance = np.average((values - average)**2, weights=sample_weight)
    return average, variance


def mean_absolute_error(y_true, y_pred, sample_weight=None):
    """Mean absolute error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    loss : float
        A positive floating point value (the best value is 0.0).

    Examples
    --------
    >>> from sklearn.metrics import mean_absolute_error
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> mean_absolute_error(y_true, y_pred)
    0.5
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> mean_absolute_error(y_true, y_pred)
    0.75

    """
    y_type, y_true, y_pred = _check_reg_targets(y_true, y_pred)
    return np.average(np.abs(y_pred - y_true).mean(axis=1),
                      weights=sample_weight)


def mean_squared_error(y_true, y_pred, sample_weight=None):
    """Mean squared error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    loss : float
        A positive floating point value (the best value is 0.0).

    Examples
    --------
    >>> from sklearn.metrics import mean_squared_error
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> mean_squared_error(y_true, y_pred)
    0.375
    >>> y_true = [[0.5, 1],[-1, 1],[7, -6]]
    >>> y_pred = [[0, 2],[-1, 2],[8, -5]]
    >>> mean_squared_error(y_true, y_pred)  # doctest: +ELLIPSIS
    0.708...

    """
    y_type, y_true, y_pred = _check_reg_targets(y_true, y_pred)
    return np.average(((y_pred - y_true) ** 2).mean(axis=1),
                      weights=sample_weight)


def explained_variance_score(y_true, y_pred, sample_weight=None):
    """Explained variance regression score function

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like
        Ground truth (correct) target values.

    y_pred : array-like
        Estimated target values.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    score : float
        The explained variance.

    Notes
    -----
    This is not a symmetric function.

    Examples
    --------
    >>> from sklearn.metrics import explained_variance_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> explained_variance_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.957...

    """
    y_type, y_true, y_pred = _check_reg_targets(y_true, y_pred)

    if y_type != "continuous":
        raise ValueError("{0} is not supported".format(y_type))

    _, numerator = _average_and_variance(y_true - y_pred, sample_weight)
    _, denominator = _average_and_variance(y_true, sample_weight)
    if denominator == 0.0:
        if numerator == 0.0:
            return 1.0
        else:
            # arbitrary set to zero to avoid -inf scores, having a constant
            # y_true is not interesting for scoring a regression anyway
            return 0.0
    return 1 - numerator / denominator


def r2_score(y_true, y_pred, sample_weight=None):
    """R^2 (coefficient of determination) regression score function.

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    z : float
        The R^2 score.

    Notes
    -----
    This is not a symmetric function.

    Unlike most other scores, R^2 score may be negative (it need not actually
    be the square of a quantity R).

    References
    ----------
    .. [1] `Wikipedia entry on the Coefficient of determination
            <http://en.wikipedia.org/wiki/Coefficient_of_determination>`_

    Examples
    --------
    >>> from sklearn.metrics import r2_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.948...
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.938...

    """
    y_type, y_true, y_pred = _check_reg_targets(y_true, y_pred)

    if sample_weight is not None:
        sample_weight = column_or_1d(sample_weight)
        weight = sample_weight[:, np.newaxis]
    else:
        weight = 1.

    numerator = (weight * (y_true - y_pred) ** 2).sum(dtype=np.float64)
    denominator = (weight * (y_true - np.average(
        y_true, axis=0, weights=sample_weight)) ** 2).sum(dtype=np.float64)

    if denominator == 0.0:
        if numerator == 0.0:
            return 1.0
        else:
            # arbitrary set to zero to avoid -inf scores, having a constant
            # y_true is not interesting for scoring a regression anyway
            return 0.0

    return 1 - numerator / denominator
