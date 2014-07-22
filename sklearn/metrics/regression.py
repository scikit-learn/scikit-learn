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
#          Manoj Kumar <manojkumarsivaraj334@gmail.com>
#          Michael Eickenberg <michael.eickenberg@gmail.com>
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


def _check_reg_targets(y_true, y_pred, output_weights):
    """Check that y_true and y_pred belong to the same regression task

    Parameters
    ----------
    y_true : array-like,

    y_pred : array-like,

    output_weights : array-like or string in ['uniform', 'variance'] or None

    Returns
    -------
    type_true : one of {'continuous', continuous-multioutput'}
        The type of the true target data, as output by
        ``utils.multiclass.type_of_target``

    y_true : array-like of shape = [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples, n_outputs]
        Estimated target values.

    output_weights : array-like of shape = [n_outputs] or in ['uniform', 'variance', None (default)]
        Custom output weights is output_weights is array-like. 'uniform' and
        'variance' indicate specific weight vectors.
        None indicates no agglomeration of scores across targets: Scores are
        returned separately per target.

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

    n_outputs = y_true.shape[1]
    output_weights_options = (None, 'uniform', 'variance')
    if output_weights not in output_weights_options:
        output_weights = check_array(output_weights, ensure_2d=False)
        if n_outputs == 1:
            raise ValueError("Custom weights are useful only in "
                             "multi-output cases.")
        elif n_outputs != len(output_weights):
            raise ValueError(("There must be equally many custom weights "
                              "(%d) as outputs (%d).") % 
                             (len(output_weights), n_outputs))
    y_type = 'continuous' if n_outputs == 1 else 'continuous-multioutput'

    return y_type, y_true, y_pred, output_weights


def _average_and_variance(values, sample_weight=None, axis=None):
    """
    Compute the (weighted) average and variance.

    Parameters
    ----------
    values : array-like of shape = [n_samples] or [n_samples, n_outputs]

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    axis : integer or None, default None
        Axis along which to calculate average and variance. Full array by
        default.

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
    n_samples, n_outputs = values.shape
    if sample_weight is not None:
        sample_weight = np.asarray(sample_weight)
        if sample_weight.ndim == 1:
            sample_weight = sample_weight.reshape((n_samples, 1))
    # if multi output but sample weight only specified in one column,
    # then we need to broadcast it over outputs
    if (sample_weight is not None and n_outputs != sample_weight.shape[1]):
        if sample_weight.shape[1] != 1:
            raise ValueError("Sample weight shape and data shape "
                             "do not correspond.")
        sample_weight = sample_weight * np.ones([1, n_outputs])

    average = np.average(values, weights=sample_weight, axis=axis)
    variance = np.average((values - average) ** 2,
                          weights=sample_weight, axis=axis)
    return average, variance


def mean_absolute_error(y_true, y_pred,
                        output_weights='uniform',
                        sample_weight=None):
    """Mean absolute error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    output_weights : string in ['uniform'] or None
                     or array-like of shape [n_outputs]
        Weights by which to average scores of outputs. Useful only if using
        multiple outputs.

        ``ndarray`` :
            array containing weights for the weighted average.

        ``'uniform'`` :
            Scores of all outputs are averaged with uniform weight.

        ``None`` :
            No averaging is performed, an array of scores is returned.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    loss : float or ndarray of shape [n_outputs]
        If output_weights is None, then mean absolute error is returned for
        each output separately.
        If output_weights is 'uniform' or an ndarray of weights, then the
        weighted average of all output scores is returned.

        MAE output is non-negative floating point. The best value is 0.0.

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
    >>> mean_absolute_error(y_true, y_pred, output_weights=None)
    array([ 0.5,  1. ])
    >>> mean_absolute_error(y_true, y_pred, output_weights=[0.3, 0.7])
    ... # doctest: +ELLIPSIS
    0.849...
    """
    y_type, y_true, y_pred, output_weights = _check_reg_targets(
        y_true, y_pred, output_weights)
    output_errors = np.average(np.abs(y_pred - y_true),
                               weights=sample_weight, axis=0)
    if output_weights is None:
        return output_errors
    elif output_weights == 'uniform':
        # pass None as weights to np.average: uniform mean
        output_weights = None

    return np.average(output_errors, weights=output_weights)


def mean_squared_error(y_true, y_pred,
                       output_weights='uniform',
                       sample_weight=None):
    """Mean squared error regression loss

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    output_weights : string in ['uniform'] or None
                     or array-like of shape [n_outputs]
        Weights by which to average scores of outputs. Useful only if using
        multiple outputs.

        ``ndarray`` :
            array containing weights for the weighted average.

        ``'uniform'`` :
            Scores of all outputs are averaged with uniform weight.

        ``None`` :
            No averaging is performed, an array of scores is returned.

    sample_weight : array-like of shape = [n_samples], optional
        Sample weights.

    Returns
    -------
    loss : float
        A non-negative floating point value (the best value is 0.0), or an 
        array of floating point values, one for each individual target.

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
    >>> mean_squared_error(y_true, y_pred, output_weights=None)
    ... # doctest: +ELLIPSIS
    array([ 0.416...,  1.        ])
    >>> mean_squared_error(y_true, y_pred, output_weights=[0.3, 0.7])
    ... # doctest: +ELLIPSIS
    0.824...

    """
    y_type, y_true, y_pred, output_weights = _check_reg_targets(
        y_true, y_pred, output_weights)
    output_errors = np.average((y_true - y_pred) ** 2, axis=0,
                                   weights=sample_weight)
    if output_weights is None:
        return output_errors
    elif output_weights == 'uniform':
        # pass None as weights to np.average: uniform mean
        output_weights = None

    return np.average(output_errors, weights=output_weights)


def explained_variance_score(y_true, y_pred,
                             output_weights='uniform',
                             sample_weight=None):
    """Explained variance regression score function

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like
        Ground truth (correct) target values.

    y_pred : array-like
        Estimated target values.

    output_weights : string in ['uniform', 'variance'] or None
                     or array-like of shape [n_outputs]
        Weights by which to average scores of outputs. Useful only if using
        multiple outputs.

        ``ndarray`` :
            array containing weights for the weighted average.

        ``'uniform'`` :
            Scores of all outputs are averaged with uniform weight.

        ``'variance'`` :
            Scores of all outputs are averaged, weighted by the variances
            of each individual output. This corresponds to a global explained

        ``None`` :
            No averaging is performed, an array of scores is returned

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
    y_type, y_true, y_pred, output_weights = _check_reg_targets(
        y_true, y_pred, output_weights)

    if y_type not in ["continuous", "continuous-multioutput"]:
        raise ValueError("{0} is not supported".format(y_type))

    _, numerator = _average_and_variance(y_true - y_pred, sample_weight,
                                         axis=0)
    _, denominator = _average_and_variance(y_true, sample_weight,
                                           axis=0)

    nonzero_numerator = numerator != 0
    nonzero_denominator = denominator != 0
    valid_score = nonzero_numerator & nonzero_denominator
    output_scores = np.ones(y_true.shape[1])

    output_scores[valid_score] = 1 - (numerator[valid_score] /
                                      denominator[valid_score])
    output_scores[nonzero_numerator & ~nonzero_denominator] = 0.
    if output_weights is None:
        # return scores individually
        return output_scores
    elif output_weights == 'uniform':
        # passing None as weights results is uniform mean
        output_weights = None
    elif output_weights == 'variance':
        output_weights = denominator

    return np.average(output_scores, weights=output_weights)


def r2_score(y_true, y_pred,
             output_weights='uniform',
             sample_weight=None):
    """R^2 (coefficient of determination) regression score function.

    Best possible score is 1.0, lower values are worse.

    Parameters
    ----------
    y_true : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Ground truth (correct) target values.

    y_pred : array-like of shape = [n_samples] or [n_samples, n_outputs]
        Estimated target values.

    output_weights : string in ['uniform', 'variance'] or None
                     or array-like of shape [n_outputs]
        Weights by which to average scores of outputs. Useful only if using
        multiple outputs.

        ``ndarray`` :
            array containing weights for the weighted average.

        ``'uniform'`` :
            Scores of all outputs are averaged with uniform weight.

        ``'variance'`` :
            Scores of all outputs are averaged, weighted by the variances
            of each individual output. This corresponds to a global explained

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
    >>> r2_score(y_true, y_pred, output_weights='variance')  # doctest: +ELLIPSIS
    0.938...

    """
    y_type, y_true, y_pred, output_weights = _check_reg_targets(
        y_true, y_pred, output_weights)

    if sample_weight is not None:
        sample_weight = column_or_1d(sample_weight)
        weight = sample_weight[:, np.newaxis]
    else:
        weight = 1.

    numerator = (weight * (y_true - y_pred) ** 2).sum(axis=0,
                                                      dtype=np.float64)
    denominator = (weight * (y_true - np.average(
        y_true, axis=0, weights=sample_weight)) ** 2).sum(axis=0,
                                                          dtype=np.float64)
    nonzero_denominator = denominator != 0
    nonzero_numerator = numerator != 0
    valid_score = nonzero_denominator & nonzero_numerator
    output_scores = np.ones([y_true.shape[1]])
    output_scores[valid_score] = 1 - (numerator[valid_score] /
                                      denominator[valid_score])
    # arbitrary set to zero to avoid -inf scores, having a constant
    # y_true is not interesting for scoring a regression anyway
    output_scores[nonzero_numerator & ~nonzero_denominator] = 0.
    if output_weights is None:
        # return scores individually
        return output_scores
    elif output_weights == 'uniform':
        # passing None as weights results is uniform mean
        output_weights = None
    elif output_weights == 'variance':
        output_weights = denominator

    return np.average(output_scores, weights=output_weights)
