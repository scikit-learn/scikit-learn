"""Metrics to assess performance on regression task.

Functions named as ``*_score`` return a scalar value to maximize: the higher
the better.

Function named as ``*_error`` or ``*_loss`` return a scalar value to minimize:
the lower the better.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Real

import numpy as np

from sklearn.exceptions import UndefinedMetricWarning
from sklearn.utils._array_api import (
    _average,
    _find_matching_floating_dtype,
    _median,
    get_namespace,
    get_namespace_and_device,
    move_to,
    size,
)
from sklearn.utils._array_api import _xlogy as xlogy
from sklearn.utils._param_validation import Interval, StrOptions, validate_params
from sklearn.utils.stats import _weighted_percentile
from sklearn.utils.validation import (
    _check_sample_weight,
    _num_samples,
    check_array,
    check_consistent_length,
    column_or_1d,
)


def mean_standardized_log_loss(
    y_true, y_pred, y_std, sample_weight=None, multioutput="uniform_average"
):
    """Mean Standardized Log Loss for uncertainty-aware regression models.

    This metric evaluates both the accuracy of predictions and the calibration
    of predicted standard deviations. It is particularly useful for evaluating
    probabilistic models that output both mean predictions and uncertainty estimates.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Ground truth (correct) target values.

    y_pred : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Estimated target values.

    y_std : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Predicted standard deviations.

    sample_weight : array-like of shape (n_samples,) or None
        Sample weights.

    multioutput : array-like of shape (n_outputs) or string in ['raw_values',
        'uniform_average', 'variance_weighted'] or None
        Custom output weights if ``multioutput`` is array-like or
        just the corresponding argument if ``multioutput`` is a
        correct keyword.

    Returns
    -------
    loss : float or ndarray of floats
        The mean standardized log loss. If multioutput is 'raw_values',
        then a numpy.ndarray of shape (n_outputs,) is returned.
        If multioutput is 'uniform_average' or 'variance_weighted', a float
        is returned.

    Notes
    -----
    Mean Standardized Log Loss is defined as:

    .. math::
        \text{MSLL} = \frac{1}{n} \sum_{i=1}^{n}
        \left(\log(\sigma_i^2) + \frac{(y_i - \mu_i)^2}{2\sigma_i^2} \right)

    where :math:`\mu_i` is the predicted mean, :math:`\sigma_i` is the predicted
    standard deviation, and :math:`y_i` is the true value.

    References
    ----------
    .. [1] G. E. H. (2019). "Evaluating Uncertainty in Regression Models".
           arXiv:1904.03228.
    """
    # Validate inputs
    y_type, y_true, y_pred, sample_weight, multioutput = _check_reg_targets(
        y_true, y_pred, sample_weight, multioutput
    )
    
    # Check if y_std is provided and valid
    if y_std is None:
        raise ValueError("y_std must be provided for mean_standardized_log_loss")
    
    y_std = check_array(y_std, ensure_2d=False, dtype="numeric")
    if y_true.shape != y_std.shape:
        raise ValueError(
            "y_true and y_std must have the same shape. "
            f"Got {y_true.shape} and {y_std.shape}"
        )
    
    # Ensure y_std is positive
    if np.any(y_std <= 0):
        raise ValueError("All predicted standard deviations must be positive")
        
    # Calculate MSLL
    xp, _, device = get_namespace_and_device(y_pred)
    y_true, y_pred, y_std = move_to(y_true, y_pred, y_std, xp=xp, device=device)
    
    # Compute log(σ²) + (y - μ)² / (2σ²) for each sample
    squared_error = (y_true - y_pred)**2
    variance = y_std**2
    
    # Use array API for computation
    msll = _average(
        xp.log(variance) + squared_error / (2 * variance),
        axis=0,
        weights=sample_weight,
        multioutput=multioutput
    )
    
    return msll


__ALL__.extend(["mean_standardized_log_loss"])

__ALL__ = [
    "max_error",
    "mean_absolute_error",
    "mean_squared_error",
    "mean_squared_log_error",
    "median_absolute_error",
    "mean_absolute_percentage_error",
    "mean_pinball_loss",
    "r2_score",
    "root_mean_squared_log_error",
    "root_mean_squared_error",
    "explained_variance_score",
    "mean_tweedie_deviance",
    "mean_poisson_deviance",
    "mean_gamma_deviance",
    "d2_tweedie_score",
    "d2_pinball_score",
    "d2_absolute_error_score",
]


def _check_reg_targets(
    y_true, y_pred, sample_weight, multioutput, dtype="numeric", xp=None, device=None
):
    """Check that y_true, y_pred and sample_weight belong to the same regression task.

    To reduce redundancy when calling `_find_matching_floating_dtype`,
    please use `_check_reg_targets_with_floating_dtype` instead.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Ground truth (correct) target values.

    y_pred : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Estimated target values.

    sample_weight : array-like of shape (n_samples,) or None
        Sample weights.

    multioutput : array-like or string in ['raw_values', uniform_average',
        'variance_weighted'] or None
        None is accepted due to backward compatibility of r2_score().

    dtype : str or list, default="numeric"
        the dtype argument passed to check_array.

    xp : module, default=None
        Precomputed array namespace module. When passed, along with `device`,
        typically from a caller that has already performed inspection of its own
        inputs, skips array namespace inspection. If `None`, the inputs are converted
        to the namespace of `y_pred`.

    device : device, default=None
        Precomputed device. When passed, along with `xp`, typically from a caller
        that has already performed inspection of its own inputs, skips device
        inspection and moves everything to that device. If `None`, the inputs are
        converted to the device of `y_pred`.

    Returns
    -------
    type_true : one of {'continuous', continuous-multioutput'}
        The type of the true target data, as output by
        'utils.multiclass.type_of_target'.

    y_true : array-like of shape (n_samples, n_outputs)
        Ground truth (correct) target values.

    y_pred : array-like of shape (n_samples, n_outputs)
        Estimated target values.

    sample_weight : array-like of shape (n_samples,) or None
        Sample weights.

    multioutput : array-like of shape (n_outputs) or string in ['raw_values',
        uniform_average', 'variance_weighted'] or None
        Custom output weights if ``multioutput`` is array-like or
        just the corresponding argument if ``multioutput`` is a
        correct keyword.
    """
    if xp is None or device is None:
        xp, _, device = get_namespace_and_device(y_pred)
    # Move to specified namespace and device, or the namespace and device of `y_pred`
    y_true, y_pred, sample_weight = move_to(
        y_true, y_pred, sample_weight, xp=xp, device=device
    )

    check_consistent_length(y_true, y_pred, sample_weight)
    y_true = check_array(y_true, ensure_2d=False, dtype=dtype)
    y_pred = check_array(y_pred, ensure_2d=False, dtype=dtype)
    if sample_weight is not None:
        sample_weight = _check_sample_weight(sample_weight, y_true, dtype=dtype)

    if y_true.ndim == 1:
        y_true = xp.reshape(y_true, (-1, 1))

    if y_pred.ndim == 1:
        y_pred = xp.reshape(y_pred, (-1, 1))

    if y_true.shape[1] != y_pred.shape[1]:
        raise ValueError(
            "y_true and y_pred have different number of output ({0}!={1})".format(
                y_true.shape[1], y_pred.shape[1]
            )
        )

    n_outputs = y_true.shape[1]
    allowed_multioutput_str = ("raw_values", "uniform_average", "variance_weighted")
    if isinstance(multioutput, str):
        if multioutput not in allowed_multioutput_str:
            raise ValueError(
                "Allowed 'multioutput' string values are {}. "
                "You provided multioutput={!r}".format(
                    allowed_multioutput_str, multioutput
                )
            )
    elif multioutput is not None:
        # Ensure `multioutput` is on the `y_pred` device and namespace.
        multioutput = move_to(multioutput, xp=xp, device=device)
        multioutput = check_array(multioutput, ensure_2d=False)
        if n_outputs == 1:
            raise ValueError("Custom weights are useful only in multi-output cases.")
        elif n_outputs != multioutput.shape[0]:
            raise ValueError(
                "There must be equally many custom weights "
                f"({multioutput.shape[0]}) as outputs ({n_outputs})."
            )
    y_type = "continuous" if n_outputs == 1 else "continuous-multioutput"

    return y_type, y_true, y_pred, sample_weight, multioutput


def _check_reg_targets_with_floating_dtype(
    y_true, y_pred, sample_weight, multioutput, xp=None, device=None
):
    """Ensures y_true, y_pred, and sample_weight correspond to same regression task.

    Extends `_check_reg_targets` by automatically selecting a suitable floating-point
    data type. If `y_pred` has a floating point dtype, it is used. Otherwise, the
    dtype is determined by `_find_matching_floating_dtype` based on `y_pred`.

    Use this private method only when converting inputs to array API-compatibles.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Ground truth (correct) target values.

    y_pred : array-like of shape (n_samples,) or (n_samples, n_outputs)
        Estimated target values.

    sample_weight : array-like of shape (n_samples,)

    multioutput : array-like or string in ['raw_values', 'uniform_average', \
        'variance_weighted'] or None
        None is accepted due to backward compatibility of r2_score().

    xp : module, default=None
        Precomputed array namespace module. When passed, along with `device`,
        typically from a caller that has already performed inspection of its own
        inputs, skips array namespace inspection. If `None`, the inputs are converted
        to the namespace of `y_pred`.

    device : device, default=None
        Precomputed device. When passed, along with `xp`, typically from a caller
        that has already performed inspection of its own inputs, skips device
        inspection and moves everything to that device. If `None`, the inputs are
        converted to the device of `y_pred`.

    Returns
    -------
    type_true : one of {'continuous', 'continuous-multioutput'}
        The type of the true target data, as output by
        'utils.multiclass.type_of_target'.

    y_true : array-like of shape (n_samples, n_outputs)
        Ground truth (correct) target values.

    y_pred : array-like of shape (n_samples, n_outputs)
        Estimated target values.

    sample_weight : array-like of shape (