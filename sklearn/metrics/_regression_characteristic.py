"""Regression Error Characteristic curve"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numbers  # For type checking Python scalars

from ..utils import check_array, check_consistent_length
from ..utils._array_api import get_namespace_and_device  # For array_api support
from ..utils._param_validation import StrOptions, validate_params


@validate_params(
    {
        "y_true": ["array-like"],
        "y_pred": ["array-like", numbers.Real],  # Allow scalar or array-like
        "loss": [StrOptions({"absolute", "squared"})],
    },
    prefer_skip_nested_validation=True,
)
def rec_curve(y_true, y_pred, *, loss="absolute"):
    """Compute Regression Error Characteristic (REC) curve.

    The REC curve evaluates regression models by plotting the error tolerance
    (deviation) on the x-axis against the percentage of data points predicted
    within that tolerance (accuracy) on the y-axis. It is the empirical
    Cumulative Distribution Function (CDF) of the error.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,)
        True target values.

    y_pred : array-like of shape (n_samples,) or scalar
        Estimated target values. If a scalar is provided, it is treated as
        a constant prediction for all samples.

    loss : {'absolute', 'squared'}, default='absolute'
        The type of loss to use for calculating deviations.
        - 'absolute': Uses absolute deviations |y_true - y_pred|.
        - 'squared': Uses squared deviations (y_true - y_pred)^2.

    Returns
    -------
    deviations : ndarray
        Sorted unique error tolerance values. These are the x-coordinates
        for the REC curve. The array will start with 0.0 if the smallest
        calculated error is greater than 0, representing zero tolerance.

    accuracy : ndarray
        The corresponding accuracy (fraction of samples with error less than
        or equal to the deviation). These are the y-coordinates for the REC
        curve. The array will start with 0.0 if `deviations` starts with an
        explicit 0.0.

    See Also
    --------
    roc_curve : Compute Receiver Operating Characteristic (ROC) curve.
    det_curve : Compute Detection Error Tradeoff (DET) curve.

    References
    ----------
    .. [1] Bi, J., & Bennett, K. P. (2003). Regression error characteristic
           curves. In Proceedings of the 20th International Conference on
           Machine Learning (ICML-03) (pp. 43-50).

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.metrics import rec_curve # Assuming function is in sklearn.metrics
    >>> y_true = np.array([1, 2, 3, 4, 5, 6])
    >>> y_pred_model1 = np.array([1.1, 2.2, 2.8, 4.3, 4.8, 6.5])
    >>> deviations, accuracy = rec_curve(y_true, y_pred_model1, loss='absolute')
    >>> deviations
    array([0.  , 0.1 , 0.2 , 0.3 , 0.5 ])
    >>> accuracy
    array([0.        , 0.16666667, 0.66666667, 0.83333333, 1.        ])

    >>> # Example with a scalar prediction (constant model)
    >>> y_pred_scalar = 3.5
    >>> dev_scalar, acc_scalar = rec_curve(y_true, y_pred_scalar)
    >>> dev_scalar
    array([0. , 0.5, 1.5, 2.5])
    >>> acc_scalar
    array([0.        , 0.33333333, 0.66666667, 1.        ])

    >>> # Example with squared loss
    >>> dev_sq, acc_sq = rec_curve(y_true, y_pred_model1, loss='squared')
    >>> dev_sq # These are squared errors
    array([0.    , 0.01  , 0.04  , 0.09  , 0.25  ])
    >>> acc_sq
    array([0.        , 0.16666667, 0.66666667, 0.83333333, 1.        ])

    >>> # For plotting with matplotlib:
    >>> # import matplotlib.pyplot as plt
    >>> # plt.figure()
    >>> # plt.plot(deviations, accuracy, marker='.', label='Model 1 (Absolute Loss)')
    >>> # plt.plot(dev_scalar, acc_scalar, marker='.', label='Constant (Absolute Loss)')
    >>> # plt.xlabel("Error Tolerance (Deviation)")
    >>> # plt.ylabel("Accuracy (Fraction of samples within tolerance)")
    >>> # plt.title("Regression Error Characteristic (REC) Curve")
    >>> # plt.legend()
    >>> # plt.grid(True)
    >>> # plt.show()
    """
    # Validate y_true and get the array namespace (xp)
    y_true_array = check_array(
        y_true, ensure_2d=False, dtype="numeric", ensure_all_finite=True
    )
    xp, _, device = get_namespace_and_device(y_true_array)

    # Handle Python native scalars (int, float)
    if isinstance(y_pred, numbers.Number):  # numbers.Real covers int, float
        y_pred_scalar_val = float(y_pred)
        y_pred = xp.full(
            y_true_array.shape,
            fill_value=y_pred_scalar_val,
            dtype=y_true_array.dtype,  # Match y_true's dtype for consistency
            device=device,
        )

    # Handle array-like with a single prediction
    if y_pred.size == 1:
        y_pred = xp.squeeze(y_pred)
        y_pred = xp.tile(y_pred, y_true_array.shape)

    y_pred_array = check_array(
        y_pred, ensure_2d=False, dtype="numeric", ensure_all_finite=True
    )
    check_consistent_length(y_true_array, y_pred_array)

    # Calculate deviations based on the chosen loss
    differences = y_true_array - y_pred_array
    if loss == "absolute":
        errors = xp.abs(differences)
    else:  # loss == "squared"
        errors = xp.square(differences)

    n_samples = y_true_array.shape[0]

    # Handle empty input (no samples)
    if n_samples == 0:
        empty_float_array = xp.asarray(
            [], dtype=xp.float64, device=xp.device(y_true_array)
        )
        return empty_float_array, empty_float_array

    # OPTIMIZED CDF CALCULATION:
    # Get unique sorted error values (deviations_calc) and their counts.
    # xp.unique_counts returns sorted unique values.
    deviations_calc, counts = xp.unique_counts(errors)

    # Calculate cumulative accuracy
    cumulative_counts = xp.cumsum(counts)
    accuracy_values = xp.astype(cumulative_counts, xp.float64) / float(n_samples)

    # Prepare output deviations and accuracy
    # Prepend (0,0) if the smallest error (first element of deviations_calc) is > 0.0,
    # ensuring the curve starts from the origin of the plot unless
    # there are samples with exactly zero error.
    if deviations_calc[0] > 0.0:
        zero_dev = xp.asarray([0.0], dtype=deviations_calc.dtype, device=device)
        zero_acc = xp.asarray([0.0], dtype=accuracy_values.dtype, device=device)

        deviations_out = xp.concatenate((zero_dev, deviations_calc))
        accuracy_out = xp.concatenate((zero_acc, accuracy_values))
    else:
        deviations_out = deviations_calc
        accuracy_out = accuracy_values

    return deviations_out, accuracy_out
