# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""
robust_minmax: Helps scales data to the [0,1] range using robust percentile.

This method is more stable and handle outliers excellently.
"""

import numpy as np
from sklearn.utils.validation import check_array


def robust_minmax_scale(X, *, quantile_range=(25.0, 75.0)):
    """
    Scale features to [0, 1] using robust quantiles instead of min/max.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Input data.
    quantile_range : tuple (lower, upper)
        Percentile range used for robust scaling.

    Returns
    -------
    X_scaled : ndarray
        Robustly min-max scaled array.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import robust_minmax_scale
    >>> X = np.array([[1], [5], [100]])
    >>> robust_minmax_scale(X)
    array([[0.  ],
           [0.04],
           [1.  ]])
    """
    X = check_array(X, accept_sparse=False)

    lower, upper = quantile_range
    q_low = np.percentile(X, lower, axis=0)
    q_high = np.percentile(X, upper, axis=0)

    scale = q_high - q_low
    scale[scale == 0] = 1.0  # avoid divide-by-zero

    X_scaled = (X - q_low) / scale
    X_scaled = np.clip(X_scaled, 0, 1)

    return X_scaled
