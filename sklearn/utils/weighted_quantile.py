"""
authors: Jasper Roebroek <roebroek.jasper@gmail.com>

The calculation is roughly 10 times as slow as np.quantile (with high number of samples), which
is not terrible as the data needs to be copied and sorted.
"""

import numpy as np
from numpy.lib.function_base import _quantile_is_valid


def weighted_quantile(a, q, weights=None, axis=None, overwrite_input=False, interpolation='linear',
                      keepdims=False, sorted=False):
    """
    Compute the q-th weighted quantile of the data along the specified axis.

    Parameters
    ----------
    a : array-like
        Input array or object that can be converted to an array.
    q : array-like of float
        Quantile or sequence of quantiles to compute, which must be between
        0 and 1 inclusive.
    weights: array-like, optional
        Weights corresponding to a.
    axis : {int, None}, optional
        Axis along which the quantiles are computed. The default is to compute
        the quantile(s) along a flattened version of the array.
    overwrite_input : bool, optional
        If True, then allow the input array `a` to be modified by intermediate
        calculations, to save memory. In this case, the contents of the input
        `a` after this function completes is undefined.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        This optional parameter specifies the interpolation method to
        use when the desired quantile lies between two data points
        ``i < j``:

            * linear: ``i + (j - i) * fraction``, where ``fraction``
              is the fractional part of the index surrounded by ``i``
              and ``j``.
            * lower: ``i``.
            * higher: ``j``.
            * nearest: ``i`` or ``j``, whichever is nearest.
            * midpoint: ``(i + j) / 2``.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left in
        the result as dimensions with size one. With this option, the
        result will broadcast correctly against the original array `a`.
    sorted : bool, optional
        If the `a` is already sorted along the given axis this can be set to
        True, to avoid the sorting step.

    Returns
    -------
    quantile : scalar or ndarray
        If `q` is a single quantile and `axis=None`, then the result
        is a scalar. If multiple quantiles are given, first axis of
        the result corresponds to the quantiles. The other axes are
        the axes that remain after the reduction of `a`. The output
        dtype is ``float64``.

    References
    ----------
    1. https://en.wikipedia.org/wiki/Percentile#The_Weighted_Percentile_method
    """
    q = np.atleast_1d(q)
    if not _quantile_is_valid(q):
        raise ValueError("Quantiles must be in the range [0, 1]")

    if q.ndim > 2:
        raise ValueError("q must be a scalar or 1D")

    if weights is None:
        return np.quantile(a, q, axis=-1, keepdims=keepdims)
    else:
        # a needs to be able to store NaN-values, thus it needs to be casted to float
        a = np.asarray(a)
        weights = np.asarray(weights)
        if a.shape[:-1] == weights.shape:
            return np.quantile(a, q, axis=axis, keepdims=keepdims)
        elif a.shape != weights.shape:
            raise IndexError("the data and weights need to be of the same shape")

    a = a.astype(np.float64, copy=not overwrite_input)
    if axis is None:
        a = a.ravel()
        weights = weights.ravel()
    elif isinstance(axis, (tuple, list)):
        raise NotImplementedError("Several axes are currently not supported.")
    else:
        a = np.moveaxis(a, axis, 0)
        weights = np.moveaxis(weights, axis, 0)

    q = np.expand_dims(q, axis=list(np.arange(1, a.ndim+1)))

    zeros = weights == 0
    a[zeros] = np.nan
    zeros_count = zeros.sum(axis=0, keepdims=True)

    if not sorted:
        # NaN-values will be sorted to the last places along the axis
        idx_sorted = np.argsort(a, axis=0)
        a_sorted = np.take_along_axis(a, idx_sorted, axis=0)
        weights_sorted = np.take_along_axis(weights, idx_sorted, axis=0)
    else:
        a_sorted = a
        weights_sorted = weights

    weights_cum = np.cumsum(weights_sorted, axis=0)
    weights_total = np.expand_dims(np.take(weights_cum, -1, axis=0), axis=0)

    weights_norm = (weights_cum - 0.5 * weights_sorted) / weights_total
    indices = np.sum(weights_norm < q, axis=1, keepdims=True) - 1

    idx_low = (indices == -1)
    high = a.shape[0] - zeros_count - 1
    idx_high = (indices == high)

    indices = np.clip(indices, 0, high - 1)

    left_weight = np.take_along_axis(weights_norm[np.newaxis, ...], indices, axis=1)
    right_weight = np.take_along_axis(weights_norm[np.newaxis, ...], indices + 1, axis=1)
    left_value = np.take_along_axis(a_sorted[np.newaxis, ...], indices, axis=1)
    right_value = np.take_along_axis(a_sorted[np.newaxis, ...], indices + 1, axis=1)

    if interpolation == 'linear':
        fraction = (q - left_weight) / (right_weight - left_weight)
    elif interpolation == 'lower':
        fraction = 0
    elif interpolation == 'higher':
        fraction = 1
    elif interpolation == 'midpoint':
        fraction = 0.5
    elif interpolation == 'nearest':
        fraction = (np.abs(left_weight - q) > np.abs(right_weight - q))
    else:
        raise ValueError("interpolation should be one of: {'linear', 'lower', 'higher', 'midpoint', 'nearest'}")

    quantiles = left_value + fraction * (right_value - left_value)

    if idx_low.sum() > 0:
        quantiles[idx_low] = np.take(a_sorted, 0, axis=0).flatten()
    if idx_high.sum() > 0:
        quantiles[idx_high] = np.take_along_axis(a_sorted, high, axis=0).flatten()

    if q.size == 1:
        quantiles = quantiles[0]

    if keepdims:
        return quantiles
    else:
        if quantiles.size == 1:
            return quantiles.item()
        else:
            return quantiles.squeeze()
