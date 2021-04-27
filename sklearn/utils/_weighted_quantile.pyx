# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# distutils: language = c

cimport numpy as np
import numpy as np
from numpy.lib.function_base import _quantile_is_valid
from libc.math cimport isnan
from libc.stdlib cimport malloc, free


cdef extern from "stdlib.h":
    ctypedef void const_void "const void"
    void qsort(void *base, int nmemb, int size,
            int(*compar)(const_void *, const_void *)) nogil


cdef struct IndexedElement:
    np.ulong_t index
    np.float32_t value


cdef int _compare(const_void *a, const_void *b):
    cdef np.float32_t v
    if isnan((<IndexedElement*> a).value): return 1
    if isnan((<IndexedElement*> b).value): return -1
    v = (<IndexedElement*> a).value-(<IndexedElement*> b).value
    if v < 0: return -1
    if v >= 0: return 1


cdef long[:] argsort(float[:] data) nogil:
    """source: https://github.com/jcrudy/cython-argsort/blob/master/cyargsort/argsort.pyx"""
    cdef np.ulong_t i
    cdef np.ulong_t n = data.shape[0]
    cdef long[:] order

    with gil:
        order = np.empty(n, dtype=np.int_)

    # Allocate index tracking array.
    cdef IndexedElement *order_struct = <IndexedElement *> malloc(n * sizeof(IndexedElement))

    # Copy data into index tracking array.
    for i in range(n):
        order_struct[i].index = i
        order_struct[i].value = data[i]

    # Sort index tracking array.
    qsort(<void *> order_struct, n, sizeof(IndexedElement), _compare)

    # Copy indices from index tracking array to output array.
    for i in range(n):
        order[i] = order_struct[i].index

    # Free index tracking array.
    free(order_struct)

    return order


cdef int _searchsorted1D(float[:] A, float x) nogil:
    """
    source: https://github.com/gesellkammer/numpyx/blob/master/numpyx.pyx
    """
    cdef:
        int imin = 0
        int imax = A.shape[0]
        int imid
    while imin < imax:
        imid = imin + ((imax - imin) / 2)
        if A[imid] < x:
            imin = imid + 1
        else:
            imax = imid
    return imin


cdef void _weighted_quantile_presorted_1D(float[:] a,
                                          float[:] q,
                                          float[:] weights,
                                          float[:] quantiles,
                                          Interpolation interpolation) nogil:
    """
    Weighted quantile (1D) on presorted data. 
    Note: the weights data will be changed
    """
    cdef long q_idx
    cdef float weights_total, weights_cum, frac
    cdef int i

    cdef int n_samples = a.shape[0]
    cdef int n_q = q.shape[0]

    weights_total = 0
    for i in range(n_samples):
        weights_total += weights[i]

    weights_cum = weights[0]
    weights[0] = 0.5 * weights[0] / weights_total
    for i in range(1, n_samples):
        weights_cum += weights[i]
        weights[i] = (weights_cum - 0.5 * weights[i]) / weights_total

    for i in range(n_q):
        q_idx = _searchsorted1D(weights, q[i]) - 1

        if q_idx == -1:
            quantiles[i] = a[0]
        elif q_idx == n_samples - 1:
            quantiles[i] = a[n_samples - 1]
        else:
            quantiles[i] = a[q_idx]
            if interpolation == linear:
                frac = (q[i] - weights[q_idx]) / (weights[q_idx + 1] - weights[q_idx])
            elif interpolation == lower:
                frac = 0
            elif interpolation == higher:
                frac = 1
            elif interpolation == midpoint:
                frac = 0.5
            elif interpolation == nearest:
                frac = (q[i] - weights[q_idx]) / (weights[q_idx + 1] - weights[q_idx])
                if frac < 0.5:
                    frac = 0
                else:
                    frac = 1

            quantiles[i] = a[q_idx] + frac * (a[q_idx + 1] - a[q_idx])


cdef void _weighted_quantile_unchecked_1D(float[:] a,
                                          float[:] q,
                                          float[:] weights,
                                          float[:] quantiles,
                                          Interpolation interpolation) nogil:
    """
    Weighted quantile (1D)
    Note: the data is not guaranteed to not be changed within this function
    """
    cdef long[:] sort_idx
    cdef int n_samples = a.shape[0]
    cdef long count_samples = 0
    cdef float[:] a_processed
    cdef float[:] weights_processed
    cdef int i

    for i in range(n_samples):
        if isnan(a[i]):
            continue
        elif weights[i] == 0:
            continue
        else:
            a[count_samples] = a[i]
            weights[count_samples] = weights[i]
            count_samples += 1

    sort_idx = argsort(a[:count_samples])

    with gil:
        a_processed = np.empty(count_samples, dtype=np.float32)
        weights_processed = np.empty(count_samples, dtype=np.float32)

    for i in range(count_samples):
        a_processed[i] = a[sort_idx[i]]
        weights_processed[i] = weights[sort_idx[i]]

    _weighted_quantile_presorted_1D(a_processed[:count_samples], q, weights_processed[:count_samples],
                                    quantiles, interpolation)


def _weighted_quantile_unchecked(a, q, weights, axis, overwrite_input=False, interpolation='linear',
                                 keepdims=False):
    """
    Numpy implementation
    Axis should not be none and a should have more than 1 dimension
    This implementation is faster than doing it manually in cython (as the
    looping currently happens with the GIL)
    """
    a = np.asarray(a, dtype=np.float32)
    weights = np.asarray(weights, dtype=np.float32)
    q = np.asarray(q, dtype=np.float32)

    a = np.moveaxis(a, axis, 0)
    weights = np.moveaxis(weights, axis, 0)

    q = np.expand_dims(q, axis=list(np.arange(1, a.ndim+1)))

    zeros = weights == 0
    a[zeros] = np.nan
    zeros_count = zeros.sum(axis=0, keepdims=True)

    idx_sorted = np.argsort(a, axis=0)
    a_sorted = np.take_along_axis(a, idx_sorted, axis=0)
    weights_sorted = np.take_along_axis(weights, idx_sorted, axis=0)

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

    return quantiles


def weighted_quantile(a, q, weights=None, axis=None, overwrite_input=False, interpolation='linear',
                      keepdims=False):
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
        return np.quantile(a, q, axis=-1, keepdims=keepdims, overwrite_input=overwrite_input,
                           interpolation=interpolation)
    else:
        a = np.asarray(a, dtype=np.float32)
        weights = np.asarray(weights, dtype=np.float32)

        if not overwrite_input:
            a = a.copy()
            weights = weights.copy()

        if a.shape != weights.shape:
            raise IndexError("the data and weights need to be of the same shape")

        q = q.astype(np.float32)

    if interpolation == 'linear':
        c_interpolation = linear
    elif interpolation == 'lower':
        c_interpolation = lower
    elif interpolation == 'higher':
        c_interpolation = higher
    elif interpolation == 'midpoint':
        c_interpolation = midpoint
    elif interpolation == 'nearest':
        c_interpolation = nearest
    else:
        raise ValueError("interpolation should be one of: {'linear', 'lower', 'higher', 'midpoint', 'nearest'}")

    if isinstance(axis, (tuple, list)):
        raise NotImplementedError("Several axes are currently not supported.")

    elif axis is not None and a.ndim > 1:
        quantiles = _weighted_quantile_unchecked(a, q, weights, axis, interpolation=interpolation,
                                                 keepdims=keepdims)

    else:
        a = a.ravel()
        weights = weights.ravel()
        quantiles = np.empty(q.size, dtype=np.float32)
        _weighted_quantile_unchecked_1D(a, q, weights, quantiles, c_interpolation)

    if q.size == 1:
        quantiles = quantiles[0]
        start_axis = 0
    else:
        start_axis = 1

    if keepdims:
        if a.ndim > 1:
            quantiles = np.moveaxis(quantiles, 0, axis)
        return quantiles
    else:
        if quantiles.size == 1:
            return quantiles.item()
        else:
            if quantiles.ndim == 1:
                return quantiles
            else:
                return np.take(quantiles, 0, axis=start_axis)
