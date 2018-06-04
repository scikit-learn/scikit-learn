"""
The arraycrop module contains functions to crop values from the edges of an
n-dimensional array.
"""
from __future__ import division, absolute_import, print_function

import numpy as np
from numpy.lib.arraypad import _validate_lengths


__all__ = ['crop']


# The below functions are retained in comments in case the NumPy architecture
# changes and we need copies of these helper functions for `crop`.
# These are identical to functions in numpy.lib.arraypad.py as of NumPy v1.11

# def _normalize_shape(ndarray, shape, cast_to_int=True):
#     """
#     Private function which does some checks and normalizes the possibly
#     much simpler representations of 'pad_width', 'stat_length',
#     'constant_values', 'end_values'.
#
#     Parameters
#     ----------
#     narray : ndarray
#         Input ndarray
#     shape : {sequence, array_like, float, int}, optional
#         The width of padding (pad_width), the number of elements on the
#         edge of the narray used for statistics (stat_length), the constant
#         value(s) to use when filling padded regions (constant_values), or the
#         endpoint target(s) for linear ramps (end_values).
#         ((before_1, after_1), ... (before_N, after_N)) unique number of
#         elements for each axis where `N` is rank of `narray`.
#         ((before, after),) yields same before and after constants for each
#         axis.
#         (constant,) or val is a shortcut for before = after = constant for
#         all axes.
#     cast_to_int : bool, optional
#         Controls if values in ``shape`` will be rounded and cast to int
#         before being returned.
#
#     Returns
#     -------
#     normalized_shape : tuple of tuples
#         val                               => ((val, val), (val, val), ...)
#         [[val1, val2], [val3, val4], ...] => ((val1, val2), (val3, val4), ...)
#         ((val1, val2), (val3, val4), ...) => no change
#         [[val1, val2], ]                  => ((val1, val2), (val1, val2), ...)
#         ((val1, val2), )                  => ((val1, val2), (val1, val2), ...)
#         [[val ,     ], ]                  => ((val, val), (val, val), ...)
#         ((val ,     ), )                  => ((val, val), (val, val), ...)
#
#     """
#     ndims = ndarray.ndim
#
#     # Shortcut shape=None
#     if shape is None:
#         return ((None, None), ) * ndims
#
#     # Convert any input `info` to a NumPy array
#     arr = np.asarray(shape)
#
#     # Switch based on what input looks like
#     if arr.ndim <= 1:
#         if arr.shape == () or arr.shape == (1,):
#             # Single scalar input
#             #   Create new array of ones, multiply by the scalar
#             arr = np.ones((ndims, 2), dtype=ndarray.dtype) * arr
#         elif arr.shape == (2,):
#             # Apply padding (before, after) each axis
#             #   Create new axis 0, repeat along it for every axis
#             arr = arr[np.newaxis, :].repeat(ndims, axis=0)
#         else:
#             fmt = "Unable to create correctly shaped tuple from %s"
#             raise ValueError(fmt % (shape,))
#
#     elif arr.ndim == 2:
#         if arr.shape[1] == 1 and arr.shape[0] == ndims:
#             # Padded before and after by the same amount
#             arr = arr.repeat(2, axis=1)
#         elif arr.shape[0] == ndims:
#             # Input correctly formatted, pass it on as `arr`
#             arr = shape
#         else:
#             fmt = "Unable to create correctly shaped tuple from %s"
#             raise ValueError(fmt % (shape,))
#
#     else:
#         fmt = "Unable to create correctly shaped tuple from %s"
#         raise ValueError(fmt % (shape,))
#
#     # Cast if necessary
#     if cast_to_int is True:
#         arr = np.round(arr).astype(int)
#
#     # Convert list of lists to tuple of tuples
#     return tuple(tuple(axis) for axis in arr.tolist())


# def _validate_lengths(narray, number_elements):
#     """
#     Private function which does some checks and reformats pad_width and
#     stat_length using _normalize_shape.
#
#     Parameters
#     ----------
#     narray : ndarray
#         Input ndarray
#     number_elements : {sequence, int}, optional
#         The width of padding (pad_width) or the number of elements on the edge
#         of the narray used for statistics (stat_length).
#         ((before_1, after_1), ... (before_N, after_N)) unique number of
#         elements for each axis.
#         ((before, after),) yields same before and after constants for each
#         axis.
#         (constant,) or int is a shortcut for before = after = constant for all
#         axes.
#
#     Returns
#     -------
#     _validate_lengths : tuple of tuples
#         int                               => ((int, int), (int, int), ...)
#         [[int1, int2], [int3, int4], ...] => ((int1, int2), (int3, int4), ...)
#         ((int1, int2), (int3, int4), ...) => no change
#         [[int1, int2], ]                  => ((int1, int2), (int1, int2), ...)
#         ((int1, int2), )                  => ((int1, int2), (int1, int2), ...)
#         [[int ,     ], ]                  => ((int, int), (int, int), ...)
#         ((int ,     ), )                  => ((int, int), (int, int), ...)
#
#     """
#     normshp = _normalize_shape(narray, number_elements)
#     for i in normshp:
#         chk = [1 if x is None else x for x in i]
#         chk = [1 if x >= 0 else -1 for x in chk]
#         if (chk[0] < 0) or (chk[1] < 0):
#             fmt = "%s cannot contain negative values."
#             raise ValueError(fmt % (number_elements,))
#     return normshp


def crop(ar, crop_width, copy=False, order='K'):
    """Crop array `ar` by `crop_width` along each dimension.

    Parameters
    ----------
    ar : array-like of rank N
        Input array.
    crop_width : {sequence, int}
        Number of values to remove from the edges of each axis.
        ``((before_1, after_1),`` ... ``(before_N, after_N))`` specifies
        unique crop widths at the start and end of each axis.
        ``((before, after),)`` specifies a fixed start and end crop
        for every axis.
        ``(n,)`` or ``n`` for integer ``n`` is a shortcut for
        before = after = ``n`` for all axes.
    copy : bool, optional
        If `True`, ensure the returned array is a contiguous copy. Normally,
        a crop operation will return a discontiguous view of the underlying
        input array.
    order : {'C', 'F', 'A', 'K'}, optional
        If ``copy==True``, control the memory layout of the copy. See
        ``np.copy``.

    Returns
    -------
    cropped : array
        The cropped array. If ``copy=False`` (default), this is a sliced
        view of the input array.
    """
    ar = np.array(ar, copy=False)
    crops = _validate_lengths(ar, crop_width)
    slices = [slice(a, ar.shape[i] - b) for i, (a, b) in enumerate(crops)]
    if copy:
        cropped = np.array(ar[slices], order=order, copy=True)
    else:
        cropped = ar[slices]
    return cropped
