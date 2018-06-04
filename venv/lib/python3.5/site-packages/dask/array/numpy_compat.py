from __future__ import absolute_import, division, print_function

from distutils.version import LooseVersion

import numpy as np
import warnings

# Taken from scikit-learn:
# https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/utils/fixes.py#L84
try:
    with warnings.catch_warnings():
        if (not np.allclose(np.divide(.4, 1, casting="unsafe"),
                            np.divide(.4, 1, casting="unsafe", dtype=np.float)) or
            not np.allclose(np.divide(1, .5, dtype='i8'), 2) or
           not np.allclose(np.divide(.4, 1), .4)):
            raise TypeError('Divide not working with dtype: '
                            'https://github.com/numpy/numpy/issues/3484')
        divide = np.divide
        ma_divide = np.ma.divide

except TypeError:
    # Divide with dtype doesn't work on Python 3
    def divide(x1, x2, out=None, dtype=None):
        """Implementation of numpy.divide that works with dtype kwarg.

        Temporary compatibility fix for a bug in numpy's version. See
        https://github.com/numpy/numpy/issues/3484 for the relevant issue."""
        x = np.divide(x1, x2, out)
        if dtype is not None:
            x = x.astype(dtype)
        return x

    ma_divide = np.ma.core._DomainedBinaryOperation(divide,
                                                    np.ma.core._DomainSafeDivide(),
                                                    0, 1)

if LooseVersion(np.__version__) < '1.12.0':
    # These functions were added in numpy 1.12.0. For previous versions they
    # are duplicated here

    # See https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory
    def _replace_nan(a, val):
        """
        If `a` is of inexact type, make a copy of `a`, replace NaNs with
        the `val` value, and return the copy together with a boolean mask
        marking the locations where NaNs were present. If `a` is not of
        inexact type, do nothing and return `a` together with a mask of None.

        Note that scalars will end up as array scalars, which is important
        for using the result as the value of the out argument in some
        operations.

        Parameters
        ----------
        a : array-like
            Input array.
        val : float
            NaN values are set to val before doing the operation.

        Returns
        -------
        y : ndarray
            If `a` is of inexact type, return a copy of `a` with the NaNs
            replaced by the fill value, otherwise return `a`.
        mask: {bool, None}
            If `a` is of inexact type, return a boolean mask marking locations of
            NaNs, otherwise return None.

        """
        is_new = not isinstance(a, np.ndarray)
        if is_new:
            a = np.array(a)
        if not issubclass(a.dtype.type, np.inexact):
            return a, None
        if not is_new:
            # need copy
            a = np.array(a, subok=True)

        mask = np.isnan(a)
        np.copyto(a, val, where=mask)
        return a, mask

    def nancumsum(a, axis=None, dtype=None, out=None):
        """
        Return the cumulative sum of array elements over a given axis treating Not a
        Numbers (NaNs) as zero.  The cumulative sum does not change when NaNs are
        encountered and leading NaNs are replaced by zeros.

        Zeros are returned for slices that are all-NaN or empty.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        a : array_like
            Input array.
        axis : int, optional
            Axis along which the cumulative sum is computed. The default
            (None) is to compute the cumsum over the flattened array.
        dtype : dtype, optional
            Type of the returned array and of the accumulator in which the
            elements are summed.  If `dtype` is not specified, it defaults
            to the dtype of `a`, unless `a` has an integer dtype with a
            precision less than that of the default platform integer.  In
            that case, the default platform integer is used.
        out : ndarray, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type will be cast if necessary. See `doc.ufuncs`
            (Section "Output arguments") for more details.

        Returns
        -------
        nancumsum : ndarray.
            A new array holding the result is returned unless `out` is
            specified, in which it is returned. The result has the same
            size as `a`, and the same shape as `a` if `axis` is not None
            or `a` is a 1-d array.

        See Also
        --------
        :func:`numpy.cumsum` : Cumulative sum across array propagating NaNs.
        isnan : Show which elements are NaN.

        Examples
        --------
        >>> np.nancumsum(1) #doctest: +SKIP
        array([1])
        >>> np.nancumsum([1]) #doctest: +SKIP
        array([1])
        >>> np.nancumsum([1, np.nan]) #doctest: +SKIP
        array([ 1.,  1.])
        >>> a = np.array([[1, 2], [3, np.nan]])
        >>> np.nancumsum(a) #doctest: +SKIP
        array([ 1.,  3.,  6.,  6.])
        >>> np.nancumsum(a, axis=0) #doctest: +SKIP
        array([[ 1.,  2.],
               [ 4.,  2.]])
        >>> np.nancumsum(a, axis=1) #doctest: +SKIP
        array([[ 1.,  3.],
               [ 3.,  3.]])

        """
        a, mask = _replace_nan(a, 0)
        return np.cumsum(a, axis=axis, dtype=dtype, out=out)

    def nancumprod(a, axis=None, dtype=None, out=None):
        """
        Return the cumulative product of array elements over a given axis treating Not a
        Numbers (NaNs) as one.  The cumulative product does not change when NaNs are
        encountered and leading NaNs are replaced by ones.

        Ones are returned for slices that are all-NaN or empty.

        .. versionadded:: 1.12.0

        Parameters
        ----------
        a : array_like
            Input array.
        axis : int, optional
            Axis along which the cumulative product is computed.  By default
            the input is flattened.
        dtype : dtype, optional
            Type of the returned array, as well as of the accumulator in which
            the elements are multiplied.  If *dtype* is not specified, it
            defaults to the dtype of `a`, unless `a` has an integer dtype with
            a precision less than that of the default platform integer.  In
            that case, the default platform integer is used instead.
        out : ndarray, optional
            Alternative output array in which to place the result. It must
            have the same shape and buffer length as the expected output
            but the type of the resulting values will be cast if necessary.

        Returns
        -------
        nancumprod : ndarray
            A new array holding the result is returned unless `out` is
            specified, in which case it is returned.

        See Also
        --------
        :func:`numpy.cumprod` : Cumulative product across array propagating NaNs.
        isnan : Show which elements are NaN.

        Examples
        --------
        >>> np.nancumprod(1) #doctest: +SKIP
        array([1])
        >>> np.nancumprod([1]) #doctest: +SKIP
        array([1])
        >>> np.nancumprod([1, np.nan]) #doctest: +SKIP
        array([ 1.,  1.])
        >>> a = np.array([[1, 2], [3, np.nan]])
        >>> np.nancumprod(a) #doctest: +SKIP
        array([ 1.,  2.,  6.,  6.])
        >>> np.nancumprod(a, axis=0) #doctest: +SKIP
        array([[ 1.,  2.],
               [ 3.,  2.]])
        >>> np.nancumprod(a, axis=1) #doctest: +SKIP
        array([[ 1.,  2.],
               [ 3.,  3.]])

        """
        a, mask = _replace_nan(a, 1)
        return np.cumprod(a, axis=axis, dtype=dtype, out=out)


def _make_sliced_dtype(dtype, index):
    if LooseVersion(np.__version__) == LooseVersion("1.14.0"):
        return _make_sliced_dtype_np_ge_14(dtype, index)
    else:
        return _make_sliced_dtype_np_lt_14(dtype, index)


def _make_sliced_dtype_np_ge_14(dtype, index):
    # For https://github.com/numpy/numpy/pull/6053, NumPy >= 1.14
    new = {
        'names': index,
        'formats': [dtype.fields[name][0] for name in index],
        'offsets': [dtype.fields[name][1] for name in index],
        'itemsize': dtype.itemsize,
    }
    return np.dtype(new)


def _make_sliced_dtype_np_lt_14(dtype, index):
    # For numpy < 1.14
    dt = np.dtype([(name, dtype[name]) for name in index])
    return dt


class _Recurser(object):
    """
    Utility class for recursing over nested iterables
    """

    # This was copied almost verbatim from numpy.core.shape_base._Recurser
    # See numpy license at https://github.com/numpy/numpy/blob/master/LICENSE.txt
    # or NUMPY_LICENSE.txt within this directory

    def __init__(self, recurse_if):
        self.recurse_if = recurse_if

    def map_reduce(self,
                   x,
                   f_map=lambda x, **kwargs: x,
                   f_reduce=lambda x, **kwargs: x,
                   f_kwargs=lambda **kwargs: kwargs,
                   **kwargs):
        """
        Iterate over the nested list, applying:
        * ``f_map`` (T -> U) to items
        * ``f_reduce`` (Iterable[U] -> U) to mapped items

        For instance, ``map_reduce([[1, 2], 3, 4])`` is::

            f_reduce([
              f_reduce([
                f_map(1),
                f_map(2)
              ]),
              f_map(3),
              f_map(4)
            ]])


        State can be passed down through the calls with `f_kwargs`,
        to iterables of mapped items. When kwargs are passed, as in
        ``map_reduce([[1, 2], 3, 4], **kw)``, this becomes::

            kw1 = f_kwargs(**kw)
            kw2 = f_kwargs(**kw1)
            f_reduce([
              f_reduce([
                f_map(1), **kw2)
                f_map(2,  **kw2)
              ],      **kw1),
              f_map(3, **kw1),
              f_map(4, **kw1)
            ]],     **kw)
        """
        def f(x, **kwargs):
            if not self.recurse_if(x):
                return f_map(x, **kwargs)
            else:
                next_kwargs = f_kwargs(**kwargs)
                return f_reduce((
                    f(xi, **next_kwargs)
                    for xi in x
                ), **kwargs)
        return f(x, **kwargs)

    def walk(self, x, index=()):
        """
        Iterate over x, yielding (index, value, entering), where

        * ``index``: a tuple of indices up to this point
        * ``value``: equal to ``x[index[0]][...][index[-1]]``. On the first iteration, is
                     ``x`` itself
        * ``entering``: bool. The result of ``recurse_if(value)``
        """
        do_recurse = self.recurse_if(x)
        yield index, x, do_recurse

        if not do_recurse:
            return
        for i, xi in enumerate(x):
            # yield from ...
            for v in self.walk(xi, index + (i,)):
                yield v
