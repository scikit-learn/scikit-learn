"""Compatibility fixes for older version of python, numpy and scipy

If you add content to this file, please give the version of the package
at which the fixe is no longer needed.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Fabian Pedregosa <fpedregosa@acm.org>
#          Lars Buitinck
#
# License: BSD 3 clause

import inspect
import warnings

import numpy as np
import scipy.sparse as sp

np_version = []
for x in np.__version__.split('.'):
    try:
        np_version.append(int(x))
    except ValueError:
        # x may be of the form dev-1ea1592
        np_version.append(x)
np_version = tuple(np_version)


try:
    from scipy.special import expit     # SciPy >= 0.10
    with np.errstate(invalid='ignore', over='ignore'):
        if np.isnan(expit(1000)):       # SciPy < 0.14
            raise ImportError("no stable expit in scipy.special")
except ImportError:
    def expit(x, out=None):
        """Logistic sigmoid function, ``1 / (1 + exp(-x))``.

        See sklearn.utils.extmath.log_logistic for the log of this function.
        """
        if out is None:
            out = np.empty(np.atleast_1d(x).shape, dtype=np.float64)
        out[:] = x

        # 1 / (1 + exp(-x)) = (1 + tanh(x / 2)) / 2
        # This way of computing the logistic is both fast and stable.
        out *= .5
        np.tanh(out, out)
        out += 1
        out *= .5

        return out.reshape(np.shape(x))


# little danse to see if np.copy has an 'order' keyword argument
if 'order' in inspect.getargspec(np.copy)[0]:
    def safe_copy(X):
        # Copy, but keep the order
        return np.copy(X, order='K')
else:
    # Before an 'order' argument was introduced, numpy wouldn't muck with
    # the ordering
    safe_copy = np.copy

try:
    if (not np.allclose(np.divide(.4, 1, casting="unsafe"),
                        np.divide(.4, 1, casting="unsafe", dtype=np.float))
            or not np.allclose(np.divide(.4, 1), .4)):
        raise TypeError('Divide not working with dtype: '
                        'https://github.com/numpy/numpy/issues/3484')
    divide = np.divide

except TypeError:
    # Compat for old versions of np.divide that do not provide support for
    # the dtype args
    def divide(x1, x2, out=None, dtype=None):
        out_orig = out
        if out is None:
            out = np.asarray(x1, dtype=dtype)
            if out is x1:
                out = x1.copy()
        else:
            if out is not x1:
                out[:] = x1
        if dtype is not None and out.dtype != dtype:
            out = out.astype(dtype)
        out /= x2
        if out_orig is None and np.isscalar(x1):
            out = np.asscalar(out)
        return out


try:
    np.array(5).astype(float, copy=False)
except TypeError:
    # Compat where astype accepted no copy argument
    def astype(array, dtype, copy=True):
        if array.dtype == dtype:
            return array
        return array.astype(dtype)
else:
    astype = np.ndarray.astype


try:
    with warnings.catch_warnings(record=True):
        # Don't raise the numpy deprecation warnings that appear in
        # 1.9, but avoid Python bug due to simplefilter('ignore')
        warnings.simplefilter('always')
        sp.csr_matrix([1.0, 2.0, 3.0]).max(axis=0)
except (TypeError, AttributeError):
    # in scipy < 14.0, sparse matrix min/max doesn't accept an `axis` argument
    # the following code is taken from the scipy 0.14 codebase

    def _minor_reduce(X, ufunc):
        major_index = np.flatnonzero(np.diff(X.indptr))
        if X.data.size == 0 and major_index.size == 0:
            # Numpy < 1.8.0 don't handle empty arrays in reduceat
            value = np.zeros_like(X.data)
        else:
            value = ufunc.reduceat(X.data, X.indptr[major_index])
        return major_index, value

    def _min_or_max_axis(X, axis, min_or_max):
        N = X.shape[axis]
        if N == 0:
            raise ValueError("zero-size array to reduction operation")
        M = X.shape[1 - axis]
        mat = X.tocsc() if axis == 0 else X.tocsr()
        mat.sum_duplicates()
        major_index, value = _minor_reduce(mat, min_or_max)
        not_full = np.diff(mat.indptr)[major_index] < N
        value[not_full] = min_or_max(value[not_full], 0)
        mask = value != 0
        major_index = np.compress(mask, major_index)
        value = np.compress(mask, value)

        from scipy.sparse import coo_matrix
        if axis == 0:
            res = coo_matrix((value, (np.zeros(len(value)), major_index)),
                             dtype=X.dtype, shape=(1, M))
        else:
            res = coo_matrix((value, (major_index, np.zeros(len(value)))),
                             dtype=X.dtype, shape=(M, 1))
        return res.A.ravel()

    def _sparse_min_or_max(X, axis, min_or_max):
        if axis is None:
            if 0 in X.shape:
                raise ValueError("zero-size array to reduction operation")
            zero = X.dtype.type(0)
            if X.nnz == 0:
                return zero
            m = min_or_max.reduce(X.data.ravel())
            if X.nnz != np.product(X.shape):
                m = min_or_max(zero, m)
            return m
        if axis < 0:
            axis += 2
        if (axis == 0) or (axis == 1):
            return _min_or_max_axis(X, axis, min_or_max)
        else:
            raise ValueError("invalid axis, use 0 for rows, or 1 for columns")

    def sparse_min_max(X, axis):
        return (_sparse_min_or_max(X, axis, np.minimum),
                _sparse_min_or_max(X, axis, np.maximum))

else:
    def sparse_min_max(X, axis):
        return (X.min(axis=axis).toarray().ravel(),
                X.max(axis=axis).toarray().ravel())


try:
    from numpy import argpartition
except ImportError:
    # numpy.argpartition was introduced in v 1.8.0
    def argpartition(a, kth, axis=-1, kind='introselect', order=None):
        return np.argsort(a, axis=axis, order=order)


try:
    from itertools import combinations_with_replacement
except ImportError:
    # Backport of itertools.combinations_with_replacement for Python 2.6,
    # from Python 3.4 documentation (http://tinyurl.com/comb-w-r), copyright
    # Python Software Foundation (https://docs.python.org/3/license.html)
    def combinations_with_replacement(iterable, r):
        # combinations_with_replacement('ABC', 2) --> AA AB AC BB BC CC
        pool = tuple(iterable)
        n = len(pool)
        if not n and r:
            return
        indices = [0] * r
        yield tuple(pool[i] for i in indices)
        while True:
            for i in reversed(range(r)):
                if indices[i] != n - 1:
                    break
            else:
                return
            indices[i:] = [indices[i] + 1] * (r - i)
            yield tuple(pool[i] for i in indices)


try:
    from numpy import isclose
except ImportError:
    def isclose(a, b, rtol=1.e-5, atol=1.e-8, equal_nan=False):
        """
        Returns a boolean array where two arrays are element-wise equal within
        a tolerance.

        This function was added to numpy v1.7.0, and the version you are
        running has been backported from numpy v1.8.1. See its documentation
        for more details.
        """
        def within_tol(x, y, atol, rtol):
            with np.errstate(invalid='ignore'):
                result = np.less_equal(abs(x-y), atol + rtol * abs(y))
            if np.isscalar(a) and np.isscalar(b):
                result = bool(result)
            return result

        x = np.array(a, copy=False, subok=True, ndmin=1)
        y = np.array(b, copy=False, subok=True, ndmin=1)
        xfin = np.isfinite(x)
        yfin = np.isfinite(y)
        if all(xfin) and all(yfin):
            return within_tol(x, y, atol, rtol)
        else:
            finite = xfin & yfin
            cond = np.zeros_like(finite, subok=True)
            # Since we're using boolean indexing, x & y must be the same shape.
            # Ideally, we'd just do x, y = broadcast_arrays(x, y). It's in
            # lib.stride_tricks, though, so we can't import it here.
            x = x * np.ones_like(cond)
            y = y * np.ones_like(cond)
            # Avoid subtraction with infinite/nan values...
            cond[finite] = within_tol(x[finite], y[finite], atol, rtol)
            # Check for equality of infinite values...
            cond[~finite] = (x[~finite] == y[~finite])
            if equal_nan:
                # Make NaN == NaN
                cond[np.isnan(x) & np.isnan(y)] = True
            return cond


if np_version < (1, 7):
    # Prior to 1.7.0, np.frombuffer wouldn't work for empty first arg.
    def frombuffer_empty(buf, dtype):
        if len(buf) == 0:
            return np.empty(0, dtype=dtype)
        else:
            return np.frombuffer(buf, dtype=dtype)
else:
    frombuffer_empty = np.frombuffer


if np_version < (1, 8):
    def in1d(ar1, ar2, assume_unique=False, invert=False):
        # Backport of numpy function in1d 1.8.1 to support numpy 1.6.2
        # Ravel both arrays, behavior for the first array could be different
        ar1 = np.asarray(ar1).ravel()
        ar2 = np.asarray(ar2).ravel()

        # This code is significantly faster when the condition is satisfied.
        if len(ar2) < 10 * len(ar1) ** 0.145:
            if invert:
                mask = np.ones(len(ar1), dtype=np.bool)
                for a in ar2:
                    mask &= (ar1 != a)
            else:
                mask = np.zeros(len(ar1), dtype=np.bool)
                for a in ar2:
                    mask |= (ar1 == a)
            return mask

        # Otherwise use sorting
        if not assume_unique:
            ar1, rev_idx = np.unique(ar1, return_inverse=True)
            ar2 = np.unique(ar2)

        ar = np.concatenate((ar1, ar2))
        # We need this to be a stable sort, so always use 'mergesort'
        # here. The values from the first array should always come before
        # the values from the second array.
        order = ar.argsort(kind='mergesort')
        sar = ar[order]
        if invert:
            bool_ar = (sar[1:] != sar[:-1])
        else:
            bool_ar = (sar[1:] == sar[:-1])
        flag = np.concatenate((bool_ar, [invert]))
        indx = order.argsort(kind='mergesort')[:len(ar1)]

        if assume_unique:
            return flag[indx]
        else:
            return flag[indx][rev_idx]
else:
    from numpy import in1d
