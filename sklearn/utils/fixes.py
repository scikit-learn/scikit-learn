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

import warnings
import sys
import functools
import os
import errno

import numpy as np
import scipy.sparse as sp
import scipy

try:
    from inspect import signature
except ImportError:
    from ..externals.funcsigs import signature


def _parse_version(version_string):
    version = []
    for x in version_string.split('.'):
        try:
            version.append(int(x))
        except ValueError:
            # x may be of the form dev-1ea1592
            version.append(x)
    return tuple(version)

euler_gamma = getattr(np, 'euler_gamma',
                      0.577215664901532860606512090082402431)

np_version = _parse_version(np.__version__)
sp_version = _parse_version(scipy.__version__)


# little danse to see if np.copy has an 'order' keyword argument
# Supported since numpy 1.7.0
if 'order' in signature(np.copy).parameters:
    def safe_copy(X):
        # Copy, but keep the order
        return np.copy(X, order='K')
else:
    # Before an 'order' argument was introduced, numpy wouldn't muck with
    # the ordering
    safe_copy = np.copy

try:
    if (not np.allclose(np.divide(.4, 1, casting="unsafe"),
                        np.divide(.4, 1, casting="unsafe", dtype=np.float64))
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
    # Compat where astype accepted no copy argument (numpy < 1.7.0)
    def astype(array, dtype, copy=True):
        if not copy and array.dtype == dtype:
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
    from numpy import partition
except ImportError:
    warnings.warn('Using `sort` instead of partition.'
                  'Upgrade numpy to 1.8 for better performace on large number'
                  'of clusters')
    def partition(a, kth, axis=-1, kind='introselect', order=None):
        return np.sort(a, axis=axis, order=order)


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


if sp_version < (0, 15):
    # Backport fix for scikit-learn/scikit-learn#2986 / scipy/scipy#4142
    from ._scipy_sparse_lsqr_backport import lsqr as sparse_lsqr
else:
    from scipy.sparse.linalg import lsqr as sparse_lsqr


def parallel_helper(obj, methodname, *args, **kwargs):
    """Helper to workaround Python 2 limitations of pickling instance methods"""
    return getattr(obj, methodname)(*args, **kwargs)


if np_version < (1, 6, 2):
    # Allow bincount to accept empty arrays
    # https://github.com/numpy/numpy/commit/40f0844846a9d7665616b142407a3d74cb65a040
    def bincount(x, weights=None, minlength=None):
        if len(x) > 0:
            return np.bincount(x, weights, minlength)
        else:
            if minlength is None:
                minlength = 0
            minlength = np.asscalar(np.asarray(minlength, dtype=np.intp))
            return np.zeros(minlength, dtype=np.intp)

else:
    from numpy import bincount


if 'exist_ok' in signature(os.makedirs).parameters:
    makedirs = os.makedirs
else:
    def makedirs(name, mode=0o777, exist_ok=False):
        """makedirs(name [, mode=0o777][, exist_ok=False])

        Super-mkdir; create a leaf directory and all intermediate ones.  Works
        like mkdir, except that any intermediate path segment (not just the
        rightmost) will be created if it does not exist. If the target
        directory already exists, raise an OSError if exist_ok is False.
        Otherwise no exception is raised.  This is recursive.

        """

        try:
            os.makedirs(name, mode=mode)
        except OSError as e:
            if (not exist_ok or e.errno != errno.EEXIST
                    or not os.path.isdir(name)):
                raise


if np_version < (1, 8, 1):
    def array_equal(a1, a2):
        # copy-paste from numpy 1.8.1
        try:
            a1, a2 = np.asarray(a1), np.asarray(a2)
        except:
            return False
        if a1.shape != a2.shape:
            return False
        return bool(np.asarray(a1 == a2).all())
else:
    from numpy import array_equal


if np_version < (1, 12):
    class MaskedArray(np.ma.MaskedArray):
        # Before numpy 1.12, np.ma.MaskedArray object is not picklable
        # This fix is needed to make our model_selection.GridSearchCV
        # picklable as the ``cv_results_`` param uses MaskedArray
        def __getstate__(self):
            """Return the internal state of the masked array, for pickling
            purposes.

            """
            cf = 'CF'[self.flags.fnc]
            data_state = super(np.ma.MaskedArray, self).__reduce__()[2]
            return data_state + (np.ma.getmaskarray(self).tostring(cf),
                                 self._fill_value)
else:
    from numpy.ma import MaskedArray    # noqa

if 'axis' not in signature(np.linalg.norm).parameters:

    def norm(X, ord=None, axis=None):
        """
        Handles the axis parameter for the norm function
        in old versions of numpy (useless for numpy >= 1.8).
        """

        if axis is None or X.ndim == 1:
            result = np.linalg.norm(X, ord=ord)
            return result

        if axis not in (0, 1):
            raise NotImplementedError("""
            The fix that adds axis parameter to the old numpy
            norm only works for 1D or 2D arrays.
            """)

        if axis == 0:
            X = X.T

        result = np.zeros(X.shape[0])
        for i in range(len(result)):
            result[i] = np.linalg.norm(X[i], ord=ord)

        return result

else:
    norm = np.linalg.norm
