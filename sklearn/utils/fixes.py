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


# Remove when minimum required NumPy >= 1.10
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


if sp_version < (0, 15):
    # Backport fix for scikit-learn/scikit-learn#2986 / scipy/scipy#4142
    from ._scipy_sparse_lsqr_backport import lsqr as sparse_lsqr
else:
    from scipy.sparse.linalg import lsqr as sparse_lsqr  # noqa


try:  # SciPy >= 0.19
    from scipy.special import comb, logsumexp
except ImportError:
    from scipy.misc import comb, logsumexp  # noqa


if sp_version >= (0, 19):
    def _argmax(arr_or_spmatrix, axis=None):
        return arr_or_spmatrix.argmax(axis=axis)
else:
    # Backport of argmax functionality from scipy 0.19.1, can be removed
    # once support for scipy 0.18 and below is dropped

    def _find_missing_index(ind, n):
        for k, a in enumerate(ind):
            if k != a:
                return k

        k += 1
        if k < n:
            return k
        else:
            return -1

    def _arg_min_or_max_axis(self, axis, op, compare):
        if self.shape[axis] == 0:
            raise ValueError("Can't apply the operation along a zero-sized "
                             "dimension.")

        if axis < 0:
            axis += 2

        zero = self.dtype.type(0)

        mat = self.tocsc() if axis == 0 else self.tocsr()
        mat.sum_duplicates()

        ret_size, line_size = mat._swap(mat.shape)
        ret = np.zeros(ret_size, dtype=int)

        nz_lines, = np.nonzero(np.diff(mat.indptr))
        for i in nz_lines:
            p, q = mat.indptr[i:i + 2]
            data = mat.data[p:q]
            indices = mat.indices[p:q]
            am = op(data)
            m = data[am]
            if compare(m, zero) or q - p == line_size:
                ret[i] = indices[am]
            else:
                zero_ind = _find_missing_index(indices, line_size)
                if m == zero:
                    ret[i] = min(am, zero_ind)
                else:
                    ret[i] = zero_ind

        if axis == 1:
            ret = ret.reshape(-1, 1)

        return np.asmatrix(ret)

    def _arg_min_or_max(self, axis, out, op, compare):
        if out is not None:
            raise ValueError("Sparse matrices do not support "
                             "an 'out' parameter.")

        # validateaxis(axis)

        if axis is None:
            if 0 in self.shape:
                raise ValueError("Can't apply the operation to "
                                 "an empty matrix.")

            if self.nnz == 0:
                return 0
            else:
                zero = self.dtype.type(0)
                mat = self.tocoo()
                mat.sum_duplicates()
                am = op(mat.data)
                m = mat.data[am]

                if compare(m, zero):
                    return mat.row[am] * mat.shape[1] + mat.col[am]
                else:
                    size = np.product(mat.shape)
                    if size == mat.nnz:
                        return am
                    else:
                        ind = mat.row * mat.shape[1] + mat.col
                        zero_ind = _find_missing_index(ind, size)
                        if m == zero:
                            return min(zero_ind, am)
                        else:
                            return zero_ind

        return _arg_min_or_max_axis(self, axis, op, compare)

    def _sparse_argmax(self, axis=None, out=None):
        return _arg_min_or_max(self, axis, out, np.argmax, np.greater)

    def _argmax(arr_or_matrix, axis=None):
        if sp.issparse(arr_or_matrix):
            return _sparse_argmax(arr_or_matrix, axis=axis)
        else:
            return arr_or_matrix.argmax(axis=axis)


def parallel_helper(obj, methodname, *args, **kwargs):
    """Workaround for Python 2 limitations of pickling instance methods"""
    return getattr(obj, methodname)(*args, **kwargs)


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
