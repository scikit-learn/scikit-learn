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

from distutils.version import LooseVersion

import numpy as np
import scipy.sparse as sp
import scipy
from scipy.sparse.linalg import lsqr as sparse_lsqr  # noqa


def _parse_version(version_string):
    version = []
    for x in version_string.split('.'):
        try:
            version.append(int(x))
        except ValueError:
            # x may be of the form dev-1ea1592
            version.append(x)
    return tuple(version)


np_version = _parse_version(np.__version__)
sp_version = _parse_version(scipy.__version__)


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
    """Workaround for Python 2 limitations of pickling instance methods

    Parameters
    ----------
    obj
    methodname
    *args
    **kwargs

    """
    return getattr(obj, methodname)(*args, **kwargs)


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


if np_version < (1, 9, 0):
    # Allow unique to return counts
    # https://github.com/numpy/numpy/commit/09fb4205a1d56090e13257a181f23514684f532b

    def unique(ar, return_index=False, return_inverse=False,
               return_counts=False):
        ar = np.asanyarray(ar).flatten()

        optional_indices = return_index or return_inverse
        optional_returns = optional_indices or return_counts

        if ar.size == 0:
            if not optional_returns:
                ret = ar
            else:
                ret = (ar,)
                if return_index:
                    ret += (np.empty(0, np.bool),)
                if return_inverse:
                    ret += (np.empty(0, np.bool),)
                if return_counts:
                    ret += (np.empty(0, np.intp),)
            return ret

        if optional_indices:
            perm = ar.argsort(
                kind='mergesort' if return_index else 'quicksort')
            aux = ar[perm]
        else:
            ar.sort()
            aux = ar
        flag = np.concatenate(([True], aux[1:] != aux[:-1]))

        if not optional_returns:
            ret = aux[flag]
        else:
            ret = (aux[flag],)
            if return_index:
                ret += (perm[flag],)
            if return_inverse:
                iflag = np.cumsum(flag) - 1
                iperm = perm.argsort()
                ret += (np.take(iflag, iperm),)
            if return_counts:
                idx = np.concatenate(np.nonzero(flag) + ([ar.size],))
                ret += (np.diff(idx),)
        return ret
else:
    from numpy import unique    # noqa


if np_version < (1, 11):
    def nanpercentile(a, q):
        """
        Compute the qth percentile of the data along the specified axis,
        while ignoring nan values.

        Returns the qth percentile(s) of the array elements.

        Parameters
        ----------
        a : array_like
            Input array or object that can be converted to an array.
        q : float in range of [0,100] (or sequence of floats)
            Percentile to compute, which must be between 0 and 100
            inclusive.

        Returns
        -------
        percentile : scalar or ndarray
            If `q` is a single percentile and `axis=None`, then the result
            is a scalar. If multiple percentiles are given, first axis of
            the result corresponds to the percentiles. The other axes are
            the axes that remain after the reduction of `a`. If the input
            contains integers or floats smaller than ``float64``, the output
            data-type is ``float64``. Otherwise, the output data-type is the
            same as that of the input. If `out` is specified, that array is
            returned instead.

        """
        data = np.compress(~np.isnan(a), a)
        if data.size:
            return np.percentile(data, q)
        else:
            size_q = 1 if np.isscalar(q) else len(q)
            return np.array([np.nan] * size_q)
else:
    from numpy import nanpercentile # noqa


# Fix for behavior inconsistency on numpy.equal for object dtypes.
# For numpy versions < 1.13, numpy.equal tests element-wise identity of objects
# instead of equality. This fix returns the mask of NaNs in an array of
# numerical or object values for all numpy versions.
if np_version < (1, 13):
    def _object_dtype_isnan(X):
        return np.frompyfunc(lambda x: x != x, 1, 1)(X).astype(bool)
else:
    def _object_dtype_isnan(X):
        return X != X


def _joblib_parallel_args(**kwargs):
    """Set joblib.Parallel arguments in a compatible way for 0.11 and 0.12+

    For joblib 0.11 this maps both ``prefer`` and ``require`` parameters to
    a specific ``backend``.

    Parameters
    ----------

    prefer : str in {'processes', 'threads'} or None
        Soft hint to choose the default backend if no specific backend
        was selected with the parallel_backend context manager.

    require : 'sharedmem' or None
        Hard condstraint to select the backend. If set to 'sharedmem',
        the selected backend will be single-host and thread-based even
        if the user asked for a non-thread based backend with
        parallel_backend.

    See joblib.Parallel documentation for more details
    """
    from . import _joblib

    if _joblib.__version__ >= LooseVersion('0.12'):
        return kwargs

    extra_args = set(kwargs.keys()).difference({'prefer', 'require'})
    if extra_args:
        raise NotImplementedError('unhandled arguments %s with joblib %s'
                                  % (list(extra_args), _joblib.__version__))
    args = {}
    if 'prefer' in kwargs:
        prefer = kwargs['prefer']
        if prefer not in ['threads', 'processes', None]:
            raise ValueError('prefer=%s is not supported' % prefer)
        args['backend'] = {'threads': 'threading',
                           'processes': 'multiprocessing',
                           None: None}[prefer]

    if 'require' in kwargs:
        require = kwargs['require']
        if require not in [None, 'sharedmem']:
            raise ValueError('require=%s is not supported' % require)
        if require == 'sharedmem':
            args['backend'] = 'threading'
    return args
