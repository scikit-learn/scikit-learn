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

if sp_version >= (1, 4):
    from scipy.sparse.linalg import lobpcg
else:
    # Backport of lobpcg functionality from scipy 1.4.0, can be removed
    # once support for sp_version < (1, 4) is dropped
    from ..externals._lobpcg import lobpcg  # noqa

if sp_version >= (1, 3):
    # Preserves earlier default choice of pinvh cutoff `cond` value.
    # Can be removed once issue #14055 is fully addressed.
    from ..externals._scipy_linalg import pinvh
else:
    from scipy.linalg import pinvh # noqa

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


# TODO: replace by copy=False, when only scipy > 1.1 is supported.
def _astype_copy_false(X):
    """Returns the copy=False parameter for
    {ndarray, csr_matrix, csc_matrix}.astype when possible,
    otherwise don't specify
    """
    if sp_version >= (1, 1) or not sp.issparse(X):
        return {'copy': False}
    else:
        return {}


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
    import joblib

    if joblib.__version__ >= LooseVersion('0.12'):
        return kwargs

    extra_args = set(kwargs.keys()).difference({'prefer', 'require'})
    if extra_args:
        raise NotImplementedError('unhandled arguments %s with joblib %s'
                                  % (list(extra_args), joblib.__version__))
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
