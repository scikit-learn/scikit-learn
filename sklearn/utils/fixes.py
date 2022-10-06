"""Compatibility fixes for older version of python, numpy and scipy

If you add content to this file, please give the version of the package
at which the fix is no longer needed.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Fabian Pedregosa <fpedregosa@acm.org>
#          Lars Buitinck
#
# License: BSD 3 clause

from functools import update_wrapper
import functools

import sklearn
import numpy as np
import scipy
import scipy.stats
import threadpoolctl
from .._config import config_context, get_config
from ..externals._packaging.version import parse as parse_version


np_version = parse_version(np.__version__)
sp_version = parse_version(scipy.__version__)


if sp_version >= parse_version("1.4"):
    from scipy.sparse.linalg import lobpcg
else:
    # Backport of lobpcg functionality from scipy 1.4.0, can be removed
    # once support for sp_version < parse_version('1.4') is dropped
    # mypy error: Name 'lobpcg' already defined (possibly by an import)
    from ..externals._lobpcg import lobpcg  # type: ignore  # noqa

try:
    from scipy.optimize._linesearch import line_search_wolfe2, line_search_wolfe1
except ImportError:  # SciPy < 1.8
    from scipy.optimize.linesearch import line_search_wolfe2, line_search_wolfe1  # type: ignore  # noqa


def _object_dtype_isnan(X):
    return X != X


class loguniform(scipy.stats.reciprocal):
    """A class supporting log-uniform random variables.

    Parameters
    ----------
    low : float
        The minimum value
    high : float
        The maximum value

    Methods
    -------
    rvs(self, size=None, random_state=None)
        Generate log-uniform random variables

    The most useful method for Scikit-learn usage is highlighted here.
    For a full list, see
    `scipy.stats.reciprocal
    <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.reciprocal.html>`_.
    This list includes all functions of ``scipy.stats`` continuous
    distributions such as ``pdf``.

    Notes
    -----
    This class generates values between ``low`` and ``high`` or

        low <= loguniform(low, high).rvs() <= high

    The logarithmic probability density function (PDF) is uniform. When
    ``x`` is a uniformly distributed random variable between 0 and 1, ``10**x``
    are random variables that are equally likely to be returned.

    This class is an alias to ``scipy.stats.reciprocal``, which uses the
    reciprocal distribution:
    https://en.wikipedia.org/wiki/Reciprocal_distribution

    Examples
    --------

    >>> from sklearn.utils.fixes import loguniform
    >>> rv = loguniform(1e-3, 1e1)
    >>> rvs = rv.rvs(random_state=42, size=1000)
    >>> rvs.min()  # doctest: +SKIP
    0.0010435856341129003
    >>> rvs.max()  # doctest: +SKIP
    9.97403052786026
    """


# remove when https://github.com/joblib/joblib/issues/1071 is fixed
def delayed(function):
    """Decorator used to capture the arguments of a function."""

    @functools.wraps(function)
    def delayed_function(*args, **kwargs):
        return _FuncWrapper(function), args, kwargs

    return delayed_function


class _FuncWrapper:
    """ "Load the global configuration before calling the function."""

    def __init__(self, function):
        self.function = function
        self.config = get_config()
        update_wrapper(self, self.function)

    def __call__(self, *args, **kwargs):
        with config_context(**self.config):
            return self.function(*args, **kwargs)


# Rename the `method` kwarg to `interpolation` for NumPy < 1.22, because
# `interpolation` kwarg was deprecated in favor of `method` in NumPy >= 1.22.
def _percentile(a, q, *, method="linear", **kwargs):
    return np.percentile(a, q, interpolation=method, **kwargs)


if np_version < parse_version("1.22"):
    percentile = _percentile
else:  # >= 1.22
    from numpy import percentile  # type: ignore  # noqa


# compatibility fix for threadpoolctl >= 3.0.0
# since version 3 it's possible to setup a global threadpool controller to avoid
# looping through all loaded shared libraries each time.
# the global controller is created during the first call to threadpoolctl.
def _get_threadpool_controller():
    if not hasattr(threadpoolctl, "ThreadpoolController"):
        return None

    if not hasattr(sklearn, "_sklearn_threadpool_controller"):
        sklearn._sklearn_threadpool_controller = threadpoolctl.ThreadpoolController()

    return sklearn._sklearn_threadpool_controller


def threadpool_limits(limits=None, user_api=None):
    controller = _get_threadpool_controller()
    if controller is not None:
        return controller.limit(limits=limits, user_api=user_api)
    else:
        return threadpoolctl.threadpool_limits(limits=limits, user_api=user_api)


threadpool_limits.__doc__ = threadpoolctl.threadpool_limits.__doc__


def threadpool_info():
    controller = _get_threadpool_controller()
    if controller is not None:
        return controller.info()
    else:
        return threadpoolctl.threadpool_info()


threadpool_info.__doc__ = threadpoolctl.threadpool_info.__doc__


# TODO: Remove when SciPy 1.9 is the minimum supported version
def _mode(a, axis=0):
    if sp_version >= parse_version("1.9.0"):
        return scipy.stats.mode(a, axis=axis, keepdims=True)
    return scipy.stats.mode(a, axis=axis)


# TODO: Remove when SciPy 1.10 is the minimum supported version
if sp_version >= parse_version("1.9.2"):
    csr_hstack = scipy.sparse.hstack
else:

    def csr_hstack(columns, dtype=np.float64):
        """
        Parameters
        ----------
        columns : list
            List of `CSR` matrices to horizontally (column-wise) stack. All
            matrices must have the same number of rows.
        dtype : dtype, default=np.float64
            The type of feature values. The output of the stacking operation is
            cast to this type if this type is wider (i.e. `float32`-->`float64`).
        Returns
        -------
        X : CSR matrix of shape (`n_rows`, `n_features_out_`)
            A CSR sparse matrix that is the result of horizontally (column-wise)
            stacking the matrices contained in `columns`.
        """
        from ._csr_hstack import _csr_hstack as cython_csr_hstack

        n_blocks = len(columns)
        if n_blocks == 0:
            raise ValueError("No matrices were provided to stack")
        if n_blocks == 1:
            return columns[0]
        other_axis_dims = set(mat.shape[0] for mat in columns)
        if len(other_axis_dims) > 1:
            raise ValueError(
                f"Mismatching dimensions along axis {0}: {other_axis_dims}"
            )
        (constant_dim,) = other_axis_dims

        # Do the stacking
        indptr_list = [mat.indptr for mat in columns]
        data_cat = np.concatenate([mat.data for mat in columns])

        # Need to check if any indices/indptr, would be too large post-
        # concatenation for np.int32. We must sum the dimensions along the axis we
        # use to stack since even empty matrices will contribute to a large post-
        # stack index value. At this point the matrices contained in `columns` may
        # be a mix of {32,64}bit integers, but this covers the case that each are
        # only 32bit while their concatenation would need to be 64bit.
        max_output_index = 0
        max_indptr = 0
        for mat in columns[:-1]:
            max_output_index += mat.shape[1]
            max_indptr = max(max_indptr, mat.indptr.max())
        if columns[-1].indices.size > 0:
            max_output_index += int(columns[-1].indices.max())
            max_indptr = max(max_indptr, columns[-1].indptr.max())
        max_int32 = np.iinfo(np.int32).max
        needs_64bit = max(max_output_index, max_indptr) > max_int32
        idx_dtype = np.int64 if needs_64bit else np.int32

        stack_dim_cat = np.array([mat.shape[1] for mat in columns], dtype=np.int64)
        if data_cat.size > 0:
            indptr_cat = np.concatenate(indptr_list).astype(idx_dtype)
            indices_cat = np.concatenate([mat.indices for mat in columns]).astype(
                idx_dtype
            )
            indptr = np.empty(constant_dim + 1, dtype=idx_dtype)
            indices = np.empty_like(indices_cat)
            data = np.empty_like(data_cat)
            cython_csr_hstack(
                n_blocks,
                constant_dim,
                stack_dim_cat,
                indptr_cat,
                indices_cat,
                data_cat,
                indptr,
                indices,
                data,
            )
        else:
            indptr = np.zeros(constant_dim + 1, dtype=idx_dtype)
            indices = np.empty(0, dtype=idx_dtype)
            data = np.empty(0, dtype=data_cat.dtype)
        sum_dim = stack_dim_cat.sum()
        return scipy.sparse.csr_matrix(
            (data.astype(dtype), indices, indptr), shape=(constant_dim, sum_dim)
        )
