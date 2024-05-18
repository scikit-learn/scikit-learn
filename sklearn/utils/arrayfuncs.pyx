"""
The :mod:`sklearn.utils.arrayfuncs` module includes a small collection of auxiliary
functions that operate on arrays.
"""

from cython cimport floating
from cython.parallel cimport prange
from libc.math cimport fabs
from libc.float cimport DBL_MAX, FLT_MAX

from ._cython_blas cimport _copy, _rotg, _rot
from ._typedefs cimport float64_t


ctypedef fused real_numeric:
    short
    int
    long
    long long
    float
    double


def min_pos(const floating[:] X):
    """Find the minimum value of an array over positive values.

    Returns the maximum representable value of the input dtype if none of the
    values are positive.

    Parameters
    ----------
    X : ndarray of shape (n,)
        Input array.

    Returns
    -------
    min_val : float
        The smallest positive value in the array, or the maximum representable value
         of the input dtype if no positive values are found.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.arrayfuncs import min_pos
    >>> X = np.array([0, -1, 2, 3, -4, 5])
    >>> min_pos(X)
    2.0
    """
    cdef Py_ssize_t i
    cdef floating min_val = FLT_MAX if floating is float else DBL_MAX
    for i in range(X.size):
        if 0. < X[i] < min_val:
            min_val = X[i]
    return min_val


def _all_with_any_reduction_axis_1(real_numeric[:, :] array, real_numeric value):
    """Check whether any row contains all values equal to `value`.

    It is equivalent to `np.any(np.all(X == value, axis=1))`, but it avoids to
    materialize the temporary boolean matrices in memory.

    Parameters
    ----------
    array: array-like
        The array to be checked.
    value: short, int, long, float, or double
        The value to use for the comparison.

    Returns
    -------
    any_all_equal: bool
        Whether or not any rows contains all values equal to `value`.
    """
    cdef Py_ssize_t i, j

    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            if array[i, j] != value:
                break
        else:  # no break
            return True
    return False


# General Cholesky Delete.
# Remove an element from the cholesky factorization
# m = columns
# n = rows
#
# TODO: put transpose as an option
def cholesky_delete(floating[:, :] L, int go_out):
    cdef:
        int n = L.shape[0]
        int m = L.strides[0]
        floating c, s
        floating *L1
        int i

    if floating is float:
        m /= sizeof(float)
    else:
        m /= sizeof(double)

    # delete row go_out
    L1 = &L[0, 0] + (go_out * m)
    for i in range(go_out, n-1):
        _copy(i + 2, L1 + m, 1, L1, 1)
        L1 += m

    L1 = &L[0, 0] + (go_out * m)
    for i in range(go_out, n-1):
        _rotg(L1 + i, L1 + i + 1, &c, &s)
        if L1[i] < 0:
            # Diagonals cannot be negative
            L1[i] = fabs(L1[i])
            c = -c
            s = -s

        L1[i + 1] = 0.  # just for cleanup
        L1 += m

        _rot(n - i - 2, L1 + i, m, L1 + i + 1, m, c, s)


def sum_parallel(const floating [:] array, int n_threads):
    """Parallel sum, always using float64 internally."""
    cdef:
        float64_t out = 0.
        int i = 0

    for i in prange(
        array.shape[0], schedule='static', nogil=True, num_threads=n_threads
    ):
        out += array[i]

    return out
