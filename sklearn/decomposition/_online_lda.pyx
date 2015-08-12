cimport cython
cimport numpy as np
import numpy as np

np.import_array()

from libc.math cimport exp, fabs
from scipy.special import psi


@cython.boundscheck(False)
@cython.wraparound(False)
def mean_change(np.ndarray[ndim=1, dtype=np.float64_t] arr_1,
                np.ndarray[ndim=1, dtype=np.float64_t] arr_2):
    """Calculate the mean difference between two arrays.

    Equivalent to np.abs(arr_1 - arr2).mean().
    """

    cdef np.float64_t total, diff
    cdef np.npy_intp i, size

    size = arr_1.shape[0]
    total = 0.0
    for i in range(size):
        diff = fabs(arr_1[i] - arr_2[i])
        total += diff

    return total / size


@cython.boundscheck(False)
@cython.wraparound(False)
def _dirichlet_expectation_1d(np.ndarray[ndim=1, dtype=np.float64_t] arr):
    """Dirichlet expectation for a single sample:
    exp(E[log(theta)]) for theta ~ Dir(arr).

    Equivalent to np.exp(psi(arr) - psi(np.sum(arr))).
    """

    cdef np.float64_t total, psi_total
    cdef np.ndarray[ndim=1, dtype=np.float64_t] d_exp
    cdef np.npy_intp i, size

    size = arr.shape[0]
    d_exp = psi(arr)

    total = 0.0
    for i in range(size):
        total += arr[i]

    psi_total = psi(total)
    for i in range(size):
        d_exp[i] = exp(d_exp[i] - psi_total)

    return d_exp


@cython.boundscheck(False)
@cython.wraparound(False)
def _dirichlet_expectation_2d(np.ndarray[ndim=2, dtype=np.float64_t] arr):
    """Dirichlet expectation for multiple samples:
    E[log(theta)] for theta ~ Dir(arr).

    Equivalent to psi(arr) - psi(np.sum(arr, axis=1))[:, np.newaxis].

    Note that unlike _dirichlet_expectation_1d, this function doesn't compute
    the exp.
    """

    cdef np.ndarray[ndim=2, dtype=np.float64_t] d_exp
    cdef np.ndarray[ndim=1, dtype=np.float64_t] psi_row
    cdef np.npy_intp i, j, n_rows, n_cols

    n_rows = arr.shape[0]
    n_cols = arr.shape[1]

    psi_row = np.zeros(n_rows, dtype=np.float64)
    for i in range(n_rows):
        for j in range(n_cols):
            psi_row[i] += arr[i, j]
    psi_row = psi(psi_row, out=psi_row)

    d_exp = psi(arr)
    for i in range(n_rows):
        for j in range(n_cols):
            d_exp[i, j] -= psi_row[i]

    return d_exp
