cimport cython
cimport numpy as np
import numpy as np

np.import_array()

from libc.math cimport exp, fabs, log
from numpy.math cimport EULER


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
def _dirichlet_expectation_1d(np.ndarray[ndim=1, dtype=np.float64_t] doc_topic,
                              double doc_topic_prior,
                              np.ndarray[ndim=1, dtype=np.float64_t] out):
    """Dirichlet expectation for a single sample:
        exp(E[log(theta)]) for theta ~ Dir(doc_topic)
    after adding doc_topic_prior to doc_topic, in-place.

    Equivalent to
        doc_topic += doc_topic_prior
        out[:] = np.exp(psi(doc_topic) - psi(np.sum(doc_topic)))
    """

    cdef np.float64_t dt, psi_total, total
    cdef np.npy_intp i, size

    size = doc_topic.shape[0]

    total = 0.0
    for i in range(size):
        dt = doc_topic[i] + doc_topic_prior
        doc_topic[i] = dt
        total += dt
    psi_total = psi(total)

    for i in range(size):
        out[i] = exp(psi(doc_topic[i]) - psi_total)


@cython.boundscheck(False)
@cython.wraparound(False)
def _dirichlet_expectation_2d(np.ndarray[ndim=2, dtype=np.float64_t] arr):
    """Dirichlet expectation for multiple samples:
    E[log(theta)] for theta ~ Dir(arr).

    Equivalent to psi(arr) - psi(np.sum(arr, axis=1))[:, np.newaxis].

    Note that unlike _dirichlet_expectation_1d, this function doesn't compute
    the exp and doesn't add in the prior.
    """
    cdef np.float64_t row_total, psi_row_total
    cdef np.ndarray[ndim=2, dtype=np.float64_t] d_exp
    cdef np.npy_intp i, j, n_rows, n_cols

    n_rows = arr.shape[0]
    n_cols = arr.shape[1]

    d_exp = np.empty_like(arr)
    for i in range(n_rows):
        row_total = 0
        for j in range(n_cols):
            row_total += arr[i, j]
        psi_row_total = psi(row_total)

        for j in range(n_cols):
            d_exp[i, j] = psi(arr[i, j]) - psi_row_total

    return d_exp


# Psi function for positive arguments. Optimized for speed, not accuracy.
#
# After: J. Bernardo (1976). Algorithm AS 103: Psi (Digamma) Function.
# https://www.uv.es/~bernardo/1976AppStatist.pdf
@cython.cdivision(True)
cdef double psi(double x) nogil:
    if x <= 1e-6:
        # psi(x) = -EULER - 1/x + O(x)
        return -EULER - 1. / x

    cdef double r, result = 0

    # psi(x + 1) = psi(x) + 1/x
    while x < 6:
        result -= 1. / x
        x += 1

    # psi(x) = log(x) - 1/(2x) - 1/(12x**2) + 1/(120x**4) - 1/(252x**6)
    #          + O(1/x**8)
    r = 1. / x
    result += log(x) - .5 * r
    r = r * r
    result -= r * ((1./12.) - r * ((1./120.) - r * (1./252.)))
    return result;
