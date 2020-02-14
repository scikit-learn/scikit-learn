# cython: cdivision=True, boundscheck=False, wraparound=False
#
# Author: Mathieu Blondel, Tom Dupre la Tour
# License: BSD 3 clause

from cython cimport floating
from libc.math cimport fabs
from cython.parallel import prange

from ..utils._cython_blas cimport _dot
from ..utils._openmp_helpers import _openmp_effective_n_threads


def _update_cdnmf_fast(floating[:, ::1] W, floating[:, ::1] HHt,
                       floating[:, ::1] XHt, Py_ssize_t[::1] permutation):
    cdef:
        floating violation = 0
        Py_ssize_t n_components = W.shape[1]
        Py_ssize_t n_samples = W.shape[0]  # n_features for H update
        Py_ssize_t i, s, t
        int num_threads = _openmp_effective_n_threads()

    with nogil:
        for s in range(n_components):
            t = permutation[s]

            for i in prange(n_samples, num_threads=num_threads):
                violation += _update_cdnmf_sample(
                    n_components, &HHt[t, 0], &W[i, 0], XHt[i, t], t)

    return violation


cdef floating _update_cdnmf_sample(Py_ssize_t n_components,
                                   floating* HHt,
                                   floating* W,
                                   floating xht,
                                   Py_ssize_t t) nogil:
    cdef:
        floating hess = HHt[t]
        floating grad, pg, violation

    # np.dot(W[i, :], HHt[t, :]) - XHt[i, t]
    grad = _dot(n_components, HHt, 1, W, 1) - xht

    pg = min(grad, 0) if W[t] == 0 else grad
    violation = fabs(pg)

    if hess != 0:
        W[t] = max(W[t] - grad / hess, 0)

    return violation
