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
        floating pg, grad, hess
        Py_ssize_t n_components = W.shape[1]
        Py_ssize_t n_samples = W.shape[0]  # n_features for H update
        Py_ssize_t i, s, t
        int num_threads = _openmp_effective_n_threads()

    with nogil:
        for s in range(n_components):
            t = permutation[s]

            for i in prange(n_samples, num_threads=num_threads):
                # np.dot(W[i, :], HHt[t, :]) - XHt[i, t]
                grad = _dot(n_components, &HHt[t, 0], 1, &W[i, 0], 1) - XHt[i, t]
                # grad = grad - XHt[i, t]  # "-=" is interpreted as reduction

                pg = min(grad, 0) if W[i, t] == 0 else grad
                violation += fabs(pg)

                hess = HHt[t, t]
                if hess != 0:
                    W[i, t] = max(W[i, t] - grad / hess, 0)

    return violation
