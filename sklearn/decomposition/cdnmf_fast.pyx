# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Mathieu Blondel, Tom Dupre la Tour
# License: BSD 3 clause

import numpy as np
cimport cython
cimport numpy as np
from libc.math cimport fabs


cdef inline double fmax(double x, double y) nogil:
    if x > y:
        return x
    return y


cdef inline double fmin(double x, double y) nogil:
    if x < y:
        return x
    return y


def _update_cdnmf_fast(double[:, ::1] W, double[:, :] HHt, double[:, :] XHt,
                       bint shuffle, int seed):

    cdef double violation = 0
    cdef int n_components = W.shape[1]
    cdef int n_samples = W.shape[0]  # n_features for H update
    cdef double grad, pg, hess
    cdef int j, i, t, r
    cdef np.ndarray[long, ndim=1] permutation_array
    cdef long* permutation = NULL

    if shuffle:
        rng = np.random.RandomState(seed)
        permutation_array = rng.permutation(n_samples)
    else:
        permutation_array = np.arange(n_samples)
    permutation = <long*> permutation_array.data

    with nogil:
        for j in range(n_samples):
            i = permutation[j]
            for t in range(n_components):

                # gradient = GW[t, i] where GW = np.dot(W, HHt) - XHt
                grad = -XHt[i, t]

                for r in range(n_components):
                    grad += HHt[t, r] * W[i, r]

                # projected gradient
                pg = fmin(0., grad) if W[i, t] == 0 else grad
                violation += fabs(pg)

                # Hessian
                hess = HHt[t, t]

                if hess != 0:
                    W[i, t] = fmax(W[i, t] - grad / hess, 0.)
                
    return violation
