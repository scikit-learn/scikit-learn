# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Author: Mathieu Blondel, Tom Dupre la Tour
# License: BSD 3 clause

cimport cython
from libc.math cimport fabs
from cython.parallel import prange


cdef inline double _update_cdnmf_samples(unsigned n_components,double xht, double [] HHt,double[] W,double hess,unsigned t) nogil:
    cdef double violation = 0
    cdef double grad = 0 - xht
    cdef double pg
    for r in range(n_components):
    # for(int r =0;r<n_components;r++)
        grad += HHt[r] * W[r]
    pg = grad
    if (W[t] == 0) and (pg > 0):
        pg = 0
    violation = fabs(pg)
    if hess != 0:
        W[t] -= grad/hess
        if W[t] < 0 :
            W[t] = 0
    return violation


def _update_cdnmf_fast(double[:, ::1] W, double[:, :] HHt, double[:, :] XHt,
                       Py_ssize_t[::1] permutation):
    cdef double violation = 0
    cdef Py_ssize_t n_components = W.shape[1]
    cdef Py_ssize_t n_samples = W.shape[0]  # n_features for H update
    cdef double pg, hess
    cdef Py_ssize_t i, r, s, t

    # print "in nogil!!"
    with nogil:
        for s in range(n_components):
            t = permutation[s]
            # Hessian
            hess = HHt[t, t]

            for i in prange(n_samples):
                violation += _update_cdnmf_samples(n_components,XHt[i, t],&HHt[t,0],&W[i,0],hess,t)
                # gradient = GW[t, i] where GW = np.dot(W, HHt) - XHt
                # grad = -XHt[i, t]

                # for r in range(n_components):
                #     grad += HHt[t, r] * W[i, r]

                # projected gradient
                # pg = min(0., grad) if W[i, t] == 0 else grad
                # violation += fabs(pg)

                # Hessian
                # hess = HHt[t, t]

                # if hess != 0:
                #     W[i, t] = max(W[i, t] - grad / hess, 0.)
                
    # print "out nogil!!"
    return violation
