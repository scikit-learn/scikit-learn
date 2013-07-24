#cython: boundscheck=False
#cython: cdivision=True

from libc.math cimport log, exp

import numpy as np
cimport numpy as np

ctypedef np.float64_t DTYPE_t


cdef DTYPE_t _inner_logistic_sigmoid(DTYPE_t x):
    """Logistic sigmoid function 1 / (1 + e ** -x)"""

    if x > 0:
        return 1. / (1 + exp(-x))
    else:
        exp_x = exp(x)
        return exp_x / (1. + exp_x)


cdef DTYPE_t _inner_log_logistic_sigmoid(DTYPE_t x):
    """Log of the logistic sigmoid function log(1 / (1 + e ** -x))"""
    if x > 0:
        return -log(1 + exp(-x))
    else:
        return x - log(1 + exp(x))


def _logistic_sigmoid(int n_samples, int n_features, 
                     np.ndarray[DTYPE_t, ndim=2] X,
                     np.ndarray[DTYPE_t, ndim=2] out):
    for i in range(n_samples):
        for j in range(n_features):
            out[i, j] = _inner_logistic_sigmoid(X[i, j])
    return out


def _log_logistic_sigmoid(int n_samples, int n_features, 
                           np.ndarray[DTYPE_t, ndim=2] X,
                           np.ndarray[DTYPE_t, ndim=2] out):
    for i in range(n_samples):
        for j in range(n_features):
            out[i, j] = _inner_log_logistic_sigmoid(X[i, j])
    return out
