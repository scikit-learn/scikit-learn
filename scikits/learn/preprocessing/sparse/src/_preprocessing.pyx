# Author: Mathieu Blondel
#
# License: BSD Style.

cimport numpy as np
import numpy as np
import numpy.linalg as linalg
cimport cython

cdef extern from "math.h":
    double fabs(double f)
    double sqrt(double f)

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def normalize_axis1_sparse(X):
    """Inplace row normalize using the l1 norm"""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[INTEGER, ndim=1] X_indices = X.indices
    cdef np.ndarray[INTEGER, ndim=1] X_indptr = X.indptr

    # the column indices for row i are stored in indices[indptr[i]:indices[i+1]]
    # and their corresponding values are stored in data[indptr[i]:indptr[i+1]]
    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += fabs(X_data[j])

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def normalize_length_axis1_sparse(X):
    """Inplace row normalize using the l2 norm"""
    cdef unsigned int n_samples = X.shape[0]
    cdef unsigned int n_features = X.shape[1]

    cdef np.ndarray[DOUBLE, ndim=1] X_data = X.data
    cdef np.ndarray[INTEGER, ndim=1] X_indices = X.indices
    cdef np.ndarray[INTEGER, ndim=1] X_indptr = X.indptr

    cdef unsigned int i
    cdef unsigned int j
    cdef double sum_

    for i in xrange(n_samples):
        sum_ = 0.0

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            sum_ += (X_data[j] * X_data[j])

        sum_ = sqrt(sum_)

        for j in xrange(X_indptr[i], X_indptr[i + 1]):
            X_data[j] /= sum_

