"""Efficient (dense) parameter vector implementation for linear models. """

cimport numpy as np


cdef extern from "math.h":
    cdef extern double sqrt(double x)


ctypedef np.float64_t DTYPE
ctypedef np.int32_t INTEGER


cdef class WeightVector(object):
    cdef public np.ndarray w
    cdef np.float64_t *w_data_ptr
    cdef double w_scale
    cdef public np.ndarray intercept
    cdef np.float64_t *intercept_data_ptr
    cdef Py_ssize_t n_features
    cdef Py_ssize_t K
    cdef double sq_norm
    cdef int fit_intercept
    cdef double intercept_decay

    cdef void add(self,  DTYPE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, int k, double c)
    cdef void dot(self, DTYPE *x_data_ptr, INTEGER *x_ind_ptr,
                    int xnnz, double *out)
    cdef void scale(self, double c)
    cdef void reset_scale(self)
    cdef double norm(self)