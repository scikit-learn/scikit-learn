"""Efficient (dense) parameter vector implementation for linear models. """

cimport numpy as np


cdef extern from "math.h":
    cdef extern double sqrt(double x)


ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


cdef class WeightVector(object):
    cdef np.ndarray w
    cdef double *w_data_ptr
    cdef double wscale
    cdef np.ndarray intercept
    cdef double *intercept_data_ptr
    cdef Py_ssize_t n_features
    cdef Py_ssize_t K
    cdef double sq_norm

    cdef void add(self,  DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, int k, double c)
    cdef double *dot(self, DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                    int xnnz)
    cdef void scale(self, double c)
    cdef void reset_wscale(self)
    cdef double norm(self)