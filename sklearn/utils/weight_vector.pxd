"""Efficient (dense) parameter vector implementation for linear models. """

cimport numpy as np


cdef extern from "math.h":
    cdef extern double sqrt(double x)


cdef class WeightVector(object):
    cdef np.ndarray w
    cdef np.ndarray aw
    cdef double *w_data_ptr
    cdef double *aw_data_ptr
    cdef double wscale
    cdef double *wscale_ptr
    cdef double average_a
    cdef double average_b
    cdef int n_features
    cdef int n_etas
    cdef double sq_norm

    cdef void add(self,  double *x_data_ptr, int *x_ind_ptr,
                  int xnnz, double *c, double update) nogil
    cdef void add_average(self,  double *x_data_ptr, int *x_ind_ptr,
                          int xnnz, double *c, double update, double num_iter) nogil
    cdef double dot(self, double *x_data_ptr, int *x_ind_ptr,
                    int xnnz) nogil
    cdef void scale(self, double c) nogil
    cdef void scale_vector(self, double l1_ratio, double* eta_ptr, double alpha) nogil
    cdef void reset_wscale(self) nogil
    cdef double norm(self) nogil
