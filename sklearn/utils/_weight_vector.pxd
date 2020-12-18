"""Efficient (dense) parameter vector implementation for linear models. """

cdef class WeightVector(object):
    cdef double *w_data_ptr
    cdef double *aw_data_ptr
    cdef double wscale
    cdef double average_a
    cdef double average_b
    cdef int n_features
    cdef double sq_norm

    cdef void add(self,  double *x_data_ptr, int *x_ind_ptr,
                  int xnnz, double c) nogil
    cdef void add_average(self,  double *x_data_ptr, int *x_ind_ptr,
                          int xnnz, double c, double num_iter) nogil
    cdef double dot(self, double *x_data_ptr, int *x_ind_ptr,
                    int xnnz) nogil
    cdef void scale(self, double c) nogil
    cdef void reset_wscale(self) nogil
    cdef double norm(self) nogil
