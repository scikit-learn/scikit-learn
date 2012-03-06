"""Utils and computational tricks for large-scale algorithms. """

cimport numpy as np


cdef extern from "math.h":
    cdef extern double sqrt(double x)


ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


################################################################################
# Efficient (dense) parameter vector implementation for linear models.
#

cdef class WeightVector:
    cdef np.ndarray w
    cdef DOUBLE *w_data_ptr
    cdef double wscale
    cdef Py_ssize_t n_features
    cdef double sq_norm

    cdef void add(self,  DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                  int xnnz, double c)
    cdef double dot(self, DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                    int xnnz)
    cdef void scale(self, double c)
    cdef void reset_wscale(self)
    cdef double norm(self)



################################################################################
# Dataset abstractions for sequential data access
#

cdef class Dataset:
    cdef Py_ssize_t n_samples

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight)
    cdef void shuffle(self, seed)


cdef class ArrayDataset(Dataset):
    cdef Py_ssize_t n_features
    cdef int current_index
    cdef int stride
    cdef DOUBLE *X_data_ptr
    cdef DOUBLE *Y_data_ptr
    cdef np.ndarray feature_indices
    cdef INTEGER *feature_indices_ptr
    cdef np.ndarray index
    cdef INTEGER *index_data_ptr
    cdef DOUBLE *sample_weight_data

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight)
    cdef void shuffle(self, seed)


cdef class CSRDataset(Dataset):
    cdef int current_index
    cdef int stride
    cdef DOUBLE *X_data_ptr
    cdef INTEGER *X_indptr_ptr
    cdef INTEGER *X_indices_ptr
    cdef DOUBLE *Y_data_ptr
    cdef np.ndarray feature_indices
    cdef INTEGER *feature_indices_ptr
    cdef np.ndarray index
    cdef INTEGER *index_data_ptr
    cdef DOUBLE *sample_weight_data

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight)
    cdef void shuffle(self, seed)
