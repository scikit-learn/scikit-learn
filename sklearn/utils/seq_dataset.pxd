"""Dataset abstractions for sequential data access. """

cimport numpy as np


cdef class SequentialDataset:
    cdef Py_ssize_t n_samples

    cdef void next(self, double **x_data_ptr, int **x_ind_ptr,
                   int *nnz, double *y, double *sample_weight) nogil
    cdef void shuffle(self, seed)


cdef class ArrayDataset(SequentialDataset):
    cdef Py_ssize_t n_features
    cdef int current_index
    cdef int stride
    cdef double *X_data_ptr
    cdef double *Y_data_ptr
    cdef np.ndarray feature_indices
    cdef int *feature_indices_ptr
    cdef np.ndarray index
    cdef int *index_data_ptr
    cdef double *sample_weight_data

    cdef void next(self, double **x_data_ptr, int **x_ind_ptr,
                   int *nnz, double *y, double *sample_weight) nogil
    cdef void shuffle(self, seed)


cdef class CSRDataset(SequentialDataset):
    cdef int current_index
    cdef int stride
    cdef double *X_data_ptr
    cdef int *X_indptr_ptr
    cdef int *X_indices_ptr
    cdef double *Y_data_ptr
    cdef np.ndarray feature_indices
    cdef int *feature_indices_ptr
    cdef np.ndarray index
    cdef int *index_data_ptr
    cdef double *sample_weight_data

    cdef void next(self, double **x_data_ptr, int **x_ind_ptr,
                   int *nnz, double *y, double *sample_weight) nogil
    cdef void shuffle(self, seed)
