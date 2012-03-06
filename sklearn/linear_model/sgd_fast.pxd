# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.

cimport numpy as np

cdef extern from "math.h":
    cdef extern double exp(double x)
    cdef extern double log(double x)
    cdef extern double sqrt(double x)
    cdef extern double pow(double x, double y)

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

# -----------------------------------------
# Headers for Loss Function Extension Types
# -----------------------------------------

cdef class LossFunction:
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Regression(LossFunction):
    cpdef double loss(self,double p, double y)
    cpdef double dloss(self,double p, double y)

cdef class Classification(LossFunction):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class ModifiedHuber(Classification):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Hinge(Classification):
    cdef double threshold
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Log(Classification):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class SquaredLoss(Regression):
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

cdef class Huber(Regression):
    cdef double c
    cpdef double loss(self, double p, double y)
    cpdef double dloss(self, double p, double y)

# -----------------------------------------
# Headers for Weight Vector Extension Types
# -----------------------------------------

cdef class WeightVector:
    cdef np.ndarray w
    cdef DOUBLE *w_data_ptr
    cdef double wscale
    cdef Py_ssize_t n_features
    cdef double add(self,  DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                    int xnnz, double c)
    cdef double dot(self, DOUBLE *x_data_ptr, INTEGER *x_ind_ptr,
                    int xnnz)
    cdef void scale(self, double c)
    cdef void reset_wscale(self)
    cdef double norm(self)

# -----------------------------------------
# Headers for Dataset Abstractions
# -----------------------------------------

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