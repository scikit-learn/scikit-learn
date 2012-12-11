"""Dataset abstractions for sequential data access. """

cimport numpy as np

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER


cdef class SequentialDataset:
    cdef Py_ssize_t n_samples

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight)
    cdef void shuffle(self, seed)


cdef class ArrayDataset(SequentialDataset):
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


cdef class CSRDataset(SequentialDataset):
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


cdef class PairwiseDataset:
    cdef Py_ssize_t n_samples

    cdef void next_pair(self, DOUBLE **a_data_ptr, DOUBLE **b_data_ptr, 
                   INTEGER **a_ind_ptr, INTEGER **b_ind_ptr,
                   int *nnz_a, int *nnz_b, DOUBLE *y_a, DOUBLE *y_b,
                   DOUBLE *sample_weight_a, DOUBLE *sample_weight_b)


cdef class PairwiseRocDataset(PairwiseDataset):
    cdef DOUBLE *Y_data_ptr
    cdef INTEGER[::1] pos_index
    cdef INTEGER[::1] neg_index
    cdef int n_pos_samples
    cdef int n_neg_samples

    cdef void init_roc_index(self)
    cdef void draw_roc_sample(self, INTEGER *a_idx, INTEGER *b_idx,
                          DOUBLE *y_a, DOUBLE *y_b)


cdef class PairwiseArrayDatasetRoc(PairwiseRocDataset):
    cdef Py_ssize_t n_features
    cdef DOUBLE *X_data_ptr
    cdef int stride
    cdef np.ndarray feature_indices
    cdef INTEGER *feature_indices_ptr
    cdef DOUBLE *sample_weight_data


cdef class PairwiseCSRDatasetRoc(PairwiseRocDataset):
    cdef DOUBLE *X_data_ptr
    cdef INTEGER *X_indptr_ptr
    cdef INTEGER *X_indices_ptr
    cdef DOUBLE *sample_weight_data


cdef class PairwiseRankDataset(PairwiseDataset):
    cdef DOUBLE *Y_data_ptr


cdef class PairwiseArrayDatasetRank(PairwiseRankDataset):
    cdef Py_ssize_t n_features
    cdef DOUBLE *X_data_ptr
    cdef int stride
    cdef np.ndarray feature_indices
    cdef INTEGER *feature_indices_ptr
    cdef DOUBLE *sample_weight_data

    cdef dict group_id_y_to_index
    cdef dict group_id_y_to_count
    cdef DOUBLE *query_data_ptr    


cdef class PairwiseCSRDatasetRank(PairwiseRankDataset):
    cdef DOUBLE *X_data_ptr
    cdef INTEGER *X_indptr_ptr
    cdef INTEGER *X_indices_ptr
    cdef DOUBLE *sample_weight_data

    cdef dict group_id_y_to_index
    cdef dict group_id_y_to_count
    cdef DOUBLE *query_data_ptr    
