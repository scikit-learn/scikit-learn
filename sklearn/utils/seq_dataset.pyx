# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.

import numpy as np
import collections

cimport numpy as np
cimport cython

cdef extern from "stdlib.h":
    int rand()


cdef class SequentialDataset:
    """Base class for datasets with sequential data access. """

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight):
        """Get the next example ``x`` from the dataset.

        Parameters
        ----------
        x_data_ptr : np.float64**
            A pointer to the double array which holds the feature
            values of the next example.
        x_ind_ptr : np.int32**
            A pointer to the int32 array which holds the feature
            indices of the next example.
        nnz : int*
            A pointer to an int holding the number of non-zero
            values of the next example.
        y : np.float64*
            The target value of the next example.
        sample_weight : np.float64*
            The weight of the next example.
        """
        raise NotImplementedError()

    cdef void shuffle(self, seed):
        """Permutes the ordering of examples.  """
        raise NotImplementedError()


cdef class ArrayDataset(SequentialDataset):
    """Dataset backed by a two-dimensional numpy array.

    The dtype of the numpy array is expected to be ``np.float64``
    and C-style memory layout.
    """

    def __cinit__(self, np.ndarray[DOUBLE, ndim=2, mode='c'] X,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] Y,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] sample_weights):
        """A ``SequentialDataset`` backed by a two-dimensional numpy array.

        Paramters
        ---------
        X : ndarray, dtype=np.float64, ndim=2, mode='c'
            The samples; a two-dimensional c-continuous numpy array of
            dtype np.float64.
        Y : ndarray, dtype=np.float64, ndim=1, mode='c'
            The target values; a one-dimensional c-continuous numpy array of
            dtype np.float64.
        sample_weights : ndarray, dtype=np.float64, ndim=1, mode='c'
            The weight of each sample; a one-dimensional c-continuous numpy
            array of dtype np.float64.
        """
        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] feature_indices = np.arange(0, self.n_features,
                                                              dtype=np.int32)
        self.feature_indices = feature_indices
        self.feature_indices_ptr = <INTEGER *> feature_indices.data
        self.current_index = -1
        self.stride = X.strides[0] / X.strides[1]
        self.X_data_ptr = <DOUBLE *>X.data
        self.Y_data_ptr = <DOUBLE *>Y.data
        self.sample_weight_data = <DOUBLE *>sample_weights.data

        # Use index array for fast shuffling
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] index = np.arange(0, self.n_samples,
                                                    dtype=np.int32)
        self.index = index
        self.index_data_ptr = <INTEGER *> index.data

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight):
        cdef int current_index = self.current_index
        if current_index >= (self.n_samples - 1):
            current_index = -1

        current_index += 1
        cdef int sample_idx = self.index_data_ptr[current_index]
        cdef int offset = sample_idx * self.stride

        y[0] = self.Y_data_ptr[sample_idx]
        x_data_ptr[0] = self.X_data_ptr + offset
        x_ind_ptr[0] = self.feature_indices_ptr
        nnz[0] = self.n_features
        sample_weight[0] = self.sample_weight_data[sample_idx]

        self.current_index = current_index

    cdef void shuffle(self, seed):
        np.random.RandomState(seed).shuffle(self.index)

cdef class CSRDataset(SequentialDataset):
    """A ``SequentialDataset`` backed by a scipy sparse CSR matrix. """

    def __cinit__(self, np.ndarray[DOUBLE, ndim=1, mode='c'] X_data,
                  np.ndarray[INTEGER, ndim=1, mode='c'] X_indptr,
                  np.ndarray[INTEGER, ndim=1, mode='c'] X_indices,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] Y,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] sample_weight):
        """Dataset backed by a scipy sparse CSR matrix.

        The feature indices of ``x`` are given by x_ind_ptr[0:nnz].
        The corresponding feature values are given by
        x_data_ptr[0:nnz].

        Parameters
        ----------
        X_data : ndarray, dtype=np.float64, ndim=1, mode='c'
            The data array of the CSR matrix; a one-dimensional c-continuous
            numpy array of dtype np.float64.
        X_indptr : ndarray, dtype=np.int32, ndim=1, mode='c'
            The index pointer array of the CSR matrix; a one-dimensional
            c-continuous numpy array of dtype np.int32.
        X_indices : ndarray, dtype=np.int32, ndim=1, mode='c'
            The column indices array of the CSR matrix; a one-dimensional
            c-continuous numpy array of dtype np.int32.
        Y : ndarray, dtype=np.float64, ndim=1, mode='c'
            The target values; a one-dimensional c-continuous numpy array of
            dtype np.float64.
        sample_weights : ndarray, dtype=np.float64, ndim=1, mode='c'
            The weight of each sample; a one-dimensional c-continuous numpy
            array of dtype np.float64.
        """
        self.n_samples = Y.shape[0]
        self.current_index = -1
        self.X_data_ptr = <DOUBLE *>X_data.data
        self.X_indptr_ptr = <INTEGER *>X_indptr.data
        self.X_indices_ptr = <INTEGER *>X_indices.data
        self.Y_data_ptr = <DOUBLE *>Y.data
        self.sample_weight_data = <DOUBLE *> sample_weight.data
        # Use index array for fast shuffling
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] index = np.arange(0, self.n_samples,
                                                    dtype=np.int32)
        self.index = index
        self.index_data_ptr = <INTEGER *> index.data

    cdef void next(self, DOUBLE **x_data_ptr, INTEGER **x_ind_ptr,
                   int *nnz, DOUBLE *y, DOUBLE *sample_weight):
        cdef int current_index = self.current_index
        if current_index >= (self.n_samples - 1):
            current_index = -1

        current_index += 1
        cdef int sample_idx = self.index_data_ptr[current_index]
        cdef int offset = self.X_indptr_ptr[sample_idx]
        y[0] = self.Y_data_ptr[sample_idx]
        x_data_ptr[0] = self.X_data_ptr + offset
        x_ind_ptr[0] = self.X_indices_ptr + offset
        nnz[0] = self.X_indptr_ptr[sample_idx + 1] - offset
        sample_weight[0] = self.sample_weight_data[sample_idx]

        self.current_index = current_index

    cdef void shuffle(self, seed):
        np.random.RandomState(seed).shuffle(self.index)


cdef class PairwiseArrayDataset:
    """Dataset backed by a two-dimensional numpy array. 

    Calling next_pair() returns a random pair of examples with disagreeing labels.

    The dtype of the numpy array is expected to be ``np.float64``
    and C-style memory layout.
    """

    cdef void next_pair(self, DOUBLE **a_data_ptr, DOUBLE **b_data_ptr, 
                   INTEGER **x_ind_ptr, int *nnz_a, int *nnz_b, 
                   DOUBLE *y_a, DOUBLE *y_b):
        raise NotImplementedError()


cdef class PairwiseArrayDatasetRoc(PairwiseArrayDataset):
    def __cinit__(self, np.ndarray[DOUBLE, ndim=2, mode='c'] X,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] Y):
        
        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] feature_indices = np.arange(0, self.n_features,
                                                              dtype=np.int32)
        self.feature_indices = feature_indices
        self.feature_indices_ptr = <INTEGER *> feature_indices.data
        self.stride = X.strides[0] / X.strides[1]
        self.X_data_ptr = <DOUBLE *>X.data
        self.Y_data_ptr = <DOUBLE *>Y.data
        
        # Create an index of positives and negatives for fast sampling
        # of disagreeing pairs
        positives = []
        negatives = []
        cdef Py_ssize_t i
        for i in range(self.n_samples):
            if self.Y_data_ptr[i] > 0:
                positives.append(i)
            else:
                negatives.append(i)
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] pos_index = np.array(positives, dtype=np.int32)
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] neg_index = np.array(negatives, dtype=np.int32)
        self.pos_index = pos_index
        self.neg_index = neg_index        
        self.pos_index_data_ptr = <INTEGER *> pos_index.data
        self.neg_index_data_ptr = <INTEGER *> neg_index.data
        self.n_pos_samples = len(pos_index)
        self.n_neg_samples = len(neg_index)

    cdef void next_pair(self, DOUBLE **a_data_ptr, DOUBLE **b_data_ptr, 
                   INTEGER **x_ind_ptr, int *nnz_a, int *nnz_b, 
                   DOUBLE *y_a, DOUBLE *y_b):

        x_ind_ptr[0] = self.feature_indices_ptr
        nnz_a[0] = self.n_features
        nnz_b[0] = self.n_features

        current_pos_index = rand() % self.n_pos_samples
        current_neg_index = rand() % self.n_neg_samples

        # For each step, randomly sample one positive and one negative
        cdef int sample_pos_idx = self.pos_index_data_ptr[current_pos_index]
        cdef int sample_neg_idx = self.neg_index_data_ptr[current_neg_index]
        cdef int pos_offset = sample_pos_idx * self.stride
        cdef int neg_offset = sample_neg_idx * self.stride

        y_a[0] = self.Y_data_ptr[sample_pos_idx]
        y_b[0] = self.Y_data_ptr[sample_neg_idx]
        a_data_ptr[0] = self.X_data_ptr + pos_offset
        b_data_ptr[0] = self.X_data_ptr + neg_offset

cdef class PairwiseArrayDatasetRank(PairwiseArrayDataset):
    def __cinit__(self, np.ndarray[DOUBLE, ndim=2, mode='c'] X,
                  np.ndarray[DOUBLE, ndim=1, mode='c'] Y,
                  int sampling_type, np.ndarray[DOUBLE, ndim=1, mode='c'] query_id):

        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]
        cdef np.ndarray[INTEGER, ndim=1,
                        mode='c'] feature_indices = np.arange(0, self.n_features,
                                                              dtype=np.int32)
        self.feature_indices = feature_indices
        self.feature_indices_ptr = <INTEGER *> feature_indices.data
        self.stride = X.strides[0] / X.strides[1]
        self.X_data_ptr = <DOUBLE *>X.data
        self.Y_data_ptr = <DOUBLE *>Y.data
        
        self.group_id_y_to_count = <dict> collections.Counter(Y)
        cdef Py_ssize_t i
        cdef DOUBLE query_data_ptr_idx
        for i in range(self.n_samples):
            query_data_ptr_idx = self.query_data_ptr[i]
            self.group_id_y_to_index[query_data_ptr_idx][self.Y_data_ptr[i]].append(i)        
        # must declare in the cinit as closures inside cdef functions not yet supported
        self.group_id_y_to_index = <dict> collections.defaultdict(lambda: \
                                                           collections.defaultdict(list))
        self.query_data_ptr = <DOUBLE *>query_id.data

        self.sampling_type = sampling_type

    cdef void next_pair(self, DOUBLE **a_data_ptr, DOUBLE **b_data_ptr, 
                   INTEGER **x_ind_ptr, int *nnz_a, int *nnz_b, 
                   DOUBLE *y_a, DOUBLE *y_b):

        x_ind_ptr[0] = self.feature_indices_ptr
        nnz_a[0] = self.n_features
        nnz_b[0] = self.n_features

        cdef int a_idx = rand() % self.n_samples
        cdef int a_offset = a_idx * self.stride
        a_data_ptr[0] = self.X_data_ptr + a_offset
        cdef DOUBLE group_id = self.query_data_ptr[a_idx]
        y_a[0] = self.Y_data_ptr[a_idx]
        y_to_list = self.group_id_y_to_count[group_id]
        cdef int y_range = self.group_id_y_to_count[group_id] - self.group_id_y_to_index[group_id][y_a[0]].size()

        cdef unsigned int random_int = rand() % y_range
        cdef int b_idx
        cdef int b_offset
        for b_y, idx_list in y_to_list.items():
            if y_a[0] == b_y:
                continue
            if random_int < idx_list.size():
                b_idx = idx_list[random_int]
                b_offset = b_idx * self.stride
                b_data_ptr[0] = self.X_data_ptr + b_offset
                break
