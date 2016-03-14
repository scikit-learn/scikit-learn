# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# Licence: BSD 3 clause

cimport cython
from libc.limits cimport INT_MAX
cimport numpy as np
import numpy as np

np.import_array()


cdef class SequentialDataset:
    """Base class for datasets with sequential data access. """

    cdef void next(self, double **x_data_ptr, int **x_ind_ptr,
                   int *nnz, double *y, double *sample_weight) nogil:
        """Get the next example ``x`` from the dataset.

        Parameters
        ----------
        x_data_ptr : double**
            A pointer to the double array which holds the feature
            values of the next example.

        x_ind_ptr : np.intc**
            A pointer to the int array which holds the feature
            indices of the next example.

        nnz : int*
            A pointer to an int holding the number of non-zero
            values of the next example.

        y : double*
            The target value of the next example.

        sample_weight : double*
            The weight of the next example.
        """
        cdef int current_index = self._get_next_index()
        self._sample(x_data_ptr, x_ind_ptr, nnz, y, sample_weight,
                     current_index)

    cdef int random(self, double **x_data_ptr, int **x_ind_ptr,
                    int *nnz, double *y, double *sample_weight) nogil:
        """Get a random example ``x`` from the dataset.

        Parameters
        ----------
        x_data_ptr : double**
            A pointer to the double array which holds the feature
            values of the next example.

        x_ind_ptr : np.intc**
            A pointer to the int array which holds the feature
            indices of the next example.

        nnz : int*
            A pointer to an int holding the number of non-zero
            values of the next example.

        y : double*
            The target value of the next example.

        sample_weight : double*
            The weight of the next example.

        Returns
        -------
        index : int
            The index sampled
        """
        cdef int current_index = self._get_random_index()
        self._sample(x_data_ptr, x_ind_ptr, nnz, y, sample_weight,
                     current_index)
        return current_index

    cdef void shuffle(self, np.uint32_t seed) nogil:
        """Permutes the ordering of examples."""
        # Fisher-Yates shuffle
        cdef int *ind = self.index_data_ptr
        cdef int n = self.n_samples
        cdef unsigned i, j
        for i in range(n - 1):
            j = i + our_rand_r(&seed) % (n - i)
            ind[i], ind[j] = ind[j], ind[i]

    cdef int _get_next_index(self) nogil:
        cdef int current_index = self.current_index
        if current_index >= (self.n_samples - 1):
            current_index = -1

        current_index += 1
        self.current_index = current_index
        return self.current_index

    cdef int _get_random_index(self) nogil:
        cdef int n = self.n_samples
        cdef int current_index = our_rand_r(&self.seed) % n
        self.current_index = current_index
        return current_index

    cdef void _sample(self, double **x_data_ptr, int **x_ind_ptr,
                      int *nnz, double *y, double *sample_weight,
                      int current_index) nogil:
        pass

    def _shuffle_py(self, np.uint32_t seed):
        """python function used for easy testing"""
        self.shuffle(seed)

    def _next_py(self):
        """python function used for easy testing"""
        cdef int current_index = self._get_next_index()
        return self._sample_py(current_index)

    def _random_py(self):
        """python function used for easy testing"""
        cdef int current_index = self._get_random_index()
        return self._sample_py(current_index)

    def _sample_py(self, int current_index):
        """python function used for easy testing"""
        cdef double* x_data_ptr
        cdef int* x_indices_ptr
        cdef int nnz, j
        cdef double y, sample_weight

        # call _sample in cython
        self._sample(&x_data_ptr, &x_indices_ptr, &nnz, &y, &sample_weight,
                     current_index)

        # transform the pointed data in numpy CSR array
        cdef np.ndarray[double, ndim=1] x_data = np.empty(nnz)
        cdef np.ndarray[int, ndim=1] x_indices = np.empty(nnz, dtype=np.int32)
        cdef np.ndarray[int, ndim=1] x_indptr = np.asarray([0, nnz],
                                                           dtype=np.int32)

        for j in range(nnz):
            x_data[j] = x_data_ptr[j]
            x_indices[j] = x_indices_ptr[j]

        cdef int sample_idx = self.index_data_ptr[current_index]

        return (x_data, x_indices, x_indptr), y, sample_weight, sample_idx

cdef class ArrayDataset(SequentialDataset):
    """Dataset backed by a two-dimensional numpy array.

    The dtype of the numpy array is expected to be ``np.float64`` (double)
    and C-style memory layout.
    """

    def __cinit__(self, np.ndarray[double, ndim=2, mode='c'] X,
                  np.ndarray[double, ndim=1, mode='c'] Y,
                  np.ndarray[double, ndim=1, mode='c'] sample_weights,
                  np.uint32_t seed=1):
        """A ``SequentialDataset`` backed by a two-dimensional numpy array.

        Parameters
        ----------
        X : ndarray, dtype=double, ndim=2, mode='c'
            The sample array, of shape(n_samples, n_features)

        Y : ndarray, dtype=double, ndim=1, mode='c'
            The target array, of shape(n_samples, )

        sample_weights : ndarray, dtype=double, ndim=1, mode='c'
            The weight of each sample, of shape(n_samples,)
        """
        if X.shape[0] > INT_MAX or X.shape[1] > INT_MAX:
            raise ValueError("More than %d samples or features not supported;"
                             " got (%d, %d)."
                             % (INT_MAX, X.shape[0], X.shape[1]))

        # keep a reference to the data to prevent garbage collection
        self.X = X
        self.Y = Y
        self.sample_weights = sample_weights

        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]

        cdef np.ndarray[int, ndim=1, mode='c'] feature_indices = \
            np.arange(0, self.n_features, dtype=np.intc)
        self.feature_indices = feature_indices
        self.feature_indices_ptr = <int *> feature_indices.data

        self.current_index = -1
        self.X_stride = X.strides[0] / X.itemsize
        self.X_data_ptr = <double *>X.data
        self.Y_data_ptr = <double *>Y.data
        self.sample_weight_data = <double *>sample_weights.data

        # Use index array for fast shuffling
        cdef np.ndarray[int, ndim=1, mode='c'] index = \
            np.arange(0, self.n_samples, dtype=np.intc)
        self.index = index
        self.index_data_ptr = <int *>index.data
        # seed should not be 0 for our_rand_r
        self.seed = max(seed, 1)

    cdef void _sample(self, double **x_data_ptr, int **x_ind_ptr,
                      int *nnz, double *y, double *sample_weight,
                      int current_index) nogil:
        cdef long long sample_idx = self.index_data_ptr[current_index]
        cdef long long offset = sample_idx * self.X_stride

        y[0] = self.Y_data_ptr[sample_idx]
        x_data_ptr[0] = self.X_data_ptr + offset
        x_ind_ptr[0] = self.feature_indices_ptr
        nnz[0] = self.n_features
        sample_weight[0] = self.sample_weight_data[sample_idx]


cdef class CSRDataset(SequentialDataset):
    """A ``SequentialDataset`` backed by a scipy sparse CSR matrix. """

    def __cinit__(self, np.ndarray[double, ndim=1, mode='c'] X_data,
                  np.ndarray[int, ndim=1, mode='c'] X_indptr,
                  np.ndarray[int, ndim=1, mode='c'] X_indices,
                  np.ndarray[double, ndim=1, mode='c'] Y,
                  np.ndarray[double, ndim=1, mode='c'] sample_weights,
                  np.uint32_t seed=1):
        """Dataset backed by a scipy sparse CSR matrix.

        The feature indices of ``x`` are given by x_ind_ptr[0:nnz].
        The corresponding feature values are given by
        x_data_ptr[0:nnz].

        Parameters
        ----------
        X_data : ndarray, dtype=double, ndim=1, mode='c'
            The data array of the CSR features matrix.

        X_indptr : ndarray, dtype=np.intc, ndim=1, mode='c'
            The index pointer array of the CSR features matrix.

        X_indices : ndarray, dtype=np.intc, ndim=1, mode='c'
            The column indices array of the CSR features matrix.

        Y : ndarray, dtype=double, ndim=1, mode='c'
            The target values.

        sample_weights : ndarray, dtype=double, ndim=1, mode='c'
            The weight of each sample.
        """
        # keep a reference to the data to prevent garbage collection
        self.X_data = X_data
        self.X_indptr = X_indptr
        self.X_indices = X_indices
        self.Y = Y
        self.sample_weights = sample_weights

        self.n_samples = Y.shape[0]
        self.current_index = -1
        self.X_data_ptr = <double *>X_data.data
        self.X_indptr_ptr = <int *>X_indptr.data
        self.X_indices_ptr = <int *>X_indices.data

        self.Y_data_ptr = <double *>Y.data
        self.sample_weight_data = <double *>sample_weights.data

        # Use index array for fast shuffling
        cdef np.ndarray[int, ndim=1, mode='c'] idx = np.arange(self.n_samples,
                                                               dtype=np.intc)
        self.index = idx
        self.index_data_ptr = <int *>idx.data
        # seed should not be 0 for our_rand_r
        self.seed = max(seed, 1)

    cdef void _sample(self, double **x_data_ptr, int **x_ind_ptr,
                      int *nnz, double *y, double *sample_weight,
                      int current_index) nogil:
        cdef long long sample_idx = self.index_data_ptr[current_index]
        cdef long long offset = self.X_indptr_ptr[sample_idx]
        y[0] = self.Y_data_ptr[sample_idx]
        x_data_ptr[0] = self.X_data_ptr + offset
        x_ind_ptr[0] = self.X_indices_ptr + offset
        nnz[0] = self.X_indptr_ptr[sample_idx + 1] - offset
        sample_weight[0] = self.sample_weight_data[sample_idx]


cdef enum:
    RAND_R_MAX = 0x7FFFFFFF


# rand_r replacement using a 32bit XorShift generator
# See http://www.jstatsoft.org/v08/i14/paper for details
# XXX copied over from sklearn/tree/_tree.pyx, should refactor
cdef inline np.uint32_t our_rand_r(np.uint32_t* seed) nogil:
    seed[0] ^= <np.uint32_t>(seed[0] << 13)
    seed[0] ^= <np.uint32_t>(seed[0] >> 17)
    seed[0] ^= <np.uint32_t>(seed[0] << 5)

    return seed[0] % (<np.uint32_t>RAND_R_MAX + 1)
