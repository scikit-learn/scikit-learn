#cython: language_level=3
#cython: boundscheck=False, cdivision=True, wraparound=False


import numpy as np
cimport numpy as np
from cython cimport floating
from libc.math cimport fmax


np.import_array()


ctypedef fused INT:
    int
    long


def _euclidean_sparse_dense_exact(floating[::1] X_data,
                                  INT[::1] X_indices,
                                  INT[::1] X_indptr,
                                  np.ndarray[floating, ndim=2, mode='c'] Y,
                                  floating[::1] y_squared_norms):
    cdef:
        int n_samples_X = X_indptr.shape[0] - 1
        int n_samples_Y = Y.shape[0]
        int n_features = Y.shape[1]

        int i, j

        floating[:, ::1] D = np.empty((n_samples_X, n_samples_Y), Y.dtype)

    for i in range(n_samples_X):
        for j in range(n_samples_Y):
            D[i, j] = _euclidean_dense_sparse_exact_1d(
                &Y[j, 0],
                y_squared_norms[j],
                &X_data[X_indptr[i]],
                &X_indices[X_indptr[i]],
                X_indptr[i + 1] - X_indptr[i])

    return np.asarray(D)


cdef floating _euclidean_dense_sparse_exact_1d(floating *x,
                                               floating x_squared_norm,
                                               floating *y_data,
                                               INT *y_indices,
                                               int y_nnz) nogil:
    """Euclidean distance between x dense and y sparse"""
    cdef:
        int i
        floating xi
        floating tmp = 0.0
        floating result = 0.0
        floating partial_x_squared_norm = 0.0
    
    # Split the loop to avoid unsafe compiler auto optimizations
    for i in range(y_nnz):
        xi = x[y_indices[i]]
        partial_x_squared_norm += xi * xi

    for i in range(y_nnz):
        tmp = y_data[i] - x[y_indices[i]]
        result += tmp * tmp 

    result += x_squared_norm - partial_x_squared_norm

    return fmax(result, 0)