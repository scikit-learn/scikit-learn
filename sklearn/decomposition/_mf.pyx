# Authors: Nicolas Tresegnie <nicolas.tresegnie@gmail.com>,
#          Artem Sobolev
# License: BSD 3 clause

import numpy as np
cimport numpy as np
np.import_array()

from libc.math cimport sqrt

ctypedef np.float64_t DOUBLE_t
ctypedef np.int32_t INTEGER_t

def factorize_matrix(
        np.ndarray[DOUBLE_t, ndim=1] data,
        np.ndarray[INTEGER_t, ndim=1] row,
        np.ndarray[INTEGER_t, ndim=1] col,
        np.ndarray[DOUBLE_t, ndim=2] L,
        np.ndarray[DOUBLE_t, ndim=2] R,
        np.ndarray[DOUBLE_t, ndim=1] bias_samples,
        np.ndarray[DOUBLE_t, ndim=1] bias_features,
        np.npy_intp n_samples,
        np.npy_intp n_features,
        np.npy_intp n_components,
        int n_iter,
        double regularization,
        double bias_regularization,
        double learning_rate,
        object random_state,
        int verbose):

    cdef tot_err = 0
    cdef np.npy_intp i, j, k
    cdef double value, estimation, error, newL

    for current_iter in range(n_iter):
        for value_index in range(data.shape[0]):
            value = data[value_index]
            i = row[value_index]
            j = col[value_index]
            
            estimation = bias_samples[i] + bias_features[j]
            for k in range(n_components):
                estimation += L[i, k] * R[k, j]
            error = value - estimation

            for k in range(n_components):
                newL = L[i, k] - learning_rate \
                          * (error * R[k, j] + regularization * L[i, k])
                R[k, j] = R[k, j] - learning_rate \
                          * (error *  L[i, k] + regularization * R[k, j])
                L[i, k] = newL

            bias_samples[i] = bias_samples[i] + learning_rate \
                              * (error  - bias_regularization * bias_samples[i])
            bias_features[j] = bias_features[j] + learning_rate \
                              * (error - bias_regularization * bias_features[j])

    return L, R, bias_samples, bias_features


def _rmse(
        np.ndarray[DOUBLE_t, ndim=1] data,
        np.ndarray[INTEGER_t, ndim=1] row,
        np.ndarray[INTEGER_t, ndim=1] col,
        np.ndarray[DOUBLE_t, ndim=2] L,
        np.ndarray[DOUBLE_t, ndim=2] R):
    
    cdef np.npy_intp i, j, k
    cdef double value, estimation, error
    cdef total_error = 0
    for value_index in range(data.shape[0]):
        value = data[value_index]
        i = row[value_index]
        j = col[value_index]
        
        estimation = 0
        for k in range(L.shape[1]):
            estimation += L[i, k] * R[k, j]
        error = value - estimation
        
        total_error += error * error
        
    return sqrt(total_error / data.shape[0])
