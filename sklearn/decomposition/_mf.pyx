import numpy as np
cimport numpy as np
np.import_array()

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INTEGER

def factorize_matrix(
        np.ndarray[DOUBLE, ndim=1] data,
        np.ndarray[INTEGER, ndim=1] row,
        np.ndarray[INTEGER, ndim=1] col,
        np.ndarray[DOUBLE, ndim=2] L,
        np.ndarray[DOUBLE, ndim=2] R,
        np.ndarray[DOUBLE, ndim=1] bias_samples,
        np.ndarray[DOUBLE, ndim=1] bias_features,
        int n_samples,
        int n_features,
        int n_iter,
        int estimated_rank,
        double learning_rate,
        double regularization,
        double bias_learning_rate,
        double bias_regularization,
        random_state,
        int verbose):

    cdef tot_err = 0
    cdef int i, j, k
    cdef double value, estimation, error, tmp

    for current_iter in range(n_iter):
        for value_index in range(data.shape[0]):
            value = data[value_index]
            i = row[value_index]
            j = col[value_index]
            
            estimation = bias_samples[i] + bias_features[j]
            for k in range(estimated_rank):
                estimation += L[i, k] * R[k, j]
            error = value - estimation

            for k in range(estimated_rank):
                tmp = L[i, k] + learning_rate * error * R[k, j] \
                      - regularization * L[i, k]
                R[k, j] = R[k, j] + learning_rate * error *  L[i, k] \
                          - regularization * R[k, j]
                L[i, k] = tmp

            bias_samples[i] = bias_samples[i] + bias_learning_rate * error \
                              - bias_regularization * bias_samples[i]
            bias_features[j] = bias_features[j] + bias_learning_rate * error \
                               - bias_regularization * bias_features[j]

    return L, R, bias_samples, bias_features


def _rmse(
        np.ndarray[DOUBLE, ndim=1] data,
        np.ndarray[INTEGER, ndim=1] row,
        np.ndarray[INTEGER, ndim=1] col,
        np.ndarray[DOUBLE, ndim=2] L,
        np.ndarray[DOUBLE, ndim=2] R):
    
    cdef int i, j, k
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
        
    return np.sqrt(total_error / data.shape[0])
