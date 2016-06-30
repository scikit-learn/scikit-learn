# Author: Nelle Varoquaux, Andrew Tulloch

# Uses the pool adjacent violators algorithm (PAVA), with the
# enhancement of searching for the longest decreasing subsequence to
# pool at each step.

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _isotonic_regression(np.ndarray[DOUBLE, ndim=1] y,
                         np.ndarray[DOUBLE, ndim=1] weight,
                         np.ndarray[DOUBLE, ndim=1] solution):
    cdef:
        DOUBLE numerator, denominator, ratio
        Py_ssize_t i, pooled, n, k

    n = y.shape[0]
    # The algorithm proceeds by iteratively updating the solution
    # array.

    # TODO - should we just pass in a pre-copied solution
    # array and mutate that?
    for i in range(n):
        solution[i] = y[i]

    if n <= 1:
        return solution

    n -= 1
    while 1:
        # repeat until there are no more adjacent violators.
        i = 0
        pooled = 0
        while i < n:
            k = i
            while k < n and solution[k] >= solution[k + 1]:
                k += 1
            if solution[i] != solution[k]:
                # solution[i:k + 1] is a decreasing subsequence, so
                # replace each point in the subsequence with the
                # weighted average of the subsequence.

                # TODO: explore replacing each subsequence with a
                # _single_ weighted point, and reconstruct the whole
                # sequence from the sequence of collapsed points.
                # Theoretically should reduce running time, though
                # initial experiments weren't promising.
                numerator = 0.0
                denominator = 0.0
                for j in range(i, k + 1):
                    numerator += solution[j] * weight[j]
                    denominator += weight[j]
                ratio = numerator / denominator
                for j in range(i, k + 1):
                    solution[j] = ratio
                pooled = 1
            i = k + 1
        # Check for convergence
        if pooled == 0:
            break

    return solution


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _make_unique(np.ndarray[dtype=np.float64_t] X,
                  np.ndarray[dtype=np.float64_t] y,
                  np.ndarray[dtype=np.float64_t] sample_weights):
    """Average targets for duplicate X, drop duplicates.

    Aggregates duplicate X values into a single X value where
    the target y is a (sample_weighted) average of the individual
    targets.

    Assumes that X is ordered, so that all duplicates follow each other.
    """
    unique_values = len(np.unique(X))
    if unique_values == len(X):
        return X, y, sample_weights
    cdef np.ndarray[dtype=np.float64_t] y_out = np.empty(unique_values)
    cdef np.ndarray[dtype=np.float64_t] x_out = np.empty(unique_values)
    cdef np.ndarray[dtype=np.float64_t] weights_out = np.empty(unique_values)

    cdef np.float64_t current_x = X[0]
    cdef np.float64_t current_y = 0
    cdef np.float64_t current_weight = 0
    cdef np.float64_t y_old = 0
    cdef int i = 0
    cdef int current_count = 0
    cdef int j
    cdef np.float64_t x
    cdef int n_samples = len(X)
    for j in range(n_samples):
        x = X[j]
        if x != current_x:
            # next unique value
            x_out[i] = current_x
            weights_out[i] = current_weight / current_count
            y_out[i] = current_y / current_weight
            i += 1
            current_x = x
            current_weight = sample_weights[j]
            current_y = y[j] * sample_weights[j]
            current_count = 1
        else:
            current_weight += sample_weights[j]
            current_y += y[j] * sample_weights[j]
            current_count += 1

    x_out[i] = current_x
    weights_out[i] = current_weight / current_count
    y_out[i] = current_y / current_weight
    return x_out, y_out, weights_out
