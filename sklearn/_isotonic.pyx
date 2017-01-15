# Author: Nelle Varoquaux, Andrew Tulloch, Antony Lee

# Uses the pool adjacent violators algorithm (PAVA), with the
# enhancement of searching for the longest decreasing subsequence to
# pool at each step.

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _inplace_contiguous_isotonic_regression(DOUBLE[::1] y, DOUBLE[::1] w):
    cdef:
        Py_ssize_t n = y.shape[0], i, k
        DOUBLE prev_y, sum_wy, sum_w
        Py_ssize_t[::1] target = np.arange(n, dtype=np.intp)

    # target describes a list of blocks.  At any time, if [i..j] (inclusive) is
    # an active block, then target[i] := j and target[j] := i.

    # For "active" indices (block starts):
    # w[i] := sum{w_orig[j], j=[i..target[i]]}
    # y[i] := sum{y_orig[j]*w_orig[j], j=[i..target[i]]} / w[i]

    with nogil:
        i = 0
        while i < n:
            k = target[i] + 1
            if k == n:
                break
            if y[i] < y[k]:
                i = k
                continue
            sum_wy = w[i] * y[i]
            sum_w = w[i]
            while True:
                # We are within a decreasing subsequence.
                prev_y = y[k]
                sum_wy += w[k] * y[k]
                sum_w += w[k]
                k = target[k] + 1
                if k == n or prev_y < y[k]:
                    # Non-singleton decreasing subsequence is finished,
                    # update first entry.
                    y[i] = sum_wy / sum_w
                    w[i] = sum_w
                    target[i] = k - 1
                    target[k - 1] = i
                    if i > 0:
                        # Backtrack if we can.  This makes the algorithm
                        # single-pass and ensures O(n) complexity.
                        i = target[i - 1]
                    # Otherwise, restart from the same point.
                    break
        # Reconstruct the solution.
        i = 0
        while i < n:
            k = target[i] + 1
            y[i + 1 : k] = y[i]
            i = k


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
