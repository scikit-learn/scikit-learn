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
        bint pooled
        DOUBLE prev_y, sum_wy, sum_w
        Py_ssize_t[::1] next_block = np.arange(1, n + 1, dtype=np.intp)
        # Make prev_block of size n+1 to avoid complications from prev_block[n].
        Py_ssize_t[::1] prev_block = np.arange(-1, n, dtype=np.intp)

    # next_block describes a list of "constant" blocks, skipping over redundant
    # indices.
    # For "active" indices (not skipped over):
    # w[i] := sum{w_orig[j], j=[i..next_block[i])}
    # y[i] := sum{y_orig[j]*w_orig[j], j=[i..next_block[i])} / w[i]
    while True:
        # Repeat until there are no more adjacent violators.
        i = 0
        pooled = False
        while i < n:
            k = i
            prev_y = y[k]
            sum_wy = w[k] * y[k]
            sum_w = w[k]
            while True:
                k = next_block[k]
                if k == n or prev_y < y[k]:
                    if k > next_block[i]:
                        # Non-singleton decreasing subsequence is finished,
                        # update first entry.
                        y[i] = sum_wy / sum_w
                        w[i] = sum_w
                        next_block[i] = k
                        prev_block[k] = i
                        if i > 0:
                            # Backtrack if we can.  This is needed to avoid
                            # O(n^2) complexity in pathological cases.
                            i = prev_block[i]
                        else:
                            i = k
                    else:
                        i = k
                    break
                # We are within a decreasing subsequence.
                prev_y = y[k]
                sum_wy += w[k] * y[k]
                sum_w += w[k]
                pooled = True
        # Check for convergence
        if not pooled:
            break

    # Reconstruct the solution.
    i = 0
    while i < n:
        k = next_block[i]
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
