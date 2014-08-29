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
