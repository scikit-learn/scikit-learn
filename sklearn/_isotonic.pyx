# Author: Nelle Varoquaux

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
        DOUBLE numerator, denominator
        Py_ssize_t i, pooled, n, k

    n = y.shape[0]
    for i in range(n):
        solution[i] = y[i]

    if n <= 1:
        return solution

    n -= 1
    while 1:
        i = 0
        pooled = 0
        while i < n:
            k = i
            while k < n and solution[k] >= solution[k+1]:
                k += 1
            if solution[i] != solution[k]:
                numerator = 0.0
                denominator = 0.0
                for alpha in range(i, k+1):
                    numerator += solution[alpha] * weight[alpha]
                    denominator += weight[alpha]
                for alpha in range(i, k+1):
                    solution[alpha] = numerator / denominator
                pooled = 1
            i = k + 1

        if pooled == 0:
            break
    return solution
