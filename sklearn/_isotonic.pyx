import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _isotonic_regression(np.ndarray[DOUBLE, ndim=1] y,
                         np.ndarray[DOUBLE, ndim=1] weight,
                         np.ndarray[DOUBLE, ndim=1] solution):

    cdef:
        Py_ssize_t current, i
        unsigned int len_active_set
        DOUBLE v, w

    len_active_set = y.shape[0]
    active_set = [[weight[i] * y[i], weight[i], [i, ]]
                  for i in range(len_active_set)]
    current = 0

    while current < len_active_set - 1:
        while current < len_active_set -1 and \
              (active_set[current][0] * active_set[current + 1][1] <= 
               active_set[current][1] * active_set[current + 1][0]):
            current += 1

        if current == len_active_set - 1:
            break

        # merge two groups
        active_set[current][0] += active_set[current + 1][0]
        active_set[current][1] += active_set[current + 1][1]
        active_set[current][2] += active_set[current + 1][2]

        active_set.pop(current + 1)
        len_active_set -= 1
        while current > 0 and \
              (active_set[current - 1][0] * active_set[current][1] > 
               active_set[current - 1][1] * active_set[current][0]):
            current -= 1
            active_set[current][0] += active_set[current + 1][0]
            active_set[current][1] += active_set[current + 1][1]
            active_set[current][2] += active_set[current + 1][2]

            active_set.pop(current + 1)
            len_active_set -= 1

    for v, w, idx in active_set:
        solution[idx] = v / w
    return solution
