cimport numpy as np


cdef void add_row_csr(np.ndarray[np.float64_t, ndim=1],
                      np.ndarray[int, ndim=1],
                      np.ndarray[int, ndim=1],
                      int i, np.ndarray[np.float64_t, ndim=1, mode="c"])

cdef inline double double_min(double a, double b):
    return b if b < a else a

cdef inline double double_max(double a, double b):
    return b if b > a else a