cimport numpy as np


cdef void add_row_csr(np.ndarray[np.float64_t, ndim=1],
                      np.ndarray[int, ndim=1],
                      np.ndarray[int, ndim=1],
                      int i, np.ndarray[np.float64_t, ndim=1, mode="c"])

