# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

import numpy as np

cimport numpy as np
cimport cython

np.import_array()

ctypedef np.float64_t DOUBLE
ctypedef np.int64_t LONG


def compute_msr(long[:] rows,
                long[:] cols,
                double[:] row_mean,
                double[:] col_mean,
                double arr_mean,
                double[:, :] X):
    cdef long n_rows = rows.shape[0]
    cdef long n_cols = cols.shape[0]

    row_msr = np.zeros(n_rows, dtype=np.float, order="c")
    col_msr = np.zeros(n_cols, dtype=np.float, order="c")
    cdef double[:] row_msr_view = row_msr
    cdef double[:] col_msr_view = col_msr

    cdef double msr = 0.0
    cdef double val = 0.0

    cdef long i
    cdef long j

    with nogil:
        for i in range(n_rows):
            for j in range(n_cols):
                val = (X[rows[i], cols[j]] - row_mean[i] - col_mean[j] + arr_mean)
                val = val * val
                row_msr_view[i] += val
                col_msr_view[j] += val
                msr += val
        for i in range(n_rows):
            row_msr_view[i] /= n_cols
        for j in range(n_cols):
            col_msr_view[j] /= n_rows
        msr /= (n_rows * n_cols)

    return msr, row_msr, col_msr
