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


def all_msr(np.ndarray[LONG, ndim=1, mode="c"] rows,
            np.ndarray[LONG, ndim=1, mode="c"] cols,
            np.ndarray[DOUBLE, ndim=2, mode="c"] X):
    cdef int n_rows = rows.shape[0]
    cdef int n_cols = cols.shape[0]

    cdef np.ndarray[DOUBLE, ndim = 1, mode = "c"] row_mean = None
    cdef np.ndarray[DOUBLE, ndim = 1, mode = "c"] col_mean = None
    row_mean = np.zeros(n_rows, dtype=np.float, order="c")
    col_mean = np.zeros(n_cols, dtype=np.float, order="c")
    cdef double arr_mean = 0.0

    cdef double val = 0.0

    for i in range(n_rows):
        for j in range(n_cols):
            val = X[rows[i], cols[j]]
            row_mean[i] += val
            col_mean[j] += val
            arr_mean += val

    for i in range(n_rows):
        row_mean[i] /= n_cols
    for j in range(n_cols):
        col_mean[j] /= n_rows
    arr_mean /= (n_rows * n_cols)

    cdef np.ndarray[DOUBLE, ndim = 1, mode = "c"] row_msr = None
    cdef np.ndarray[DOUBLE, ndim = 1, mode = "c"] col_msr = None
    row_msr = np.zeros(n_rows, dtype=np.float, order="c")
    col_msr = np.zeros(n_cols, dtype=np.float, order="c")

    cdef double msr = 0.0

    for i in range(n_rows):
        for j in range(n_cols):
            val = (X[rows[i], cols[j]] - row_mean[i] - col_mean[j] + arr_mean)
            val = val * val
            row_msr[i] += val
            col_msr[j] += val
            msr += val
    for i in range(n_rows):
        row_msr[i] /= n_cols
    for j in range(n_cols):
        col_msr[j] /= n_rows
    msr /= (n_rows * n_cols)

    return msr, row_msr, col_msr
