# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
#
# Licence: BSD 3 clause

from libc.stdint cimport int64_t, int32_t
cimport numpy as np
import numpy as np
import cython


np.import_array()

ctypedef np.float64_t double_t
ctypedef float [:, :] float_array_2d_t
ctypedef double [:, :] double_array_2d_t

cdef fused floating_array_2d_t:
    float_array_2d_t
    double_array_2d_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _chi2_kernel_fast(floating_array_2d_t X,
                      floating_array_2d_t Y,
                      floating_array_2d_t result):
    cdef int i, j, k
    cdef int n_samples_X = X.shape[0]
    cdef int n_samples_Y = Y.shape[0]
    cdef int n_features = X.shape[1]
    cdef double res, nom, denom
    for i in xrange(n_samples_X):
        for j in xrange(n_samples_Y):
            res = 0
            for k in xrange(n_features):
                denom = (X[i, k] - Y[j, k])
                nom = (X[i, k] + Y[j, k])
                if nom != 0:
                    res  += denom * denom / nom
            result[i, j] = -res


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _manhattan_distances_coo(X, Y):
    cdef:
        int i, j, r_, r, idx
        int n_x = len(X.data)
        int n_y = len(Y.data)
        int n_samples_Y = Y.shape[0]
        int n_samples_X = X.shape[0]

        np.ndarray[double_t, ndim=1] xdata = X.data
        np.ndarray[int, ndim=1] xrow = X.row
        np.ndarray[int, ndim=1] xcol = X.col
        np.ndarray[double_t, ndim=1] ydata = Y.data
        np.ndarray[int, ndim=1] yrow = Y.row
        np.ndarray[int, ndim=1] ycol = Y.col

        int64_t nnz = X.nnz * Y.shape[0] + Y.nnz * X.shape[0]
        np.ndarray[double_t, ndim=1] data = np.empty(nnz, dtype=np.float64)
        np.ndarray[int32_t, ndim=1] rows = np.empty(nnz, dtype=np.int32)
        np.ndarray[int32_t, ndim=1] cols = np.empty(nnz, dtype=np.int32)

    idx = 0
    for i in xrange(n_x):
        r_ = xrow[i] * n_samples_Y
        for j in xrange(n_samples_Y):
            r = r_ + j
            rows[idx] = r
            cols[idx] = xcol[i]
            data[idx] = xdata[i]
            idx += 1

    for i in xrange(n_y):
        for j in xrange(n_samples_X):
            r = j * n_samples_Y + yrow[i]
            rows[idx] = r
            cols[idx] = ycol[i]
            data[idx] = -ydata[i]
            idx += 1
    return data, rows, cols


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def _manhattan_distances_csr(X, Y):
    cdef:
        int i, j, k, r_, r, idx, xi, yi
        int xstart, xstop, ystart, ystop, xcol, ycol
        int n_x = len(X.data)
        int n_y = len(Y.data)
        int n_samples_Y = Y.shape[0]
        int n_samples_X = X.shape[0]

        np.ndarray[np.float64_t, ndim=1] xdata = X.data
        np.ndarray[int, ndim=1] x_indptr = X.indptr
        np.ndarray[int, ndim=1] x_indices = X.indices
        np.ndarray[np.float64_t, ndim=1] ydata = Y.data
        np.ndarray[int, ndim=1] y_indptr = Y.indptr
        np.ndarray[int, ndim=1] y_indices = Y.indices

        np.int64_t nnz = X.nnz * Y.shape[0] + Y.nnz * X.shape[0]
        int ptr_size = X.shape[0] * Y.shape[0] + 1
        np.ndarray[np.float64_t, ndim=1] data = np.empty(nnz, dtype=np.float64)
        np.ndarray[np.int32_t, ndim=1] indices = np.empty(nnz, dtype=np.int32)
        np.ndarray[np.int32_t, ndim=1] indptr = np.empty(ptr_size,
                                                         dtype=np.int32)

    idx = 0
    indptr[0] = 0
    for i in xrange(n_samples_X):
        r_ = i * n_samples_Y
        xstart, xstop = x_indptr[i], x_indptr[i+1]
        for j in xrange(n_samples_Y):
            ystart, ystop = y_indptr[j], y_indptr[j+1]
            r = r_ + j
            nnz = (xstop - xstart) + (ystop - ystart)
            k, xi, yi = 0, xstart, ystart
            xcol , ycol = -1, -1
            # insertion sort for the indices of row r
            while k < nnz:
                if xstart <= xi < xstop:
                    xcol = x_indices[xi]

                if ystart <= yi < ystop:
                    ycol = y_indices[yi]

                if xcol == ycol:
                    data[idx] = xdata[xi] - ydata[yi]
                    indices[idx] = xcol
                    xi += 1
                    yi += 1
                    nnz -= 1  # shared column, do one less iteration
                elif (xi < xstop and xcol < ycol) or yi == ystop:
                    data[idx] = xdata[xi]
                    indices[idx] = xcol
                    xi += 1
                elif (yi < ystop and ycol < xcol) or xi == xstop:
                    data[idx] = ydata[yi]
                    indices[idx] = ycol
                    yi += 1

                k += 1
                idx += 1
            indptr[r+1] = idx

    return data, indices, indptr
