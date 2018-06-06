#!python
#cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
import numpy as np
cimport numpy as np
from numpy cimport ndarray as nd_arr
from cython.parallel cimport prange

# faster than np.sqrt
from libc.math cimport sqrt


'''
Port struct defined in _mds_pertubations.h
Contains:
   - The index k of the point that is moved
   - The step that this point is moved
   - The error after the move
'''
cdef extern from "src/_mds_pertubations.h":
    cdef struct pertub_res:
        int k
        double step
        double error

    ctypedef pertub_res pertub_res_t

'''
Get the pertubation that yields the best error. Defined in _mds_pertubations.c
'''
cdef extern pertub_res min_pertub_error(
        double* xs, double radius, double* d_current,
        double* d_goal, int ii, int x_rows, int x_cols,
        double percent, int n_jobs)


def distance_matrix(double[:, :] A):

    cdef:
        Py_ssize_t nrow = A.shape[0]
        Py_ssize_t ncol = A.shape[1]
        Py_ssize_t ii, jj, kk
        np.ndarray[np.float64_t, ndim=2] D = np.zeros((nrow, nrow), np.double)
        double tmpss, diff

    for ii in range(nrow):
        for jj in range(ii + 1, nrow):
            tmpss = 0
            for kk in range(ncol):
                diff = A[ii, kk] - A[jj, kk]
                tmpss += diff * diff
            tmpss = sqrt(tmpss)
            D[ii, jj] = tmpss
            D[jj, ii] = tmpss
    return D


cpdef double mse(nd_arr[np.float64_t, ndim=1] d_goal, nd_arr[np.float64_t, ndim=1] d):
    cdef:
        Py_ssize_t N = d.shape[0]
        Py_ssize_t ii = 0

        double s = 0, diff = 0

    for ii in range(N):
        diff = d_goal[ii] - d[ii]
        s += diff * diff
    return s


cpdef double mse2(nd_arr[np.float64_t, ndim=2] d_goal, nd_arr[np.float64_t, ndim=2] d):
    cdef:
        Py_ssize_t N = d.shape[0]
        Py_ssize_t ii = 0

        double s = 0, diff = 0

    for ii in range(N):
        for jj in range(ii + 1):
            diff = d_goal[ii, jj] - d[ii, jj]
            s += diff * diff
    return s


cpdef nd_arr[np.float64_t, ndim=2] update_distance_matrix(
        nd_arr[np.float64_t, ndim=2] xs,
        nd_arr[np.float64_t, ndim=2] d_current,
        int ii,
        double optimum_step,
        int optimum_k):
    cdef:
        Py_ssize_t N = d_current.shape[0]
        Py_ssize_t jj = 0
        double d = 0

    for jj in range(N):
        if ii != jj:
            d = sqrt(d_current[ii, jj] ** 2 -
                (xs[ii, optimum_k] - xs[jj, optimum_k]) ** 2 +
                (xs[ii, optimum_k] + optimum_step - xs[jj, optimum_k]) ** 2)
            d_current[ii, jj] = d
            d_current[jj, ii] = d
    d_current[ii, ii] = 0
    return d_current


cpdef (double, int, double) c_pertub_error(
    nd_arr[np.float64_t, ndim=2, mode="c"] xs,
    double radius,
    nd_arr[np.float64_t, ndim=2, mode="c"] d_current,
    nd_arr[np.float64_t, ndim=2, mode="c"] d_goal,
    int ii,
    double percent=.5,
    int n_jobs=1):
    cdef:
        int x_cols = xs.shape[1]
        int x_rows = xs.shape[0]
        int jj, kk, ll

        double e = 0
        double d_temp = 0

    cdef:
        pertub_res optimum = min_pertub_error(
            &xs[0,0], radius, &d_current[0, 0], &d_goal[0, 0],
            ii, x_rows, x_cols, percent, n_jobs)
        double optimum_error = optimum.error
        int optimum_k = optimum.k
        double optimum_step = optimum.step

    return optimum_error, optimum_k, optimum_step
