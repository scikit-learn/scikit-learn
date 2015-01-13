# Author: x0l <x0l@github.com>
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
cimport cython

from libc.stdlib  cimport malloc, free
from libc.math cimport log
from libc.float cimport DBL_MAX

ctypedef np.float64_t DOUBLE
ctypedef np.int32_t INT

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults

np.import_array()


@cython.cdivision(True)
cdef inline double likelihood(int n, double c) nogil:
    cdef double u = c / n
    return (n - 1) * log((n - 1) / (n - u)) - log(u)


@cython.cdivision(True)
cdef inline int offset(int n, int i, int j) nogil:
    return (((n << 1) - i - 3) * i >> 1) + j - 1


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void _gaussian_factors(int n, double* X, double damping, double* linkage):
    cdef:
        int n_active, i, j, ci, cj, idx, i_max, j_max, offset_j
        double max_delta

        int *active = <int*>malloc(n * sizeof(int))

        int *ns = <int*>malloc(n * sizeof(int))
        double *ls = <double*>malloc(n * sizeof(double))
        double *cs = <double*>malloc(n * sizeof(double))

        int NN = n * (n - 1) / 2
        double *dc = <double*>malloc(NN * sizeof(double))
        double *dl = <double*>malloc(NN * sizeof(double))

    for i in range(n):
        ns[i] = 1
        ls[i] = 0.0
        cs[i] = 1.0
        active[i] = i
        for j in range(i + 1, n):
            idx = offset(n, i, j)
            dc[idx] = 2 * X[n * i + j]
            dl[idx] = likelihood(2, 2.0 - damping + dc[idx])

    # Main loop
    n_active = n
    while n_active > 1:
        # This nested loop is the bottleneck
        max_delta = -DBL_MAX
        for ci in range(n_active):
            idx = offset(n, active[ci], 0)
            for cj in range(ci + 1, n_active):
                if max_delta < dl[idx + active[cj]]:
                    max_delta = dl[idx + active[cj]]
                    offset_j = cj
                    i_max = active[ci]
                    j_max = active[cj]

        # Update everything
        for i in range(offset_j, n_active-1):
            active[i] = active[i+1]

        idx = offset(n, i_max, j_max)

        ns[i_max] += ns[j_max]
        cs[i_max] += cs[j_max] + dc[idx]
        ls[i_max] += ls[j_max] + max_delta

        linkage[4*(n - n_active) + 0] = i_max
        linkage[4*(n - n_active) + 1] = j_max
        linkage[4*(n - n_active) + 2] = max_delta
        linkage[4*(n - n_active) + 3] = ns[i_max]

        n_active -= 1
        for i in range(n_active):
            j = active[i]
            if j == i_max:
                continue

            i0, i1 = i_max, j
            if i_max > j:
                i0, i1 = i1, i0

            j0, j1 = j_max, j
            if j_max > j:
                j0, j1 = j1, j0

            idx = offset(n, i0, i1)

            dc[idx] += dc[offset(n, j0, j1)]

            dl[idx] = likelihood(ns[i_max] + ns[j],
                                 dc[idx] + cs[i_max] + cs[j] - damping)
            dl[idx] -= ls[i_max] + ls[j]

    free(dl)
    free(dc)
    free(cs)
    free(ls)
    free(ns)
    free(active)


def gaussian_factors(np.ndarray[DOUBLE, ndim=2] X, double damping=1e-4):
    cdef int n = X.shape[0]
    cdef np.ndarray[DOUBLE, ndim=2] Z

    Z = np.zeros((n-1, 4), dtype=np.float64)

    _gaussian_factors(n, <double*>X.data, damping, <double*>Z.data)

    for i, p in enumerate(Z):
        v = Z[i + 1:, 0:2]
        v[v == p[0]] = i + n

    labels = np.arange(n)
    best = Z[:, 2].cumsum().argmax()
    for i in range(best + 1):
        labels[labels == int(Z[i, 0])] = i + n
        labels[labels == int(Z[i, 1])] = i + n

    _, labels = np.unique(labels, return_inverse=True)

    return Z, labels