cimport numpy as np
import numpy as np
cimport cython

ctypedef np.float64_t DTYPE_t
ctypedef np.int_t DTYPE

@cython.boundscheck(False)
@cython.wraparound(False)
# Checks for smallest reachability distance
# In case of tie, preserves order and returns first instance
# as sorted by distance
cpdef quick_scan(double[:] rdists, double[:] dists):
    cdef Py_ssize_t n
    cdef int idx
    cdef int i
    cdef double rdist
    cdef double dist
    rdist = np.inf
    dist = np.inf
    n = len(rdists)
    for i from 0 <= i < n:
        if rdists[i] < rdist:
            rdist = rdists[i]
            dist = dists[i]
            idx = i
        if rdists[i] == rdist:
            if dists[i] < dist:
                dist = dists[i]
                idx = i
    return idx
