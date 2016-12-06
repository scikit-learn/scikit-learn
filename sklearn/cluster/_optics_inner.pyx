
cimport numpy as np
import numpy as np
cimport cython

ctypedef np.float64_t DTYPE_t
ctypedef np.int_t DTYPE

cpdef min_heap(double[:] rdists, double[:] dists, long[:] indices):
    cdef int len_heap
    cdef long idx
    len_heap = len(rdists)
    heapify(rdists, dists, indices, len_heap)
    idx = indices[0]
    return idx

@cython.boundscheck(False)
cpdef heapify(double[:] rdists, double[:] dists, 
            long[:] indices, int len_heap):
    cdef int node
    cdef int i
    node = (len_heap - 2) // 2
    for i in np.r_[node:-1:-1]:
        bubble_down(rdists, dists, indices, i, len_heap)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bubble_down(double[:] rdists, double[:] dists, 
                long[:] indices, int index, int len_heap):
    cdef bint go
    cdef int left
    cdef int right
    cdef double Crdist
    cdef double Cdist
    cdef double Lrdist
    cdef double Ldist
    cdef double Rrdist
    cdef double Rdist
    cdef tuple C
    cdef tuple L
    cdef tuple R
    go = True
    left = index*2 + 1
    if left == len_heap - 1:
        right = left
    else:
        right = index*2 + 2
    while go:
        Crdist = rdists[index]
        Cdist = dists[index]
        C = tuple((Crdist,Cdist))
        Lrdist = rdists[left]
        Ldist = dists[left]
        L = tuple((Lrdist,Ldist))
        Rrdist = rdists[right]
        Rdist = dists[right]
        R = tuple((Rrdist,Rdist))
        if C > (L or R):
            if L <= R:
                rdists[index], rdists[left] = rdists[left], rdists[index]
                dists[index], dists[left] = dists[left], dists[index]
                indices[index], indices[left] = indices[left], indices[index]
                index = left
            else:
                rdists[index], rdists[right] = rdists[right], rdists[index]
                dists[index], dists[right] = dists[right], dists[index]
                indices[index], indices[right] = indices[right], indices[index]
                index = right
            go = index < (len_heap // 2) - 1
            left = index*2 + 1
            if left == len_heap - 1:
                right = left
            else:
                right = index*2 + 2
        else:
            go = False
