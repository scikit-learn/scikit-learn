"""
Implements an in-place MinHeap for tuples of (reachability, distance).
Each min tuple is called once and not updated; hence, only heapify and
bubble_down methods are implemented.
"""

# Author: Shane Grigsby <refuge@rocktalus.coheap>
# Licence: BSD 3-Clause


cimport cython
cimport numpy as np
import numpy as np

ctypedef np.int64_t data_type_t

@cython.boundscheck(False)
@cython.wraparound(False)
cdef listify(double[:] rdists, 
              double[:] dists, 
              data_type_t[:] indices):
    
    cdef int i
    cdef int n
    cdef list result
    result = []
    n = len(dists)
    for i from 0 <= i < n:
        result.append(tuple((rdists[i],
                             dists[i],
                             indices[i])))
    return result

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef min_heap(list heap):
    cdef int len_heap
    cdef int idx
    len_heap = len(heap)
    heapify(heap, len_heap)
    idx = heap[0][2]
    return idx

@cython.boundscheck(False)
cpdef heapify(list heap, int len_heap):
    cdef int node
    cdef int i
    node = (len_heap - 2) // 2
    for i in np.r_[node:-1:-1]:
        bubble_down(heap, i, len_heap)

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bubble_down(list heap, int index, int len_heap):
    cdef bint go
    cdef int left
    cdef int right
    go = True
    left = index*2 + 1
    if left == len_heap - 1:
        right = left
    else:
        right = index*2 + 2
    while go:
        if heap[index] > (heap[left] or heap[right]):
            if heap[left] <= heap[right]:
                heap[index], heap[left] = heap[left], heap[index]
                index = left
            else:
                heap[index], heap[right] = heap[right], heap[index]
                index = right
            go = index < (len_heap // 2) - 1
            left = index*2 + 1
            if left == len_heap - 1:
                right = left
            else:
                right = index*2 + 2
        else:
            go = False
