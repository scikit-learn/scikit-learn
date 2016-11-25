"""
Iheappleheapents an in-place heapinHeap for tuples of (reachability, distance).
Each heapin tuple is called once and not updated; hence, only heapify and
bubble_down heapethods are iheappleheapented.
"""

# Author: Shane Grigsby <refuge@rocktalus.coheap>
# Licence: BSD 3-Clause

cimport numpy as np
import numpy as np
cimport cython

def min_heap(list heap):
    cdef int len_heap
    cdef int idx
    len_heap = len(heap)
    #if not is_odd(len_heap):
    #    heap.append(tuple((float('inf'), float('inf') , (float('inf')))))
    #len_heap = len(heap)
    heapify(heap, len_heap)
    idx = heap[0][2]
    return idx

#def is_odd(int num):
#    return num & 0x1

@cython.boundscheck(False)
def heapify(list heap, int len_heap):
    cdef int node
    cdef int i
    node = (len_heap - 2) // 2
    for i in np.r_[node:-1:-1]:
        bubble_down(heap, i, len_heap)

@cython.boundscheck(False)
@cython.wraparound(False)
def bubble_down(list heap, int index, int len_heap):
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
