#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

import numpy as np
cimport numpy as np
np.import_array()  # required in order to use C-API

cimport cython
cimport numpy as np
from libc.math cimport fabs, sqrt, exp, cos, pow
from cython cimport floating, integral, numeric

from ._typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t
from ._typedefs import DTYPE, ITYPE


cdef int _simultaneous_sort(
    floating* dist,
    integral* idx,
    integral size
) nogil except -1:
    """
    Perform a recursive quicksort on the dist array, simultaneously
    performing the same swaps on the idx array.

    TODO: test if the following algorithms are better:
      - introselect via std::nth_element
      - heap-sort-like
    """
    cdef:
        integral pivot_idx, i, store_idx
        floating pivot_val

    # in the small-array case, do things efficiently
    if size <= 1:
        pass
    elif size == 2:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
    elif size == 3:
        if dist[0] > dist[1]:
            dual_swap(dist, idx, 0, 1)
        if dist[1] > dist[2]:
            dual_swap(dist, idx, 1, 2)
            if dist[0] > dist[1]:
                dual_swap(dist, idx, 0, 1)
    else:
        # Determine the pivot using the median-of-three rule.
        # The smallest of the three is moved to the beginning of the array,
        # the middle (the pivot value) is moved to the end, and the largest
        # is moved to the pivot index.
        pivot_idx = size // 2
        if dist[0] > dist[size - 1]:
            dual_swap(dist, idx, 0, size - 1)
        if dist[size - 1] > dist[pivot_idx]:
            dual_swap(dist, idx, size - 1, pivot_idx)
            if dist[0] > dist[size - 1]:
                dual_swap(dist, idx, 0, size - 1)
        pivot_val = dist[size - 1]

        # partition indices about pivot.  At the end of this operation,
        # pivot_idx will contain the pivot value, everything to the left
        # will be smaller, and everything to the right will be larger.
        store_idx = 0
        for i in range(size - 1):
            if dist[i] < pivot_val:
                dual_swap(dist, idx, i, store_idx)
                store_idx += 1
        dual_swap(dist, idx, store_idx, size - 1)
        pivot_idx = store_idx

        # recursively sort each side of the pivot
        if pivot_idx > 1:
            _simultaneous_sort(dist, idx, pivot_idx)
        if pivot_idx + 2 < size:
            _simultaneous_sort(dist + pivot_idx + 1,
                               idx + pivot_idx + 1,
                               size - pivot_idx - 1)
    return 0


cdef int _push(
    floating* dist,
    integral* idx,
    integral size,
    floating val,
    integral i_val,
) nogil except -1:
    """push (val, i_val) into the heap (dist, idx) of the given size"""
    cdef:
        integral current_idx, left_child_idx, right_child_idx, swap_idx

    # check if val should be in heap
    if val > dist[0]:
        return 0

    # insert val at position zero
    dist[0] = val
    idx[0] = i_val

    # descend the heap, swapping values until the max heap criterion is met
    current_idx = 0
    while True:
        left_child_idx = 2 * current_idx + 1
        right_child_idx = left_child_idx + 1

        if left_child_idx >= size:
            break
        elif right_child_idx >= size:
            if dist[left_child_idx] > val:
                swap_idx = left_child_idx
            else:
                break
        elif dist[left_child_idx] >= dist[right_child_idx]:
            if val < dist[left_child_idx]:
                swap_idx = left_child_idx
            else:
                break
        else:
            if val < dist[right_child_idx]:
                swap_idx = right_child_idx
            else:
                break

        dist[current_idx] = dist[swap_idx]
        idx[current_idx] = idx[swap_idx]

        current_idx = swap_idx

    dist[current_idx] = val
    idx[current_idx] = i_val

    return 0


cdef class NeighborsHeap:
    """A max-heap structure to keep track of distances/indices of neighbors

    This implements an efficient pre-allocated set of fixed-size heaps
    for chasing neighbors, holding both an index and a distance.
    When any row of the heap is full, adding an additional point will push
    the furthest point off the heap.

    Parameters
    ----------
    n_pts : int
        the number of heaps to use
    n_nbrs : int
        the size of each heap.
    """

    def __cinit__(self):
        self.distances_arr = np.zeros((1, 1), dtype=DTYPE, order='C')
        self.indices_arr = np.zeros((1, 1), dtype=ITYPE, order='C')
        self.distances = self.distances_arr
        self.indices = self.indices_arr

    def __init__(self, n_pts, n_nbrs):
        self.distances_arr = np.full((n_pts, n_nbrs), np.inf, dtype=DTYPE,
                                     order='C')
        self.indices_arr = np.zeros((n_pts, n_nbrs), dtype=ITYPE, order='C')
        self.distances = self.distances_arr
        self.indices = self.indices_arr

    def get_arrays(self, sort=True):
        """Get the arrays of distances and indices within the heap.

        If sort=True, then simultaneously sort the indices and distances,
        so the closer points are listed first.
        """
        if sort:
            self._sort()
        return self.distances_arr, self.indices_arr

    cdef inline DTYPE_t largest(self, ITYPE_t row) nogil except -1:
        """Return the largest distance in the given row"""
        return self.distances[row, 0]

    def push(self, ITYPE_t row, DTYPE_t val, ITYPE_t i_val):
        return self._push(row, val, i_val)

    cdef int _push(self, ITYPE_t row, DTYPE_t val,
                   ITYPE_t i_val) nogil except -1:
        """push (val, i_val) into the given row"""
        cdef ITYPE_t size = self.distances.shape[1]
        cdef DTYPE_t* dist_arr = &self.distances[row, 0]
        cdef ITYPE_t* ind_arr = &self.indices[row, 0]

        return _push(dist_arr, ind_arr, size, val, i_val)


    cdef int _sort(self) except -1:
        """simultaneously sort the distances and indices"""
        cdef DTYPE_t[:, ::1] distances = self.distances
        cdef ITYPE_t[:, ::1] indices = self.indices
        cdef ITYPE_t row
        for row in range(distances.shape[0]):
            _simultaneous_sort(&distances[row, 0],
                               &indices[row, 0],
                               distances.shape[1])
        return 0
