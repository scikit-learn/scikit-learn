#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

cimport cython
cimport numpy as np
from libc.math cimport fabs, sqrt, exp, cos, pow
from cython cimport floating, integral, numeric

from ._typedefs cimport DTYPE_t, ITYPE_t, DITYPE_t
from ._typedefs import DTYPE, ITYPE

cdef inline void swap(numeric* arr, integral i1, integral i2):
    """swap the values at index i1 and i2 of arr"""
    cdef numeric tmp = arr[i1]
    arr[i1] = arr[i2]
    arr[i2] = tmp

cdef inline void dual_swap(
    floating* dist,
    integral* idx,
    integral i1,
    integral i2
) nogil:
    """swap the values at index i1 and i2 of both dist and idx"""
    cdef:
        floating dtmp = dist[i1]
        integral itmp = idx[i1]

    dist[i1] = dist[i2]
    dist[i2] = dtmp

    idx[i1] = idx[i2]
    idx[i2] = itmp

cdef int _simultaneous_sort(
    floating* dist,
    integral* idx,
    integral size
) nogil except -1

cdef int _push(
    floating* dist,
    integral* idx,
    integral size,
    floating val,
    integral i_val,
) nogil except -1


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
    cdef np.ndarray distances_arr
    cdef np.ndarray indices_arr

    cdef DTYPE_t[:, ::1] distances
    cdef ITYPE_t[:, ::1] indices

    cdef inline DTYPE_t largest(self, ITYPE_t row) nogil except -1

    cdef int _push(self, ITYPE_t row, DTYPE_t val,
                   ITYPE_t i_val) nogil except -1

    cdef int _sort(self) except -1
