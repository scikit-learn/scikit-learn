# cython: language_level=3
from cython cimport floating
cimport numpy as np

cdef int intro_select(
        floating *data,
        int *indices,
        int pivot,
        int n_points) nogil except -1

cpdef np.ndarray[int, ndim=2, mode='c'] argpartition(
        np.ndarray[floating, ndim=2, mode='c'] data,
        int pivot)