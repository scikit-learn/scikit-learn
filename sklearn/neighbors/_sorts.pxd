# cython: language_level=3
cimport numpy as np

from ._typedefs cimport DTYPE_t, ITYPE_t

cdef int intro_select(
        DTYPE_t *data,
        ITYPE_t *indices,
        ITYPE_t pivot,
        ITYPE_t n_points) nogil except -1

cpdef np.ndarray[ITYPE_t, ndim=2, mode='c'] argpartition(
        np.ndarray[DTYPE_t, ndim=2, mode='c'] data,
        ITYPE_t pivot)