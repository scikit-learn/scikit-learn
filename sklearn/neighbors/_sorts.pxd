# cython: language_level=3
cimport numpy as np
from cython cimport floating, integral

cdef integral intro_select(
        floating *data,
        integral *indices,
        integral pivot,
        integral n_points,
) nogil

cpdef np.ndarray[integral, ndim=2, mode='c'] argpartition(
        np.ndarray[floating, ndim=2, mode='c'] data,
        integral pivot,
)