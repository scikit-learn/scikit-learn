cimport numpy as cnp
from ._typedefs cimport ITYPE_t

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil

cdef void sort(floating* Xf, cnp.npy_intp* samples, cnp.npy_intp n) nogil
