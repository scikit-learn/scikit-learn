cimport numpy as cnp
from ._typedefs cimport ITYPE_t

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    ITYPE_t *idx,
    ITYPE_t size,
) nogil

cdef void sort(cnp.npy_float32* Xf, cnp.npy_intp* samples, cnp.npy_intp n) nogil

cdef void swap(
    cnp.npy_float32* Xf,
    cnp.npy_intp* samples,
    cnp.npy_intp i,
    cnp.npy_intp j,
) nogil

cdef cnp.npy_float32 median3(cnp.npy_float32* Xf, cnp.npy_intp n) nogil

cdef void introsort(cnp.npy_float32* Xf, cnp.npy_intp *samples, cnp.npy_intp n, int maxd) nogil

cdef void sift_down(
    cnp.npy_float32* Xf,
    cnp.npy_intp* samples,
    cnp.npy_intp start,
    cnp.npy_intp end,
) nogil

cdef void heapsort(cnp.npy_float32* Xf, cnp.npy_intp* samples, cnp.npy_intp n) nogil
