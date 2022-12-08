cimport numpy as cnp

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    cnp.intp_t *idx,
    cnp.intp_t size,
) nogil

cdef void sort(floating* Xf, cnp.intp_t* samples, cnp.intp_t n) nogil
