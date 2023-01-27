cimport numpy as cnp

from cython cimport floating

cdef int simultaneous_quick_sort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) nogil

cdef void simultaneous_introsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) nogil
