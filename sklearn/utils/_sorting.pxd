cimport numpy as cnp

from cython cimport floating

cdef void simultaneous_quick_sort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) noexcept nogil

cdef void simultaneous_introsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) noexcept nogil

cdef void simultaneous_heapsort(
    floating* values,
    cnp.intp_t* indices,
    cnp.intp_t size,
) noexcept nogil
