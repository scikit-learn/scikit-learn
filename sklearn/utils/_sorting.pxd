from cython cimport floating

from ._typedefs cimport intp_t

cdef void simultaneous_quicksort(
    floating* values,
    intp_t* indices,
    intp_t size,
) noexcept nogil

cdef void simultaneous_introsort(
    floating* values,
    intp_t* indices,
    intp_t size,
) noexcept nogil

cdef void simultaneous_heapsort(
    floating* values,
    intp_t* indices,
    intp_t size,
) noexcept nogil
