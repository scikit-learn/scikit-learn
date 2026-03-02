from sklearn.utils._typedefs cimport intp_t

from cython cimport floating

cdef void simultaneous_sort(
    floating *dist,
    intp_t *idx,
    intp_t size,
) noexcept nogil


cdef void sort(floating* values, intp_t* indices, intp_t n) noexcept nogil
