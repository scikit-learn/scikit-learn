from sklearn.utils._typedefs cimport intp_t

from cython cimport floating

cdef void simultaneous_sort(
    floating* values,
    intp_t* indices,
    intp_t n,
    bint use_three_way_partition=*,
) noexcept nogil
