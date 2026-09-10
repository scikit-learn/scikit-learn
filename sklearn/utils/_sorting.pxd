from libc.stdint cimport int32_t

from sklearn.utils._typedefs cimport intp_t

from cython cimport floating

ctypedef fused index_t:
    intp_t
    int32_t

cdef void simultaneous_sort(
    floating* values,
    index_t* indices,
    intp_t n,
    bint use_three_way_partition=*,
) noexcept nogil
