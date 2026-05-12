from sklearn.utils._typedefs cimport intp_t

from cython cimport floating

cdef int simultaneous_sort(
    floating *dist,
    intp_t *idx,
    intp_t size,
) noexcept nogil

cdef void sort(
    floating* feature_values,
    intp_t* samples,
    intp_t n,
) noexcept nogil
