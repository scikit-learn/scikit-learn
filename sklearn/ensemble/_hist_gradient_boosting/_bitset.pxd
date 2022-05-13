from .common cimport BITSET_DTYPE_C, BITSET_INNER_DTYPE_C, X_BINNED_DTYPE_C, X_DTYPE_C


cdef void init_bitset(BITSET_DTYPE_C bitset) nogil

cdef void set_bitset(BITSET_DTYPE_C bitset, X_BINNED_DTYPE_C val) nogil

cdef unsigned char in_bitset(BITSET_DTYPE_C bitset, X_BINNED_DTYPE_C val) nogil

cpdef unsigned char in_bitset_memoryview(const BITSET_INNER_DTYPE_C[:] bitset,
                                         X_BINNED_DTYPE_C val) nogil

cdef unsigned char in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C [:, :] bitset,
    X_BINNED_DTYPE_C val,
    unsigned int row) nogil
