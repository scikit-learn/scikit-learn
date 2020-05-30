# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
from .common cimport BITSET_INNER_DTYPE_C

cdef inline void init_bitset(BITSET_DTYPE_C bitset) nogil: # OUT
    cdef:
        unsigned int i

    for i in range(8):
        bitset[i] = 0

cdef inline void set_bitset(X_BINNED_DTYPE_C val,
                            BITSET_DTYPE_C bitset) nogil: # OUT
    cdef:
        unsigned int i1 = val // 32
        unsigned int i2 = val % 32

    # It is assumed that val < 256 so that i1 < 8
    bitset[i1] |= (1 << i2)

cdef inline unsigned char in_bitset(X_BINNED_DTYPE_C val,
                                    BITSET_DTYPE_C bitset) nogil:
    cdef:
        unsigned int i1 = val // 32
        unsigned int i2 = val % 32

    return (bitset[i1] >> i2) & 1


def set_bitset_py(X_BINNED_DTYPE_C val, BITSET_INNER_DTYPE_C[:] bitset):
    cdef:
        unsigned int i1 = val // 32
        unsigned int i2 = val % 32

    # It is assumed that val < 256 or i1 < 8
    bitset[i1] |= (1 << i2)
