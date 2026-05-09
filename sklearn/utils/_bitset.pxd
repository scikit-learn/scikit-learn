# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._typedefs cimport float64_t, uint8_t, uint32_t

# Fixed-width bitset storage shared by categorical split implementations.
# e.g. in HGBT and DT implmenetations.
#
# A bitset stores BITSET_LENGTH words, with one bit per supported category code.
# BITSET_INNER_BITS is derived from the word dtype so changing
# BITSET_INNER_DTYPE_C updates bit addressing automatically.
#
# BITSET_LENGTH is intentionally kept as a literal enum value: Cython does not
# allow deriving it from sizeof(BITSET_DTYPE_C) in a cdef enum. The total number
# of representable categories is BITSET_LENGTH * BITSET_INNER_BITS.
ctypedef uint32_t BITSET_INNER_DTYPE_C
cdef enum:
    BITSET_INNER_BITS = sizeof(BITSET_INNER_DTYPE_C) * 8
    BITSET_LENGTH = 8
ctypedef BITSET_INNER_DTYPE_C[BITSET_LENGTH] BITSET_DTYPE_C


cdef void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil

cdef void set_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil

cdef uint8_t in_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil

cpdef uint8_t in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, uint8_t val
) noexcept nogil

cdef uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset, uint8_t val, unsigned int row
) noexcept nogil
