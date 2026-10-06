# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._typedefs cimport float64_t, uint8_t, uint32_t

# Bitsets are stored as BITSET_LENGTH words of type BITSET_INNER_DTYPE_C.
# BITSET_INNER_DTYPE_C is uint32_t, so BITSET_INNER_BITS currently resolves to 32.
# Here we keep arithmetic derived from BITSET_INNER_BITS instead of hard-coding
# those values at call sites.
ctypedef uint32_t BITSET_INNER_DTYPE_C

# Note: we have to use an enum because assigning to const in pxd files
# are not allowed, see https://github.com/cython/cython/issues/4369.
cdef enum:
    BITSET_LENGTH = 8
    BITSET_INNER_BITS = sizeof(BITSET_INNER_DTYPE_C) * 8  # =32
    # Total number of representable categories: 8 * 32 = 256.
    N_BITSETS = BITSET_LENGTH * BITSET_INNER_BITS

ctypedef BITSET_INNER_DTYPE_C[BITSET_LENGTH] BITSET_DTYPE_C

cdef inline void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil:
    cdef:
        unsigned int i

    for i in range(BITSET_LENGTH):
        bitset[i] = 0


cdef inline void set_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil:
    bitset[val // BITSET_INNER_BITS] |= (1 << (val % BITSET_INNER_BITS))


cdef inline uint8_t in_bitset(
    const BITSET_INNER_DTYPE_C* bitset, uint8_t val
) noexcept nogil:
    return (bitset[val // BITSET_INNER_BITS] >> (val % BITSET_INNER_BITS)) & 1


cpdef inline uint8_t in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, uint8_t val
) noexcept nogil:
    return (bitset[val // BITSET_INNER_BITS] >> (val % BITSET_INNER_BITS)) & 1


cdef inline uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset, uint8_t val, unsigned int row
) noexcept nogil:
    # Same as above but works on 2d memory views to avoid the creation of 1d
    # memory views. See https://github.com/scikit-learn/scikit-learn/issues/17299
    return (
        bitset[row, val // BITSET_INNER_BITS] >> (val % BITSET_INNER_BITS)
    ) & 1


cpdef inline void set_bitset_memoryview(BITSET_INNER_DTYPE_C[:] bitset, uint8_t val):
    bitset[val // BITSET_INNER_BITS] |= (1 << (val % BITSET_INNER_BITS))
