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

cdef void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil

cdef void set_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil

cdef uint8_t in_bitset(const BITSET_INNER_DTYPE_C* bitset, uint8_t val) noexcept nogil

cpdef uint8_t in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, uint8_t val
) noexcept nogil

cdef uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset, uint8_t val, unsigned int row
) noexcept nogil
