# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.utils._typedefs cimport float64_t, uint8_t, uint32_t

ctypedef uint32_t BITSET_INNER_DTYPE_C
ctypedef BITSET_INNER_DTYPE_C[8] BITSET_DTYPE_C

cdef void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil

cdef void set_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil

cdef uint8_t in_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil

cpdef uint8_t in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, uint8_t val
) noexcept nogil

cdef uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset, uint8_t val, unsigned int row
) noexcept nogil
