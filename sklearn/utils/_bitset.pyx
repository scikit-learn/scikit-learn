"""
A bitset is a data structure used to represent sets of integers in [0, n]. For decision
trees, we use them to represent sets of features indices (e.g. features that go to the
left child, or features that are categorical). For familiarity with bitsets and bitwise
operations see:
https://en.wikipedia.org/wiki/Bit_array
https://en.wikipedia.org/wiki/Bitwise_operation
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause


cdef inline void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil:
    cdef:
        unsigned int i

    for i in range(8):
        bitset[i] = 0


cdef inline void set_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil:
    bitset[val // 32] |= (1 << (val % 32))


cdef inline uint8_t in_bitset(BITSET_DTYPE_C bitset, uint8_t val) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


cpdef inline uint8_t in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, uint8_t val
) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


cdef inline uint8_t in_bitset_2d_memoryview(
    const BITSET_INNER_DTYPE_C[:, :] bitset, uint8_t val, unsigned int row
) noexcept nogil:
    # Same as above but works on 2d memory views to avoid the creation of 1d
    # memory views. See https://github.com/scikit-learn/scikit-learn/issues/17299
    return (bitset[row, val // 32] >> (val % 32)) & 1


cpdef inline void set_bitset_memoryview(BITSET_INNER_DTYPE_C[:] bitset, uint8_t val):
    bitset[val // 32] |= (1 << (val % 32))


def set_raw_bitset_from_binned_bitset(
    BITSET_INNER_DTYPE_C[:] raw_bitset,
    BITSET_INNER_DTYPE_C[:] binned_bitset,
    float64_t[:] categories,
):
    """Set the raw_bitset from the values of the binned bitset

    categories is a mapping from binned category value to raw category value.
    """
    cdef:
        int binned_cat_value
        float64_t raw_cat_value

    for binned_cat_value, raw_cat_value in enumerate(categories):
        if in_bitset_memoryview(binned_bitset, binned_cat_value):
            set_bitset_memoryview(raw_bitset, <uint8_t>raw_cat_value)
