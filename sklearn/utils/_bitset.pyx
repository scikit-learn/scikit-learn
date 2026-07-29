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
