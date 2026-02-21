from ._typedefs cimport uint8_t, int32_t, uint64_t


# A bitset is a data structure used to represent sets of integers in [0, n]. We
# use them to represent sets of features indices (e.g. features that go to the
# left child, or features that are categorical). For familiarity with bitsets
# and bitwise operations:
# https://en.wikipedia.org/wiki/Bit_array
# https://en.wikipedia.org/wiki/Bitwise_operation


cdef inline void init_bitset(BITSET_DTYPE_C bitset) noexcept nogil:  # OUT
    cdef:
        unsigned int i

    for i in range(8):
        bitset[i] = 0


cdef inline void set_bitset(BITSET_DTYPE_C bitset,  # OUT
                            X_BINNED_DTYPE_C val) noexcept nogil:
    bitset[val // 32] |= (1 << (val % 32))


cdef inline void set_bitset_memoryview(BITSET_INNER_DTYPE_C[:] bitset,  # OUT
                                       X_BINNED_DTYPE_C val) noexcept nogil:
    bitset[val // 32] |= (1 << (val % 32))


cdef inline uint8_t in_bitset(BITSET_DTYPE_C bitset,
                              X_BINNED_DTYPE_C val) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


cdef inline uint8_t in_bitset_memoryview(const BITSET_INNER_DTYPE_C[:] bitset,
                                         X_BINNED_DTYPE_C val) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


cpdef inline uint8_t py_in_bitset_memoryview(const BITSET_INNER_DTYPE_C[:] bitset,
                                             X_BINNED_DTYPE_C val) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


cdef inline uint8_t in_bitset_2d_memoryview(const BITSET_INNER_DTYPE_C[:, :] bitset,
                                            X_BINNED_DTYPE_C val,
                                            unsigned int row) noexcept nogil:
    # Same as above but works on 2d memory views to avoid the creation of 1d
    # memory views. See https://github.com/scikit-learn/scikit-learn/issues/17299
    return (bitset[row, val // 32] >> (val % 32)) & 1


cpdef inline void py_set_bitset_memoryview(BITSET_INNER_DTYPE_C[:] bitset,  # OUT
                                        X_BINNED_DTYPE_C val):
    bitset[val // 32] |= (1 << (val % 32))


def set_raw_bitset_from_binned_bitset(BITSET_INNER_DTYPE_C[:] raw_bitset,  # OUT
                                      BITSET_INNER_DTYPE_C[:] binned_bitset,
                                      X_DTYPE_C[:] categories):
    """Set the raw_bitset from the values of the binned bitset

    categories is a mapping from binned category value to raw category value.
    """
    cdef:
        int binned_cat_value
        X_DTYPE_C raw_cat_value

    for binned_cat_value, raw_cat_value in enumerate(categories):
        if py_in_bitset_memoryview(binned_bitset, binned_cat_value):
            py_set_bitset_memoryview(raw_bitset, <X_BINNED_DTYPE_C>raw_cat_value)


cdef inline void flip_bitset_up_to_nbits(
    BITSET_INNER_DTYPE_C[:] bitset,
    int32_t n_categories
) noexcept nogil:
    cdef int32_t full_words = n_categories // 32
    cdef int32_t bits_in_last = n_categories % 32
    cdef int i

    for i in range(full_words):
        bitset[i] = ~bitset[i]

    if bits_in_last:
        # mask off upper unused bits in last word
        bitset[full_words] = (~bitset[full_words]) & ((1 << bits_in_last) - 1)


cdef inline void bs_from_template_memoryview(
    BITSET_INNER_DTYPE_C[:] bitset,   # OUT: a length‐8 uint32_t memoryview
    BITSET_INNER_DTYPE_C   template,  # template is now a 32-bit mask over “present categories”
    intp_t[:]             cat_offs,   # e.g. [2, 5, 7] if you have 3 present categories
    intp_t                ncats_present
) noexcept nogil:
    """
    For i in [0..ncats_present-1]:
      If (template & (1<<i)) != 0, then set bit at index cat_offs[i] in the 256-bit bitset.
    """
    cdef intp_t i, target_bit, word_idx, bit_in_word

    # Zero the entire array first (important if this is called in a loop).
    bitset[:] = 0

    # Now set each bit.  Since we only know cat_offs[i] can be anywhere in [0..255],
    # we figure out which uint32_t word to hit and which bit within it:
    for i in range(ncats_present):
        if template & ((<BITSET_INNER_DTYPE_C>1) << i):
            # target_bit = cat_offs[i];   # e.g. if cat_offs[i]==70, that's bit #70 of the 256-bit mask
            target_bit = cat_offs[i]
            word_idx   = target_bit // 32  # which of the 8 uint32_t words
            bit_in_word = target_bit % 32  # which bit within that 32-bit word
            # set that single bit
            bitset[word_idx] |= (<BITSET_INNER_DTYPE_C>1 << bit_in_word)