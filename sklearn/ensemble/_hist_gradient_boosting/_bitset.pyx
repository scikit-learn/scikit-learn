from .common cimport Bitsets
from .common cimport BITSET_INNER_DTYPE_C
from .common cimport X_DTYPE_C
from .common cimport X_BINNED_DTYPE_C
from ...utils._typedefs cimport uint8_t

import numpy as np


# A bitset is a data structure used to represent sets of integers in [0, n]. We
# use them to represent sets of features indices (e.g. features that go to the
# left child, or features that are categorical). For familiarity with bitsets
# and bitwise operations:
# https://en.wikipedia.org/wiki/Bit_array
# https://en.wikipedia.org/wiki/Bitwise_operation


cdef inline void set_bitset(
    BITSET_INNER_DTYPE_C* bitset, X_BINNED_DTYPE_C val
) noexcept nogil:
    bitset[val // 32] |= (1 << (val % 32))


cdef inline unsigned char in_bitset(
    const BITSET_INNER_DTYPE_C* bitset, X_BINNED_DTYPE_C val
) noexcept nogil:
    return (bitset[val // 32] >> (val % 32)) & 1


def set_bitset_memoryview(
    BITSET_INNER_DTYPE_C[:] bitset, X_BINNED_DTYPE_C val
):
    """For testing only."""
    set_bitset(&bitset[0], val)


def in_bitset_memoryview(
    const BITSET_INNER_DTYPE_C[:] bitset, X_BINNED_DTYPE_C val
):
    """For testing only."""
    return in_bitset(&bitset[0], val)


cdef create_feature_bitset_array(Bitsets b, int bitset_idx):
    return np.copy(
        b.bitsets_view[
            b.offsets_view[bitset_idx] : b.offsets_view[bitset_idx + 1]
        ]
    )


def copyto_feature_bitset_array(
    Bitsets dst, int bitset_idx, BITSET_INNER_DTYPE_C[::1] src
):
    np.copyto(
        dst=dst.bitsets[
            dst.offsets_view[bitset_idx] : dst.offsets_view[bitset_idx + 1]
        ],
        src=src,
    )


def set_raw_bitset_from_binned_bitset(
    Bitsets raw_bitsets,  # OUT
    BITSET_INNER_DTYPE_C[:] binned_bitset,
    X_DTYPE_C[:] categories,
    unsigned int bitset_idx,
):
    """Set the raw_bitset from the values of the binned bitset

    categories is a mapping from binned category value to raw category value.
    """
    cdef:
        int binned_cat_value
        X_DTYPE_C raw_cat_value
        BITSET_INNER_DTYPE_C* raw_bitset = raw_bitsets.at(bitset_idx)

    for binned_cat_value, raw_cat_value in enumerate(categories):
        if in_bitset(&binned_bitset[0], binned_cat_value):
            set_bitset(raw_bitset, <X_BINNED_DTYPE_C>raw_cat_value)


def set_known_cat_bitset_from_known_categories(
    Bitsets known_cat_bitsets,  # OUT
    object known_categories,
    uint8_t[:] is_categorical,
):
    """Set the bitsets for known categories.

    Parameters
    ----------
    known_cat_bitsets : Bitsets
        The bitsets to set/fill.
    known_categories : list of {ndarray, None} of shape (n_features,) \
            dtype=X_DTYPE, default=none
        For each categorical feature, the array indicates the set of unique
        categorical values. These should be the possible values over all the
        data, not just the training data. For continuous features, the
        corresponding entry should be None.
    is_categorical : ndarray of shape (n_features,), dtype=np.uint8
        Indicator for categorical features.
    """
    cdef:
        unsigned int f_idx
        unsigned int n_features = is_categorical.shape[0]
        X_DTYPE_C raw_cat_value
        BITSET_INNER_DTYPE_C* raw_bitset

    # TODO: Think twice if it's worth to parallelize, which would require to wrap
    # know_categories in a single array and a 2nd offset array.
    for f_idx in range(n_features):
        if is_categorical[f_idx]:
            raw_bitset = known_cat_bitsets.at(f_idx)
            for raw_cat_value in known_categories[f_idx]:
                set_bitset(raw_bitset, <X_BINNED_DTYPE_C>raw_cat_value)
