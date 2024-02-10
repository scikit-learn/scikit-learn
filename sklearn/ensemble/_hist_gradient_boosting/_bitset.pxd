from .common cimport X_BINNED_DTYPE_FUSED_C
from .common cimport BITSET_INNER_DTYPE_C
from .common cimport X_DTYPE_C
from .common cimport Bitsets

cdef void set_bitset(
    BITSET_INNER_DTYPE_C* bitset, X_BINNED_DTYPE_FUSED_C val
) noexcept nogil

cdef unsigned char in_bitset(
    const BITSET_INNER_DTYPE_C* bitset, X_BINNED_DTYPE_FUSED_C val
) noexcept nogil

cdef create_feature_bitset_array(Bitsets b, int bitset_idx)
