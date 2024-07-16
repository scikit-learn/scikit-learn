# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _partitioner.pyx for details.
from ..utils._typedefs cimport BITSET_t, float32_t, float64_t, intp_t, int8_t, int32_t, uint32_t


# Mitigate precision differences between 32 bit and 64 bit
cdef float32_t FEATURE_THRESHOLD = 1e-7


ctypedef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.split_value.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the BITSET_INNER_DTYPE_C view (`SplitValue.cat_split``) is
    # used. It works in one of two ways, indicated by the value of its least
    # significant bit (LSB). If the LSB is 0, then cat_split acts as a bitfield
    # for up to 64 categories, sending samples left if the bit corresponding to
    # their category is 1 or right if it is 0. If the LSB is 1, then the most
    # significant 32 bits of cat_split make a random seed. To evaluate a
    # sample, use the random seed to flip a coin (category_value + 1) times and
    # send it left if the last flip gives 1; otherwise right. This second
    # method allows up to 2**31 category values, but can only be used for
    # RandomSplitter.
    float64_t threshold
    BITSET_t cat_split


cdef class BasePartitioner:
    cdef intp_t[::1] samples
    cdef float32_t[::1] feature_values
    cdef intp_t start
    cdef intp_t end
    cdef intp_t n_missing
    cdef const unsigned char[::1] missing_values_in_feature_mask
    cdef const int32_t[::1] n_categories
    cdef BITSET_t[::1] cat_cache

    cdef void sort_samples_and_feature_values(
        self, intp_t current_feature
    ) noexcept nogil
    cdef void init_node_split(
        self,
        intp_t start,
        intp_t end
    ) noexcept nogil
    cdef void find_min_max(
        self,
        intp_t current_feature,
        float32_t* min_feature_value_out,
        float32_t* max_feature_value_out,
    ) noexcept nogil
    cdef void next_p(
        self,
        intp_t* p_prev,
        intp_t* p
    ) noexcept nogil
    cdef intp_t partition_samples(
        self,
        SplitValue split_value,
        intp_t feature,
    ) noexcept nogil
    cdef void partition_samples_final(
        self,
        intp_t best_pos,
        SplitValue best_split_value,
        intp_t best_feature,
        intp_t n_missing,
    ) noexcept nogil
    cdef intp_t partition_samples_category(
        self,
        BITSET_t cat_split
    ) noexcept nogil


cdef class DensePartitioner(BasePartitioner):
    """Partitioner specialized for dense data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    cdef const float32_t[:, :] X

    cdef float32_t[:] sort_value
    cdef float32_t[:] sort_density
    cdef int32_t[:] cat_offset
    cdef intp_t[:] sorted_cat
    cdef bint breiman_shortcut

    cdef void _breiman_sort_categories(
        self,
        intp_t start,
        intp_t end,
        int32_t ncat,
        intp_t ncat_present,
        const int32_t[:] cat_offset,
        intp_t[:] sorted_cat,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
    ) noexcept nogil

cdef class SparsePartitioner(BasePartitioner):
    """Partitioner specialized for sparse CSC data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    cdef const float32_t[::1] X_data
    cdef const int32_t[::1] X_indices
    cdef const int32_t[::1] X_indptr

    cdef intp_t n_total_samples

    cdef intp_t[::1] index_to_samples
    cdef intp_t[::1] sorted_samples

    cdef intp_t start_positive
    cdef intp_t end_negative
    cdef bint is_samples_sorted

    cdef float32_t[:] sort_value
    cdef float32_t[:] sort_density
    cdef int32_t[:] cat_offset
    cdef intp_t[:] sorted_cat
    cdef bint breiman_shortcut

    cdef void extract_nnz(
        self,
        intp_t feature
    ) noexcept nogil
    cdef intp_t _partition(
        self,
        float64_t threshold,
        intp_t zero_pos
    ) noexcept nogil

    cdef void _breiman_sort_categories(
        self,
        intp_t start,
        intp_t end,
        int32_t ncat,
        intp_t ncat_present,
        const int32_t[:] cat_offset,
        intp_t[:] sorted_cat,
        const float64_t[:, ::1] y,
        const float64_t[:] sample_weight,
    ) noexcept nogil
