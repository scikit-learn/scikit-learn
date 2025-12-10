# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _partitioner.pyx for details.

from cython cimport floating

from sklearn.utils._typedefs cimport (
    float32_t, float64_t, int8_t, int32_t, intp_t, uint8_t, uint32_t
)
from sklearn.tree._splitter cimport SplitRecord


# Mitigate precision differences between 32 bit and 64 bit
cdef const float32_t FEATURE_THRESHOLD = 1e-7


# We provide here the abstract interface for a Partitioner that would be
# theoretically shared between the Dense and Sparse partitioners. However,
# we leave it commented out for now as it is not used in the current
# implementation due to the performance hit from vtable lookups when using
# inheritance based polymorphism. It is left here for future reference.
#
# Note: Instead, in `_splitter.pyx`, we define a fused type that can be used
# to represent both the dense and sparse partitioners.
#
# cdef class BasePartitioner:
#     cdef intp_t[::1] samples
#     cdef float32_t[::1] feature_values
#     cdef intp_t start
#     cdef intp_t end
#     cdef intp_t n_missing
#     cdef const uint8_t[::1] missing_values_in_feature_mask

#     cdef void sort_samples_and_feature_values(
#         self, intp_t current_feature
#     ) noexcept nogil
#     cdef void init_node_split(
#         self,
#         intp_t start,
#         intp_t end
#     ) noexcept nogil
#     cdef void find_min_max(
#         self,
#         intp_t current_feature,
#         float32_t* min_feature_value_out,
#         float32_t* max_feature_value_out,
#     ) noexcept nogil
#     cdef void next_p(
#         self,
#         intp_t* p_prev,
#         intp_t* p
#     ) noexcept nogil
#     cdef intp_t partition_samples(
#         self,
#         float64_t current_threshold
#     ) noexcept nogil
#     cdef void partition_samples_final(
#         self,
#         intp_t best_pos,
#         float64_t best_threshold,
#         intp_t best_feature,
#         intp_t n_missing,
#     ) noexcept nogil


cdef class DensePartitioner:
    """Partitioner specialized for dense data.

    Note that this partitioner is agnostic to the splitting strategy (best vs. random).
    """
    cdef const float32_t[:, :] X
    cdef intp_t[::1] samples
    cdef float32_t[::1] feature_values
    cdef intp_t start
    cdef intp_t end
    cdef intp_t n_missing
    cdef const uint8_t[::1] missing_values_in_feature_mask

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
        float64_t current_threshold
    ) noexcept nogil
    cdef void partition_samples_final(
        self,
        intp_t best_pos,
        float64_t best_threshold,
        intp_t best_feature,
        intp_t n_missing,
    ) noexcept nogil


cdef class SparsePartitioner:
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

    cdef intp_t[::1] samples
    cdef float32_t[::1] feature_values
    cdef intp_t start
    cdef intp_t end
    cdef intp_t n_missing
    cdef const uint8_t[::1] missing_values_in_feature_mask

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
        float64_t current_threshold
    ) noexcept nogil
    cdef void partition_samples_final(
        self,
        intp_t best_pos,
        float64_t best_threshold,
        intp_t best_feature,
        intp_t n_missing,
    ) noexcept nogil

    cdef void extract_nnz(
        self,
        intp_t feature
    ) noexcept nogil
    cdef intp_t _partition(
        self,
        float64_t threshold,
        intp_t zero_pos
    ) noexcept nogil


cdef void shift_missing_values_to_left_if_required(
    SplitRecord* best,
    intp_t[::1] samples,
    intp_t end,
) noexcept nogil


cdef void sort(floating* feature_values, intp_t* samples, intp_t n) noexcept nogil
