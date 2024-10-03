# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.

cimport numpy as cnp

from ..neighbors._quad_tree cimport Cell
from ..utils._typedefs cimport (BITSET_t, float32_t, float64_t, int32_t,
                                intp_t, uint8_t, uint32_t, uint64_t)


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


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    # SplitValue split_value             # Generalized threshold for categorical and
    #                                    # non-categorical features
    float64_t threshold
    BITSET_t cat_split
    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left     # Whether features have missing values


cdef struct SplitRecord:
    # Data to track sample split
    intp_t feature         # Which feature to split on.
    intp_t pos             # Split samples array at the given position,
    #                      # i.e. count of samples below threshold for feature.
    #                      # pos is >= end if the node is a leaf.
    SplitValue split_value    # Generalized threshold for categorical and
    #                         # non-categorical features to split samples.
    float64_t improvement     # Impurity improvement given parent node.
    float64_t impurity_left   # Impurity of the left split.
    float64_t impurity_right  # Impurity of the right split.
    float64_t lower_bound     # Lower bound on value of both children for monotonicity
    float64_t upper_bound     # Upper bound on value of both children for monotonicity
    uint8_t missing_go_to_left  # Controls if missing values go to the left node.
    intp_t n_missing            # Number of missing values for the feature being split on


cdef struct ParentInfo:
    # Structure to store information about the parent of a node
    # This is passed to the splitter, to provide information about the previous split

    float64_t lower_bound           # the lower bound of the parent's impurity
    float64_t upper_bound           # the upper bound of the parent's impurity
    float64_t impurity              # the impurity of the parent
    intp_t n_constant_features      # the number of constant features found in parent



cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    # It corresponds to the maximum representable value for
    # 32-bit signed integers (i.e. 2^31 - 1).
    RAND_R_MAX = 2147483647


# safe_realloc(&p, n) resizes the allocation of p to n * sizeof(*p) bytes or
# raises a MemoryError. It never calls free, since that's __dealloc__'s job.
#   cdef float32_t *p = NULL
#   safe_realloc(&p, n)
# is equivalent to p = malloc(n * sizeof(*p)) with error checking.
ctypedef fused realloc_ptr:
    # Add pointer types here as needed.
    (float32_t*)
    (intp_t*)
    (int32_t*)
    (uint8_t*)
    (uint32_t*)
    (uint64_t*)
    (WeightedPQueueRecord*)
    (float64_t*)
    (float64_t**)
    (Node*)
    (Cell*)
    (Node**)

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil


cdef cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size)

cdef cnp.ndarray int32t_ptr_to_ndarray(int32_t* data, intp_t size)

cdef intp_t rand_int(intp_t low, intp_t high,
                     uint32_t* random_state) noexcept nogil


cdef float64_t rand_uniform(float64_t low, float64_t high,
                            uint32_t* random_state) noexcept nogil


cdef float64_t log(float64_t x) noexcept nogil

# =============================================================================
# WeightedPQueue data structure
# =============================================================================

# A record stored in the WeightedPQueue
cdef struct WeightedPQueueRecord:
    float64_t data
    float64_t weight

cdef class WeightedPQueue:
    cdef intp_t capacity
    cdef intp_t array_ptr
    cdef WeightedPQueueRecord* array_

    cdef bint is_empty(self) noexcept nogil
    cdef int reset(self) except -1 nogil
    cdef intp_t size(self) noexcept nogil
    cdef int push(self, float64_t data, float64_t weight) except -1 nogil
    cdef int remove(self, float64_t data, float64_t weight) noexcept nogil
    cdef int pop(self, float64_t* data, float64_t* weight) noexcept nogil
    cdef int peek(self, float64_t* data, float64_t* weight) noexcept nogil
    cdef float64_t get_weight_from_index(self, intp_t index) noexcept nogil
    cdef float64_t get_value_from_index(self, intp_t index) noexcept nogil


# =============================================================================
# WeightedMedianCalculator data structure
# =============================================================================

cdef class WeightedMedianCalculator:
    cdef intp_t initial_capacity
    cdef WeightedPQueue samples
    cdef float64_t total_weight
    cdef intp_t k
    cdef float64_t sum_w_0_k  # represents sum(weights[0:k]) = w[0] + w[1] + ... + w[k-1]
    cdef intp_t size(self) noexcept nogil
    cdef int push(self, float64_t data, float64_t weight) except -1 nogil
    cdef int reset(self) except -1 nogil
    cdef int update_median_parameters_post_push(
        self, float64_t data, float64_t weight,
        float64_t original_median) noexcept nogil
    cdef int remove(self, float64_t data, float64_t weight) noexcept nogil
    cdef int pop(self, float64_t* data, float64_t* weight) noexcept nogil
    cdef int update_median_parameters_post_remove(
        self, float64_t data, float64_t weight,
        float64_t original_median) noexcept nogil
    cdef float64_t get_median(self) noexcept nogil


cdef void setup_cat_cache(
    BITSET_t[:] cachebits,
    BITSET_t cat_split,
    int32_t n_categories
) noexcept nogil

cdef BITSET_t bs_set(BITSET_t value, intp_t i) noexcept nogil
cdef BITSET_t bs_reset(BITSET_t value, intp_t i) noexcept nogil
cdef BITSET_t bs_flip(BITSET_t value, intp_t i) noexcept nogil
cdef BITSET_t bs_flip_all(BITSET_t value, intp_t n_low_bits) noexcept nogil
cdef bint bs_get(BITSET_t value, intp_t i) noexcept nogil
cdef BITSET_t bs_from_template(
    uint64_t template,
    int32_t[:] cat_offs,
    intp_t ncats_present
) noexcept nogil
