# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.
cimport numpy as cnp
from sklearn.neighbors._quad_tree cimport Cell
from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t, uint64_t
from sklearn.utils._bitset cimport BITSET_DTYPE_C

cdef enum:
    SPLIT_NUMERIC = 0
    SPLIT_CATEGORICAL_BITSET = 1
    SPLIT_CATEGORICAL_HASH = 2

ctypedef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.split_value.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the interpretation of categorical_bitset
    # depends on split_kind:
    # - SPLIT_CATEGORICAL_BITSET: stores the set of categories that go left.
    # - SPLIT_CATEGORICAL_HASH: categorical_bitset[0] stores the hash seed.
    float64_t threshold
    # Array size = ceil(MAX_NUM_CATEGORIES / 32).
    # Currently MAX_NUM_CATEGORIES = 256, so 256/32 = 8 words.
    # If you change MAX_NUM_CATEGORIES, update BITSET_DTYPE_C accordingly.
    BITSET_DTYPE_C categorical_bitset

cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    SplitValue split_value               # Generalized threshold for categorical and
    #                                    # non-categorical features

    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left     # Whether missing values go to the left child
    uint8_t split_kind             # Encoding used by split_value


cdef bint goes_left(
    SplitValue split_value,
    bint missing_go_to_left,
    uint8_t split_kind,
    float32_t value,
) noexcept nogil


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
    (uint8_t*)
    (uint32_t*)
    (uint64_t*)
    (float64_t*)
    (float64_t**)
    (Node*)
    (Cell*)
    (Node**)

cdef int safe_realloc(realloc_ptr* p, size_t nelems) except -1 nogil


cdef cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size)


cdef intp_t rand_int(intp_t low, intp_t high,
                     uint32_t* random_state) noexcept nogil


cdef float64_t rand_uniform(float64_t low, float64_t high,
                            uint32_t* random_state) noexcept nogil


cdef float64_t log(float64_t x) noexcept nogil


cdef class WeightedFenwickTree:
    cdef intp_t size         # number of leaves (ranks)
    cdef float64_t* tree_w   # BIT for weights
    cdef float64_t* tree_wy  # BIT for weighted targets
    cdef intp_t max_pow2     # highest power of two <= n
    cdef float64_t total_w   # running total weight
    cdef float64_t total_wy  # running total weighted target

    cdef void reset(self, intp_t size) noexcept nogil
    cdef void add(self, intp_t idx, float64_t y, float64_t w) noexcept nogil
    cdef intp_t search(
        self,
        float64_t t,
        float64_t* cw_out,
        float64_t* cwy_out,
        intp_t* prev_idx_out,
    ) noexcept nogil
