# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.
cimport numpy as cnp
from sklearn.neighbors._quad_tree cimport Cell
from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t, uint64_t


ctypedef uint64_t bitword_t

cdef enum:
    BITWORD_SHIFT = 6          # log2(64)
    BITWORD_MASK  = 63         # 64 - 1
    BITWORD_BITS  = 64

cdef enum:
    MAX_CAT_BITSET_WORDS = 4   # ceil(256 / 64)

ctypedef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.split_value.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the BITSET_INNER_DTYPE_C view (`SplitValue.categorical_split``) is
    # used. It works in one of two ways, indicated by the value of its least
    # significant bit (LSB). If the LSB is 0, then categorical_split acts as a bitfield
    # for up to 64 categories, sending samples left if the bit corresponding to
    # their category is 1 or right if it is 0. If the LSB is 1, then the most
    # significant 32 bits of categorical_split make a random seed. To evaluate a
    # sample, use the random seed to flip a coin (category_value + 1) times and
    # send it left if the last flip gives 1; otherwise right. This second
    # method allows up to 2**31 category values, but can only be used for
    # RandomSplitter.
    float64_t threshold
    # BITSET_DTYPE_C categorical_split
    # Array size = ceil(MAX_NUM_CATEGORIES / 64).
    # Currently MAX_NUM_CATEGORIES = 256, so 256/64 = 4 words.
    # If you change MAX_NUM_CATEGORIES, update this array size accordingly.
    uint64_t[MAX_CAT_BITSET_WORDS] categorical_bitset

cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    SplitValue split_value             # Generalized threshold for categorical and
                                       # non-categorical features
    # float64_t threshold                  # Threshold value at the node

    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left     # Whether features have missing values


cdef intp_t n_words_for_nbits(intp_t n_bits) noexcept

cdef void set_bit_fast(bitword_t* words, intp_t c) noexcept nogil
cdef bint in_bitset_words_fast(const bitword_t* words, intp_t c) noexcept nogil

cdef bint goes_left(
    SplitValue split_value,
    bint missing_go_to_left,
    bint is_categorical,
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


