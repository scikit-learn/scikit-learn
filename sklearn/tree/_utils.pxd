# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# See _utils.pyx for details.
from libc.math cimport isnan
from libc.math cimport log as ln
from libc.stdlib cimport realloc

cimport numpy as cnp
from sklearn.neighbors._quad_tree cimport Cell
from sklearn.utils._typedefs cimport float32_t, float64_t, intp_t, uint8_t, int32_t, uint32_t, uint64_t
from sklearn.utils._bitset cimport BITSET_DTYPE_C, BITSET_INNER_DTYPE_C, N_BITSETS, in_bitset
from sklearn.utils._random cimport our_rand_r

cdef enum:
    MAX_NUM_CATEGORIES = N_BITSETS


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    intp_t left_child                    # id of the left child of the node
    intp_t right_child                   # id of the right child of the node
    intp_t feature                       # Feature used for splitting the node
    # Threshold for numerical features splits:
    # - feature values less than or equal to the threshold go left, and values greater than the threshold go right.
    float64_t threshold
    # Threshold for categorical features splits:
    # - left_cat_bitset stores the set of categories that go to the left child.
    BITSET_DTYPE_C left_cat_bitset

    float64_t impurity                   # Impurity of the node (i.e., the value of the criterion)
    intp_t n_node_samples                # Number of samples at the node
    float64_t weighted_n_node_samples    # Weighted number of samples at the node
    uint8_t missing_go_to_left           # Whether features have missing values


cdef inline bint goes_left(
    float64_t threshold,
    const BITSET_INNER_DTYPE_C* left_cat_bitset,
    bint missing_go_to_left,
    bint is_categorical,
    float32_t value,
) noexcept nogil:
    if isnan(value):
        return missing_go_to_left
    elif is_categorical:
        return in_bitset(left_cat_bitset, <uint8_t> value)
    else:
        return value <= threshold


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


cdef inline cnp.ndarray sizet_ptr_to_ndarray(intp_t* data, intp_t size):
    """Return copied data as 1D numpy array of intp's."""
    cdef cnp.npy_intp shape[1]
    shape[0] = <cnp.npy_intp> size
    return cnp.PyArray_SimpleNewFromData(1, shape, cnp.NPY_INTP, data).copy()


cdef inline intp_t rand_int(intp_t low, intp_t high,
                            uint32_t* random_state) noexcept nogil:
    """Generate a random integer in [low; end)."""
    return low + our_rand_r(random_state) % (high - low)


cdef inline float64_t rand_uniform(float64_t low, float64_t high,
                                   uint32_t* random_state) noexcept nogil:
    """Generate a random float64_t in [low; high)."""
    return ((high - low) * <float64_t> our_rand_r(random_state) /
            <float64_t> RAND_R_MAX) + low


cdef inline float64_t log(float64_t x) noexcept nogil:
    return ln(x) / ln(2.0)

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
