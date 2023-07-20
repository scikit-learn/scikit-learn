# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

# See _utils.pyx for details.
import cython
from libcpp.vector cimport vector

cimport numpy as cnp

from sklearn.neighbors._quad_tree cimport Cell


ctypedef cnp.npy_float32 DTYPE_t          # Type of X
ctypedef cnp.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef cnp.npy_intp SIZE_t              # Type for indices and counters
ctypedef cnp.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef cnp.npy_uint32 UINT32_t          # Unsigned 32 bit integer
ctypedef cnp.npy_uint64 UINT64_t          # Unsigned 64 bit integer

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    # It corresponds to the maximum representable value for
    # 32-bit signed integers (i.e. 2^31 - 1).
    RAND_R_MAX = 2147483647


cdef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the UINT64_t view (`SplitValue.cat_split``) is
    # used. It works in one of two ways, indicated by the value of its least
    # significant bit (LSB). If the LSB is 0, then cat_split acts as a bitfield
    # for up to 64 categories, sending samples left if the bit corresponding to
    # their category is 1 or right if it is 0. If the LSB is 1, then the most
    # significant 32 bits of cat_split make a random seed. To evaluate a
    # sample, use the random seed to flip a coin (category_value + 1) times and
    # send it left if the last flip gives 1; otherwise right. This second
    # method allows up to 2**31 category values, but can only be used for
    # RandomSplitter.
    DOUBLE_t threshold
    UINT64_t cat_split

cdef struct SplitRecord:
    # Data to track sample split
    SIZE_t feature         # Which feature to split on.
    SIZE_t pos             # Split samples array at the given position,
    #                      # i.e. count of samples below threshold for feature.
    #                      # pos is >= end if the node is a leaf.
    # SplitValue split_value # Generalized threshold for categorical and
    #                      # non-categorical features
    DOUBLE_t threshold
    UINT64_t cat_split
    double improvement     # Impurity improvement given parent node.
    double impurity_left   # Impurity of the left split.
    double impurity_right  # Impurity of the right split.
    double lower_bound     # Lower bound on value of both children for monotonicity
    double upper_bound     # Upper bound on value of both children for monotonicity
    unsigned char missing_go_to_left  # Controls if missing values go to the left node.
    SIZE_t n_missing        # Number of missing values for the feature being split on

cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    SIZE_t left_child                    # id of the left child of the node
    SIZE_t right_child                   # id of the right child of the node
    SIZE_t feature                       # Feature used for splitting the node
    # SplitValue split_value             # Generalized threshold for categorical and
    #                                    # non-categorical features
    DOUBLE_t threshold
    UINT64_t cat_split
    DOUBLE_t impurity                    # Impurity of the node (i.e., the value of the criterion)
    SIZE_t n_node_samples                # Number of samples at the node
    DOUBLE_t weighted_n_node_samples     # Weighted number of samples at the node
    unsigned char missing_go_to_left     # Whether features have missing values


# safe_realloc(&p, n) resizes the allocation of p to n * sizeof(*p) bytes or
# raises a MemoryError. It never calls free, since that's __dealloc__'s job.
#   cdef DTYPE_t *p = NULL
#   safe_realloc(&p, n)
# is equivalent to p = malloc(n * sizeof(*p)) with error checking.
ctypedef fused realloc_ptr:
    # Add pointer types here as needed.
    (DTYPE_t*)
    (SIZE_t*)
    (unsigned char*)
    (WeightedPQueueRecord*)
    (DOUBLE_t*)
    (DOUBLE_t**)
    (Node*)
    (Cell*)
    (Node**)
    (void**)
    (INT32_t*)
    (UINT32_t*)
    (UINT64_t*)

cdef realloc_ptr safe_realloc(realloc_ptr* p, size_t nelems) except * nogil


cdef cnp.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size)


cdef cnp.ndarray int32_ptr_to_ndarray(INT32_t* data, SIZE_t size)


cdef SIZE_t rand_int(SIZE_t low, SIZE_t high,
                     UINT32_t* random_state) noexcept nogil


cdef double rand_uniform(double low, double high,
                         UINT32_t* random_state) noexcept nogil


cdef double log(double x) noexcept nogil


cdef void setup_cat_cache(
    vector[UINT64_t]& cachebits,
    UINT64_t cat_split,
    INT32_t n_categories
) noexcept nogil


cdef bint goes_left(
    DTYPE_t feature_value,
    # SplitValue split,
    # DOUBLE_t threshold,
    # INT32_t n_categories,
    Node* node,
    const INT32_t[:] n_categories,
    vector[UINT64_t]& cachebits
) noexcept nogil

# =============================================================================
# WeightedPQueue data structure
# =============================================================================

# A record stored in the WeightedPQueue
cdef struct WeightedPQueueRecord:
    DOUBLE_t data
    DOUBLE_t weight

cdef class WeightedPQueue:
    cdef SIZE_t capacity
    cdef SIZE_t array_ptr
    cdef WeightedPQueueRecord* array_

    cdef bint is_empty(self) noexcept nogil
    cdef int reset(self) except -1 nogil
    cdef SIZE_t size(self) noexcept nogil
    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) except -1 nogil
    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) noexcept nogil
    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) noexcept nogil
    cdef int peek(self, DOUBLE_t* data, DOUBLE_t* weight) noexcept nogil
    cdef DOUBLE_t get_weight_from_index(self, SIZE_t index) noexcept nogil
    cdef DOUBLE_t get_value_from_index(self, SIZE_t index) noexcept nogil


# =============================================================================
# WeightedMedianCalculator data structure
# =============================================================================

cdef class WeightedMedianCalculator:
    cdef SIZE_t initial_capacity
    cdef WeightedPQueue samples
    cdef DOUBLE_t total_weight
    cdef SIZE_t k
    cdef DOUBLE_t sum_w_0_k  # represents sum(weights[0:k]) = w[0] + w[1] + ... + w[k-1]
    cdef SIZE_t size(self) noexcept nogil
    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) except -1 nogil
    cdef int reset(self) except -1 nogil
    cdef int update_median_parameters_post_push(
        self, DOUBLE_t data, DOUBLE_t weight,
        DOUBLE_t original_median) noexcept nogil
    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) noexcept nogil
    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) noexcept nogil
    cdef int update_median_parameters_post_remove(
        self, DOUBLE_t data, DOUBLE_t weight,
        DOUBLE_t original_median) noexcept nogil
    cdef DOUBLE_t get_median(self) noexcept nogil


cdef UINT64_t bs_set(UINT64_t value, SIZE_t i) noexcept nogil
cdef UINT64_t bs_reset(UINT64_t value, SIZE_t i) noexcept nogil
cdef UINT64_t bs_flip(UINT64_t value, SIZE_t i) noexcept nogil
cdef UINT64_t bs_flip_all(UINT64_t value, SIZE_t n_low_bits) noexcept nogil
cdef bint bs_get(UINT64_t value, SIZE_t i) noexcept nogil
cdef UINT64_t bs_from_template(UINT64_t template,
                               INT32_t *cat_offs,
                               SIZE_t ncats_present) noexcept nogil
