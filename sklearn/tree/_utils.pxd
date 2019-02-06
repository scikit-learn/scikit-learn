# cython: language_level=3

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

# See _utils.pyx for details.

import numpy as np
cimport numpy as np

from ..neighbors.quad_tree cimport Cell

ctypedef np.npy_float32 DTYPE_t          # Type of X
ctypedef np.npy_float64 DOUBLE_t         # Type of y, sample_weight
ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_int32 INT32_t            # Signed 32 bit integer
ctypedef np.npy_uint32 UINT32_t          # Unsigned 32 bit integer
ctypedef np.npy_uint64 UINT64_t          # Unsigned 64 bit integer
ctypedef UINT64_t BITSET_t

ctypedef union SplitValue:
    # Union type to generalize the concept of a threshold to categorical
    # features. The floating point view, i.e. ``SplitValue.threshold`` is used
    # for numerical features, where feature values less than or equal to the
    # threshold go left, and values greater than the threshold go right.
    #
    # For categorical features, the BITSET_t view (`SplitValue.cat_split``) is
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
    BITSET_t cat_split


ctypedef struct SplitRecord:
    # Data to track sample split
    SIZE_t feature         # Which feature to split on.
    SIZE_t pos             # Split samples array at the given position,
                           # i.e. count of samples below threshold for feature.
                           # pos is >= end if the node is a leaf.
    SplitValue split_value # Generalized threshold for categorical and
                           # non-categorical features
    double improvement     # Impurity improvement given parent node.
    double impurity_left   # Impurity of the left split.
    double impurity_right  # Impurity of the right split.


cdef struct Node:
    # Base storage structure for the nodes in a Tree object

    SIZE_t left_child                    # id of the left child of the node
    SIZE_t right_child                   # id of the right child of the node
    SIZE_t feature                       # Feature used for splitting the node
    SplitValue split_value               # Generalized threshold for categorical and
                                         # non-categorical features
    DOUBLE_t impurity                    # Impurity of the node (i.e., the value of the criterion)
    SIZE_t n_node_samples                # Number of samples at the node
    DOUBLE_t weighted_n_node_samples     # Weighted number of samples at the node

# cdef struct Node  # Forward declaration

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF


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
    (Node*)
    (Cell*)
    (Node**)
    (StackRecord*)
    (PriorityHeapRecord*)
    (void**)
    (INT32_t*)
    (UINT32_t*)
    (BITSET_t*)

cdef realloc_ptr safe_realloc(realloc_ptr* p, size_t nelems, size_t elem_bytes) nogil except *


cdef np.ndarray sizet_ptr_to_ndarray(SIZE_t* data, SIZE_t size)

cdef np.ndarray int32_ptr_to_ndarray(INT32_t* data, SIZE_t size)

cdef SIZE_t rand_int(SIZE_t low, SIZE_t high,
                     UINT32_t* random_state) nogil


cdef double rand_uniform(double low, double high,
                         UINT32_t* random_state) nogil


cdef double log(double x) nogil


cdef void setup_cat_cache(BITSET_t* cachebits, UINT64_t cat_split,
                          INT32_t n_categories) nogil


cdef bint goes_left(DTYPE_t feature_value, SplitValue split,
                    INT32_t n_categories, BITSET_t* cachebits) nogil


# =============================================================================
# Stack data structure
# =============================================================================

# A record on the stack for depth-first tree growing
cdef struct StackRecord:
    SIZE_t start
    SIZE_t end
    SIZE_t depth
    SIZE_t parent
    bint is_left
    double impurity
    SIZE_t n_constant_features

cdef class Stack:
    cdef SIZE_t capacity
    cdef SIZE_t top
    cdef StackRecord* stack_

    cdef bint is_empty(self) nogil
    cdef int push(self, SIZE_t start, SIZE_t end, SIZE_t depth, SIZE_t parent,
                  bint is_left, double impurity,
                  SIZE_t n_constant_features) nogil except -1
    cdef int pop(self, StackRecord* res) nogil


# =============================================================================
# PriorityHeap data structure
# =============================================================================

# A record on the frontier for best-first tree growing
cdef struct PriorityHeapRecord:
    SIZE_t node_id
    SIZE_t start
    SIZE_t end
    SIZE_t pos
    SIZE_t depth
    bint is_leaf
    double impurity
    double impurity_left
    double impurity_right
    double improvement

cdef class PriorityHeap:
    cdef SIZE_t capacity
    cdef SIZE_t heap_ptr
    cdef PriorityHeapRecord* heap_

    cdef bint is_empty(self) nogil
    cdef void heapify_up(self, PriorityHeapRecord* heap, SIZE_t pos) nogil
    cdef void heapify_down(self, PriorityHeapRecord* heap, SIZE_t pos, SIZE_t heap_length) nogil
    cdef int push(self, SIZE_t node_id, SIZE_t start, SIZE_t end, SIZE_t pos,
                  SIZE_t depth, bint is_leaf, double improvement,
                  double impurity, double impurity_left,
                  double impurity_right) nogil except -1
    cdef int pop(self, PriorityHeapRecord* res) nogil

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

    cdef bint is_empty(self) nogil
    cdef int reset(self) nogil except -1
    cdef SIZE_t size(self) nogil
    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) nogil except -1
    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) nogil
    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) nogil
    cdef int peek(self, DOUBLE_t* data, DOUBLE_t* weight) nogil
    cdef DOUBLE_t get_weight_from_index(self, SIZE_t index) nogil
    cdef DOUBLE_t get_value_from_index(self, SIZE_t index) nogil


# =============================================================================
# WeightedMedianCalculator data structure
# =============================================================================

cdef class WeightedMedianCalculator:
    cdef SIZE_t initial_capacity
    cdef WeightedPQueue samples
    cdef DOUBLE_t total_weight
    cdef SIZE_t k
    cdef DOUBLE_t sum_w_0_k            # represents sum(weights[0:k])
                                       # = w[0] + w[1] + ... + w[k-1]

    cdef SIZE_t size(self) nogil
    cdef int push(self, DOUBLE_t data, DOUBLE_t weight) nogil except -1
    cdef int reset(self) nogil except -1
    cdef int update_median_parameters_post_push(
        self, DOUBLE_t data, DOUBLE_t weight,
        DOUBLE_t original_median) nogil
    cdef int remove(self, DOUBLE_t data, DOUBLE_t weight) nogil
    cdef int pop(self, DOUBLE_t* data, DOUBLE_t* weight) nogil
    cdef int update_median_parameters_post_remove(
        self, DOUBLE_t data, DOUBLE_t weight,
        DOUBLE_t original_median) nogil
    cdef DOUBLE_t get_median(self) nogil


cdef BITSET_t bs_set(BITSET_t value, SIZE_t i) nogil
cdef BITSET_t bs_reset(BITSET_t value, SIZE_t i) nogil
cdef BITSET_t bs_flip(BITSET_t value, SIZE_t i) nogil
cdef BITSET_t bs_flip_all(BITSET_t value, SIZE_t n_low_bits) nogil
cdef bint bs_get(BITSET_t value, SIZE_t i) nogil
cdef BITSET_t bs_from_template(UINT64_t template,
                               INT32_t *cat_offs,
                               SIZE_t ncats_present) nogil
