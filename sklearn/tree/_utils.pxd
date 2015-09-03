# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#
# Licence: BSD 3 clause

# See _utils.pyx for details.

import numpy as np
cimport numpy as np

ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef np.npy_float64 DOUBLE_t 

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
    DOUBLE_t impurity
    DOUBLE_t weight
    DOUBLE_t* node_value
    DOUBLE_t yw_sq_sum

cdef class Stack:
    cdef SIZE_t capacity
    cdef SIZE_t top
    cdef StackRecord* stack_

    cdef bint is_empty(self) nogil
    cdef int push(self, SIZE_t start, SIZE_t end, SIZE_t depth, SIZE_t parent,
                  bint is_left, DOUBLE_t impurity, DOUBLE_t weight, 
                  DOUBLE_t yw_sq_sum, DOUBLE_t* node_value) nogil
    cdef int pop(self, StackRecord* res) nogil


# =============================================================================
# PriorityHeap data structure
# =============================================================================

# A record on the frontier for best-first tree growing
cdef struct SplitRecord:
    # Data to track sample split
    SIZE_t parent
    SIZE_t depth
    SIZE_t is_left
    SIZE_t feature            # Which feature to split on.
    SIZE_t start
    SIZE_t end
    SIZE_t pos                # Split samples array at the given position,
                              # i.e. count of samples below threshold for feature.
                              # pos is >= end if the node is a leaf.
    DOUBLE_t threshold          # Threshold to split at.
    DOUBLE_t improvement        # Impurity improvement given parent node.
    DOUBLE_t impurity_left      # Impurity of the left split.
    DOUBLE_t impurity_right     # Impurity of the right split.
    DOUBLE_t impurity           # Impurity of the current node
    DOUBLE_t weight             # Weight of the current node
    DOUBLE_t weight_left        # Weight of the left child
    DOUBLE_t weight_right       # Weight of the right child
    DOUBLE_t yw_sq_sum
    DOUBLE_t yw_sq_sum_left
    DOUBLE_t yw_sq_sum_right
    DOUBLE_t* node_value         # Value predicted by this node
    DOUBLE_t* node_value_left    # Value predicted by the left child
    DOUBLE_t* node_value_right   # Value predicted by the right child

cdef class PriorityHeap:
    cdef SIZE_t capacity
    cdef SIZE_t heap_ptr
    cdef SplitRecord* heap_

    cdef bint is_empty(self) nogil
    cdef int push(self, SplitRecord split) nogil
    cdef int pop(self, SplitRecord* res) nogil
