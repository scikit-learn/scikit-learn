"""
The definition of the Fibonacci heap data structure, which provide a fast
find-minimum operation needed for a number of algorithms such as Dijkstra's
algorithm for shortest graph path searches.
"""

cimport numpy as np
ctypedef np.float64_t DTYPE_t


######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.
#

cdef struct FibonacciNode:
    unsigned int index, rank, state
    DTYPE_t val
    FibonacciNode *parent
    FibonacciNode *left_sibling
    FibonacciNode *right_sibling
    FibonacciNode *children

ctypedef FibonacciNode* pFibonacciNode

cdef void initialize_node(FibonacciNode* node,
                          unsigned int index,
                          DTYPE_t val=*)

cdef FibonacciNode* rightmost_sibling(FibonacciNode* node)

cdef FibonacciNode* leftmost_sibling(FibonacciNode* node)

cdef void add_child(FibonacciNode* node, FibonacciNode* new_child)

cdef void add_sibling(FibonacciNode* node, FibonacciNode* new_sibling)

cdef void remove(FibonacciNode* node)


######################################################################
# FibonacciHeap structure
#  This structure and operations on it use the FibonacciNode
#  routines to implement a Fibonacci heap

cdef struct FibonacciHeap:
    FibonacciNode* min_node
    pFibonacciNode[100] roots_by_rank  # maximum number of nodes is ~2^100.

cdef void insert_node(FibonacciHeap* heap,
                      FibonacciNode* node)

cdef void decrease_val(FibonacciHeap* heap,
                       FibonacciNode* node,
                       DTYPE_t newval)

cdef void link(FibonacciHeap* heap, FibonacciNode* node)

cdef FibonacciNode* remove_min(FibonacciHeap* heap)
