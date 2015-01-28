#!python

cimport cython
cimport numpy as np

from ..metrics._dist_metrics cimport DistanceMetric

from ..utils._typedefs cimport DTYPE_t, ITYPE_t

# Some compound datatypes used below:
cdef struct NodeHeapData_t:
    DTYPE_t val
    ITYPE_t i1
    ITYPE_t i2

cdef struct NodeData_t:
    ITYPE_t idx_start, idx_end, is_leaf
    DTYPE_t radius

cdef enum KernelType:
    GAUSSIAN_KERNEL = 1
    TOPHAT_KERNEL = 2
    EPANECHNIKOV_KERNEL = 3
    EXPONENTIAL_KERNEL = 4
    LINEAR_KERNEL = 5
    COSINE_KERNEL = 6

######################################################################
# NodeHeap : min-heap used to keep track of nodes during
#            breadth-first query
cdef class NodeHeap:
    """NodeHeap

    This is a min-heap implementation for keeping track of nodes
    during a breadth-first search.  Unlike the NeighborsHeap above,
    the NodeHeap does not have a fixed size and must be able to grow
    as elements are added.

    Internally, the data is stored in a simple binary heap which meets
    the min heap condition:

        heap[i].val < min(heap[2 * i + 1].val, heap[2 * i + 2].val)
    """
    cdef np.ndarray data_arr
    cdef NodeHeapData_t[::1] data
    cdef ITYPE_t n

    cdef int resize(self, ITYPE_t new_size) except -1

    cdef int push(self, NodeHeapData_t data) except -1

    cdef NodeHeapData_t peek(self)

    cdef NodeHeapData_t pop(self)

    cdef void clear(self)

cdef class NeighborsHeap:
    """A max-heap structure to keep track of distances/indices of neighbors

    This implements an efficient pre-allocated set of fixed-size heaps
    for chasing neighbors, holding both an index and a distance.
    When any row of the heap is full, adding an additional point will push
    the furthest point off the heap.

    Parameters
    ----------
    n_pts : int
        the number of heaps to use
    n_nbrs : int
        the size of each heap.
    """
    cdef np.ndarray distances_arr
    cdef np.ndarray indices_arr

    cdef DTYPE_t[:, ::1] distances
    cdef ITYPE_t[:, ::1] indices

    cdef inline DTYPE_t largest(self, ITYPE_t row) nogil except -1

    cdef int _push(self, ITYPE_t row, DTYPE_t val, ITYPE_t i_val) nogil except -1

    cdef int _sort(self) except -1


######################################################################
# Binary Tree class
cdef class BinaryTree:

    cdef np.ndarray data_arr
    cdef np.ndarray sample_weight_arr
    cdef np.ndarray idx_array_arr
    cdef np.ndarray node_data_arr
    cdef np.ndarray node_bounds_arr

    cdef readonly DTYPE_t[:, ::1] data
    cdef readonly const DTYPE_t[::1] sample_weight
    cdef public DTYPE_t sum_weight

    cdef public ITYPE_t[::1] idx_array
    cdef public NodeData_t[::1] node_data
    cdef public DTYPE_t[:, :, ::1] node_bounds

    cdef ITYPE_t leaf_size
    cdef ITYPE_t n_levels
    cdef ITYPE_t n_nodes

    cdef DistanceMetric dist_metric
    cdef int euclidean

    # variables to keep track of building & querying stats
    cdef int n_trims
    cdef int n_leaves
    cdef int n_splits
    cdef int n_calls

    cdef inline DTYPE_t dist(self, DTYPE_t* x1, DTYPE_t* x2,
                             ITYPE_t size) nogil except -1

    cdef inline DTYPE_t rdist(self, DTYPE_t* x1, DTYPE_t* x2,
                              ITYPE_t size) nogil except -1

    cdef int _recursive_build(self, NodeData_t[::1] node_data, ITYPE_t i_node, ITYPE_t idx_start,
                              ITYPE_t idx_end) except -1

    cdef int _query_single_depthfirst(self, ITYPE_t i_node,
                                      DTYPE_t* pt, ITYPE_t i_pt,
                                      NeighborsHeap heap,
                                      DTYPE_t reduced_dist_LB) nogil except -1

    cdef int _query_single_breadthfirst(self, DTYPE_t* pt,
                                        ITYPE_t i_pt,
                                        NeighborsHeap heap,
                                        NodeHeap nodeheap) except -1

    cdef int _query_dual_depthfirst(self, ITYPE_t i_node1,
                                    BinaryTree other, ITYPE_t i_node2,
                                    DTYPE_t[::1] bounds,
                                    NeighborsHeap heap,
                                    DTYPE_t reduced_dist_LB) except -1

    cdef int _query_dual_breadthfirst(self, BinaryTree other,
                                      NeighborsHeap heap,
                                      NodeHeap nodeheap) except -1

    cdef ITYPE_t _query_radius_single(self,
                                      ITYPE_t i_node,
                                      DTYPE_t* pt, DTYPE_t r,
                                      ITYPE_t* indices,
                                      DTYPE_t* distances,
                                      ITYPE_t count,
                                      int count_only,
                                      int return_distance) nogil except -1

    cdef DTYPE_t _kde_single_breadthfirst(self, DTYPE_t* pt,
                                          KernelType kernel, DTYPE_t h,
                                          DTYPE_t log_knorm,
                                          DTYPE_t log_atol, DTYPE_t log_rtol,
                                          NodeHeap nodeheap,
                                          DTYPE_t* node_log_min_bounds,
                                          DTYPE_t* node_log_bound_spreads)

    cdef int _kde_single_depthfirst(
                   self, ITYPE_t i_node, DTYPE_t* pt,
                   KernelType kernel, DTYPE_t h,
                   DTYPE_t log_knorm,
                   DTYPE_t log_atol, DTYPE_t log_rtol,
                   DTYPE_t local_log_min_bound,
                   DTYPE_t local_log_bound_spread,
                   DTYPE_t* global_log_min_bound,
                   DTYPE_t* global_log_bound_spread) nogil except -1

    cdef int _two_point_single(self, ITYPE_t i_node, DTYPE_t* pt, DTYPE_t* r,
                               ITYPE_t* count, ITYPE_t i_min,
                               ITYPE_t i_max) nogil except -1

    cdef int _two_point_dual(self, ITYPE_t i_node1,
                             BinaryTree other, ITYPE_t i_node2,
                             DTYPE_t* r, ITYPE_t* count,
                             ITYPE_t i_min, ITYPE_t i_max) nogil except -1

    cdef int allocate_date(self, ITYPE_t n_nodes, ITYPE_t n_features) except -1

    cdef int init_node(self, ITYPE_t i_node,
                       ITYPE_t idx_start, ITYPE_t idx_end) nogil except -1

    cdef DTYPE_t min_rdist(BinaryTree tree, ITYPE_t i_node,
                           DTYPE_t* pt) nogil except -1

    cdef DTYPE_t min_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt) nogil except -1

    cdef DTYPE_t max_rdist(BinaryTree tree, ITYPE_t i_node,
                           DTYPE_t* pt) nogil except -1

    cdef DTYPE_t max_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt) nogil except -1

    cdef int min_max_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt, DTYPE_t* min_dist,
                          DTYPE_t* max_dist) nogil except -1

    cdef DTYPE_t min_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                BinaryTree tree2, ITYPE_t i_node2) nogil except -1

    cdef DTYPE_t min_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                               BinaryTree tree2, ITYPE_t i_node2) nogil except -1

    cdef DTYPE_t max_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                BinaryTree tree2, ITYPE_t i_node2) nogil except -1

    cdef DTYPE_t max_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                               BinaryTree tree2, ITYPE_t i_node2) nogil except -1
