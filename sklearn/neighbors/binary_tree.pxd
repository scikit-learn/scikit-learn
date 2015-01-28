#!python

cimport cython
cimport numpy as np

from dist_metrics cimport DistanceMetric

from typedefs cimport DTYPE_t, ITYPE_t


cdef struct NodeData_t:
    ITYPE_t idx_start, idx_end, is_leaf
    DTYPE_t radius


######################################################################
# Binary Tree class
cdef class BinaryTree:

    cdef np.ndarray data_arr
    cdef np.ndarray idx_array_arr
    cdef np.ndarray node_data_arr
    cdef np.ndarray node_bounds_arr

    cdef readonly DTYPE_t[:, ::1] data
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
                             ITYPE_t size) except -1

    cdef inline DTYPE_t rdist(self, DTYPE_t* x1, DTYPE_t* x2,
                              ITYPE_t size) except -1

    cdef int _recursive_build(self, ITYPE_t i_node, ITYPE_t idx_start,
                              ITYPE_t idx_end) except -1

    cdef int _query_single_depthfirst(self, ITYPE_t i_node,
                                      DTYPE_t* pt, ITYPE_t i_pt,
                                      NeighborsHeap heap,
                                      DTYPE_t reduced_dist_LB) except -1

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
                                      int return_distance) except -1

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
                   DTYPE_t* global_log_bound_spread) except -1

    cdef int _two_point_single(self, ITYPE_t i_node, DTYPE_t* pt, DTYPE_t* r,
                               ITYPE_t* count, ITYPE_t i_min,
                               ITYPE_t i_max) except -1

    cdef int _two_point_dual(self, ITYPE_t i_node1,
                             BinaryTree other, ITYPE_t i_node2,
                             DTYPE_t* r, ITYPE_t* count,
                             ITYPE_t i_min, ITYPE_t i_max) except -1

    cdef int allocate_date(self, ITYPE_t n_nodes, ITYPE_t n_features) except -1

    cdef int init_node(self, ITYPE_t i_node,
                       ITYPE_t idx_start, ITYPE_t idx_end) except -1

    cdef DTYPE_t min_rdist(BinaryTree tree, ITYPE_t i_node,
                           DTYPE_t* pt) except -1

    cdef DTYPE_t min_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt) except -1

    cdef DTYPE_t max_rdist(BinaryTree tree, ITYPE_t i_node,
                           DTYPE_t* pt) except -1

    cdef DTYPE_t max_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt) except -1

    cdef int min_max_dist(BinaryTree tree, ITYPE_t i_node,
                          DTYPE_t* pt, DTYPE_t* min_dist,
                          DTYPE_t* max_dist) except -1

    cdef DTYPE_t min_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                BinaryTree tree2, ITYPE_t i_node2) except -1

    cdef DTYPE_t min_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                               BinaryTree tree2, ITYPE_t i_node2) except -1

    cdef DTYPE_t max_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                BinaryTree tree2, ITYPE_t i_node2) except -1

    cdef DTYPE_t max_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                               BinaryTree tree2, ITYPE_t i_node2) except -1
