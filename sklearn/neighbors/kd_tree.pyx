#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

# By Jake Vanderplas (2013) <jakevdp@cs.washington.edu>
# written for the scikit-learn project
# License: BSD

__all__ = ['KDTree']

DOC_DICT = {'BinaryTree': 'KDTree', 'binary_tree': 'kd_tree'}

VALID_METRICS = ['EuclideanDistance', 'ManhattanDistance',
                 'ChebyshevDistance', 'MinkowskiDistance']


#----------------------------------------------------------------------
# Here's our big hack: we can't subclass BinaryTree, because polymorphism
# doesn't work in cython.  The dual-tree queries defined in BinaryTree
# break if we try this approach.  So we use a literal include to "inherit"
# all the boiler-plate code, and assign BinaryTree to KDTree.  The
# specifics of the implementation are the functions in this module.
include "binary_tree.pxi"
KDTree = BinaryTree

#----------------------------------------------------------------------
# The functions below specialized the Binary Tree as a KD Tree
#
#   Note that these functions use the concept of "reduced distance".
#   The reduced distance, defined for some metrics, is a quantity which
#   is more efficient to compute than the distance, but preserves the
#   relative rankings of the true distance.  For example, the reduced
#   distance for the Euclidean metric is the squared-euclidean distance.
#   For some metrics, the reduced distance is simply the distance.


cdef int allocate_data(BinaryTree bt, ITYPE_t n_nodes,
                       ITYPE_t n_features) except -1:
    """Allocate arrays needed for the KD Tree"""
    bt.node_bounds = np.zeros((2, n_nodes, n_features), dtype=DTYPE)
    return 0


cdef int init_node(BinaryTree bt, ITYPE_t i_node,
                   ITYPE_t idx_start, ITYPE_t idx_end) except -1:
    """Initialize the node for the dataset stored in bt.data"""
    cdef ITYPE_t n_features = bt.data.shape[1]
    cdef ITYPE_t i, j
    cdef DTYPE_t rad = 0

    cdef DTYPE_t* lower_bounds = &bt.node_bounds[0, i_node, 0]
    cdef DTYPE_t* upper_bounds = &bt.node_bounds[1, i_node, 0]
    cdef DTYPE_t* data = &bt.data[0, 0]
    cdef ITYPE_t* idx_array = &bt.idx_array[0]

    cdef DTYPE_t* data_row

    # determine Node bounds
    for j in range(n_features):
        lower_bounds[j] = INF
        upper_bounds[j] = -INF

    # Compute the actual data range.  At build time, this is slightly
    # slower than using the previously-computed bounds of the parent node,
    # but leads to more compact trees and thus faster queries.
    for i in range(idx_start, idx_end):
        data_row = data + idx_array[i] * n_features
        for j in range(n_features):
            lower_bounds[j] = fmin(lower_bounds[j], data_row[j])
            upper_bounds[j] = fmax(upper_bounds[j], data_row[j])
        if bt.dm.p == INF:
            rad = fmax(rad, 0.5 * (upper_bounds[j] - lower_bounds[j]))
        else:
            rad += pow(0.5 * abs(upper_bounds[j] - lower_bounds[j]), bt.dm.p)

    bt.node_data[i_node].idx_start = idx_start
    bt.node_data[i_node].idx_end = idx_end

    # The radius will hold the size of the circumscribed hypersphere measured
    # with the specified metric: in querying, this is used as a measure of the
    # size of each node when deciding which nodes to split.
    bt.node_data[i_node].radius = pow(rad, 1. / bt.dm.p)
    return 0


cdef DTYPE_t min_rdist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt) except -1:
    """Compute the minimum reduced-distance between a point and a node"""
    cdef ITYPE_t n_features = bt.data.shape[1]
    cdef DTYPE_t d, d_lo, d_hi, rdist=0.0
    cdef ITYPE_t j

    if bt.dm.p == INF:
        for j in range(n_features):
            d_lo = bt.node_bounds[0, i_node, j] - pt[j]
            d_hi = pt[j] - bt.node_bounds[1, i_node, j]
            d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
            rdist = fmax(rdist, 0.5 * d)
    else:
        # here we'll use the fact that x + abs(x) = 2 * max(x, 0)
        for j in range(n_features):
            d_lo = bt.node_bounds[0, i_node, j] - pt[j]
            d_hi = pt[j] - bt.node_bounds[1, i_node, j]
            d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
            rdist += pow(0.5 * d, bt.dm.p)

    return rdist


cdef DTYPE_t min_dist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt) except -1:
    """Compute the minimum distance between a point and a node"""
    if bt.dm.p == INF:
        return min_rdist(bt, i_node, pt)
    else:
        return pow(min_rdist(bt, i_node, pt), 1. / bt.dm.p)


cdef DTYPE_t max_rdist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt) except -1:
    """Compute the maximum reduced-distance between a point and a node"""
    cdef ITYPE_t n_features = bt.data.shape[1]

    cdef DTYPE_t d, d_lo, d_hi, rdist=0.0
    cdef ITYPE_t j

    if bt.dm.p == INF:
        for j in range(n_features):
            rdist = fmax(rdist, fabs(pt[j] - bt.node_bounds[0, i_node, j]))
            rdist = fmax(rdist, fabs(pt[j] - bt.node_bounds[1, i_node, j]))
    else:
        for j in range(n_features):
            d_lo = fabs(pt[j] - bt.node_bounds[0, i_node, j])
            d_hi = fabs(pt[j] - bt.node_bounds[1, i_node, j])
            rdist += pow(fmax(d_lo, d_hi), bt.dm.p)

    return rdist


cdef DTYPE_t max_dist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt) except -1:
    """Compute the maximum distance between a point and a node"""
    if bt.dm.p == INF:
        return max_rdist(bt, i_node, pt)
    else:
        return pow(max_rdist(bt, i_node, pt), 1. / bt.dm.p)


cdef inline int min_max_dist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt,
                             DTYPE_t* min_dist, DTYPE_t* max_dist) except -1:
    """Compute the minimum and maximum distance between a point and a node"""
    cdef ITYPE_t n_features = bt.data.shape[1]

    cdef DTYPE_t d, d_lo, d_hi
    cdef ITYPE_t j

    min_dist[0] = 0.0
    max_dist[0] = 0.0

    if bt.dm.p == INF:
        for j in range(n_features):
            d_lo = bt.node_bounds[0, i_node, j] - pt[j]
            d_hi = pt[j] - bt.node_bounds[1, i_node, j]
            d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
            min_dist[0] = fmax(min_dist[0], 0.5 * d)
            max_dist[0] = fmax(max_dist[0],
                               fabs(pt[j] - bt.node_bounds[0, i_node, j]))
            max_dist[0] = fmax(max_dist[0],
                               fabs(pt[j] - bt.node_bounds[1, i_node, j]))
    else:
        # as above, use the fact that x + abs(x) = 2 * max(x, 0)
        for j in range(n_features):
            d_lo = bt.node_bounds[0, i_node, j] - pt[j]
            d_hi = pt[j] - bt.node_bounds[1, i_node, j]
            d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
            min_dist[0] += pow(0.5 * d, bt.dm.p)
            max_dist[0] += pow(fmax(fabs(d_lo), fabs(d_hi)), bt.dm.p)

        min_dist[0] = pow(min_dist[0], 1. / bt.dm.p)
        max_dist[0] = pow(max_dist[0], 1. / bt.dm.p)

    return 0


cdef inline DTYPE_t min_rdist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                   BinaryTree bt2, ITYPE_t i_node2) except -1:
    """Compute the minimum reduced distance between two nodes"""
    cdef ITYPE_t n_features = bt1.data.shape[1]

    cdef DTYPE_t d, d1, d2, rdist=0.0
    cdef DTYPE_t zero = 0.0
    cdef ITYPE_t j

    if bt1.dm.p == INF:
        for j in range(n_features):
            d1 = (bt1.node_bounds[0, i_node1, j]
                  - bt2.node_bounds[1, i_node2, j])
            d2 = (bt2.node_bounds[0, i_node2, j]
                  - bt1.node_bounds[1, i_node1, j])
            d = (d1 + fabs(d1)) + (d2 + fabs(d2))

            rdist = fmax(rdist, 0.5 * d)
    else:
        # here we'll use the fact that x + abs(x) = 2 * max(x, 0)
        for j in range(n_features):
            d1 = (bt1.node_bounds[0, i_node1, j]
                  - bt2.node_bounds[1, i_node2, j])
            d2 = (bt2.node_bounds[0, i_node2, j]
                  - bt1.node_bounds[1, i_node1, j])
            d = (d1 + fabs(d1)) + (d2 + fabs(d2))

            rdist += pow(0.5 * d, bt1.dm.p)

    return rdist


cdef inline DTYPE_t min_dist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                  BinaryTree bt2, ITYPE_t i_node2) except -1:
    """Compute the minimum distance between two nodes"""
    return bt1.dm._rdist_to_dist(min_rdist_dual(bt1, i_node1,
                                                bt2, i_node2))


cdef inline DTYPE_t max_rdist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                   BinaryTree bt2, ITYPE_t i_node2) except -1:
    """Compute the maximum reduced distance between two nodes"""
    cdef ITYPE_t n_features = bt1.data.shape[1]

    cdef DTYPE_t d, d1, d2, rdist=0.0
    cdef DTYPE_t zero = 0.0
    cdef ITYPE_t j

    if bt1.dm.p == INF:
        for j in range(n_features):
            rdist = fmax(rdist, fabs(bt1.node_bounds[0, i_node1, j]
                                     - bt2.node_bounds[1, i_node2, j]))
            rdist = fmax(rdist, fabs(bt1.node_bounds[1, i_node1, j]
                                     - bt2.node_bounds[0, i_node2, j]))
    else:
        for j in range(n_features):
            d1 = fabs(bt1.node_bounds[0, i_node1, j]
                      - bt2.node_bounds[1, i_node2, j])
            d2 = fabs(bt1.node_bounds[1, i_node1, j]
                      - bt2.node_bounds[0, i_node2, j])
            rdist += pow(fmax(d1, d2), bt1.dm.p)

    return rdist


cdef inline DTYPE_t max_dist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                  BinaryTree bt2, ITYPE_t i_node2) except -1:
    """Compute the maximum distance between two nodes"""
    return bt1.dm._rdist_to_dist(max_rdist_dual(bt1, i_node1,
                                                bt2, i_node2))
