#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

# By Jake Vanderplas (2013) <jakevdp@cs.washington.edu>
# written for the scikit-learn project
# License: BSD

__all__ = ['BallTree']

DOC_DICT = {'BinaryTree':'BallTree', 'binary_tree':'ball_tree'}

VALID_METRICS = ['EuclideanDistance', 'SEuclideanDistance',
                 'ManhattanDistance', 'ChebyshevDistance',
                 'MinkowskiDistance', 'WMinkowskiDistance',
                 'MahalanobisDistance', 'HammingDistance',
                 'CanberraDistance', 'BrayCurtisDistance',
                 'JaccardDistance', 'MatchingDistance',
                 'DiceDistance', 'KulsinskiDistance',
                 'RogerStanimotoDistance', 'RussellRaoDistance',
                 'SokalMichenerDistance', 'SokalSneathDistance',
                 'PyFuncDistance', 'HaversineDistance']

include "binary_tree.pxi"

BallTree = BinaryTree

#----------------------------------------------------------------------
# These functions make  the BinaryTree a BallTree
#
cdef int allocate_data(BinaryTree bt, ITYPE_t n_nodes,
                       ITYPE_t n_features) except -1:
    bt.node_bounds = np.zeros((1, n_nodes, n_features), dtype=DTYPE)
    return 0

cdef int init_node(BinaryTree bt, ITYPE_t i_node,
                   ITYPE_t idx_start, ITYPE_t idx_end) except -1:
    cdef ITYPE_t n_features = bt.data.shape[1]
    cdef ITYPE_t n_points = idx_end - idx_start

    cdef ITYPE_t i, j
    cdef DTYPE_t radius
    cdef DTYPE_t *this_pt

    cdef ITYPE_t* idx_array = &bt.idx_array[0]
    cdef DTYPE_t* data = &bt.data[0, 0]
    cdef DTYPE_t* centroid = &bt.node_bounds[0, i_node, 0]

    # determine Node centroid
    for j in range(n_features):
        centroid[j] = 0

    for i in range(idx_start, idx_end):
        this_pt = data + n_features * idx_array[i]
        for j from 0 <= j < n_features:
            centroid[j] += this_pt[j]

    for j in range(n_features):
        centroid[j] /= n_points

    # determine Node radius
    radius = 0
    for i in range(idx_start, idx_end):
        radius = fmax(radius,
                      bt.rdist(centroid,
                                 data + n_features * idx_array[i],
                                 n_features))

    bt.node_data[i_node].radius = bt.dm.rdist_to_dist(radius)
    bt.node_data[i_node].idx_start = idx_start
    bt.node_data[i_node].idx_end = idx_end
    return 0

cdef inline DTYPE_t min_dist(BinaryTree bt, ITYPE_t i_node,
                             DTYPE_t* pt) except -1:
    cdef DTYPE_t dist_pt = bt.dist(pt, &bt.node_bounds[0, i_node, 0],
                                     bt.data.shape[1])
    return fmax(0, dist_pt - bt.node_data[i_node].radius)

cdef inline DTYPE_t max_dist(BinaryTree bt, ITYPE_t i_node,
                             DTYPE_t* pt) except -1:
    cdef DTYPE_t dist_pt = bt.dist(pt, &bt.node_bounds[0, i_node, 0],
                                     bt.data.shape[1])
    return dist_pt + bt.node_data[i_node].radius

cdef inline int min_max_dist(BinaryTree bt, ITYPE_t i_node, DTYPE_t* pt,
                             DTYPE_t* min_dist, DTYPE_t* max_dist) except -1:
    cdef DTYPE_t dist_pt = bt.dist(pt, &bt.node_bounds[0, i_node, 0],
                                     bt.data.shape[1])
    cdef DTYPE_t rad = bt.node_data[i_node].radius
    min_dist[0] = fmax(0, dist_pt - rad)
    max_dist[0] = dist_pt + rad
    return 0

cdef inline DTYPE_t min_rdist(BinaryTree bt, ITYPE_t i_node,
                              DTYPE_t* pt) except -1:
    if bt.euclidean:
        return euclidean_dist_to_rdist(min_dist(bt, i_node, pt))
    else:
        return bt.dm.dist_to_rdist(min_dist(bt, i_node, pt))

cdef inline DTYPE_t max_rdist(BinaryTree bt, ITYPE_t i_node,
                              DTYPE_t* pt) except -1:
    if bt.euclidean:
        return euclidean_dist_to_rdist(max_dist(bt, i_node, pt))
    else:
        return bt.dm.dist_to_rdist(max_dist(bt, i_node, pt))

cdef inline DTYPE_t min_dist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                  BinaryTree bt2, ITYPE_t i_node2) except -1:
    cdef DTYPE_t dist_pt = bt1.dist(&bt2.node_bounds[0, i_node2, 0],
                                    &bt1.node_bounds[0, i_node1, 0],
                                    bt1.data.shape[1])
    return fmax(0, (dist_pt - bt1.node_data[i_node1].radius
                    - bt2.node_data[i_node2].radius))

cdef inline DTYPE_t max_dist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                  BinaryTree bt2, ITYPE_t i_node2) except -1:
    cdef DTYPE_t dist_pt = bt1.dist(&bt2.node_bounds[0, i_node2, 0],
                                    &bt1.node_bounds[0, i_node1, 0],
                                    bt1.data.shape[1])
    return (dist_pt + bt1.node_data[i_node1].radius
            + bt2.node_data[i_node2].radius)

cdef inline DTYPE_t min_rdist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                   BinaryTree bt2, ITYPE_t i_node2) except -1:
    if bt1.euclidean:
        return euclidean_dist_to_rdist(min_dist_dual(bt1, i_node1,
                                                     bt2, i_node2))
    else:
        return bt1.dm.dist_to_rdist(min_dist_dual(bt1, i_node1,
                                                  bt2, i_node2))

cdef inline DTYPE_t max_rdist_dual(BinaryTree bt1, ITYPE_t i_node1,
                                   BinaryTree bt2, ITYPE_t i_node2) except -1:
    if bt1.euclidean:
        return euclidean_dist_to_rdist(max_dist_dual(bt1, i_node1,
                                                     bt2, i_node2))
    else:
        return bt1.dm.dist_to_rdist(max_dist_dual(bt1, i_node1,
                                                  bt2, i_node2))
