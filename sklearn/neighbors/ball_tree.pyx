#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

# Author: Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD 3 clause

__all__ = ['BallTree']

DOC_DICT = {'BinaryTree': 'BallTree', 'binary_tree': 'ball_tree'}

VALID_METRICS = ['EuclideanDistance', 'SEuclideanDistance',
                 'ManhattanDistance', 'ChebyshevDistance',
                 'MinkowskiDistance', 'WMinkowskiDistance',
                 'MahalanobisDistance', 'HammingDistance',
                 'CanberraDistance', 'BrayCurtisDistance',
                 'JaccardDistance', 'MatchingDistance',
                 'DiceDistance', 'KulsinskiDistance',
                 'RogersTanimotoDistance', 'RussellRaoDistance',
                 'SokalMichenerDistance', 'SokalSneathDistance',
                 'PyFuncDistance', 'HaversineDistance']


include "binary_tree.pxi"

# Inherit BallTree from BinaryTree
cdef class BallTree(BinaryTree):
    __doc__ = CLASS_DOC.format(**DOC_DICT)
    pass


#----------------------------------------------------------------------
# The functions below specialized the Binary Tree as a Ball Tree
#
#   Note that these functions use the concept of "reduced distance".
#   The reduced distance, defined for some metrics, is a quantity which
#   is more efficient to compute than the distance, but preserves the
#   relative rankings of the true distance.  For example, the reduced
#   distance for the Euclidean metric is the squared-euclidean distance.
#   For some metrics, the reduced distance is simply the distance.

cdef int allocate_data(BinaryTree tree, ITYPE_t n_nodes,
                       ITYPE_t n_features) except -1:
    """Allocate arrays needed for the KD Tree"""
    tree.node_bounds_arr = np.zeros((1, n_nodes, n_features), dtype=DTYPE)
    tree.node_bounds = get_memview_DTYPE_3D(tree.node_bounds_arr)
    return 0


cdef int init_node(BinaryTree tree, ITYPE_t i_node,
                   ITYPE_t idx_start, ITYPE_t idx_end) except -1:
    """Initialize the node for the dataset stored in tree.data"""
    cdef ITYPE_t n_features = tree.data.shape[1]
    cdef ITYPE_t n_points = idx_end - idx_start

    cdef ITYPE_t i, j
    cdef DTYPE_t radius
    cdef DTYPE_t *this_pt

    cdef ITYPE_t* idx_array = &tree.idx_array[0]
    cdef DTYPE_t* data = &tree.data[0, 0]
    cdef DTYPE_t* centroid = &tree.node_bounds[0, i_node, 0]

    cdef bint with_sample_weight = tree.sample_weight is not None
    cdef DTYPE_t* sample_weight
    cdef DTYPE_t sum_weight_node
    if with_sample_weight:
        sample_weight = &tree.sample_weight[0]

    # determine Node centroid
    for j in range(n_features):
        centroid[j] = 0

    if with_sample_weight:
        sum_weight_node = 0
        for i in range(idx_start, idx_end):
            sum_weight_node += sample_weight[idx_array[i]]
            this_pt = data + n_features * idx_array[i]
            for j from 0 <= j < n_features:
                centroid[j] += this_pt[j] * sample_weight[idx_array[i]]

        for j in range(n_features):
            centroid[j] /= sum_weight_node
    else:
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
                      tree.rdist(centroid,
                                 data + n_features * idx_array[i],
                                 n_features))

    tree.node_data[i_node].radius = tree.dist_metric._rdist_to_dist(radius)
    tree.node_data[i_node].idx_start = idx_start
    tree.node_data[i_node].idx_end = idx_end
    return 0


cdef inline DTYPE_t min_dist(BinaryTree tree, ITYPE_t i_node,
                             DTYPE_t* pt) nogil except -1:
    """Compute the minimum distance between a point and a node"""
    cdef DTYPE_t dist_pt = tree.dist(pt, &tree.node_bounds[0, i_node, 0],
                                     tree.data.shape[1])
    return fmax(0, dist_pt - tree.node_data[i_node].radius)


cdef inline DTYPE_t max_dist(BinaryTree tree, ITYPE_t i_node,
                             DTYPE_t* pt) except -1:
    """Compute the maximum distance between a point and a node"""
    cdef DTYPE_t dist_pt = tree.dist(pt, &tree.node_bounds[0, i_node, 0],
                                     tree.data.shape[1])
    return dist_pt + tree.node_data[i_node].radius


cdef inline int min_max_dist(BinaryTree tree, ITYPE_t i_node, DTYPE_t* pt,
                             DTYPE_t* min_dist, DTYPE_t* max_dist) nogil except -1:
    """Compute the minimum and maximum distance between a point and a node"""
    cdef DTYPE_t dist_pt = tree.dist(pt, &tree.node_bounds[0, i_node, 0],
                                     tree.data.shape[1])
    cdef DTYPE_t rad = tree.node_data[i_node].radius
    min_dist[0] = fmax(0, dist_pt - rad)
    max_dist[0] = dist_pt + rad
    return 0


cdef inline DTYPE_t min_rdist(BinaryTree tree, ITYPE_t i_node,
                              DTYPE_t* pt) nogil except -1:
    """Compute the minimum reduced-distance between a point and a node"""
    if tree.euclidean:
        return euclidean_dist_to_rdist(min_dist(tree, i_node, pt))
    else:
        return tree.dist_metric._dist_to_rdist(min_dist(tree, i_node, pt))


cdef inline DTYPE_t max_rdist(BinaryTree tree, ITYPE_t i_node,
                              DTYPE_t* pt) except -1:
    """Compute the maximum reduced-distance between a point and a node"""
    if tree.euclidean:
        return euclidean_dist_to_rdist(max_dist(tree, i_node, pt))
    else:
        return tree.dist_metric._dist_to_rdist(max_dist(tree, i_node, pt))


cdef inline DTYPE_t min_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                  BinaryTree tree2, ITYPE_t i_node2) except -1:
    """compute the minimum distance between two nodes"""
    cdef DTYPE_t dist_pt = tree1.dist(&tree2.node_bounds[0, i_node2, 0],
                                      &tree1.node_bounds[0, i_node1, 0],
                                      tree1.data.shape[1])
    return fmax(0, (dist_pt - tree1.node_data[i_node1].radius
                    - tree2.node_data[i_node2].radius))


cdef inline DTYPE_t max_dist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                  BinaryTree tree2, ITYPE_t i_node2) except -1:
    """compute the maximum distance between two nodes"""
    cdef DTYPE_t dist_pt = tree1.dist(&tree2.node_bounds[0, i_node2, 0],
                                      &tree1.node_bounds[0, i_node1, 0],
                                      tree1.data.shape[1])
    return (dist_pt + tree1.node_data[i_node1].radius
            + tree2.node_data[i_node2].radius)


cdef inline DTYPE_t min_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                   BinaryTree tree2, ITYPE_t i_node2) except -1:
    """compute the minimum reduced distance between two nodes"""
    if tree1.euclidean:
        return euclidean_dist_to_rdist(min_dist_dual(tree1, i_node1,
                                                     tree2, i_node2))
    else:
        return tree1.dist_metric._dist_to_rdist(min_dist_dual(tree1, i_node1,
                                                              tree2, i_node2))


cdef inline DTYPE_t max_rdist_dual(BinaryTree tree1, ITYPE_t i_node1,
                                   BinaryTree tree2, ITYPE_t i_node2) except -1:
    """compute the maximum reduced distance between two nodes"""
    if tree1.euclidean:
        return euclidean_dist_to_rdist(max_dist_dual(tree1, i_node1,
                                                     tree2, i_node2))
    else:
        return tree1.dist_metric._dist_to_rdist(max_dist_dual(tree1, i_node1,
                                                              tree2, i_node2))
