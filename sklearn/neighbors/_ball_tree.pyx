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

cimport numpy as np

from libc.math cimport fmax, fmin

from ._binary_tree cimport BinaryTree
from ..utils._typedefs cimport DTYPE_t, ITYPE_t
from ..metrics._dist_metrics cimport get_valid_metric_ids

from ._binary_tree import CLASS_DOC
from ..utils._typedefs import DTYPE, ITYPE

import numpy as np

np.import_array()

cdef class BallTree(BinaryTree):
    __doc__ = CLASS_DOC.format(**DOC_DICT)

    valid_metrics = get_valid_metric_ids(VALID_METRICS)

    # Implementations of abstract methods.
    #
    # Note that these functions use the concept of "reduced distance".
    # The reduced distance, defined for some metrics, is a quantity which
    # is more efficient to compute than the distance, but preserves the
    # relative rankings of the true distance.  For example, the reduced
    # distance for the Euclidean metric is the squared-euclidean distance.
    # For some metrics, the reduced distance is simply the distance.

    cdef int allocate_data(self, ITYPE_t n_nodes, ITYPE_t n_features) except -1:
        self.node_bounds_arr = np.zeros((1, n_nodes, n_features), dtype=DTYPE)
        self.node_bounds = self.node_bounds_arr
        return 0

    cdef int init_node(self, ITYPE_t i_node,
                       ITYPE_t idx_start, ITYPE_t idx_end) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t n_points = idx_end - idx_start
        cdef ITYPE_t i, j
        cdef DTYPE_t radius
        cdef DTYPE_t *this_pt

        cdef ITYPE_t* idx_array = &self.idx_array[0]
        cdef DTYPE_t* data = &self.data[0, 0]
        cdef DTYPE_t* centroid = &self.node_bounds[0, i_node, 0]

        # determine Node centroid
        for j in range(n_features):
            centroid[j] = 0

        for i in range(idx_start, idx_end):
            this_pt = data + n_features * idx_array[i]
            for j in range(n_features):
                centroid[j] += this_pt[j]

        for j in range(n_features):
            centroid[j] /= n_points

        # determine Node radius
        radius = 0
        for i in range(idx_start, idx_end):
            radius = fmax(radius,
                          self.rdist(centroid,
                                     data + n_features * idx_array[i],
                                     n_features))

        self.node_data[i_node].radius = self.dist_metric._rdist_to_dist(radius)
        self.node_data[i_node].idx_start = idx_start
        self.node_data[i_node].idx_end = idx_end
        return 0

    cdef DTYPE_t min_dist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef DTYPE_t dist_pt = self.dist(pt, &self.node_bounds[0, i_node, 0],
                                         self.data.shape[1])
        return fmax(0, dist_pt - self.node_data[i_node].radius)

    cdef DTYPE_t max_dist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef DTYPE_t dist_pt = self.dist(pt, &self.node_bounds[0, i_node, 0],
                                         self.data.shape[1])
        return dist_pt + self.node_data[i_node].radius

    cdef int min_max_dist(self, ITYPE_t i_node, DTYPE_t* pt,
                          DTYPE_t* min_dist, DTYPE_t* max_dist) nogil except -1:
        cdef DTYPE_t dist_pt = self.dist(pt, &self.node_bounds[0, i_node, 0],
                                         self.data.shape[1])
        cdef DTYPE_t rad = self.node_data[i_node].radius
        min_dist[0] = fmax(0, dist_pt - rad)
        max_dist[0] = dist_pt + rad
        return 0

    cdef DTYPE_t min_rdist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef DTYPE_t d
        if self.euclidean:
            d = self.min_dist(i_node, pt)
            return d * d
        else:
            return self.dist_metric._dist_to_rdist(self.min_dist(i_node, pt))

    cdef DTYPE_t max_rdist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef DTYPE_t d
        if self.euclidean:
            d = self.max_dist(i_node, pt)
            return d * d
        else:
            return self.dist_metric._dist_to_rdist(self.max_dist(i_node, pt))

    cdef DTYPE_t min_dist_dual(self, ITYPE_t i_node1,
                               BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t dist_pt = self.dist(&other.node_bounds[0, i_node2, 0],
                                         &self.node_bounds[0, i_node1, 0],
                                         self.data.shape[1])
        return fmax(0, (dist_pt - self.node_data[i_node1].radius
                        - other.node_data[i_node2].radius))

    cdef DTYPE_t max_dist_dual(self, ITYPE_t i_node1,
                               BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t dist_pt = self.dist(&other.node_bounds[0, i_node2, 0],
                                         &self.node_bounds[0, i_node1, 0],
                                         self.data.shape[1])
        return (dist_pt + self.node_data[i_node1].radius
                + other.node_data[i_node2].radius)

    cdef DTYPE_t min_rdist_dual(self, ITYPE_t i_node1,
                                BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t d = self.min_dist_dual(i_node1, other, i_node2)
        if self.euclidean:
            return d * d
        else:
            return self.dist_metric._dist_to_rdist(d)

    cdef DTYPE_t max_rdist_dual(self, ITYPE_t i_node1,
                                BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t d = self.max_dist_dual(i_node1, other, i_node2)
        if self.euclidean:
            return d * d
        else:
            return self.dist_metric._dist_to_rdist(d)
