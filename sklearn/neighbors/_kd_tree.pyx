# By Jake Vanderplas (2013) <jakevdp@cs.washington.edu>
# written for the scikit-learn project
# License: BSD

__all__ = ['KDTree']

DOC_DICT = {'BinaryTree': 'KDTree', 'binary_tree': 'kd_tree'}

VALID_METRICS = ['EuclideanDistance', 'ManhattanDistance',
                 'ChebyshevDistance', 'MinkowskiDistance']

cimport numpy as np

from numpy.math cimport INFINITY, isinf
from libc.math cimport fabs, pow

from ._binary_tree cimport BinaryTree
from ..metrics._dist_metrics cimport get_valid_metric_ids
from ..utils._typedefs cimport DTYPE_t, ITYPE_t

from ._binary_tree import CLASS_DOC
from ..utils._typedefs import DTYPE, ITYPE


import numpy as np

np.import_array()


cdef class KDTree(BinaryTree):
    __doc__ = CLASS_DOC.format(**DOC_DICT)

    valid_metrics = get_valid_metric_ids(VALID_METRICS)

    # Implement abstract methods.
    #
    # Note that these functions use the concept of "reduced distance".
    # The reduced distance, defined for some metrics, is a quantity which
    # is more efficient to compute than the distance, but preserves the
    # relative rankings of the true distance.  For example, the reduced
    # distance for the Euclidean metric is the squared-euclidean distance.
    # For some metrics, the reduced distance is simply the distance.

    cdef int allocate_data(self, ITYPE_t n_nodes, ITYPE_t n_features) except -1:
        self.node_bounds_arr = np.zeros((2, n_nodes, n_features), dtype=DTYPE)
        self.node_bounds = self.node_bounds_arr
        return 0

    cdef int init_node(self, ITYPE_t i_node,
                       ITYPE_t idx_start, ITYPE_t idx_end) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef ITYPE_t i, j
        cdef DTYPE_t rad = 0

        cdef DTYPE_t* lower_bounds = &self.node_bounds[0, i_node, 0]
        cdef DTYPE_t* upper_bounds = &self.node_bounds[1, i_node, 0]
        cdef DTYPE_t* data = &self.data[0, 0]
        cdef ITYPE_t* idx_array = &self.idx_array[0]

        cdef DTYPE_t* data_row

        # determine Node bounds
        for j in range(n_features):
            lower_bounds[j] = INFINITY
            upper_bounds[j] = -INFINITY

        # Compute the actual data range.  At build time, this is slightly
        # slower than using the previously-computed bounds of the parent node,
        # but leads to more compact trees and thus faster queries.
        for i in range(idx_start, idx_end):
            data_row = data + idx_array[i] * n_features
            for j in range(n_features):
                lower_bounds[j] = min(lower_bounds[j], data_row[j])
                upper_bounds[j] = max(upper_bounds[j], data_row[j])
            if isinf(self.dist_metric.p) == 1:
                rad = max(rad, 0.5 * (upper_bounds[j] - lower_bounds[j]))
            else:
                rad += pow(0.5 * fabs(upper_bounds[j] - lower_bounds[j]),
                           self.dist_metric.p)

        self.node_data[i_node].idx_start = idx_start
        self.node_data[i_node].idx_end = idx_end

        # The radius will hold the size of the circumscribed hypersphere
        # measured with the specified metric: in querying, this is used as a
        # measure of the size of each node when deciding which nodes to split.
        self.node_data[i_node].radius = pow(rad, 1. / self.dist_metric.p)
        return 0

    cdef DTYPE_t min_rdist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]
        cdef DTYPE_t d, d_lo, d_hi, rdist=0.0
        cdef ITYPE_t j

        if isinf(self.dist_metric.p) == 1:
            for j in range(n_features):
                d_lo = self.node_bounds[0, i_node, j] - pt[j]
                d_hi = pt[j] - self.node_bounds[1, i_node, j]
                d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
                rdist = max(rdist, 0.5 * d)
        else:
            # here we'll use the fact that x + abs(x) = 2 * max(x, 0)
            for j in range(n_features):
                d_lo = self.node_bounds[0, i_node, j] - pt[j]
                d_hi = pt[j] - self.node_bounds[1, i_node, j]
                d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
                rdist += pow(0.5 * d, self.dist_metric.p)

        return rdist

    cdef DTYPE_t min_dist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        if isinf(self.dist_metric.p) == 1:
            return self.min_rdist(i_node, pt)
        else:
            return pow(self.min_rdist(i_node, pt), 1. / self.dist_metric.p)

    cdef DTYPE_t max_rdist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t d, d_lo, d_hi, rdist=0.0
        cdef ITYPE_t j

        if isinf(self.dist_metric.p) == 1:
            for j in range(n_features):
                rdist = max(rdist, fabs(pt[j] - self.node_bounds[0, i_node, j]))
                rdist = max(rdist, fabs(pt[j] - self.node_bounds[1, i_node, j]))
        else:
            for j in range(n_features):
                d_lo = fabs(pt[j] - self.node_bounds[0, i_node, j])
                d_hi = fabs(pt[j] - self.node_bounds[1, i_node, j])
                rdist += pow(max(d_lo, d_hi), self.dist_metric.p)

        return rdist

    cdef DTYPE_t max_dist(self, ITYPE_t i_node, DTYPE_t* pt) nogil except -1:
        if isinf(self.dist_metric.p) == 1:
            return self.max_rdist(i_node, pt)
        else:
            return pow(self.max_rdist(i_node, pt), 1. / self.dist_metric.p)

    cdef int min_max_dist(self, ITYPE_t i_node, DTYPE_t* pt,
                          DTYPE_t* min_dist, DTYPE_t* max_dist) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t d, d_lo, d_hi
        cdef ITYPE_t j

        min_dist[0] = 0.0
        max_dist[0] = 0.0

        if isinf(self.dist_metric.p) == 1:
            for j in range(n_features):
                d_lo = self.node_bounds[0, i_node, j] - pt[j]
                d_hi = pt[j] - self.node_bounds[1, i_node, j]
                d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
                min_dist[0] = max(min_dist[0], 0.5 * d)
                max_dist[0] = max(max_dist[0],
                                  fabs(pt[j] - self.node_bounds[0, i_node, j]))
                max_dist[0] = max(max_dist[0],
                                  fabs(pt[j] - self.node_bounds[1, i_node, j]))
        else:
            # as above, use the fact that x + abs(x) = 2 * max(x, 0)
            for j in range(n_features):
                d_lo = self.node_bounds[0, i_node, j] - pt[j]
                d_hi = pt[j] - self.node_bounds[1, i_node, j]
                d = (d_lo + fabs(d_lo)) + (d_hi + fabs(d_hi))
                min_dist[0] += pow(0.5 * d, self.dist_metric.p)
                max_dist[0] += pow(max(fabs(d_lo), fabs(d_hi)),
                                   self.dist_metric.p)

            min_dist[0] = pow(min_dist[0], 1. / self.dist_metric.p)
            max_dist[0] = pow(max_dist[0], 1. / self.dist_metric.p)

        return 0

    cdef DTYPE_t min_rdist_dual(self, ITYPE_t i_node1,
                                BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t d, d1, d2, rdist=0.0
        cdef ITYPE_t j

        if isinf(self.dist_metric.p) == 1:
            for j in range(n_features):
                d1 = (self.node_bounds[0, i_node1, j]
                     - other.node_bounds[1, i_node2, j])
                d2 = (other.node_bounds[0, i_node2, j]
                     - self.node_bounds[1, i_node1, j])
                d = (d1 + fabs(d1)) + (d2 + fabs(d2))

                rdist = max(rdist, 0.5 * d)
        else:
            # here we'll use the fact that x + abs(x) = 2 * max(x, 0)
            for j in range(n_features):
                d1 = (self.node_bounds[0, i_node1, j]
                      - other.node_bounds[1, i_node2, j])
                d2 = (other.node_bounds[0, i_node2, j]
                      - self.node_bounds[1, i_node1, j])
                d = (d1 + fabs(d1)) + (d2 + fabs(d2))

                rdist += pow(0.5 * d, self.dist_metric.p)

        return rdist

    cdef DTYPE_t min_dist_dual(self, ITYPE_t i_node1,
                               BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t rd = self.min_rdist_dual(i_node1, other, i_node2)
        return self.dist_metric._rdist_to_dist(rd)

    cdef DTYPE_t max_rdist_dual(self, ITYPE_t i_node1,
                                BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef ITYPE_t n_features = self.data.shape[1]

        cdef DTYPE_t d, d1, d2, rdist=0.0
        cdef ITYPE_t j

        if isinf(self.dist_metric.p) == 1:
            for j in range(n_features):
                rdist = max(rdist, fabs(self.node_bounds[0, i_node1, j]
                                        - other.node_bounds[1, i_node2, j]))
                rdist = max(rdist, fabs(self.node_bounds[1, i_node1, j]
                                        - other.node_bounds[0, i_node2, j]))
        else:
            for j in range(n_features):
                d1 = fabs(self.node_bounds[0, i_node1, j]
                          - other.node_bounds[1, i_node2, j])
                d2 = fabs(self.node_bounds[1, i_node1, j]
                          - other.node_bounds[0, i_node2, j])
                rdist += pow(max(d1, d2), self.dist_metric.p)

        return rdist

    cdef DTYPE_t max_dist_dual(self, ITYPE_t i_node1,
                               BinaryTree other, ITYPE_t i_node2) nogil except -1:
        cdef DTYPE_t rd = self.max_rdist_dual(i_node1, other, i_node2)
        return self.dist_metric._rdist_to_dist(rd)
