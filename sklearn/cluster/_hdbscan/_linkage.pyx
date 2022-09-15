# Minimum spanning tree single linkage implementation for hdbscan
# Authors: Leland McInnes, Steve Astels
# License: 3-clause BSD

import numpy as np
cimport numpy as cnp
import cython

from libc.float cimport DBL_MAX

from ...metrics._dist_metrics cimport DistanceMetric


cpdef cnp.ndarray[cnp.double_t, ndim=2] mst_from_distance_matrix(
    cnp.ndarray[cnp.double_t, ndim=2] distance_matrix
):

    cdef:
        cnp.ndarray[cnp.intp_t, ndim=1] node_labels
        cnp.ndarray[cnp.intp_t, ndim=1] current_labels
        cnp.ndarray[cnp.double_t, ndim=1] current_distances
        cnp.ndarray[cnp.double_t, ndim=1] left
        cnp.ndarray[cnp.double_t, ndim=1] right
        cnp.ndarray[cnp.double_t, ndim=2] result

        cnp.ndarray label_filter

        cnp.intp_t current_node
        cnp.intp_t new_node_index
        cnp.intp_t new_node
        cnp.intp_t i

    result = np.zeros((distance_matrix.shape[0] - 1, 3))
    node_labels = np.arange(distance_matrix.shape[0], dtype=np.intp)
    current_node = 0
    current_distances = np.infty * np.ones(distance_matrix.shape[0])
    current_labels = node_labels
    for i in range(1, node_labels.shape[0]):
        label_filter = current_labels != current_node
        current_labels = current_labels[label_filter]
        left = current_distances[label_filter]
        right = distance_matrix[current_node][current_labels]
        current_distances = np.where(left < right, left, right)

        new_node_index = np.argmin(current_distances)
        new_node = current_labels[new_node_index]
        result[i - 1, 0] = <double> current_node
        result[i - 1, 1] = <double> new_node
        result[i - 1, 2] = current_distances[new_node_index]
        current_node = new_node

    return result


cpdef cnp.ndarray[cnp.double_t, ndim=2] mst_from_data_matrix(
    cnp.ndarray[cnp.double_t, ndim=2, mode='c'] raw_data,
    cnp.ndarray[cnp.double_t, ndim=1, mode='c'] core_distances,
    DistanceMetric dist_metric,
    cnp.double_t alpha=1.0
):

    cdef:
        cnp.ndarray[cnp.double_t, ndim=1] current_distances_arr
        cnp.ndarray[cnp.double_t, ndim=1] current_sources_arr
        cnp.ndarray[cnp.int8_t, ndim=1] in_tree_arr
        cnp.ndarray[cnp.double_t, ndim=2] result_arr

        cnp.double_t * current_distances
        cnp.double_t * current_sources
        cnp.double_t * current_core_distances
        cnp.double_t * raw_data_ptr
        cnp.int8_t * in_tree
        cnp.double_t[:, ::1] raw_data_view
        cnp.double_t[:, ::1] result

        cnp.ndarray label_filter

        cnp.intp_t current_node
        cnp.intp_t source_node
        cnp.intp_t right_node
        cnp.intp_t left_node
        cnp.intp_t new_node
        cnp.intp_t i
        cnp.intp_t j
        cnp.intp_t dim
        cnp.intp_t num_features

        double current_node_core_distance
        double right_value
        double left_value
        double core_value
        double new_distance

    dim = raw_data.shape[0]
    num_features = raw_data.shape[1]

    raw_data_view = (<cnp.double_t[:raw_data.shape[0], :raw_data.shape[1]:1]> (
        <cnp.double_t *> raw_data.data))
    raw_data_ptr = (<cnp.double_t *> &raw_data_view[0, 0])

    result_arr = np.zeros((dim - 1, 3))
    in_tree_arr = np.zeros(dim, dtype=np.int8)
    current_node = 0
    current_distances_arr = np.infty * np.ones(dim)
    current_sources_arr = np.ones(dim)

    result = (<cnp.double_t[:dim - 1, :3:1]> (<cnp.double_t *> result_arr.data))
    in_tree = (<cnp.int8_t *> in_tree_arr.data)
    current_distances = (<cnp.double_t *> current_distances_arr.data)
    current_sources = (<cnp.double_t *> current_sources_arr.data)
    current_core_distances = (<cnp.double_t *> core_distances.data)

    for i in range(1, dim):

        in_tree[current_node] = 1

        current_node_core_distance = current_core_distances[current_node]

        new_distance = DBL_MAX
        source_node = 0
        new_node = 0

        for j in range(dim):
            if in_tree[j]:
                continue

            right_value = current_distances[j]
            right_source = current_sources[j]

            left_value = dist_metric.dist(&raw_data_ptr[num_features *
                                                        current_node],
                                          &raw_data_ptr[num_features * j],
                                          num_features)
            left_source = current_node

            if alpha != 1.0:
                left_value /= alpha

            core_value = core_distances[j]
            if (current_node_core_distance > right_value or
                    core_value > right_value or
                    left_value > right_value):
                if right_value < new_distance:
                    new_distance = right_value
                    source_node = right_source
                    new_node = j
                continue

            if core_value > current_node_core_distance:
                if core_value > left_value:
                    left_value = core_value
            else:
                if current_node_core_distance > left_value:
                    left_value = current_node_core_distance

            if left_value < right_value:
                current_distances[j] = left_value
                current_sources[j] = left_source
                if left_value < new_distance:
                    new_distance = left_value
                    source_node = left_source
                    new_node = j
            else:
                if right_value < new_distance:
                    new_distance = right_value
                    source_node = right_source
                    new_node = j

        result[i - 1, 0] = <double> source_node
        result[i - 1, 1] = <double> new_node
        result[i - 1, 2] = new_distance
        current_node = new_node

    return result_arr


cdef class UnionFind (object):

    cdef:
        cnp.ndarray parent_arr
        cnp.ndarray size_arr
        cnp.intp_t next_label
        cnp.intp_t *parent
        cnp.intp_t *size

    def __init__(self, N):
        self.parent_arr = -1 * np.ones(2 * N - 1, dtype=np.intp, order='C')
        self.next_label = N
        self.size_arr = np.hstack((np.ones(N, dtype=np.intp),
                                   np.zeros(N-1, dtype=np.intp)))
        self.parent = (<cnp.intp_t *> self.parent_arr.data)
        self.size = (<cnp.intp_t *> self.size_arr.data)

    cdef void union(self, cnp.intp_t m, cnp.intp_t n):
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.parent[m] = self.next_label
        self.parent[n] = self.next_label
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.next_label += 1

        return

    @cython.wraparound(True)
    cdef cnp.intp_t fast_find(self, cnp.intp_t n):
        cdef cnp.intp_t p
        p = n
        while self.parent_arr[n] != -1:
            n = self.parent_arr[n]
        # label up to the root
        while self.parent_arr[p] != n:
            p, self.parent_arr[p] = self.parent_arr[p], n
        return n

@cython.wraparound(True)
cpdef cnp.ndarray[cnp.double_t, ndim=2] label(cnp.double_t[:,:] L):

    cdef:
        cnp.ndarray[cnp.double_t, ndim=2] result_arr
        cnp.double_t[:, ::1] result

        cnp.intp_t N, a, aa, b, bb, index
        cnp.double_t delta

    result_arr = np.zeros((L.shape[0], L.shape[1] + 1))
    result = (<cnp.double_t[:L.shape[0], :4:1]> (
        <cnp.double_t *> result_arr.data))
    N = L.shape[0] + 1
    U = UnionFind(N)

    for index in range(L.shape[0]):

        a = <cnp.intp_t> L[index, 0]
        b = <cnp.intp_t> L[index, 1]
        delta = L[index, 2]

        aa, bb = U.fast_find(a), U.fast_find(b)

        result[index][0] = aa
        result[index][1] = bb
        result[index][2] = delta
        result[index][3] = U.size[aa] + U.size[bb]

        U.union(aa, bb)

    return result_arr
