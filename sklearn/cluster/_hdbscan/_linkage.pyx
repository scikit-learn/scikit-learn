# Minimum spanning tree single linkage implementation for hdbscan
# Authors: Leland McInnes <leland.mcinnes@gmail.com>
#          Steve Astels <sastels@gmail.com>
# License: 3-clause BSD

cimport numpy as cnp

import numpy as np
import cython

from libc.float cimport DBL_MAX

from ...metrics._dist_metrics cimport DistanceMetric
from ...cluster._hierarchical_fast cimport UnionFind
from ...utils._typedefs cimport ITYPE_t, DTYPE_t
from ...utils._typedefs import ITYPE, DTYPE

cpdef cnp.ndarray[cnp.float64_t, ndim=2] mst_from_distance_matrix(
    cnp.ndarray[cnp.float64_t, ndim=2] distance_matrix
):

    cdef:
        cnp.ndarray[cnp.intp_t, ndim=1] node_labels
        cnp.ndarray[cnp.intp_t, ndim=1] current_labels
        cnp.ndarray[cnp.float64_t, ndim=1] current_distances, left, right
        cnp.ndarray[cnp.float64_t, ndim=2] result

        cnp.ndarray label_filter

        cnp.intp_t n_samples = distance_matrix.shape[0]
        cnp.intp_t current_node, new_node_index, new_node, i

    result = np.zeros((n_samples - 1, 3))
    node_labels = np.arange(n_samples, dtype=np.intp)
    current_node = 0
    current_distances = np.infty * np.ones(n_samples)
    current_labels = node_labels
    for i in range(1, n_samples):
        label_filter = current_labels != current_node
        current_labels = current_labels[label_filter]
        left = current_distances[label_filter]
        right = distance_matrix[current_node][current_labels]
        current_distances = np.minimum(left, right)

        new_node_index = np.argmin(current_distances)
        new_node = current_labels[new_node_index]
        result[i - 1, 0] = <cnp.float64_t> current_node
        result[i - 1, 1] = <cnp.float64_t> new_node
        result[i - 1, 2] = current_distances[new_node_index]
        current_node = new_node

    return result


cpdef cnp.ndarray[cnp.float64_t, ndim=2] mst_from_data_matrix(
    cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] raw_data,
    cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] core_distances,
    DistanceMetric dist_metric,
    cnp.float64_t alpha=1.0
):

    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1] current_distances_arr
        cnp.ndarray[cnp.float64_t, ndim=1] current_sources_arr
        cnp.ndarray[cnp.int8_t, ndim=1] in_tree_arr
        cnp.ndarray[cnp.float64_t, ndim=2] result_arr

        cnp.float64_t * current_distances
        cnp.float64_t * current_sources
        cnp.float64_t * current_core_distances
        cnp.float64_t * raw_data_ptr
        cnp.int8_t * in_tree
        cnp.float64_t[:, ::1] raw_data_view
        cnp.float64_t[:, ::1] result

        cnp.ndarray label_filter

        cnp.intp_t current_node, source_node, right_node, left_node, new_node
        cnp.intp_t i, j, n_samples, num_features

        cnp.float64_t current_node_core_distance, new_distance
        cnp.float64_t right_value, left_value, core_value

    n_samples = raw_data.shape[0]
    num_features = raw_data.shape[1]

    raw_data_view = (<cnp.float64_t[:n_samples, :num_features:1]> (
        <cnp.float64_t *> raw_data.data))
    raw_data_ptr = (<cnp.float64_t *> &raw_data_view[0, 0])

    result_arr = np.zeros((n_samples - 1, 3))
    in_tree_arr = np.zeros(n_samples, dtype=np.int8)
    current_node = 0
    current_distances_arr = np.infty * np.ones(n_samples)
    current_sources_arr = np.ones(n_samples)

    result = (<cnp.float64_t[:n_samples - 1, :3:1]> (<cnp.float64_t *> result_arr.data))
    in_tree = (<cnp.int8_t *> in_tree_arr.data)
    current_distances = (<cnp.float64_t *> current_distances_arr.data)
    current_sources = (<cnp.float64_t *> current_sources_arr.data)
    current_core_distances = (<cnp.float64_t *> core_distances.data)

    for i in range(1, n_samples):

        in_tree[current_node] = 1

        current_node_core_distance = current_core_distances[current_node]

        new_distance = DBL_MAX
        source_node = 0
        new_node = 0

        for j in range(n_samples):
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

        result[i - 1, 0] = <cnp.float64_t> source_node
        result[i - 1, 1] = <cnp.float64_t> new_node
        result[i - 1, 2] = new_distance
        current_node = new_node

    return result_arr

@cython.wraparound(True)
cpdef cnp.ndarray[cnp.float64_t, ndim=2] label(cnp.float64_t[:,:] L):

    cdef:
        cnp.ndarray[cnp.float64_t, ndim=2] result_arr
        cnp.float64_t[:, ::1] result

        cnp.intp_t N, a, aa, b, bb, index
        cnp.float64_t delta

    result_arr = np.zeros((L.shape[0], L.shape[1] + 1))
    result = (<cnp.float64_t[:L.shape[0], :4:1]> (
        <cnp.float64_t *> result_arr.data))
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
