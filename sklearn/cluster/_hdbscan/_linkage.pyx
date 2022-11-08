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

# Numpy structured dtype representing a single ordered edge in Prim's algorithm
MST_edge_dtype = np.dtype([
    ("current_node", np.intp),
    ("next_node", np.intp),
    ("distance", np.float64),
])

# Packed shouldn't make a difference since they're all 8-byte quantities,
# but it's included just to be safe.
ctypedef packed struct MST_edge_t:
    cnp.intp_t current_node
    cnp.intp_t next_node
    cnp.float64_t distance

cpdef cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst_from_mutual_reachability(
    cnp.ndarray[cnp.float64_t, ndim=2] mutual_reachability
):
    cdef:
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] node_labels
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] current_labels
        cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] min_reachability, left, right
        cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst

        cnp.ndarray[cnp.uint8_t, mode='c'] label_filter

        cnp.intp_t n_samples = mutual_reachability.shape[0]
        cnp.intp_t current_node, new_node_index, new_node, i

    mst = np.empty(n_samples - 1, dtype=MST_edge_dtype)
    node_labels = np.arange(n_samples, dtype=np.intp)
    current_node = 0
    min_reachability = np.infty * np.ones(n_samples)
    current_labels = node_labels
    for i in range(1, n_samples):
        label_filter = current_labels != current_node
        current_labels = current_labels[label_filter]
        left = min_reachability[label_filter]
        right = mutual_reachability[current_node][current_labels]
        min_reachability = np.minimum(left, right)

        new_node_index = np.argmin(min_reachability)
        new_node = current_labels[new_node_index]
        mst[i - 1].current_node = current_node
        mst[i - 1].next_node = new_node
        mst[i - 1].distance = min_reachability[new_node_index]
        current_node = new_node

    return mst


cpdef cnp.ndarray[MST_edge_t, ndim=1] mst_from_data_matrix(
    cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] raw_data,
    cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] core_distances,
    DistanceMetric dist_metric,
    cnp.float64_t alpha=1.0
):

    cdef:
        cnp.int8_t[::1] in_tree
        cnp.float64_t[::1] min_reachability, current_sources
        cnp.float64_t[::1] current_core_distances = core_distances
        cnp.float64_t[:, ::1] raw_data_view = raw_data
        cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst
        cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] mst_arr

        cnp.ndarray[cnp.uint8_t, mode='c'] label_filter

        cnp.intp_t current_node, source_node, right_node, left_node, new_node
        cnp.intp_t i, j, n_samples, num_features

        cnp.float64_t current_node_core_dist, new_reachability, mutual_reachability_distance
        cnp.float64_t next_node_min_reach, pair_distance, next_node_core_dist

    n_samples = raw_data.shape[0]
    num_features = raw_data.shape[1]

    mst = np.empty(n_samples - 1, dtype=MST_edge_dtype)

    in_tree = np.zeros(n_samples, dtype=np.int8)
    min_reachability = np.infty * np.ones(n_samples)
    current_sources = np.ones(n_samples)

    current_node = 0

    for i in range(1, n_samples):

        in_tree[current_node] = 1

        current_node_core_dist = current_core_distances[current_node]

        new_reachability = DBL_MAX
        source_node = 0
        new_node = 0

        for j in range(n_samples):
            if in_tree[j]:
                continue

            next_node_min_reach = min_reachability[j]
            next_node_source = current_sources[j]

            pair_distance = dist_metric.dist(
                &raw_data_view[current_node, 0],
                &raw_data_view[j, 0],
                num_features
            )

            if alpha != 1.0:
                pair_distance /= alpha

            next_node_core_dist = core_distances[j]
            mutual_reachability_distance = max(
                current_node_core_dist,
                next_node_core_dist,
                pair_distance
            )
            if mutual_reachability_distance > next_node_min_reach:
                if next_node_min_reach < new_reachability:
                    new_reachability = next_node_min_reach
                    source_node = next_node_source
                    new_node = j
                continue

            if mutual_reachability_distance < next_node_min_reach:
                min_reachability[j] = mutual_reachability_distance
                current_sources[j] = current_node
                if mutual_reachability_distance < new_reachability:
                    new_reachability = mutual_reachability_distance
                    source_node = current_node
                    new_node = j
            else:
                if next_node_min_reach < new_reachability:
                    new_reachability = next_node_min_reach
                    source_node = next_node_source
                    new_node = j

        mst[i - 1].current_node = source_node
        mst[i - 1].next_node = new_node
        mst[i - 1].distance = new_reachability
        current_node = new_node

    return mst

@cython.wraparound(True)
cpdef cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] label(MST_edge_t[::1] mst):

    cdef:
        cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] single_linkage

        # Note mst.shape[0] is one fewer than the number of samples
        cnp.intp_t n_samples = mst.shape[0] + 1
        cnp.intp_t current_node_ancestor, next_node_ancestor
        cnp.intp_t current_node, next_node, index
        cnp.float64_t distance

    single_linkage = np.zeros((n_samples - 1, 4))
    U = UnionFind(n_samples)

    for i in range(n_samples - 1):

        current_node = mst[i].current_node
        next_node = mst[i].next_node
        distance = mst[i].distance

        current_node_ancestor, next_node_ancestor = (
            U.fast_find(current_node),
            U.fast_find(next_node)
        )

        single_linkage[i][0] = current_node_ancestor
        single_linkage[i][1] = next_node_ancestor
        single_linkage[i][2] = distance
        single_linkage[i][3] = U.size[current_node_ancestor] + U.size[next_node_ancestor]

        U.union(current_node_ancestor, next_node_ancestor)

    return single_linkage
