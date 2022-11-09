# Minimum spanning tree single linkage implementation for hdbscan
# Authors: Leland McInnes <leland.mcinnes@gmail.com>
#          Steve Astels <sastels@gmail.com>
#          Meekail Zain <zainmeekail@gmail.com>
# License: 3-clause BSD

cimport numpy as cnp
from libc.float cimport DBL_MAX

import numpy as np
from ...metrics._dist_metrics cimport DistanceMetric
from ...cluster._hierarchical_fast cimport UnionFind
from ...utils._typedefs cimport ITYPE_t, DTYPE_t
from ...utils._typedefs import ITYPE, DTYPE

# Numpy structured dtype representing a single ordered edge in Prim's algorithm
MST_edge_dtype = np.dtype([
    ("current_node", np.int64),
    ("next_node", np.int64),
    ("distance", np.float64),
])

# Packed shouldn't make a difference since they're all 8-byte quantities,
# but it's included just to be safe.
ctypedef packed struct MST_edge_t:
    cnp.int64_t current_node
    cnp.int64_t next_node
    cnp.float64_t distance

cpdef cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst_from_mutual_reachability(
    cnp.ndarray[cnp.float64_t, ndim=2] mutual_reachability
):
    """Compute the Minimum Spanning Tree (MST) representation of the mutual-
    reachability graph using Prim's algorithm.

    Parameters
    ----------
    mutual_reachability : ndarray of shape (n_samples, n_samples)
        Array of mutual-reachabilities between samples.

    Returns
    -------
    mst : ndarray of shape (n_samples - 1,)
        The MST representation of the mutual-reahability graph. The MST is
        represented as a collecteion of edges. Each edge is an instance of a
        custom dtype `MST_edge_dtype` with the following specification:

        MST_edge_dtype = np.dtype([
            ("current_node", np.int64),
            ("next_node", np.int64),
            ("distance", np.float64),
        ])
    """
    cdef:
        cnp.ndarray[cnp.int64_t, ndim=1, mode='c'] current_labels
        cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] min_reachability, left, right
        cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst

        cnp.ndarray[cnp.uint8_t, mode='c'] label_filter

        cnp.int64_t n_samples = mutual_reachability.shape[0]
        cnp.int64_t current_node, new_node_index, new_node, i

    mst = np.empty(n_samples - 1, dtype=MST_edge_dtype)
    current_labels = np.arange(n_samples, dtype=np.int64)
    current_node = 0
    min_reachability = np.infty * np.ones(n_samples, dtype=np.float64)
    for i in range(0, n_samples - 1):
        label_filter = current_labels != current_node
        current_labels = current_labels[label_filter]
        left = min_reachability[label_filter]
        right = mutual_reachability[current_node][current_labels]
        min_reachability = np.minimum(left, right)

        new_node_index = np.argmin(min_reachability)
        new_node = current_labels[new_node_index]
        mst[i].current_node = current_node
        mst[i].next_node = new_node
        mst[i].distance = min_reachability[new_node_index]
        current_node = new_node

    return mst


cpdef cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst_from_data_matrix(
    const cnp.float64_t[:, ::1] raw_data,
    const cnp.float64_t[::1] core_distances,
    DistanceMetric dist_metric,
):
    """Compute the Minimum Spanning Tree (MST) representation of the mutual-
    reachability graph generated from the provided `raw_data` and
    `core_distances` using Prim's algorithm.

    Parameters
    ----------
    raw_data : ndarray of shape (n_samples, n_features)
        Input array of data samples.

    core_distances : ndarray of shape (n_samples,)
        An array containing the core-distance calculated for each corresponding
        sample.

    dist_metric : DistanceMetric
        The distance metric to use when calculating pairwise distances for
        determining mutual-reachability.

    Returns
    -------
    mst : ndarray of shape (n_samples - 1,)
        The MST representation of the mutual-reahability graph. The MST is
        represented as a collecteion of edges. Each edge is an instance of a
        custom dtype `MST_edge_dtype` with the following specification:

        MST_edge_dtype = np.dtype([
            ("current_node", np.int64),
            ("next_node", np.int64),
            ("distance", np.float64),
        ])
    """

    cdef:
        cnp.int8_t[::1] in_tree
        cnp.float64_t[::1] min_reachability
        cnp.int64_t[::1] current_sources
        cnp.ndarray[MST_edge_t, ndim=1, mode='c'] mst

        cnp.int64_t current_node, source_node, right_node, left_node, new_node, next_node_source
        cnp.int64_t i, j, n_samples, num_features

        cnp.float64_t current_node_core_dist, new_reachability, mutual_reachability_distance
        cnp.float64_t next_node_min_reach, pair_distance, next_node_core_dist

    n_samples = raw_data.shape[0]
    num_features = raw_data.shape[1]

    mst = np.empty(n_samples - 1, dtype=MST_edge_dtype)

    in_tree = np.zeros(n_samples, dtype=np.int8)
    min_reachability = np.infty * np.ones(n_samples, dtype=np.float64)
    current_sources = np.ones(n_samples, dtype=np.int64)

    current_node = 0

    for i in range(0, n_samples - 1):

        in_tree[current_node] = 1

        current_node_core_dist = core_distances[current_node]

        new_reachability = DBL_MAX
        source_node = 0
        new_node = 0

        for j in range(n_samples):
            if in_tree[j]:
                continue

            next_node_min_reach = min_reachability[j]
            next_node_source = current_sources[j]

            pair_distance = dist_metric.dist(
                &raw_data[current_node, 0],
                &raw_data[j, 0],
                num_features
            )

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

        mst[i].current_node = source_node
        mst[i].next_node = new_node
        mst[i].distance = new_reachability
        current_node = new_node

    return mst

cpdef cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] make_single_linkage(const MST_edge_t[::1] mst):
    """Construct a single-linkage tree from an MST.

    Parameters
    ----------
    mst : ndarray of shape (n_samples - 1,)
        The MST representation of the mutual-reahability graph. The MST is
        represented as a collecteion of edges. Each edge is an instance of a
        custom dtype `MST_edge_dtype` with the following specification:

        MST_edge_dtype = np.dtype([
            ("current_node", np.int64),
            ("next_node", np.int64),
            ("distance", np.float64),
        ])

    Returns
    -------
    single_linkage : ndarray of shape (n_samples - 1, 4)
        The single-linkage tree tree (dendrogram) built from the MST. Each
        of the array represents the following:

        - left node/cluster
        - right node/cluster
        - distance
        - new cluster size
    """
    cdef:
        cnp.ndarray[cnp.float64_t, ndim=2, mode='c'] single_linkage

        # Note mst.shape[0] is one fewer than the number of samples
        cnp.int64_t n_samples = mst.shape[0] + 1
        cnp.int64_t current_node_cluster, next_node_cluster
        cnp.int64_t current_node, next_node, index
        cnp.float64_t distance
        UnionFind U = UnionFind(n_samples)

    single_linkage = np.zeros((n_samples - 1, 4))

    for i in range(n_samples - 1):

        current_node = mst[i].current_node
        next_node = mst[i].next_node
        distance = mst[i].distance

        current_node_cluster, next_node_cluster = (
            U.fast_find(current_node),
            U.fast_find(next_node)
        )

        # TODO: Update this to an array of structs (AoS).
        # Should be done simultaneously in _tree.pyx to ensure compatability.
        single_linkage[i][0] = <cnp.float64_t> current_node_cluster
        single_linkage[i][1] = <cnp.float64_t> next_node_cluster
        single_linkage[i][2] = distance
        single_linkage[i][3] = U.size[current_node_cluster] + U.size[next_node_cluster]

        U.union(current_node_cluster, next_node_cluster)

    return single_linkage
