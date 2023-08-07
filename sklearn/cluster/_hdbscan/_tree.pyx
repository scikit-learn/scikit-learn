# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# Copyright (c) 2015, Leland McInnes
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


cimport numpy as cnp
from libc.math cimport isinf
import cython

import numpy as np

cdef extern from "numpy/arrayobject.h":
    intp_t * PyArray_SHAPE(cnp.PyArrayObject *)

cdef cnp.float64_t INFTY = np.inf
cdef cnp.intp_t NOISE = -1

HIERARCHY_dtype = np.dtype([
    ("left_node", np.intp),
    ("right_node", np.intp),
    ("value", np.float64),
    ("cluster_size", np.intp),
])

CONDENSED_dtype = np.dtype([
    ("parent", np.intp),
    ("child", np.intp),
    ("value", np.float64),
    ("cluster_size", np.intp),
])

cpdef tuple tree_to_labels(
    const HIERARCHY_t[::1] single_linkage_tree,
    cnp.intp_t min_cluster_size=10,
    cluster_selection_method="eom",
    bint allow_single_cluster=False,
    cnp.float64_t cluster_selection_epsilon=0.0,
    max_cluster_size=None,
):
    cdef:
        cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] labels
        cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] probabilities

    condensed_tree = _condense_tree(single_linkage_tree, min_cluster_size)
    labels, probabilities = _get_clusters(
        condensed_tree,
        _compute_stability(condensed_tree),
        cluster_selection_method,
        allow_single_cluster,
        cluster_selection_epsilon,
        max_cluster_size,
    )

    return (labels, probabilities)

cdef list bfs_from_hierarchy(
    const HIERARCHY_t[::1] hierarchy,
    cnp.intp_t bfs_root
):
    """
    Perform a breadth first search on a tree in scipy hclust format.
    """

    cdef list process_queue, next_queue, result
    cdef cnp.intp_t n_samples = hierarchy.shape[0] + 1
    cdef cnp.intp_t node
    process_queue = [bfs_root]
    result = []

    while process_queue:
        result.extend(process_queue)
        # By construction, node i is formed by the union of nodes
        # hierarchy[i - n_samples, 0] and hierarchy[i - n_samples, 1]
        process_queue = [
            x - n_samples
            for x in process_queue
            if x >= n_samples
        ]
        if process_queue:
            next_queue = []
            for node in process_queue:
                next_queue.extend(
                    [
                        hierarchy[node].left_node,
                        hierarchy[node].right_node,
                    ]
                )
            process_queue = next_queue
    return result


cpdef cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] _condense_tree(
    const HIERARCHY_t[::1] hierarchy,
    cnp.intp_t min_cluster_size=10
):
    """Condense a tree according to a minimum cluster size. This is akin
    to the runt pruning procedure of Stuetzle. The result is a much simpler
    tree that is easier to visualize. We include extra information on the
    lambda value at which individual points depart clusters for later
    analysis and computation.

    Parameters
    ----------
    hierarchy : ndarray of shape (n_samples,), dtype=HIERARCHY_dtype
        A single linkage hierarchy in scipy.cluster.hierarchy format.

    min_cluster_size : int, optional (default 10)
        The minimum size of clusters to consider. Clusters smaller than this
        are pruned from the tree.

    Returns
    -------
    condensed_tree : ndarray of shape (n_samples,), dtype=CONDENSED_dtype
        Effectively an edgelist encoding a parent/child pair, along with a
        value and the corresponding cluster_size in each row providing a tree
        structure.
    """

    cdef:
        cnp.intp_t root = 2 * hierarchy.shape[0]
        cnp.intp_t n_samples = hierarchy.shape[0] + 1
        cnp.intp_t next_label = n_samples + 1
        list result_list, node_list = bfs_from_hierarchy(hierarchy, root)

        cnp.intp_t[::1] relabel
        cnp.uint8_t[::1] ignore

        cnp.intp_t node, sub_node, left, right
        cnp.float64_t lambda_value, distance
        cnp.intp_t left_count, right_count
        HIERARCHY_t children

    relabel = np.empty(root + 1, dtype=np.intp)
    relabel[root] = n_samples
    result_list = []
    ignore = np.zeros(len(node_list), dtype=bool)

    for node in node_list:
        if ignore[node] or node < n_samples:
            continue

        children = hierarchy[node - n_samples]
        left = children.left_node
        right = children.right_node
        distance = children.value
        if distance > 0.0:
            lambda_value = 1.0 / distance
        else:
            lambda_value = INFTY

        if left >= n_samples:
            left_count = hierarchy[left - n_samples].cluster_size
        else:
            left_count = 1

        if right >= n_samples:
            right_count = <cnp.intp_t> hierarchy[right - n_samples].cluster_size
        else:
            right_count = 1

        if left_count >= min_cluster_size and right_count >= min_cluster_size:
            relabel[left] = next_label
            next_label += 1
            result_list.append(
                (relabel[node], relabel[left], lambda_value, left_count)
            )

            relabel[right] = next_label
            next_label += 1
            result_list.append(
                (relabel[node], relabel[right], lambda_value, right_count)
            )

        elif left_count < min_cluster_size and right_count < min_cluster_size:
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < n_samples:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < n_samples:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

        elif left_count < min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < n_samples:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

        else:
            relabel[left] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < n_samples:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

    return np.array(result_list, dtype=CONDENSED_dtype)


cdef dict _compute_stability(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree
):

    cdef:
        cnp.float64_t[::1] result, births
        cnp.intp_t[:] parents = condensed_tree['parent']

        cnp.intp_t parent, cluster_size, result_index, idx
        cnp.float64_t lambda_val
        CONDENSED_t condensed_node
        cnp.intp_t largest_child = condensed_tree['child'].max()
        cnp.intp_t smallest_cluster = np.min(parents)
        cnp.intp_t num_clusters = np.max(parents) - smallest_cluster + 1
        dict stability_dict = {}

    largest_child = max(largest_child, smallest_cluster)
    births = np.full(largest_child + 1, np.nan, dtype=np.float64)

    for idx in range(PyArray_SHAPE(<cnp.PyArrayObject*> condensed_tree)[0]):
        condensed_node = condensed_tree[idx]
        births[condensed_node.child] = condensed_node.value

    births[smallest_cluster] = 0.0

    result = np.zeros(num_clusters, dtype=np.float64)
    for idx in range(PyArray_SHAPE(<cnp.PyArrayObject*> condensed_tree)[0]):
        condensed_node = condensed_tree[idx]
        parent = condensed_node.parent
        lambda_val = condensed_node.value
        cluster_size = condensed_node.cluster_size

        result_index = parent - smallest_cluster
        result[result_index] += (lambda_val - births[parent]) * cluster_size

    for idx in range(num_clusters):
        stability_dict[idx + smallest_cluster] = result[idx]

    return stability_dict


cdef list bfs_from_cluster_tree(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree,
    cnp.intp_t bfs_root
):

    cdef:
        list result = []
        cnp.ndarray[cnp.intp_t, ndim=1] process_queue = (
            np.array([bfs_root], dtype=np.intp)
        )
        cnp.ndarray[cnp.intp_t, ndim=1] children = condensed_tree['child']
        cnp.intp_t[:] parents = condensed_tree['parent']

    while len(process_queue) > 0:
        result.extend(process_queue.tolist())
        process_queue = children[np.isin(parents, process_queue)]

    return result


cdef cnp.float64_t[::1] max_lambdas(cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree):

    cdef:
        cnp.intp_t parent, current_parent, idx
        cnp.float64_t lambda_val, max_lambda
        cnp.float64_t[::1] deaths
        cnp.intp_t largest_parent = condensed_tree['parent'].max()

    deaths = np.zeros(largest_parent + 1, dtype=np.float64)
    current_parent = condensed_tree[0].parent
    max_lambda = condensed_tree[0].value

    for idx in range(1, PyArray_SHAPE(<cnp.PyArrayObject*> condensed_tree)[0]):
        parent = condensed_tree[idx].parent
        lambda_val = condensed_tree[idx].value

        if parent == current_parent:
            max_lambda = max(max_lambda, lambda_val)
        else:
            deaths[current_parent] = max_lambda
            current_parent = parent
            max_lambda = lambda_val

    deaths[current_parent] = max_lambda  # value for last parent
    return deaths


@cython.final
cdef class TreeUnionFind:

    cdef cnp.intp_t[:, ::1] data
    cdef cnp.uint8_t[::1] is_component

    def __init__(self, size):
        cdef cnp.intp_t idx
        self.data = np.zeros((size, 2), dtype=np.intp)
        for idx in range(size):
            self.data[idx, 0] = idx
        self.is_component = np.ones(size, dtype=np.uint8)

    cdef void union(self, cnp.intp_t x, cnp.intp_t y):
        cdef cnp.intp_t x_root = self.find(x)
        cdef cnp.intp_t y_root = self.find(y)

        if self.data[x_root, 1] < self.data[y_root, 1]:
            self.data[x_root, 0] = y_root
        elif self.data[x_root, 1] > self.data[y_root, 1]:
            self.data[y_root, 0] = x_root
        else:
            self.data[y_root, 0] = x_root
            self.data[x_root, 1] += 1
        return

    cdef cnp.intp_t find(self, cnp.intp_t x):
        if self.data[x, 0] != x:
            self.data[x, 0] = self.find(self.data[x, 0])
            self.is_component[x] = False
        return self.data[x, 0]


cpdef cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] labelling_at_cut(
        const HIERARCHY_t[::1] linkage,
        cnp.float64_t cut,
        cnp.intp_t min_cluster_size
):
    """Given a single linkage tree and a cut value, return the
    vector of cluster labels at that cut value. This is useful
    for Robust Single Linkage, and extracting DBSCAN results
    from a single HDBSCAN run.

    Parameters
    ----------
    linkage : ndarray of shape (n_samples,), dtype=HIERARCHY_dtype
        The single linkage tree in scipy.cluster.hierarchy format.

    cut : double
        The cut value at which to find clusters.

    min_cluster_size : int
        The minimum cluster size; clusters below this size at
        the cut will be considered noise.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        The cluster labels for each point in the data set;
        a label of -1 denotes a noise assignment.
    """

    cdef:
        cnp.intp_t n, cluster, root, n_samples, cluster_label
        cnp.intp_t[::1] unique_labels, cluster_size
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] result
        TreeUnionFind union_find
        dict cluster_label_map
        HIERARCHY_t node

    root = 2 * linkage.shape[0]
    n_samples = root // 2 + 1
    result = np.empty(n_samples, dtype=np.intp)
    union_find = TreeUnionFind(root + 1)

    cluster = n_samples
    for node in linkage:
        if node.value < cut:
            union_find.union(node.left_node, cluster)
            union_find.union(node.right_node, cluster)
        cluster += 1

    cluster_size = np.zeros(cluster, dtype=np.intp)
    for n in range(n_samples):
        cluster = union_find.find(n)
        cluster_size[cluster] += 1
        result[n] = cluster

    cluster_label_map = {-1: NOISE}
    cluster_label = 0
    unique_labels = np.unique(result)

    for cluster in unique_labels:
        if cluster_size[cluster] < min_cluster_size:
            cluster_label_map[cluster] = NOISE
        else:
            cluster_label_map[cluster] = cluster_label
            cluster_label += 1

    for n in range(n_samples):
        result[n] = cluster_label_map[result[n]]

    return result


cpdef cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] _do_labelling(
        cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree,
        set clusters,
        dict cluster_label_map,
        cnp.intp_t allow_single_cluster,
        cnp.float64_t cluster_selection_epsilon
):
    """Given a condensed tree, clusters and a labeling map for the clusters,
    return an array containing the labels of each point based on cluster
    membership. Note that this is where points may be marked as noisy
    outliers. The determination of some points as noise is in large, single-
    cluster datasets is controlled by the `allow_single_cluster` and
    `cluster_selection_epsilon` parameters.

    Parameters
    ----------
    condensed_tree : ndarray of shape (n_samples,), dtype=CONDENSED_dtype
        Effectively an edgelist encoding a parent/child pair, along with a
        value and the corresponding cluster_size in each row providing a tree
        structure.

    clusters : set
        The set of nodes corresponding to identified clusters. These node
        values should be the same as those present in `condensed_tree`.

    cluster_label_map : dict
        A mapping from the node values present in `clusters` to the labels
        which will be returned.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        The cluster labels for each point in the data set;
        a label of -1 denotes a noise assignment.
    """

    cdef:
        cnp.intp_t root_cluster
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] result
        cnp.ndarray[cnp.intp_t, ndim=1] parent_array, child_array
        cnp.ndarray[cnp.float64_t, ndim=1] lambda_array
        TreeUnionFind union_find
        cnp.intp_t n, parent, child, cluster
        cnp.float64_t threshold

    child_array = condensed_tree['child']
    parent_array = condensed_tree['parent']
    lambda_array = condensed_tree['value']

    root_cluster = np.min(parent_array)
    result = np.empty(root_cluster, dtype=np.intp)
    union_find = TreeUnionFind(np.max(parent_array) + 1)

    for n in range(PyArray_SHAPE(<cnp.PyArrayObject*> condensed_tree)[0]):
        child = child_array[n]
        parent = parent_array[n]
        if child not in clusters:
            union_find.union(parent, child)

    for n in range(root_cluster):
        cluster = union_find.find(n)
        label = NOISE
        if cluster != root_cluster:
            label = cluster_label_map[cluster]
        elif len(clusters) == 1 and allow_single_cluster:
            # There can only be one edge with this particular child hence this
            # expression extracts a unique, scalar lambda value.
            parent_lambda = lambda_array[child_array == n]
            if cluster_selection_epsilon != 0.0:
                threshold = 1 / cluster_selection_epsilon
            else:
                # The threshold should be calculated per-sample based on the
                # largest lambda of any simbling node.
                threshold = lambda_array[parent_array == cluster].max()
            if parent_lambda >= threshold:
                label = cluster_label_map[cluster]

        result[n] = label

    return result


cdef cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] get_probabilities(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree,
    dict cluster_map,
    cnp.intp_t[::1] labels
):

    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] result
        cnp.float64_t[:] lambda_array
        cnp.float64_t[::1] deaths
        cnp.intp_t[:] child_array, parent_array
        cnp.intp_t root_cluster, n, point, cluster_num, cluster
        cnp.float64_t max_lambda, lambda_val

    child_array = condensed_tree['child']
    parent_array = condensed_tree['parent']
    lambda_array = condensed_tree['value']

    result = np.zeros(labels.shape[0])
    deaths = max_lambdas(condensed_tree)
    root_cluster = np.min(parent_array)

    for n in range(PyArray_SHAPE(<cnp.PyArrayObject*> condensed_tree)[0]):
        point = child_array[n]
        if point >= root_cluster:
            continue

        cluster_num = labels[point]
        if cluster_num == -1:
            continue

        cluster = cluster_map[cluster_num]
        max_lambda = deaths[cluster]
        if max_lambda == 0.0 or isinf(lambda_array[n]):
            result[point] = 1.0
        else:
            lambda_val = min(lambda_array[n], max_lambda)
            result[point] = lambda_val / max_lambda

    return result


cpdef list recurse_leaf_dfs(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] cluster_tree,
    cnp.intp_t current_node
):
    cdef cnp.intp_t[:] children
    cdef cnp.intp_t child

    children = cluster_tree[cluster_tree['parent'] == current_node]['child']
    if children.shape[0] == 0:
        return [current_node,]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


cpdef list get_cluster_tree_leaves(cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] cluster_tree):
    cdef cnp.intp_t root
    if PyArray_SHAPE(<cnp.PyArrayObject*> cluster_tree)[0] == 0:
        return []
    root = cluster_tree['parent'].min()
    return recurse_leaf_dfs(cluster_tree, root)

cdef cnp.intp_t traverse_upwards(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] cluster_tree,
    cnp.float64_t cluster_selection_epsilon,
    cnp.intp_t leaf,
    cnp.intp_t allow_single_cluster
):
    cdef cnp.intp_t root, parent
    cdef cnp.float64_t parent_eps

    root = cluster_tree['parent'].min()
    parent = cluster_tree[cluster_tree['child'] == leaf]['parent']
    if parent == root:
        if allow_single_cluster:
            return parent
        else:
            return leaf  # return node closest to root

    parent_eps = 1 / cluster_tree[cluster_tree['child'] == parent]['value']
    if parent_eps > cluster_selection_epsilon:
        return parent
    else:
        return traverse_upwards(
            cluster_tree,
            cluster_selection_epsilon,
            parent,
            allow_single_cluster
        )

cdef set epsilon_search(
    set leaves,
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] cluster_tree,
    cnp.float64_t cluster_selection_epsilon,
    cnp.intp_t allow_single_cluster
):
    cdef:
        list selected_clusters = list()
        list processed = list()
        cnp.intp_t leaf, epsilon_child, sub_node
        cnp.float64_t eps
        cnp.uint8_t[:] leaf_nodes
        cnp.ndarray[cnp.intp_t, ndim=1] children = cluster_tree['child']
        cnp.ndarray[cnp.float64_t, ndim=1] distances = cluster_tree['value']

    for leaf in leaves:
        leaf_nodes = children == leaf
        eps = 1 / distances[leaf_nodes][0]
        if eps < cluster_selection_epsilon:
            if leaf not in processed:
                epsilon_child = traverse_upwards(
                    cluster_tree,
                    cluster_selection_epsilon,
                    leaf,
                    allow_single_cluster
                )
                selected_clusters.append(epsilon_child)

                for sub_node in bfs_from_cluster_tree(cluster_tree, epsilon_child):
                    if sub_node != epsilon_child:
                        processed.append(sub_node)
        else:
            selected_clusters.append(leaf)

    return set(selected_clusters)


@cython.wraparound(True)
cdef tuple _get_clusters(
    cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] condensed_tree,
    dict stability,
    cluster_selection_method='eom',
    cnp.uint8_t allow_single_cluster=False,
    cnp.float64_t cluster_selection_epsilon=0.0,
    max_cluster_size=None
):
    """Given a tree and stability dict, produce the cluster labels
    (and probabilities) for a flat clustering based on the chosen
    cluster selection method.

    Parameters
    ----------
    condensed_tree : ndarray of shape (n_samples,), dtype=CONDENSED_dtype
        Effectively an edgelist encoding a parent/child pair, along with a
        value and the corresponding cluster_size in each row providing a tree
        structure.

    stability : dict
        A dictionary mapping cluster_ids to stability values

    cluster_selection_method : string, optional (default 'eom')
        The method of selecting clusters. The default is the
        Excess of Mass algorithm specified by 'eom'. The alternate
        option is 'leaf'.

    allow_single_cluster : boolean, optional (default False)
        Whether to allow a single cluster to be selected by the
        Excess of Mass algorithm.

    cluster_selection_epsilon: double, optional (default 0.0)
        A distance threshold for cluster splits.

    max_cluster_size: int, default=None
        The maximum size for clusters located by the EOM clusterer. Can
        be overridden by the cluster_selection_epsilon parameter in
        rare cases.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        An integer array of cluster labels, with -1 denoting noise.

    probabilities : ndarray (n_samples,)
        The cluster membership strength of each sample.

    stabilities : ndarray (n_clusters,)
        The cluster coherence strengths of each cluster.
    """
    cdef:
        list node_list
        cnp.ndarray[CONDENSED_t, ndim=1, mode='c'] cluster_tree
        cnp.uint8_t[::1] child_selection
        cnp.ndarray[cnp.intp_t, ndim=1, mode='c'] labels
        dict is_cluster, cluster_sizes
        cnp.float64_t subtree_stability
        cnp.intp_t node, sub_node, cluster, n_samples
        cnp.ndarray[cnp.float64_t, ndim=1, mode='c'] probs

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        node_list = sorted(stability.keys(), reverse=True)[:-1]
        # (exclude root)

    cluster_tree = condensed_tree[condensed_tree['cluster_size'] > 1]
    is_cluster = {cluster: True for cluster in node_list}
    n_samples = np.max(condensed_tree[condensed_tree['cluster_size'] == 1]['child']) + 1

    if max_cluster_size is None:
        max_cluster_size = n_samples + 1  # Set to a value that will never be triggered
    cluster_sizes = {
        child: cluster_size for child, cluster_size
        in zip(cluster_tree['child'], cluster_tree['cluster_size'])
    }
    if allow_single_cluster:
        # Compute cluster size for the root node
        cluster_sizes[node_list[-1]] = np.sum(
            cluster_tree[cluster_tree['parent'] == node_list[-1]]['cluster_size'])

    if cluster_selection_method == 'eom':
        for node in node_list:
            child_selection = (cluster_tree['parent'] == node)
            subtree_stability = np.sum([
                stability[child] for
                child in cluster_tree['child'][child_selection]])
            if subtree_stability > stability[node] or cluster_sizes[node] > max_cluster_size:
                is_cluster[node] = False
                stability[node] = subtree_stability
            else:
                for sub_node in bfs_from_cluster_tree(cluster_tree, node):
                    if sub_node != node:
                        is_cluster[sub_node] = False

        if cluster_selection_epsilon != 0.0 and PyArray_SHAPE(<cnp.PyArrayObject*> cluster_tree)[0] > 0:
            eom_clusters = [c for c in is_cluster if is_cluster[c]]
            selected_clusters = []
            # first check if eom_clusters only has root node, which skips epsilon check.
            if (len(eom_clusters) == 1 and eom_clusters[0] == cluster_tree['parent'].min()):
                if allow_single_cluster:
                    selected_clusters = eom_clusters
            else:
                selected_clusters = epsilon_search(
                    set(eom_clusters),
                    cluster_tree,
                    cluster_selection_epsilon,
                    allow_single_cluster
                )
            for c in is_cluster:
                if c in selected_clusters:
                    is_cluster[c] = True
                else:
                    is_cluster[c] = False

    elif cluster_selection_method == 'leaf':
        leaves = set(get_cluster_tree_leaves(cluster_tree))
        if len(leaves) == 0:
            for c in is_cluster:
                is_cluster[c] = False
            is_cluster[condensed_tree['parent'].min()] = True

        if cluster_selection_epsilon != 0.0:
            selected_clusters = epsilon_search(
                leaves,
                cluster_tree,
                cluster_selection_epsilon,
                allow_single_cluster
            )
        else:
            selected_clusters = leaves

        for c in is_cluster:
            if c in selected_clusters:
                is_cluster[c] = True
            else:
                is_cluster[c] = False

    clusters = set([c for c in is_cluster if is_cluster[c]])
    cluster_map = {c: n for n, c in enumerate(sorted(list(clusters)))}
    reverse_cluster_map = {n: c for c, n in cluster_map.items()}

    labels = _do_labelling(
        condensed_tree,
        clusters,
        cluster_map,
        allow_single_cluster,
        cluster_selection_epsilon
    )
    probs = get_probabilities(condensed_tree, reverse_cluster_map, labels)

    return (labels, probs)
