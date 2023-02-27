# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD


cimport numpy as cnp

import cython
import numpy as np

cdef cnp.float64_t INFTY = np.inf
cdef cnp.intp_t NOISE = -1


HIERARCHY_dtype = np.dtype([
    ("left_node", np.intp),
    ("right_node", np.intp),
    ("value", np.float64),
    ("cluster_size", np.intp),
])

cdef list bfs_from_hierarchy(
    cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy,
    cnp.intp_t bfs_root
):
    """
    Perform a breadth first search on a single linkage hierarchy in
    scipy.cluster.hierarchy format.
    """
    # NOTE: We keep `process_queue` as a list rather than a memory-view to
    # retain semantics which make the below runtime algorithm convenient, such
    # as not being required to preallocate space, and allowing easy checks to
    # determine whether the list is emptey or not.
    cdef list process_queue, next_queue
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

cpdef cnp.ndarray condense_tree(
    cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy,
    cnp.intp_t min_cluster_size=10
):
    """Condense an MST according to a minimum cluster size. This is akin
    to the runt pruning procedure of Stuetzle. The result is a much simpler
    tree that is easier to visualize. We include extra information on the
    lambda value at which individual points depart clusters for later
    analysis and computation.

    Parameters
    ----------
    hierarchy : ndarray of shape (n_samples,), dtype=HIERARCHY_dtype
        A single linkage hierarchy in scipy.cluster.hierarchy format.

    min_cluster_size : int, optional (default 10)
        The minimum size of clusters to consider. Clusters smaler than this
        are pruned from the tree.

    Returns
    -------
    condensed_tree : ndarray of shape (n_samples,), dtype=HIERARCHY_dtype
        Effectively an edgelist with a parent, child, lambda_val
        and cluster_size in each row providing a tree structure.
    """

    cdef:
        cnp.intp_t root = 2 * hierarchy.shape[0]
        cnp.intp_t num_points = hierarchy.shape[0] + 1
        cnp.intp_t next_label = num_points + 1
        list result_list, node_list = bfs_from_hierarchy(hierarchy, root)

        cnp.intp_t[::1] relabel
        cnp.uint8_t[::1] ignore
        HIERARCHY_t children

        cnp.intp_t node, sub_node, left, right
        cnp.float64_t lambda_value, distance
        cnp.intp_t left_count, right_count

    relabel = np.empty(root + 1, dtype=np.intp)
    relabel[root] = num_points
    result_list = []
    ignore = np.zeros(len(node_list), dtype=np.uint8)

    for node in node_list:
        # Guarantee that node is a cluster node of interest
        if ignore[node] or node < num_points:
            continue

        children = hierarchy[node - num_points]
        left = children.left_node
        right = children.right_node
        distance = children.value
        if distance > 0.0:
            lambda_value = 1.0 / distance
        else:
            lambda_value = INFTY

        # Guarantee that left is a cluster node
        if left >= num_points:
            left_count = hierarchy[left - num_points].cluster_size
        else:
            left_count = 1

        # Guarantee that right is a cluster node
        if right >= num_points:
            right_count = hierarchy[right - num_points].cluster_size
        else:
            right_count = 1

        # Each child is regarded as a proper cluster
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

        # Each child is regarded as a collection of single-sample clusters
        elif left_count < min_cluster_size and right_count < min_cluster_size:
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

        # One child is a collection of single-sample clusters, while the other
        # is a persistance of the parent node cluster
        elif left_count < min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

        # One child is a collection of single-sample clusters, while the other
        # is a persistance of the parent node cluster
        else:
            relabel[left] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append(
                        (relabel[node], sub_node, lambda_value, 1)
                    )
                ignore[sub_node] = True

    return np.array(result_list, dtype=HIERARCHY_dtype)


cpdef dict compute_stability(cnp.ndarray[HIERARCHY_t, ndim=1] condensed_tree):

    cdef:
        cnp.float64_t[::1] result, births_arr
        cnp.ndarray[cnp.intp_t, ndim=1] parents

        cnp.intp_t parent, cluster_size, result_index
        cnp.float64_t lambda_val
        HIERARCHY_t condensed_node
        cnp.float64_t[:, :] result_pre_dict


    parents = condensed_tree['left_node']
    cdef cnp.intp_t largest_child = condensed_tree['right_node'].max()
    cdef cnp.intp_t smallest_cluster = parents.min()
    cdef cnp.intp_t num_clusters = parents.max() - smallest_cluster + 1

    largest_child = max(largest_child, smallest_cluster)

    births_arr = np.full(largest_child + 1, np.nan, dtype=np.float64)

    for idx in range(condensed_tree.shape[0]):
        condensed_node = condensed_tree[idx]
        births_arr[condensed_node.right_node] = condensed_node.value

    births_arr[smallest_cluster] = 0.0

    result = np.zeros(num_clusters, dtype=np.float64)
    for condensed_node in condensed_tree:
        parent = condensed_node.left_node
        lambda_val = condensed_node.value
        cluster_size = condensed_node.cluster_size

        result_index = parent - smallest_cluster
        result[result_index] += (lambda_val - births_arr[parent]) * cluster_size

    result_pre_dict = np.vstack(
        (
            np.arange(smallest_cluster, parents.max() + 1),
            result
        )
    ).T

    return dict(result_pre_dict)


cdef list bfs_from_cluster_tree(
    cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy,
    cnp.intp_t bfs_root,
):

    cdef list result
    cdef cnp.ndarray[cnp.intp_t, ndim=1] process_queue, children = hierarchy['right_node']
    cdef cnp.intp_t[:] parents = hierarchy['left_node']

    result = []
    process_queue = np.array([bfs_root], dtype=np.intp)

    while process_queue.shape[0] > 0:
        result.extend(process_queue.tolist())
        process_queue = children[np.isin(parents, process_queue)]

    return result


cdef cnp.float64_t[::1] max_lambdas(cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy):

    cdef cnp.intp_t parent, current_parent, idx
    cdef cnp.float64_t lambda_val, max_lambda

    cdef cnp.float64_t[::1] deaths

    cdef cnp.intp_t largest_parent = hierarchy['left_node'].max()

    deaths = np.zeros(largest_parent + 1, dtype=np.float64)

    current_parent = hierarchy[0].left_node
    max_lambda = hierarchy[0].value

    for idx in range(1, hierarchy.shape[0]):
        parent = hierarchy[idx].left_node
        lambda_val = hierarchy[idx].value

        if parent == current_parent:
            max_lambda = max(max_lambda, lambda_val)
        else:
            deaths[current_parent] = max_lambda
            current_parent = parent
            max_lambda = lambda_val

    deaths[current_parent] = max_lambda # value for last parent
    return deaths


cdef class TreeUnionFind:

    cdef cnp.intp_t[:, ::1] data
    cdef cnp.uint8_t[::1] is_component

    def __init__(self, size):
        cdef cnp.ndarray[cnp.intp_t, ndim=2] data_arr
        data_arr = np.zeros((size, 2), dtype=np.intp)
        data_arr.T[0] = np.arange(size)
        self.data = data_arr
        self.is_component = np.ones(size, dtype=bool)

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

cpdef cnp.ndarray[cnp.intp_t, ndim=1] labelling_at_cut(
        cnp.ndarray[HIERARCHY_t, ndim=1] linkage,
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

    cut : float
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
        cnp.intp_t root, num_points
        cnp.intp_t[::1] unique_labels, cluster_size
        cnp.ndarray[cnp.intp_t, ndim=1] result
        TreeUnionFind union_find
        cnp.intp_t n, cluster, cluster_id, cluster_label
        dict cluster_label_map
        HIERARCHY_t node

    root = 2 * linkage.shape[0]
    num_points = root // 2 + 1

    result = np.empty(num_points, dtype=np.intp)

    union_find = TreeUnionFind(root + 1)

    cluster = num_points
    for node in linkage:
        if node.value < cut:
            union_find.union(node.left_node, cluster)
            union_find.union(node.right_node, cluster)
        cluster += 1

    cluster_size = np.zeros(cluster, dtype=np.intp)
    for n in range(num_points):
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

    for n in range(num_points):
        result[n] = cluster_label_map[result[n]]

    return result


cdef cnp.ndarray[cnp.intp_t, ndim=1] do_labelling(
        cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy,
        set clusters,
        dict cluster_label_map,
        cnp.uint8_t allow_single_cluster,
        cnp.float64_t cluster_selection_epsilon):
    cdef:
        cnp.ndarray[cnp.intp_t, ndim=1] result, parents, children
        cnp.ndarray[cnp.float64_t, ndim=1] lambdas
        TreeUnionFind union_find
        cnp.intp_t root_cluster, parent, child, cluster, n
        cnp.int64_t cluster_label, label

    children = hierarchy['right_node']
    parents = hierarchy['left_node']
    lambdas = hierarchy['value']

    root_cluster = parents.min()
    result = np.empty(root_cluster, dtype=np.intp)

    union_find = TreeUnionFind(parents.max() + 1)

    for n in range(hierarchy.shape[0]):
        child = children[n]
        parent = parents[n]
        if child not in clusters:
            union_find.union(parent, child)

    for n in range(root_cluster):
        cluster = union_find.find(n)
        label = NOISE
        if cluster != root_cluster:
            label = cluster_label_map[cluster]
            result[n] = label
            continue
        if len(clusters) == 1 and allow_single_cluster:
            if cluster_selection_epsilon != 0.0:
                if lambdas[children == n] >= 1 / cluster_selection_epsilon :
                    label = cluster_label_map[cluster]
            elif lambdas[children == n] >= lambdas[parents == cluster].max():
                label = cluster_label_map[cluster]
        result[n] = label
    return result


cdef cnp.ndarray[cnp.float64_t, ndim=1] get_probabilities(
    cnp.ndarray[HIERARCHY_t, ndim=1] hierarchy,
    dict cluster_map,
    cnp.intp_t[:] labels
):
    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1] result
        cnp.float64_t[:] lambdas
        cnp.float64_t[::1] deaths
        cnp.ndarray[cnp.intp_t, ndim=1] children, parents
        cnp.intp_t n, point, cluster_num, cluster, root_cluster
        cnp.float64_t max_lambda, lambda_val

    children = hierarchy['right_node']
    parents = hierarchy['left_node']
    lambdas = hierarchy['value']

    result = np.zeros(labels.shape[0])
    deaths = max_lambdas(hierarchy)
    root_cluster = parents.min()

    for n in range(hierarchy.shape[0]):
        point = children[n]
        if point >= root_cluster:
            continue

        cluster_num = labels[point]

        if cluster_num == -1:
            continue

        cluster = cluster_map[cluster_num]
        max_lambda = deaths[cluster]
        if max_lambda == 0.0 or not np.isfinite(lambdas[n]):
            result[point] = 1.0
        else:
            lambda_val = min(lambdas[n], max_lambda)
            result[point] = lambda_val / max_lambda

    return result


cdef list recurse_leaf_dfs(
    cnp.ndarray[HIERARCHY_t, ndim=1] cluster_tree,
    cnp.intp_t current_node,
):
    cdef cnp.intp_t[:] children
    cdef cnp.intp_t child

    children = cluster_tree[cluster_tree['left_node'] == current_node]['right_node']
    if children.shape[0] == 0:
        return [current_node,]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


cdef list get_cluster_tree_leaves(cnp.ndarray[HIERARCHY_t, ndim=1] cluster_tree):
    cdef cnp.intp_t root
    if cluster_tree.shape[0] == 0:
        return []
    root = cluster_tree['left_node'].min()
    return recurse_leaf_dfs(cluster_tree, root)

cdef cnp.intp_t traverse_upwards(
    cnp.ndarray[HIERARCHY_t, ndim=1] cluster_tree,
    cnp.float64_t cluster_selection_epsilon,
    cnp.intp_t leaf,
    cnp.intp_t allow_single_cluster
):

    cdef cnp.intp_t root, parent
    cdef cnp.float64_t parent_eps

    root = cluster_tree['left_node'].min()
    parent = cluster_tree[cluster_tree['right_node'] == leaf]['left_node']
    if parent == root:
        if allow_single_cluster:
            return parent
        else:
            return leaf #return node closest to root

    parent_eps = 1 / cluster_tree[cluster_tree['right_node'] == parent]['value']
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
    cnp.ndarray[HIERARCHY_t, ndim=1] cluster_tree,
    cnp.float64_t cluster_selection_epsilon,
    cnp.intp_t allow_single_cluster
):

    cdef list selected_clusters = list()
    cdef list processed = list()
    cdef cnp.intp_t leaf, epsilon_child, sub_node
    cdef cnp.float64_t eps

    for leaf in leaves:
        eps = 1 / cluster_tree['value'][cluster_tree['right_node'] == leaf][0]
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
cpdef tuple get_clusters(
    cnp.ndarray hierarchy,
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
    hierarchy : ndarray of shape (n_samples,), dtype=HIERARCHY_dtype
        The condensed tree to extract flat clusters from

    stability : dict
        A dictionary mapping cluster_ids to stability values

    cluster_selection_method : string, optional (default 'eom')
        The method of selecting clusters. The default is the
        Excess of Mass algorithm specified by 'eom'. The alternate
        option is 'leaf'.

    allow_single_cluster : boolean, optional (default False)
        Whether to allow a single cluster to be selected by the
        Excess of Mass algorithm.

    cluster_selection_epsilon: float, optional (default 0.0)
        A distance threshold for cluster splits.

    max_cluster_size: int, default=None
        The maximum size for clusters located by the EOM clusterer. Can
        be overridden by the cluster_selection_epsilon parameter in
        rare cases.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        An integer array of cluster labels, with -1 denoting noise.

    probabilities : ndarray of shape (n_samples,)
        The cluster membership strength of each sample.

    stabilities : ndarray (n_clusters,)
        The cluster coherence strengths of each cluster.
    """
    cdef:
        list node_list
        cnp.ndarray[HIERARCHY_t, ndim=1] cluster_tree
        cnp.uint8_t[:] child_selection
        cnp.ndarray[cnp.intp_t, ndim=1] labels
        dict is_cluster, cluster_sizes
        cnp.float64_t subtree_stability, max_lambda
        cnp.intp_t node, sub_node, cluster, num_points
        cnp.ndarray[cnp.float64_t, ndim=1] probs

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        # exclude root
        node_list = sorted(stability.keys(), reverse=True)[:-1]

    cluster_tree = hierarchy[hierarchy['cluster_size'] > 1]
    is_cluster = {cluster: True for cluster in node_list}
    num_points = np.max(hierarchy[hierarchy['cluster_size'] == 1]['right_node']) + 1
    max_lambda = np.max(hierarchy['value'])

    if max_cluster_size is None:
        # Set to a value that will never be triggered
        max_cluster_size = num_points + 1
    cluster_sizes = {
        child: cluster_size for child, cluster_size
        in zip(cluster_tree['right_node'], cluster_tree['cluster_size'])
    }
    if allow_single_cluster:
        # Compute cluster size for the root node
        cluster_sizes[node_list[-1]] = np.sum(
            cluster_tree[cluster_tree['left_node'] == node_list[-1]]['cluster_size'])

    if cluster_selection_method == 'eom':
        for node in node_list:
            child_selection = (cluster_tree['left_node'] == node)
            subtree_stability = np.sum([
                stability[child] for
                child in cluster_tree['right_node'][child_selection]])
            if subtree_stability > stability[node] or cluster_sizes[node] > max_cluster_size:
                is_cluster[node] = False
                stability[node] = subtree_stability
            else:
                for sub_node in bfs_from_cluster_tree(cluster_tree, node):
                    if sub_node != node:
                        is_cluster[sub_node] = False

        if cluster_selection_epsilon != 0.0 and cluster_tree.shape[0] > 0:
            eom_clusters = [c for c in is_cluster if is_cluster[c]]
            selected_clusters = []
            # first check if eom_clusters only has root node, which skips epsilon check.
            if (len(eom_clusters) == 1 and eom_clusters[0] == cluster_tree['left_node'].min()):
                if allow_single_cluster:
                    selected_clusters = eom_clusters
            else:
                selected_clusters = epsilon_search(set(eom_clusters), cluster_tree, cluster_selection_epsilon, allow_single_cluster)
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
            is_cluster[hierarchy['left_node'].min()] = True

        if cluster_selection_epsilon != 0.0:
            selected_clusters = epsilon_search(leaves, cluster_tree, cluster_selection_epsilon, allow_single_cluster)
        else:
            selected_clusters = leaves

        for c in is_cluster:
                if c in selected_clusters:
                    is_cluster[c] = True
                else:
                    is_cluster[c] = False
    else:
        raise ValueError('Invalid Cluster Selection Method: %s\n'
                         'Should be one of: "eom", "leaf"\n')

    clusters = set([c for c in is_cluster if is_cluster[c]])
    cluster_map = {c: n for n, c in enumerate(sorted(list(clusters)))}
    reverse_cluster_map = {n: c for c, n in cluster_map.items()}

    labels = do_labelling(
        hierarchy,
        clusters,
        cluster_map,
        allow_single_cluster,
        cluster_selection_epsilon
    )
    probs = get_probabilities(hierarchy, reverse_cluster_map, labels)

    return (labels, probs)
