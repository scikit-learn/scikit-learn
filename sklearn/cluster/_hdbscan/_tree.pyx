# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np

cimport numpy as cnp

import cython


cdef cnp.float64_t INFTY = np.inf
cdef cnp.intp_t NOISE = -1


cdef list bfs_from_hierarchy(
    cnp.ndarray[cnp.float64_t, ndim=2] hierarchy,
    cnp.intp_t bfs_root
):
    """
    Perform a breadth first search on a tree in scipy hclust format.
    """

    cdef list process_queue, next_queue
    cdef cnp.intp_t n_samples = hierarchy.shape[0] + 1
    cdef cnp.intp_t node
    process_queue = [bfs_root]
    result = []

    while process_queue:
        result.extend(process_queue)
        process_queue = [
            x - n_samples
            for x in process_queue
            if x >= n_samples
        ]
        if process_queue:
            process_queue = (
                hierarchy[process_queue, :2]
                .flatten()
                .astype(np.intp)
                .tolist()
            )

    return result


cpdef cnp.ndarray condense_tree(
    cnp.ndarray[cnp.float64_t, ndim=2] hierarchy,
    cnp.intp_t min_cluster_size=10
):
    """Condense a tree according to a minimum cluster size. This is akin
    to the runt pruning procedure of Stuetzle. The result is a much simpler
    tree that is easier to visualize. We include extra information on the
    lambda value at which individual points depart clusters for later
    analysis and computation.

    Parameters
    ----------
    hierarchy : ndarray (n_samples - 1, 4)
        A single linkage hierarchy in scipy.cluster.hierarchy format.

    min_cluster_size : int, optional (default 10)
        The minimum size of clusters to consider. Clusters smaler than this
        are pruned from the tree.

    Returns
    -------
    condensed_tree : numpy recarray
        Effectively an edgelist with a parent, child, lambda_val
        and child_size in each row providing a tree structure.
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
    relabel = np.empty(root + 1, dtype=np.intp)
    relabel[root] = n_samples
    result_list = []
    ignore = np.zeros(len(node_list), dtype=bool)

    for node in node_list:
        if ignore[node] or node < n_samples:
            continue

        children = hierarchy[node - n_samples]
        left = <cnp.intp_t> children[0]
        right = <cnp.intp_t> children[1]
        distance = children[2]
        if distance > 0.0:
            lambda_value = 1.0 / distance
        else:
            lambda_value = INFTY

        if left >= n_samples:
            left_count = <cnp.intp_t> hierarchy[left - n_samples][3]
        else:
            left_count = 1

        if right >= n_samples:
            right_count = <cnp.intp_t> hierarchy[right - n_samples][3]
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

    return np.array(result_list, dtype=[('parent', np.intp),
                                        ('child', np.intp),
                                        ('lambda_val', np.float64),
                                        ('child_size', np.intp)])


cpdef dict compute_stability(cnp.ndarray condensed_tree):

    cdef:
        cnp.float64_t[::1] result, births
        cnp.ndarray condensed_node
        cnp.intp_t[:] parents = condensed_tree['parent']
        cnp.float64_t[:] lambdas = condensed_tree['lambda_val']
        cnp.intp_t[:] sizes = condensed_tree['child_size']

        cnp.intp_t parent, cluster_size, result_index
        cnp.float64_t lambda_val
        cnp.float64_t[:, :] result_pre_dict
        cnp.intp_t largest_child = condensed_tree['child'].max()
        cnp.intp_t smallest_cluster = np.min(parents)
        cnp.intp_t num_clusters = np.max(parents) - smallest_cluster + 1
        cnp.ndarray sorted_child_data = np.sort(condensed_tree[['child', 'lambda_val']], axis=0)
        cnp.intp_t[:] sorted_children = sorted_child_data['child'].copy()
        cnp.float64_t[:] sorted_lambdas = sorted_child_data['lambda_val'].copy()

    largest_child = max(largest_child, smallest_cluster)
    births = np.full(largest_child + 1, np.nan, dtype=np.float64)

    if largest_child < smallest_cluster:
        largest_child = smallest_cluster

    births = np.nan * np.ones(largest_child + 1, dtype=np.float64)
    current_child = -1
    min_lambda = 0
    for idx in range(condensed_tree.shape[0]):
        child = sorted_children[idx]
        lambda_val = sorted_lambdas[idx]

        if child == current_child:
            min_lambda = min(min_lambda, lambda_val)
        elif current_child != -1:
            births[current_child] = min_lambda
            current_child = child
            min_lambda = lambda_val
        else:
            # Initialize
            current_child = child
            min_lambda = lambda_val

    if current_child != -1:
        births[current_child] = min_lambda
    births[smallest_cluster] = 0.0

    result = np.zeros(num_clusters, dtype=np.float64)
    for idx in range(condensed_tree.shape[0]):
        parent = parents[idx]
        lambda_val = lambdas[idx]
        child_size = sizes[idx]

        result_index = parent - smallest_cluster
        result[result_index] += (lambda_val - births[parent]) * child_size

    result_pre_dict = np.vstack(
        (
            np.arange(smallest_cluster, np.max(parents) + 1),
            result
        )
    ).T

    return dict(result_pre_dict)


cdef list bfs_from_cluster_tree(cnp.ndarray hierarchy, cnp.intp_t bfs_root):

    cdef list result
    cdef cnp.ndarray[cnp.intp_t, ndim=1] to_process

    result = []
    to_process = np.array([bfs_root], dtype=np.intp)

    while to_process.shape[0] > 0:
        result.extend(to_process.tolist())
        to_process = hierarchy['child'][np.in1d(hierarchy['parent'], to_process)]

    return result


cdef max_lambdas(cnp.ndarray hierarchy):

    cdef:
        cnp.ndarray sorted_parent_data
        cnp.intp_t[:] sorted_parents
        cnp.float64_t[:] sorted_lambdas, deaths
        cnp.intp_t parent, current_parent
        cnp.float64_t lambda_val, max_lambda
        cnp.intp_t largest_parent = hierarchy['parent'].max()

    sorted_parent_data = np.sort(hierarchy[['parent', 'lambda_val']], axis=0)
    deaths = np.zeros(largest_parent + 1, dtype=np.float64)
    sorted_parents = sorted_parent_data['parent']
    sorted_lambdas = sorted_parent_data['lambda_val']

    current_parent = -1
    max_lambda = 0

    for row in range(sorted_parent_data.shape[0]):
        parent = sorted_parents[row]
        lambda_val = sorted_lambdas[row]

        if parent == current_parent:
            max_lambda = max(max_lambda, lambda_val)
        elif current_parent != -1:
            deaths[current_parent] = max_lambda
            current_parent = parent
            max_lambda = lambda_val
        else:
            # Initialize
            current_parent = parent
            max_lambda = lambda_val

    deaths[current_parent] = max_lambda # value for last parent

    return deaths


cdef class TreeUnionFind (object):

    cdef cnp.ndarray _data_arr
    cdef cnp.intp_t[:, ::1] _data
    cdef cnp.ndarray is_component

    def __init__(self, size):
        self._data_arr = np.zeros((size, 2), dtype=np.intp)
        self._data_arr.T[0] = np.arange(size)
        self._data = self._data_arr
        self.is_component = np.ones(size, dtype=bool)

    cdef union_(self, cnp.intp_t x, cnp.intp_t y):
        cdef cnp.intp_t x_root = self.find(x)
        cdef cnp.intp_t y_root = self.find(y)

        if self._data[x_root, 1] < self._data[y_root, 1]:
            self._data[x_root, 0] = y_root
        elif self._data[x_root, 1] > self._data[y_root, 1]:
            self._data[y_root, 0] = x_root
        else:
            self._data[y_root, 0] = x_root
            self._data[x_root, 1] += 1

        return

    cdef find(self, cnp.intp_t x):
        if self._data[x, 0] != x:
            self._data[x, 0] = self.find(self._data[x, 0])
            self.is_component[x] = False
        return self._data[x, 0]

    cdef cnp.ndarray[cnp.intp_t, ndim=1] components(self):
        return self.is_component.nonzero()[0]


cpdef cnp.ndarray[cnp.intp_t, ndim=1] labelling_at_cut(
        cnp.ndarray linkage,
        cnp.float64_t cut,
        cnp.intp_t min_cluster_size
):
    """Given a single linkage tree and a cut value, return the
    vector of cluster labels at that cut value. This is useful
    for Robust Single Linkage, and extracting DBSCAN results
    from a single HDBSCAN run.

    Parameters
    ----------
    linkage : ndarray (n_samples - 1, 4)
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
        cnp.intp_t n, cluster, cluster_id, root, n_samples
        cnp.ndarray[cnp.intp_t, ndim=1] result
        cnp.intp_t[:] unique_labels, cluster_size
        TreeUnionFind union_find

    root = 2 * linkage.shape[0]
    n_samples = root // 2 + 1
    result = np.empty(n_samples, dtype=np.intp)
    union_find = TreeUnionFind(<cnp.intp_t> root + 1)

    cluster = n_samples
    for row in linkage:
        if row[2] < cut:
            union_find.union_(<cnp.intp_t> row[0], cluster)
            union_find.union_(<cnp.intp_t> row[1], cluster)
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


cdef cnp.ndarray[cnp.intp_t, ndim=1] do_labelling(
        cnp.ndarray hierarchy,
        set clusters,
        dict cluster_label_map,
        cnp.intp_t allow_single_cluster,
        cnp.float64_t cluster_selection_epsilon
):

    cdef:
        cnp.intp_t root_cluster
        cnp.ndarray[cnp.intp_t, ndim=1] result
        cnp.intp_t[:] parent_array, child_array
        cnp.float64_t[:] lambda_array
        TreeUnionFind union_find
        cnp.intp_t n, parent, child, cluster

    child_array = hierarchy['child']
    parent_array = hierarchy['parent']
    lambda_array = hierarchy['lambda_val']

    root_cluster = np.min(parent_array)
    result = np.empty(root_cluster, dtype=np.intp)
    union_find = TreeUnionFind(np.max(parent_array) + 1)

    for n in range(hierarchy.shape[0]):
        child = child_array[n]
        parent = parent_array[n]
        if child not in clusters:
            union_find.union_(parent, child)

    for n in range(root_cluster):
        cluster = union_find.find(n)
        if cluster < root_cluster:
            result[n] = NOISE
        elif cluster == root_cluster:
            if len(clusters) == 1 and allow_single_cluster:
                if cluster_selection_epsilon != 0.0:
                    if hierarchy['lambda_val'][hierarchy['child'] == n] >= 1 / cluster_selection_epsilon :
                        result[n] = cluster_label_map[cluster]
                    else:
                        result[n] = NOISE
                elif hierarchy['lambda_val'][hierarchy['child'] == n] >= \
                     hierarchy['lambda_val'][hierarchy['parent'] == cluster].max():
                    result[n] = cluster_label_map[cluster]
                else:
                    result[n] = NOISE
            else:
                result[n] = NOISE
        else:
            result[n] = cluster_label_map[cluster]

    return result


cdef get_probabilities(cnp.ndarray hierarchy, dict cluster_map, cnp.ndarray labels):

    cdef:
        cnp.ndarray[cnp.float64_t, ndim=1] result
        cnp.float64_t[:] lambda_array
        cnp.float64_t[::1] deaths
        cnp.intp_t[:] child_array, parent_array
        cnp.intp_t root_cluster, n, point, cluster_num, cluster
        cnp.float64_t max_lambda, lambda_val

    child_array = hierarchy['child']
    parent_array = hierarchy['parent']
    lambda_array = hierarchy['lambda_val']

    result = np.zeros(labels.shape[0])
    deaths = max_lambdas(hierarchy)
    root_cluster = np.min(parent_array)

    for n in range(hierarchy.shape[0]):
        point = child_array[n]
        if point >= root_cluster:
            continue

        cluster_num = labels[point]
        if cluster_num == -1:
            continue

        cluster = cluster_map[cluster_num]
        max_lambda = deaths[cluster]
        if max_lambda == 0.0 or not np.isfinite(lambda_array[n]):
            result[point] = 1.0
        else:
            lambda_val = min(lambda_array[n], max_lambda)
            result[point] = lambda_val / max_lambda

    return result


cpdef list recurse_leaf_dfs(cnp.ndarray cluster_tree, cnp.intp_t current_node):
    cdef cnp.intp_t[:] children
    cdef cnp.intp_t child

    children = cluster_tree[cluster_tree['parent'] == current_node]['child']
    if len(children) == 0:
        return [current_node,]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


cpdef list get_cluster_tree_leaves(cnp.ndarray cluster_tree):
    cdef cnp.intp_t root
    if cluster_tree.shape[0] == 0:
        return []
    root = cluster_tree['parent'].min()
    return recurse_leaf_dfs(cluster_tree, root)

cdef cnp.intp_t traverse_upwards(
    cnp.ndarray cluster_tree,
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
            return leaf #return node closest to root

    parent_eps = 1/cluster_tree[cluster_tree['child'] == parent]['lambda_val']
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
    cnp.ndarray cluster_tree,
    cnp.float64_t cluster_selection_epsilon,
    cnp.intp_t allow_single_cluster
):
    cdef:
        list selected_clusters = list()
        list processed = list()
        cnp.intp_t leaf, epsilon_child, sub_node
        cnp.float64_t eps

    for leaf in leaves:
        eps = 1/cluster_tree['lambda_val'][cluster_tree['child'] == leaf][0]
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
    tree : numpy recarray
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
        cnp.ndarray cluster_tree
        cnp.uint8_t[:] child_selection
        cnp.ndarray[cnp.intp_t, ndim=1] labels
        dict is_cluster, cluster_sizes
        cnp.float64_t subtree_stability, max_lambda
        cnp.intp_t node, sub_node, cluster, n_samples
        cnp.ndarray[cnp.float64_t, ndim=1] probs

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        node_list = sorted(stability.keys(), reverse=True)[:-1]
        # (exclude root)

    cluster_tree = hierarchy[hierarchy['child_size'] > 1]
    is_cluster = {cluster: True for cluster in node_list}
    n_samples = np.max(hierarchy[hierarchy['child_size'] == 1]['child']) + 1
    max_lambda = np.max(hierarchy['lambda_val'])

    if max_cluster_size is None:
        max_cluster_size = n_samples + 1  # Set to a value that will never be triggered
    cluster_sizes = {child: child_size for child, child_size
                 in zip(cluster_tree['child'], cluster_tree['child_size'])}
    if allow_single_cluster:
        # Compute cluster size for the root node
        cluster_sizes[node_list[-1]] = np.sum(
            cluster_tree[cluster_tree['parent'] == node_list[-1]]['child_size'])

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

        if cluster_selection_epsilon != 0.0 and cluster_tree.shape[0] > 0:
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
            is_cluster[hierarchy['parent'].min()] = True

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

    labels = do_labelling(
        hierarchy,
        clusters,
        cluster_map,
        allow_single_cluster,
        cluster_selection_epsilon
    )
    probs = get_probabilities(hierarchy, reverse_cluster_map, labels)

    return (labels, probs)
