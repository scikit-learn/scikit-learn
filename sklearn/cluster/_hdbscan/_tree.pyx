# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD


cimport numpy as cnp

import cython
import numpy as np


cdef cnp.float64_t INFTY = np.inf

CONDENSED_dtype = np.dtype([
    ("parent", np.intp),
    ("child", np.intp),
    ("lambda_val", np.float64),
    ("cluster_size", np.intp),
])

HIERARCHY_dtype = np.dtype([
    ("left_child", np.intp),
    ("right_child", np.intp),
    ("lambda_val", np.float32),
    ("cluster_size", np.intp),
])

ctypedef packed struct CONDENSED_t:
    cnp.intp_t parent
    cnp.intp_t child
    cnp.float64_t lambda_val
    cnp.intp_t cluster_size

ctypedef packed struct HIERARCHY_t:
    cnp.intp_t left_child
    cnp.intp_t right_child
    cnp.float_t lambda_val
    cnp.intp_t cluster_size


cdef list bfs_from_hierarchy(
    cnp.ndarray[cnp.float64_t, ndim=2] hierarchy,
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
    cdef list process_queue
    cdef cnp.intp_t n_samples = hierarchy.shape[0] + 1

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
            process_queue = (
                hierarchy[process_queue, :2]
                .flatten()
                .astype(np.intp)
                .tolist()
            )
    return result

cdef inline void _shatter_cluster(
    cnp.ndarray[cnp.float64_t, ndim=2] hierarchy,
    cnp.uint8_t[::1] ignore,
    cnp.intp_t[::1] relabel,
    cnp.intp_t root,
    cnp.intp_t node,
    cnp.intp_t num_points,
    cnp.float64_t lambda_value,
    list result_list,
):
    """Take what would be a cluster but has too few points and assign each of
    its constituent points as single-sample clusters, removing them from
    consideration.
    """
    cdef cnp.intp_t sub_node
    for sub_node in bfs_from_hierarchy(hierarchy, root):
        if sub_node < num_points:
            result_list.append(
                (relabel[node], sub_node, lambda_value, 1)
            )
        ignore[sub_node] = True

cpdef cnp.ndarray condense_tree(
    cnp.ndarray[cnp.float64_t, ndim=2] hierarchy,
    cnp.intp_t min_cluster_size=10
):
    """Condense an MST according to a minimum cluster size. This is akin
    to the runt pruning procedure of Stuetzle. The result is a much simpler
    tree that is easier to visualize. We include extra information on the
    lambda value at which individual points depart clusters for later
    analysis and computation.

    Parameters
    ----------
    hierarchy : ndarray (n_samples, 4)
        A single linkage hierarchy in scipy.cluster.hierarchy format.

    min_cluster_size : int, optional (default 10)
        The minimum size of clusters to consider. Smaller "runt"
        clusters are pruned from the tree.

    Returns
    -------
    condensed_tree : numpy recarray
        Effectively an edgelist with a parent, child, lambda_val
        and cluster_size in each row providing a tree structure.
    """

    cdef cnp.intp_t root = 2 * hierarchy.shape[0]
    cdef cnp.intp_t num_points = hierarchy.shape[0] + 1
    cdef cnp.intp_t next_label = num_points + 1
    cdef list node_list = bfs_from_hierarchy(hierarchy, root)
    cdef list result_list

    cdef cnp.intp_t[::1] relabel
    cdef cnp.uint8_t[::1] ignore
    cdef cnp.float64_t[::1] children

    cdef cnp.intp_t node
    cdef cnp.intp_t sub_node
    cdef cnp.intp_t left
    cdef cnp.intp_t right
    cdef cnp.float64_t lambda_value
    cdef cnp.float64_t distance
    cdef cnp.intp_t left_count
    cdef cnp.intp_t right_count

    relabel = np.empty(root + 1, dtype=np.intp)
    relabel[root] = num_points
    result_list = []
    ignore = np.zeros(len(node_list), dtype=np.uint8)

    for node in node_list:
        # Guarantee that node is a cluster node of interest
        if ignore[node] or node < num_points:
            continue

        children = hierarchy[node - num_points]
        left = <cnp.intp_t> children[0]
        right = <cnp.intp_t> children[1]
        distance = children[2]
        if distance > 0.0:
            lambda_value = 1.0 / distance
        else:
            lambda_value = INFTY

        # Guarantee that left is a cluster node
        if left >= num_points:
            left_count = <cnp.intp_t> hierarchy[left - num_points][3]
        else:
            left_count = 1

        # Guarantee that right is a cluster node
        if right >= num_points:
            right_count = <cnp.intp_t> hierarchy[right - num_points][3]
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
            _shatter_cluster(
                hierarchy=hierarchy,
                ignore=ignore,
                relabel=relabel,
                root=left,
                node=node,
                num_points=num_points,
                lambda_value=lambda_value,
                result_list=result_list,
            )
            _shatter_cluster(
                hierarchy=hierarchy,
                ignore=ignore,
                relabel=relabel,
                root=right,
                node=node,
                num_points=num_points,
                lambda_value=lambda_value,
                result_list=result_list,
            )

        # One child is a collection of single-sample clusters, while the other
        # is a persistance of the parent node cluster
        elif left_count < min_cluster_size:
            relabel[right] = relabel[node]
            _shatter_cluster(
                hierarchy=hierarchy,
                ignore=ignore,
                relabel=relabel,
                root=left,
                node=node,
                num_points=num_points,
                lambda_value=lambda_value,
                result_list=result_list,
            )

        # One child is a collection of single-sample clusters, while the other
        # is a persistance of the parent node cluster
        else:
            relabel[left] = relabel[node]
            _shatter_cluster(
                hierarchy=hierarchy,
                ignore=ignore,
                relabel=relabel,
                root=right,
                node=node,
                num_points=num_points,
                lambda_value=lambda_value,
                result_list=result_list,
            )

    return np.array(result_list, dtype=CONDENSED_dtype)


cpdef dict compute_stability(cnp.ndarray[CONDENSED_t, ndim=1] condensed_tree):

    cdef cnp.float64_t[::1] result_arr
    cdef cnp.ndarray[cnp.intp_t, ndim=1] parents
    parents = condensed_tree['parent']

    cdef cnp.intp_t parent
    cdef cnp.intp_t cluster_size
    cdef cnp.intp_t result_index
    cdef cnp.float64_t lambda_val

    cdef cnp.float64_t[::1] births_arr
    cdef CONDENSED_t condensed_node


    cdef cnp.intp_t largest_child = condensed_tree['child'].max()
    cdef cnp.intp_t smallest_cluster = parents.min()
    cdef cnp.intp_t num_clusters = parents.max() - smallest_cluster + 1

    largest_child = max(largest_child, smallest_cluster)

    births_arr = np.full(largest_child + 1, np.nan, dtype=np.float64)

    for idx in range(condensed_tree.shape[0]):
        condensed_node = condensed_tree[idx]
        births_arr[condensed_node.child] = condensed_node.lambda_val

    births_arr[smallest_cluster] = 0.0

    result_arr = np.zeros(num_clusters, dtype=np.float64)
    for condensed_node in condensed_tree:
        parent = condensed_node.parent
        lambda_val = condensed_node.lambda_val
        cluster_size = condensed_node.cluster_size

        result_index = parent - smallest_cluster
        result_arr[result_index] += (lambda_val - births_arr[parent]) * cluster_size

    result_pre_dict = np.vstack(
        (
            np.arange(smallest_cluster, parents.max() + 1),
            result_arr
        )
    ).T

    return dict(result_pre_dict)


cdef list bfs_from_cluster_tree(cnp.ndarray[CONDENSED_t, ndim=1] tree, cnp.intp_t bfs_root):

    cdef list result
    cdef cnp.ndarray[cnp.intp_t, ndim=1] process_queue

    result = []
    process_queue = np.array([bfs_root], dtype=np.intp)

    while process_queue.shape[0] > 0:
        result.extend(process_queue.tolist())
        process_queue = tree['child'][np.isin(tree['parent'], process_queue)]

    return result


cdef cnp.float64_t[::1] max_lambdas(cnp.ndarray[CONDENSED_t, ndim=1] tree):

    cdef cnp.intp_t parent
    cdef cnp.intp_t current_parent
    cdef cnp.float64_t lambda_
    cdef cnp.float64_t max_lambda

    cdef cnp.float64_t[::1] deaths

    cdef cnp.intp_t largest_parent = tree['parent'].max()

    deaths = np.zeros(largest_parent + 1, dtype=np.float64)

    current_parent = tree[0].parent
    max_lambda = tree[0].lambda_val

    for idx in range(1, tree.shape[0]):
        parent = tree[idx].parent
        lambda_ = tree[idx].lambda_val

        if parent == current_parent:
            max_lambda = max(max_lambda, lambda_)
        else:
            deaths[current_parent] = max_lambda
            current_parent = parent
            max_lambda = lambda_

    deaths[current_parent] = max_lambda # value for last parent
    return deaths


cdef class TreeUnionFind (object):

    cdef cnp.ndarray _data_arr
    cdef cnp.intp_t[:, ::1] _data
    cdef cnp.ndarray is_component

    def __init__(self, size):
        self._data_arr = np.zeros((size, 2), dtype=np.intp)
        self._data_arr.T[0] = np.arange(size)
        self._data = (<cnp.intp_t[:size, :2:1]> (
            <cnp.intp_t *> self._data_arr.data))
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
    linkage : ndarray (n_samples, 4)
        The single linkage tree in scipy.cluster.hierarchy format.

    cut : float
        The cut value at which to find clusters.

    min_cluster_size : int
        The minimum cluster size; clusters below this size at
        the cut will be considered noise.

    Returns
    -------
    labels : ndarray (n_samples,)
        The cluster labels for each point in the data set;
        a label of -1 denotes a noise assignment.
    """

    cdef cnp.intp_t root
    cdef cnp.intp_t num_points
    cdef cnp.ndarray[cnp.intp_t, ndim=1] result_arr
    cdef cnp.ndarray[cnp.intp_t, ndim=1] unique_labels
    cdef cnp.ndarray[cnp.intp_t, ndim=1] cluster_size
    cdef cnp.intp_t *result
    cdef TreeUnionFind union_find
    cdef cnp.intp_t n
    cdef cnp.intp_t cluster
    cdef cnp.intp_t cluster_id

    root = 2 * linkage.shape[0]
    num_points = root // 2 + 1

    result_arr = np.empty(num_points, dtype=np.intp)
    result = (<cnp.intp_t *> result_arr.data)

    union_find = TreeUnionFind(<cnp.intp_t> root + 1)

    cluster = num_points
    for node in linkage:
        if node[2] < cut:
            union_find.union_(<cnp.intp_t> node[0], cluster)
            union_find.union_(<cnp.intp_t> node[1], cluster)
        cluster += 1

    cluster_size = np.zeros(cluster, dtype=np.intp)
    for n in range(num_points):
        cluster = union_find.find(n)
        cluster_size[cluster] += 1
        result[n] = cluster

    cluster_label_map = {-1: -1}
    cluster_label = 0
    unique_labels = np.unique(result_arr)

    for cluster in unique_labels:
        if cluster_size[cluster] < min_cluster_size:
            cluster_label_map[cluster] = -1
        else:
            cluster_label_map[cluster] = cluster_label
            cluster_label += 1

    for n in range(num_points):
        result[n] = cluster_label_map[result[n]]

    return result_arr


cdef cnp.ndarray[cnp.intp_t, ndim=1] do_labelling(
        cnp.ndarray tree,
        set clusters,
        dict cluster_label_map,
        cnp.intp_t allow_single_cluster,
        cnp.float64_t cluster_selection_epsilon):

    cdef cnp.intp_t root_cluster
    cdef cnp.ndarray[cnp.intp_t, ndim=1] result_arr
    cdef cnp.ndarray[cnp.intp_t, ndim=1] parent_array
    cdef cnp.ndarray[cnp.intp_t, ndim=1] child_array
    cdef cnp.ndarray[cnp.float64_t, ndim=1] lambda_array
    cdef cnp.intp_t *result
    cdef TreeUnionFind union_find
    cdef cnp.intp_t parent
    cdef cnp.intp_t child
    cdef cnp.intp_t n
    cdef cnp.intp_t cluster

    child_array = tree['child']
    parent_array = tree['parent']
    lambda_array = tree['lambda_val']

    root_cluster = parent_array.min()
    result_arr = np.empty(root_cluster, dtype=np.intp)
    result = (<cnp.intp_t *> result_arr.data)

    union_find = TreeUnionFind(parent_array.max() + 1)

    for n in range(tree.shape[0]):
        child = child_array[n]
        parent = parent_array[n]
        if child not in clusters:
            union_find.union_(parent, child)

    for n in range(root_cluster):
        cluster = union_find.find(n)
        if cluster < root_cluster:
            result[n] = -1
        elif cluster == root_cluster:
            if len(clusters) == 1 and allow_single_cluster:
                if cluster_selection_epsilon != 0.0:
                    if tree['lambda_val'][tree['child'] == n] >= 1 / cluster_selection_epsilon :
                        result[n] = cluster_label_map[cluster]
                    else:
                        result[n] = -1
                elif tree['lambda_val'][tree['child'] == n] >= \
                     tree['lambda_val'][tree['parent'] == cluster].max():
                    result[n] = cluster_label_map[cluster]
                else:
                    result[n] = -1
            else:
                result[n] = -1
        else:
            result[n] = cluster_label_map[cluster]

    return result_arr


cdef get_probabilities(cnp.ndarray[CONDENSED_t, ndim=1] tree, dict cluster_map, cnp.ndarray labels):

    cdef cnp.ndarray[cnp.float64_t, ndim=1] result
    cdef cnp.float64_t[::1] deaths
    cdef cnp.ndarray[cnp.float64_t, ndim=1] lambda_array
    cdef cnp.ndarray[cnp.intp_t, ndim=1] child_array
    cdef cnp.ndarray[cnp.intp_t, ndim=1] parent_array
    cdef cnp.intp_t root_cluster
    cdef cnp.intp_t n
    cdef cnp.intp_t point
    cdef cnp.intp_t cluster_num
    cdef cnp.intp_t cluster
    cdef cnp.float64_t max_lambda
    cdef cnp.float64_t lambda_

    child_array = tree['child']
    parent_array = tree['parent']
    lambda_array = tree['lambda_val']

    result = np.zeros(labels.shape[0])
    deaths = max_lambdas(tree)
    root_cluster = parent_array.min()

    for n in range(tree.shape[0]):
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
            lambda_ = min(lambda_array[n], max_lambda)
            result[point] = lambda_ / max_lambda

    return result


cpdef list recurse_leaf_dfs(cnp.ndarray cluster_tree, cnp.intp_t current_node):
    children = cluster_tree[cluster_tree['parent'] == current_node]['child']
    if len(children) == 0:
        return [current_node,]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


cpdef list get_cluster_tree_leaves(cnp.ndarray cluster_tree):
    if cluster_tree.shape[0] == 0:
        return []
    root = cluster_tree['parent'].min()
    return recurse_leaf_dfs(cluster_tree, root)

cpdef cnp.intp_t traverse_upwards(cnp.ndarray cluster_tree, cnp.float64_t cluster_selection_epsilon, cnp.intp_t leaf, cnp.intp_t allow_single_cluster):

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
        return traverse_upwards(cluster_tree, cluster_selection_epsilon, parent, allow_single_cluster)

cpdef set epsilon_search(set leaves, cnp.ndarray cluster_tree, cnp.float64_t cluster_selection_epsilon, cnp.intp_t allow_single_cluster):

    selected_clusters = list()
    processed = list()

    for leaf in leaves:
        eps = 1/cluster_tree['lambda_val'][cluster_tree['child'] == leaf][0]
        if eps < cluster_selection_epsilon:
            if leaf not in processed:
                epsilon_child = traverse_upwards(cluster_tree, cluster_selection_epsilon, leaf, allow_single_cluster)
                selected_clusters.append(epsilon_child)

                for sub_node in bfs_from_cluster_tree(cluster_tree, epsilon_child):
                    if sub_node != epsilon_child:
                        processed.append(sub_node)
        else:
            selected_clusters.append(leaf)

    return set(selected_clusters)

@cython.wraparound(True)
cpdef tuple get_clusters(cnp.ndarray tree, dict stability,
                         cluster_selection_method='eom',
                         allow_single_cluster=False,
                         cluster_selection_epsilon=0.0,
                         max_cluster_size=None):
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

    cluster_selection_epsilon: float, optional (default 0.0)
        A distance threshold for cluster splits.

    max_cluster_size: int, default=None
        The maximum size for clusters located by the EOM clusterer. Can
        be overridden by the cluster_selection_epsilon parameter in
        rare cases.

    Returns
    -------
    labels : ndarray (n_samples,)
        An integer array of cluster labels, with -1 denoting noise.

    probabilities : ndarray (n_samples,)
        The cluster membership strength of each sample.

    stabilities : ndarray (n_clusters,)
        The cluster coherence strengths of each cluster.
    """
    cdef list node_list
    cdef cnp.ndarray cluster_tree
    cdef cnp.ndarray child_selection
    cdef dict is_cluster
    cdef dict cluster_sizes
    cdef float subtree_stability
    cdef cnp.intp_t node
    cdef cnp.intp_t sub_node
    cdef cnp.intp_t cluster
    cdef cnp.intp_t num_points
    cdef cnp.ndarray labels
    cdef cnp.float64_t max_lambda

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        node_list = sorted(stability.keys(), reverse=True)[:-1]
        # (exclude root)

    cluster_tree = tree[tree['cluster_size'] > 1]
    is_cluster = {cluster: True for cluster in node_list}
    num_points = np.max(tree[tree['cluster_size'] == 1]['child']) + 1
    max_lambda = np.max(tree['lambda_val'])

    if max_cluster_size is None:
        max_cluster_size = num_points + 1  # Set to a value that will never be triggered
    cluster_sizes = {child: cluster_size for child, cluster_size
                 in zip(cluster_tree['child'], cluster_tree['cluster_size'])}
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

        if cluster_selection_epsilon != 0.0 and cluster_tree.shape[0] > 0:
            eom_clusters = [c for c in is_cluster if is_cluster[c]]
            selected_clusters = []
            # first check if eom_clusters only has root node, which skips epsilon check.
            if (len(eom_clusters) == 1 and eom_clusters[0] == cluster_tree['parent'].min()):
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
            is_cluster[tree['parent'].min()] = True

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

    labels = do_labelling(tree, clusters, cluster_map,
                          allow_single_cluster, cluster_selection_epsilon)
    probs = get_probabilities(tree, reverse_cluster_map, labels)

    return (labels, probs)
