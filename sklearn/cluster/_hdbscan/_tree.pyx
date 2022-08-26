# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np

cimport numpy as np

import cython


cdef np.double_t INFTY = np.inf


cdef list bfs_from_hierarchy(np.ndarray[np.double_t, ndim=2] hierarchy,
                             np.intp_t bfs_root):
    """
    Perform a breadth first search on a tree in scipy hclust format.
    """

    cdef list to_process
    cdef np.intp_t max_node
    cdef np.intp_t num_points
    cdef np.intp_t dim

    dim = hierarchy.shape[0]
    max_node = 2 * dim
    num_points = max_node - dim + 1

    to_process = [bfs_root]
    result = []

    while to_process:
        result.extend(to_process)
        to_process = [x - num_points for x in
                      to_process if x >= num_points]
        if to_process:
            to_process = hierarchy[to_process,
                                   :2].flatten().astype(np.intp).tolist()

    return result


cpdef np.ndarray condense_tree(np.ndarray[np.double_t, ndim=2] hierarchy,
                               np.intp_t min_cluster_size=10):
    """Condense a tree according to a minimum cluster size. This is akin
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
        and child_size in each row providing a tree structure.
    """

    cdef np.intp_t root
    cdef np.intp_t num_points
    cdef np.intp_t next_label
    cdef list node_list
    cdef list result_list

    cdef np.ndarray[np.intp_t, ndim=1] relabel
    cdef np.ndarray[np.int_t, ndim=1] ignore
    cdef np.ndarray[np.double_t, ndim=1] children

    cdef np.intp_t node
    cdef np.intp_t sub_node
    cdef np.intp_t left
    cdef np.intp_t right
    cdef double lambda_value
    cdef np.intp_t left_count
    cdef np.intp_t right_count

    root = 2 * hierarchy.shape[0]
    num_points = root // 2 + 1
    next_label = num_points + 1

    node_list = bfs_from_hierarchy(hierarchy, root)

    relabel = np.empty(root + 1, dtype=np.intp)
    relabel[root] = num_points
    result_list = []
    ignore = np.zeros(len(node_list), dtype=int)

    for node in node_list:
        if ignore[node] or node < num_points:
            continue

        children = hierarchy[node - num_points]
        left = <np.intp_t> children[0]
        right = <np.intp_t> children[1]
        if children[2] > 0.0:
            lambda_value = 1.0 / children[2]
        else:
            lambda_value = INFTY

        if left >= num_points:
            left_count = <np.intp_t> hierarchy[left - num_points][3]
        else:
            left_count = 1

        if right >= num_points:
            right_count = <np.intp_t> hierarchy[right - num_points][3]
        else:
            right_count = 1

        if left_count >= min_cluster_size and right_count >= min_cluster_size:
            relabel[left] = next_label
            next_label += 1
            result_list.append((relabel[node], relabel[left], lambda_value,
                                left_count))

            relabel[right] = next_label
            next_label += 1
            result_list.append((relabel[node], relabel[right], lambda_value,
                                right_count))

        elif left_count < min_cluster_size and right_count < min_cluster_size:
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node,
                                        lambda_value, 1))
                ignore[sub_node] = True

            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node,
                                        lambda_value, 1))
                ignore[sub_node] = True

        elif left_count < min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node,
                                        lambda_value, 1))
                ignore[sub_node] = True

        else:
            relabel[left] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node,
                                        lambda_value, 1))
                ignore[sub_node] = True

    return np.array(result_list, dtype=[('parent', np.intp),
                                        ('child', np.intp),
                                        ('lambda_val', float),
                                        ('child_size', np.intp)])


cpdef dict compute_stability(np.ndarray condensed_tree):

    cdef np.ndarray[np.double_t, ndim=1] result_arr
    cdef np.ndarray sorted_child_data
    cdef np.ndarray[np.intp_t, ndim=1] sorted_children
    cdef np.ndarray[np.double_t, ndim=1] sorted_lambdas

    cdef np.ndarray[np.intp_t, ndim=1] parents
    cdef np.ndarray[np.intp_t, ndim=1] sizes
    cdef np.ndarray[np.double_t, ndim=1] lambdas

    cdef np.intp_t child
    cdef np.intp_t parent
    cdef np.intp_t child_size
    cdef np.intp_t result_index
    cdef np.intp_t current_child
    cdef np.float64_t lambda_
    cdef np.float64_t min_lambda

    cdef np.ndarray[np.double_t, ndim=1] births_arr
    cdef np.double_t *births

    cdef np.intp_t largest_child = condensed_tree['child'].max()
    cdef np.intp_t smallest_cluster = condensed_tree['parent'].min()
    cdef np.intp_t num_clusters = (condensed_tree['parent'].max() -
                                   smallest_cluster + 1)

    if largest_child < smallest_cluster:
        largest_child = smallest_cluster

    sorted_child_data = np.sort(condensed_tree[['child', 'lambda_val']],
                                axis=0)
    births_arr = np.nan * np.ones(largest_child + 1, dtype=np.double)
    births = (<np.double_t *> births_arr.data)
    sorted_children = sorted_child_data['child'].copy()
    sorted_lambdas = sorted_child_data['lambda_val'].copy()

    parents = condensed_tree['parent']
    sizes = condensed_tree['child_size']
    lambdas = condensed_tree['lambda_val']

    current_child = -1
    min_lambda = 0

    for row in range(sorted_child_data.shape[0]):
        child = <np.intp_t> sorted_children[row]
        lambda_ = sorted_lambdas[row]

        if child == current_child:
            min_lambda = min(min_lambda, lambda_)
        elif current_child != -1:
            births[current_child] = min_lambda
            current_child = child
            min_lambda = lambda_
        else:
            # Initialize
            current_child = child
            min_lambda = lambda_

    if current_child != -1:
        births[current_child] = min_lambda
    births[smallest_cluster] = 0.0

    result_arr = np.zeros(num_clusters, dtype=np.double)

    for i in range(condensed_tree.shape[0]):
        parent = parents[i]
        lambda_ = lambdas[i]
        child_size = sizes[i]
        result_index = parent - smallest_cluster

        result_arr[result_index] += (lambda_ - births[parent]) * child_size

    result_pre_dict = np.vstack((np.arange(smallest_cluster,
                                           condensed_tree['parent'].max() + 1),
                                 result_arr)).T

    return dict(result_pre_dict)


cdef list bfs_from_cluster_tree(np.ndarray tree, np.intp_t bfs_root):

    cdef list result
    cdef np.ndarray[np.intp_t, ndim=1] to_process

    result = []
    to_process = np.array([bfs_root], dtype=np.intp)

    while to_process.shape[0] > 0:
        result.extend(to_process.tolist())
        to_process = tree['child'][np.in1d(tree['parent'], to_process)]

    return result


cdef max_lambdas(np.ndarray tree):

    cdef np.ndarray sorted_parent_data
    cdef np.ndarray[np.intp_t, ndim=1] sorted_parents
    cdef np.ndarray[np.double_t, ndim=1] sorted_lambdas

    cdef np.intp_t parent
    cdef np.intp_t current_parent
    cdef np.float64_t lambda_
    cdef np.float64_t max_lambda

    cdef np.ndarray[np.double_t, ndim=1] deaths_arr
    cdef np.double_t *deaths

    cdef np.intp_t largest_parent = tree['parent'].max()

    sorted_parent_data = np.sort(tree[['parent', 'lambda_val']], axis=0)
    deaths_arr = np.zeros(largest_parent + 1, dtype=np.double)
    deaths = (<np.double_t *> deaths_arr.data)
    sorted_parents = sorted_parent_data['parent']
    sorted_lambdas = sorted_parent_data['lambda_val']

    current_parent = -1
    max_lambda = 0

    for row in range(sorted_parent_data.shape[0]):
        parent = <np.intp_t> sorted_parents[row]
        lambda_ = sorted_lambdas[row]

        if parent == current_parent:
            max_lambda = max(max_lambda, lambda_)
        elif current_parent != -1:
            deaths[current_parent] = max_lambda
            current_parent = parent
            max_lambda = lambda_
        else:
            # Initialize
            current_parent = parent
            max_lambda = lambda_

    deaths[current_parent] = max_lambda # value for last parent

    return deaths_arr


cdef class TreeUnionFind (object):

    cdef np.ndarray _data_arr
    cdef np.intp_t[:, ::1] _data
    cdef np.ndarray is_component

    def __init__(self, size):
        self._data_arr = np.zeros((size, 2), dtype=np.intp)
        self._data_arr.T[0] = np.arange(size)
        self._data = (<np.intp_t[:size, :2:1]> (
            <np.intp_t *> self._data_arr.data))
        self.is_component = np.ones(size, dtype=bool)

    cdef union_(self, np.intp_t x, np.intp_t y):
        cdef np.intp_t x_root = self.find(x)
        cdef np.intp_t y_root = self.find(y)

        if self._data[x_root, 1] < self._data[y_root, 1]:
            self._data[x_root, 0] = y_root
        elif self._data[x_root, 1] > self._data[y_root, 1]:
            self._data[y_root, 0] = x_root
        else:
            self._data[y_root, 0] = x_root
            self._data[x_root, 1] += 1

        return

    cdef find(self, np.intp_t x):
        if self._data[x, 0] != x:
            self._data[x, 0] = self.find(self._data[x, 0])
            self.is_component[x] = False
        return self._data[x, 0]

    cdef np.ndarray[np.intp_t, ndim=1] components(self):
        return self.is_component.nonzero()[0]


cpdef np.ndarray[np.intp_t, ndim=1] labelling_at_cut(
        np.ndarray linkage,
        np.double_t cut,
        np.intp_t min_cluster_size):
    """Given a single linkage tree and a cut value, return the
    vector of cluster labels at that cut value. This is useful
    for Robust Single Linkage, and extracting DBSCAN results
    from a single HDBSCAN run.

    Parameters
    ----------
    linkage : ndarray (n_samples, 4)
        The single linkage tree in scipy.cluster.hierarchy format.

    cut : double
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

    cdef np.intp_t root
    cdef np.intp_t num_points
    cdef np.ndarray[np.intp_t, ndim=1] result_arr
    cdef np.ndarray[np.intp_t, ndim=1] unique_labels
    cdef np.ndarray[np.intp_t, ndim=1] cluster_size
    cdef np.intp_t *result
    cdef TreeUnionFind union_find
    cdef np.intp_t n
    cdef np.intp_t cluster
    cdef np.intp_t cluster_id

    root = 2 * linkage.shape[0]
    num_points = root // 2 + 1

    result_arr = np.empty(num_points, dtype=np.intp)
    result = (<np.intp_t *> result_arr.data)

    union_find = TreeUnionFind(<np.intp_t> root + 1)

    cluster = num_points
    for row in linkage:
        if row[2] < cut:
            union_find.union_(<np.intp_t> row[0], cluster)
            union_find.union_(<np.intp_t> row[1], cluster)
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


cdef np.ndarray[np.intp_t, ndim=1] do_labelling(
        np.ndarray tree,
        set clusters,
        dict cluster_label_map,
        np.intp_t allow_single_cluster,
        np.double_t cluster_selection_epsilon):

    cdef np.intp_t root_cluster
    cdef np.ndarray[np.intp_t, ndim=1] result_arr
    cdef np.ndarray[np.intp_t, ndim=1] parent_array
    cdef np.ndarray[np.intp_t, ndim=1] child_array
    cdef np.ndarray[np.double_t, ndim=1] lambda_array
    cdef np.intp_t *result
    cdef TreeUnionFind union_find
    cdef np.intp_t parent
    cdef np.intp_t child
    cdef np.intp_t n
    cdef np.intp_t cluster

    child_array = tree['child']
    parent_array = tree['parent']
    lambda_array = tree['lambda_val']

    root_cluster = parent_array.min()
    result_arr = np.empty(root_cluster, dtype=np.intp)
    result = (<np.intp_t *> result_arr.data)

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


cdef get_probabilities(np.ndarray tree, dict cluster_map, np.ndarray labels):

    cdef np.ndarray[np.double_t, ndim=1] result
    cdef np.ndarray[np.double_t, ndim=1] deaths
    cdef np.ndarray[np.double_t, ndim=1] lambda_array
    cdef np.ndarray[np.intp_t, ndim=1] child_array
    cdef np.ndarray[np.intp_t, ndim=1] parent_array
    cdef np.intp_t root_cluster
    cdef np.intp_t n
    cdef np.intp_t point
    cdef np.intp_t cluster_num
    cdef np.intp_t cluster
    cdef np.double_t max_lambda
    cdef np.double_t lambda_

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


cpdef list recurse_leaf_dfs(np.ndarray cluster_tree, np.intp_t current_node):
    children = cluster_tree[cluster_tree['parent'] == current_node]['child']
    if len(children) == 0:
        return [current_node,]
    else:
        return sum([recurse_leaf_dfs(cluster_tree, child) for child in children], [])


cpdef list get_cluster_tree_leaves(np.ndarray cluster_tree):
    if cluster_tree.shape[0] == 0:
        return []
    root = cluster_tree['parent'].min()
    return recurse_leaf_dfs(cluster_tree, root)

cpdef np.intp_t traverse_upwards(np.ndarray cluster_tree, np.double_t cluster_selection_epsilon, np.intp_t leaf, np.intp_t allow_single_cluster):

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

cpdef set epsilon_search(set leaves, np.ndarray cluster_tree, np.double_t cluster_selection_epsilon, np.intp_t allow_single_cluster):

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
cpdef tuple get_clusters(np.ndarray tree, dict stability,
                         cluster_selection_method='eom',
                         allow_single_cluster=False,
                         cluster_selection_epsilon=0.0,
                         max_cluster_size=0):
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

    max_cluster_size: int, optional (default 0)
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
    cdef np.ndarray cluster_tree
    cdef np.ndarray child_selection
    cdef dict is_cluster
    cdef dict cluster_sizes
    cdef float subtree_stability
    cdef np.intp_t node
    cdef np.intp_t sub_node
    cdef np.intp_t cluster
    cdef np.intp_t num_points
    cdef np.ndarray labels
    cdef np.double_t max_lambda

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    if allow_single_cluster:
        node_list = sorted(stability.keys(), reverse=True)
    else:
        node_list = sorted(stability.keys(), reverse=True)[:-1]
        # (exclude root)

    cluster_tree = tree[tree['child_size'] > 1]
    is_cluster = {cluster: True for cluster in node_list}
    num_points = np.max(tree[tree['child_size'] == 1]['child']) + 1
    max_lambda = np.max(tree['lambda_val'])

    if max_cluster_size <= 0:
        max_cluster_size = num_points + 1  # Set to a value that will never be triggered
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
