#cython: boundscheck=False, nonecheck=False, initializedcheck=False
# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes
# License: 3-clause BSD

import numpy as np
cimport numpy as np

import itertools

cpdef np.ndarray get_points(np.ndarray[np.double_t, ndim=2] hierarchy):
    """
    Extract original point information from the single linkage hierarchy
    providing a list of points for each node in the hierarchy tree.
    """

    cdef int max_node
    cdef int num_points
    cdef int node
    cdef int parent
    cdef np.ndarray[np.double_t, ndim=1] row
    cdef np.ndarray result
    cdef int left_child
    cdef int right_child

    max_node = 2 * hierarchy.shape[0]
    num_points = max_node - (hierarchy.shape[0] - 1)

    result = np.empty(max_node + 1, dtype=object)

    for node in range(num_points):
        result[node] = [node]

    for parent, row in enumerate(hierarchy, num_points):
        left_child = <int> row[0]
        right_child = <int> row[1]
        result[parent] = result[left_child] + result[right_child]

    return result

cdef list bfs_from_hierarchy(np.ndarray[np.double_t, ndim=2] hierarchy, int bfs_root):

    """
    Perform a breadth first search on a tree in scipy hclust format.
    """
    
    cdef list to_process
    cdef int max_node
    cdef int num_points
    cdef int dim

    dim = hierarchy.shape[0]
    max_node = 2 * dim
    num_points = max_node - dim + 1

    to_process = [bfs_root]
    result = []

    while to_process:
        result.extend(to_process)
        to_process = [int(x - num_points) for x in 
                          to_process if x >= num_points]
        if to_process:
            to_process = hierarchy[to_process,:2].flatten().tolist()

    return result
        
cpdef tuple condense_tree(np.ndarray[np.double_t, ndim=2] hierarchy,
                          np.ndarray points,
                          int min_cluster_size=10):

    cdef int root
    cdef int num_points
    cdef int next_label
    cdef list node_list
    cdef list result_list
    cdef list new_points

    cdef np.ndarray[np.int64_t, ndim=1] relabel
    cdef np.ndarray[np.int_t, ndim=1] ignore
    cdef np.ndarray[np.double_t, ndim=1] children

    cdef int node
    cdef int sub_node
    cdef int left
    cdef int right
    cdef double lambda_value
    cdef int left_count
    cdef int right_count
    
    root = 2 * hierarchy.shape[0]
    num_points = root // 2 + 1
    next_label = num_points + 1
    
    node_list = bfs_from_hierarchy(hierarchy, root)
    
    relabel = np.empty(len(node_list), dtype=np.int64)
    relabel[root] = num_points
    result_list = []
    new_points = [[x] for x in range(num_points)]
    new_points.append(points[root])
    ignore = np.zeros(len(node_list), dtype=np.int)
    
    for node in node_list:
        if ignore[node] or node < num_points:
            continue
            
        children = hierarchy[node - num_points]
        left = <int> children[0]
        right = <int> children[1]
        lambda_value = 1.0 / children[2] if children[2] > 0.0 else np.inf
        left_count = <int> (hierarchy[left - num_points][3] 
                            if left >= num_points else 1)
        right_count = <int> (hierarchy[right - num_points][3] 
                             if right >= num_points else 1)
        
        if left_count > min_cluster_size and right_count > min_cluster_size:
            new_points.append(points[left])
            relabel[left] = next_label
            next_label += 1
            result_list.append((relabel[node], relabel[left], lambda_value, left_count))
            
            new_points.append(points[right])
            relabel[right] = next_label
            next_label += 1
            result_list.append((relabel[node], relabel[right], lambda_value, right_count))
            
        elif left_count <= min_cluster_size and right_count <= min_cluster_size:
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node, lambda_value, 1))
                ignore[sub_node] = True
                
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node, lambda_value, 1))
                ignore[sub_node] = True
                
        elif left_count <= min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node, lambda_value, 1))
                ignore[sub_node] = True
                
        else:
            relabel[left] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result_list.append((relabel[node], sub_node, lambda_value, 1))
                ignore[sub_node] = True
                
    return np.array(result_list, dtype=[
                                        ('parent', int),
                                        ('child', int),
                                        ('lambda', float),
                                        ('child_size', int)
                                       ]), new_points


cdef int keyfunc(row):
    return <int> row[0]
 
cpdef dict compute_stability(np.ndarray condensed_tree):

    cdef dict result
    cdef np.ndarray sorted_child_data
    cdef np.ndarray selection
    cdef list birth_pairs
    cdef dict births
    cdef int from_i
    cdef int to_i

    result = {}
    sorted_child_data = np.sort(condensed_tree[['child', 'lambda']], axis=0)
    birth_pairs = [min(a[1]) for a in itertools.groupby(sorted_child_data,
                                                        key=keyfunc)]
    births = dict(birth_pairs)
    from_i = condensed_tree['child'][condensed_tree['child_size'] == 1].max()
    from_i += 1
    to_i = condensed_tree['parent'].max() + 1
    for i in range(from_i, to_i):
        selection = (condensed_tree['parent'] == i)
        result[i] = np.sum((condensed_tree['lambda'][selection] - 
                            births.get(i, np.nan)) *
                           condensed_tree['child_size'][selection])
    return result

cdef list bfs_from_cluster_tree(np.ndarray tree, int bfs_root):

    cdef list result
    cdef list to_process

    result = []
    to_process = [bfs_root]
    
    while to_process:
        result.extend(to_process)
        to_process = tree['child'][np.in1d(tree['parent'], to_process)].tolist()

    return result

cpdef list get_clusters(np.ndarray tree, dict stability, list points):
    """
    The tree is assumed to have numeric node ids such that a reverse numeric
    sort is equivalent to a topological sort.
    """
    cdef list node_list
    cdef np.ndarray cluster_tree
    cdef np.ndarray child_selection
    cdef dict is_cluster
    cdef float subtree_stability
    cdef int node
    cdef int sub_node
    cdef int cluster

    # Assume clusters are ordered by numeric id equivalent to
    # a topological sort of the tree; This is valid given the
    # current implementation above, so don't change that ... or
    # if you do, change this accordingly!
    node_list = sorted(stability.keys(), reverse=True)[:-1] # (exclude root)
    cluster_tree = tree[tree['child_size'] > 1]
    is_cluster = {cluster:True for cluster in node_list}

    for node in node_list:
        child_selection = (cluster_tree['parent'] == node)
        subtree_stability = np.sum([stability[child] for 
                                    child in cluster_tree['child'][child_selection]])
        if subtree_stability > stability[node]:
            is_cluster[node] = False
            stability[node] = subtree_stability
        else:
            for sub_node in bfs_from_cluster_tree(cluster_tree, node):
                if sub_node != node:
                    is_cluster[sub_node] = False

    return [points[cluster] for cluster in is_cluster if is_cluster[cluster]]

    
    
    
       
