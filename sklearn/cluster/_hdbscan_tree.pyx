#cython: boundscheck=False, wraparound=False
# Tree handling (condensing, finding stable clusters) for hdbscan
# Authors: Leland McInnes, David Thomson, John Healy
# License: 3-clause BSD

import numpy as np
cimport numpy as np

import igraph
import pandas as pd

vint = np.vectorize(int)

cpdef object igraph_to_tree(np.ndarray[np.double_t, ndim=2] hierarchy):

    cdef np.ndarray[np.long_t, ndim=2] left
    cdef np.ndarray[np.long_t, ndim=2] right
    cdef np.ndarray[np.long_t, ndim=2] edge_list
 
    cdef np.ndarray[np.double_t, ndim=1] row
    
    cdef int max_node
    cdef int num_points
    cdef int node
    cdef int parent
    cdef int left_child
    cdef int right_child
    cdef int dim
    
    dim = hierarchy.shape[0]
    
    left = np.vstack((np.arange(dim + 1, dim * 2 + 1),
                      hierarchy[:,0].astype(np.int))).T
    right = np.vstack((np.arange(dim + 1, dim * 2 + 1),
                       hierarchy[:,1].astype(np.int))).T
    edge_list = np.vstack((left, right))
    
    max_node = 2 * dim
    
    node_attr_dict = {
        'dist' : np.zeros(max_node + 1, dtype=np.float64),
        'count' : np.zeros(max_node + 1, dtype=np.int64),
        'points' : np.zeros(max_node + 1, dtype=object)
    }
 
    num_points = max_node - dim + 1
    
    for node in range(num_points):
        node_attr_dict['dist'][node] = 0.0
        node_attr_dict['count'][node] = 1
        node_attr_dict['points'][node] = [node]
        
    for parent, row in enumerate(hierarchy, num_points):
        left_child, right_child = int(row[0]), int(row[1])
        node_attr_dict['dist'][parent] = row[2]
        node_attr_dict['count'][parent] = int(row[3])
        node_attr_dict['points'][parent] = \
                                      node_attr_dict['points'][left_child] + \
                                      node_attr_dict['points'][right_child]

    return igraph.Graph(n=max_node, edges=vint(edge_list).tolist(), 
                        directed=True, vertex_attrs=node_attr_dict)

cpdef np.ndarray get_points(np.ndarray[np.double_t, ndim=2] hierarchy):

    cdef np.ndarray[np.double_t, ndim=1] row
    cdef np.ndarray result
    
    cdef int max_node
    cdef int num_points
    cdef int node
    cdef int parent
    cdef int left_child
    cdef int right_child
    cdef int dim

    dim = hierarchy.shape[0]
    max_node = 2 * dim
    num_points = max_node - dim + 1
    result = np.empty(dim, dtype=np.object)

    for node in range(num_points):
        result[node] = [node]
        
    for parent, row in enumerate(hierarchy, num_points):
        left_child, right_child = <int> row[0], <int> row[1]
        result[parent] = result[left_child] + result[right_child]

    return result

cdef list bfs_from_hierarchy(np.ndarray[np.double_t, ndim=2] hierarchy, int start):

    cdef list for_processing
    cdef int max_node
    cdef int num_points
    cdef int dim

    dim = hierarchy.shape[0]
    max_node = 2 * dim
    num_points = max_node - dim + 1

    for_processing = [start]
    result = []

    while len(for_processing) > 0:
        result.extend(for_processing)
        for_processing = [(x - num_points) for x in 
                          for_processing if x >= num_points]
        if for_processing:
            for_processing = for_processing.flatten().tolist()

    return result
        
cpdef tuple condense_tree(np.ndarray[np.double_t, ndim=2] hierarchy,
                          np.ndarray points,
                          int min_cluster_size=10):

    cdef int root
    cdef int num_points
    cdef int next_label
    cdef list node_list
    cdef np.ndarray[np.int64_t, ndim=1] relabel
    cdef list new_points
    cdef np.ndarray[np.int64_t, ndim=1] ignore
    cdef np.ndarray children
    cdef int node
    cdef int sub_node
    cdef int left
    cdef int right
    cdef double lambda_value
    cdef int left_count
    cdef int right_count
    cdef int edge_index
    cdef np.ndarray result
    
    root = 2 * hierarchy.shape[0]
    num_points = root // 2 + 1
    next_label = num_points + 1
    
    node_list = bfs_from_hierarchy(hierarchy, root)
    
    result = np.empty(len(node_list), dtype=[
        ('parent', int),
        ('child', int),
        ('lambda', float),
        ('size', int)
    ])
    edge_index = 0
    relabel = np.empty(len(node_list))
    relabel[root] = num_points
    edge_list = []
    edge_attrs = {'lambda':[], 'size':[]}
    new_points = [[x] for x in range(num_points)]
    new_points.append(points[root])
    ignore = np.zeros(len(node_list), dtype=np.int)
    
    for node in node_list:
        if ignore[node] or node < num_points:
            continue
            
        children = hierarchy[node - num_points]
        left, right = <int> children[0], <int> children[1]
        lambda_value = 1.0 / children[2]
        left_count = hierarchy[left - num_points][3] if left >= num_points else 1
        right_count = hierarchy[right - num_points][3] if right >= num_points else 1
        
        if left_count > min_cluster_size and right_count > min_cluster_size:
            new_points.append(points[left])
            relabel[left] = next_label
            next_label += 1
            result[edge_index] = (relabel[node], relabel[left], lambda_value, left_count)
            edge_index += 1
            
            new_points.append(points[right])
            relabel[right] = next_label
            next_label += 1
            result[edge_index] = (relabel[node], relabel[right], lambda_value, right_count)
            edge_index += 1
            
        elif left_count <= min_cluster_size and right_count <= min_cluster_size:
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result[edge_index] = (relabel[node], sub_node, lambda_value, 1)
                    edge_index += 1
                ignore[sub_node] = True
                
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result[edge_index] = (relabel[node], sub_node, lambda_value, 1)
                    edge_index += 1
                ignore[sub_node] = True
                
        elif left_count <= min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, left):
                if sub_node < num_points:
                    result[edge_index] = (relabel[node], sub_node, lambda_value, 1)
                    edge_index += 1
                ignore[sub_node] = True
                
        else:
            relabel[left] = relabel[node]
            for sub_node in bfs_from_hierarchy(hierarchy, right):
                if sub_node < num_points:
                    result[edge_index] = (relabel[node], sub_node, lambda_value, 1)
                    edge_index += 1
                ignore[sub_node] = True
                
    return result, np.ndarray(new_points, dtype=object)
                        
cpdef object igraph_condense_tree(object tree, int min_cluster_size=10):

    cdef int root
    cdef int num_points
    cdef int next_label
    cdef list node_list
    cdef dict relabel
    cdef list edge_list
    cdef dict edge_attrs
    cdef dict node_attrs
    cdef dict ignore
    cdef int node
    cdef int left
    cdef int right
    cdef double lambda_value
    cdef int left_count
    cdef int right_count
    
    root = max(tree.vs).index
    num_points = root // 2 + 1
    next_label = num_points + 1
    
    node_list = tree.bfs(root)[0]
    
    relabel = {root:num_points}
    edge_list = []
    edge_attrs = {'lambda':[], 'size':[]}
    node_attrs = {'points':[[x] for x in range(num_points)]}
    node_attrs['points'].append(tree.vs[root]['points'])
    ignore = {x:False for x in node_list}
    
    for node in node_list:
        if ignore[node] or node < num_points:
            continue
            
        left, right = tree.successors(node)
        lambda_value = 1.0 / tree.vs[node]['dist']
        left_count = tree.vs[left]['count']
        right_count = tree.vs[right]['count']
        
        if left_count > min_cluster_size and right_count > min_cluster_size:
            node_attrs['points'].append(tree.vs[left]['points'])
            relabel[left] = next_label
            next_label += 1
            edge_list.append([relabel[node],relabel[left]])
            edge_attrs['lambda'].append(lambda_value)
            edge_attrs['size'].append(left_count)
            
            node_attrs['points'].append(tree.vs[right]['points'])
            relabel[right] = next_label
            next_label += 1
            edge_list.append([relabel[node],relabel[right]])
            edge_attrs['lambda'].append(lambda_value)
            edge_attrs['size'].append(right_count)
            
        elif left_count <= min_cluster_size and right_count <= min_cluster_size:
            for sub_node in tree.bfsiter(left):
                if sub_node.index < num_points:
                    edge_list.append([relabel[node], sub_node.index])
                    edge_attrs['lambda'].append(lambda_value)
                    edge_attrs['size'].append(1)
                ignore[sub_node.index] = True
                
            for sub_node in tree.bfsiter(right):
                if sub_node.index < num_points:
                    edge_list.append([relabel[node], sub_node.index])
                    edge_attrs['lambda'].append(lambda_value)
                    edge_attrs['size'].append(1)
                ignore[sub_node.index] = True
                
        elif left_count <= min_cluster_size:
            relabel[right] = relabel[node]
            for sub_node in tree.bfsiter(left):
                if sub_node.index < num_points:
                    edge_list.append([relabel[node], sub_node.index])
                    edge_attrs['lambda'].append(lambda_value)
                    edge_attrs['size'].append(1)
                ignore[sub_node.index] = True
                
        else:
            relabel[left] = relabel[node]
            for sub_node in tree.bfsiter(right):
                if sub_node.index < num_points:
                    edge_list.append([relabel[node], sub_node.index])
                    edge_attrs['lambda'].append(lambda_value)
                    edge_attrs['size'].append(1)
                ignore[sub_node.index] = True
                
    return igraph.Graph(len(node_attrs['points']), edges=edge_list,
                        directed=True, vertex_attrs=node_attrs,
                        edge_attrs=edge_attrs)
                        
def igraph_tree_to_dataframe(tree):
    return pd.DataFrame([(edge.source, edge.target, edge['lambda'], 
                          edge['size']) for edge in tree.es],
                        columns=('parent', 'child', 'lambda', 'child_size'))

cpdef list igraph_get_clusters(object tree, dict stability):

    cdef int root
    cdef int node
    cdef float subtree_stability
    cdef int x
    
    root = tree.topological_sorting(mode='OUT')[0]
    
    # Set default states
    tree.vs.set_attribute_values('is_cluster', True)
    tree.vs.set_attribute_values('stability', 0.0)
    
    # Set stability for nodes that have it
    tree.vs[stability.keys()].set_attrbute_values('stability', 
                                                  stability.values())
    
    # We ignore the root for stability
    tree.vs[root]['is_cluster'] = False
    tree.vs[root]['stability'] = 0.0
    
    for node in tree.topological_sorting(mode='IN'):
        # Nodes less than the root are points and we can ignore them
        if node <= root:
            tree.vs[node]['is_cluster'] = False
            continue
            
        subtree_stability = sum([tree.vs[child]['stability'] for child in
                                 tree.successors(node)])
        if tree.vs[node]['stability'] < subtree_stability:
            tree.vs[node]['is_cluster'] = False
            tree.vs[node]['stability'] = subtree_stability            
        else:
            for sub_node in tree.bfsiter(node):
                if sub_node.index != node:
                    sub_node['is_cluster'] = False
                    
    return [tree.vs[x]['points'] for x in range(len(tree.vs)) if 
            tree.vs[x]['is_cluster']]

            
            
    
    
    
       
