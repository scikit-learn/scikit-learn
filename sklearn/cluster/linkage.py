
# TODO: Consider nearest-neighbor chain algorithm ?
# https://secure.wikimedia.org/wikipedia/en/wiki/Nearest-neighbor_chain_algorithm

from collections import defaultdict
import itertools
from heapq import heapify, heappop, heappush

import numpy as np

from . import _inertia


class Linkage(object):
    
    def __init__(self):
        self.A = []
        
    def has_more_candidates(self):
        return len(self.distances) > 0
    
    def fetch_candidate(self):
        return heappop(self.distances)
    
    def get_nodes_connected_to_clusters(self, cluster_roots):
        return set.union(*[set(self.A[cluster_root]) 
                                for cluster_root in cluster_roots])
        
    def get_nodes_of_clusters(self, cluster_roots, children, n_samples):
        nodes = set()
        for cluster_root in cluster_roots:
            if cluster_root < n_samples:
                nodes.add(cluster_root)
            else:
                nodes = nodes.union(
                    self.get_nodes_of_clusters(children[cluster_root - n_samples],
                                               children, n_samples))
        return list(nodes)
        

class WardsLinkage(Linkage):

    def __init__(self, X, connectivity, n_nodes, n_features, n_samples):
        super(WardsLinkage, self).__init__()
        # build moments as a list
        self.moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features))]
        self.moments[0][:n_samples] = 1
        self.moments[1][:n_samples] = X

        # create distances matrix
        coord_row = []
        coord_col = []
        for ind, row in enumerate(connectivity.rows):
            self.A.append(row)
            # We keep only the upper triangular for the moments
            # Generator expressions are faster than arrays on the following
            row = [i for i in row if i < ind]
            coord_row.extend(len(row) * [ind, ])
            coord_col.extend(row)
    
        coord_row = np.array(coord_row, dtype=np.int)
        coord_col = np.array(coord_col, dtype=np.int)
        
        distances = np.empty(len(coord_row), dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   coord_row, coord_col, distances)
    
        self.distances = zip(distances, coord_row, coord_col)
        heapify(self.distances)

    def update(self, child_node1, child_node2, parent_node, parent):
        # update the moments
        for p in range(2):
            self.moments[p][parent_node] = \
                    self.moments[p][child_node1] + self.moments[p][child_node2]
        # update the structure matrix A and the inertia matrix
        coord_col = []
        visited = np.empty(parent.shape[0], dtype=bool)
        visited[:] = False
        visited[parent_node] = True     
        for l in set(self.A[child_node1]).union(self.A[child_node2]):
            while parent[l] != l:
                l = parent[l]
            if not visited[l]:
                visited[l] = True
                coord_col.append(l)
                self.A[l].append(parent_node)
        self.A.append(coord_col)
        coord_col = np.array(coord_col, dtype=np.int)
        coord_row = np.empty_like(coord_col)
        coord_row.fill(parent_node)
        
        distances = np.empty(len(coord_row), dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   coord_row, coord_col, distances)
        
        for tupl in itertools.izip(distances, coord_row, coord_col):
            heappush(self.distances, tupl)
    
    def compute_distance(self, i, j):
        # Compute inertia of the cluster obtained when merging i and j
        distances = np.empty(1, dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   np.array([i]), np.array([j]), distances)
        return distances[0]


class ConstrainedCompleteLinkage(Linkage):
    
    def __init__(self, X, connectivity, *args):
        super(ConstrainedCompleteLinkage, self).__init__()
        
        self.X = X
        # Heap of possible cluster merges, sorted by heir distances
        self.distances = [] 
        # Distances between two clusters, represented by their 
        # root node indices 
        self.distance_dict = {} 
        # Mapping from a pair of clusters to the set of node pairs which
        # connect these two clusters
        self.cluster_connections = defaultdict(set)
        # Mapping from a a cluster to a mapping from nodes in this cluster
        # to their distance to the most distant node in the cluster
        self.max_dist_of_node_in_cluster = defaultdict(dict)  
        # Determine distances between all connected nodes 
        for ind1, row in enumerate(connectivity.rows):
            self.A.append(row)
            self.max_dist_of_node_in_cluster[ind1][ind1] = 0.0
            for ind2 in row:
                self.cluster_connections[(ind1, ind2)] = set([(ind1, ind2)])
                self.cluster_connections[(ind2, ind1)] = set([(ind1, ind2)])
                # Compute distance between two connected nodes
                dist = np.linalg.norm(X[ind1] - X[ind2])
                self.distances.append((dist, ind1, ind2))
                self.distance_dict[(ind1, ind2)] = dist
                self.distance_dict[(ind2, ind1)] = dist
                
        # Enforce symmetry of A
        for ind1, row in enumerate(connectivity.rows):
            for ind2 in row:
                if ind1 not in self.A[ind2]:
                    self.A[ind2].append(ind1)
        
        heapify(self.distances)

    def update(self, child_node1, child_node2, parent_node, parent):
        # Update max_dist_of_node_in_cluster for new cluster "parent_node" 
        for node in self.max_dist_of_node_in_cluster[child_node1].keys():
            min_inter_cluster_dist = np.inf
            for (node1, node2) in self.cluster_connections[(child_node1, 
                                                            child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if node1 not in self.max_dist_of_node_in_cluster[child_node1]:
                    node1, node2 = node2, node1
                    
                # TODO: graph distance instead of euclidean
                inter_cluster_dist = \
                    self.max_dist_of_node_in_cluster[child_node2][node2] \
                        + self.distance_dict[(node1, node2)] \
                        + np.linalg.norm(self.X[node] - self.X[node1])  
                min_inter_cluster_dist = \
                        min(min_inter_cluster_dist, inter_cluster_dist)
            
            # Take maximum of intra- und inter-cluster distance          
            self.max_dist_of_node_in_cluster[parent_node][node] = \
                    max(self.max_dist_of_node_in_cluster[child_node1][node], 
                        min_inter_cluster_dist)
        for node in self.max_dist_of_node_in_cluster[child_node2].keys():
            min_inter_cluster_dist = np.inf
            for (node1, node2) in self.cluster_connections[(child_node1, 
                                                            child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if node1 not in self.max_dist_of_node_in_cluster[child_node1]:
                    node1, node2 = node2, node1
                    
                 # TODO: graph distance instead of euclidean
                inter_cluster_dist = \
                    self.max_dist_of_node_in_cluster[child_node1][node1] \
                        + self.distance_dict[(node1, node2)] \
                        + np.linalg.norm(self.X[node] - self.X[node2])
                min_inter_cluster_dist = \
                        min(min_inter_cluster_dist, inter_cluster_dist)
                      
            # Take maximum of intra- und inter-cluster distance   
            self.max_dist_of_node_in_cluster[parent_node][node] = \
                    max(self.max_dist_of_node_in_cluster[child_node2][node],
                        min_inter_cluster_dist)
        
        # Cleaning up
        self.max_dist_of_node_in_cluster.pop(child_node1)
        self.max_dist_of_node_in_cluster.pop(child_node2)
        
        # Determine all other clusters that are connected to one of the child
        # clusters. These cluster will also be connected to the parent cluster
        coord_col = []
        visited = np.empty(parent_node + 1, dtype=bool)
        visited[:] = False
        visited[parent_node] = True
        for l in self.get_nodes_connected_to_clusters([child_node1, 
                                                       child_node2]):
            while parent[l] != l:
                l = parent[l]
            if not visited[l]:
                visited[l] = True
                coord_col.append(l)
                self.A[l].append(parent_node)
        self.A.append(coord_col)
        
        # Determine for all connected clusters the distance to the newly formed
        # cluster
        for l in coord_col:
            self.cluster_connections[(parent_node, l)] = \
                        set.union(self.cluster_connections[(child_node1, l)],
                                  self.cluster_connections[(child_node2, l)])
            self.cluster_connections[(l, parent_node)] = \
                            self.cluster_connections[(parent_node, l)]
            # Find the distance between pair of nodes that are most distant 
            # in clusters rooted at l and parent_node
            inter_cluster_dist = np.inf
            for (node1, node2) in self.cluster_connections[(l, parent_node)]:
                # Canonical ordering (connecting node 1 in l)
                if node1 not in self.max_dist_of_node_in_cluster[l]:
                    node1, node2 = node2, node1
                dist = self.max_dist_of_node_in_cluster[l][node1] \
                        + self.distance_dict[(node1, node2)] \
                        + self.max_dist_of_node_in_cluster[parent_node][node2]
                inter_cluster_dist = min(inter_cluster_dist, dist)
            
            self.distance_dict[(l, parent_node)] = inter_cluster_dist
            self.distance_dict[(parent_node, l)] = inter_cluster_dist
            heappush(self.distances, (inter_cluster_dist, parent_node, l))
            
            # Cleaning up
            self.cluster_connections.pop((child_node1, l), None)
            self.cluster_connections.pop((child_node2, l), None)
            self.cluster_connections.pop((l, child_node1), None)
            self.cluster_connections.pop((l, child_node2), None)
    
    def compute_distance(self, i, j):
        return self.distance_dict[(i, j)]
    
    def get_nontrivial_clusters(self):
        return [cluster_root 
                    for cluster_root in self.max_dist_of_node_in_cluster
                        if cluster_root > self.X.shape[0]]
