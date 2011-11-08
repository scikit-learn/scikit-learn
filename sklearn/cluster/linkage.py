
# TODO: Consider nearest-neighbor chain algorithm ?
# https://secure.wikimedia.org/wikipedia/en/wiki/Nearest-neighbor_chain_algorithm

import itertools
from heapq import heapify, heappop, heappush

import numpy as np

from . import _inertia

class Linkage(object):
    
    def __init__(self, X):
        self.X = X
        
        self.A = []
        
        # Heap of possible cluster merges, sorted by their distances
        self.distances = []
        
        # Dictionary mapping node in tree to indices of all data points
        # which belong to the cluster rooted at this node
        self.cluster_nodes = dict(zip(xrange(X.shape[0]), 
                                      [[i] for i in xrange(X.shape[0])])) 
        
    def update(self, child_node1, child_node2, parent_node):
        self.cluster_nodes[parent_node] = \
            self.cluster_nodes[child_node1] + self.cluster_nodes[child_node2]
        
    def has_more_candidates(self):
        return len(self.distances) > 0
    
    def fetch_candidate(self):
        return heappop(self.distances)
    
    def get_nodes_connected_to_clusters(self, cluster_roots):
        return set.union(*[set(self.A[cluster_root]) 
                                for cluster_root in cluster_roots])
        
    def get_nodes_of_cluster(self, cluster_root):
        return self.cluster_nodes[cluster_root]
        

class WardsLinkage(Linkage):

    def __init__(self, X, connectivity):
        super(WardsLinkage, self).__init__(X)
        
        n_samples, n_features = X.shape
        n_nodes = n_samples * 2 -1
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
        super(WardsLinkage, self).update(child_node1, child_node2, parent_node)
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



class CompleteLinkage(Linkage):
    
    def __init__(self, X, connectivity, 
                 base_distance=lambda x,y: np.linalg.norm(x - y), *args):
        super(CompleteLinkage, self).__init__(X)
        
        self.base_distance = base_distance
        # Distances between two clusters, represented by their 
        # root node indices 
        self.distance_dict = {}  
        # Determine distances between all connected nodes 
        for ind1, row in enumerate(connectivity.rows):
            self.A.append(row)
            for ind2 in row:
                if ind1 == ind2: continue
                # Compute distance between two connected nodes
                dist = self.base_distance(X[ind1], X[ind2])
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
        super(CompleteLinkage, self).update(child_node1, child_node2, parent_node)       
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
        for connected_cluster in coord_col:
            # Find the distance between pair of nodes that are most distant 
            # in clusters rooted at connected_cluster and parent_node.
            max_dist = 0.0
            # Do a brute force search
            for node1 in self.get_nodes_of_cluster(parent_node):
                for node2 in self.get_nodes_of_cluster(connected_cluster):
                    max_dist = max(max_dist, self.base_distance(self.X[node1],
                                                                self.X[node2]))
            
            self.distance_dict[(connected_cluster, parent_node)] = max_dist
            self.distance_dict[(parent_node, connected_cluster)] = max_dist
            heappush(self.distances, (max_dist, parent_node, connected_cluster))
                
    def compute_distance(self, i, j):
        return self.distance_dict[(i, j)]
    