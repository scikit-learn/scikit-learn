""" Linkage criteria for hierarchical clustering

Different linkage criteria that can be used by the hierarchical clustering
methods in hierarchical.py. Currently implemented are
  * Ward's linkage criterion
  * Complete linkage (farthest neighbor) linkage criterion

Authors : Jan Hendrik Metzen
License: BSD 3 clause

"""

import itertools
from heapq import heapify, heappop, heappush

import numpy as np
from scipy.spatial.distance import pdist, squareform

from . import _inertia


class Linkage(object):
    """ Base-class for different linkage criterion for hierarchical clustering.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    Methods
    -------
    update:
        Merge two clusters into new cluster and update internal structures.

    get_cluster_distance:
        Returns the distance of two  clusters.

    fetch_candidate:
        Returns the next merge candidate (the one with minimal distance).
        
    get_nodes_connected_to_clusters:
        Returns all clusters that are connected to the given clusters.
        
    get_nodes_of_cluster:
        Returns all nodes belonging to cluster with index *cluster_root*.

    Attributes
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered
    
    A : list of list of int
        The i-th entry in A contains the indices of all clusters that are
        connected to the cluster with index i
    
    distances : list (heapified) of tuples (float, int, int)
       The min-heap of merge candidates where each merge candidate consists
       of the associated distance and the indices of the two clusters to be 
       merged  
    
    cluster_nodes : dict (int to list of int)
       Mapping cluster index to indices of all data points
       that belong to the cluster
    """
    
    def __init__(self, X):
        self.X = X
        
        self.A = []
        
        # Heap of possible cluster merges, sorted by their distances
        self.distances = []
        
        # Dictionary mapping node in tree to indices of all data points
        # that belong to the cluster rooted at this node
        self.cluster_nodes = dict(zip(xrange(X.shape[0]), 
                                      [[i] for i in xrange(X.shape[0])])) 
        
    def update(self, child_node1, child_node2, parent_node):
        """ Merge two clusters into new cluster and update internal structures.
        
        This method should be called when two cluster with indices 
        *child_node1* and *child_node2* are merged into a new cluster with
        index *parent_node*. The method determines which other clusters
        are connected to parent_node and compute their distances. This
        distance is computed according to the complete-linkage criterion.
        
        Parameters
        ----------
        child_node1 : int
            Index of the first child cluster that is to be merged
        
        child_node2 : int
            Index of the second child cluster that is to be merged
        
        parent_node : int
            Index of the newly formed cluster after merging
        """
        self.cluster_nodes[parent_node] = \
            self.cluster_nodes[child_node1] + self.cluster_nodes[child_node2]
        
    def has_more_candidates(self):
        """ Returns whether there are possible merges of clusters left."""
        return len(self.distances) > 0
    
    def fetch_candidate(self):
        """ Returns the next merge candidate (the one with minimal distance).
        
        Return
        ------
        merge_candidate : tuple (float, int, int)
            A tuple specifying the next merge candidate. The tuple consists
            of the distance of the two clusters to be merged and the indices
            if these two clusters.
        """
        return heappop(self.distances)
    
    def get_nodes_connected_to_clusters(self, cluster_indices):
        """ Returns all clusters that are connected to the given clusters. 
        
        Parameters
        ----------
        cluster_indices : iterable of ints
            The indices of the clusters for which all connected clusters
            are determined 
        
        Return
        ------
        nodes : set
            The set of all nodes (i.e. cluster indices) that are connected 
            to the clusters with the given indices.
        """
        return set.union(*[set(self.A[cluster_index]) 
                                    for cluster_index in cluster_indices])
        
    def get_nodes_of_cluster(self, cluster_root):
        """ Returns all nodes belonging to cluster with index *cluster_root*.
        
        Returns all nodes (i.e. indices of datapoints in X) belonging to 
        the cluster with index *cluster_root*.
        
        Parameters
        ----------
        cluster_root : int
            The index of the cluster for which all datapoints belonging to it
            are determined
        
        Return
        ------
        nodes : list of int
            A list of datapoint indices of the cluster.
        """
        return self.cluster_nodes[cluster_root]
    
    def get_nontrivial_clusters(self):
        """ Returns the root nodes of all non-trivial clusters.
        
        Return
        ------
        _ : list of int
            A list of tree indices that correspond to non-trivial clusters, i.e.
            clusters that contain more than one datapoint.        
        """
        return [root for root, nodes in self.cluster_nodes.iteritems()
                            if len(nodes) > 1]
        

class WardsLinkage(Linkage):
    """ Ward's linkage criterion for hierarchical clustering.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered
        
    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neighboring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.

    Methods
    -------
    update:
        Merge two clusters into new cluster and update internal structures.
        
    get_cluster_distance:
        Returns the distance of two  clusters.
        
    fetch_candidate:
        Returns the next merge candidate (the one with minimal distance).
        
    get_nodes_connected_to_clusters:
        Returns all clusters that are connected to the given clusters. 
        
    get_nodes_of_cluster:
        Returns all nodes belonging to cluster with index *cluster_root*.

    Attributes
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered
    
    A : list of list of int
        The i-th entry in A contains the indices of all clusters that are 
        connected to the cluster with index i
    
    distances : list (heapified) of tuples (float, int, int)
       The min-heap of merge candidates where each merge candidate consists
       of the associated distance and the indices of the two clusters to be 
       merged  
    
    cluster_nodes : dict (int to list of int)
       Mapping cluster index to indices of all data points
       that belong to the cluster
       
    moments : list of arrays
       A list whose first entry is the array containing as i-the entry the
       number of datapoints contained in cluster with index i. The second entry
       contains as i-th entry the sum of all datapoints contained in cluster 
       with index i.       
    """
    
    def __init__(self, X, connectivity):
        super(WardsLinkage, self).__init__(X)
        
        n_samples, n_features = X.shape
        n_nodes = n_samples * 2 - 1
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
        """ Merge two clusters into new cluster and update internal structures.
        
        This method should be called when two cluster with indices 
        *child_node1* and *child_node2* are merged into a new cluster with
        index *parent_node*. The method determines which other clusters
        are connected to parent_node and compute their distances. This
        distance is computed according to the Wards criterion.
        
        Parameters
        ----------
        child_node1 : int
            Index of the first child cluster that is to be merged
        
        child_node2 : int
            Index of the second child cluster that is to be merged
        
        parent_node : int
            Index of the newly formed cluster after merging
            
        parent: array with 2 * n_samples - 1 entries
            A mapping from a cluster's index to the index of its parent
            cluster (i.e. the cluster it has been merged into). Unmerged
            clusters have themselves as parent.
        """
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
                # cluster l is now connected to parent_node instead of 
                # the two child nodes
                self.A[l].append(parent_node)
                if child_node1 in self.A[l]:
                    self.A[l].remove(child_node1)
                if child_node2 in self.A[l]:
                    self.A[l].remove(child_node2)
        self.A.append(coord_col)
        coord_col = np.array(coord_col, dtype=np.int)
        coord_row = np.empty_like(coord_col)
        coord_row.fill(parent_node)
        
        distances = np.empty(len(coord_row), dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   coord_row, coord_col, distances)
        
        for tupl in itertools.izip(distances, coord_row, coord_col):
            heappush(self.distances, tupl)
    
    def get_cluster_distance(self, i, j):
        """ Returns the distance of the clusters with indices *i* and *j*."""
        # Compute inertia of the cluster obtained when merging i and j
        distances = np.empty(1, dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   np.array([i]), np.array([j]), distances)
        return distances[0]


class CompleteLinkage(Linkage):
    """ Complete Linkage criterion for hierarchical clustering

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered
        
    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neighboring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.
            
    base_distance : function or str
        A metric string that is understood by scipy.spatial.distance.pdist
        or a binary function that returns for two datapoints their distance.
        Defaults to "euclidean".
        
        Note: it is typically much faster to pass a string than a function
              when using precompute_distances.
              
    distance_matrix: array of shape (n_samples, n_samples) or None
        A matrix which stores at index (i,j) the distance between X[i,:] and
        X[j,:]. If this matrix is given, the passed base_distance is ignored.
        Defaults to None.
              
    precompute_distances : bool
        If True, the distance between two datapoints are computed all initially
        and then stored in distance_matrix. This is typically faster then
        computing distances only when required but may require memory that 
        grows quadratically with the number of samples. Defaults to True.
        
        Note: this parameter is ignored if a non-None distance_matrix is 
              passed.

    Methods
    -------
    update:
        Merge two clusters into new cluster and update internal structures.
        
    get_cluster_distance:
        Returns the distance of two  clusters.
        
    fetch_candidate:
        Returns the next merge candidate (the one with minimal distance).
        
    get_nodes_connected_to_clusters:
        Returns all clusters that are connected to the given clusters. 
        
    get_nodes_of_cluster:
        Returns all nodes belonging to cluster with index *cluster_root*.

    Attributes
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered
    
    A : list of list of int
        The i-th entry in A contains the indices of all clusters that are 
        connected to the cluster with index i
    
    distances : list (heapified) of tuples (float, int, int)
       The min-heap of merge candidates where each merge candidate consists
       of the associated distance and the indices of the two clusters to be 
       merged  
    
    cluster_nodes : dict (int to list of int)
       Mapping cluster index to indices of all data points
       that belong to the cluster
       
    distance_matrix: array of shape (n_samples, n_samples) or None
        A matrix which stores at index (i,j) the distance between X[i,:] and
        X[j,:]. 
               
    distance_dict: dict of (int, int) -> float
       Mapping from indices of two clusters to their distance
    """

    def __init__(self, X, connectivity, base_distance="euclidean",
                 distance_matrix=None, precompute_distances=True,
                 *args):
        super(CompleteLinkage, self).__init__(X)
        
        # Set up distance matrix or base distance function
        if distance_matrix is not None:
            assert X.shape[0] == distance_matrix.shape[0] \
                    and X.shape[0] == distance_matrix.shape[1], \
                "distance_matrix must be a square matrix of size "\
                "n_samples x n_samples."
            self.distance_matrix = distance_matrix
        elif precompute_distances:
            self.distance_matrix = squareform(pdist(X, base_distance))
        else:
            self.distance_matrix = None
            if isinstance(base_distance, basestring):               
                self.base_distance = lambda i, j: pdist(X[[i, j], :],
                                                        base_distance)
            else:
                assert isinstance(base_distance, type(lambda : None)), \
                   "base_distance must be either a string that is understood "\
                   "by pdist or a binary function that returns the distance "\
                   "of two points."
                    
                self.base_distance = lambda i, j: base_distance(self.X[i, :],
                                                                self.X[j, :])
            
        # Distances between two clusters, represented by their 
        # root node indices 
        self.distance_dict = {}  
        
        # Determine distances between all connected nodes and create matrix A  
        for ind1, row in enumerate(connectivity.rows):
            self.A.append(row)
            for ind2 in row:
                if ind1 == ind2:
                    continue
                dist = self.get_cluster_distance(ind1, ind2)
                self.distances.append((dist, ind1, ind2))
                
        # Enforce symmetry of A
        for ind1, row in enumerate(connectivity.rows):
            for ind2 in row:
                if ind1 not in self.A[ind2]:
                    self.A[ind2].append(ind1)
        
        heapify(self.distances)

    def update(self, child_node1, child_node2, parent_node, parent):
        """ Merge two clusters into new cluster and update internal structures.
        
        This method should be called when two cluster with indices 
        *child_node1* and *child_node2* are merged into a new cluster with
        index *parent_node*. The method determines which other clusters
        are connected to parent_node and compute their distances. This
        distance is computed according to the complete-linkage criterion.
        
        Parameters
        ----------
        child_node1 : int
            Index of the first child cluster that is to be merged
        
        child_node2 : int
            Index of the second child cluster that is to be merged
        
        parent_node : int
            Index of the newly formed cluster after merging
            
        parent: array with 2 * n_samples - 1 entries
            A mapping from a cluster's index to the index of its parent
            cluster (i.e. the cluster it has been merged into). Unmerged
            clusters have themselves as parent.
        """
        super(CompleteLinkage, self).update(child_node1, child_node2, 
                                            parent_node)       
        # Determine all other clusters that are connected to one of the child
        # clusters. These cluster will also be connected to the parent cluster
        coord_col = []
        visited = np.empty(parent_node + 1, dtype=bool)
        visited[:] = False
        visited[parent_node] = True
        for l in self.get_nodes_connected_to_clusters([child_node1, 
                                                       child_node2]):
            if l in [child_node1, child_node2]: 
                continue
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
            if connected_cluster in [child_node1, child_node2]:
                continue           
            # The distance of the connected cluster to the the newly formed 
            # parent cluster is the maximum of the distances of the connected
            # cluster to the two child clusters
            max_dist = max(self.get_cluster_distance(connected_cluster, 
                                                     child_node1),
                           self.get_cluster_distance(connected_cluster, 
                                                     child_node2))
            
            # Update internal structures with the distance of the connected
            # cluster to the parent cluster
            self.distance_dict[(connected_cluster, parent_node)] = max_dist
            self.distance_dict[(parent_node, connected_cluster)] = max_dist
            heappush(self.distances, 
                     (max_dist, parent_node, connected_cluster))
            
            # Clean up
            self.distance_dict.pop((connected_cluster, child_node1), None)
            self.distance_dict.pop((connected_cluster, child_node2), None)
            self.distance_dict.pop((child_node1, connected_cluster), None)
            self.distance_dict.pop((child_node2, connected_cluster), None)
                
    def get_cluster_distance(self, i, j):
        """ Returns the distance of the clusters with indices *i* and *j*."""
        if i < self.X.shape[0] and j < self.X.shape[0]:
            # Trivial clusters consisting each of 1 datapoint -> return 
            # distance of datapoints
            if self.distance_matrix is not None:
                return self.distance_matrix[i, j]
            else:
                return self.base_distance(i, j)
        if not (i, j) in self.distance_dict:
            # Distance of clusters not known -> compute it
            if self.distance_matrix is not None:
                dist = \
                  np.max(self.distance_matrix[self.get_nodes_of_cluster(i), :]
                                             [:, self.get_nodes_of_cluster(j)])
            else:
                dist = max(self.base_distance(k, l)
                                for k in self.get_nodes_of_cluster(i)
                                    for l in self.get_nodes_of_cluster(j))
            self.distance_dict[(i, j)] = dist 
            self.distance_dict[(j, i)] = dist
            
        return self.distance_dict[(i, j)]
