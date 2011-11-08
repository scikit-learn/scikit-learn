"""Hierarchical Agglomerative Clustering

These routines perform some hierarchical agglomerative clustering of some
input data. Currently, only complete-linkage and ward's criterion are
implemented as linkage criteria.

Authors : Vincent Michel, Bertrand Thirion, Alexandre Gramfort,
          Gael Varoquaux, Jan Hendrik Metzen
License: BSD 3 clause
"""

from collections import defaultdict
from heapq import heapify, heappop, heappush
import itertools
import warnings

import numpy as np
from scipy import sparse
from scipy.cluster import hierarchy

from ..base import BaseEstimator
from ..utils._csgraph import cs_graph_components
from ..externals.joblib import Memory

from . import _inertia
from ._feature_agglomeration import AgglomerationTransform


###############################################################################

class WardsLinkage(object):

    def __init__(self, X, connectivity, n_nodes, n_features, n_samples):
        # build moments as a list
        self.moments = [np.zeros(n_nodes), np.zeros((n_nodes, n_features))]
        self.moments[0][:n_samples] = 1
        self.moments[1][:n_samples] = X
    
        self.A = []
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
        visited = np.empty(len(self.A), dtype=bool)
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
            
    def has_more_candidates(self):
        return len(self.distances) > 0
    
    def fetch_candidate(self):
        return heappop(self.distances)
    
    def compute_distance(self, i, j):
        # Compute inertia of the cluster obtained when merging i and j
        distances = np.empty(1, dtype=np.float)
        _inertia.compute_ward_dist(self.moments[0], self.moments[1],
                                   np.array([i]), np.array([j]), distances)
        return distances[0]


class CompleteLinkage(object):
    
    def __init__(self, X, connectivity, n_nodes, n_features, n_samples):
        self.X = X
        self.A = [] 
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
            minInterClusterDist = np.inf
            for (node1, node2) in self.cluster_connections[(child_node1, 
                                                            child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if node1 not in self.max_dist_of_node_in_cluster[child_node1]:
                    node1, node2 = node2, node1
                    
                # TODO: graph distance instead of euclidean
                interClusterDist = \
                    self.max_dist_of_node_in_cluster[child_node2][node2] \
                        + self.distance_dict[(node1, node2)] \
                        + np.linalg.norm(self.X[node] - self.X[node1])  
                minInterClusterDist = \
                        min(minInterClusterDist, interClusterDist)
            
            # Take maximum of intra- und inter-cluster distance          
            self.max_dist_of_node_in_cluster[parent_node][node] = \
                    max(self.max_dist_of_node_in_cluster[child_node1][node], 
                        minInterClusterDist)
        for node in self.max_dist_of_node_in_cluster[child_node2].keys():
            minInterClusterDist = np.inf
            for (node1, node2) in self.cluster_connections[(child_node1, 
                                                            child_node2)]:
                # Canonical ordering (connecting node 1 in cluster 1)
                if node1 not in self.max_dist_of_node_in_cluster[child_node1]:
                    node1, node2 = node2, node1
                    
                 # TODO: graph distance instead of euclidean
                interClusterDist = \
                    self.max_dist_of_node_in_cluster[child_node1][node1] \
                        + self.distance_dict[(node1, node2)] \
                        + np.linalg.norm(self.X[node] - self.X[node2])
                minInterClusterDist = \
                        min(minInterClusterDist, interClusterDist)
                      
            # Take maximum of intra- und inter-cluster distance   
            self.max_dist_of_node_in_cluster[parent_node][node] = \
                    max(self.max_dist_of_node_in_cluster[child_node2][node],
                        minInterClusterDist)
        
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
            # Find the distance between pair of nodeWs that are most distant 
            # in clusters rooted at l and parent_node
            interClusterDist = np.inf
            for (node1, node2) in self.cluster_connections[(l, parent_node)]:
                # Canonical ordering (connecting node 1 in l)
                if node1 not in self.max_dist_of_node_in_cluster[l]:
                    node1, node2 = node2, node1
                dist = self.max_dist_of_node_in_cluster[l][node1] \
                        + self.distance_dict[(node1, node2)] \
                        + self.max_dist_of_node_in_cluster[parent_node][node2]
                interClusterDist = min(interClusterDist, dist)
            
            self.distance_dict[(l, parent_node)] = interClusterDist
            self.distance_dict[(parent_node, l)] = interClusterDist
            heappush(self.distances, (interClusterDist, parent_node, l))
            
            # Cleaning up
            self.cluster_connections.pop((child_node1, l), None)
            self.cluster_connections.pop((child_node2, l), None)
            self.cluster_connections.pop((l, child_node1), None)
            self.cluster_connections.pop((l, child_node2), None)
            
    def has_more_candidates(self):
        return len(self.distances) > 0
    
    def fetch_candidate(self):
        return heappop(self.distances)
    
    def compute_distance(self, i, j):
        return self.distance_dict[(i, j)]
    
    def get_nodes_connected_to_clusters(self, cluster_roots):
        return set.union(*[set(self.A[cluster_root]) 
                                for cluster_root in cluster_roots])
    
    def get_nontrivial_clusters(self):
        return [cluster_root 
                    for cluster_root in self.max_dist_of_node_in_cluster
                        if cluster_root > self.X.shape[0]]
    
            
def dendrogram(X, connectivity=None, n_components=None, return_inertias=False, 
              merge_replay=[], linkage_criteria='CompleteLinkage', 
              copy=True):
    """Hierarchical clustering based on a Feature matrix.

    The inertia matrix uses a Heapq-based representation.

    This is the structured version, that takes into account a some topological
    structure between samples.

    Parameters
    ----------
    X : array of shape (n_samples, n_features)
        feature matrix  representing n_samples samples to be clustered

    connectivity : sparse matrix.
        connectivity matrix. Defines for each sample the neigbhoring samples
        following a given structure of the data. The matrix is assumed to
        be symmetric and only the upper triangular half is used.
        Default is None, i.e, the Ward algorithm is unstructured.

    n_components : int (optional)
        Number of connected components. If None the number of connected
        components is estimated from the connectivity matrix.
        
    return_inertias : bool (optional)
        If true, the inertias are returned additionally

    copy : bool (optional)
        Make a copy of connectivity or work inplace. If connectivity
        is not of LIL type there will be a copy in any case.

    Returns
    -------
    children : list of pairs. Length of n_nodes
               list of the children of each nodes.
               Leaves of the tree have empty list of children.

    n_components : sparse matrix.
        The number of connected components in the graph.

    n_leaves : int
        The number of leaves in the tree
        
    heights :  list of floats. Length n_nodes
        The inertias associated to the tree's nodes. Only returned, if 
        return_inertias is True
    """
    X = np.asarray(X)
    n_samples, n_features = X.shape
    if X.ndim == 1:
        X = np.reshape(X, (-1, 1))

    try:
        linkage_criteria = eval(linkage_criteria)
    except:
        raise Exception("Unknown distance class %s" % linkage_criteria)

    # Compute the number of nodes
    if connectivity is not None:
        if n_components is None:
            n_components, _ = cs_graph_components(connectivity)
        if n_components > 1:
            warnings.warn("the number of connected components of the"
            " connectivity matrix is %d > 1. The tree will be stopped early."
            % n_components)
    else:
        out = hierarchy.ward(X)
        children_ = out[:, :2].astype(np.int)
        return children_, 1, n_samples, None

    n_nodes = 2 * n_samples - 1

    if (connectivity.shape[0] != n_samples or
        connectivity.shape[1] != n_samples):
        raise ValueError('Wrong shape for connectivity matrix: %s '
                         'when X is %s' % (connectivity.shape, X.shape))
    # convert connectivity matrix to LIL eventually with a copy
    if sparse.isspmatrix_lil(connectivity) and copy:
        connectivity = connectivity.copy()
    else:
        connectivity = connectivity.tolil()

    # Remove diagonal from connectivity matrix
    connectivity.setdiag(np.zeros(connectivity.shape[0]))        
    # Compute distances between connected nodes
    linkage = linkage_criteria(X, connectivity, n_nodes, n_features, 
                                      n_samples)

    # prepare the main fields
    parent = np.arange(n_nodes, dtype=np.int)
    heights = np.zeros(n_nodes)
    open_nodes = np.ones(n_nodes, dtype=bool)
    children = []
    merges = []
    
    # The indices which are contained in the replayed merges
    reserved_indices = set([])
    for (i_, j_), k_ in merge_replay:
        reserved_indices.add(i_)
        reserved_indices.add(j_)
    for (i_, j_), k_ in merge_replay:
        if k_ in reserved_indices:
             # Add k_ later since we don't know the new index yet
            reserved_indices.remove(k_)
    
    index_mapping = dict(zip(reserved_indices, reserved_indices))
    
    augmentation = None
    # recursive merge loop
    replaying = len(merge_replay) > 0
    for k in xrange(n_samples, n_nodes):
        # Fetch merge that will be reapplied next (if any)
        if len(merge_replay) > 0:
            (i_, j_), k_ = merge_replay[0]
            merge_replay = merge_replay[1:]
             # Associate indices
            i = index_mapping[i_]
            j = index_mapping[j_]
            # Compute merge distance
            merge_distance = linkage.compute_distance(i, j)
            index_mapping[k_] = k
#            print "Reapply", merge_distance, i, j, k
        else:  # No merges to be reapplied left
            if augmentation: 
                # There are nodes that can be added to cluster created during
                # replaying 
                merge_distance, i, j = augmentation
                augmentation = None
#                print "Augmentation", merge_distance, i, j, k
            else:
                # Identify the merge that will be applied next using the 
                # standard method
                merge_distance = np.inf
                while linkage.has_more_candidates():
                    merge_distance, i, j = linkage.fetch_candidate()
                    if open_nodes[i] and open_nodes[j]:
                        break
                    if not linkage.has_more_candidates():
                        merge_distance = np.inf
                        break
                if not linkage.has_more_candidates() \
                                        and merge_distance == np.inf: 
                    # Merge unconnected components with height infinity
                    i = k - 1
                    desc = set(_hc_get_descendent([i], np.array(children), 
                                                  n_samples, 
                                                  add_intermediate_nodes=True))
                    for j in xrange(k - 2, 0, -1):
                        if j not in desc:
                            break      
#                print "Novel", merge_distance, i, j, k     
            # Update index mapping
            if i in index_mapping.values():
                i_ = index_mapping.keys()[index_mapping.values().index(i)]
                index_mapping[i_] = k
            if j in index_mapping.values():
                j_ = index_mapping.keys()[index_mapping.values().index(j)]
                index_mapping[j_] = k  

        parent[i], parent[j], heights[k] = k, k, merge_distance
        merges.append(((i, j), k))
        children.append([i, j])
        open_nodes[i], open_nodes[j] = False, False
        
        # Add new possible merges and their distances to heap
        linkage.update(i, j, k, parent)
        
        # If we are finished with replaying
        if replaying and len(merge_replay) == 0:
            # Check for possible augmentations of the clusters created during 
            # replaying
            possible_augmentations = []
            for cluster_root in linkage.get_nontrivial_clusters():
                for node in linkage.get_nodes_connected_to_clusters(
                                                            [cluster_root]):
                    if not open_nodes[node]:
                        continue
                    dist = linkage.compute_distance(node, cluster_root)
                    if dist < heights[cluster_root]:
                        possible_augmentations.append(
                                                (dist - heights[cluster_root],
                                                 heights[cluster_root],
                                                 node, cluster_root))
            if len(possible_augmentations) > 0:
                augmentation = min(possible_augmentations)[1:]
            else:
                replaying = False

    # Separate leaves in children (empty lists up to now)
    n_leaves = n_samples
    children = np.array(children)  # return numpy array for efficient caching

    if return_inertias:
        return children, n_components, n_leaves, merges, heights
    else:
        return children, n_components, n_leaves, merges


###############################################################################
# Functions for cutting  hierarchical clustering tree

def _hc_get_descendent(ind, children, n_leaves, add_intermediate_nodes=False):
    """Function returning all the descendent leaves of a set of nodes.

    Parameters
    ----------
    ind : list of int
        A list that indicates the nodes for which we want the descendents.

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        This is not defined for leaves.

    n_leaves : int
        Number of leaves.

    Return
    ------
    descendent : list of int
    """
    descendent = []
    while len(ind) != 0:
        i = ind.pop()
        if i < n_leaves:
            descendent.append(i)
        else:
            if add_intermediate_nodes:
                descendent.append(i)
            ci = children[i - n_leaves]
            ind.extend((ci[0], ci[1]))
    return descendent


def _hc_cut(n_clusters, children, n_leaves):
    """Function cutting the ward tree for a given number of clusters.

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters to form.

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        Leaves have empty list of children and are not stored.

    n_leaves : int
        Number of leaves of the tree.

    Return
    ------
    labels : array [n_points]
        cluster labels for each point

    """
    nodes = [np.max(children[-1]) + 1]
    for i in xrange(n_clusters - 1):
        nodes.extend(children[np.max(nodes) - n_leaves])
        nodes.remove(np.max(nodes))
    labels = np.zeros(n_leaves, dtype=np.int)
    for i, node in enumerate(nodes):
        labels[_hc_get_descendent([node], children, n_leaves)] = i
    return labels


def _hc_cut_inertia(max_inertia, children, n_leaves, inertias):
    """Function cutting the ward tree for a given maximal inertia.

    Parameters
    ----------
    max_inertia : float
        The maximal inertia a cluster is allowed to have. Determines indirectly
        how many clusters are formed. 

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        Leaves have empty list of children and are not stored.

    n_leaves : int
        Number of leaves of the tree.
        
    inertias :  list of floats. Length of n_nodes
        The inertias associated to the tree's nodes.

    Return
    ------
    labels : array
        cluster labels for each point

    """
    open_nodes = [len(inertias) - 1]  # root of the tree
    cluster_roots = [] 
    while open_nodes != []:
        node = open_nodes[0]
        open_nodes = open_nodes[1:]
        if inertias[node] <= max_inertia: 
            # This tree node is the root of a cluster
            cluster_roots.append(node)
        else:
            # Tree node induces subtree with too large inertia; split it
            open_nodes.extend(children[node - n_leaves])
            
    labels = np.zeros(n_leaves, dtype=np.int)
    for i, node in enumerate(cluster_roots):
        labels[_hc_get_descendent([node], children, n_leaves)] = i
    return labels


###############################################################################
# Class for Hierarchical clustering

class HierarchicalClustering(BaseEstimator):
    """Hierarchical clustering: constructs a tree and cuts it.

    Parameters
    ----------
    n_clusters : int or ndarray or None
        The number of clusters to find.
        
    max_inertia : float or None
        If this value is not None, the max_inertia is used to determine the
        number of returned clusters. In this case, the n_clusters parameter
        is agnored and may be None.

    connectivity : sparse matrix.
        Connectivity matrix. Defines for each sample the neigbhoring
        samples following a given structure of the data.
        Default is None, i.e, the hiearchical clustering algorithm is
        unstructured.

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    copy : bool
        Copy the connectivity matrix or work inplace.

    n_components : int (optional)
        The number of connected components in the graph defined by the
        connectivity matrix. If not set, it is estimated.

    Methods
    -------
    fit:
        Compute the clustering

    Attributes
    ----------
    children_ : array-like, shape = [n_nodes, 2]
        List of the children of each nodes.
        Leaves of the tree do not appear.

    labels_ : array [n_points]
        cluster labels for each point

    n_leaves_ : int
        Number of leaves in the hiearchical tree.

    """

    def __init__(self, n_clusters=2, max_inertia=None,
                 memory=Memory(cachedir=None, verbose=0), connectivity=None,
                 copy=True, n_components=None):
        self.n_clusters = n_clusters
        self.max_inertia = max_inertia
        self.memory = memory
        self.copy = copy
        self.n_components = n_components
        self.connectivity = connectivity
        self.merges = []

    def fit(self, X):
        """Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The samples a.k.a. observations.

        Returns
        -------
        self
        """
        memory = self.memory
        if isinstance(memory, basestring):
            memory = Memory(cachedir=memory)

        # Construct the tree
        if self.max_inertia is None:
            children, n_components, n_leaves, self.merges = \
                    memory.cache(dendrogram)(X, self.connectivity,
                                            n_components=self.n_components,
                                            merge_replay=self.merges,
                                            copy=self.copy)
            # Cut the tree based on number of desired clusters
            self.labels_ = _hc_cut(self.n_clusters, children, n_leaves)
        else:
            children, n_components, n_leaves, self.merges, inertias = \
                    memory.cache(dendrogram)(X, self.connectivity,
                                            n_components=self.n_components,
                                            return_inertias=True, 
                                            merge_replay=self.merges,
                                            copy=self.copy)
            # Cut the tree based on maximally allowed inertia
            self.labels_ = _hc_cut_inertia(self.max_inertia, children, 
                                           n_leaves, inertias)
            
        # Undo merges in the uppermost part of the tree that have been 
        # "cut off" by _hc_cut_inertia or _hc_cut
        if len(np.unique(self.labels_)) > 1:
            self.merges = self.merges[:-len(np.unique(self.labels_)) + 1]
            
        return self


###############################################################################
class Ward(HierarchicalClustering):
    pass


# Ward-based feature agglomeration
class WardAgglomeration(AgglomerationTransform, Ward):
    """Feature agglomeration based on Ward hierarchical clustering

    Parameters
    ----------
    n_clusters : int or ndarray
        The number of clusters.

    connectivity : sparse matrix
        connectivity matrix. Defines for each feature the neigbhoring
        features following a given structure of the data.
        Default is None, i.e, the hiearchical agglomeration algorithm is
        unstructured.

    memory : Instance of joblib.Memory or string
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    copy : bool
        Copy the connectivity matrix or work inplace.

    n_components : int (optional)
        The number of connected components in the graph defined by the
        connectivity matrix. If not set, it is estimated.

    Methods
    -------
    fit:
        Compute the clustering of features

    Attributes
    ----------
    children_ : array-like, shape = [n_nodes, 2]
        List of the children of each nodes.
        Leaves of the tree do not appear.

    labels_ : array [n_points]
        cluster labels for each point

    n_leaves_ : int
        Number of leaves in the hiearchical tree.

    """

    def fit(self, X, y=None, **params):
        """Fit the hierarchical clustering on the data

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The data

        Returns
        -------
        self
        """
        return Ward.fit(self, X.T, **params)
